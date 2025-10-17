#!/usr/bin/env python3

import numpy as np
import gzip
import json
import sys

def parse_deeptools_matrix(matrix_file):
    """Parse deepTools computeMatrix output"""
    with gzip.open(matrix_file, 'rt') as f:
        # Parse JSON metadata from first line
        metadata_line = f.readline().strip()
        metadata = json.loads(metadata_line[1:])  # Remove '@' prefix

        # Extract sample information
        samples = metadata['sample_labels']

        # Read data lines
        data_lines = []
        for line in f:
            if line.strip():
                data_lines.append(line.strip().split('\t'))

        # Extract matrix values (skip first 6 genomic coordinate columns)
        matrix_values = []
        for line in data_lines:
            values = [float(x) for x in line[6:]]  # Skip chr,start,end...
            matrix_values.append(values)

        # Convert to numpy array - shape: (n_regions, n_samples * n_bins)
        data_matrix = np.array(matrix_values)

    return samples, data_matrix

def calculate_tss_enrichment(matrix_file, output_file):
    """Calculate TSS enrichment scores"""
    samples, data = parse_deeptools_matrix(matrix_file)

    if len(samples) == 0:
        raise ValueError("No samples found in matrix file")

    if data.size == 0:
        raise ValueError("No data found in matrix file")

    n_samples = len(samples)
    n_bins = data.shape[1] // n_samples

    print(f"Processing {n_samples} samples with {n_bins} bins each")
    print(f"Sample names: {samples}")

    # Validate that data dimensions make sense
    if data.shape[1] % n_samples != 0:
        raise ValueError(f"Data columns ({data.shape[1]}) not evenly divisible by number of samples ({n_samples})")

    # Validate matrix dimensions
    if n_bins < 50:
        raise ValueError(f"Matrix has too few bins ({n_bins}). Expected at least 50 bins for TSS enrichment calculation.")

    enrichments = []

    center = n_bins // 2
    window_size = 8       # ±400 bp if 50 bp bins (broader TSS capture)
    flank_distance = 25   # 1.25kb from center to start of flanks
    flank_size = 10       # ±500 bp flanks
    max_flank = flank_distance + flank_size  # Total distance from center

    # Validate parameters against matrix size
    if center - max_flank < 0 or center + max_flank >= n_bins:
        # Adjust parameters for smaller matrices
        max_flank = min(max_flank, center, n_bins - center - 1)
        flank_distance = max(window_size + 5, max_flank - flank_size)
        flank_size = min(flank_size, max_flank - flank_distance)

    for i, sample in enumerate(samples):
        sample_data = data[:, i*n_bins:(i+1)*n_bins]

        # Define safe bounds with explicit validation
        left_start = max(center - flank_distance - flank_size, 0)
        left_end   = max(center - flank_distance, 0)
        right_start = min(center + flank_distance, n_bins)
        right_end   = min(center + flank_distance + flank_size, n_bins)

        center_start = max(center - window_size, 0)
        center_end   = min(center + window_size, n_bins)

        # Additional bounds validation
        if left_start >= left_end or right_start >= right_end or center_start >= center_end:
            print(f"Warning: Invalid region bounds for sample {sample}, using fallback calculation")
            enrichment = np.nan
        else:
            # Compute means only on valid slices
            tss_region = np.mean(sample_data[:, center_start:center_end], axis=1)
            flank_left = np.mean(sample_data[:, left_start:left_end], axis=1)
            flank_right = np.mean(sample_data[:, right_start:right_end], axis=1)
            flanking = (flank_left + flank_right) / 2

            valid = flanking > 1e-3
            enrichment = np.mean(tss_region[valid] / flanking[valid]) if np.any(valid) else np.nan


        enrichments.append((sample, enrichment))

    # Write results
    with open(output_file, 'w') as f:
        f.write('sample_id\tTSS_Enrichment\n')
        for sample, enrichment in enrichments:
            f.write(f'{sample}\t{enrichment:.3f}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_tss_enrichment.py <matrix_file> <output_file>")
        sys.exit(1)

    matrix_file = sys.argv[1]
    output_file = sys.argv[2]

    calculate_tss_enrichment(matrix_file, output_file)
