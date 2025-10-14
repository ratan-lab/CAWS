#!/usr/bin/env python3

import numpy as np
import gzip
import json
import sys

def parse_deeptools_matrix(matrix_file):
    """Parse deepTools computeMatrix output"""
    with gzip.open(matrix_file, 'rt') as f:
        lines = f.readlines()

    samples = []
    matrices = []
    data_start = False

    for line in lines:
        if line.startswith('@') and 'sample_labels' in line:
            # Extract sample names from the JSON-like format
            # Example line: @sample_labels ["sample1","sample2","sample3"]
            try:
                # Find the JSON array in the line
                start_bracket = line.find('[')
                end_bracket = line.rfind(']')
                if start_bracket != -1 and end_bracket != -1:
                    samples_str = line[start_bracket:end_bracket+1]
                    samples = json.loads(samples_str)
                else:
                    # Fallback: try to extract from quoted strings
                    import re
                    quoted_samples = re.findall(r'"([^"]+)"', line)
                    samples = quoted_samples if quoted_samples else ["unknown_sample"]
            except (json.JSONDecodeError, ValueError) as e:
                print(f"Warning: Could not parse sample labels from line: {line.strip()}")
                print(f"Error: {e}")
                # Use fallback sample names
                samples = [f"sample_{i+1}" for i in range(1)]  # Will be adjusted later based on data
        elif not line.startswith('@'):
            if not data_start:
                data_start = True
            # Skip first 6 columns (chr, start, end, name, score, strand)
            values = [float(x) for x in line.strip().split()[6:]]
            matrices.append(values)

    data_matrix = np.array(matrices)

    # If we couldn't parse sample names properly, infer from data dimensions
    if len(samples) == 1 and samples[0].startswith('sample_'):
        # Infer number of samples from data columns
        if data_matrix.size > 0:
            total_cols = data_matrix.shape[1]
            # Assume bins per sample is consistent, try common values
            for bins_per_sample in [200, 201, 160, 161, 400, 401]:
                if total_cols % bins_per_sample == 0:
                    n_samples = total_cols // bins_per_sample
                    samples = [f"sample_{i+1}" for i in range(n_samples)]
                    print(f"Inferred {n_samples} samples with {bins_per_sample} bins each")
                    break
            else:
                # Fallback: assume single sample
                samples = ["sample_1"]
                print(f"Warning: Could not infer sample count from {total_cols} columns, assuming single sample")

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
    if n_bins < 200:
        raise ValueError(f"Matrix has too few bins ({n_bins}). Expected at least 200 bins for TSS enrichment calculation.")

    enrichments = []

    center = n_bins // 2
    window_size = 4   # ±200 bp if 50 bp bins
    flank_size = 20   # 1000 bp if 50 bp bins
    max_flank = 100   # 5000 bp if 50 bp bins

    # Validate parameters against matrix size
    if center - max_flank < 0 or center + max_flank >= n_bins:
        # Adjust parameters for smaller matrices
        max_flank = min(max_flank, center, n_bins - center - 1)
        flank_size = min(flank_size, max_flank - window_size - 1)

    for i, sample in enumerate(samples):
        sample_data = data[:, i*n_bins:(i+1)*n_bins]

        # Define safe bounds with explicit validation
        left_start = max(center - max_flank, 0)
        left_end   = max(center - flank_size, 0)
        right_start = min(center + flank_size, n_bins)
        right_end   = min(center + max_flank, n_bins)

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
