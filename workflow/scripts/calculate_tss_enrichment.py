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
            samples_str = line.split('sample_labels')[1].strip()
            samples_str = samples_str.split("]")[0] + ']'
            samples_str = samples_str[2:]
            samples = json.loads(samples_str)
        elif not line.startswith('@'):
            if not data_start:
                data_start = True
            # Skip first 6 columns (chr, start, end, name, score, strand)
            values = [float(x) for x in line.strip().split()[6:]]
            matrices.append(values)

    return samples, np.array(matrices)

def calculate_tss_enrichment(matrix_file, output_file):
    """Calculate TSS enrichment scores"""
    samples, data = parse_deeptools_matrix(matrix_file)

    n_samples = len(samples)
    n_bins = data.shape[1] // n_samples

    enrichments = []

    center = n_bins // 2
    window_size = 4   # ±200 bp if 50 bp bins
    flank_size = 20   # 1000 bp if 50 bp bins
    max_flank = 100   # 5000 bp if 50 bp bins
    
    for i, sample in enumerate(samples):
        sample_data = data[:, i*n_bins:(i+1)*n_bins]
    
        # Define safe bounds
        left_start = max(center - max_flank, 0)
        left_end   = max(center - flank_size, 0)
        right_start = min(center + flank_size, n_bins)
        right_end   = min(center + max_flank, n_bins)
    
        center_start = max(center - window_size, 0)
        center_end   = min(center + window_size, n_bins)
    
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
        f.write('Sample\tTSS_Enrichment\n')
        for sample, enrichment in enrichments:
            f.write(f'{sample}\t{enrichment:.3f}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_tss_enrichment.py <matrix_file> <output_file>")
        sys.exit(1)

    matrix_file = sys.argv[1]
    output_file = sys.argv[2]

    calculate_tss_enrichment(matrix_file, output_file)
