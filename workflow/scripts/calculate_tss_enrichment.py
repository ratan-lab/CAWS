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

    for i, sample in enumerate(samples):
        # Extract data for this sample
        sample_data = data[:, i*n_bins:(i+1)*n_bins]

        # Calculate mean signal in central region (±200bp around TSS)
        center = n_bins // 2
        window_size = 10  # ±200bp with 20bp bins = 10 bins
        tss_region = np.mean(sample_data[:, center-window_size:center+window_size], axis=1)

        # Calculate mean signal in flanking regions (±1000-2000bp)
        flank_size = 50  # 1000bp with 20bp bins = 50 bins
        flank_left = np.mean(sample_data[:, center-100:center-flank_size], axis=1)
        flank_right = np.mean(sample_data[:, center+flank_size:center+100], axis=1)
        flanking = (flank_left + flank_right) / 2

        # Calculate enrichment score (avoid division by zero)
        valid_regions = (flanking > 0.001)
        if np.sum(valid_regions) > 0:
            enrichment = np.mean(tss_region[valid_regions] / flanking[valid_regions])
        else:
            enrichment = 0.0

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