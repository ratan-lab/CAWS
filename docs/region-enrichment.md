# Region enrichment

CAWS measures enrichment of every sample over region sets: the transcription start sites of protein-coding genes (region set `TSS`, when `gtf_file` is set; GENCODE `gene_type` or Ensembl `gene_biotype`) and any BED files listed in `enrichment_regions`, for example predicted G-quadruplex sites (PQS) for BG4, or ATAC-seq peaks of the same cell line to ask whether peaks of a protein lie in open chromatin.

```json
"enrichment_regions": {"PQS": "/path/to/pqs.bed", "ATAC_K562": "/path/to/K562_atac.narrowPeak"}
```

Region files can be BED, BED.gz or narrowPeak; only regions on the chromosomes in `chromosome_file` are used, and chromosome names must match the reference. Each region has a **reference point**: the summit for narrowPeak files, the site itself for single-base regions (such as TSSs, strand-aware), and the center otherwise.

**Three measures, for every sample including the IgG controls:**
1. **Signal profile** – signal ±2 kb around the reference points (deepTools `computeMatrix`, 50 bp bins), plotted as a profile and a heatmap. The profile score is the mean signal within ±400 bp of the reference point divided by the mean signal 1.25–1.75 kb away, averaged over regions. It is computed from the normalized bigWig tracks, but as a ratio within each sample it does not depend on the normalization. The profile plots and heatmaps show normalized signal, so they can be compared between samples.
2. **Reads in regions** – the percentage of filtered fragments overlapping the regions, divided by the percentage of the analyzed genome (listed chromosomes minus blacklist) the regions cover. 1 means no enrichment.
3. **Peak overlap** – the percentage of a sample's SEACR and MACS3 peaks, and of the consensus peaks, overlapping the regions, compared with the same peaks shuffled 10 times within their chromosomes outside the blacklist (fixed random seeds, so results are reproducible). The fold enrichment is observed / mean shuffled overlap.

**Compare with IgG.** pA-Tn5 inserts preferentially into accessible chromatin, especially for low-abundance targets and IgG, so most CUT&Tag samples appear enriched at open-chromatin regions (ATAC peaks, promoters, many PQS). Enrichment of a sample by itself therefore does not separate binding from accessibility bias; the `*_vs_igg` columns divide each measure by that of the sample's matched IgG control (the control of its `group`). Values clearly above 1 indicate enrichment beyond the accessibility background.

**Outputs:**
- `stats/region_enrichment.txt` — one row per region set and sample (and per consensus peak set): profile score, reads-in-regions percentage and fold, peak-overlap percentage and fold for SEACR and MACS3, and the IgG ratios; also shown on the report's Region Enrichment page
- `stats/enrichment/{region}_profile.pdf`, `_heatmap.png`, `_matrix.gz` — profile plot, heatmap and deepTools matrix
- `stats/enrichment/{region}_profile_score.txt`, `_reads.txt`, `_peaks.txt` — the three measures in detail
- `annotations/tss.bed`, `annotations/regions/` — TSSs, reference points and merged regions
