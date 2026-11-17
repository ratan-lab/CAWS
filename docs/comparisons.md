# Comparing conditions

To ask whether a perturbation (a knockout, a knockdown, a drug) creates new peaks or changes existing ones, CAWS can prepare the input for a differential analysis with DESeq2, edgeR or DiffBind: one reference region set for the conditions being compared, and the number of fragments of every sample in every region.

```json
"comparisons": {
    "DHX36_KD": ["GM13069_cells_BG4_siNC", "GM13069_cells_BG4_siDHX36"],
    "CPT_in_KO": ["Jurkat_nuclei_TOP2B_KO1_UT", "Jurkat_nuclei_TOP2B_KO1_CPT"]
}
```

Each entry names a comparison and lists two or more conditions from the samplesheet `condition` column (not `Control`). Comparison names may contain letters, digits and `. _ -`. Leave `comparisons` out, or set it to `{}`, to skip this step.

## How the regions are made

1. **Pooled peak call.** The filtered fragments of all samples of the comparison are pooled, after downsampling each sample to the depth of the smallest one, so that no library or condition dominates the pool. MACS3 calls peaks on the pool with a relaxed threshold (`comparison_qvalue`, default 0.05) and `--call-summits`, so that neighboring peaks keep their own summits. The IgG controls of the conditions' groups are pooled as control when `igg_control` is true.
2. **Replicate support.** A pooled peak is supported by a replicate if it overlaps one of that replicate's MACS3 peaks by at least 50% of either peak (the rule used for [consensus peaks](consensus-peaks.md)). A pooled peak is kept if it is supported by at least `consensus_min_replicates` replicates (capped at the number of replicates) of **at least one** condition. A peak present only in the knockout is therefore kept.
3. **Summit-centered regions.** Each kept peak becomes a region of ±`comparison_window` bp (default 250) around its summit. Where regions overlap, the more significant one is kept, so regions are of equal width and never merged. With `macs3_broad: true` there are no summits, and the pooled broad peaks are used as they are.
4. **Counts.** Every sample's filtered fragments (all of them, not the downsampled pool) are counted in the region that contains their midpoint, so each fragment is counted at most once.

The relaxed pooled call finds peaks that are weak in each sample on its own. The replicate-support filter removes the false positives that the relaxed threshold lets through.

## Outputs

`comparisons/{name}/`:

| File | Content |
|------|---------|
| `regions.tsv` | one row per region: coordinates, `region_id`, summit, −log10 q of the pooled peak, and for each condition the number of supporting replicates (e.g. `2/2`) |
| `regions.bed` | the regions as BED (chrom, start, end, `region_id`), e.g. for IGV |
| `counts.tsv` | fragments per region (rows) and sample (columns, by `sampleID`) |
| `samples.tsv` | one row per sample: `sampleID`, `label`, `condition`, `filtered_fragments` (as used for CPM) and, with `spikein_genome`, `spikein_fragments` |
| `pooled_peaks.narrowPeak` | the pooled MACS3 call before the support filter (`pooled_peaks.broadPeak` with `macs3_broad: true`) |

A region supported only in one condition (e.g. `0/2` and `2/2`) is a candidate for a gained or lost peak. Whether it differs is for the statistical test to decide: peak calls depend on thresholds and depth, so "called in one condition only" is not evidence on its own.

## Differential analysis

With two replicates per condition DESeq2, edgeR and DiffBind work: they estimate variability across all regions, not from the two replicates of a single region. Power is lower than with three or more replicates, but the tests are valid. With one replicate per condition there is no variance estimate, and only fold changes can be ranked.

```r
library(DESeq2)
counts  <- read.delim("comparisons/DHX36_KD/counts.tsv", row.names = 1, check.names = FALSE)
samples <- read.delim("comparisons/DHX36_KD/samples.tsv", row.names = 1)
samples$condition <- factor(samples$condition)
dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ condition)

# spike-in normalization (if the spike-in counts are high enough); omit for DESeq2's default
sizeFactors(dds) <- samples$spikein_fragments / mean(samples$spikein_fragments)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "GM13069_cells_BG4_siDHX36", "GM13069_cells_BG4_siNC"))
```

**Normalization decides the answer when the perturbation changes binding globally.** DESeq2's and edgeR's default normalizations assume that most regions do not change, so a global gain or loss is normalized away. Spike-in factors ([Normalization](normalization.md)) or `filtered_fragments` (library size) do not make that assumption. Running the test with two normalizations and reporting the regions that hold up under both is a robust choice. If they disagree a lot, the perturbation probably changed binding globally.

**Checks with few replicates.** Replicates should cluster by condition in a PCA of the counts (e.g. DESeq2's `plotPCA(vst(dds))`), and the MA plot should be centered near zero unless a global change is expected.

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `comparisons` | `{}` | map of comparison name → list of conditions |
| `comparison_window` | `250` | half-width (bp) of the regions around the summits; ±250 bp suits transcription factors, chromatin proteins and DNA structures (BG4, S9.6) |
| `comparison_qvalue` | `0.05` | MACS3 q-value cutoff for the pooled call |

Adding a comparison to a finished run rebuilds the filtered fragment files of the samples involved (they are temporary), which takes about as long as the original filtering step. The replicate support uses `consensus_min_replicates`. For broad marks (H3K27me3, H3K36me3), set `macs3_broad: true` to use the pooled broad peaks as regions.
