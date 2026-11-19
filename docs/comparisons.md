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

`counts.tsv` is the count matrix and `samples.tsv` the sample table for DESeq2 or edgeR. Use the raw counts in `counts.tsv`, never CPM or bigWig values. The statistical analysis is left to the user. The steps below are a recommended starting point for DESeq2.

With two replicates per condition DESeq2, edgeR and DiffBind work: they estimate variability across all regions, not from the two replicates of a single region. Power is lower than with three or more replicates, but the tests are valid. With one replicate per condition there is no variance estimate, and only fold changes can be ranked.

**1. Read the data.** The columns of `counts.tsv` are the sampleIDs and match the rows of `samples.tsv`. Conditions are listed in the order of the comparison in the configuration. Setting the factor levels in that order makes the first listed condition the reference.

```r
library(DESeq2)
counts  <- read.delim("comparisons/DHX36_KD/counts.tsv", row.names = 1, check.names = FALSE)
samples <- read.delim("comparisons/DHX36_KD/samples.tsv", row.names = 1)
samples$condition <- factor(samples$condition, levels = unique(samples$condition))
dds <- DESeqDataSetFromMatrix(counts, samples[colnames(counts), ], design = ~ condition)
```

**2. Pre-filter lightly.** Every region already passed the replicate-support filter, so drop only regions with almost no fragments. `results()` also filters low-count regions on its own (independent filtering), so stricter pre-filtering is not needed.

```r
smallest_group <- min(table(samples$condition))
dds <- dds[rowSums(counts(dds) >= 10) >= smallest_group, ]
```

**3. Choose the normalization.** DESeq2's and edgeR's default normalizations assume that most regions do not change. Counted over peaks only, this assumption is stronger than for RNA-seq, and a perturbation that changes binding globally is normalized away. Library size (`filtered_fragments`) and spike-in counts ([Normalization](normalization.md)) do not make that assumption.

```r
geo <- function(x) x / exp(mean(log(x)))
dds_default <- estimateSizeFactors(dds)
dds_libsize <- dds; sizeFactors(dds_libsize) <- geo(samples[colnames(dds), "filtered_fragments"])
dds_spikein <- dds; sizeFactors(dds_spikein) <- geo(samples[colnames(dds), "spikein_fragments"])  # with spikein_genome
```

If the normalizations give similar results, the default is fine. If they disagree a lot, the perturbation probably changed binding globally. The spike-in factors are then the most trustworthy, provided each sample has thousands of spike-in fragments. Report which normalization was used, and preferably the regions that hold up under two of them.

**4. Account for how replicates were processed.** If the replicates were processed in pairs (for example replicate 1 of every condition on the same day or with the same pA-Tn5 batch), add the pairing to the design. Add a `replicate` column to `samples` first. With two replicates per condition this costs one degree of freedom but often removes much of the noise.

```r
design(dds) <- ~ replicate + condition
```

**5. Check the samples.** Replicates should cluster by condition in a PCA, and the dispersion estimates should decrease smoothly with the mean. With two replicates, one outlier dominates the results, so check before trusting the list of regions. `vst()` needs at least 1,000 regions. For smaller region sets use `varianceStabilizingTransformation(dds, blind = TRUE)`.

```r
plotPCA(vst(dds, blind = TRUE), intgroup = "condition")
dds <- DESeq(dds)
plotDispEsts(dds)
```

**6. Test for a minimum effect and shrink fold changes.** Testing against a fold-change threshold (`lfcThreshold`, here 1.4-fold) gives valid p-values for "changed by more than", which filtering by fold change after the test does not. Shrunken fold changes keep regions with few fragments from dominating rankings, plots and heatmaps. Use them for those, and the unshrunken results for significance.

```r
res <- results(dds, contrast = c("condition", "GM13069_cells_BG4_siDHX36", "GM13069_cells_BG4_siNC"),
               lfcThreshold = 0.5, alpha = 0.05)
res_shrunk <- lfcShrink(dds, contrast = c("condition", "GM13069_cells_BG4_siDHX36", "GM13069_cells_BG4_siNC"),
                        type = "ashr")
plotMA(res_shrunk)
```

The MA plot should be centered near zero unless a global change is expected.

**7. Interpret the regions.** Join the results to `regions.tsv` by `region_id` for the coordinates, summits and replicate support. Regions supported only in one condition (e.g. `0/2` and `2/2`) and significantly changed are the clearest gained or lost peaks. Inspect the top regions in IGV with the normalized tracks, and annotate them, for example with ChIPseeker.

```r
res_df <- as.data.frame(res) |> tibble::rownames_to_column("region_id") |>
  dplyr::left_join(read.delim("comparisons/DHX36_KD/regions.tsv", check.names = FALSE), by = "region_id")
```

**Comparisons of more than two conditions.** A comparison can list more than two conditions. All of them then share one region set and one count matrix, and DESeq2 estimates variability from all samples. Test for any difference with `DESeq(dds, test = "LRT", reduced = ~ 1)`, then follow up with specific contrasts. For two factors, such as genotype and treatment, split `condition` into one column per factor and use a design such as `~ genotype + treatment + genotype:treatment`.

**Duplicate removal.** Handle duplicates the same way in all samples of a comparison (`dedup`). Otherwise differences in duplicate removal appear as differences in binding.

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `comparisons` | `{}` | map of comparison name → list of conditions |
| `comparison_window` | `250` | half-width (bp) of the regions around the summits. ±250 bp suits transcription factors, chromatin proteins and DNA structures (BG4, S9.6) |
| `comparison_qvalue` | `0.05` | MACS3 q-value cutoff for the pooled call |

Adding a comparison to a finished run rebuilds the filtered fragment files of the samples involved (they are temporary), which takes about as long as the original filtering step. The replicate support uses `consensus_min_replicates`. For broad marks (H3K27me3, H3K36me3), set `macs3_broad: true` to use the pooled broad peaks as regions.
