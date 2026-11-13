# Consensus peaks

For every experimental condition with at least two replicates, CAWS builds one consensus peak set per peak caller, following the ENCODE approach for histone marks:

1. **Pool the replicates.** The filtered fragments of all replicates of the condition are combined. Paired-end samples contribute their sequenced fragments; single-end reads are extended to 147 bp along their strand. The IgG controls of all groups that the condition's samples belong to are pooled the same way when `igg_control` is `true`.
2. **Call peaks on the pool.** SEACR runs on a bedgraph of the pooled fragments (with the pooled IgG bedgraph, or `seacr_qvalue` without IgG); MACS3 runs on the pooled fragments as BEDPE (with the pooled IgG as control, or `--nolambda` without), in broad mode if `macs3_broad` is set.
3. **Keep peaks supported by the replicates.** A pooled peak is supported by a replicate if it overlaps one of that replicate's peaks by at least 50% of the length of either peak. Pooled peaks supported by at least `consensus_min_replicates` replicates (default 2, capped at the number of replicates) form the consensus set.

**Design decisions:**
- **Peak boundaries come from the pooled peak call.** Merging overlapping replicate peaks instead (e.g. with `bedtools merge`) chains peaks together, so consensus regions can become wider than any peak that was called, especially for broad marks and SEACR domains. Calling peaks on the pool keeps natural boundaries and needs no merge distance.
- **50% of either peak.** Requiring the overlap to cover half of the pooled peak *or* half of the replicate peak ignores small overlaps at peak edges, while a narrow replicate peak lying inside a broad pooled domain (common for single-end replicates of broad marks) still counts as support.
- **Mixed read layouts.** Extending single-end reads to 147 bp (one nucleosome, the same assumption MACS3 makes for single-end samples) lets single-end and paired-end replicates be pooled in one condition.
- **Pooled controls.** A condition can span several groups, each with its own IgG; pooling those controls gives a background that matches the pooled sample.
- **Single-replicate conditions are skipped.** One sample cannot show agreement between replicates; CAWS does not create pseudo-replicates.

**Outputs** (`peaks/consensus/`):
- `{seacr,macs3}/{condition}.consensus.bed` — consensus peaks: chromosome, start, end, number of supporting replicates, supporting sample names
- `seacr/{condition}.pooled.stringent.bed`, `macs3/{condition}.pooled_peaks.narrowPeak` (or `.broadPeak`) — all peaks called on the pooled replicates
- `stats/consensus_peaks.txt` — per condition and caller: number of replicates, support threshold, pooled peaks, consensus peaks and median width (also shown on the report's Sample Details page)
