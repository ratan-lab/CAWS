# Duplicate removal

By default (`dedup: "controls"`), CAWS removes PCR duplicates from IgG controls only and keeps them in experimental samples.

**Why targets are not deduplicated by default.** CUT&Tag libraries are PCR-amplified, but pA-Tn5 is tethered at its target and inserts preferentially into accessible DNA nearby, so independent cells often produce fragments with identical ends. In high-quality libraries most apparent duplicates at target sites are therefore real signal, and removing them caps the signal exactly where it is strongest. This is especially true for single-end data, where a "duplicate" only needs the same start position and strand. IgG controls, in contrast, come from little DNA, have low complexity and are dominated by true PCR duplicates, which would otherwise inflate the background used for peak calling. Without unique molecular identifiers (UMIs) the two cannot be told apart, so this is a trade-off; nf-core/cutandrun uses the same default.

**When to deduplicate target samples.** The default suits abundant targets such as histone marks, which are spread over the genome, amplify like the background and do not run out. Low-abundance targets and DNA/RNA structures behave differently: few tethered pA-Tn5 complexes and few target molecules give low-complexity libraries that need more PCR cycles, so most duplicates are PCR copies; structured, GC-rich target fragments (e.g. G-quadruplexes) amplify less efficiently than background fragments, so PCR over-represents the background; and when targets are scarce, Tn5 cuts the same accessible sites repeatedly.

| Target | Suggested setting |
|---|---|
| Abundant histone marks (H3K27me3, H3K4me3, H3K36me3, …) | `"controls"` (default) |
| Transcription factors, other chromatin proteins | deduplicate (`true` for the condition) |
| DNA/RNA structures (BG4, S9.6, Z-DNA antibodies, …) | deduplicate (`true` for the condition) |
| Any low-input or heavily amplified library | deduplicate that sample (samplesheet `dedup` = `yes`) |

The library complexity metrics help decide: a PBC1 below 0.8 (moderate or severe PCR bottleneck), a low NRF, or a high `DuplicationRate` with a small `EstimatedLibrarySize` compared with other libraries of the same mark indicate a library dominated by PCR duplicates (see [QC Metrics](qc-metrics.md)). For abundant marks, real Tn5 "duplicates" also lower these values somewhat, which is why CAWS reports them rather than deciding automatically.

**Configuration.** `dedup` in `config.json` is `"controls"`, `true` or `false`, or a mapping per condition, e.g.

```json
"dedup": {"default": "controls", "BG4": true, "S9.6": true}
```

The samplesheet `dedup` column (`yes`/`no`) overrides the condition setting for a sample, which overrides the default. CAWS warns at startup when replicates of one condition are treated differently, because their consensus peaks then pool deduplicated and non-deduplicated data.

**How it works.** The runs of a sample are merged first, so duplicates are found across all sequencing of a library. Picard MarkDuplicates runs on every sample, so duplication metrics are always reported. Where duplicates are removed, downstream steps use the deduplicated alignments; otherwise they use the alignments unchanged. Both are stored in `dedupalignments/`, and the `Deduplicated` column of `stats/stats.txt` records the choice for each sample. MACS3 is run with `--keep-dup all` (as in nf-core/cutandrun, ePeak and CUT&RUNTools), so it does not remove duplicates on its own and this setting alone decides how duplicates are handled.
