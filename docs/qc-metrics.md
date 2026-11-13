# QC metrics

**`stats/stats.txt`** (one row per sample; also the table on the report's Alignment page). A *fragment* is a read pair for paired-end samples and a read for single-end samples. Alignment counts are taken from the merged alignments of each sample (all runs), the same way for FASTQ and BAM input.

| Column | Meaning |
|---|---|
| `SequencingDepth` | Fragments in the alignments, aligned or not (after trimming for FASTQ input) |
| `MappedFragNum` | Aligned fragments (paired-end: properly paired) |
| `UniqueMappedFragNum` | Aligned fragments without an alternative alignment (no `XS` tag on the read or either mate); matches Bowtie2's "aligned exactly 1 time" |
| `AlignmentRate` | `MappedFragNum` / `SequencingDepth` |
| `UniqueAlignmentRate` | `UniqueMappedFragNum` / `SequencingDepth` |
| `MtNumReads`, `MtAlignmentRate` | Fragments on `mt_chrom` (counted once per fragment, before duplicate removal) and their share of `SequencingDepth`. High values usually indicate open-chromatin background rather than targeted signal |
| `EstimatedLibrarySize` | Picard's estimate of the number of distinct molecules in the library (paired-end only; `NA` for single-end) |
| `DuplicationRate` | Percentage of aligned fragments marked as duplicates by Picard. Calculated for every sample, whether or not duplicates are removed |
| `Deduplicated` | Whether duplicates were removed for this sample (see [Duplicate Removal](duplicate-removal.md)) |
| `NRF`, `PBC1`, `PBC2` | ENCODE library complexity before duplicate removal, after the MAPQ and chromosome filters: NRF = distinct fragments / fragments; PBC1 = positions with exactly one fragment / distinct positions; PBC2 = positions with one / positions with two fragments. A position is a fragment's start and end (paired-end) or a read's 5′ end and strand (single-end) |
| `PCRBottleneck` | ENCODE grading of PBC1: severe (< 0.5), moderate (0.5–0.8), mild (0.8–0.9), none (≥ 0.9) |
| `UnalignedFragNum` | Unaligned fragments in the alignments; `0` usually means a BAM that was filtered to aligned reads, in which case alignment rate and contamination cannot be judged |
| `SpikeinFragNum`, `SpikeinRate` | Only when `spikein_genome` is set: fragments that did not align to the host and align to the spike-in genome (MAPQ ≥ `minquality`, deduplicated when the sample is), and their percentage of `SequencingDepth`. See [Normalization](normalization.md) |
| `UniqueFragNum` | Aligned fragments that are not duplicates, i.e. the number remaining after deduplication (for samples that are not deduplicated, the number that *would* remain) |
| `<Genome>OnlyRate` | One column per contaminant genome: percentage of screened reads (a subsample of `contaminant_subset` read 1 sequences) that map to that genome and not to the host or any other screened genome (FastQ Screen "one genome" hits) |

**`stats/peaks_consolidated.txt`** (one row per experimental sample; `seacr_` and `macs3_` columns for each caller).

| Column | Meaning |
|---|---|
| `peakN` | Number of peaks |
| `peakReprodNum`, `peakReprodRate` | Peaks that overlap a peak of any other replicate of the same condition, as a count and as a percentage of `peakN` (`NA` for conditions with one sample) |
| `FRiP` | Fraction of reads in peaks, in percent: filtered fragments overlapping a peak / all filtered fragments |
| `median_width`, `mean_width` | Peak width in bp |
| `pct_narrow`, `pct_broad` | Percentage of peaks ≤ 500 bp and > 2000 bp |
| `median_gc` | Median GC fraction of the peak sequences |

Consensus peak statistics (`stats/consensus_peaks.txt`) are described under [Consensus Peaks](consensus-peaks.md) and enrichment statistics (`stats/region_enrichment.txt`) under [Region Enrichment](region-enrichment.md).
