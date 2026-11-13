# Configuration

Start from the template `config.json` in the repository. Keys marked *optional* can be left out and take the default shown; all others must be present.

| Parameter | Description |
|------------|-------------|
| `samplesheet` | Path to TSV file describing samples (columns: `sampleID`, `condition`, `group`, `read1`, `read2`) |
| `reference_fa` | Reference genome FASTA file |
| `bt2_idx` | Bowtie2 index prefix for the reference genome |
| `reference_fai` | FASTA index file |
| `chromosome_file` | Text file with one chromosome name per line; only reads on these chromosomes are used for peak calling, FRiP, coverage correlation and coverage tracks |
| `blacklist` | BED file of genomic regions to exclude |
| `contaminant_genomes` | Map of contaminant name → Bowtie2 index prefix, screened with FastQ Screen against the host genome (e.g., `{"Ecoli": "...", "Bovine": "..."}`). Names must be valid identifiers (letters, digits, underscore; not starting with a digit); `Host` is reserved |
| `contaminant_subset` | Number of read 1 sequences sampled per library for FastQ Screen (the template uses `1000000`; `0` screens all reads) |
| `spikein_genome` | *Optional*, default `""`: name of the `contaminant_genomes` entry that serves as spike-in, usually `"Ecoli"` (carry-over DNA from pA-Tn5). When set, spike-in fragments are counted for every sample and reported in `stats/stats.txt`. See [Normalization](normalization.md) |
| `normalization` | *Optional*, default `"cpm"`: scaling of the coverage tracks (TDF, bigWig) and region signal profiles. `"cpm"` scales to fragments per million filtered fragments. `"spikein"` scales to fragments per 10,000 spike-in fragments and needs `spikein_genome`. Peak calling is not affected. See [Normalization](normalization.md) |
| `trim_adapters` | `true`/`false` — enable adapter trimming via Trim Galore. Also selects the Bowtie2 mode: `--local` with trimming, `--end-to-end` without |
| `igg_control` | `true`/`false` — use IgG controls as background for peak calling (SEACR control bedgraph, MACS3 `-c`). Requires exactly one `Control` per samplesheet `group` |
| `mt_chrom` | Name of mitochondrial chromosome (e.g., `"chrM"`) |
| `dedup` | *Optional*, default `"controls"`: remove PCR duplicates (Picard MarkDuplicates) from IgG controls only. `true` removes them from all samples, `false` from none. A mapping sets it per condition, e.g. `{"default": "controls", "BG4": true}`; the samplesheet `dedup` column overrides both per sample. See [Duplicate Removal](duplicate-removal.md) |
| `minquality` | Minimum mapping quality (MAPQ) for reads used in peak calling, FRiP, coverage correlation, coverage tracks and spike-in counts |
| `fragment_min`, `fragment_max` | Minimum/maximum fragment length for paired-end alignment (Bowtie2 `-I`/`-X`). Fragments longer than 1000 bp are additionally removed before peak calling |
| `genome_size` | Genome size for MACS3 peak calling: `"hs"` (human), `"mm"` (mouse), `"ce"` (C. elegans), `"dm"` (Drosophila), or numeric value (e.g., `2.7e9`) |
| `seacr_qvalue` | Used only when `igg_control` is `false`: SEACR's numeric threshold, i.e. the top fraction of peaks kept by total signal (e.g. `0.01` = top 1%). With IgG, SEACR uses the control bedgraph instead |
| `macs3_qvalue_with_control`, `macs3_qvalue_no_control` | MACS3 q-value cutoffs with and without IgG control. Single-end samples are called with `--nomodel --extsize 147` (Tn5 cuts lack the strand bimodality MACS3 models); paired-end samples use actual fragments (`-f BAMPE`) |
| `macs3_broad`, `macs3_broad_cutoff` | *Optional*, defaults `false` and `0.1`: call broad peaks with MACS3 (`--broad --broad-cutoff`), e.g. for H3K27me3 or H3K36me3. Output is then `*_peaks.broadPeak` instead of `*_peaks.narrowPeak` |
| `consensus_min_replicates` | *Optional*, default `2`: minimum number of replicates that must support a pooled peak for it to enter the consensus peak set (capped at the number of replicates in the condition) |
| `heatmap_window` | Distance (bp) on each side of SEACR peak summits for the peak heatmaps |
| `igv_genome` | Genome for IGV TDF generation and the IGV session: an igvtools genome ID (e.g., `"hg38"`, `"mm10"`) or a path to a chrom.sizes file (the session then uses `reference_fa` as the genome) |
| `gtf_file` | Path to GTF file; its protein-coding TSSs become the built-in `TSS` region set (optional; `""` to skip) |
| `enrichment_regions` | *Optional*, default `{}`: map of name → BED file (BED, BED.gz or narrowPeak) of regions to test for enrichment, e.g. `{"PQS": "/path/pqs.bed", "ATAC": "/path/atac_peaks.narrowPeak"}`. See [Region Enrichment](region-enrichment.md) |
| `outdir` | Output directory for analysis results. Snakemake runs inside this directory |

**Example `config.json`:**

```json
{
    "samplesheet": "/data/CAWS/samplesheet.tsv",
    "reference_fa": "/share/gatk_bundle/hg38/hs38DH.fa",
    "bt2_idx": "/share/gatk_bundle/hg38/hs38DH",
    "reference_fai": "/share/gatk_bundle/hg38/hs38DH.fa.fai",
    "chromosome_file": "/share/gatk_bundle/hg38/hs38DH.primary.txt",
    "blacklist": "/share/gatk_bundle/hg38/blacklist.bed.gz",
    "contaminant_genomes": {
        "Ecoli": "/share/gatk_bundle/ecoli/MG1655",
        "Bovine": "/share/gatk_bundle/bovine/ARS-UCD1.2"
    },
    "contaminant_subset": 1000000,
    "trim_adapters": true,
    "igg_control": false,
    "mt_chrom": "chrM",
    "dedup": "controls",
    "minquality": 3,
    "fragment_min": 10,
    "fragment_max": 700,
    "genome_size": "hs",
    "seacr_qvalue": 0.01,
    "macs3_qvalue_with_control": 0.01,
    "macs3_qvalue_no_control": 0.001,
    "heatmap_window": 3000,
    "gtf_file": "",
    "igv_genome": "hg38",
    "outdir": "/data/CAWS/20230324"
}
```
