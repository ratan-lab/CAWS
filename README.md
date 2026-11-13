# CAWS: CUT&Tag Analysis Workflow System  
## *A Snakemake-based pipeline for comprehensive CUT&Tag data analysis*  

[![test](https://github.com/ratan-lab/CAWS/actions/workflows/test.yml/badge.svg)](https://github.com/ratan-lab/CAWS/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-≥3.11-blue.svg)](https://python.org)

---

## Overview  

**CAWS (CUT&Tag Analysis Workflow System)** is a comprehensive Snakemake pipeline for analyzing CUT&Tag (Cleavage Under Targets and Tagmentation) sequencing data. CUT&Tag is an epigenomic profiling method that maps protein–DNA interactions and histone modifications in situ, with lower background and fewer cells than ChIP-seq (Kaya-Okur et al., 2019).  

---

## Why Use CAWS?

- **CUT&Tag Optimized** – Purpose-built for CUT&Tag data with specialized QC and peak-calling strategies
- **Flexible Input** – Supports single-end, paired-end, and mixed datasets in the same analysis
- **Comprehensive Analysis** – End-to-end workflow from raw reads to interactive visualizations
- **Dual Peak Calling** – Generates peak calls using both SEACR and MACS3 for robust peak detection
- **Interactive Reports** – HTML dashboards with Plotly-based visualizations
- **Reproducible** – Per-rule conda environments with pinned tool and R package versions (`workflow/envs/`)

---

## Target Audience  

CAWS is designed for:  
- Bioinformaticians analyzing CUT&Tag datasets  
- Researchers studying chromatin biology or epigenomics  
- Laboratories seeking standardized workflows  
- Users requiring comprehensive QC and visualization  

---

## Prerequisites

### Required Software
- **Snakemake** ≥9.0.0 (CI tests the latest release; also used with 9.8.1)
- **Conda** ≥24.7.1
- **Python** ≥3.11 (required by Snakemake 9)

All other tools, including R and its packages, are installed automatically by Snakemake into per-rule conda environments (`workflow/envs/`). `envs/snakemake-runner.yaml` provides an environment with Snakemake, the SLURM executor plugin and pandas.

**Important Notes**:
- Snakemake 9.0+ requires Conda ≥24.7.1 for environment management
- If you encounter conda version errors, see Troubleshooting section below  

### System Requirements  
- **Storage**: ~50 GB free space (typical dataset)  
- **Memory**: ≥8 GB RAM (16 GB+ recommended)  
- **CPU**: ≥4 cores  
- **OS**: Linux (x86-64); HPC cluster execution recommended. macOS and Windows are not supported because several conda environments pin Linux (linux-64) builds.  

### Input Data Requirements
- **FASTQ** files (gzipped) or **BAM** alignments — **single-end**, **paired-end** or both, and several sequencing runs per library
- Reference genome (**FASTA**, its `.fai` index, and a Bowtie2 index)
- **Chromosome list** (plain text, one chromosome name per line) — only reads on these chromosomes are used for peak calling and coverage
- **Blacklist** regions (**BED**)
- **Contaminant genomes** (Bowtie2 indexes, at least one, e.g. *E. coli*). *E. coli* also serves as the spike-in for optional spike-in normalization
- **Annotation** file (**GTF**, optional, for TSS enrichment) and other **region sets** (BED, optional) for enrichment analysis
- **Sample metadata** (**TSV**)

**Use absolute paths** in `config.json` and the samplesheet. The pipeline runs inside `outdir`, so relative paths are resolved from there, not from where you launch Snakemake.

---

## Quick Start  

### 1. Clone and Set Up

```bash
# Clone the repository
git clone https://github.com/ratan-lab/CAWS.git
cd CAWS
```

Create the environment that runs Snakemake (once):

```bash
conda env create -f envs/snakemake-runner.yaml
conda activate caws-snakemake
```

The tools used by each step are installed automatically by Snakemake (`--use-conda`).

---

### 2. Configure SLURM Profile (for HPC execution)

Edit `profiles/slurm/config.yaml` to set your SLURM allocation and partition:

```yaml
default-resources:
  - runtime=180m        # always give a unit: a bare number is read as seconds
  - mem_mb=9000
  - threads=1
  - slurm_account=your-allocation-name  # e.g., ratan
  - slurm_partition=your-partition      # e.g., standard
```

**Skip this step** if running locally (not on a cluster).

---

### 3. Prepare Analysis Configuration Files  

**Sample Sheet (`samplesheet.tsv`)**

The pipeline supports **paired-end**, **single-end**, and **mixed** datasets. Specify single-end samples by using `-`, `NA`, or similar values in the `read2` column:

```tsv
sampleID	condition	group	read1	read2
sample1_PE	Treatment	group1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
sample2_SE	Treatment	group2	/path/to/sample2.fastq.gz	-
sample3_PE	Control	group1	/path/to/sample3_R1.fastq.gz	/path/to/sample3_R2.fastq.gz
sample4_SE	Control	group2	/path/to/sample4.fastq.gz	NA
```

- `condition`: `Control` (exact, case-sensitive) marks IgG controls; any other value is an experimental condition. Peaks are called only for experimental samples. Samples with the same condition are treated as replicates for peak reproducibility.
- `group`: links experimental samples to their IgG control. With `igg_control: true`, each group must contain exactly one `Control` and at least one experimental sample, and experimental samples and their control must have the same read type (both SE or both PE).
- **Several rows with the same `sampleID`** are sequencing runs (lanes, re-sequencing) of one library, i.e. technical replicates. Each run is trimmed and aligned separately and the alignments are merged before duplicate removal, so duplicates are found across all runs and the library counts as **one** replicate. The rows must agree on `condition`, `group` and read layout. **Every distinct `sampleID` is treated as a biological replicate.**
- **BAM input**: if `read1` ends in `.bam`, the sample is taken as aligned (leave `read2` empty), e.g. alignments made with settings you need. Single-/paired-end is detected from the BAM, which must be aligned to the configured reference (chromosome names and lengths are checked). Trimming and Bowtie2 are skipped; FastQC runs on the BAM and the contaminant screen on read 1 extracted from it. FASTQ and BAM rows can be mixed. Note that MAPQ scales differ between aligners (`minquality`), and that a BAM without unaligned reads gives a 100% alignment rate and underestimates contamination.
- `dedup` (optional column): `yes` or `no` overrides the global `dedup` setting for that sample, e.g. to remove duplicates from a low-input library expected to be over-amplified. Leave empty (or omit the column) to use the global setting.

**Single-end indicators**: Use any of these values in the `read2` column to mark a sample as single-end: `-`, `NA`, `na`, `n/a`, `null`, `none`, `nan`, `nil`, `empty`, `0`, `false`, `missing`, `absent`, `single`, or `se`.

**Configuration File (`config.json`)**

Edit `config.json` with your file paths and parameters. See the Configuration Details section below for a complete description of all parameters.

---

### 4. Run the Analysis

**Local Execution (small datasets)**
```bash
snakemake --cores 8 --use-conda --configfile config.json
```

**SLURM Cluster Execution (recommended for large datasets)**

**Option 1: Simple Execution with Profile**

Submit as a SLURM job:
```bash
# Create a submission script
cat > run_caws.sh << 'EOF'
#!/bin/bash
#SBATCH -A your-allocation      # Change to your allocation name
#SBATCH -p standard             # Partition for controller job
#SBATCH -c 1                    # Controller only needs 1 core
#SBATCH -t 48:00:00             # Max runtime for controller
#SBATCH --mem=4G                # Memory for controller
#SBATCH -o caws_%j.out          # Output log
#SBATCH -e caws_%j.err          # Error log
#SBATCH -J caws_pipeline        # Job name

# Activate the Snakemake environment (created once with: conda env create -f envs/snakemake-runner.yaml)
module load miniforge   # if conda is provided as a module on your cluster
conda activate caws-snakemake

# Run pipeline with profile
snakemake \
  --profile profiles/slurm \
  --configfile config.json
EOF

# Submit the job
sbatch run_caws.sh
```

**Option 2: Advanced Execution with Custom Settings**

For more control over execution and conda environment caching:
```bash
#!/bin/bash
#SBATCH -A your-allocation
#SBATCH -p standard
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH --mem=4G
#SBATCH -o caws_%j.out
#SBATCH -e caws_%j.err
#SBATCH -J caws_pipeline

module load miniforge   # if conda is provided as a module on your cluster
conda activate caws-snakemake

# Run pipeline with explicit parameters
snakemake \
  --profile /path/to/CAWS/profiles/slurm \
  --snakefile /path/to/CAWS/workflow/Snakefile \
  --configfile config.json \
  --cores 1 \
  --use-conda \
  --conda-prefix /scratch/your-user/conda-envs \
  --latency-wait 60 \
  --rerun-incomplete \
  --keep-going \
  --jobs 100
```

**Interactive Execution:**
```bash
module load miniforge   # if conda is provided as a module on your cluster
conda activate caws-snakemake

# Run in current session (use screen or tmux recommended)
snakemake --profile profiles/slurm --configfile config.json
```

**Notes:**
- `--conda-prefix` caches conda environments for faster subsequent runs
- `--jobs` overrides the profile's limit on concurrent SLURM jobs (100)
- `--notemp` can be added to preserve temporary files for debugging
- Use absolute paths if running from a different directory
- Temporary files (sorting, Picard) go to `$TMPDIR` as set where `snakemake` is launched, and that path is passed to every job. Point it at a directory on shared scratch that the compute nodes can see, e.g. `export TMPDIR=/scratch/$USER/tmp` (after `mkdir -p`) before the `snakemake` command, or set it with `--default-resources tmpdir=/scratch/$USER/tmp`. The node-local `/tmp` works but is often small; a directory that exists only where `snakemake` runs (e.g. a per-job directory of the controller job) makes jobs fail

---

### 5. View Results  

Open `<outdir>/reports/cutntag_analysis_report.html` in a web browser.

---

## Testing the Installation

The repository includes a small test dataset in `.test/data` (about 7 MB): H3K27me3 and IgG CUT&Tag reads from K562 cells (Kaya-Okur et al., 2019), subset to a 4 Mb window of chr22, with paired-end and single-end samples. See `.test/data/PROVENANCE.txt` for how it was made. To run it:

```bash
conda env create -f envs/snakemake-runner.yaml
conda create -n caws-setup -c conda-forge -c bioconda bowtie2=2.5.1 samtools

# Build indexes and write the test samplesheet and config into .test/run/
conda run -n caws-setup bash .test/setup.sh

conda activate caws-snakemake
snakemake -s workflow/Snakefile --configfile .test/run/config.json --use-conda --cores 4
bash .test/check_outputs.sh
```

Results are written to `.test/run/out/`. The first run takes longer because Snakemake creates the conda environments for each tool. GitHub Actions (`.github/workflows/test.yml`) runs this test on every push to `main` and on pull requests, together with dry runs of all `igg_control`/`dedup` combinations.

---

## Example Analysis with Public Data

This example demonstrates a complete CAWS analysis using publicly available CUT&Tag data. We'll analyze H3K4me3 (active promoter mark) samples from Kaya-Okur et al. (2019), *Nature Communications* 10:1930.

### Dataset Information

- **Study**: CUT&Tag profiling of H3K4me3 in human K562 cells
- **Data Source**: NCBI GEO [GSE124557](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124557)
- **Sample Size**: 4 samples (2 H3K4me3 + 2 IgG controls)

| Sample | GEO | SRA runs |
|---|---|---|
| H3K4me3_rep1 | GSM3536518 | SRR8383516, SRR8383517 |
| H3K4me3_rep2 | GSM3680226 | SRR8754610 |
| IgG_rep1 | GSM3560264 | SRR8435051, SRR8435052 |
| IgG_rep2 | GSM3680227 | SRR8754611, SRR8754612 |

- **Sequencing**: Paired-end Illumina
- **Reference**: Human (hg38)

### Step 1: Download Example Data

```bash
# Create working directory
mkdir -p cutntag_example && cd cutntag_example

# Download FASTQ files from SRA using SRA Toolkit
# Install SRA Toolkit if needed: conda install -c bioconda sra-tools

for run in SRR8383516 SRR8383517 SRR8754610 SRR8435051 SRR8435052 SRR8754611 SRR8754612; do
    fasterq-dump --split-files "$run"
    gzip "${run}"_1.fastq "${run}"_2.fastq
done
```

Some libraries were sequenced in two runs; they are listed as two rows with the same `sampleID` in the samplesheet (Step 3), and CAWS merges them.

### Step 2: Prepare Reference Files

```bash
# Download human hg38 reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build Bowtie2 index
bowtie2-build hg38.fa hg38_index

# Index FASTA file
samtools faidx hg38.fa

# Download blacklist regions
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz

# List of chromosomes to analyze, one name per line (here: primary chromosomes, without chrM)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
cut -f1 hg38.chrom.sizes | grep -vE '_|chrM|chrEBV' > hg38.chromosomes.txt

# Download GTF annotation for TSS enrichment (optional)
# GENCODE uses UCSC-style chromosome names (chr1), matching the hg38 reference above.
# Ensembl GTFs also work, but use "1" instead of "chr1" and must be paired with an Ensembl-style reference.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
gunzip gencode.v47.annotation.gtf.gz

# Contaminant genomes for FastQ Screen (QC only)
# E. coli (pA-Tn5 carry-over)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
bowtie2-build GCF_000005845.2_ASM584v2_genomic.fna ecoli_index
# Bovine (BSA / serum carry-over); building this index takes a few hours
wget https://ftp.ensembl.org/pub/release-110/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
bowtie2-build --threads 8 Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz bovine_index
```

### Step 3: Create Sample Sheet

Create `samplesheet.tsv`:

```tsv
sampleID	condition	group	read1	read2
H3K4me3_rep1	H3K4me3	rep1	/path/to/SRR8383516_1.fastq.gz	/path/to/SRR8383516_2.fastq.gz
H3K4me3_rep1	H3K4me3	rep1	/path/to/SRR8383517_1.fastq.gz	/path/to/SRR8383517_2.fastq.gz
H3K4me3_rep2	H3K4me3	rep2	/path/to/SRR8754610_1.fastq.gz	/path/to/SRR8754610_2.fastq.gz
IgG_rep1	Control	rep1	/path/to/SRR8435051_1.fastq.gz	/path/to/SRR8435051_2.fastq.gz
IgG_rep1	Control	rep1	/path/to/SRR8435052_1.fastq.gz	/path/to/SRR8435052_2.fastq.gz
IgG_rep2	Control	rep2	/path/to/SRR8754611_1.fastq.gz	/path/to/SRR8754611_2.fastq.gz
IgG_rep2	Control	rep2	/path/to/SRR8754612_1.fastq.gz	/path/to/SRR8754612_2.fastq.gz
```

Rows with the same `sampleID` are the runs of one library and are merged. Each H3K4me3 replicate is paired (`group`) with the IgG control from the same series in GSE124557 (same experimenter prefix, HK or SH), and the two H3K4me3 samples share a `condition`, so they are compared for reproducibility.

**Important**: Replace `/path/to/` with absolute paths to your downloaded FASTQ files.

### Step 4: Create Configuration File

Create `config.json`:

```json
{
    "samplesheet": "/path/to/samplesheet.tsv",
    "reference_fa": "/path/to/hg38.fa",
    "bt2_idx": "/path/to/hg38_index",
    "reference_fai": "/path/to/hg38.fa.fai",
    "chromosome_file": "/path/to/hg38.chromosomes.txt",
    "blacklist": "/path/to/hg38-blacklist.v2.bed",
    "contaminant_genomes": {
        "Ecoli": "/path/to/ecoli_index",
        "Bovine": "/path/to/bovine_index"
    },
    "contaminant_subset": 1000000,
    "spikein_genome": "Ecoli",
    "normalization": "cpm",
    "gtf_file": "/path/to/gencode.v47.annotation.gtf",
    "trim_adapters": true,
    "igg_control": true,
    "mt_chrom": "chrM",
    "dedup": "controls",
    "minquality": 10,
    "fragment_min": 0,
    "fragment_max": 700,
    "genome_size": "hs",
    "seacr_qvalue": 0.01,
    "macs3_qvalue_with_control": 0.05,
    "macs3_qvalue_no_control": 0.01,
    "heatmap_window": 3000,
    "igv_genome": "hg38",
    "outdir": "/path/to/cutntag_example/results"
}
```

### Step 5: Run the Pipeline

**Local execution** (for testing with small datasets):

```bash
# Clone CAWS
git clone https://github.com/ratan-lab/CAWS.git

# Run pipeline
snakemake \
  --snakefile CAWS/workflow/Snakefile \
  --configfile config.json \
  --cores 8 \
  --use-conda
```

**SLURM cluster execution** (recommended):

```bash
# Edit CAWS/profiles/slurm/config.yaml to set your account and partition

# Submit as batch job
cat > run_example.sh << 'EOF'
#!/bin/bash
#SBATCH -A your-allocation
#SBATCH -p standard
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --mem=4G
#SBATCH -o example_%j.out
#SBATCH -J caws_example

module load miniforge   # if conda is provided as a module on your cluster
conda activate caws-snakemake

snakemake \
  --snakefile CAWS/workflow/Snakefile \
  --profile CAWS/profiles/slurm \
  --configfile config.json
EOF

sbatch run_example.sh
```

### Step 6: Expected Results

See [Output Structure](#output-structure) for the full layout. **Key Outputs**:
1. **Peak files**: `results/peaks/seacr/H3K4me3_rep1.peaks.stringent.bed`, `results/peaks/macs3/H3K4me3_rep1_peaks.narrowPeak`
2. **Coverage tracks**: `results/tdf/H3K4me3_rep1.tdf` (for IGV) and `results/bigwig/H3K4me3_rep1.bw` (for other genome browsers and deepTools), both normalized as set by `normalization` (CPM by default). See [Normalization](#normalization)
3. **QC metrics**: `results/stats/stats.txt` (alignment, duplication, mtDNA, contaminants, *E. coli* spike-in fragments) and `results/stats/peaks_consolidated.txt` (peak counts, reproducibility, FRiP, widths, GC)
4. **Final report**: `results/reports/cutntag_analysis_report.html`

**What to expect for H3K4me3** (rough guides, not thresholds enforced by the pipeline):
- **Alignment rate**: typically >80% for good CUT&Tag libraries
- **FRiP**: high, as H3K4me3 signal is concentrated at promoters
- **TSS enrichment** (`stats/region_enrichment.txt`, region `TSS`): profile score well above 1, and above the IgG controls
- **Fragment size**: peak at ~150–200 bp (mono-nucleosomal)

### Step 7: Explore Results

Open the HTML report in a web browser:

```bash
# On local machine
open results/reports/cutntag_analysis_report.html

# Or copy from HPC to local
scp user@hpc:~/cutntag_example/results/reports/cutntag_analysis_report.html .
```

The interactive report includes:
- Per-sample peak summary (peak counts, FRiP, peak GC, TSS profile score)
- Alignment, duplication, mtDNA and contaminant rates
- Fragment size distributions (paired-end samples)
- Correlation of binned coverage between samples
- Peak counts, FRiP, widths and between-sample peak overlap for SEACR and MACS3

- Enrichment over TSSs and any other region sets

Peak heatmaps (`stats/seacr/*.heatmap.png`) and region profile plots and heatmaps (`stats/enrichment/`) are written as separate files.

---

## Configuration Details  

<details>
<summary><b>Click to expand configuration documentation</b></summary>

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
| `spikein_genome` | *Optional*, default `""`: name of the `contaminant_genomes` entry that serves as spike-in, usually `"Ecoli"` (carry-over DNA from pA-Tn5). When set, spike-in fragments are counted for every sample and reported in `stats/stats.txt`. See [Normalization](#normalization) |
| `normalization` | *Optional*, default `"cpm"`: scaling of the coverage tracks (TDF, bigWig) and region signal profiles. `"cpm"` scales to fragments per million filtered fragments. `"spikein"` scales to fragments per 10,000 spike-in fragments and needs `spikein_genome`. Peak calling is not affected. See [Normalization](#normalization) |
| `trim_adapters` | `true`/`false` — enable adapter trimming via Trim Galore. Also selects the Bowtie2 mode: `--local` with trimming, `--end-to-end` without |
| `igg_control` | `true`/`false` — use IgG controls as background for peak calling (SEACR control bedgraph, MACS3 `-c`). Requires exactly one `Control` per samplesheet `group` |
| `mt_chrom` | Name of mitochondrial chromosome (e.g., `"chrM"`) |
| `dedup` | *Optional*, default `"controls"`: remove PCR duplicates (Picard MarkDuplicates) from IgG controls only. `true` removes them from all samples, `false` from none. A mapping sets it per condition, e.g. `{"default": "controls", "BG4": true}`; the samplesheet `dedup` column overrides both per sample. See [Duplicate Removal](#duplicate-removal) |
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
| `enrichment_regions` | *Optional*, default `{}`: map of name → BED file (BED, BED.gz or narrowPeak) of regions to test for enrichment, e.g. `{"PQS": "/path/pqs.bed", "ATAC": "/path/atac_peaks.narrowPeak"}`. See [Region Enrichment](#region-enrichment) |
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

</details>

---

## Running on SLURM

<details>
<summary><b>Click to view SLURM execution details</b></summary>

### Snakemake 9.x (Profile-based execution)

**Step 1**: Configure the SLURM profile

Edit `profiles/slurm/config.yaml` to set your allocation and partition in the `default-resources` section:
```yaml
default-resources:
  - runtime=180m        # always give a unit: a bare number is read as seconds
  - mem_mb=9000
  - threads=1
  - slurm_account=your-allocation-name  # e.g., ratan
  - slurm_partition=your-partition      # e.g., standard
```

**Step 2**: Activate the Snakemake environment

```bash
module load miniforge   # if conda is provided as a module on your cluster
conda activate caws-snakemake
```

**Step 3**: Execute the pipeline

**Simple execution:**
```bash
snakemake --profile profiles/slurm --configfile config.json
```

**Advanced execution with conda caching:**
```bash
snakemake \
  --profile profiles/slurm \
  --configfile config.json \
  --conda-prefix /scratch/$USER/conda-envs \
  --jobs 100
```

**As a SLURM job:**
```bash
#!/bin/bash
#SBATCH -A your-allocation
#SBATCH -p standard
#SBATCH -c 1
#SBATCH -t 48:00:00
#SBATCH --mem=4G
#SBATCH -o caws_%j.out
#SBATCH -e caws_%j.err

module load miniforge   # if conda is provided as a module on your cluster
conda activate caws-snakemake

snakemake \
  --profile profiles/slurm \
  --configfile config.json \
  --conda-prefix /scratch/$USER/conda-envs
```

### Profile Configuration Details

The profile (`profiles/slurm/config.yaml`) handles:
- Job submission to SLURM with executor plugin
- Resource allocation per rule (memory, CPU, runtime)
- Conda environments (`use-conda`), `--rerun-incomplete` and `--keep-going`
- Concurrent job limits (default: 100, override with `--jobs`)

### Useful Options

- `--conda-prefix /path/to/envs`: Cache conda environments (recommended for repeated runs)
- `--jobs N`: Override maximum concurrent SLURM jobs (default: 100)
- `--notemp`: Keep temporary files for debugging
- `--rerun-incomplete`: Rerun incomplete jobs (already enabled in profile)
- `--dryrun` or `-n`: Preview execution without running
- `--printshellcmds`: Show shell commands (already enabled in profile)

### Running Several Analyses at Once

Give every analysis its own config with its own `outdir`, and launch each as a separate controller job. Snakemake keeps its state (lock, job metadata, SLURM logs) in `outdir/.snakemake`, and all CAWS results and temporary directories are written below `outdir`, so runs do not overwrite each other. If two runs point at the same `outdir` by mistake, the second stops with "Directory cannot be locked".

These can be shared safely between runs:
- The CAWS code, the SLURM profile, references, indexes and region files, which are only read.
- The conda environments in `--conda-prefix`. To avoid two runs building the same environment at the same time, create them once before launching several runs: `snakemake --profile profiles/slurm --configfile config.json --conda-create-envs-only`.
- `$TMPDIR`. Sorting and Picard write uniquely named files there. Several runs that sort large files at the same time share its space, so use a directory on scratch.

Each controller submits up to `jobs` (100) SLURM jobs, so several runs together can reach your account's job limits.

### Monitoring

Monitor jobs with:
```bash
squeue -u $USER           # Check job status
tail -f caws_*.out        # Watch controller output
tail -f logs/rulename/*   # Watch specific rule logs
```

</details>

---

## Features

### Data Input Flexibility
- **Single-End & Paired-End Support** – Process SE, PE, or mixed datasets seamlessly
- **Automatic Detection** – Pipeline automatically identifies sequencing type from samplesheet
- **SE/PE Compatibility Validation** – Ensures experimental samples and controls have matching sequencing types
- **Optimized Processing** – Separate optimized workflows for SE and PE data with appropriate biological handling

- **BAM Input and Multiple Runs** – Aligned BAM files can be used instead of FASTQ; several sequencing runs of a library are merged before duplicate removal
- **Input Validation** – Checks the samplesheet (files exist, consistent rows per sample, one control per group, matching read types), BAM references, region files and the contaminant and dedup configuration at startup; checks that GTF chromosome names match the reference when extracting TSSs

### Core Pipeline
- Adapter trimming (**Trim Galore**) and read QC (**FastQC**, on trimmed reads when trimming is enabled)
- Alignment (**Bowtie2**; paired-end in concordant-pair mode)
- Duplicate removal (**Picard**) from IgG controls by default, configurable for all samples or per sample; duplication rate and estimated library size (paired-end) are reported for every sample
- Filtering by mapping quality, chromosome list and blacklist
- Peak calling with both **SEACR** (stringent) and **MACS3** (narrow or broad peaks); IgG control used by both when enabled. Paired-end samples use actual fragments (`BAMPE`); single-end samples use a fixed 147 bp extension
- Consensus peaks per condition (ENCODE-style): peaks are called on the pooled replicates of each condition and kept if they overlap peaks in at least `consensus_min_replicates` replicates by at least 50% of either peak. Pooling uses fragments, with single-end reads extended to 147 bp, so conditions can mix single-end and paired-end replicates; IgG controls of the condition's groups are pooled the same way
- Contamination screening (**FastQ Screen**: *E. coli*, bovine, or any configured genome) and mitochondrial read fraction
- Coverage tracks normalized to counts per million (CPM, default) or to *E. coli* carry-over spike-in fragments (optional), with spike-in fragments counted for every sample when a spike-in genome is set

### QC Metrics and Visualization
- Fragment size distributions (paired-end)
- Coverage correlation between samples (Pearson, log2 counts in 500 bp bins; 50,000 bins sampled)
- FRiP, peak reproducibility between replicates, and pairwise peak overlap between samples (shared / total merged peak intervals)
- Peak width distributions (narrow/broad categories) and median peak GC content
- Enrichment over TSSs and user-supplied region sets (e.g. PQS, ATAC peaks): signal profile, reads in regions and peak overlap vs. shuffled peaks, each relative to the matched IgG
- Library complexity (ENCODE NRF, PBC1, PBC2) with PCR bottleneck grading
- Signal heatmaps around SEACR peak summits (**deepTools**)
- Normalized coverage tracks from the filtered alignments as TDF (IGV) and bigWig (fragment coverage for paired-end, read coverage for single-end), plus unfiltered, unnormalized bigWig tracks used for the peak heatmaps
- Interactive HTML dashboard (**flexdashboard**/**Plotly**) with exportable tables
- **MultiQC** report of FastQC, Trim Galore, Bowtie2, Picard and FastQ Screen results
- IGV session file that loads all coverage tracks and peak sets
- Record of the configuration, samplesheet and pinned software versions used for the run

---

## Output Structure

<details>
<summary><b>Click to expand output directory layout</b></summary>

```
output_directory/
├── qc/                  # FastQC reports; fastq_screen/ contaminant screens
├── trimmed/             # Adapter-trimmed reads
├── alignments/          # alignment statistics, complexity, spike-in counts, fragment lengths, mtDNA counts, unfiltered bigWigs; units/: Bowtie2 logs per run
├── dedupalignments/     # BAMs used downstream (deduplicated where dedup applies), Picard metrics, filtered BAMs
├── bedalignments/       # Raw-count bedgraphs (SEACR input) and 500 bp binned counts
├── peaks/
│   ├── seacr/           # SEACR peaks
│   ├── macs3/           # MACS3 peaks
│   └── consensus/       # pooled peak calls and consensus peaks per condition (seacr/, macs3/)
├── stats/               # stats.txt, samples.tsv, peaks_consolidated.txt, consensus_peaks.txt, region_enrichment.txt
│   ├── enrichment/      # per region set: profile plot, heatmap, matrix and the three enrichment measures
│   ├── seacr/           # SEACR metrics, peak GC, peak heatmaps
│   └── macs3/           # MACS3 metrics, peak GC
├── tdf/                 # IGV TDF coverage tracks, normalized (one per sample), igv_session.xml
├── bigwig/              # the same normalized coverage tracks as bigWig (one per sample)
├── reports/             # HTML report, MultiQC report, pipeline_info/ (config, samplesheet, software versions)
├── annotations/         # tss.bed; regions/: reference points and merged regions per region set
└── logs/                # Log files
```

**Note on alignments**: every sample has its BAM in `dedupalignments/`: deduplicated where `dedup` applies, otherwise an unchanged copy of the alignments. The raw `alignments/{sample}.sorted.bam` (and index) are temporary and removed once all steps that use them have finished; run with `--notemp` to keep them.

**Note on IGV sessions**: open `tdf/igv_session.xml` in IGV (File → Open Session) to load the TDF tracks of all samples (experimental samples grouped by condition, then controls), the consensus peaks and the per-sample peak calls. Paths are relative, so keep the output directory layout when copying it.

**Note on consensus peaks**: `peaks/consensus/{seacr,macs3}/{condition}.consensus.bed` lists, for each consensus peak, the number and names of supporting replicates. Conditions with a single replicate have no consensus peaks.

**Note on coverage tracks**: the TDF and bigWig tracks hold the same coverage of the filtered fragments, scaled to counts per million filtered fragments (CPM, excluding fragments on `mt_chrom` from the total) or, with `normalization: "spikein"`, to fragments per 10,000 spike-in fragments. See [Normalization](#normalization). The raw-count bedgraphs in `bedalignments/` are kept for SEACR. The unfiltered `alignments/{sample}_raw.bw` tracks are not normalized and are used for the peak heatmaps only.

</details>

---

## Interactive HTML Dashboard  

**Report:** `reports/cutntag_analysis_report.html`  

### Dashboard Pages  
- **Sample Details:** per-sample peak counts, FRiP and peak GC for SEACR and MACS3, plus TSS profile score (if a GTF is provided), and consensus peak counts per condition  
- **Alignment:** alignment statistics table (including whether each sample was deduplicated) and rates (alignment, duplication, mtDNA, contaminants)  
- **Fragment Length:** fragment size distributions (paired-end samples)  
- **Coverage Similarity:** correlation of binned coverage between all samples  
- **SEACR / MACS3:** peak counts, FRiP, width distributions and width categories by condition  
- **SEACR Peaks / MACS3 Peaks:** pairwise peak overlap heatmaps (shared / total merged peak intervals)  
- **Region Enrichment:** enrichment over TSSs and user region sets, with IgG ratios  

**Interactive Features:**  
- Plotly-based zoom/pan/hover  
- Sortable, searchable tables with export buttons  

A separate **MultiQC** report (`reports/multiqc_report.html`) collects the FastQC, Trim Galore, Bowtie2, Picard and FastQ Screen results in the standard MultiQC format.

---

## QC Metrics  

<details>
<summary><b>Click to expand the definitions of the reported metrics</b></summary>

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
| `Deduplicated` | Whether duplicates were removed for this sample (see [Duplicate Removal](#duplicate-removal)) |
| `NRF`, `PBC1`, `PBC2` | ENCODE library complexity before duplicate removal, after the MAPQ and chromosome filters: NRF = distinct fragments / fragments; PBC1 = positions with exactly one fragment / distinct positions; PBC2 = positions with one / positions with two fragments. A position is a fragment's start and end (paired-end) or a read's 5′ end and strand (single-end) |
| `PCRBottleneck` | ENCODE grading of PBC1: severe (< 0.5), moderate (0.5–0.8), mild (0.8–0.9), none (≥ 0.9) |
| `UnalignedFragNum` | Unaligned fragments in the alignments; `0` usually means a BAM that was filtered to aligned reads, in which case alignment rate and contamination cannot be judged |
| `SpikeinFragNum`, `SpikeinRate` | Only when `spikein_genome` is set: fragments that did not align to the host and align to the spike-in genome (MAPQ ≥ `minquality`, deduplicated when the sample is), and their percentage of `SequencingDepth`. See [Normalization](#normalization) |
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

Consensus peak statistics (`stats/consensus_peaks.txt`) are described under [Consensus Peaks](#consensus-peaks) and enrichment statistics (`stats/region_enrichment.txt`) under [Region Enrichment](#region-enrichment).
</details>

---

## Normalization

<details>
<summary><b>Click to expand details</b></summary>

CAWS scales the coverage tracks (`tdf/`, `bigwig/`) and the region signal profiles in one of two ways, set by `normalization`. Peak calling, FRiP, consensus peaks and the other statistics use the unscaled data and are the same in both modes.

**CPM (default).** Coverage is divided by the number of filtered fragments of the sample (excluding `mt_chrom`) and multiplied by one million. This corrects for sequencing depth. It assumes that every library contains about the same total amount of signal, so it shows *where* signal is, relative to the rest of the library. It cannot show a change in the *total* amount of signal. Libraries are amplified and sequenced to similar depths, so a sample in which the target is twice as abundant everywhere looks the same as the untreated sample after CPM.

**Spike-in.** CUT&Tag reactions contain carry-over *E. coli* DNA from the production of pA-Tn5. Its amount depends on the enzyme, not on the sample, so it can serve as a fixed reference. With `normalization: "spikein"`, coverage is multiplied by 10,000 / the number of spike-in fragments of the sample (fragments per 10,000 spike-in fragments, as in the CUT&Tag protocol of Zheng, Ahmad and Henikoff). This keeps global changes. For example, if a G4-stabilizing ligand doubles BG4 signal at all G4s, spike-in scaled tracks show the doubling, while CPM hides most of it and shows an apparent loss of background.

**When to use spike-in normalization:**
- Treatments or genotypes expected to change the total amount of the target, e.g. G4 ligands with BG4, inhibitors of a histone methyltransferase, or knockdown of the target protein.
- Only when the samples to be compared had the same number of cells and used the **same batch of pA-Tn5**. Carry-over differs between batches, and different cell numbers change the ratio of human to *E. coli* DNA for technical reasons.
- Only with enough spike-in fragments. Check `SpikeinFragNum` in `stats/stats.txt` first. Fewer than about 1,000 fragments give noisy scale factors, and CAWS writes a warning to `logs/count_spikein/` for such samples.

Use CPM for questions about *where* the target binds (peak landscapes, region enrichment) and for histone marks without expected global changes. With spike-in normalization, keep in mind that strong competition for pA-Tn5 could also change how efficiently the carry-over DNA is tagmented. An exogenous spike-in, such as a fixed amount of *Drosophila* chromatin added to each reaction, avoids this. Its genome can be given in `contaminant_genomes` and selected as `spikein_genome` in the same way.

**Configuration:**

```json
"contaminant_genomes": {"Ecoli": "/path/to/ecoli_index"},
"spikein_genome": "Ecoli",
"normalization": "spikein"
```

Setting `spikein_genome` with `"normalization": "cpm"` counts spike-in fragments for every sample (columns `SpikeinFragNum` and `SpikeinRate` in `stats/stats.txt`) without changing the tracks. This shows whether spike-in normalization would be usable before switching it on.

**How spike-in fragments are counted.** The reads of a sample that did not align to the host genome (pairs with neither mate aligned, for paired-end data) are aligned to the spike-in genome with the same Bowtie2 settings as the host. Fragments with MAPQ ≥ `minquality` are counted (paired-end: concordant pairs). Duplicates are removed from the count when the sample is deduplicated, so the spike-in is counted the same way as the host signal it scales. Because only reads that did not align to the host are used, sequences shared by the two genomes are not counted as spike-in. BAM input must contain the unaligned reads, otherwise no spike-in fragments are found, and spike-in normalization stops with an error for that sample.

**Outputs:**
- `bigwig/{sample}.bw` and `tdf/{sample}.tdf`: normalized coverage tracks (the scale factor used is in `logs/bam_to_bed_pe/` or `logs/bam_to_bed_se/`)
- `alignments/{sample}.spikein.txt`: spike-in fragment count
- `stats/stats.txt`: `SpikeinFragNum` and `SpikeinRate`
</details>

---

## Duplicate Removal  

<details>
<summary><b>Click to expand details</b></summary>

By default (`dedup: "controls"`), CAWS removes PCR duplicates from IgG controls only and keeps them in experimental samples.

**Why targets are not deduplicated by default.** CUT&Tag libraries are PCR-amplified, but pA-Tn5 is tethered at its target and inserts preferentially into accessible DNA nearby, so independent cells often produce fragments with identical ends. In high-quality libraries most apparent duplicates at target sites are therefore real signal, and removing them caps the signal exactly where it is strongest. This is especially true for single-end data, where a "duplicate" only needs the same start position and strand. IgG controls, in contrast, come from little DNA, have low complexity and are dominated by true PCR duplicates, which would otherwise inflate the background used for peak calling. Without unique molecular identifiers (UMIs) the two cannot be told apart, so this is a trade-off; nf-core/cutandrun uses the same default.

**When to deduplicate target samples.** The default suits abundant targets such as histone marks, which are spread over the genome, amplify like the background and do not run out. Low-abundance targets and DNA/RNA structures behave differently: few tethered pA-Tn5 complexes and few target molecules give low-complexity libraries that need more PCR cycles, so most duplicates are PCR copies; structured, GC-rich target fragments (e.g. G-quadruplexes) amplify less efficiently than background fragments, so PCR over-represents the background; and when targets are scarce, Tn5 cuts the same accessible sites repeatedly.

| Target | Suggested setting |
|---|---|
| Abundant histone marks (H3K27me3, H3K4me3, H3K36me3, …) | `"controls"` (default) |
| Transcription factors, other chromatin proteins | deduplicate (`true` for the condition) |
| DNA/RNA structures (BG4, S9.6, Z-DNA antibodies, …) | deduplicate (`true` for the condition) |
| Any low-input or heavily amplified library | deduplicate that sample (samplesheet `dedup` = `yes`) |

The library complexity metrics help decide: a PBC1 below 0.8 (moderate or severe PCR bottleneck), a low NRF, or a high `DuplicationRate` with a small `EstimatedLibrarySize` compared with other libraries of the same mark indicate a library dominated by PCR duplicates (see [QC Metrics](#qc-metrics)). For abundant marks, real Tn5 "duplicates" also lower these values somewhat, which is why CAWS reports them rather than deciding automatically.

**Configuration.** `dedup` in `config.json` is `"controls"`, `true` or `false`, or a mapping per condition, e.g.

```json
"dedup": {"default": "controls", "BG4": true, "S9.6": true}
```

The samplesheet `dedup` column (`yes`/`no`) overrides the condition setting for a sample, which overrides the default. CAWS warns at startup when replicates of one condition are treated differently, because their consensus peaks then pool deduplicated and non-deduplicated data.

**How it works.** The runs of a sample are merged first, so duplicates are found across all sequencing of a library. Picard MarkDuplicates runs on every sample, so duplication metrics are always reported. Where duplicates are removed, downstream steps use the deduplicated alignments; otherwise they use the alignments unchanged. Both are stored in `dedupalignments/`, and the `Deduplicated` column of `stats/stats.txt` records the choice for each sample. MACS3 is run with `--keep-dup all` (as in nf-core/cutandrun, ePeak and CUT&RUNTools), so it does not remove duplicates on its own and this setting alone decides how duplicates are handled.
</details>

---

## Region Enrichment

<details>
<summary><b>Click to expand details</b></summary>

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
</details>

---

## Consensus Peaks  

<details>
<summary><b>Click to expand details</b></summary>

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
</details>

---

## Peak Width Analysis  

<details>
<summary><b>Click to expand details</b></summary>

**Metrics:**  
- Summary statistics (mean, median, quartiles)  
- Peak categories:  
  - Narrow ≤ 500 bp (typical of point-source or promoter peaks)  
  - Broad > 2000 bp (typical of broad domains, e.g. H3K27me3)  
- Reported for both SEACR and MACS3  

**Visualizations:**  
- Histograms and category proportions  
- Sample-level statistics tables  
- Interactive Plotly charts  
</details>

---

### Key Dependencies  
- **Snakemake** — Workflow management  
- **Bowtie2** — Alignment  
- **SEACR** — CUT&RUN/CUT&Tag peak calling  
- **MACS3** — Model-based ChIP-seq analysis  
- **FastQC**, **Trim Galore**, **Picard**, **FastQ Screen**, **samtools**, **bedtools**, **deepTools**, **igvtools**, **bedGraphToBigWig** (UCSC)  
- **R packages** — `tidyverse`, `GenomicRanges`, `Rsamtools`, `plotly`, `flexdashboard`, `DT`  

---

## Troubleshooting

<details>
<summary><b>Click to view common issues and solutions</b></summary>

### Conda Version Error

**Error**: `CreateCondaEnvironmentException: Conda must be version 24.7.1 or later, found version X.X.X`

**Solution 1: Update Conda (Recommended)**
```bash
# Update conda in your base environment
conda update -n base conda

# Verify version
conda --version  # Should be ≥24.7.1
```

**Solution 2: Install Conda in Snakemake Environment**
```bash
# Activate your snakemake environment
conda activate caws-snakemake

# Install updated conda
conda install -c conda-forge conda>=24.7.1

# Verify
conda --version
```

### Module Load Errors

**Error**: `module: command not found` or modules don't load

**Solution**: Ensure you're on a system with environment modules. If running locally without modules:
```bash
# Skip module loading, just activate conda directly
conda activate your-snakemake-env
snakemake --cores 8 --use-conda --configfile config.json
```

### SLURM Job Failures

**Error**: Jobs fail with memory or time limits

**Solution**: Adjust resource allocations in `profiles/slurm/config.yaml`:
```yaml
set-resources:
  rule_name:
    mem_mb: 40000      # Increase memory
    runtime: "600m"    # Increase time (always give a unit, e.g. "600m" or "10h")
```

### Empty Peak Files

**Issue**: MACS3 produces empty peak files for some samples

**Cause**: No regions passed the threshold (common for low-depth or low-quality libraries), or the peak caller exited with an error

**Expected Behavior**: The pipeline continues and writes empty peak files, which are skipped in the statistics. Check the MACS3 log:
```bash
cat logs/call_macs_peaks_pe/sample_name.log   # or logs/call_macs_peaks_se/
```

If this affects many samples, consider a less stringent `macs3_qvalue_with_control` / `macs3_qvalue_no_control`. SEACR failures are also tolerated (empty `peaks/seacr/*.peaks.stringent.bed`); rerun the SEACR command shown by `snakemake -p` to see its error.

### Conda Environment Creation Hangs

**Issue**: Environment creation takes very long or appears stuck

**Solution**: Snakemake 9 uses conda's built-in libmamba solver. Make sure it is active:
```bash
conda config --show solver   # should report: libmamba
conda config --set solver libmamba
```

### Permission Denied Errors

**Issue**: Cannot write to output directory or conda-prefix

**Solution**: Ensure you have write permissions:
```bash
# Check output directory
ls -ld /path/to/outdir

# Use a directory in your scratch space for conda-prefix:
snakemake --conda-prefix /scratch/$USER/conda-envs --configfile config.json
```

### File Not Found Errors

**Issue**: Input files (FASTQ, reference, etc.) not found

**Solution**:
1. Check all paths in `config.json` are absolute paths
2. Verify files exist and are readable:
```bash
# Test file access
ls -lh /path/to/file
```
3. Ensure paths don't contain spaces or special characters

</details>

---

## Getting Help and Contributing

- **Questions and bug reports**: open an issue at https://github.com/ratan-lab/CAWS/issues. Please include your `config.json`, the samplesheet (paths can be redacted), the failing rule's log from `logs/`, and the Snakemake and conda versions.
- **Contributing**: see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## Use of AI Assistance

Since January 2026, development of CAWS has been assisted by generative AI tools, primarily Claude Code (Anthropic). AI assistance has been used for writing and modifying workflow code, scripts, and documentation. The architecture of the pipeline, the choice of tools and methods, and the overall design of the system are the result of decisions made by the human developers. All AI-assisted changes are reviewed by a human developer before being committed, and are checked by running the workflow (including dry runs and tests on small datasets) before release.

---

## License  

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.  

---

## Acknowledgements  

We thank the developers of all tools integrated into CAWS and the wider bioinformatics community for providing the foundational software that makes this workflow possible.  

