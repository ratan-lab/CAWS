# CAWS: CUT&Tag Analysis Workflow System  
## *A Snakemake-based pipeline for comprehensive CUT&Tag data analysis*  

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-≥3.7-blue.svg)](https://python.org)
[![R](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://r-project.org)

---

## Overview  

**CAWS (CUT&Tag Analysis Workflow System)** is a comprehensive Snakemake pipeline for analyzing CUT&Tag (Cleavage Under Targets and Tagmentation) sequencing data. CUT&Tag is an advanced epigenomic profiling method that maps protein–DNA interactions and chromatin modifications with greater sensitivity and specificity than traditional ChIP-seq approaches.  

---

## Why Use CAWS?

- **CUT&Tag Optimized** – Purpose-built for CUT&Tag data with specialized QC and peak-calling strategies
- **Flexible Input** – Supports single-end, paired-end, and mixed datasets in the same analysis
- **Comprehensive Analysis** – End-to-end workflow from raw reads to interactive visualizations
- **Dual Peak Calling** – Generates peak calls using both SEACR and MACS3 for robust peak detection
- **Interactive Reports** – HTML dashboards with Plotly-based visualizations
- **Reproducible** – Conda-based environment management for consistent results

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
- **Snakemake** ≥9.0.0 (tested with 9.8.1)
- **Conda** or **Mamba**
- **Python** ≥3.7
- **R** ≥4.0

**Note**: CAWS requires Snakemake 9.0 or later for profile-based SLURM execution with executor plugins.  

### System Requirements  
- **Storage**: ~50 GB free space (typical dataset)  
- **Memory**: ≥8 GB RAM (16 GB+ recommended)  
- **CPU**: ≥4 cores  
- **OS**: Linux or macOS (cluster execution recommended)  

### Input Data Requirements
- **FASTQ** files (gzipped) — supports both **single-end** and **paired-end** sequencing data
- Reference genome (**FASTA** with Bowtie2 index)
- **Annotation** file (**GTF**, optional for TSS enrichment)
- **Blacklist** regions (**BED**)
- **Sample metadata** (**TSV**)

---

## Quick Start  

### 1. Clone and Set Up

```bash
# Clone the repository
git clone https://github.com/ratan-lab/CAWS.git
cd CAWS
```

**Note**: CAWS uses conda environments managed automatically by Snakemake. No need to create a master environment.

---

### 2. Configure SLURM Profile (for HPC execution)

Edit `profiles/slurm/config.yaml` to set your SLURM allocation and partition:

```yaml
slurm:
  account: your-allocation-name  # e.g., ratan-lab
  partition: your-partition      # e.g., standard
  extra: "--parsable"
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

# Load required modules
module load miniforge
module load snakemake/9.8.1

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

# Load required modules
module load miniforge
module load snakemake/9.8.1

# Run pipeline with explicit parameters
snakemake \
  --profile /path/to/CAWS/profiles/slurm \
  --snakefile /path/to/CAWS/workflow/Snakefile \
  --configfile config.json \
  --cores 1 \
  --use-conda \
  --conda-frontend mamba \
  --conda-prefix /scratch/your-user/conda-envs \
  --latency-wait 60 \
  --rerun-incomplete \
  --keep-going \
  --jobs 100
```

**Interactive Execution:**
```bash
# Load modules
module load miniforge
module load snakemake/9.8.1

# Run in current session (use screen or tmux recommended)
snakemake --profile profiles/slurm --configfile config.json
```

**Notes:**
- `--conda-prefix` caches conda environments for faster subsequent runs
- `--jobs 100` overrides the profile setting (default: 50)
- `--notemp` can be added to preserve temporary files for debugging
- Use absolute paths if running from a different directory

---

### 5. View Results  

```bash
reports/cutntag_analysis_report.html
```

---

## Configuration Details  

<details>
<summary><b>Click to expand configuration documentation</b></summary>

| Parameter | Description |
|------------|-------------|
| `samplesheet` | Path to TSV file describing samples (columns: `sampleID`, `condition`, `group`, `read1`, `read2`) |
| `reference_fa` | Reference genome FASTA file |
| `bt2_idx` | Bowtie2 index prefix for the reference genome |
| `reference_fai` | FASTA index file |
| `chromosome_file` | List of chromosomes to include |
| `blacklist` | BED file of genomic regions to exclude |
| `ecoli_reference` | E. coli reference genome (for contamination assessment) |
| `ecoli_bt2_idx` | Bowtie2 index prefix for E. coli |
| `trim_adapters` | `true`/`false` — enable adapter trimming via Trim Galore |
| `igg_control` | `true`/`false` — use IgG controls as background for peak calling |
| `mt_chrom` | Name of mitochondrial chromosome (e.g., `"chrM"`) |
| `dedup` | `true`/`false` — remove PCR duplicates before peak calling |
| `minquality` | Minimum mapping quality threshold |
| `fragment_min`, `fragment_max` | Fragment size range for alignments |
| `genome_size` | Genome size for MACS3 peak calling: `"hs"` (human), `"mm"` (mouse), `"ce"` (C. elegans), `"dm"` (Drosophila), or numeric value (e.g., `2.7e9`) |
| `seacr_qvalue`, `macs3_qvalue_with_control`, `macs3_qvalue_no_control` | Statistical thresholds for SEACR/MACS3 |
| `heatmap_window` | Window size (bp) around peaks for heatmaps |
| `gtf_file` | Path to GTF file for TSS enrichment (optional) |
| `outdir` | Output directory for analysis results |

**Example `config.json`:**

```json
{
    "samplesheet": "samplesheet.tsv",
    "reference_fa": "/share/gatk_bundle/hg38/hs38DH.fa",
    "bt2_idx": "/share/gatk_bundle/hg38/hs38DH",
    "reference_fai": "/share/gatk_bundle/hg38/hs38DH.fa.fai",
    "chromosome_file": "/share/gatk_bundle/hg38/hs38DH.primary.txt",
    "blacklist": "/share/gatk_bundle/hg38/blacklist.bed.gz",
    "ecoli_reference": "/share/gatk_bundle/ecoli/MG1655.fa",
    "ecoli_bt2_idx": "/share/gatk_bundle/ecoli/MG1655",
    "trim_adapters": true,
    "igg_control": false,
    "mt_chrom": "chrM",
    "dedup": true,
    "minquality": 3,
    "fragment_min": 10,
    "fragment_max": 700,
    "genome_size": "hs",
    "seacr_qvalue": 0.01,
    "macs3_qvalue_with_control": 0.01,
    "macs3_qvalue_no_control": 0.001,
    "heatmap_window": 3000,
    "gtf_file": "",
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

Edit `profiles/slurm/config.yaml` to set your allocation and partition:
```yaml
slurm:
  account: your-allocation-name  # e.g., ratan-lab
  partition: your-partition      # e.g., standard
```

**Step 2**: Load required modules

```bash
module load miniforge
module load snakemake/9.8.1
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

module load miniforge
module load snakemake/9.8.1

snakemake \
  --profile profiles/slurm \
  --configfile config.json \
  --conda-prefix /scratch/$USER/conda-envs
```

### Profile Configuration Details

The profile (`profiles/slurm/config.yaml`) handles:
- Job submission to SLURM with executor plugin
- Resource allocation per rule (memory, CPU, runtime)
- Conda environment management with mamba
- Job monitoring and resubmission on failure
- Concurrent job limits (default: 50, override with `--jobs`)

### Useful Options

- `--conda-prefix /path/to/envs`: Cache conda environments (recommended for repeated runs)
- `--jobs N`: Override maximum concurrent SLURM jobs (default: 50)
- `--notemp`: Keep temporary files for debugging
- `--rerun-incomplete`: Rerun incomplete jobs (useful after interruptions)
- `--dryrun` or `-n`: Preview execution without running
- `--printshellcmds`: Show shell commands (already enabled in profile)

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

### Core Pipeline
- Quality control (**FastQC**)
- Adapter trimming (**Trim Galore**)
- Alignment (**Bowtie2**)
- Duplicate removal (**Picard**)
- Peak calling (**SEACR** & **MACS3**)
- Contamination detection (mtDNA, *E. coli*)
- Fragment size & correlation analysis (PE only)
- FRiP score and reproducibility metrics

### Enhanced Analysis (v2.0)  
- **TSS enrichment profiling**  
- **Peak width and shape distribution**  
- **Interactive HTML dashboards (Plotly)**  
- **SEACR vs MACS3 method comparison**  

---

## Output Structure

<details>
<summary><b>Click to expand output directory layout</b></summary>

```
output_directory/
├── qc/                  # FastQC reports
├── trimmed/             # Adapter-trimmed reads
├── alignments/          # BAM alignments (deleted on successful completion)
├── dedupalignments/     # Deduplicated BAMs (quality-filtered versions preserved)
├── bedalignments/       # BED files for peak calling
├── peaks/
│   ├── seacr/           # SEACR peaks
│   └── macs3/           # MACS3 peaks
├── stats/
│   ├── seacr/           # SEACR metrics
│   └── macs3/           # MACS3 metrics
├── reports/             # Interactive HTML reports
├── annotations/         # GTF-based TSS data
└── logs/                # Log files
```

**Note on alignment cleanup**: Upon successful pipeline completion, the `alignments/` directory is automatically deleted to save disk space. Quality-filtered BAM files used for peak calling are preserved in `dedupalignments/` (if deduplication is enabled). If you need to retain the original alignment files, interrupt the pipeline before completion or modify the cleanup behavior in the Snakefile (lines 965-966).

</details>

---

## Interactive HTML Dashboard  

**Report:** `reports/cutntag_analysis_report.html`  

### Dashboard Sections  
- **Overview:** Sample summaries, alignment metrics, and FRiP scores  
- **Peak Analysis:** Width distributions, peak counts, and statistics  
- **Quality Control:** Fragment size profiles, reproducibility metrics, TSS enrichment  
- **Method Comparison:** SEACR vs MACS3 correlation and overlap  

**Interactive Features:**  
- Plotly-based zoom/pan/hover  
- Sortable, searchable data tables  
- Conditional TSS panels (shown only if GTF is provided)  

---

## TSS Enrichment Analysis  

<details>
<summary><b>Click to expand details</b></summary>

Enable by specifying a GTF file in your configuration:

```json
"gtf_file": "/path/to/genes.gtf"
```

**Features:**  
- Extracts TSS sites from GTF (`.gtf` or `.gtf.gz`)  
- Calculates enrichment ±2 kb around TSS  
- Generates TSS enrichment plots and scores  
- Integrates results into the dashboard  

**Outputs:**  
- `annotations/tss.bed` — Extracted TSS coordinates  
- `stats/tss_enrichment.txt` — Per-sample enrichment scores  
</details>

---

## Peak Width Analysis  

<details>
<summary><b>Click to expand details</b></summary>

**Metrics:**  
- Summary statistics (mean, median, quartiles)  
- Peak categories:  
  - Narrow ≤ 500 bp (high specificity)  
  - Broad > 2000 bp (open chromatin regions)  
- Comparison between SEACR and MACS3  

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
- **FastQC**, **Trim Galore**, **Picard**, **deepTools**  
- **R packages** — `ggplot2`, `plotly`, `flexdashboard`  

---

## License  

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.  

---

## Acknowledgements  

We thank the developers of all tools integrated into CAWS and the wider bioinformatics community for providing the foundational software that makes this workflow possible.  

