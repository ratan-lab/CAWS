# CAWS: CUT&Tag Analysis Workflow System  
## *A Snakemake-based pipeline for comprehensive CUT&Tag data analysis*  

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Python](https://img.shields.io/badge/python-≥3.7-blue.svg)](https://python.org)
[![R](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://r-project.org)

---

## Overview  

**CAWS (CUT&Tag Analysis Workflow System)** is a comprehensive Snakemake pipeline for analyzing CUT&Tag (Cleavage Under Targets and Tagmentation) sequencing data. CUT&Tag is an advanced epigenomic profiling method that maps protein–DNA interactions and chromatin modifications with greater sensitivity and specificity than traditional ChIP-seq approaches.  

---

## Why Use CAWS?  

- **CUT&Tag Optimized** – Purpose-built for CUT&Tag data with specialized QC and peak-calling strategies  
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
- **Snakemake** ≥6.0.0  
- **Conda** or **Mamba**  
- **Python** ≥3.7  
- **R** ≥4.0  

### System Requirements  
- **Storage**: ~50 GB free space (typical dataset)  
- **Memory**: ≥8 GB RAM (16 GB+ recommended)  
- **CPU**: ≥4 cores  
- **OS**: Linux or macOS (cluster execution recommended)  

### Input Data Requirements  
- Paired-end **FASTQ** files (gzipped)  
- Reference genome (**FASTA** with Bowtie2 index)  
- **Annotation** file (**GTF**, optional for TSS enrichment)  
- **Blacklist** regions (**BED**)  
- **Sample metadata** (**TSV**)  

---

## Quick Start  

### 1. Clone and Set Up  

```bash
# Clone the repository
git clone https://github.com/your-username/CAWS.git
cd CAWS

# Create and activate Conda environment
conda env create -f environment.yml
conda activate caws
```

---

### 2. Prepare Configuration Files  

**Sample Sheet (`samplesheet.tsv`)**  

```tsv
sampleID	condition	group	read1	read2
sample1	Treatment	group1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
sample2	Control	group1	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz
```

**Configuration File (`config.json`)**  

```bash
cp config.example.json config.json
# Edit config.json with your file paths and parameters
```

---

### 3. Run the Analysis  

**Local Execution (small datasets)**  
```bash
snakemake --cores 8 --use-conda --configfile config.json
```

**SLURM Cluster Execution (recommended for large datasets)**  
```bash
snakemake \
    --cores 1 \
    --configfile config.json \
    --use-conda \
    --jobs 50 \
    --cluster-config clusterconfig.json \
    --cluster "sbatch --partition {cluster.queue} --cpus-per-task {cluster.nCPUs} --mem {cluster.memory}"
```

---

### 4. View Results  

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
<summary><b>Click to view SLURM command example</b></summary>

```bash
module load snakemake

snakemake \
    --cores 1 \
    --configfile config.json \
    --conda-frontend mamba \
    --latency-wait 60 \
    --notemp \
    --rerun-incomplete \
    --reason \
    --keep-going \
    --jobs 100 \
    --use-conda \
    --conda-prefix /data/CAWS/environment \
    --verbose \
    --cluster-config clusterconfig.json \
    --cluster "sbatch --partition {cluster.queue} -J {cluster.name} --cpus-per-task {cluster.nCPUs} --mem {cluster.memory} --time {cluster.maxTime} -o '{cluster.output}' -e '{cluster.error}' --mail-type=None --parsable -A {cluster.account}"
```
</details>

---

## Features  

### Core Pipeline  
- Quality control (**FastQC**)  
- Adapter trimming (**Trim Galore**)  
- Alignment (**Bowtie2**)  
- Duplicate removal (**Picard**)  
- Peak calling (**SEACR** & **MACS3**)  
- Contamination detection (mtDNA, *E. coli*)  
- Fragment size & correlation analysis  
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
├── alignments/          # BAM alignments
├── dedupalignments/     # Deduplicated BAMs
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

