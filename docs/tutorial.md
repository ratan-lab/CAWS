# Tutorial: H3K4me3 in K562 cells

This example demonstrates a complete CAWS analysis using publicly available CUT&Tag data. We'll analyze H3K4me3 (active promoter mark) samples from Kaya-Okur et al. (2019), *Nature Communications* 10:1930.

## Dataset Information

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

## Step 1: Download Example Data

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

## Step 2: Prepare Reference Files

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

## Step 3: Create Sample Sheet

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

## Step 4: Create Configuration File

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

## Step 5: Run the Pipeline

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

## Step 6: Expected Results

See [Output Structure](outputs.md#output-structure) for the full layout. **Key Outputs**:
1. **Peak files**: `results/peaks/seacr/H3K4me3_rep1.peaks.stringent.bed`, `results/peaks/macs3/H3K4me3_rep1_peaks.narrowPeak`
2. **Coverage tracks**: `results/tdf/H3K4me3_rep1.tdf` (for IGV) and `results/bigwig/H3K4me3_rep1.bw` (for other genome browsers and deepTools), both normalized as set by `normalization` (CPM by default). See [Normalization](normalization.md)
3. **QC metrics**: `results/stats/stats.txt` (alignment, duplication, mtDNA, contaminants, *E. coli* spike-in fragments) and `results/stats/peaks_consolidated.txt` (peak counts, reproducibility, FRiP, widths, GC)
4. **Final report**: `results/reports/cutntag_analysis_report.html`

**What to expect for H3K4me3** (rough guides, not thresholds enforced by the pipeline):
- **Alignment rate**: typically >80% for good CUT&Tag libraries
- **FRiP**: high, as H3K4me3 signal is concentrated at promoters
- **TSS enrichment** (`stats/region_enrichment.txt`, region `TSS`): profile score well above 1, and above the IgG controls
- **Fragment size**: peak at ~150–200 bp (mono-nucleosomal)

## Step 7: Explore Results

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
