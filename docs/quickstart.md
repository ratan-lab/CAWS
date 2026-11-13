# Quick start

Complete the [installation](installation.md) first.

## Input data
- **FASTQ** files (gzipped) or **BAM** alignments — **single-end**, **paired-end** or both, and several sequencing runs per library
- Reference genome (**FASTA**, its `.fai` index, and a Bowtie2 index)
- **Chromosome list** (plain text, one chromosome name per line) — only reads on these chromosomes are used for peak calling and coverage
- **Blacklist** regions (**BED**)
- **Contaminant genomes** (Bowtie2 indexes, at least one, e.g. *E. coli*). *E. coli* also serves as the spike-in for optional spike-in normalization
- **Annotation** file (**GTF**, optional, for TSS enrichment) and other **region sets** (BED, optional) for enrichment analysis
- **Sample metadata** (**TSV**)

**Use absolute paths** in `config.json` and the samplesheet. The pipeline runs inside `outdir`, so relative paths are resolved from there, not from where you launch Snakemake.

## Prepare the samplesheet and config  

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

Edit `config.json` with your file paths and parameters. See [Configuration](configuration.md) for a complete description of all parameters.

## Run the analysis

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

## View the results  

Open `<outdir>/reports/cutntag_analysis_report.html` in a web browser.
