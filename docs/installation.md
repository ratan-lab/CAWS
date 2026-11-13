# Installation

## Required Software
- **Snakemake** ≥9.0.0 (CI tests the latest release; also used with 9.8.1)
- **Conda** ≥24.7.1
- **Python** ≥3.11 (required by Snakemake 9)

All other tools, including R and its packages, are installed automatically by Snakemake into per-rule conda environments (`workflow/envs/`). `envs/snakemake-runner.yaml` provides an environment with Snakemake, the SLURM executor plugin and pandas.

**Important Notes**:
- Snakemake 9.0+ requires Conda ≥24.7.1 for environment management
- If you encounter conda version errors, see [Troubleshooting](troubleshooting.md)  

## System Requirements  
- **Storage**: ~50 GB free space (typical dataset)  
- **Memory**: ≥8 GB RAM (16 GB+ recommended)  
- **CPU**: ≥4 cores  
- **OS**: Linux (x86-64); HPC cluster execution recommended. macOS and Windows are not supported because several conda environments pin Linux (linux-64) builds.

## Clone and set up

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

## Configure the SLURM profile (for HPC execution)

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

## Testing the installation

The repository includes a small test dataset in `.test/data` (about 7 MB): H3K27me3 and IgG CUT&Tag reads from K562 cells (Kaya-Okur et al., 2019), subset to a 4 Mb window of chr22, with paired-end and single-end samples. See `.test/data/PROVENANCE.txt` for how it was made. To run it:

```bash
conda create -n caws-setup -c conda-forge -c bioconda bowtie2=2.5.1 samtools

# Build indexes and write the test samplesheet and config into .test/run/
conda run -n caws-setup bash .test/setup.sh

conda activate caws-snakemake
snakemake -s workflow/Snakefile --configfile .test/run/config.json --use-conda --cores 4
bash .test/check_outputs.sh
```

Results are written to `.test/run/out/`. The first run takes longer because Snakemake creates the conda environments for each tool. GitHub Actions (`.github/workflows/test.yml`) runs this test on every push to `main` and on pull requests, together with dry runs of all `igg_control`/`dedup` combinations.
