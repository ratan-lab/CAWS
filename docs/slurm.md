# Running on SLURM

## Snakemake 9.x (Profile-based execution)

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

## Profile Configuration Details

The profile (`profiles/slurm/config.yaml`) handles:
- Job submission to SLURM with executor plugin
- Resource allocation per rule (memory, CPU, runtime)
- Conda environments (`use-conda`), `--rerun-incomplete` and `--keep-going`
- Concurrent job limits (default: 100, override with `--jobs`)

## Useful Options

- `--conda-prefix /path/to/envs`: Cache conda environments (recommended for repeated runs)
- `--jobs N`: Override maximum concurrent SLURM jobs (default: 100)
- `--notemp`: Keep temporary files for debugging
- `--rerun-incomplete`: Rerun incomplete jobs (already enabled in profile)
- `--dryrun` or `-n`: Preview execution without running
- `--printshellcmds`: Show shell commands (already enabled in profile)

## Running Several Analyses at Once

Give every analysis its own config with its own `outdir`, and launch each as a separate controller job. Snakemake keeps its state (lock, job metadata, SLURM logs) in `outdir/.snakemake`, and all CAWS results and temporary directories are written below `outdir`, so runs do not overwrite each other. If two runs point at the same `outdir` by mistake, the second stops with "Directory cannot be locked".

These can be shared safely between runs:
- The CAWS code, the SLURM profile, references, indexes and region files, which are only read.
- The conda environments in `--conda-prefix`. To avoid two runs building the same environment at the same time, create them once before launching several runs: `snakemake --profile profiles/slurm --configfile config.json --conda-create-envs-only`.
- `$TMPDIR`. Sorting and Picard write uniquely named files there. Several runs that sort large files at the same time share its space, so use a directory on scratch.

Each controller submits up to `jobs` (100) SLURM jobs, so several runs together can reach your account's job limits.

## Monitoring

Monitor jobs with:
```bash
squeue -u $USER           # Check job status
tail -f caws_*.out        # Watch controller output
tail -f logs/rulename/*   # Watch specific rule logs
```
