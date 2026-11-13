# Troubleshooting

## Conda Version Error

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

## Module Load Errors

**Error**: `module: command not found` or modules don't load

**Solution**: Ensure you're on a system with environment modules. If running locally without modules:
```bash
# Skip module loading, just activate conda directly
conda activate your-snakemake-env
snakemake --cores 8 --use-conda --configfile config.json
```

## SLURM Job Failures

**Error**: Jobs fail with memory or time limits

**Solution**: Adjust resource allocations in `profiles/slurm/config.yaml`:
```yaml
set-resources:
  rule_name:
    mem_mb: 40000      # Increase memory
    runtime: "600m"    # Increase time (always give a unit, e.g. "600m" or "10h")
```

## Empty Peak Files

**Issue**: MACS3 produces empty peak files for some samples

**Cause**: No regions passed the threshold (common for low-depth or low-quality libraries), or the peak caller exited with an error

**Expected Behavior**: The pipeline continues and writes empty peak files, which are skipped in the statistics. Check the MACS3 log:
```bash
cat logs/call_macs_peaks_pe/sample_name.log   # or logs/call_macs_peaks_se/
```

If this affects many samples, consider a less stringent `macs3_qvalue_with_control` / `macs3_qvalue_no_control`. SEACR failures are also tolerated (empty `peaks/seacr/*.peaks.stringent.bed`); rerun the SEACR command shown by `snakemake -p` to see its error.

## Conda Environment Creation Hangs

**Issue**: Environment creation takes very long or appears stuck

**Solution**: Snakemake 9 uses conda's built-in libmamba solver. Make sure it is active:
```bash
conda config --show solver   # should report: libmamba
conda config --set solver libmamba
```

## Permission Denied Errors

**Issue**: Cannot write to output directory or conda-prefix

**Solution**: Ensure you have write permissions:
```bash
# Check output directory
ls -ld /path/to/outdir

# Use a directory in your scratch space for conda-prefix:
snakemake --conda-prefix /scratch/$USER/conda-envs --configfile config.json
```

## File Not Found Errors

**Issue**: Input files (FASTQ, reference, etc.) not found

**Solution**:
1. Check all paths in `config.json` are absolute paths
2. Verify files exist and are readable:
```bash
# Test file access
ls -lh /path/to/file
```
3. Ensure paths don't contain spaces or special characters
