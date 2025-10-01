# CUT&Tag Analysis With Snakemake

After cloning the repository, we need to fill in the name of the user account in `clusterconfig.json` that will be used to provide SUs for analysis. The `config.json` should be filled with the details for analysis. 

Here are the details of the various entries to be filled:
- samplesheet: Absolute path for a TSV file with the following columns
    - sampleID: A unique identifier for each sample
    - condition: A identifier shared by biological replicates
    - group: Currently unused, but will be used to specify controls in future
    - read1: Absolute path of the zipped fastq file of read1
    - read2: Absolute path of the zipped fastq file of read2
- reference_fa: Absolute path of the reference genome to be used for alignment
- bt2_idx: Absolute path for the bowtie2 index of the reference genome
- reference_fai: Absolute path for the fasta index of the reference genome
- chromosome_file: Absolute path of a file with chromosomes to include in the analysis. The remaining contigs and chromosomes are not used in peak calls
- blacklist: Absolute path of a file with regions that should be removed prior to peak calls. The blacklist for humans are available here : https://sites.google.com/site/anshulkundaje/projects/blacklists
- ecoli_reference: Absolute path of the E.Coli reference genome
- ecoli_bt2_idx: Absolute path for the bowtie2 index of the E.Coli genome
- trim_adapters: true/false. Should adapters be trimmed? We use trim_galore which identifies the adapters before trimming them
- igg_control: true/false. Should IgG controls be used for peak calling? When true, control samples (marked as "Control" in the condition column) are used as background for both SEACR and MACS3. Controls are matched to experiments via the 'group' field. The pipeline validates that each group has exactly one Control sample and at least one experimental sample when this option is enabled
- mt_chrom: The name of the mitochondrial genome in the reference. Reads aligning to the mtDNA are counted in stats, and removed prior to peak calls
- dedup: true/false. Should the PCR duplicates be removed prior to peak calls?
- minquality: What is the minimum mapping quality we should include in the counts for calling peaks
- fragment_min: Minimum fragment size for alignment (default: 10)
- fragment_max: Maximum fragment size for alignment (default: 700)
- seacr_qvalue: Statistical threshold for SEACR peak calling when no control is used (default: 0.01)
- macs3_qvalue_with_control: Q-value threshold for MACS3 peak calling with IgG control (default: 0.01)
- macs3_qvalue_no_control: Q-value threshold for MACS3 peak calling without control (default: 0.001)
- heatmap_window: Window size (bp) around peak centers for heatmap generation (default: 3000)
- gtf_file: Absolute path to GTF annotation file for TSS enrichment analysis (optional - leave empty to skip)
- outdir: Absolute path of the output directory of the analysis

Example config file

```
{
    "samplesheet": "samplesheet.tsv",
    "reference_fa": "/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg38dh/hs38DH.fa",
    "bt2_idx": "/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg38dh/hs38DH",
    "reference_fai": "/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg38dh/hs38DH.fa.fai",
    "chromosome_file": "/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg38dh/hs38DH.primary.txt", 
    "blacklist": "/share/gatk_bundle/ftp.broadinstitute.org/bundle/2.8/hg38dh/blacklist.bed.gz",
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
    "outdir": "/nv/vol169/cphg_ratan/ar7jq/CAWS/20230324"
}
```

### Running the analysis on SLURM cluster

```
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
    --conda-prefix /nv/vol169/cphg_ratan/ar7jq/CAWS/environment \
    --verbose \
    --cluster-config clusterconfig.json \
    --cluster " sbatch --partition {cluster.queue} -J {cluster.name} --cpus-per-task {cluster.nCPUs} --mem {cluster.memory} --time {cluster.maxTime} -o \"{cluster.output}\" -e \"{cluster.error}\" --mail-type=None --parsable -A {cluster.account} "
```

## Features

### Core Pipeline
- Quality control with FastQC
- Adapter trimming with trim_galore
- Alignment with Bowtie2
- Duplicate removal with Picard
- Peak calling with SEACR and MACS3
- Contamination detection (mtDNA, E.coli)
- Fragment size analysis
- Replicate correlation analysis
- Peak reproducibility assessment
- FRiP score calculation

### Enhanced Analysis & QC (v2.0)
- **TSS enrichment calculation** - Signal enrichment around transcription start sites
- **Peak width distribution analysis** - Detailed peak characteristics and categorization
- **Interactive HTML reports** - Comprehensive dashboard with plotly visualizations
- **Peak quality metrics** - Narrow vs broad peak classification
- **Method comparison** - Side-by-side SEACR vs MACS3 analysis

## Output Structure

```
output_directory/
├── qc/                   # FastQC quality control reports
├── trimmed/              # Adapter-trimmed reads (if enabled)
├── alignments/           # Alignment files and temporary data
├── dedupalignments/      # Deduplicated alignments (if enabled)
├── bedalignments/        # Processed BED files for peak calling
├── peaks/                # Peak calling results
│   ├── seacr/           # SEACR peak files
│   └── macs3/           # MACS3 peak files
├── stats/                # Summary statistics and plots
│   ├── seacr/           # SEACR-specific analysis
│   └── macs3/           # MACS3-specific analysis
├── reports/              # Interactive HTML reports
├── annotations/          # TSS annotations (if GTF provided)
└── logs/                 # Log files
```

## Key Output Files

### Summary Reports
- `reports/cutntag_analysis_report.html` - **Interactive HTML dashboard** with all analysis results
- `stats/stats.txt` - Alignment and QC statistics
- `stats/peaks_consolidated.txt` - Combined peak analysis results with width metrics

### Static Visualizations
- `stats/fragments.pdf` - Fragment size distribution
- `stats/replicates.pdf` - Replicate correlation analysis
- `stats/seacr/peaks_fig.pdf` - SEACR peak analysis plots (multi-page with width analysis)
- `stats/macs3/peaks_fig.pdf` - MACS3 peak analysis plots (multi-page with width analysis)
- `stats/seacr/peaks_width_analysis.pdf` - Dedicated SEACR peak width analysis
- `stats/macs3/peaks_width_analysis.pdf` - Dedicated MACS3 peak width analysis

### Enhanced Analysis Files
- `stats/tss_enrichment.txt` - TSS enrichment scores (if GTF provided)
- `stats/tss_enrichment.pdf` - TSS enrichment profile plots (if GTF provided)
- `stats/seacr/peaks_width_stats.txt` - Detailed SEACR peak width statistics
- `stats/macs3/peaks_width_stats.txt` - Detailed MACS3 peak width statistics
- `annotations/tss.bed` - Extracted transcription start sites (if GTF provided)

## Interactive HTML Report

The pipeline generates a comprehensive interactive HTML dashboard (`reports/cutntag_analysis_report.html`) featuring:

### Multi-page Dashboard
- **Overview**: Sample summaries, alignment stats, peak count/FRiP comparisons
- **Peak Analysis**: Interactive peak width distributions, statistics tables, categorization
- **Quality Control**: Reproducibility analysis, fragment distributions, TSS enrichment
- **Method Comparison**: SEACR vs MACS3 correlation plots and comparisons

### Interactive Features
- **Plotly integration**: Hover tooltips, zoom, pan functionality
- **Responsive tables**: Sortable, searchable data tables
- **Conditional content**: TSS analysis shown only when GTF file provided
- **Modern design**: Professional dashboard layout

## TSS Enrichment Analysis

### Configuration
To enable TSS enrichment analysis, provide a GTF file path in your config:

```json
{
    "gtf_file": "/path/to/your/genes.gtf"
}
```

### Analysis Features
- Extracts transcription start sites from GTF annotations
- Calculates signal enrichment around TSS (±2kb window)
- Generates TSS enrichment profile plots
- Computes enrichment scores (TSS signal vs flanking regions)
- Integrates results into HTML dashboard

### Output Files
- `annotations/tss.bed` - Extracted TSS coordinates
- `stats/tss_enrichment.txt` - Per-sample enrichment scores
- `stats/tss_enrichment.pdf` - TSS profile visualization

## Peak Width Analysis

### Metrics Calculated
- **Summary statistics**: median, mean, quartiles, min/max widths
- **Peak categorization**:
  - Narrow peaks: ≤500bp (indicate high specificity)
  - Broad peaks: >2000bp (indicate open chromatin regions)
- **Method comparison**: SEACR vs MACS3 peak characteristics

### Visualizations
- Peak width distribution histograms by condition
- Peak category percentages (narrow vs broad)
- Sample-level width statistics tables

## Citation

If you use this pipeline, please cite the relevant tools:

- Snakemake
- Bowtie2
- SEACR
- MACS3
- FastQC
- trim_galore
- Picard
- deepTools (for TSS analysis)
- R/ggplot2, plotly, flexdashboard (for visualizations)
