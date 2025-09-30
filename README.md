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
