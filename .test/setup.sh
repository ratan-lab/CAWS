#!/usr/bin/env bash
# Prepare a CAWS test run from the small dataset in .test/data (see .test/data/PROVENANCE.txt):
# decompress references, build Bowtie2 indexes, write samplesheet and config.
# Requires bowtie2 and samtools on PATH. Then run from the repository root:
#   snakemake -s workflow/Snakefile --configfile .test/run/config.json --use-conda --cores 4
set -euo pipefail

TEST=$(cd "$(dirname "$0")" && pwd)
DATA=$TEST/data
RUN=$TEST/run
mkdir -p "$RUN/ref" "$RUN/out"

gzip -dc "$DATA/ref/genome.fa.gz" > "$RUN/ref/genome.fa"
samtools faidx "$RUN/ref/genome.fa"
bowtie2-build -q --threads 2 "$RUN/ref/genome.fa" "$RUN/ref/genome"
bowtie2-build -q "$DATA/ref/ecoli.fa.gz" "$RUN/ref/ecoli"

F=$DATA/fastq
{
    printf "sampleID\tcondition\tgroup\tread1\tread2\n"
    printf "K27me3_rep1\tH3K27me3\tPE\t%s\t%s\n" "$F/K27me3_rep1_R1.fq.gz" "$F/K27me3_rep1_R2.fq.gz"
    printf "K27me3_rep2\tH3K27me3\tPE\t%s\t%s\n" "$F/K27me3_rep2_R1.fq.gz" "$F/K27me3_rep2_R2.fq.gz"
    printf "IgG_PE\tControl\tPE\t%s\t%s\n" "$F/IgG_R1.fq.gz" "$F/IgG_R2.fq.gz"
    printf "K27me3_SE\tH3K27me3\tSE\t%s\t-\n" "$F/K27me3_SE_R1.fq.gz"
    printf "IgG_SE\tControl\tSE\t%s\t-\n" "$F/IgG_R1.fq.gz"
} > "$RUN/samplesheet.tsv"

cat > "$RUN/config.json" <<EOF
{
    "samplesheet": "$RUN/samplesheet.tsv",
    "reference_fa": "$RUN/ref/genome.fa",
    "bt2_idx": "$RUN/ref/genome",
    "reference_fai": "$RUN/ref/genome.fa.fai",
    "chromosome_file": "$DATA/ref/chromosomes.txt",
    "blacklist": "$DATA/ref/blacklist.bed",
    "contaminant_genomes": {
        "Ecoli": "$RUN/ref/ecoli"
    },
    "contaminant_subset": 1000000,
    "trim_adapters": true,
    "igg_control": true,
    "mt_chrom": "chrM",
    "dedup": true,
    "minquality": 3,
    "fragment_min": 10,
    "fragment_max": 700,
    "genome_size": 4000000,
    "seacr_qvalue": 0.01,
    "macs3_qvalue_with_control": 0.01,
    "macs3_qvalue_no_control": 0.001,
    "heatmap_window": 3000,
    "gtf_file": "$DATA/ref/genes.gtf.gz",
    "igv_genome": "$DATA/ref/genome.chrom.sizes",
    "outdir": "$RUN/out"
}
EOF

echo "Test run prepared in $RUN"
