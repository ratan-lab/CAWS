#!/usr/bin/env bash
# Prepare a CAWS test run from the small dataset in .test/data (see .test/data/PROVENANCE.txt):
# decompress references, build Bowtie2 indexes, write samplesheet and config.
# Requires bowtie2 and samtools on PATH. Then run from the repository root:
#   snakemake -s workflow/Snakefile --configfile .test/run/config.json --use-conda --cores 4
set -euo pipefail

TEST=$(cd "$(dirname "$0")" && pwd)
DATA=$TEST/data
RUN=$TEST/run
mkdir -p "$RUN/ref" "$RUN/out" "$RUN/fastq" "$RUN/bam"

gzip -dc "$DATA/ref/genome.fa.gz" > "$RUN/ref/genome.fa"
samtools faidx "$RUN/ref/genome.fa"
bowtie2-build -q --threads 2 "$RUN/ref/genome.fa" "$RUN/ref/genome"
bowtie2-build -q "$DATA/ref/ecoli.fa.gz" "$RUN/ref/ecoli"

F=$DATA/fastq

# K27me3_rep1 as two sequencing runs of one library (technical replicates, merged by CAWS)
lines=$(gzip -dc "$F/K27me3_rep1_R1.fq.gz" | wc -l)
half=$(( lines / 8 * 4 ))
for r in R1 R2; do
    gzip -dc "$F/K27me3_rep1_$r.fq.gz" | awk -v n="$half" 'NR <= n' | gzip > "$RUN/fastq/K27me3_rep1_runA_$r.fq.gz"
    gzip -dc "$F/K27me3_rep1_$r.fq.gz" | tail -n +"$(( half + 1 ))" | gzip > "$RUN/fastq/K27me3_rep1_runB_$r.fq.gz"
done

# K27me3_rep2 as BAM input, aligned with different settings than CAWS uses (unsorted)
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -X 700 -p 2 -x "$RUN/ref/genome" \
    -1 "$F/K27me3_rep2_R1.fq.gz" -2 "$F/K27me3_rep2_R2.fq.gz" 2> "$RUN/bam/K27me3_rep2.bowtie2.txt" \
    | samtools view -b -o "$RUN/bam/K27me3_rep2.bam" -

# a user-supplied region set for enrichment: gene bodies from the test GTF
gzip -dc "$DATA/ref/genes.gtf.gz" | awk -F'\t' -v OFS='\t' '$3 == "gene" {
    match($9, /gene_id "[^"]+"/); id = substr($9, RSTART + 9, RLENGTH - 10)
    print $1, $4 - 1, $5, id, 0, $7}' > "$RUN/ref/gene_bodies.bed"

{
    # K27me3_rep1: two runs (technical replicates); K27me3_rep2: BAM input
    # dedup column: K27me3_rep2 overrides the default ("controls") to exercise per-sample deduplication
    printf "sampleID\tcondition\tgroup\tread1\tread2\tdedup\n"
    printf "K27me3_rep1\tH3K27me3\tPE\t%s\t%s\t\n" "$RUN/fastq/K27me3_rep1_runA_R1.fq.gz" "$RUN/fastq/K27me3_rep1_runA_R2.fq.gz"
    printf "K27me3_rep1\tH3K27me3\tPE\t%s\t%s\t\n" "$RUN/fastq/K27me3_rep1_runB_R1.fq.gz" "$RUN/fastq/K27me3_rep1_runB_R2.fq.gz"
    printf "K27me3_rep2\tH3K27me3\tPE\t%s\t\tyes\n" "$RUN/bam/K27me3_rep2.bam"
    printf "IgG_PE\tControl\tPE\t%s\t%s\t\n" "$F/IgG_R1.fq.gz" "$F/IgG_R2.fq.gz"
    printf "K27me3_SE\tH3K27me3\tSE\t%s\t-\t\n" "$F/K27me3_SE_R1.fq.gz"
    printf "IgG_SE\tControl\tSE\t%s\t-\t\n" "$F/IgG_R1.fq.gz"
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
    "dedup": "controls",
    "minquality": 3,
    "fragment_min": 10,
    "fragment_max": 700,
    "genome_size": 4000000,
    "seacr_qvalue": 0.01,
    "macs3_qvalue_with_control": 0.01,
    "macs3_qvalue_no_control": 0.001,
    "macs3_broad": false,
    "macs3_broad_cutoff": 0.1,
    "consensus_min_replicates": 2,
    "enrichment_regions": {
        "GeneBodies": "$RUN/ref/gene_bodies.bed"
    },
    "heatmap_window": 3000,
    "gtf_file": "$DATA/ref/genes.gtf.gz",
    "igv_genome": "$DATA/ref/genome.chrom.sizes",
    "outdir": "$RUN/out"
}
EOF

echo "Test run prepared in $RUN"
