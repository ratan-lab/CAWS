#!/bin/bash
#SBATCH -J caws_testdata
#SBATCH -A your-allocation        # TODO: set your allocation
#SBATCH -p standard
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o caws_testdata_%j.out
#
# Build the CAWS test dataset from Kaya-Okur et al. 2019 (GEO GSE124557, K562).
#
# Reads are aligned genome-wide to hg38; the test dataset keeps every read pair
# whose fragment lies entirely inside one window, plus E. coli, chrM and
# unaligned pairs sampled at the same fraction (window pairs / genome-wide pairs),
# so QC rates stay close to those of the full library (with a floor of MIN_PAIRS
# pairs per category, so every QC path has data).
#
# Usage (activate the env from environment.yaml first; sbatch exports it):
#   export HG38_FA=/path/to/hg38.fa HG38_IDX=/path/to/hg38_bowtie2_prefix
#   sbatch build_test_data.sh select                         # download, align, rank windows
#   sbatch build_test_data.sh extract chr22:START-END        # write the test dataset
set -euo pipefail

MODE=${1:?usage: build_test_data.sh select | extract chrN:START-END}
HG38_FA=${HG38_FA:?set HG38_FA to an hg38 FASTA with UCSC names (chr22, chrM)}
HG38_IDX=${HG38_IDX:?set HG38_IDX to the matching Bowtie2 index prefix}
WORK=${WORK:-$PWD/caws_testdata_work}
OUT=${OUT:-$PWD/caws_testdata}
THREADS=${SLURM_CPUS_PER_TASK:-8}
SEED=42
MIN_PAIRS=200   # floor for sampled chrM / unaligned / E. coli pairs so every QC path has data
WINDOW_CHROM=chr22
WINDOW_BP=4000000
BIN_BP=500000

ECOLI_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
GTF_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
BLACKLIST_URL=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz

# test sample -> SRA runs (all paired-end, 2 x 25 bp)
declare -A RUNS=(
  [K27me3_rep1]="SRR8383510 SRR8383511"   # GSM3536515 K562_K27me3 (2018-09-19)
  [K27me3_rep2]="SRR8435047"              # GSM3560261 K562_H3K27me3 (2019-01-04)
  [IgG]="SRR8435051 SRR8435052"           # GSM3560264 K562_IgG (2018-09-19)
  [K27me3_SE]="SRR8754597"                # GSM3680219 K562_H3K27me3 fresh 1/1000; R1 used as single-end
)
SAMPLES="K27me3_rep1 K27me3_rep2 IgG K27me3_SE"

mkdir -p "$WORK"/{sra,fastq,bam,names,ref,tmp}

fetch_references() {
    [ -s "$WORK/ref/ecoli.1.bt2" ] && return
    wget -q -O "$WORK/ref/ecoli.fa.gz" "$ECOLI_URL"
    wget -q -O "$WORK/ref/gencode.gtf.gz" "$GTF_URL"
    wget -q -O "$WORK/ref/blacklist.bed.gz" "$BLACKLIST_URL"
    bowtie2-build -q --threads "$THREADS" "$WORK/ref/ecoli.fa.gz" "$WORK/ref/ecoli"
}

download_and_align() {
    local s=$1 r m url md5
    [ -f "$WORK/names/$s.ecoli.txt" ] && return   # written last, only on success
    # original FASTQs (full quality scores) from ENA, checked against ENA's MD5s
    for r in ${RUNS[$s]}; do
        wget -q -O "$WORK/sra/$r.tsv" \
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$r&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
        for m in 1 2; do
            url=$(tail -1 "$WORK/sra/$r.tsv" | cut -f2 | cut -d';' -f$m)
            md5=$(tail -1 "$WORK/sra/$r.tsv" | cut -f3 | cut -d';' -f$m)
            wget -q -O "$WORK/sra/${r}_$m.fastq.gz" "https://$url"
            echo "$md5  $WORK/sra/${r}_$m.fastq.gz" | md5sum -c --quiet -
        done
    done
    for m in 1 2; do
        for r in ${RUNS[$s]}; do cat "$WORK/sra/${r}_$m.fastq.gz"; done > "$WORK/fastq/${s}_$m.fq.gz"
    done
    rm -f "$WORK/sra/"*

    # same settings as the CAWS align_reads_pe rule
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -p "$THREADS" -x "$HG38_IDX" -1 "$WORK/fastq/${s}_1.fq.gz" -2 "$WORK/fastq/${s}_2.fq.gz" \
        2> "$WORK/bam/$s.hg38.log" \
        | samtools sort -@ 4 -o "$WORK/bam/$s.bam"
    samtools index "$WORK/bam/$s.bam"

    # concordant E. coli pairs (read names), same settings as the old CAWS E. coli check
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --no-unal --phred33 -I 10 -X 700 \
        -p "$THREADS" -x "$WORK/ref/ecoli" -1 "$WORK/fastq/${s}_1.fq.gz" -2 "$WORK/fastq/${s}_2.fq.gz" \
        2> "$WORK/bam/$s.ecoli.log" \
        | samtools view -f 0x42 - | cut -f1 > "$WORK/names/$s.ecoli.tmp"
    mv "$WORK/names/$s.ecoli.tmp" "$WORK/names/$s.ecoli.txt"
}

select_window() {
    fetch_references
    for s in $SAMPLES; do download_and_align "$s"; done

    # proper-pair fragments per bin on WINDOW_CHROM, and protein-coding TSSs per bin
    for s in $SAMPLES; do
        samtools view -f 0x42 -F 0x900 -q 10 "$WORK/bam/$s.bam" "$WINDOW_CHROM" \
            | awk -v b=$BIN_BP '{print int($4/b)}' | sort -n | uniq -c \
            | awk -v s="$s" '{print s"\t"$2"\t"$1}'
    done > "$WORK/bins.tsv"
    zcat "$WORK/ref/gencode.gtf.gz" \
        | awk -v c="$WINDOW_CHROM" -v b=$BIN_BP -F'\t' \
            '$1==c && $3=="gene" && $9~/gene_type "protein_coding"/ {print "genes\t"int(($7=="+"?$4:$5)/b)"\t1"}' \
        >> "$WORK/bins.tsv"

    python3 - "$WORK/bins.tsv" "$WINDOW_CHROM" "$WINDOW_BP" "$BIN_BP" <<'EOF'
import sys, collections
path, chrom, win, binbp = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
c = collections.defaultdict(lambda: collections.Counter())
for line in open(path):
    key, b, n = line.split("\t")
    c[key][int(b)] += int(n)
nb = win // binbp
last = max(max(v) for v in c.values())
rows = []
for start in range(0, last - nb + 2):
    tot = {k: sum(c[k][b] for b in range(start, start + nb)) for k in c}
    treat = tot.get("K27me3_rep1", 0) + tot.get("K27me3_rep2", 0)
    rows.append((treat, start, tot))
rows.sort(reverse=True)
print(f"Top {chrom} windows of {win/1e6:.1f} Mb by H3K27me3 fragments (rep1+rep2):")
print("region\trep1\trep2\tSE\tIgG\tprotein_coding_genes")
for treat, start, tot in rows[:10]:
    s, e = start * binbp, start * binbp + win
    print(f"{chrom}:{s}-{e}\t{tot.get('K27me3_rep1',0)}\t{tot.get('K27me3_rep2',0)}\t"
          f"{tot.get('K27me3_SE',0)}\t{tot.get('IgG',0)}\t{tot.get('genes',0)}")
EOF
}

sample_names() {  # subsample names on stdin at fraction $1, keeping at least MIN_PAIRS
    gawk -v f="$1" -v min=$MIN_PAIRS -v seed=$SEED \
        '{n[NR] = $0} END {if (NR * f < min) f = (NR ? min / NR : 0); srand(seed); for (i = 1; i <= NR; i++) if (rand() < f) print n[i]}'
}

extract_window() {
    local region=$1 chrom start end
    chrom=${region%%:*}; start=${region#*:}; start=${start%-*}; end=${region##*-}
    mkdir -p "$OUT"/{fastq,ref}
    : > "$WORK/summary.tsv"

    for s in $SAMPLES; do
        local bam="$WORK/bam/$s.bam" total win f
        # read pairs whose whole fragment lies inside the window
        samtools view -f 0x42 -F 0x900 "$bam" "$region" \
            | awk -v s="$start" -v e="$end" \
                '{l = ($4 < $8 ? $4 : $8); r = l + ($9 < 0 ? -$9 : $9); if (l > s && r < e) print $1}' \
            | sort -u > "$WORK/names/$s.window.txt"
        total=$(samtools view -c -f 0x42 -F 0x900 "$bam")
        win=$(wc -l < "$WORK/names/$s.window.txt")
        f=$(awk -v w="$win" -v t="$total" 'BEGIN{printf "%.8f", w/t}')

        {
            cat "$WORK/names/$s.window.txt"
            samtools view -f 0x42 -F 0x900 "$bam" chrM | cut -f1 | sample_names "$f"
            samtools view -f 77 "$bam" | cut -f1 | sample_names "$f"   # pairs unaligned to hg38
            sample_names "$f" < "$WORK/names/$s.ecoli.txt"
        } | sort -u > "$WORK/names/$s.keep.txt"

        local mates="1 2"
        [ "$s" = "K27me3_SE" ] && mates="1"
        for m in $mates; do
            seqkit grep -j "$THREADS" -f "$WORK/names/$s.keep.txt" "$WORK/fastq/${s}_$m.fq.gz" \
                > "$WORK/tmp/${s}_$m.fq"
        done
        if [ "$mates" = "1 2" ]; then  # mates must stay paired record-for-record (original IDs)
            cmp -s <(seqkit seq -n -i "$WORK/tmp/${s}_1.fq") <(seqkit seq -n -i "$WORK/tmp/${s}_2.fq") \
                || { echo "ERROR: R1/R2 out of sync for $s" >&2; exit 1; }
        fi
        for m in $mates; do
            seqkit replace -p '^.+$' -r "${s}.{nr}" "$WORK/tmp/${s}_$m.fq" \
                | gzip -9 > "$OUT/fastq/${s}_R$m.fq.gz"
            rm "$WORK/tmp/${s}_$m.fq"
        done
        printf "%s\t%s\t%s\t%s\t%s\n" "$s" "${RUNS[$s]}" "$total" "$win" \
            "$(wc -l < "$WORK/names/$s.keep.txt")" >> "$WORK/summary.tsv"
    done

    # host reference: chrom with everything outside the window masked to N (keeps hg38
    # coordinates, so GTF/blacklist/IGV need no liftover), plus chrM
    samtools faidx "$HG38_FA" "$chrom" > "$WORK/ref/chrom.fa"
    local len
    len=$(awk 'NR>1{n+=length($0)} END{print n}' "$WORK/ref/chrom.fa")
    printf "%s\t0\t%s\n" "$chrom" "$start" > "$WORK/ref/outside.bed"
    [ "$end" -lt "$len" ] && printf "%s\t%s\t%s\n" "$chrom" "$end" "$len" >> "$WORK/ref/outside.bed"
    bedtools maskfasta -fi "$WORK/ref/chrom.fa" -bed "$WORK/ref/outside.bed" -fo "$WORK/ref/masked.fa"
    { cat "$WORK/ref/masked.fa"; samtools faidx "$HG38_FA" chrM; } | gzip -9 > "$OUT/ref/genome.fa.gz"
    printf "%s\t%s\n" "$chrom" "$len" > "$OUT/ref/genome.chrom.sizes"
    samtools faidx "$HG38_FA" chrM | awk 'NR>1{n+=length($0)} END{print "chrM\t"n}' >> "$OUT/ref/genome.chrom.sizes"
    echo "$chrom" > "$OUT/ref/chromosomes.txt"

    cp "$WORK/ref/ecoli.fa.gz" "$OUT/ref/ecoli.fa.gz"
    zcat "$WORK/ref/gencode.gtf.gz" \
        | awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" '/^#/ || ($1==c && $4>s && $5<e)' \
        | gzip -9 > "$OUT/ref/genes.gtf.gz"
    zcat "$WORK/ref/blacklist.bed.gz" \
        | awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" '$1==c && $3>s && $2<e' \
        > "$OUT/ref/blacklist.bed"

    {
        echo "CAWS test dataset"
        echo "Source: Kaya-Okur et al. 2019, GEO GSE124557 (K562); built $(date +%F) with build_test_data.sh"
        echo "Window: $region (hg38); host reference = $chrom masked outside the window, plus chrM"
        echo "Kept: all proper pairs fully inside the window, plus chrM, hg38-unaligned and E. coli"
        echo "pairs sampled at fraction = window pairs / genome-wide proper pairs, with at least"
        echo "$MIN_PAIRS pairs per category (gawk srand($SEED))"
        echo
        printf "sample\truns\tgenome_wide_pairs\twindow_pairs\tpairs_kept\n"
        cat "$WORK/summary.tsv"
    } > "$OUT/PROVENANCE.txt"

    cat "$OUT/PROVENANCE.txt"
    du -sh "$OUT" "$OUT"/*
}

case "$MODE" in
    select)  select_window ;;
    extract) extract_window "${2:?extract needs a region, e.g. chr22:20000000-24000000}" ;;
    *)       echo "unknown mode: $MODE" >&2; exit 1 ;;
esac
