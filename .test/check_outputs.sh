#!/usr/bin/env bash
# Sanity checks on the output of the .test run (see .test/setup.sh).
# Thresholds are loose on purpose: tool versions shift exact numbers.
set -euo pipefail

OUT=$(cd "$(dirname "$0")" && pwd)/run/out
fail=0
check() { if eval "$2"; then echo "ok   - $1"; else echo "FAIL - $1"; fail=1; fi; }

check "stats.txt has 5 samples" '[ $(tail -n +2 "$OUT/stats/stats.txt" | wc -l) -eq 5 ]'
check "alignment rate > 70% for every sample" \
    'tail -n +2 "$OUT/stats/stats.txt" | awk -F"\t" "{sub(/%/, \"\", \$5); if (\$5 + 0 <= 70) exit 1}"'
check "peaks_consolidated.txt has one row per H3K27me3 sample" \
    '[ $(tail -n +2 "$OUT/stats/peaks_consolidated.txt" | wc -l) -eq 3 ]'
for s in K27me3_rep1 K27me3_rep2 K27me3_SE; do
    check "$s SEACR peaks >= 10" '[ $(wc -l < "$OUT/peaks/seacr/$s.peaks.stringent.bed") -ge 10 ]'
done
for s in K27me3_rep1 K27me3_rep2; do
    check "$s MACS3 peaks >= 10" '[ $(wc -l < "$OUT/peaks/macs3/${s}_peaks.narrowPeak") -ge 10 ]'
done
check "K27me3_SE MACS3 peaks >= 1" '[ -s "$OUT/peaks/macs3/K27me3_SE_peaks.narrowPeak" ]'
check "TSS enrichment for 5 samples" '[ $(tail -n +2 "$OUT/stats/tss_enrichment.txt" | wc -l) -eq 5 ]'
for s in K27me3_rep1 K27me3_rep2 IgG_PE K27me3_SE IgG_SE; do
    check "$s TDF exists" '[ -s "$OUT/tdf/$s.tdf" ]'
done
check "HTML report exists" '[ -s "$OUT/reports/cutntag_analysis_report.html" ]'
check "deduplication follows dedup setting and samplesheet override" \
    '[ "$(tail -n +2 "$OUT/stats/stats.txt" | cut -f 1,11 | LC_ALL=C sort | tr "\t\n" ":," )" = "IgG_PE:yes,IgG_SE:yes,K27me3_SE:no,K27me3_rep1:no,K27me3_rep2:yes," ]'
for caller in seacr macs3; do
    check "$caller consensus peaks >= 1" '[ -s "$OUT/peaks/consensus/$caller/H3K27me3.consensus.bed" ]'
done
check "consensus summary has seacr and macs3 rows" '[ $(tail -n +2 "$OUT/stats/consensus_peaks.txt" | wc -l) -eq 2 ]'
check "MultiQC report exists" '[ -s "$OUT/reports/multiqc_report.html" ]'
check "software versions recorded" 'grep -q "^bowtie" "$OUT/reports/pipeline_info/software_versions.tsv"'
check "IGV session lists all TDF tracks" '[ $(grep -c "\.tdf\"" "$OUT/tdf/igv_session.xml") -eq 5 ]'

exit $fail
