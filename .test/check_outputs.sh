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

exit $fail
