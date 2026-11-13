# One key per fragment from SAM records, for duplicate counting:
#   paired-end (-v paired=yes, read 1 of proper pairs): chromosome, leftmost position, end
#   single-end: chromosome, 5' end, strand
BEGIN { OFS = "\t" }
paired == "yes" { s = ($9 > 0) ? $4 : $8; l = ($9 < 0) ? -$9 : $9; print $3, s, s + l; next }
int($2 / 16) % 2 == 0 { print $3, $4, "+"; next }
{
    c = $6; len = 0
    while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
        n = substr(c, 1, RLENGTH - 1); op = substr(c, RLENGTH, 1)
        if (op ~ /[MDN=X]/) len += n
        c = substr(c, RLENGTH + 1)
    }
    print $3, $4 + len - 1, "-"
}
