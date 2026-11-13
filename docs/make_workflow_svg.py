"""Generate the CAWS workflow figure (docs/workflow.svg).

Usage: python3 docs/make_workflow_svg.py docs/workflow.svg
For a PDF (e.g. for the manuscript): rsvg-convert -f pdf docs/workflow.svg -o workflow.pdf
"""
import sys

FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"

# Okabe-Ito: (stroke, card fill, panel tint, dark text)
GREY = ("#8C8C8C", "#F0F0F0", "#FAFAFA", "#4D4D4D")
BLUE = ("#0072B2", "#D6E9F5", "#F2F8FC", "#005A8C")
ORANGE = ("#E69F00", "#FCEBC4", None, "#8F6300")
PINK = ("#CC79A7", "#F5E1EC", None, "#8E4A70")
GREEN = ("#009E73", "#CDEFE3", "#EDF8F4", "#00704F")
DARK = ("#2F3640", "#2F3640", "#F4F5F6", "#2F3640")
NEUTRAL_PANEL = "#F7F7F7"

COLS = {1: 30, 2: 400, 3: 770, 4: 1140}  # column x; width 320
W = COLS[4] + 320 + 30
CW, CARD_PAD, GAP, TOP = 320, 15, 20, 112
CARD_W = CW - 2 * CARD_PAD

out = []
def add(s): out.append(s)

def esc(t):
    t = t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return t.replace("E. coli", '<tspan font-style="italic">E. coli</tspan>')

def card_height(lines): return 62 + 24 * len(lines)

def card(col, y, title, lines, pal, fill=None):
    x, h = COLS[col] + CARD_PAD, card_height(lines)
    add(f'<rect x="{x}" y="{y}" width="{CARD_W}" height="{h}" rx="12" fill="{fill or pal[1]}" stroke="{pal[0]}" stroke-width="2"/>')
    add(f'<text x="{x + 18}" y="{y + 35}" font-size="24" font-weight="bold" fill="#1A1A1A">{esc(title)}</text>')
    for i, ln in enumerate(lines):
        add(f'<text x="{x + 18}" y="{y + 63 + 24 * i}" font-size="20" fill="#3D3D3D">{esc(ln)}</text>')
    return x, y, CARD_W, h

def stack(col, cards, pal_default=None):
    """Cards of a column from the top, GAP apart; each entry is (title, lines[, palette])."""
    boxes, y = [], TOP
    for entry in cards:
        title, lines = entry[0], entry[1]
        pal = entry[2] if len(entry) > 2 else pal_default
        boxes.append(card(col, y, title, lines, pal))
        y += card_height(lines) + GAP
    return boxes

def header(col, text, pal):
    x = COLS[col]
    add(f'<text x="{x + 4}" y="58" font-size="27" font-weight="bold" fill="{pal[3]}">{esc(text)}</text>')
    add(f'<rect x="{x}" y="70" width="{CW}" height="5" rx="2.5" fill="{pal[0]}"/>')

def arrow(d, color="#4A4A4A", dashed=False, width=2.4):
    dash = ' stroke-dasharray="7 6"' if dashed else ""
    mk = "arrG" if color == GREEN[0] else "arr"
    add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}"{dash} '
        f'stroke-linejoin="round" marker-end="url(#{mk})"/>')

def line(d, color="#4A4A4A", width=2.4):
    add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round"/>')

mid = lambda b: b[1] + b[3] / 2
bottom = lambda b: b[1] + b[3]

# ---------------- content
INPUTS = [("FASTQ or BAM", ["SE, PE or mixed, several", "runs per library"]),
          ("Samplesheet", ["conditions, IgG controls,", "dedup per sample"]),
          ("Reference genome", ["Bowtie2 index, chromosome", "list, blacklist"]),
          ("Region sets", ["TSSs from a GTF, BED or", "narrowPeak files (optional)"]),
          ("Contaminant genomes", ["for screening and spike-in,", "e.g. E. coli, bovine"])]
PROCESSING = [("Input validation", ["samplesheet, BAMs, regions"]),
              ("Trim Galore, Bowtie2", ["FastQC, per sequencing run,", "BAM input imported"]),
              ("Merge runs", ["one library per sample"]),
              ("Duplicate removal", ["Picard, IgG by default,", "per condition or sample"]),
              ("Filtering", ["MAPQ, chromosome list,", "blacklist"])]
TRACKS_PEAKS = [("Coverage tracks", ["CPM or E. coli spike-in,", "bigWig and TDF"], PINK),
                ("SEACR", ["stringent, ± IgG control"], ORANGE),
                ("MACS3", ["narrow or broad, ± IgG,", "SE: 147 bp extension"], ORANGE),
                ("Consensus peaks", ["SEACR and MACS3 on", "pooled replicates, kept", "if replicates agree"], ORANGE)]
RESULTS = ["peaks, consensus peaks", "tracks, IGV session", "filtered BAMs, QC tables", "MultiQC, run record"]
PAGES = ["Sample details", "Alignment", "Fragment length", "Coverage similarity", "SEACR and MACS3",
         "Peak overlap", "Region enrichment"]
QC_ROW1 = [("Contaminant screen", ["FastQ Screen, spike-in", "fragment counts"]),
           ("Alignment QC", ["alignment, duplication,", "mtDNA, complexity (PBC),", "fragment sizes (PE)"]),
           ("Peak QC", ["FRiP, replicate", "reproducibility, width,", "GC, peak overlap"])]
QC_ROW2 = [("Sample correlation", ["binned coverage"]),
           ("Region enrichment", ["profile, reads and peaks", "in region sets, vs IgG"])]

# ---------------- layout from content
def column_end(cards): return TOP + sum(card_height(c[1]) for c in cards) + GAP * (len(cards) - 1)
panel_bottom = max(column_end(c) for c in (INPUTS, PROCESSING, TRACKS_PEAKS)) + 22
band_y = panel_bottom + 40
row1 = band_y + 20
row2 = row1 + max(card_height(c[1]) for c in QC_ROW1) + 14
band_bottom = row2 + max(card_height(c[1]) for c in QC_ROW2) + 20
H = band_bottom + 30

add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
    f'font-family="{FONT}">')
add('<defs>'
    '<marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5.5" markerHeight="5.5" orient="auto-start-reverse">'
    '<path d="M0,0 L10,5 L0,10 z" fill="#4A4A4A"/></marker>'
    '<marker id="arrG" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5.5" markerHeight="5.5" orient="auto-start-reverse">'
    f'<path d="M0,0 L10,5 L0,10 z" fill="{GREEN[0]}"/></marker>'
    '</defs>')
add(f'<rect width="{W}" height="{H}" fill="white"/>')

# ---------------- panels
ph = panel_bottom - 90
add(f'<rect x="{COLS[1]}" y="90" width="{CW}" height="{ph}" rx="16" fill="{GREY[2]}"/>')
add(f'<rect x="{COLS[2]}" y="90" width="{CW}" height="{ph}" rx="16" fill="{BLUE[2]}"/>')
add(f'<rect x="{COLS[3]}" y="90" width="{CW}" height="{ph}" rx="16" fill="{NEUTRAL_PANEL}"/>')
add(f'<rect x="{COLS[4]}" y="90" width="{CW}" height="{band_bottom - 90}" rx="16" fill="{DARK[2]}"/>')
add(f'<rect x="{COLS[1]}" y="{band_y}" width="{COLS[3] + CW - COLS[1]}" height="{band_bottom - band_y}" rx="16" '
    f'fill="{GREEN[2]}" stroke="{GREEN[0]}" stroke-width="2"/>')

header(1, "Inputs", GREY)
header(2, "Processing", BLUE)
header(3, "Tracks and peaks", (ORANGE[0], None, None, "#2F3640"))
header(4, "Outputs", DARK)
# two-color underline for the mixed column
add(f'<rect x="{COLS[3]}" y="70" width="{CW / 2}" height="5" rx="2.5" fill="{PINK[0]}"/>')
add(f'<rect x="{COLS[3] + CW / 2}" y="70" width="{CW / 2}" height="5" rx="2.5" fill="{ORANGE[0]}"/>')

# ---------------- columns
c1 = stack(1, INPUTS, GREY)
c2 = stack(2, PROCESSING, BLUE)
c3 = stack(3, TRACKS_PEAKS)

# processing chain
cx2 = COLS[2] + CW / 2
for b, nb in zip(c2, c2[1:]):
    arrow(f"M{cx2},{bottom(b)} V{nb[1] - 2}")

# inputs bracket -> validation
r1 = COLS[1] + CARD_PAD + CARD_W
bx = COLS[2] - 32
for b in c1:
    line(f"M{r1},{mid(b)} H{bx}")
line(f"M{bx},{mid(c1[0])} V{mid(c1[-1])}")
arrow(f"M{bx},{mid(c2[0])} H{c2[0][0] - 2}")

# filtering -> tracks, SEACR, MACS3 (bus in the gutter)
filt = c2[-1]
r2 = COLS[2] + CARD_PAD + CARD_W
l3 = COLS[3] + CARD_PAD
lane = r2 + 35
line(f"M{r2},{mid(filt)} H{lane}")
line(f"M{lane},{min(mid(c3[0]), mid(filt))} V{max(mid(c3[-1]), mid(filt))}")
for b in c3:  # tracks, SEACR, MACS3, and consensus peaks (called on the pooled filtered fragments)
    arrow(f"M{lane},{mid(b)} H{l3 - 2}")
cx3 = COLS[3] + CW / 2

# column 3 bracket -> result files
r3 = COLS[3] + CARD_PAD + CARD_W
bx3 = COLS[3] + CW + 22
for b in c3:
    line(f"M{r3},{mid(b)} H{bx3}")
line(f"M{bx3},{mid(c3[0])} V{mid(c3[-1])}")

# ---------------- column 4: outputs
res = card(4, TOP, "Result files", RESULTS, DARK, fill="white")
arrow(f"M{bx3},{mid(res)} H{res[0] - 2}")
rx = COLS[4] + CARD_PAD
rep_y = bottom(res) + 30
rep_h = band_bottom - 20 - rep_y
add(f'<rect x="{rx}" y="{rep_y}" width="{CARD_W}" height="{rep_h}" rx="12" fill="{DARK[1]}"/>')
add(f'<text x="{rx + 18}" y="{rep_y + 40}" font-size="25" font-weight="bold" fill="white">Interactive</text>')
add(f'<text x="{rx + 18}" y="{rep_y + 70}" font-size="25" font-weight="bold" fill="white">HTML report</text>')
add(f'<text x="{rx + 18}" y="{rep_y + 98}" font-size="19" fill="#C9CED6">flexdashboard, Plotly</text>')
step = (rep_h - 135) / len(PAGES)
for i, p in enumerate(PAGES):
    yy = rep_y + 160 + step * i
    add(f'<rect x="{rx + 18}" y="{yy - 28}" width="{CARD_W - 36}" height="44" rx="8" fill="#424B57"/>')
    add(f'<text x="{rx + 32}" y="{yy}" font-size="19" fill="white">{esc(p)}</text>')

# ---------------- QC band
qc1 = [card(i + 1, row1, t, l, GREEN) for i, (t, l) in enumerate(QC_ROW1)]
add(f'<text x="{COLS[1] + 25}" y="{row2 + 42}" font-size="28" font-weight="bold" fill="{GREEN[3]}">Quality control</text>')
add(f'<text x="{COLS[1] + 25}" y="{row2 + 72}" font-size="20" fill="{GREEN[3]}">metrics collected in</text>')
add(f'<text x="{COLS[1] + 25}" y="{row2 + 96}" font-size="20" fill="{GREEN[3]}">the HTML report</text>')
qc2 = [card(i + 2, row2, t, l, GREEN) for i, (t, l) in enumerate(QC_ROW2)]

# connectors into QC
cx1 = COLS[1] + CW / 2
arrow(f"M{cx1},{bottom(c1[-1])} V{row1 - 2}", color=GREEN[0], dashed=True)
arrow(f"M{cx2},{bottom(filt)} V{row1 - 2}", color=GREEN[0], dashed=True)
arrow(f"M{cx3},{bottom(c3[-1])} V{row1 - 2}", color=GREEN[0], dashed=True)
# QC band -> report
band_r = COLS[3] + CW
arrow(f"M{band_r},{(band_y + band_bottom) / 2} H{rx - 2}", color=GREEN[0], width=2.6)

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out))
