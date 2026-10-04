#!/usr/bin/env python3
"""Generate the nf-core-style metro map for hybrid-metatranscriptomics as a standalone SVG.

    python3 extra/make_metro_map.py [assets/metro_map.svg]

Coordinates are hand-placed: after adding or renaming a step, render it (e.g. rsvg-convert) and
check that no label collides.
"""
import sys
from html import escape

W, H = 1800, 930
FONT = "Helvetica, Arial, sans-serif"
MONO = "Menlo, Consolas, 'DejaVu Sans Mono', monospace"

C = dict(
    long="#e6550d",      # long reads
    short="#2171b5",     # short reads
    asm="#6baed6",       # assembly-only reads
    ref="#737373",       # genome / reference databases
    tx="#756bb1",        # transcriptome / reference ladder
    ann="#31a354",       # annotation
    de="#d6266b",        # quantification & DE
    gmap="#17a2b8",      # gene maps
    dtu="#b07aa1",       # isoform switching
    rep="#8c6d31",       # reporting
)
INK = "#222222"
SUB = "#555555"
LW = 8

out = []
def add(s): out.append(s)

# ------------------------------------------------------------------ primitives
def band(x0, x1, y0, y1, n, title, shade, right=False):
    add(f'<rect x="{x0}" y="{y0}" width="{x1-x0}" height="{y1-y0}" rx="10" fill="{shade}"/>')
    tx, ta = (right, "end") if right else (x0 + 12, "start")
    add(f'<text x="{tx}" y="{y0+24}" text-anchor="{ta}" font-family="{FONT}" font-size="15" font-weight="700" fill="#444">'
        f'<tspan fill="#999">{n}</tspan>  {escape(title)}</text>')

def line(pts, color, dash=None, width=LW):
    d = "M " + " L ".join(f"{x},{y}" for x, y in pts)
    da = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" '
        f'stroke-linejoin="round" stroke-linecap="round"{da}/>')

def station(x, y, optional=False):
    da = ' stroke-dasharray="4 3"' if optional else ""
    add(f'<circle cx="{x}" cy="{y}" r="8" fill="#fff" stroke="{INK}" stroke-width="3"{da}/>')

def hub(x, y0, y1, optional=False):
    """interchange spanning several parallel tracks"""
    da = ' stroke-dasharray="4 3"' if optional else ""
    add(f'<rect x="{x-11}" y="{y0-11}" width="22" height="{y1-y0+22}" rx="11" '
        f'fill="#fff" stroke="{INK}" stroke-width="3"{da}/>')

def terminus(x, y, color):
    add(f'<rect x="{x-10}" y="{y-10}" width="20" height="20" rx="4" fill="{color}" stroke="{INK}" stroke-width="2.5"/>')

def pill_width(t):
    return len(t) * 6.3 + 10

def label(x, y, title, sub=None, flags=(), anchor="middle"):
    """title at baseline y, then sub, then option pills"""
    add(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="{FONT}" font-size="13" '
        f'font-weight="700" fill="{INK}">{escape(title)}</text>')
    yy = y
    if sub:
        yy += 14
        add(f'<text x="{x}" y="{yy}" text-anchor="{anchor}" font-family="{FONT}" font-size="11" '
            f'fill="{SUB}">{escape(sub)}</text>')
    if flags:
        yy += 6
        ws = [pill_width(f) for f in flags]
        total = sum(ws) + 4 * (len(flags) - 1)
        x0 = {"middle": x - total / 2, "start": x, "end": x - total}[anchor]
        for f, w in zip(flags, ws):
            add(f'<rect x="{x0:.1f}" y="{yy}" width="{w:.1f}" height="15" rx="7.5" fill="#eef1f4" stroke="#c9d1d9"/>')
            add(f'<text x="{x0 + w/2:.1f}" y="{yy+11}" text-anchor="middle" font-family="{MONO}" '
                f'font-size="10" fill="#333">{escape(f)}</text>')
            x0 += w + 4

def above(x, row, title, sub=None, flags=(), anchor="middle"):
    h = 13 + (14 if sub else 0) + (21 if flags else 0)
    label(x, row - 16 - h + 11, title, sub, flags, anchor)

def below(x, row, title, sub=None, flags=(), anchor="middle"):
    label(x, row + 30, title, sub, flags, anchor)

# ------------------------------------------------------------------ layout
add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
add(f'<text x="24" y="30" font-family="{FONT}" font-size="20" font-weight="700" fill="{INK}">'
    f'hybrid-metatranscriptomics</text>')
add(f'<text x="300" y="30" font-family="{FONT}" font-size="13" fill="{SUB}">host-only reference '
    f'transcriptome from long + short reads, contaminant reads removed before assembly</text>')

T1, B1 = 44, 500      # tier 1 band
T2, B2 = 512, 916     # tier 2 band
band(20, 170, T1, B1, 1, "Inputs", "#f4f6f8")
band(176, 644, T1, B1, 2, "Read preprocessing", "#f9fafb")
band(650, 790, T1, B1, 3, "Decontamination", "#f4f6f8")
band(796, 1264, T1, B1, 4, "Four assembly routes", "#f9fafb")
band(1270, 1780, T1, B1, 5, "Reference ladder", "#f4f6f8")
band(1222, 1780, T2, B2, 6, "Quantification & annotation", "#f9fafb", right=1686)
band(796, 1216, T2, B2, 7, "Differential expression", "#f4f6f8")
band(20, 790, T2, B2, 8, "Reporting", "#f9fafb")

# rows
yA, yL, yB = 150, 180, 215          # route A, long reads, route B
yS, yS2 = 300, 314                  # short reads, assembly-only reads (parallel)
yD, yD2 = 360, 374                  # route D pair
yG = 440                            # genome / reference databases
yT = 262                            # transcriptome after the ladder
yU, yM, yN, yAn = 570, 640, 700, 800   # DTU, main quant/DE, gene maps, annotation

# ---- lines (drawn first, stations on top)
# references
line([(60, yG), (860, yG)], C["ref"])
# long reads -> A and B
line([(60, yL), (750, yL), (780, yA), (1320, yA)], C["long"])
line([(750, yL), (785, yB), (1320, yB)], C["long"])
# short reads -> C and D
line([(60, yS), (1320, yS)], C["short"])
line([(800, yS), (860, yD), (1320, yD)], C["short"])
# assembly-only reads, parallel to the short reads, into C (to StringTie) and D (Trinity)
line([(60, yD), (106, yS2), (1110, yS2)], C["asm"], dash="10 6")
line([(800, yS2), (860, yD2), (1000, yD2)], C["asm"], dash="10 6")
# transcriptome after the ladder
line([(1320, yT), (1700, yT), (1700, yM)], C["tx"])
# annotation branch
line([(1700, yM), (1700, yAn), (1000, yAn)], C["ann"])
# quantification & DE main line
line([(1700, yM), (860, yM)], C["de"])
# gene maps
line([(1500, yM), (1440, yN), (1160, yN)], C["gmap"])
# isoform switching (optional)
line([(1500, yM), (1430, yU), (1330, yU)], C["dtu"], dash="10 6")
# reporting
line([(860, yM), (420, yM)], C["rep"])

# ---- 1 inputs
terminus(60, yL, C["long"]);  above(40, yL, "Long reads", "Iso-Seq HQ isoforms", anchor="start")
terminus(60, yS, C["short"]); above(40, yS, "Short reads", "paired-end", anchor="start")
terminus(60, yD, C["asm"]);   below(40, yD, "Assembly-only reads", "optional, not tested", anchor="start")
terminus(60, yG, C["ref"]);   below(40, yG, "Genome + contaminants", "FASTA or accession", anchor="start")

# ---- 2 preprocessing
station(230, yL);     above(230, yL, "To FASTA", "seqkit fq2fa if FASTQ")
hub(230, yS, yS2);    above(230, yS, "Merge lanes", "cat")
hub(350, yS, yS2);    above(350, yS, "fastp", "trim, filter", ["--skip_fastp"])
hub(470, yL, yS2);    above(470, yL, "rRNA depletion", "BBDuk · SILVA/Rfam", ["--skip_rrna"])
hub(590, yS, yS2);    above(590, yS, "FastQC", "read QC")
station(350, yG);     below(350, yG, "Fetch references", "datasets · UniProt · SortMeRNA")
station(590, yG);     below(590, yG, "Combined index", "minimap2 host + contam.")

# ---- 3 decontamination
hub(720, yL, yG);     above(720, yL, "Read decontamination", "minimap2 · best hit wins", ["--skip_decontam"])
add(f'<text x="735" y="{yG+30}" font-family="{FONT}" font-size="13" font-weight="700" fill="{INK}">read-level</text>')
add(f'<text x="735" y="{yG+44}" font-family="{FONT}" font-size="11" fill="{SUB}">before assembly</text>')

# ---- 4 assembly routes
add(f'<text x="806" y="{yA-6}" font-family="{FONT}" font-size="12" font-weight="700" fill="{C["long"]}">A</text>')
add(f'<text x="806" y="{yB+24}" font-family="{FONT}" font-size="12" font-weight="700" fill="{C["long"]}">B</text>')
add(f'<text x="806" y="{yS-8}" font-family="{FONT}" font-size="12" font-weight="700" fill="{C["short"]}">C</text>')
add(f'<text x="840" y="{yD2+22}" font-family="{FONT}" font-size="12" font-weight="700" fill="{C["short"]}">D</text>')
station(860, yA);     above(860, yA, "pbmm2", "align to genome")
station(965, yA);     above(965, yA, "isoseq collapse", "genome loci")
station(1110, yA);    above(1110, yA, "Genome sequence", "models from genome", ["--longread_genome_sequence"])
station(1220, yA, optional=True); below(1220, yA, "SQANTI3", "QC report", ["--skip_sqanti"])
station(880, yB);     below(880, yB, "vsearch", "reference-free", ["--skip_reffree"])
station(860, yG);     below(860, yG, "STAR index", "")
hub(1000, yS, yS2);   above(1000, yS, "STAR", "align")
hub(1110, yS, yS2);   above(1110, yS, "StringTie", "per sample", ["--skip_stringtie"])
station(1220, yS);    above(1220, yS, "Merge", "gffcompare: new loci")
hub(1000, yD, yD2);   below(1000, yD2, "Trinity", "de novo", ["--skip_denovo"])
station(1150, yD);    below(1150, yD, "Novel loci", "vsearch + minimap2")

# ---- 5 ladder
hub(1320, yA, yD);    above(1320, yA, "Reference ladder", "A → +C → +D → +B → dedup", ["--dge_reference"])
station(1450, yT);    above(1450, yT, "BUSCO", "every rung", ["--skip_busco"])
station(1570, yT);    above(1570, yT, "Contribution", "per source, per rung")
terminus(1700, yT, C["tx"])
below(1688, yT, "Reference transcriptome", "selected rung", anchor="end")

# ---- 6 quantification & annotation
station(1620, yM);    below(1620, yM, "Expression filter", "TPM ≥ 1 in ≥ 3", ["--expr_filter_tpm"])
station(1500, yM);    above(1512, yM, "salmon", "decoy-aware quant", anchor="start")
station(1320, yN);    below(1320, yN, "Gene maps", "locus · tx2gene · Corset",
                            ["--dge_primary_map", "--skip_corset"])
station(1330, yU, optional=True); above(1330, yU, "IsoformSwitchAnalyzeR", "DTU, off by default",
                                        ["--skip_dtu false"])
station(1600, yAn);   below(1600, yAn, "TransDecoder", "ORFs → peptides", ["--skip_annotation"])
station(1480, yAn);   below(1480, yAn, "DIAMOND", "Swiss-Prot + Tox-Prot")
station(1360, yAn);   below(1360, yAn, "eggNOG-mapper", "GO / KEGG", ["--skip_eggnog"])
station(1250, yAn);   below(1250, yAn, "InterProScan", "Pfam, SMART, ...", ["--skip_interproscan"])

# ---- 7 differential expression
station(1110, yAn);   below(1110, yAn, "Annotation table", "per transcript & gene")
hub(1160, yM, yN);    above(1160, yM, "DESeq2", "tximport · apeglm", ["--skip_dge"])
hub(1000, yM, yAn);   above(1000, yM, "Enrichment", "clusterProfiler ORA/GSEA", ["--skip_enrichment"])
station(860, yM);     above(860, yM, "Figures", "volcano, PCA, heatmaps", ["--skip_figures"])

# ---- 8 reporting
station(700, yM);     above(700, yM, "MultiQC", "all QC", ["--skip_multiqc"])
station(560, yM);     above(560, yM, "Dashboard", "interactive HTML", ["--skip_dashboard"])
terminus(420, yM, C["rep"])
above(420, yM, "Results", "--outdir")

# ---- legend
lx, ly = 40, 700
add(f'<text x="{lx}" y="{ly}" font-family="{FONT}" font-size="13" font-weight="700" fill="{INK}">Legend</text>')
items = [("long", "Long reads", None), ("short", "Short reads", None), ("asm", "Assembly-only reads", "10 6"),
         ("ref", "Genome & databases", None), ("tx", "Transcriptome", None), ("ann", "Annotation", None),
         ("de", "Quantification & DE", None), ("gmap", "Gene maps", None), ("dtu", "Isoform switching", "10 6"),
         ("rep", "Reporting", None)]
for i, (k, name, dash) in enumerate(items):
    cx = lx + (i % 3) * 190
    cy = ly + 22 + (i // 3) * 24
    line([(cx, cy), (cx + 36, cy)], C[k], dash=dash, width=7)
    add(f'<text x="{cx+46}" y="{cy+4}" font-family="{FONT}" font-size="12" fill="{INK}">{escape(name)}</text>')
sx, sy = lx + 590, ly + 22
station(sx, sy); add(f'<text x="{sx+16}" y="{sy+4}" font-family="{FONT}" font-size="12" fill="{INK}">step</text>')
station(sx, sy + 24, optional=True)
add(f'<text x="{sx+16}" y="{sy+28}" font-family="{FONT}" font-size="12" fill="{INK}">optional / QC-only step</text>')
hub(sx, sy + 48, sy + 56)
add(f'<text x="{sx+16}" y="{sy+56}" font-family="{FONT}" font-size="12" fill="{INK}">shared step (several inputs)</text>')
label(sx - 8, sy + 92, "", None, ["--option"], anchor="start")
add(f'<text x="{sx+56}" y="{sy+105}" font-family="{FONT}" font-size="12" fill="{INK}">parameter that switches or sets the step</text>')

add('</svg>')
open(sys.argv[1] if len(sys.argv) > 1 else "assets/metro_map.svg", "w").write("\n".join(out) + "\n")
