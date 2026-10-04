#!/usr/bin/env python3
"""Build the interactive run dashboard: one self-contained HTML file, no external assets.

Reads whatever the run published (a directory of staged outputs) plus the run manifest, and
emits results AND methods side by side: every stage card carries the arguments it actually ran
with, the reason it exists, and the numbers it produced, so a figure can be traced back to the
step that made it without leaving the page.

Missing inputs are not errors — a skipped stage is reported as skipped, and a section with no
data is dropped. That keeps the dashboard honest for partial runs.

Usage: build_dashboard.py --dir <staged inputs> [--manifest run_manifest.json] --out index.html
"""
import argparse
import base64
import csv
import glob
import json
import math
import os
import re
import sys
from collections import OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="directory of staged pipeline outputs")
ap.add_argument("--manifest", default=None, help="run_manifest.json written by the workflow")
ap.add_argument("--out", required=True)
ap.add_argument("--title", default="Transcriptome run")
ap.add_argument("--max-points", type=int, default=2500,
                help="background (non-significant) volcano points kept per contrast")
ap.add_argument("--sig-ceiling", type=float, default=0.25,
                help="every gene with padj below this is embedded in full, so the volcano's "
                     "interactive thresholds give exact counts up to it")
ap.add_argument("--max-genes", type=int, default=300, help="gene-table rows kept per contrast")
args = ap.parse_args()

D = args.dir


# ------------------------------------------------------------------ small readers
def path(*parts):
    return os.path.join(D, *parts)


def exists(p):
    return p and os.path.exists(p)


def read_tsv(name, required_cols=None):
    """A TSV with a header row -> list of dicts. Comment lines (#) are dropped."""
    p = path(name)
    if not exists(p):
        return []
    with open(p, newline="") as fh:
        rows = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
    if required_cols and rows and not all(c in rows[0] for c in required_cols):
        return []
    return rows


def read_csv_rows(p):
    if not exists(p):
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def num(v, default=0.0):
    try:
        if v is None or v == "" or v == "NA":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def i(v, default=0):
    return int(num(v, default))


def first_glob(pattern):
    hits = sorted(glob.glob(path(pattern)))
    return hits[0] if hits else None


def one_dir(pattern):
    hits = [h for h in sorted(glob.glob(path(pattern))) if os.path.isdir(h)]
    return hits[0] if hits else None


# ------------------------------------------------------------------ manifest / params
manifest = {}
if args.manifest and os.path.exists(args.manifest):
    try:
        manifest = json.load(open(args.manifest))
    except (ValueError, OSError) as e:
        sys.stderr.write(f"[dashboard] could not read manifest: {e}\n")

params = manifest.get("params", {}) or {}
wf = manifest.get("workflow", {}) or {}


def P(key, default=None):
    v = params.get(key, default)
    return default if v in (None, "null", "") else v


def on(key):
    """A step ran when its skip flag is false."""
    v = params.get(key)
    return not (v is True or str(v).lower() == "true")


# ------------------------------------------------------------------ per-sample QC
rrna = {r["sample"]: r for r in read_tsv("rrna_per_sample.tsv")}
decon = {r["sample"]: r for r in read_tsv("decontam_per_sample.tsv")}
mapping = {r["sample"]: r for r in read_tsv("mapping_rates.tsv")}

sample_ids = sorted(set(rrna) | set(decon) | set(mapping))
samples = []
for s in sample_ids:
    r, d, m = rrna.get(s, {}), decon.get(s, {}), mapping.get(s, {})
    reads_in = i(r.get("reads_in")) or i(d.get("pairs_in"))
    kept = i(d.get("pairs_kept")) or i(r.get("reads_kept")) or reads_in
    samples.append({
        "sample": s,
        "condition": r.get("condition") or d.get("condition") or "",
        "reads_in": reads_in,
        "rrna_removed": i(r.get("rrna_removed")),
        "contam_removed": i(d.get("contaminant_removed")),
        "kept": kept,
        "rrna_pct": num(r.get("rrna_pct")),
        "contam_pct": num(d.get("pct_removed")),
        "mapped_pct": num(m.get("percent_mapped")),
        "fragments": i(m.get("fragments_processed")),
        "library_type": m.get("library_type", ""),
    })

# ------------------------------------------------------- assembly-only short reads
# Libraries given through --assembly_samplesheet: preprocessed like the experiment's, assembled
# (StringTie + Trinity), never quantified or tested. Kept apart from `samples` so every per-sample
# QC number, the accounting check and the DE tabs still describe only the experiment.
_ar = {r["sample"]: r for r in read_tsv("rrna_assembly_only.tsv")}
_ad = {r["sample"]: r for r in read_tsv("decontam_assembly_only.tsv")}
assembly_reads = []
for s in sorted(set(_ar) | set(_ad)):
    r, d = _ar.get(s, {}), _ad.get(s, {})
    assembly_reads.append({"sample": s,
                           "reads_in": i(r.get("reads_in")) or i(d.get("pairs_in")),
                           "rrna_removed": i(r.get("rrna_removed")),
                           "contam_removed": i(d.get("contaminant_removed")),
                           "kept": i(d.get("pairs_kept")) or i(r.get("reads_kept"))})


def lanes_per_sample(sheet):
    """Rows per sample in a sample sheet: the lane count when it is the same for every sample."""
    try:
        with open(str(sheet), newline="") as fh:
            n = defaultdict(int)
            for row in csv.DictReader(fh):
                n[row.get("sample", "")] += 1
        vals = set(n.values())
        return vals.pop() if len(vals) == 1 else None
    except (OSError, TypeError):
        return None


lanes = lanes_per_sample(P("samplesheet"))

# ------------------------------------------------------------------ long reads
lr_rrna = read_tsv("longread_rrna.tsv")
lr_samples = read_tsv("longread_per_sample.tsv")
lr_note = ""
_p = path("longread_per_sample.tsv")
if exists(_p):
    for line in open(_p):
        if line.startswith("#"):
            lr_note = line.lstrip("#").strip()
for r in lr_samples:
    for k in ("isoforms_in", "rrna_removed", "host", "no_alignment", "contaminant",
              "ambiguous", "contaminant_removed", "kept"):
        r[k] = i(r.get(k))
    r["pct_removed"] = num(r.get("pct_removed"))
lr_decon = []
p = path("longread_decontam.tsv")
if exists(p):
    with open(p) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 2 and f[0] not in ("class",):
                lr_decon.append({"class": f[0], "n": i(f[1])})

# ------------------------------------------------------------------ reference ladder + BUSCO
BUSCO_RE = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%")


def parse_busco(text):
    m = BUSCO_RE.search(text or "")
    if not m:
        return None
    c, s, d, f, mi = (float(x) for x in m.groups())
    return {"complete": c, "single": s, "dup": d, "frag": f, "missing": mi}


ladder = []
for r in read_tsv("route_contribution.tsv"):
    row = {
        "reference": r.get("reference", ""),
        "basis": r.get("basis", ""),
        "transcripts": i(r.get("n_transcripts")),
        "bp": i(r.get("total_bp")),
        "n50": i(r.get("n50")),
        "genes": i(r.get("n_genes")) if r.get("n_genes") not in ("-", "") else None,
        "busco": parse_busco(r.get("BUSCO", "")),
    }
    if row["busco"] is None:
        bp = path(f"{row['reference']}.busco.txt")
        if exists(bp):
            row["busco"] = parse_busco(open(bp).read())
    ladder.append(row)

# The expression-filtered set is a rung too, and the most consequential one: it is the reference
# every DE number is computed against. route_contribution.tsv is written before the filter runs, so
# it cannot carry it -- but SEQ_STATS_EXPR and BUSCO_EXPR both publish for it. Appending it here is
# what makes the assembly and the DE reference comparable on one line each, which is the whole
# point of reporting completeness on the assembly and testing on the subset.
for _st in sorted(glob.glob(path("*_expressed.stats.tsv"))):
    _rows = read_tsv(os.path.basename(_st))
    if not _rows:
        continue
    _name = os.path.basename(_st)[:-len(".stats.tsv")]
    if any(r["reference"] == _name for r in ladder):
        continue
    _bp = path(f"{_name}.busco.txt")
    ladder.append({
        "reference": _name,
        "basis": "expression-filtered (the reference DE and enrichment run against)",
        "transcripts": i(_rows[0].get("num_seqs")),
        "bp": i(_rows[0].get("sum_len")),
        "n50": None,
        "genes": None,
        "busco": parse_busco(open(_bp).read()) if exists(_bp) else None,
    })

# ------------------------------------------------------------------ annotation
annot_summary = []
sp = first_glob("annotation.*summary*.txt")
if sp:
    for line in open(sp):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        m = re.match(r"^(.*?)\s+([\d,]+)$", line)
        if m:
            annot_summary.append({"label": m.group(1).strip(), "value": int(m.group(2).replace(",", ""))})

annot_path = None
for cand in sorted(glob.glob(path("annotation.*.tsv"))):
    if "locus" not in os.path.basename(cand):
        annot_path = cand
        break
    annot_path = annot_path or cand

annot_index = {}      # gene -> {desc, flag}
annot_evidence = []
interest_label = P("interest_label", "family of interest")
if annot_path:
    with open(annot_path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames or []
        counts = OrderedDict((k, 0) for k in
                             ["ORF (coding)", "protein hit", "InterPro", "GO", "KEGG KO", interest_label])
        total = 0
        for r in rd:
            total += 1
            if r.get("coding_status") == "coding":
                counts["ORF (coding)"] += 1
            for key, col in (("protein hit", "protein_hit"), ("InterPro", "InterPro"),
                             ("GO", "GOs"), ("KEGG KO", "KEGG_ko")):
                if r.get(col):
                    counts[key] += 1
            if r.get(interest_label):
                counts[interest_label] += 1
            g = r.get("gene")
            if g and g not in annot_index:
                desc = r.get("protein_desc") or r.get("eggnog_desc") or ""
                annot_index[g] = {"desc": desc.split(" OS=")[0][:90],
                                  "flag": bool(r.get(interest_label))}
    annot_evidence = [{"label": k, "n": v, "total": total} for k, v in counts.items() if v or k in counts]
    annot_evidence = [e for e in annot_evidence if e["n"] > 0]

# ------------------------------------------------------------------ differential expression
def load_dge(dirpath):
    if not dirpath:
        return None
    out = {"dir": os.path.basename(dirpath), "contrasts": [], "pca": []}
    for r in read_csv_rows(os.path.join(dirpath, "pca_data.csv")):
        grp = r.get("condition") or r.get("group") or ""
        out["pca"].append({"x": num(r.get("PC1")), "y": num(r.get("PC2")),
                           "group": grp, "name": r.get("name") or r.get("", "")})
    summary = {r["contrast"]: r for r in read_csv_rows(os.path.join(dirpath, "DE_summary.csv"))
               if r.get("contrast")}
    for f in sorted(glob.glob(os.path.join(dirpath, "DE_*.csv"))):
        if os.path.basename(f) == "DE_summary.csv":
            continue
        name = os.path.basename(f)[3:-4]
        rows = read_csv_rows(f)
        if not rows:
            continue
        key = list(rows[0].keys())[0]          # the gene id column has an empty header
        below, background, top = [], [], []
        for r in rows:
            padj, lfc = r.get("padj"), num(r.get("log2FoldChange"))
            if padj in (None, "", "NA"):
                continue
            pa = num(padj, 1.0)
            # 5 decimals on -log10(padj): at 3 the rounding flips genes sitting at padj ~ 0.0499,
            # and the volcano's live counts then disagree with DESeq2's own summary by a gene or two
            pt = [round(lfc, 4), round(-math.log10(max(pa, 1e-300)), 5)]
            # Everything below the ceiling is embedded in FULL. The volcano's padj filter is
            # interactive, so a thinned point near the threshold would make its live count wrong;
            # only the clearly non-significant cloud is sampled, and it is context, never counted.
            (below if pa < args.sig_ceiling else background).append(pt)
        step = max(1, len(background) // max(1, args.max_points))
        points = below + background[::step]
        for r in rows[:args.max_genes]:
            if r.get("padj") in (None, "", "NA"):
                continue
            g = r[key]
            a = annot_index.get(g, {})
            top.append({"gene": g, "base": round(num(r.get("baseMean")), 1),
                        "lfc": round(num(r.get("log2FoldChange")), 2),
                        "padj": float(f"{num(r.get('padj'), 1.0):.3g}"),
                        "desc": a.get("desc", ""), "flag": bool(a.get("flag"))})
        s = summary.get(name, {})
        out["contrasts"].append({
            "id": name,
            "tested": i(s.get("tested")) or len([1 for r in rows if r.get("padj") not in ("", "NA", None)]),
            # two different questions, kept apart: sig_padj asks only about padj, sig also
            # applies the fold-change cut and is what the bars and the volcano colour by
            "sig_padj": i(s.get("sig_padj")),
            "sig": i(s.get("sig")) or (i(s.get("up")) + i(s.get("down"))),
            "up": i(s.get("up")), "down": i(s.get("down")),
            "points": points, "n_points_total": len(rows), "top": top,
            "ceiling": args.sig_ceiling, "n_exact": len(below),
            "n_background": len(background), "n_background_shown": len(background[::step]),
        })
    return out if out["contrasts"] or out["pca"] else None


MAP_KINDS = ("string", "corset", "locus")


def map_kind_of(label):
    """Which gene map a <prefix>_<label> directory belongs to, from the label's suffix.

    The suffixes are the labels the processes emit: LOCUS_TX2GENE -> <ref>_locus, CORSET ->
    <ref>_corset; the string-stripped map carries the bare reference label.
    """
    return "locus" if label.endswith("_locus") else "corset" if label.endswith("_corset") else "string"


def map_dir(prefix, kind):
    """The <prefix>_<label> directory of one gene map, e.g. map_dir("dge", "locus")."""
    for c in sorted(glob.glob(path(prefix + "_*"))):
        if os.path.isdir(c) and map_kind_of(os.path.basename(c)[len(prefix) + 1:]) == kind:
            return c
    return None


dge = load_dge(map_dir("dge", "string"))
dge_corset = load_dge(map_dir("dge", "corset"))
dge_locus = load_dge(map_dir("dge", "locus"))
dge_by_map = {"string": dge, "corset": dge_corset, "locus": dge_locus}

# Every headline DE number, the diagnostics and the enrichment default to this map. Falling back
# to the string map only when the named one did not run keeps a skip_corset / no-locus run working.
primary_map = str(P("dge_primary_map", "string") or "string")
if not dge_by_map.get(primary_map):
    primary_map = next((k for k in MAP_KINDS if dge_by_map.get(k)), primary_map)
dge_primary = dge_by_map.get(primary_map)


def map_outputs(prefix):
    """The published folders of every map that ran, as conf/modules.config names them."""
    dirs = [f"{prefix}/" if k == "string" else f"{prefix}_{k}/"
            for k in MAP_KINDS if map_dir(prefix, k)]
    return " · ".join(dirs) or f"{prefix}/"

# ------------------------------------------------------------------ enrichment
def load_enrichment(dirpath):
    if not dirpath:
        return {}
    out = {}
    for f in sorted(glob.glob(os.path.join(dirpath, "*.csv"))):
        base = os.path.basename(f)[:-4]
        if base.endswith("_summary") or base == "enrichment_summary":
            continue
        m = re.match(r"^(.*)_(GO|KEGG)_(ORA|GSEA)$", base)
        if not m:
            continue
        contrast, db, method = m.groups()
        rows = []
        for r in read_csv_rows(f):
            padj = num(r.get("p.adjust"), 1.0)      # read p.adjust BY NAME - see README
            rows.append({"id": r.get("ID", ""), "desc": (r.get("Description") or r.get("ID") or "")[:70],
                         "padj": padj, "nes": num(r.get("NES"), 0.0),
                         "size": i(r.get("setSize") or r.get("Count"))})
        rows.sort(key=lambda x: x["padj"])
        out.setdefault(contrast, {})[f"{db} {method}"] = rows[:20]
    return out


enrichment_by_map = {k: e for k in MAP_KINDS for e in [load_enrichment(map_dir("enrichment", k))] if e}


def count_enrichment(dirpath, alpha=0.05):
    """Significant terms per analysis, counted over the FULL csv.

    load_enrichment keeps only the top 20 rows per analysis for plotting, so counting from its
    payload would silently cap any analysis that returned more. The methods diagram quotes these
    numbers next to each gene map, so they have to be the real ones.
    """
    out = {}
    if not dirpath:
        return out
    for f in sorted(glob.glob(os.path.join(dirpath, "*.csv"))):
        m = re.match(r"^(.*)_(GO|KEGG)_(ORA|GSEA)$", os.path.basename(f)[:-4])
        if not m:
            continue
        _, db, method = m.groups()
        out[f"{db} {method}"] = out.get(f"{db} {method}", 0) + \
            sum(1 for r in read_csv_rows(f) if num(r.get("p.adjust"), 1.0) < alpha)
    return out


# ------------------------------------------------------------------ locus rebuild
locus = None
lp = path("locus_rebuild_report.txt")
if exists(lp):
    text = open(lp).read()
    def grab(pattern, cast=int):
        m = re.search(pattern, text)
        if not m:
            return None
        try:
            return cast(m.group(1).replace(",", ""))
        except ValueError:
            return None
    locus = {
        "before": grab(r"gene count BEFORE.*?:\s*([\d,]+)"),
        "after": grab(r"gene count AFTER.*?:\s*([\d,]+)"),
        "placed": grab(r"placed on the genome.*?:\s*([\d,]+)"),
        "unplaced": grab(r"unplaced.*?:\s*([\d,]+)"),
        "ratio_before": grab(r"BEFORE:.*?([\d.]+) genes/accession", float),
        "ratio_after": grab(r"AFTER :.*?([\d.]+) genes/accession", float),
        "overmerge_before": grab(r"BEFORE: [\d,]+ / [\d,]+ annotated genes \(([\d.]+)%\)", float),
        "overmerge_after": grab(r"AFTER : [\d,]+ / [\d,]+ annotated loci  \(([\d.]+)%\)", float),
        "warning": "WARNING: over-merging rose sharply" in text,
        "report": text,
    }


# ------------------------------------------------- where the three gene maps agree
# The maps cannot be compared on gene ids — each invents its own (PB.1 / LOCUS_000123 /
# Cluster-7.2). The common denominator is the TRANSCRIPT: for each map, take the genes it calls
# significant and expand them to the transcripts they contain. Two maps agree about a transcript
# when both implicate it, whatever they call the gene it belongs to.
def gene_to_transcripts(t2g_path):
    g2t = {}
    try:
        with open(t2g_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    g2t.setdefault(parts[1], []).append(parts[0])
    except OSError:
        return {}
    return g2t


def significant_transcripts(dge_dir, g2t, alpha, lfc_cut):
    """contrast -> set of transcripts belonging to a significant gene."""
    out = {}
    for f in sorted(glob.glob(os.path.join(dge_dir, "DE_*.csv"))):
        base = os.path.basename(f)
        if base.startswith("DE_summary"):
            continue
        contrast = base[3:-4]
        tx = set()
        try:
            with open(f) as fh:
                rd = csv.reader(fh)
                header = next(rd, None)
                if not header:
                    continue
                idx = {name: k for k, name in enumerate(header)}
                gi = 0
                pi, li = idx.get("padj"), idx.get("log2FoldChange")
                if pi is None or li is None:
                    continue
                for row in rd:
                    if len(row) <= max(pi, li):
                        continue
                    try:
                        padj, lfc = float(row[pi]), float(row[li])
                    except ValueError:
                        continue
                    if padj < alpha and abs(lfc) >= lfc_cut:
                        tx.update(g2t.get(row[gi], ()))
        except OSError:
            continue
        out[contrast] = tx
    return out


venn = {}
_alpha = num(P("de_alpha", 0.05), 0.05)
_lfc = num(P("de_lfc", 1), 1.0)
_per_map = {}
for d in sorted(glob.glob(path("dge_*"))):
    if not os.path.isdir(d):
        continue
    label = os.path.basename(d)[4:]                      # dge_<label>
    t2g = path(f"tx2gene.{label}.tsv")
    if not exists(t2g):
        continue
    _per_map[map_kind_of(label)] = significant_transcripts(d, gene_to_transcripts(t2g), _alpha, _lfc)

if len(_per_map) >= 2:
    kinds = [k for k in ("string", "corset", "locus") if k in _per_map]
    contrasts = sorted(set().union(*[set(v) for v in _per_map.values()]))
    NAMES = {"string": "String-stripped", "corset": "Corset", "locus": "Locus rebuild"}
    for c in contrasts:
        sets = {k: _per_map[k].get(c, set()) for k in kinds}
        regions, total = {}, set().union(*sets.values()) if sets else set()
        if len(kinds) == 3:
            a, b, cc = (sets[kinds[0]], sets[kinds[1]], sets[kinds[2]])
            regions = {
                "A": len(a - b - cc), "B": len(b - a - cc), "C": len(cc - a - b),
                "AB": len((a & b) - cc), "AC": len((a & cc) - b), "BC": len((b & cc) - a),
                "ABC": len(a & b & cc),
            }
        elif len(kinds) == 2:
            a, b = sets[kinds[0]], sets[kinds[1]]
            regions = {"A": len(a - b), "B": len(b - a), "AB": len(a & b)}
        venn[c] = {"labels": [NAMES[k] for k in kinds],
                   "sizes": [len(sets[k]) for k in kinds],
                   "regions": regions, "union": len(total)}

# ------------------------------------------------------- diagnostics embedded in the page
# The figure set already exists as PNGs next to the dashboard, but a reader who opens the
# dashboard should not have to go and find them: the panels that say whether the DE results can
# be believed at all (dispersion fit, p-value shape, MA, sample clustering) are embedded here.
FIGURES_SHOWN = [
    ("F11_dispersion", "Dispersion fit",
     "DESeq2 shrinks each gene's dispersion toward a fitted trend. A cloud that does not follow "
     "the red line means the model is wrong and every p-value below is suspect."),
    ("F12_pvalue_histograms", "p-value distribution",
     "A flat bulk with a spike near zero is what a healthy test looks like. A hill in the middle "
     "means the null is mis-specified; a spike at one means over-conservative filtering."),
    ("F10_sample_distances", "Sample-to-sample distance",
     "Euclidean distance on variance-stabilised counts. If the treatment were the dominant axis "
     "of variation, replicates would block together here."),
]

def load_figures(d):
    """One map's diagnostic PNGs -> data URIs, so the page stays a single file that works offline."""
    out = {}
    if not d:
        return out
    for f in sorted(glob.glob(os.path.join(d, "*.png"))):
        name = os.path.splitext(os.path.basename(f))[0]
        keep = any(name == n for n, _, _ in FIGURES_SHOWN) or name.startswith("F13_MA_")
        if not keep:
            continue
        try:
            with open(f, "rb") as fh:
                out[name] = "data:image/png;base64," + base64.b64encode(fh.read()).decode()
        except OSError:
            pass
    return out


# Per map: each figure set is drawn from that map's DESeq2 fit, so the Expression tab shows the
# set matching the map its switcher is on.
figures_by_map = {k: f for k in MAP_KINDS for f in [load_figures(map_dir("figures", k))] if f}
_fig_names = set().union(*figures_by_map.values()) if figures_by_map else set()
figure_meta = [{"name": n, "title": t, "why": w} for n, t, w in FIGURES_SHOWN if n in _fig_names]
ma_figures = sorted(n for n in _fig_names if n.startswith("F13_MA_"))

# ---- the three gene maps, side by side rather than three separate claims
gene_maps = []
corset_summary = first_glob("*_corset.summary.txt")
corset = None
if corset_summary:
    ctext = open(corset_summary).read()

    def cgrab(pattern):
        m = re.search(pattern, ctext)
        return i(m.group(1).replace(",", "")) if m else None
    corset = {"transcripts": cgrab(r"transcripts clustered\s*:\s*([\d,]+)"),
              "clusters": cgrab(r"clusters \(genes\)\s*:\s*([\d,]+)"),
              "stripped": cgrab(r"string-stripped genes\s*:\s*([\d,]+)")}

def de_total(d):
    """Significant genes summed over contrasts, so the maps compare on results not just counts."""
    if not d or not d.get("contrasts"):
        return None
    return sum(i(c.get("sig")) for c in d["contrasts"])


if corset and corset.get("stripped"):
    gene_maps.append({"primary": primary_map == "string", "name": "String-stripped ids", "basis": "per-source naming rules",
                      "genes": corset["stripped"], "sig": de_total(dge),
                      "note": "the default map; cannot reconcile identity ACROSS assembly sources"})
if corset and corset.get("clusters"):
    gene_maps.append({"primary": primary_map == "corset", "name": "Corset", "basis": "shared reads + expression",
                      "genes": corset["clusters"], "sig": de_total(dge_corset),
                      "note": "clusters contigs by the reads they share, not by sequence identity"})
if locus and locus.get("after"):
    gene_maps.append({"primary": primary_map == "locus", "name": "Locus rebuild", "basis": "genomic coordinates",
                      "genes": locus["after"], "sig": de_total(dge_locus),
                      "note": "groups transcripts sharing a splice junction on the genome"})

# ---- the three DE branches, for the methods diagram
# The Genes tab already switches between the maps, but the diagram has to show all three side by
# side: they are parallel branches of the same run, not alternatives a reader picks between. Each
# branch carries its own DESeq2 result and its own enrichment, and which one is primary is a
# parameter (dge_primary_map), not a property of the method.
_enr_alpha = float(P("de_alpha", 0.05) or 0.05)
map_flow = []
for _key, _label, _basis, _d, _enrdir in [
    ("string", "String-stripped ids", "per-source naming rules", dge, map_dir("enrichment", "string")),
    ("corset", "Corset", "shared reads + expression", dge_corset, map_dir("enrichment", "corset")),
    ("locus", "Locus rebuild", "genomic coordinates", dge_locus, map_dir("enrichment", "locus")),
]:
    if not _d or not _d.get("contrasts"):
        continue
    _cs = _d["contrasts"]
    map_flow.append({
        "key": _key, "label": _label, "basis": _basis,
        "primary": primary_map == _key,
        "tested": max([i(c.get("tested")) for c in _cs] or [0]),
        "sig": sum(i(c.get("sig")) for c in _cs),
        "up": sum(i(c.get("up")) for c in _cs),
        "down": sum(i(c.get("down")) for c in _cs),
        "contrasts": len(_cs),
        "per_contrast": [{"id": c["id"], "sig": i(c.get("sig"))} for c in _cs],
        "terms": count_enrichment(_enrdir, _enr_alpha),
        "enr_alpha": _enr_alpha,
    })

# ---- what the DE filters removed, per contrast
def load_filtering(d):
    fs = os.path.join(d, "filtering_summary.csv") if d else None
    if not fs or not exists(fs):
        return []
    return [{"contrast": r.get("contrast", ""),
             "in_model": i(r.get("genes_in_model")),
             "tested": i(r.get("tested")),
             "threshold": r.get("filter_threshold_basemean", ""),
             "filtered": i(r.get("dropped_independent_filter")),
             "cooks": i(r.get("dropped_cooks_outlier")),
             "shrinkage": r.get("shrinkage", "")} for r in read_csv_rows(fs)]


de_filtering_by_map = {k: f for k in MAP_KINDS for f in [load_filtering(map_dir("dge", k))] if f}

# ------------------------------------------------------------------ headline numbers
def fmt(n):
    if n is None:
        return "—"
    if isinstance(n, float) and not n.is_integer():
        return f"{n:,.2f}"
    return f"{int(n):,}"


ref_set = P("dge_reference", "all_sources_dedup")
ref_row = next((r for r in ladder if r["reference"] == ref_set), None) or (ladder[-1] if ladder else None)

total_in = sum(s["reads_in"] for s in samples)
total_kept = sum(s["kept"] for s in samples)
total_contam = sum(s["contam_removed"] for s in samples)
total_rrna = sum(s["rrna_removed"] for s in samples)
mapped = [s["mapped_pct"] for s in samples if s["mapped_pct"]]
# the headline number uses the SAME rule the bars and volcanoes use, so the page cannot quote two
# different totals for "differentially expressed"
de_total = sum(c["sig"] for c in dge_primary["contrasts"]) if dge_primary else 0
de_total_padj_only = sum(c["sig_padj"] for c in dge_primary["contrasts"]) if dge_primary else 0

# ------------------------------------------------------------- the datasets this run touched
# One block per dataset, built only from what is present. Nothing is listed unconditionally, so a
# short-read-only run, a run without annotation and a run of another species all describe
# themselves correctly without the page knowing anything about this project.
datasets = []


def pct_of(part, whole, digits=2):
    return f"{100.0 * part / whole:.{digits}f}%" if whole else "—"


def block(bid, title, kind, subtitle, rows, tab, funnel=None, note=""):
    rows = [r for r in rows if r and r.get("value") not in (None, "", "—")]
    if not rows and not funnel:
        return
    datasets.append({"id": bid, "title": title, "kind": kind, "subtitle": subtitle,
                     "rows": rows, "tab": tab, "funnel": funnel or [], "note": note})


def R(label, value, note=""):
    return {"label": label, "value": value, "note": note}


def step(label, value, note=""):
    return {"label": label, "value": value, "note": note}


# ---- short reads -------------------------------------------------------------------
if samples:
    by_cond = defaultdict(int)
    for x in samples:
        by_cond[x["condition"] or "—"] += 1
    cond_txt = ", ".join(f"{k} {v}" for k, v in sorted(by_cond.items()))
    unit = (rrna[samples[0]["sample"]].get("unit") if rrna.get(samples[0]["sample"]) else None) or "pairs"
    sr_funnel = []
    if total_in:
        sr_funnel.append(step(f"{unit} sequenced", total_in))
    if total_rrna:
        sr_funnel.append(step("after rRNA depletion", total_in - total_rrna))
    if total_contam:
        sr_funnel.append(step("after contaminant removal", total_kept))
    elif total_kept and not total_rrna:
        sr_funnel.append(step("read pairs kept", total_kept))
    if frag_total := sum(x["fragments"] for x in samples):
        sr_funnel.append(step("quantified", frag_total))
    block("short_reads", "Short reads", "reads",
          f"{len(samples)} samples · {len(by_cond)} condition{'s' if len(by_cond) != 1 else ''} ({cond_txt})",
          [R("rRNA removed", f"{fmt(total_rrna)} {unit}" if total_rrna else "",
             pct_of(total_rrna, total_in) + " of sequenced" if total_rrna else ""),
           R("contaminant removed", f"{fmt(total_contam)} {unit}" if total_contam else "",
             pct_of(total_contam, total_contam + total_kept) + " of what reached it" if total_contam else ""),
           R("kept for assembly", fmt(total_kept) if total_kept else "", f"host {unit}"),
           R("mean mapping rate", f"{sum(mapped) / len(mapped):.2f}%" if mapped else "",
             f"{min(mapped):.2f}–{max(mapped):.2f}% across samples" if mapped else ""),
           R("library type", samples[0]["library_type"] or "", "auto-detected by salmon")],
          "qc", sr_funnel)

# ---- assembly-only short reads ------------------------------------------------------
if assembly_reads:
    a_in = sum(x["reads_in"] for x in assembly_reads)
    a_kept = sum(x["kept"] for x in assembly_reads)
    a_rrna = sum(x["rrna_removed"] for x in assembly_reads)
    a_contam = sum(x["contam_removed"] for x in assembly_reads)
    block("assembly_reads", "Short reads, assembly only", "reads",
          f"{len(assembly_reads)} librar{'ies' if len(assembly_reads) != 1 else 'y'} · "
          f"{', '.join(x['sample'] for x in assembly_reads)} · never quantified or tested",
          [R("rRNA removed", f"{fmt(a_rrna)} pairs" if a_rrna else "", pct_of(a_rrna, a_in) + " of sequenced" if a_rrna else ""),
           R("contaminant removed", f"{fmt(a_contam)} pairs" if a_contam else "",
             pct_of(a_contam, a_contam + a_kept) + " of what reached it" if a_contam else ""),
           R("kept for assembly", fmt(a_kept) if a_kept else "", "host pairs, into StringTie and Trinity only"),
           R("sample sheet", os.path.basename(str(P("assembly_samplesheet", ""))), "--assembly_samplesheet")],
          "methods",
          [step("pairs sequenced", a_in)] + ([step("after rRNA depletion", a_in - a_rrna)] if a_rrna else [])
          + ([step("after contaminant removal", a_kept)] if a_contam else []),
          "Same preprocessing as the experiment's libraries. They add assembly evidence only: salmon, "
          "the expression filter and DESeq2 see the experiment's libraries alone.")

# ---- long reads --------------------------------------------------------------------
lr_first = lr_rrna[0] if lr_rrna else None
lr_cls = {str(r["class"]).split()[0].lower(): r["n"] for r in lr_decon}
lr_kept, lr_removed = lr_cls.get("kept"), lr_cls.get("removed")
lr_pre = i(lr_first.get("kept")) if lr_first else (
    (lr_kept + lr_removed) if (lr_kept is not None and lr_removed is not None) else None)
lr_in = i(lr_first.get("isoforms_in")) if lr_first else lr_pre
if lr_in or lr_samples:
    rung = {r["reference"]: r for r in ladder}
    lr_funnel = []
    if lr_in:
        lr_funnel.append(step("sequences in", lr_in))
    if lr_first and i(lr_first.get("rrna_removed")):
        lr_funnel.append(step("after rRNA depletion", lr_pre))
    if lr_removed:
        lr_funnel.append(step("after contaminant removal", lr_kept))
    if "isoseq" in rung:
        lr_funnel.append(step("collapsed models", rung["isoseq"]["transcripts"]))
    if "isoseq_stringtie_lrfree" in rung and "isoseq_stringtie" in rung:
        lr_funnel.append(step("added reference-free",
                              rung["isoseq_stringtie_lrfree"]["transcripts"] - rung["isoseq_stringtie"]["transcripts"]))
    n_lr = len(lr_samples) if lr_samples else (len(lr_rrna) or 1)
    a_rung = rung.get("isoseq")
    block("long_reads", "Long reads", "reads",
          f"{n_lr} sample{'s' if n_lr != 1 else ''} · pooled for reference building",
          [R("rRNA removed", fmt(i(lr_first.get("rrna_removed"))) if lr_first and i(lr_first.get("rrna_removed")) else "",
             (lr_first.get("pct_removed", "") + "% of input") if lr_first else ""),
           R("contaminant removed", fmt(lr_removed) if lr_removed else "",
             pct_of(lr_removed, lr_pre) + " of what reached it" if lr_removed and lr_pre else ""),
           R("kept for assembly", fmt(lr_kept) if lr_kept else ""),
           R("BUSCO of this source alone",
             f"{a_rung['busco']['complete']:.1f}%" if a_rung and a_rung.get("busco") else "",
             "before the short-read routes are layered on"),
           R("per-sample detail", f"{len(lr_samples)} samples in Reads & QC" if lr_samples else "",
             "recovered from the membership map" if lr_samples else "")],
          "qc", lr_funnel,
          "Long reads are the only source of full-length isoform structure, so they set the "
          "isoform boundaries every short-read route is measured against.")

# ---- references --------------------------------------------------------------------
decon_scaffolds = {}
_dl = path("decon_index.log")
if exists(_dl):
    for line in open(_dl):
        m = re.match(r"scaffolds (\w+):\s+(\d+)", line.strip())
        if m:
            decon_scaffolds[m.group(1)] = int(m.group(2))

genome_src = P("genome") or P("genome_accession")
if genome_src:
    total_scaf = decon_scaffolds.get("total")
    contam_scaf = decon_scaffolds.get("contaminant")
    block("genome", "Reference genome", "reference",
          str(P("genome_name", "") or os.path.basename(str(genome_src))),
          [R("source", "downloaded by accession" if not P("genome") else "local FASTA",
             str(P("genome_accession", "")) if not P("genome") else os.path.basename(str(P("genome")))),
           R("scaffolds", fmt(total_scaf - contam_scaf) if total_scaf and contam_scaf is not None else "",
             "counted in the combined decontamination index"),
           R("STAR index", "prebuilt" if P("star_index") else ("built by the run" if on("skip_stringtie") else ""),
             "no annotation in the index: junctions are discovered de novo"),
           R("used as", "assembly target, salmon decoys, locus grouping", "")],
          "methods")

contam_src = P("contaminant_fasta") or P("contaminant_dir") or P("contaminant_accessions")
if contam_src and on("skip_decontam"):
    block("contaminants", "Contaminant genomes", "reference",
          f"removed at the read level, before assembly · prefix {P('contaminant_prefix', 'CONTAM_')}",
          [R("source", "accessions, downloaded" if P("contaminant_accessions") and not P("contaminant_fasta")
             else ("directory of FASTAs" if P("contaminant_dir") and not P("contaminant_fasta") else "local FASTA"),
             str(P("contaminant_accessions", "")) if not P("contaminant_fasta") else os.path.basename(str(P("contaminant_fasta")))),
           R("scaffolds", fmt(decon_scaffolds.get("contaminant")) if decon_scaffolds.get("contaminant") else ""),
           R("short-read pairs removed", fmt(total_contam) if total_contam else "",
             pct_of(total_contam, total_contam + total_kept) + " pooled" if total_contam else ""),
           R("long-read sequences removed", fmt(lr_removed) if lr_removed else ""),
           R("ambiguous placements", "kept" if P("keep_ambiguous_long") else "dropped",
             "a sequence that scores alike on both references")],
          "methods")

# ---- the reference transcriptome this run built ------------------------------------
if ref_row:
    genes = ref_row.get("genes")
    locus_after = locus.get("after") if locus else None
    block("transcriptome", "Reference transcriptome", "product",
          f"{ref_set} · {len(ladder)} variant{'s' if len(ladder) != 1 else ''} built and scored",
          [R("transcripts", fmt(ref_row["transcripts"])),
           R("genes", fmt(genes) if genes else "",
             f"upper bound; {fmt(locus_after)} by genomic locus" if locus_after else "from string-stripped ids"),
           R("N50", fmt(ref_row["n50"]) if ref_row.get("n50") else ""),
           R("BUSCO complete",
             f"{ref_row['busco']['complete']:.1f}%" if ref_row.get("busco") else "",
             f"{P('busco_lineage', '')} · duplicated {ref_row['busco']['dup']:.1f}%" if ref_row.get("busco") else ""),
           R("sources combined", str(len(ladder)), "each layered on and measured separately")],
          "reference")

# ---- annotation --------------------------------------------------------------------
if annot_evidence:
    tot = annot_evidence[0]["total"]
    ev = {e["label"]: e["n"] for e in annot_evidence}
    block("annotation", "Annotation", "product", f"evidence on {fmt(tot)} transcripts",
          [R(k, f"{fmt(v)}", pct_of(v, tot, 1)) for k, v in ev.items()],
          "annotation")

# ---- expression --------------------------------------------------------------------
if dge_primary:
    block("expression", "Expression", "product",
          f"{P('design', '~ condition')} · baseline {P('reference_level', 'first level')}",
          [R("samples quantified", fmt(len(mapping)) if mapping else fmt(len(samples))),
           R("contrasts tested", fmt(len(dge_primary["contrasts"]))),
           R("genes tested", fmt(max((c["tested"] for c in dge_primary["contrasts"]), default=0)),
             f"per contrast, at most, on the {primary_map} map"),
           R("differentially expressed", fmt(de_total),
             f"padj < {P('de_alpha', 0.05)} AND |log2FC| > {P('de_lfc', 1)}, summed over "
             f"{len(dge_primary['contrasts'])} contrasts, on the {primary_map} map "
             f"({fmt(de_total_padj_only)} pass padj alone)"),
           R("gene map", f"{len(gene_maps)} maps compared" if len(gene_maps) > 1 else "string-stripped ids",
             "see the Genes tab" if dge_locus else "")],
          "expression")

# ---- the headline figure: whatever this run's main product is ----------------------
if ref_row:
    hero = {"value": fmt(ref_row["transcripts"]), "label": "transcripts in the assembled transcriptome",
            "sub": f"{ref_set} · {fmt(ref_row['genes']) if ref_row.get('genes') else '—'} genes"
                   + (f" · BUSCO {ref_row['busco']['complete']:.1f}% complete" if ref_row.get("busco") else "")}
elif dge_primary:
    hero = {"value": fmt(de_total), "label": "differentially expressed genes",
            "sub": f"padj < 0.05 across {len(dge_primary['contrasts'])} contrasts, {primary_map} map"}
elif total_kept:
    hero = {"value": fmt(total_kept), "label": "host read pairs kept",
            "sub": f"{len(samples)} samples"}
else:
    hero = {"value": fmt(len(samples)) if samples else "—", "label": "samples processed", "sub": ""}

# an end-to-end accounting check the reader can see rather than take on trust
frag = sum(s["fragments"] for s in samples)
accounting = None
if frag and total_kept:
    accounting = {"kept": total_kept, "processed": frag, "ok": frag == total_kept}


# ------------------------------------------------------------------ method stages
def M(label, value):
    return {"label": label, "value": value}


stages = []


# ------------------------------------------------------------------ expression filter
# The filter log is the only place the expressed/rescued split and the DE-eligible count are
# written down. Parsed rather than recomputed so the dashboard cannot disagree with the pipeline.
expr_filter = {}
_ef = first_glob("transcriptome/*_expressed.filter.txt") or first_glob("*_expressed.filter.txt")
if _ef:
    _txt = open(_ef).read()
    def _ef_num(pat):
        m = re.search(pat, _txt)
        return int(m.group(1).replace(",", "")) if m else None
    expr_filter = {
        "rule":      (re.search(r"expression filter:\s*(.+)", _txt).group(1).strip()
                      if re.search(r"expression filter:\s*(.+)", _txt) else ""),
        "kept":      _ef_num(r"transcripts kept\s*:\s*([\d,]+)"),
        "expressed": _ef_num(r"expressed\s*:\s*([\d,]+)"),
        "rescued":   _ef_num(r"rescued by ORF\+hom:\s*([\d,]+)"),
        "dropped":   _ef_num(r"transcripts dropped\s*:\s*([\d,]+)"),
        "eligible":  _ef_num(r"DE-eligible transcripts \(expressed only\):\s*([\d,]+)"),
    }
    _g = re.search(r"genes\s+([\d,]+)\s*->\s*([\d,]+)", _txt)
    if _g:
        expr_filter["genes_before"] = int(_g.group(1).replace(",", ""))
        expr_filter["genes_after"] = int(_g.group(2).replace(",", ""))


def stage(sid, name, group, tool, ran, why, arglist, metrics, out=None):
    stages.append({"id": sid, "name": name, "group": group, "tool": tool,
                   "status": "ran" if ran else "skipped", "why": why,
                   "args": [a for a in arglist if a], "metrics": [m for m in metrics if m["value"] != "—"],
                   "output": out or ""})


stage("qc", "Read trimming", "Reads", "fastp", on("skip_fastp"),
      "Adapter and quality trimming. Sample sheet rows sharing a sample name are lanes of one "
      "library and are concatenated first — gzip members concatenate legally, so no re-compression.",
      [f"args = {P('ext_fastp', '--detect_adapter_for_pe')}"],
      [M("samples", fmt(len(samples))), M("read pairs in", fmt(total_in) if total_in else "—")],
      "qc/fastp/")

stage("rrna", "rRNA depletion", "Reads", "bbduk", on("skip_rrna"),
      "k-mer match against an rRNA database; a pair is dropped when either mate matches. bbduk "
      "rather than SortMeRNA for native paired output and a per-database %rRNA report.",
      [f"k = {P('rrna_kmer', 31)}",
       f"database = {os.path.basename(str(P('rrna_fasta', 'downloaded SortMeRNA set')))}"],
      [M("read pairs removed", fmt(total_rrna) if total_rrna else "—"),
       M("mean % rRNA", f"{sum(s['rrna_pct'] for s in samples) / len(samples):.2f}%" if samples and total_rrna else "—")],
      "qc/rrna/")

stage("decon", "Contaminant removal", "Reads", "minimap2 + seqkit", on("skip_decontam"),
      "Host and contaminant genomes go into ONE index, so each read's primary alignment is its best "
      "placement genome-wide: a conserved host read that also aligns to a contaminant still scores "
      "best on the host and is kept, and reads aligning nowhere are kept too. Only a positive "
      "contaminant placement removes a pair. Done before assembly, because an assembly built from "
      "mixed reads cannot be cleanly separated afterwards.",
      [f"index = host + {P('contaminant_prefix', 'CONTAM_')}-prefixed contaminant genomes",
       "minimap2 -ax sr, primary non-supplementary records only (-F 0x900)",
       f"long-read classifier margin = {P('contaminant_margin', 1.1)}"],
      [M("read pairs removed", fmt(total_contam) if total_contam else "—"),
       M("pooled % removed", f"{100 * total_contam / max(1, total_contam + total_kept):.2f}%" if total_contam else "—"),
       M("long-read isoforms removed",
         fmt(next((r["n"] for r in lr_decon if r["class"] == "removed"), None)) if lr_decon else "—")],
      "decontam/")

def rung(name):
    return next((r for r in ladder if r["reference"] == name), None)


def added(later, earlier):
    a, b = rung(later), rung(earlier)
    return fmt(a["transcripts"] - b["transcripts"]) if a and b else "—"


lr_row = rung("isoseq")
stage("routea", "Long reads on the genome", "Assembly", "pbmm2 + isoseq collapse",
      on("skip_longread") and on("skip_genome_guided"),
      "Spliced alignment of the long reads to the reference genome, then collapse to non-redundant "
      "isoform models. This is the only source that yields real genomic loci and isoform structure — "
      "everything else is layered on top of it.",
      ["pbmm2 align --preset ISOSEQ --sort", "isoseq collapse"],
      [M("transcripts", fmt(lr_row["transcripts"]) if lr_row else "—"),
       M("BUSCO complete", f"{lr_row['busco']['complete']:.1f}%" if lr_row and lr_row["busco"] else "—")],
      "assembly/isoseq/")

stage("routeb", "Long reads, reference-free", "Assembly", "vsearch",
      on("skip_longread") and on("skip_reffree"),
      "Collapses near-identical isoforms without consulting the genome, so loci missing from a draft "
      "assembly survive. vsearch rather than cd-hit-est: the identity definition is stated rather "
      "than implied, where cd-hit measures matches over the shorter sequence without saying so. "
      "The params are still named cdhit_* for continuity with earlier runs.",
      [f"vsearch --cluster_fast --id {P('cdhit_longread_id', 0.99)}", f"prefix = {P('prefix_lrfree', 'PBR_')}"],
      [M("added to the reference", added("isoseq_stringtie_lrfree", "isoseq_stringtie"))],
      "assembly/clustering/")

boost_row = rung("isoseq_stringtie")
stage("stringtie", "Short reads on the genome", "Assembly", "STAR + StringTie", on("skip_stringtie"),
      "Two-pass alignment with no annotation in the index — junctions are discovered de novo, because "
      "a draft genome's annotation would bias exactly the discovery this step exists for. The merged "
      "assembly is compared to the long-read models and only novel loci are added.",
      ["STAR --twopassMode Basic --outSAMstrandField intronMotif",
       f"stringtie {'--rf' if P('strandedness') == 'reverse' else ''}".strip(),
       f"gffcompare class codes kept = {P('stringtie_class_codes', 'u|x|i|y')}",
       f"minimum novel transcript length = {P('min_transcript_len', 300)} bp"],
      [M("reference after this step", fmt(boost_row["transcripts"]) if boost_row else "—"),
       M("added", added("isoseq_stringtie", "isoseq"))],
      "assembly/stringtie/")

fin_row = rung("isoseq_stringtie_trinity")
stage("trinity", "Short reads de novo", "Assembly", "Trinity", on("skip_denovo"),
      "Genome-free assembly across every sample — plus any assembly-only libraries — so loci absent "
      "from the draft genome are recovered. "
      "Its output is audited for residual contamination, collapsed against itself, then reduced to "
      "what the genome-guided reference does not already contain.",
      [f"--SS_lib_type {'RF' if P('strandedness') == 'reverse' else 'FR'}",
       f"cd-hit-est / -est-2d -c {P('cdhit_dedup_id', 0.95)}",
       f"prefix = {P('prefix_denovo', 'TR_')}"],
      [M("added", added("isoseq_stringtie_trinity", "isoseq_stringtie")),
       M("assembly-only libraries", fmt(len(assembly_reads)) if assembly_reads else "—")],
      "assembly/denovo/")

stage("dedup", "Cross-source deduplication", "Assembly", "vsearch", True,
      "The layered set counts one biological locus once per source that recovered it. Deduplication "
      "reduces that; it does not eliminate it — which is what the locus rebuild below addresses.",
      [f"vsearch --cluster_fast --id {P('cdhit_dedup_id', 0.95)}"],
      [M("reference used downstream", fmt(ref_row["transcripts"]) if ref_row else "—"),
       M("BUSCO complete", f"{ref_row['busco']['complete']:.1f}%" if ref_row and ref_row["busco"] else "—")],
      "transcriptome/")

stage("filter", "Expression filter", "Expression", "salmon pass 1 + filter_by_expression.py", True,
      "The ladder applies no expression evidence at any rung, so the deduplicated set carries every "
      "unexpressed fragment the four sources produced. A first salmon pass over that set measures "
      "what the experiment actually saw, and the reference is cut to the survivors. What this cuts "
      "is the reference for QUANTIFICATION AND TESTING — not the assembly. The annotation and the "
      "BUSCO completeness quoted for this transcriptome are both measured upstream, on the full "
      "deduplicated set, so a dropped transcript keeps its annotation row and still counts toward "
      "reported completeness; it is simply never quantified and never tested.",
      [f"expressed: TPM >= {P('expr_filter_tpm', 1)} in >= {P('expr_filter_samples', 3)} libraries",
       "rescue (ORF AND homology): ON — those transcripts are quantified but never tested"
       if P("expr_filter_rescue", False) else
       "rescue (ORF AND homology): OFF — measured on this data to add 0 DE genes and 0.13 pts of "
       "mapping, while costing 4.2 pts of BUSCO duplication",
       "a presence rule, not a mean: a transcript on in one group and off elsewhere survives"],
      [M("DE reference", fmt(expr_filter["kept"]) if expr_filter.get("kept") else "—"),
       M("of those, expressed", fmt(expr_filter["expressed"]) if expr_filter.get("expressed") else "—"),
       M("rescued by ORF+homology", fmt(expr_filter["rescued"]) if expr_filter.get("rescued") else "—"),
       M("dropped", fmt(expr_filter["dropped"]) if expr_filter.get("dropped") else "—")],
      "transcriptome/")

stage("annot", "Functional annotation", "Annotation",
      "TransDecoder · DIAMOND · eggNOG-mapper · InterProScan", on("skip_annotation"),
      "ORFs first, then homology and domains. Column has_ORF is yes/no for every row — the ORF count "
      "is coding_status == \"coding\", not the number of non-empty cells.",
      ["TransDecoder.Predict --single_best_only",
       "diamond blastp -k 1 -e 1e-5 --more-sensitive",
       "eggnog-mapper -m diamond" if on("skip_eggnog") else "eggNOG-mapper: skipped",
       f"interproscan -appl {P('interproscan_appl', '')}" if on("skip_interproscan") else "InterProScan: skipped"],
      [M(a["label"], fmt(a["value"])) for a in annot_summary[:6]],
      "annotation/")

stage("quant", "Quantification", "Expression", "salmon", True,
      "Selective alignment against a decoy-aware index (the genome is the decoy), which stops reads "
      "from genomic or unassembled regions being forced onto a transcript. Library type is "
      "auto-detected and reported, not assumed.",
      ["salmon index -d <genome decoys> -k 31",
       "salmon quant -l A --gcBias --seqBias --dumpEq --numBootstraps 30",
       "run TWICE: once on the deduplicated set to measure expression, once on the filtered set"],
      [M("mean mapping", f"{sum(mapped) / len(mapped):.2f}%" if mapped else "—"),
       M("library type", samples[0]["library_type"] if samples and samples[0]["library_type"] else "—"),
       M("fragments processed", fmt(frag) if frag else "—")],
      "quantification/")

stage("dge", "Differential expression", "Expression", "tximport + DESeq2", on("skip_dge"),
      "Gene-level counts aggregated with tximport, then DESeq2 on the configured design. Only genes "
      "with at least one expressed transcript are tested; where the ORF+homology rescue is enabled "
      "those transcripts stay in the quantification and their reads still count toward their gene, "
      "without themselves being testable. Fold changes are "
      "SHRUNKEN: apeglm for every contrast, releveling and refitting the Wald test where the "
      "contrast is not already a model coefficient, with ashr only as a fallback if that fails. "
      "The raw estimate is kept as log2FoldChange_unshrunken. sig applies both the padj and the "
      "fold-change cut; sig_padj drops the fold-change cut, so it is the larger number.",
      [f"design = {P('design', '~ condition')}",
       f"baseline level = {P('reference_level', 'first level')}",
       f"lfcShrink type = {P('de_shrink', 'apeglm')}",
       f"DE-eligible: genes with >=1 transcript at TPM >= {P('expr_filter_tpm', 1)} in >= {P('expr_filter_samples', 3)} libraries",
       f"counts >= {P('de_min_count', 5)} in >= {P('de_min_samples', 3)} samples",
       f"significance: padj < {P('de_alpha', 0.05)}, |log2FC| > {P('de_lfc', 1)}"],
      [M("contrasts", fmt(len(dge_primary["contrasts"])) if dge_primary else "—"),
       M(f"DE genes, all contrasts ({primary_map} map)", fmt(de_total) if dge_primary else "—"),
       M("padj alone, no fold-change cut", fmt(de_total_padj_only) if dge_primary else "—")],
      map_outputs("dge"))

stage("enrich", "Enrichment", "Expression", "clusterProfiler", on("skip_enrichment"),
      "GO/KEGG over-representation and GSEA against this run's own gene→term map, so no model-organism "
      "database is needed. GO terms come from eggNOG and InterProScan together; KEGG is tested at "
      "pathway level, not on orthologue ids. fgsea's multilevel algorithm is randomised: the seed is "
      "fixed, otherwise the same inputs return a different set of significant terms on every run.",
      [f"seed = {P('random_seed', 42)} (set.seed + clusterProfiler seed=TRUE)",
       "universe = every gene in the model (the DE-eligible set), shared by all contrasts",
       "GSEA ranked on the Wald statistic, minGSSize 10",
       f"tables written in full; significance = p.adjust < {P('de_alpha', 0.05)}"],
      [M(f"contrasts with terms ({primary_map} map)",
         fmt(len(enrichment_by_map[primary_map])) if enrichment_by_map.get(primary_map) else "—")],
      map_outputs("enrichment"))

stage("locus", "Gene identity from genomic loci", "Genes", "minimap2 + bedtools",
      bool(P("run_locus_rebuild", True)) and locus is not None,
      "String-stripped transcript ids apply a different rule per assembly source and reconcile nothing "
      "ACROSS sources, so one locus recovered by three sources counts three times. This places every "
      "transcript on the genome and groups connected components sharing a splice junction. It halves "
      "the split-gene problem and roughly doubles gene fusion — both directions are measured below.",
      [f"grouping = {P('locus_groupby', 'junction')}",
       f"minimum transcript coverage = {P('locus_mincov', 0.8)}",
       f"single-exon attach threshold = {P('locus_se_overlap', 0.5)}"],
      [M("genes before", fmt(locus["before"]) if locus else "—"),
       M("genes after", fmt(locus["after"]) if locus else "—")],
      "transcriptome/locus_rebuild_report.txt")


# ------------------------------------------------------------------ file inventory
inventory = []
for p in sorted(glob.glob(path("*"))):
    base = os.path.basename(p)
    if base.startswith("."):
        continue
    if os.path.isdir(p):
        n = len(glob.glob(os.path.join(p, "*")))
        inventory.append({"name": base + "/", "size": f"{n} files", "kind": "directory"})
    else:
        inventory.append({"name": base, "size": f"{os.path.getsize(p) / 1024:.0f} KB", "kind": "file"})

DATA = {
    "meta": {
        "title": args.title,
        "reference_set": ref_set,
        "run_name": wf.get("run_name", ""),
        "start": wf.get("start", ""),
        "profile": wf.get("profile", ""),
        "revision": wf.get("revision", ""),
        "nextflow": wf.get("nextflow", ""),
        "command_line": wf.get("command_line", ""),
        "outdir": wf.get("outdir", P("outdir", "")),
        "interest_label": interest_label,
        "primary_map": primary_map,
    },
    "hero": hero,
    "datasets": datasets,
    "accounting": accounting,
    "stages": stages,
    "expr_filter": expr_filter,
    "samples": samples,
    "lanes": lanes,
    "assembly_reads": assembly_reads,
    "longread": {"rrna": lr_rrna, "decontam": lr_decon, "samples": lr_samples, "note": lr_note},
    "ladder": ladder,
    "annotation": {"summary": annot_summary, "evidence": annot_evidence,
                   "reference": ref_set,
                   "note": ("Built on the selected ladder rung BEFORE the expression filter, so "
                            "this describes more transcripts than the DE tabs — deliberately. "
                            "Annotation and BUSCO completeness are properties of the ASSEMBLY and "
                            "are measured on the full deduplicated set; the expression filter "
                            "produces the reference for quantification and testing and has no say "
                            "over what is annotated. A transcript with a protein hit, an InterPro "
                            "domain or a KEGG assignment therefore keeps its row here even when "
                            "the experiment never saw it above threshold: that is a statement "
                            "about what this species encodes, not a result. Gene ids here are the "
                            "string-stripped map; the locus and Corset variants are re-keyed "
                            "copies of this same table.")},
    "dge": dge,
    "dge_locus": dge_locus,
    "dge_corset": dge_corset,
    "enrichment_by_map": enrichment_by_map,
    "locus": locus,
    "figures_by_map": figures_by_map,
    "figure_meta": figure_meta,
    "ma_figures": ma_figures,
    "gene_maps": gene_maps,
    "map_flow": map_flow,
    "venn": venn,
    "de_filtering_by_map": de_filtering_by_map,
    "related": [r for r in [
        {"label": "MultiQC report", "href": "../qc/multiqc/multiqc_report.html",
         "note": "per-sample read QC: fastp, FastQC, rRNA, STAR, salmon"} if on("skip_multiqc") else None,
        {"label": f"Figure set, {primary_map} map (PDF + 300 dpi PNG)",
         "href": f"../{os.path.basename(map_dir('figures', primary_map))}/",
         "note": "the publication figures for the primary gene map; the other maps' sets sit beside it"}
        if on("skip_figures") and map_dir("figures", primary_map) else None,
        {"label": "run_summary.txt", "href": "../run_summary.txt",
         "note": "the same headline numbers as plain text, with the file each came from"},
        {"label": "Reference transcriptome", "href": f"../transcriptome/{ref_set}.fasta",
         "note": "the FASTA everything downstream was quantified against"},
    ] if r],
    "params": {k: params[k] for k in sorted(params)} if params else {},
    "inventory": inventory,
}

STYLE = r"""
*,*::before,*::after{box-sizing:border-box}
:root{
  color-scheme:light;
  --plane:#f9f9f7; --surface:#fcfcfb; --raise:#ffffff;
  --ink:#0b0b0b; --ink-2:#52514e; --muted:#898781;
  --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,.10);
  --s1:#2a78d6; --s2:#eb6834; --s3:#1baf7a; --s4:#eda100;
  --up:#e34948; --down:#2a78d6;
  --ramp-1:#2a78d6; --ramp-2:#86b6ef;
  --good:#0ca30c; --warning:#fab219; --serious:#ec835a; --critical:#d03b3b;
  --wash:rgba(42,120,214,.08);
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    color-scheme:dark;
    --plane:#0d0d0d; --surface:#1a1a19; --raise:#222220;
    --ink:#ffffff; --ink-2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,.10);
    --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#c98500;
    --up:#e66767; --down:#3987e5;
    --ramp-1:#3987e5; --ramp-2:#256abf;
    --wash:rgba(57,135,229,.14);
  }
}
:root[data-theme="dark"]{
  color-scheme:dark;
  --plane:#0d0d0d; --surface:#1a1a19; --raise:#222220;
  --ink:#ffffff; --ink-2:#c3c2b7; --muted:#898781;
  --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,.10);
  --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#c98500;
  --up:#e66767; --down:#3987e5;
  --ramp-1:#3987e5; --ramp-2:#256abf;
  --wash:rgba(57,135,229,.14);
}
html,body{margin:0;padding:0}
body{background:var(--plane);color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:15px;line-height:1.5}
a{color:inherit}
.wrap{max-width:1180px;margin:0 auto;padding:0 24px 64px}
header.top{position:sticky;top:0;z-index:20;background:var(--plane);border-bottom:1px solid var(--border)}
header.top .wrap{padding-top:18px;padding-bottom:0}
.title-row{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap}
h1{font-size:20px;font-weight:650;margin:0;letter-spacing:-.01em}
.runmeta{color:var(--muted);font-size:13px}
.spacer{flex:1}
button{font:inherit;color:inherit;background:none;border:0;cursor:pointer}
.ghost{border:1px solid var(--border);border-radius:7px;padding:5px 11px;font-size:13px;color:var(--ink-2)}
.ghost:hover{background:var(--wash)}
.ghost[aria-pressed="true"]{background:var(--wash);color:var(--ink);border-color:var(--s1)}
nav.tabs{display:flex;gap:2px;margin-top:14px;overflow-x:auto}
nav.tabs button{padding:9px 14px;font-size:14px;color:var(--ink-2);border-bottom:2px solid transparent;white-space:nowrap}
nav.tabs button:hover{color:var(--ink)}
nav.tabs button[aria-selected="true"]{color:var(--ink);border-bottom-color:var(--s1);font-weight:600}
section[role="tabpanel"]{padding-top:26px}
section[hidden]{display:none}
h2{font-size:16px;font-weight:650;margin:32px 0 4px;letter-spacing:-.005em}
h2:first-child{margin-top:0}
.sub{color:var(--ink-2);font-size:13.5px;margin:0 0 16px;max-width:78ch}
.card{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:18px 20px;margin-bottom:16px}
.grid{display:grid;gap:14px}
.g2{grid-template-columns:repeat(auto-fit,minmax(330px,1fr))}
.g3{grid-template-columns:repeat(auto-fit,minmax(215px,1fr))}
.hero{display:flex;flex-wrap:wrap;align-items:baseline;gap:16px;padding:22px 24px}
.hero .fig{font-size:52px;font-weight:640;line-height:1;letter-spacing:-.02em}
.hero .lab{font-size:15px;color:var(--ink-2);max-width:32ch}
.hero .sub2{color:var(--muted);font-size:13px;width:100%}
/* dataset panels on the overview */
.dataset{display:flex;flex-direction:column;gap:2px}
.dshead{display:flex;align-items:baseline;gap:10px}
.dshead h3{font-size:15px;font-weight:640;margin:0}
.dshead .kind{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.06em}
.dssub{color:var(--ink-2);font-size:13px;margin:2px 0 10px}
.barrows{margin:2px 0 12px}
.barrow{display:grid;grid-template-columns:minmax(96px,34%) 1fr auto;gap:10px;align-items:center;padding:3px 0}
.barrow .k{color:var(--ink-2);font-size:12.5px}
.barrow .v{font-size:12.5px;font-variant-numeric:tabular-nums;color:var(--ink)}
.barrow .track{height:10px;background:var(--grid);border-radius:5px;overflow:hidden}
.barrow .fill{height:100%;border-radius:0 5px 5px 0;min-width:2px}
dl.facts{display:grid;grid-template-columns:minmax(110px,40%) 1fr;gap:2px 12px;margin:0}
dl.facts dt{color:var(--ink-2);font-size:12.5px;padding:3px 0}
dl.facts dd{margin:0;padding:3px 0;font-size:12.5px}
dl.facts .fv{font-variant-numeric:tabular-nums;font-weight:600}
dl.facts .fn{color:var(--muted);margin-left:7px}
.golink{align-self:flex-start;margin-top:12px}
.tile .label{color:var(--ink-2);font-size:12.5px}
.tile .value{font-size:26px;font-weight:640;letter-spacing:-.01em;margin-top:3px}
.tile .note{color:var(--muted);font-size:12px;margin-top:3px}
.note-line{color:var(--muted);font-size:12.5px;margin-top:10px}
.check{display:flex;gap:9px;align-items:flex-start;font-size:13.5px;color:var(--ink-2)}
.dot{width:9px;height:9px;border-radius:50%;flex:0 0 auto;margin-top:5px}
.figcap{display:flex;align-items:baseline;gap:10px;margin-bottom:2px}
.figcap h3{font-size:14.5px;font-weight:620;margin:0}
.figcap .why{color:var(--muted);font-size:12.5px}
.legend{display:flex;gap:14px;flex-wrap:wrap;margin:8px 0 6px;font-size:12.5px;color:var(--ink-2)}
.legend span{display:inline-flex;align-items:center;gap:6px}
.key{width:11px;height:11px;border-radius:3px;display:inline-block}
.key.line{height:3px;width:14px;border-radius:2px}
svg{display:block;width:100%;overflow:visible}
svg text{font-family:inherit}
.tick{fill:var(--muted);font-size:11px;font-variant-numeric:tabular-nums}
.axlab{fill:var(--ink-2);font-size:11.5px}
.vallab{fill:var(--ink-2);font-size:11.5px;font-variant-numeric:tabular-nums}
.mark:focus{outline:2px solid var(--s1);outline-offset:2px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:6px 10px;border-bottom:1px solid var(--border);vertical-align:top}
th{color:var(--ink-2);font-weight:600;font-size:12px;white-space:nowrap}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.scroll{overflow-x:auto;max-height:460px;overflow-y:auto}
.stage{border:1px solid var(--border);border-radius:11px;background:var(--surface);margin-bottom:10px}
.stage>summary{list-style:none;cursor:pointer;padding:13px 16px;display:flex;gap:12px;align-items:center}
.stage>summary::-webkit-details-marker{display:none}
.stage summary:hover{background:var(--wash)}
.stage .nm{font-weight:600;font-size:14.5px}
.stage .tool{color:var(--muted);font-size:12.5px}
.stage .body{padding:0 16px 16px 16px;border-top:1px solid var(--border);margin-top:2px}
.stage .why{color:var(--ink-2);font-size:13.5px;margin:12px 0;max-width:80ch}
.badge{font-size:11px;padding:2px 8px;border-radius:20px;border:1px solid var(--border);color:var(--ink-2)}
.badge.ran{border-color:color-mix(in srgb,var(--good) 45%,transparent)}
.badge.skipped{color:var(--muted)}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px}
.args{margin:0;padding:0;list-style:none}
.args li{padding:5px 9px;background:var(--plane);border:1px solid var(--border);border-radius:6px;
  margin-bottom:5px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:var(--ink-2);
  overflow-x:auto;white-space:pre}
.metrics{display:flex;flex-wrap:wrap;gap:18px;margin-top:10px}
.metrics div{font-size:12.5px;color:var(--muted)}
.metrics b{display:block;font-size:16px;color:var(--ink);font-weight:620}
.rail{display:flex;gap:6px;flex-wrap:wrap;margin:6px 0 2px}
.rail button{border:1px solid var(--border);border-radius:8px;padding:8px 11px;font-size:12.5px;
  color:var(--ink-2);background:var(--surface);display:flex;gap:8px;align-items:center}
.rail button:hover{border-color:var(--s1);color:var(--ink)}
.rail .arrow{color:var(--muted)}
.rail .n{font-variant-numeric:tabular-nums;color:var(--ink);font-weight:600}
.controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:14px}
select,input[type=search],input[type=number]{font:inherit;font-size:13px;color:var(--ink);background:var(--surface);
  border:1px solid var(--border);border-radius:7px;padding:6px 9px}
input[type=search]{min-width:230px}
input[type=number]{width:82px;font-variant-numeric:tabular-nums}
.controls .why{color:var(--ink-2)}
.tip{position:fixed;z-index:60;pointer-events:none;background:var(--raise);color:var(--ink);
  border:1px solid var(--border);border-radius:9px;padding:9px 11px;font-size:12.5px;
  box-shadow:0 6px 22px rgba(0,0,0,.16);max-width:290px;opacity:0;transition:opacity .09s}
.tip .tv{font-weight:640;font-size:14px}
.tip .tr{display:flex;gap:8px;align-items:baseline;margin-top:3px}
.tip .tk{width:12px;height:3px;border-radius:2px;flex:0 0 auto}
.tip .tn{color:var(--ink-2)}
.caveat{border-left:3px solid var(--warning);padding:2px 0 2px 13px;color:var(--ink-2);font-size:13.5px;margin:10px 0}
.warnbox{border-left:3px solid var(--critical)}
.empty{color:var(--muted);font-size:13.5px;padding:10px 0}
@media (max-width:640px){.hero .fig{font-size:40px}.wrap{padding:0 14px 48px}}
@media print{nav.tabs,.ghost{display:none}section[hidden]{display:block!important}}
"""

SCRIPT = r"""
const DATA = JSON.parse(document.getElementById('run-data').textContent);
const SVGNS = 'http://www.w3.org/2000/svg';
const W = 760;                      // chart coordinate width; the SVG scales to its card

/* ---------- tiny DOM helpers. Labels come from tool output: always textContent. ---------- */
function h(tag, attrs, kids) {
  const n = document.createElement(tag);
  for (const k in (attrs || {})) {
    if (k === 'text') n.textContent = attrs[k];
    else if (k === 'class') n.className = attrs[k];
    else if (k.startsWith('on')) n.addEventListener(k.slice(2), attrs[k]);
    else if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
  }
  (kids || []).forEach(c => c && n.appendChild(c));
  return n;
}
function s(tag, attrs, text) {
  const n = document.createElementNS(SVGNS, tag);
  for (const k in (attrs || {})) if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
  if (text !== undefined) n.textContent = text;
  return n;
}
const padjOf = p => Math.pow(10, -p[1]);
const isSig = (p, padjT, lfcT) => padjOf(p) < padjT && Math.abs(p[0]) >= lfcT;
const fmt = n => (n === null || n === undefined || Number.isNaN(n)) ? '—'
  : (Math.abs(n) >= 1000 ? Math.round(n).toLocaleString() : (Number.isInteger(n) ? String(n) : n.toFixed(2)));
const pct = n => (n === null || n === undefined) ? '—' : n.toFixed(2) + '%';
const clean = t => String(t === null || t === undefined ? '' : t);

/* rounded on the data end only, square at the baseline */
function barPath(x, y, w, hgt, r) {
  r = Math.max(0, Math.min(r, w, hgt / 2));
  if (w <= 0) return `M${x} ${y}h0`;
  return `M${x} ${y}h${w - r}a${r} ${r} 0 0 1 ${r} ${r}v${hgt - 2 * r}a${r} ${r} 0 0 1 ${-r} ${r}H${x}z`;
}
/* a label never overflows its gutter: measure-by-proxy, truncate, keep the full text in the tooltip */
const CHAR_W = 5.6;                 /* ~11.5px system sans, averaged */
function fit(text, px) {
  const t = String(text === null || text === undefined ? '' : text);
  const max = Math.max(3, Math.floor(px / CHAR_W));
  return t.length <= max ? t : t.slice(0, max - 1) + '\u2026';
}
function ticks(max, n) {
  if (max <= 0) return [0];
  const raw = max / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(v => v >= raw) || mag * 10;
  const out = []; for (let v = 0; v <= max * 1.0001; v += step) out.push(v);
  return out;
}

/* ---------- one shared tooltip: enhances, never gates (every value is in the table view) ---------- */
const tip = h('div', { class: 'tip', role: 'status' });
document.body.appendChild(tip);
function showTip(evt, title, rows) {
  tip.textContent = '';
  tip.appendChild(h('div', { class: 'tv', text: clean(title) }));
  (rows || []).forEach(r => {
    const line = h('div', { class: 'tr' });
    if (r.color) { const k = h('span', { class: 'tk' }); k.style.background = r.color; line.appendChild(k); }
    line.appendChild(h('span', { class: 'tn', text: clean(r.label) }));
    line.appendChild(h('b', { text: clean(r.value) }));
    tip.appendChild(line);
  });
  tip.style.opacity = '1';
  moveTip(evt);
}
function moveTip(evt) {
  const pad = 14, r = tip.getBoundingClientRect();
  let x = (evt.clientX || 0) + pad, y = (evt.clientY || 0) + pad;
  if (x + r.width > innerWidth - 8) x = (evt.clientX || 0) - r.width - pad;
  if (y + r.height > innerHeight - 8) y = (evt.clientY || 0) - r.height - pad;
  tip.style.left = Math.max(8, x) + 'px'; tip.style.top = Math.max(8, y) + 'px';
}
function hideTip() { tip.style.opacity = '0'; }
function hoverable(node, title, rows) {
  node.classList.add('mark');
  node.setAttribute('tabindex', '0');
  node.addEventListener('pointerenter', e => showTip(e, title, rows));
  node.addEventListener('pointermove', moveTip);
  node.addEventListener('pointerleave', hideTip);
  node.addEventListener('focus', e => {
    const b = node.getBoundingClientRect();
    showTip({ clientX: b.left + b.width / 2, clientY: b.top }, title, rows);
  });
  node.addEventListener('blur', hideTip);
}

/* ---------- figure card: chart + its table twin ---------- */
function card(opts) {
  const body = h('div');
  const chart = h('div');
  // a card may be table-only: appendChild(undefined) throws, and the whole tab then fails to build
  const hasChart = !!opts.chart;
  if (hasChart) chart.appendChild(opts.chart);
  const table = h('div', { class: 'scroll' });
  // opts.table is a spec for buildTable; opts.tableNode is an already-built element, for a table
  // assembled from several sources rather than from one row set
  const hasTable = !!(opts.table || opts.tableNode);
  if (opts.tableNode) table.appendChild(opts.tableNode);
  else if (opts.table) table.appendChild(buildTable(opts.table));
  // the table hides behind a toggle only when there is a chart to toggle against
  if (hasTable && hasChart) table.setAttribute('hidden', '');
  const toggle = h('button', {
    class: 'ghost', 'aria-pressed': 'false', text: 'Table',
    onclick: () => {
      const showTable = chart.hasAttribute('hidden');
      if (showTable) { chart.removeAttribute('hidden'); table.setAttribute('hidden', ''); toggle.setAttribute('aria-pressed', 'false'); }
      else { chart.setAttribute('hidden', ''); table.removeAttribute('hidden'); toggle.setAttribute('aria-pressed', 'true'); }
    }
  });
  const cap = h('div', { class: 'figcap' }, [
    h('h3', { text: opts.title }),
    opts.why ? h('span', { class: 'why', text: opts.why }) : null,
    h('span', { class: 'spacer' }),
    (hasTable && hasChart) ? toggle : null
  ]);
  body.appendChild(cap);
  if (opts.legend) body.appendChild(legendRow(opts.legend));
  if (hasChart) body.appendChild(chart);
  if (hasTable) body.appendChild(table);
  if (opts.note) body.appendChild(h('p', { class: 'note-line', text: opts.note }));
  return h('div', { class: 'card' }, [body]);
}
function legendRow(items) {
  const row = h('div', { class: 'legend' });
  items.forEach(it => {
    const key = h('span', { class: it.line ? 'key line' : 'key' });
    key.style.background = it.color;
    row.appendChild(h('span', {}, [key, h('span', { text: it.label })]));
  });
  return row;
}
function buildTable(spec) {
  const t = h('table');
  const tr = h('tr');
  spec.cols.forEach(c => tr.appendChild(h('th', { class: c.num ? 'num' : '', text: c.label })));
  t.appendChild(h('thead', {}, [tr]));
  const tb = h('tbody');
  spec.rows.forEach(r => {
    const row = h('tr');
    spec.cols.forEach(c => row.appendChild(h('td', { class: c.num ? 'num' : '', text: clean(r[c.key]) })));
    tb.appendChild(row);
  });
  t.appendChild(tb);
  return t;
}

/* ---------- horizontal bars, one series (magnitude) ---------- */
function barH(spec) {
  const rows = spec.rows, gutter = spec.gutter || 150, rowH = spec.rowH || 26, barH_ = 16;
  const right = spec.right || 74, plot = W - gutter - right;
  const max = spec.max || Math.max(1, ...rows.map(r => r.value));
  const H = rows.length * rowH + 34;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': spec.title || '' });
  ticks(max, 4).forEach(v => {
    const x = gutter + plot * v / max;
    svg.appendChild(s('line', { x1: x, x2: x, y1: 0, y2: rows.length * rowH + 4, style: 'stroke:var(--grid);stroke-width:1' }));
    const tk = s('text', { x: x, y: rows.length * rowH + 20, 'text-anchor': 'middle', class: 'tick' },
      spec.tickFmt ? spec.tickFmt(v) : fmt(v));
    svg.appendChild(tk);
  });
  svg.appendChild(s('line', { x1: gutter, x2: gutter, y1: 0, y2: rows.length * rowH + 4, style: 'stroke:var(--axis);stroke-width:1' }));
  rows.forEach((r, idx) => {
    const y = idx * rowH + (rowH - barH_) / 2;
    const w = Math.max(0, plot * (r.value / max));
    const g = s('g', {});
    const shown = fit(r.label, gutter - 12);
    g.appendChild(s('text', { x: gutter - 8, y: y + barH_ - 3, 'text-anchor': 'end', class: 'axlab' }, shown));
    const p = s('path', { d: barPath(gutter, y, w, barH_, 4), style: `fill:${r.color || spec.color || 'var(--s1)'}` });
    g.appendChild(p);
    g.appendChild(s('text', { x: gutter + w + 7, y: y + barH_ - 3, class: 'vallab' },
      spec.valFmt ? spec.valFmt(r) : fmt(r.value)));
    const rows_ = (r.tip || [{ label: spec.measure || 'value', value: spec.valFmt ? spec.valFmt(r) : fmt(r.value) }]);
    hoverable(g, r.label, rows_);
    svg.appendChild(g);
  });
  return svg;
}

/* ---------- horizontal stacked bars (part-to-whole per sample) ---------- */
function stackedBarH(spec) {
  const rows = spec.rows, gutter = spec.gutter || 150, rowH = 24, barH_ = 15, right = 66;
  const plot = W - gutter - right;
  const max = Math.max(1, ...rows.map(r => r.parts.reduce((a, p) => a + p.value, 0)));
  const H = rows.length * rowH + 34;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': spec.title || '' });
  ticks(max, 4).forEach(v => {
    const x = gutter + plot * v / max;
    svg.appendChild(s('line', { x1: x, x2: x, y1: 0, y2: rows.length * rowH + 4, style: 'stroke:var(--grid);stroke-width:1' }));
    svg.appendChild(s('text', { x: x, y: rows.length * rowH + 20, 'text-anchor': 'middle', class: 'tick' },
      spec.tickFmt ? spec.tickFmt(v) : fmt(v)));
  });
  svg.appendChild(s('line', { x1: gutter, x2: gutter, y1: 0, y2: rows.length * rowH + 4, style: 'stroke:var(--axis);stroke-width:1' }));
  rows.forEach((r, idx) => {
    const y = idx * rowH + (rowH - barH_) / 2;
    let x = gutter;
    const total = r.parts.reduce((a, p) => a + p.value, 0);
    const g = s('g', {});
    g.appendChild(s('text', { x: gutter - 8, y: y + barH_ - 3, 'text-anchor': 'end', class: 'axlab' }, fit(r.label, gutter - 12)));
    r.parts.forEach((p, k) => {
      const w = plot * (p.value / max);
      if (w > 0.4) {
        const last = k === r.parts.length - 1;
        const d = last ? barPath(x, y, w, barH_, 4) : `M${x} ${y}h${w}v${barH_}h${-w}z`;
        g.appendChild(s('path', { d: d, style: `fill:${p.color}` }));
      }
      x += w + 2;                 /* 2px surface gap does the separating, not a stroke */
    });
    g.appendChild(s('text', { x: gutter + plot * (total / max) + 8, y: y + barH_ - 3, class: 'vallab' },
      spec.endLabel ? spec.endLabel(r) : ''));
    hoverable(g, r.label, r.parts.map(p => ({ color: p.color, label: p.label, value: fmt(p.value) }))
      .concat(r.extraTip || []));
    svg.appendChild(g);
  });
  return svg;
}

/* ---------- meter rows: a filled share against a same-ramp track ---------- */
function meterRows(spec) {
  const rows = spec.rows, gutter = spec.gutter || 150, rowH = 26, barH_ = 15, right = 70;
  const plot = W - gutter - right, H = rows.length * rowH + 10;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': spec.title || '' });
  rows.forEach((r, idx) => {
    const y = idx * rowH + (rowH - barH_) / 2;
    const g = s('g', {});
    g.appendChild(s('text', { x: gutter - 8, y: y + barH_ - 3, 'text-anchor': 'end', class: 'axlab' }, fit(r.label, gutter - 12)));
    g.appendChild(s('path', { d: barPath(gutter, y, plot, barH_, 4), style: 'fill:var(--grid)' }));
    let x = gutter;
    r.parts.forEach(p => {
      const w = plot * p.value / 100;
      if (w > 0.4) g.appendChild(s('path', { d: barPath(x, y, w, barH_, 4), style: `fill:${p.color}` }));
      x += w + 2;
    });
    g.appendChild(s('text', { x: gutter + plot + 8, y: y + barH_ - 3, class: 'vallab' }, r.endLabel));
    hoverable(g, r.label, r.tip);
    svg.appendChild(g);
  });
  return svg;
}

/* ---------- scatter with a nearest-point hover layer (dots are too small to aim at) ---------- */
function scatter(spec) {
  const H = spec.height || 330, padL = 54, padR = 16, padT = 12, padB = 38;
  const pw = W - padL - padR, ph = H - padT - padB;
  const xs = spec.points.map(p => p[0]), ys = spec.points.map(p => p[1]);
  const x0 = spec.xMin !== undefined ? spec.xMin : Math.min(...xs), x1 = spec.xMax !== undefined ? spec.xMax : Math.max(...xs);
  const y0 = spec.yMin !== undefined ? spec.yMin : Math.min(...ys), y1 = spec.yMax !== undefined ? spec.yMax : Math.max(...ys);
  const sx = v => padL + pw * (v - x0) / ((x1 - x0) || 1);
  const sy = v => padT + ph - ph * (v - y0) / ((y1 - y0) || 1);
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': spec.title || '' });
  ticks(y1, 4).forEach(v => {
    if (v < y0) return;
    svg.appendChild(s('line', { x1: padL, x2: W - padR, y1: sy(v), y2: sy(v), style: 'stroke:var(--grid);stroke-width:1' }));
    svg.appendChild(s('text', { x: padL - 8, y: sy(v) + 4, 'text-anchor': 'end', class: 'tick' }, fmt(v)));
  });
  const xt = ticks(Math.max(Math.abs(x0), Math.abs(x1)), 3);
  xt.concat(xt.filter(v => v > 0).map(v => -v)).forEach(v => {
    if (v < x0 || v > x1) return;
    svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: padT, y2: padT + ph, style: 'stroke:var(--grid);stroke-width:1' }));
    svg.appendChild(s('text', { x: sx(v), y: H - 16, 'text-anchor': 'middle', class: 'tick' }, fmt(v)));
  });
  (spec.vlines || []).forEach(v => {
    if (v.x < x0 || v.x > x1) return;
    svg.appendChild(s('line', { x1: sx(v.x), x2: sx(v.x), y1: padT, y2: padT + ph, style: 'stroke:var(--axis);stroke-width:1' }));
  });
  (spec.hlines || []).forEach(v => {
    svg.appendChild(s('line', { x1: padL, x2: W - padR, y1: sy(v.y), y2: sy(v.y), style: 'stroke:var(--axis);stroke-width:1' }));
    if (v.label) svg.appendChild(s('text', { x: W - padR, y: sy(v.y) - 5, 'text-anchor': 'end', class: 'tick' }, v.label));
  });
  svg.appendChild(s('text', { x: padL + pw / 2, y: H - 2, 'text-anchor': 'middle', class: 'axlab' }, spec.xLabel || ''));
  svg.appendChild(s('text', { x: 12, y: padT + ph / 2, class: 'axlab', transform: `rotate(-90 12 ${padT + ph / 2})`, 'text-anchor': 'middle' }, spec.yLabel || ''));
  const pts = [];
  spec.points.forEach(p => {
    const cls = spec.classOf(p);
    const c = s('circle', { cx: sx(p[0]), cy: sy(p[1]), r: cls.r || 2.4, style: `fill:${cls.color};opacity:${cls.opacity || 1}` });
    if (cls.ring) c.setAttribute('style', `fill:${cls.color};stroke:var(--surface);stroke-width:2`);
    svg.appendChild(c);
    pts.push({ x: sx(p[0]), y: sy(p[1]), p: p, cls: cls });
  });
  const hit = s('rect', { x: padL, y: padT, width: pw, height: ph, style: 'fill:transparent' });
  const focus = s('circle', { r: 5.5, style: 'fill:none;stroke:var(--ink);stroke-width:1.5;opacity:0' });
  svg.appendChild(hit); svg.appendChild(focus);
  hit.addEventListener('pointermove', evt => {
    const box = svg.getBoundingClientRect(), k = W / box.width;
    const mx = (evt.clientX - box.left) * k, my = (evt.clientY - box.top) * k;
    let best = null, bd = 1e9;
    pts.forEach(q => { const d = (q.x - mx) ** 2 + (q.y - my) ** 2; if (d < bd) { bd = d; best = q; } });
    if (best && bd < 900) {
      focus.setAttribute('cx', best.x); focus.setAttribute('cy', best.y); focus.style.opacity = '1';
      showTip(evt, spec.tipTitle(best.p), spec.tipRows(best.p));
    } else { focus.style.opacity = '0'; hideTip(); }
  });
  hit.addEventListener('pointerleave', () => { focus.style.opacity = '0'; hideTip(); });
  return svg;
}
"""

SCRIPT2 = r"""
/* ================================ sections ================================ */
const S1 = 'var(--s1)', S2 = 'var(--s2)', S3 = 'var(--s3)', S4 = 'var(--s4)';
const UP = 'var(--up)', DOWN = 'var(--down)', MUTED = 'var(--muted)';

/* magnitude rows in plain HTML, so the text stays at reading size inside a narrow card */
function barRows(steps, opts) {
  opts = opts || {};
  const max = Math.max(1, ...steps.map(x => x.value));
  const wrap = h('div', { class: 'barrows' });
  steps.forEach(x => {
    const track = h('div', { class: 'track' });
    const fill = h('div', { class: 'fill' });
    fill.style.width = (100 * x.value / max).toFixed(2) + '%';
    fill.style.background = opts.color || 'var(--s1)';
    track.appendChild(fill);
    const row = h('div', { class: 'barrow' }, [
      h('span', { class: 'k', text: x.label }),
      track,
      h('span', { class: 'v', text: fmt(x.value) })
    ]);
    if (x.note) hoverable(row, x.label, [{ label: x.note, value: fmt(x.value) }]);
    wrap.appendChild(row);
  });
  return wrap;
}

function datasetCard(d) {
  const card = h('div', { class: 'card dataset' });
  card.appendChild(h('div', { class: 'dshead' }, [
    h('h3', { text: d.title }),
    h('span', { class: 'kind', text: d.kind })
  ]));
  if (d.subtitle) card.appendChild(h('p', { class: 'dssub', text: d.subtitle }));
  if ((d.funnel || []).length) card.appendChild(barRows(d.funnel));
  if ((d.rows || []).length) {
    const dl = h('dl', { class: 'facts' });
    d.rows.forEach(r => {
      dl.appendChild(h('dt', { text: r.label }));
      dl.appendChild(h('dd', {}, [
        h('span', { class: 'fv', text: r.value }),
        r.note ? h('span', { class: 'fn', text: r.note }) : null
      ]));
    });
    card.appendChild(dl);
  }
  if (d.note) card.appendChild(h('p', { class: 'note-line', text: d.note }));
  if (d.tab) {
    const label = (TABS.find(t => t[0] === d.tab) || [, d.tab])[1];
    card.appendChild(h('button', { class: 'ghost golink', text: label + ' →', onclick: () => showTab(d.tab) }));
  }
  return card;
}

function sectionOverview() {
  const out = h('div');
  const hero = DATA.hero || {};
  out.appendChild(h('div', { class: 'card hero' }, [
    h('div', { class: 'fig', text: hero.value || '—' }),
    h('div', { class: 'lab', text: hero.label || '' }),
    h('div', { class: 'sub2', text: hero.sub || '' })
  ]));

  /* one block per dataset the run touched - inputs first, then what was built from them */
  const D = DATA.datasets || [];
  if (D.length) {
    const inputs = D.filter(d => d.kind !== 'product');
    const products = D.filter(d => d.kind === 'product');
    if (inputs.length) {
      out.appendChild(h('h2', { text: 'Data in' }));
      out.appendChild(h('p', { class: 'sub', text: 'What the run consumed, and what survived each filter. Only what this run actually had is listed.' }));
      const g = h('div', { class: 'grid g2' });
      inputs.forEach(d => g.appendChild(datasetCard(d)));
      out.appendChild(g);
    }
    if (products.length) {
      out.appendChild(h('h2', { text: 'Results out' }));
      out.appendChild(h('p', { class: 'sub', text: 'What the run produced. Each panel links to the tab that shows how it was measured.' }));
      const g = h('div', { class: 'grid g2' });
      products.forEach(d => g.appendChild(datasetCard(d)));
      out.appendChild(g);
    }
  }

  const a = DATA.accounting;
  if (a) {
    const dot = h('span', { class: 'dot' });
    dot.style.background = a.ok ? 'var(--good)' : 'var(--critical)';
    out.appendChild(h('div', { class: 'card' }, [
      h('div', { class: 'figcap' }, [h('h3', { text: 'End-to-end read accounting' })]),
      h('div', { class: 'check' }, [dot, h('span', {
        text: a.ok
          ? `${fmt(a.processed)} fragments quantified = ${fmt(a.kept)} read pairs kept after contaminant removal. Nothing was dropped or double-counted between the two ends of the run.`
          : `${fmt(a.processed)} fragments quantified vs ${fmt(a.kept)} read pairs kept after contaminant removal — these should match. A difference means samples were dropped or duplicated in the fan-out.`
      })])
    ]));
  }

  /* the pipeline as a rail of its own numbers - click a stage to read its method */
  const rail = h('div', { class: 'rail' });
  (DATA.stages || []).forEach((st, idx) => {
    if (idx) rail.appendChild(h('span', { class: 'arrow', text: '→' }));
    const m = st.metrics && st.metrics[0];
    const b = h('button', { title: st.tool, onclick: () => { showTab('methods'); openStage(st.id); } }, [
      h('span', { text: st.name }),
      m ? h('span', { class: 'n', text: m.value }) : null
    ]);
    if (st.status === 'skipped') b.style.opacity = '.5';
    rail.appendChild(b);
  });
  out.appendChild(h('h2', { text: 'How it ran' }));
  out.appendChild(h('div', { class: 'card' }, [
    h('div', { class: 'figcap' }, [h('h3', { text: 'The run, stage by stage' }),
      h('span', { class: 'why', text: 'each stage shows what it produced — open one for the method behind it' })]),
    rail
  ]));

  const caveats = [];
  if (DATA.dge) caveats.push('Fold changes are shrunken with apeglm, so the volcano x-axis reflects evidence rather than ratio size. The unshrunken estimate is kept in the CSV as log2FoldChange_unshrunken.');
  if (DATA.dge) caveats.push('The significant-gene count carries no fold-change filter, so it exceeds up + down; DESeq2 independent filtering leaves padj = NA for genes that were never tested.');
  if (DATA.gene_maps && DATA.gene_maps.length > 1)
    caveats.push(`${DATA.gene_maps.length} gene maps were built and none is unambiguously right — they agree on a minority of the transcripts they implicate (see "Where the maps agree" in the Genes tab). Quote a range and name the map, not a single gene count.`);
  else if (DATA.locus) caveats.push('Two gene maps exist and neither is unambiguously right — see the Genes tab.');
  const used = (DATA.ladder || []).find(r => r.reference === DATA.meta.reference_set);
  if (used && used.busco && used.busco.dup > 20)
    caveats.push(`BUSCO duplication is ${used.busco.dup.toFixed(1)}% on ${used.reference}, the reference everything was quantified against: on a de novo catalogue that usually means one locus is represented several times, not that the genome is duplicated.`);
  if ((DATA.longread.samples || []).length && DATA.longread.note) caveats.push(DATA.longread.note);
  if (caveats.length) {
    const box = h('div', { class: 'card' }, [h('div', { class: 'figcap' }, [h('h3', { text: 'Read the numbers with these in mind' })])]);
    caveats.forEach(c => box.appendChild(h('p', { class: 'caveat', text: c })));
    out.appendChild(box);
  }
  return out;
}


/* The shape of the run: what feeds what, and where each reported number is produced.
   Drawn rather than described because the things readers get wrong here are all structural: the
   annotation and the reported BUSCO completeness are measured on the FULL assembly, above the
   expression filter, which cuts only the reference for quantification and testing; and the three
   gene maps are three re-aggregations of ONE quantification, not three analyses, so they are drawn
   as parallel branches with their own DESeq2 and their own enrichment rather than folded into a
   single box that implies one answer. Boxes are numbered in execution order so prose can point at
   one. */
// step numbers of the boxes the prose refers to; set when the diagram is drawn (an optional input
// box shifts them)
let FLOW_N = { annot: 8, filt: 10 };
function flowDiagram() {
  const W = 780, H = 706;
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label',
    'Flow: reads are trimmed, rRNA-depleted and decontaminated, then assembled by four routes ' +
    '(any assembly-only libraries join the two short-read routes and nothing else) ' +
    'into a reference ladder; the ladder is deduplicated, and that full set is what gets ' +
    'annotated and scored for completeness. A first salmon pass measures expression, the ' +
    'expression filter cuts the reference down to what the experiment saw, and a second salmon ' +
    'pass quantifies it. The run then splits into three parallel branches — one per gene map — ' +
    'each with its own DESeq2 test and its own GO/KEGG enrichment; one is marked primary and the ' +
    'other two are sensitivity arms. Boxes are numbered in execution order.');

  const mk = (tag, attrs, text) => {
    const e = document.createElementNS(NS, tag);
    Object.entries(attrs).forEach(([k, v]) => { if (v !== null) e.setAttribute(k, v); });
    if (text !== undefined) e.textContent = text;
    return e;
  };
  const defs = mk('defs', {});
  const marker = mk('marker', { id: 'fl-arrow', viewBox: '0 0 10 10', refX: 9, refY: 5,
                                markerWidth: 6, markerHeight: 6, orient: 'auto-start-reverse' });
  marker.appendChild(mk('path', { d: 'M 0 0 L 10 5 L 0 10 z', fill: 'currentColor', 'fill-opacity': 0.55 }));
  defs.appendChild(marker);
  svg.appendChild(defs);

  /* Boxes are numbered as they are created, which is execution order, so prose can name "step 11"
     instead of describing a box. opts.lines takes several sub-lines — the three DE branches carry
     their own counts — and opts.primary draws the heavier border on the map dge_primary_map names. */
  let step = 0;
  const box = (x, y, w, hgt, title, colour, opts) => {
    const o = opts || {};
    const lines = (o.lines || []).filter(t => t);
    const g = mk('g', {});
    g.appendChild(mk('rect', { x, y, width: w, height: hgt, rx: 5,
      fill: colour, 'fill-opacity': o.dashed ? 0.06 : ((o.primary || o.emphasis) ? 0.2 : 0.14),
      stroke: colour, 'stroke-width': o.dashed ? 1 : ((o.primary || o.emphasis) ? 2.6 : 1.5),
      'stroke-dasharray': o.dashed ? '5 3' : null }));
    step += 1;
    g.appendChild(mk('circle', { cx: x + 14, cy: y + 14, r: 9, fill: colour, 'fill-opacity': 0.92 }));
    g.appendChild(mk('text', { x: x + 14, y: y + 17.6, 'text-anchor': 'middle', 'font-size': 10,
      'font-weight': 700, fill: '#fff' }, String(step)));
    // the title is centred in the space right of the step circle (x+24 .. x+w), not over it
    g.appendChild(mk('text', { x: x + (w + 24) / 2, y: y + (lines.length ? 19 : hgt / 2 + 4),
      'text-anchor': 'middle', 'font-size': o.titleSize || 12.5, 'font-weight': 600, fill: 'currentColor' }, title));
    lines.forEach((t, k) => g.appendChild(mk('text', { x: x + w / 2, y: y + 35 + k * 13.5,
      'text-anchor': 'middle', 'font-size': 10.5, fill: 'currentColor',
      'fill-opacity': k === 0 ? 0.72 : 0.58 }, t)));
    if (o.primary) g.appendChild(mk('text', { x: x + w - 9, y: y + 17.5, 'text-anchor': 'end',
      'font-size': 8.5, 'font-weight': 700, fill: colour, 'letter-spacing': 0.6 }, 'PRIMARY'));
    svg.appendChild(g);
    return { x, y, w, h: hgt, cx: x + w / 2, cy: y + hgt / 2, bottom: y + hgt, right: x + w, n: step };
  };
  const arrow = (a, b, label) => {
    svg.appendChild(mk('path', {
      d: `M ${a.cx} ${a.bottom} L ${b.cx} ${b.y - 2}`,
      stroke: 'currentColor', 'stroke-opacity': 0.45, 'stroke-width': 1.4,
      fill: 'none', 'marker-end': 'url(#fl-arrow)' }));
    if (label) svg.appendChild(mk('text', { x: (a.cx + b.cx) / 2 + 6, y: (a.bottom + b.y) / 2 + 3,
      'font-size': 10, fill: 'currentColor', 'fill-opacity': 0.65 }, label));
  };
  // midOverride: a long edge that would otherwise route its horizontal leg straight through an
  // unrelated box is dropped into the gap between them instead.
  const elbow = (a, b, label, midOverride) => {
    const midY = midOverride !== undefined ? midOverride : (a.bottom + b.y) / 2;
    svg.appendChild(mk('path', {
      d: `M ${a.cx} ${a.bottom} V ${midY} H ${b.cx} V ${b.y - 2}`,
      stroke: 'currentColor', 'stroke-opacity': 0.45, 'stroke-width': 1.4,
      fill: 'none', 'marker-end': 'url(#fl-arrow)' }));
    if (label) svg.appendChild(mk('text', { x: b.cx + 8, y: midY - 4, 'font-size': 10,
      fill: 'currentColor', 'fill-opacity': 0.65 }, label));
  };
  /* The annotation supplies terms to all three enrichment boxes. Three long diagonals across the
     middle of the figure would be unreadable, so it runs down the left gutter as one dashed rail
     and drops into each branch — dashed because it carries labels, not counts. */
  const rail = (a, targets, gutterX, dropY, railY, label) => {
    if (!targets.length) return;
    const far = Math.max(...targets.map(t => t.cx));
    svg.appendChild(mk('path', {
      d: `M ${a.cx} ${a.bottom} V ${dropY} H ${gutterX} V ${railY} H ${far}`,
      stroke: 'currentColor', 'stroke-opacity': 0.4, 'stroke-width': 1.4,
      'stroke-dasharray': '4 3', fill: 'none' }));
    targets.forEach(t => svg.appendChild(mk('path', {
      d: `M ${t.cx} ${railY} V ${t.y - 2}`,
      stroke: 'currentColor', 'stroke-opacity': 0.4, 'stroke-width': 1.4,
      'stroke-dasharray': '4 3', fill: 'none', 'marker-end': 'url(#fl-arrow)' })));
    if (label) svg.appendChild(mk('text', { x: gutterX + 6, y: railY - 5, 'font-size': 10,
      fill: 'currentColor', 'fill-opacity': 0.65 }, label));
  };

  const fmtN = v => (v === null || v === undefined) ? null : Number(v).toLocaleString();
  const F = DATA.expr_filter || {};
  const refName = (DATA.annotation || {}).reference || 'all_sources_dedup';
  const rung = nm => (DATA.ladder || []).find(r => r.reference === nm) || null;
  const refRow = rung(refName), deRow = rung(refName + '_expressed');
  const busco = r => (r && r.busco) ? `BUSCO C:${r.busco.complete}% D:${r.busco.dup}%` : null;

  /* Assembly-only libraries (--assembly_samplesheet) get their own input box, dashed like every
     box that feeds evidence rather than counts. They join the two short-read assembly routes and
     nothing else, so the quantification half of the figure is unchanged. */
  const AR = DATA.assembly_reads || [];
  const nS = (DATA.samples || []).length;
  const illum = `Illumina · ${nS} sample${nS === 1 ? '' : 's'}` + (DATA.lanes > 1 ? ` × ${DATA.lanes} lanes` : '');
  let reads, extra = null, long;
  if (AR.length) {
    // PacBio | Illumina | assembly-only, so each input's edges run on their own level and the
    // assembly-only edges stay on the right, over the two routes they feed
    long  = box(20, 10, 236, 44, 'PacBio Iso-Seq', S1, { lines: ['pooled HQ · rRNA → decontam'] });
    reads = box(272, 10, 236, 44, `Illumina · ${nS} samples`, S1,
                { lines: [(DATA.lanes > 1 ? `× ${DATA.lanes} lanes · ` : '') + 'fastp → rRNA → decontam'] });
    extra = box(524, 10, 236, 44, `Public · ${AR.length} librar${AR.length === 1 ? 'y' : 'ies'}`, S1,
                { dashed: true, lines: ['assembly only · not quantified'] });
  } else {
    long  = box(40, 10, 300, 44, 'PacBio Iso-Seq · pooled HQ', S1, { lines: ['rRNA → decontaminate'] });
    reads = box(440, 10, 300, 44, illum, S1, { lines: ['merge → fastp → rRNA → decontaminate'] });
  }

  // four routes span the input row (20 .. 760); PacBio feeds A and B, short reads C and D
  const RT = { titleSize: 12 };
  const rA = box(20, 96, 176, 46, 'A  Iso-Seq on genome', S2, { ...RT, lines: ['pbmm2 + collapse'] });
  const rB = box(208, 96, 176, 46, 'B  reference-free', S2, { ...RT, lines: ['vsearch cluster'] });
  const rC = box(396, 96, 176, 46, 'C  STAR + StringTie', S2, { ...RT, lines: ['novel loci only'] });
  const rD = box(584, 96, 176, 46, 'D  Trinity de novo', S2, { ...RT, lines: ['genome-absent loci'] });

  const ladderB = box(218, 186, 340, 46, 'Reference ladder', S3,
                      { lines: ['isoseq → isoseq_stringtie → …_trinity_lrfree → ' + refName] });

  // The assembly, and the numbers quoted for it: annotation and completeness are BOTH measured
  // here, on the full deduplicated set, not on whatever the expression filter leaves behind.
  const annot = box(30, 262, 310, 62, 'The assembly  ·  ' + refName, MUTED,
                    { lines: ['TransDecoder → DIAMOND · eggNOG · InterProScan',
                              [refRow && refRow.transcripts ? fmtN(refRow.transcripts) + ' transcripts' : null,
                               busco(refRow)].filter(Boolean).join('  ·  ') || 'annotated before the filter'] });
  const q1    = box(420, 262, 310, 62, 'salmon pass 1  ·  ' + refName, MUTED,
                    { lines: ['measures what the experiment saw', 'TPM only — nothing here is reported'] });

  const filt  = box(218, 356, 512, 62, 'Expression filter  →  the DE reference', S3,
                    { lines: [F.rule || 'TPM ≥ 1 in ≥ 3 libraries',
                              [F.kept ? fmtN(F.kept) + ' kept, ' + fmtN(F.dropped) + ' dropped' : null,
                               busco(deRow)].filter(Boolean).join('  ·  ')
                              || 'cuts the reference for testing, not the assembly'] });
  const q2    = box(380, 450, 350, 46, 'salmon pass 2  ·  one quantification', S1,
                    { lines: ['decoy-aware, --dumpEq, 30 bootstraps — the reported counts'] });

  /* Three parallel branches, not three analyses. Same counts, same DESeq2 code, same thresholds;
     only the transcript→gene map differs, which is why a difference between them is a difference
     in gene identity and nothing else. Which one is primary is a parameter, so it is marked
     rather than assumed. */
  const MF = (DATA.map_flow || []).slice(0, 3);
  const cols = [[112, 210], [335, 210], [558, 210]];
  const deBoxes = [], enrBoxes = [];
  (MF.length ? MF : [{ label: 'gene map', basis: '', primary: false }]).forEach((m, k) => {
    const [x, w] = cols[k] || cols[2];
    const col = m.primary ? S3 : S2;
    deBoxes.push(box(x, 528, w, 62, m.label, col, {
      primary: !!m.primary,
      lines: ['DESeq2 · ' + (m.basis || 'gene map'),
              (m.sig === undefined || m.sig === null) ? 'no results in this run'
                : `${m.contrasts} contrasts · ${fmtN(m.sig)} sig (${fmtN(m.up)}↑/${fmtN(m.down)}↓)`]
    }));
  });
  (MF.length ? MF : [{}]).forEach((m, k) => {
    const [x, w] = cols[k] || cols[2];
    const t = m.terms || {};
    const keys = Object.keys(t);
    const sum = pat => keys.filter(kk => kk.endsWith(pat)).reduce((a, kk) => a + t[kk], 0);
    const tot = keys.reduce((a, kk) => a + t[kk], 0);
    enrBoxes.push(box(x, 626, w, 58, 'Enrichment', m.primary ? S3 : MUTED, {
      emphasis: !!m.primary, dashed: !m.primary,
      lines: [keys.length ? `${fmtN(tot)} terms at padj < ${m.enr_alpha || 0.05}` : 'GO / KEGG · ORA + GSEA',
              keys.length ? `ORA ${fmtN(sum('ORA'))}  ·  GSEA ${fmtN(sum('GSEA'))}` : 'not run']
    }));
  });

  // Long reads feed only the two long-read routes (A, B) and short reads only the two short-read
  // routes (C, D), so the input edges never cross. The assembly-only edges leave from the right
  // of their box, one level higher, and join C and D beside the Illumina arrows.
  elbow(long, rA, undefined, 72); elbow(long, rB, undefined, 72);
  elbow(reads, rC, undefined, 72); elbow(reads, rD, undefined, 72);
  if (extra) {
    const sx = extra.right - 40;
    [rC, rD].forEach(b => svg.appendChild(mk('path', {
      d: `M ${sx} ${extra.bottom} V 64 H ${b.cx + 12} V ${b.y - 2}`,
      stroke: 'currentColor', 'stroke-opacity': 0.45, 'stroke-width': 1.4, 'stroke-dasharray': '4 3',
      fill: 'none', 'marker-end': 'url(#fl-arrow)' })));
  }
  elbow(rA, ladderB); elbow(rB, ladderB); elbow(rC, ladderB); elbow(rD, ladderB);
  elbow(ladderB, annot, 'annotated in full'); elbow(ladderB, q1);
  elbow(q1, filt, 'TPM');
  // only drawn when the run actually used the ORF+homology rescue: with it off (the default) the
  // filter takes no input from the annotation at all, and an edge here would be a lie
  if (F.rescued) elbow(annot, filt, 'ORF + homology');
  FLOW_N = { annot: annot.n, filt: filt.n };
  arrow(filt, q2);
  deBoxes.forEach(b => elbow(q2, b, undefined, 512));
  deBoxes.forEach((b, k) => { if (enrBoxes[k]) arrow(b, enrBoxes[k]); });
  rail(annot, enrBoxes, 66, 338, 608, 'terms for enrichment');

  return svg;
}

/* Every threshold the run applied, in the order it applied them, with what each one removed.
   These are scattered across the parameter dump, the filter log and the DESeq2 filtering summary;
   a reader checking whether a number is defensible should not have to assemble them by hand. */
function thresholdTable() {
  const P = DATA.params || {};
  const F = DATA.expr_filter || {};
  const n = v => (v === null || v === undefined) ? '—' : Number(v).toLocaleString();
  const p = (k, d) => (P[k] !== undefined && P[k] !== null) ? String(P[k]) : String(d);
  const rows = [];
  const add = (stage, rule, level, effect) => rows.push({ stage, rule, level, effect });

  add('rRNA depletion', `bbduk k-mer match, k = ${p('rrna_kmer', 31)}`, 'read pair',
      'a pair is dropped when EITHER mate matches');
  add('Contaminant removal', 'best primary alignment is to a contaminant genome', 'read pair',
      'reads aligning nowhere are KEPT, not dropped');
  add('Long-read clustering (route B)', `vsearch --id ${p('cdhit_longread_id', 0.99)}`, 'transcript',
      'collapses near-identical long-read isoforms');
  add('Cross-source deduplication', `vsearch --id ${p('cdhit_dedup_id', 0.95)}`, 'transcript',
      F.genes_before ? `${n(F.genes_before)} genes enter the expression filter` : 'produces all_sources_dedup');
  add('Novel-vs-reference call', 'minimap2 -x asm20 -c --secondary=no', 'transcript',
      'exact match counts, not the approximate ones -x alone gives');

  // The ORF+homology rescue is off by default, and when it is off the "kept" and "DE-eligible"
  // sets are the same transcripts. Showing two rows that differ only in wording would invent a
  // distinction this run did not make, so the split is only spelled out when it is real.
  const rescued = Number(F.rescued || 0);
  add('Expression filter — the DE reference', F.rule ||
      `TPM >= ${p('expr_filter_tpm', 1)} in >= ${p('expr_filter_samples', 3)} libraries`,
      'transcript',
      F.kept ? `${n(F.kept)} kept, ${n(F.dropped)} dropped${rescued ? ` (${n(F.expressed)} expressed + ${n(rescued)} rescued)` : ''}` : '—');
  add('→ what it feeds', 'salmon pass 2, and every count the DE analysis uses', 'transcript',
      'NOT the annotation and NOT the reported BUSCO — both are measured before this cut');
  add('→ what it does NOT touch', 'the assembly, its annotation, its completeness', 'transcript',
      'dropped transcripts keep their annotation row; they are simply never quantified or tested');

  if (rescued) {
    add('DE-eligible set', `gene has >= 1 transcript at TPM >= ${p('expr_filter_tpm', 1)} in >= ${p('expr_filter_samples', 3)} libraries`,
        'gene',
        F.eligible ? `${n(F.eligible)} expressed transcripts define it` : 'rescued-only genes are not tested');
    add('→ what it feeds', 'DESeq2, and therefore the GO/KEGG ORA background', 'gene',
        'a tested gene keeps the counts of ALL its transcripts, rescued ones included');
  }

  add('DESeq2 pre-filter', `>= ${p('de_min_count', 5)} counts in >= ${p('de_min_samples', 3)} samples`, 'gene',
      'a presence rule, not a total: a sum passes a gene carried by one sample');
  add('Independent filtering', 'baseMean threshold chosen per contrast by DESeq2', 'gene',
      'see the per-contrast table below; it differs between contrasts');
  add("Cook's distance", 'default DESeq2 cut', 'gene',
      'with n = 5 per group DESeq2 cannot replace outliers, so these genes leave the test');

  add('Significance', `padj < ${p('de_alpha', 0.05)} AND |log2FC| >= ${p('de_lfc', 1)}`, 'gene',
      'the fold change is SHRUNKEN — this is the headline "sig" number');
  add('Fold-change shrinkage', `lfcShrink type = ${p('de_shrink', 'apeglm')}`, 'gene',
      'apeglm for every contrast; it shrinks sparse high-dispersion genes hard, by design');
  add('Enrichment', `p.adjust < ${p('de_alpha', 0.05)}, GSEA minGSSize 10, seed = ${p('random_seed', 42)}`, 'term',
      'tables are written in full; significance is applied once, here');
  add('Locus rebuild', `coverage >= ${p('locus_mincov', 0.9)}, group by ${p('locus_groupby', 'junction')}`, 'transcript',
      'transcripts below the coverage cut stay as their own gene (LOCUSU_*)');

  return buildTable({
    cols: [{ key: 'stage', label: 'Step' }, { key: 'rule', label: 'Rule' },
           { key: 'level', label: 'Applied to' }, { key: 'effect', label: 'Effect on this run' }],
    rows: rows
  });
}

function sectionMethods() {
  const out = h('div');
  out.appendChild(card({
    title: 'How the run fits together',
    why: 'what feeds what, and where each reported number is made',
    chart: flowDiagram(),
    note: ((DATA.assembly_reads || []).length
          ? 'The dashed input box is assembly-only: those libraries go through the same preprocessing ' +
            'and into StringTie and Trinity (dashed edges), but are never quantified or tested. ' : '') +
          'Boxes are numbered in execution order. Four things the diagram exists to make plain. ' +
          'Salmon runs TWICE — the first pass exists only to measure expression for the filter, ' +
          'and every number reported as a quantification comes from the second. The assembly and ' +
          'the DE reference are different sets and are measured separately: annotation and BUSCO ' +
          `completeness are both reported on the full deduplicated assembly (box ${FLOW_N.annot}), while only ` +
          `the expression-filtered subset (box ${FLOW_N.filt}) is indexed, quantified and tested — so the ` +
          'Annotation tab legitimately describes more transcripts than the DE tabs. The run then ' +
          'splits into three PARALLEL branches, one per gene map, each with its own DESeq2 and ' +
          'its own enrichment; they are three re-aggregations of ONE quantification, so their ' +
          'results are not independent evidence and agreement between them is not replication. ' +
          'The branch marked PRIMARY is the one dge_primary_map names and the one the headline ' +
          'numbers come from; the other two are published beside it as sensitivity arms. The ' +
          'dashed rail is the annotation supplying GO and KEGG terms to each enrichment — the ' +
          'ORA background is not the annotation but the genes that survived into each DESeq2 ' +
          'results table, which is why the three branches can return different term counts from ' +
          'the same term source.'
  }));
  out.appendChild(card({
    title: 'Filters and thresholds',
    why: 'every cut this run applied, in the order it applied them',
    tableNode: thresholdTable(),
    note: 'Read top to bottom: each row narrows what the row below it sees. The row most often ' +
          'misread is the expression filter — it decides what is QUANTIFIED AND TESTED, and ' +
          'nothing else. The transcript count, the annotation and the BUSCO completeness quoted ' +
          'for this transcriptome are all measured above it, on the full deduplicated assembly, ' +
          'so a transcript can be absent from the DE tables and still present in the annotation.'
  }));
  out.appendChild(h('p', { class: 'sub', text: 'Every stage as it actually ran: the arguments come from this run’s parameters, the numbers from its outputs. A skipped stage stays listed so the shape of the run is visible.' }));
  const groups = [];
  (DATA.stages || []).forEach(st => {
    let g = groups.find(x => x.name === st.group);
    if (!g) { g = { name: st.group, items: [] }; groups.push(g); }
    g.items.push(st);
  });
  groups.forEach(g => {
    out.appendChild(h('h2', { text: g.name }));
    g.items.forEach(st => {
      const body = h('div', { class: 'body' });
      body.appendChild(h('p', { class: 'why', text: st.why }));
      if (st.args.length) {
        const ul = h('ul', { class: 'args' });
        st.args.forEach(a => ul.appendChild(h('li', { text: a })));
        body.appendChild(ul);
      }
      if (st.metrics.length) {
        const m = h('div', { class: 'metrics' });
        st.metrics.forEach(x => m.appendChild(h('div', {}, [h('b', { text: x.value }), h('span', { text: x.label })])));
        body.appendChild(m);
      }
      if (st.output) body.appendChild(h('p', { class: 'note-line', text: 'output → ' + st.output }));
      const d = h('details', { class: 'stage', id: 'stage-' + st.id }, [
        h('summary', {}, [
          h('span', { class: 'nm', text: st.name }),
          h('span', { class: 'tool', text: st.tool }),
          h('span', { class: 'spacer' }),
          h('span', { class: 'badge ' + st.status, text: st.status })
        ]),
        body
      ]);
      out.appendChild(d);
    });
  });
  return out;
}

function sectionQC() {
  const out = h('div');
  const S = DATA.samples || [];
  if (!S.length) return h('p', { class: 'empty', text: 'No per-sample QC tables in this run.' });

  const rows = S.map(x => ({
    label: x.sample,
    parts: [
      { label: 'kept (host)', value: x.kept, color: S1 },
      { label: 'rRNA removed', value: x.rrna_removed, color: S2 },
      { label: 'contaminant removed', value: x.contam_removed, color: S3 }
    ],
    extraTip: [{ label: 'condition', value: x.condition }]
  }));
  out.appendChild(card({
    title: 'What survived per sample',
    why: 'read pairs, in the order they were removed',
    legend: [{ color: S1, label: 'kept (host)' }, { color: S2, label: 'rRNA removed' }, { color: S3, label: 'contaminant removed' }],
    chart: stackedBarH({
      rows: rows, gutter: 160,
      tickFmt: v => v >= 1e6 ? (v / 1e6).toFixed(0) + 'M' : fmt(v),
      endLabel: r => {
        const t = r.parts.reduce((a, p) => a + p.value, 0);
        return t ? (100 * r.parts[0].value / t).toFixed(0) + '% kept' : '';
      }
    }),
    table: {
      cols: [{ key: 'sample', label: 'Sample' }, { key: 'condition', label: 'Condition' },
      { key: 'reads_in', label: 'Pairs in', num: true }, { key: 'rrna_removed', label: 'rRNA removed', num: true },
      { key: 'contam_removed', label: 'Contaminant removed', num: true }, { key: 'kept', label: 'Kept', num: true },
      { key: 'contam_pct', label: '% contaminant', num: true }],
      rows: S.map(x => Object.assign({}, x, {
        reads_in: fmt(x.reads_in), rrna_removed: fmt(x.rrna_removed),
        contam_removed: fmt(x.contam_removed), kept: fmt(x.kept), contam_pct: pct(x.contam_pct)
      }))
    },
    note: 'Contamination is removed at the read level, before assembly: an assembly built from mixed reads cannot be cleanly separated afterwards. Reads that align to neither reference are kept — on a fragmented draft genome, unmapped is not evidence of contamination.'
  }));

  const mapped = S.filter(x => x.mapped_pct);
  if (mapped.length) {
    out.appendChild(card({
      title: 'Quantification mapping rate',
      why: 'salmon, selective alignment against a decoy-aware index',
      chart: barH({
        rows: mapped.map(x => ({ label: x.sample, value: x.mapped_pct, tip: [{ label: 'mapped', value: pct(x.mapped_pct) }, { label: 'fragments', value: fmt(x.fragments) }] })),
        max: 100, gutter: 160, color: S1, measure: 'mapped',
        valFmt: r => r.value.toFixed(1) + '%', tickFmt: v => v + '%'
      }),
      table: {
        cols: [{ key: 'sample', label: 'Sample' }, { key: 'library_type', label: 'Library' },
        { key: 'fragments', label: 'Fragments', num: true }, { key: 'mapped_pct', label: '% mapped', num: true }],
        rows: mapped.map(x => ({ sample: x.sample, library_type: x.library_type, fragments: fmt(x.fragments), mapped_pct: pct(x.mapped_pct) }))
      }
    }));
  }

  const lr = DATA.longread || {};
  const LS = lr.samples || [];
  if (LS.length || (lr.decontam || []).length || (lr.rrna || []).length) {
    out.appendChild(h('h2', { text: 'Long reads' }));
    out.appendChild(h('p', { class: 'sub', text: 'Long reads are pooled before they are mapped — one pass over the contaminant reference instead of one per sample — so these per-sample numbers are recovered from the membership map afterwards, not from separate runs.' }));
  }

  if (LS.length) {
    /* the same composition chart as the short reads, in the same colours */
    out.appendChild(card({
      title: 'Long-read isoforms per sample',
      why: 'in the order they were removed',
      legend: [{ color: S1, label: 'kept (host)' }, { color: S2, label: 'rRNA removed' }, { color: S3, label: 'contaminant removed' }],
      chart: stackedBarH({
        rows: LS.map(x => ({
          label: x.sample,
          parts: [
            { label: 'kept (host)', value: x.kept, color: S1 },
            { label: 'rRNA removed', value: x.rrna_removed, color: S2 },
            { label: 'contaminant removed', value: x.contaminant_removed, color: S3 }
          ],
          extraTip: [{ label: 'condition', value: x.condition }, { label: 'isoforms in', value: fmt(x.isoforms_in) }]
        })),
        gutter: 160,
        tickFmt: v => v >= 1e6 ? (v / 1e6).toFixed(1) + 'M' : (v >= 1000 ? (v / 1000).toFixed(0) + 'k' : fmt(v)),
        endLabel: r => {
          const t = r.parts.reduce((a, p) => a + p.value, 0);
          return t ? (100 * r.parts[0].value / t).toFixed(0) + '% kept' : '';
        }
      }),
      table: {
        cols: [{ key: 'sample', label: 'Sample' }, { key: 'condition', label: 'Condition' },
        { key: 'isoforms_in', label: 'Isoforms in', num: true }, { key: 'rrna_removed', label: 'rRNA removed', num: true },
        { key: 'contaminant_removed', label: 'Contaminant removed', num: true }, { key: 'kept', label: 'Kept', num: true },
        { key: 'pct_removed', label: '% contaminant', num: true }],
        rows: LS.map(x => ({
          sample: x.sample, condition: x.condition, isoforms_in: fmt(x.isoforms_in),
          rrna_removed: fmt(x.rrna_removed), contaminant_removed: fmt(x.contaminant_removed),
          kept: fmt(x.kept), pct_removed: pct(x.pct_removed)
        }))
      },
      note: lr.note || ''
    }));

    /* where each sample's isoforms placed - the long-read-specific question */
    out.appendChild(card({
      title: 'Where each sample’s isoforms placed',
      why: 'best alignment across the host and contaminant references',
      legend: [{ color: S1, label: 'host' }, { color: S2, label: 'no alignment' },
      { color: S3, label: 'contaminant' }, { color: S4, label: 'ambiguous' }],
      chart: stackedBarH({
        rows: LS.map(x => ({
          label: x.sample,
          parts: [
            { label: 'host', value: x.host, color: S1 },
            { label: 'no alignment', value: x.no_alignment, color: S2 },
            { label: 'contaminant', value: x.contaminant, color: S3 },
            { label: 'ambiguous', value: x.ambiguous, color: S4 }
          ],
          extraTip: [{ label: 'condition', value: x.condition }]
        })),
        gutter: 160,
        tickFmt: v => v >= 1e6 ? (v / 1e6).toFixed(1) + 'M' : (v >= 1000 ? (v / 1000).toFixed(0) + 'k' : fmt(v)),
        endLabel: r => {
          const t = r.parts.reduce((a, p) => a + p.value, 0);
          return t ? (100 * r.parts[0].value / t).toFixed(0) + '% host' : '';
        }
      }),
      table: {
        cols: [{ key: 'sample', label: 'Sample' }, { key: 'host', label: 'Host', num: true },
        { key: 'no_alignment', label: 'No alignment', num: true }, { key: 'contaminant', label: 'Contaminant', num: true },
        { key: 'ambiguous', label: 'Ambiguous', num: true }],
        rows: LS.map(x => ({
          sample: x.sample, host: fmt(x.host), no_alignment: fmt(x.no_alignment),
          contaminant: fmt(x.contaminant), ambiguous: fmt(x.ambiguous)
        }))
      },
      note: 'Contaminant and ambiguous are dropped; host and no-alignment are kept. Only sequences seen in at least one alignment are classified, so "no alignment" is the complement — on a fragmented draft genome, unmapped is not evidence of contamination.'
    }));
  }

  if ((lr.decontam || []).length) {
    const c = h('div', { class: 'card' }, [h('div', { class: 'figcap' }, [
      h('h3', { text: 'Pooled long-read classification' }),
      h('span', { class: 'why', text: 'the set the reference was actually built from' })])]);
    c.appendChild(buildTable({
      cols: [{ key: 'class', label: 'Class' }, { key: 'n', label: 'Sequences', num: true }],
      rows: lr.decontam.map(r => ({ class: r.class, n: fmt(r.n) }))
    }));
    if (LS.length) c.appendChild(h('p', { class: 'note-line', text: 'A sequence detected in several samples is one sequence here and one row per sample above, so the per-sample columns sum to more than these totals.' }));
    out.appendChild(c);
  }
  return out;
}

function sectionReference() {
  const L = DATA.ladder || [];
  if (!L.length) return h('p', { class: 'empty', text: 'No reference ladder in this run.' });
  const out = h('div');
  out.appendChild(h('p', { class: 'sub', text: 'Each assembly source is layered onto the previous one and only its novel contribution is added, so the table says what every source was worth instead of hiding it in one merge.' }));

  let measure = 'transcripts';
  const holder = h('div');
  function draw() {
    holder.textContent = '';
    const rows = L.filter(r => measure === 'transcripts' || r.genes).map(r => ({
      label: r.reference, value: measure === 'transcripts' ? r.transcripts : r.genes,
      tip: [{ label: 'transcripts', value: fmt(r.transcripts) }, { label: 'genes', value: r.genes ? fmt(r.genes) : '—' },
      { label: 'N50', value: r.n50 ? fmt(r.n50) : '—' }, { label: 'basis', value: r.basis }]
    }));
    holder.appendChild(card({
      title: measure === 'transcripts' ? 'Transcripts per reference variant' : 'Genes per reference variant',
      why: 'the assembly is ' + (DATA.meta.reference_set || '') + '; the expression-filtered rung below it is what DE runs against',
      chart: barH({ rows: rows, gutter: 150, color: S1, measure: measure }),
      table: {
        cols: [{ key: 'reference', label: 'Reference' }, { key: 'basis', label: 'Basis' },
        { key: 'transcripts', label: 'Transcripts', num: true }, { key: 'genes', label: 'Genes', num: true },
        { key: 'n50', label: 'N50', num: true }, { key: 'busco', label: 'BUSCO complete', num: true }],
        rows: L.map(r => ({
          reference: r.reference, basis: r.basis, transcripts: fmt(r.transcripts),
          genes: r.genes ? fmt(r.genes) : '—', n50: r.n50 ? fmt(r.n50) : '—',
          busco: r.busco ? r.busco.complete.toFixed(1) + '%' : '—'
        }))
      },
      note: measure === 'genes' ? 'These gene counts come from string-stripped transcript ids and are an upper bound on distinct loci — see the Genes tab.' : ''
    }));
  }
  const controls = h('div', { class: 'controls' }, [
    h('span', { class: 'why', text: 'Measure' }),
    h('button', { class: 'ghost', 'aria-pressed': 'true', text: 'Transcripts', onclick: e => { measure = 'transcripts'; setPair(e.target); draw(); } }),
    h('button', { class: 'ghost', 'aria-pressed': 'false', text: 'Genes', onclick: e => { measure = 'genes'; setPair(e.target); draw(); } })
  ]);
  function setPair(active) {
    [...controls.querySelectorAll('button')].forEach(b => b.setAttribute('aria-pressed', String(b === active)));
  }
  out.appendChild(controls);
  out.appendChild(holder);
  draw();

  const withB = L.filter(r => r.busco);
  if (withB.length) {
    out.appendChild(card({
      title: 'BUSCO completeness',
      why: 'complete = single + duplicated; the track is the full lineage set',
      legend: [{ color: 'var(--ramp-1)', label: 'complete, single copy' }, { color: 'var(--ramp-2)', label: 'complete, duplicated' }],
      chart: meterRows({
        rows: withB.map(r => ({
          label: r.reference,
          parts: [{ value: r.busco.single, color: 'var(--ramp-1)' }, { value: r.busco.dup, color: 'var(--ramp-2)' }],
          endLabel: r.busco.complete.toFixed(1) + '%',
          tip: [{ color: 'var(--ramp-1)', label: 'single', value: r.busco.single.toFixed(1) + '%' },
          { color: 'var(--ramp-2)', label: 'duplicated', value: r.busco.dup.toFixed(1) + '%' },
          { label: 'fragmented', value: r.busco.frag.toFixed(1) + '%' },
          { label: 'missing', value: r.busco.missing.toFixed(1) + '%' }]
        })), gutter: 150
      }),
      table: {
        cols: [{ key: 'reference', label: 'Reference' }, { key: 'c', label: 'Complete', num: true }, { key: 's', label: 'Single', num: true },
        { key: 'd', label: 'Duplicated', num: true }, { key: 'f', label: 'Fragmented', num: true }, { key: 'm', label: 'Missing', num: true }],
        rows: withB.map(r => ({
          reference: r.reference, c: r.busco.complete.toFixed(1) + '%', s: r.busco.single.toFixed(1) + '%',
          d: r.busco.dup.toFixed(1) + '%', f: r.busco.frag.toFixed(1) + '%', m: r.busco.missing.toFixed(1) + '%'
        }))
      },
      note: 'A rising duplicated fraction as sources are layered on is the same phenomenon as the inflated gene count: one locus recovered by several assemblers.'
    }));
  }
  return out;
}

function sectionAnnotation() {
  const A = DATA.annotation || {};
  if (!(A.evidence || []).length && !(A.summary || []).length)
    return h('p', { class: 'empty', text: 'Annotation was skipped in this run.' });
  const out = h('div');
  // where these numbers come from: readers otherwise assume they describe the DE reference
  if (A.note) out.appendChild(h('p', { class: 'caveat', text:
    'Source: the ' + (A.reference || 'selected') + ' transcript set. ' + A.note }));
  if ((A.evidence || []).length) {
    const total = A.evidence[0].total;
    out.appendChild(card({
      title: 'Annotation evidence',
      why: `transcripts carrying each line of evidence, of ${fmt(total)}`,
      chart: barH({
        rows: A.evidence.map(e => ({
          label: e.label, value: e.n,
          tip: [{ label: 'transcripts', value: fmt(e.n) }, { label: 'of total', value: (100 * e.n / total).toFixed(1) + '%' }]
        })), gutter: 175, color: S1, measure: 'transcripts'
      }),
      table: {
        cols: [{ key: 'label', label: 'Evidence' }, { key: 'n', label: 'Transcripts', num: true }, { key: 'p', label: '% of set', num: true }],
        rows: A.evidence.map(e => ({ label: e.label, n: fmt(e.n), p: (100 * e.n / total).toFixed(1) + '%' }))
      },
      note: 'has_ORF is yes/no on every row of the annotation table, so counting non-empty values there returns the whole table — the ORF count is coding_status == "coding".'
    }));
  }
  if ((A.summary || []).length) {
    out.appendChild(h('div', { class: 'card' }, [
      h('div', { class: 'figcap' }, [h('h3', { text: 'Annotation summary' })]),
      buildTable({
        cols: [{ key: 'label', label: 'Metric' }, { key: 'value', label: 'Count', num: true }],
        rows: A.summary.map(x => ({ label: x.label, value: fmt(x.value) }))
      })
    ]));
  }
  return out;
}
"""

SCRIPT3 = r"""

/* An embedded diagnostic figure. The PNG is a data URI, so this still works offline. */
function figureCard(name, title, why, note) {
  const src = ((DATA.figures_by_map || {})[activeMapKey()] || {})[name];
  if (!src) return null;
  const img = h('img', { src: src, alt: title, loading: 'lazy',
                         style: 'width:100%;height:auto;display:block;border-radius:4px' });
  return card({ title: title, why: why, chart: img, note: note || '' });
}


/* Which gene map the DE plots are showing. The three maps are separate answers to "what is a
   gene", so every plot drawn from one of them says which, and can be switched to the others
   instead of leaving the reader to guess. */
const MAP_LABELS = { string: 'String-stripped ids', corset: 'Corset', locus: 'Locus rebuild' };
function geneMapSets() {
  const out = [];
  if (DATA.dge) out.push({ key: 'string', data: DATA.dge });
  if (DATA.dge_corset) out.push({ key: 'corset', data: DATA.dge_corset });
  if (DATA.dge_locus) out.push({ key: 'locus', data: DATA.dge_locus });
  return out.map(m => ({
    key: m.key, data: m.data,
    label: MAP_LABELS[m.key] + (m.key === (DATA.meta || {}).primary_map ? ' (primary)' : '')
  }));
}
/* Opens on the primary map, so the first thing a reader sees is the result the headline quotes. */
let ACTIVE_MAP = (DATA.meta || {}).primary_map || null;
function activeMapKey() {
  const sets = geneMapSets();
  const hit = sets.find(m => m.key === ACTIVE_MAP);
  return (hit || sets[0] || { key: null }).key;
}
function activeDge() {
  const sets = geneMapSets();
  if (!sets.length) return null;
  const hit = sets.find(m => m.key === ACTIVE_MAP);
  return (hit || sets[0]).data;
}
function activeMapLabel() {
  const sets = geneMapSets();
  const hit = sets.find(m => m.key === ACTIVE_MAP);
  return (hit || sets[0] || { label: '—' }).label;
}
/* The switcher. Changing the map rebuilds the tab it sits in, so every plot below follows. */
function mapSwitcher(tabId) {
  const sets = geneMapSets();
  if (sets.length < 2) return null;
  const sel = h('select', {
    'aria-label': 'gene map used for the plots below',
    onchange: e => {
      ACTIVE_MAP = e.target.value;
      const panel = document.getElementById('panel-' + tabId);
      panel.textContent = '';
      built[tabId] = false;
      showTab(tabId);
    }
  }, sets.map(m => h('option', { value: m.key, text: m.label,
                                 selected: m.key === (ACTIVE_MAP || sets[0].key) ? '' : null })));
  return h('div', { class: 'card' }, [
    h('div', { class: 'figcap' }, [
      h('h3', { text: 'Gene map for the plots below' }),
      h('span', { class: 'why', text: 'the same counts, re-aggregated under a different definition of a gene' })
    ]),
    h('div', { class: 'filters' }, [h('label', { text: 'Map' }), sel]),
    h('p', { class: 'note-line', text:
      'Every DE figure in this tab is drawn from the map selected here. The three disagree — see ' +
      '"Where the maps agree" in the Genes tab — so a gene list is only meaningful alongside the map it came from.' })
  ]);
}

function sectionExpression() {
  const D = activeDge();
  if (!D || !D.contrasts.length) return h('p', { class: 'empty', text: 'Differential expression was skipped in this run.' });
  const out = h('div');
  const sw = mapSwitcher('expression');
  if (sw) out.appendChild(sw);

  /* summary first - it is about every contrast, so it sits above the filter row */
  const rows = [];
  D.contrasts.forEach(c => {
    rows.push({ label: c.id + ' · up', value: c.up, color: UP, tip: [{ color: UP, label: 'up (padj<0.05, log2FC>1)', value: fmt(c.up) }, { label: 'genes tested', value: fmt(c.tested) }] });
    rows.push({ label: c.id + ' · down', value: c.down, color: DOWN, tip: [{ color: DOWN, label: 'down (padj<0.05, log2FC<-1)', value: fmt(c.down) }, { label: 'genes tested', value: fmt(c.tested) }] });
  });
  out.appendChild(card({
    title: 'Differentially expressed genes — ' + activeMapLabel(),
    why: 'padj < 0.05 and |log2FC| > 1, counted on the gene map named above',
    legend: [{ color: UP, label: 'up' }, { color: DOWN, label: 'down' }],
    chart: barH({ rows: rows, gutter: 210, measure: 'genes' }),
    table: {
      cols: [{ key: 'id', label: 'Contrast' }, { key: 'tested', label: 'Genes tested', num: true },
      { key: 'sigp', label: 'padj only', num: true }, { key: 'sig', label: 'padj + |log2FC|', num: true },
      { key: 'up', label: 'Up', num: true }, { key: 'down', label: 'Down', num: true }],
      rows: D.contrasts.map(c => ({ id: c.id, tested: fmt(c.tested), sigp: fmt(c.sig_padj), sig: fmt(c.sig), up: fmt(c.up), down: fmt(c.down) }))
    },
    note: 'The padj < 0.05 column has no fold-change filter, so it exceeds up + down. Genes DESeq2 filtered out independently carry padj = NA and are not in "tested".'
  }));

  if ((D.pca || []).length) {
    const groups = [...new Set(D.pca.map(p => p.group))];
    const palette = [S1, S2, S3];
    const colorOf = g => groups.indexOf(g) < 3 ? palette[groups.indexOf(g)] : MUTED;
    const pts = D.pca.map(p => [p.x, p.y, p.group, p.name]);
    out.appendChild(card({
      title: 'Sample similarity (PCA of variance-stabilised counts)',
      why: 'hover a point for the sample',
      legend: groups.slice(0, 3).map(g => ({ color: colorOf(g), label: g }))
        .concat(groups.length > 3 ? [{ color: MUTED, label: 'other conditions (see table)' }] : []),
      chart: scatter({
        points: pts, height: 300, xLabel: 'PC1', yLabel: 'PC2',
        classOf: p => ({ color: colorOf(p[2]), r: 5, ring: true }),
        tipTitle: p => p[3] || p[2],
        tipRows: p => [{ color: colorOf(p[2]), label: 'condition', value: p[2] }, { label: 'PC1', value: p[0].toFixed(2) }, { label: 'PC2', value: p[1].toFixed(2) }]
      }),
      table: {
        cols: [{ key: 'name', label: 'Sample' }, { key: 'group', label: 'Condition' }, { key: 'x', label: 'PC1', num: true }, { key: 'y', label: 'PC2', num: true }],
        rows: D.pca.map(p => ({ name: p.name, group: p.group, x: p.x.toFixed(2), y: p.y.toFixed(2) }))
      },
      note: 'If the conditions do not separate here, the design is not the dominant axis of variation: the DE statistics still hold, but single-gene claims need independent confirmation.'
    }));
  }

  /* one filter row, scoping everything below it */
  let contrast = D.contrasts[0].id;
  const sel = h('select', { 'aria-label': 'Contrast', onchange: e => { contrast = e.target.value; drawScoped(); } });
  D.contrasts.forEach(c => sel.appendChild(h('option', { value: c.id, text: c.id })));
  out.appendChild(h('div', { class: 'controls' }, [h('span', { class: 'why', text: 'Contrast' }), sel]));
  const scoped = h('div');
  out.appendChild(scoped);

  function drawScoped() {
    scoped.textContent = '';
    const c = D.contrasts.find(x => x.id === contrast);
    const maxY = Math.max(4, ...c.points.map(p => p[1]));
    scoped.appendChild(card({
      title: 'Volcano — ' + c.id + '  ·  ' + activeMapLabel(),
      why: `${fmt(c.n_points_total)} genes tested; every gene with padj < ${c.ceiling} is plotted, the rest sampled for context`,
      legend: [{ color: UP, label: 'up' }, { color: DOWN, label: 'down' }, { color: MUTED, label: 'not significant' }],
      chart: scatter({
        points: c.points, height: 340,
        xLabel: 'log2 fold change (apeglm-shrunken)', yLabel: '−log10 padj',
        yMin: 0, yMax: maxY,
        vlines: [{ x: -1 }, { x: 1 }], hlines: [{ y: -Math.log10(0.05), label: 'padj 0.05' }],
        classOf: p => isSig(p, 0.05, 1) ? { color: p[0] > 0 ? UP : DOWN, r: 2.8 } : { color: MUTED, r: 2, opacity: .5 },
        tipTitle: p => (isSig(p, 0.05, 1) ? 'significant' : 'not significant'),
        tipRows: p => [{ label: 'log2FC', value: p[0].toFixed(2) }, { label: 'padj', value: Math.pow(10, -p[1]).toExponential(2) }]
      }),
      table: {
        cols: [{ key: 'gene', label: 'Gene' }, { key: 'base', label: 'baseMean', num: true },
        { key: 'lfc', label: 'log2FC', num: true }, { key: 'padj', label: 'padj', num: true }, { key: 'desc', label: 'Annotation' }],
        rows: c.top.slice(0, 100).map(g => ({ gene: g.gene, base: fmt(g.base), lfc: g.lfc.toFixed(2), padj: String(g.padj), desc: g.desc }))
      },
      note: 'Fold changes are shrunken (apeglm), so the x-axis reflects evidence rather than ratio size; the raw estimate is kept as log2FoldChange_unshrunken in the CSV. The table lists the top 100 by padj; the Genes tab has a volcano with adjustable thresholds.'
    }));

    /* Diagnostics: whether the numbers above can be believed at all. */
    (DATA.figure_meta || []).forEach(f => {
      const c = figureCard(f.name, f.title + ' — ' + activeMapLabel(), f.why);
      if (c) out.appendChild(c);
    });
    (DATA.ma_figures || []).forEach(n => {
      const nm = n.replace('F13_MA_', '').replace(/_/g, ' ');
      const c = figureCard(n, 'MA — ' + nm + ' · ' + activeMapLabel(),
        'where the significant genes sit in expression',
        'Significance is applied to the shrunken fold change, so low-count genes no longer fly to the edges.');
      if (c) out.appendChild(c);
    });

    const FILT = (DATA.de_filtering_by_map || {})[activeMapKey()] || [];
    if (FILT.length) {
      out.appendChild(card({
        title: 'What the filters removed — ' + activeMapLabel(),
        why: 'per contrast, before the test was read',
        table: {
          cols: [{ key: 'contrast', label: 'Contrast' }, { key: 'in_model', label: 'Genes in model', num: true },
                 { key: 'tested', label: 'Tested', num: true }, { key: 'threshold', label: 'Filter baseMean', num: true },
                 { key: 'filtered', label: 'Independent filter', num: true },
                 { key: 'cooks', label: "Cook's outliers", num: true }, { key: 'shrinkage', label: 'Shrinkage' }],
          rows: FILT.map(r => ({ contrast: r.contrast, in_model: fmt(r.in_model),
            tested: fmt(r.tested), threshold: r.threshold, filtered: fmt(r.filtered),
            cooks: fmt(r.cooks), shrinkage: r.shrinkage }))
        },
        note: "Independent filtering trims low-count genes at a per-contrast baseMean threshold. Cook's outliers leave the test entirely — with 5 replicates DESeq2 cannot replace them, so they carry padj = NA."
      }));
    }

    if ((DATA.gene_maps || []).length) {
      out.appendChild(card({
        title: 'Three gene maps, same reference',
        why: 'none is unambiguously right — compare, do not average',
        chart: barH({ rows: (DATA.gene_maps || []).map((m, idx) => ({
            label: m.name, value: m.genes, color: [S1, S2, S3][idx % 3],
            tip: [{ label: m.basis, value: fmt(m.genes) + ' genes' }] })), gutter: 190, measure: 'genes' }),
        table: {
          cols: [{ key: 'name', label: 'Map' }, { key: 'basis', label: 'Built from' },
                 { key: 'genes', label: 'Genes', num: true },
                 { key: 'sig', label: 'DE genes (all contrasts)', num: true },
                 { key: 'note', label: 'Caveat' }],
          rows: (DATA.gene_maps || []).map(m => ({ name: m.name + (m.primary ? '  (primary)' : ''), basis: m.basis, genes: fmt(m.genes),
                 sig: (m.sig === null || m.sig === undefined) ? '—' : fmt(m.sig), note: m.note }))
        },
        note: 'Every map runs through the same quantifications and the same DESeq2 code, so a difference in DE genes is a difference in gene identity and nothing else. The headline numbers come from the map named by dge_primary_map; the others are published beside them as sensitivity analyses.'
      }));
    }

    const E = ((DATA.enrichment_by_map || {})[activeMapKey()] || {})[contrast];
    if (E && Object.keys(E).length) {
      let analysis = Object.keys(E)[0];
      const holder = h('div');
      const asel = h('select', { 'aria-label': 'Analysis', onchange: e => { analysis = e.target.value; drawEnr(); } });
      Object.keys(E).forEach(k => asel.appendChild(h('option', { value: k, text: k })));
      function drawEnr() {
        holder.textContent = '';
        const terms = (E[analysis] || []).filter(t => t.padj < 0.25).slice(0, 12);
        if (!terms.length) { holder.appendChild(h('p', { class: 'empty', text: 'No terms returned for this analysis.' })); return; }
        holder.appendChild(card({
          title: analysis + ' — ' + contrast + ' · ' + activeMapLabel(),
          why: 'terms ranked by adjusted p-value',
          chart: barH({
            rows: terms.map(t => ({
              label: t.desc || t.id, value: -Math.log10(Math.max(t.padj, 1e-300)),
              color: t.padj < 0.05 ? S1 : MUTED,
              tip: [{ label: 'padj', value: t.padj.toExponential(2) }, { label: 'set size', value: fmt(t.size) },
              { label: 'NES', value: t.nes ? t.nes.toFixed(2) : '—' }, { label: 'term', value: t.id },
              { label: 'description', value: t.desc || '' }]
            })),
            gutter: 300, right: 60, measure: '−log10 padj', valFmt: r => r.value.toFixed(1)
          }),
          table: {
            cols: [{ key: 'id', label: 'Term' }, { key: 'desc', label: 'Description' }, { key: 'padj', label: 'padj', num: true },
            { key: 'nes', label: 'NES', num: true }, { key: 'size', label: 'Size', num: true }],
            rows: (E[analysis] || []).map(t => ({ id: t.id, desc: t.desc, padj: t.padj.toExponential(2), nes: t.nes ? t.nes.toFixed(2) : '—', size: fmt(t.size) }))
          },
          note: 'Bars below the padj 0.05 line are shown in grey — they are returned by the test, not significant. GSEA is seeded, so this list is reproducible; unseeded it is not.'
        }));
      }
      scoped.appendChild(h('div', { class: 'controls' }, [h('span', { class: 'why', text: 'Analysis' }), asel]));
      scoped.appendChild(holder);
      drawEnr();
    }
  }
  drawScoped();
  return out;
}


/* A three-set Venn of the transcripts each gene map calls differentially expressed.
   Drawn to scale in position but not in area: the counts are the point, and area-proportional
   three-circle Venns cannot represent most count combinations honestly anyway. */
function vennSvg(v) {
  const W = 420, H = 300, r = 88;
  const cx = [W / 2 - 52, W / 2 + 52, W / 2], cy = [H / 2 - 34, H / 2 - 34, H / 2 + 52];
  const cols = [S1, S2, S3];
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('role', 'img');
  const labs = v.labels || [];
  svg.setAttribute('aria-label',
    'Venn of differentially expressed transcripts shared between the gene maps: ' +
    labs.map((l, i) => `${l} ${fmt((v.sizes || [])[i] || 0)}`).join(', '));

  const mk = (tag, attrs) => {
    const e = document.createElementNS('http://www.w3.org/2000/svg', tag);
    Object.entries(attrs).forEach(([k, val]) => e.setAttribute(k, val));
    return e;
  };
  const n = Math.min(labs.length, 3);
  for (let i = 0; i < n; i++) {
    svg.appendChild(mk('circle', { cx: cx[i], cy: cy[i], r: r, fill: cols[i],
                                   'fill-opacity': 0.22, stroke: cols[i], 'stroke-width': 2 }));
  }
  const R = v.regions || {};
  const put = (x, y, val, strong) => {
    if (val === undefined) return;
    const t = mk('text', { x: x, y: y, 'text-anchor': 'middle', 'font-size': strong ? 15 : 13,
                           'font-weight': strong ? 600 : 400, fill: 'currentColor' });
    t.textContent = fmt(val);
    svg.appendChild(t);
  };
  if (n === 3) {
    put(cx[0] - 46, cy[0] - 6, R.A); put(cx[1] + 46, cy[1] - 6, R.B); put(cx[2], cy[2] + 58, R.C);
    put(W / 2, cy[0] - 18, R.AB); put(cx[0] - 20, cy[2] + 8, R.AC); put(cx[1] + 20, cy[2] + 8, R.BC);
    put(W / 2, cy[0] + 44, R.ABC, true);
  } else if (n === 2) {
    put(cx[0] - 44, cy[0], R.A); put(cx[1] + 44, cy[1], R.B); put(W / 2, cy[0], R.AB, true);
  }
  for (let i = 0; i < n; i++) {
    const t = mk('text', { x: cx[i], y: i === 2 ? H - 6 : 20, 'text-anchor': 'middle',
                           'font-size': 12, fill: cols[i], 'font-weight': 600 });
    t.textContent = labs[i];
    svg.appendChild(t);
  }
  return svg;
}

function sectionGenes() {
  const out = h('div');

  /* What each gene map makes of the same experiment, and where they agree. */
  const V = DATA.venn || {};
  const vContrasts = Object.keys(V);
  if (vContrasts.length) {
    const maps = (DATA.gene_maps || []);
    if (maps.length) {
      out.appendChild(card({
        title: 'The same experiment under three definitions of a gene',
        why: 'identical quantifications, identical DESeq2 code — only the transcript-to-gene map differs',
        table: {
          cols: [{ key: 'name', label: 'Map' }, { key: 'basis', label: 'Built from' },
                 { key: 'genes', label: 'Genes', num: true },
                 { key: 'sig', label: 'DE genes (all contrasts)', num: true }],
          rows: maps.map(m => ({ name: m.name + (m.primary ? '  (primary)' : ''), basis: m.basis,
                                 genes: fmt(m.genes),
                                 sig: (m.sig === null || m.sig === undefined) ? '—' : fmt(m.sig) }))
        },
        note: 'A map that calls more genes differentially expressed is not the better map. Read this against the over- and under-merging measured below.'
      }));
    }

    let vc = vContrasts[0];
    const vennBox = h('div');
    const draw = () => {
      vennBox.textContent = '';
      const v = V[vc];
      vennBox.appendChild(vennSvg(v));
      vennBox.appendChild(h('p', { class: 'note-line', text:
        `${fmt(v.union)} transcripts are implicated by at least one map in ${vc}; ` +
        `${fmt((v.regions || {}).ABC || 0)} by all three. Counted on TRANSCRIPTS, because the maps ` +
        `do not share gene ids — a transcript is counted when it belongs to a gene that map calls ` +
        `significant. A map that merges more transcripts into each gene therefore expands to more ` +
        `transcripts, so the circle sizes are not normalised for that.` }));
    };
    const picker = h('div', { class: 'filters' }, [
      h('label', { text: 'Contrast' }),
      h('select', {
        onchange: e => { vc = e.target.value; draw(); }
      }, vContrasts.map(c => h('option', { value: c, text: c })))
    ]);
    const vcard = h('div', { class: 'card' }, [
      h('div', { class: 'figcap' }, [
        h('h3', { text: 'Where the maps agree' }),
        h('span', { class: 'why', text: 'shared differentially expressed transcripts' })
      ]),
      picker, vennBox
    ]);
    draw();
    out.appendChild(vcard);
  }

  const L = DATA.locus;
  if (L && L.before && L.after) {
    const c = h('div', { class: 'card' }, [
      h('div', { class: 'figcap' }, [h('h3', { text: 'How many genes? A range, not one number' })]),
      h('p', { class: 'sub', text: 'String-stripped transcript ids apply a different rule per assembly source and reconcile nothing across sources, so a locus recovered by three sources counts three times. Grouping by genomic locus fixes that and introduces the opposite error: tandem duplicates and nested genes merge on a fragmented draft. Corset sits between them, clustering on shared reads rather than on sequence or position, and is run label-free so the gene definition never sees the conditions being tested. The two error directions are measured below; the third map is compared in the panels above.' })
    ]);
    c.appendChild(barH({
      rows: [
        { label: 'string-stripped ids', value: L.before, color: S2, tip: [{ label: 'genes', value: fmt(L.before) }, { label: 'error mode', value: 'under-merged (upper bound)' }] },
        { label: 'genomic loci', value: L.after, color: S1, tip: [{ label: 'genes', value: fmt(L.after) }, { label: 'error mode', value: 'over-merged (lower bound)' }] }
      ], gutter: 175, measure: 'genes'
    }));
    const met = [];
    if (L.ratio_before !== null && L.ratio_before !== undefined)
      met.push({ k: 'Genes per protein accession', a: L.ratio_before, b: L.ratio_after, good: 1.0, hint: 'above 1.00 = one protein split across several genes (under-merged)' });
    if (L.overmerge_before !== null && L.overmerge_before !== undefined)
      met.push({ k: 'Loci carrying >1 distinct accession', a: L.overmerge_before, b: L.overmerge_after, suffix: '%', hint: 'rising = distinct genes fused into one locus (over-merged)' });
    if (met.length) {
      c.appendChild(buildTable({
        cols: [{ key: 'k', label: 'Validation metric' }, { key: 'a', label: 'String-stripped', num: true }, { key: 'b', label: 'Genomic loci', num: true }, { key: 'hint', label: 'Reads as' }],
        rows: met.map(m => ({ k: m.k, a: m.a.toFixed(2) + (m.suffix || ''), b: m.b.toFixed(2) + (m.suffix || ''), hint: m.hint }))
      }));
    }
    if (L.warning) c.appendChild(h('p', { class: 'caveat warnbox', text: 'The rebuild reported a sharp rise in over-merging. Raise the locus coverage threshold, or keep the string-stripped map for anything that depends on gene identity.' }));
    c.appendChild(h('p', { class: 'note-line', text: `${fmt(L.placed)} transcripts placed on the genome, ${fmt(L.unplaced)} unplaced (each kept as its own gene).` }));
    const det = h('details', {}, [h('summary', { class: 'ghost', text: 'Full validation report' })]);
    det.appendChild(h('pre', { class: 'mono scroll', text: L.report }));
    c.appendChild(det);
    out.appendChild(c);
  }

  const sw2 = mapSwitcher('genes');
  if (sw2) out.appendChild(sw2);
  const D = activeDge();
  if (!D) {
    if (!L) out.appendChild(h('p', { class: 'empty', text: 'No gene-level results in this run: differential expression and the locus rebuild were both skipped.' }));
    return out;
  }
  const all = [];
  D.contrasts.forEach(c => c.top.forEach(g => all.push(Object.assign({ contrast: c.id }, g))));
  if (!all.length) return out;

  const ceiling = D.contrasts[0].ceiling || 0.25;
  let f = { contrast: 'all', dir: 'all', q: '', padj: 0.05, lfc: 1 };
  let sortKey = 'padj', sortDir = 1;
  const body = h('div', { class: 'scroll' });
  const volcanoes = h('div');
  const count = h('span', { class: 'why' });

  const csel = h('select', { 'aria-label': 'Contrast', onchange: e => { f.contrast = e.target.value; draw(); } });
  csel.appendChild(h('option', { value: 'all', text: 'all contrasts' }));
  D.contrasts.forEach(c => csel.appendChild(h('option', { value: c.id, text: c.id })));
  const dsel = h('select', { 'aria-label': 'Direction', onchange: e => { f.dir = e.target.value; draw(); } });
  [['all', 'up and down'], ['up', 'up only'], ['down', 'down only'], ['flag', DATA.meta.interest_label || 'family of interest']]
    .forEach(([v, t]) => dsel.appendChild(h('option', { value: v, text: t })));
  const padjIn = h('input', {
    type: 'number', min: '0.0001', max: String(ceiling), step: '0.01', value: '0.05',
    'aria-label': 'adjusted p-value threshold',
    title: `only genes below padj ${ceiling} are embedded, so the threshold is capped there`,
    oninput: e => {
      const v = parseFloat(e.target.value);
      if (!(v > 0)) return;
      f.padj = Math.min(v, ceiling);
      if (v > ceiling) e.target.value = String(ceiling);   // show the value that is actually applied
      draw();
    }
  });
  const lfcIn = h('input', {
    type: 'number', min: '0', step: '0.5', value: '1', 'aria-label': 'minimum absolute log2 fold change',
    oninput: e => { const v = parseFloat(e.target.value); if (v >= 0) { f.lfc = v; draw(); } }
  });
  const q = h('input', { type: 'search', placeholder: 'search gene id or annotation', oninput: e => { f.q = e.target.value.toLowerCase(); draw(); } });

  function passes(g) {
    return (f.contrast === 'all' || g.contrast === f.contrast)
      && (f.dir === 'all' || (f.dir === 'up' && g.lfc > 0) || (f.dir === 'down' && g.lfc < 0) || (f.dir === 'flag' && g.flag))
      && g.padj < f.padj && Math.abs(g.lfc) >= f.lfc
      && (!f.q || g.gene.toLowerCase().includes(f.q) || (g.desc || '').toLowerCase().includes(f.q));
  }

  function drawTable() {
    let rows = all.filter(passes);
    rows.sort((a, b) => (a[sortKey] > b[sortKey] ? 1 : a[sortKey] < b[sortKey] ? -1 : 0) * sortDir);
    count.textContent = `${rows.length} of ${all.length} rows`;
    body.textContent = '';
    const cols = [{ key: 'gene', label: 'Gene' }, { key: 'contrast', label: 'Contrast' },
    { key: 'base', label: 'baseMean', num: true }, { key: 'lfc', label: 'log2FC', num: true },
    { key: 'padj', label: 'padj', num: true }, { key: 'desc', label: 'Annotation' }];
    const t = h('table');
    const tr = h('tr');
    cols.forEach(c => {
      const th = h('th', { class: c.num ? 'num' : '', text: c.label + (sortKey === c.key ? (sortDir > 0 ? ' ↑' : ' ↓') : '') });
      th.style.cursor = 'pointer';
      th.addEventListener('click', () => { if (sortKey === c.key) sortDir *= -1; else { sortKey = c.key; sortDir = 1; } drawTable(); });
      tr.appendChild(th);
    });
    t.appendChild(h('thead', {}, [tr]));
    const tb = h('tbody');
    rows.slice(0, 400).forEach(g => {
      const row = h('tr');
      row.appendChild(h('td', { text: g.gene }));
      row.appendChild(h('td', { text: g.contrast }));
      row.appendChild(h('td', { class: 'num', text: fmt(g.base) }));
      const lfc = h('td', { class: 'num', text: g.lfc.toFixed(2) });
      const dot = h('span', { class: 'key' });
      dot.style.background = g.lfc > 0 ? UP : DOWN; dot.style.marginRight = '6px';
      lfc.prepend(dot);
      row.appendChild(lfc);
      row.appendChild(h('td', { class: 'num', text: String(g.padj) }));
      row.appendChild(h('td', { text: (g.flag ? '★ ' : '') + (g.desc || '') }));
      tb.appendChild(row);
    });
    t.appendChild(tb);
    body.appendChild(t);
    if (!rows.length) body.appendChild(h('p', { class: 'empty', text: 'No gene passes the current filters.' }));
    else if (rows.length > 400) body.appendChild(h('p', { class: 'note-line', text: 'Showing the first 400 rows of the current filter.' }));
  }

  /* one volcano per contrast, sharing the filter row above - small multiples when "all" is chosen */
  function volcanoFor(c, small) {
    const ctx = { color: MUTED, r: small ? 1.7 : 2, opacity: .45 };
    const cls = p => {
      if (!isSig(p, f.padj, f.lfc)) return ctx;
      if (f.dir === 'up' && p[0] <= 0) return ctx;
      if (f.dir === 'down' && p[0] >= 0) return ctx;
      return { color: p[0] > 0 ? UP : DOWN, r: small ? 2.2 : 2.8 };
    };
    let up = 0, down = 0;
    c.points.forEach(p => { if (isSig(p, f.padj, f.lfc)) { p[0] > 0 ? up++ : down++; } });
    const maxY = Math.max(-Math.log10(f.padj) * 1.3, ...c.points.map(p => p[1]));
    return card({
      title: c.id,
      why: `${fmt(up)} up · ${fmt(down)} down at padj < ${f.padj} and |log2FC| ≥ ${f.lfc}`,
      chart: scatter({
        points: c.points, height: small ? 260 : 360,
        xLabel: 'log2 fold change (apeglm-shrunken)', yLabel: '−log10 padj',
        yMin: 0, yMax: maxY,
        vlines: f.lfc > 0 ? [{ x: -f.lfc }, { x: f.lfc }] : [],
        hlines: [{ y: -Math.log10(f.padj), label: 'padj ' + f.padj }],
        classOf: cls,
        tipTitle: p => (isSig(p, f.padj, f.lfc) ? (p[0] > 0 ? 'up at these thresholds' : 'down at these thresholds') : 'below the thresholds'),
        tipRows: p => [{ label: 'log2FC', value: p[0].toFixed(2) }, { label: 'padj', value: padjOf(p).toExponential(2) }]
      })
    });
  }

  function drawVolcano() {
    volcanoes.textContent = '';
    volcanoes.appendChild(legendRow([{ color: UP, label: 'up' }, { color: DOWN, label: 'down' },
    { color: MUTED, label: 'below the thresholds' }]));
    const list = f.contrast === 'all' ? D.contrasts : D.contrasts.filter(c => c.id === f.contrast);
    if (list.length > 1) {
      const grid = h('div', { class: 'grid g2' });
      list.forEach(c => grid.appendChild(volcanoFor(c, true)));
      volcanoes.appendChild(grid);
    } else if (list.length) {
      volcanoes.appendChild(volcanoFor(list[0], false));
    }
    const c0 = D.contrasts[0];
    volcanoes.appendChild(h('p', { class: 'note-line', text:
      `Every gene with padj < ${ceiling} is plotted, so the counts above are exact at any threshold up to it; `
      + `the ${fmt(c0.n_background)} genes further from significance are sampled down to ${fmt(c0.n_background_shown)} points for context and never counted. `
      + `The contrast, direction and threshold fields scope both the table and these plots; the text search narrows the table only.` }));
  }

  function draw() { drawTable(); drawVolcano(); }

  out.appendChild(h('h2', { text: 'Top differentially expressed genes — ' + activeMapLabel() }));
  out.appendChild(h('p', { class: 'sub', text: 'The strongest genes per contrast by adjusted p-value, joined to their annotation. Click a column to sort; ★ marks the configured family of interest. The thresholds below scope this table and the volcano plots underneath it.' }));
  out.appendChild(h('div', { class: 'controls' }, [
    csel, dsel,
    h('span', { class: 'why', text: 'padj <' }), padjIn,
    h('span', { class: 'why', text: '|log2FC| ≥' }), lfcIn,
    q, count
  ]));
  out.appendChild(h('div', { class: 'card' }, [body]));
  out.appendChild(h('h2', { text: 'Volcano' }));
  out.appendChild(h('p', { class: 'sub', text: 'The same slice as the table, plotted. Points outside the thresholds stay on the plot in grey — a volcano with its context removed hides how ordinary the significant genes are.' }));
  out.appendChild(volcanoes);
  draw();
  return out;
}

function sectionProvenance() {
  const out = h('div');
  const M = DATA.meta || {};
  out.appendChild(h('div', { class: 'card' }, [
    h('div', { class: 'figcap' }, [h('h3', { text: 'This run' })]),
    buildTable({
      cols: [{ key: 'k', label: 'Field' }, { key: 'v', label: 'Value' }],
      rows: [['run name', M.run_name], ['started', M.start], ['profile', M.profile], ['revision', M.revision],
      ['Nextflow', M.nextflow], ['results directory', M.outdir], ['reference set', M.reference_set],
      ['command line', M.command_line]].filter(r => r[1]).map(r => ({ k: r[0], v: r[1] }))
    })
  ]));

  const P = DATA.params || {};
  const keys = Object.keys(P);
  if (keys.length) {
    const body = h('div', { class: 'scroll' });
    const q = h('input', { type: 'search', placeholder: 'filter parameters', oninput: e => draw(e.target.value.toLowerCase()) });
    function draw(term) {
      body.textContent = '';
      const rows = keys.filter(k => !term || k.toLowerCase().includes(term) || String(P[k]).toLowerCase().includes(term))
        .map(k => ({ k: k, v: P[k] === null ? '—' : (typeof P[k] === 'object' ? JSON.stringify(P[k]) : String(P[k])) }));
      body.appendChild(buildTable({ cols: [{ key: 'k', label: 'Parameter' }, { key: 'v', label: 'Value' }], rows: rows }));
    }
    draw('');
    out.appendChild(h('h2', { text: 'Parameters' }));
    out.appendChild(h('p', { class: 'sub', text: 'Every value this run used. Re-running with the same parameter file reproduces it; the arguments quoted on the Methods tab come from here.' }));
    out.appendChild(h('div', { class: 'controls' }, [q]));
    out.appendChild(h('div', { class: 'card' }, [body]));
  }

  if ((DATA.related || []).length) {
    const c = h('div', { class: 'card' }, [h('div', { class: 'figcap' }, [
      h('span', { class: 'why', text: 'paths relative to this page' })])]);
    const t = h('table');
    const tb = h('tbody');
    DATA.related.forEach(r => {
      const a = h('a', { href: r.href, text: r.label });
      tb.appendChild(h('tr', {}, [h('td', {}, [a]), h('td', { text: r.note })]));
    });
    t.appendChild(tb); c.appendChild(t);
    out.appendChild(h('h2', { text: 'Related outputs' }));
    out.appendChild(c);
  }

  if ((DATA.inventory || []).length) {
    out.appendChild(h('h2', { text: 'Files behind this page' }));
    out.appendChild(h('div', { class: 'card' }, [buildTable({
      cols: [{ key: 'name', label: 'Output' }, { key: 'size', label: 'Size' }],
      rows: DATA.inventory
    })]));
  }
  return out;
}

/* ================================ page ================================ */
const TABS = [
  ['overview', 'Overview', sectionOverview],
  ['methods', 'Methods', sectionMethods],
  ['qc', 'Reads & QC', sectionQC],
  ['reference', 'Reference', sectionReference],
  ['annotation', 'Annotation', sectionAnnotation],
  ['expression', 'Expression', sectionExpression],
  ['genes', 'Genes', sectionGenes],
  ['provenance', 'Provenance', sectionProvenance],
];
const built = {};
function showTab(id) {
  TABS.forEach(([tid]) => {
    const panel = document.getElementById('panel-' + tid), btn = document.getElementById('tab-' + tid);
    const on = tid === id;
    btn.setAttribute('aria-selected', String(on));
    if (on) { panel.removeAttribute('hidden'); } else { panel.setAttribute('hidden', ''); }
  });
  if (!built[id]) {
    const [, , render] = TABS.find(t => t[0] === id);
    let node;
    try { node = render(); }
    catch (e) { node = h('p', { class: 'empty', text: 'This section could not be built: ' + e.message }); }
    document.getElementById('panel-' + id).appendChild(node);
    built[id] = true;
  }
  try { history.replaceState(null, '', '#' + id); } catch (e) { /* file:// */ }
}
function openStage(sid) {
  const d = document.getElementById('stage-' + sid);
  if (d) { d.open = true; if (d.scrollIntoView) d.scrollIntoView({ block: 'center', behavior: 'smooth' }); }
}

(function boot() {
  const nav = document.getElementById('tabs');
  TABS.forEach(([id, label]) => nav.appendChild(h('button', {
    id: 'tab-' + id, role: 'tab', 'aria-selected': 'false', 'aria-controls': 'panel-' + id,
    text: label, onclick: () => showTab(id)
  })));
  const main = document.getElementById('panels');
  TABS.forEach(([id]) => main.appendChild(h('section', { id: 'panel-' + id, role: 'tabpanel', 'aria-labelledby': 'tab-' + id, hidden: '' })));

  const themeBtn = document.getElementById('theme');
  let theme = null;
  try { theme = localStorage.getItem('dash-theme'); } catch (e) { }
  applyTheme(theme);
  themeBtn.addEventListener('click', () => {
    theme = theme === 'dark' ? 'light' : theme === 'light' ? null : 'dark';
    applyTheme(theme);
    try { theme ? localStorage.setItem('dash-theme', theme) : localStorage.removeItem('dash-theme'); } catch (e) { }
  });
  function applyTheme(t) {
    if (t) document.documentElement.setAttribute('data-theme', t);
    else document.documentElement.removeAttribute('data-theme');
    themeBtn.textContent = t === 'dark' ? 'Dark' : t === 'light' ? 'Light' : 'Auto';
  }

  const start = (location.hash || '').replace('#', '');
  showTab(TABS.some(t => t[0] === start) ? start : 'overview');
})();
"""

PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>__STYLE__</style>
</head>
<body>
<header class="top">
  <div class="wrap">
    <div class="title-row">
      <h1>__TITLE__</h1>
      <span class="runmeta">__RUNMETA__</span>
      <span class="spacer"></span>
      <button id="theme" class="ghost" title="Light / dark / follow the system">Auto</button>
    </div>
    <nav class="tabs" id="tabs" role="tablist"></nav>
  </div>
</header>
<main class="wrap" id="panels"></main>
<script type="application/json" id="run-data">__DATA__</script>
<script>__SCRIPT__</script>
</body>
</html>
"""

meta_bits = [b for b in [DATA["meta"]["reference_set"] and f"reference {DATA['meta']['reference_set']}",
                         DATA["meta"]["run_name"], DATA["meta"]["start"][:19].replace("T", " ") if DATA["meta"]["start"] else ""]
             if b]

data_json = json.dumps(DATA, separators=(",", ":"), default=str).replace("<", "\\u003c").replace("&", "\\u0026")
page = (PAGE
        .replace("__TITLE__", args.title)
        .replace("__RUNMETA__", " · ".join(meta_bits))
        .replace("__STYLE__", STYLE)
        .replace("__DATA__", data_json)
        .replace("__SCRIPT__", SCRIPT + SCRIPT2 + SCRIPT3))

os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
with open(args.out, "w") as fh:
    fh.write(page)

print(f"[dashboard] {args.out}  ({len(page) / 1024:.0f} KB)")
print(f"[dashboard] stages {len(stages)} · samples {len(samples)} · reference variants {len(ladder)}"
      f" · contrasts {len(dge_primary['contrasts']) if dge_primary else 0}"
      f" · primary map {primary_map} · maps with DE {sorted(k for k, v in dge_by_map.items() if v)}"
      f" · enrichment {sorted(enrichment_by_map)} · figures {sorted(figures_by_map)}")
