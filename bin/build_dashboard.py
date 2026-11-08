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
import array
import base64
import csv
import glob
import gzip
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

# The expression-filtered set is a rung too, and the most consequential one: it is the set that is
# annotated, quantified and tested. route_contribution.tsv is written before the filter runs, so it
# cannot carry it -- but SEQ_STATS_EXPR and BUSCO_EXPR both publish for it. Appending it here is
# what makes the assembly and the analysed set comparable on one line each, which is the whole
# point of reporting completeness on the assembly and analysing the subset.
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
        "basis": "expression-filtered (annotated, quantified and tested)",
        "transcripts": i(_rows[0].get("num_seqs")),
        "bp": i(_rows[0].get("sum_len")),
        "n50": None,
        "genes": None,
        "busco": parse_busco(open(_bp).read()) if exists(_bp) else None,
    })

# ------------------------------------------------- gene-map evaluation (optional, if staged)
# Which transcript->gene map is primary is a parameter, and the pipeline cannot settle it from its
# own outputs: deciding needs INDEPENDENT truth (marker orthologs, a reference genome's gene
# models) and a label-permutation null. Those come from an external evaluation, so the dashboard
# renders this card only when that evaluation is staged beside the run's own tables. Any project
# that supplies the same two files gets the card; nothing here is specific to one species.
#
#   genemap_eval.tsv : one row per (map, truth)  -- map, truth, F1, ARI, split_pct, merge_pct,
#                      genes_per_label, fused_pct  (extra columns are ignored)
#   genemap_null.tsv : one row per map           -- map, false_calls_mean, n_perms
genemap_eval = []
_ge = first_glob("genemap_eval.tsv") or first_glob("*/genemap_eval.tsv")
if _ge:
    for r in read_tsv(os.path.basename(_ge)) or []:
        def _f(k):
            v = (r.get(k) or "").strip()
            try:
                return float(v)
            except ValueError:
                return None
        if r.get("map") and r.get("truth"):
            genemap_eval.append({
                "map": r["map"].strip(), "truth": r["truth"].strip(),
                "f1": _f("F1"), "ari": _f("ARI"), "split": _f("split_pct"),
                "merge": _f("merge_pct"), "gpl": _f("genes_per_label"),
                "fused": _f("fused_pct"),
            })

genemap_null = {}
_gn = first_glob("genemap_null.tsv") or first_glob("*/genemap_null.tsv")
if _gn:
    for r in read_tsv(os.path.basename(_gn)) or []:
        if r.get("map"):
            try:
                genemap_null[r["map"].strip()] = {
                    "false_calls": float(r["false_calls_mean"]),
                    "n_perms": int(float(r.get("n_perms") or 0)),
                }
            except (ValueError, KeyError):
                pass

genemap_verdict = ""
_gv = first_glob("genemap_verdict.txt") or first_glob("*/genemap_verdict.txt")
if _gv:
    genemap_verdict = open(_gv).read().strip()

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

# ------------------------------------------- the expression matrices, for the Explorer tab
# Everything above reports numbers the pipeline already computed. This block ships the matrices
# THEMSELVES so the page can answer a question nobody wrote a figure for: what does gene X do
# across the 15 libraries, which transcripts of it carry the signal, and does a call survive
# dropping the fold-change shrinkage.
#
# Three decisions keep that from bloating the page:
#   1. Quantise. Expression goes in as uint16 log2(value + 1) * 1000 — 0.001 of a log2 unit, far
#      finer than the measurement — and DE statistics as int16 at 1/100, with -32768 for "not
#      tested". Raw doubles as JSON would be ~8x this and carry no extra information.
#   2. Dictionary-encode the annotation. 44,461 transcripts share ~12,500 distinct descriptions,
#      so each row stores an index into one string table.
#   3. gzip + base64 every array, and decompress in the browser on first use. The page stays one
#      file that works from file://, and the Explorer's cost is only paid if the tab is opened.
#
# Row order is the canonical index for a map: normalized_counts.csv's order, which DESeq2 also
# wrote the DE tables in, so one gene is one row across counts and every contrast without a join.
EXPR_SCALE = 1000          # uint16 counts = log2(x + 1) * EXPR_SCALE
# The DE statistics are int32, not int16, and the reason is the decision boundary. The page
# recomputes "significant" live from these arrays, under the same rule deseq2.R applied
# (padj < alpha AND |log2FC| > lfc_cut). At 1/100 of a log2 unit a gene whose fold change is
# really 1.0049 is stored as 1.00 and quietly stops being significant, and the page's count then
# disagrees with DE_summary.csv by a gene or two — measured on this run, not hypothetical. 1e-4 on
# the fold change and 1e-5 on -log10(padj) put the rounding far below the width of the decision.
LFC_SCALE = 10000          # int32: log2 fold change, shrunken and unshrunken
MLP_SCALE = 100000         # int32: -log10(padj), the same 5 decimals load_dge() settled on
BASE_SCALE = 100           # int32: log2(baseMean + 1), only ever displayed
DE_NA = -2147483648        # "DESeq2 did not test this gene in this contrast"
_I32_MAX = 2147483000      # leaves DE_NA unambiguous


def _gz_b64(raw):
    return base64.b64encode(gzip.compress(raw, 9)).decode("ascii")


# The browser reads these as Uint16Array / Int16Array / Int32Array / Uint8Array, so the widths
# have to be exactly 2/2/4/1 bytes. CPython's `array` sizes 'i' and 'I' by the C int, so assert
# rather than assume: a platform where that differs would produce a silently mis-parsed page.
_WIDTH = {"H": 2, "h": 2, "i": 4, "I": 4, "B": 1}


def _enc(typecode, values):
    """A typed array as gzip+base64, always little-endian (what every browser reads)."""
    a = array.array(typecode, values)
    if a.itemsize != _WIDTH[typecode]:
        raise SystemExit(f"[ERROR] array '{typecode}' is {a.itemsize} bytes here, "
                         f"not {_WIDTH[typecode]}: the dashboard's typed arrays would misparse")
    if sys.byteorder == "big":
        a.byteswap()
    return _gz_b64(a.tobytes())


def _enc_text(lines):
    return _gz_b64(("\n".join(lines)).encode("utf-8"))


def _q_expr(v):
    q = int(round(math.log2(max(v, 0.0) + 1.0) * EXPR_SCALE))
    return 0 if q < 0 else (65535 if q > 65535 else q)


def _q_de(v, scale):
    q = int(round(v * scale))
    return -_I32_MAX if q < -_I32_MAX else (_I32_MAX if q > _I32_MAX else q)


# ---- the string table every row's description points into
desc_list, desc_id = [], {}


def _desc_ix(text):
    if not text:
        return -1
    ix = desc_id.get(text)
    if ix is None:
        ix = desc_id[text] = len(desc_list)
        desc_list.append(text)
    return ix


# transcript -> (description index, interest flag), read from the transcript-level annotation
tx_desc, tx_flag = {}, {}
if annot_path:
    with open(annot_path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = r.get("transcript")
            if not t:
                continue
            d = (r.get("protein_desc") or r.get("eggnog_desc") or "").split(" OS=")[0][:90]
            if d:
                tx_desc[t] = _desc_ix(d)
            if r.get(interest_label):
                tx_flag[t] = 1


def _read_matrix(csv_path):
    """normalized_counts.csv / vst_counts.csv -> (sample names, row ids, quantised values)."""
    ids, vals = [], array.array("H")
    with open(csv_path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd, None)
        if not header:
            return None
        cols = [c for c in header[1:]]
        for row in rd:
            if not row or len(row) < 2:
                continue
            ids.append(row[0])
            for v in row[1:len(cols) + 1]:
                vals.append(_q_expr(num(v)))
    return cols, ids, vals


def _load_de_arrays(dirpath, index_of):
    """Per contrast: shrunken LFC, the unshrunken LFC beside it, -log10(padj) and baseMean,
    every array in the map's canonical row order so the page can switch between them."""
    out = {}
    n = len(index_of)
    for f in sorted(glob.glob(os.path.join(dirpath, "DE_*.csv"))):
        if os.path.basename(f) == "DE_summary.csv":
            continue
        name = os.path.basename(f)[3:-4]
        lfc = array.array("i", [DE_NA]) * n
        raw = array.array("i", [DE_NA]) * n
        mlp = array.array("i", [DE_NA]) * n
        base = array.array("i", [DE_NA]) * n
        seen = 0
        with open(f, newline="") as fh:
            rd = csv.DictReader(fh)
            key = (rd.fieldnames or [""])[0]
            has_raw = "log2FoldChange_unshrunken" in (rd.fieldnames or [])
            for r in rd:
                ix = index_of.get(r.get(key))
                if ix is None:
                    continue
                padj = r.get("padj")
                if padj in (None, "", "NA"):
                    continue          # left at DE_NA: independent filtering or a Cook's outlier
                seen += 1
                pa = max(num(padj, 1.0), 1e-300)
                lfc[ix] = _q_de(num(r.get("log2FoldChange")), LFC_SCALE)
                raw[ix] = (_q_de(num(r.get("log2FoldChange_unshrunken")), LFC_SCALE)
                           if has_raw else DE_NA)
                mlp[ix] = _q_de(-math.log10(pa), MLP_SCALE)
                base[ix] = _q_de(math.log2(num(r.get("baseMean")) + 1.0), BASE_SCALE)
        if not seen:
            continue
        out[name] = {
            "lfc": _enc("i", lfc), "mlp": _enc("i", mlp), "base": _enc("i", base),
            "raw": _enc("i", raw) if has_raw else None, "tested": seen,
        }
    return out


explorer = None
_expl_maps, _expl_samples, _expl_cols = {}, [], []
for _kind in MAP_KINDS:
    _dd = map_dir("dge", _kind)
    if not _dd:
        continue
    _m = _read_matrix(os.path.join(_dd, "normalized_counts.csv"))
    if not _m:
        continue
    _cols, _ids, _vals = _m
    if not _expl_cols:
        _expl_cols = _cols
    elif _cols != _expl_cols:
        # every map is quantified on the same libraries; if that ever stops being true the
        # matrices cannot share one column order, and shipping them anyway would mislabel columns
        print(f"[dashboard] explorer: {_kind} map has different libraries, skipped", file=sys.stderr)
        continue
    _index_of = {g: n for n, g in enumerate(_ids)}
    # transcript -> gene row, from the map this DESeq2 run actually used
    _t2g = os.path.join(_dd, "tx2gene.used.tsv")
    _gene_desc = [-1] * len(_ids)
    _gene_flag = [0] * len(_ids)
    _tx_of_gene = defaultdict(list)
    if os.path.exists(_t2g):
        with open(_t2g) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                ix = _index_of.get(p[1])
                if ix is None:
                    continue
                _tx_of_gene[ix].append(p[0])
                if _gene_desc[ix] < 0 and p[0] in tx_desc:
                    _gene_desc[ix] = tx_desc[p[0]]
                if p[0] in tx_flag:
                    _gene_flag[ix] = 1
    _expl_maps[_kind] = {
        "label": {"string": "String-stripped ids", "corset": "Corset",
                  "locus": "Locus rebuild"}[_kind],
        "n": len(_ids),
        "ids": _enc_text(_ids),
        "counts": _enc("H", _vals),
        "desc": _enc("i", _gene_desc),
        "flag": _enc("B", _gene_flag),
        "ntx": _enc("I", [len(_tx_of_gene.get(n, ())) for n in range(len(_ids))]),
        "de": _load_de_arrays(_dd, _index_of),
        "t2g_rows": sum(len(v) for v in _tx_of_gene.values()),
    }
    # the transcript -> gene row map, for the isoform view. Written once per map, against the
    # transcript order established below.
    _expl_maps[_kind]["_tx_of_gene"] = _tx_of_gene

# ---- transcript level: salmon's own TPM, the quantity the expression filter was applied to
_quant_dirs = sorted(
    d for d in glob.glob(path("*"))
    if os.path.isdir(d) and os.path.exists(os.path.join(d, "quant.sf"))
)
_tx_payload = None
if _quant_dirs and _expl_cols:
    _by_name = {os.path.basename(d): d for d in _quant_dirs}
    _missing = [c for c in _expl_cols if c not in _by_name]
    if _missing:
        print(f"[dashboard] explorer: no quant.sf for {', '.join(_missing[:3])}"
              f"{'...' if len(_missing) > 3 else ''}; transcript level omitted", file=sys.stderr)
    else:
        _tx_ids, _tx_ix = [], {}
        _per_sample = []
        for _c in _expl_cols:
            _vals_s = {}
            with open(os.path.join(_by_name[_c], "quant.sf"), newline="") as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                for r in rd:
                    nm = r["Name"]
                    if nm not in _tx_ix:
                        _tx_ix[nm] = len(_tx_ids)
                        _tx_ids.append(nm)
                    _vals_s[nm] = num(r.get("TPM"))
            _per_sample.append(_vals_s)
        _ntx, _ns = len(_tx_ids), len(_expl_cols)
        _tpm = array.array("H", [0]) * (_ntx * _ns)
        for _j, _vals_s in enumerate(_per_sample):
            for _nm, _v in _vals_s.items():
                _tpm[_tx_ix[_nm] * _ns + _j] = _q_expr(_v)
        _tx_payload = {
            "n": _ntx,
            "ids": _enc_text(_tx_ids),
            "tpm": _enc("H", _tpm),
            "desc": _enc("i", [tx_desc.get(t, -1) for t in _tx_ids]),
            "flag": _enc("B", [tx_flag.get(t, 0) for t in _tx_ids]),
        }
        # per map, the gene row each transcript belongs to (-1 = not in that map's model)
        for _kind, _mp in _expl_maps.items():
            _owner = [-1] * _ntx
            for _grow, _txs in _mp.pop("_tx_of_gene").items():
                for _t in _txs:
                    _jx = _tx_ix.get(_t)
                    if _jx is not None:
                        _owner[_jx] = _grow
            _mp["tx_gene"] = _enc("i", _owner)
for _mp in _expl_maps.values():
    _mp.pop("_tx_of_gene", None)

if _expl_maps:
    _cond_of = {}
    for _kind in MAP_KINDS:
        _d = dge_by_map.get(_kind)
        if _d:
            for _p in _d.get("pca", []):
                if _p.get("name") and _p.get("group"):
                    _cond_of.setdefault(_p["name"], _p["group"])
    _expl_samples = [{"name": c, "condition": _cond_of.get(c, "")} for c in _expl_cols]
    _conds = []
    for _s in _expl_samples:
        if _s["condition"] and _s["condition"] not in _conds:
            _conds.append(_s["condition"])
    explorer = {
        "scale": EXPR_SCALE, "na": DE_NA, "lfc_scale": LFC_SCALE,
        "mlp_scale": MLP_SCALE, "base_scale": BASE_SCALE,
        # the rule deseq2.R applied, so the page's live counts start where DE_summary.csv does
        "alpha": num(P("de_alpha", 0.05), 0.05) or 0.05,
        "lfc_cut": num(P("de_lfc", 1.0), 1.0),
        "samples": _expl_samples, "conditions": _conds,
        "desc_table": _enc_text(desc_list),
        "interest_label": interest_label,
        "maps": _expl_maps,
        "tx": _tx_payload,
        "primary": primary_map,
    }
    print(f"[dashboard] explorer: maps {sorted(_expl_maps)} · "
          f"genes {[(k, v['n']) for k, v in sorted(_expl_maps.items())]} · "
          f"transcripts {_tx_payload['n'] if _tx_payload else 0} · "
          f"descriptions {len(desc_list)}")


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

# ---- the gene-level FASTA each map ships, keyed by map
# A gene map is a table of ids: without this, a "gene" has counts and no sequence. GENE_FASTA
# writes one representative per gene, and its gene_rep.tsv is what says how many of them carry a
# peptide — the number that decides whether a protein-level use of the set is possible at all.
gene_seqs = {}
for _gr in sorted(glob.glob(path("*.gene_rep.tsv"))):
    try:
        with open(_gr) as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh][1:]
    except OSError:
        continue
    if not rows:
        continue
    _label = os.path.basename(_gr)[:-len(".gene_rep.tsv")]
    gene_seqs[map_kind_of(_label)] = {
        "file": os.path.basename(_gr).replace(".gene_rep.tsv", ".genes.fasta"),
        "genes": len(rows),
        "pep": sum(1 for r in rows if len(r) > 4 and r[4] not in ("NA", "")),
    }

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
    gene_maps.append({"key": "string", "primary": primary_map == "string", "name": "String-stripped ids", "basis": "per-source naming rules",
                      "genes": corset["stripped"], "sig": de_total(dge), "seqs": gene_seqs.get("string"),
                      "note": "the default map; cannot reconcile identity ACROSS assembly sources"})
if corset and corset.get("clusters"):
    gene_maps.append({"key": "corset", "primary": primary_map == "corset", "name": "Corset", "basis": "shared reads + expression",
                      "genes": corset["clusters"], "sig": de_total(dge_corset), "seqs": gene_seqs.get("corset"),
                      "note": "clusters contigs by the reads they share, not by sequence identity"})
if locus and locus.get("after"):
    gene_maps.append({"key": "locus", "primary": primary_map == "locus", "name": "Locus rebuild", "basis": "genomic coordinates",
                      "genes": locus["after"], "sig": de_total(dge_locus), "seqs": gene_seqs.get("locus"),
                      "note": "groups transcripts sharing a splice junction on the genome"})

# The gene count the headline quotes has to be the PRIMARY map's. The expression filter writes a
# gene count in its own log, but it counts the map it was handed — the string-stripped one — so
# quoting that next to a locus-primary run overstates the gene count by every split the locus
# rebuild reconciles (5,144 genes on this dataset, 34,853 against 29,709).
_primary_genes = next((m["genes"] for m in gene_maps if m.get("primary")), None)
_primary_name = next((m["name"] for m in gene_maps if m.get("primary")), None)

# The annotation module counts genes from the tx2gene it was HANDED — the string-stripped one,
# because annotation runs before the locus rebuild and the rebuild needs the genome placement.
# Every other row of that summary is a per-transcript count and is therefore identical under any
# gene map; the gene row is the only one that moves. So label it for what it is and put the primary
# map's own count beside it, instead of leaving a reader to assume it matches the headline.
if _primary_genes and _primary_name:
    for _k, _row in enumerate(annot_summary):
        if _row["label"].strip().lower() != "genes":
            continue
        _row["label"] = "genes (string-stripped)"
        if _row["value"] != _primary_genes:
            annot_summary.insert(_k + 1, {"label": f"genes ({_primary_name.lower()}, primary)",
                                          "value": _primary_genes})
        break

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
# the analysed set: the same rung after the expression filter, which is what is annotated and
# tested. Absent when the filter is off (expr_filter_tpm = 0), in which case ref_row IS the
# analysed set.
de_row = next((r for r in ladder if r["reference"] == f"{ref_set}_expressed"), None)

# ------------------------------------------------------------------ expression filter
# The filter log is the only place the kept/dropped accounting and the DE-eligible count are
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
        "dropped":   _ef_num(r"transcripts dropped\s*:\s*([\d,]+)"),
        "eligible":  _ef_num(r"DE-eligible transcripts\s*:\s*([\d,]+)"),
    }
    _g = re.search(r"genes\s+([\d,]+)\s*->\s*([\d,]+)", _txt)
    if _g:
        expr_filter["genes_before"] = int(_g.group(1).replace(",", ""))
        expr_filter["genes_after"] = int(_g.group(2).replace(",", ""))

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
    # BUSCO is quoted for the analysed set, falling back to the assembly only when no filter ran
    # (in which case they are the same set)
    _busco_row = de_row if (de_row and de_row.get("busco")) else ref_row
    _busco = _busco_row.get("busco")
    # lead with the primary map, like the headline; the string-stripped count follows as the
    # upper bound it is
    _an_genes = _primary_genes or expr_filter.get("genes_after") or genes
    _an_genes_src = _primary_name if _primary_genes else "string-stripped ids"
    block("transcriptome", "Reference transcriptome", "product",
          f"{ref_set} · {len(ladder)} variant{'s' if len(ladder) != 1 else ''} built and scored",
          [R("transcripts analysed",
             fmt(de_row["transcripts"]) if de_row else fmt(ref_row["transcripts"]),
             "after the expression filter: annotated, quantified and tested" if de_row
             else "no expression filter ran, so the assembly is the analysed set"),
           R("genes", fmt(_an_genes) if _an_genes else "",
             " · ".join(x for x in [
                 f"by {_an_genes_src.lower()}, the primary map" if _primary_genes else _an_genes_src,
                 (f"{fmt(expr_filter['genes_after'])} by string-stripped ids, an upper bound"
                  if expr_filter.get("genes_after")
                  and expr_filter.get("genes_after") != _an_genes else ""),
             ] if x)),
           # completeness for the set that was actually analysed; the assembly's own score sits on
           # the row below, so the two can never be read as one number
           R("BUSCO complete",
             f"{_busco['complete']:.1f}%" if _busco else "",
             (f"{P('busco_lineage', '')} · duplicated {_busco['dup']:.1f}% · measured on "
              f"{_busco_row['reference']}") if _busco else ""),
           R("the assembly, before the filter",
             fmt(ref_row["transcripts"]) if de_row else "",
             (" · ".join(x for x in [
                 f"N50 {fmt(ref_row['n50'])}" if ref_row.get("n50") else "",
                 (f"BUSCO {ref_row['busco']['complete']:.1f}% complete"
                  if ref_row.get("busco") else ""),
             ] if x)) if de_row else ""),
           R("N50", fmt(ref_row["n50"]) if (ref_row.get("n50") and not de_row) else ""),
           R("sources combined", str(len(ladder)), "each layered on and measured separately")],
          "reference")

# ---- annotation --------------------------------------------------------------------
if annot_evidence:
    tot = annot_evidence[0]["total"]
    ev = {e["label"]: e["n"] for e in annot_evidence}
    block("annotation", "Annotation", "product",
          f"evidence on {fmt(tot)} transcripts — the expression-filtered set",
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
# The analysed set leads, not the assembly. It is what was annotated, quantified and tested, so it
# is the denominator every other number on this page shares; heading the page with the assembly's
# larger count invites the reader to scale the DE and annotation counts against a set that was
# never tested. The assembly keeps its own row on the Reference tab, and is named here as the set
# this one was filtered out of.
if de_row:
    # the primary map's gene count, named, falling back to the filter's own (string-stripped) count
    # only when no gene map reported one — and saying which it is either way, because the two
    # differ by thousands and an unlabelled gene count beside a named primary map reads as that
    # map's
    _hero_genes = _primary_genes if _primary_genes else expr_filter.get("genes_after")
    _hero_genes_src = _primary_name if _primary_genes else "string-stripped ids"
    hero = {"value": fmt(de_row["transcripts"]),
            "label": "transcripts analysed — annotated, quantified and tested",
            "sub": " · ".join(x for x in [
                de_row["reference"],
                f"{fmt(_hero_genes)} genes ({_hero_genes_src})" if _hero_genes else "",
                (f"BUSCO {de_row['busco']['complete']:.1f}% complete"
                 if de_row.get("busco") else ""),
                (f"filtered from {fmt(ref_row['transcripts'])} assembled" if ref_row else ""),
            ] if x)}
elif ref_row:
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

# The clustering identity thresholds. Renamed cluster_* in 2.0.1dev because they had been named
# for CD-HIT, which this pipeline does not run; the cdhit_* fallback keeps a results directory
# written by an earlier run readable, since this page is rebuilt from the params json THAT run
# recorded, not from the current config.
_CLUSTER_LR_ID    = P("cluster_longread_id", P("cdhit_longread_id", 0.99))
_CLUSTER_DEDUP_ID = P("cluster_dedup_id", P("cdhit_dedup_id", 0.95))

stage("routeb", "Long reads, reference-free", "Assembly", "vsearch",
      on("skip_longread") and on("skip_reffree"),
      "Collapses near-identical isoforms without consulting the genome, so loci missing from a draft "
      "assembly survive. vsearch rather than cd-hit-est: the identity definition is stated rather "
      "than implied, where cd-hit measures matches over the shorter sequence without saying so.",
      [f"vsearch --cluster_fast --id {_CLUSTER_LR_ID}", f"prefix = {P('prefix_lrfree', 'PBR_')}"],
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
       f"vsearch --cluster_fast --id {_CLUSTER_DEDUP_ID} to collapse it against itself",
       f"then minimap2 -x asm20: already represented when an alignment covers "
       f"{P('novelty_min_cov', 0.9)} of a transcript at identity {_CLUSTER_DEDUP_ID}",
       f"prefix = {P('prefix_denovo', 'TR_')}"],
      [M("added", added("isoseq_stringtie_trinity", "isoseq_stringtie")),
       M("assembly-only libraries", fmt(len(assembly_reads)) if assembly_reads else "—")],
      "assembly/denovo/")

stage("dedup", "Cross-source deduplication", "Assembly", "vsearch", True,
      "The layered set counts one biological locus once per source that recovered it. Deduplication "
      "reduces that; it does not eliminate it — which is what the locus rebuild below addresses.",
      [f"vsearch --cluster_fast --id {_CLUSTER_DEDUP_ID}"],
      [M("the assembly (before the filter)", fmt(ref_row["transcripts"]) if ref_row else "—"),
       M("BUSCO complete", f"{ref_row['busco']['complete']:.1f}%" if ref_row and ref_row["busco"] else "—")],
      "transcriptome/")

stage("filter", "Expression filter", "Expression", "salmon pass 1 + filter_by_expression.py", True,
      "The ladder applies no expression evidence at any rung, so the deduplicated set carries every "
      "unexpressed fragment the four sources produced. A first salmon pass over that set measures "
      "what the experiment actually saw, and the reference is cut to the survivors. Everything after "
      "this point runs on that one set: it is what gets annotated, re-quantified, tested and drawn, "
      "so no number on this dashboard has to be qualified with which reference it came from. What "
      "the filter does NOT touch is the assembly — every ladder rung keeps its own FASTA, stats and "
      "BUSCO summary, and the assembly's own count and completeness stay on the Reference tab. The "
      "headline count and the headline BUSCO are the analysed set's, because that is the set every "
      "other number on this page describes.",
      [f"expressed: TPM >= {P('expr_filter_tpm', 1)} in >= {P('expr_filter_samples', 3)} libraries",
       "a presence rule, not a mean: a transcript on in one group and off elsewhere survives",
       "nothing is rescued on homology: measured on this data that added 0 DE genes and 0.13 pts "
       "of mapping, while costing 4.2 pts of BUSCO duplication"],
      [M("the analysed set", fmt(expr_filter["kept"]) if expr_filter.get("kept") else "—"),
       M("dropped", fmt(expr_filter["dropped"]) if expr_filter.get("dropped") else "—"),
       # The filter is handed the string-stripped map, so its log reports that map's count and not
       # the primary map's — and the two differ by thousands. The set the filter KEEPS is the
       # analysed set, which is exactly the set the locus rebuild is defined on, so the primary
       # map's count for it is available and is shown beside it rather than left to be inferred.
       M("genes (string-stripped)",
         fmt(expr_filter["genes_after"]) if expr_filter.get("genes_after") else "—")]
      + ([M(f"genes ({_primary_name.lower()}, primary)", fmt(_primary_genes))]
         if _primary_genes and _primary_name
         and _primary_genes != expr_filter.get("genes_after") else [])
      + [M("BUSCO complete", f"{de_row['busco']['complete']:.1f}%" if de_row and de_row["busco"] else "—")],
      "transcriptome/")

stage("annot", "Functional annotation", "Annotation",
      "TransDecoder · DIAMOND · eggNOG-mapper · InterProScan", on("skip_annotation"),
      "ORFs first, then homology and domains. It runs on the expression-filtered set, so this table "
      "and the DE tables describe the same transcripts: an annotation row here is a transcript the "
      "experiment saw. Column has_ORF is yes/no for every row — the ORF count is "
      "coding_status == \"coding\", not the number of non-empty cells.",
      ["TransDecoder.Predict --single_best_only",
       "diamond blastp -k 1 -e 1e-5 --more-sensitive",
       "eggnog-mapper -m diamond" if on("skip_eggnog") else "eggNOG-mapper: skipped",
       f"interproscan -appl {P('interproscan_appl', '')}" if on("skip_interproscan") else "InterProScan: skipped"],
      [M(a["label"], fmt(a["value"])) for a in annot_summary[:6]],
      "annotation/")

stage("quant", "Quantification", "Expression", "salmon", True,
      "Selective alignment against a decoy-aware index (the genome is the decoy), which stops reads "
      "from genomic or unassembled regions being forced onto a transcript. Library type is "
      "auto-detected and reported, not assumed. "
      "NOT bit-reproducible, and this is the limit on how exactly any number downstream of it can "
      "be repeated: `salmon index -p 16` builds pufferfish's position table in parallel, so from a "
      "byte-identical reference FASTA and an identical command it produces a different seq.bin / "
      "pos.bin / ctable.bin (every summary figure in its own log identical). Quantification then "
      "redistributes multi-mapping reads between near-identical transcripts. Measured across two "
      "runs of this reference: total assigned reads per library exactly equal, 27% of transcripts "
      "given a different count, but only 0.5% moved by more than 1% and 0.04% by more than half. "
      "A handful move enough to carry a gene over the significance cut, so a re-run reproduces the "
      "DE counts to within a gene or two rather than exactly.",
      ["salmon index -d <genome decoys> -k 31",
       "salmon quant -l A --gcBias --seqBias --dumpEq --numBootstraps 30",
       "the reported pass; the earlier pass on the deduplicated set only fed the filter"],
      [M("mean mapping", f"{sum(mapped) / len(mapped):.2f}%" if mapped else "—"),
       M("library type", samples[0]["library_type"] if samples and samples[0]["library_type"] else "—"),
       M("fragments processed", fmt(frag) if frag else "—")],
      "quantification/")

stage("dge", "Differential expression", "Expression", "tximport + DESeq2", on("skip_dge"),
      "Gene-level counts aggregated with tximport, then DESeq2 on the configured design. Every gene "
      "in the reference is testable, because the reference is already the expressed set — the "
      "DE-eligible list is passed through and applied all the same, so a change to how the "
      "reference is built cannot silently widen what gets tested. Fold changes are "
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
      "the split-gene problem and roughly doubles gene fusion — both directions are measured below. "
      "A single-exon transcript has no junction to share and so cannot be grouped that way: it is "
      "attached to the multi-exon locus it sits inside, and single-exon transcripts with no such "
      "host are then clustered among themselves by overlap. A transcript that fails the coverage "
      "cut, or that does not align at all, is kept as a gene of its own rather than dropped.",
      [f"grouping = {P('locus_groupby', 'junction')} — link two transcripts that share an identical "
       "(scaffold, strand, donor, acceptor) junction, then take connected components",
       f"minimum transcript coverage = {P('locus_mincov', 0.8)} — fraction of the transcript that "
       "must align before it is placed at all",
       f"single-exon attach threshold = {P('locus_se_overlap', 0.5)} — a junctionless transcript "
       "joins a locus when at least this fraction of ITS OWN genomic length is covered by that "
       "locus's exons on the same strand. Lower it and fragments get swept into neighbouring "
       "genes; raise it and each becomes its own gene"],
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
        # the set that was annotated, quantified and tested: the same rung after the
        # expression filter, or the rung itself when the filter is off
        "analysed_set": de_row["reference"] if de_row else ref_set,
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
                   "note": ("Built on the expression-filtered reference, so this table and the DE "
                            "tabs describe the SAME transcripts: every row is a transcript the "
                            "experiment saw at TPM >= the filter threshold, and an unexpressed "
                            "gene model from the assembly has no row here. That is deliberate — it "
                            "makes an annotation count and a DE count comparable without a "
                            "footnote — but it also means this table is not a catalogue of what "
                            "the species encodes: the assembly holds more gene models than this. "
                            "BUSCO is reported for this set, with the assembly's own completeness "
                            "kept beside it on the Reference tab. Every row below except the gene "
                            "count is a count of TRANSCRIPTS, so it is the same number under any "
                            "gene map. The gene row is the exception: annotation runs before the "
                            "locus rebuild and is keyed on the string-stripped map, so the "
                            "module's own gene figure is the string-stripped one and the primary "
                            "map's count is listed beside it. The locus and Corset variants are "
                            "re-keyed copies of this same table.")},
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
    "explorer": explorer,
    "genemap_eval": genemap_eval,
    "genemap_null": genemap_null,
    "genemap_verdict": genemap_verdict,
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
.linkish{color:var(--s1);font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;
  text-align:left;padding:0;text-decoration:underline;text-decoration-style:dotted}
.linkish:hover{text-decoration-style:solid}
.controls label.why{display:inline-flex;align-items:center;gap:5px}
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
  const used = (DATA.ladder || []).find(r => r.reference === DATA.meta.analysed_set)
            || (DATA.ladder || []).find(r => r.reference === DATA.meta.reference_set);
  if (used && used.busco && used.busco.dup > 20)
    caveats.push(`BUSCO duplication is ${used.busco.dup.toFixed(1)}% on ${used.reference}, the set that was annotated, quantified and tested: on a de novo catalogue that usually means one locus is represented several times, not that the genome is duplicated.`);
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
   expression filter sits ABOVE the annotation, so one set of transcripts is annotated, quantified
   and tested while the completeness quoted for the assembly belongs to the larger set the filter
   cut from; and the three gene maps are three re-aggregations of ONE quantification, not three
   analyses, so they are drawn as parallel branches with their own DESeq2 and their own enrichment
   rather than folded into a single box that implies one answer. Boxes are numbered in execution
   order so prose can point at one. */
// step numbers of the boxes the prose refers to; set when the diagram is drawn (an optional input
// box shifts them)
let FLOW_N = { annot: 10, filt: 9 };
function flowDiagram() {
  const W = 780, H = 772;
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label',
    'Flow: reads are trimmed, rRNA-depleted and decontaminated, then assembled by four routes ' +
    '(any assembly-only libraries join the two short-read routes and nothing else) ' +
    'into a reference ladder; the ladder is deduplicated, and that full set is what is scored for ' +
    'completeness. A first salmon pass measures expression, the expression filter cuts the ' +
    'reference down to what the experiment saw, and that filtered set is then annotated and ' +
    'quantified by a second salmon pass. The run then splits into three parallel branches — one ' +
    'per gene map — each with its own DESeq2 test and its own GO/KEGG enrichment; one is marked ' +
    'primary and the other two are sensitivity arms. Boxes are numbered in execution order.');

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

  /* The assembly is now drawn as what it is: the set the completeness number belongs to, and
     nothing else. Dashed and muted, because no reported annotation or expression number comes
     from it — the filter below it decides what the rest of the run sees. */
  const asm   = box(218, 252, 512, 46, 'The assembly  ·  ' + refName, MUTED,
                    { dashed: true,
                      lines: [[refRow && refRow.transcripts ? fmtN(refRow.transcripts) + ' transcripts' : null,
                               busco(refRow)].filter(Boolean).join('  ·  ')
                              || 'scored for completeness; not annotated'] });
  const q1    = box(218, 330, 512, 46, 'salmon pass 1  ·  ' + refName, MUTED,
                    { lines: ['measures what the experiment saw — TPM only, nothing here is reported'] });

  const filt  = box(218, 408, 512, 62, 'Expression filter  →  the analysed set', S3,
                    { lines: [F.rule || 'TPM ≥ 1 in ≥ 3 libraries',
                              [F.kept ? fmtN(F.kept) + ' kept, ' + fmtN(F.dropped) + ' dropped' : null,
                               busco(deRow)].filter(Boolean).join('  ·  ')
                              || 'everything below runs on what it keeps'] });
  const annot = box(30, 500, 310, 62, 'Annotation  ·  ' + (DATA.meta.analysed_set || refName), S2,
                    { lines: ['TransDecoder → DIAMOND · eggNOG · InterProScan',
                              F.kept ? fmtN(F.kept) + ' transcripts — the same set as the DE tables'
                                     : 'the filtered set, not the assembly'] });
  const q2    = box(380, 500, 350, 62, 'salmon pass 2  ·  one quantification', S1,
                    { lines: ['decoy-aware, --dumpEq, 30 bootstraps', 'the counts every DE number comes from'] });

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
    deBoxes.push(box(x, 598, w, 62, m.label, col, {
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
    enrBoxes.push(box(x, 692, w, 58, 'Enrichment', m.primary ? S3 : MUTED, {
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
  // one chain, top to bottom: the assembly is quantified once to measure expression, the filter
  // cuts it to what the experiment saw, and that set is what is annotated AND what is quantified
  // for the record — which is the ordering the rest of the page assumes
  elbow(ladderB, asm);
  elbow(asm, q1);
  elbow(q1, filt, 'TPM');
  elbow(filt, annot, 'the filtered set');
  arrow(filt, q2);
  FLOW_N = { annot: annot.n, filt: filt.n };
  deBoxes.forEach(b => elbow(q2, b, undefined, 580));
  deBoxes.forEach((b, k) => { if (enrBoxes[k]) arrow(b, enrBoxes[k]); });
  rail(annot, enrBoxes, 66, 576, 676, 'terms for enrichment');

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
  add('Long-read clustering (route B)', `vsearch --id ${p('cluster_longread_id', p('cdhit_longread_id', 0.99))}`, 'transcript',
      'collapses near-identical long-read isoforms');
  add('Cross-source deduplication', `vsearch --id ${p('cluster_dedup_id', p('cdhit_dedup_id', 0.95))}`, 'transcript',
      F.genes_before ? `${n(F.genes_before)} genes (string-stripped: the locus map needs genome placement and is built for the analysed set only) enter the expression filter`
                     : 'produces all_sources_dedup');
  add('Novel-vs-reference call', 'minimap2 -x asm20 -c --secondary=no', 'transcript',
      'exact match counts, not the approximate ones -x alone gives');

  // There is one way into the analysed set and it is expression, so "kept" and "DE-eligible" are
  // the same transcripts. The two rows below say what that set feeds and what it leaves alone,
  // which is the pair readers get wrong.
  add('Expression filter — the analysed set', F.rule ||
      `TPM >= ${p('expr_filter_tpm', 1)} in >= ${p('expr_filter_samples', 3)} libraries`,
      'transcript',
      F.kept ? `${n(F.kept)} kept, ${n(F.dropped)} dropped` : '—');
  add('→ what it feeds', 'the annotation, salmon pass 2, DESeq2, enrichment and the figures',
      'transcript', 'everything reported past this row describes the kept transcripts and nothing else');
  add('→ what it does NOT touch', 'the assembly: it keeps its own FASTA, count and BUSCO',
      'transcript', 'every rung is still scored separately — see the Reference tab');

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
          'and every number reported as a quantification comes from the second. The expression ' +
          `filter (box ${FLOW_N.filt}) sits ABOVE the annotation (box ${FLOW_N.annot}), so one set of ` +
          'transcripts is annotated, quantified and tested and the Annotation tab describes the ' +
          'same transcripts as the DE tabs; what the assembly keeps that this set does not is a ' +
          'larger transcript count and its own completeness score, both kept on the Reference ' +
          'tab beside the analysed set the headline quotes. The run then ' +
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
    note: 'Read top to bottom: each row narrows what the row below it sees. The row that decides ' +
          'the most is the expression filter: everything below it — the annotation, the reported ' +
          'quantification, DESeq2, enrichment — sees only what it kept, and the headline count ' +
          'and BUSCO are measured on that kept set. The full deduplicated assembly is scored ' +
          'above the filter and kept beside it on the Reference tab, where its larger transcript ' +
          'count describes more transcripts than any other number on this page.'
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
      why: 'the assembly is ' + (DATA.meta.reference_set || '') + '; the expression-filtered rung below it is what the annotation and every DE number run on',
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
    'Source: the ' + (DATA.meta.analysed_set || A.reference || 'selected') + ' transcript set. ' + A.note }));
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

SCRIPT4 = r"""
/* ========================= Explorer: the matrices themselves =========================
 * Every other tab reports a number the pipeline computed. This one ships the expression
 * matrix and the DE statistics and lets the reader ask their own question of them: what
 * does this gene do across the 15 libraries, which of its transcripts carries that, and
 * does the call survive dropping the fold-change shrinkage.
 *
 * The arrays arrive gzip+base64 as quantised typed arrays and are inflated on first use,
 * so a reader who never opens this tab pays nothing for it.
 */
const EXPL = DATA.explorer || null;
const XC = {};                                   /* inflate-once cache, keyed by caller */

function b64bytes(b64) {
  const bin = atob(b64), u = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
  return u;
}
async function gunzip(u8) {
  const stream = new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}
/* A typed view over the inflated bytes. Response.arrayBuffer() hands back a buffer this
 * array owns outright, so the view can start at 0 and span the whole thing. */
async function xarr(key, b64, Ctor) {
  if (XC[key]) return XC[key];
  const raw = await gunzip(b64bytes(b64));
  return (XC[key] = new Ctor(raw.buffer, 0, raw.byteLength / Ctor.BYTES_PER_ELEMENT));
}
async function xlines(key, b64) {
  if (XC[key]) return XC[key];
  const raw = await gunzip(b64bytes(b64));
  return (XC[key] = new TextDecoder().decode(raw).split('\n'));
}

const NS = EXPL ? EXPL.samples.length : 0;
const CONDS = EXPL ? (EXPL.conditions.length ? EXPL.conditions : ['all']) : [];
/* which columns belong to each condition, in matrix column order */
const COND_COLS = EXPL ? CONDS.map(c => EXPL.samples
  .map((s, i) => [s, i]).filter(([s]) => (s.condition || 'all') === c).map(([, i]) => i)) : [];
const COND_COLOR = [S1, S2, S3, S4];

/* log2(x+1), the unit everything is stored and plotted in */
const lg = (M, row, j) => M.counts[row * NS + j] / EXPL.scale;
const unlg = v => Math.pow(2, v) - 1;
/* Each DE array carries its own scale; DE_NA means DESeq2 returned no p-value for that gene. */
const LFCS = EXPL ? EXPL.lfc_scale : 1, MLPS = EXPL ? EXPL.mlp_scale : 1;
const deV = (arr, row, scale) => (!arr || arr[row] === EXPL.na) ? null : arr[row] / scale;
const deL = (arr, row) => deV(arr, row, LFCS);          /* a log2 fold change */
const deP = (arr, row) => deV(arr, row, MLPS);          /* -log10(adjusted p) */
/* The fold-change cut is STRICT, because deseq2.R counts `log2FoldChange > lfc_cut`. Using >=
 * here would make this page disagree with DE_summary.csv for any gene sitting exactly on it. */
const overCut = (v, cut) => v !== null && Math.abs(v) > cut;

async function loadLevel(level, mapKey) {
  const key = level === 'tx' ? 'tx' : 'g:' + mapKey;
  if (XC['LVL:' + key]) return XC['LVL:' + key];
  const out = { level: level, key: key };
  if (level === 'tx') {
    const T = EXPL.tx;
    out.n = T.n;
    out.ids = await xlines(key + ':ids', T.ids);
    out.counts = await xarr(key + ':c', T.tpm, Uint16Array);
    out.desc = await xarr(key + ':d', T.desc, Int32Array);
    out.flag = await xarr(key + ':f', T.flag, Uint8Array);
    out.de = {};
    out.unit = 'TPM';
  } else {
    const M = EXPL.maps[mapKey];
    out.n = M.n;
    out.ids = await xlines(key + ':ids', M.ids);
    out.counts = await xarr(key + ':c', M.counts, Uint16Array);
    out.desc = await xarr(key + ':d', M.desc, Int32Array);
    out.flag = await xarr(key + ':f', M.flag, Uint8Array);
    out.ntx = await xarr(key + ':n', M.ntx, Uint32Array);
    out.txGene = M.tx_gene ? await xarr(key + ':tg', M.tx_gene, Int32Array) : null;
    out.de = {};
    for (const c in (M.de || {})) {
      const d = M.de[c];
      out.de[c] = {
        lfc: await xarr(key + ':' + c + ':l', d.lfc, Int32Array),
        raw: d.raw ? await xarr(key + ':' + c + ':r', d.raw, Int32Array) : null,
        mlp: await xarr(key + ':' + c + ':p', d.mlp, Int32Array),
        base: await xarr(key + ':' + c + ':b', d.base, Int32Array),
        tested: d.tested
      };
    }
    out.unit = 'normalised count';
  }
  return (XC['LVL:' + key] = out);
}
async function descTable() { return xlines('DESCS', EXPL.desc_table); }

/* ---------- getting a plot out of the page ---------- */
/* The charts paint with CSS custom properties, which mean nothing in a standalone file, so
 * every var() is resolved against the live document before the SVG leaves the page. */
function inlineVars(node) {
  const cs = getComputedStyle(document.documentElement);
  const resolve = v => String(v).replace(/var\((--[\w-]+)\)/g,
    (_, n) => (cs.getPropertyValue(n) || '').trim() || '#888');
  node.querySelectorAll('*').forEach(el => {
    ['style', 'fill', 'stroke'].forEach(a => {
      const v = el.getAttribute(a);
      if (v && v.indexOf('var(') >= 0) el.setAttribute(a, resolve(v));
    });
  });
  node.querySelectorAll('text').forEach(el => {
    const c = getComputedStyle(el);
    if (!el.getAttribute('fill')) el.setAttribute('fill', c.fill);
    el.setAttribute('font-size', c.fontSize);
    el.setAttribute('font-family', 'system-ui, -apple-system, Segoe UI, Roboto, sans-serif');
  });
  return node;
}
function svgSource(svg) {
  const c = svg.cloneNode(true);
  c.setAttribute('xmlns', SVGNS);
  const vb = (c.getAttribute('viewBox') || '0 0 760 400').split(/\s+/);
  c.setAttribute('width', vb[2]); c.setAttribute('height', vb[3]);
  inlineVars(c);
  const bg = s('rect', { x: 0, y: 0, width: vb[2], height: vb[3],
    fill: (getComputedStyle(document.body).backgroundColor || '#fff') });
  c.insertBefore(bg, c.firstChild);
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(c);
}
function saveBlob(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
}
function downloadSVG(svg, name) { saveBlob(new Blob([svgSource(svg)], { type: 'image/svg+xml' }), name + '.svg'); }
function downloadPNG(svg, name, scale) {
  const src = svgSource(svg), vb = (svg.getAttribute('viewBox') || '0 0 760 400').split(/\s+/);
  const k = scale || 3, img = new Image();
  img.onload = () => {
    const cv = document.createElement('canvas');
    cv.width = vb[2] * k; cv.height = vb[3] * k;
    const g = cv.getContext('2d');
    g.scale(k, k); g.drawImage(img, 0, 0);
    cv.toBlob(b => b && saveBlob(b, name + '.png'), 'image/png');
  };
  img.onerror = () => downloadSVG(svg, name);
  img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(src)));
}
function exportRow(svgHolder, name, tsv) {
  return h('div', { class: 'controls', style: 'margin:10px 0 0' }, [
    h('button', { class: 'ghost', text: 'Download PNG',
      onclick: () => { const g = svgHolder.querySelector('svg'); if (g) downloadPNG(g, name, 3); } }),
    h('button', { class: 'ghost', text: 'Download SVG',
      onclick: () => { const g = svgHolder.querySelector('svg'); if (g) downloadSVG(g, name); } }),
    tsv ? h('button', { class: 'ghost', text: 'Download data (TSV)',
      onclick: () => saveBlob(new Blob([tsv()], { type: 'text/tab-separated-values' }), name + '.tsv') }) : null,
    tsv ? h('button', {
      class: 'ghost', text: 'Copy data',
      onclick: e => {
        const t = tsv(), btn = e.target;
        const done = () => { btn.textContent = 'Copied'; setTimeout(() => { btn.textContent = 'Copy data'; }, 1400); };
        if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(t).then(done, () => { });
        else { const ta = h('textarea', { }); ta.value = t; document.body.appendChild(ta); ta.select();
               try { document.execCommand('copy'); done(); } catch (err) { } ta.remove(); }
      }
    }) : null
  ]);
}

/* A CSS custom property as literal rgb, and a blend of two: the heatmap has to carry its own
 * colours out of the page, where neither var() nor color-mix() means anything. */
function cssRGB(name) {
  const raw = (getComputedStyle(document.documentElement).getPropertyValue(name) || '').trim();
  let m = raw.match(/^#([0-9a-f]{3}|[0-9a-f]{6})$/i);
  if (m) {
    const x = m[1].length === 3 ? m[1].split('').map(c => c + c).join('') : m[1];
    return [parseInt(x.slice(0, 2), 16), parseInt(x.slice(2, 4), 16), parseInt(x.slice(4, 6), 16)];
  }
  m = raw.match(/rgba?\(([^)]+)\)/);
  if (m) { const p = m[1].split(/[,\s/]+/).map(Number); return [p[0] | 0, p[1] | 0, p[2] | 0]; }
  return [136, 136, 136];
}
function mixRGB(a, b, t) {
  const k = Math.max(0, Math.min(1, t));
  return 'rgb(' + [0, 1, 2].map(i => Math.round(a[i] + (b[i] - a[i]) * k)).join(',') + ')';
}

/* ---------- the charts the selection can be drawn as ---------- */
function statsOf(M, row, linear) {
  return COND_COLS.map(cols => {
    const v = cols.map(j => linear ? unlg(lg(M, row, j)) : lg(M, row, j));
    const m = v.reduce((a, b) => a + b, 0) / (v.length || 1);
    const sd = v.length > 1
      ? Math.sqrt(v.reduce((a, b) => a + (b - m) * (b - m), 0) / (v.length - 1)) : 0;
    return { mean: m, sd: sd, v: v };
  });
}

/* feature x condition, mean with a standard-deviation whisker */
function chartConditionMeans(M, rows, ids, linear) {
  const H = 60 + rows.length * (18 * CONDS.length + 16), padL = 190, padR = 20, padT = 10;
  const pw = W - padL - padR;
  const st = rows.map(r => statsOf(M, r, linear));
  const max = Math.max(1e-9, ...st.flat().map(x => x.mean + x.sd));
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'condition means' });
  const sx = v => padL + pw * v / max;
  ticks(max, 4).forEach(v => {
    svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: padT, y2: H - 30, style: 'stroke:var(--grid)' }));
    svg.appendChild(s('text', { x: sx(v), y: H - 14, 'text-anchor': 'middle', class: 'tick' }, fmt(v)));
  });
  let y = padT + 6;
  rows.forEach((r, ri) => {
    svg.appendChild(s('text', { x: padL - 10, y: y + 11, 'text-anchor': 'end', class: 'axlab' },
      fit(ids[r], padL - 16)));
    CONDS.forEach((c, ci) => {
      const d = st[ri][ci], bh = 13;
      const bar = s('path', { d: barPath(padL, y, Math.max(0, sx(d.mean) - padL), bh, 3),
        style: `fill:${COND_COLOR[ci % 4]}` });
      svg.appendChild(bar);
      hoverable(bar, ids[r], [
        { label: c + ' mean', value: fmt(d.mean), color: COND_COLOR[ci % 4] },
        { label: 'SD', value: fmt(d.sd) },
        { label: 'libraries', value: d.v.map(x => fmt(x)).join(', ') }]);
      if (d.sd > 0) {
        svg.appendChild(s('line', { x1: sx(Math.max(0, d.mean - d.sd)), x2: sx(d.mean + d.sd),
          y1: y + bh / 2, y2: y + bh / 2, style: 'stroke:var(--ink-2);stroke-width:1.2' }));
        [Math.max(0, d.mean - d.sd), d.mean + d.sd].forEach(v =>
          svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: y + 2, y2: y + bh - 2,
            style: 'stroke:var(--ink-2);stroke-width:1.2' })));
      }
      y += 18;
    });
    y += 16;
  });
  return svg;
}

/* every library as its own dot, so n and spread are visible rather than summarised */
function chartStrip(M, rows, ids, linear) {
  const rowH = 46, H = 40 + rows.length * rowH, padL = 190, padR = 20, padT = 14;
  const pw = W - padL - padR;
  const all = [];
  rows.forEach(r => COND_COLS.forEach(cols => cols.forEach(j =>
    all.push(linear ? unlg(lg(M, r, j)) : lg(M, r, j)))));
  const max = Math.max(1e-9, ...all);
  const sx = v => padL + pw * v / max;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'per-library values' });
  ticks(max, 4).forEach(v => {
    svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: padT, y2: H - 26, style: 'stroke:var(--grid)' }));
    svg.appendChild(s('text', { x: sx(v), y: H - 10, 'text-anchor': 'middle', class: 'tick' }, fmt(v)));
  });
  rows.forEach((r, ri) => {
    const y0 = padT + ri * rowH;
    svg.appendChild(s('text', { x: padL - 10, y: y0 + 18, 'text-anchor': 'end', class: 'axlab' },
      fit(ids[r], padL - 16)));
    CONDS.forEach((c, ci) => {
      const yy = y0 + 8 + ci * 11;
      const vals = COND_COLS[ci].map(j => ({ v: linear ? unlg(lg(M, r, j)) : lg(M, r, j), j: j }));
      const m = vals.reduce((a, b) => a + b.v, 0) / (vals.length || 1);
      svg.appendChild(s('line', { x1: sx(m), x2: sx(m), y1: yy - 4, y2: yy + 4,
        style: `stroke:${COND_COLOR[ci % 4]};stroke-width:2` }));
      vals.forEach(o => {
        const c2 = s('circle', { cx: sx(o.v), cy: yy, r: 3,
          style: `fill:${COND_COLOR[ci % 4]};opacity:.72` });
        svg.appendChild(c2);
        hoverable(c2, ids[r], [{ label: EXPL.samples[o.j].name, value: fmt(o.v), color: COND_COLOR[ci % 4] },
                               { label: c + ' mean', value: fmt(m) }]);
      });
    });
  });
  return svg;
}

/* one bar per library, grouped and coloured by condition: the raw shape of the data */
function chartSamples(M, rows, ids, linear) {
  const panelH = 120, H = rows.length * (panelH + 26) + 26, padL = 60, padR = 14;
  const pw = W - padL - padR, bw = pw / NS;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'per-library bars' });
  const condOf = j => Math.max(0, CONDS.indexOf(EXPL.samples[j].condition || 'all'));
  rows.forEach((r, ri) => {
    const top = ri * (panelH + 26) + 18;
    const vals = [];
    for (let j = 0; j < NS; j++) vals.push(linear ? unlg(lg(M, r, j)) : lg(M, r, j));
    const max = Math.max(1e-9, ...vals);
    svg.appendChild(s('text', { x: padL, y: top - 5, class: 'axlab' }, fit(ids[r], pw)));
    ticks(max, 3).forEach(v => {
      const y = top + panelH - panelH * v / max;
      svg.appendChild(s('line', { x1: padL, x2: W - padR, y1: y, y2: y, style: 'stroke:var(--grid)' }));
      svg.appendChild(s('text', { x: padL - 6, y: y + 4, 'text-anchor': 'end', class: 'tick' }, fmt(v)));
    });
    vals.forEach((v, j) => {
      const hh = panelH * v / max, ci = condOf(j);
      const bar = s('rect', { x: padL + j * bw + 1.5, y: top + panelH - hh,
        width: Math.max(1, bw - 3), height: Math.max(0, hh), rx: 2,
        style: `fill:${COND_COLOR[ci % 4]}` });
      svg.appendChild(bar);
      hoverable(bar, EXPL.samples[j].name, [
        { label: ids[r], value: fmt(v), color: COND_COLOR[ci % 4] },
        { label: 'condition', value: EXPL.samples[j].condition || '—' }]);
    });
    if (ri === rows.length - 1) {
      for (let j = 0; j < NS; j++) {
        const t = s('text', { x: padL + j * bw + bw / 2, y: top + panelH + 14, 'text-anchor': 'end',
          class: 'tick', transform: `rotate(-38 ${padL + j * bw + bw / 2} ${top + panelH + 14})` },
          fit(EXPL.samples[j].name, 74));
        svg.appendChild(t);
      }
    }
  });
  return svg;
}

/* features x libraries. Row z-scores by default: a heatmap of absolute values just restates
 * which genes are abundant, which the table already says. */
function chartHeatmap(M, rows, ids, linear, zscore) {
  const cellH = Math.max(7, Math.min(16, Math.floor(620 / Math.max(1, rows.length))));
  const padL = 190, padR = 58, padT = 56;
  const H = padT + rows.length * cellH + 16, pw = W - padL - padR, cw = pw / NS;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'expression heatmap' });
  const mat = rows.map(r => {
    const v = []; for (let j = 0; j < NS; j++) v.push(linear ? unlg(lg(M, r, j)) : lg(M, r, j));
    if (!zscore) return v;
    const m = v.reduce((a, b) => a + b, 0) / NS;
    const sd = Math.sqrt(v.reduce((a, b) => a + (b - m) * (b - m), 0) / Math.max(1, NS - 1));
    return v.map(x => sd > 0 ? (x - m) / sd : 0);
  });
  const flat = mat.flat();
  const lim = zscore ? Math.max(0.5, ...flat.map(Math.abs)) : Math.max(1e-9, ...flat);
  /* Resolved to literal rgb() rather than left as color-mix(var(...)): an exported SVG is read
   * outside this document, where neither the custom properties nor color-mix would resolve. */
  const UPC = cssRGB('--up'), DNC = cssRGB('--down'), RMP = cssRGB('--ramp-1'), BGC = cssRGB('--plane');
  const colr = v => {
    if (zscore) {
      const t = Math.max(-1, Math.min(1, v / lim));
      return mixRGB(BGC, t >= 0 ? UPC : DNC, Math.abs(t));
    }
    return mixRGB(BGC, RMP, Math.max(0, Math.min(1, v / lim)));
  };
  for (let j = 0; j < NS; j++) {
    const x = padL + j * cw + cw / 2;
    svg.appendChild(s('text', { x: x, y: padT - 8, class: 'tick', 'text-anchor': 'start',
      transform: `rotate(-54 ${x} ${padT - 8})` }, fit(EXPL.samples[j].name, 70)));
    const ci = Math.max(0, CONDS.indexOf(EXPL.samples[j].condition || 'all'));
    svg.appendChild(s('rect', { x: padL + j * cw, y: padT - 5, width: cw - 1, height: 3,
      style: `fill:${COND_COLOR[ci % 4]}` }));
  }
  mat.forEach((v, ri) => {
    const y = padT + ri * cellH;
    if (cellH >= 9) svg.appendChild(s('text', { x: padL - 8, y: y + cellH - 2, 'text-anchor': 'end',
      class: 'tick' }, fit(ids[rows[ri]], padL - 14)));
    v.forEach((x, j) => {
      const cell = s('rect', { x: padL + j * cw, y: y, width: Math.max(1, cw - 1),
        height: Math.max(1, cellH - 1), style: `fill:${colr(x)}` });
      svg.appendChild(cell);
      hoverable(cell, ids[rows[ri]], [
        { label: EXPL.samples[j].name, value: zscore ? x.toFixed(2) + ' z' : fmt(x) },
        { label: 'condition', value: EXPL.samples[j].condition || '—' }]);
    });
  });
  /* the scale, so a colour can be read as a number */
  const gw = 12, gx = W - padR + 12, gy = padT, gh = Math.min(170, rows.length * cellH);
  for (let k = 0; k < 24; k++) {
    const t = 1 - k / 23, v = zscore ? (t * 2 - 1) * lim : t * lim;
    svg.appendChild(s('rect', { x: gx, y: gy + k * (gh / 24), width: gw, height: gh / 24 + 0.6,
      style: `fill:${colr(v)}` }));
  }
  svg.appendChild(s('text', { x: gx + gw + 4, y: gy + 8, class: 'tick' },
    zscore ? '+' + lim.toFixed(1) : fmt(lim)));
  svg.appendChild(s('text', { x: gx + gw + 4, y: gy + gh, class: 'tick' },
    zscore ? '-' + lim.toFixed(1) : '0'));
  return svg;
}

/* every feature as a line across the libraries: co-regulation, or the lack of it */
function chartProfile(M, rows, ids, linear) {
  const H = 330, padL = 58, padR = 16, padT = 12, padB = 62;
  const pw = W - padL - padR, ph = H - padT - padB;
  const series = rows.map(r => {
    const v = []; for (let j = 0; j < NS; j++) v.push(linear ? unlg(lg(M, r, j)) : lg(M, r, j));
    return v;
  });
  const max = Math.max(1e-9, ...series.flat());
  const sx = j => padL + (NS > 1 ? pw * j / (NS - 1) : pw / 2);
  const sy = v => padT + ph - ph * v / max;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'expression profiles' });
  ticks(max, 4).forEach(v => {
    svg.appendChild(s('line', { x1: padL, x2: W - padR, y1: sy(v), y2: sy(v), style: 'stroke:var(--grid)' }));
    svg.appendChild(s('text', { x: padL - 6, y: sy(v) + 4, 'text-anchor': 'end', class: 'tick' }, fmt(v)));
  });
  /* condition bands, so a step at a group boundary is visible */
  CONDS.forEach((c, ci) => {
    const cols = COND_COLS[ci]; if (!cols.length) return;
    const x0 = sx(Math.min(...cols)), x1 = sx(Math.max(...cols));
    svg.appendChild(s('rect', { x: x0 - 3, y: padT + ph + 6, width: Math.max(4, x1 - x0 + 6), height: 4,
      style: `fill:${COND_COLOR[ci % 4]}` }));
    svg.appendChild(s('text', { x: (x0 + x1) / 2, y: padT + ph + 22, 'text-anchor': 'middle',
      class: 'tick' }, c));
  });
  for (let j = 0; j < NS; j++)
    svg.appendChild(s('text', { x: sx(j), y: H - 8, 'text-anchor': 'end', class: 'tick',
      transform: `rotate(-38 ${sx(j)} ${H - 8})` }, fit(EXPL.samples[j].name, 70)));
  series.forEach((v, ri) => {
    const col = COND_COLOR[ri % 4];
    const d = v.map((x, j) => `${j ? 'L' : 'M'}${sx(j).toFixed(1)} ${sy(x).toFixed(1)}`).join('');
    const line = s('path', { d: d, style: `fill:none;stroke:${col};stroke-width:1.8;opacity:.85` });
    svg.appendChild(line);
    hoverable(line, ids[rows[ri]], v.map((x, j) =>
      ({ label: EXPL.samples[j].name, value: fmt(x) })).slice(0, 15));
    v.forEach((x, j) => svg.appendChild(s('circle', { cx: sx(j), cy: sy(x), r: 2.2,
      style: `fill:${col}` })));
  });
  return svg;
}

/* ---------- shrinkage: the same contrast with and without apeglm ---------- */
/* Both axes are fold changes of the same genes, so the frame is square and symmetric and the
 * identity line is drawn: shrinkage shows up as points pulled towards y = 0. */
function chartShrinkScatter(M, contrast, padjT, lfcT, pts) {
  const H = 380, pad = 54, padT = 12, padR = 14;
  const pw = W - pad - padR, ph = H - padT - 40;
  const lim = Math.max(1, ...pts.map(p => Math.max(Math.abs(p.raw), Math.abs(p.lfc)))) * 1.04;
  const sx = v => pad + pw * (v + lim) / (2 * lim);
  const sy = v => padT + ph - ph * (v + lim) / (2 * lim);
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img',
    'aria-label': 'shrunken against unshrunken log2 fold change' });
  const tk = ticks(lim, 3); const axv = tk.concat(tk.filter(v => v > 0).map(v => -v));
  axv.forEach(v => {
    svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: padT, y2: padT + ph, style: 'stroke:var(--grid)' }));
    svg.appendChild(s('line', { x1: pad, x2: W - padR, y1: sy(v), y2: sy(v), style: 'stroke:var(--grid)' }));
    svg.appendChild(s('text', { x: sx(v), y: H - 22, 'text-anchor': 'middle', class: 'tick' }, fmt(v)));
    svg.appendChild(s('text', { x: pad - 7, y: sy(v) + 4, 'text-anchor': 'end', class: 'tick' }, fmt(v)));
  });
  svg.appendChild(s('line', { x1: sx(-lim), y1: sy(-lim), x2: sx(lim), y2: sy(lim),
    style: 'stroke:var(--axis);stroke-width:1;stroke-dasharray:4 3' }));
  svg.appendChild(s('line', { x1: pad, x2: W - padR, y1: sy(0), y2: sy(0), style: 'stroke:var(--axis)' }));
  svg.appendChild(s('text', { x: pad + pw / 2, y: H - 5, 'text-anchor': 'middle', class: 'axlab' },
    'unshrunken log2 fold change (DESeq2 MLE)'));
  svg.appendChild(s('text', { x: 12, y: padT + ph / 2, class: 'axlab',
    transform: `rotate(-90 12 ${padT + ph / 2})`, 'text-anchor': 'middle' },
    'shrunken log2 fold change (reported)'));
  pts.forEach(p => {
    const sigS = p.padj < padjT && overCut(p.lfc, lfcT);
    const sigR = p.padj < padjT && overCut(p.raw, lfcT);
    const col = (sigS && sigR) ? S3 : sigR ? S2 : sigS ? S1 : MUTED;
    const c = s('circle', { cx: sx(p.raw), cy: sy(p.lfc), r: (sigS || sigR) ? 2.8 : 1.8,
      style: `fill:${col};opacity:${(sigS || sigR) ? .9 : .32}` });
    svg.appendChild(c);
    if (sigS || sigR) hoverable(c, p.id, [
      { label: 'shrunken log2FC', value: p.lfc.toFixed(2), color: S1 },
      { label: 'unshrunken log2FC', value: p.raw.toFixed(2), color: S2 },
      { label: 'padj', value: p.padj.toExponential(2) },
      { label: 'called by', value: (sigS && sigR) ? 'both' : sigR ? 'unshrunken only' : 'shrunken only' }]);
  });
  return svg;
}

function chartVolcano(pts, padjT, lfcT, useRaw, label) {
  const H = 360, pad = 54, padT = 12, padR = 14;
  const pw = W - pad - padR, ph = H - padT - 42;
  const xv = p => useRaw ? p.raw : p.lfc;
  const lim = Math.max(1, ...pts.map(p => Math.abs(xv(p)))) * 1.04;
  const ymax = Math.max(2, ...pts.map(p => p.mlp)) * 1.04;
  const sx = v => pad + pw * (v + lim) / (2 * lim);
  const sy = v => padT + ph - ph * v / ymax;
  const svg = s('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'volcano plot' });
  ticks(ymax, 4).forEach(v => {
    svg.appendChild(s('line', { x1: pad, x2: W - padR, y1: sy(v), y2: sy(v), style: 'stroke:var(--grid)' }));
    svg.appendChild(s('text', { x: pad - 7, y: sy(v) + 4, 'text-anchor': 'end', class: 'tick' }, fmt(v)));
  });
  const tk = ticks(lim, 3); tk.concat(tk.filter(v => v > 0).map(v => -v)).forEach(v => {
    svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: padT, y2: padT + ph, style: 'stroke:var(--grid)' }));
    svg.appendChild(s('text', { x: sx(v), y: H - 24, 'text-anchor': 'middle', class: 'tick' }, fmt(v)));
  });
  [-lfcT, lfcT].forEach(v => svg.appendChild(s('line', { x1: sx(v), x2: sx(v), y1: padT, y2: padT + ph,
    style: 'stroke:var(--axis);stroke-dasharray:4 3' })));
  const yT = -Math.log10(padjT);
  if (yT <= ymax) svg.appendChild(s('line', { x1: pad, x2: W - padR, y1: sy(yT), y2: sy(yT),
    style: 'stroke:var(--axis);stroke-dasharray:4 3' }));
  svg.appendChild(s('text', { x: pad + pw / 2, y: H - 6, 'text-anchor': 'middle', class: 'axlab' }, label));
  svg.appendChild(s('text', { x: 12, y: padT + ph / 2, class: 'axlab',
    transform: `rotate(-90 12 ${padT + ph / 2})`, 'text-anchor': 'middle' }, '-log10 adjusted p'));
  pts.forEach(p => {
    const x = xv(p), sig = p.padj < padjT && overCut(x, lfcT);
    const col = sig ? (x > 0 ? UP : DOWN) : MUTED;
    const c = s('circle', { cx: sx(x), cy: sy(p.mlp), r: sig ? 2.8 : 1.7,
      style: `fill:${col};opacity:${sig ? .88 : .3}` });
    svg.appendChild(c);
    if (sig) hoverable(c, p.id, [
      { label: 'log2FC ' + (useRaw ? '(unshrunken)' : '(shrunken)'), value: x.toFixed(2), color: col },
      { label: 'the other estimate', value: (useRaw ? p.lfc : p.raw).toFixed(2) },
      { label: 'padj', value: p.padj.toExponential(2) }]);
  });
  return svg;
}

/* ================================ the Explorer tab ================================ */
function sectionExplorer() {
  const out = h('div');
  if (!EXPL) {
    out.appendChild(h('p', { class: 'empty', text:
      'No expression matrices were staged for this run, so there is nothing to explore here.' }));
    return out;
  }
  if (typeof DecompressionStream === 'undefined') {
    out.appendChild(h('p', { class: 'caveat', text:
      'This browser has no DecompressionStream, which the page uses to unpack the expression ' +
      'matrices. Every other tab works; open this one in a current Firefox, Chrome or Safari.' }));
    return out;
  }
  out.appendChild(h('p', { class: 'sub', text:
    'The matrices themselves, not a summary of them. Pick a level and a gene map, filter down to ' +
    'the features you want, then plot or export exactly that selection. Gene-level values are ' +
    'DESeq2 median-of-ratios normalised counts; transcript-level values are salmon TPM, the ' +
    'quantity the expression filter was applied to. Fold changes are shown as reported ' +
    '(apeglm-shrunken) and as DESeq2 estimated them before shrinkage.' }));

  const X = {
    level: 'gene', map: (EXPL.maps[EXPL.primary] ? EXPL.primary : Object.keys(EXPL.maps)[0]),
    q: '', flagOnly: false, contrast: '', sigOnly: true, lfcMode: 'lfc',
    padj: EXPL.alpha, lfc: EXPL.lfc_cut, dir: 'all', sort: 'expr', limit: 100, linear: false,
    chart: 'cond', zscore: true, sel: new Set(), focus: null, M: null, descs: null
  };
  const descOf = r => {
    const ix = X.M.desc[r];
    return (ix >= 0 && X.descs && X.descs[ix]) ? X.descs[ix] : '';
  };
  const rowMean = r => {
    let t = 0; for (let j = 0; j < NS; j++) t += lg(X.M, r, j);
    return t / NS;
  };
  const contrastsOf = () => Object.keys((X.M && X.M.de) || {});

  /* ---- the three cards, filled in by draw() once the arrays are in ---- */
  const browseBox = h('div'), plotBox = h('div'), shrinkBox = h('div'), isoBox = h('div');
  const status = h('p', { class: 'empty', text: 'Unpacking the expression matrices…' });
  out.appendChild(status);
  out.appendChild(browseBox); out.appendChild(plotBox);
  out.appendChild(shrinkBox); out.appendChild(isoBox);

  /* ---------------- filtering ---------------- */
  function rowsNow() {
    const M = X.M, de = X.contrast ? M.de[X.contrast] : null;
    const lfcArr = de ? ((X.lfcMode === 'raw' && de.raw) ? de.raw : de.lfc) : null;
    const q = X.q.toLowerCase();
    const keep = [];
    for (let r = 0; r < M.n; r++) {
      if (X.flagOnly && !M.flag[r]) continue;
      if (q && M.ids[r].toLowerCase().indexOf(q) < 0 && descOf(r).toLowerCase().indexOf(q) < 0) continue;
      if (de) {
        const mlp = deP(de.mlp, r);
        if (mlp === null) continue;                /* not tested in this contrast */
        const L = deL(lfcArr, r);
        if (L === null) continue;
        if (X.dir === 'up' && !(L > 0)) continue;
        if (X.dir === 'down' && !(L < 0)) continue;
        if (X.sigOnly && !(Math.pow(10, -mlp) < X.padj && overCut(L, X.lfc))) continue;
      }
      keep.push(r);
    }
    const de2 = de, la = lfcArr;
    if (X.sort === 'id') keep.sort((a, b) => M.ids[a] < M.ids[b] ? -1 : M.ids[a] > M.ids[b] ? 1 : 0);
    else if (X.sort === 'padj' && de2) keep.sort((a, b) => (deP(de2.mlp, b) || 0) - (deP(de2.mlp, a) || 0));
    else if (X.sort === 'lfc' && de2) keep.sort((a, b) => Math.abs(deL(la, b) || 0) - Math.abs(deL(la, a) || 0));
    else { const mm = new Map(); keep.forEach(r => mm.set(r, rowMean(r))); keep.sort((a, b) => mm.get(b) - mm.get(a)); }
    return keep;
  }

  /* ---------------- card A: browse ---------------- */
  const tableHolder = h('div', { class: 'scroll' });
  const counter = h('span', { class: 'why' });
  const selCount = h('span', { class: 'why' });

  function numIn(val, step, min, max, onset, title) {
    return h('input', { type: 'number', value: String(val), step: String(step), min: String(min),
      max: max === null ? null : String(max), title: title || null,
      oninput: e => { const v = parseFloat(e.target.value); if (!Number.isNaN(v)) { onset(v); draw(); } } });
  }
  function sel(options, val, onset, aria) {
    const n = h('select', { 'aria-label': aria, onchange: e => { onset(e.target.value); draw(); } },
      options.map(([v, t]) => h('option', { value: v, text: t })));
    n.value = val;
    return n;
  }

  function controlsA() {
    const levelOpts = [['gene', 'Genes']];
    if (EXPL.tx) levelOpts.push(['tx', 'Transcripts']);
    const row1 = h('div', { class: 'controls' }, [
      h('label', { class: 'why', text: 'Level' }),
      sel(levelOpts, X.level, v => { if (v !== X.level) { X.level = v; X.sel.clear(); X.focus = null; X.M = null; } }, 'Level'),
      X.level === 'gene' && Object.keys(EXPL.maps).length > 1
        ? h('label', { class: 'why', text: 'Gene map' }) : null,
      X.level === 'gene' && Object.keys(EXPL.maps).length > 1
        ? sel(Object.keys(EXPL.maps).map(k => [k, EXPL.maps[k].label + (k === EXPL.primary ? ' (primary)' : '')]),
              X.map, v => { X.map = v; X.sel.clear(); X.focus = null; X.M = null; }, 'Gene map') : null,
      h('input', { type: 'search', placeholder: 'search id or annotation', value: X.q,
        oninput: e => { X.q = e.target.value; draw(); } }),
      (() => {
        const w = h('label', { class: 'why' });
        const cb = h('input', { type: 'checkbox', onchange: e => { X.flagOnly = e.target.checked; draw(); } });
        cb.checked = X.flagOnly;
        w.appendChild(cb); w.appendChild(document.createTextNode(' ' + (EXPL.interest_label || 'flagged') + ' only'));
        return w;
      })()
    ]);
    const cs = contrastsOf();
    const row2 = cs.length ? h('div', { class: 'controls' }, [
      h('label', { class: 'why', text: 'Contrast' }),
      sel([['', '— none (expression only)']].concat(cs.map(c => [c, c])), X.contrast,
          v => { X.contrast = v; if (X.sort !== 'expr' && !v) X.sort = 'expr'; }, 'Contrast'),
      X.contrast ? h('label', { class: 'why', text: 'Fold change' }) : null,
      X.contrast ? sel([['lfc', 'shrunken (reported)'], ['raw', 'unshrunken (MLE)']], X.lfcMode,
                       v => { X.lfcMode = v; }, 'Fold-change estimate') : null,
      X.contrast ? h('label', { class: 'why', text: 'padj <' }) : null,
      X.contrast ? numIn(X.padj, 0.01, 0.0001, 1, v => { X.padj = v; }) : null,
      X.contrast ? h('label', { class: 'why', text: '|log2FC| >' }) : null,
      X.contrast ? numIn(X.lfc, 0.5, 0, null, v => { X.lfc = v; }) : null,
      X.contrast ? sel([['all', 'up and down'], ['up', 'up only'], ['down', 'down only']], X.dir,
                       v => { X.dir = v; }, 'Direction') : null,
      X.contrast ? (() => {
        const w = h('label', { class: 'why' });
        const cb = h('input', { type: 'checkbox', onchange: e => { X.sigOnly = e.target.checked; draw(); } });
        cb.checked = X.sigOnly;
        w.appendChild(cb); w.appendChild(document.createTextNode(' significant only'));
        return w;
      })() : null
    ]) : null;
    const sortOpts = [['expr', 'mean expression'], ['id', 'id']];
    if (X.contrast) { sortOpts.push(['padj', 'adjusted p']); sortOpts.push(['lfc', '|log2FC|']); }
    const row3 = h('div', { class: 'controls' }, [
      h('label', { class: 'why', text: 'Sort by' }),
      sel(sortOpts, X.sort, v => { X.sort = v; }, 'Sort'),
      h('label', { class: 'why', text: 'Show' }),
      sel([['50', '50 rows'], ['100', '100 rows'], ['250', '250 rows'], ['500', '500 rows']],
          String(X.limit), v => { X.limit = parseInt(v, 10); }, 'Rows shown'),
      counter, h('span', { class: 'spacer' }), selCount
    ]);
    return [row1, row2, row3].filter(Boolean);
  }

  function drawTable(rows) {
    const shown = rows.slice(0, X.limit);
    counter.textContent = `${fmt(rows.length)} of ${fmt(X.M.n)} ${X.level === 'tx' ? 'transcripts' : 'genes'}`
      + (rows.length > X.limit ? ` · showing the first ${X.limit}` : '');
    selCount.textContent = `${X.sel.size} selected`;
    const de = X.contrast ? X.M.de[X.contrast] : null;
    const lfcArr = de ? ((X.lfcMode === 'raw' && de.raw) ? de.raw : de.lfc) : null;
    const t = h('table');
    const hr = h('tr');
    ['', X.level === 'tx' ? 'Transcript' : 'Gene', 'Annotation'].forEach(l => hr.appendChild(h('th', { text: l })));
    if (X.level === 'gene') hr.appendChild(h('th', { class: 'num', text: 'tx' }));
    CONDS.forEach(c => hr.appendChild(h('th', { class: 'num', text: c })));
    if (de) ['log2FC', 'other est.', 'padj'].forEach(l => hr.appendChild(h('th', { class: 'num', text: l })));
    t.appendChild(h('thead', {}, [hr]));
    const tb = h('tbody');
    shown.forEach(r => {
      const tr = h('tr');
      const cb = h('input', { type: 'checkbox', 'aria-label': 'select ' + X.M.ids[r],
        onchange: e => {
          if (e.target.checked) X.sel.add(r); else X.sel.delete(r);
          selCount.textContent = `${X.sel.size} selected`;
          drawPlot(); drawIso();
        } });
      cb.checked = X.sel.has(r);
      tr.appendChild(h('td', {}, [cb]));
      const idCell = h('td', {}, [h('button', { class: 'linkish', text: X.M.ids[r],
        title: 'focus this feature', onclick: () => { X.focus = r; drawIso(); } })]);
      tr.appendChild(idCell);
      tr.appendChild(h('td', { text: descOf(r) || '—' }));
      if (X.level === 'gene') tr.appendChild(h('td', { class: 'num', text: X.M.ntx ? fmt(X.M.ntx[r]) : '—' }));
      COND_COLS.forEach(cols => {
        const m = cols.reduce((a, j) => a + lg(X.M, r, j), 0) / (cols.length || 1);
        tr.appendChild(h('td', { class: 'num', text: X.linear ? fmt(unlg(m)) : m.toFixed(2) }));
      });
      if (de) {
        const L = deL(lfcArr, r), O = deL((X.lfcMode === 'raw' && de.raw) ? de.lfc : de.raw, r);
        const mlp = deP(de.mlp, r);
        tr.appendChild(h('td', { class: 'num', text: L === null ? '—' : L.toFixed(2) }));
        tr.appendChild(h('td', { class: 'num', text: O === null ? '—' : O.toFixed(2) }));
        tr.appendChild(h('td', { class: 'num',
          text: mlp === null ? '—' : Math.pow(10, -mlp).toExponential(1) }));
      }
      tb.appendChild(tr);
    });
    t.appendChild(tb);
    tableHolder.textContent = '';
    tableHolder.appendChild(t);
    return shown;
  }

  /* ---------------- card B: plot whatever is selected ---------------- */
  const CHARTS = [
    ['cond', 'Condition means ± SD', 12],
    ['strip', 'Every library as a dot', 12],
    ['samples', 'Bars per library', 9],
    ['profile', 'Profiles across libraries', 25],
    ['heat', 'Heatmap', 200]
  ];
  const plotHolder = h('div');

  function selectionTSV() {
    const rows = [...X.sel];
    const head = [X.level === 'tx' ? 'transcript' : 'gene', 'annotation']
      .concat(EXPL.samples.map(s => s.name));
    const lines = [head.join('\t')];
    rows.forEach(r => {
      const v = [];
      for (let j = 0; j < NS; j++) { const q = lg(X.M, r, j); v.push((X.linear ? unlg(q) : q).toFixed(4)); }
      lines.push([X.M.ids[r], descOf(r).replace(/\t/g, ' ')].concat(v).join('\t'));
    });
    lines.push('');
    lines.push('# values: ' + (X.linear ? X.M.unit : 'log2(' + X.M.unit + ' + 1)'));
    lines.push('# conditions: ' + EXPL.samples.map(s => s.name + '=' + (s.condition || '?')).join(' '));
    lines.push('# stored to ~0.04% of value; for exact numbers use the run\'s '
      + (X.level === 'tx' ? 'quantification/<sample>/quant.sf' : 'normalized_counts.csv'));
    return lines.join('\n');
  }

  function drawPlot() {
    const spec = CHARTS.find(c => c[0] === X.chart) || CHARTS[0];
    const cap = spec[2];
    const rows = [...X.sel].slice(0, cap);
    const over = X.sel.size - rows.length;
    plotHolder.textContent = '';
    if (!X.sel.size) {
      plotHolder.appendChild(h('p', { class: 'empty', text:
        'Nothing selected yet. Tick rows in the table above — or filter down to what you want ' +
        'and use "select all shown" — and the plot is drawn from exactly that set.' }));
      return;
    }
    const chartBox = h('div');
    let svg;
    if (X.chart === 'cond') svg = chartConditionMeans(X.M, rows, X.M.ids, X.linear);
    else if (X.chart === 'strip') svg = chartStrip(X.M, rows, X.M.ids, X.linear);
    else if (X.chart === 'samples') svg = chartSamples(X.M, rows, X.M.ids, X.linear);
    else if (X.chart === 'profile') svg = chartProfile(X.M, rows, X.M.ids, X.linear);
    else svg = chartHeatmap(X.M, rows, X.M.ids, X.linear, X.zscore);
    chartBox.appendChild(svg);
    plotHolder.appendChild(chartBox);
    plotHolder.appendChild(legendRow(CONDS.map((c, i) => ({ label: c, color: COND_COLOR[i % 4] }))));
    if (over > 0) plotHolder.appendChild(h('p', { class: 'caveat', text:
      `${fmt(X.sel.size)} features are selected and this chart draws the first ${cap}; ` +
      `${fmt(over)} are not shown. The heatmap takes up to 200, and the TSV export takes all of them.` }));
    plotHolder.appendChild(exportRow(chartBox,
      `${X.level === 'tx' ? 'transcripts' : 'genes'}_${X.chart}`, selectionTSV));
  }

  /* ---------------- card C: with and without shrinkage ---------------- */
  const shrinkHolder = h('div');
  let shrinkView = 'scatter';

  function shrinkPoints(contrast) {
    const de = X.M.de[contrast];
    if (!de || !de.raw) return null;
    const sig = [], bg = [];
    for (let r = 0; r < X.M.n; r++) {
      const mlp = deP(de.mlp, r); if (mlp === null) continue;
      const L = deL(de.lfc, r), R = deL(de.raw, r);
      if (L === null || R === null) continue;
      const padj = Math.pow(10, -mlp);
      const o = { id: X.M.ids[r], lfc: L, raw: R, padj: padj, mlp: mlp };
      /* Everything that either estimate calls significant is kept in full, so the live counts
       * are exact; only the indifferent cloud is thinned, and it is never counted. */
      ((padj < X.padj && (overCut(L, X.lfc) || overCut(R, X.lfc))) ? sig : bg).push(o);
    }
    const step = Math.max(1, Math.floor(bg.length / 3500));
    return { all: sig.concat(bg.filter((_, i) => i % step === 0)), sig: sig, nbg: bg.length,
             shown_bg: bg.filter((_, i) => i % step === 0).length };
  }

  function drawShrink() {
    shrinkHolder.textContent = '';
    const cs = contrastsOf();
    if (!cs.length) return;
    const c = X.contrast || cs[0];
    const de = X.M.de[c];
    if (!de) return;
    if (!de.raw) {
      shrinkHolder.appendChild(h('p', { class: 'empty', text:
        'This run’s DE tables carry no unshrunken fold change, so there is nothing to compare.' }));
      return;
    }
    /* exact counts over every tested gene, not over the thinned scatter */
    let nS = 0, nR = 0, both = 0, onlyS = 0, onlyR = 0, tested = 0;
    for (let r = 0; r < X.M.n; r++) {
      const mlp = deP(de.mlp, r); if (mlp === null) continue;
      tested++;
      const padj = Math.pow(10, -mlp), L = deL(de.lfc, r), R = deL(de.raw, r);
      if (L === null || R === null) continue;
      const a = padj < X.padj && overCut(L, X.lfc), b = padj < X.padj && overCut(R, X.lfc);
      if (a) nS++; if (b) nR++;
      if (a && b) both++; else if (a) onlyS++; else if (b) onlyR++;
    }
    const pts = shrinkPoints(c);
    shrinkHolder.appendChild(h('div', { class: 'controls' }, [
      h('label', { class: 'why', text: 'Contrast' }),
      sel(cs.map(x => [x, x]), c, v => { X.contrast = v; }, 'Contrast'),
      h('label', { class: 'why', text: 'padj <' }),
      numIn(X.padj, 0.01, 0.0001, 1, v => { X.padj = v; }),
      h('label', { class: 'why', text: '|log2FC| >' }),
      numIn(X.lfc, 0.5, 0, null, v => { X.lfc = v; }),
      h('span', { class: 'spacer' }),
      sel([['scatter', 'one against the other'], ['vshrunk', 'volcano, shrunken'],
           ['vraw', 'volcano, unshrunken']], shrinkView, v => { shrinkView = v; }, 'View')
    ]));
    shrinkHolder.appendChild(buildTable({
      cols: [{ key: 'k', label: 'At padj < ' + X.padj + ' and |log2FC| > ' + X.lfc },
             { key: 'v', label: 'Genes', num: true }, { key: 'n', label: 'Reads as' }],
      rows: [
        { k: 'Tested in this contrast', v: fmt(tested), n: 'genes DESeq2 returned a p-value for' },
        { k: 'Significant, shrunken log2FC (reported)', v: fmt(nS), n: 'what every other tab counts' },
        { k: 'Significant, unshrunken log2FC', v: fmt(nR), n: 'the MLE, before apeglm' },
        { k: 'Called by both', v: fmt(both), n: 'robust to the choice' },
        { k: 'Only with shrinkage', v: fmt(onlyS), n: 'shrinkage moved it past the cut — rare' },
        { k: 'Only without shrinkage', v: fmt(onlyR), n: 'a large but badly-estimated fold change' }
      ]
    }));
    const box = h('div');
    if (shrinkView === 'scatter') {
      box.appendChild(chartShrinkScatter(X.M, c, X.padj, X.lfc, pts.all));
      box.appendChild(legendRow([
        { label: 'called by both', color: S3 }, { label: 'unshrunken only', color: S2 },
        { label: 'shrunken only', color: S1 }, { label: 'neither', color: MUTED }]));
    } else {
      const useRaw = shrinkView === 'vraw';
      box.appendChild(chartVolcano(pts.all, X.padj, X.lfc, useRaw,
        (useRaw ? 'unshrunken' : 'shrunken') + ' log2 fold change — ' + c));
      box.appendChild(legendRow([{ label: 'up', color: UP }, { label: 'down', color: DOWN },
                                 { label: 'not significant', color: MUTED }]));
    }
    shrinkHolder.appendChild(box);
    shrinkHolder.appendChild(exportRow(box, 'shrinkage_' + c, null));
    shrinkHolder.appendChild(h('p', { class: 'note-line', text:
      `Every gene either estimate calls significant is plotted; of the ${fmt(pts.nbg)} indifferent ` +
      `genes, ${fmt(pts.shown_bg)} are drawn as context. The counts in the table are exact over all ` +
      `${fmt(tested)} tested genes.` }));
    shrinkHolder.appendChild(h('p', { class: 'note-line', text:
      'Why the shrunken estimate is the one reported: a log2 fold change of 6 from four reads ' +
      'against one is not evidence of a six-fold change, and apeglm pulls exactly those towards ' +
      'zero while leaving well-measured genes alone. The adjusted p-value is the same in both ' +
      'columns — shrinkage changes the effect size, not the test — so a gene that appears ' +
      'only in the unshrunken count is one whose significance never depended on the |log2FC| cut ' +
      'being met by a reliable estimate.' }));
  }

  /* ---------------- card D: the transcripts behind one gene ---------------- */
  const isoHolder = h('div');
  let TXL = null, txWanted = false;

  function drawIso() {
    isoHolder.textContent = '';
    if (X.level !== 'gene') {
      isoHolder.appendChild(h('p', { class: 'empty', text:
        'Showing transcripts already — switch the level to Genes to group them.' }));
      return;
    }
    if (!EXPL.tx || !X.M.txGene) {
      isoHolder.appendChild(h('p', { class: 'empty', text:
        'No transcript-level quantification was staged for this run.' }));
      return;
    }
    const focus = X.focus !== null ? X.focus : (X.sel.size ? [...X.sel][0] : null);
    if (focus === null) {
      isoHolder.appendChild(h('p', { class: 'empty', text:
        'Click a gene id in the table above and its transcripts are broken out here.' }));
      return;
    }
    if (!TXL) {
      isoHolder.appendChild(h('p', { class: 'empty', text: 'Unpacking the transcript matrix…' }));
      if (!txWanted) {
        txWanted = true;
        loadLevel('tx', null).then(t => { TXL = t; drawIso(); },
          e => { isoHolder.textContent = ''; isoHolder.appendChild(
            h('p', { class: 'empty', text: 'Transcript matrix could not be read: ' + e.message })); });
      }
      return;
    }
    const tg = X.M.txGene, rows = [];
    for (let j = 0; j < tg.length; j++) if (tg[j] === focus) rows.push(j);
    const gid = X.M.ids[focus];
    isoHolder.appendChild(h('div', { class: 'controls' }, [
      h('span', { class: 'why', text: `${gid} — ${fmt(rows.length)} transcript` +
        (rows.length === 1 ? '' : 's') + ' · ' + (descOf(focus) || 'no annotation') }),
      h('span', { class: 'spacer' }),
      h('button', { class: 'ghost', text: 'Clear', onclick: () => { X.focus = null; drawIso(); } })
    ]));
    if (!rows.length) {
      isoHolder.appendChild(h('p', { class: 'empty', text:
        'No transcript of this gene is in the quantification — which should not happen; ' +
        'the gene was built from transcripts.' }));
      return;
    }
    const box = h('div');
    box.appendChild(rows.length === 1
      ? chartSamples(TXL, rows, TXL.ids, X.linear)
      : chartConditionMeans(TXL, rows.slice(0, 12), TXL.ids, X.linear));
    isoHolder.appendChild(box);
    isoHolder.appendChild(legendRow(CONDS.map((c, i) => ({ label: c, color: COND_COLOR[i % 4] }))));
    isoHolder.appendChild(buildTable({
      cols: [{ key: 't', label: 'Transcript' }, { key: 'a', label: 'Annotation' }]
        .concat(CONDS.map(c => ({ key: 'c_' + c, label: c + ' (TPM)', num: true })))
        .concat([{ key: 'sh', label: '% of gene TPM', num: true }]),
      rows: (() => {
        const tot = CONDS.map((_, ci) => rows.reduce((a, r) =>
          a + COND_COLS[ci].reduce((b, j) => b + unlg(lg(TXL, r, j)), 0) / (COND_COLS[ci].length || 1), 0));
        const totAll = tot.reduce((a, b) => a + b, 0) || 1;
        return rows.slice(0, 40).map(r => {
          const o = { t: TXL.ids[r] };
          const ix = TXL.desc[r];
          o.a = (ix >= 0 && X.descs && X.descs[ix]) ? X.descs[ix] : '—';
          let mine = 0;
          CONDS.forEach((c, ci) => {
            const m = COND_COLS[ci].reduce((b, j) => b + unlg(lg(TXL, r, j)), 0) / (COND_COLS[ci].length || 1);
            o['c_' + c] = fmt(m); mine += m;
          });
          o.sh = (100 * mine / totAll).toFixed(1) + '%';
          return o;
        });
      })()
    }));
    isoHolder.appendChild(exportRow(box, 'isoforms_' + gid, () => {
      const head = ['transcript', 'annotation'].concat(EXPL.samples.map(s => s.name));
      const lines = [head.join('\t')];
      rows.forEach(r => {
        const v = []; for (let j = 0; j < NS; j++) v.push(unlg(lg(TXL, r, j)).toFixed(4));
        const ix = TXL.desc[r];
        lines.push([TXL.ids[r], (ix >= 0 && X.descs ? X.descs[ix] : '')].concat(v).join('\t'));
      });
      return lines.join('\n');
    }));
    isoHolder.appendChild(h('p', { class: 'note-line', text:
      'TPM, so the shares are comparable within a library. A gene whose differential expression ' +
      'comes from one transcript while the others hold still is a candidate for isoform switching ' +
      'rather than for a change in gene output.' }));
  }

  /* ---------------- assembly ---------------- */
  function draw() {
    const rows = rowsNow();
    drawTable(rows);
    browseBox.textContent = '';
    const cA = h('div', { class: 'card' });
    cA.appendChild(h('div', { class: 'figcap' }, [
      h('h3', { text: 'Browse the expression matrix' }),
      h('span', { class: 'why', text: X.level === 'tx' ? 'salmon TPM per transcript'
                                                       : 'normalised counts per gene, ' + EXPL.maps[X.map].label }),
      h('span', { class: 'spacer' }),
      (() => {
        const w = h('label', { class: 'why' });
        const cb = h('input', { type: 'checkbox',
          onchange: e => { X.linear = e.target.checked; draw(); } });
        cb.checked = X.linear;
        w.appendChild(cb);
        w.appendChild(document.createTextNode(' linear scale (otherwise log2)'));
        return w;
      })()
    ]));
    controlsA().forEach(n => cA.appendChild(n));
    cA.appendChild(h('div', { class: 'controls' }, [
      h('button', { class: 'ghost', text: `Select all shown (${Math.min(rows.length, X.limit)})`,
        onclick: () => { rows.slice(0, X.limit).forEach(r => X.sel.add(r)); draw(); } }),
      h('button', { class: 'ghost', text: 'Select the top 50 of this filter',
        onclick: () => { rows.slice(0, 50).forEach(r => X.sel.add(r)); draw(); } }),
      h('button', { class: 'ghost', text: 'Clear selection',
        onclick: () => { X.sel.clear(); draw(); } })
    ]));
    cA.appendChild(tableHolder);
    cA.appendChild(h('p', { class: 'note-line', text:
      'Condition columns are the mean over that condition’s libraries. "other est." is the same ' +
      'gene’s fold change under the estimate not currently chosen, so the two can be read ' +
      'against each other without switching. Rows are ranked by the sort, and the number shown is ' +
      'capped for the browser’s sake — the filter always reports the full count.' }));
    cA.appendChild(h('p', { class: 'note-line', text:
      'Precision: expression is stored to about 0.04% of its value, which is far below the ' +
      'measurement but means the last digit or two here can differ from normalized_counts.csv / ' +
      'quant.sf — use those files when an exact count matters. Fold changes and adjusted ' +
      'p-values are stored finely enough to reproduce DE_summary.csv exactly, so the significance ' +
      'counts on this tab are not approximations.' }));
    browseBox.appendChild(cA);

    plotBox.textContent = '';
    const cB = h('div', { class: 'card' });
    cB.appendChild(h('div', { class: 'figcap' }, [
      h('h3', { text: 'Plot the selection' }),
      h('span', { class: 'why', text: 'drawn from exactly the rows ticked above' })
    ]));
    cB.appendChild(h('div', { class: 'controls' }, [
      h('label', { class: 'why', text: 'Chart' }),
      sel(CHARTS.map(c => [c[0], c[1]]), X.chart, v => { X.chart = v; }, 'Chart type'),
      X.chart === 'heat' ? (() => {
        const w = h('label', { class: 'why' });
        const cb = h('input', { type: 'checkbox', onchange: e => { X.zscore = e.target.checked; draw(); } });
        cb.checked = X.zscore;
        w.appendChild(cb); w.appendChild(document.createTextNode(' row z-score'));
        return w;
      })() : null
    ]));
    cB.appendChild(plotHolder);
    plotBox.appendChild(cB);
    drawPlot();

    shrinkBox.textContent = '';
    if (contrastsOf().length) {
      const cC = h('div', { class: 'card' });
      cC.appendChild(h('div', { class: 'figcap' }, [
        h('h3', { text: 'Fold changes with and without shrinkage' }),
        h('span', { class: 'why', text: 'apeglm-shrunken against the raw DESeq2 estimate' })
      ]));
      cC.appendChild(shrinkHolder);
      shrinkBox.appendChild(cC);
      drawShrink();
    }

    isoBox.textContent = '';
    const cD = h('div', { class: 'card' });
    cD.appendChild(h('div', { class: 'figcap' }, [
      h('h3', { text: 'The transcripts behind one gene' }),
      h('span', { class: 'why', text: 'which isoform carries the gene-level signal' })
    ]));
    cD.appendChild(isoHolder);
    isoBox.appendChild(cD);
    drawIso();
  }

  /* arrays first, then the whole tab; a failure says so rather than rendering half a page */
  (async () => {
    try {
      X.descs = await descTable();
      X.M = await loadLevel(X.level, X.map);
      const cs = contrastsOf();
      if (!X.contrast && cs.length) X.contrast = cs[0];
      status.remove();
      draw();
    } catch (e) {
      status.className = 'caveat';
      status.textContent = 'The expression matrices could not be unpacked: ' + e.message;
    }
  })();

  /* the level / map switches drop X.M, so every redraw after one has to reload first */
  const origDraw = draw;
  draw = function () {
    if (!X.M) {
      status.className = 'empty';
      status.textContent = 'Unpacking…';
      browseBox.before(status);
      loadLevel(X.level, X.map).then(m => {
        X.M = m;
        const cs = contrastsOf();
        if (X.contrast && !cs.includes(X.contrast)) X.contrast = cs.length ? cs[0] : '';
        if (!X.contrast && cs.length) X.contrast = cs[0];
        status.remove();
        origDraw();
      });
      return;
    }
    origDraw();
  };
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
      /* The reader arrives here from a headline gene count and finds a smaller "tested" number,
         so the card has to account for the whole gap, not just the part DESeq2 is responsible
         for. Three steps, and the first happens BEFORE DESeq2 is called: a presence pre-filter on
         raw counts, which is why genes-in-model is already below the map's gene count. */
      /* the ACTIVE map's gene count, not the primary one: this card follows the map switcher,
         and a primary-map figure here would print a negative pre-filter on the other two */
      const mapGenes = ((DATA.gene_maps || []).find(m => m.key === activeMapKey()) || {}).genes;
      const inModel = FILT[0] && FILT[0].in_model;
      const preDrop = (mapGenes && inModel) ? mapGenes - inModel : null;
      const chain = preDrop === null ? '' :
        `Chain from the headline gene count: ${fmt(mapGenes)} genes in the ${activeMapLabel()} map ` +
        `→ ${fmt(preDrop)} removed by the presence pre-filter, which is applied before DESeq2 is ` +
        `called (a gene needs the minimum count in the minimum number of libraries) → ` +
        `${fmt(inModel)} genes in the model, the figure all contrasts share → the per-contrast ` +
        `Tested column. Only the last step is DESeq2's own. `;
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
        note: chain + "Independent filtering trims low-count genes at a per-contrast baseMean threshold, so it varies between contrasts. Cook's outliers leave the test entirely — with 5 replicates DESeq2 cannot replace them, so they carry pvalue = NA — and because the cutoff is a property of the design rather than of a contrast, the same genes are dropped from every contrast."
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
                 { key: 'seqs', label: 'Gene FASTA', num: true },
                 { key: 'note', label: 'Caveat' }],
          rows: (DATA.gene_maps || []).map(m => ({ name: m.name + (m.primary ? '  (primary)' : ''), basis: m.basis, genes: fmt(m.genes),
                 sig: (m.sig === null || m.sig === undefined) ? '—' : fmt(m.sig),
                 seqs: m.seqs ? fmt(m.seqs.genes) + ' nt · ' + fmt(m.seqs.pep) + ' aa' : '—',
                 note: m.note }))
        },
        note: 'Every map runs through the same quantifications and the same DESeq2 code, so a difference in DE genes is a difference in gene identity and nothing else. The headline numbers come from the map named by dge_primary_map; the others are published beside them as sensitivity analyses. Each map also ships <code>&lt;set&gt;.genes.fasta</code> and <code>.genes.pep</code> in transcriptome/ — one representative per gene (its longest transcript, and the longest ORF on any of its transcripts), so a gene set can be BLASTed or scored in protein mode; the aa count is how many genes have an ORF at all.'
      }));
    }

    /* Why one map is primary. The card above compares the maps on what the run itself can see
       (gene counts, DE counts) -- which cannot settle the choice, because a map that calls more DE
       genes is not thereby better. This one shows the two things that can settle it: agreement
       with INDEPENDENT truth, and how many DE calls each map invents when the condition labels are
       shuffled. Rendered only when an external evaluation is staged. */
    if ((DATA.genemap_eval || []).length) {
      const EV = DATA.genemap_eval, NUL = DATA.genemap_null || {};
      const truths = [...new Set(EV.map(r => r.truth))];
      const maps   = [...new Set(EV.map(r => r.map))];
      const get = (m, t) => EV.find(r => r.map === m && r.truth === t) || null;
      const pm = (DATA.meta || {}).primary_map || '';
      // the primary map's own row, so every claim below is anchored to a number in this table
      const best = t => {
        const vals = maps.map(m => (get(m, t) || {}).f1).filter(v => v !== null && v !== undefined);
        return vals.length ? Math.max(...vals) : null;
      };
      const nulls = maps.map(m => (NUL[m] || {}).false_calls).filter(v => v !== null && v !== undefined);
      const bestNull = nulls.length ? Math.min(...nulls) : null;
      const isPrimary = m => m === pm || m.toLowerCase() === pm.toLowerCase();

      const cols = [{ key: 'map', label: 'Map' }];
      truths.forEach(t => cols.push({ key: 'f1_' + t, label: 'Matches ' + t + ' truth (F1)', num: true }));
      cols.push({ key: 'split', label: 'Split (under-merged) %', num: true },
                { key: 'fused', label: 'Fused (over-merged) %', num: true });
      if (nulls.length) cols.push({ key: 'nullcalls', label: 'False DE calls / shuffle', num: true });

      const rows = maps.map(m => {
        const r = { map: m + (isPrimary(m) ? '  (primary)' : '') };
        truths.forEach(t => {
          const e = get(m, t);
          const v = e && e.f1 !== null && e.f1 !== undefined ? e.f1 : null;
          r['f1_' + t] = v === null ? '\u2014'
                       : v.toFixed(3) + (best(t) !== null && Math.abs(v - best(t)) < 1e-9 ? '  \u2190 best' : '');
        });
        const any = get(m, truths[0]);
        r.split = any && any.split !== null && any.split !== undefined ? any.split.toFixed(1) + '%' : '\u2014';
        r.fused = any && any.fused !== null && any.fused !== undefined ? any.fused.toFixed(1) + '%' : '\u2014';
        if (nulls.length) {
          const n = (NUL[m] || {}).false_calls;
          r.nullcalls = (n === null || n === undefined) ? '\u2014'
                      : n.toFixed(1) + (bestNull !== null && Math.abs(n - bestNull) < 1e-9 ? '  \u2190 fewest' : '');
        }
        return r;
      });

      const chartTruth = truths[0];
      out.appendChild(card({
        title: 'Why ' + (pm || 'the primary map') + ' is the primary map',
        why: 'decided on independent truth and a label-permutation null \u2014 not on which map calls more genes',
        chart: barH({
          rows: maps.map((m, idx) => {
            const e = get(m, chartTruth);
            return { label: m + (isPrimary(m) ? ' \u2605' : ''),
                     value: e && e.f1 ? e.f1 : 0,
                     color: isPrimary(m) ? S1 : MUTED,
                     tip: truths.map(t => {
                       const x = get(m, t);
                       return { label: 'F1 vs ' + t, value: x && x.f1 ? x.f1.toFixed(3) : '\u2014' };
                     }).concat((NUL[m] ? [{ label: 'false calls / shuffle', value: NUL[m].false_calls.toFixed(1) }] : [])) };
          }), gutter: 190, measure: 'F1 vs ' + chartTruth + ' truth'
        }),
        table: { cols: cols, rows: rows },
        note: (DATA.genemap_verdict || '') +
          ' How to read it. F1 is agreement with a grouping the pipeline never saw, so it is the ' +
          'only column here that speaks to whether a map\u2019s genes ARE genes; 1.000 would be ' +
          'perfect agreement. Split% and Fused% are the two ways a map can be wrong and they pull ' +
          'against each other \u2014 splitting one gene across several ids (under-merging) versus ' +
          'fusing distinct genes into one id (over-merging) \u2014 so no map can minimise both and ' +
          'the right map is the one with the better total error, not a zero in either column. ' +
          'False DE calls per shuffle comes from re-running the whole DE test on permuted condition ' +
          'labels: with the labels destroyed every call is false, so a map that reports more of them ' +
          'is manufacturing signal. A map built using the condition labels will score badly here ' +
          'however good its F1, because its gene definition has already seen what DESeq2 then tests.'
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
  ['explorer', 'Explorer', sectionExplorer],
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
        .replace("__SCRIPT__", SCRIPT + SCRIPT2 + SCRIPT4 + SCRIPT3))

os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
with open(args.out, "w") as fh:
    fh.write(page)

print(f"[dashboard] {args.out}  ({len(page) / 1024:.0f} KB)")
print(f"[dashboard] stages {len(stages)} · samples {len(samples)} · reference variants {len(ladder)}"
      f" · contrasts {len(dge_primary['contrasts']) if dge_primary else 0}"
      f" · primary map {primary_map} · maps with DE {sorted(k for k, v in dge_by_map.items() if v)}"
      f" · enrichment {sorted(enrichment_by_map)} · figures {sorted(figures_by_map)}")
