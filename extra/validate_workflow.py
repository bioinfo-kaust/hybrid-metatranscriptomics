#!/usr/bin/env python3
"""Static checks for the failure modes this pipeline has actually hit.

A stub run proves the DAG connects and the scripts parse. It does not prove that a stub declares
every output it promises, that a referenced parameter exists, or that a tool named in params.mod
is loadable. Each check below exists because the corresponding mistake reached the cluster:

  stub-outputs   CORSET gained an output whose stub never created it. CORSET is error_ignore, so
                 the task failed silently, its channels emitted nothing, and two downstream
                 processes vanished from a run that still reported success.
  params         a params.X that is never defined resolves to null and is interpolated into a
                 command line as the empty string — a silently different command, not an error.
  scripts        bin/ scripts are staged and run per task; a syntax error surfaces only when that
                 task runs, which can be days into a run.
  modules        an environment module that does not exist fails the task at load time.

Usage:  extra/validate_workflow.py [--pipeline DIR]
Exit code is non-zero if any check fails, so it can gate a submission.
"""
import argparse
import ast
import glob
import os
import re
import shutil
import subprocess
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--pipeline", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
args = ap.parse_args()
ROOT = args.pipeline

FAIL, WARN = [], []


def report(ok, check, detail):
    print(f"  [{'ok' if ok else 'FAIL'}] {check}: {detail}")


# ---------------------------------------------------------------- 1. stub completeness
def process_blocks(text):
    """Yield (name, body) for each `process NAME { ... }`, brace-matched."""
    for m in re.finditer(r"^process\s+(\w+)\s*\{", text, re.M):
        i, depth = m.end() - 1, 0
        for j in range(i, len(text)):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    yield m.group(1), text[i:j]
                    break


def expand_braces(name):
    """foo_R{1,2}.fq -> [foo_R{1,2}.fq, foo_R1.fq, foo_R2.fq]; a stub may write either form."""
    out = [name]
    m = re.search(r"\{([^{}]*,[^{}]*)\}", name)
    if m:
        for alt in m.group(1).split(","):
            out.extend(expand_braces(name[:m.start()] + alt + name[m.end():]))
    return out


def declared_outputs(body):
    """Output names a process promises. Optional outputs are excluded: a stub is not required to
    produce what the real task may legitimately omit."""
    m = re.search(r"\n\s*output:(.*?)\n\s*(script|shell|exec|stub):", body, re.S)
    if not m:
        return []
    names = []
    for line in m.group(1).splitlines():
        if "optional" in line and "true" in line:
            continue
        for raw in re.findall(r"""path[\s(]+["']([^"']+)["']""", line):
            names.append(raw)
    return names


def stub_body(body):
    m = re.search(r"\n\s*stub:(.*)$", body, re.S)
    return m.group(1) if m else None


print("\n[1] every declared output is created by the stub")
for nf in sorted(glob.glob(os.path.join(ROOT, "modules", "local", "*.nf"))):
    text = open(nf).read()
    for name, body in process_blocks(text):
        stub = stub_body(body)
        outs = declared_outputs(body)
        if stub is None:
            if outs:
                FAIL.append(f"{name}: declares outputs but has no stub block")
                report(False, name, "no stub block, so -stub-run cannot exercise it")
            continue
        missing = []
        for o in outs:
            if o.startswith("versions"):
                continue
            hit = False
            for variant in expand_braces(o):
                frags = [f for f in re.split(r"\$\{[^}]*\}", variant) if f and f not in ("/", "*")]
                core = max(frags, key=len).strip("*/") if frags else ""
                if core and core in stub:
                    hit = True
                    break
                # e.g. path("${meta.id}"): no literal to match, so compare the interpolation itself
                if not frags and variant.strip("*/") in stub:
                    hit = True
                    break
            if not hit:
                missing.append(o)
        if missing:
            FAIL.append(f"{name}: stub does not create {', '.join(missing)}")
            report(False, name, f"stub is missing {', '.join(missing)}")
if not any(f.split(":")[0] for f in FAIL):
    report(True, "stub outputs", "every process stub creates what its outputs declare")


# ---------------------------------------------------------------- 2. params are defined
print("\n[2] every params.X referenced is defined in nextflow.config")
cfg = open(os.path.join(ROOT, "nextflow.config")).read()
defined = set(re.findall(r"^\s*([a-z_][a-z0-9_]*)\s*=", cfg, re.M | re.I))
for extra_cfg in glob.glob(os.path.join(ROOT, "conf", "*.config")):
    defined |= set(re.findall(r"^\s*([a-z_][a-z0-9_]*)\s*=", open(extra_cfg).read(), re.M | re.I))
# params.mod.<tool> / params.conda_pkg.<tool> are map members, not top-level names
defined |= {"mod", "conda_pkg"}

referenced = set()
for f in (glob.glob(os.path.join(ROOT, "**", "*.nf"), recursive=True)
          + glob.glob(os.path.join(ROOT, "conf", "*.config"))):
    text_f = open(f).read()
    for m in re.finditer(r"params\.([a-z_][a-z0-9_]*)\s*(.?)", text_f, re.I):
        # params.collectEntries { } / params.each { } are method calls on the params map itself
        if m.group(2) in ("(", "{"):
            continue
        referenced.add(m.group(1))

undefined = sorted(referenced - defined)
if undefined:
    FAIL.append(f"undefined params referenced: {', '.join(undefined)}")
    report(False, "params", f"referenced but never defined: {', '.join(undefined)}")
else:
    report(True, "params", f"{len(referenced)} referenced, all defined")

unused = sorted(d for d in defined - referenced
                if d not in {"mod", "conda_pkg"} and not d.startswith(("max_", "skip_")))
if unused:
    WARN.append(f"defined but never referenced: {', '.join(unused[:12])}")


# ---------------------------------------------------------------- 3. bin/ scripts parse
print("\n[3] every bin/ script parses")
bad = []
skipped_R = []
for f in sorted(glob.glob(os.path.join(ROOT, "bin", "*"))):
    if f.endswith(".py"):
        try:
            ast.parse(open(f).read())
        except SyntaxError as e:
            bad.append(f"{os.path.basename(f)}: {e}")
    elif f.endswith(".R"):
        if shutil.which("Rscript") is None:
            skipped_R.append(os.path.basename(f))     # no R on PATH: report it, do not crash
            continue
        r = subprocess.run(["Rscript", "-e", f'invisible(parse("{f}"))'],
                           capture_output=True, text=True)
        if r.returncode != 0:
            bad.append(f"{os.path.basename(f)}: {r.stderr.strip().splitlines()[-1] if r.stderr else 'parse error'}")
if bad:
    FAIL.extend(bad)
    for b in bad:
        report(False, "script", b)
else:
    n = len(glob.glob(os.path.join(ROOT, "bin", "*.py"))) + len(glob.glob(os.path.join(ROOT, "bin", "*.R")))
    note = f" ({len(skipped_R)} R not checked: no Rscript on PATH)" if skipped_R else ""
    report(True, "scripts", f"{n - len(skipped_R)} of {n} python/R scripts parse{note}")


# ---------------------------------------------------------------- 4. bin/ scripts are executable
print("\n[4] bin/ scripts are executable (Nextflow puts bin/ on PATH)")
notx = [os.path.basename(f) for f in glob.glob(os.path.join(ROOT, "bin", "*"))
        if os.path.isfile(f) and not os.access(f, os.X_OK)]
if notx:
    FAIL.append(f"not executable: {', '.join(notx)}")
    report(False, "bin", f"not executable: {', '.join(notx)}")
else:
    report(True, "bin", "all executable")


# ---------------------------------------------------------------- 5. environment modules exist
print("\n[5] every module named in params.mod is loadable")
mods = dict(re.findall(r"^\s*(\w+)\s*:\s*'([^']+)'", re.search(r"mod\s*=\s*\[(.*?)\]", cfg, re.S).group(1), re.M))
missing_mods = []
for tool, mod in sorted(mods.items()):
    r = subprocess.run(["bash", "-lc", f"module avail {mod} 2>&1"], capture_output=True, text=True)
    if mod.split("/")[0] not in (r.stdout + r.stderr):
        missing_mods.append(f"{tool} -> {mod}")
if missing_mods:
    WARN.append(f"modules not found: {', '.join(missing_mods)}")
    report(False, "modules", f"not found: {', '.join(missing_mods)}")
else:
    report(True, "modules", f"{len(mods)} modules resolve")


# ---------------------------------------------------------------- verdict
print("\n" + "=" * 70)
for w in WARN:
    print(f"  warn: {w}")
if FAIL:
    print(f"  {len(FAIL)} FAILED check(s):")
    for f in FAIL:
        print(f"    - {f}")
    sys.exit(1)
print("  all checks passed")
