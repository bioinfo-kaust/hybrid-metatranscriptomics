#!/usr/bin/env bash
# Write this task's versions.yml.   Usage: emit_versions.sh "<process>" [tool ...]
#
# Provenance must never be the thing that breaks a run, so every lookup is best-effort: a tool
# that is absent, or that refuses to report a version, is recorded as "n/a" rather than failing
# the task. The file is always written, so the output declaration can stay non-optional and the
# collector downstream never has to special-case a missing file.
set -uo pipefail
proc="${1:-unknown}"; shift || true

ver () {
  local t=$1 v=""
  command -v "$t" >/dev/null 2>&1 || { echo "n/a (not on PATH)"; return; }
  case "$t" in
    salmon|vsearch|stringtie|fastp|multiqc|busco|corset)
                      v=$("$t" --version 2>&1 | head -1) ;;
    minimap2|seqkit|samtools|pbmm2|diamond|gffread|gffcompare|bowtie2|hisat2)
                      v=$("$t" --version 2>&1 | head -1) ;;
    STAR)             v=$(STAR --version 2>&1 | head -1) ;;
    Trinity)          v=$(Trinity --version 2>&1 | grep -m1 -i version) ;;
    TransDecoder.LongOrfs|TransDecoder.Predict)
                      v=$("$t" --version 2>&1 | head -1) ;;
    emapper.py)       v=$(emapper.py --version 2>&1 | head -1) ;;
    interproscan.sh)  v=$(interproscan.sh --version 2>&1 | grep -m1 -i version) ;;
    fastqc)           v=$(fastqc --version 2>&1 | head -1) ;;
    bbduk.sh)         v=$(bbduk.sh --version 2>&1 | grep -m1 -i bbmap) ;;
    sqanti3_qc.py)    v=$(sqanti3_qc.py --version 2>&1 | head -1) ;;
    python3)          v=$(python3 --version 2>&1 | head -1) ;;
    Rscript)          v=$(Rscript --version 2>&1 | head -1) ;;
    *)                v=$("$t" --version 2>&1 | head -1) ;;
  esac
  # strip the tool's own name and any leading "version", collapse whitespace
  v=$(printf '%s' "$v" | sed -E "s/^[Vv]ersion[: ]*//; s/^$t[ ,v]*//I; s/[[:space:]]+/ /g" | head -c 200)
  [ -n "$v" ] && printf '%s' "$v" || printf 'n/a (no version output)'
}

{
  printf '"%s":\n' "$proc"
  for t in "$@"; do printf '    %s: "%s"\n' "$t" "$(ver "$t")"; done
  # a process with no external tool still gets an entry, so every task yields a parseable file
  if [ "$#" -eq 0 ]; then printf '    bash: "%s"\n' "${BASH_VERSION:-n/a}"; fi
} > versions.yml
# never let provenance be the thing that fails a task: callers run under `bash -ue`
exit 0
