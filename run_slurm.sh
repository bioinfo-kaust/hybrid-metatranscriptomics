#!/usr/bin/env bash
#SBATCH --job-name=nf_driver
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#
# Nextflow driver job for SLURM: submits the pipeline tasks and waits for them.
#
# This exists because the pipeline had been launched from interactive shells, which carry a
# 10-hour wall limit here. TRINITY alone is a 9 h 17 m task (measured), so a from-scratch run
# cannot fit in one: the driver would be killed mid-run and every task with it. The driver itself
# is tiny — it only submits and waits — so 2 cpus is plenty; the 3-day limit is what matters.
#
# Site settings come from the environment, so nothing cluster-specific lives in the repository:
#   PROJ         project directory holding nextflow/ (this repo), nextflow_work/ and logs/
#   PARAMS_FILE  the run's params file (required)
#   SITE_CONFIG  optional config with site paths and the account, if kept out of PARAMS_FILE
#   NF_PROFILE   default: hpc,envmodules
# Pass the account and log paths to sbatch, e.g.
#   PROJ=/path/to/project PARAMS_FILE=my_params.yml sbatch --account=<alloc> --output=<logs>/nf_driver.%j.out \
#       --error=<logs>/nf_driver.%j.err run_slurm.sh [-resume <session-id>] [-c overlay.config]
#
# Anything after the script name is passed straight to `nextflow run`.
set -uo pipefail
PROJ=${PROJ:?set PROJ to the project directory}
PARAMS_FILE=${PARAMS_FILE:?set PARAMS_FILE to the params file for this run}
SITE_CONFIG=${SITE_CONFIG:-}
NF_PROFILE=${NF_PROFILE:-hpc,envmodules}
mkdir -p "$PROJ/logs/slurm"
cd "$PROJ/nextflow"

# The execution report, timeline and trace are per-RUN, not per-pipeline: a second run would
# otherwise overwrite the first one's, and those files are the only record of what was cached, what
# re-ran and how long each task took. Stamped with the SLURM job id so a report can be matched to
# the driver's .out/.err, or with a timestamp when this is run outside SLURM.
STAMP="${SLURM_JOB_ID:-$(date +%Y%m%d-%H%M%S)}"
REPORTS="$PROJ/logs/slurm/nf_run.${STAMP}"

module load nextflow 2>/dev/null || true

# Keep the driver's own scratch off the shared filesystem's metadata hot path, and let a task
# that loses its node be retried rather than taking the run down with it.
export NXF_OPTS='-Xms1g -Xmx8g'
export NXF_ANSI_LOG=false

echo "[$(date)] driver starting on $(hostname)"
echo "[$(date)] reports: ${REPORTS}.{report.html,timeline.html,trace.txt}"
echo "[$(date)] extra args: $*"
echo "[$(date)] nextflow: $(nextflow -version 2>&1 | grep -i version | head -1)"

# The run comes from PARAMS_FILE, plus SITE_CONFIG when the paths are kept separate. A
# command-line --param still overrides either.
[ -s "$PARAMS_FILE" ] || { echo "[ERROR] params file not found: $PARAMS_FILE" >&2; exit 1; }
SITE_ARGS=()
if [ -n "$SITE_CONFIG" ]; then
    [ -s "$SITE_CONFIG" ] || { echo "[ERROR] site config not found: $SITE_CONFIG" >&2; exit 1; }
    SITE_ARGS=(-c "$SITE_CONFIG")
fi
echo "[$(date)] profile: $NF_PROFILE   params: $PARAMS_FILE   site config: ${SITE_CONFIG:-none}"
nextflow run ./main.nf \
    -profile "$NF_PROFILE" \
    -params-file "$PARAMS_FILE" \
    "${SITE_ARGS[@]}" \
    -work-dir "$PROJ/nextflow_work" \
    -with-report "${REPORTS}.report.html" \
    -with-timeline "${REPORTS}.timeline.html" \
    -with-trace "${REPORTS}.trace.txt" \
    "$@"
rc=$?

echo "[$(date)] driver finished, exit ${rc}"
exit $rc
