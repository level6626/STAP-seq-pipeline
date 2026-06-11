#!/usr/bin/env bash
set -Eeuo pipefail

SAMPLE="${SAMPLE:-TAPS_27ac_rep1_S7}"
OUTDIR="${OUTDIR:-results/taps_pipeline/${SAMPLE}}"
CONDA_ENV_BIN="${CONDA_ENV_BIN:-/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin}"

BAM="${BAM:-${OUTDIR}/${SAMPLE}.star.r3.sorted.bam}"
REFERENCE_FASTA="${REFERENCE_FASTA:-../data/hg38/hg38.fa}"
OUT_SUMMARY="${OUT_SUMMARY:-${OUTDIR}/${SAMPLE}.taps_noncpg_context_summary.tsv}"
MIN_MAPQ="${MIN_MAPQ:-0}"
INCLUDE_CPG="${INCLUDE_CPG:-0}"
MAX_RECORDS="${MAX_RECORDS:-0}"
PROGRESS_EVERY="${PROGRESS_EVERY:-1000000}"

PYTHON="${PYTHON:-${CONDA_ENV_BIN}/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${LOGDIR}"

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || { echo "Required file is missing or empty: ${path}" >&2; exit 1; }
}

require_file "${BAM}"
require_file "${REFERENCE_FASTA}"

ARGS=(
    "${PYTHON}" "${SCRIPT_DIR}/count_taps_noncpg_conversion.py"
    --bam "${BAM}"
    --reference-fasta "${REFERENCE_FASTA}"
    --out-summary "${OUT_SUMMARY}"
    --min-mapq "${MIN_MAPQ}"
    --max-records "${MAX_RECORDS}"
    --progress-every "${PROGRESS_EVERY}"
)
if [[ "${INCLUDE_CPG}" == "1" || "${INCLUDE_CPG}" == "true" ]]; then
    ARGS+=(--include-cpg)
fi

{
    printf '[%(%F %T)T] command:' -1
    printf ' %q' "${ARGS[@]}"
    printf '\n'
    "${ARGS[@]}"
    printf '[%(%F %T)T] done\n' -1
} >"${LOGDIR}/count_taps_noncpg_conversion.log" 2>&1

ls -lh "${OUT_SUMMARY}"
