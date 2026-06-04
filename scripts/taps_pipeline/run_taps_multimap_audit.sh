#!/usr/bin/env bash
set -Eeuo pipefail

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-TAPS_27ac_rep1_S7}"
SOURCE_OUTDIR="${SOURCE_OUTDIR:-results/taps_pipeline/${SAMPLE}}"
OUTDIR="${OUTDIR:-results/taps_pipeline/${SAMPLE}_multimap_audit}"
CONDA_ENV_BIN="${CONDA_ENV_BIN:-/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin}"

TAGGED_TRIM_R3="${TAGGED_TRIM_R3:-${SOURCE_OUTDIR}/work/${SAMPLE}.tagged.trim.R3.fastq.gz}"
PREPARE_STATS="${PREPARE_STATS:-${SOURCE_OUTDIR}/${SAMPLE}.prepare_taps.stats.tsv}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24}"

MAX_MULTIMAP="${MAX_MULTIMAP:-100}"
MIN_MAPQ="${MIN_MAPQ:-0}"
STAR_EXTRA="${STAR_EXTRA:---outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3}"

PYTHON="${PYTHON:-${CONDA_ENV_BIN}/python}"
STAR="${STAR:-${CONDA_ENV_BIN}/STAR}"
SAMTOOLS="${SAMTOOLS:-${CONDA_ENV_BIN}/samtools}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
WORKDIR="${OUTDIR}/work"
mkdir -p "${LOGDIR}" "${WORKDIR}"

trap 'echo "ERROR at line ${LINENO}. See ${LOGDIR} for logs." >&2' ERR

log_msg() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOGDIR}/multimap_audit.progress.log" >&2
}

run_logged() {
    local name="$1"
    shift
    log_msg "START ${name}"
    {
        printf '[%(%F %T)T] command:' -1
        printf ' %q' "$@"
        printf '\n'
        "$@"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    log_msg "DONE  ${name}"
}

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || { echo "Required file is missing or empty: ${path}" >&2; exit 1; }
}

require_file "${TAGGED_TRIM_R3}"
require_file "${PREPARE_STATS}"
require_file "${STAR_INDEX_DIR}/Genome"

STAR_PREFIX="${WORKDIR}/${SAMPLE}.star_multimap."
rm -rf \
    "${STAR_PREFIX}_STARtmp" \
    "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" \
    "${STAR_PREFIX}Log.final.out" \
    "${STAR_PREFIX}Log.out" \
    "${STAR_PREFIX}Log.progress.out" \
    "${STAR_PREFIX}SJ.out.tab"

SORTED_BAM="${OUTDIR}/${SAMPLE}.star.r3.multimap${MAX_MULTIMAP}.sorted.bam"
run_logged star_multimap_align \
    "${STAR}" \
    --runMode alignReads \
    --genomeLoad NoSharedMemory \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX_DIR}" \
    --readFilesCommand zcat \
    --readFilesIn "${TAGGED_TRIM_R3}" \
    --outFileNamePrefix "${STAR_PREFIX}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMreadID Standard \
    --outFilterMultimapNmax "${MAX_MULTIMAP}" \
    --outSAMmultNmax "${MAX_MULTIMAP}" \
    ${STAR_EXTRA}

mv "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" "${SORTED_BAM}"
cp "${STAR_PREFIX}Log.final.out" "${LOGDIR}/${SAMPLE}.STAR.multimap.Log.final.out"
run_logged index_multimap_bam "${SAMTOOLS}" index "${SORTED_BAM}"
run_logged flagstat_multimap_bam "${SAMTOOLS}" flagstat "${SORTED_BAM}"

MAPPING_BY_CODE="${OUTDIR}/${SAMPLE}.mapping_by_meth_code.multimap${MAX_MULTIMAP}.tsv"
run_logged quantify_mapping_by_meth_code \
    "${PYTHON}" "${SCRIPT_DIR}/quantify_mapping_by_meth_code.py" \
    --bam "${SORTED_BAM}" \
    --prepare-stats "${PREPARE_STATS}" \
    --star-final-log "${LOGDIR}/${SAMPLE}.STAR.multimap.Log.final.out" \
    --out "${MAPPING_BY_CODE}" \
    --min-mapq "${MIN_MAPQ}"

log_msg "Multimap audit complete for ${SAMPLE}"
ls -lh "${SORTED_BAM}" "${MAPPING_BY_CODE}" "${LOGDIR}/${SAMPLE}.STAR.multimap.Log.final.out" \
    | tee -a "${LOGDIR}/multimap_audit.progress.log" >&2
