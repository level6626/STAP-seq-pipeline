#!/usr/bin/env bash
set -Eeuo pipefail

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-TAPS_27ac_rep1_S7}"
RUN_DIR="${RUN_DIR:-../data/Run202}"
OUTDIR="${OUTDIR:-results/taps_pipeline/${SAMPLE}}"
CONDA_ENV_BIN="${CONDA_ENV_BIN:-/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin}"

R1="${R1:-${RUN_DIR}/${SAMPLE}_R1_001.fastq.gz}"
R2="${R2:-${RUN_DIR}/${SAMPLE}_R2_001.fastq.gz}"
R3="${R3:-${RUN_DIR}/${SAMPLE}_R3_001.fastq.gz}"
REFERENCE_FASTA="${REFERENCE_FASTA:-../data/hg38/hg38.fa}"

ALIGNER="${ALIGNER:-star}"
ALIGN_READS="${ALIGN_READS:-r3}" # r3 or paired
STAR_INDEX_DIR="${STAR_INDEX_DIR:-/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24}"
STAR_EXTRA="${STAR_EXTRA:---outFilterMultimapNmax 1 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3}"
BOWTIE2_INDEX_PREFIX="${BOWTIE2_INDEX_PREFIX:-${OUTDIR}/reference/bowtie2_index/reference}"
BOWTIE2_EXTRA="${BOWTIE2_EXTRA:---end-to-end --very-sensitive}"
BUILD_BOWTIE2_INDEX="${BUILD_BOWTIE2_INDEX:-auto}"

MAX_READS="${MAX_READS:-0}"
R2_ORIENTATION="${R2_ORIENTATION:-forward}"
FILTER_R2_MOTIF="${FILTER_R2_MOTIF:-}"
TRIM_R1_UMI="${TRIM_R1_UMI:-0}"
FASTP_EXTRA="${FASTP_EXTRA:---dont_eval_duplication}"
MIN_MAPQ="${MIN_MAPQ:-0}"
DEDUP_BY_UMI_IN_COUNTER="${DEDUP_BY_UMI_IN_COUNTER:-0}"

PYTHON="${PYTHON:-${CONDA_ENV_BIN}/python}"
FASTP="${FASTP:-${CONDA_ENV_BIN}/fastp}"
STAR="${STAR:-${CONDA_ENV_BIN}/STAR}"
BOWTIE2="${BOWTIE2:-${CONDA_ENV_BIN}/bowtie2}"
BOWTIE2_BUILD="${BOWTIE2_BUILD:-${CONDA_ENV_BIN}/bowtie2-build}"
SAMTOOLS="${SAMTOOLS:-${CONDA_ENV_BIN}/samtools}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
WORKDIR="${OUTDIR}/work"
mkdir -p "${LOGDIR}" "${WORKDIR}" "$(dirname "${BOWTIE2_INDEX_PREFIX}")"

trap 'echo "ERROR at line ${LINENO}. See ${LOGDIR} for logs." >&2' ERR

log_msg() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOGDIR}/pipeline.progress.log" >&2
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

run_bash_logged() {
    local name="$1"
    local command="$2"
    log_msg "START ${name}"
    {
        printf '[%(%F %T)T] command: %s\n' -1 "${command}"
        bash -o pipefail -c "${command}"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    log_msg "DONE  ${name}"
}

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || { echo "Required file is missing or empty: ${path}" >&2; exit 1; }
}

require_file "${R1}"
require_file "${R2}"
require_file "${R3}"
require_file "${REFERENCE_FASTA}"

if [[ ! -s "${REFERENCE_FASTA}.fai" ]]; then
    run_logged reference_faidx "${SAMTOOLS}" faidx "${REFERENCE_FASTA}"
fi

TAGGED_R1="${WORKDIR}/${SAMPLE}.tagged.R1.fastq.gz"
TAGGED_R3="${WORKDIR}/${SAMPLE}.tagged.R3.fastq.gz"
PREP_STATS="${OUTDIR}/${SAMPLE}.prepare_taps.stats.tsv"

PREP_ARGS=(
    "${PYTHON}" "${SCRIPT_DIR}/prepare_taps_fastqs.py"
    --r1 "${R1}"
    --r2 "${R2}"
    --r3 "${R3}"
    --out-r1 "${TAGGED_R1}"
    --out-r3 "${TAGGED_R3}"
    --stats "${PREP_STATS}"
    --r2-orientation "${R2_ORIENTATION}"
)
if [[ "${MAX_READS}" != "0" ]]; then
    PREP_ARGS+=(--max-reads "${MAX_READS}")
fi
if [[ -n "${FILTER_R2_MOTIF}" ]]; then
    PREP_ARGS+=(--filter-r2-motif "${FILTER_R2_MOTIF}")
fi
if [[ "${TRIM_R1_UMI}" == "1" || "${TRIM_R1_UMI}" == "true" ]]; then
    PREP_ARGS+=(--trim-r1-umi)
fi
run_logged prepare_taps_fastqs "${PREP_ARGS[@]}"

TRIM_R1="${WORKDIR}/${SAMPLE}.tagged.trim.R1.fastq.gz"
TRIM_R3="${WORKDIR}/${SAMPLE}.tagged.trim.R3.fastq.gz"
run_logged fastp_trim \
    "${FASTP}" \
    --thread "${THREADS}" \
    --in1 "${TAGGED_R1}" \
    --in2 "${TAGGED_R3}" \
    --out1 "${TRIM_R1}" \
    --out2 "${TRIM_R3}" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 20 \
    --html "${LOGDIR}/${SAMPLE}.fastp.html" \
    --json "${LOGDIR}/${SAMPLE}.fastp.json" \
    ${FASTP_EXTRA}

SORTED_BAM="${OUTDIR}/${SAMPLE}.${ALIGNER}.${ALIGN_READS}.sorted.bam"
if [[ "${ALIGNER}" == "star" ]]; then
    [[ "${ALIGN_READS}" == "r3" ]] || { echo "STAR mode currently supports ALIGN_READS=r3." >&2; exit 1; }
    require_file "${STAR_INDEX_DIR}/Genome"
    STAR_PREFIX="${WORKDIR}/${SAMPLE}.star."
    rm -rf \
        "${STAR_PREFIX}_STARtmp" \
        "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" \
        "${STAR_PREFIX}Log.final.out" \
        "${STAR_PREFIX}Log.out" \
        "${STAR_PREFIX}Log.progress.out" \
        "${STAR_PREFIX}SJ.out.tab"
    run_logged star_align \
        "${STAR}" \
        --runMode alignReads \
        --genomeLoad NoSharedMemory \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX_DIR}" \
        --readFilesCommand zcat \
        --readFilesIn "${TRIM_R3}" \
        --outFileNamePrefix "${STAR_PREFIX}" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMreadID Standard \
        ${STAR_EXTRA}
    mv "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" "${SORTED_BAM}"
    cp "${STAR_PREFIX}Log.final.out" "${LOGDIR}/${SAMPLE}.STAR.Log.final.out"
elif [[ "${ALIGNER}" == "bowtie2" ]]; then
    if [[ "${BUILD_BOWTIE2_INDEX}" == "1" || "${BUILD_BOWTIE2_INDEX}" == "true" || "${BUILD_BOWTIE2_INDEX}" == "auto" ]]; then
        if [[ ! -s "${BOWTIE2_INDEX_PREFIX}.1.bt2" && ! -s "${BOWTIE2_INDEX_PREFIX}.1.bt2l" ]]; then
            run_logged bowtie2_build "${BOWTIE2_BUILD}" "${REFERENCE_FASTA}" "${BOWTIE2_INDEX_PREFIX}"
        fi
    fi
    if [[ "${ALIGN_READS}" == "paired" ]]; then
        run_bash_logged bowtie2_align_sort \
            "'${BOWTIE2}' ${BOWTIE2_EXTRA} -p ${THREADS} -x '${BOWTIE2_INDEX_PREFIX}' -1 '${TRIM_R1}' -2 '${TRIM_R3}' 2>'${LOGDIR}/bowtie2.stderr.log' | '${SAMTOOLS}' sort -@ ${THREADS} -o '${SORTED_BAM}' -"
    elif [[ "${ALIGN_READS}" == "r3" ]]; then
        run_bash_logged bowtie2_align_sort \
            "'${BOWTIE2}' ${BOWTIE2_EXTRA} -p ${THREADS} -x '${BOWTIE2_INDEX_PREFIX}' -U '${TRIM_R3}' 2>'${LOGDIR}/bowtie2.stderr.log' | '${SAMTOOLS}' sort -@ ${THREADS} -o '${SORTED_BAM}' -"
    else
        echo "ALIGN_READS must be r3 or paired." >&2
        exit 1
    fi
else
    echo "ALIGNER must be star or bowtie2." >&2
    exit 1
fi

run_logged index_sorted "${SAMTOOLS}" index "${SORTED_BAM}"
run_logged flagstat_sorted "${SAMTOOLS}" flagstat "${SORTED_BAM}"

SITE_COUNTS="${OUTDIR}/${SAMPLE}.taps_cpg_sites.tsv"
SUMMARY_COUNTS="${OUTDIR}/${SAMPLE}.taps_meth_code_summary.tsv"
COUNT_ARGS=(
    "${PYTHON}" "${SCRIPT_DIR}/count_taps_cpg_conversion.py"
    --bam "${SORTED_BAM}"
    --reference-fasta "${REFERENCE_FASTA}"
    --out-sites "${SITE_COUNTS}"
    --out-summary "${SUMMARY_COUNTS}"
    --min-mapq "${MIN_MAPQ}"
)
if [[ "${DEDUP_BY_UMI_IN_COUNTER}" == "1" || "${DEDUP_BY_UMI_IN_COUNTER}" == "true" ]]; then
    COUNT_ARGS+=(--dedup-by-umi)
fi
run_logged count_taps_cpg_conversion "${COUNT_ARGS[@]}"

log_msg "Pipeline complete for ${SAMPLE}"
ls -lh "${PREP_STATS}" "${SORTED_BAM}" "${SITE_COUNTS}" "${SUMMARY_COUNTS}" \
    | tee -a "${LOGDIR}/pipeline.progress.log" >&2
