#!/usr/bin/env bash
set -Eeuo pipefail

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-STAP_TSS_500_oligos_S2}"
RUN_DIR="${RUN_DIR:-../data/Run202}"
OUTDIR="${OUTDIR:-results/barcode_pipeline/${SAMPLE}}"

R1="${R1:-${RUN_DIR}/${SAMPLE}_R1_001.fastq.gz}"
R2="${R2:-${RUN_DIR}/${SAMPLE}_R2_001.fastq.gz}"
R3="${R3:-${RUN_DIR}/${SAMPLE}_R3_001.fastq.gz}"
OLIGO_XLSX="${OLIGO_XLSX:-../data/meta/STAP_Seq_oligos.xlsx}"

REFERENCE_FASTA="${REFERENCE_FASTA:-}"
BOWTIE2_INDEX_PREFIX="${BOWTIE2_INDEX_PREFIX:-${OUTDIR}/reference/bowtie2_index/reference}"
BUILD_BOWTIE2_INDEX="${BUILD_BOWTIE2_INDEX:-auto}"
BOWTIE2_EXTRA="${BOWTIE2_EXTRA:---end-to-end --very-sensitive}"

MAX_READS="${MAX_READS:-0}"
BARCODE_SEARCH_BASES="${BARCODE_SEARCH_BASES:-0}"
BARCODE_ORIENTATION="${BARCODE_ORIENTATION:-both}"
MAX_BARCODE_MISMATCHES="${MAX_BARCODE_MISMATCHES:-1}"
KEEP_R3_BARCODE="${KEEP_R3_BARCODE:-1}"

MIN_MAPQ="${MIN_MAPQ:-0}"
DEDUP_WRITE_STATS="${DEDUP_WRITE_STATS:-0}"

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
require_file "${OLIGO_XLSX}"

DEMUX_R1="${WORKDIR}/${SAMPLE}.demux.R1.fastq.gz"
DEMUX_R3="${WORKDIR}/${SAMPLE}.demux.R3.fastq.gz"
DEMUX_STATS="${OUTDIR}/${SAMPLE}.demux.stats.tsv"
OLIGO_METADATA="${OUTDIR}/${SAMPLE}.oligo_metadata.tsv"

DEMUX_ARGS=(
    python "${SCRIPT_DIR}/demux_extract_stap_triplets.py"
    --r1 "${R1}"
    --r2 "${R2}"
    --r3 "${R3}"
    --oligo-xlsx "${OLIGO_XLSX}"
    --out-r1 "${DEMUX_R1}"
    --out-r3 "${DEMUX_R3}"
    --stats "${DEMUX_STATS}"
    --metadata-out "${OLIGO_METADATA}"
    --max-barcode-mismatches "${MAX_BARCODE_MISMATCHES}"
    --barcode-search-bases "${BARCODE_SEARCH_BASES}"
    --barcode-orientation "${BARCODE_ORIENTATION}"
)
if [[ "${MAX_READS}" != "0" ]]; then
    DEMUX_ARGS+=(--max-reads "${MAX_READS}")
fi
if [[ "${KEEP_R3_BARCODE}" == "1" || "${KEEP_R3_BARCODE}" == "true" ]]; then
    DEMUX_ARGS+=(--keep-r3-barcode)
else
    DEMUX_ARGS+=(--trim-r3-barcode)
fi
run_logged demux_extract "${DEMUX_ARGS[@]}"

TRIM_R1="${WORKDIR}/${SAMPLE}.demux.trim.R1.fastq.gz"
TRIM_R3="${WORKDIR}/${SAMPLE}.demux.trim.R3.fastq.gz"
run_logged fastp_trim \
    fastp \
    --thread "${THREADS}" \
    --in1 "${DEMUX_R1}" \
    --in2 "${DEMUX_R3}" \
    --out1 "${TRIM_R1}" \
    --out2 "${TRIM_R3}" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 20 \
    --html "${LOGDIR}/${SAMPLE}.fastp.html" \
    --json "${LOGDIR}/${SAMPLE}.fastp.json"

if [[ "${BUILD_BOWTIE2_INDEX}" == "1" || "${BUILD_BOWTIE2_INDEX}" == "true" || "${BUILD_BOWTIE2_INDEX}" == "auto" ]]; then
    if [[ ! -s "${BOWTIE2_INDEX_PREFIX}.1.bt2" && ! -s "${BOWTIE2_INDEX_PREFIX}.1.bt2l" ]]; then
        [[ -n "${REFERENCE_FASTA}" ]] || { echo "REFERENCE_FASTA is required to build a Bowtie2 index." >&2; exit 1; }
        run_logged bowtie2_build bowtie2-build "${REFERENCE_FASTA}" "${BOWTIE2_INDEX_PREFIX}"
    fi
fi

SORTED_BAM="${OUTDIR}/${SAMPLE}.bowtie2.sorted.bam"
run_bash_logged bowtie2_align_sort \
    "bowtie2 ${BOWTIE2_EXTRA} -p ${THREADS} -x '${BOWTIE2_INDEX_PREFIX}' -1 '${TRIM_R1}' -2 '${TRIM_R3}' 2>'${LOGDIR}/bowtie2.stderr.log' | samtools sort -@ ${THREADS} -o '${SORTED_BAM}' -"
run_logged index_sorted samtools index "${SORTED_BAM}"
run_logged flagstat_sorted samtools flagstat "${SORTED_BAM}"

DEDUP_BAM="${OUTDIR}/${SAMPLE}.bowtie2.dedup.bam"
DEDUP_ARGS=(
    umi_tools dedup
    --paired
    --stdin="${SORTED_BAM}"
    --stdout="${DEDUP_BAM}"
)
if [[ "${DEDUP_WRITE_STATS}" == "1" || "${DEDUP_WRITE_STATS}" == "true" ]]; then
    DEDUP_ARGS+=(--output-stats="${OUTDIR}/${SAMPLE}.dedup")
fi
run_logged umi_dedup "${DEDUP_ARGS[@]}"
run_logged index_dedup samtools index "${DEDUP_BAM}"
run_logged flagstat_dedup samtools flagstat "${DEDUP_BAM}"

TSS_COUNTS="${OUTDIR}/${SAMPLE}.tss_by_oligo_meth.tsv"
run_logged tss_count_by_oligo_meth \
    python "${SCRIPT_DIR}/count_tss_by_oligo_meth.py" \
    --bam "${DEDUP_BAM}" \
    --out "${TSS_COUNTS}" \
    --min-mapq "${MIN_MAPQ}"

log_msg "Pipeline complete for ${SAMPLE}"
ls -lh "${DEMUX_STATS}" "${OLIGO_METADATA}" "${SORTED_BAM}" "${DEDUP_BAM}" "${TSS_COUNTS}" \
    | tee -a "${LOGDIR}/pipeline.progress.log" >&2
