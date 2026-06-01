#!/usr/bin/env bash
set -Eeuo pipefail

# Standard-tools comparator for STAP_TSS triplet libraries.
# Run a tiny smoke test with MAX_READS set; unset or set MAX_READS=0 for full data.

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-STAP_TSS_500_oligos_S2}"
RUN_DIR="${RUN_DIR:-../data/Run202}"
OUTDIR="${OUTDIR:-results/standard_tools/${SAMPLE}}"

R1="${R1:-${RUN_DIR}/${SAMPLE}_R1_001.fastq.gz}"
R2="${R2:-${RUN_DIR}/${SAMPLE}_R2_001.fastq.gz}"
R3="${R3:-${RUN_DIR}/${SAMPLE}_R3_001.fastq.gz}"

ALIGNER="${ALIGNER:-bowtie2}" # bowtie2 or star
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
BOWTIE2_INDEX_PREFIX="${BOWTIE2_INDEX_PREFIX:-${OUTDIR}/reference/bowtie2_index/reference}"
BUILD_BOWTIE2_INDEX="${BUILD_BOWTIE2_INDEX:-auto}"
BOWTIE2_EXTRA="${BOWTIE2_EXTRA:---local --very-sensitive-local -X 1000}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-}"
STAR_EXTRA="${STAR_EXTRA:---outFilterMultimapNmax 20 --alignEndsType Local}"

CHROM_SIZES="${CHROM_SIZES:-${OUTDIR}/reference/chrom.sizes}"
CANDIDATE_WINDOWS="${CANDIDATE_WINDOWS:-}"
MAX_READS="${MAX_READS:-0}"
KEEP_INTERMEDIATES="${KEEP_INTERMEDIATES:-0}"
DEDUP_WRITE_STATS="${DEDUP_WRITE_STATS:-0}"

# R1 starts with 8 nt RNA UMI. R2 is the 17 nt molecule/plasmid barcode.
R1_BC_PATTERN="${R1_BC_PATTERN:-NNNNNNNN}"
R2_BC_PATTERN="${R2_BC_PATTERN:-NNNNNNNNNNNNNNNNN}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
WORKDIR="${OUTDIR}/work"
mkdir -p "${LOGDIR}" "${WORKDIR}" "${OUTDIR}/reference"
mkdir -p "$(dirname "${BOWTIE2_INDEX_PREFIX}")"
mkdir -p "$(dirname "${CHROM_SIZES}")"

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
    if [[ ! -s "${path}" ]]; then
        echo "Required file is missing or empty: ${path}" >&2
        exit 1
    fi
}

file_report() {
    local path
    for path in "$@"; do
        [[ -e "${path}" ]] && ls -lh "${path}" | tee -a "${LOGDIR}/pipeline.progress.log" >&2
    done
}

require_file "${R1}"
require_file "${R2}"
require_file "${R3}"

INPUT_R1="${R1}"
INPUT_R2="${R2}"
INPUT_R3="${R3}"

if [[ "${MAX_READS}" != "0" ]]; then
    log_msg "Using first ${MAX_READS} reads only for smoke testing"
    INPUT_R1="${WORKDIR}/${SAMPLE}.subset.R1.fastq.gz"
    INPUT_R2="${WORKDIR}/${SAMPLE}.subset.R2.fastq.gz"
    INPUT_R3="${WORKDIR}/${SAMPLE}.subset.R3.fastq.gz"
    run_logged subset_r1 seqkit head -n "${MAX_READS}" "${R1}" -o "${INPUT_R1}"
    run_logged subset_r2 seqkit head -n "${MAX_READS}" "${R2}" -o "${INPUT_R2}"
    run_logged subset_r3 seqkit head -n "${MAX_READS}" "${R3}" -o "${INPUT_R3}"
    file_report "${INPUT_R1}" "${INPUT_R2}" "${INPUT_R3}"
fi

UMI_R1="${WORKDIR}/${SAMPLE}.R1.umi.fastq.gz"
UMI_R2_DISCARD="${WORKDIR}/${SAMPLE}.R2.extracted.discard.fastq.gz"
UMI_R3="${WORKDIR}/${SAMPLE}.R3.umi.fastq.gz"

run_logged umi_extract_r1_r2 \
    umi_tools extract \
    --extract-method=string \
    --bc-pattern="${R1_BC_PATTERN}" \
    --bc-pattern2="${R2_BC_PATTERN}" \
    -I "${INPUT_R1}" \
    -S "${UMI_R1}" \
    --read2-in="${INPUT_R2}" \
    --read2-out="${UMI_R2_DISCARD}"

run_logged sync_r3_umi_header \
    python "${SCRIPT_DIR}/sync_fastq_headers.py" \
    --template "${UMI_R1}" \
    --sequences "${INPUT_R3}" \
    --output "${UMI_R3}" \
    --keep-sequence-comment

rm -f "${UMI_R2_DISCARD}"
file_report "${UMI_R1}" "${UMI_R3}"

TRIM_R1="${WORKDIR}/${SAMPLE}.R1.umi.trim.fastq.gz"
TRIM_R3="${WORKDIR}/${SAMPLE}.R3.umi.trim.fastq.gz"

run_logged fastp_trim \
    fastp \
    --thread "${THREADS}" \
    --in1 "${UMI_R1}" \
    --in2 "${UMI_R3}" \
    --out1 "${TRIM_R1}" \
    --out2 "${TRIM_R3}" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 20 \
    --html "${LOGDIR}/${SAMPLE}.fastp.html" \
    --json "${LOGDIR}/${SAMPLE}.fastp.json"

file_report "${TRIM_R1}" "${TRIM_R3}"

if [[ ! -s "${CHROM_SIZES}" ]]; then
    if [[ -z "${REFERENCE_FASTA}" ]]; then
        echo "CHROM_SIZES does not exist and REFERENCE_FASTA is not set." >&2
        exit 1
    fi
    run_bash_logged make_chrom_sizes "seqkit fx2tab -n -l '${REFERENCE_FASTA}' > '${CHROM_SIZES}'"
fi

SORTED_BAM="${OUTDIR}/${SAMPLE}.${ALIGNER}.sorted.bam"

case "${ALIGNER}" in
    bowtie2)
        if [[ "${BUILD_BOWTIE2_INDEX}" == "1" || "${BUILD_BOWTIE2_INDEX}" == "true" || "${BUILD_BOWTIE2_INDEX}" == "auto" ]]; then
            if [[ ! -s "${BOWTIE2_INDEX_PREFIX}.1.bt2" && ! -s "${BOWTIE2_INDEX_PREFIX}.1.bt2l" ]]; then
                [[ -n "${REFERENCE_FASTA}" ]] || { echo "REFERENCE_FASTA is required to build a Bowtie2 index." >&2; exit 1; }
                run_logged bowtie2_build bowtie2-build "${REFERENCE_FASTA}" "${BOWTIE2_INDEX_PREFIX}"
            fi
        fi
        run_bash_logged bowtie2_align_sort \
            "bowtie2 ${BOWTIE2_EXTRA} -p ${THREADS} -x '${BOWTIE2_INDEX_PREFIX}' -1 '${TRIM_R1}' -2 '${TRIM_R3}' 2>'${LOGDIR}/bowtie2.stderr.log' | samtools sort -@ ${THREADS} -o '${SORTED_BAM}' -"
        ;;
    star)
        [[ -n "${STAR_INDEX_DIR}" ]] || { echo "STAR_INDEX_DIR is required for ALIGNER=star." >&2; exit 1; }
        STAR_PREFIX="${WORKDIR}/${SAMPLE}.STAR."
        run_logged star_align \
            STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX_DIR}" \
            --readFilesIn "${TRIM_R1}" "${TRIM_R3}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${STAR_PREFIX}" \
            --outSAMtype BAM SortedByCoordinate \
            ${STAR_EXTRA}
        mv "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" "${SORTED_BAM}"
        ;;
    *)
        echo "Unknown ALIGNER=${ALIGNER}; expected bowtie2 or star." >&2
        exit 1
        ;;
esac

run_logged index_sorted samtools index "${SORTED_BAM}"
run_logged flagstat_sorted samtools flagstat "${SORTED_BAM}"
file_report "${SORTED_BAM}" "${SORTED_BAM}.bai"

DEDUP_BAM="${OUTDIR}/${SAMPLE}.${ALIGNER}.dedup.bam"
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
file_report "${DEDUP_BAM}" "${DEDUP_BAM}.bai"

TSS_BED="${OUTDIR}/${SAMPLE}.read1_tss.bed"
TSS_BEDGRAPH="${OUTDIR}/${SAMPLE}.read1_tss.bedGraph"
TSS_BW="${OUTDIR}/${SAMPLE}.read1_tss.bw"

run_bash_logged tss_bed \
    "samtools view -f 64 -F 2820 '${DEDUP_BAM}' | gawk -f '${SCRIPT_DIR}/sam_read1_tss_to_bed.awk' | sort -k1,1 -k2,2n > '${TSS_BED}'"

run_bash_logged tss_bedgraph \
    "bedtools genomecov -bg -i '${TSS_BED}' -g '${CHROM_SIZES}' | sort -k1,1 -k2,2n > '${TSS_BEDGRAPH}'"

if command -v bedGraphToBigWig >/dev/null 2>&1; then
    run_logged tss_bigwig bedGraphToBigWig "${TSS_BEDGRAPH}" "${CHROM_SIZES}" "${TSS_BW}"
fi

if [[ -n "${CANDIDATE_WINDOWS}" ]]; then
    require_file "${CANDIDATE_WINDOWS}"
    run_bash_logged candidate_window_counts \
        "bedtools coverage -counts -a '${CANDIDATE_WINDOWS}' -b '${TSS_BED}' > '${OUTDIR}/${SAMPLE}.candidate_window_tss_counts.tsv'"
fi

run_bash_logged tss_count_summary \
    "printf 'tss_events\\t'; wc -l < '${TSS_BED}'"

if [[ "${KEEP_INTERMEDIATES}" == "0" && "${MAX_READS}" == "0" ]]; then
    log_msg "KEEP_INTERMEDIATES=0, but retaining extracted/trimmed FASTQs for auditability; remove ${WORKDIR} manually when satisfied."
fi

log_msg "Pipeline complete for ${SAMPLE}"
file_report "${TSS_BED}" "${TSS_BEDGRAPH}" "${TSS_BW:-}"
