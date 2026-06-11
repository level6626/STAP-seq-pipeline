#!/usr/bin/env bash
set -Eeuo pipefail

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-EM_pSTAP_cell_27ac_S14}"
RUN_DIR="${RUN_DIR:-../data/Run188}"
OUTDIR="${OUTDIR:-results/emseq_pipeline/${SAMPLE}}"
CONDA_ENV_BIN="${CONDA_ENV_BIN:-/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin}"

R1="${R1:-${RUN_DIR}/${SAMPLE}_R1_001.fastq.gz}"
R2="${R2:-${RUN_DIR}/${SAMPLE}_R2_001.fastq.gz}"
R3="${R3:-${RUN_DIR}/${SAMPLE}_R3_001.fastq.gz}"

DEFAULT_REFERENCE_FASTA="../data/hg38/hg38.fa"
if [[ "${SAMPLE}" == *500* ]]; then
    DEFAULT_REFERENCE_FASTA="/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa"
fi
REFERENCE_FASTA="${REFERENCE_FASTA:-${DEFAULT_REFERENCE_FASTA}}"

PYTHON="${PYTHON:-${CONDA_ENV_BIN}/python}"
FASTP="${FASTP:-${CONDA_ENV_BIN}/fastp}"
SAMTOOLS="${SAMTOOLS:-${CONDA_ENV_BIN}/samtools}"
UMI_TOOLS="${UMI_TOOLS:-${CONDA_ENV_BIN}/umi_tools}"
BISMARK="${BISMARK:-${CONDA_ENV_BIN}/bismark}"
BISMARK_GENOME_PREPARATION="${BISMARK_GENOME_PREPARATION:-${CONDA_ENV_BIN}/bismark_genome_preparation}"
BISMARK_METHYLATION_EXTRACTOR="${BISMARK_METHYLATION_EXTRACTOR:-${CONDA_ENV_BIN}/bismark_methylation_extractor}"
BOWTIE2="${BOWTIE2:-${CONDA_ENV_BIN}/bowtie2}"

MAX_READS="${MAX_READS:-0}"
R2_ORIENTATION="${R2_ORIENTATION:-forward}"
R1_UMI_LENGTH="${R1_UMI_LENGTH:-8}"
EXPECTED_R1_LENGTH="${EXPECTED_R1_LENGTH:-55}"
UMI_SOURCE="${UMI_SOURCE:-r1+r2_barcode}"
KEEP_R1_UMI="${KEEP_R1_UMI:-0}"
FASTP_EXTRA="${FASTP_EXTRA:---dont_eval_duplication}"
BISMARK_EXTRA="${BISMARK_EXTRA:-}"
BISMARK_EXTRACTOR_EXTRA="${BISMARK_EXTRACTOR_EXTRA:---bedGraph --CX_context --comprehensive --gzip}"
BISMARK_PARALLEL="${BISMARK_PARALLEL:-1}"
EXTRACTOR_PARALLEL="${EXTRACTOR_PARALLEL:-1}"
PAD_OLIGO_REFERENCE="${PAD_OLIGO_REFERENCE:-auto}"
OLIGO_REFERENCE_RIGHT_PAD="${OLIGO_REFERENCE_RIGHT_PAD:-100}"
OLIGO_REFERENCE_PAD_BASE="${OLIGO_REFERENCE_PAD_BASE:-N}"
BUILD_BISMARK_INDEX="${BUILD_BISMARK_INDEX:-auto}"
RUN_BISMARK="${RUN_BISMARK:-1}"
RUN_DEDUP="${RUN_DEDUP:-1}"
RUN_METHYLATION_EXTRACTOR="${RUN_METHYLATION_EXTRACTOR:-1}"
DEDUP_WRITE_STATS="${DEDUP_WRITE_STATS:-0}"
REUSE_TRIMMED_FASTQS="${REUSE_TRIMMED_FASTQS:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
WORKDIR="${OUTDIR}/work"
mkdir -p "${LOGDIR}" "${WORKDIR}"

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

require_executable() {
    local path="$1"
    [[ -x "${path}" ]] || { echo "Required executable is missing: ${path}" >&2; exit 1; }
}

require_file "${R1}"
require_file "${R2}"
require_file "${R3}"
require_file "${REFERENCE_FASTA}"
require_executable "${PYTHON}"
require_executable "${FASTP}"
require_executable "${SAMTOOLS}"
require_executable "${UMI_TOOLS}"

TAGGED_R1="${WORKDIR}/${SAMPLE}.emseq.tagged.R1.fastq.gz"
TAGGED_R3="${WORKDIR}/${SAMPLE}.emseq.tagged.R3.fastq.gz"
PREP_STATS="${OUTDIR}/${SAMPLE}.prepare_emseq.stats.tsv"
TRIM_R1="${WORKDIR}/${SAMPLE}.emseq.tagged.trim.R1.fastq.gz"
TRIM_R3="${WORKDIR}/${SAMPLE}.emseq.tagged.trim.R3.fastq.gz"

if [[ "${REUSE_TRIMMED_FASTQS}" == "1" || "${REUSE_TRIMMED_FASTQS}" == "true" ]]; then
    require_file "${PREP_STATS}"
    require_file "${TRIM_R1}"
    require_file "${TRIM_R3}"
    log_msg "REUSE_TRIMMED_FASTQS=1; reusing ${TRIM_R1} and ${TRIM_R3}"
else
    PREP_ARGS=(
        "${PYTHON}" "${SCRIPT_DIR}/prepare_emseq_fastqs.py"
        --r1 "${R1}"
        --r2 "${R2}"
        --r3 "${R3}"
        --out-r1 "${TAGGED_R1}"
        --out-r3 "${TAGGED_R3}"
        --stats "${PREP_STATS}"
        --r2-orientation "${R2_ORIENTATION}"
        --r1-umi-length "${R1_UMI_LENGTH}"
        --expected-r1-length "${EXPECTED_R1_LENGTH}"
        --umi-source "${UMI_SOURCE}"
    )
    if [[ "${MAX_READS}" != "0" ]]; then
        PREP_ARGS+=(--max-reads "${MAX_READS}")
    fi
    if [[ "${KEEP_R1_UMI}" == "1" || "${KEEP_R1_UMI}" == "true" ]]; then
        PREP_ARGS+=(--keep-r1-umi)
    fi
    run_logged prepare_emseq_fastqs "${PREP_ARGS[@]}"

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
fi

if [[ "${RUN_BISMARK}" != "1" && "${RUN_BISMARK}" != "true" ]]; then
    log_msg "RUN_BISMARK=0; stopping after FASTQ preparation and trimming"
    ls -lh "${PREP_STATS}" "${TRIM_R1}" "${TRIM_R3}" | tee -a "${LOGDIR}/pipeline.progress.log" >&2
    exit 0
fi

require_executable "${BISMARK}"
require_executable "${BISMARK_GENOME_PREPARATION}"
require_executable "${BISMARK_METHYLATION_EXTRACTOR}"
require_executable "${BOWTIE2}"

REF_BASENAME="$(basename "${REFERENCE_FASTA}")"
USE_PADDED_REFERENCE=0
if [[ "${PAD_OLIGO_REFERENCE}" == "1" || "${PAD_OLIGO_REFERENCE}" == "true" ]]; then
    USE_PADDED_REFERENCE=1
elif [[ "${PAD_OLIGO_REFERENCE}" == "auto" && ( "${REF_BASENAME}" == data_500.fa || "${REFERENCE_FASTA}" == */data_500/* ) ]]; then
    USE_PADDED_REFERENCE=1
fi
if [[ "${USE_PADDED_REFERENCE}" == "1" ]]; then
    PADDED_REFERENCE_DIR="${OUTDIR}/reference/right_padded_oligo"
    PADDED_REFERENCE_FASTA="${PADDED_REFERENCE_DIR}/${REF_BASENAME%.fa}.rightpad${OLIGO_REFERENCE_RIGHT_PAD}${OLIGO_REFERENCE_PAD_BASE}.fa"
    run_logged write_right_padded_oligo_fasta \
        "${PYTHON}" "${SCRIPT_DIR}/write_right_padded_fasta.py" \
        --input "${REFERENCE_FASTA}" \
        --out "${PADDED_REFERENCE_FASTA}" \
        --right-pad "${OLIGO_REFERENCE_RIGHT_PAD}" \
        --pad-base "${OLIGO_REFERENCE_PAD_BASE}"
    REFERENCE_FASTA="${PADDED_REFERENCE_FASTA}"
    BISMARK_GENOME_DIR="${PADDED_REFERENCE_DIR}"
    log_msg "Using right-padded oligo reference for Bismark: ${REFERENCE_FASTA}"
fi

if [[ -z "${BISMARK_GENOME_DIR:-}" ]]; then
    ref_dir="$(cd "$(dirname "${REFERENCE_FASTA}")" && pwd)"
    ref_base="$(basename "${REFERENCE_FASTA}")"
    if [[ -w "${ref_dir}" ]]; then
        BISMARK_GENOME_DIR="${ref_dir}"
    else
        BISMARK_GENOME_DIR="${OUTDIR}/reference/bismark_genome"
        mkdir -p "${BISMARK_GENOME_DIR}"
        ln -sf "$(cd "$(dirname "${REFERENCE_FASTA}")" && pwd)/${ref_base}" "${BISMARK_GENOME_DIR}/${ref_base}"
    fi
fi
mkdir -p "${BISMARK_GENOME_DIR}"

if [[ ! -e "${BISMARK_GENOME_DIR}/$(basename "${REFERENCE_FASTA}")" ]]; then
    ln -sf "$(cd "$(dirname "${REFERENCE_FASTA}")" && pwd)/$(basename "${REFERENCE_FASTA}")" \
        "${BISMARK_GENOME_DIR}/$(basename "${REFERENCE_FASTA}")"
fi

if [[ "${BUILD_BISMARK_INDEX}" == "1" || "${BUILD_BISMARK_INDEX}" == "true" || "${BUILD_BISMARK_INDEX}" == "auto" ]]; then
    if [[ ! -s "${BISMARK_GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2" && ! -s "${BISMARK_GENOME_DIR}/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2l" ]]; then
        run_logged bismark_genome_preparation \
            "${BISMARK_GENOME_PREPARATION}" \
            --bowtie2 \
            --path_to_aligner "$(dirname "${BOWTIE2}")" \
            "${BISMARK_GENOME_DIR}"
    fi
fi

BISMARK_DIR="${OUTDIR}/bismark"
mkdir -p "${BISMARK_DIR}"
BISMARK_ARGS=(
    "${BISMARK}"
    --genome "${BISMARK_GENOME_DIR}"
    --output_dir "${BISMARK_DIR}"
    --path_to_bowtie2 "$(dirname "${BOWTIE2}")"
    --samtools_path "$(dirname "${SAMTOOLS}")"
)
if [[ "${BISMARK_PARALLEL}" == "1" ]]; then
    BISMARK_ARGS+=(--basename "${SAMPLE}")
else
    log_msg "BISMARK_PARALLEL=${BISMARK_PARALLEL}; omitting --basename because Bismark rejects --basename with multicore mode"
    BISMARK_ARGS+=(--parallel "${BISMARK_PARALLEL}")
fi
BISMARK_ARGS+=(
    -1 "${TRIM_R1}"
    -2 "${TRIM_R3}"
)
if [[ -n "${BISMARK_EXTRA}" ]]; then
    # shellcheck disable=SC2206
    BISMARK_EXTRA_ARGS=(${BISMARK_EXTRA})
    BISMARK_ARGS+=("${BISMARK_EXTRA_ARGS[@]}")
fi
run_logged bismark_align "${BISMARK_ARGS[@]}"

BISMARK_BAM_CANDIDATES=("${BISMARK_DIR}"/*_bismark_bt2_pe.bam)
if [[ ${#BISMARK_BAM_CANDIDATES[@]} -ne 1 || ! -s "${BISMARK_BAM_CANDIDATES[0]}" ]]; then
    echo "Expected exactly one Bismark paired-end BAM in ${BISMARK_DIR}, found ${#BISMARK_BAM_CANDIDATES[@]}" >&2
    printf 'Candidate: %s\n' "${BISMARK_BAM_CANDIDATES[@]}" >&2
    exit 1
fi
BISMARK_BAM="${BISMARK_BAM_CANDIDATES[0]}"
SORTED_BAM="${OUTDIR}/${SAMPLE}.bismark.sorted.bam"
run_logged sort_bismark_bam "${SAMTOOLS}" sort -@ "${THREADS}" -o "${SORTED_BAM}" "${BISMARK_BAM}"
run_logged index_sorted "${SAMTOOLS}" index "${SORTED_BAM}"
run_logged flagstat_sorted "${SAMTOOLS}" flagstat "${SORTED_BAM}"

DEDUP_BAM="${OUTDIR}/${SAMPLE}.bismark.dedup.bam"
if [[ "${RUN_DEDUP}" == "1" || "${RUN_DEDUP}" == "true" ]]; then
    DEDUP_ARGS=(
        "${UMI_TOOLS}" dedup
        --paired
        --extract-umi-method=read_id
        --umi-separator=_
        --stdin "${SORTED_BAM}"
        --stdout "${DEDUP_BAM}"
        --log "${LOGDIR}/${SAMPLE}.umi_tools.dedup.log"
    )
    if [[ "${DEDUP_WRITE_STATS}" == "1" || "${DEDUP_WRITE_STATS}" == "true" ]]; then
        DEDUP_ARGS+=(--output-stats "${OUTDIR}/${SAMPLE}.umi_tools.dedup")
    fi
    run_logged umi_tools_dedup "${DEDUP_ARGS[@]}"
    run_logged index_dedup "${SAMTOOLS}" index "${DEDUP_BAM}"
else
    DEDUP_BAM="${SORTED_BAM}"
fi

if [[ "${RUN_METHYLATION_EXTRACTOR}" == "1" || "${RUN_METHYLATION_EXTRACTOR}" == "true" ]]; then
    EXTRACT_DIR="${OUTDIR}/methylation_extractor"
    EXTRACT_BAM="${WORKDIR}/${SAMPLE}.bismark.extract.name_sorted.bam"
    mkdir -p "${EXTRACT_DIR}"
    run_logged name_sort_for_methylation_extractor \
        "${SAMTOOLS}" sort -n -@ "${THREADS}" -o "${EXTRACT_BAM}" "${DEDUP_BAM}"
    run_logged bismark_methylation_extractor \
        "${BISMARK_METHYLATION_EXTRACTOR}" \
        --paired-end \
        --parallel "${EXTRACTOR_PARALLEL}" \
        --genome_folder "${BISMARK_GENOME_DIR}" \
        --samtools_path "$(dirname "${SAMTOOLS}")" \
        --output "${EXTRACT_DIR}" \
        ${BISMARK_EXTRACTOR_EXTRA} \
        "${EXTRACT_BAM}"
fi

log_msg "Pipeline complete for ${SAMPLE}"
ls -lh "${PREP_STATS}" "${SORTED_BAM}" "${DEDUP_BAM}" | tee -a "${LOGDIR}/pipeline.progress.log" >&2
