#!/usr/bin/env bash
set -Eeuo pipefail

SAMPLE="${SAMPLE:-EM_pSTAP_cell_27ac_S14}"
OUTDIR="${OUTDIR:-results/emseq_pipeline/${SAMPLE}}"
CONDA_ENV_BIN="${CONDA_ENV_BIN:-/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin}"
MIN_MAPQ="${MIN_MAPQ:-0}"
OLIGO_BARCODE_LENGTH="${OLIGO_BARCODE_LENGTH:-10}"
RUN_OLIGO_BARCODE_QC="${RUN_OLIGO_BARCODE_QC:-auto}"
OLIGO_REFERENCE_FASTA="${OLIGO_REFERENCE_FASTA:-}"

PYTHON="${PYTHON:-${CONDA_ENV_BIN}/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${LOGDIR}"

SORTED_BAM="${BAM:-${OUTDIR}/${SAMPLE}.bismark.sorted.bam}"
DEDUP_BAM="${DEDUP_BAM:-${OUTDIR}/${SAMPLE}.bismark.dedup.bam}"
TRIM_R1="${TRIM_R1:-${OUTDIR}/work/${SAMPLE}.emseq.tagged.trim.R1.fastq.gz}"

run_logged() {
    local name="$1"
    shift
    printf '[%(%F %T)T] START %s\n' -1 "${name}" | tee -a "${LOGDIR}/emseq_qc.progress.log" >&2
    {
        printf '[%(%F %T)T] command:' -1
        printf ' %q' "$@"
        printf '\n'
        "$@"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    printf '[%(%F %T)T] DONE  %s\n' -1 "${name}" | tee -a "${LOGDIR}/emseq_qc.progress.log" >&2
}

[[ -s "${SORTED_BAM}" ]] || { echo "Missing BAM: ${SORTED_BAM}" >&2; exit 1; }
[[ -s "${DEDUP_BAM}" ]] || { echo "Missing dedup BAM: ${DEDUP_BAM}" >&2; exit 1; }
[[ -s "${TRIM_R1}" ]] || { echo "Missing trimmed R1 FASTQ: ${TRIM_R1}" >&2; exit 1; }

run_logged quantify_emseq_mapping_by_meth_code \
    "${PYTHON}" "${SCRIPT_DIR}/quantify_emseq_mapping_by_meth_code.py" \
    --bam "${SORTED_BAM}" \
    --trimmed-r1 "${TRIM_R1}" \
    --out "${OUTDIR}/${SAMPLE}.emseq_mapping_by_meth_code.tsv" \
    --min-mapq "${MIN_MAPQ}"

run_logged count_emseq_xm_conversion \
    "${PYTHON}" "${SCRIPT_DIR}/count_emseq_xm_conversion.py" \
    --bam "${DEDUP_BAM}" \
    --out-cpg-summary "${OUTDIR}/${SAMPLE}.emseq_cpg_conversion_by_meth_code.tsv" \
    --out-noncpg-summary "${OUTDIR}/${SAMPLE}.emseq_noncpg_conversion_by_meth_code.tsv" \
    --min-mapq "${MIN_MAPQ}"

if [[ -z "${OLIGO_REFERENCE_FASTA}" && "${SAMPLE}" == *500* ]]; then
    OLIGO_REFERENCE_FASTA="/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa"
fi
if [[ "${RUN_OLIGO_BARCODE_QC}" == "1" || "${RUN_OLIGO_BARCODE_QC}" == "true" || ( "${RUN_OLIGO_BARCODE_QC}" == "auto" && -n "${OLIGO_REFERENCE_FASTA}" ) ]]; then
    [[ -s "${OLIGO_REFERENCE_FASTA}" ]] || { echo "Missing oligo reference FASTA: ${OLIGO_REFERENCE_FASTA}" >&2; exit 1; }
    run_logged quantify_oligo_barcode_mapping_accuracy \
        "${PYTHON}" "${SCRIPT_DIR}/quantify_oligo_barcode_mapping_accuracy.py" \
        --bam "${SORTED_BAM}" \
        --oligo-fasta "${OLIGO_REFERENCE_FASTA}" \
        --out-summary "${OUTDIR}/${SAMPLE}.emseq_oligo_barcode_mapping_accuracy.tsv" \
        --barcode-length "${OLIGO_BARCODE_LENGTH}" \
        --min-mapq "${MIN_MAPQ}"
fi

printf '[%(%F %T)T] EM-seq QC complete for %s\n' -1 "${SAMPLE}" | tee -a "${LOGDIR}/emseq_qc.progress.log" >&2
