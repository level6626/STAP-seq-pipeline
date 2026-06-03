#!/usr/bin/env bash
set -Eeuo pipefail

# Download UCSC hg38 and build a STAR index compatible with the active STAR.

THREADS="${THREADS:-12}"
REFERENCE_DIR="${REFERENCE_DIR:-../data/hg38}"
FASTA_GZ="${FASTA_GZ:-${REFERENCE_DIR}/hg38.fa.gz}"
FASTA="${FASTA:-${REFERENCE_DIR}/hg38.fa}"
FASTA_URL="${FASTA_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz}"
GTF="${GTF:-/gpfs/data/zhou-lab/dcai/data/hg38/gencode.v24.annotation.gtf}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-${REFERENCE_DIR}/STAR_index_2.7.11b_gencode_v24}"
SJDB_OVERHANG="${SJDB_OVERHANG:-100}"
LOGDIR="${LOGDIR:-${REFERENCE_DIR}/logs}"

mkdir -p "${REFERENCE_DIR}" "${STAR_INDEX_DIR}" "${LOGDIR}"

trap 'echo "ERROR at line ${LINENO}. See ${LOGDIR} for logs." >&2' ERR

log_msg() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOGDIR}/build_hg38_star_index.progress.log" >&2
}

require_file() {
    local path="$1"
    [[ -s "${path}" ]] || { echo "Required file is missing or empty: ${path}" >&2; exit 1; }
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || { echo "Required command is not available: $1" >&2; exit 1; }
}

require_command STAR
require_command gzip
if command -v curl >/dev/null 2>&1; then
    DOWNLOAD=(curl --fail --location --retry 3 --continue-at - --output)
elif command -v wget >/dev/null 2>&1; then
    DOWNLOAD=(wget --continue --output-document)
else
    echo "Either curl or wget is required to download hg38.fa.gz." >&2
    exit 1
fi

require_file "${GTF}"

if [[ ! -s "${FASTA_GZ}" ]]; then
    log_msg "Downloading ${FASTA_URL}"
    "${DOWNLOAD[@]}" "${FASTA_GZ}" "${FASTA_URL}" >"${LOGDIR}/download_hg38.log" 2>&1
else
    log_msg "Reusing existing ${FASTA_GZ}"
fi

log_msg "Validating ${FASTA_GZ}"
gzip -t "${FASTA_GZ}"

if [[ ! -s "${FASTA}" ]]; then
    log_msg "Decompressing ${FASTA_GZ}"
    gzip -dc "${FASTA_GZ}" >"${FASTA}.tmp"
    mv -f "${FASTA}.tmp" "${FASTA}"
else
    log_msg "Reusing existing ${FASTA}"
fi

require_file "${FASTA}"

if [[ -s "${STAR_INDEX_DIR}/Genome" && -s "${STAR_INDEX_DIR}/SA" && -s "${STAR_INDEX_DIR}/genomeParameters.txt" ]]; then
    log_msg "STAR index already exists at ${STAR_INDEX_DIR}; remove it manually to rebuild"
    exit 0
fi

log_msg "Building STAR index with $(STAR --version)"
STAR \
    --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX_DIR}" \
    --genomeFastaFiles "${FASTA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    >"${LOGDIR}/STAR_genomeGenerate.log" 2>&1

require_file "${STAR_INDEX_DIR}/Genome"
require_file "${STAR_INDEX_DIR}/SA"
require_file "${STAR_INDEX_DIR}/genomeParameters.txt"
log_msg "STAR index complete: ${STAR_INDEX_DIR}"
