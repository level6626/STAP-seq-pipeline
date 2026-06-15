#!/usr/bin/env bash
set -Eeuo pipefail

# Split the existing aligned STAP_TSS_27ac_rep1_S3 Bowtie2 deduplicated BAM by
# R2 methylation code and create per-code read1 TSS BigWig tracks.

PROJECT_DIR="${PROJECT_DIR:-/gpfs/data/zhou-lab/yczhang/methylation}"
PIPELINE_DIR="${PIPELINE_DIR:-${PROJECT_DIR}/STAP-seq-pipeline}"

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-STAP_TSS_27ac_rep1_S3}"
RUN_LABEL="${RUN_LABEL:-STAP_TSS_27ac_rep1_S3_bowtie}"
BAM="${BAM:-${PIPELINE_DIR}/results/standard_tools/${RUN_LABEL}/${SAMPLE}.bowtie2.dedup.bam}"
OUTDIR="${OUTDIR:-${PIPELINE_DIR}/results/standard_tools/${RUN_LABEL}/methylation_split}"
CHROM_SIZES="${CHROM_SIZES:-${PROJECT_DIR}/data/hg38/hg38.chrom.sizes}"

export THREADS SAMPLE BAM OUTDIR CHROM_SIZES

cd "${PIPELINE_DIR}"
scripts/standard_tools/split_methylation_from_bam.sh
