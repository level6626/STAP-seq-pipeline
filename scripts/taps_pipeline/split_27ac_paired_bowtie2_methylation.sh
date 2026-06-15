#!/usr/bin/env bash
set -Eeuo pipefail

# Split the existing TAPS_27ac_rep1_S7 paired Bowtie2 BAM by R2 methylation code
# and create per-code alignment-coverage BigWig tracks.

PROJECT_DIR="${PROJECT_DIR:-/gpfs/data/zhou-lab/yczhang/methylation}"
PIPELINE_DIR="${PIPELINE_DIR:-${PROJECT_DIR}/STAP-seq-pipeline}"

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-TAPS_27ac_rep1_S7}"
RUN_LABEL="${RUN_LABEL:-TAPS_27ac_rep1_S7_paired_bowtie2}"
BAM="${BAM:-${PIPELINE_DIR}/results/taps_pipeline/${RUN_LABEL}/${SAMPLE}.bowtie2.paired.sorted.bam}"
OUTDIR="${OUTDIR:-${PIPELINE_DIR}/results/taps_pipeline/${RUN_LABEL}/methylation_split}"
CHROM_SIZES="${CHROM_SIZES:-${PROJECT_DIR}/data/hg38/hg38.chrom.sizes}"

export THREADS SAMPLE BAM OUTDIR CHROM_SIZES

cd "${PIPELINE_DIR}"
scripts/taps_pipeline/split_methylation_from_bam.sh
