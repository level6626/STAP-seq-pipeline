cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=12 \
MAX_READS=0 \
SAMPLE=TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7_bismark_paired \
ALIGNER=bismark \
ALIGN_READS=paired \
REFERENCE_FASTA=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/hg38.fa \
BISMARK_GENOME_DIR=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38 \
BUILD_BISMARK_INDEX=0 \
BISMARK_PARALLEL=2 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh