THREADS=12 \
SAMPLE=TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7_paired_bowtie2 \
ALIGNER=bowtie2 \
ALIGN_READS=paired \
REFERENCE_FASTA=../data/hg38/hg38.fa \
BOWTIE2_INDEX_PREFIX=../data/hg38/bowtie2_index/hg38 \
BUILD_BOWTIE2_INDEX=0 \
TRIM_R1_UMI=1 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh