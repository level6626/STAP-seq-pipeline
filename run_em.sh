cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=12 \
BISMARK_PARALLEL=2 \
EXTRACTOR_PARALLEL=4 \
REUSE_TRIMMED_FASTQS=1 \
BUILD_BISMARK_INDEX=0 \
SAMPLE=EM_pSTAP_cell_27ac_S14 \
OUTDIR=results/emseq_pipeline/EM_pSTAP_cell_27ac_S14 \
REFERENCE_FASTA=../data/hg38/hg38.fa \
scripts/emseq_pipeline/run_emseq_pipeline.sh