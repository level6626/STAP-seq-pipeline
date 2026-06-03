THREADS=12 \
MAX_READS=0 \
SAMPLE=STAP_TSS_27ac_rep1_S3 \
OUTDIR=results/standard_tools/STAP_TSS_27ac_rep1_S3 \
ALIGNER=star \
REUSE_UMI_FASTQS=1 \
FASTP_EXTRA=--dont_eval_duplication \
STAR_INDEX_DIR=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24 \
CHROM_SIZES=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/hg38.chrom.sizes \
scripts/standard_tools/run_stap_standard_tools.sh
