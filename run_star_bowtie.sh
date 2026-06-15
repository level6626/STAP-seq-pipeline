THREADS=12 \
MAX_READS=0 \
SAMPLE=STAP_TSS_BasPromoter_S6 \
OUTDIR=results/standard_tools/STAP_TSS_BasPromoter_S6_bowtie \
ALIGNER=bowtie2 \
REFERENCE_FASTA=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/hg38.fa \
BOWTIE2_INDEX_PREFIX=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/bowtie2_index \
CHROM_SIZES=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/hg38.chrom.sizes \
scripts/standard_tools/run_stap_standard_tools.sh
