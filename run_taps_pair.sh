THREADS=12 \
SAMPLE=TAPS_BasPromoter_S10 \
R1=../data/Run202/TAPS_BasPromoter_S10_R1_001.filter.fastq.gz \
R2=../data/Run202/TAPS_BasPromoter_S10_R2_001.filter.fastq.gz \
R3=../data/Run202/TAPS_BasPromoter_S10_R3_001.filter.fastq.gz \
OUTDIR=results/taps_pipeline/TAPS_BasPromoter_S10_paired_bowtie2 \
ALIGNER=bowtie2 \
ALIGN_READS=paired \
REFERENCE_FASTA=../data/hg38/hg38.fa \
BOWTIE2_INDEX_PREFIX=../data/hg38/bowtie2_index/hg38 \
BUILD_BOWTIE2_INDEX=0 \
TRIM_R1_UMI=1 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh