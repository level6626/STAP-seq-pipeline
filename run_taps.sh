THREADS=12 \
MAX_READS=0 \
SAMPLE=TAPS_27ac_rep2_S8 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep2_S8 \
R1=../data/Run202/TAPS_27ac_rep2_S8_R1_001.filter.fastq.gz \
R2=../data/Run202/TAPS_27ac_rep2_S8_R2_001.filter.fastq.gz \
R3=../data/Run202/TAPS_27ac_rep2_S8_R3_001.filter.fastq.gz \
ALIGNER=star \
ALIGN_READS=r3 \
REFERENCE_FASTA=../data/hg38/hg38.fa \
STAR_INDEX_DIR=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh