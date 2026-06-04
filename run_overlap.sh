cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

python \
  scripts/compare_raw_r2_barcodes.py \
  --stap-r2 ../data/Run202/STAP_TSS_27ac_rep1_S3_R2_001.fastq.gz \
  --taps-r2 ../data/Run202/TAPS_27ac_rep1_S7_R2_001.fastq.gz \
  --outdir results/raw_barcode_overlap/STAP_TSS_27ac_rep1_S3__TAPS_27ac_rep1_S7