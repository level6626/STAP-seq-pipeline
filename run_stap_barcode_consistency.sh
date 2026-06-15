cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

python \
  scripts/standard_tools/quantify_stap_barcode_mapping_consistency.py \
  --bam results/standard_tools/STAP_TSS_27ac_rep1_S3_bowtie/STAP_TSS_27ac_rep1_S3.bowtie2.dedup.bam \
  --out-family results/standard_tools/STAP_TSS_27ac_rep1_S3_bowtie/STAP_TSS_27ac_rep1_S3.barcode_mapping_consistency.mapq40.w50.families.tsv \
  --out-summary results/standard_tools/STAP_TSS_27ac_rep1_S3_bowtie/STAP_TSS_27ac_rep1_S3.barcode_mapping_consistency.mapq40.w50.summary.tsv \
  --window-bp 50 \
  --min-mapq 40 \
  --require-known-meth-code