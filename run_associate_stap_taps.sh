python scripts/associate_stap_taps_reads.py \
  --stap-bam results/standard_tools/STAP_TSS_BasPromoter_S6_bowtie/STAP_TSS_BasPromoter_S6.bowtie2.dedup.bam \
  --taps-bam results/taps_pipeline/TAPS_BasPromoter_S10_paired_bowtie2/TAPS_BasPromoter_S10.bowtie2.paired.sorted.bam \
  --outdir results/stap_taps_association/STAP_TSS_BasPromoter_S6_bowtie_TAPS_BasPromoter_S10_paired_bowtie2 \
  --min-mapq 20