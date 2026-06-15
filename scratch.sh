READ='VH01570:206:AAJ2NL3M5:1:2102:62332:39905|METH_CODE=AAA|METH_LABEL=0%|METH_EXPECTED=0|R2=AAATGAAGGGAGGGAAT|UMI=CTCACAACAAATGAAGGGAGGGAAT'
CORE="${READ%%|*}"

# for fq in results/taps_pipeline/TAPS_27ac_rep1_S7_paired_bowtie2/work/TAPS_27ac_rep1_S7.tagged*.fastq.gz
# do
#   echo "$fq"
#   zcat "$fq" | awk -v id="@$READ" '
#     $0==id {print; getline; print; getline; print; getline; print; exit}
#   '
# done

for fq in /gpfs/data/zhou-lab/yczhang/methylation/data/Run202/TAPS_27ac_rep1_S7_R{1,2,3}_001.fastq.gz
do
  echo "$fq"
  zcat "$fq" | awk -v id="@$CORE" '
    index($0,id)==1 {print; getline; print; getline; print; getline; print; exit}
  '
done