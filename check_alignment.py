from pathlib import Path
import pysam
import re
from collections import Counter

prep = Path("results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.prepare_taps.stats.tsv")
bam_path = "results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.star.r3.sorted.bam"

meth_re = re.compile(r"(?:^|\|)METH_CODE=([^|]+)")

input_counts = {}
for line in prep.read_text().splitlines()[1:]:
    metric, count = line.rstrip("\n").split("\t")
    if metric.startswith("meth_code_"):
        # e.g. meth_code_TTT_100%
        code = metric.split("_")[2]
        input_counts[code] = int(count)

mapped = Counter()
nm_sum = Counter()
aligned_bases = Counter()

with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue

        m = meth_re.search(read.query_name)
        if not m:
            continue
        code = m.group(1)
        if code not in {"AAA", "TTT"}:
            continue

        mapped[code] += 1

        try:
            nm = read.get_tag("nM")  # STAR
        except KeyError:
            try:
                nm = read.get_tag("NM")
            except KeyError:
                nm = None

        if nm is not None:
            nm_sum[code] += nm
            aligned_bases[code] += read.query_alignment_length or 0

print("code\tinput_tagged\tmapped\tmapped_fraction\tavg_NM\tNM_per_100bp")
for code in ["AAA", "TTT"]:
    inp = input_counts.get(code, 0)
    mp = mapped.get(code, 0)
    frac = mp / inp if inp else float("nan")
    avg_nm = nm_sum[code] / mp if mp else float("nan")
    nm100 = 100 * nm_sum[code] / aligned_bases[code] if aligned_bases[code] else float("nan")
    print(f"{code}\t{inp}\t{mp}\t{frac:.6f}\t{avg_nm:.6f}\t{nm100:.6f}")