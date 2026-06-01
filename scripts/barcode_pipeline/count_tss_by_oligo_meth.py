#!/usr/bin/env python3
"""Count deduplicated read1 TSS positions grouped by OLIGO and METH tags."""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

import pysam


METH_RE = re.compile(r"(?:^|\|)METH=([^|]+)")
OLIGO_UMI_RE = re.compile(r"\|OLIGO=(.+)_([ACGTNacgtn]{22})$")


def parse_tags(query_name: str) -> tuple[str, str]:
    meth_match = METH_RE.search(query_name)
    oligo_match = OLIGO_UMI_RE.search(query_name)
    if not meth_match or not oligo_match:
        return "UNKNOWN", "UNKNOWN"
    return oligo_match.group(1), meth_match.group(1)


def read1_tss_1based(read: pysam.AlignedSegment) -> int:
    if read.is_reverse:
        return int(read.reference_end)
    return int(read.reference_start) + 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    args = parser.parse_args()

    counts: Counter[tuple[str, str, str, int]] = Counter()
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_read1:
                continue
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_qcfail or read.is_duplicate:
                continue
            if read.mapping_quality < args.min_mapq:
                continue
            oligo, meth = parse_tags(read.query_name)
            chrom = bam.get_reference_name(read.reference_id)
            counts[(oligo, meth, chrom, read1_tss_1based(read))] += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        out.write("Oligo_ID\tMeth_State\tChromosome\tTSS_Position\tCount\n")
        for (oligo, meth, chrom, pos), count in sorted(counts.items()):
            out.write(f"{oligo}\t{meth}\t{chrom}\t{pos}\t{count}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
