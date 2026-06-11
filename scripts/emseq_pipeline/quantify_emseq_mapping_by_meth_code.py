#!/usr/bin/env python3
"""Quantify Bismark unique alignment rate by EM-seq methylation code."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter
from pathlib import Path
from typing import BinaryIO

import pysam


METHYLATION = {
    "TT": ("0_pct", 0.0),
    "AA": ("100_pct", 1.0),
    "GG": ("40_pct", 0.4),
    "CC": ("10_pct", 0.1),
    "AT": ("1_pct", 0.01),
    "TA": ("0.1_pct", 0.001),
}

METH_RE = re.compile(r"(?:^|\|)METH_CODE=([^|/]+)")
METH_RE_BYTES = re.compile(rb"(?:^|\|)METH_CODE=([^|/\s]+)")


def open_binary(path: Path) -> BinaryIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")  # type: ignore[return-value]
    return path.open("rb")


def parse_meth(query_name: str) -> str:
    match = METH_RE.search(query_name)
    return match.group(1) if match else "UNKNOWN"


def parse_meth_header(header: bytes) -> str:
    match = METH_RE_BYTES.search(header)
    if not match:
        return "UNKNOWN"
    return match.group(1).decode("ascii", "replace")


def count_trimmed_fastq_codes(path: Path, progress_every: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    with open_binary(path) as handle:
        record_idx = 0
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ record in {path}")
            record_idx += 1
            counts[parse_meth_header(header.rstrip(b"\r\n"))] += 1
            if progress_every and record_idx % progress_every == 0:
                print(f"counted {record_idx:,} FASTQ records", file=sys.stderr, flush=True)
    return counts


def count_bam_codes(path: Path, min_mapq: int, progress_every: int) -> tuple[Counter[str], Counter[str]]:
    aligned_pairs: Counter[str] = Counter()
    record_stats: Counter[str] = Counter()
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            record_stats["bam_records"] += 1
            if progress_every and record_stats["bam_records"] % progress_every == 0:
                print(f"processed {record_stats['bam_records']:,} BAM records", file=sys.stderr, flush=True)
            if read.is_secondary or read.is_supplementary:
                record_stats["skip_secondary_or_supplementary"] += 1
                continue
            if read.is_read2:
                record_stats["skip_read2"] += 1
                continue
            code = parse_meth(read.query_name)
            if read.is_unmapped:
                record_stats["read1_unmapped_in_bam"] += 1
                continue
            if read.is_qcfail or read.mapping_quality < min_mapq:
                record_stats["read1_filtered"] += 1
                continue
            aligned_pairs[code] += 1
            record_stats["read1_aligned_counted"] += 1
    return aligned_pairs, record_stats


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--trimmed-r1", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()

    input_counts = count_trimmed_fastq_codes(args.trimmed_r1, args.progress_every)
    aligned_counts, record_stats = count_bam_codes(args.bam, args.min_mapq, args.progress_every)

    codes = sorted(set(input_counts) | set(aligned_counts) | set(METHYLATION))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_methylation\ttrimmed_input_pairs\t"
            "unique_aligned_pairs\tmissing_or_nonunique_pairs\tunique_alignment_rate\n"
        )
        for code in codes:
            input_count = input_counts[code]
            aligned = aligned_counts[code]
            missing = input_count - aligned
            rate = aligned / input_count if input_count else float("nan")
            label, expected = METHYLATION.get(code, ("Unknown", float("nan")))
            out.write(
                f"{code}\t{label}\t{expected:g}\t{input_count}\t{aligned}\t"
                f"{missing}\t{rate:.6g}\n"
            )
        out.write("\nmetric\tcount\n")
        out.write(f"bam_path\t{args.bam}\n")
        out.write(f"trimmed_r1\t{args.trimmed_r1}\n")
        out.write(f"min_mapq\t{args.min_mapq}\n")
        for key, value in sorted(record_stats.items()):
            out.write(f"{key}\t{value}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
