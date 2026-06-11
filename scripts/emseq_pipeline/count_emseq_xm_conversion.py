#!/usr/bin/env python3
"""Count EM-seq methylation/conversion calls by code from Bismark XM tags."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

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
XM_CONTEXT = {
    "Z": "CpG",
    "z": "CpG",
    "X": "CHG",
    "x": "CHG",
    "H": "CHH",
    "h": "CHH",
    "U": "Unknown",
    "u": "Unknown",
}


def parse_meth(query_name: str) -> str:
    match = METH_RE.search(query_name)
    return match.group(1) if match else "UNKNOWN"


def is_usable(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or read.mapping_quality < min_mapq
    )


def add_call(
    counts: Counter[tuple[str, str, str]],
    code: str,
    context: str,
    methylated: bool,
) -> None:
    outcome = "methylated" if methylated else "converted"
    counts[(code, context, outcome)] += 1
    if context in {"CHG", "CHH"}:
        counts[(code, "CpH", outcome)] += 1


def write_summary(
    path: Path,
    counts: Counter[tuple[str, str, str]],
    contexts: list[str],
    read_stats: Counter[str],
    args: argparse.Namespace,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    codes = sorted({key[0] for key in counts} | set(METHYLATION))
    with path.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_methylation\tcontext\t"
            "methylated_count\tconverted_count\tcallable_count\t"
            "methylation_rate\tconversion_rate\n"
        )
        for code in codes:
            label, expected = METHYLATION.get(code, ("Unknown", float("nan")))
            for context in contexts:
                methylated = counts[(code, context, "methylated")]
                converted = counts[(code, context, "converted")]
                callable_count = methylated + converted
                if callable_count == 0:
                    continue
                methylation_rate = methylated / callable_count
                conversion_rate = converted / callable_count
                out.write(
                    f"{code}\t{label}\t{expected:g}\t{context}\t"
                    f"{methylated}\t{converted}\t{callable_count}\t"
                    f"{methylation_rate:.6g}\t{conversion_rate:.6g}\n"
                )
        out.write("\nmetric\tcount\n")
        out.write(f"bam_path\t{args.bam}\n")
        out.write(f"min_mapq\t{args.min_mapq}\n")
        out.write(f"count_read1\t{int(args.count_read1)}\n")
        out.write(f"count_read2\t{int(args.count_read2)}\n")
        out.write(
            "note\tRates are counted from Bismark XM tags. methylation_rate is retained/unconverted C; "
            "conversion_rate is converted/unmethylated C.\n"
        )
        for key, value in sorted(read_stats.items()):
            out.write(f"{key}\t{value}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--out-cpg-summary", required=True, type=Path)
    parser.add_argument("--out-noncpg-summary", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    parser.add_argument("--count-read1", action="store_true", default=True)
    parser.add_argument("--count-read2", action="store_true", default=True)
    parser.add_argument(
        "--include-cpg-in-noncpg-summary",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include CpG rows in the non-CpG summary as a baseline.",
    )
    args = parser.parse_args()

    counts: Counter[tuple[str, str, str]] = Counter()
    read_stats: Counter[str] = Counter()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            read_stats["bam_records"] += 1
            if args.progress_every and read_stats["bam_records"] % args.progress_every == 0:
                print(f"processed {read_stats['bam_records']:,} BAM records", file=sys.stderr, flush=True)
            if not is_usable(read, args.min_mapq):
                read_stats["skipped_records"] += 1
                continue
            if read.is_read1 and not args.count_read1:
                read_stats["skipped_read1"] += 1
                continue
            if read.is_read2 and not args.count_read2:
                read_stats["skipped_read2"] += 1
                continue
            if not read.has_tag("XM"):
                read_stats["skipped_no_xm"] += 1
                continue

            code = parse_meth(read.query_name)
            xm = read.get_tag("XM")
            for char in xm:
                context = XM_CONTEXT.get(char)
                if context is None:
                    continue
                add_call(counts, code, context, char.isupper())
                read_stats[f"{context}_calls"] += 1
            read_stats["usable_records"] += 1

    write_summary(args.out_cpg_summary, counts, ["CpG"], read_stats, args)
    noncpg_contexts = ["CpH", "CHG", "CHH"]
    if args.include_cpg_in_noncpg_summary:
        noncpg_contexts.append("CpG")
    write_summary(args.out_noncpg_summary, counts, noncpg_contexts, read_stats, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
