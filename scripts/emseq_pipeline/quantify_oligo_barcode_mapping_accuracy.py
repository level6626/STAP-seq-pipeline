#!/usr/bin/env python3
"""Quantify oligo mapping accuracy using terminal oligo barcodes in R3/read2."""

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
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def parse_meth(query_name: str) -> str:
    match = METH_RE.search(query_name)
    return match.group(1) if match else "UNKNOWN"


def fasta_records(path: Path):
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def load_terminal_barcodes(path: Path, barcode_length: int) -> dict[str, str]:
    barcodes: dict[str, str] = {}
    for name, seq in fasta_records(path):
        if len(seq) < barcode_length:
            continue
        barcodes[name] = seq[-barcode_length:]
    return barcodes


def classify_read2(seq: str, barcode: str) -> tuple[str, int]:
    seq = seq.upper()
    barcode = barcode.upper()
    expected_r3_prefix = reverse_complement(barcode)
    prefix = seq[: len(barcode)]
    suffix = seq[-len(barcode) :]
    prefix_rc = reverse_complement(prefix)
    suffix_rc = reverse_complement(suffix)

    if prefix == expected_r3_prefix:
        return "expected_r3_prefix_exact", 0
    if prefix == barcode:
        return "r3_prefix_forward_exact", 0
    if suffix == expected_r3_prefix:
        return "r3_suffix_expected_exact", 0
    if suffix == barcode:
        return "r3_suffix_forward_exact", 0
    if prefix_rc == barcode:
        return "r3_prefix_rc_forward_exact", 0
    if suffix_rc == barcode:
        return "r3_suffix_rc_forward_exact", 0

    distances = [
        hamming(prefix, expected_r3_prefix),
        hamming(prefix, barcode),
        hamming(suffix, expected_r3_prefix),
        hamming(suffix, barcode),
        hamming(prefix_rc, barcode),
        hamming(suffix_rc, barcode),
    ]
    best = min(distances)
    if best == 1:
        return "one_mismatch_any_end_or_orientation", best
    return "mismatch", best


def is_usable(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or read.mapping_quality < min_mapq
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--oligo-fasta", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--barcode-length", type=int, default=10)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()

    barcodes = load_terminal_barcodes(args.oligo_fasta, args.barcode_length)
    if not barcodes:
        raise ValueError(f"No terminal barcodes loaded from {args.oligo_fasta}")

    counts: Counter[tuple[str, str]] = Counter()
    distance_sums: Counter[str] = Counter()
    stats: Counter[str] = Counter()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            stats["bam_records"] += 1
            if args.progress_every and stats["bam_records"] % args.progress_every == 0:
                print(f"processed {stats['bam_records']:,} BAM records", file=sys.stderr, flush=True)
            if not read.is_read2:
                stats["skip_not_read2"] += 1
                continue
            if not is_usable(read, args.min_mapq):
                stats["skip_unusable_read2"] += 1
                continue
            chrom = bam.get_reference_name(read.reference_id)
            if chrom is None:
                stats["skip_no_reference"] += 1
                continue
            barcode = barcodes.get(chrom)
            if barcode is None:
                stats["skip_reference_without_barcode"] += 1
                continue
            seq = read.query_sequence
            if not seq or len(seq) < args.barcode_length:
                stats["skip_short_or_missing_query"] += 1
                continue

            code = parse_meth(read.query_name)
            status, distance = classify_read2(seq, barcode)
            counts[(code, status)] += 1
            counts[(code, "total_checked")] += 1
            distance_sums[code] += distance
            stats["checked_read2"] += 1

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    statuses = [
        "expected_r3_prefix_exact",
        "r3_prefix_forward_exact",
        "r3_suffix_expected_exact",
        "r3_suffix_forward_exact",
        "r3_prefix_rc_forward_exact",
        "r3_suffix_rc_forward_exact",
        "one_mismatch_any_end_or_orientation",
        "mismatch",
    ]
    exact_statuses = set(statuses[:6])
    with args.out_summary.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_methylation\tchecked_read2\t"
            "expected_r3_prefix_exact\tany_exact\tone_mismatch_any\tmismatch\t"
            "expected_r3_prefix_exact_rate\tany_exact_rate\tone_mismatch_or_exact_rate\t"
            "mean_best_hamming_distance\n"
        )
        codes = sorted({code for code, _status in counts} | set(METHYLATION))
        for code in codes:
            total = counts[(code, "total_checked")]
            if total == 0:
                continue
            exact = sum(counts[(code, status)] for status in exact_statuses)
            one_mm = counts[(code, "one_mismatch_any_end_or_orientation")]
            mismatch = counts[(code, "mismatch")]
            expected_prefix = counts[(code, "expected_r3_prefix_exact")]
            label, expected = METHYLATION.get(code, ("Unknown", float("nan")))
            out.write(
                f"{code}\t{label}\t{expected:g}\t{total}\t{expected_prefix}\t"
                f"{exact}\t{one_mm}\t{mismatch}\t"
                f"{expected_prefix / total:.6g}\t{exact / total:.6g}\t"
                f"{(exact + one_mm) / total:.6g}\t{distance_sums[code] / total:.6g}\n"
            )

        out.write("\nstatus\tcount\n")
        for status in statuses:
            out.write(f"{status}\t{sum(counts[(code, status)] for code in codes)}\n")

        out.write("\nmetric\tcount\n")
        out.write(f"bam_path\t{args.bam}\n")
        out.write(f"oligo_fasta\t{args.oligo_fasta}\n")
        out.write(f"barcode_length\t{args.barcode_length}\n")
        out.write(f"terminal_barcodes_loaded\t{len(barcodes)}\n")
        out.write(
            "note\tExpected orientation is R3 prefix equals reverse-complement of the mapped oligo terminal barcode.\n"
        )
        for key, value in sorted(stats.items()):
            out.write(f"{key}\t{value}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
