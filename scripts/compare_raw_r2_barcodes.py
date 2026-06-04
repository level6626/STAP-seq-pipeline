#!/usr/bin/env python3
"""Compare plasmid R2 barcode overlap directly from raw FASTQ/FASTA files."""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path
from typing import Iterator, TextIO


METHYLATION = {
    "TTT": ("100%", 1.0),
    "AAA": ("0%", 0.0),
    "CAT": ("60%", 0.6),
    "AGT": ("40%", 0.4),
    "TGA": ("20%", 0.2),
    "TAG": ("10%", 0.1),
    "CTA": ("1%", 0.01),
    "ATG": ("0.1%", 0.001),
}

RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")  # type: ignore[return-value]
    return path.open("r")


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def fasta_records(handle: TextIO, first_header: str) -> Iterator[tuple[str, str]]:
    header = first_header.rstrip("\n")
    seq_parts: list[str] = []
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            yield header, "".join(seq_parts).upper()
            header = line
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    yield header, "".join(seq_parts).upper()


def fastq_records(handle: TextIO, first_header: str) -> Iterator[tuple[str, str]]:
    header = first_header.rstrip("\n")
    while header:
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("Truncated FASTQ record")
        yield header, seq.rstrip("\n").upper()
        header = handle.readline().rstrip("\n")


def sequence_records(path: Path) -> Iterator[tuple[str, str]]:
    with open_text(path) as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith("@"):
            yield from fastq_records(handle, first)
        elif first.startswith(">"):
            yield from fasta_records(handle, first)
        else:
            raise ValueError(f"Could not detect FASTQ/FASTA format for {path}")


def orient_r2(seq: str, orientation: str) -> tuple[str, str]:
    seq = seq.upper()
    if orientation == "forward":
        return seq, "forward"
    if orientation == "reverse-complement":
        return reverse_complement(seq), "reverse-complement"

    rc = reverse_complement(seq)
    forward_known = seq[:3] in METHYLATION
    rc_known = rc[:3] in METHYLATION
    if forward_known and not rc_known:
        return seq, "forward"
    if rc_known and not forward_known:
        return rc, "reverse-complement"
    return seq, "forward"


def count_raw_barcodes(
    path: Path,
    orientation: str,
    r2_length: int,
    max_records: int,
    require_known_code: bool,
) -> tuple[Counter[str], Counter[str]]:
    counts: Counter[str] = Counter()
    stats: Counter[str] = Counter()

    for _header, seq in sequence_records(path):
        if max_records and stats["records_seen"] >= max_records:
            stats["stopped_at_max_records"] = max_records
            break
        stats["records_seen"] += 1
        if len(seq) < r2_length:
            stats["skipped_too_short"] += 1
            continue
        raw_r2 = seq[:r2_length]
        r2, used_orientation = orient_r2(raw_r2, orientation)
        stats[f"orientation_{used_orientation}"] += 1
        if require_known_code and r2[:3] not in METHYLATION:
            stats["skipped_unknown_meth_code"] += 1
            continue
        counts[r2] += 1
        stats["usable_records"] += 1

        if stats["records_seen"] % 5_000_000 == 0:
            print(
                f"{path.name}: records={stats['records_seen']} unique_barcodes={len(counts)}",
                file=sys.stderr,
                flush=True,
            )

    stats["unique_barcodes"] = len(counts)
    return counts, stats


def meth_info(r2: str) -> tuple[str, str]:
    meth_code = r2[:3]
    meth_label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
    return meth_label, f"{expected:g}"


def write_counts(path: Path, assay: str, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write(
            f"r2_barcode\tmeth_code\tmeth_label\texpected_methylation\trandom14\t"
            f"{assay}_raw_reads\n"
        )
        for r2 in sorted(counts):
            meth_label, expected = meth_info(r2)
            out.write(f"{r2}\t{r2[:3]}\t{meth_label}\t{expected}\t{r2[3:]}\t{counts[r2]}\n")


def write_overlap(path: Path, stap_counts: Counter[str], taps_counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write(
            "r2_barcode\tmeth_code\tmeth_label\texpected_methylation\trandom14\t"
            "in_stap_raw\tin_taps_raw\tstap_raw_reads\ttaps_raw_reads\n"
        )
        for r2 in sorted(set(stap_counts) | set(taps_counts)):
            meth_label, expected = meth_info(r2)
            out.write(
                f"{r2}\t{r2[:3]}\t{meth_label}\t{expected}\t{r2[3:]}\t"
                f"{int(stap_counts[r2] > 0)}\t{int(taps_counts[r2] > 0)}\t"
                f"{stap_counts[r2]}\t{taps_counts[r2]}\n"
            )


def write_summary(
    path: Path,
    stap_counts: Counter[str],
    taps_counts: Counter[str],
    stap_stats: Counter[str],
    taps_stats: Counter[str],
    min_stap_reads: int,
    min_taps_reads: int,
) -> None:
    all_barcodes = set(stap_counts) | set(taps_counts)
    meth_codes = sorted({r2[:3] for r2 in all_barcodes} | set(METHYLATION))

    def stap_present(r2: str) -> bool:
        return stap_counts[r2] >= min_stap_reads

    def taps_present(r2: str) -> bool:
        return taps_counts[r2] >= min_taps_reads

    with path.open("w") as out:
        out.write(
            "scope\tmeth_code\tmeth_label\texpected_methylation\t"
            "stap_raw_barcodes\ttaps_raw_barcodes\toverlap_raw_barcodes\t"
            "fraction_stap_raw_in_taps_raw\tfraction_taps_raw_in_stap_raw\t"
            "stap_raw_reads\tstap_raw_reads_in_overlap_barcodes\t"
            "taps_raw_reads\ttaps_raw_reads_in_overlap_barcodes\n"
        )
        for scope, code_filter in [("all", None), *[("meth_code", code) for code in meth_codes]]:
            if code_filter is None:
                barcodes = all_barcodes
                meth_code = "ALL"
                meth_label = "ALL"
                expected = "nan"
            else:
                barcodes = {r2 for r2 in all_barcodes if r2[:3] == code_filter}
                meth_code = code_filter
                meth_label, expected = meth_info(code_filter + "N" * 14)

            stap_set = {r2 for r2 in barcodes if stap_present(r2)}
            taps_set = {r2 for r2 in barcodes if taps_present(r2)}
            overlap = stap_set & taps_set
            frac_stap = len(overlap) / len(stap_set) if stap_set else float("nan")
            frac_taps = len(overlap) / len(taps_set) if taps_set else float("nan")
            stap_reads = sum(stap_counts[r2] for r2 in stap_set)
            taps_reads = sum(taps_counts[r2] for r2 in taps_set)
            stap_overlap_reads = sum(stap_counts[r2] for r2 in overlap)
            taps_overlap_reads = sum(taps_counts[r2] for r2 in overlap)
            out.write(
                f"{scope}\t{meth_code}\t{meth_label}\t{expected}\t"
                f"{len(stap_set)}\t{len(taps_set)}\t{len(overlap)}\t"
                f"{frac_stap:.6g}\t{frac_taps:.6g}\t"
                f"{stap_reads}\t{stap_overlap_reads}\t{taps_reads}\t{taps_overlap_reads}\n"
            )

        out.write("\nmetric\tassay\tcount\n")
        for assay, stats in [("stap_raw", stap_stats), ("taps_raw", taps_stats)]:
            for key, value in sorted(stats.items()):
                out.write(f"{key}\t{assay}\t{value}\n")
        out.write(f"min_stap_reads\tparameter\t{min_stap_reads}\n")
        out.write(f"min_taps_reads\tparameter\t{min_taps_reads}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stap-r2", required=True, type=Path)
    parser.add_argument("--taps-r2", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--r2-length", type=int, default=17)
    parser.add_argument(
        "--stap-r2-orientation",
        choices=["forward", "reverse-complement", "auto"],
        default="forward",
    )
    parser.add_argument(
        "--taps-r2-orientation",
        choices=["forward", "reverse-complement", "auto"],
        default="forward",
    )
    parser.add_argument("--max-stap-records", type=int, default=0)
    parser.add_argument("--max-taps-records", type=int, default=0)
    parser.add_argument("--min-stap-reads", type=int, default=1)
    parser.add_argument("--min-taps-reads", type=int, default=1)
    parser.add_argument(
        "--require-known-code",
        action="store_true",
        help="Discard R2 barcodes whose first 3 bp are not one of the designed methylation codes.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    stap_counts, stap_stats = count_raw_barcodes(
        args.stap_r2,
        args.stap_r2_orientation,
        args.r2_length,
        args.max_stap_records,
        args.require_known_code,
    )
    taps_counts, taps_stats = count_raw_barcodes(
        args.taps_r2,
        args.taps_r2_orientation,
        args.r2_length,
        args.max_taps_records,
        args.require_known_code,
    )

    write_counts(args.outdir / "raw_barcode_counts_stap.tsv", "stap", stap_counts)
    write_counts(args.outdir / "raw_barcode_counts_taps.tsv", "taps", taps_counts)
    write_overlap(args.outdir / "raw_barcode_overlap.tsv", stap_counts, taps_counts)
    write_summary(
        args.outdir / "raw_barcode_overlap_summary.tsv",
        stap_counts,
        taps_counts,
        stap_stats,
        taps_stats,
        args.min_stap_reads,
        args.min_taps_reads,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
