#!/usr/bin/env python3
"""Prepare TAPS FASTQs by carrying R2 methylation/barcode tags into read names."""

from __future__ import annotations

import argparse
import gzip
import re
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


def open_text(path: Path, mode: str) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t")  # type: ignore[return-value]
    return path.open(mode)


def fastq_records(handle: TextIO) -> Iterator[tuple[str, str, str, str]]:
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("Truncated FASTQ record")
        yield header.rstrip("\n"), seq.rstrip("\n").upper(), plus.rstrip("\n"), qual.rstrip("\n")


def read_token(header: str) -> str:
    return header[1:].split()[0] if header.startswith("@") else header.split()[0]


def normalize_read_id(header: str) -> str:
    token = read_token(header)
    return re.sub(r"([/._ -][123])$", "", token)


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


def orient_r2(seq: str, orientation: str) -> tuple[str | None, str]:
    candidates: list[tuple[str, str]] = []
    if orientation in {"forward", "both"}:
        candidates.append((seq.upper(), "forward"))
    if orientation in {"reverse-complement", "both"}:
        candidates.append((reverse_complement(seq), "reverse-complement"))

    hits = [(oriented, label) for oriented, label in candidates if oriented[:3] in METHYLATION]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        labels = {label for _oriented, label in hits}
        if len(labels) == 1:
            return hits[0]
        return None, "ambiguous"
    return None, "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--r3", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r3", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--r1-umi-length", type=int, default=8)
    parser.add_argument("--r2-length", type=int, default=17)
    parser.add_argument(
        "--r2-orientation",
        choices=["forward", "reverse-complement", "both"],
        default="forward",
        help="Orientation used before parsing R2[:3] as the methylation code.",
    )
    parser.add_argument(
        "--filter-r2-motif",
        default="",
        help="Discard triplets whose raw R2 contains this motif. Empty disables filtering.",
    )
    parser.add_argument(
        "--trim-r1-umi",
        action="store_true",
        help="Remove the leading R1 UMI from the emitted R1 sequence.",
    )
    args = parser.parse_args()

    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r3.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)

    stats: Counter[str] = Counter()
    with open_text(args.r1, "r") as r1_fh, open_text(args.r2, "r") as r2_fh, open_text(
        args.r3, "r"
    ) as r3_fh, open_text(args.out_r1, "w") as out_r1, open_text(args.out_r3, "w") as out_r3:
        for idx, (rec1, rec2, rec3) in enumerate(
            zip(fastq_records(r1_fh), fastq_records(r2_fh), fastq_records(r3_fh)), start=1
        ):
            if args.max_reads and idx > args.max_reads:
                break
            stats["total_triplets"] += 1

            h1, s1, _p1, q1 = rec1
            h2, s2, _p2, _q2 = rec2
            h3, s3, _p3, q3 = rec3

            if normalize_read_id(h1) != normalize_read_id(h2) or normalize_read_id(h1) != normalize_read_id(h3):
                stats["discard_read_name_mismatch"] += 1
                continue
            if len(s1) < args.r1_umi_length or len(q1) < args.r1_umi_length:
                stats["discard_r1_too_short"] += 1
                continue
            if len(s2) < args.r2_length:
                stats["discard_r2_too_short"] += 1
                continue
            if args.filter_r2_motif and args.filter_r2_motif.upper() in s2:
                stats["discard_r2_filter_motif"] += 1
                continue

            oriented_r2, orientation = orient_r2(s2[: args.r2_length], args.r2_orientation)
            if oriented_r2 is None:
                stats[f"discard_r2_{orientation}_meth_code"] += 1
                continue

            meth_code = oriented_r2[:3]
            meth_label, meth_expected = METHYLATION[meth_code]
            r1_umi = s1[: args.r1_umi_length]
            molecule_umi = f"{r1_umi}{oriented_r2}"
            original = read_token(h1)
            new_header = (
                f"@{original}|METH_CODE={meth_code}|METH_LABEL={meth_label}|"
                f"METH_EXPECTED={meth_expected:g}|R2={oriented_r2}|UMI={molecule_umi}"
            )

            out_s1 = s1[args.r1_umi_length :] if args.trim_r1_umi else s1
            out_q1 = q1[args.r1_umi_length :] if args.trim_r1_umi else q1
            if not out_s1 or not s3:
                stats["discard_empty_after_trim"] += 1
                continue

            out_r1.write(f"{new_header}\n{out_s1}\n+\n{out_q1}\n")
            out_r3.write(f"{new_header}\n{s3}\n+\n{q3}\n")
            stats["written_triplets"] += 1
            stats[f"meth_code_{meth_code}_{meth_label}"] += 1
            stats[f"r2_orientation_{orientation}"] += 1

            if idx % 1_000_000 == 0:
                print(
                    f"processed={idx} written={stats['written_triplets']} "
                    f"unknown={stats['discard_r2_unknown_meth_code']}",
                    file=sys.stderr,
                    flush=True,
                )

    with args.stats.open("w") as out:
        out.write("metric\tcount\n")
        for key, value in sorted(stats.items()):
            out.write(f"{key}\t{value}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
