#!/usr/bin/env python3
"""Prepare custom 3-read EM-seq FASTQs for Bismark and UMI-tools."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter
from pathlib import Path
from typing import BinaryIO, Iterator


METHYLATION = {
    b"TT": "0_pct",
    b"AA": "100_pct",
    b"GG": "40_pct",
    b"CC": "10_pct",
    b"AT": "1_pct",
    b"TA": "0.1_pct",
}

RC_TABLE = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def open_binary(path: Path, mode: str) -> BinaryIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return path.open(mode)


def fastq_records(handle: BinaryIO) -> Iterator[tuple[bytes, bytes, bytes, bytes]]:
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("Truncated FASTQ record")
        yield (
            header.rstrip(b"\r\n"),
            seq.rstrip(b"\r\n").upper(),
            plus.rstrip(b"\r\n"),
            qual.rstrip(b"\r\n"),
        )


def read_token(header: bytes) -> bytes:
    token = header[1:].split()[0] if header.startswith(b"@") else header.split()[0]
    return token


def normalize_read_id(header: bytes) -> bytes:
    token = read_token(header)
    return re.sub(rb"([/._ -][123])$", b"", token)


def reverse_complement(seq: bytes) -> bytes:
    return seq.translate(RC_TABLE)[::-1].upper()


def orient_r2(seq: bytes, orientation: str) -> tuple[bytes, str]:
    if orientation == "forward":
        return seq.upper(), "forward"
    if orientation == "reverse-complement":
        return reverse_complement(seq), "reverse-complement"

    candidates = [(seq.upper(), "forward"), (reverse_complement(seq), "reverse-complement")]
    hits = [(oriented, label) for oriented, label in candidates if oriented[:2] in METHYLATION]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1 and hits[0][0] == hits[1][0]:
        return hits[0]
    if len(hits) > 1:
        return seq.upper(), "ambiguous"
    return seq.upper(), "unknown"


def sanitize_tag(value: str) -> str:
    value = re.sub(r"\s+", "-", value.strip())
    value = re.sub(r"[^A-Za-z0-9_.:+-]", "-", value)
    return value.strip("-") or "NA"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stream synchronized R1/R2/R3 EM-seq FASTQs, trim the leading R1 UMI, "
            "and carry R2 methylation/plasmid-barcode metadata into paired FASTQ "
            "read names for Bismark and UMI-tools."
        )
    )
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--r3", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r3", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--r1-umi-length", type=int, default=8)
    parser.add_argument("--r2-meth-length", type=int, default=2)
    parser.add_argument("--expected-r1-length", type=int, default=55)
    parser.add_argument(
        "--r2-orientation",
        choices=["forward", "reverse-complement", "both"],
        default="forward",
        help="Orientation used before parsing R2[:2] as the methylation code.",
    )
    parser.add_argument(
        "--keep-r1-umi",
        action="store_true",
        help="Keep the leading R1 UMI in emitted R1 instead of trimming it.",
    )
    parser.add_argument(
        "--umi-source",
        choices=["r1", "r1+r2", "r1+r2_barcode"],
        default="r1+r2_barcode",
        help=(
            "Final read-name UMI suffix used by umi_tools. r1+r2_barcode is the "
            "default molecule key: R1 UMI plus R2 bases 3-end."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r3.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)

    stats: Counter[str] = Counter()
    with open_binary(args.r1, "rb") as r1_fh, open_binary(args.r2, "rb") as r2_fh, open_binary(
        args.r3, "rb"
    ) as r3_fh, open_binary(args.out_r1, "wb") as out_r1, open_binary(args.out_r3, "wb") as out_r3:
        triplets = zip(fastq_records(r1_fh), fastq_records(r2_fh), fastq_records(r3_fh))
        for idx, (rec1, rec2, rec3) in enumerate(triplets, start=1):
            if args.max_reads and idx > args.max_reads:
                break
            stats["total_triplets"] += 1

            h1, s1, _p1, q1 = rec1
            h2, s2, _p2, _q2 = rec2
            h3, s3, _p3, q3 = rec3

            if normalize_read_id(h1) != normalize_read_id(h2) or normalize_read_id(h1) != normalize_read_id(h3):
                stats["discard_read_name_mismatch"] += 1
                continue
            if len(s1) != len(q1) or len(s3) != len(q3):
                stats["discard_sequence_quality_length_mismatch"] += 1
                continue
            if len(s1) < args.r1_umi_length + 1 or len(q1) < args.r1_umi_length + 1:
                stats["discard_r1_too_short"] += 1
                continue
            if len(s2) < args.r2_meth_length:
                stats["discard_r2_too_short"] += 1
                continue

            stats[f"r1_length_{len(s1)}"] += 1
            stats[f"r2_length_{len(s2)}"] += 1
            stats[f"r3_length_{len(s3)}"] += 1
            if args.expected_r1_length and len(s1) != args.expected_r1_length:
                stats[f"r1_length_not_expected_{args.expected_r1_length}"] += 1

            r1_umi = s1[: args.r1_umi_length]
            oriented_r2, orientation = orient_r2(s2, args.r2_orientation)
            meth_code = oriented_r2[: args.r2_meth_length]
            plasmid_barcode = oriented_r2[args.r2_meth_length :]
            meth_label = METHYLATION.get(meth_code, "Unknown")

            if args.umi_source == "r1":
                molecule_umi = r1_umi
            elif args.umi_source == "r1+r2":
                molecule_umi = r1_umi + oriented_r2
            else:
                molecule_umi = r1_umi + plasmid_barcode

            if not args.keep_r1_umi:
                out_s1 = s1[args.r1_umi_length :]
                out_q1 = q1[args.r1_umi_length :]
            else:
                out_s1 = s1
                out_q1 = q1
            if not out_s1 or not s3:
                stats["discard_empty_after_trim"] += 1
                continue

            original = read_token(h1).decode("ascii", "replace")
            tag = (
                f"{original}|METH_CODE={meth_code.decode('ascii', 'replace')}|"
                f"METH_LABEL={sanitize_tag(meth_label)}|"
                f"R2={oriented_r2.decode('ascii', 'replace')}|"
                f"R2_BARCODE={plasmid_barcode.decode('ascii', 'replace')}|"
                f"UMI={molecule_umi.decode('ascii', 'replace')}"
            ).encode()
            new_header = b"@" + tag + b"|UMITOOLS_" + molecule_umi

            out_r1.write(new_header + b"\n" + out_s1 + b"\n+\n" + out_q1 + b"\n")
            out_r3.write(new_header + b"\n" + s3 + b"\n+\n" + q3 + b"\n")

            stats["written_triplets"] += 1
            stats[f"meth_code_{meth_code.decode('ascii', 'replace')}_{sanitize_tag(meth_label)}"] += 1
            stats[f"r2_orientation_{orientation}"] += 1
            if meth_label == "Unknown":
                stats["meth_label_Unknown"] += 1

            if idx % 1_000_000 == 0:
                print(
                    f"processed={idx} written={stats['written_triplets']} "
                    f"unknown={stats['meth_label_Unknown']}",
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
