#!/usr/bin/env python3
"""Copy UMI-tagged FASTQ headers from a template FASTQ onto another FASTQ."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Iterator, TextIO


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
        yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def first_token(header: str) -> str:
    return header[1:].split()[0] if header.startswith("@") else header.split()[0]


def strip_appended_umi(read_id: str) -> str:
    return read_id.rsplit("_", 1)[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--sequences", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--keep-sequence-comment",
        action="store_true",
        help="Keep the text after the first whitespace from the sequence FASTQ header.",
    )
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(args.template, "r") as template_fh, open_text(
        args.sequences, "r"
    ) as seq_fh, open_text(args.output, "w") as out_fh:
        template_iter = fastq_records(template_fh)
        seq_iter = fastq_records(seq_fh)
        while True:
            try:
                template = next(template_iter)
            except StopIteration:
                try:
                    next(seq_iter)
                except StopIteration:
                    break
                raise ValueError("Sequence FASTQ has more records than template FASTQ")

            try:
                sequence = next(seq_iter)
            except StopIteration as exc:
                raise ValueError("Template FASTQ has more records than sequence FASTQ") from exc

            template_id = strip_appended_umi(first_token(template[0]))
            sequence_id = first_token(sequence[0])
            if template_id != sequence_id:
                raise ValueError(
                    f"Read order mismatch at record {count + 1}: "
                    f"{template_id!r} != {sequence_id!r}"
                )

            header = template[0]
            if args.keep_sequence_comment and " " in sequence[0]:
                header = f"@{first_token(template[0])} {sequence[0].split(' ', 1)[1]}"

            out_fh.write(f"{header}\n{sequence[1]}\n+\n{sequence[3]}\n")
            count += 1

    print(f"wrote {count} records to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
