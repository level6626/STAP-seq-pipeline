#!/usr/bin/env python3
"""Write a FASTA with right-side padding on every record."""

from __future__ import annotations

import argparse
from pathlib import Path


def fasta_records(path: Path):
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def wrap(seq: str, width: int):
    for start in range(0, len(seq), width):
        yield seq[start : start + width]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--right-pad", type=int, default=100)
    parser.add_argument("--pad-base", default="N")
    parser.add_argument("--line-width", type=int, default=80)
    args = parser.parse_args()

    if args.right_pad < 0:
        raise ValueError("--right-pad must be non-negative")
    pad_base = args.pad_base.upper()
    if len(pad_base) != 1 or pad_base not in {"A", "C", "G", "T", "N"}:
        raise ValueError("--pad-base must be one of A/C/G/T/N")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pad = pad_base * args.right_pad
    with args.out.open("w") as out:
        for name, seq in fasta_records(args.input):
            out.write(f">{name}\n")
            for line in wrap(seq + pad, args.line_width):
                out.write(f"{line}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
