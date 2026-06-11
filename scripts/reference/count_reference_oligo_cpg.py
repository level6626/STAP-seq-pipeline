#!/usr/bin/env python3
"""Count CpG dinucleotides in each reference oligo sequence."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


DEFAULT_REFERENCE_FASTA = Path(
    "/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa"
)
DEFAULT_OUTDIR = Path("data/meta")


def read_fasta(path: Path):
    name = None
    seq_parts: list[str] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"Missing FASTA record name at {path}:{line_number}")
                seq_parts = []
            else:
                if name is None:
                    raise ValueError(f"Found sequence before first FASTA header at {path}:{line_number}")
                seq_parts.append(line)
    if name is not None:
        yield name, "".join(seq_parts).upper()


def count_cpg(seq: str) -> int:
    return sum(1 for i in range(len(seq) - 1) if seq[i : i + 2] == "CG")


def default_output_path(reference_fasta: Path, outdir: Path) -> Path:
    stem = reference_fasta.name
    if stem.endswith(".fa"):
        stem = stem[:-3]
    elif stem.endswith(".fasta"):
        stem = stem[:-6]
    else:
        stem = reference_fasta.stem
    return outdir / f"{stem}_cpg_counts.tsv"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-fasta",
        type=Path,
        default=DEFAULT_REFERENCE_FASTA,
        help=f"Reference oligo FASTA. Default: {DEFAULT_REFERENCE_FASTA}",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Output folder used when --out is not set. Default: {DEFAULT_OUTDIR}",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output TSV path. Default: <outdir>/<reference>_cpg_counts.tsv",
    )
    args = parser.parse_args()

    if not args.reference_fasta.exists():
        raise SystemExit(f"Missing reference FASTA: {args.reference_fasta}")

    out_path = args.out if args.out is not None else default_output_path(args.reference_fasta, args.outdir)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fields = ["sequence_id", "length_bp", "cpg_count", "cpg_per_100bp"]
    records_written = 0
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for sequence_id, seq in read_fasta(args.reference_fasta):
            cpg_count = count_cpg(seq)
            length = len(seq)
            writer.writerow(
                {
                    "sequence_id": sequence_id,
                    "length_bp": length,
                    "cpg_count": cpg_count,
                    "cpg_per_100bp": f"{(cpg_count / length * 100) if length else 0:.6g}",
                }
            )
            records_written += 1

    print(f"wrote CpG counts for {records_written} sequences to {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
