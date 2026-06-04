#!/usr/bin/env python3
"""Build a reporter FASTA from the STAP oligo metadata workbook."""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
BARCODE_SCRIPT = SCRIPT_DIR.parent / "barcode_pipeline" / "demux_extract_stap_triplets.py"


def load_barcode_module():
    spec = importlib.util.spec_from_file_location("demux_extract_stap_triplets", BARCODE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {BARCODE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def clean_seq(seq: str) -> str:
    return re.sub(r"[^ACGTN]", "", seq.upper())


def fasta_wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--oligo-xlsx", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--include-sheets",
        default="",
        help="Comma-separated sheet names to keep. Empty keeps all sheets.",
    )
    parser.add_argument(
        "--sequence-columns",
        default=(
            "oligo sequence+barcode,oligo_sequence+barcode,oligo sequence,"
            "oligo_sequence,oligomer sequence"
        ),
    )
    parser.add_argument("--min-length", type=int, default=20)
    args = parser.parse_args()

    demux = load_barcode_module()
    records = demux.load_oligos(args.oligo_xlsx)
    keep_sheets = {x.strip() for x in args.include_sheets.split(",") if x.strip()}
    sequence_columns = [x.strip() for x in args.sequence_columns.split(",") if x.strip()]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    seen_ids: set[str] = set()
    with args.out.open("w") as out:
        for record in records:
            if keep_sheets and record.sheet not in keep_sheets:
                continue
            seq = ""
            for column in sequence_columns:
                if column in record.metadata:
                    seq = clean_seq(record.metadata[column])
                    if len(seq) >= args.min_length:
                        break
            if len(seq) < args.min_length:
                continue
            name = record.oligo_id
            if name in seen_ids:
                suffix = 2
                while f"{name}.{suffix}" in seen_ids:
                    suffix += 1
                name = f"{name}.{suffix}"
            seen_ids.add(name)
            out.write(f">{name}\n{fasta_wrap(seq)}\n")
            written += 1
    print(f"wrote {written} reporter sequences to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
