#!/usr/bin/env python3
"""Count TAPS-like conversion at non-CpG cytosines, grouped by R2 methylation code."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pysam


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

METH_RE = re.compile(r"(?:^|\|)METH_CODE=([^|]+)")
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
CONTEXTS = ("CpA", "CpC", "CpT")
STRANDS = ("plus_c", "minus_c")
OUTCOMES = ("converted", "unconverted", "other")


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


def classify_plus_c_context(ref_span: str, ref_idx: int) -> str | None:
    """Return CpN context for a plus-strand reference C."""
    if ref_idx + 1 >= len(ref_span):
        return None
    next_base = ref_span[ref_idx + 1]
    if next_base not in COMPLEMENT:
        return None
    return f"Cp{next_base}"


def classify_minus_c_context(ref_span: str, ref_idx: int) -> str | None:
    """Return CpN context for a minus-strand C represented as a plus-strand G."""
    if ref_idx - 1 < 0:
        return None
    prev_plus_base = ref_span[ref_idx - 1]
    if prev_plus_base not in COMPLEMENT:
        return None
    minus_next_base = COMPLEMENT[prev_plus_base]
    return f"Cp{minus_next_base}"


def add_count(
    counts: Counter[tuple[str, str, str, str]],
    meth_code: str,
    context: str,
    observed_strand: str,
    outcome: str,
) -> None:
    counts[(meth_code, context, observed_strand, outcome)] += 1
    counts[(meth_code, context, "both", outcome)] += 1
    counts[(meth_code, "CpH", observed_strand, outcome)] += 1
    counts[(meth_code, "CpH", "both", outcome)] += 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a tagged TAPS BAM and count background TAPS-like conversion "
            "at non-CpG cytosines. Plus-strand non-CpG C>T and opposite-strand "
            "non-CpG G>A are counted as converted."
        )
    )
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument(
        "--include-cpg",
        action="store_true",
        help="Also report CpG conversion rows as a sanity-check baseline.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=0,
        help="Stop after this many BAM records. Use 0 for all records.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress to stderr every N BAM records. Use 0 to disable.",
    )
    args = parser.parse_args()

    counts: Counter[tuple[str, str, str, str]] = Counter()
    read_stats: Counter[str] = Counter()

    fasta = pysam.FastaFile(str(args.reference_fasta))
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            read_stats["bam_records"] += 1
            if args.max_records and read_stats["bam_records"] > args.max_records:
                read_stats["stopped_at_max_records"] = args.max_records
                break
            if args.progress_every and read_stats["bam_records"] % args.progress_every == 0:
                print(
                    f"processed {read_stats['bam_records']:,} BAM records",
                    file=sys.stderr,
                    flush=True,
                )
            if not is_usable(read, args.min_mapq):
                read_stats["skipped_records"] += 1
                continue

            chrom = bam.get_reference_name(read.reference_id)
            if chrom is None:
                read_stats["skipped_no_reference_name"] += 1
                continue
            start = read.reference_start
            end = read.reference_end
            if start is None or end is None or end <= start:
                read_stats["skipped_no_reference_span"] += 1
                continue

            ref_span = fasta.fetch(chrom, max(0, start - 1), end + 1).upper()
            offset = 1 if start > 0 else 0
            query_seq = read.query_sequence
            if not query_seq:
                read_stats["skipped_no_query_sequence"] += 1
                continue

            read_stats["usable_records"] += 1
            meth_code = parse_meth(read.query_name)

            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                ref_idx = rpos - start + offset
                if ref_idx < 0 or ref_idx >= len(ref_span):
                    continue

                ref_base = ref_span[ref_idx]
                read_base = query_seq[qpos].upper()
                if ref_base == "C":
                    context = classify_plus_c_context(ref_span, ref_idx)
                    if context is None:
                        read_stats["skipped_unknown_context"] += 1
                        continue
                    if context == "CpG" and not args.include_cpg:
                        continue
                    if context != "CpG":
                        read_stats["non_cpg_cytosine_observations"] += 1
                    outcome = "converted" if read_base == "T" else "unconverted" if read_base == "C" else "other"
                    if context == "CpG":
                        counts[(meth_code, context, "plus_c", outcome)] += 1
                        counts[(meth_code, context, "both", outcome)] += 1
                    else:
                        add_count(counts, meth_code, context, "plus_c", outcome)
                elif ref_base == "G":
                    context = classify_minus_c_context(ref_span, ref_idx)
                    if context is None:
                        read_stats["skipped_unknown_context"] += 1
                        continue
                    if context == "CpG" and not args.include_cpg:
                        continue
                    if context != "CpG":
                        read_stats["non_cpg_cytosine_observations"] += 1
                    outcome = "converted" if read_base == "A" else "unconverted" if read_base == "G" else "other"
                    if context == "CpG":
                        counts[(meth_code, context, "minus_c", outcome)] += 1
                        counts[(meth_code, context, "both", outcome)] += 1
                    else:
                        add_count(counts, meth_code, context, "minus_c", outcome)

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    codes = sorted({key[0] for key in counts} | set(METHYLATION))
    contexts = ["CpH", *CONTEXTS]
    if args.include_cpg:
        contexts.append("CpG")

    with args.out_summary.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_conversion\tcontext\tobserved_strand\t"
            "converted_count\tunconverted_count\tother_count\tcallable_count\t"
            "conversion_rate\n"
        )
        for meth_code in codes:
            label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
            for context in contexts:
                for observed_strand in (*STRANDS, "both"):
                    converted = counts[(meth_code, context, observed_strand, "converted")]
                    unconverted = counts[(meth_code, context, observed_strand, "unconverted")]
                    other = counts[(meth_code, context, observed_strand, "other")]
                    if converted == 0 and unconverted == 0 and other == 0:
                        continue
                    callable_count = converted + unconverted
                    rate = converted / callable_count if callable_count else float("nan")
                    out.write(
                        f"{meth_code}\t{label}\t{expected:g}\t{context}\t{observed_strand}\t"
                        f"{converted}\t{unconverted}\t{other}\t{callable_count}\t{rate:.6g}\n"
                    )

        out.write("\nmetric\tcount\n")
        out.write(f"bam_path\t{args.bam}\n")
        out.write(f"reference_fasta\t{args.reference_fasta}\n")
        out.write(f"min_mapq\t{args.min_mapq}\n")
        out.write(f"include_cpg\t{int(args.include_cpg)}\n")
        for key, value in sorted(read_stats.items()):
            out.write(f"{key}\t{value}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
