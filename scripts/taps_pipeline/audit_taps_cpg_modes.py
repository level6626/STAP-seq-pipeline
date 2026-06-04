#!/usr/bin/env python3
"""Audit TAPS CpG conversion calls with several strand/base counting modes."""

from __future__ import annotations

import argparse
import re
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
    mode: str,
    meth_code: str,
    ref_base: str,
    read_base: str,
) -> None:
    if ref_base == "C":
        if read_base == "T":
            outcome = "converted"
        elif read_base == "C":
            outcome = "unconverted"
        else:
            outcome = "other"
    elif ref_base == "G":
        if read_base == "A":
            outcome = "converted"
        elif read_base == "G":
            outcome = "unconverted"
        else:
            outcome = "other"
    else:
        return
    counts[(mode, meth_code, outcome)] += 1


def write_summary(
    path: Path,
    counts: Counter[tuple[str, str, str]],
    stats: Counter[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    modes = sorted({mode for mode, _code, _outcome in counts})
    codes = sorted({code for _mode, code, _outcome in counts} | set(METHYLATION))
    with path.open("w") as out:
        out.write(
            "mode\tmeth_code\tmeth_label\texpected_conversion\tconverted_count\t"
            "unconverted_count\tother_count\tcallable_count\tconversion_rate\n"
        )
        for mode in modes:
            for code in codes:
                converted = counts[(mode, code, "converted")]
                unconverted = counts[(mode, code, "unconverted")]
                other = counts[(mode, code, "other")]
                callable_count = converted + unconverted
                rate = converted / callable_count if callable_count else float("nan")
                label, expected = METHYLATION.get(code, ("UNKNOWN", float("nan")))
                out.write(
                    f"{mode}\t{code}\t{label}\t{expected:g}\t{converted}\t"
                    f"{unconverted}\t{other}\t{callable_count}\t{rate:.6g}\n"
                )
        out.write("\nmetric\tcount\n")
        for key, value in sorted(stats.items()):
            out.write(f"{key}\t{value}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--max-records", type=int, default=0)
    args = parser.parse_args()

    counts: Counter[tuple[str, str, str]] = Counter()
    stats: Counter[str] = Counter()

    fasta = pysam.FastaFile(str(args.reference_fasta))
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if args.max_records and stats["bam_records_seen"] >= args.max_records:
                stats["stopped_at_max_records"] = args.max_records
                break
            stats["bam_records_seen"] += 1
            if not is_usable(read, args.min_mapq):
                stats["skipped_unusable_alignment"] += 1
                continue

            chrom = bam.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end
            query_seq = read.query_sequence
            if chrom is None or start is None or end is None or end <= start or not query_seq:
                stats["skipped_missing_alignment_data"] += 1
                continue

            meth_code = parse_meth(read.query_name)
            ref_span = fasta.fetch(chrom, max(0, start - 1), end + 1).upper()
            offset = 1 if start > 0 else 0
            stats["usable_alignments"] += 1

            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                ref_idx = rpos - start + offset
                if ref_idx < 0 or ref_idx >= len(ref_span):
                    continue
                ref_base = ref_span[ref_idx]
                read_base = query_seq[qpos].upper()

                is_cpg_c = ref_base == "C" and ref_idx + 1 < len(ref_span) and ref_span[ref_idx + 1] == "G"
                is_cpg_g = ref_base == "G" and ref_idx - 1 >= 0 and ref_span[ref_idx - 1] == "C"
                if not (is_cpg_c or is_cpg_g):
                    continue

                stats["cpg_base_observations"] += 1

                # Current production mode: count both CpG bases when observed.
                if is_cpg_c or is_cpg_g:
                    add_call(counts, "both_cpg_bases", meth_code, ref_base, read_base)

                # Only count reference CpG C observations.
                if is_cpg_c:
                    add_call(counts, "ref_c_only_all_reads", meth_code, ref_base, read_base)

                # Only count reference CpG G observations.
                if is_cpg_g:
                    add_call(counts, "ref_g_only_all_reads", meth_code, ref_base, read_base)

                # Strand-aware candidate for single-end directional reads:
                # forward alignments count plus-strand C; reverse alignments count reference G.
                if (not read.is_reverse and is_cpg_c) or (read.is_reverse and is_cpg_g):
                    add_call(counts, "read_strand_expected", meth_code, ref_base, read_base)

                # Opposite of the above. If this shows high TTT conversion, read orientation
                # assumptions are reversed.
                if (not read.is_reverse and is_cpg_g) or (read.is_reverse and is_cpg_c):
                    add_call(counts, "read_strand_opposite", meth_code, ref_base, read_base)

    write_summary(args.out, counts, stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
