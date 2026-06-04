#!/usr/bin/env python3
"""Count TAPS C-to-T/G-to-A conversion at CpGs, grouped by R2 methylation code."""

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
UMI_RE = re.compile(r"(?:^|\|)UMI=([^|]+)")


def parse_meth(query_name: str) -> str:
    match = METH_RE.search(query_name)
    if not match:
        return "UNKNOWN"
    return match.group(1)


def parse_umi(query_name: str) -> str:
    match = UMI_RE.search(query_name)
    if not match:
        return query_name
    return match.group(1)


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
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--out-sites", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument(
        "--dedup-by-umi",
        action="store_true",
        help="Within each CpG/strand/methylation bin, count a UMI only once. Best for small references.",
    )
    args = parser.parse_args()

    site_counts: Counter[tuple[str, str, int, int, str, str]] = Counter()
    summary_counts: Counter[tuple[str, str]] = Counter()
    read_stats: Counter[str] = Counter()
    seen: set[tuple[str, str, int, str, str]] = set()

    fasta = pysam.FastaFile(str(args.reference_fasta))
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            read_stats["bam_records"] += 1
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
            meth_code = parse_meth(read.query_name)
            umi = parse_umi(read.query_name)

            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                ref_idx = rpos - start + offset
                if ref_idx < 0 or ref_idx >= len(ref_span):
                    continue
                ref_base = ref_span[ref_idx]
                read_base = query_seq[qpos].upper()

                if ref_base == "C" and ref_idx + 1 < len(ref_span) and ref_span[ref_idx + 1] == "G":
                    if read_base not in {"C", "T"}:
                        outcome = "other"
                    else:
                        outcome = "converted" if read_base == "T" else "unconverted"
                    cpg_c_pos = rpos + 1
                    cpg_g_pos = rpos + 2
                    observed_strand = "plus_c"
                elif ref_base == "G" and ref_idx - 1 >= 0 and ref_span[ref_idx - 1] == "C":
                    if read_base not in {"G", "A"}:
                        outcome = "other"
                    else:
                        outcome = "converted" if read_base == "A" else "unconverted"
                    cpg_c_pos = rpos
                    cpg_g_pos = rpos + 1
                    observed_strand = "minus_c"
                else:
                    continue

                if args.dedup_by_umi:
                    dedup_key = (meth_code, chrom, cpg_c_pos, observed_strand, umi)
                    if dedup_key in seen:
                        continue
                    seen.add(dedup_key)

                key = (meth_code, chrom, cpg_c_pos, cpg_g_pos, observed_strand, outcome)
                site_counts[key] += 1
                summary_counts[(meth_code, outcome)] += 1
                read_stats["cpg_observations"] += 1

    args.out_sites.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    with args.out_sites.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_conversion\tchrom\tcpg_c_pos_1based\t"
            "cpg_g_pos_1based\tobserved_strand\tconverted_count\tunconverted_count\t"
            "other_count\tcallable_count\tconversion_rate\n"
        )
        grouped: dict[tuple[str, str, int, int, str], Counter[str]] = {}
        for (meth_code, chrom, c_pos, g_pos, strand, outcome), count in site_counts.items():
            grouped.setdefault((meth_code, chrom, c_pos, g_pos, strand), Counter())[outcome] += count
        for meth_code, chrom, c_pos, g_pos, strand in sorted(grouped):
            counts = grouped[(meth_code, chrom, c_pos, g_pos, strand)]
            converted = counts["converted"]
            unconverted = counts["unconverted"]
            other = counts["other"]
            callable_count = converted + unconverted
            rate = converted / callable_count if callable_count else float("nan")
            label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
            out.write(
                f"{meth_code}\t{label}\t{expected:g}\t{chrom}\t{c_pos}\t{g_pos}\t{strand}\t"
                f"{converted}\t{unconverted}\t{other}\t{callable_count}\t{rate:.6g}\n"
            )

    with args.out_summary.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_conversion\tconverted_count\tunconverted_count\t"
            "other_count\tcallable_count\tconversion_rate\n"
        )
        for meth_code in sorted({key[0] for key in summary_counts} | set(METHYLATION)):
            converted = summary_counts[(meth_code, "converted")]
            unconverted = summary_counts[(meth_code, "unconverted")]
            other = summary_counts[(meth_code, "other")]
            callable_count = converted + unconverted
            rate = converted / callable_count if callable_count else float("nan")
            label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
            out.write(
                f"{meth_code}\t{label}\t{expected:g}\t{converted}\t{unconverted}\t"
                f"{other}\t{callable_count}\t{rate:.6g}\n"
            )
        out.write("\nmetric\tcount\n")
        for key, value in sorted(read_stats.items()):
            out.write(f"{key}\t{value}\n")
        if args.dedup_by_umi:
            out.write(f"dedup_keys_seen\t{len(seen)}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
