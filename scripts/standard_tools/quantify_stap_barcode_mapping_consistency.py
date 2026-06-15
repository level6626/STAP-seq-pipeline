#!/usr/bin/env python3
"""Quantify mapping consistency within STAP plasmid/methylation barcode families."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

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

STAP_SUFFIX_RE = re.compile(r"_([ACGTNacgtn]{25})$")


@dataclass
class BarcodeFamily:
    r2_barcode: str
    meth_code: str
    plasmid_barcode: str
    reads: int = 0
    umis: set[str] = field(default_factory=set)
    locus_counts: Counter[tuple[str, str, int]] = field(default_factory=Counter)


@dataclass(frozen=True)
class FamilyMetrics:
    key: tuple[str, str]
    r2_barcode: str
    meth_code: str
    plasmid_barcode: str
    meth_label: str
    expected_methylation: float
    n_reads: int
    n_unique_umis: int
    n_loci: int
    major_locus: tuple[str, str, int]
    major_locus_reads: int
    exact_consistency: float
    exact_discordant_reads: int
    window_locus: tuple[str, str, int]
    window_reads: int
    window_consistency: float
    window_discordant_reads: int


def parse_stap_barcode(query_name: str) -> tuple[str | None, str | None, str | None]:
    match = STAP_SUFFIX_RE.search(query_name)
    if not match:
        return None, None, None
    combined = match.group(1).upper()
    r1_umi = combined[:8]
    r2 = combined[8:]
    return r2, r1_umi, combined


def is_usable(read: pysam.AlignedSegment, min_mapq: int, include_duplicates: bool) -> bool:
    if (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.mapping_quality < min_mapq
    ):
        return False
    if read.is_duplicate and not include_duplicates:
        return False
    return True


def keep_read_number(read: pysam.AlignedSegment, selection: str) -> bool:
    if selection == "both":
        return True
    if selection == "read1":
        return read.is_read1
    if selection == "read2":
        return read.is_read2
    raise ValueError(f"Unknown read selection: {selection}")


def read_tss_1based(read: pysam.AlignedSegment) -> int | None:
    if read.reference_start < 0:
        return None
    if read.is_reverse:
        return read.reference_end
    return read.reference_start + 1


def locus_string(locus: tuple[str, str, int]) -> str:
    chrom, strand, tss = locus
    return f"{chrom}:{strand}:{tss}"


def dense_window(
    locus_counts: Counter[tuple[str, str, int]],
    window_bp: int,
) -> tuple[tuple[str, str, int], int]:
    if not locus_counts:
        raise ValueError("Cannot choose a window from an empty locus counter")
    if window_bp == 0:
        locus, count = max(locus_counts.items(), key=lambda item: (item[1], item[0]))
        return locus, count

    by_chrom_strand: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for (chrom, strand, tss), count in locus_counts.items():
        by_chrom_strand[(chrom, strand)].append((tss, count))

    best_locus = max(locus_counts.items(), key=lambda item: (item[1], item[0]))[0]
    best_count = 0
    for (chrom, strand), positions in by_chrom_strand.items():
        positions.sort()
        left = 0
        right = 0
        running = 0
        for center_pos, _center_count in positions:
            while right < len(positions) and positions[right][0] <= center_pos + window_bp:
                running += positions[right][1]
                right += 1
            while left < len(positions) and positions[left][0] < center_pos - window_bp:
                running -= positions[left][1]
                left += 1
            candidate = (chrom, strand, center_pos)
            if running > best_count or (running == best_count and candidate < best_locus):
                best_count = running
                best_locus = candidate
    return best_locus, best_count


def median(values: list[float]) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def format_float(value: float) -> str:
    return f"{value:.6g}"


def iter_family_metrics(
    families: dict[tuple[str, str], BarcodeFamily],
    window_bp: int,
) -> Iterable[FamilyMetrics]:
    for key in sorted(families):
        family = families[key]
        if not family.locus_counts:
            continue
        major_locus, major_count = max(
            family.locus_counts.items(),
            key=lambda item: (item[1], item[0]),
        )
        window_locus, window_count = dense_window(family.locus_counts, window_bp)
        meth_label, expected = METHYLATION.get(family.meth_code, ("UNKNOWN", float("nan")))
        yield FamilyMetrics(
            key=key,
            r2_barcode=family.r2_barcode,
            meth_code=family.meth_code,
            plasmid_barcode=family.plasmid_barcode,
            meth_label=meth_label,
            expected_methylation=expected,
            n_reads=family.reads,
            n_unique_umis=len(family.umis),
            n_loci=len(family.locus_counts),
            major_locus=major_locus,
            major_locus_reads=major_count,
            exact_consistency=major_count / family.reads,
            exact_discordant_reads=family.reads - major_count,
            window_locus=window_locus,
            window_reads=window_count,
            window_consistency=window_count / family.reads,
            window_discordant_reads=family.reads - window_count,
        )


def write_family_table(
    path: Path,
    metrics: list[FamilyMetrics],
    families: dict[tuple[str, str], BarcodeFamily],
    max_loci_report: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write(
            "meth_code\tmeth_label\texpected_methylation\tplasmid_barcode\tr2_barcode\t"
            "n_reads\tn_unique_r1_umis\tn_loci\tmajor_locus\tmajor_locus_reads\t"
            "exact_consistency\texact_discordant_reads\twindow_locus\twindow_reads\t"
            "window_consistency\twindow_discordant_reads\tall_loci_counts\n"
        )
        for row in metrics:
            family = families[row.key]
            loci = sorted(
                family.locus_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
            locus_items = [f"{locus_string(locus)}={count}" for locus, count in loci[:max_loci_report]]
            if len(loci) > max_loci_report:
                locus_items.append(f"...{len(loci) - max_loci_report}_more")
            out.write(
                f"{row.meth_code}\t{row.meth_label}\t{format_float(row.expected_methylation)}\t"
                f"{row.plasmid_barcode}\t{row.r2_barcode}\t{row.n_reads}\t"
                f"{row.n_unique_umis}\t{row.n_loci}\t{locus_string(row.major_locus)}\t"
                f"{row.major_locus_reads}\t{format_float(row.exact_consistency)}\t"
                f"{row.exact_discordant_reads}\t{locus_string(row.window_locus)}\t"
                f"{row.window_reads}\t{format_float(row.window_consistency)}\t"
                f"{row.window_discordant_reads}\t{';'.join(locus_items)}\n"
            )


def passes_summary_filter(row: FamilyMetrics, min_reads: int, min_unique_umis: int) -> bool:
    return row.n_reads >= min_reads and row.n_unique_umis >= min_unique_umis


def summarize_rows(
    rows: list[FamilyMetrics],
    consistency_threshold: float,
) -> dict[str, str]:
    families = len(rows)
    reads = sum(row.n_reads for row in rows)
    unique_umis = sum(row.n_unique_umis for row in rows)
    exact_discordant = sum(row.exact_discordant_reads for row in rows)
    window_discordant = sum(row.window_discordant_reads for row in rows)
    multilocus = sum(row.n_loci > 1 for row in rows)
    exact_pass = sum(row.exact_consistency >= consistency_threshold for row in rows)
    window_pass = sum(row.window_consistency >= consistency_threshold for row in rows)
    return {
        "families": str(families),
        "reads": str(reads),
        "unique_r1_umis": str(unique_umis),
        "multilocus_families": str(multilocus),
        "fraction_multilocus_families": format_float(multilocus / families) if families else "nan",
        "median_exact_consistency": format_float(median([row.exact_consistency for row in rows])),
        "families_exact_consistency_ge_threshold": str(exact_pass),
        "fraction_exact_consistency_ge_threshold": format_float(exact_pass / families) if families else "nan",
        "read_weighted_exact_discordance": format_float(exact_discordant / reads) if reads else "nan",
        "median_window_consistency": format_float(median([row.window_consistency for row in rows])),
        "families_window_consistency_ge_threshold": str(window_pass),
        "fraction_window_consistency_ge_threshold": format_float(window_pass / families) if families else "nan",
        "read_weighted_window_discordance": format_float(window_discordant / reads) if reads else "nan",
    }


def write_summary_table(
    path: Path,
    metrics: list[FamilyMetrics],
    stats: Counter[str],
    args: argparse.Namespace,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    filtered = [
        row
        for row in metrics
        if passes_summary_filter(row, args.min_family_reads, args.min_family_umis)
    ]
    meth_codes = sorted({row.meth_code for row in filtered} | set(METHYLATION))
    columns = [
        "scope",
        "meth_code",
        "meth_label",
        "expected_methylation",
        "families",
        "reads",
        "unique_r1_umis",
        "multilocus_families",
        "fraction_multilocus_families",
        "median_exact_consistency",
        "families_exact_consistency_ge_threshold",
        "fraction_exact_consistency_ge_threshold",
        "read_weighted_exact_discordance",
        "median_window_consistency",
        "families_window_consistency_ge_threshold",
        "fraction_window_consistency_ge_threshold",
        "read_weighted_window_discordance",
    ]

    with path.open("w") as out:
        out.write("\t".join(columns) + "\n")
        for meth_code in ["ALL", *meth_codes]:
            if meth_code == "ALL":
                rows = filtered
                meth_label = "ALL"
                expected = "nan"
                scope = "all"
            else:
                rows = [row for row in filtered if row.meth_code == meth_code]
                meth_label, expected_value = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
                expected = format_float(expected_value)
                scope = "meth_code"
            summary = summarize_rows(rows, args.consistency_threshold)
            out.write(
                "\t".join(
                    [
                        scope,
                        meth_code,
                        meth_label,
                        expected,
                        *[summary[column] for column in columns[4:]],
                    ]
                )
                + "\n"
            )

        out.write("\nmetric\tvalue\n")
        out.write(f"bam_path\t{args.bam}\n")
        out.write(f"read_selection\t{args.read_selection}\n")
        out.write(f"min_mapq\t{args.min_mapq}\n")
        out.write(f"include_duplicates\t{int(args.include_duplicates)}\n")
        out.write(f"require_known_meth_code\t{int(args.require_known_meth_code)}\n")
        out.write(f"window_bp\t{args.window_bp}\n")
        out.write(f"consistency_threshold\t{format_float(args.consistency_threshold)}\n")
        out.write(f"min_family_reads\t{args.min_family_reads}\n")
        out.write(f"min_family_umis\t{args.min_family_umis}\n")
        out.write(f"families_total_before_summary_filter\t{len(metrics)}\n")
        out.write(f"families_total_after_summary_filter\t{len(filtered)}\n")
        for key, value in sorted(stats.items()):
            out.write(f"{key}\t{value}\n")


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--out-family", required=True, type=Path)
    parser.add_argument("--out-summary", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument(
        "--read-selection",
        choices=["read1", "read2", "both"],
        default="read1",
        help="Which aligned read records to use when defining barcode-family loci.",
    )
    parser.add_argument(
        "--include-duplicates",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Include alignments flagged duplicate. The deduplicated BAM normally does not need this.",
    )
    parser.add_argument(
        "--require-known-meth-code",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Drop R2 barcodes whose first 3 bp are not one of the expected methylation codes.",
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=5,
        help="Window size for tolerant consistency. A value of 5 counts loci within +/-5 bp.",
    )
    parser.add_argument(
        "--consistency-threshold",
        type=float,
        default=0.9,
        help="Threshold used for pass/fail counts in the summary table.",
    )
    parser.add_argument(
        "--min-family-reads",
        type=int,
        default=2,
        help="Minimum family reads included in summary metrics. Per-family output keeps all families.",
    )
    parser.add_argument(
        "--min-family-umis",
        type=int,
        default=0,
        help="Minimum unique R1 UMIs included in summary metrics.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=0,
        help="Stop after this many BAM records. Use 0 for the full BAM.",
    )
    parser.add_argument("--max-loci-report", type=int, default=50)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    if args.window_bp < 0:
        raise ValueError("--window-bp must be >= 0")
    if not 0 <= args.consistency_threshold <= 1:
        raise ValueError("--consistency-threshold must be between 0 and 1")
    if args.min_family_reads < 1:
        raise ValueError("--min-family-reads must be >= 1")
    if args.min_family_umis < 0:
        raise ValueError("--min-family-umis must be >= 0")
    if args.max_loci_report < 1:
        raise ValueError("--max-loci-report must be >= 1")

    families: dict[tuple[str, str], BarcodeFamily] = {}
    stats: Counter[str] = Counter()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if args.max_records and stats["bam_records_seen"] >= args.max_records:
                stats["stopped_at_max_records"] = args.max_records
                break
            stats["bam_records_seen"] += 1
            if args.progress_every and stats["bam_records_seen"] % args.progress_every == 0:
                print(f"processed {stats['bam_records_seen']:,} BAM records", file=sys.stderr, flush=True)
            if not keep_read_number(read, args.read_selection):
                stats["skipped_read_selection"] += 1
                continue
            if not is_usable(read, args.min_mapq, args.include_duplicates):
                stats["skipped_unusable_alignment"] += 1
                continue
            r2_barcode, r1_umi, _combined = parse_stap_barcode(read.query_name)
            if r2_barcode is None or r1_umi is None:
                stats["skipped_unparsed_barcode"] += 1
                continue
            if len(r2_barcode) != 17:
                stats["skipped_bad_r2_length"] += 1
                continue
            meth_code = r2_barcode[:3]
            if meth_code not in METHYLATION:
                stats["unknown_meth_code_records"] += 1
                if args.require_known_meth_code:
                    stats["skipped_unknown_meth_code"] += 1
                    continue
            tss = read_tss_1based(read)
            if tss is None:
                stats["skipped_missing_tss"] += 1
                continue
            chrom = read.reference_name
            if chrom is None:
                stats["skipped_missing_reference_name"] += 1
                continue
            strand = "-" if read.is_reverse else "+"
            plasmid_barcode = r2_barcode[3:]
            key = (meth_code, plasmid_barcode)
            family = families.get(key)
            if family is None:
                family = BarcodeFamily(
                    r2_barcode=r2_barcode,
                    meth_code=meth_code,
                    plasmid_barcode=plasmid_barcode,
                )
                families[key] = family
            family.reads += 1
            family.umis.add(r1_umi)
            family.locus_counts[(chrom, strand, tss)] += 1
            stats["usable_records"] += 1

    metrics = list(iter_family_metrics(families, args.window_bp))
    write_family_table(args.out_family, metrics, families, args.max_loci_report)
    write_summary_table(args.out_summary, metrics, stats, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
