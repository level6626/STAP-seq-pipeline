#!/usr/bin/env python3
"""Compare plasmid R2 barcode overlap between STAP-seq and TAPS-seq BAMs."""

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

STAP_SUFFIX_RE = re.compile(r"_([ACGTNacgtn]{25})$")
TAPS_R2_RE = re.compile(r"(?:^|\|)R2=([ACGTNacgtn]{17})(?:\||$)")
TAPS_UMI_RE = re.compile(r"(?:^|\|)UMI=([^|]+)(?:\||$)")


def is_usable(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or read.mapping_quality < min_mapq
    )


def parse_stap_r2(query_name: str) -> tuple[str | None, str | None]:
    match = STAP_SUFFIX_RE.search(query_name)
    if not match:
        return None, None
    combined = match.group(1).upper()
    return combined[-17:], combined


def parse_taps_r2(query_name: str) -> tuple[str | None, str | None]:
    match = TAPS_R2_RE.search(query_name)
    if not match:
        return None, None
    r2 = match.group(1).upper()
    umi_match = TAPS_UMI_RE.search(query_name)
    molecule = umi_match.group(1) if umi_match else query_name
    return r2, molecule


def count_barcodes(
    bam_path: Path,
    assay: str,
    min_mapq: int,
    max_records: int,
) -> tuple[Counter[str], dict[str, set[str]], Counter[str]]:
    records: Counter[str] = Counter()
    molecules: dict[str, set[str]] = {}
    stats: Counter[str] = Counter()

    parser = parse_stap_r2 if assay == "stap" else parse_taps_r2
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            stats["bam_records_seen"] += 1
            if max_records and stats["bam_records_seen"] > max_records:
                stats["stopped_at_max_records"] = max_records
                break
            if not is_usable(read, min_mapq):
                stats["skipped_unusable_alignment"] += 1
                continue
            r2, molecule = parser(read.query_name)
            if r2 is None or molecule is None:
                stats["skipped_unparsed_barcode"] += 1
                continue
            if len(r2) != 17:
                stats["skipped_bad_r2_length"] += 1
                continue
            records[r2] += 1
            molecules.setdefault(r2, set()).add(molecule)
            stats["usable_records"] += 1

    return records, molecules, stats


def write_counts(
    path: Path,
    assay: str,
    records: Counter[str],
    molecules: dict[str, set[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write(
            f"r2_barcode\tmeth_code\tmeth_label\texpected_methylation\trandom14\t"
            f"{assay}_records\t{assay}_molecules\n"
        )
        for r2 in sorted(records):
            meth_code = r2[:3]
            meth_label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
            out.write(
                f"{r2}\t{meth_code}\t{meth_label}\t{expected:g}\t{r2[3:]}\t"
                f"{records[r2]}\t{len(molecules.get(r2, set()))}\n"
            )


def write_overlap(
    path: Path,
    stap_records: Counter[str],
    stap_molecules: dict[str, set[str]],
    taps_records: Counter[str],
    taps_molecules: dict[str, set[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    all_barcodes = sorted(set(stap_records) | set(taps_records))
    with path.open("w") as out:
        out.write(
            "r2_barcode\tmeth_code\tmeth_label\texpected_methylation\trandom14\t"
            "in_stap\tin_taps\tstap_records\tstap_molecules\t"
            "taps_records\ttaps_molecules\n"
        )
        for r2 in all_barcodes:
            meth_code = r2[:3]
            meth_label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
            stap_molecule_count = len(stap_molecules.get(r2, set()))
            taps_molecule_count = len(taps_molecules.get(r2, set()))
            out.write(
                f"{r2}\t{meth_code}\t{meth_label}\t{expected:g}\t{r2[3:]}\t"
                f"{int(stap_molecule_count > 0)}\t{int(taps_molecule_count > 0)}\t"
                f"{stap_records[r2]}\t{stap_molecule_count}\t"
                f"{taps_records[r2]}\t{taps_molecule_count}\n"
            )


def write_summary(
    path: Path,
    stap_records: Counter[str],
    stap_molecules: dict[str, set[str]],
    taps_records: Counter[str],
    taps_molecules: dict[str, set[str]],
    stap_stats: Counter[str],
    taps_stats: Counter[str],
    min_stap_molecules: int,
    min_taps_molecules: int,
) -> None:
    all_barcodes = set(stap_records) | set(taps_records)
    meth_codes = sorted({r2[:3] for r2 in all_barcodes} | set(METHYLATION))

    def present_stap(r2: str) -> bool:
        return len(stap_molecules.get(r2, set())) >= min_stap_molecules

    def present_taps(r2: str) -> bool:
        return len(taps_molecules.get(r2, set())) >= min_taps_molecules

    with path.open("w") as out:
        out.write(
            "scope\tmeth_code\tmeth_label\texpected_methylation\t"
            "stap_barcodes\ttaps_barcodes\toverlap_barcodes\t"
            "fraction_stap_in_taps\tfraction_taps_in_stap\n"
        )
        for scope, code_filter in [("all", None), *[("meth_code", code) for code in meth_codes]]:
            if code_filter is None:
                barcodes = all_barcodes
                meth_code = "ALL"
                meth_label = "ALL"
                expected = "nan"
            else:
                barcodes = {r2 for r2 in all_barcodes if r2[:3] == code_filter}
                meth_code = code_filter
                meth_label, expected_value = METHYLATION.get(code_filter, ("UNKNOWN", float("nan")))
                expected = f"{expected_value:g}"
            stap_set = {r2 for r2 in barcodes if present_stap(r2)}
            taps_set = {r2 for r2 in barcodes if present_taps(r2)}
            overlap = stap_set & taps_set
            frac_stap = len(overlap) / len(stap_set) if stap_set else float("nan")
            frac_taps = len(overlap) / len(taps_set) if taps_set else float("nan")
            out.write(
                f"{scope}\t{meth_code}\t{meth_label}\t{expected}\t"
                f"{len(stap_set)}\t{len(taps_set)}\t{len(overlap)}\t"
                f"{frac_stap:.6g}\t{frac_taps:.6g}\n"
            )

        out.write("\nmetric\tassay\tcount\n")
        for assay, stats in [("stap", stap_stats), ("taps", taps_stats)]:
            for key, value in sorted(stats.items()):
                out.write(f"{key}\t{assay}\t{value}\n")
        out.write(f"min_stap_molecules\tparameter\t{min_stap_molecules}\n")
        out.write(f"min_taps_molecules\tparameter\t{min_taps_molecules}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stap-bam", required=True, type=Path)
    parser.add_argument("--taps-bam", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--max-stap-records", type=int, default=0)
    parser.add_argument("--max-taps-records", type=int, default=0)
    parser.add_argument("--min-stap-molecules", type=int, default=1)
    parser.add_argument("--min-taps-molecules", type=int, default=1)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    stap_records, stap_molecules, stap_stats = count_barcodes(
        args.stap_bam, "stap", args.min_mapq, args.max_stap_records
    )
    taps_records, taps_molecules, taps_stats = count_barcodes(
        args.taps_bam, "taps", args.min_mapq, args.max_taps_records
    )

    write_counts(args.outdir / "barcode_counts_stap.tsv", "stap", stap_records, stap_molecules)
    write_counts(args.outdir / "barcode_counts_taps.tsv", "taps", taps_records, taps_molecules)
    write_overlap(
        args.outdir / "barcode_overlap.tsv",
        stap_records,
        stap_molecules,
        taps_records,
        taps_molecules,
    )
    write_summary(
        args.outdir / "barcode_overlap_summary.tsv",
        stap_records,
        stap_molecules,
        taps_records,
        taps_molecules,
        stap_stats,
        taps_stats,
        args.min_stap_molecules,
        args.min_taps_molecules,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
