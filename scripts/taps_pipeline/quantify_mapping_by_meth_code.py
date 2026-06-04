#!/usr/bin/env python3
"""Quantify TAPS mapping status by R2 methylation code from a tagged BAM."""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

import pysam


METHYLATION = {
    "TTT": "100%",
    "AAA": "0%",
    "CAT": "60%",
    "AGT": "40%",
    "TGA": "20%",
    "TAG": "10%",
    "CTA": "1%",
    "ATG": "0.1%",
}

METH_RE = re.compile(r"(?:^|\|)METH_CODE=([^|]+)")


def parse_meth(query_name: str) -> str:
    match = METH_RE.search(query_name)
    return match.group(1) if match else "UNKNOWN"


def load_prepare_counts(path: Path | None) -> dict[str, int]:
    if path is None:
        return {}
    counts: dict[str, int] = {}
    with path.open() as handle:
        header = handle.readline()
        for line in handle:
            metric, count = line.rstrip("\n").split("\t")
            if metric.startswith("meth_code_"):
                code = metric.split("_")[2]
                counts[code] = int(count)
    return counts


def parse_star_final_log(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    metrics: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = line.rstrip("\n").split("|", 1)
            metrics[key.strip()] = value.strip()
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--prepare-stats",
        type=Path,
        help="Optional prepare_taps stats file for per-code input denominators.",
    )
    parser.add_argument(
        "--star-final-log",
        type=Path,
        help="Optional STAR Log.final.out for global multimapping metrics.",
    )
    parser.add_argument("--min-mapq", type=int, default=0)
    args = parser.parse_args()

    input_counts = load_prepare_counts(args.prepare_stats)
    star_metrics = parse_star_final_log(args.star_final_log)

    counts: Counter[tuple[str, str]] = Counter()
    nh_hist: Counter[int] = Counter()
    total_primary_records = Counter()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            code = parse_meth(read.query_name)
            total_primary_records[code] += 1
            if read.is_unmapped:
                counts[(code, "unmapped_in_bam")] += 1
                continue
            if read.is_qcfail or read.is_duplicate or read.mapping_quality < args.min_mapq:
                counts[(code, "mapped_filtered")] += 1
                continue
            nh = read.get_tag("NH") if read.has_tag("NH") else 0
            nh_hist[nh] += 1
            if nh == 1:
                counts[(code, "unique_in_bam")] += 1
            elif nh > 1:
                counts[(code, "multi_in_bam")] += 1
            else:
                counts[(code, "mapped_unknown_nh")] += 1

    all_codes = sorted(
        set(input_counts)
        | {code for code, _status in counts}
        | set(METHYLATION)
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        out.write(
            "meth_code\tmeth_label\tinput_tagged\tprimary_bam_records\t"
            "unique_in_bam\tmulti_in_bam\tmapped_unknown_nh\tmapped_filtered\t"
            "unmapped_in_bam\tbam_mapped_records\tmapped_fraction_vs_input\t"
            "unique_fraction_vs_input\tmulti_fraction_among_bam_mapped\n"
        )
        for code in all_codes:
            input_count = input_counts.get(code, 0)
            unique = counts[(code, "unique_in_bam")]
            multi = counts[(code, "multi_in_bam")]
            unknown = counts[(code, "mapped_unknown_nh")]
            filtered = counts[(code, "mapped_filtered")]
            unmapped = counts[(code, "unmapped_in_bam")]
            bam_mapped = unique + multi + unknown + filtered
            mapped_fraction = bam_mapped / input_count if input_count else float("nan")
            unique_fraction = unique / input_count if input_count else float("nan")
            multi_fraction = multi / (unique + multi + unknown) if (unique + multi + unknown) else float("nan")
            out.write(
                f"{code}\t{METHYLATION.get(code, 'UNKNOWN')}\t{input_count}\t"
                f"{total_primary_records[code]}\t{unique}\t{multi}\t{unknown}\t"
                f"{filtered}\t{unmapped}\t{bam_mapped}\t{mapped_fraction:.6g}\t"
                f"{unique_fraction:.6g}\t{multi_fraction:.6g}\n"
            )

        out.write("\nmetric\tvalue\n")
        out.write(f"bam_path\t{args.bam}\n")
        out.write(f"min_mapq\t{args.min_mapq}\n")
        for nh, count in sorted(nh_hist.items()):
            out.write(f"NH_{nh}_primary_mapped_records\t{count}\n")
        if nh_hist and sum(count for nh, count in nh_hist.items() if nh > 1) == 0:
            out.write(
                "warning\tNo NH>1 primary alignments were observed in this BAM. "
                "If STAR was run with --outFilterMultimapNmax 1, per-code "
                "multimapper rates cannot be recovered from this BAM; rerun STAR "
                "with multimappers retained.\n"
            )
        if star_metrics:
            out.write("\nstar_final_metric\tvalue\n")
            for key in (
                "Number of input reads",
                "Uniquely mapped reads number",
                "Uniquely mapped reads %",
                "Number of reads mapped to multiple loci",
                "% of reads mapped to multiple loci",
                "Number of reads mapped to too many loci",
                "% of reads mapped to too many loci",
            ):
                if key in star_metrics:
                    out.write(f"{key}\t{star_metrics[key]}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
