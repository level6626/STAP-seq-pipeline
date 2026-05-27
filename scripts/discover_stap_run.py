#!/usr/bin/env python3
"""Create a STAP-only sample manifest from an Illumina run folder."""

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path


FASTQ_RE = re.compile(r"^(STAP_.+)_R([123])_001\.fastq\.gz$")


def first_read_length(path):
    with gzip.open(path, "rt") as handle:
        handle.readline()
        seq = handle.readline().strip()
    return len(seq)


def classify_sample(sample, reads):
    if sample.startswith("STAP_TSS_") and reads == {"R1", "R2", "R3"}:
        if "500_oligos" in sample or "6_oligos" in sample:
            return "tss_oligo_control"
        return "tss_genomic"
    if sample.startswith("STAP_scTSSV2C_") and reads == {"R1", "R2"}:
        return "sc_tssv2c_unprocessed"
    return "unknown_stap_layout"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="../data/Run202", help="Run folder")
    parser.add_argument("--out", required=True, help="Output TSV")
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    samples = defaultdict(dict)
    for fastq in sorted(run_dir.glob("STAP_*.fastq.gz")):
        match = FASTQ_RE.match(fastq.name)
        if not match:
            continue
        sample, read_num = match.groups()
        samples[sample][f"R{read_num}"] = fastq

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        fields = [
            "sample",
            "assay_branch",
            "read_count",
            "R1",
            "R2",
            "R3",
            "R1_len",
            "R2_len",
            "R3_len",
        ]
        handle.write("\t".join(fields) + "\n")
        for sample in sorted(samples):
            reads = samples[sample]
            read_keys = set(reads)
            row = {
                "sample": sample,
                "assay_branch": classify_sample(sample, read_keys),
                "read_count": str(len(read_keys)),
                "R1": str(reads.get("R1", "")),
                "R2": str(reads.get("R2", "")),
                "R3": str(reads.get("R3", "")),
                "R1_len": str(first_read_length(reads["R1"])) if "R1" in reads else "",
                "R2_len": str(first_read_length(reads["R2"])) if "R2" in reads else "",
                "R3_len": str(first_read_length(reads["R3"])) if "R3" in reads else "",
            }
            handle.write("\t".join(row[field] for field in fields) + "\n")


if __name__ == "__main__":
    main()
