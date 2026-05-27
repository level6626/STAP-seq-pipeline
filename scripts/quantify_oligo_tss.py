#!/usr/bin/env python3
"""Dependency-light STAP TSS oligo quantification.

This is for the Run202 oligo-control triplets. It assigns each read pair to a
design sequence by exact k-mer matching against R1/R3, parses R1/R2 UMIs and
methylation barcodes, and reports raw and UMI-collapsed counts per design
position and methylation code.
"""

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
METHYLATION_CODES = {
    "TTT": "100",
    "AAA": "0",
    "CAT": "60",
    "AGT": "40",
    "TGA": "20",
    "TAG": "10",
    "CTA": "1",
    "ATG": "0.1",
}


def revcomp(seq):
    return seq.translate(DNA_COMP)[::-1].upper()


def read_fasta(path):
    records = {}
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records[name] = "".join(seq).upper()
    return records


def build_kmer_index(records, k):
    index = {}
    ambiguous = set()
    for ref_id, seq in records.items():
        for pos in range(0, len(seq) - k + 1):
            kmer = seq[pos : pos + k]
            if "N" in kmer:
                continue
            hits = [(kmer, "+", pos), (revcomp(kmer), "-", pos)]
            for key, strand, ref_pos in hits:
                value = (ref_id, strand, ref_pos)
                if key in index and index[key] != value:
                    ambiguous.add(key)
                else:
                    index[key] = value
    for key in ambiguous:
        index.pop(key, None)
    return index


def fastq_iter(path):
    with gzip.open(path, "rt") as handle:
        while True:
            name = handle.readline().rstrip()
            if not name:
                return
            seq = handle.readline().rstrip().upper()
            handle.readline()
            qual = handle.readline().rstrip()
            yield name.split()[0][1:], seq, qual


def best_hit(seq, index, k):
    for read_pos in range(0, len(seq) - k + 1):
        hit = index.get(seq[read_pos : read_pos + k])
        if hit:
            ref_id, strand, ref_pos = hit
            if strand == "+":
                tss_pos = ref_pos - read_pos + 1
            else:
                tss_pos = ref_pos + k + read_pos
            return ref_id, strand, tss_pos, read_pos
    return None


def parse_r2(seq, orientation):
    oriented = revcomp(seq) if orientation == "reverse-complement" else seq
    code = oriented[:3]
    return {
        "r2_oriented": oriented,
        "methylation_code": code,
        "methylation_percent": METHYLATION_CODES.get(code, "unknown"),
        "r2_random_barcode": oriented[3:],
    }


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample prefix without _R1/_R2/_R3")
    parser.add_argument("--run-dir", default="../data/Run202")
    parser.add_argument("--reference", required=True, help="Design FASTA")
    parser.add_argument("--outdir", default="results/oligo_tss")
    parser.add_argument("--kmer", type=int, default=25)
    parser.add_argument("--r1-umi-length", type=int, default=8)
    parser.add_argument(
        "--r2-orientation",
        choices=["forward", "reverse-complement"],
        default="forward",
        help="Orientation used to parse the R2 methylation code and random barcode",
    )
    parser.add_argument("--max-reads", type=int, default=0, help="Debug limit; 0 means all reads")
    args = parser.parse_args()

    if not args.sample.startswith("STAP_"):
        raise SystemExit("Refusing to process non-STAP sample")

    run_dir = Path(args.run_dir)
    r1 = run_dir / f"{args.sample}_R1_001.fastq.gz"
    r2 = run_dir / f"{args.sample}_R2_001.fastq.gz"
    r3 = run_dir / f"{args.sample}_R3_001.fastq.gz"
    for path in (r1, r2, r3):
        if not path.exists():
            raise SystemExit(f"Missing input FASTQ: {path}")

    records = read_fasta(args.reference)
    index = build_kmer_index(records, args.kmer)
    outdir = Path(args.outdir) / args.sample
    outdir.mkdir(parents=True, exist_ok=True)

    raw_counts = Counter()
    umi_counts = defaultdict(set)
    summary = Counter()

    assignment_path = outdir / "assignments.tsv.gz"
    with gzip.open(assignment_path, "wt", newline="") as out_handle:
        fields = [
            "read_id",
            "oligo_id",
            "position",
            "strand",
            "support_read",
            "r1_umi",
            "r2_raw",
            "r2_oriented",
            "methylation_code",
            "methylation_percent",
            "r2_random_barcode",
            "dedup_key",
        ]
        writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for i, ((id1, seq1, _), (id2, seq2, _), (id3, seq3, _)) in enumerate(
            zip(fastq_iter(r1), fastq_iter(r2), fastq_iter(r3)), start=1
        ):
            if args.max_reads and i > args.max_reads:
                break
            summary["total_reads"] += 1
            if id1 != id2 or id1 != id3:
                summary["read_id_mismatch"] += 1
                continue
            hit = best_hit(seq1, index, args.kmer)
            support = "R1"
            if hit is None:
                hit = best_hit(seq3, index, args.kmer)
                support = "R3"
            if hit is None:
                summary["unassigned"] += 1
                continue
            oligo_id, strand, pos, _ = hit
            r1_umi = seq1[: args.r1_umi_length]
            r2_info = parse_r2(seq2, args.r2_orientation)
            methylation_code = r2_info["methylation_code"]
            if r2_info["methylation_percent"] == "unknown":
                summary["unknown_methylation_code"] += 1
            key = (oligo_id, pos, strand, methylation_code)
            dedup_key = f"{r1_umi}:{r2_info['r2_oriented']}"
            raw_counts[key] += 1
            umi_counts[key].add(dedup_key)
            summary["assigned"] += 1
            writer.writerow(
                {
                    "read_id": id1,
                    "oligo_id": oligo_id,
                    "position": pos,
                    "strand": strand,
                    "support_read": support,
                    "r1_umi": r1_umi,
                    "r2_raw": seq2,
                    "r2_oriented": r2_info["r2_oriented"],
                    "methylation_code": methylation_code,
                    "methylation_percent": r2_info["methylation_percent"],
                    "r2_random_barcode": r2_info["r2_random_barcode"],
                    "dedup_key": dedup_key,
                }
            )

    count_rows = []
    for key in sorted(raw_counts):
        oligo_id, pos, strand, methylation_code = key
        count_rows.append(
            {
                "sample": args.sample,
                "oligo_id": oligo_id,
                "position": pos,
                "strand": strand,
                "methylation_code": methylation_code,
                "methylation_percent": METHYLATION_CODES.get(methylation_code, "unknown"),
                "raw_count": raw_counts[key],
                "dedup_count": len(umi_counts[key]),
            }
        )
    write_tsv(
        outdir / "tss_counts.tsv",
        count_rows,
        [
            "sample",
            "oligo_id",
            "position",
            "strand",
            "methylation_code",
            "methylation_percent",
            "raw_count",
            "dedup_count",
        ],
    )

    summary_rows = [{"metric": key, "value": value} for key, value in sorted(summary.items())]
    summary_rows.extend(
        [
            {"metric": "reference_records", "value": len(records)},
            {"metric": "non_ambiguous_kmers", "value": len(index)},
            {"metric": "assignment_file", "value": assignment_path},
        ]
    )
    write_tsv(outdir / "summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
