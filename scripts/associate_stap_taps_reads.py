#!/usr/bin/env python3
"""Associate aligned STAP and TAPS reads with matching R2 barcodes."""

from __future__ import annotations

import argparse
import gzip
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO

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


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")  # type: ignore[return-value]
    return path.open("w")


def is_usable(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or read.mapping_quality < min_mapq
    )


def keep_read_number(read: pysam.AlignedSegment, selection: str) -> bool:
    if selection == "both":
        return True
    if selection == "read1":
        return read.is_read1
    if selection == "read2":
        return read.is_read2
    raise ValueError(f"Unknown read selection: {selection}")


def read_number(read: pysam.AlignedSegment) -> str:
    if read.is_read1:
        return "read1"
    if read.is_read2:
        return "read2"
    return "unpaired"


def parse_stap_barcode(query_name: str) -> tuple[str | None, str | None, str | None]:
    match = STAP_SUFFIX_RE.search(query_name)
    if not match:
        return None, None, None
    combined = match.group(1).upper()
    return combined[8:], combined[:8], combined


def parse_taps_barcode(query_name: str) -> tuple[str | None, str | None, str | None]:
    match = TAPS_R2_RE.search(query_name)
    if not match:
        return None, None, None
    r2 = match.group(1).upper()
    umi_match = TAPS_UMI_RE.search(query_name)
    molecule_umi = umi_match.group(1).upper() if umi_match else ""
    r1_umi = molecule_umi[: -len(r2)] if molecule_umi.endswith(r2) else ""
    return r2, r1_umi, molecule_umi


def alignment_row(
    read: pysam.AlignedSegment,
    r2: str,
    r1_umi: str,
    molecule_umi: str,
) -> tuple[str, str, str, str, str, str, str, int, int, str, int, str, str]:
    return (
        r2,
        r2[:3],
        r2[3:],
        read.query_name,
        r1_umi,
        molecule_umi,
        read.reference_name,
        read.reference_start + 1,
        read.reference_end or 0,
        "-" if read.is_reverse else "+",
        read.mapping_quality,
        read.cigarstring or "",
        read_number(read),
    )


def configure_db(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute(
        """
        CREATE TABLE stap (
            r2 TEXT NOT NULL,
            meth_code TEXT NOT NULL,
            plasmid_barcode TEXT NOT NULL,
            read_id TEXT NOT NULL,
            r1_umi TEXT NOT NULL,
            molecule_umi TEXT NOT NULL,
            chrom TEXT NOT NULL,
            start_1based INTEGER NOT NULL,
            end_1based INTEGER NOT NULL,
            strand TEXT NOT NULL,
            mapq INTEGER NOT NULL,
            cigar TEXT NOT NULL,
            read_number TEXT NOT NULL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE taps (
            r2 TEXT NOT NULL,
            meth_code TEXT NOT NULL,
            plasmid_barcode TEXT NOT NULL,
            read_id TEXT NOT NULL,
            r1_umi TEXT NOT NULL,
            molecule_umi TEXT NOT NULL,
            chrom TEXT NOT NULL,
            start_1based INTEGER NOT NULL,
            end_1based INTEGER NOT NULL,
            strand TEXT NOT NULL,
            mapq INTEGER NOT NULL,
            cigar TEXT NOT NULL,
            read_number TEXT NOT NULL
        )
        """
    )


def insert_rows(
    conn: sqlite3.Connection,
    table: str,
    rows: list[tuple[str, str, str, str, str, str, str, int, int, str, int, str, str]],
) -> None:
    if not rows:
        return
    conn.executemany(
        f"""
        INSERT INTO {table} (
            r2, meth_code, plasmid_barcode, read_id, r1_umi, molecule_umi,
            chrom, start_1based, end_1based, strand, mapq, cigar, read_number
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )


def load_bam(
    conn: sqlite3.Connection,
    table: str,
    bam_path: Path,
    assay: str,
    min_mapq: int,
    read_selection: str,
    require_known_meth_code: bool,
    max_records: int,
    batch_size: int,
) -> Counter[str]:
    parser = parse_stap_barcode if assay == "stap" else parse_taps_barcode
    stats: Counter[str] = Counter()
    batch: list[tuple[str, str, str, str, str, str, str, int, int, str, int, str, str]] = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if max_records and stats["bam_records_seen"] >= max_records:
                stats["stopped_at_max_records"] = max_records
                break
            stats["bam_records_seen"] += 1
            if not is_usable(read, min_mapq):
                stats["skipped_unusable_alignment"] += 1
                continue
            if not keep_read_number(read, read_selection):
                stats["skipped_read_selection"] += 1
                continue
            r2, r1_umi, molecule_umi = parser(read.query_name)
            if r2 is None or r1_umi is None or molecule_umi is None:
                stats["skipped_unparsed_barcode"] += 1
                continue
            if len(r2) != 17:
                stats["skipped_bad_r2_length"] += 1
                continue
            meth_code = r2[:3]
            if meth_code not in METHYLATION:
                stats["skipped_unknown_meth_code"] += 1
                if require_known_meth_code:
                    continue

            batch.append(alignment_row(read, r2, r1_umi, molecule_umi))
            stats["loaded_records"] += 1
            if len(batch) >= batch_size:
                insert_rows(conn, table, batch)
                batch.clear()

    insert_rows(conn, table, batch)
    conn.commit()
    return stats


def create_indexes(conn: sqlite3.Connection) -> None:
    conn.execute("CREATE INDEX idx_stap_r2 ON stap(r2)")
    conn.execute("CREATE INDEX idx_taps_r2 ON taps(r2)")
    conn.commit()


def write_associations(
    conn: sqlite3.Connection,
    path: Path,
    max_associations: int,
    max_pairs_per_barcode: int,
) -> Counter[str]:
    stats: Counter[str] = Counter()
    per_barcode: Counter[str] = Counter()
    query = """
        SELECT
            s.r2, s.meth_code, s.plasmid_barcode,
            s.read_id, s.r1_umi, s.molecule_umi, s.chrom, s.start_1based,
            s.end_1based, s.strand, s.mapq, s.cigar, s.read_number,
            t.read_id, t.r1_umi, t.molecule_umi, t.chrom, t.start_1based,
            t.end_1based, t.strand, t.mapq, t.cigar, t.read_number
        FROM stap AS s
        JOIN taps AS t ON s.r2 = t.r2
    """
    header = (
        "r2_barcode\tmeth_code\tplasmid_barcode\t"
        "stap_read_id\tstap_r1_umi\tstap_molecule_umi\tstap_chrom\t"
        "stap_start_1based\tstap_end_1based\tstap_strand\tstap_mapq\t"
        "stap_cigar\tstap_read_number\t"
        "taps_read_id\ttaps_r1_umi\ttaps_molecule_umi\ttaps_chrom\t"
        "taps_start_1based\ttaps_end_1based\ttaps_strand\ttaps_mapq\t"
        "taps_cigar\ttaps_read_number\n"
    )

    with open_output(path) as out:
        out.write(header)
        for row in conn.execute(query):
            r2 = row[0]
            stats["association_rows_considered"] += 1
            if max_pairs_per_barcode and per_barcode[r2] >= max_pairs_per_barcode:
                stats["skipped_max_pairs_per_barcode"] += 1
                continue
            out.write("\t".join(map(str, row)) + "\n")
            per_barcode[r2] += 1
            stats["association_rows_written"] += 1
            if max_associations and stats["association_rows_written"] >= max_associations:
                stats["stopped_at_max_associations"] = max_associations
                break

    stats["associated_barcodes_written"] = len(per_barcode)
    return stats


def write_barcode_summary(conn: sqlite3.Connection, path: Path) -> None:
    query = """
        WITH
        stap_counts AS (
            SELECT r2, meth_code, plasmid_barcode, COUNT(*) AS stap_records
            FROM stap
            GROUP BY r2
        ),
        taps_counts AS (
            SELECT r2, meth_code, plasmid_barcode, COUNT(*) AS taps_records
            FROM taps
            GROUP BY r2
        ),
        keys AS (
            SELECT r2 FROM stap_counts
            UNION
            SELECT r2 FROM taps_counts
        )
        SELECT
            keys.r2,
            COALESCE(stap_counts.meth_code, taps_counts.meth_code) AS meth_code,
            COALESCE(stap_counts.plasmid_barcode, taps_counts.plasmid_barcode) AS plasmid_barcode,
            COALESCE(stap_counts.stap_records, 0) AS stap_records,
            COALESCE(taps_counts.taps_records, 0) AS taps_records,
            COALESCE(stap_counts.stap_records, 0) * COALESCE(taps_counts.taps_records, 0)
                AS association_rows
        FROM keys
        LEFT JOIN stap_counts ON keys.r2 = stap_counts.r2
        LEFT JOIN taps_counts ON keys.r2 = taps_counts.r2
        ORDER BY keys.r2
    """
    with open_output(path) as out:
        out.write(
            "r2_barcode\tmeth_code\tmeth_label\texpected_methylation\t"
            "plasmid_barcode\tin_stap\tin_taps\tstap_records\t"
            "taps_records\tassociation_rows\n"
        )
        for r2, meth_code, plasmid_barcode, stap_records, taps_records, associations in conn.execute(
            query
        ):
            meth_label, expected = METHYLATION.get(meth_code, ("UNKNOWN", float("nan")))
            out.write(
                f"{r2}\t{meth_code}\t{meth_label}\t{expected:g}\t"
                f"{plasmid_barcode}\t{int(stap_records > 0)}\t{int(taps_records > 0)}\t"
                f"{stap_records}\t{taps_records}\t{associations}\n"
            )


def scalar(conn: sqlite3.Connection, query: str) -> int:
    value = conn.execute(query).fetchone()[0]
    return int(value or 0)


def write_run_summary(
    conn: sqlite3.Connection,
    path: Path,
    stap_stats: Counter[str],
    taps_stats: Counter[str],
    association_stats: Counter[str],
    args: argparse.Namespace,
) -> None:
    metrics = Counter()
    metrics["stap_loaded_records"] = scalar(conn, "SELECT COUNT(*) FROM stap")
    metrics["taps_loaded_records"] = scalar(conn, "SELECT COUNT(*) FROM taps")
    metrics["stap_unique_barcodes"] = scalar(conn, "SELECT COUNT(DISTINCT r2) FROM stap")
    metrics["taps_unique_barcodes"] = scalar(conn, "SELECT COUNT(DISTINCT r2) FROM taps")
    metrics["overlap_barcodes"] = scalar(
        conn,
        """
        SELECT COUNT(*) FROM (
            SELECT DISTINCT stap.r2
            FROM stap
            INNER JOIN taps ON stap.r2 = taps.r2
        )
        """,
    )
    metrics["association_rows_possible"] = scalar(
        conn,
        """
        WITH
        s AS (SELECT r2, COUNT(*) n FROM stap GROUP BY r2),
        t AS (SELECT r2, COUNT(*) n FROM taps GROUP BY r2)
        SELECT SUM(s.n * t.n)
        FROM s INNER JOIN t ON s.r2 = t.r2
        """,
    )

    with open_output(path) as out:
        out.write("metric\tscope\tvalue\n")
        for key, value in sorted(metrics.items()):
            out.write(f"{key}\tderived\t{value}\n")
        for key, value in sorted(stap_stats.items()):
            out.write(f"{key}\tstap\t{value}\n")
        for key, value in sorted(taps_stats.items()):
            out.write(f"{key}\ttaps\t{value}\n")
        for key, value in sorted(association_stats.items()):
            out.write(f"{key}\tassociations\t{value}\n")

        parameter_items = {
            "stap_bam": args.stap_bam,
            "taps_bam": args.taps_bam,
            "min_mapq": args.min_mapq,
            "read_selection": args.read_selection,
            "require_known_meth_code": int(not args.include_unknown_meth_code),
            "max_stap_records": args.max_stap_records,
            "max_taps_records": args.max_taps_records,
            "max_associations": args.max_associations,
            "max_pairs_per_barcode": args.max_pairs_per_barcode,
        }
        for key, value in parameter_items.items():
            out.write(f"{key}\tparameter\t{value}\n")


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be >= 0")
    return parsed


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Associate aligned STAP RNA reads with aligned TAPS DNA reads by the "
            "shared 17-bp R2 methylation/plasmid barcode."
        )
    )
    parser.add_argument("--stap-bam", required=True, type=Path)
    parser.add_argument("--taps-bam", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument(
        "--read-selection",
        choices=["read1", "read2", "both"],
        default="read1",
        help="Which aligned mates to associate. Default: read1.",
    )
    parser.add_argument(
        "--include-unknown-meth-code",
        action="store_true",
        help="Keep reads whose R2 first 3 bp are not one of the known methylation codes.",
    )
    parser.add_argument("--max-stap-records", type=positive_int, default=0)
    parser.add_argument("--max-taps-records", type=positive_int, default=0)
    parser.add_argument(
        "--max-associations",
        type=positive_int,
        default=0,
        help="Stop writing association rows after this many rows. 0 means no limit.",
    )
    parser.add_argument(
        "--max-pairs-per-barcode",
        type=positive_int,
        default=0,
        help="Maximum association rows to write per barcode. 0 means no limit.",
    )
    parser.add_argument("--batch-size", type=positive_int, default=10000)
    parser.add_argument(
        "--keep-sqlite",
        action="store_true",
        help="Keep the intermediate SQLite database in the output directory.",
    )
    parser.add_argument(
        "--sqlite-db",
        type=Path,
        default=None,
        help="Optional path for the intermediate SQLite database.",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    db_path = args.sqlite_db or args.outdir / "associate_stap_taps_reads.sqlite"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()

    conn = sqlite3.connect(db_path)
    try:
        configure_db(conn)
        require_known = not args.include_unknown_meth_code
        stap_stats = load_bam(
            conn,
            "stap",
            args.stap_bam,
            "stap",
            args.min_mapq,
            args.read_selection,
            require_known,
            args.max_stap_records,
            args.batch_size,
        )
        taps_stats = load_bam(
            conn,
            "taps",
            args.taps_bam,
            "taps",
            args.min_mapq,
            args.read_selection,
            require_known,
            args.max_taps_records,
            args.batch_size,
        )
        create_indexes(conn)

        association_stats = write_associations(
            conn,
            args.outdir / "stap_taps_read_associations.tsv.gz",
            args.max_associations,
            args.max_pairs_per_barcode,
        )
        write_barcode_summary(conn, args.outdir / "barcode_association_summary.tsv")
        write_run_summary(
            conn,
            args.outdir / "associate_stap_taps_reads.summary.tsv",
            stap_stats,
            taps_stats,
            association_stats,
            args,
        )
    finally:
        conn.close()

    if not args.keep_sqlite and db_path.exists():
        db_path.unlink()

    print(f"Wrote outputs to {args.outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
