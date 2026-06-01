#!/usr/bin/env python3
"""Demultiplex STAP triplet FASTQs by R2 methylation and R3 oligo barcode."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


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

DNA_ALPHABET = "ACGT"
XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
XLSX_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKG_REL_NS = "{http://schemas.openxmlformats.org/package/2006/relationships}"


@dataclass(frozen=True)
class OligoRecord:
    oligo_id: str
    sheet: str
    row_number: int
    barcode: str
    metadata: dict[str, str]


@dataclass(frozen=True)
class BarcodeHit:
    record: OligoRecord
    mismatches: int
    offset: int
    orientation: str


def open_text(path: Path, mode: str) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t")  # type: ignore[return-value]
    return path.open(mode)


def fastq_records(handle: TextIO) -> Iterator[tuple[str, str, str, str]]:
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("Truncated FASTQ record")
        yield header.rstrip("\n"), seq.rstrip("\n").upper(), plus.rstrip("\n"), qual.rstrip("\n")


def read_token(header: str) -> str:
    return header[1:].split()[0] if header.startswith("@") else header.split()[0]


def normalize_read_id(header: str) -> str:
    token = read_token(header)
    return re.sub(r"([/._ -][123])$", "", token)


def sanitize(value: str) -> str:
    value = value.strip()
    value = re.sub(r"\s+", "-", value)
    value = re.sub(r"[^A-Za-z0-9_.:+-]", "-", value)
    return value.strip("-") or "NA"


def column_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    idx = 0
    for ch in letters:
        idx = idx * 26 + ord(ch.upper()) - 64
    return idx - 1


def xlsx_rows(path: Path) -> Iterator[tuple[str, int, list[str]]]:
    with zipfile.ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{XLSX_NS}si"):
                shared_strings.append("".join(t.text or "" for t in item.iter(f"{XLSX_NS}t")))

        rel_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rels = {
            rel.attrib["Id"]: rel.attrib["Target"].lstrip("/")
            for rel in rel_root.findall(f"{PKG_REL_NS}Relationship")
        }
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))

        for sheet in workbook.findall(f".//{XLSX_NS}sheet"):
            sheet_name = sheet.attrib["name"]
            target = rels[sheet.attrib[f"{XLSX_REL_NS}id"]]
            sheet_path = target if target.startswith("xl/") else f"xl/{target}"
            worksheet = ET.fromstring(zf.read(sheet_path))

            for row in worksheet.findall(f".//{XLSX_NS}sheetData/{XLSX_NS}row"):
                values: dict[int, str] = {}
                max_col = -1
                for cell in row.findall(f"{XLSX_NS}c"):
                    idx = column_index(cell.attrib.get("r", "A1"))
                    cell_type = cell.attrib.get("t")
                    value = ""
                    v = cell.find(f"{XLSX_NS}v")
                    if v is not None and v.text is not None:
                        value = v.text
                        if cell_type == "s":
                            value = shared_strings[int(value)]
                    else:
                        inline = cell.find(f"{XLSX_NS}is")
                        if inline is not None:
                            value = "".join(t.text or "" for t in inline.iter(f"{XLSX_NS}t"))
                    values[idx] = value.strip()
                    max_col = max(max_col, idx)
                if max_col >= 0:
                    yield sheet_name, int(row.attrib.get("r", "0")), [
                        values.get(i, "") for i in range(max_col + 1)
                    ]


def load_oligos(path: Path) -> list[OligoRecord]:
    by_sheet: dict[str, list[tuple[int, list[str]]]] = defaultdict(list)
    for sheet, row_number, values in xlsx_rows(path):
        by_sheet[sheet].append((row_number, values))

    records: list[OligoRecord] = []
    seen_barcodes: Counter[str] = Counter()
    for sheet, rows in by_sheet.items():
        if not rows:
            continue
        header = [v.strip() for v in rows[0][1]]
        header_lc = [h.lower() for h in header]
        if "barcode" not in header_lc:
            continue
        barcode_idx = header_lc.index("barcode")

        for row_number, values in rows[1:]:
            if barcode_idx >= len(values):
                continue
            barcode = re.sub(r"[^ACGTN]", "", values[barcode_idx].upper())
            if not barcode:
                continue
            metadata = {
                header[i] or f"col{i + 1}": values[i]
                for i in range(min(len(header), len(values)))
                if values[i] != ""
            }
            source_id = ""
            for key in ("oligo_id", "name", "TSS#", "Synth TSS inds"):
                if key in metadata:
                    source_id = metadata[key]
                    break
            if not source_id:
                source_id = f"row{row_number}"
            oligo_id = sanitize(f"{sheet}:{source_id}:{barcode}")
            records.append(OligoRecord(oligo_id, sheet, row_number, barcode, metadata))
            seen_barcodes[barcode] += 1

    duplicated = [bc for bc, n in seen_barcodes.items() if n > 1]
    if duplicated:
        joined = ", ".join(duplicated[:10])
        raise ValueError(f"Duplicate oligo barcodes in metadata; first duplicates: {joined}")
    return records


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def neighbors_one_mismatch(seq: str) -> Iterator[str]:
    for i, base in enumerate(seq):
        for alt in DNA_ALPHABET:
            if alt != base:
                yield f"{seq[:i]}{alt}{seq[i + 1:]}"


def build_barcode_maps(
    records: list[OligoRecord], max_mismatches: int, orientations: set[str]
) -> dict[int, dict[str, BarcodeHit | None]]:
    maps: dict[int, dict[str, BarcodeHit | None]] = defaultdict(dict)
    for record in records:
        oriented = []
        if "forward" in orientations:
            oriented.append((record.barcode, "forward"))
        if "reverse-complement" in orientations:
            oriented.append((reverse_complement(record.barcode), "reverse-complement"))
        for barcode, orientation in oriented:
            length = len(barcode)
            hit = BarcodeHit(record, 0, 0, orientation)
            if barcode in maps[length] and maps[length][barcode] != hit:
                maps[length][barcode] = None
            else:
                maps[length][barcode] = hit
            if max_mismatches >= 1:
                for neighbor in neighbors_one_mismatch(barcode):
                    mm_hit = BarcodeHit(record, 1, 0, orientation)
                    if neighbor in maps[length] and maps[length][neighbor] != mm_hit:
                        maps[length][neighbor] = None
                    else:
                        maps[length][neighbor] = mm_hit
    return maps


def find_barcode(
    seq: str,
    barcode_maps: dict[int, dict[str, BarcodeHit | None]],
    search_bases: int,
) -> BarcodeHit | None | str:
    hits: list[BarcodeHit] = []
    for length, lookup in barcode_maps.items():
        last_start = min(search_bases, max(0, len(seq) - length))
        for offset in range(last_start + 1):
            query = seq[offset : offset + length]
            hit = lookup.get(query)
            if hit is None and query in lookup:
                return "ambiguous"
            if hit is not None:
                hits.append(BarcodeHit(hit.record, hit.mismatches, offset, hit.orientation))
        if hits:
            break
    if not hits:
        return None
    hits.sort(key=lambda h: (h.offset, h.mismatches, h.record.oligo_id))
    best = hits[0]
    ties = [h for h in hits if (h.offset, h.mismatches) == (best.offset, best.mismatches)]
    if len({(h.record.oligo_id, h.orientation) for h in ties}) > 1:
        return "ambiguous"
    return best


def write_metadata(records: list[OligoRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({key for record in records for key in record.metadata})
    with path.open("w") as out:
        out.write("\t".join(["oligo_id", "sheet", "row_number", "barcode", *keys]) + "\n")
        for record in records:
            row = [
                record.oligo_id,
                record.sheet,
                str(record.row_number),
                record.barcode,
                *[record.metadata.get(key, "") for key in keys],
            ]
            out.write("\t".join(row) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--r3", required=True, type=Path)
    parser.add_argument("--oligo-xlsx", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r3", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--metadata-out", required=True, type=Path)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--max-barcode-mismatches", type=int, default=1, choices=[0, 1])
    parser.add_argument("--barcode-search-bases", type=int, default=0)
    parser.add_argument(
        "--barcode-orientation",
        choices=["forward", "reverse-complement", "both"],
        default="both",
    )
    r3_barcode_group = parser.add_mutually_exclusive_group()
    r3_barcode_group.add_argument(
        "--keep-r3-barcode",
        dest="keep_r3_barcode",
        action="store_true",
        default=True,
        help="Keep the matched barcode in R3. This is the default.",
    )
    r3_barcode_group.add_argument(
        "--trim-r3-barcode",
        dest="keep_r3_barcode",
        action="store_false",
        help="Trim the matched barcode and any leading offset from R3.",
    )
    args = parser.parse_args()

    records = load_oligos(args.oligo_xlsx)
    if not records:
        raise ValueError(f"No oligo barcode records loaded from {args.oligo_xlsx}")
    write_metadata(records, args.metadata_out)

    orientations = (
        {"forward", "reverse-complement"}
        if args.barcode_orientation == "both"
        else {args.barcode_orientation}
    )
    barcode_maps = build_barcode_maps(records, args.max_barcode_mismatches, orientations)

    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r3.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)

    stats: Counter[str] = Counter()
    with open_text(args.r1, "r") as r1_fh, open_text(args.r2, "r") as r2_fh, open_text(
        args.r3, "r"
    ) as r3_fh, open_text(args.out_r1, "w") as out_r1, open_text(args.out_r3, "w") as out_r3:
        for idx, (rec1, rec2, rec3) in enumerate(
            zip(fastq_records(r1_fh), fastq_records(r2_fh), fastq_records(r3_fh)), start=1
        ):
            if args.max_reads and idx > args.max_reads:
                break
            stats["total_triplets"] += 1

            h1, s1, _p1, q1 = rec1
            h2, s2, _p2, _q2 = rec2
            h3, s3, _p3, q3 = rec3

            if normalize_read_id(h1) != normalize_read_id(h2) or normalize_read_id(h1) != normalize_read_id(h3):
                stats["discard_read_name_mismatch"] += 1
                continue
            if len(s1) < 9 or len(q1) < 9:
                stats["discard_r1_too_short"] += 1
                continue
            if len(s2) < 17:
                stats["discard_r2_too_short"] += 1
                continue

            r1_umi = s1[:8]
            meth_code = s2[:3]
            r2_umi = s2[3:17]
            meth = METHYLATION.get(meth_code)
            if meth is None:
                stats[f"discard_unknown_meth_{meth_code}"] += 1
                continue

            barcode_hit = find_barcode(s3, barcode_maps, args.barcode_search_bases)
            if barcode_hit is None:
                stats["discard_no_oligo_barcode"] += 1
                continue
            if barcode_hit == "ambiguous":
                stats["discard_ambiguous_oligo_barcode"] += 1
                continue
            assert isinstance(barcode_hit, BarcodeHit)

            combined_umi = f"{r1_umi}{r2_umi}"
            original = read_token(h1)
            new_header = f"@{original}|METH={meth}|OLIGO={barcode_hit.record.oligo_id}_{combined_umi}"

            trim_r3 = 0 if args.keep_r3_barcode else barcode_hit.offset + len(barcode_hit.record.barcode)
            trimmed_r1_seq = s1[8:]
            trimmed_r1_qual = q1[8:]
            trimmed_r3_seq = s3[trim_r3:]
            trimmed_r3_qual = q3[trim_r3:]
            if not trimmed_r1_seq or not trimmed_r3_seq:
                stats["discard_empty_after_trim"] += 1
                continue

            out_r1.write(f"{new_header}\n{trimmed_r1_seq}\n+\n{trimmed_r1_qual}\n")
            out_r3.write(f"{new_header}\n{trimmed_r3_seq}\n+\n{trimmed_r3_qual}\n")

            stats["written_triplets"] += 1
            stats[f"meth_{meth_code}_{meth}"] += 1
            stats[f"barcode_orientation_{barcode_hit.orientation}"] += 1
            stats[f"barcode_mismatches_{barcode_hit.mismatches}"] += 1
            stats[f"barcode_offset_{barcode_hit.offset}"] += 1

            if idx % 1_000_000 == 0:
                print(
                    f"processed={idx} written={stats['written_triplets']} "
                    f"no_barcode={stats['discard_no_oligo_barcode']}",
                    file=sys.stderr,
                    flush=True,
                )

    with args.stats.open("w") as out:
        out.write("metric\tcount\n")
        out.write(f"oligo_records_loaded\t{len(records)}\n")
        for key, value in sorted(stats.items()):
            out.write(f"{key}\t{value}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
