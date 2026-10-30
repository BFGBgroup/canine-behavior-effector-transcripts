#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence


CLASSIFICATION_REQUIRED_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "class_code",
    "novelty_category",
]

BED_COLUMNS = [
    ("1", "chr", "GTF seqname"),
    ("2", "exon_start", "0-based BED start"),
    ("3", "exon_end", "0-based BED end, half-open"),
    ("4", "query_transcript_id", "query transcript_id from gffcompare annotated GTF"),
    ("5", "query_gene_id", "query gene_id joined from gffcompare classification"),
    ("6", "strand", "GTF strand"),
    ("7", "class_code", "gffcompare class_code from tmap"),
    ("8", "novelty_category", "workflow-level class-code interpretation"),
]

QC_FIELDS = ["metric", "value", "note"]
MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]
COLUMNS_FIELDS = ["column_number", "field_name", "description"]


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return text if text else "NA"


def fail(message: str) -> None:
    raise RuntimeError(message)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract transcript exon BED from gffcompare annotated GTF and transcript classification."
    )
    parser.add_argument("--annotated-gtf", required=True, type=Path)
    parser.add_argument("--classification-tsv", required=True, type=Path)
    parser.add_argument("--out-bed", required=True, type=Path)
    parser.add_argument("--qc-summary", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--columns-tsv", required=True, type=Path)
    return parser.parse_args(argv)


def require_readable(path: Path, label: str) -> None:
    if not path.exists():
        fail(f"{label} does not exist: {path}")
    if not path.is_file():
        fail(f"{label} is not a regular file: {path}")
    if path.stat().st_size == 0:
        fail(f"{label} is empty: {path}")


def require_fields(path: Path, fieldnames: Sequence[str], required: Sequence[str]) -> None:
    missing = [field for field in required if field not in fieldnames]
    if missing:
        fail(f"{path} missing required fields: {','.join(missing)}")


def parse_gtf_attributes(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for key, quoted_value in re.findall(r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"', text):
        attrs[key] = quoted_value
    for chunk in text.rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk or " " not in chunk:
            continue
        key, value = chunk.split(None, 1)
        attrs.setdefault(key, value.strip().strip('"'))
    return attrs


def read_classification(path: Path) -> Dict[str, Dict[str, str]]:
    logging.info("read transcript classification: %s", path)
    mapping: Dict[str, Dict[str, str]] = {}
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames or [], CLASSIFICATION_REQUIRED_FIELDS)
        for row_number, row in enumerate(reader, start=1):
            transcript_id = clean(row.get("query_transcript_id"))
            if transcript_id == "NA":
                fail(f"{path} row {row_number} query_transcript_id is empty")
            if transcript_id in mapping:
                fail(f"{path} in query_transcript_id duplicate: {transcript_id}")
            mapping[transcript_id] = {
                "query_gene_id": clean(row.get("query_gene_id")),
                "class_code": clean(row.get("class_code")),
                "novelty_category": clean(row.get("novelty_category")),
            }
            if row_number % 50000 == 0:
                logging.info("read classification rows: %s", row_number)
    logging.info("classification read completed: %s transcripts", len(mapping))
    return mapping


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean(row.get(field)) for field in fieldnames})


def manifest_row(item: str, path: Path, note: str) -> Dict[str, str]:
    return {
        "item": item,
        "path": str(path),
        "exists": "1" if path.exists() else "0",
        "size_bytes": str(path.stat().st_size) if path.exists() and path.is_file() else "NA",
        "note": note,
    }


def extract_exons(
    annotated_gtf: Path,
    classification: Mapping[str, Mapping[str, str]],
    out_bed: Path,
) -> Dict[str, int]:
    logging.info("from annotated GTF extract exon BED: %s", annotated_gtf)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    tmp_bed = out_bed.with_suffix(out_bed.suffix + ".tmp")

    exon_rows = 0
    invalid_exon_rows = 0
    missing_classification_transcripts = set()
    transcripts_with_exon = set()

    with annotated_gtf.open("r", errors="replace") as input_handle, tmp_bed.open("w", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        for line_number, line in enumerate(input_handle, start=1):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom, _, _, start_text, end_text, _, strand, _, attr_text = parts[:9]
            try:
                gtf_start = int(start_text)
                gtf_end = int(end_text)
            except ValueError:
                invalid_exon_rows += 1
                continue
            if gtf_start < 1 or gtf_end < gtf_start:
                invalid_exon_rows += 1
                continue
            attrs = parse_gtf_attributes(attr_text)
            transcript_id = clean(attrs.get("transcript_id"))
            if transcript_id == "NA":
                invalid_exon_rows += 1
                continue
            info = classification.get(transcript_id)
            if info is None:
                missing_classification_transcripts.add(transcript_id)
                continue
            bed_start = gtf_start - 1
            bed_end = gtf_end
            writer.writerow(
                [
                    chrom,
                    bed_start,
                    bed_end,
                    transcript_id,
                    info["query_gene_id"],
                    strand,
                    info["class_code"],
                    info["novelty_category"],
                ]
            )
            exon_rows += 1
            transcripts_with_exon.add(transcript_id)
            if exon_rows % 250000 == 0:
                logging.info("written exon BED rows: %s", exon_rows)

    if exon_rows == 0:
        tmp_bed.unlink(missing_ok=True)
        fail("not written any exon BED row, please check annotated GTF whether contains exon feature")
    if invalid_exon_rows > 0:
        tmp_bed.unlink(missing_ok=True)
        fail(f"The annotated GTF contains invalid exon rows: {invalid_exon_rows}")
    if missing_classification_transcripts:
        tmp_bed.unlink(missing_ok=True)
        examples = ",".join(sorted(missing_classification_transcripts)[:10])
        fail(
            "Some transcript_id values in the annotated GTF have no classification; "
            f"missing_count={len(missing_classification_transcripts)} examples={examples}"
        )

    tmp_bed.replace(out_bed)
    return {
        "exon_rows_written": exon_rows,
        "transcripts_with_exon_written": len(transcripts_with_exon),
        "invalid_exon_rows": invalid_exon_rows,
        "missing_classification_transcript_count": len(missing_classification_transcripts),
    }


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    setup_logging()
    logging.info("Step 1.6 transcript exon BED extract start")

    require_readable(args.annotated_gtf, "gffcompare annotated GTF")
    require_readable(args.classification_tsv, "transcript classification TSV")

    classification = read_classification(args.classification_tsv)
    metrics = extract_exons(args.annotated_gtf, classification, args.out_bed)

    args.qc_summary.parent.mkdir(parents=True, exist_ok=True)
    qc_rows = [
        {"metric": "classification_transcript_count", "value": str(len(classification)), "note": str(args.classification_tsv)},
        {"metric": "exon_rows_written", "value": str(metrics["exon_rows_written"]), "note": str(args.out_bed)},
        {
            "metric": "transcripts_with_exon_written",
            "value": str(metrics["transcripts_with_exon_written"]),
            "note": "unique transcript_id values in BED",
        },
        {
            "metric": "missing_classification_transcript_count",
            "value": str(metrics["missing_classification_transcript_count"]),
            "note": "strictly required to be 0",
        },
        {"metric": "invalid_exon_rows", "value": str(metrics["invalid_exon_rows"]), "note": "strictly required to be 0"},
        {"metric": "bed_coordinate_system", "value": "0based_half_open", "note": "BED start = GTF start - 1, BED end = GTF end"},
        {"metric": "pass_status", "value": "PASS", "note": "PASS requires non-empty BED and complete classification join"},
    ]
    write_tsv(args.qc_summary, QC_FIELDS, qc_rows)

    column_rows = [
        {"column_number": number, "field_name": field, "description": description}
        for number, field, description in BED_COLUMNS
    ]
    write_tsv(args.columns_tsv, COLUMNS_FIELDS, column_rows)

    manifest_rows = [
        manifest_row("annotated_gtf_input", args.annotated_gtf, "gffcompare annotated/combined GTF"),
        manifest_row("classification_tsv_input", args.classification_tsv, "Step 1.6 transcript-level classification"),
        manifest_row("transcript_exon_bed", args.out_bed, "headerless 8-column BED for Step 1.7 anchor overlap"),
        manifest_row("qc_summary_tsv", args.qc_summary, "transcript exon BED QC summary"),
        manifest_row("columns_tsv", args.columns_tsv, "transcript exon BED column definitions"),
        {
            "item": "manifest_tsv",
            "path": str(args.manifest),
            "exists": "1",
            "size_bytes": "NA",
            "note": "this manifest; size is finalized after this row is written",
        },
    ]
    write_tsv(args.manifest, MANIFEST_FIELDS, manifest_rows)

    logging.info("Step 1.6 transcript exon BED extract completed: %s", args.out_bed)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except Exception as exc:
        logging.exception("Step 1.6 transcript exon BED extract failed: %s", exc)
        raise SystemExit(1)
