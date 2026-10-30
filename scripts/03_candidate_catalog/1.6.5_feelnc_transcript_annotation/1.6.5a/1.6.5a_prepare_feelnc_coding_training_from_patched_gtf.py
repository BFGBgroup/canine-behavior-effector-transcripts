#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple


QC_FIELDS = ["metric", "value", "note"]
MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%F %T",
        stream=sys.stderr,
    )


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return text if text else "NA"


def fail(message: str) -> None:
    raise RuntimeError(message)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare FEELnc protein-coding training GTF from patched CanFam3.1 GTF."
    )
    parser.add_argument("--patched-gtf", required=True, type=Path)
    parser.add_argument("--out-training-gtf", required=True, type=Path)
    parser.add_argument("--out-transcript-ids", required=True, type=Path)
    parser.add_argument("--qc-summary", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--min-coding-transcripts", type=int, default=1000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def require_readable(path: Path, label: str) -> None:
    if not path.exists():
        fail(f"{label} does not exist: {path}")
    if not path.is_file():
        fail(f"{label} is not a regular file: {path}")
    if path.stat().st_size == 0:
        fail(f"{label} is empty: {path}")


def prepare_outputs(paths: Sequence[Path], overwrite: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        joined = "\n".join(str(path) for path in existing)
        fail("output already exists; set --overwrite only after confirmation:\n" + joined)
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)


def parse_gtf_attributes(attribute_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for raw_part in attribute_text.strip().rstrip(";").split(";"):
        part = raw_part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, value = part.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def parse_gtf_line(line: str, row_number: int, path: Path) -> Tuple[List[str], Dict[str, str]]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 9:
        fail(f"{path} line {row_number} is not a 9-column GTF row")
    attrs = parse_gtf_attributes(fields[8])
    return fields, attrs


def scan_cds_transcript_ids(path: Path) -> Tuple[Set[str], Counter]:
    logging.info("stage 1/4 scan patched GTF in CDS transcript_id: %s", path)
    coding_transcript_ids: Set[str] = set()
    metrics: Counter = Counter()

    with path.open("r", errors="replace") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields, attrs = parse_gtf_line(line, row_number, path)
            feature = fields[2]
            metrics["input_feature_rows"] += 1
            metrics[f"input_{feature}_rows"] += 1
            if "gene_id" not in attrs:
                metrics["input_missing_gene_id_rows"] += 1
            transcript_id = clean(attrs.get("transcript_id"))
            if transcript_id == "NA":
                metrics["input_missing_transcript_id_rows"] += 1
            if feature == "CDS" and transcript_id != "NA":
                coding_transcript_ids.add(transcript_id)
            if row_number % 250000 == 0:
                logging.info("stage 1/4 scanned GTF row count: %s", row_number)

    metrics["coding_transcript_id_count"] = len(coding_transcript_ids)
    logging.info("stage 1/4 completed, CDS-containing transcript count: %s", len(coding_transcript_ids))
    return coding_transcript_ids, metrics


def write_transcript_ids(path: Path, transcript_ids: Iterable[str]) -> int:
    logging.info("stage 2/4 write protein-coding transcript ID list: %s", path)
    count = 0
    with path.open("w", newline="") as handle:
        for transcript_id in sorted(transcript_ids):
            handle.write(transcript_id + "\n")
            count += 1
    return count


def seqname_group(seqname: str) -> str:
    if seqname.startswith("chr"):
        return "chr_like"
    if seqname.startswith("NW_"):
        return "NW"
    return "other"


def extract_training_gtf(path: Path, output_gtf: Path, coding_transcript_ids: Set[str]) -> Counter:
    logging.info("stage 3/4 extract protein-coding training GTF: %s", output_gtf)
    metrics: Counter = Counter()

    with path.open("r", errors="replace") as input_handle, output_gtf.open("w", newline="") as output_handle:
        for row_number, line in enumerate(input_handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields, attrs = parse_gtf_line(line, row_number, path)
            transcript_id = clean(attrs.get("transcript_id"))
            if transcript_id not in coding_transcript_ids:
                continue

            output_handle.write(line)
            feature = fields[2]
            group = seqname_group(fields[0])
            metrics["output_rows"] += 1
            metrics[f"output_{feature}_rows"] += 1
            metrics[f"output_seqname_{group}_rows"] += 1
            if "gene_id" not in attrs:
                metrics["output_missing_gene_id_rows"] += 1
            if row_number % 250000 == 0:
                logging.info("stage 3/4 scanned GTF row count: %s; written rows: %s", row_number, metrics["output_rows"])

    logging.info("stage 3/4 completed, training GTF rows: %s", metrics["output_rows"])
    return metrics


def manifest_row(item: str, path: Path, note: str) -> Dict[str, str]:
    exists = path.exists()
    return {
        "item": item,
        "path": str(path),
        "exists": "1" if exists else "0",
        "size_bytes": str(path.stat().st_size) if exists and path.is_file() else "NA",
        "note": note,
    }


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_qc_rows(
    args: argparse.Namespace,
    scan_metrics: Counter,
    output_metrics: Counter,
    id_rows_written: int,
    start_time: datetime,
) -> List[Dict[str, str]]:
    pass_status = (
        scan_metrics["coding_transcript_id_count"] >= args.min_coding_transcripts
        and id_rows_written == scan_metrics["coding_transcript_id_count"]
        and output_metrics["output_rows"] > 0
        and output_metrics["output_transcript_rows"] == scan_metrics["coding_transcript_id_count"]
        and output_metrics["output_exon_rows"] > 0
        and output_metrics["output_CDS_rows"] > 0
        and output_metrics["output_missing_gene_id_rows"] == 0
    )
    elapsed_sec = int((datetime.now() - start_time).total_seconds())

    rows = [
        {"metric": "patched_gtf_input", "value": str(args.patched_gtf), "note": "source GTF with patched gene_id"},
        {"metric": "input_feature_rows", "value": str(scan_metrics["input_feature_rows"]), "note": "non-comment GTF rows scanned"},
        {"metric": "input_transcript_rows", "value": str(scan_metrics["input_transcript_rows"]), "note": "feature == transcript"},
        {"metric": "input_exon_rows", "value": str(scan_metrics["input_exon_rows"]), "note": "feature == exon"},
        {"metric": "input_CDS_rows", "value": str(scan_metrics["input_CDS_rows"]), "note": "feature == CDS"},
        {"metric": "input_missing_gene_id_rows", "value": str(scan_metrics["input_missing_gene_id_rows"]), "note": "should be 0 for patched GTF"},
        {"metric": "input_missing_transcript_id_rows", "value": str(scan_metrics["input_missing_transcript_id_rows"]), "note": "rows without transcript_id"},
        {"metric": "coding_transcript_id_count", "value": str(scan_metrics["coding_transcript_id_count"]), "note": "unique transcript_id values observed on CDS rows"},
        {"metric": "transcript_id_rows_written", "value": str(id_rows_written), "note": str(args.out_transcript_ids)},
        {"metric": "training_gtf_rows_written", "value": str(output_metrics["output_rows"]), "note": str(args.out_training_gtf)},
        {"metric": "training_gtf_transcript_rows", "value": str(output_metrics["output_transcript_rows"]), "note": "feature == transcript"},
        {"metric": "training_gtf_exon_rows", "value": str(output_metrics["output_exon_rows"]), "note": "feature == exon"},
        {"metric": "training_gtf_CDS_rows", "value": str(output_metrics["output_CDS_rows"]), "note": "feature == CDS"},
        {"metric": "training_gtf_missing_gene_id_rows", "value": str(output_metrics["output_missing_gene_id_rows"]), "note": "should be 0"},
        {"metric": "training_gtf_chr_like_rows", "value": str(output_metrics["output_seqname_chr_like_rows"]), "note": "seqname starts with chr"},
        {"metric": "training_gtf_NW_rows", "value": str(output_metrics["output_seqname_NW_rows"]), "note": "seqname starts with NW_"},
        {"metric": "training_gtf_other_seqname_rows", "value": str(output_metrics["output_seqname_other_rows"]), "note": "seqname is neither chr* nor NW_*"},
        {"metric": "runtime_sec", "value": str(elapsed_sec), "note": "wall-clock runtime of this script"},
        {"metric": "pass_status", "value": "PASS" if pass_status else "FAIL", "note": "PASS requires non-empty CDS-derived training set and no missing gene_id"},
    ]
    return rows


def build_manifest_rows(args: argparse.Namespace) -> List[Dict[str, str]]:
    return [
        manifest_row("patched_gtf_input", args.patched_gtf, "CanFam3.1 chr/NW GTF with gene_id patched"),
        manifest_row("protein_coding_transcript_ids_output", args.out_transcript_ids, "CDS-containing transcript IDs for FEELnc coding training"),
        manifest_row("protein_coding_training_gtf_output", args.out_training_gtf, "FEELnc known mRNA training GTF from patched GTF"),
        manifest_row("qc_summary_output", args.qc_summary, "QC metrics for patched-GTF-derived training set"),
        manifest_row("manifest_output", args.manifest, "File manifest for this preparation step"),
    ]


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    args = parse_args(argv)
    start_time = datetime.now()

    logging.info("Step 1.6.5a FEELnc coding training set preparation start")
    require_readable(args.patched_gtf, "patched GTF")
    if args.min_coding_transcripts < 1:
        fail("--min-coding-transcripts must be positive")
    prepare_outputs(
        [args.out_training_gtf, args.out_transcript_ids, args.qc_summary, args.manifest],
        overwrite=args.overwrite,
    )

    coding_transcript_ids, scan_metrics = scan_cds_transcript_ids(args.patched_gtf)
    if len(coding_transcript_ids) < args.min_coding_transcripts:
        fail(
            f"too few CDS-containing transcripts: {len(coding_transcript_ids)} "
            f"< {args.min_coding_transcripts}"
        )

    id_rows_written = write_transcript_ids(args.out_transcript_ids, coding_transcript_ids)
    output_metrics = extract_training_gtf(args.patched_gtf, args.out_training_gtf, coding_transcript_ids)

    logging.info("stage 4/4 write QC summary and manifest")
    write_tsv(args.qc_summary, QC_FIELDS, build_qc_rows(args, scan_metrics, output_metrics, id_rows_written, start_time))
    write_tsv(args.manifest, MANIFEST_FIELDS, build_manifest_rows(args))

    logging.info("Step 1.6.5a FEELnc coding training set preparation completed")
    logging.info("training GTF: %s", args.out_training_gtf)
    logging.info("transcript IDs: %s", args.out_transcript_ids)
    logging.info("QC summary: %s", args.qc_summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        logging.exception("Step 1.6.5a FEELnc coding training set preparation failed: %s", exc)
        raise SystemExit(1)
