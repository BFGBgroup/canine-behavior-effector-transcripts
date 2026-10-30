#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


BASIC_INFO_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "seqname",
    "multi_seqname_flag",
    "source",
    "strand",
    "multi_strand_flag",
    "transcript_start",
    "transcript_end",
    "exon_count",
    "transcript_length",
    "single_exon_flag",
    "feature_rows",
    "transcript_row_count",
    "gtf_class_code",
    "class_code",
    "reference_gene_id",
    "reference_transcript_id",
    "reference_gene_name",
    "novelty_category",
    "candidate_relevance",
    "initial_priority_for_anchor_overlap",
    "classification_num_exons",
    "classification_length",
    "FPKM",
    "TPM",
    "cov",
    "major_iso_id",
    "ref_match_len",
    "seqname_group",
    "consistency_flags",
]
CLASS_SUMMARY_FIELDS = ["class_code", "transcript_count", "single_exon_count", "multi_exon_count"]
ISSUE_EXAMPLE_FIELDS = [
    "issue_type",
    "query_transcript_id",
    "row_number",
    "query_gene_id",
    "class_code",
    "gtf_value",
    "classification_value",
    "detail",
]
QC_FIELDS = ["metric", "value", "note"]
MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]
EXPECTED_GTF_FEATURES = {"transcript", "exon"}
REQUIRED_CLASSIFICATION_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "reference_gene_id",
    "reference_transcript_id",
    "class_code",
    "novelty_category",
    "candidate_relevance",
    "initial_priority_for_anchor_overlap",
    "num_exons",
    "FPKM",
    "TPM",
    "cov",
    "length",
    "major_iso_id",
    "ref_match_len",
]


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
        description="Prepare full gffcompare query transcripts for FEELnc annotation."
    )
    parser.add_argument("--annotated-gtf", required=True, type=Path)
    parser.add_argument("--classification-tsv", required=True, type=Path)
    parser.add_argument("--out-query-gtf", required=True, type=Path)
    parser.add_argument("--out-basic-info", required=True, type=Path)
    parser.add_argument("--out-class-summary", required=True, type=Path)
    parser.add_argument("--out-issue-examples", required=True, type=Path)
    parser.add_argument("--qc-summary", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--max-issue-examples", type=int, default=200)
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
        if not part or " " not in part:
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


def seqname_group(seqname: str) -> str:
    if seqname.startswith("chr"):
        return "chr_like"
    if seqname.startswith("NW_"):
        return "NW"
    return "other"


def parse_int(text: str, label: str) -> int | None:
    if text in {"", "-", "NA"}:
        return None
    try:
        return int(text)
    except ValueError:
        logging.warning("cannot parse integer %s=%s", label, text)
        return None


def load_classification(path: Path, max_issue_examples: int) -> Tuple[Dict[str, Dict[str, str]], Counter, List[Dict[str, str]]]:
    logging.info("stage 1/5 read gffcompare transcript classification: %s", path)
    metrics: Counter = Counter()
    rows_by_tid: Dict[str, Dict[str, str]] = {}
    issue_examples: List[Dict[str, str]] = []

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"classification table has no header: {path}")
        missing_fields = [field for field in REQUIRED_CLASSIFICATION_FIELDS if field not in reader.fieldnames]
        if missing_fields:
            fail("classification table missing required fields: " + ",".join(missing_fields))

        for row_number, row in enumerate(reader, start=2):
            tid = clean(row.get("query_transcript_id"))
            if tid == "NA":
                metrics["classification_missing_transcript_id_rows"] += 1
                continue
            if tid in rows_by_tid:
                metrics["classification_duplicate_transcript_id_rows"] += 1
                if len(issue_examples) < max_issue_examples:
                    issue_examples.append(
                        {
                            "issue_type": "duplicate_classification_transcript_id",
                            "query_transcript_id": tid,
                            "row_number": str(row_number),
                            "query_gene_id": clean(row.get("query_gene_id")),
                            "class_code": clean(row.get("class_code")),
                            "gtf_value": "NA",
                            "classification_value": clean(row.get("query_transcript_id")),
                            "detail": "classification TSV contains a repeated query_transcript_id; first row is kept",
                        }
                    )
                continue
            rows_by_tid[tid] = {key: clean(value) for key, value in row.items()}
            metrics["classification_rows"] += 1
            if row_number % 25000 == 0:
                logging.info("stage 1/5 read classification row count: %s", row_number - 1)

    metrics["classification_unique_transcripts"] = len(rows_by_tid)
    logging.info("stage 1/5 completed, classification transcript count: %s", len(rows_by_tid))
    return rows_by_tid, metrics, issue_examples


def empty_info() -> MutableMapping[str, object]:
    return {
        "gene_ids": set(),
        "gene_names": set(),
        "seqnames": set(),
        "sources": set(),
        "strands": set(),
        "transcript_start": None,
        "transcript_end": None,
        "transcript_row_numbers": [],
        "exon_count": 0,
        "transcript_length": 0,
        "feature_rows": 0,
        "transcript_row_count": 0,
        "gtf_class_codes": set(),
        "ref_gene_ids": set(),
        "cmp_refs": set(),
    }


def update_transcript_info(
    info: MutableMapping[str, object],
    fields: Sequence[str],
    attrs: Mapping[str, str],
    row_number: int,
) -> None:
    feature = fields[2]
    start = int(fields[3])
    end = int(fields[4])
    if start > end:
        fail(f"GTF coordinate start > end for transcript_id={attrs.get('transcript_id', 'NA')}")

    info["feature_rows"] = int(info["feature_rows"]) + 1
    info["seqnames"].add(fields[0])
    info["sources"].add(fields[1])
    info["strands"].add(fields[6])
    info["gene_ids"].add(clean(attrs.get("gene_id")))
    if attrs.get("gene_name"):
        info["gene_names"].add(clean(attrs.get("gene_name")))
    if attrs.get("ref_gene_id"):
        info["ref_gene_ids"].add(clean(attrs.get("ref_gene_id")))
    if attrs.get("cmp_ref"):
        info["cmp_refs"].add(clean(attrs.get("cmp_ref")))

    if feature == "transcript":
        info["transcript_row_count"] = int(info["transcript_row_count"]) + 1
        info["transcript_row_numbers"].append(row_number)
        info["transcript_start"] = start
        info["transcript_end"] = end
        if attrs.get("class_code"):
            info["gtf_class_codes"].add(clean(attrs.get("class_code")))
    elif feature == "exon":
        info["exon_count"] = int(info["exon_count"]) + 1
        info["transcript_length"] = int(info["transcript_length"]) + (end - start + 1)


def scan_annotated_gtf(path: Path) -> Tuple[Dict[str, MutableMapping[str, object]], List[str], Counter]:
    logging.info("stage 2/5 scan annotated GTF and summary transcript basic information: %s", path)
    metrics: Counter = Counter()
    info_by_tid: Dict[str, MutableMapping[str, object]] = defaultdict(empty_info)
    transcript_order: List[str] = []
    seen_tids: Set[str] = set()

    with path.open("r", errors="replace") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields, attrs = parse_gtf_line(line, row_number, path)
            feature = fields[2]
            metrics["gtf_feature_rows"] += 1
            metrics[f"gtf_{feature}_rows"] += 1
            if feature not in EXPECTED_GTF_FEATURES:
                metrics["gtf_unexpected_feature_rows"] += 1
                metrics[f"gtf_unexpected_{feature}_rows"] += 1

            if "transcript_id" not in attrs:
                metrics["gtf_missing_transcript_id_rows"] += 1
                continue
            if "gene_id" not in attrs:
                metrics["gtf_missing_gene_id_rows"] += 1

            tid = clean(attrs.get("transcript_id"))
            if tid == "NA":
                metrics["gtf_missing_transcript_id_rows"] += 1
                continue
            if tid not in seen_tids:
                transcript_order.append(tid)
                seen_tids.add(tid)
            update_transcript_info(info_by_tid[tid], fields, attrs, row_number)

            if row_number % 250000 == 0:
                logging.info("stage 2/5 scanned GTF row count: %s", row_number)

    metrics["gtf_unique_transcripts"] = len(info_by_tid)
    logging.info("stage 2/5 completed, GTF transcript count: %s", len(info_by_tid))
    return info_by_tid, transcript_order, metrics


def joined(values: Iterable[object]) -> str:
    cleaned = sorted(clean(value) for value in values if clean(value) != "NA")
    return ";".join(cleaned) if cleaned else "NA"


def build_basic_rows(
    transcript_order: Sequence[str],
    info_by_tid: Mapping[str, MutableMapping[str, object]],
    classification: Mapping[str, Mapping[str, str]],
    max_issue_examples: int,
) -> Tuple[List[Dict[str, str]], Counter, List[Dict[str, str]]]:
    logging.info("stage 3/5 merge GTF basic information and classification annotation")
    rows: List[Dict[str, str]] = []
    metrics: Counter = Counter()
    issue_examples: List[Dict[str, str]] = []

    for tid in transcript_order:
        info = info_by_tid[tid]
        class_row = classification.get(tid)
        flags: List[str] = []

        if class_row is None:
            metrics["gtf_transcripts_missing_classification"] += 1
            flags.append("missing_classification")
            class_code = "NA"
            reference_gene_id = joined(info["ref_gene_ids"])
            reference_transcript_id = joined(info["cmp_refs"])
        else:
            class_code = clean(class_row.get("class_code"))
            reference_gene_id = clean(class_row.get("reference_gene_id"))
            reference_transcript_id = clean(class_row.get("reference_transcript_id"))

        gtf_class_code = joined(info["gtf_class_codes"])
        if gtf_class_code != "NA" and class_code != "NA" and gtf_class_code != class_code:
            metrics["class_code_mismatch_transcripts"] += 1
            flags.append("class_code_mismatch")

        exon_count = int(info["exon_count"])
        transcript_length = int(info["transcript_length"])
        seqname_count = len(info["seqnames"])
        strand_count = len(info["strands"])
        multi_seqname_flag = "1" if seqname_count > 1 else "0"
        multi_strand_flag = "1" if strand_count > 1 else "0"
        if multi_seqname_flag == "1":
            metrics["multi_seqname_transcripts"] += 1
            flags.append("multi_seqname")
        if multi_strand_flag == "1":
            metrics["multi_strand_transcripts"] += 1
            flags.append("multi_strand")
        if int(info["transcript_row_count"]) != 1:
            metrics["transcripts_without_single_transcript_row"] += 1
            flags.append("transcript_row_count_not_1")
            if len(issue_examples) < max_issue_examples:
                issue_examples.append(
                    {
                        "issue_type": "duplicate_or_missing_gtf_transcript_feature",
                        "query_transcript_id": tid,
                        "row_number": joined(info["transcript_row_numbers"]),
                        "query_gene_id": clean(class_row.get("query_gene_id")) if class_row else joined(info["gene_ids"]),
                        "class_code": clean(class_row.get("class_code")) if class_row else "NA",
                        "gtf_value": clean(info["transcript_row_count"]),
                        "classification_value": "1",
                        "detail": "annotated GTF should contain exactly one transcript feature row for each query_transcript_id",
                    }
                )
        if exon_count < 1:
            metrics["transcripts_without_exon"] += 1
            flags.append("no_exon")

        if class_row is not None:
            class_num_exons = parse_int(clean(class_row.get("num_exons")), "num_exons")
            class_length = parse_int(clean(class_row.get("length")), "length")
            if class_num_exons is not None and class_num_exons != exon_count:
                metrics["exon_count_mismatch_transcripts"] += 1
                flags.append("exon_count_mismatch")
            if class_length is not None and class_length != transcript_length:
                metrics["length_mismatch_transcripts"] += 1
                flags.append("length_mismatch")
                if len(issue_examples) < max_issue_examples:
                    issue_examples.append(
                        {
                            "issue_type": "length_mismatch",
                            "query_transcript_id": tid,
                            "row_number": clean(class_row.get("tmap_row_number")),
                            "query_gene_id": clean(class_row.get("query_gene_id")),
                            "class_code": clean(class_row.get("class_code")),
                            "gtf_value": str(transcript_length),
                            "classification_value": clean(class_row.get("length")),
                            "detail": "exon-summed GTF transcript_length differs from classification length",
                        }
                    )

        seqname = joined(info["seqnames"])
        seqname_group_value = "multi" if multi_seqname_flag == "1" else seqname_group(seqname)
        row = {
            "query_transcript_id": tid,
            "query_gene_id": clean(class_row.get("query_gene_id")) if class_row else joined(info["gene_ids"]),
            "seqname": seqname,
            "multi_seqname_flag": multi_seqname_flag,
            "source": joined(info["sources"]),
            "strand": joined(info["strands"]),
            "multi_strand_flag": multi_strand_flag,
            "transcript_start": clean(info["transcript_start"]),
            "transcript_end": clean(info["transcript_end"]),
            "exon_count": str(exon_count),
            "transcript_length": str(transcript_length),
            "single_exon_flag": "1" if exon_count == 1 else "0",
            "feature_rows": str(info["feature_rows"]),
            "transcript_row_count": str(info["transcript_row_count"]),
            "gtf_class_code": gtf_class_code,
            "class_code": class_code,
            "reference_gene_id": reference_gene_id,
            "reference_transcript_id": reference_transcript_id,
            "reference_gene_name": joined(info["gene_names"]),
            "novelty_category": clean(class_row.get("novelty_category")) if class_row else "NA",
            "candidate_relevance": clean(class_row.get("candidate_relevance")) if class_row else "NA",
            "initial_priority_for_anchor_overlap": clean(class_row.get("initial_priority_for_anchor_overlap")) if class_row else "NA",
            "classification_num_exons": clean(class_row.get("num_exons")) if class_row else "NA",
            "classification_length": clean(class_row.get("length")) if class_row else "NA",
            "FPKM": clean(class_row.get("FPKM")) if class_row else "NA",
            "TPM": clean(class_row.get("TPM")) if class_row else "NA",
            "cov": clean(class_row.get("cov")) if class_row else "NA",
            "major_iso_id": clean(class_row.get("major_iso_id")) if class_row else "NA",
            "ref_match_len": clean(class_row.get("ref_match_len")) if class_row else "NA",
            "seqname_group": seqname_group_value,
            "consistency_flags": ";".join(flags) if flags else "ok",
        }
        rows.append(row)

    classification_missing = sorted(set(classification) - set(info_by_tid))
    metrics["classification_transcripts_missing_gtf"] = len(classification_missing)
    metrics["basic_info_rows"] = len(rows)
    return rows, metrics, issue_examples


def copy_gtf(input_path: Path, output_path: Path) -> Counter:
    logging.info("stage 4/5 write all records FEELnc query GTF: %s", output_path)
    metrics: Counter = Counter()
    with input_path.open("r", errors="replace") as input_handle, output_path.open("w", newline="") as output_handle:
        for row_number, line in enumerate(input_handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields, attrs = parse_gtf_line(line, row_number, input_path)
            if clean(attrs.get("transcript_id")) == "NA":
                continue
            output_handle.write(line)
            metrics["output_query_gtf_rows"] += 1
            metrics[f"output_query_gtf_{fields[2]}_rows"] += 1
            if fields[2] not in EXPECTED_GTF_FEATURES:
                metrics["output_query_gtf_unexpected_feature_rows"] += 1
            if row_number % 250000 == 0:
                logging.info("stage 4/5 scanned GTF row count: %s; written rows: %s", row_number, metrics["output_query_gtf_rows"])
    return metrics


def class_summary_rows(basic_rows: Sequence[Mapping[str, str]]) -> List[Dict[str, str]]:
    summary: Dict[str, Counter] = defaultdict(Counter)
    for row in basic_rows:
        class_code = row["class_code"]
        summary[class_code]["transcript_count"] += 1
        if row["single_exon_flag"] == "1":
            summary[class_code]["single_exon_count"] += 1
        else:
            summary[class_code]["multi_exon_count"] += 1

    return [
        {
            "class_code": class_code,
            "transcript_count": str(counts["transcript_count"]),
            "single_exon_count": str(counts["single_exon_count"]),
            "multi_exon_count": str(counts["multi_exon_count"]),
        }
        for class_code, counts in sorted(summary.items(), key=lambda item: (-item[1]["transcript_count"], item[0]))
    ]


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def manifest_row(item: str, path: Path, note: str) -> Dict[str, str]:
    exists = path.exists()
    return {
        "item": item,
        "path": str(path),
        "exists": "1" if exists else "0",
        "size_bytes": str(path.stat().st_size) if exists and path.is_file() else "NA",
        "note": note,
    }


def build_qc_rows(
    args: argparse.Namespace,
    classification_metrics: Counter,
    gtf_metrics: Counter,
    merge_metrics: Counter,
    copy_metrics: Counter,
    start_time: datetime,
) -> List[Dict[str, str]]:
    pass_status = (
        classification_metrics["classification_duplicate_transcript_id_rows"] == 0
        and classification_metrics["classification_missing_transcript_id_rows"] == 0
        and gtf_metrics["gtf_missing_transcript_id_rows"] == 0
        and gtf_metrics["gtf_missing_gene_id_rows"] == 0
        and gtf_metrics["gtf_unexpected_feature_rows"] == 0
        and gtf_metrics["gtf_unique_transcripts"] == classification_metrics["classification_unique_transcripts"]
        and merge_metrics["gtf_transcripts_missing_classification"] == 0
        and merge_metrics["classification_transcripts_missing_gtf"] == 0
        and merge_metrics["class_code_mismatch_transcripts"] == 0
        and merge_metrics["exon_count_mismatch_transcripts"] == 0
        and merge_metrics["length_mismatch_transcripts"] == 0
        and merge_metrics["multi_seqname_transcripts"] == 0
        and merge_metrics["multi_strand_transcripts"] == 0
        and copy_metrics["output_query_gtf_rows"] == gtf_metrics["gtf_feature_rows"]
        and copy_metrics["output_query_gtf_unexpected_feature_rows"] == 0
    )
    elapsed_sec = int((datetime.now() - start_time).total_seconds())

    return [
        {"metric": "annotated_gtf_input", "value": str(args.annotated_gtf), "note": "full gffcompare annotated query GTF"},
        {"metric": "classification_tsv_input", "value": str(args.classification_tsv), "note": "Step 1.6 transcript classification table"},
        {"metric": "classification_rows", "value": str(classification_metrics["classification_rows"]), "note": "non-duplicate rows loaded"},
        {"metric": "classification_unique_transcripts", "value": str(classification_metrics["classification_unique_transcripts"]), "note": "unique query_transcript_id values"},
        {"metric": "classification_duplicate_transcript_id_rows", "value": str(classification_metrics["classification_duplicate_transcript_id_rows"]), "note": "should be 0"},
        {"metric": "gtf_feature_rows", "value": str(gtf_metrics["gtf_feature_rows"]), "note": "non-comment GTF rows scanned"},
        {"metric": "gtf_transcript_rows", "value": str(gtf_metrics["gtf_transcript_rows"]), "note": "feature == transcript"},
        {"metric": "gtf_exon_rows", "value": str(gtf_metrics["gtf_exon_rows"]), "note": "feature == exon"},
        {"metric": "gtf_unexpected_feature_rows", "value": str(gtf_metrics["gtf_unexpected_feature_rows"]), "note": "rows whose feature is not transcript/exon"},
        {"metric": "unexpected_feature_count", "value": str(gtf_metrics["gtf_unexpected_feature_rows"]), "note": "alias of gtf_unexpected_feature_rows; should be 0"},
        {"metric": "gtf_unique_transcripts", "value": str(gtf_metrics["gtf_unique_transcripts"]), "note": "unique transcript_id values in annotated GTF"},
        {"metric": "gtf_missing_transcript_id_rows", "value": str(gtf_metrics["gtf_missing_transcript_id_rows"]), "note": "should be 0"},
        {"metric": "gtf_missing_gene_id_rows", "value": str(gtf_metrics["gtf_missing_gene_id_rows"]), "note": "should be 0"},
        {"metric": "gtf_transcripts_missing_classification", "value": str(merge_metrics["gtf_transcripts_missing_classification"]), "note": "GTF transcript IDs absent from classification table"},
        {"metric": "classification_transcripts_missing_gtf", "value": str(merge_metrics["classification_transcripts_missing_gtf"]), "note": "classification transcript IDs absent from annotated GTF"},
        {"metric": "class_code_mismatch_transcripts", "value": str(merge_metrics["class_code_mismatch_transcripts"]), "note": "class_code disagreement between GTF attributes and classification TSV"},
        {"metric": "exon_count_mismatch_transcripts", "value": str(merge_metrics["exon_count_mismatch_transcripts"]), "note": "num_exons disagreement between GTF and classification TSV"},
        {"metric": "length_mismatch_transcripts", "value": str(merge_metrics["length_mismatch_transcripts"]), "note": "length disagreement between exon-summed GTF length and classification TSV"},
        {"metric": "multi_seqname_transcripts", "value": str(merge_metrics["multi_seqname_transcripts"]), "note": "transcript_id observed on more than one seqname"},
        {"metric": "multi_strand_transcripts", "value": str(merge_metrics["multi_strand_transcripts"]), "note": "transcript_id observed on more than one strand"},
        {"metric": "transcripts_without_single_transcript_row", "value": str(merge_metrics["transcripts_without_single_transcript_row"]), "note": "should be 0"},
        {"metric": "transcripts_without_exon", "value": str(merge_metrics["transcripts_without_exon"]), "note": "should be 0"},
        {"metric": "basic_info_rows", "value": str(merge_metrics["basic_info_rows"]), "note": str(args.out_basic_info)},
        {"metric": "output_query_gtf_rows", "value": str(copy_metrics["output_query_gtf_rows"]), "note": str(args.out_query_gtf)},
        {"metric": "output_query_gtf_transcript_rows", "value": str(copy_metrics["output_query_gtf_transcript_rows"]), "note": "feature == transcript in output query GTF"},
        {"metric": "output_query_gtf_exon_rows", "value": str(copy_metrics["output_query_gtf_exon_rows"]), "note": "feature == exon in output query GTF"},
        {"metric": "output_query_gtf_unexpected_feature_rows", "value": str(copy_metrics["output_query_gtf_unexpected_feature_rows"]), "note": "output query GTF rows whose feature is not transcript/exon"},
        {"metric": "max_issue_examples", "value": str(args.max_issue_examples), "note": str(args.out_issue_examples)},
        {"metric": "runtime_sec", "value": str(elapsed_sec), "note": "wall-clock runtime of this script"},
        {"metric": "pass_status", "value": "PASS" if pass_status else "FAIL", "note": "PASS requires full GTF/classification coverage and no key mismatches"},
    ]


def build_manifest_rows(args: argparse.Namespace) -> List[Dict[str, str]]:
    rows = [
        manifest_row("annotated_gtf_input", args.annotated_gtf, "Step 1.6 full gffcompare annotated query GTF"),
        manifest_row("classification_tsv_input", args.classification_tsv, "Step 1.6 transcript classification table"),
        manifest_row("query_gtf_output", args.out_query_gtf, "full query GTF for downstream FEELnc steps"),
        manifest_row("query_basic_info_output", args.out_basic_info, "transcript-level GTF/classification merged metadata"),
        manifest_row("class_code_summary_output", args.out_class_summary, "class_code counts for the full query set"),
        manifest_row("issue_examples_output", args.out_issue_examples, "representative duplicate transcript and length-mismatch examples"),
        manifest_row("qc_summary_output", args.qc_summary, "QC metrics for query transcript preparation"),
    ]
    rows.append(
        {
            "item": "manifest_output",
            "path": str(args.manifest),
            "exists": "1",
            "size_bytes": "NA",
            "note": "this manifest file; size is not known before writing",
        }
    )
    return rows


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    args = parse_args(argv)
    start_time = datetime.now()

    logging.info("Step 1.6.5b FEELnc query transcript preparation start")
    require_readable(args.annotated_gtf, "annotated GTF")
    require_readable(args.classification_tsv, "classification TSV")
    prepare_outputs(
        [
            args.out_query_gtf,
            args.out_basic_info,
            args.out_class_summary,
            args.out_issue_examples,
            args.qc_summary,
            args.manifest,
        ],
        overwrite=args.overwrite,
    )

    classification, classification_metrics, classification_issue_examples = load_classification(
        args.classification_tsv,
        args.max_issue_examples,
    )
    info_by_tid, transcript_order, gtf_metrics = scan_annotated_gtf(args.annotated_gtf)
    basic_rows, merge_metrics, merge_issue_examples = build_basic_rows(
        transcript_order,
        info_by_tid,
        classification,
        args.max_issue_examples,
    )

    logging.info("stage 4/5 write transcript basic info, class_code summary and issue examples")
    write_tsv(args.out_basic_info, BASIC_INFO_FIELDS, basic_rows)
    write_tsv(args.out_class_summary, CLASS_SUMMARY_FIELDS, class_summary_rows(basic_rows))
    write_tsv(args.out_issue_examples, ISSUE_EXAMPLE_FIELDS, classification_issue_examples + merge_issue_examples)
    copy_metrics = copy_gtf(args.annotated_gtf, args.out_query_gtf)

    logging.info("stage 5/5 write QC summary and manifest")
    write_tsv(
        args.qc_summary,
        QC_FIELDS,
        build_qc_rows(args, classification_metrics, gtf_metrics, merge_metrics, copy_metrics, start_time),
    )
    write_tsv(args.manifest, MANIFEST_FIELDS, build_manifest_rows(args))

    logging.info("Step 1.6.5b FEELnc query transcript preparation completed")
    logging.info("query GTF: %s", args.out_query_gtf)
    logging.info("basic info: %s", args.out_basic_info)
    logging.info("QC summary: %s", args.qc_summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        logging.exception("Step 1.6.5b FEELnc query transcript preparation failed: %s", exc)
        raise SystemExit(1)
