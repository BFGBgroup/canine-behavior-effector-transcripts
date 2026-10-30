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
from typing import Dict, List, Mapping, Sequence, Set, Tuple


STATUS_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "class_code",
    "transcript_length",
    "exon_count",
    "single_exon_flag",
    "FEELnc_filter_status",
    "FEELnc_CPS",
    "FEELnc_RF_label",
    "FEELnc_codpot_category",
    "FEELnc_noORF_flag",
    "FEELnc_biotype",
    "FEELnc_classifier_best_isBest",
    "FEELnc_lncRNA_class",
    "FEELnc_lncRNA_subtype",
    "FEELnc_lncRNA_location",
    "FEELnc_partner_gene",
    "FEELnc_partner_transcript",
    "FEELnc_partner_direction",
    "FEELnc_partner_distance",
    "FEELnc_consistency_flags",
]
SUMMARY_FIELDS = ["metric", "value", "note"]
MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]
REQUIRED_BASIC_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "class_code",
    "transcript_length",
    "exon_count",
    "single_exon_flag",
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
        description="Collect FEELnc module outputs for all-query transcript annotation."
    )
    parser.add_argument("--query-basic-info", required=True, type=Path)
    parser.add_argument("--query-gtf", required=True, type=Path)
    parser.add_argument("--filter-gtf", required=True, type=Path)
    parser.add_argument("--filter-log", required=True, type=Path)
    parser.add_argument("--codpot-rf", required=True, type=Path)
    parser.add_argument("--codpot-rf-summary", required=True, type=Path)
    parser.add_argument("--codpot-lncrna-gtf", required=True, type=Path)
    parser.add_argument("--codpot-mrna-gtf", required=True, type=Path)
    parser.add_argument("--codpot-tucp-gtf", required=True, type=Path)
    parser.add_argument("--codpot-noorf-gtf", required=True, type=Path)
    parser.add_argument("--classifier-input-gtf", required=True, type=Path)
    parser.add_argument("--classifier-classes", required=True, type=Path)
    parser.add_argument("--classifier-log", required=True, type=Path)
    parser.add_argument("--status-tsv", required=True, type=Path)
    parser.add_argument("--summary-tsv", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
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


def iter_gtf_transcript_ids(path: Path) -> Tuple[Set[str], Counter]:
    metrics: Counter = Counter()
    transcript_ids: Set[str] = set()
    if not path.exists() or path.stat().st_size == 0:
        metrics["missing_or_empty_gtf"] = 1
        return transcript_ids, metrics

    with path.open("r", errors="replace") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                fail(f"{path} line {row_number} is not a 9-column GTF row")
            metrics["gtf_rows"] += 1
            metrics[f"gtf_{fields[2]}_rows"] += 1
            attrs = parse_gtf_attributes(fields[8])
            tid = clean(attrs.get("transcript_id"))
            if tid == "NA":
                metrics["missing_transcript_id_rows"] += 1
                continue
            transcript_ids.add(tid)
    metrics["unique_transcripts"] = len(transcript_ids)
    return transcript_ids, metrics


def read_basic_info(path: Path) -> List[Dict[str, str]]:
    logging.info("stage 1/5 read 1.6.5b transcript basic info: %s", path)
    rows: List[Dict[str, str]] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"query basic info has no header: {path}")
        missing = [field for field in REQUIRED_BASIC_FIELDS if field not in reader.fieldnames]
        if missing:
            fail("query basic info missing required fields: " + ",".join(missing))
        for row_number, row in enumerate(reader, start=2):
            tid = clean(row.get("query_transcript_id"))
            if tid == "NA":
                fail(f"query basic info line {row_number} has empty query_transcript_id")
            rows.append({key: clean(value) for key, value in row.items()})
            if row_number % 25000 == 0:
                logging.info("stage 1/5 read basic info row count: %s", row_number - 1)
    logging.info("stage 1/5 completed, basic info transcript count: %s", len(rows))
    return rows


def read_rf_scores(path: Path) -> Tuple[Dict[str, Dict[str, str]], Counter]:
    logging.info("stage 2/5 read FEELnc_codpot RF score: %s", path)
    scores: Dict[str, Dict[str, str]] = {}
    metrics: Counter = Counter()
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"codpot RF table has no header: {path}")
        for required in ["name", "coding_potential", "label"]:
            if required not in reader.fieldnames:
                fail(f"codpot RF table missing required field {required}: {path}")
        for row_number, row in enumerate(reader, start=2):
            tid = clean(row.get("name"))
            if tid == "NA":
                metrics["rf_missing_name_rows"] += 1
                continue
            if tid in scores:
                metrics["rf_duplicate_name_rows"] += 1
                continue
            scores[tid] = {key: clean(value) for key, value in row.items()}
            metrics["rf_rows"] += 1
            if row_number % 25000 == 0:
                logging.info("stage 2/5 read RF row count: %s", row_number - 1)
    metrics["rf_unique_transcripts"] = len(scores)
    return scores, metrics


def read_classifier_best(path: Path) -> Tuple[Dict[str, Dict[str, str]], Counter]:
    logging.info("stage 3/5 read FEELnc_classifier classes: %s", path)
    best_by_tid: Dict[str, Dict[str, str]] = {}
    metrics: Counter = Counter()
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"classifier classes table has no header: {path}")
        for required in ["isBest", "lncRNA_transcript", "partnerRNA_gene", "partnerRNA_transcript", "direction", "type", "distance", "subtype", "location"]:
            if required not in reader.fieldnames:
                fail(f"classifier classes table missing required field {required}: {path}")
        for row in reader:
            metrics["classifier_interaction_rows"] += 1
            if clean(row.get("isBest")) != "1":
                continue
            tid = clean(row.get("lncRNA_transcript"))
            if tid == "NA":
                metrics["classifier_best_missing_transcript_rows"] += 1
                continue
            if tid in best_by_tid:
                metrics["classifier_duplicate_best_transcripts"] += 1
                continue
            best_by_tid[tid] = {key: clean(value) for key, value in row.items()}
    metrics["classifier_best_transcripts"] = len(best_by_tid)
    return best_by_tid, metrics


def assign_codpot_categories(category_sets: Mapping[str, Set[str]]) -> Tuple[Dict[str, str], Counter]:
    category_by_tid: Dict[str, str] = {}
    metrics: Counter = Counter()
    for category, ids in category_sets.items():
        metrics[f"codpot_{category}_transcripts"] = len(ids)
        for tid in ids:
            if tid in category_by_tid and category_by_tid[tid] != category:
                metrics["codpot_category_overlap_transcripts"] += 1
                category_by_tid[tid] = category_by_tid[tid] + ";" + category
            else:
                category_by_tid[tid] = category
    return category_by_tid, metrics


def biotype_from_category(category: str) -> str:
    if category in {"lncRNA", "noORF"}:
        return "lncRNA"
    if category == "mRNA":
        return "mRNA"
    if category == "TUCp":
        return "TUCp"
    return "unknown"


def build_status_rows(
    basic_rows: Sequence[Mapping[str, str]],
    filter_pass_ids: Set[str],
    rf_scores: Mapping[str, Mapping[str, str]],
    codpot_category_by_tid: Mapping[str, str],
    classifier_best: Mapping[str, Mapping[str, str]],
) -> Tuple[List[Dict[str, str]], Counter]:
    logging.info("stage 4/5 merge FEELnc module status to transcript-level table")
    rows: List[Dict[str, str]] = []
    metrics: Counter = Counter()
    seen_tids: Set[str] = set()

    for basic in basic_rows:
        tid = clean(basic.get("query_transcript_id"))
        if tid in seen_tids:
            metrics["basic_duplicate_transcript_rows"] += 1
        seen_tids.add(tid)

        category = clean(codpot_category_by_tid.get(tid))
        rf_row = rf_scores.get(tid, {})
        best = classifier_best.get(tid, {})
        flags: List[str] = []

        filter_status = "pass" if tid in filter_pass_ids else "filtered_or_not_lncRNA_like"
        if filter_status != "pass":
            flags.append("FEELnc_filter_removed")
            metrics["filter_removed_or_not_lncRNA_like"] += 1
        else:
            metrics["filter_pass"] += 1

        if category == "NA":
            category = "missing_codpot_category"
            flags.append("missing_codpot_category")
            metrics["missing_codpot_category"] += 1

        cps = clean(rf_row.get("coding_potential"))
        label = clean(rf_row.get("label"))
        if cps == "NA":
            flags.append("missing_RF_score")
            metrics["missing_RF_score_transcripts"] += 1

        biotype = biotype_from_category(category)
        if biotype == "lncRNA" and not best:
            flags.append("no_classifier_best_partner")
            metrics["lncRNA_without_classifier_best_partner"] += 1

        rows.append(
            {
                "query_transcript_id": tid,
                "query_gene_id": clean(basic.get("query_gene_id")),
                "class_code": clean(basic.get("class_code")),
                "transcript_length": clean(basic.get("transcript_length")),
                "exon_count": clean(basic.get("exon_count")),
                "single_exon_flag": clean(basic.get("single_exon_flag")),
                "FEELnc_filter_status": filter_status,
                "FEELnc_CPS": cps,
                "FEELnc_RF_label": label,
                "FEELnc_codpot_category": category,
                "FEELnc_noORF_flag": "1" if category == "noORF" else "0",
                "FEELnc_biotype": biotype,
                "FEELnc_classifier_best_isBest": clean(best.get("isBest")),
                "FEELnc_lncRNA_class": clean(best.get("type")),
                "FEELnc_lncRNA_subtype": clean(best.get("subtype")),
                "FEELnc_lncRNA_location": clean(best.get("location")),
                "FEELnc_partner_gene": clean(best.get("partnerRNA_gene")),
                "FEELnc_partner_transcript": clean(best.get("partnerRNA_transcript")),
                "FEELnc_partner_direction": clean(best.get("direction")),
                "FEELnc_partner_distance": clean(best.get("distance")),
                "FEELnc_consistency_flags": ";".join(flags) if flags else "ok",
            }
        )

    metrics["status_rows"] = len(rows)
    return rows, metrics


def parse_codpot_summary(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not path.exists() or path.stat().st_size == 0:
        return values
    with path.open("r", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#") or ":" not in text:
                continue
            key, value = text.split(":", 1)
            values[key.strip().lstrip("-")] = value.strip()
    return values


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summary_row(metric: str, value: object, note: str) -> Dict[str, str]:
    return {"metric": metric, "value": clean(value), "note": note}


def build_summary_rows(
    args: argparse.Namespace,
    basic_count: int,
    gtf_metrics: Mapping[str, Counter],
    rf_metrics: Counter,
    codpot_metrics: Counter,
    classifier_metrics: Counter,
    status_metrics: Counter,
    codpot_summary: Mapping[str, str],
    start_time: datetime,
) -> List[Dict[str, str]]:
    codpot_total_unique = (
        codpot_metrics["codpot_lncRNA_transcripts"]
        + codpot_metrics["codpot_mRNA_transcripts"]
        + codpot_metrics["codpot_TUCp_transcripts"]
        + codpot_metrics["codpot_noORF_transcripts"]
    )
    pass_status = (
        basic_count > 0
        and gtf_metrics["query"]["unique_transcripts"] == basic_count
        and gtf_metrics["filter"]["unique_transcripts"] > 0
        and rf_metrics["rf_unique_transcripts"] > 0
        and codpot_total_unique == basic_count
        and status_metrics["missing_codpot_category"] == 0
        and codpot_metrics["codpot_category_overlap_transcripts"] == 0
        and classifier_metrics["classifier_interaction_rows"] >= classifier_metrics["classifier_best_transcripts"]
        and status_metrics["status_rows"] == basic_count
        and status_metrics["basic_duplicate_transcript_rows"] == 0
    )
    elapsed_sec = int((datetime.now() - start_time).total_seconds())

    rows = [
        summary_row("query_basic_info_input", args.query_basic_info, "1.6.5b transcript metadata"),
        summary_row("query_gtf_input", args.query_gtf, "1.6.5b full query GTF used by codpot"),
        summary_row("filter_gtf_input", args.filter_gtf, "FEELnc_filter pass GTF used only as filter evidence"),
        summary_row("codpot_rf_input", args.codpot_rf, "FEELnc_codpot RF score table"),
        summary_row("classifier_classes_input", args.classifier_classes, "FEELnc_classifier interaction table"),
        summary_row("query_basic_info_transcripts", basic_count, "all transcript rows in 1.6.5b basic info"),
        summary_row("query_gtf_unique_transcripts", gtf_metrics["query"]["unique_transcripts"], "unique transcript_id in query GTF"),
        summary_row("filter_pass_transcripts", gtf_metrics["filter"]["unique_transcripts"], "unique transcript_id in FEELnc_filter output"),
        summary_row("filter_removed_or_not_lncRNA_like", status_metrics["filter_removed_or_not_lncRNA_like"], "query transcripts absent from filter pass GTF"),
        summary_row("codpot_RF_rows", rf_metrics["rf_rows"], "rows in FEELnc RF score table"),
        summary_row("codpot_RF_unique_transcripts", rf_metrics["rf_unique_transcripts"], "unique transcript names in FEELnc RF table"),
        summary_row("codpot_lncRNA_transcripts", codpot_metrics["codpot_lncRNA_transcripts"], "transcripts in codpot .lncRNA.gtf"),
        summary_row("codpot_mRNA_transcripts", codpot_metrics["codpot_mRNA_transcripts"], "transcripts in codpot .mRNA.gtf"),
        summary_row("codpot_TUCp_transcripts", codpot_metrics["codpot_TUCp_transcripts"], "transcripts in codpot .TUCp.gtf"),
        summary_row("codpot_noORF_transcripts", codpot_metrics["codpot_noORF_transcripts"], "transcripts in codpot .noORF.gtf"),
        summary_row("codpot_total_unique_transcripts", codpot_total_unique, "lncRNA + mRNA + TUCp + noORF unique transcript categories"),
        summary_row("codpot_category_overlap_transcripts", codpot_metrics["codpot_category_overlap_transcripts"], "should be 0"),
        summary_row("missing_codpot_category", status_metrics["missing_codpot_category"], "should be 0; every query transcript should be assigned by FEELnc_codpot"),
        summary_row("missing_RF_score_transcripts", status_metrics["missing_RF_score_transcripts"], "may include noORF depending on FEELnc output"),
        summary_row("classifier_interaction_rows", classifier_metrics["classifier_interaction_rows"], "all classifier interaction rows"),
        summary_row("classifier_best_transcripts", classifier_metrics["classifier_best_transcripts"], "transcripts with isBest == 1"),
        summary_row("classifier_duplicate_best_transcripts", classifier_metrics["classifier_duplicate_best_transcripts"], "should be 0"),
        summary_row("lncRNA_without_classifier_best_partner", status_metrics["lncRNA_without_classifier_best_partner"], "lncRNA/noORF transcripts with no best partner row"),
        summary_row("status_rows", status_metrics["status_rows"], str(args.status_tsv)),
        summary_row("runtime_sec", elapsed_sec, "wall-clock runtime of this collector"),
        summary_row("pass_status", "PASS" if pass_status else "FAIL", "PASS checks output consistency, not biological correctness"),
    ]
    for key, value in sorted(codpot_summary.items()):
        rows.append(summary_row("codpot_summary_" + key, value, str(args.codpot_rf_summary)))
    return rows


def manifest_row(item: str, path: Path, note: str) -> Dict[str, str]:
    exists = path.exists()
    return {
        "item": item,
        "path": str(path),
        "exists": "1" if exists else "0",
        "size_bytes": str(path.stat().st_size) if exists and path.is_file() else "NA",
        "note": note,
    }


def build_manifest_rows(args: argparse.Namespace) -> List[Dict[str, str]]:
    return [
        manifest_row("query_basic_info_input", args.query_basic_info, "1.6.5b transcript metadata"),
        manifest_row("query_gtf_input", args.query_gtf, "1.6.5b full query GTF"),
        manifest_row("filter_gtf_input", args.filter_gtf, "FEELnc_filter pass GTF"),
        manifest_row("filter_log_input", args.filter_log, "FEELnc_filter log"),
        manifest_row("codpot_rf_input", args.codpot_rf, "FEELnc_codpot RF table"),
        manifest_row("codpot_rf_summary_input", args.codpot_rf_summary, "FEELnc_codpot summary"),
        manifest_row("codpot_lncrna_gtf_input", args.codpot_lncrna_gtf, "FEELnc_codpot lncRNA GTF"),
        manifest_row("codpot_mrna_gtf_input", args.codpot_mrna_gtf, "FEELnc_codpot mRNA GTF"),
        manifest_row("codpot_tucp_gtf_input", args.codpot_tucp_gtf, "optional FEELnc_codpot TUCp GTF"),
        manifest_row("codpot_noorf_gtf_input", args.codpot_noorf_gtf, "optional FEELnc_codpot noORF GTF"),
        manifest_row("classifier_input_gtf", args.classifier_input_gtf, "lncRNA/noORF GTF passed to FEELnc_classifier"),
        manifest_row("classifier_classes_input", args.classifier_classes, "FEELnc_classifier classes table"),
        manifest_row("classifier_log_input", args.classifier_log, "FEELnc_classifier log"),
        manifest_row("status_output", args.status_tsv, "transcript-level FEELnc module status"),
        manifest_row("summary_output", args.summary_tsv, "collector QC summary"),
        {
            "item": "manifest_output",
            "path": str(args.manifest),
            "exists": "1",
            "size_bytes": "NA",
            "note": "this manifest file; size is not known before writing",
        },
    ]


def main(argv: Sequence[str] | None = None) -> int:
    configure_logging()
    args = parse_args(argv)
    start_time = datetime.now()

    logging.info("Step 1.6.5c FEELnc module output collection start")
    for path, label in [
        (args.query_basic_info, "query basic info"),
        (args.query_gtf, "query GTF"),
        (args.filter_gtf, "FEELnc filter GTF"),
        (args.filter_log, "FEELnc filter log"),
        (args.codpot_rf, "FEELnc codpot RF"),
        (args.codpot_rf_summary, "FEELnc codpot RF summary"),
        (args.codpot_lncrna_gtf, "FEELnc codpot lncRNA GTF"),
        (args.codpot_mrna_gtf, "FEELnc codpot mRNA GTF"),
        (args.classifier_input_gtf, "FEELnc classifier input GTF"),
        (args.classifier_classes, "FEELnc classifier classes"),
        (args.classifier_log, "FEELnc classifier log"),
    ]:
        require_readable(path, label)

    prepare_outputs([args.status_tsv, args.summary_tsv, args.manifest], overwrite=args.overwrite)

    basic_rows = read_basic_info(args.query_basic_info)
    _query_ids, query_gtf_metrics = iter_gtf_transcript_ids(args.query_gtf)
    filter_pass_ids, filter_gtf_metrics = iter_gtf_transcript_ids(args.filter_gtf)
    rf_scores, rf_metrics = read_rf_scores(args.codpot_rf)
    lnc_ids, lnc_metrics = iter_gtf_transcript_ids(args.codpot_lncrna_gtf)
    mrna_ids, mrna_metrics = iter_gtf_transcript_ids(args.codpot_mrna_gtf)
    tucp_ids, tucp_metrics = iter_gtf_transcript_ids(args.codpot_tucp_gtf)
    noorf_ids, noorf_metrics = iter_gtf_transcript_ids(args.codpot_noorf_gtf)
    classifier_best, classifier_metrics = read_classifier_best(args.classifier_classes)

    codpot_category_by_tid, codpot_metrics = assign_codpot_categories(
        {"lncRNA": lnc_ids, "mRNA": mrna_ids, "TUCp": tucp_ids, "noORF": noorf_ids}
    )
    status_rows, status_metrics = build_status_rows(
        basic_rows,
        filter_pass_ids,
        rf_scores,
        codpot_category_by_tid,
        classifier_best,
    )
    codpot_summary = parse_codpot_summary(args.codpot_rf_summary)

    logging.info("stage 5/5 write status, summary and manifest")
    write_tsv(args.status_tsv, STATUS_FIELDS, status_rows)
    write_tsv(
        args.summary_tsv,
        SUMMARY_FIELDS,
        build_summary_rows(
            args,
            len(basic_rows),
            {
                "query": query_gtf_metrics,
                "filter": filter_gtf_metrics,
                "lncRNA": lnc_metrics,
                "mRNA": mrna_metrics,
                "TUCp": tucp_metrics,
                "noORF": noorf_metrics,
            },
            rf_metrics,
            codpot_metrics,
            classifier_metrics,
            status_metrics,
            codpot_summary,
            start_time,
        ),
    )
    write_tsv(args.manifest, MANIFEST_FIELDS, build_manifest_rows(args))

    logging.info("Step 1.6.5c FEELnc module output collection completed")
    logging.info("status TSV: %s", args.status_tsv)
    logging.info("summary TSV: %s", args.summary_tsv)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        logging.exception("Step 1.6.5c FEELnc module output collection failed: %s", exc)
        raise SystemExit(1)
