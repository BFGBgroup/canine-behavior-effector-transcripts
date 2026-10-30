#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


CLASS_CODE_CATEGORY = {
    "=": "known_transcript_match",
    "c": "contained_in_reference",
    "j": "known_gene_novel_isoform",
    "e": "single_exon_overlap_reference",
    "i": "intronic_novel_transcript",
    "o": "generic_exonic_overlap",
    "p": "possible_polymerase_run_on_or_readthrough",
    "r": "repeat_overlap",
    "s": "opposite_strand_intron_match",
    "u": "intergenic_novel_transcript",
    "x": "antisense_novel_transcript",
    "y": "contains_reference_within_introns",
    "k": "contains_reference_transcript",
    "m": "retained_intron_or_fragment",
    "n": "retained_intron_or_fragment",
    ".": "unclassified",
}

CLASS_CODE_RELEVANCE = {
    "u": "candidate_novel_gene_or_lncRNA_priority",
    "x": "candidate_antisense_lncRNA_priority",
    "i": "candidate_intronic_transcript_caution",
    "j": "known_gene_novel_isoform_not_new_gene",
    "p": "readthrough_or_runon_caution",
    "=": "known_transcript_low_priority",
    "c": "known_or_contained_low_priority",
    "e": "known_or_exonic_overlap_low_priority",
    "o": "exonic_overlap_caution",
    "r": "repeat_overlap_artifact_risk",
}

CLASS_CODE_PRIORITY = {
    "u": "high",
    "x": "high",
    "i": "medium",
    "j": "medium",
    "p": "low",
    "=": "low",
    "c": "low",
    "e": "low",
    "o": "low",
    "r": "excluded_candidate",
}

TMAP_REQUIRED_FIELDS = ["class_code", "qry_id", "qry_gene_id"]

TRANSCRIPT_CLASSIFICATION_FIELDS = [
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
    "tmap_row_number",
]

CLASS_SUMMARY_FIELDS = [
    "class_code",
    "transcript_count",
    "novelty_category",
    "candidate_relevance",
    "initial_priority_for_anchor_overlap",
]

MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]
QC_FIELDS = ["metric", "value", "note"]


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
        description="Summarize gffcompare class-code output for five-dogs StringTie3 discovery."
    )
    parser.add_argument("--tmap", required=True, type=Path)
    parser.add_argument("--annotated-gtf", required=True, type=Path)
    parser.add_argument("--tracking", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--loci", required=True, type=Path)
    parser.add_argument("--refmap", required=True, type=str)
    parser.add_argument("--merged-gtf", required=True, type=Path)
    parser.add_argument("--reference-gtf", required=True, type=Path)
    parser.add_argument("--step15-qc", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--gffcompare-prefix", required=True)
    parser.add_argument("--gffcompare-command", required=True)
    return parser.parse_args(argv)


def require_readable(path: Path, label: str) -> None:
    if str(path) == "NA":
        fail(f"{label} path is NA")
    if not path.exists():
        fail(f"{label} does not exist: {path}")
    if not path.is_file():
        fail(f"{label} is not a regular file: {path}")
    if path.stat().st_size == 0:
        fail(f"{label} is empty: {path}")


def maybe_path(value: str) -> Path | None:
    if clean(value) == "NA":
        return None
    return Path(value)


def require_fields(path: Path, fieldnames: Sequence[str], required: Sequence[str]) -> None:
    missing = [field for field in required if field not in fieldnames]
    if missing:
        fail(f"{path} missing required fields: {','.join(missing)}")


def class_category(code: str) -> str:
    return CLASS_CODE_CATEGORY.get(code, "other_or_unclassified")


def class_relevance(code: str) -> str:
    return CLASS_CODE_RELEVANCE.get(code, "other_caution")


def class_priority(code: str) -> str:
    return CLASS_CODE_PRIORITY.get(code, "medium")


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean(row.get(field)) for field in fieldnames})


def parse_tmap(tmap: Path) -> Tuple[List[Dict[str, str]], Counter]:
    logging.info("read gffcompare tmap: %s", tmap)
    rows: List[Dict[str, str]] = []
    code_counter: Counter = Counter()
    with tmap.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(tmap, reader.fieldnames or [], TMAP_REQUIRED_FIELDS)
        for row_number, row in enumerate(reader, start=1):
            code = clean(row.get("class_code"))
            code_counter[code] += 1
            rows.append(
                {
                    "query_transcript_id": clean(row.get("qry_id")),
                    "query_gene_id": clean(row.get("qry_gene_id")),
                    "reference_gene_id": clean(row.get("ref_gene_id")),
                    "reference_transcript_id": clean(row.get("ref_id")),
                    "class_code": code,
                    "novelty_category": class_category(code),
                    "candidate_relevance": class_relevance(code),
                    "initial_priority_for_anchor_overlap": class_priority(code),
                    "num_exons": clean(row.get("num_exons")),
                    "FPKM": clean(row.get("FPKM")),
                    "TPM": clean(row.get("TPM")),
                    "cov": clean(row.get("cov")),
                    "length": clean(row.get("len")),
                    "major_iso_id": clean(row.get("major_iso_id")),
                    "ref_match_len": clean(row.get("ref_match_len")),
                    "tmap_row_number": str(row_number),
                }
            )
            if row_number % 50000 == 0:
                logging.info("read tmap rows: %s", row_number)
    logging.info("tmap read completed: %s rows", len(rows))
    return rows, code_counter


def parse_tracking_reference(value: str) -> Tuple[str, str]:
    text = clean(value)
    if text in {"NA", "-"}:
        return "NA", "NA"
    parts = text.split("|", 1)
    ref_gene_id = clean(parts[0])
    ref_id = clean(parts[1] if len(parts) > 1 else "NA")
    return ref_gene_id, ref_id


def parse_tracking_query(value: str, tracking: Path, row_number: int) -> Dict[str, str]:
    text = clean(value)
    if ":" not in text:
        fail(f"{tracking} row {row_number} query field is missing the qN: prefix: {text}")
    _, payload = text.split(":", 1)
    parts = payload.split("|")
    if len(parts) < 7:
        fail(f"{tracking} row {row_number} query field has too few columns: {text}")
    return {
        "query_gene_id": clean(parts[0]),
        "query_transcript_id": clean(parts[1]),
        "num_exons": clean(parts[2]),
        "FPKM": clean(parts[3]),
        "TPM": clean(parts[4]),
        "cov": clean(parts[5]),
        "length": clean(parts[6]),
    }


def parse_tracking(tracking: Path) -> Tuple[List[Dict[str, str]], Counter]:
    logging.info("not provided tmap, read instead from gffcompare tracking read class-code: %s", tracking)
    rows: List[Dict[str, str]] = []
    code_counter: Counter = Counter()
    with tracking.open("r", errors="replace") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                fail(f"{tracking}: row {row_number} has too few columns to parse the tracking record")
            ref_gene_id, ref_id = parse_tracking_reference(parts[2])
            code = clean(parts[3])
            query = parse_tracking_query(parts[4], tracking, row_number)
            code_counter[code] += 1
            rows.append(
                {
                    "query_transcript_id": query["query_transcript_id"],
                    "query_gene_id": query["query_gene_id"],
                    "reference_gene_id": ref_gene_id,
                    "reference_transcript_id": ref_id,
                    "class_code": code,
                    "novelty_category": class_category(code),
                    "candidate_relevance": class_relevance(code),
                    "initial_priority_for_anchor_overlap": class_priority(code),
                    "num_exons": query["num_exons"],
                    "FPKM": query["FPKM"],
                    "TPM": query["TPM"],
                    "cov": query["cov"],
                    "length": query["length"],
                    "major_iso_id": "NA",
                    "ref_match_len": "NA",
                    "tmap_row_number": str(row_number),
                }
            )
            if row_number % 50000 == 0:
                logging.info("read tracking rows: %s", row_number)
    logging.info("tracking read completed: %s rows", len(rows))
    return rows, code_counter


def count_gtf_features(path: Path) -> Counter:
    counts: Counter = Counter()
    if not path.exists() or path.stat().st_size == 0:
        return counts
    with path.open("r", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            feature = parts[2]
            counts[feature] += 1
            if feature == "transcript":
                seqname = parts[0]
                if seqname.startswith("chr"):
                    counts["transcript_on_chr_like_contig"] += 1
                else:
                    counts["transcript_on_non_chr_contig"] += 1
    return counts


def read_step15_qc(path: Path) -> Dict[str, str]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "metric" not in reader.fieldnames or "value" not in reader.fieldnames:
            return {}
        return {clean(row.get("metric")): clean(row.get("value")) for row in reader}


def file_manifest_row(item: str, path: Path | None, note: str) -> Dict[str, str]:
    if path is None:
        return {
            "item": item,
            "path": "NA",
            "exists": "0",
            "size_bytes": "NA",
            "note": note,
        }
    return {
        "item": item,
        "path": str(path),
        "exists": "1" if path.exists() else "0",
        "size_bytes": str(path.stat().st_size) if path.exists() and path.is_file() else "NA",
        "note": note,
    }


def build_class_summary(code_counter: Counter) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for code, count in code_counter.most_common():
        rows.append(
            {
                "class_code": code,
                "transcript_count": str(count),
                "novelty_category": class_category(code),
                "candidate_relevance": class_relevance(code),
                "initial_priority_for_anchor_overlap": class_priority(code),
            }
        )
    return rows


def build_qc_rows(
    args: argparse.Namespace,
    rows: Sequence[Mapping[str, str]],
    code_counter: Counter,
    merged_counts: Counter,
    reference_counts: Counter,
    annotated_counts: Counter,
    step15_qc: Mapping[str, str],
    tmap_path: Path | None,
    refmap_path: Path | None,
    classification_source: str,
) -> List[Dict[str, str]]:
    required_outputs = [args.tracking, args.stats, args.loci]
    if tmap_path is not None:
        required_outputs.append(tmap_path)
    if clean(str(args.annotated_gtf)) != "NA":
        required_outputs.append(args.annotated_gtf)
    required_ok = all(path.exists() and path.stat().st_size > 0 for path in required_outputs)
    pass_status = "PASS" if required_ok and len(rows) > 0 and merged_counts.get("transcript", 0) > 0 else "FAIL"
    return [
        {"metric": "step15_pass_status", "value": step15_qc.get("pass_status", "NA"), "note": "from Step 1.5 QC"},
        {"metric": "gffcompare_classification_rows", "value": str(len(rows)), "note": "query transcripts parsed from tmap or tracking"},
        {"metric": "gffcompare_classification_source", "value": classification_source, "note": "tmap when available, otherwise tracking"},
        {"metric": "merged_gtf_transcript_count", "value": str(merged_counts.get("transcript", 0)), "note": str(args.merged_gtf)},
        {"metric": "merged_gtf_exon_count", "value": str(merged_counts.get("exon", 0)), "note": str(args.merged_gtf)},
        {
            "metric": "merged_gtf_transcripts_on_chr_like_contig",
            "value": str(merged_counts.get("transcript_on_chr_like_contig", 0)),
            "note": "seqname starts with chr",
        },
        {
            "metric": "merged_gtf_transcripts_on_non_chr_contig",
            "value": str(merged_counts.get("transcript_on_non_chr_contig", 0)),
            "note": "seqname does not start with chr",
        },
        {"metric": "reference_gtf_transcript_count", "value": str(reference_counts.get("transcript", 0)), "note": str(args.reference_gtf)},
        {"metric": "annotated_gtf_transcript_count", "value": str(annotated_counts.get("transcript", 0)), "note": str(args.annotated_gtf)},
        {"metric": "class_code_equal_count", "value": str(code_counter.get("=", 0)), "note": "known transcript match"},
        {"metric": "class_code_j_count", "value": str(code_counter.get("j", 0)), "note": "known gene novel isoform"},
        {"metric": "class_code_u_count", "value": str(code_counter.get("u", 0)), "note": "intergenic novel transcript priority"},
        {"metric": "class_code_x_count", "value": str(code_counter.get("x", 0)), "note": "antisense novel transcript priority"},
        {"metric": "class_code_i_count", "value": str(code_counter.get("i", 0)), "note": "intronic novel transcript caution"},
        {"metric": "class_code_p_count", "value": str(code_counter.get("p", 0)), "note": "possible readthrough/run-on caution"},
        {
            "metric": "refmap_exists",
            "value": "1" if refmap_path is not None and refmap_path.exists() else "0",
            "note": str(refmap_path) if refmap_path is not None else "NA",
        },
        {"metric": "pass_status", "value": pass_status, "note": "PASS requires non-empty gffcompare outputs and non-empty merged GTF"},
    ]


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    setup_logging()

    logging.info("Step 1.6 gffcompare summary start")
    args.outdir.mkdir(parents=True, exist_ok=True)

    tmap_path = maybe_path(str(args.tmap))
    if tmap_path is not None:
        require_readable(tmap_path, "gffcompare tmap")
    require_readable(args.tracking, "gffcompare tracking")
    require_readable(args.stats, "gffcompare stats")
    require_readable(args.loci, "gffcompare loci")
    require_readable(args.merged_gtf, "merged StringTie3 GTF")
    require_readable(args.reference_gtf, "reference dog GTF")
    require_readable(args.step15_qc, "Step 1.5 QC")
    if clean(str(args.annotated_gtf)) != "NA":
        require_readable(args.annotated_gtf, "gffcompare annotated/combined GTF")

    refmap_path = maybe_path(args.refmap)
    if refmap_path is not None and (not refmap_path.exists() or refmap_path.stat().st_size == 0):
        logging.warning("refmap does not exist or is empty, will be manifest/QC record in: %s", refmap_path)

    if tmap_path is not None:
        classification_source = "tmap"
        classification_rows, code_counter = parse_tmap(tmap_path)
    else:
        classification_source = "tracking"
        classification_rows, code_counter = parse_tracking(args.tracking)
    class_summary_rows = build_class_summary(code_counter)

    logging.info("statistics merged/reference/annotated GTF feature counts")
    merged_counts = count_gtf_features(args.merged_gtf)
    reference_counts = count_gtf_features(args.reference_gtf)
    annotated_counts = count_gtf_features(args.annotated_gtf) if clean(str(args.annotated_gtf)) != "NA" else Counter()
    step15_qc = read_step15_qc(args.step15_qc)

    classification_path = args.outdir / "1.6_gffcompare_transcript_classification.tsv"
    class_summary_path = args.outdir / "1.6_gffcompare_class_code_summary.tsv"
    qc_path = args.outdir / "1.6_gffcompare_qc_summary.tsv"
    manifest_path = args.outdir / "1.6_gffcompare_manifest.tsv"

    logging.info("write transcript classification: %s", classification_path)
    write_tsv(classification_path, TRANSCRIPT_CLASSIFICATION_FIELDS, classification_rows)

    logging.info("write class-code summary: %s", class_summary_path)
    write_tsv(class_summary_path, CLASS_SUMMARY_FIELDS, class_summary_rows)

    qc_rows = build_qc_rows(
        args=args,
        rows=classification_rows,
        code_counter=code_counter,
        merged_counts=merged_counts,
        reference_counts=reference_counts,
        annotated_counts=annotated_counts,
        step15_qc=step15_qc,
        tmap_path=tmap_path,
        refmap_path=refmap_path,
        classification_source=classification_source,
    )

    logging.info("write QC summary: %s", qc_path)
    write_tsv(qc_path, QC_FIELDS, qc_rows)

    manifest_rows = [
        file_manifest_row("merged_gtf_input", args.merged_gtf, "Step 1.5 merged StringTie3 transcriptome"),
        file_manifest_row("reference_gtf_input", args.reference_gtf, "CanFam3.1 dog reference annotation"),
        file_manifest_row("step15_qc_input", args.step15_qc, "Step 1.5 QC summary"),
        file_manifest_row("gffcompare_tmap", tmap_path, "per-query transcript class-code mapping when produced by gffcompare"),
        file_manifest_row("gffcompare_refmap", refmap_path, "reference-to-query mapping when produced by gffcompare"),
        file_manifest_row("gffcompare_tracking", args.tracking, "gffcompare tracking output"),
        file_manifest_row("gffcompare_stats", args.stats, "gffcompare stats output"),
        file_manifest_row("gffcompare_loci", args.loci, "gffcompare loci output"),
        file_manifest_row("gffcompare_annotated_or_combined_gtf", args.annotated_gtf, "annotated/combined GTF output"),
        file_manifest_row("transcript_classification_tsv", classification_path, "parsed transcript-level classification"),
        file_manifest_row("class_code_summary_tsv", class_summary_path, "class-code counts"),
        file_manifest_row("qc_summary_tsv", qc_path, "Step 1.6 QC summary"),
        {
            "item": "manifest_tsv",
            "path": str(manifest_path),
            "exists": "1",
            "size_bytes": "NA",
            "note": "Step 1.6 manifest; size is finalized after this row is written",
        },
        {
            "item": "gffcompare_prefix",
            "path": args.gffcompare_prefix,
            "exists": "NA",
            "size_bytes": "NA",
            "note": "output prefix passed to gffcompare",
        },
        {
            "item": "gffcompare_command",
            "path": args.gffcompare_command,
            "exists": "NA",
            "size_bytes": "NA",
            "note": "command recorded by batch wrapper",
        },
    ]

    logging.info("write manifest: %s", manifest_path)
    write_tsv(manifest_path, MANIFEST_FIELDS, manifest_rows)

    pass_row = next((row for row in qc_rows if row["metric"] == "pass_status"), None)
    if pass_row and pass_row["value"] != "PASS":
        fail(f"Step 1.6 QC failed: {pass_row}")

    logging.info("Step 1.6 gffcompare summary completed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except Exception as exc:
        logging.exception("Step 1.6 gffcompare summary failed: %s", exc)
        raise SystemExit(1)
