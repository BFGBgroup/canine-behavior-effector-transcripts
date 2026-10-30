#!/usr/bin/python3.10


from __future__ import annotations

import csv
import gzip
import math
import os
import shutil
import sys
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, TextIO, Tuple


SOURCE_RUN_ID = os.environ.get("SOURCE_RUN_ID", "20260827_201055")
CATALOGUE_ROOT = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.5_no_gwas_build_expression_ase_regulatory_catalogue'
)
CATALOGUE_DIR = CATALOGUE_ROOT / SOURCE_RUN_ID
VARIANT_CATALOGUE = (
    CATALOGUE_DIR / "no_gwas_expression_ase_regulatory_variant_catalogue.tsv.gz"
)
VALID_BLOCK_CATALOGUE = (
    CATALOGUE_DIR / "no_gwas_expression_ase_valid_block_catalogue.tsv"
)
INVALID_BLOCK_QC = CATALOGUE_DIR / "no_gwas_expression_ase_invalid_block_qc.tsv"
CATALOGUE_SUMMARY = (
    CATALOGUE_DIR / "no_gwas_expression_ase_catalogue_summary.tsv"
)
SOURCE_COMPLETION_MARKER = CATALOGUE_DIR / "execution_completion.marker"

OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.6_no_gwas_export_expression_ase_catalogue_deliverables'
)
RUN_ID = os.environ.get("RUN_ID", datetime.now().strftime("%Y%m%d_%H%M%S"))
RUN_DIR = OUTPUT_ROOT / RUN_ID
COMPLETION_MARKER = RUN_DIR / "execution_completion.marker"

SCRIPT_DIR = Path(
    '/path/to/repository/scripts/10_molecular_finemapping/8.6_no_gwas_export_expression_ase_catalogue_deliverables'
)
LOG_DIR = SCRIPT_DIR / "log"
MASTER_LOG = LOG_DIR / (
    "8.6.export_no_gwas_expression_ase_catalogue_deliverables.{0}.log".format(
        RUN_ID
    )
)

CONCURRENT_TASKS = 1
THREADS_PER_TASK = 1
MASTER_LOG_HANDLE: Optional[TextIO] = None

ALL_VARIANTS_TSV = (
    RUN_DIR / "no_gwas_expression_ase_all_interpretable_variants.tsv.gz"
)
PRIORITY_VARIANTS_TSV = (
    RUN_DIR / "no_gwas_expression_ase_priority_variants_ppa_ge_0.1.tsv.gz"
)
HIGH_CONFIDENCE_TSV = (
    RUN_DIR / "no_gwas_expression_ase_high_confidence_variants_ppa_ge_0.5.tsv.gz"
)
STRONG_TSV = RUN_DIR / "no_gwas_expression_ase_strong_variants_ppa_ge_0.8.tsv.gz"
VERY_HIGH_TSV = (
    RUN_DIR / "no_gwas_expression_ase_very_high_variants_ppa_ge_0.95.tsv.gz"
)
TOP_VARIANTS_TSV = RUN_DIR / "no_gwas_expression_ase_top_variants.tsv.gz"
CREDIBLE_SET_TSV = (
    RUN_DIR / "no_gwas_expression_ase_credible_set_variants.tsv.gz"
)
PRIORITY_MARKER_RECURRENCE_TSV = (
    RUN_DIR
    / "no_gwas_expression_ase_priority_marker_recurrence_ppa_ge_0.1.tsv.gz"
)
EXPORT_SUMMARY_TSV = RUN_DIR / "no_gwas_expression_ase_export_summary.tsv"

ALL_VARIANTS_BED = (
    RUN_DIR / "no_gwas_expression_ase_all_interpretable_variants.bed.gz"
)
PRIORITY_BED = (
    RUN_DIR / "no_gwas_expression_ase_priority_variants_ppa_ge_0.1.bed.gz"
)
HIGH_CONFIDENCE_BED = (
    RUN_DIR / "no_gwas_expression_ase_high_confidence_variants_ppa_ge_0.5.bed.gz"
)
STRONG_BED = RUN_DIR / "no_gwas_expression_ase_strong_variants_ppa_ge_0.8.bed.gz"
VERY_HIGH_BED = (
    RUN_DIR / "no_gwas_expression_ase_very_high_variants_ppa_ge_0.95.bed.gz"
)
TOP_VARIANTS_BED = RUN_DIR / "no_gwas_expression_ase_top_variants.bed.gz"
CREDIBLE_SET_BED = (
    RUN_DIR / "no_gwas_expression_ase_credible_set_variants.bed.gz"
)

PRIORITY_THRESHOLD = 0.1
HIGH_CONFIDENCE_THRESHOLD = 0.5
STRONG_THRESHOLD = 0.8
VERY_HIGH_THRESHOLD = 0.95
PROGRESS_EVERY = 1_000_000
ALLOWED_PPA_STATUS = {"finite", "finite_boundary_adjusted"}
ALLOWED_FEATURE_ORIGINS = {
    "reference_annotation_gene",
    "dog5_manual_new_gene",
}
EXPECTED_EXPRESSION_MATRIX_SCOPE = (
    "reference_annotation_plus_dog5_manual_new_gene"
)
EXPECTED_FEATURE_ORIGIN_RULE = (
    "gene-=>reference_annotation_gene;DOG5_NEWGENE_=>dog5_manual_new_gene"
)
EXPECTED_CATALOGUE_ROLE_BY_FLAGS = {
    ("0", "0"): "finite_background",
    ("1", "0"): "top_variant",
    ("0", "1"): "credible_set",
    ("1", "1"): "top_and_credible_set",
}

VARIANT_REQUIRED = [
    "block_id",
    "tissue",
    "feature_id",
    "feature_origin",
    "expression_matrix_scope",
    "marker_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "ppa",
    "ppa_tier",
    "catalogue_role",
    "is_top_variant",
    "is_credible_set_member",
    "result_status",
    "ppa_status",
    "result_interpretable",
]
VALID_BLOCK_REQUIRED = ["block_id", "result_status", "result_interpretable"]
INVALID_BLOCK_REQUIRED = ["block_id", "result_status", "result_interpretable"]
SUMMARY_REQUIRED = [
    "run_id",
    "source_run_id",
    "expression_matrix_scope",
    "feature_origin_rule",
    "interpretable_variant_rows",
    "top_variant_rows",
    "credible_set_variant_rows",
    "valid_blocks",
    "invalid_blocks_qc_only",
    "ppa_ge_0.95_variant_rows",
    "ppa_ge_0.8_variant_rows",
    "ppa_ge_0.5_variant_rows",
    "ppa_ge_0.1_variant_rows",
    "association_fdr_screening",
    "null_causal_model",
    "interpretation_scope",
    "not_interpretable_as",
]


def log(message: str) -> None:
    print(
        "[{0} ] {1}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message),
        flush=True,
    )


def initialize_execution_logging() -> None:
    global MASTER_LOG_HANDLE
    validate_run_id(RUN_ID, "RUN_ID")
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MASTER_LOG_HANDLE = MASTER_LOG.open("x", buffering=1)
    os.dup2(MASTER_LOG_HANDLE.fileno(), sys.stdout.fileno())
    os.dup2(MASTER_LOG_HANDLE.fileno(), sys.stderr.fileno())


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def safe_float(value: object) -> float:
    try:
        return float(clean(value))
    except (TypeError, ValueError):
        return float("nan")


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(clean(value))
    except (TypeError, ValueError):
        return default


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return "{0:.12g}".format(value)


def validate_run_id(value: str, label: str) -> None:
    if not value or Path(value).name != value or value in {".", ".."}:
        raise ValueError("{0} must be one directory name: {1}".format(label, value))


def required_columns(
    path: Path, observed: Sequence[str], required: Sequence[str]
) -> None:
    missing = sorted(set(required) - set(observed))
    if missing:
        raise ValueError(
            "{0} missing required columns: {1}".format(path, ",".join(missing))
        )


def read_key_value_tsv(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns(path, list(reader.fieldnames or []), ["key", "value"])
        for row in reader:
            key = clean(row.get("key"))
            if not key:
                raise ValueError("{0} contains an empty key".format(path))
            if key in values:
                raise ValueError("{0} contains duplicate key: {1}".format(path, key))
            values[key] = clean(row.get("value"))
    return values


def read_completion_marker(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            key, separator, value = text.partition("=")
            if not separator or not key:
                raise ValueError("invalid completion marker line: {0}".format(text))
            if key in values:
                raise ValueError("duplicate completion marker key: {0}".format(key))
            values[key] = value
    return values


def inspect_block_table(
    path: Path, required: Sequence[str], expected_status: str
) -> int:
    seen = set()
    count = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns(path, list(reader.fieldnames or []), required)
        for row in reader:
            block_id = clean(row.get("block_id"))
            if not block_id:
                raise ValueError("{0} contains an empty block_id".format(path))
            if block_id in seen:
                raise ValueError(
                    "{0} contains duplicate block_id: {1}".format(path, block_id)
                )
            seen.add(block_id)
            status = clean(row.get("result_status"))
            interpretable = clean(row.get("result_interpretable"))
            if expected_status == "valid":
                if status != "valid" or interpretable != "1":
                    raise ValueError(
                        "{0} contains a non-interpretable valid block: {1}".format(
                            path, block_id
                        )
                    )
            elif status == "valid" or interpretable == "1":
                raise ValueError(
                    "{0} contains an interpretable invalid block: {1}".format(
                        path, block_id
                    )
                )
            count += 1
    return count


def read_variant_header(path: Path) -> List[str]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
    required_columns(path, fields, VARIANT_REQUIRED)
    return fields


def validate_inputs() -> Tuple[List[str], Dict[str, str], int, int]:
    validate_run_id(SOURCE_RUN_ID, "SOURCE_RUN_ID")
    for path in [
        VARIANT_CATALOGUE,
        VALID_BLOCK_CATALOGUE,
        INVALID_BLOCK_QC,
        CATALOGUE_SUMMARY,
        SOURCE_COMPLETION_MARKER,
    ]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("missing or empty 8.5 input: {0}".format(path))
    if RUN_DIR.exists():
        raise FileExistsError(
            "refusing to reuse existing 8.6 run directory: {0}".format(RUN_DIR)
        )

    variant_fields = read_variant_header(VARIANT_CATALOGUE)
    summary = read_key_value_tsv(CATALOGUE_SUMMARY)
    source_completion = read_completion_marker(SOURCE_COMPLETION_MARKER)
    missing_summary = sorted(set(SUMMARY_REQUIRED) - set(summary))
    if missing_summary:
        raise ValueError(
            "{0} missing required keys: {1}".format(
                CATALOGUE_SUMMARY, ",".join(missing_summary)
            )
        )
    if summary["run_id"] != SOURCE_RUN_ID:
        raise ValueError(
            "8.5 summary run_id does not match SOURCE_RUN_ID: {0}!={1}".format(
                summary["run_id"], SOURCE_RUN_ID
            )
        )
    if summary["expression_matrix_scope"] != EXPECTED_EXPRESSION_MATRIX_SCOPE:
        raise ValueError(
            "unexpected 8.5 expression_matrix_scope: {0}".format(
                summary["expression_matrix_scope"]
            )
        )
    if summary["feature_origin_rule"] != EXPECTED_FEATURE_ORIGIN_RULE:
        raise ValueError(
            "unexpected 8.5 feature_origin_rule: {0}".format(
                summary["feature_origin_rule"]
            )
        )
    if source_completion.get("run_id") != SOURCE_RUN_ID:
        raise ValueError(
            "8.5 completion run_id does not match SOURCE_RUN_ID: {0}!={1}".format(
                source_completion.get("run_id", ""), SOURCE_RUN_ID
            )
        )
    if source_completion.get("execution_status") != "success":
        raise ValueError(
            "8.5 completion marker is not successful: {0}".format(
                source_completion.get("execution_status", "")
            )
        )

    valid_count = inspect_block_table(
        VALID_BLOCK_CATALOGUE, VALID_BLOCK_REQUIRED, "valid"
    )
    invalid_count = inspect_block_table(
        INVALID_BLOCK_QC, INVALID_BLOCK_REQUIRED, "invalid"
    )
    if valid_count != safe_int(summary["valid_blocks"], -1):
        raise ValueError(
            "valid block count mismatch: table={0}; summary={1}".format(
                valid_count, summary["valid_blocks"]
            )
        )
    if invalid_count != safe_int(summary["invalid_blocks_qc_only"], -1):
        raise ValueError(
            "invalid block count mismatch: table={0}; summary={1}".format(
                invalid_count, summary["invalid_blocks_qc_only"]
            )
        )
    return variant_fields, summary, valid_count, invalid_count


def open_gzip_tsv_writer(
    path: Path, fieldnames: Sequence[str]
) -> Tuple[TextIO, csv.DictWriter]:
    handle: TextIO = gzip.open(path, "wt", newline="")
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=list(fieldnames),
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    return handle, writer


def write_bed_row(handle: TextIO, row: Mapping[str, object], ppa: float) -> None:
    chrom = clean(row.get("chrom"))
    pos = safe_int(row.get("pos"), -1)
    if not chrom or pos <= 0:
        raise ValueError(
            "invalid BED coordinate: chrom={0}; pos={1}; marker={2}".format(
                chrom, row.get("pos", ""), row.get("marker_id", "")
            )
        )
    name = "{0}|{1}|{2}|{3}".format(
        clean(row.get("marker_id")),
        clean(row.get("tissue")),
        clean(row.get("feature_id")),
        clean(row.get("catalogue_role")),
    )
    score = max(0, min(1000, int(round(ppa * 1000))))
    fields = [
        chrom,
        str(pos - 1),
        str(pos),
        name,
        str(score),
        ".",
        format_float(ppa),
        clean(row.get("marker_id")),
        clean(row.get("ref")),
        clean(row.get("alt")),
        clean(row.get("tissue")),
        clean(row.get("feature_id")),
        clean(row.get("block_id")),
        clean(row.get("ppa_tier")),
        clean(row.get("catalogue_role")),
    ]
    handle.write("\t".join(fields) + "\n")


def update_priority_recurrence(
    grouped: Dict[Tuple[str, str, str, str, str], Dict[str, object]],
    row: Mapping[str, str],
    ppa: float,
) -> None:
    key = (
        clean(row.get("marker_id")),
        clean(row.get("chrom")),
        clean(row.get("pos")),
        clean(row.get("ref")),
        clean(row.get("alt")),
    )
    if not all(key):
        raise ValueError("priority row has an incomplete marker key: {0}".format(key))
    bucket = grouped.get(key)
    if bucket is None:
        bucket = {
            "catalogue_row_count": 0,
            "block_ids": set(),
            "feature_ids": set(),
            "feature_origins": set(),
            "tissues": set(),
            "top_block_ids": set(),
            "credible_block_ids": set(),
            "max_ppa": ppa,
            "max_ppa_tier": clean(row.get("ppa_tier")),
            "max_ppa_tissue": clean(row.get("tissue")),
            "max_ppa_feature_id": clean(row.get("feature_id")),
        }
        grouped[key] = bucket
    bucket["catalogue_row_count"] = int(bucket["catalogue_row_count"]) + 1
    bucket["block_ids"].add(clean(row.get("block_id")))
    bucket["feature_ids"].add(clean(row.get("feature_id")))
    bucket["feature_origins"].add(clean(row.get("feature_origin")))
    bucket["tissues"].add(clean(row.get("tissue")))
    if clean(row.get("is_top_variant")) == "1":
        bucket["top_block_ids"].add(clean(row.get("block_id")))
    if clean(row.get("is_credible_set_member")) == "1":
        bucket["credible_block_ids"].add(clean(row.get("block_id")))
    if ppa > float(bucket["max_ppa"]):
        bucket["max_ppa"] = ppa
        bucket["max_ppa_tier"] = clean(row.get("ppa_tier"))
        bucket["max_ppa_tissue"] = clean(row.get("tissue"))
        bucket["max_ppa_feature_id"] = clean(row.get("feature_id"))


def join_sorted(values: Iterable[str]) -> str:
    return ";".join(sorted({clean(value) for value in values if clean(value)}))


def build_priority_recurrence_rows(
    grouped: Mapping[Tuple[str, str, str, str, str], Mapping[str, object]]
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for key, bucket in grouped.items():
        marker_id, chrom, pos, ref, alt = key
        rows.append(
            {
                "marker_id": marker_id,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "catalogue_row_count": bucket["catalogue_row_count"],
                "block_count": len(bucket["block_ids"]),
                "feature_count": len(bucket["feature_ids"]),
                "tissue_count": len(bucket["tissues"]),
                "top_variant_block_count": len(bucket["top_block_ids"]),
                "credible_set_block_count": len(bucket["credible_block_ids"]),
                "max_ppa": format_float(float(bucket["max_ppa"])),
                "max_ppa_tier": bucket["max_ppa_tier"],
                "max_ppa_tissue": bucket["max_ppa_tissue"],
                "max_ppa_feature_id": bucket["max_ppa_feature_id"],
                "tissues": join_sorted(bucket["tissues"]),
                "feature_ids": join_sorted(bucket["feature_ids"]),
                "feature_origins": join_sorted(bucket["feature_origins"]),
                "block_ids": join_sorted(bucket["block_ids"]),
                "interpretation_note": (
                    "marker recurrence among PPA >= 0.1 exploratory no-GWAS "
                    "expression+ASE candidates; not a phenotype causal variant"
                ),
            }
        )
    rows.sort(
        key=lambda row: (
            -safe_float(row["max_ppa"]),
            -safe_int(row["block_count"]),
            clean(row["marker_id"]),
        )
    )
    return rows


def write_priority_recurrence(
    grouped: Mapping[Tuple[str, str, str, str, str], Mapping[str, object]]
) -> int:
    fields = [
        "marker_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "catalogue_row_count",
        "block_count",
        "feature_count",
        "tissue_count",
        "top_variant_block_count",
        "credible_set_block_count",
        "max_ppa",
        "max_ppa_tier",
        "max_ppa_tissue",
        "max_ppa_feature_id",
        "tissues",
        "feature_ids",
        "feature_origins",
        "block_ids",
        "interpretation_note",
    ]
    rows = build_priority_recurrence_rows(grouped)
    handle, writer = open_gzip_tsv_writer(PRIORITY_MARKER_RECURRENCE_TSV, fields)
    try:
        writer.writerows(rows)
    finally:
        handle.close()
    return len(rows)


def write_export_summary(
    source_summary: Mapping[str, str],
    counts: Mapping[str, int],
    unique_priority_markers: int,
    valid_blocks: int,
    invalid_blocks: int,
) -> None:
    rows = [
        {"key": "run_id", "value": RUN_ID},
        {"key": "source_run_id", "value": SOURCE_RUN_ID},
        {"key": "source_variant_catalogue", "value": str(VARIANT_CATALOGUE)},
        {"key": "source_valid_block_catalogue", "value": str(VALID_BLOCK_CATALOGUE)},
        {"key": "source_invalid_block_qc", "value": str(INVALID_BLOCK_QC)},
        {"key": "source_catalogue_summary", "value": str(CATALOGUE_SUMMARY)},
        {
            "key": "source_execution_completion_marker",
            "value": str(SOURCE_COMPLETION_MARKER),
        },
        {
            "key": "expression_matrix_scope",
            "value": source_summary["expression_matrix_scope"],
        },
        {
            "key": "feature_origin_rule",
            "value": source_summary["feature_origin_rule"],
        },
        {"key": "all_interpretable_variant_rows", "value": counts["all"]},
        {"key": "priority_variant_rows_ppa_ge_0.1", "value": counts["priority"]},
        {"key": "high_confidence_variant_rows_ppa_ge_0.5", "value": counts["high"]},
        {"key": "strong_variant_rows_ppa_ge_0.8", "value": counts["strong"]},
        {"key": "very_high_variant_rows_ppa_ge_0.95", "value": counts["very_high"]},
        {"key": "top_variant_rows", "value": counts["top"]},
        {"key": "credible_set_variant_rows", "value": counts["credible"]},
        {
            "key": "unique_priority_markers_ppa_ge_0.1",
            "value": unique_priority_markers,
        },
        {"key": "valid_blocks", "value": valid_blocks},
        {"key": "invalid_blocks_qc_only", "value": invalid_blocks},
        {
            "key": "association_fdr_screening",
            "value": source_summary["association_fdr_screening"],
        },
        {
            "key": "null_causal_model",
            "value": source_summary["null_causal_model"],
        },
        {
            "key": "interpretation_scope",
            "value": source_summary["interpretation_scope"],
        },
        {
            "key": "not_interpretable_as",
            "value": source_summary["not_interpretable_as"],
        },
    ]
    with EXPORT_SUMMARY_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["key", "value"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def compare_count(
    observed: int, summary: Mapping[str, str], summary_key: str, label: str
) -> None:
    expected = safe_int(summary.get(summary_key), -1)
    if observed != expected:
        raise ValueError(
            "{0} count mismatch: observed={1}; 8.5_summary={2}".format(
                label, observed, expected
            )
        )


def main() -> None:
    log("[start] export five dogs no-GWAS expression + ASE catalogue delivery file")
    log("[config] SOURCE_RUN_ID={0}".format(SOURCE_RUN_ID))
    log("[config] VARIANT_CATALOGUE={0}".format(VARIANT_CATALOGUE))
    log("[config] RUN_DIR={0}".format(RUN_DIR))
    log("[config] MASTER_LOG={0}".format(MASTER_LOG))
    log(
        "[config] concurrent_tasks={0}; threads_per_task={1}; "
        "resource_baseline=preserved_serial_script".format(
            CONCURRENT_TASKS, THREADS_PER_TASK
        )
    )

    variant_fields, source_summary, valid_blocks, invalid_blocks = validate_inputs()
    log(
        "[stage:validate] input check completed; valid_blocks={0}; invalid_blocks={1}".format(
            valid_blocks, invalid_blocks
        )
    )

    RUN_DIR.mkdir(parents=True)
    tsv_paths = {
        "priority": PRIORITY_VARIANTS_TSV,
        "high": HIGH_CONFIDENCE_TSV,
        "strong": STRONG_TSV,
        "very_high": VERY_HIGH_TSV,
        "top": TOP_VARIANTS_TSV,
        "credible": CREDIBLE_SET_TSV,
    }
    tsv_handles: Dict[str, TextIO] = {}
    tsv_writers: Dict[str, csv.DictWriter] = {}
    for label, path in tsv_paths.items():
        handle, writer = open_gzip_tsv_writer(path, variant_fields)
        tsv_handles[label] = handle
        tsv_writers[label] = writer

    bed_paths = {
        "all": ALL_VARIANTS_BED,
        "priority": PRIORITY_BED,
        "high": HIGH_CONFIDENCE_BED,
        "strong": STRONG_BED,
        "very_high": VERY_HIGH_BED,
        "top": TOP_VARIANTS_BED,
        "credible": CREDIBLE_SET_BED,
    }
    bed_handles = {
        label: gzip.open(path, "wt", newline="") for label, path in bed_paths.items()
    }

    counts = {
        "all": 0,
        "priority": 0,
        "high": 0,
        "strong": 0,
        "very_high": 0,
        "top": 0,
        "credible": 0,
    }
    priority_recurrence: Dict[
        Tuple[str, str, str, str, str], Dict[str, object]
    ] = {}

    log("[stage:stream] start streaming read 8.5 variant catalogue")
    try:
        with gzip.open(VARIANT_CATALOGUE, "rt", newline="") as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            for row in reader:
                counts["all"] += 1
                ppa = safe_float(row.get("ppa"))
                if (
                    clean(row.get("result_status")) != "valid"
                    or clean(row.get("result_interpretable")) != "1"
                    or clean(row.get("ppa_status")) not in ALLOWED_PPA_STATUS
                    or clean(row.get("feature_origin"))
                    not in ALLOWED_FEATURE_ORIGINS
                    or clean(row.get("expression_matrix_scope"))
                    != EXPECTED_EXPRESSION_MATRIX_SCOPE
                    or not math.isfinite(ppa)
                    or not 0.0 <= ppa <= 1.0
                ):
                    raise ValueError(
                        "8.5 variant catalogue contains an invalid export row "
                        "at data row {0}: block_id={1}; marker_id={2}; "
                        "result_status={3}; result_interpretable={4}; "
                        "ppa_status={5}; feature_origin={6}; "
                        "expression_matrix_scope={7}; ppa={8}".format(
                            counts["all"],
                            clean(row.get("block_id")),
                            clean(row.get("marker_id")),
                            clean(row.get("result_status")),
                            clean(row.get("result_interpretable")),
                            clean(row.get("ppa_status")),
                            clean(row.get("feature_origin")),
                            clean(row.get("expression_matrix_scope")),
                            clean(row.get("ppa")),
                        )
                    )

                top_token = clean(row.get("is_top_variant"))
                credible_token = clean(row.get("is_credible_set_member"))
                observed_role = clean(row.get("catalogue_role"))
                expected_role = EXPECTED_CATALOGUE_ROLE_BY_FLAGS.get(
                    (top_token, credible_token)
                )
                if expected_role is None or observed_role != expected_role:
                    raise ValueError(
                        "8.5 variant catalogue contains inconsistent top/credible "
                        "flags at data row {0}: block_id={1}; marker_id={2}; "
                        "is_top_variant={3}; is_credible_set_member={4}; "
                        "catalogue_role={5}; expected_role={6}".format(
                            counts["all"],
                            clean(row.get("block_id")),
                            clean(row.get("marker_id")),
                            top_token,
                            credible_token,
                            observed_role,
                            expected_role or "binary_0_or_1_flags_required",
                        )
                    )

                write_bed_row(bed_handles["all"], row, ppa)

                is_top = top_token == "1"
                is_credible = credible_token == "1"
                if is_top:
                    counts["top"] += 1
                    tsv_writers["top"].writerow(row)
                    write_bed_row(bed_handles["top"], row, ppa)
                if is_credible:
                    counts["credible"] += 1
                    tsv_writers["credible"].writerow(row)
                    write_bed_row(bed_handles["credible"], row, ppa)

                if ppa >= PRIORITY_THRESHOLD:
                    counts["priority"] += 1
                    tsv_writers["priority"].writerow(row)
                    write_bed_row(bed_handles["priority"], row, ppa)
                    update_priority_recurrence(priority_recurrence, row, ppa)
                if ppa >= HIGH_CONFIDENCE_THRESHOLD:
                    counts["high"] += 1
                    tsv_writers["high"].writerow(row)
                    write_bed_row(bed_handles["high"], row, ppa)
                if ppa >= STRONG_THRESHOLD:
                    counts["strong"] += 1
                    tsv_writers["strong"].writerow(row)
                    write_bed_row(bed_handles["strong"], row, ppa)
                if ppa >= VERY_HIGH_THRESHOLD:
                    counts["very_high"] += 1
                    tsv_writers["very_high"].writerow(row)
                    write_bed_row(bed_handles["very_high"], row, ppa)

                if counts["all"] % PROGRESS_EVERY == 0:
                    log(
                        "[stage:stream] processed={0}; priority={1}; "
                        "top={2}; credible={3}".format(
                            counts["all"],
                            counts["priority"],
                            counts["top"],
                            counts["credible"],
                        )
                    )
    finally:
        for handle in tsv_handles.values():
            handle.close()
        for handle in bed_handles.values():
            handle.close()

    log("[stage:validate] check 8.6 row count and 8.5 summary")
    compare_count(
        counts["all"], source_summary, "interpretable_variant_rows", "all variants"
    )
    compare_count(counts["top"], source_summary, "top_variant_rows", "top variants")
    compare_count(
        counts["credible"],
        source_summary,
        "credible_set_variant_rows",
        "credible-set variants",
    )
    compare_count(
        counts["priority"],
        source_summary,
        "ppa_ge_0.1_variant_rows",
        "PPA >= 0.1 variants",
    )
    compare_count(
        counts["high"],
        source_summary,
        "ppa_ge_0.5_variant_rows",
        "PPA >= 0.5 variants",
    )
    compare_count(
        counts["strong"],
        source_summary,
        "ppa_ge_0.8_variant_rows",
        "PPA >= 0.8 variants",
    )
    compare_count(
        counts["very_high"],
        source_summary,
        "ppa_ge_0.95_variant_rows",
        "PPA >= 0.95 variants",
    )

    log("[stage:write] copy all records gzip catalogue and write priority recurrence")
    shutil.copyfile(VARIANT_CATALOGUE, ALL_VARIANTS_TSV)
    unique_priority_markers = write_priority_recurrence(priority_recurrence)
    write_export_summary(
        source_summary,
        counts,
        unique_priority_markers,
        valid_blocks,
        invalid_blocks,
    )

    log("[summary] all_interpretable_variant_rows={0}".format(counts["all"]))
    log("[summary] priority_ppa_ge_0.1={0}".format(counts["priority"]))
    log("[summary] high_confidence_ppa_ge_0.5={0}".format(counts["high"]))
    log("[summary] strong_ppa_ge_0.8={0}".format(counts["strong"]))
    log("[summary] very_high_ppa_ge_0.95={0}".format(counts["very_high"]))
    log("[summary] top_variant_rows={0}".format(counts["top"]))
    log("[summary] credible_set_variant_rows={0}".format(counts["credible"]))
    log("[summary] unique_priority_markers={0}".format(unique_priority_markers))
    log("[output] {0}".format(RUN_DIR))
    COMPLETION_MARKER.write_text(
        "run_id={0}\nexecution_status=success\nfinished_at={1}\n".format(
            RUN_ID, datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ),
        encoding="utf-8",
    )
    log(
        "[done] five dogs no-GWAS expression + ASE catalogue delivery file export completed; "
        "task completed / done"
    )
    log("execution_status=success")
    log(
        "[reminder] execution completion does not establish result completeness, format, counts or biological quality have passed review; "
        "please separately check output."
    )


if __name__ == "__main__":
    try:
        initialize_execution_logging()
        main()
    except BaseException as exc:
        if MASTER_LOG_HANDLE is not None:
            log("[ERROR] {0}: {1}".format(type(exc).__name__, clean(exc)))
            traceback.print_exc()
            log("execution_status=failed")
        else:
            traceback.print_exc()
        raise
    finally:
        if MASTER_LOG_HANDLE is not None:
            MASTER_LOG_HANDLE.close()
