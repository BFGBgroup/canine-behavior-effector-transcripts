#!/usr/bin/python3.10


from __future__ import annotations

import csv
import gzip
import math
import os
import sys
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, TextIO, Tuple


SOURCE_RUN_ID = os.environ.get("SOURCE_RUN_ID", "20260827_224912")
EXPORT_ROOT = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.6_no_gwas_export_expression_ase_catalogue_deliverables'
)
EXPORT_DIR = EXPORT_ROOT / SOURCE_RUN_ID
SOURCE_COMPLETION_MARKER = EXPORT_DIR / "execution_completion.marker"

OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.7_no_gwas_final_catalogue_qc_report'
)
RUN_ID = os.environ.get("RUN_ID", datetime.now().strftime("%Y%m%d_%H%M%S"))
RUN_DIR = OUTPUT_ROOT / RUN_ID
COMPLETION_MARKER = RUN_DIR / "execution_completion.marker"

SCRIPT_DIR = Path(
    '/path/to/repository/scripts/10_molecular_finemapping/8.7_no_gwas_final_catalogue_qc_report'
)
LOG_DIR = SCRIPT_DIR / "log"
MASTER_LOG = LOG_DIR / (
    "8.7.final_qc_no_gwas_expression_ase_catalogue.{0}.log".format(RUN_ID)
)

CONCURRENT_TASKS = 1
THREADS_PER_TASK = 1
MASTER_LOG_HANDLE: Optional[TextIO] = None

ALL_VARIANTS_TSV = (
    EXPORT_DIR / "no_gwas_expression_ase_all_interpretable_variants.tsv.gz"
)
PRIORITY_VARIANTS_TSV = (
    EXPORT_DIR / "no_gwas_expression_ase_priority_variants_ppa_ge_0.1.tsv.gz"
)
HIGH_CONFIDENCE_TSV = (
    EXPORT_DIR / "no_gwas_expression_ase_high_confidence_variants_ppa_ge_0.5.tsv.gz"
)
STRONG_TSV = (
    EXPORT_DIR / "no_gwas_expression_ase_strong_variants_ppa_ge_0.8.tsv.gz"
)
VERY_HIGH_TSV = (
    EXPORT_DIR / "no_gwas_expression_ase_very_high_variants_ppa_ge_0.95.tsv.gz"
)
TOP_VARIANTS_TSV = EXPORT_DIR / "no_gwas_expression_ase_top_variants.tsv.gz"
CREDIBLE_SET_TSV = (
    EXPORT_DIR / "no_gwas_expression_ase_credible_set_variants.tsv.gz"
)
PRIORITY_MARKER_RECURRENCE_TSV = (
    EXPORT_DIR
    / "no_gwas_expression_ase_priority_marker_recurrence_ppa_ge_0.1.tsv.gz"
)
EXPORT_SUMMARY_TSV = EXPORT_DIR / "no_gwas_expression_ase_export_summary.tsv"

ALL_VARIANTS_BED = (
    EXPORT_DIR / "no_gwas_expression_ase_all_interpretable_variants.bed.gz"
)
PRIORITY_BED = (
    EXPORT_DIR / "no_gwas_expression_ase_priority_variants_ppa_ge_0.1.bed.gz"
)
HIGH_CONFIDENCE_BED = (
    EXPORT_DIR
    / "no_gwas_expression_ase_high_confidence_variants_ppa_ge_0.5.bed.gz"
)
STRONG_BED = (
    EXPORT_DIR / "no_gwas_expression_ase_strong_variants_ppa_ge_0.8.bed.gz"
)
VERY_HIGH_BED = (
    EXPORT_DIR / "no_gwas_expression_ase_very_high_variants_ppa_ge_0.95.bed.gz"
)
TOP_VARIANTS_BED = EXPORT_DIR / "no_gwas_expression_ase_top_variants.bed.gz"
CREDIBLE_SET_BED = (
    EXPORT_DIR / "no_gwas_expression_ase_credible_set_variants.bed.gz"
)

QC_SUMMARY = RUN_DIR / "no_gwas_expression_ase_final_qc_summary.tsv"
QC_ISSUES = RUN_DIR / "no_gwas_expression_ase_final_qc_issues.tsv"
FILE_MANIFEST = RUN_DIR / "no_gwas_expression_ase_final_file_manifest.tsv"
TISSUE_SUMMARY = RUN_DIR / "no_gwas_expression_ase_final_tissue_summary.tsv"
FINAL_REPORT_MD = RUN_DIR / "no_gwas_expression_ase_final_report.md"

PROGRESS_EVERY = 5_000_000
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
    "catalogue_role",
    "block_id",
    "tissue",
    "feature_id",
    "feature_origin",
    "expression_matrix_scope",
    "marker_id",
    "chrom",
    "pos",
    "ppa",
    "ppa_tier",
    "result_status",
    "ppa_status",
    "result_interpretable",
    "is_top_variant",
    "is_credible_set_member",
]
RECURRENCE_REQUIRED = [
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
SUMMARY_REQUIRED = [
    "run_id",
    "source_run_id",
    "source_execution_completion_marker",
    "expression_matrix_scope",
    "feature_origin_rule",
    "all_interpretable_variant_rows",
    "priority_variant_rows_ppa_ge_0.1",
    "high_confidence_variant_rows_ppa_ge_0.5",
    "strong_variant_rows_ppa_ge_0.8",
    "very_high_variant_rows_ppa_ge_0.95",
    "top_variant_rows",
    "credible_set_variant_rows",
    "unique_priority_markers_ppa_ge_0.1",
    "valid_blocks",
    "invalid_blocks_qc_only",
    "association_fdr_screening",
    "null_causal_model",
    "interpretation_scope",
    "not_interpretable_as",
]

TSV_EXPORTS: List[Tuple[str, Path, str, str]] = [
    (
        "all_interpretable_variants",
        ALL_VARIANTS_TSV,
        "all_interpretable_variant_rows",
        "all interpretable variant-block catalogue rows",
    ),
    (
        "priority_ppa_ge_0.1",
        PRIORITY_VARIANTS_TSV,
        "priority_variant_rows_ppa_ge_0.1",
        "variant-block rows with PPA >= 0.1",
    ),
    (
        "high_confidence_ppa_ge_0.5",
        HIGH_CONFIDENCE_TSV,
        "high_confidence_variant_rows_ppa_ge_0.5",
        "variant-block rows with PPA >= 0.5",
    ),
    (
        "strong_ppa_ge_0.8",
        STRONG_TSV,
        "strong_variant_rows_ppa_ge_0.8",
        "variant-block rows with PPA >= 0.8",
    ),
    (
        "very_high_ppa_ge_0.95",
        VERY_HIGH_TSV,
        "very_high_variant_rows_ppa_ge_0.95",
        "variant-block rows with PPA >= 0.95",
    ),
    (
        "top_variants",
        TOP_VARIANTS_TSV,
        "top_variant_rows",
        "all tied top-ranked variant-block rows",
    ),
    (
        "credible_set_variants",
        CREDIBLE_SET_TSV,
        "credible_set_variant_rows",
        "95% credible-set variant-block rows",
    ),
    (
        "priority_marker_recurrence_ppa_ge_0.1",
        PRIORITY_MARKER_RECURRENCE_TSV,
        "unique_priority_markers_ppa_ge_0.1",
        "unique marker recurrence among PPA >= 0.1 rows",
    ),
]

BED_EXPORTS: List[Tuple[str, Path, Path]] = [
    ("all_interpretable_variants", ALL_VARIANTS_BED, ALL_VARIANTS_TSV),
    ("priority_ppa_ge_0.1", PRIORITY_BED, PRIORITY_VARIANTS_TSV),
    (
        "high_confidence_ppa_ge_0.5",
        HIGH_CONFIDENCE_BED,
        HIGH_CONFIDENCE_TSV,
    ),
    ("strong_ppa_ge_0.8", STRONG_BED, STRONG_TSV),
    ("very_high_ppa_ge_0.95", VERY_HIGH_BED, VERY_HIGH_TSV),
    ("top_variants", TOP_VARIANTS_BED, TOP_VARIANTS_TSV),
    ("credible_set_variants", CREDIBLE_SET_BED, CREDIBLE_SET_TSV),
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


def write_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(fieldnames),
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def add_issue(
    issues: List[Dict[str, object]],
    check_name: str,
    severity: str,
    detail: str,
) -> None:
    issues.append(
        {"check_name": check_name, "severity": severity, "detail": detail}
    )


def load_export_summary(path: Path) -> Dict[str, str]:
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
    missing = sorted(set(SUMMARY_REQUIRED) - set(values))
    if missing:
        raise ValueError(
            "{0} missing required keys: {1}".format(path, ",".join(missing))
        )
    return values


def load_completion_marker(path: Path) -> Dict[str, str]:
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


def validate_inputs() -> Tuple[Dict[str, str], Dict[str, str]]:
    validate_run_id(SOURCE_RUN_ID, "SOURCE_RUN_ID")
    validate_run_id(RUN_ID, "RUN_ID")
    input_paths = [
        path for _label, path, _key, _purpose in TSV_EXPORTS
    ] + [EXPORT_SUMMARY_TSV, SOURCE_COMPLETION_MARKER]
    input_paths.extend(path for _label, path, _tsv_path in BED_EXPORTS)
    for path in input_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("missing or empty 8.6 input: {0}".format(path))
    if RUN_DIR.exists():
        raise FileExistsError(
            "refusing to reuse existing 8.7 run directory: {0}".format(RUN_DIR)
        )
    summary = load_export_summary(EXPORT_SUMMARY_TSV)
    source_completion = load_completion_marker(SOURCE_COMPLETION_MARKER)
    missing_completion = sorted(
        {"run_id", "execution_status", "finished_at"} - set(source_completion)
    )
    if missing_completion:
        raise ValueError(
            "{0} missing required keys: {1}".format(
                SOURCE_COMPLETION_MARKER, ",".join(missing_completion)
            )
        )
    if summary["run_id"] != SOURCE_RUN_ID:
        raise ValueError(
            "8.6 summary run_id does not match SOURCE_RUN_ID: {0}!={1}".format(
                summary["run_id"], SOURCE_RUN_ID
            )
        )
    if summary["expression_matrix_scope"] != EXPECTED_EXPRESSION_MATRIX_SCOPE:
        raise ValueError(
            "unexpected 8.6 expression_matrix_scope: {0}".format(
                summary["expression_matrix_scope"]
            )
        )
    if summary["feature_origin_rule"] != EXPECTED_FEATURE_ORIGIN_RULE:
        raise ValueError(
            "unexpected 8.6 feature_origin_rule: {0}".format(
                summary["feature_origin_rule"]
            )
        )
    if summary["association_fdr_screening"] != "not_performed":
        raise ValueError(
            "unexpected 8.6 association_fdr_screening: {0}".format(
                summary["association_fdr_screening"]
            )
        )
    if summary["null_causal_model"] != "not_included":
        raise ValueError(
            "unexpected 8.6 null_causal_model: {0}".format(
                summary["null_causal_model"]
            )
        )
    if source_completion.get("run_id") != SOURCE_RUN_ID:
        raise ValueError(
            "8.6 completion run_id does not match SOURCE_RUN_ID: {0}!={1}".format(
                source_completion.get("run_id", ""), SOURCE_RUN_ID
            )
        )
    if source_completion.get("execution_status") != "success":
        raise ValueError(
            "8.6 completion marker is not successful: {0}".format(
                source_completion.get("execution_status", "")
            )
        )
    return summary, source_completion


def expected_ppa_tier(ppa: float) -> str:
    if ppa >= 0.95:
        return "ppa_ge_0.95"
    if ppa >= 0.8:
        return "ppa_0.8_to_lt_0.95"
    if ppa >= 0.5:
        return "ppa_0.5_to_lt_0.8"
    if ppa >= 0.1:
        return "ppa_0.1_to_lt_0.5"
    return "ppa_lt_0.1"


def expected_feature_origin(feature_id: str) -> str:
    if feature_id.startswith("DOG5_NEWGENE_"):
        return "dog5_manual_new_gene"
    if feature_id.startswith("gene-"):
        return "reference_annotation_gene"
    return ""


def init_tissue_row(tissue: str) -> Dict[str, object]:
    return {
        "tissue": tissue,
        "all_interpretable_variants": 0,
        "top_variants": 0,
        "credible_set_variants": 0,
        "ppa_ge_0.1": 0,
        "ppa_ge_0.5": 0,
        "ppa_ge_0.8": 0,
        "ppa_ge_0.95": 0,
    }


def inspect_variant_export(
    label: str,
    path: Path,
    issues: List[Dict[str, object]],
    minimum_ppa: float | None = None,
    required_flag: str | None = None,
    tissue_rows: Dict[str, Dict[str, object]] | None = None,
) -> int:
    count = 0
    bad_status = 0
    bad_ppa = 0
    bad_tier = 0
    bad_flags = 0
    bad_role = 0
    bad_origin = 0
    bad_scope = 0
    below_threshold = 0
    bad_required_flag = 0

    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns(path, list(reader.fieldnames or []), VARIANT_REQUIRED)
        for row in reader:
            count += 1
            ppa = safe_float(row.get("ppa"))
            if (
                clean(row.get("result_status")) != "valid"
                or clean(row.get("result_interpretable")) != "1"
                or clean(row.get("ppa_status")) not in ALLOWED_PPA_STATUS
            ):
                bad_status += 1
            if not math.isfinite(ppa) or not 0.0 <= ppa <= 1.0:
                bad_ppa += 1
            elif clean(row.get("ppa_tier")) != expected_ppa_tier(ppa):
                bad_tier += 1
            if (
                clean(row.get("is_top_variant")) not in {"0", "1"}
                or clean(row.get("is_credible_set_member")) not in {"0", "1"}
            ):
                bad_flags += 1
            else:
                expected_role = EXPECTED_CATALOGUE_ROLE_BY_FLAGS[
                    (
                        clean(row.get("is_top_variant")),
                        clean(row.get("is_credible_set_member")),
                    )
                ]
                if clean(row.get("catalogue_role")) != expected_role:
                    bad_role += 1
            expected_origin = expected_feature_origin(clean(row.get("feature_id")))
            if (
                not expected_origin
                or clean(row.get("feature_origin")) != expected_origin
            ):
                bad_origin += 1
            if (
                clean(row.get("expression_matrix_scope"))
                != EXPECTED_EXPRESSION_MATRIX_SCOPE
            ):
                bad_scope += 1
            if minimum_ppa is not None and (
                not math.isfinite(ppa) or ppa < minimum_ppa
            ):
                below_threshold += 1
            if required_flag is not None and clean(row.get(required_flag)) != "1":
                bad_required_flag += 1

            if tissue_rows is not None and math.isfinite(ppa):
                tissue = clean(row.get("tissue")) or "NA"
                bucket = tissue_rows.setdefault(tissue, init_tissue_row(tissue))
                bucket["all_interpretable_variants"] = (
                    safe_int(bucket["all_interpretable_variants"]) + 1
                )
                if clean(row.get("is_top_variant")) == "1":
                    bucket["top_variants"] = safe_int(bucket["top_variants"]) + 1
                if clean(row.get("is_credible_set_member")) == "1":
                    bucket["credible_set_variants"] = (
                        safe_int(bucket["credible_set_variants"]) + 1
                    )
                if ppa >= 0.1:
                    bucket["ppa_ge_0.1"] = safe_int(bucket["ppa_ge_0.1"]) + 1
                if ppa >= 0.5:
                    bucket["ppa_ge_0.5"] = safe_int(bucket["ppa_ge_0.5"]) + 1
                if ppa >= 0.8:
                    bucket["ppa_ge_0.8"] = safe_int(bucket["ppa_ge_0.8"]) + 1
                if ppa >= 0.95:
                    bucket["ppa_ge_0.95"] = safe_int(bucket["ppa_ge_0.95"]) + 1

            if count % PROGRESS_EVERY == 0:
                log("[stage:tsv] {0}: processed={1}".format(label, count))

    if bad_status:
        add_issue(
            issues,
            "variant_interpretable_status",
            "ERROR",
            "{0}: invalid status rows={1}".format(label, bad_status),
        )
    if bad_ppa:
        add_issue(
            issues,
            "variant_ppa_range",
            "ERROR",
            "{0}: non-finite or out-of-range PPA rows={1}".format(label, bad_ppa),
        )
    if bad_tier:
        add_issue(
            issues,
            "variant_ppa_tier",
            "ERROR",
            "{0}: inconsistent ppa_tier rows={1}".format(label, bad_tier),
        )
    if bad_flags:
        add_issue(
            issues,
            "catalogue_binary_flags",
            "ERROR",
            "{0}: non-binary top/credible flags={1}".format(label, bad_flags),
        )
    if bad_role:
        add_issue(
            issues,
            "catalogue_role",
            "ERROR",
            "{0}: catalogue_role inconsistent with top/credible flags={1}".format(
                label, bad_role
            ),
        )
    if bad_origin:
        add_issue(
            issues,
            "feature_origin",
            "ERROR",
            "{0}: feature_id and feature_origin mismatches={1}".format(
                label, bad_origin
            ),
        )
    if bad_scope:
        add_issue(
            issues,
            "expression_matrix_scope",
            "ERROR",
            "{0}: unexpected expression_matrix_scope rows={1}".format(
                label, bad_scope
            ),
        )
    if below_threshold:
        add_issue(
            issues,
            "threshold_membership",
            "ERROR",
            "{0}: rows below threshold {1}={2}".format(
                label, minimum_ppa, below_threshold
            ),
        )
    if bad_required_flag:
        add_issue(
            issues,
            "catalogue_flag_membership",
            "ERROR",
            "{0}: rows without {1}=1: {2}".format(
                label, required_flag, bad_required_flag
            ),
        )
    return count


def inspect_recurrence_export(
    path: Path, issues: List[Dict[str, object]]
) -> int:
    count = 0
    bad_marker = 0
    bad_ppa = 0
    bad_tier = 0
    bad_counts = 0
    bad_origins = 0
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns(path, list(reader.fieldnames or []), RECURRENCE_REQUIRED)
        for row in reader:
            count += 1
            max_ppa = safe_float(row.get("max_ppa"))
            catalogue_count = safe_int(row.get("catalogue_row_count"), -1)
            block_count = safe_int(row.get("block_count"), -1)
            feature_count = safe_int(row.get("feature_count"), -1)
            tissue_count = safe_int(row.get("tissue_count"), -1)
            top_count = safe_int(row.get("top_variant_block_count"), -1)
            credible_count = safe_int(row.get("credible_set_block_count"), -1)
            origins = {
                token
                for token in clean(row.get("feature_origins")).split(";")
                if token
            }
            if not clean(row.get("marker_id")):
                bad_marker += 1
            if not math.isfinite(max_ppa) or not 0.1 <= max_ppa <= 1.0:
                bad_ppa += 1
            elif clean(row.get("max_ppa_tier")) != expected_ppa_tier(max_ppa):
                bad_tier += 1
            if not origins or not origins.issubset(ALLOWED_FEATURE_ORIGINS):
                bad_origins += 1
            if (
                catalogue_count < 1
                or block_count < 1
                or feature_count < 1
                or tissue_count < 1
                or block_count > catalogue_count
                or feature_count > catalogue_count
                or tissue_count > catalogue_count
                or top_count < 0
                or credible_count < 0
                or top_count > block_count
                or credible_count > block_count
            ):
                bad_counts += 1
    if bad_marker:
        add_issue(
            issues,
            "recurrence_marker_id",
            "ERROR",
            "priority recurrence rows with empty marker_id: {0}".format(bad_marker),
        )
    if bad_ppa:
        add_issue(
            issues,
            "recurrence_max_ppa",
            "ERROR",
            "priority recurrence rows with invalid max_ppa: {0}".format(bad_ppa),
        )
    if bad_tier:
        add_issue(
            issues,
            "recurrence_max_ppa_tier",
            "ERROR",
            "priority recurrence rows with inconsistent max_ppa_tier: {0}".format(
                bad_tier
            ),
        )
    if bad_origins:
        add_issue(
            issues,
            "recurrence_feature_origins",
            "ERROR",
            "priority recurrence rows with unsupported feature_origins: {0}".format(
                bad_origins
            ),
        )
    if bad_counts:
        add_issue(
            issues,
            "recurrence_counts",
            "ERROR",
            "priority recurrence rows with inconsistent counts: {0}".format(
                bad_counts
            ),
        )
    return count


def inspect_bed_export(
    label: str, path: Path, issues: List[Dict[str, object]]
) -> int:
    count = 0
    bad_columns = 0
    bad_coords = 0
    bad_score = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            count += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 15:
                bad_columns += 1
                continue
            start = safe_int(fields[1], -1)
            end = safe_int(fields[2], -1)
            score = safe_int(fields[4], -1)
            if start < 0 or end != start + 1:
                bad_coords += 1
            if score < 0 or score > 1000:
                bad_score += 1
            if count % PROGRESS_EVERY == 0:
                log("[stage:bed] {0}: processed={1}".format(label, count))
    if bad_columns or bad_coords or bad_score:
        add_issue(
            issues,
            "bed_format",
            "ERROR",
            "{0}: bad_columns={1}; bad_coords={2}; bad_score={3}".format(
                label, bad_columns, bad_coords, bad_score
            ),
        )
    return count


def compare_count(
    label: str,
    observed: int,
    expected: int,
    check_name: str,
    issues: List[Dict[str, object]],
) -> None:
    if observed != expected:
        add_issue(
            issues,
            check_name,
            "ERROR",
            "{0}: observed={1}; expected={2}".format(label, observed, expected),
        )


def build_file_manifest(
    row_counts: Mapping[Path, int],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for label, path, summary_key, purpose in TSV_EXPORTS:
        rows.append(
            {
                "file_label": label,
                "file_type": "tsv.gz",
                "path": str(path),
                "data_rows": row_counts[path],
                "bytes": path.stat().st_size,
                "summary_key": summary_key,
                "purpose": purpose,
            }
        )
    rows.append(
        {
            "file_label": "export_summary",
            "file_type": "tsv",
            "path": str(EXPORT_SUMMARY_TSV),
            "data_rows": row_counts[EXPORT_SUMMARY_TSV],
            "bytes": EXPORT_SUMMARY_TSV.stat().st_size,
            "summary_key": "NA",
            "purpose": "8.6 export count summary and interpretation boundary",
        }
    )
    rows.append(
        {
            "file_label": "source_execution_completion_marker",
            "file_type": "marker",
            "path": str(SOURCE_COMPLETION_MARKER),
            "data_rows": row_counts[SOURCE_COMPLETION_MARKER],
            "bytes": SOURCE_COMPLETION_MARKER.stat().st_size,
            "summary_key": "NA",
            "purpose": "8.6 execution completion provenance",
        }
    )
    for label, path, _tsv_path in BED_EXPORTS:
        rows.append(
            {
                "file_label": label,
                "file_type": "bed.gz",
                "path": str(path),
                "data_rows": row_counts[path],
                "bytes": path.stat().st_size,
                "summary_key": "NA",
                "purpose": "single-base BED export for {0}".format(label),
            }
        )
    return rows


def write_markdown_report(
    source_summary: Mapping[str, str],
    issue_count: int,
    tissue_count: int,
) -> None:
    status = "PASS" if issue_count == 0 else "FAIL"
    lines = [
        "# No-GWAS Expression + ASE Catalogue Final QC Report",
        "",
        "## Status",
        "",
        "- final_qc_status: `{0}`".format(status),
        "- issue_count: `{0}`".format(issue_count),
        "- source_8.6_run_id: `{0}`".format(SOURCE_RUN_ID),
        "- tissue_count: `{0}`".format(tissue_count),
        "- expression_matrix_scope: `{0}`".format(
            source_summary["expression_matrix_scope"]
        ),
        "- feature_origin_rule: `{0}`".format(
            source_summary["feature_origin_rule"]
        ),
        "",
        "## Core Counts",
        "",
        "- all_interpretable_variant_rows: `{0}`".format(
            source_summary["all_interpretable_variant_rows"]
        ),
        "- priority_variant_rows_ppa_ge_0.1: `{0}`".format(
            source_summary["priority_variant_rows_ppa_ge_0.1"]
        ),
        "- high_confidence_variant_rows_ppa_ge_0.5: `{0}`".format(
            source_summary["high_confidence_variant_rows_ppa_ge_0.5"]
        ),
        "- strong_variant_rows_ppa_ge_0.8: `{0}`".format(
            source_summary["strong_variant_rows_ppa_ge_0.8"]
        ),
        "- very_high_variant_rows_ppa_ge_0.95: `{0}`".format(
            source_summary["very_high_variant_rows_ppa_ge_0.95"]
        ),
        "- top_variant_rows: `{0}`".format(source_summary["top_variant_rows"]),
        "- credible_set_variant_rows: `{0}`".format(
            source_summary["credible_set_variant_rows"]
        ),
        "- unique_priority_markers_ppa_ge_0.1: `{0}`".format(
            source_summary["unique_priority_markers_ppa_ge_0.1"]
        ),
        "",
        "## Interpretation Boundary",
        "",
        "This report freezes an exploratory no-GWAS expression+ASE fine-mapping catalogue.",
        "Association/FDR screening and a null-causal model were not applied.",
        "A PASS result confirms the current file structure, counts, and internal consistency only.",
        "It does not establish association-significant or phenotype-causal variants or genes.",
        "",
    ]
    FINAL_REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    log("[start] five dogs no-GWAS expression + ASE catalogue final QC")
    log("[config] SOURCE_RUN_ID={0}".format(SOURCE_RUN_ID))
    log("[config] EXPORT_DIR={0}".format(EXPORT_DIR))
    log("[config] RUN_DIR={0}".format(RUN_DIR))
    log("[config] MASTER_LOG={0}".format(MASTER_LOG))
    log(
        "[config] concurrent_tasks={0}; threads_per_task={1}; "
        "resource_baseline=preserved_serial_script".format(
            CONCURRENT_TASKS, THREADS_PER_TASK
        )
    )

    source_summary, source_completion = validate_inputs()
    log("[stage:validate] 8.6 input file, summary and completion marker contract check completed")

    issues: List[Dict[str, object]] = []
    row_counts: Dict[Path, int] = {}
    tissue_buckets: Dict[str, Dict[str, object]] = {}

    variant_jobs = [
        (
            "all_interpretable_variants",
            ALL_VARIANTS_TSV,
            None,
            None,
            tissue_buckets,
        ),
        ("priority_ppa_ge_0.1", PRIORITY_VARIANTS_TSV, 0.1, None, None),
        (
            "high_confidence_ppa_ge_0.5",
            HIGH_CONFIDENCE_TSV,
            0.5,
            None,
            None,
        ),
        ("strong_ppa_ge_0.8", STRONG_TSV, 0.8, None, None),
        ("very_high_ppa_ge_0.95", VERY_HIGH_TSV, 0.95, None, None),
        ("top_variants", TOP_VARIANTS_TSV, None, "is_top_variant", None),
        (
            "credible_set_variants",
            CREDIBLE_SET_TSV,
            None,
            "is_credible_set_member",
            None,
        ),
    ]
    log("[stage:tsv] start streaming check variant TSV exports")
    for label, path, threshold, required_flag, tissue_rows in variant_jobs:
        log("[stage:tsv] checking {0}".format(path.name))
        row_counts[path] = inspect_variant_export(
            label,
            path,
            issues,
            minimum_ppa=threshold,
            required_flag=required_flag,
            tissue_rows=tissue_rows,
        )

    log("[stage:tsv] checking priority marker recurrence")
    row_counts[PRIORITY_MARKER_RECURRENCE_TSV] = inspect_recurrence_export(
        PRIORITY_MARKER_RECURRENCE_TSV, issues
    )
    row_counts[EXPORT_SUMMARY_TSV] = len(source_summary)
    row_counts[SOURCE_COMPLETION_MARKER] = len(source_completion)

    for label, path, summary_key, _purpose in TSV_EXPORTS:
        compare_count(
            label,
            row_counts[path],
            safe_int(source_summary[summary_key], -1),
            "tsv_count_match",
            issues,
        )

    log("[stage:bed] start streaming check BED exports")
    for label, bed_path, tsv_path in BED_EXPORTS:
        log("[stage:bed] checking {0}".format(bed_path.name))
        row_counts[bed_path] = inspect_bed_export(label, bed_path, issues)
        compare_count(
            label,
            row_counts[bed_path],
            row_counts[tsv_path],
            "bed_count_match",
            issues,
        )

    log("[stage:manifest] record file path, row count and size")
    manifest_rows = build_file_manifest(row_counts)
    tissue_rows = [
        tissue_buckets[tissue] for tissue in sorted(tissue_buckets)
    ]
    final_status = "PASS" if not issues else "FAIL"

    summary_rows: List[Dict[str, object]] = [
        {"key": "final_qc_status", "value": final_status},
        {"key": "issue_count", "value": len(issues)},
        {"key": "run_id", "value": RUN_ID},
        {"key": "source_export_dir", "value": str(EXPORT_DIR)},
        {
            "key": "source_8.6_execution_completion_marker",
            "value": str(SOURCE_COMPLETION_MARKER),
        },
        {
            "key": "source_8.6_execution_status",
            "value": source_completion["execution_status"],
        },
        {"key": "final_report_dir", "value": str(RUN_DIR)},
    ]
    summary_rows.extend(
        {
            "key": "source_8.6_{0}".format(key),
            "value": value,
        }
        for key, value in source_summary.items()
    )
    summary_rows.extend(
        {
            "key": "observed_rows_{0}".format(label),
            "value": row_counts[path],
        }
        for label, path, _summary_key, _purpose in TSV_EXPORTS
    )
    summary_rows.extend(
        {
            "key": "observed_bed_rows_{0}".format(label),
            "value": row_counts[path],
        }
        for label, path, _tsv_path in BED_EXPORTS
    )

    log("[stage:write] write final QC, manifest, tissue summary and report")
    RUN_DIR.mkdir(parents=True)
    write_tsv(QC_SUMMARY, ["key", "value"], summary_rows)
    write_tsv(QC_ISSUES, ["check_name", "severity", "detail"], issues)
    write_tsv(
        FILE_MANIFEST,
        [
            "file_label",
            "file_type",
            "path",
            "data_rows",
            "bytes",
            "summary_key",
            "purpose",
        ],
        manifest_rows,
    )
    write_tsv(
        TISSUE_SUMMARY,
        [
            "tissue",
            "all_interpretable_variants",
            "top_variants",
            "credible_set_variants",
            "ppa_ge_0.1",
            "ppa_ge_0.5",
            "ppa_ge_0.8",
            "ppa_ge_0.95",
        ],
        tissue_rows,
    )
    write_markdown_report(source_summary, len(issues), len(tissue_rows))

    log("[summary] final_qc_status={0}".format(final_status))
    log("[summary] issue_count={0}".format(len(issues)))
    log("[summary] manifest_files={0}".format(len(manifest_rows)))
    log("[summary] tissues={0}".format(len(tissue_rows)))
    log("[output] {0}".format(RUN_DIR))
    COMPLETION_MARKER.write_text(
        "run_id={0}\nexecution_status=success\nfinal_qc_status={1}\n"
        "finished_at={2}\n".format(
            RUN_ID,
            final_status,
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        ),
        encoding="utf-8",
    )
    log(
        "[done] five dogs no-GWAS expression + ASE catalogue final QC completed; "
        "task completed / done"
    )
    log("execution_status=success")
    log("final_qc_status={0}".format(final_status))
    log(
        "[reminder] execution completion does not establish result completeness, format, counts or biological quality have passed review; "
        "consider together with final_qc_status and issues table separately check."
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
