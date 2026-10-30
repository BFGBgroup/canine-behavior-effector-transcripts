#!/usr/bin/python3.10


from __future__ import annotations

import csv
import gzip
import heapq
import math
import os
import shutil
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, TextIO, Tuple


PREPARE_RUN_ID = "20260821_181104"
PLASMA_RUN_ID = "20260822_233245"
PREPARE_OUT_DIR = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs/20260821_181104'
)
RUN_OUT_DIR = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.3_no_gwas_run_expression_ase_plasma_joint/20260822_233245'
)
SUMMARY_OUT_DIR = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.4_no_gwas_summarize_expression_ase_plasma_joint'
)
SCRIPT_DIR = Path(
    '/path/to/repository/scripts/10_molecular_finemapping/8.4_no_gwas_summarize_expression_ase_plasma_joint'
)

MANIFEST = PREPARE_OUT_DIR / "no_gwas_plasma_candidate_blocks.tsv"
PLASMA_INPUT_ROOT = PREPARE_OUT_DIR / "plasma_inputs"
PLASMA_RESULT_ROOT = RUN_OUT_DIR / "plasma_results"
PREPARE_COMPLETION_MARKER = PREPARE_OUT_DIR / "execution_completion.marker"
RUN_LOG_DIR = Path(
    '/path/to/repository/scripts/10_molecular_finemapping/8.3_no_gwas_run_expression_ase_plasma_joint/log'
)
SUCCESS_LOG = RUN_LOG_DIR / (
    "8.3.no_gwas_expression_ase_plasma_joint.{0}.success.tsv".format(PLASMA_RUN_ID)
)
FAILED_LOG = RUN_LOG_DIR / (
    "8.3.no_gwas_expression_ase_plasma_joint.{0}.failed.tsv".format(PLASMA_RUN_ID)
)
SKIPPED_LOG = RUN_LOG_DIR / (
    "8.3.no_gwas_expression_ase_plasma_joint.{0}.skipped.tsv".format(PLASMA_RUN_ID)
)
RUN_SUMMARY = RUN_OUT_DIR / "8.3.no_gwas_expression_ase_plasma_joint.run_summary.tsv"

PARALLEL_JOBS = 15
THREADS_PER_TASK = 1
RUN_ID = os.environ.get("RUN_ID", datetime.now().strftime("%Y%m%d_%H%M%S"))
SUMMARY_DIR = SUMMARY_OUT_DIR / RUN_ID
PARTS_DIR = SUMMARY_DIR / ".parts"
BLOCK_SUMMARY = SUMMARY_DIR / "no_gwas_plasma_joint_block_summary.tsv"
VALID_BLOCKS = SUMMARY_DIR / "no_gwas_plasma_joint_block_valid.tsv"
INVALID_BLOCKS = SUMMARY_DIR / "no_gwas_plasma_joint_block_invalid.tsv"
TOP_VARIANTS = SUMMARY_DIR / "no_gwas_plasma_joint_top_variants.tsv.gz"
CREDIBLE_SET_VARIANTS = SUMMARY_DIR / "no_gwas_plasma_joint_credible_set_variants.tsv.gz"
VARIANT_SUMMARY = SUMMARY_DIR / "no_gwas_plasma_joint_variant_ppa_summary.tsv.gz"
SUMMARY_QC = SUMMARY_DIR / "no_gwas_plasma_joint_run_qc.tsv"
COMPLETION_MARKER = SUMMARY_DIR / "execution_completion.marker"

LOG_DIR = SCRIPT_DIR / "log"
MASTER_LOG = LOG_DIR / (
    "8.4.summarize_no_gwas_expression_ase_plasma_joint.{0}.log".format(RUN_ID)
)
SHARD_SUCCESS_LOG = LOG_DIR / (
    "8.4.summarize_no_gwas_expression_ase_plasma_joint.{0}.success.tsv".format(RUN_ID)
)
SHARD_FAILED_LOG = LOG_DIR / (
    "8.4.summarize_no_gwas_expression_ase_plasma_joint.{0}.failed.tsv".format(RUN_ID)
)

PPA_BOUND_TOLERANCE = 1e-12
PPA_SUM_TOLERANCE = 1e-6
CREDIBLE_SET_CONFIDENCE = 0.95
CREDIBLE_SET_SUM_TOLERANCE = 1e-6


for thread_variable in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]:
    os.environ[thread_variable] = str(THREADS_PER_TASK)

MASTER_LOG_HANDLE: Optional[TextIO] = None

MANIFEST_REQUIRED = [
    "block_id",
    "tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "expression_gene_id",
    "expression_source",
    "chrom",
    "selected_sample_count",
    "common_marker_count",
    "total_reads",
    "plasma_input_allowed",
    "input_dir",
]
MARKER_REQUIRED = [
    "chrom",
    "pos",
    "marker_id",
    "ref",
    "alt",
    "phase_scope",
    "local_ps_by_sample",
]
SAMPLE_REQUIRED = [
    "sample",
    "tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "expression_gene_id",
    "expression_source",
    "total_exp",
    "chrom",
    "phase_scope",
    "count_source",
    "orientation_status",
    "hap1_count",
    "hap2_count",
    "total_count",
]

BLOCK_FIELDS = [
    "block_id",
    "tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "expression_gene_id",
    "expression_source",
    "chrom",
    "selected_sample_count",
    "common_marker_count",
    "total_reads",
    "observed_sample_count",
    "observed_marker_count",
    "runner_status",
    "runner_exit_code",
    "runner_reason",
    "result_status",
    "result_interpretable",
    "qc_category",
    "ppa_status",
    "finite_ppa_count",
    "nonfinite_ppa_count",
    "ppa_boundary_adjusted_count",
    "ppa_sum",
    "credible_set_size",
    "credible_set_ppa_sum",
    "max_ppa",
    "top_ppa_tie_count",
    "top_marker_is_unique",
    "top_marker_id",
    "top_marker_pos",
    "qc_reason",
    "runner_stderr_last_line",
    "input_dir",
    "result_dir",
]
SHARD_STATUS_FIELDS = [
    "run_id",
    "shard_id",
    "status",
    "exit_code",
    "started_at",
    "finished_at",
    "output_path",
    "task_log",
    "reason",
]
VARIANT_FIELDS = [
    "block_id",
    "tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "expression_gene_id",
    "expression_source",
    "selected_sample_count",
    "common_marker_count",
    "total_reads",
    "chrom",
    "pos",
    "marker_id",
    "ref",
    "alt",
    "phase_scope",
    "local_ps_by_sample",
    "ppa",
    "marker_rank_by_ppa",
    "credible_set_member",
    "result_status",
    "result_interpretable",
]
TOP_FIELDS = [
    "block_id",
    "tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "expression_gene_id",
    "expression_source",
    "selected_sample_count",
    "common_marker_count",
    "total_reads",
    "chrom",
    "pos",
    "marker_id",
    "ref",
    "alt",
    "phase_scope",
    "local_ps_by_sample",
    "ppa",
    "top_ppa_tie_count",
    "top_marker_is_unique",
]


def log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S ")
    print("[{0}] {1}".format(timestamp, message), flush=True)


def current_time() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def initialize_execution_logging() -> None:

    global MASTER_LOG_HANDLE
    if not RUN_ID or Path(RUN_ID).name != RUN_ID or RUN_ID in {".", ".."}:
        raise ValueError("RUN_ID must be one directory name: {0}".format(RUN_ID))
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MASTER_LOG_HANDLE = MASTER_LOG.open("x", buffering=1)
    os.dup2(MASTER_LOG_HANDLE.fileno(), sys.stdout.fileno())
    os.dup2(MASTER_LOG_HANDLE.fileno(), sys.stderr.fileno())
    for path in [SHARD_SUCCESS_LOG, SHARD_FAILED_LOG]:
        with path.open("x", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                fieldnames=SHARD_STATUS_FIELDS,
                lineterminator="\n",
            )
            writer.writeheader()


def append_shard_status(path: Path, row: Mapping[str, object]) -> None:

    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=SHARD_STATUS_FIELDS,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writerow({field: row.get(field, "") for field in SHARD_STATUS_FIELDS})


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(clean(value))
    except (TypeError, ValueError):
        return default


def format_float(value: float) -> str:
    return "NA" if not math.isfinite(value) else "{0:.12g}".format(value)


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError("missing or empty TSV: {0}".format(path))
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def required_columns(path: Path, observed: Sequence[str], required: Sequence[str]) -> None:
    missing = [field for field in required if field not in observed]
    if missing:
        raise ValueError("{0} missing required columns: {1}".format(path, ",".join(missing)))


def unique_rows_by_key(
    path: Path,
    required: Sequence[str],
    key: str = "block_id",
) -> Dict[str, Dict[str, str]]:
    fields, rows = read_tsv(path)
    required_columns(path, fields, list(required) + [key])
    indexed: Dict[str, Dict[str, str]] = {}
    for row in rows:
        item_key = clean(row.get(key))
        if not item_key:
            raise ValueError("{0} contains an empty {1}".format(path, key))
        if item_key in indexed:
            raise ValueError("{0} contains duplicate {1}: {2}".format(path, key, item_key))
        indexed[item_key] = row
    return indexed


def load_ready_manifest(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError("missing or empty TSV: {0}".format(path))
    ready: List[Dict[str, str]] = []
    seen = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns(path, list(reader.fieldnames or []), MANIFEST_REQUIRED)
        for row in reader:
            if clean(row.get("plasma_input_allowed")) != "1":
                continue
            block_id = clean(row.get("block_id"))
            if not block_id:
                raise ValueError("{0} contains an empty block_id".format(path))
            if block_id in seen:
                raise ValueError(
                    "{0} contains duplicate ready block_id: {1}".format(path, block_id)
                )
            seen.add(block_id)
            ready.append({field: clean(row.get(field)) for field in MANIFEST_REQUIRED})
    return ready


def read_numeric_vector(path: Path) -> List[float]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError("missing or empty numeric vector: {0}".format(path))
    values: List[float] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            token = line.strip()
            if not token:
                raise ValueError("{0}:{1} is empty".format(path, line_number))
            try:
                values.append(float(token))
            except ValueError as exc:
                raise ValueError(
                    "{0}:{1} is not numeric: {2}".format(path, line_number, token)
                ) from exc
    return values


def open_tsv_writer(
    path: Path,
    fieldnames: Sequence[str],
    compressed: bool = False,
    write_header: bool = True,
) -> Tuple[TextIO, csv.DictWriter]:
    if compressed:
        handle: TextIO = gzip.open(path, "wt", newline="")
    else:
        handle = path.open("w", newline="")
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=list(fieldnames),
        lineterminator="\n",
        extrasaction="ignore",
    )
    if write_header:
        writer.writeheader()
    return handle, writer


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(fieldnames),
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def last_nonempty_line(path: Path) -> str:
    if not path.is_file():
        return ""
    last = ""
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.strip():
                last = line.strip()
    return last


def common_block_values(row: Mapping[str, str]) -> Dict[str, object]:
    return {
        "block_id": clean(row.get("block_id")),
        "tissue": clean(row.get("tissue")),
        "feature_id": clean(row.get("feature_id")),
        "feature_primary_id": clean(row.get("feature_primary_id")),
        "feature_label": clean(row.get("feature_label")),
        "feature_class": clean(row.get("feature_class")),
        "expression_gene_id": clean(row.get("expression_gene_id")),
        "expression_source": clean(row.get("expression_source")),
        "chrom": clean(row.get("chrom")),
        "selected_sample_count": clean(row.get("selected_sample_count")),
        "common_marker_count": clean(row.get("common_marker_count")),
        "total_reads": clean(row.get("total_reads")),
    }


def make_unprocessed_block(
    row: Mapping[str, str],
    runner_status: str,
    runner_exit_code: str,
    runner_reason: str,
    qc_category: str,
    qc_reason: str,
) -> Dict[str, object]:
    block_id = clean(row.get("block_id"))
    input_dir = Path(clean(row.get("input_dir"))) if clean(row.get("input_dir")) else PLASMA_INPUT_ROOT / block_id
    result_dir = PLASMA_RESULT_ROOT / block_id
    block = common_block_values(row)
    block.update(
        {
            "observed_sample_count": "NA",
            "observed_marker_count": "NA",
            "runner_status": runner_status,
            "runner_exit_code": runner_exit_code or "NA",
            "runner_reason": runner_reason,
            "result_status": "invalid",
            "result_interpretable": 0,
            "qc_category": qc_category,
            "ppa_status": "not_evaluated",
            "finite_ppa_count": "NA",
            "nonfinite_ppa_count": "NA",
            "ppa_boundary_adjusted_count": "NA",
            "ppa_sum": "NA",
            "credible_set_size": "NA",
            "credible_set_ppa_sum": "NA",
            "max_ppa": "NA",
            "top_ppa_tie_count": "NA",
            "top_marker_is_unique": "NA",
            "top_marker_id": "NA",
            "top_marker_pos": "NA",
            "qc_reason": qc_reason,
            "runner_stderr_last_line": last_nonempty_line(result_dir / "run.stderr.log"),
            "input_dir": str(input_dir),
            "result_dir": str(result_dir),
        }
    )
    return block


def validate_status_logs(
    ready_rows: Sequence[Mapping[str, str]],
    success_rows: Mapping[str, Mapping[str, str]],
    failed_rows: Mapping[str, Mapping[str, str]],
    skipped_rows: Mapping[str, Mapping[str, str]],
) -> None:
    ready_by_id = {clean(row.get("block_id")): row for row in ready_rows}
    ready = set(ready_by_id)
    overlap = (set(success_rows) & set(failed_rows)) | (set(success_rows) & set(skipped_rows))
    overlap |= set(failed_rows) & set(skipped_rows)
    if overlap:
        raise ValueError("8.3 status logs overlap for block_id: {0}".format(sorted(overlap)[0]))
    terminal = set(success_rows) | set(failed_rows) | set(skipped_rows)
    unknown = terminal - ready
    if unknown:
        raise ValueError("8.3 status logs contain block_id absent from ready manifest: {0}".format(sorted(unknown)[0]))
    missing = ready - terminal
    if missing:
        raise ValueError("ready manifest block absent from 8.3 status logs: {0}".format(sorted(missing)[0]))

    status_groups = [
        ("success", success_rows, "joint_expression_ase", 0),
        ("failed", failed_rows, "joint_expression_ase", None),
        ("skipped", skipped_rows, "not_run", 0),
    ]
    for expected_status, rows, expected_mode, expected_exit_code in status_groups:
        for block_id, status_row in rows.items():
            manifest_row = ready_by_id[block_id]
            if clean(status_row.get("run_id")) != PLASMA_RUN_ID:
                raise ValueError("8.3 run_id mismatch for block_id: {0}".format(block_id))
            if clean(status_row.get("status")) != expected_status:
                raise ValueError("8.3 status mismatch for block_id: {0}".format(block_id))
            if clean(status_row.get("analysis_mode")) != expected_mode:
                raise ValueError("8.3 analysis_mode mismatch for block_id: {0}".format(block_id))
            if expected_exit_code is not None and safe_int(status_row.get("exit_code"), -1) != expected_exit_code:
                raise ValueError("8.3 exit_code mismatch for block_id: {0}".format(block_id))
            if expected_status == "failed" and safe_int(status_row.get("exit_code"), 0) == 0:
                raise ValueError("failed 8.3 block has zero exit_code: {0}".format(block_id))
            if clean(status_row.get("input_dir")) != clean(manifest_row.get("input_dir")):
                raise ValueError("8.3 input_dir mismatch for block_id: {0}".format(block_id))
            if safe_int(status_row.get("selected_sample_count"), -1) != safe_int(
                manifest_row.get("selected_sample_count"), -2
            ):
                raise ValueError("8.3 selected_sample_count mismatch for block_id: {0}".format(block_id))
            if safe_int(status_row.get("marker_count"), -1) != safe_int(
                manifest_row.get("common_marker_count"), -2
            ):
                raise ValueError("8.3 marker_count mismatch for block_id: {0}".format(block_id))
            expected_output = str(PLASMA_RESULT_ROOT / block_id)
            observed_output = clean(status_row.get("output_path"))
            if expected_status == "skipped":
                if observed_output:
                    raise ValueError("skipped 8.3 block unexpectedly has output_path: {0}".format(block_id))
                if safe_int(status_row.get("total_exp_unique_count"), -1) != 1:
                    raise ValueError("skipped 8.3 block is not constant-expression: {0}".format(block_id))
                if clean(status_row.get("reason")) != "constant_total_expression_not_joint_estimable":
                    raise ValueError("unexpected skipped reason for block_id: {0}".format(block_id))
            elif observed_output != expected_output:
                raise ValueError("8.3 output_path mismatch for block_id: {0}".format(block_id))


def ppa_ranks(ppas: Sequence[float]) -> List[int]:
    rank_by_value: Dict[float, int] = {}
    for index, value in enumerate(sorted(ppas, reverse=True), 1):
        rank_by_value.setdefault(value, index)
    return [rank_by_value[value] for value in ppas]


def process_completed_block(
    row: Mapping[str, str],
    runner_row: Mapping[str, str],
    variant_writer: csv.DictWriter,
    credible_writer: csv.DictWriter,
    top_writer: csv.DictWriter,
) -> Tuple[Dict[str, object], int, int]:
    block_id = clean(row.get("block_id"))
    input_dir = Path(clean(row.get("input_dir"))) if clean(row.get("input_dir")) else PLASMA_INPUT_ROOT / block_id
    result_dir = PLASMA_RESULT_ROOT / block_id
    markers_path = input_dir / "markers.tsv"
    samples_path = input_dir / "samples.tsv"
    ppas_path = result_dir / "ppas.txt"
    cset_path = result_dir / "cset.txt"

    block = common_block_values(row)
    block.update(
        {
            "runner_status": clean(runner_row.get("status")) or "success",
            "runner_exit_code": 0,
            "runner_reason": "",
            "runner_stderr_last_line": last_nonempty_line(result_dir / "run.stderr.log"),
            "input_dir": str(input_dir),
            "result_dir": str(result_dir),
        }
    )

    missing = [
        name
        for name, path in [
            ("markers.tsv", markers_path),
            ("samples.tsv", samples_path),
            ("ppas.txt", ppas_path),
            ("cset.txt", cset_path),
        ]
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        block.update(
            {
                "result_status": "invalid",
                "result_interpretable": 0,
                "qc_category": "invalid_missing_result_file",
                "ppa_status": "not_evaluated",
                "qc_reason": "missing_or_empty_files={0}".format(",".join(missing)),
            }
        )
        return block, 0, 0

    marker_fields, marker_rows = read_tsv(markers_path)
    sample_fields, sample_rows = read_tsv(samples_path)
    required_columns(markers_path, marker_fields, MARKER_REQUIRED)
    required_columns(samples_path, sample_fields, SAMPLE_REQUIRED)

    marker_count = len(marker_rows)
    sample_count = len(sample_rows)
    expected_markers = safe_int(row.get("common_marker_count"), -1)
    expected_samples = safe_int(row.get("selected_sample_count"), -1)
    block["observed_marker_count"] = marker_count
    block["observed_sample_count"] = sample_count

    structural_reasons: List[str] = []
    if marker_count != expected_markers:
        structural_reasons.append(
            "manifest_marker_count_mismatch={0}!={1}".format(expected_markers, marker_count)
        )
    if sample_count != expected_samples:
        structural_reasons.append(
            "manifest_sample_count_mismatch={0}!={1}".format(expected_samples, sample_count)
        )
    runner_markers = safe_int(runner_row.get("marker_count"), -1)
    runner_samples = safe_int(runner_row.get("selected_sample_count"), -1)
    if runner_markers != expected_markers:
        structural_reasons.append(
            "runner_marker_count_mismatch={0}!={1}".format(runner_markers, expected_markers)
        )
    if runner_samples != expected_samples:
        structural_reasons.append(
            "runner_sample_count_mismatch={0}!={1}".format(runner_samples, expected_samples)
        )
    if any(clean(marker.get("chrom")) != clean(row.get("chrom")) for marker in marker_rows):
        structural_reasons.append("marker_chrom_mismatch")
    if any(clean(marker.get("phase_scope")) != clean(marker.get("chrom")) for marker in marker_rows):
        structural_reasons.append("phase_scope_not_chromosome")
    if len({clean(sample.get("sample")) for sample in sample_rows}) != sample_count:
        structural_reasons.append("duplicate_sample_rows")
    if any(clean(sample.get("tissue")) != clean(row.get("tissue")) for sample in sample_rows):
        structural_reasons.append("sample_tissue_mismatch")
    if any(clean(sample.get("feature_id")) != clean(row.get("feature_id")) for sample in sample_rows):
        structural_reasons.append("sample_feature_id_mismatch")
    if any(clean(sample.get("chrom")) != clean(row.get("chrom")) for sample in sample_rows):
        structural_reasons.append("sample_chrom_mismatch")
    if any(clean(sample.get("phase_scope")) != clean(row.get("chrom")) for sample in sample_rows):
        structural_reasons.append("sample_phase_scope_not_chromosome")
    if any(clean(sample.get("orientation_status")) != "consistent" for sample in sample_rows):
        structural_reasons.append("sample_orientation_not_consistent")

    try:
        ppas = read_numeric_vector(ppas_path)
        cset = read_numeric_vector(cset_path)
    except (FileNotFoundError, ValueError) as exc:
        block.update(
            {
                "result_status": "invalid",
                "result_interpretable": 0,
                "qc_category": "invalid_result_parse",
                "ppa_status": "parse_error",
                "qc_reason": clean(exc),
            }
        )
        return block, 0, 0

    if len(ppas) != marker_count:
        structural_reasons.append("ppa_length_mismatch={0}!={1}".format(len(ppas), marker_count))
    if len(cset) != marker_count:
        structural_reasons.append("cset_length_mismatch={0}!={1}".format(len(cset), marker_count))
    if structural_reasons:
        block.update(
            {
                "result_status": "invalid",
                "result_interpretable": 0,
                "qc_category": "invalid_result_structure",
                "ppa_status": "structural_error",
                "finite_ppa_count": sum(math.isfinite(value) for value in ppas),
                "nonfinite_ppa_count": sum(not math.isfinite(value) for value in ppas),
                "ppa_boundary_adjusted_count": "NA",
                "qc_reason": ";".join(structural_reasons),
            }
        )
        return block, 0, 0

    finite_ppa_count = sum(math.isfinite(value) for value in ppas)
    nonfinite_ppa_count = marker_count - finite_ppa_count
    ppa_boundary_adjusted_count = sum(
        math.isfinite(value)
        and (
            -PPA_BOUND_TOLERANCE <= value < 0.0
            or 1.0 < value <= 1.0 + PPA_BOUND_TOLERANCE
        )
        for value in ppas
    )
    ppas = [
        min(1.0, max(0.0, value))
        if math.isfinite(value)
        and -PPA_BOUND_TOLERANCE <= value <= 1.0 + PPA_BOUND_TOLERANCE
        else value
        for value in ppas
    ]
    if finite_ppa_count == marker_count:
        ppa_status = (
            "finite_boundary_adjusted"
            if ppa_boundary_adjusted_count
            else "finite"
        )
    elif finite_ppa_count == 0:
        ppa_status = "all_nonfinite"
    else:
        ppa_status = "partial_nonfinite"

    numerical_reasons: List[str] = []
    if nonfinite_ppa_count:
        numerical_reasons.append("nonfinite_ppa")
    if any(
        math.isfinite(value)
        and (
            value < -PPA_BOUND_TOLERANCE
            or value > 1.0 + PPA_BOUND_TOLERANCE
        )
        for value in ppas
    ):
        numerical_reasons.append("ppa_outside_0_1")
    ppa_sum = sum(value for value in ppas if math.isfinite(value))
    if finite_ppa_count == marker_count and abs(ppa_sum - 1.0) > PPA_SUM_TOLERANCE:
        numerical_reasons.append("ppa_sum_not_1")
    if any(not math.isfinite(value) for value in cset):
        numerical_reasons.append("nonfinite_cset")
    if any(math.isfinite(value) and value not in (0.0, 1.0) for value in cset):
        numerical_reasons.append("cset_not_binary")
    cset_is_binary = all(math.isfinite(value) and value in (0.0, 1.0) for value in cset)
    credible_indices_for_qc = [
        index for index, value in enumerate(cset) if value == 1.0
    ]
    if finite_ppa_count == marker_count and cset_is_binary:
        if not credible_indices_for_qc:
            numerical_reasons.append("empty_credible_set")
        else:
            credible_ppa_for_qc = sum(ppas[index] for index in credible_indices_for_qc)
            if (
                credible_ppa_for_qc
                < CREDIBLE_SET_CONFIDENCE - CREDIBLE_SET_SUM_TOLERANCE
            ):
                numerical_reasons.append("credible_set_ppa_below_confidence")

    interpretable = 0 if numerical_reasons else 1
    result_status = "valid" if interpretable else "invalid"
    block.update(
        {
            "result_status": result_status,
            "result_interpretable": interpretable,
            "qc_category": "valid" if interpretable else "invalid_result_numerical",
            "ppa_status": ppa_status,
            "finite_ppa_count": finite_ppa_count,
            "nonfinite_ppa_count": nonfinite_ppa_count,
            "ppa_boundary_adjusted_count": ppa_boundary_adjusted_count,
            "ppa_sum": format_float(ppa_sum) if finite_ppa_count else "NA",
            "qc_reason": ";".join(numerical_reasons),
        }
    )

    ranks: List[object]
    if interpretable:
        ranks = ppa_ranks(ppas)
        max_ppa = max(ppas)
        top_indices = [index for index, value in enumerate(ppas) if value == max_ppa]
        credible_indices = credible_indices_for_qc
        credible_set_ppa_sum = sum(ppas[index] for index in credible_indices)
        top_idx = top_indices[0]
        top_marker = marker_rows[top_idx]
        block.update(
            {
                "credible_set_size": len(credible_indices),
                "credible_set_ppa_sum": format_float(credible_set_ppa_sum),
                "max_ppa": format_float(max_ppa),
                "top_ppa_tie_count": len(top_indices),
                "top_marker_is_unique": 1 if len(top_indices) == 1 else 0,
                "top_marker_id": clean(top_marker.get("marker_id")),
                "top_marker_pos": clean(top_marker.get("pos")),
            }
        )
    else:
        ranks = ["NA"] * marker_count
        top_indices = []
        credible_indices = []
        block.update(
            {
                "credible_set_size": "NA",
                "credible_set_ppa_sum": "NA",
                "max_ppa": "NA",
                "top_ppa_tie_count": "NA",
                "top_marker_is_unique": "NA",
                "top_marker_id": "NA",
                "top_marker_pos": "NA",
            }
        )

    common = common_block_values(row)
    credible_index_set = set(credible_indices)
    for index, marker in enumerate(marker_rows):
        variant = {
            **common,
            **marker,
            "ppa": format_float(ppas[index]),
            "marker_rank_by_ppa": ranks[index],
            "credible_set_member": 1 if index in credible_index_set else 0 if interpretable else "NA",
            "result_status": result_status,
            "result_interpretable": interpretable,
        }
        variant_writer.writerow({field: variant.get(field, "") for field in VARIANT_FIELDS})
        if index in credible_index_set:
            credible_writer.writerow({field: variant.get(field, "") for field in VARIANT_FIELDS})

    for index in top_indices:
        top = {
            **common,
            **marker_rows[index],
            "ppa": format_float(ppas[index]),
            "top_ppa_tie_count": len(top_indices),
            "top_marker_is_unique": 1 if len(top_indices) == 1 else 0,
        }
        top_writer.writerow({field: top.get(field, "") for field in TOP_FIELDS})

    return block, marker_count, len(credible_indices)


def partition_ready_rows(
    ready_rows: Sequence[Dict[str, str]],
) -> List[List[Dict[str, str]]]:

    shards: List[List[Dict[str, str]]] = [[] for _ in range(PARALLEL_JOBS)]
    loads = [(0, shard_id) for shard_id in range(PARALLEL_JOBS)]
    heapq.heapify(loads)
    indexed_rows = list(enumerate(ready_rows))
    indexed_rows.sort(
        key=lambda item: (-safe_int(item[1].get("common_marker_count")), item[0])
    )

    for manifest_order, row in indexed_rows:
        current_load, shard_id = heapq.heappop(loads)
        shard_row = dict(row)
        shard_row["_summary_order"] = str(manifest_order)
        shards[shard_id].append(shard_row)
        new_load = current_load + max(0, safe_int(row.get("common_marker_count")))
        heapq.heappush(loads, (new_load, shard_id))

    for shard in shards:
        shard.sort(key=lambda row: safe_int(row.get("_summary_order")))
    return shards


def shard_paths(parts_dir: Path, shard_id: int) -> Dict[str, Path]:
    prefix = parts_dir / "part-{0:02d}".format(shard_id)
    return {
        "block": Path(str(prefix) + ".block.tsv"),
        "top": Path(str(prefix) + ".top.tsv.gz"),
        "credible": Path(str(prefix) + ".credible.tsv.gz"),
        "variant": Path(str(prefix) + ".variant.tsv.gz"),
    }


def process_shard(
    shard_id: int,
    rows: Sequence[Dict[str, str]],
    success_rows: Mapping[str, Mapping[str, str]],
    failed_rows: Mapping[str, Mapping[str, str]],
    skipped_rows: Mapping[str, Mapping[str, str]],
    parts_dir: Path,
) -> Dict[str, object]:

    paths = shard_paths(parts_dir, shard_id)
    block_handle, block_writer = open_tsv_writer(
        paths["block"], BLOCK_FIELDS, write_header=False
    )
    top_handle, top_writer = open_tsv_writer(
        paths["top"], TOP_FIELDS, compressed=True, write_header=False
    )
    credible_handle, credible_writer = open_tsv_writer(
        paths["credible"], VARIANT_FIELDS, compressed=True, write_header=False
    )
    variant_handle, variant_writer = open_tsv_writer(
        paths["variant"], VARIANT_FIELDS, compressed=True, write_header=False
    )

    valid_count = 0
    invalid_count = 0
    variant_count = 0
    credible_variant_count = 0
    try:
        for row in rows:
            block_id = clean(row.get("block_id"))
            if block_id in failed_rows:
                failed = failed_rows[block_id]
                stderr_last_line = last_nonempty_line(
                    PLASMA_RESULT_ROOT / block_id / "run.stderr.log"
                )
                non_positive_definite = (
                    stderr_last_line
                    == "numpy.linalg.LinAlgError: Matrix is not positive definite"
                )
                block = make_unprocessed_block(
                    row,
                    "failed",
                    clean(failed.get("exit_code")),
                    clean(failed.get("reason")),
                    (
                        "not_estimable_non_positive_definite"
                        if non_positive_definite
                        else "runner_execution_failure"
                    ),
                    (
                        "joint_model_non_positive_definite"
                        if non_positive_definite
                        else "runner_failed"
                    ),
                )
            elif block_id in skipped_rows:
                skipped = skipped_rows[block_id]
                block = make_unprocessed_block(
                    row,
                    "skipped",
                    clean(skipped.get("exit_code")),
                    clean(skipped.get("reason")),
                    "not_estimable_constant_expression",
                    "joint_not_estimable_constant_expression",
                )
            elif block_id not in success_rows:
                block = make_unprocessed_block(
                    row,
                    "missing",
                    "",
                    "",
                    "missing_terminal_status",
                    "block_absent_from_8.3_status_logs",
                )
            else:
                block, block_variant_count, block_credible_count = process_completed_block(
                    row,
                    success_rows[block_id],
                    variant_writer,
                    credible_writer,
                    top_writer,
                )
                variant_count += block_variant_count
                credible_variant_count += block_credible_count

            block_writer.writerow({field: block.get(field, "") for field in BLOCK_FIELDS})
            if clean(block.get("result_status")) == "valid":
                valid_count += 1
            else:
                invalid_count += 1
    finally:
        block_handle.close()
        top_handle.close()
        credible_handle.close()
        variant_handle.close()

    return {
        "shard_id": shard_id,
        "block_count": len(rows),
        "valid_count": valid_count,
        "invalid_count": invalid_count,
        "variant_count": variant_count,
        "credible_variant_count": credible_variant_count,
        "paths": {key: str(value) for key, value in paths.items()},
    }


def merge_text_parts(
    destination: Path,
    fieldnames: Sequence[str],
    part_paths: Sequence[Path],
) -> None:
    with destination.open("w", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            delimiter="\t",
            fieldnames=list(fieldnames),
            lineterminator="\n",
        )
        writer.writeheader()
        for part_path in part_paths:
            with part_path.open() as part_handle:
                shutil.copyfileobj(part_handle, output_handle)


def merge_gzip_parts(
    destination: Path,
    fieldnames: Sequence[str],
    part_paths: Sequence[Path],
) -> None:

    with gzip.open(destination, "wt", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            delimiter="\t",
            fieldnames=list(fieldnames),
            lineterminator="\n",
        )
        writer.writeheader()
    with destination.open("ab") as output_binary:
        for part_path in part_paths:
            with part_path.open("rb") as part_binary:
                shutil.copyfileobj(part_binary, output_binary)


def check_input_contract() -> None:
    required_files = [
        MANIFEST,
        PREPARE_COMPLETION_MARKER,
        SUCCESS_LOG,
        FAILED_LOG,
        SKIPPED_LOG,
        RUN_SUMMARY,
    ]
    missing_files = [str(path) for path in required_files if not path.is_file()]
    if missing_files:
        raise FileNotFoundError(
            "missing required 8.4 input files: {0}".format(",".join(missing_files))
        )
    required_directories = [PLASMA_INPUT_ROOT, PLASMA_RESULT_ROOT]
    missing_directories = [
        str(path) for path in required_directories if not path.is_dir()
    ]
    if missing_directories:
        raise FileNotFoundError(
            "missing required 8.4 input directories: {0}".format(
                ",".join(missing_directories)
            )
        )


def check_output_contract() -> None:
    if not RUN_ID or Path(RUN_ID).name != RUN_ID or RUN_ID in {".", ".."}:
        raise ValueError("RUN_ID must be one directory name: {0}".format(RUN_ID))
    if SUMMARY_DIR.exists():
        raise FileExistsError(
            "refusing to reuse existing run directory: {0}".format(SUMMARY_DIR)
        )
    outputs = [
        BLOCK_SUMMARY,
        VALID_BLOCKS,
        INVALID_BLOCKS,
        TOP_VARIANTS,
        CREDIBLE_SET_VARIANTS,
        VARIANT_SUMMARY,
        SUMMARY_QC,
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError(
            "refusing to overwrite existing 8.4 outputs: {0}".format(",".join(existing))
        )


def main() -> None:
    log("[start] summary five dogs no-GWAS PLASMA expression + ASE joint result")
    log("[config] PREPARE_RUN_ID={0}".format(PREPARE_RUN_ID))
    log("[config] PLASMA_RUN_ID={0}".format(PLASMA_RUN_ID))
    log("[config] MANIFEST={0}".format(MANIFEST))
    log("[config] SUCCESS_LOG={0}".format(SUCCESS_LOG))
    log("[config] FAILED_LOG={0}".format(FAILED_LOG))
    log("[config] SUMMARY_DIR={0}".format(SUMMARY_DIR))
    log(
        "[config] parallel_jobs={0}; threads_per_task={1}; resource_source=preserved_existing".format(
            PARALLEL_JOBS, THREADS_PER_TASK
        )
    )

    check_input_contract()
    ready_rows = load_ready_manifest(MANIFEST)
    status_required = [
        "run_id",
        "input_dir",
        "output_path",
        "selected_sample_count",
        "marker_count",
        "total_exp_unique_count",
        "analysis_mode",
        "status",
        "exit_code",
        "reason",
    ]
    success_rows = unique_rows_by_key(
        SUCCESS_LOG,
        status_required,
    )
    failed_rows = unique_rows_by_key(
        FAILED_LOG,
        status_required,
    )
    skipped_rows = unique_rows_by_key(
        SKIPPED_LOG,
        status_required,
    )
    run_summary_fields, run_summary_rows = read_tsv(RUN_SUMMARY)
    required_columns(RUN_SUMMARY, run_summary_fields, ["key", "value"])
    run_summary_values = {clean(row.get("key")): clean(row.get("value")) for row in run_summary_rows}

    validate_status_logs(ready_rows, success_rows, failed_rows, skipped_rows)
    expected_log_counts = {
        "tasks": len(ready_rows),
        "success": len(success_rows),
        "failed": len(failed_rows),
        "skipped": len(skipped_rows),
        "terminal_total": len(success_rows) + len(failed_rows) + len(skipped_rows),
    }
    for key, observed_count in expected_log_counts.items():
        summary_count = safe_int(run_summary_values.get(key), -1)
        if summary_count != observed_count:
            raise ValueError(
                "{0} count mismatch: run_summary={1}; observed={2}".format(
                    key, summary_count, observed_count
                )
            )
    if run_summary_values.get("run_id") != PLASMA_RUN_ID:
        raise ValueError("8.3 run_summary run_id does not match PLASMA_RUN_ID")
    if run_summary_values.get("prepare_manifest") != str(MANIFEST):
        raise ValueError("8.3 run_summary prepare_manifest does not match 8.4 MANIFEST")
    if run_summary_values.get("result_dir") != str(PLASMA_RESULT_ROOT):
        raise ValueError("8.3 run_summary result_dir does not match PLASMA_RESULT_ROOT")
    if safe_int(run_summary_values.get("max_causal"), -1) != 1:
        raise ValueError("8.3 max_causal is not 1")
    try:
        observed_confidence = float(run_summary_values.get("confidence", ""))
    except ValueError as exc:
        raise ValueError("8.3 run_summary confidence is not numeric") from exc
    if abs(observed_confidence - CREDIBLE_SET_CONFIDENCE) > 1e-12:
        raise ValueError(
            "8.3 confidence mismatch: {0}!={1}".format(
                observed_confidence, CREDIBLE_SET_CONFIDENCE
            )
        )
    log(
        "[stage:validate] ready={0}; success={1}; failed={2}; skipped={3}".format(
            len(ready_rows), len(success_rows), len(failed_rows), len(skipped_rows)
        )
    )

    check_output_contract()
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    PARTS_DIR.mkdir()
    shards = partition_ready_rows(ready_rows)
    shard_loads = [
        sum(max(0, safe_int(row.get("common_marker_count"))) for row in shard)
        for shard in shards
    ]
    log(
        "[stage:partition] workers={0}; shard_blocks={1}-{2}; shard_markers={3}-{4}".format(
            PARALLEL_JOBS,
            min(len(shard) for shard in shards),
            max(len(shard) for shard in shards),
            min(shard_loads),
            max(shard_loads),
        )
    )

    log("[stage:parallel] submit {0} items independent shard".format(PARALLEL_JOBS))
    shard_results: Dict[int, Dict[str, object]] = {}
    failed_shard_count = 0
    with ProcessPoolExecutor(max_workers=PARALLEL_JOBS) as executor:
        futures: Dict[object, Tuple[int, str]] = {}
        for shard_id, shard in enumerate(shards):
            shard_ids = {clean(row.get("block_id")) for row in shard}
            shard_success = {
                block_id: row for block_id, row in success_rows.items() if block_id in shard_ids
            }
            shard_failed = {
                block_id: row for block_id, row in failed_rows.items() if block_id in shard_ids
            }
            shard_skipped = {
                block_id: row for block_id, row in skipped_rows.items() if block_id in shard_ids
            }
            future = executor.submit(
                process_shard,
                shard_id,
                shard,
                shard_success,
                shard_failed,
                shard_skipped,
                PARTS_DIR,
            )
            futures[future] = (shard_id, current_time())

        for completed_count, future in enumerate(as_completed(futures), 1):
            submitted_shard_id, started_at = futures[future]
            finished_at = current_time()
            try:
                result = future.result()
            except Exception as exc:
                failed_shard_count += 1
                reason = "{0}: {1}".format(type(exc).__name__, clean(exc))
                reason = reason.replace("\t", " ").replace("\n", " ")
                append_shard_status(
                    SHARD_FAILED_LOG,
                    {
                        "run_id": RUN_ID,
                        "shard_id": submitted_shard_id,
                        "status": "failed",
                        "exit_code": 1,
                        "started_at": started_at,
                        "finished_at": finished_at,
                        "output_path": str(PARTS_DIR),
                        "task_log": str(MASTER_LOG),
                        "reason": reason,
                    },
                )
                log(
                    "[stage:parallel] completed_shards={0}/{1}; shard={2}; status=failed; reason={3}".format(
                        completed_count,
                        PARALLEL_JOBS,
                        submitted_shard_id,
                        reason,
                    )
                )
                continue
            shard_id = safe_int(result.get("shard_id"), -1)
            shard_results[shard_id] = result
            append_shard_status(
                SHARD_SUCCESS_LOG,
                {
                    "run_id": RUN_ID,
                    "shard_id": shard_id,
                    "status": "success",
                    "exit_code": 0,
                    "started_at": started_at,
                    "finished_at": finished_at,
                    "output_path": str(PARTS_DIR),
                    "task_log": str(MASTER_LOG),
                    "reason": "",
                },
            )
            log(
                "[stage:parallel] completed_shards={0}/{1}; shard={2}; status=success; blocks={3}".format(
                    completed_count,
                    PARALLEL_JOBS,
                    shard_id,
                    result.get("block_count", ""),
                )
            )

    if failed_shard_count:
        raise RuntimeError("{0} 8.4 shards failed".format(failed_shard_count))
    if sorted(shard_results) != list(range(PARALLEL_JOBS)):
        raise RuntimeError(
            "not all {0} shards returned a result".format(PARALLEL_JOBS)
        )
    if sum(safe_int(result.get("block_count")) for result in shard_results.values()) != len(ready_rows):
        raise RuntimeError("merged shard block count does not match ready manifest")

    ordered_results = [shard_results[shard_id] for shard_id in range(PARALLEL_JOBS)]
    block_parts = [Path(result["paths"]["block"]) for result in ordered_results]
    top_parts = [Path(result["paths"]["top"]) for result in ordered_results]
    credible_parts = [Path(result["paths"]["credible"]) for result in ordered_results]
    variant_parts = [Path(result["paths"]["variant"]) for result in ordered_results]

    log("[stage:merge] merge {0} items shard result".format(PARALLEL_JOBS))
    merge_text_parts(BLOCK_SUMMARY, BLOCK_FIELDS, block_parts)
    merge_gzip_parts(TOP_VARIANTS, TOP_FIELDS, top_parts)
    merge_gzip_parts(CREDIBLE_SET_VARIANTS, VARIANT_FIELDS, credible_parts)
    merge_gzip_parts(VARIANT_SUMMARY, VARIANT_FIELDS, variant_parts)

    block_fields, block_rows = read_tsv(BLOCK_SUMMARY)
    required_columns(BLOCK_SUMMARY, block_fields, BLOCK_FIELDS)
    if len(block_rows) != len(ready_rows):
        raise RuntimeError(
            "block summary row count mismatch: {0}!={1}".format(
                len(block_rows), len(ready_rows)
            )
        )
    valid_rows = [
        row for row in block_rows if clean(row.get("result_status")) == "valid"
    ]
    invalid_rows = [
        row for row in block_rows if clean(row.get("result_status")) != "valid"
    ]
    write_tsv(VALID_BLOCKS, BLOCK_FIELDS, valid_rows)
    write_tsv(INVALID_BLOCKS, BLOCK_FIELDS, invalid_rows)

    variant_count = sum(
        safe_int(result.get("variant_count")) for result in ordered_results
    )
    credible_variant_count = sum(
        safe_int(result.get("credible_variant_count")) for result in ordered_results
    )
    ppa_boundary_adjusted_blocks = sum(
        safe_int(row.get("ppa_boundary_adjusted_count")) > 0 for row in block_rows
    )
    ppa_boundary_adjusted_values = sum(
        safe_int(row.get("ppa_boundary_adjusted_count")) for row in block_rows
    )
    non_positive_definite_blocks = sum(
        clean(row.get("qc_category")) == "not_estimable_non_positive_definite"
        for row in block_rows
    )
    constant_expression_blocks = sum(
        clean(row.get("qc_category")) == "not_estimable_constant_expression"
        for row in block_rows
    )
    qc_rows = [
        {"key": "run_id", "value": RUN_ID},
        {"key": "prepare_run_id", "value": PREPARE_RUN_ID},
        {"key": "plasma_run_id", "value": PLASMA_RUN_ID},
        {"key": "prepare_manifest", "value": MANIFEST},
        {"key": "plasma_result_root", "value": PLASMA_RESULT_ROOT},
        {"key": "runner_success_log", "value": SUCCESS_LOG},
        {"key": "runner_failed_log", "value": FAILED_LOG},
        {"key": "runner_skipped_log", "value": SKIPPED_LOG},
        {"key": "parallel_jobs", "value": PARALLEL_JOBS},
        {"key": "threads_per_task", "value": THREADS_PER_TASK},
        {"key": "resource_source", "value": "preserved_existing"},
        {"key": "manifest_ready_blocks", "value": len(ready_rows)},
        {"key": "runner_success_blocks", "value": len(success_rows)},
        {"key": "runner_failed_blocks", "value": len(failed_rows)},
        {"key": "runner_skipped_blocks", "value": len(skipped_rows)},
        {
            "key": "not_estimable_non_positive_definite_blocks",
            "value": non_positive_definite_blocks,
        },
        {
            "key": "not_estimable_constant_expression_blocks",
            "value": constant_expression_blocks,
        },
        {"key": "valid_blocks", "value": len(valid_rows)},
        {"key": "invalid_blocks", "value": len(invalid_rows)},
        {"key": "variant_rows_written", "value": variant_count},
        {"key": "credible_set_variant_rows_written", "value": credible_variant_count},
        {"key": "ppa_boundary_tolerance", "value": PPA_BOUND_TOLERANCE},
        {"key": "ppa_boundary_adjusted_blocks", "value": ppa_boundary_adjusted_blocks},
        {"key": "ppa_boundary_adjusted_values", "value": ppa_boundary_adjusted_values},
        {"key": "ppa_sum_tolerance", "value": PPA_SUM_TOLERANCE},
        {"key": "credible_set_confidence", "value": CREDIBLE_SET_CONFIDENCE},
        {
            "key": "credible_set_sum_tolerance",
            "value": CREDIBLE_SET_SUM_TOLERANCE,
        },
    ]
    for key in [
        "search_mode",
        "max_causal",
        "confidence",
        "overdispersion_global",
        "plasma_python",
        "plasma_script",
    ]:
        qc_rows.append({"key": "8.3_{0}".format(key), "value": run_summary_values.get(key, "")})
    write_tsv(SUMMARY_QC, ["key", "value"], qc_rows)

    log(
        "[summary] block_result_status=valid={0};invalid={1}".format(
            len(valid_rows), len(invalid_rows)
        )
    )
    log("[summary] variant_rows={0}; credible_set_variant_rows={1}".format(
        variant_count, credible_variant_count
    ))
    log("[output] {0}".format(BLOCK_SUMMARY))
    log("[output] {0}".format(TOP_VARIANTS))
    log("[output] {0}".format(VARIANT_SUMMARY))
    log("[output] {0}".format(CREDIBLE_SET_VARIANTS))
    COMPLETION_MARKER.write_text(
        "run_id={0}\nexecution_status=success\nfinished_at={1}\n".format(
            RUN_ID, current_time()
        )
    )
    log("[done] five dogs 8.4 PLASMA summary completed; task completed / done")
    log("execution_status=success")
    log(
        "[reminder] execution completion does not establish result completeness, format, counts or biological quality have passed review; please separately check output."
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
