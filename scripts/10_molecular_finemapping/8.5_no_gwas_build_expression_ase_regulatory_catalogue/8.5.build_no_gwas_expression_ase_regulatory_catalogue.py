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


SOURCE_RUN_ID = os.environ.get("SOURCE_RUN_ID", "20260827_143014")
SUMMARY_ROOT = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.4_no_gwas_summarize_expression_ase_plasma_joint'
)
SUMMARY_DIR = SUMMARY_ROOT / SOURCE_RUN_ID
VARIANT_SUMMARY = SUMMARY_DIR / "no_gwas_plasma_joint_variant_ppa_summary.tsv.gz"
VALID_BLOCKS = SUMMARY_DIR / "no_gwas_plasma_joint_block_valid.tsv"
INVALID_BLOCKS = SUMMARY_DIR / "no_gwas_plasma_joint_block_invalid.tsv"
SOURCE_QC = SUMMARY_DIR / "no_gwas_plasma_joint_run_qc.tsv"
SOURCE_COMPLETION_MARKER = SUMMARY_DIR / "execution_completion.marker"

OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.5_no_gwas_build_expression_ase_regulatory_catalogue'
)
RUN_ID = os.environ.get("RUN_ID", datetime.now().strftime("%Y%m%d_%H%M%S"))
OUTPUT_DIR = OUTPUT_ROOT / RUN_ID
COMPLETION_MARKER = OUTPUT_DIR / "execution_completion.marker"

SCRIPT_DIR = Path(
    '/path/to/repository/scripts/10_molecular_finemapping/8.5_no_gwas_build_expression_ase_regulatory_catalogue'
)
LOG_DIR = SCRIPT_DIR / "log"
MASTER_LOG = LOG_DIR / (
    "8.5.build_no_gwas_expression_ase_regulatory_catalogue.{0}.log".format(RUN_ID)
)

CONCURRENT_TASKS = 1
THREADS_PER_TASK = 1
EXPRESSION_MATRIX_SCOPE = "reference_annotation_plus_dog5_manual_new_gene"
MASTER_LOG_HANDLE: Optional[TextIO] = None

VARIANT_CATALOGUE = (
    OUTPUT_DIR / "no_gwas_expression_ase_regulatory_variant_catalogue.tsv.gz"
)
TOP_VARIANT_CATALOGUE = (
    OUTPUT_DIR / "no_gwas_expression_ase_top_variant_catalogue.tsv.gz"
)
CREDIBLE_SET_CATALOGUE = (
    OUTPUT_DIR / "no_gwas_expression_ase_credible_set_catalogue.tsv.gz"
)
VALID_BLOCK_CATALOGUE = (
    OUTPUT_DIR / "no_gwas_expression_ase_valid_block_catalogue.tsv"
)
INVALID_BLOCK_QC = OUTPUT_DIR / "no_gwas_expression_ase_invalid_block_qc.tsv"
TISSUE_SUMMARY = (
    OUTPUT_DIR / "no_gwas_expression_ase_catalogue_tissue_summary.tsv"
)
CATALOGUE_SUMMARY = (
    OUTPUT_DIR / "no_gwas_expression_ase_catalogue_summary.tsv"
)

PROGRESS_EVERY = 1_000_000

VARIANT_REQUIRED = [
    "block_id",
    "tissue",
    "feature_id",
    "chrom",
    "pos",
    "marker_id",
    "ref",
    "alt",
    "ppa",
    "marker_rank_by_ppa",
    "credible_set_member",
    "result_status",
    "result_interpretable",
]
BLOCK_REQUIRED = [
    "block_id",
    "tissue",
    "feature_id",
    "chrom",
    "observed_marker_count",
    "result_status",
    "result_interpretable",
    "ppa_status",
    "ppa_boundary_adjusted_count",
    "max_ppa",
    "top_ppa_tie_count",
    "top_marker_is_unique",
    "credible_set_size",
    "credible_set_ppa_sum",
    "qc_reason",
]
VARIANT_DERIVED_FIELDS = [
    "catalogue_layer",
    "catalogue_role",
    "feature_origin",
    "expression_matrix_scope",
    "is_top_variant",
    "is_credible_set_member",
    "ppa_tier",
    "ppa_status",
    "ppa_boundary_adjusted_count",
    "block_max_ppa",
    "block_top_ppa_tie_count",
    "block_top_marker_is_unique",
    "block_credible_set_size",
    "block_credible_set_ppa_sum",
    "interpretation_note",
]
BLOCK_DERIVED_FIELDS = [
    "catalogue_layer",
    "feature_origin",
    "expression_matrix_scope",
    "interpretation_note",
]
TISSUE_FIELDS = [
    "tissue",
    "valid_blocks",
    "invalid_blocks",
    "interpretable_variant_rows",
    "top_variant_rows",
    "credible_set_variant_rows",
    "ppa_ge_0.95_variant_rows",
    "ppa_ge_0.8_variant_rows",
    "ppa_ge_0.5_variant_rows",
    "ppa_ge_0.1_variant_rows",
]


def log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S ")
    print("[{0}] {1}".format(timestamp, message), flush=True)


def initialize_execution_logging() -> None:
    global MASTER_LOG_HANDLE
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    MASTER_LOG_HANDLE = MASTER_LOG.open("x", buffering=1)
    os.dup2(MASTER_LOG_HANDLE.fileno(), sys.stdout.fileno())
    os.dup2(MASTER_LOG_HANDLE.fileno(), sys.stderr.fileno())


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(clean(value))
    except (TypeError, ValueError):
        return default


def safe_float(value: object) -> float:
    try:
        return float(clean(value))
    except (TypeError, ValueError):
        return float("nan")


def required_columns(
    path: Path,
    observed: Sequence[str],
    required: Sequence[str],
) -> None:
    missing = [field for field in required if field not in observed]
    if missing:
        raise ValueError(
            "{0} missing required columns: {1}".format(path, ",".join(missing))
        )


def merge_fields(preferred: Sequence[str], observed: Sequence[str]) -> List[str]:
    merged: List[str] = []
    for field in list(preferred) + list(observed):
        if field not in merged:
            merged.append(field)
    return merged


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError("missing or empty TSV: {0}".format(path))
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(row) for row in reader]


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
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def open_gzip_writer(
    path: Path,
    fieldnames: Sequence[str],
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


def unique_by_block_id(
    path: Path,
    fields: Sequence[str],
    rows: Sequence[Dict[str, str]],
) -> Dict[str, Dict[str, str]]:
    required_columns(path, fields, BLOCK_REQUIRED)
    indexed: Dict[str, Dict[str, str]] = {}
    for row in rows:
        block_id = clean(row.get("block_id"))
        if not block_id:
            raise ValueError("{0} contains an empty block_id".format(path))
        if block_id in indexed:
            raise ValueError(
                "{0} contains duplicate block_id: {1}".format(path, block_id)
            )
        indexed[block_id] = row
    return indexed


def read_source_qc(path: Path) -> Dict[str, str]:
    fields, rows = read_tsv(path)
    required_columns(path, fields, ["key", "value"])
    values: Dict[str, str] = {}
    for row in rows:
        key = clean(row.get("key"))
        if not key:
            raise ValueError("{0} contains an empty QC key".format(path))
        if key in values:
            raise ValueError("{0} contains duplicate QC key: {1}".format(path, key))
        values[key] = clean(row.get("value"))
    return values


def read_completion_marker(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = raw_line.partition("=")
        if not separator or not key.strip():
            raise ValueError("invalid completion marker line: {0}".format(raw_line))
        key = key.strip()
        if key in values:
            raise ValueError("duplicate completion marker key: {0}".format(key))
        values[key] = value.strip()
    return values


def check_paths() -> None:
    if (
        not SOURCE_RUN_ID
        or Path(SOURCE_RUN_ID).name != SOURCE_RUN_ID
        or SOURCE_RUN_ID in {".", ".."}
    ):
        raise ValueError(
            "SOURCE_RUN_ID must be one directory name: {0}".format(SOURCE_RUN_ID)
        )
    if not RUN_ID or Path(RUN_ID).name != RUN_ID or RUN_ID in {".", ".."}:
        raise ValueError("RUN_ID must be one directory name: {0}".format(RUN_ID))
    for path in [
        VARIANT_SUMMARY,
        VALID_BLOCKS,
        INVALID_BLOCKS,
        SOURCE_QC,
        SOURCE_COMPLETION_MARKER,
    ]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("missing or empty 8.4 input: {0}".format(path))
    if OUTPUT_DIR.exists():
        raise FileExistsError(
            "refusing to reuse existing 8.5 run directory: {0}".format(OUTPUT_DIR)
        )


def validate_block_partition(
    valid_blocks: Mapping[str, Mapping[str, str]],
    invalid_blocks: Mapping[str, Mapping[str, str]],
    source_qc: Mapping[str, str],
) -> None:
    overlap = set(valid_blocks) & set(invalid_blocks)
    if overlap:
        raise ValueError(
            "valid and invalid block tables overlap: {0}".format(sorted(overlap)[0])
        )
    if any(
        clean(row.get("result_status")) != "valid"
        or clean(row.get("result_interpretable")) != "1"
        for row in valid_blocks.values()
    ):
        raise ValueError("valid block table contains a non-interpretable block")
    if any(
        clean(row.get("result_status")) == "valid"
        or clean(row.get("result_interpretable")) == "1"
        for row in invalid_blocks.values()
    ):
        raise ValueError("invalid block table contains an interpretable block")

    expected_valid = safe_int(source_qc.get("valid_blocks"), -1)
    expected_invalid = safe_int(source_qc.get("invalid_blocks"), -1)
    expected_total = safe_int(source_qc.get("manifest_ready_blocks"), -1)
    if len(valid_blocks) != expected_valid:
        raise ValueError(
            "valid block count mismatch: table={0}; QC={1}".format(
                len(valid_blocks), expected_valid
            )
        )
    if len(invalid_blocks) != expected_invalid:
        raise ValueError(
            "invalid block count mismatch: table={0}; QC={1}".format(
                len(invalid_blocks), expected_invalid
            )
        )
    if len(valid_blocks) + len(invalid_blocks) != expected_total:
        raise ValueError(
            "block partition does not match manifest_ready_blocks: {0}!={1}".format(
                len(valid_blocks) + len(invalid_blocks), expected_total
            )
        )


def ppa_tier(ppa: float) -> str:
    if ppa >= 0.95:
        return "ppa_ge_0.95"
    if ppa >= 0.8:
        return "ppa_0.8_to_lt_0.95"
    if ppa >= 0.5:
        return "ppa_0.5_to_lt_0.8"
    if ppa >= 0.1:
        return "ppa_0.1_to_lt_0.5"
    return "ppa_lt_0.1"


def catalogue_role(is_top: bool, is_credible: bool) -> str:
    if is_top and is_credible:
        return "top_and_credible_set"
    if is_top:
        return "top_variant"
    if is_credible:
        return "credible_set"
    return "finite_background"


def feature_origin(feature_id: str) -> str:

    if feature_id.startswith("DOG5_NEWGENE_"):
        return "dog5_manual_new_gene"
    if feature_id.startswith("gene-"):
        return "reference_annotation_gene"
    raise ValueError(
        "unsupported feature_id prefix for feature_origin: {0}".format(feature_id)
    )


def add_feature_metadata(row: Mapping[str, str]) -> Dict[str, str]:
    annotated = dict(row)
    annotated.update(
        {
            "feature_origin": feature_origin(clean(row.get("feature_id"))),
            "expression_matrix_scope": EXPRESSION_MATRIX_SCOPE,
        }
    )
    return annotated


def init_tissue_bucket() -> Dict[str, int]:
    return {
        "valid_blocks": 0,
        "invalid_blocks": 0,
        "interpretable_variant_rows": 0,
        "top_variant_rows": 0,
        "credible_set_variant_rows": 0,
        "ppa_ge_0.95_variant_rows": 0,
        "ppa_ge_0.8_variant_rows": 0,
        "ppa_ge_0.5_variant_rows": 0,
        "ppa_ge_0.1_variant_rows": 0,
    }


def annotate_block(
    row: Mapping[str, str],
    layer: str,
    note: str,
) -> Dict[str, object]:
    annotated: Dict[str, object] = dict(row)
    annotated.update(
        {
            "catalogue_layer": layer,
            "interpretation_note": note,
        }
    )
    return annotated


def update_threshold_counts(
    bucket: Dict[str, int],
    ppa: float,
) -> None:
    if ppa >= 0.95:
        bucket["ppa_ge_0.95_variant_rows"] += 1
    if ppa >= 0.8:
        bucket["ppa_ge_0.8_variant_rows"] += 1
    if ppa >= 0.5:
        bucket["ppa_ge_0.5_variant_rows"] += 1
    if ppa >= 0.1:
        bucket["ppa_ge_0.1_variant_rows"] += 1


def main() -> None:
    log("[start] build five dogs no-GWAS expression + ASE regulatory catalogue")
    log("[config] SOURCE_RUN_ID={0}".format(SOURCE_RUN_ID))
    log("[config] VARIANT_SUMMARY={0}".format(VARIANT_SUMMARY))
    log("[config] OUTPUT_DIR={0}".format(OUTPUT_DIR))
    log(
        "[config] concurrent_tasks={0}; threads_per_task={1}; resource_source=preserved_serial".format(
            CONCURRENT_TASKS, THREADS_PER_TASK
        )
    )
    log("[config] expression_matrix_scope={0}".format(EXPRESSION_MATRIX_SCOPE))

    check_paths()
    valid_fields, valid_rows = read_tsv(VALID_BLOCKS)
    invalid_fields, invalid_rows = read_tsv(INVALID_BLOCKS)
    valid_rows = [add_feature_metadata(row) for row in valid_rows]
    invalid_rows = [add_feature_metadata(row) for row in invalid_rows]
    valid_by_id = unique_by_block_id(VALID_BLOCKS, valid_fields, valid_rows)
    invalid_by_id = unique_by_block_id(INVALID_BLOCKS, invalid_fields, invalid_rows)
    source_qc = read_source_qc(SOURCE_QC)
    source_completion = read_completion_marker(SOURCE_COMPLETION_MARKER)
    validate_block_partition(valid_by_id, invalid_by_id, source_qc)
    if clean(source_qc.get("run_id")) != SOURCE_RUN_ID:
        raise ValueError(
            "8.4 QC run_id does not match SOURCE_RUN_ID: {0}!={1}".format(
                source_qc.get("run_id", ""), SOURCE_RUN_ID
            )
        )
    if clean(source_completion.get("run_id")) != SOURCE_RUN_ID:
        raise ValueError(
            "8.4 completion run_id does not match SOURCE_RUN_ID: {0}!={1}".format(
                source_completion.get("run_id", ""), SOURCE_RUN_ID
            )
        )
    if clean(source_completion.get("execution_status")) != "success":
        raise ValueError("8.4 completion marker is not successful")
    log(
        "[stage:validate] valid_blocks={0}; invalid_blocks={1}".format(
            len(valid_by_id), len(invalid_by_id)
        )
    )

    with gzip.open(VARIANT_SUMMARY, "rt", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        variant_fields = list(reader.fieldnames or [])
        required_columns(VARIANT_SUMMARY, variant_fields, VARIANT_REQUIRED)

        catalogue_fields = merge_fields(VARIANT_DERIVED_FIELDS, variant_fields)
        OUTPUT_DIR.mkdir(parents=True)
        variant_handle, variant_writer = open_gzip_writer(
            VARIANT_CATALOGUE, catalogue_fields
        )
        top_handle, top_writer = open_gzip_writer(
            TOP_VARIANT_CATALOGUE, catalogue_fields
        )
        credible_handle, credible_writer = open_gzip_writer(
            CREDIBLE_SET_CATALOGUE, catalogue_fields
        )

        source_variant_rows = 0
        interpretable_variant_rows = 0
        skipped_invalid_variant_rows = 0
        top_variant_rows = 0
        credible_set_variant_rows = 0
        threshold_counts = init_tissue_bucket()
        tissue_counts: Dict[str, Dict[str, int]] = {}
        for row in valid_rows:
            tissue = clean(row.get("tissue")) or "NA"
            tissue_counts.setdefault(tissue, init_tissue_bucket())[
                "valid_blocks"
            ] += 1
        for row in invalid_rows:
            tissue = clean(row.get("tissue")) or "NA"
            tissue_counts.setdefault(tissue, init_tissue_bucket())[
                "invalid_blocks"
            ] += 1

        try:
            for row in reader:
                source_variant_rows += 1
                block_id = clean(row.get("block_id"))
                if block_id in invalid_by_id:
                    if (
                        clean(row.get("result_status")) == "valid"
                        or clean(row.get("result_interpretable")) == "1"
                    ):
                        raise ValueError(
                            "invalid block has an interpretable variant row: {0}".format(
                                block_id
                            )
                        )
                    skipped_invalid_variant_rows += 1
                    continue
                if block_id not in valid_by_id:
                    raise ValueError(
                        "variant row references an unknown block_id: {0}".format(
                            block_id
                        )
                    )
                if (
                    clean(row.get("result_status")) != "valid"
                    or clean(row.get("result_interpretable")) != "1"
                ):
                    raise ValueError(
                        "valid block has a non-interpretable variant row: {0}".format(
                            block_id
                        )
                    )

                block = valid_by_id[block_id]
                if clean(row.get("chrom")) != clean(block.get("chrom")):
                    raise ValueError(
                        "variant/block chromosome mismatch for {0}".format(block_id)
                    )
                ppa = safe_float(row.get("ppa"))
                rank = safe_int(row.get("marker_rank_by_ppa"), -1)
                credible_token = clean(row.get("credible_set_member"))
                if not math.isfinite(ppa) or not 0.0 <= ppa <= 1.0:
                    raise ValueError(
                        "valid block contains an invalid PPA: {0}".format(block_id)
                    )
                if rank < 1:
                    raise ValueError(
                        "valid block contains an invalid marker rank: {0}".format(
                            block_id
                        )
                    )
                if credible_token not in {"0", "1"}:
                    raise ValueError(
                        "valid block contains a non-binary credible_set_member: {0}".format(
                            block_id
                        )
                    )

                is_top = rank == 1
                is_credible = credible_token == "1"
                annotated: Dict[str, object] = dict(row)
                annotated.update(
                    {
                        "catalogue_layer": (
                            "no_gwas_expression_ase_regulatory_variant"
                        ),
                        "catalogue_role": catalogue_role(is_top, is_credible),
                        "feature_origin": clean(block.get("feature_origin")),
                        "expression_matrix_scope": clean(
                            block.get("expression_matrix_scope")
                        ),
                        "is_top_variant": int(is_top),
                        "is_credible_set_member": int(is_credible),
                        "ppa_tier": ppa_tier(ppa),
                        "ppa_status": clean(block.get("ppa_status")),
                        "ppa_boundary_adjusted_count": clean(
                            block.get("ppa_boundary_adjusted_count")
                        ),
                        "block_max_ppa": clean(block.get("max_ppa")),
                        "block_top_ppa_tie_count": clean(
                            block.get("top_ppa_tie_count")
                        ),
                        "block_top_marker_is_unique": clean(
                            block.get("top_marker_is_unique")
                        ),
                        "block_credible_set_size": clean(
                            block.get("credible_set_size")
                        ),
                        "block_credible_set_ppa_sum": clean(
                            block.get("credible_set_ppa_sum")
                        ),
                        "interpretation_note": (
                            "exploratory no-GWAS expression+ASE joint fine-mapping "
                            "candidate; association/FDR screening and a null-causal "
                            "model were not applied"
                        ),
                    }
                )
                output_row = {
                    field: annotated.get(field, "") for field in catalogue_fields
                }
                variant_writer.writerow(output_row)
                if is_top:
                    top_writer.writerow(output_row)
                    top_variant_rows += 1
                if is_credible:
                    credible_writer.writerow(output_row)
                    credible_set_variant_rows += 1

                interpretable_variant_rows += 1
                tissue = clean(row.get("tissue")) or "NA"
                tissue_bucket = tissue_counts.setdefault(
                    tissue, init_tissue_bucket()
                )
                tissue_bucket["interpretable_variant_rows"] += 1
                if is_top:
                    tissue_bucket["top_variant_rows"] += 1
                if is_credible:
                    tissue_bucket["credible_set_variant_rows"] += 1
                update_threshold_counts(tissue_bucket, ppa)
                update_threshold_counts(threshold_counts, ppa)

                if source_variant_rows % PROGRESS_EVERY == 0:
                    log(
                        "[stage:stream] source_rows={0}; catalogue_rows={1}; "
                        "top_rows={2}; credible_rows={3}".format(
                            source_variant_rows,
                            interpretable_variant_rows,
                            top_variant_rows,
                            credible_set_variant_rows,
                        )
                    )
        finally:
            variant_handle.close()
            top_handle.close()
            credible_handle.close()

    valid_note = (
        "numerically interpretable PLASMA block; association/FDR screening and "
        "a null-causal model were not applied"
    )
    invalid_note = "PLASMA block failed runner, structural, or numerical QC"
    valid_catalogue = [
        annotate_block(
            row,
            "no_gwas_expression_ase_valid_block",
            valid_note,
        )
        for row in valid_rows
    ]
    invalid_catalogue = [
        annotate_block(
            row,
            "no_gwas_expression_ase_invalid_block_qc",
            invalid_note,
        )
        for row in invalid_rows
    ]
    block_fields = merge_fields(BLOCK_DERIVED_FIELDS, valid_fields)
    invalid_block_fields = merge_fields(BLOCK_DERIVED_FIELDS, invalid_fields)
    tissue_rows = [
        {"tissue": tissue, **tissue_counts[tissue]}
        for tissue in sorted(tissue_counts)
    ]
    unique_top_blocks = sum(
        clean(row.get("top_marker_is_unique")) == "1" for row in valid_rows
    )
    tied_top_blocks = sum(
        clean(row.get("top_marker_is_unique")) == "0" for row in valid_rows
    )

    log("[stage:write] write block, tissue and global QC table")
    write_tsv(VALID_BLOCK_CATALOGUE, block_fields, valid_catalogue)
    write_tsv(INVALID_BLOCK_QC, invalid_block_fields, invalid_catalogue)
    write_tsv(TISSUE_SUMMARY, TISSUE_FIELDS, tissue_rows)

    summary_rows = [
        {"key": "run_id", "value": RUN_ID},
        {"key": "source_run_id", "value": SOURCE_RUN_ID},
        {"key": "source_variant_summary", "value": str(VARIANT_SUMMARY)},
        {"key": "source_valid_blocks", "value": str(VALID_BLOCKS)},
        {"key": "source_invalid_blocks", "value": str(INVALID_BLOCKS)},
        {"key": "source_qc", "value": str(SOURCE_QC)},
        {
            "key": "source_execution_completion_marker",
            "value": str(SOURCE_COMPLETION_MARKER),
        },
        {
            "key": "expression_matrix_scope",
            "value": EXPRESSION_MATRIX_SCOPE,
        },
        {
            "key": "feature_origin_rule",
            "value": (
                "gene-=>reference_annotation_gene;"
                "DOG5_NEWGENE_=>dog5_manual_new_gene"
            ),
        },
        {"key": "source_variant_rows", "value": source_variant_rows},
        {
            "key": "skipped_invalid_variant_rows",
            "value": skipped_invalid_variant_rows,
        },
        {
            "key": "interpretable_variant_rows",
            "value": interpretable_variant_rows,
        },
        {"key": "top_variant_rows", "value": top_variant_rows},
        {
            "key": "credible_set_variant_rows",
            "value": credible_set_variant_rows,
        },
        {"key": "valid_blocks", "value": len(valid_by_id)},
        {"key": "invalid_blocks_qc_only", "value": len(invalid_by_id)},
        {"key": "unique_top_blocks", "value": unique_top_blocks},
        {"key": "tied_top_blocks", "value": tied_top_blocks},
        {
            "key": "ppa_ge_0.95_variant_rows",
            "value": threshold_counts["ppa_ge_0.95_variant_rows"],
        },
        {
            "key": "ppa_ge_0.8_variant_rows",
            "value": threshold_counts["ppa_ge_0.8_variant_rows"],
        },
        {
            "key": "ppa_ge_0.5_variant_rows",
            "value": threshold_counts["ppa_ge_0.5_variant_rows"],
        },
        {
            "key": "ppa_ge_0.1_variant_rows",
            "value": threshold_counts["ppa_ge_0.1_variant_rows"],
        },
        {
            "key": "8.3_search_mode",
            "value": source_qc.get("8.3_search_mode", ""),
        },
        {
            "key": "8.3_max_causal",
            "value": source_qc.get("8.3_max_causal", ""),
        },
        {
            "key": "8.3_confidence",
            "value": source_qc.get("8.3_confidence", ""),
        },
        {
            "key": "association_fdr_screening",
            "value": "not_performed",
        },
        {
            "key": "null_causal_model",
            "value": "not_included",
        },
        {
            "key": "interpretation_scope",
            "value": (
                "exploratory no-GWAS expression+ASE joint fine-mapping candidates"
            ),
        },
        {
            "key": "not_interpretable_as",
            "value": (
                "association-significant regulatory variants, phenotype causal "
                "variants, or causal genes"
            ),
        },
    ]
    write_tsv(CATALOGUE_SUMMARY, ["key", "value"], summary_rows)

    log(
        "[summary] source_rows={0}; catalogue_rows={1}; skipped_invalid_rows={2}".format(
            source_variant_rows,
            interpretable_variant_rows,
            skipped_invalid_variant_rows,
        )
    )
    log(
        "[summary] valid_blocks={0}; invalid_blocks={1}; top_rows={2}; "
        "credible_rows={3}".format(
            len(valid_by_id),
            len(invalid_by_id),
            top_variant_rows,
            credible_set_variant_rows,
        )
    )
    log("[output] {0}".format(VARIANT_CATALOGUE))
    log("[output] {0}".format(TOP_VARIANT_CATALOGUE))
    log("[output] {0}".format(CREDIBLE_SET_CATALOGUE))
    log("[output] {0}".format(CATALOGUE_SUMMARY))
    COMPLETION_MARKER.write_text(
        "run_id={0}\nexecution_status=success\nfinished_at={1}\n".format(
            RUN_ID, datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ),
        encoding="utf-8",
    )
    log("[done] five dogs 8.5 regulatory catalogue build completed; task completed / done")
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
