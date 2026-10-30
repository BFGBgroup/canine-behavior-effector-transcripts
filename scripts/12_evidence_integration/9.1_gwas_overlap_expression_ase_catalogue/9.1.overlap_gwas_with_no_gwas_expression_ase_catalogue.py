#!/usr/bin/python3.10


from __future__ import annotations

import csv
import gzip
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


BASE_DIR = Path("/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal")
SCRIPT_DIR = BASE_DIR / "scripts/9.1_gwas_overlap_expression_ase_catalogue"
GWAS_SUMMARY_DIR = Path("/path/to/project/legacy/reference/gwas/summary_statistics")

SOURCE_EXPORT_RUN_ID = "20260827_224912"
SOURCE_QC_RUN_ID = "20260828_122551"
NO_GWAS_EXPORT_DIR = (
    BASE_DIR / "8.6_no_gwas_export_expression_ase_catalogue_deliverables" / SOURCE_EXPORT_RUN_ID
)
NO_GWAS_CATALOGUE = (
    NO_GWAS_EXPORT_DIR / "no_gwas_expression_ase_priority_variants_ppa_ge_0.1.tsv.gz"
)
NO_GWAS_FINAL_QC_DIR = (
    BASE_DIR
    / "8.7_no_gwas_final_catalogue_qc_report"
    / SOURCE_QC_RUN_ID
)
NO_GWAS_FINAL_QC_SUMMARY = (
    NO_GWAS_FINAL_QC_DIR / "no_gwas_expression_ase_final_qc_summary.tsv"
)
NO_GWAS_FINAL_QC_COMPLETION_MARKER = (
    NO_GWAS_FINAL_QC_DIR / "execution_completion.marker"
)

OUTPUT_ROOT = BASE_DIR / "9.1_gwas_overlap_expression_ase_catalogue"
RUN_ID = os.environ.get("RUN_ID", datetime.now().strftime("%Y%m%d_%H%M%S"))
RUN_DIR = OUTPUT_ROOT / RUN_ID

P_VALUE_COLUMN = "p_lrt"
GENOMEWIDE_P_THRESHOLD = 2.97e-6
SUGGESTIVE_P_THRESHOLD = 1e-5
EXPLORATORY_POINT_P_THRESHOLD = 5e-5
GWAS_CANDIDATE_P_THRESHOLD = 1e-4
TRAIT_FILE_SUFFIX = "_GWAS_summary_statistics.tsv.gz"
EXPECTED_EXPRESSION_MATRIX_SCOPE = (
    "reference_annotation_plus_dog5_manual_new_gene"
)
EXPECTED_FEATURE_ORIGIN_RULE = (
    "gene-=>reference_annotation_gene;DOG5_NEWGENE_=>dog5_manual_new_gene"
)

GWAS_TRAIT_INVENTORY_TSV = RUN_DIR / "gwas_trait_marker_inventory.tsv"
GWAS_CANDIDATES_TSV = RUN_DIR / "gwas_candidate_markers_p_le_1e-4.tsv"
EXACT_ALL_P_TSV = RUN_DIR / "gwas_catalogue_exact_variant_overlaps_all_p.tsv"
EXACT_CANDIDATE_TSV = RUN_DIR / "gwas_catalogue_exact_variant_overlaps_p_le_1e-4.tsv"
ALLELE_MISMATCH_QC_TSV = RUN_DIR / "gwas_catalogue_position_allele_mismatch_qc.tsv"
TRAIT_FEATURE_SUMMARY_TSV = RUN_DIR / "gwas_catalogue_trait_feature_summary.tsv"
RUN_SUMMARY_TSV = RUN_DIR / "gwas_catalogue_overlap_run_summary.tsv"
README_MD = RUN_DIR / "README.md"


GWAS_FINEMAP_RUN = Path(
    '/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352'
)
GWAS_LOCUS_MANIFEST = GWAS_FINEMAP_RUN / "locus_run_manifest.tsv"
POINT_CANDIDATE_TSV = RUN_DIR / "gwas_catalogue_exact_variant_overlaps_p_lt_5e-5.tsv"
REGIONAL_TSV = RUN_DIR / "gwas_catalogue_main_modelled_region_overlaps.tsv"
REGIONAL_SUMMARY_TSV = RUN_DIR / "gwas_catalogue_main_region_scope_summary.tsv"
P_VALUE_ISSUES_TSV = RUN_DIR / "gwas_p_value_issues.tsv"
REGION_COLUMNS = [
    "gwas_locus_id", "gwas_analysis_scope", "gwas_locus_lead_id",
    "gwas_locus_lead_p_lrt", "gwas_locus_lead_p_tier",
    "gwas_membership_source_branch", "gwas_modelled_variant_id",
    "gwas_locus_status", "gwas_locus_model_status",
    "gwas_membership_source_file",
]


CATALOGUE_REQUIRED_COLUMNS = [
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

GWAS_REQUIRED_COLUMNS = [
    "chr",
    "rs",
    "ps",
    "allele1",
    "allele0",
    "af",
    "beta",
    "se",
    "p_wald",
    "p_lrt",
    "p_score",
]


TRAIT_FEATURE_SUMMARY_COLUMNS = [
    "trait", "tissue", "feature_id", "feature_label", "block_count",
    "exact_all_p_overlap_rows", "exact_all_p_marker_count",
    "exact_candidate_overlap_rows", "exact_candidate_marker_count",
    "max_exact_overlap_ppa", "min_exact_all_p_gwas_p", "min_exact_candidate_gwas_p",
]

CATALOGUE_OUTPUT_COLUMNS = list(CATALOGUE_REQUIRED_COLUMNS)

GWAS_OUTPUT_COLUMNS = [
    "trait",
    "gwas_marker_id",
    "gwas_chrom",
    "gwas_pos",
    "gwas_allele1",
    "gwas_allele0",
    "gwas_af",
    "gwas_beta",
    "gwas_se",
    "gwas_p_wald",
    "gwas_p_lrt",
    "gwas_p_score",
    "gwas_or",
    "gwas_or_95ci_lower",
    "gwas_or_95ci_upper",
    "gwas_p_value_column",
    "gwas_p_value",
    "gwas_variant_p_tier",
    "gwas_p_value_status",
]

EXACT_OUTPUT_COLUMNS = CATALOGUE_OUTPUT_COLUMNS + GWAS_OUTPUT_COLUMNS + [
    "normalized_variant_key",
    "overlap_type",
    "allele_match_status",
    "gwas_allele1_catalogue_relation",
    "gwas_allele0_catalogue_relation",
    "allele_order_status",
    "beta_alignment_status",
    "distance_to_gwas_marker_bp",
    "gwas_modelled_locus_count",
    "gwas_region_match_status",
]


def log(message: str) -> None:
    print("[{0} ] {1}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message), flush=True)


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def safe_int(value: object, default: int = -1) -> int:
    try:
        return int(float(clean(value)))
    except (TypeError, ValueError):
        return default


def safe_float(value: object) -> float:
    try:
        return float(clean(value))
    except (TypeError, ValueError):
        return float("nan")


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    return "{0:.8g}".format(value)


def normalize_chrom(value: object) -> str:
    chrom = clean(value)
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom


def chrom_for_output(value: object) -> str:
    chrom = normalize_chrom(value)
    return "chr{0}".format(chrom) if chrom else ""


def allele_set(*values: object) -> Tuple[str, ...]:
    alleles = sorted(
        {
            clean(value).upper()
            for value in values
            if clean(value) and clean(value).upper() != "NA"
        }
    )
    return tuple(alleles)


def exact_allele_audit(
    catalogue_row: Mapping[str, object],
    gwas_row: Mapping[str, object],
) -> Dict[str, str]:
    ref = clean(catalogue_row.get("ref")).upper()
    alt = clean(catalogue_row.get("alt")).upper()
    allele1 = clean(gwas_row.get("allele1")).upper()
    allele0 = clean(gwas_row.get("allele0")).upper()
    chrom = chrom_for_output(catalogue_row.get("chrom"))
    pos = clean(catalogue_row.get("pos"))

    if allele1 == alt and allele0 == ref:
        order_status = "allele1_ALT_allele0_REF"
        allele1_relation = "ALT"
        allele0_relation = "REF"
    elif allele1 == ref and allele0 == alt:
        order_status = "allele1_REF_allele0_ALT"
        allele1_relation = "REF"
        allele0_relation = "ALT"
    else:
        order_status = "unexpected_after_allele_set_match"
        allele1_relation = "UNRESOLVED"
        allele0_relation = "UNRESOLVED"

    return {
        "normalized_variant_key": "{0}:{1}:{2}".format(
            chrom,
            pos,
            "/".join(allele_set(ref, alt)),
        ),
        "gwas_allele1_catalogue_relation": allele1_relation,
        "gwas_allele0_catalogue_relation": allele0_relation,
        "allele_order_status": order_status,
        "beta_alignment_status": "not_aligned_effect_allele_semantics_unverified",
    }


def p_value_status(value: object) -> str:

    token = clean(value)
    if token.lower() in {"", "na", "nan", ".", "none"}:
        return "missing"
    number = safe_float(token)
    return "valid" if math.isfinite(number) and 0 <= number <= 1 else "invalid"


def gwas_variant_p_tier(p_value: float) -> str:
    if not math.isfinite(p_value) or not 0 <= p_value <= 1:
        return "unavailable"
    if p_value < GENOMEWIDE_P_THRESHOLD:
        return "genomewide"
    if p_value < SUGGESTIVE_P_THRESHOLD:
        return "suggestive"
    if p_value < EXPLORATORY_POINT_P_THRESHOLD:
        return "exploratory_point"
    if p_value <= GWAS_CANDIDATE_P_THRESHOLD:
        return "extended_exploratory_point"
    return "background"


def load_region_memberships() -> Dict[Tuple, List[Dict[str, str]]]:
    fields, manifest = read_tsv(GWAS_LOCUS_MANIFEST)
    required_columns(GWAS_LOCUS_MANIFEST, fields, [
        "branch", "locus_id", "trait", "analysis_scope", "gwas_lead_variant_id",
        "gwas_lead_p_lrt", "gwas_lead_p_lrt_tier", "variant_pip", "status", "model_status",
    ])
    memberships: Dict[Tuple, List[Dict[str, str]]] = {}
    seen_loci = set()
    for locus in manifest:
        if locus["branch"] != "main":
            continue
        locus_key = (locus["trait"], locus["locus_id"])
        if locus_key in seen_loci:
            raise ValueError("duplicate main region: {0}".format(locus_key))
        seen_loci.add(locus_key)
        if locus["analysis_scope"] not in {"primary", "exploratory"}:
            raise ValueError("unknown region span: {0}".format(locus["analysis_scope"]))
        path = Path(locus["variant_pip"])
        fields, variants = read_tsv(path)
        required_columns(path, fields, [
            "locus_id", "trait", "analysis_scope", "chrom", "pos", "ld_a1", "ld_a2", "ld_variant_id",
        ])
        seen_variants = set()
        for variant in variants:
            if any(variant[k] != locus[k] for k in ("locus_id", "trait", "analysis_scope")):
                raise ValueError("region manifest and variant_pip identity mismatch: {0}".format(path))
            alleles = allele_set(variant["ld_a1"], variant["ld_a2"])
            position = safe_int(variant["pos"])
            if position <= 0 or len(alleles) != 2:
                raise ValueError("model variant identity invalid: {0}".format(path))
            key = (variant["trait"], normalize_chrom(variant["chrom"]), position, alleles)

            if key in seen_variants:
                raise ValueError("the same region duplicate variant identity: {0}: {1}".format(path, key))
            seen_variants.add(key)
            memberships.setdefault(key, []).append({
                "gwas_locus_id": locus["locus_id"],
                "gwas_analysis_scope": locus["analysis_scope"],
                "gwas_locus_lead_id": locus["gwas_lead_variant_id"],
                "gwas_locus_lead_p_lrt": locus["gwas_lead_p_lrt"],
                "gwas_locus_lead_p_tier": locus["gwas_lead_p_lrt_tier"],
                "gwas_membership_source_branch": "main",
                "gwas_modelled_variant_id": variant["ld_variant_id"],
                "gwas_locus_status": locus["status"],
                "gwas_locus_model_status": locus["model_status"],
                "gwas_membership_source_file": str(path),
            })
    if not seen_loci:
        raise ValueError("region list no main branch")
    log("[region] main_loci={0}; variant_trait_keys={1}".format(len(seen_loci), len(memberships)))
    return memberships


def annotate_regions(rows: Sequence[Dict[str, object]], memberships: Mapping) -> List[Dict[str, object]]:
    regional = []
    for row in rows:
        key = (row["trait"], normalize_chrom(row["gwas_chrom"]), safe_int(row["gwas_pos"]),
               allele_set(row["gwas_allele1"], row["gwas_allele0"]))
        matches = memberships.get(key, [])
        row["gwas_modelled_locus_count"] = len(matches)
        row["gwas_region_match_status"] = (
            "matched_main_modelled_variant" if matches else "not_found_in_main_modelled_variants"
        )
        for match in matches:

            regional.append({**row, **match,
                             "overlap_type": "exact_variant_main_modelled_region_allele_set_match"})
    return regional


def trait_from_path(path: Path) -> str:
    if path.name.endswith(TRAIT_FILE_SUFFIX):
        return path.name[: -len(TRAIT_FILE_SUFFIX)]
    return path.name.replace(".tsv.gz", "")


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError("missing or empty TSV: {0}".format(path))
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", newline="")
    else:
        handle = path.open(newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = list(reader.fieldnames or [])
    return fields, rows


def write_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def required_columns(path: Path, observed: Sequence[str], required: Sequence[str]) -> None:
    missing = sorted(set(required) - set(observed))
    if missing:
        raise ValueError("{0} missing required columns: {1}".format(path, ",".join(missing)))


def validate_run_id(name: str, value: str) -> None:
    if len(value) != 15 or value[8] != "_" or not value.replace("_", "").isdigit():
        raise ValueError("{0} must use YYYYMMDD_HHMMSS format: {1}".format(name, value))


def validate_inputs() -> None:
    validate_run_id("SOURCE_EXPORT_RUN_ID", SOURCE_EXPORT_RUN_ID)
    validate_run_id("SOURCE_QC_RUN_ID", SOURCE_QC_RUN_ID)
    validate_run_id("RUN_ID", RUN_ID)
    for path in (
        NO_GWAS_CATALOGUE,
        GWAS_LOCUS_MANIFEST,
        NO_GWAS_FINAL_QC_SUMMARY,
        NO_GWAS_FINAL_QC_COMPLETION_MARKER,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("missing or empty input: {0}".format(path))
    if not GWAS_SUMMARY_DIR.is_dir():
        raise FileNotFoundError("missing GWAS summary directory: {0}".format(GWAS_SUMMARY_DIR))


    if RUN_DIR.exists():
        unexpected = list(RUN_DIR.iterdir())
        if unexpected:
            raise FileExistsError("run directory already contains outputs: {0}".format(RUN_DIR))


def load_key_value_tsv(path: Path) -> Dict[str, str]:
    fields, rows = read_tsv(path)
    required_columns(path, fields, ["key", "value"])
    return {clean(row.get("key")): clean(row.get("value")) for row in rows}


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


def require_final_qc_pass(
    summary_path: Path,
    completion_marker_path: Path,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    summary = load_key_value_tsv(summary_path)
    completion = load_completion_marker(completion_marker_path)
    status = clean(summary.get("final_qc_status"))
    issue_count = safe_int(summary.get("issue_count"), default=-1)
    if status != "PASS" or issue_count != 0:
        raise ValueError(
            "no-GWAS final QC did not pass: status={0}; issue_count={1}".format(
                status, issue_count
            )
        )
    if clean(summary.get("source_8.6_run_id")) != SOURCE_EXPORT_RUN_ID:
        raise ValueError("8.7 QC does not point to the configured 8.6 run")
    if clean(summary.get("run_id")) != SOURCE_QC_RUN_ID:
        raise ValueError("8.7 QC summary does not match the configured 8.7 run")
    if (
        clean(summary.get("source_8.6_expression_matrix_scope"))
        != EXPECTED_EXPRESSION_MATRIX_SCOPE
    ):
        raise ValueError("8.7 QC contains an unexpected expression matrix scope")
    if (
        clean(summary.get("source_8.6_feature_origin_rule"))
        != EXPECTED_FEATURE_ORIGIN_RULE
    ):
        raise ValueError("8.7 QC contains an unexpected feature origin rule")
    if clean(completion.get("run_id")) != SOURCE_QC_RUN_ID:
        raise ValueError("8.7 completion marker does not match the configured run")
    if clean(completion.get("execution_status")) != "success":
        raise ValueError("8.7 completion marker is not successful")
    if clean(completion.get("final_qc_status")) != "PASS":
        raise ValueError("8.7 completion marker final QC status is not PASS")
    return summary, completion


def expected_feature_origin(feature_id: str) -> str:
    if feature_id.startswith("DOG5_NEWGENE_"):
        return "dog5_manual_new_gene"
    if feature_id.startswith("gene-"):
        return "reference_annotation_gene"
    return ""


def build_catalogue_payload(row: Mapping[str, str]) -> Dict[str, object]:
    payload = {column: clean(row.get(column)) for column in CATALOGUE_OUTPUT_COLUMNS}
    payload["chrom"] = chrom_for_output(row.get("chrom"))
    payload["ref"] = clean(row.get("ref")).upper()
    payload["alt"] = clean(row.get("alt")).upper()
    return payload


def build_gwas_payload(
    trait: str,
    row: Mapping[str, str],
    p_value: float,
) -> Dict[str, object]:
    return {
        "trait": trait,
        "gwas_marker_id": clean(row.get("rs")),
        "gwas_chrom": chrom_for_output(row.get("chr")),
        "gwas_pos": clean(row.get("ps")),
        "gwas_allele1": clean(row.get("allele1")).upper(),
        "gwas_allele0": clean(row.get("allele0")).upper(),
        "gwas_af": clean(row.get("af")),
        "gwas_beta": clean(row.get("beta")),
        "gwas_se": clean(row.get("se")),
        "gwas_p_wald": clean(row.get("p_wald")),
        "gwas_p_lrt": clean(row.get("p_lrt")),
        "gwas_p_score": clean(row.get("p_score")),
        "gwas_or": clean(row.get("OR")),
        "gwas_or_95ci_lower": clean(row.get("OR_95CI_LOWER")),
        "gwas_or_95ci_upper": clean(row.get("OR_95CI_UPPER")),
        "gwas_p_value_column": P_VALUE_COLUMN,
        "gwas_p_value": format_float(p_value),
        "gwas_variant_p_tier": (gwas_variant_p_tier(p_value) if p_value_status(row.get(P_VALUE_COLUMN)) == "valid" else p_value_status(row.get(P_VALUE_COLUMN)) + "_p"),
        "gwas_p_value_status": p_value_status(row.get(P_VALUE_COLUMN)),
    }


def load_catalogue(
    path: Path,
) -> Tuple[List[Dict[str, str]], Dict[Tuple[str, int], List[Dict[str, str]]]]:
    fields, rows = read_tsv(path)
    required_columns(path, fields, CATALOGUE_REQUIRED_COLUMNS)
    kept_rows: List[Dict[str, str]] = []
    by_position: Dict[Tuple[str, int], List[Dict[str, str]]] = {}
    invalid_lines: List[int] = []

    for line_number, row in enumerate(rows, start=2):
        chrom = normalize_chrom(row.get("chrom"))
        pos = safe_int(row.get("pos"))
        ppa = safe_float(row.get("ppa"))
        feature_id = clean(row.get("feature_id"))
        feature_origin = expected_feature_origin(feature_id)
        valid = (
            chrom
            and pos > 0
            and len(allele_set(row.get("ref"), row.get("alt"))) == 2
            and feature_origin
            and clean(row.get("feature_origin")) == feature_origin
            and clean(row.get("expression_matrix_scope"))
            == EXPECTED_EXPRESSION_MATRIX_SCOPE
            and clean(row.get("result_status")) == "valid"
            and clean(row.get("result_interpretable")) == "1"
            and math.isfinite(ppa)
            and 0.1 <= ppa <= 1.0 + 1e-8
        )
        if not valid:
            invalid_lines.append(line_number)
            continue
        kept_rows.append(row)
        by_position.setdefault((chrom, pos), []).append(row)

    if invalid_lines:
        raise ValueError(
            "priority catalogue contains {0} invalid rows; first lines: {1}".format(
                len(invalid_lines),
                ",".join(str(value) for value in invalid_lines[:10]),
            )
        )
    return kept_rows, by_position


def discover_gwas_files(path: Path) -> List[Path]:
    files = sorted(path.glob("*_GWAS_summary_statistics.tsv.gz"))
    if not files:
        raise FileNotFoundError("no GWAS summary files found in {0}".format(path))
    return files


def stream_gwas(
    gwas_files: Sequence[Path],
    catalogue_by_position: Mapping[Tuple[str, int], Sequence[Dict[str, str]]],
    p_issue_writer: csv.DictWriter,
) -> Tuple[
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
]:
    inventory_rows: List[Dict[str, object]] = []
    candidate_rows: List[Dict[str, object]] = []
    exact_all_rows: List[Dict[str, object]] = []
    exact_candidate_rows: List[Dict[str, object]] = []
    mismatch_rows: List[Dict[str, object]] = []

    for file_index, path in enumerate(gwas_files, start=1):
        trait = trait_from_path(path)
        log("[1/3] scanning GWAS {0}/{1}: {2}".format(file_index, len(gwas_files), path))
        counts = {
            "row_count": 0,
            "finite_p_count": 0,
            "genomewide_count": 0,
            "strong_count": 0,
            "loose_count": 0,
            "candidate_count": 0,
            "position_hit_rows": 0,
            "exact_all_rows": 0,
            "exact_candidate_rows": 0,
            "mismatch_all_rows": 0,
        }
        min_p = float("inf")

        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            required_columns(path, fields, GWAS_REQUIRED_COLUMNS)

            for line_number, row in enumerate(reader, start=2):
                counts["row_count"] += 1
                chrom = normalize_chrom(row.get("chr"))
                pos = safe_int(row.get("ps"))
                p_value = safe_float(row.get(P_VALUE_COLUMN))
                status = p_value_status(row.get(P_VALUE_COLUMN))
                finite_p = status == "valid"
                if not finite_p:
                    p_issue_writer.writerow({"source_file": str(path), "line_number": line_number,
                                             "trait": trait, "gwas_marker_id": row.get("rs", ""),
                                             "raw_p_lrt": row.get(P_VALUE_COLUMN, ""), "status": status})
                    p_value = float("nan")

                if finite_p:
                    counts["finite_p_count"] += 1
                    min_p = min(min_p, p_value)
                    if p_value < GENOMEWIDE_P_THRESHOLD:
                        counts["genomewide_count"] += 1
                    if p_value < SUGGESTIVE_P_THRESHOLD:
                        counts["strong_count"] += 1
                    if p_value < EXPLORATORY_POINT_P_THRESHOLD:
                        counts["loose_count"] += 1
                    if p_value <= GWAS_CANDIDATE_P_THRESHOLD:
                        counts["candidate_count"] += 1
                        candidate_rows.append(build_gwas_payload(trait, row, p_value))

                if not chrom or pos <= 0:
                    continue
                catalogue_rows = catalogue_by_position.get((chrom, pos))
                if not catalogue_rows:
                    continue

                counts["position_hit_rows"] += len(catalogue_rows)
                gwas_alleles = allele_set(row.get("allele1"), row.get("allele0"))
                for catalogue_row in catalogue_rows:
                    catalogue_alleles = allele_set(
                        catalogue_row.get("ref"),
                        catalogue_row.get("alt"),
                    )
                    if len(gwas_alleles) == 2 and gwas_alleles == catalogue_alleles:
                        exact_row = {
                            **build_catalogue_payload(catalogue_row),
                            **build_gwas_payload(trait, row, p_value),
                            **exact_allele_audit(catalogue_row, row),
                            "overlap_type": "exact_variant_all_p_allele_set_match",
                            "allele_match_status": "allele_set_match",
                            "distance_to_gwas_marker_bp": 0,
                        }
                        exact_all_rows.append(exact_row)
                        counts["exact_all_rows"] += 1
                        if finite_p and p_value <= GWAS_CANDIDATE_P_THRESHOLD:
                            candidate_exact_row = dict(exact_row)
                            candidate_exact_row["overlap_type"] = (
                                "exact_variant_p_le_1e-4_allele_set_match"
                            )
                            exact_candidate_rows.append(candidate_exact_row)
                            counts["exact_candidate_rows"] += 1
                    else:
                        mismatch_rows.append(
                            {
                                "trait": trait,
                                "catalogue_marker_id": clean(catalogue_row.get("marker_id")),
                                "catalogue_chrom": chrom_for_output(catalogue_row.get("chrom")),
                                "catalogue_pos": clean(catalogue_row.get("pos")),
                                "catalogue_ref": clean(catalogue_row.get("ref")).upper(),
                                "catalogue_alt": clean(catalogue_row.get("alt")).upper(),
                                "gwas_marker_id": clean(row.get("rs")),
                                "gwas_chrom": chrom_for_output(row.get("chr")),
                                "gwas_pos": clean(row.get("ps")),
                                "gwas_allele1": clean(row.get("allele1")).upper(),
                                "gwas_allele0": clean(row.get("allele0")).upper(),
                                "gwas_p_value_column": P_VALUE_COLUMN,
                                "gwas_p_value": format_float(p_value),
                                "gwas_variant_p_tier": gwas_variant_p_tier(p_value),
                                "allele_match_status": "position_match_allele_set_mismatch",
                            }
                        )
                        counts["mismatch_all_rows"] += 1

        inventory_rows.append(
            {
                "trait": trait,
                "gwas_summary_file": str(path),
                "row_count": counts["row_count"],
                "valid_p_value_count": counts["finite_p_count"],
                "p_value_column": P_VALUE_COLUMN,
                "min_p_value": format_float(min_p),
                "genomewide_p_threshold": GENOMEWIDE_P_THRESHOLD,
                "genomewide_p_marker_count": counts["genomewide_count"],
                "suggestive_p_threshold": SUGGESTIVE_P_THRESHOLD,
                "cumulative_p_lt_1e-5_marker_count": counts["strong_count"],
                "exploratory_point_p_threshold": EXPLORATORY_POINT_P_THRESHOLD,
                "cumulative_p_lt_5e-5_marker_count": counts["loose_count"],
                "candidate_p_threshold": GWAS_CANDIDATE_P_THRESHOLD,
                "candidate_marker_count": counts["candidate_count"],
                "exact_position_hit_catalogue_row_count_all_p": counts["position_hit_rows"],
                "exact_allele_match_catalogue_row_count_all_p": counts["exact_all_rows"],
                "exact_allele_match_catalogue_row_count_p_le_1e-4": counts[
                    "exact_candidate_rows"
                ],
                "position_allele_mismatch_catalogue_row_count_all_p": counts[
                    "mismatch_all_rows"
                ],
            }
        )

    return (
        inventory_rows,
        candidate_rows,
        exact_all_rows,
        exact_candidate_rows,
        mismatch_rows,
    )


def max_finite(current: object, value: float) -> float:
    current_value = safe_float(current)
    if not math.isfinite(current_value):
        return value
    if not math.isfinite(value):
        return current_value
    return max(current_value, value)


def min_finite(current: object, value: float) -> float:
    current_value = safe_float(current)
    if not math.isfinite(current_value):
        return value
    if not math.isfinite(value):
        return current_value
    return min(current_value, value)


def summarize_trait_feature_overlaps(
    exact_all_rows: Sequence[Mapping[str, object]],
    exact_candidate_rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    grouped: MutableMapping[Tuple[str, str, str, str], Dict[str, object]] = {}

    def ensure(row: Mapping[str, object]) -> Dict[str, object]:
        key = (
            clean(row.get("trait")),
            clean(row.get("tissue")),
            clean(row.get("feature_id")),
            clean(row.get("feature_label")),
        )
        if key not in grouped:
            grouped[key] = {
                "trait": key[0],
                "tissue": key[1],
                "feature_id": key[2],
                "feature_label": key[3],
                "block_ids": set(),
                "exact_all_rows": 0,
                "exact_all_markers": set(),
                "exact_candidate_rows": 0,
                "exact_candidate_markers": set(),
                "max_exact_ppa": float("nan"),
                "min_exact_all_p": float("nan"),
                "min_exact_candidate_p": float("nan"),
            }
        return grouped[key]

    for row in exact_all_rows:
        item = ensure(row)
        item["block_ids"].add(clean(row.get("block_id")))
        item["exact_all_rows"] += 1
        item["exact_all_markers"].add(clean(row.get("marker_id")))
        item["max_exact_ppa"] = max_finite(item["max_exact_ppa"], safe_float(row.get("ppa")))
        item["min_exact_all_p"] = min_finite(
            item["min_exact_all_p"],
            (safe_float(row.get("gwas_p_lrt"))
             if row.get("gwas_p_value_status") == "valid" else float("nan")),
        )

    for row in exact_candidate_rows:
        item = ensure(row)
        item["exact_candidate_rows"] += 1
        item["exact_candidate_markers"].add(clean(row.get("marker_id")))
        item["min_exact_candidate_p"] = min_finite(
            item["min_exact_candidate_p"],
            (safe_float(row.get("gwas_p_lrt"))
             if row.get("gwas_p_value_status") == "valid" else float("nan")),
        )

    summary_rows: List[Dict[str, object]] = []
    for item in grouped.values():
        summary_rows.append(
            {
                "trait": item["trait"],
                "tissue": item["tissue"],
                "feature_id": item["feature_id"],
                "feature_label": item["feature_label"],
                "block_count": len(item["block_ids"]),
                "exact_all_p_overlap_rows": item["exact_all_rows"],
                "exact_all_p_marker_count": len(item["exact_all_markers"]),
                "exact_candidate_overlap_rows": item["exact_candidate_rows"],
                "exact_candidate_marker_count": len(item["exact_candidate_markers"]),
                "max_exact_overlap_ppa": format_float(item["max_exact_ppa"]),
                "min_exact_all_p_gwas_p": format_float(item["min_exact_all_p"]),
                "min_exact_candidate_gwas_p": format_float(item["min_exact_candidate_p"]),
            }
        )
    return sorted(
        summary_rows,
        key=lambda row: (
            clean(row.get("trait")),
            clean(row.get("tissue")),
            clean(row.get("feature_id")),
        ),
    )


def write_readme(path: Path, summary: Mapping[str, object]) -> None:
    text = """# Step 9.1 GWAS overlap with expression+ASE catalogue

This step reports exact variant overlaps between the fine-map catalogue and GWAS, using all
GWAS summary statistics to find variants at the same coordinate with the same unordered allele set. The `p_lrt <= 1e-4`
candidate list contains exploratory GWAS markers only; no physical-distance window overlap is constructed.

Interpretation limits:

- The all-p exact-match table indicates catalogue variants with GWAS summary records at the same coordinate and with the same allele set; P values may be missing or invalid and must be interpreted using gwas_p_value_status.
- Results with `p_lrt <= 1e-4` are exploratory candidates, not necessarily genome-wide significant.
- Input is restricted to molecular fine-mapping PPA>=0.1; all_p means no restriction on GWAS P values, not inclusion of all molecular variants.
- Single-variant categories: P<2.97e-6 is genomewide; up to <1e-5 is suggestive; up to <5e-5 is exploratory_point; up to <=1e-4 is extended_exploratory_point; all others are background.
- cumulative_p_lt_1e-5 and cumulative_p_lt_5e-5 in inventory are cumulative counts that include more significant sites and must not be added as mutually exclusive categories.
- exact_candidate and min_exact_candidate in trait_feature_summary use the extended P<=1e-4 tier, not the regional primary tier.
- gwas_variant_p_tier and gwas_analysis_scope are independent; the latter is inherited from the main manifest and is not recalculated from the overlapping SNP P value.
- gwas_catalogue_exact_variant_overlaps_p_lt_5e-5.tsv contains single-variant candidates; p_le_1e-4 is the extended exploratory tier.
- gwas_catalogue_main_modelled_region_overlaps.tsv joins the main variant table by trait, coordinate, and allele set; one overlap may belong to multiple regions, so expanded regional row counts are not unique-variant counts.
- The regional layer applies no single-variant P threshold or CS filter and does not represent shared credible-set membership or colocalization.
- not_found_in_main_modelled_variants in all_p does not imply a physical location outside regions; this step does not infer interval membership.
- Missing or invalid P values are reported in gwas_p_value_issues.tsv and excluded from single-variant candidates; exact identity and model membership may still be retained.
- Main-branch manifest: {gwas_manifest}; regional summaries include primary and exploratory regions and must not be combined into genome-wide significant results.
- Exact matching allows different allele orders but requires identical unordered allele sets.
- Outputs record whether allele1/allele0 corresponds to REF or ALT; GWAS beta is not direction-aligned.

Key settings:

- GWAS p-value column: `{p_value_column}`
- Candidate p threshold: `{candidate_p_threshold}`
- Source 8.6 run: `{source_export_run_id}`
- Source 8.7 run: `{source_qc_run_id}`
- Expression matrix scope: `{expression_matrix_scope}`
- Genome build assumption: CanFam3.1

Run summary:

- GWAS trait files: `{gwas_file_count}`
- Catalogue rows: `{catalogue_row_count}`
- GWAS candidate markers: `{candidate_rows}`
- Exact overlap rows, all p: `{exact_all_rows}`
- Exact overlap rows, p<=1e-4: `{exact_candidate_rows}`
- Position/allele mismatch rows, all p: `{mismatch_rows}`
""".format(
        gwas_manifest=GWAS_LOCUS_MANIFEST,
        p_value_column=P_VALUE_COLUMN,
        candidate_p_threshold=GWAS_CANDIDATE_P_THRESHOLD,
        source_export_run_id=SOURCE_EXPORT_RUN_ID,
        source_qc_run_id=SOURCE_QC_RUN_ID,
        expression_matrix_scope=EXPECTED_EXPRESSION_MATRIX_SCOPE,
        gwas_file_count=summary.get("gwas_file_count", 0),
        catalogue_row_count=summary.get("catalogue_row_count", 0),
        candidate_rows=summary.get("candidate_rows", 0),
        exact_all_rows=summary.get("exact_all_rows", 0),
        exact_candidate_rows=summary.get("exact_candidate_rows", 0),
        mismatch_rows=summary.get("mismatch_rows", 0),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def main() -> None:
    log("Step 9.1 GWAS overlap started")
    validate_inputs()
    log("script_dir={0}".format(SCRIPT_DIR))
    log("source_export_run_id={0}".format(SOURCE_EXPORT_RUN_ID))
    log("source_qc_run_id={0}".format(SOURCE_QC_RUN_ID))
    log("run_id={0}".format(RUN_ID))
    log("gwas_summary_dir={0}".format(GWAS_SUMMARY_DIR))
    log("no_gwas_catalogue={0}".format(NO_GWAS_CATALOGUE))
    log("no_gwas_final_qc_summary={0}".format(NO_GWAS_FINAL_QC_SUMMARY))
    log(
        "no_gwas_final_qc_completion_marker={0}".format(
            NO_GWAS_FINAL_QC_COMPLETION_MARKER
        )
    )
    log("output_dir={0}".format(RUN_DIR))

    log("[0/3] validating frozen catalogue and loading inputs")
    final_qc_summary, final_qc_completion = require_final_qc_pass(
        NO_GWAS_FINAL_QC_SUMMARY,
        NO_GWAS_FINAL_QC_COMPLETION_MARKER,
    )
    catalogue_rows, catalogue_by_position = load_catalogue(NO_GWAS_CATALOGUE)
    expected_rows = safe_int(
        final_qc_summary.get("observed_rows_priority_ppa_ge_0.1"),
        default=-1,
    )
    if expected_rows < 0 or len(catalogue_rows) != expected_rows:
        raise ValueError(
            "priority catalogue row count mismatch: observed={0}; expected={1}".format(
                len(catalogue_rows),
                expected_rows,
            )
        )
    gwas_files = discover_gwas_files(GWAS_SUMMARY_DIR)
    log(
        "[0/3] catalogue_rows={0}; unique_positions={1}; gwas_files={2}".format(
            len(catalogue_rows),
            len(catalogue_by_position),
            len(gwas_files),
        )
    )

    memberships = load_region_memberships()
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    with P_VALUE_ISSUES_TSV.open("x", newline="") as issue_handle:
        p_issue_writer = csv.DictWriter(issue_handle, delimiter="\t", lineterminator="\n",
                                       fieldnames=["source_file", "line_number", "trait", "gwas_marker_id", "raw_p_lrt", "status"])
        p_issue_writer.writeheader()
        inventory_rows, candidate_rows, exact_all_rows, exact_candidate_rows, mismatch_rows = stream_gwas(
            gwas_files, catalogue_by_position, p_issue_writer)
    regional_rows = annotate_regions(exact_all_rows, memberships)
    annotate_regions(exact_candidate_rows, memberships)

    point_rows = [{**row, "overlap_type": "exact_variant_p_lt_5e-5_allele_set_match"}
                  for row in exact_all_rows if row["gwas_p_value_status"] == "valid"
                  and safe_float(row["gwas_p_lrt"]) < EXPLORATORY_POINT_P_THRESHOLD]
    write_tsv(POINT_CANDIDATE_TSV, EXACT_OUTPUT_COLUMNS, point_rows)
    write_tsv(REGIONAL_TSV, EXACT_OUTPUT_COLUMNS + REGION_COLUMNS, regional_rows)
    scope_rows = []
    for scope in ("primary", "exploratory"):
        selected = [row for row in regional_rows if row["gwas_analysis_scope"] == scope]
        scope_rows.append({
            "gwas_analysis_scope": scope, "region_expanded_rows": len(selected),
            "unique_variants": len({r["normalized_variant_key"] for r in selected}),
            "unique_trait_variants": len({(r["trait"], r["normalized_variant_key"]) for r in selected}),
            "unique_loci": len({(r["trait"], r["gwas_locus_id"]) for r in selected}),
            "unique_molecular_blocks": len({r["block_id"] for r in selected}),
        })
    write_tsv(REGIONAL_SUMMARY_TSV, list(scope_rows[0]), scope_rows)

    log("[2/3] summarizing exact variant overlaps")
    summary_rows = summarize_trait_feature_overlaps(
        exact_all_rows,
        exact_candidate_rows,
    )

    log("[2/3] writing outputs")
    write_tsv(
        GWAS_TRAIT_INVENTORY_TSV,
        list(inventory_rows[0].keys()) if inventory_rows else [],
        inventory_rows,
    )
    write_tsv(GWAS_CANDIDATES_TSV, GWAS_OUTPUT_COLUMNS, candidate_rows)
    write_tsv(EXACT_ALL_P_TSV, EXACT_OUTPUT_COLUMNS, exact_all_rows)
    write_tsv(EXACT_CANDIDATE_TSV, EXACT_OUTPUT_COLUMNS, exact_candidate_rows)
    write_tsv(
        ALLELE_MISMATCH_QC_TSV,
        [
            "trait",
            "catalogue_marker_id",
            "catalogue_chrom",
            "catalogue_pos",
            "catalogue_ref",
            "catalogue_alt",
            "gwas_marker_id",
            "gwas_chrom",
            "gwas_pos",
            "gwas_allele1",
            "gwas_allele0",
            "gwas_p_value_column",
            "gwas_p_value",
            "gwas_variant_p_tier",
            "allele_match_status",
        ],
        mismatch_rows,
    )
    write_tsv(
        TRAIT_FEATURE_SUMMARY_TSV,
        TRAIT_FEATURE_SUMMARY_COLUMNS,
        summary_rows,
    )

    run_summary_rows = [
        {"key": "gwas_locus_manifest", "value": str(GWAS_LOCUS_MANIFEST)},
        {"key": "region_membership_branch", "value": "main"},
        {"key": "region_join_rule", "value": "trait_CHROM_POS_unordered_allele_set;no_point_p_filter;not_CS_filter"},
        {"key": "point_candidates_p_lt_5e-5_rows", "value": len(point_rows)},
        {"key": "region_expanded_rows", "value": len(regional_rows)},
        {"key": "exact_rows_without_main_model_membership", "value": sum(r["gwas_modelled_locus_count"] == 0 for r in exact_all_rows)},
        {"key": "exact_unique_variants", "value": len({r["normalized_variant_key"] for r in exact_all_rows})},
        {"key": "exact_unique_trait_variants", "value": len({(r["trait"], r["normalized_variant_key"]) for r in exact_all_rows})},
        {"key": "exact_unique_variant_gene_tissue", "value": len({(r["normalized_variant_key"], r["feature_id"], r["tissue"]) for r in exact_all_rows})},
        {"key": "run_id", "value": RUN_ID},
        {"key": "source_8.6_run_id", "value": SOURCE_EXPORT_RUN_ID},
        {"key": "source_8.7_run_id", "value": SOURCE_QC_RUN_ID},
        {
            "key": "expression_matrix_scope",
            "value": EXPECTED_EXPRESSION_MATRIX_SCOPE,
        },
        {"key": "feature_origin_rule", "value": EXPECTED_FEATURE_ORIGIN_RULE},
        {"key": "gwas_summary_dir", "value": str(GWAS_SUMMARY_DIR)},
        {"key": "no_gwas_catalogue", "value": str(NO_GWAS_CATALOGUE)},
        {"key": "no_gwas_final_qc_summary", "value": str(NO_GWAS_FINAL_QC_SUMMARY)},
        {
            "key": "no_gwas_final_qc_completion_marker",
            "value": str(NO_GWAS_FINAL_QC_COMPLETION_MARKER),
        },
        {
            "key": "no_gwas_final_qc_execution_status",
            "value": final_qc_completion.get("execution_status", ""),
        },
        {"key": "no_gwas_final_qc_status", "value": final_qc_summary.get("final_qc_status", "")},
        {"key": "no_gwas_final_qc_issue_count", "value": final_qc_summary.get("issue_count", "")},
        {"key": "p_value_column", "value": P_VALUE_COLUMN},
        {"key": "genomewide_p_threshold", "value": GENOMEWIDE_P_THRESHOLD},
        {"key": "suggestive_p_threshold", "value": SUGGESTIVE_P_THRESHOLD},
        {"key": "exploratory_point_p_threshold", "value": EXPLORATORY_POINT_P_THRESHOLD},
        {"key": "gwas_candidate_p_threshold", "value": GWAS_CANDIDATE_P_THRESHOLD},
        {"key": "proximity_analysis_applied", "value": "0"},
        {"key": "primary_result_type", "value": "exact_variant_allele_set_match"},
        {"key": "allele_matching_rule", "value": "CHROM_POS_unordered_allele_set"},
        {
            "key": "beta_alignment_status",
            "value": "not_aligned_effect_allele_semantics_unverified",
        },
        {"key": "gwas_file_count", "value": len(gwas_files)},
        {"key": "catalogue_row_count", "value": len(catalogue_rows)},
        {"key": "catalogue_unique_position_count", "value": len(catalogue_by_position)},
        {"key": "gwas_candidate_marker_rows", "value": len(candidate_rows)},
        {"key": "exact_overlap_rows_all_p", "value": len(exact_all_rows)},
        {"key": "exact_overlap_rows_p_le_1e-4", "value": len(exact_candidate_rows)},
        {"key": "allele_mismatch_rows_all_p", "value": len(mismatch_rows)},
        {"key": "trait_feature_summary_rows", "value": len(summary_rows)},
        {
            "key": "interpretation_boundary",
            "value": "only exact coordinate plus unordered allele-set matches are reported; physical proximity alone is not treated as overlap",
        },
    ]
    write_tsv(RUN_SUMMARY_TSV, ["key", "value"], run_summary_rows)
    write_readme(
        README_MD,
        {
            "gwas_file_count": len(gwas_files),
            "catalogue_row_count": len(catalogue_rows),
            "candidate_rows": len(candidate_rows),
            "exact_all_rows": len(exact_all_rows),
            "exact_candidate_rows": len(exact_candidate_rows),
            "mismatch_rows": len(mismatch_rows),
        },
    )

    log("[3/3] gwas_candidate_rows={0}".format(len(candidate_rows)))
    log("[3/3] exact_overlap_rows_all_p={0}".format(len(exact_all_rows)))
    log("[3/3] exact_overlap_rows_p_le_1e-4={0}".format(len(exact_candidate_rows)))
    log("[3/3] outputs={0}".format(RUN_DIR))
    log("[done] Step 9.1 GWAS overlap finished")


if __name__ == "__main__":
    main()
