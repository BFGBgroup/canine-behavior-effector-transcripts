#!/usr/bin/python3.10


from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, MutableMapping, Sequence, Set, TextIO, Tuple


DEFAULT_EQTL_EXPORT_DIR = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.6_no_gwas_export_expression_ase_catalogue_deliverables/20260827_224912'
)
DEFAULT_EQTL_FINAL_QC = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.7_no_gwas_final_catalogue_qc_report/20260828_122551/no_gwas_expression_ase_final_qc_summary.tsv'
)
DEFAULT_GWAS_RUN_DIR = Path(
    '/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352'
)

EQTL_REQUIRED_COLUMNS = (
    "catalogue_layer", "catalogue_role", "feature_origin",
    "expression_matrix_scope", "is_top_variant", "is_credible_set_member",
    "ppa_tier", "ppa_status", "ppa_boundary_adjusted_count", "block_max_ppa",
    "block_top_ppa_tie_count", "block_top_marker_is_unique",
    "block_credible_set_size", "block_credible_set_ppa_sum",
    "interpretation_note", "block_id", "tissue",
    "feature_id", "feature_primary_id", "feature_label", "feature_class",
    "expression_gene_id", "expression_source", "selected_sample_count",
    "common_marker_count", "total_reads", "chrom", "pos", "marker_id",
    "ref", "alt", "phase_scope", "local_ps_by_sample", "ppa",
    "marker_rank_by_ppa", "credible_set_member", "result_status",
    "result_interpretable",
)
GWAS_MANIFEST_REQUIRED_COLUMNS = (
    "locus_id", "trait", "analysis_scope", "branch", "status", "model_status",
    "gwas_lead_variant_id", "gwas_lead_pos", "gwas_lead_p_lrt",
    "gwas_lead_p_lrt_tier", "gwas_lead_in_harmonized_variants", "variant_count",
    "reference_sample_count", "n_locus_median", "min_observed_p_lrt",
    "converged", "niter", "L", "coverage", "credible_set_count", "max_pip",
    "output_dir", "variant_pip", "credible_sets", "credible_set_members",
)
GWAS_VARIANT_REQUIRED_COLUMNS = (
    "matrix_variant_index", "locus_id", "trait", "analysis_scope", "chrom",
    "pos", "ld_variant_id", "ld_a1", "ld_a2", "beta_for_gwas_allele1",
    "beta_for_ld_a1", "se", "z_for_ld_a1", "p_lrt", "n_obs",
    "ambiguous_at_cg", "harmonization_status", "pip", "pip_rank",
    "credible_sets", "in_credible_set",
)
GWAS_MEMBER_REQUIRED_COLUMNS = (
    "locus_id", "trait", "analysis_scope", "branch", "cs_index", "cs_id",
    "matrix_variant_index", "chrom", "pos", "ld_variant_id", "ld_a1",
    "ld_a2", "pip", "z_for_ld_a1", "p_lrt",
)
GWAS_CS_REQUIRED_COLUMNS = (
    "locus_id", "trait", "analysis_scope", "branch", "cs_index", "cs_id",
    "cs_size", "coverage",
)

OVERLAP_COLUMNS = (
    "result_layer", "normalized_variant_key", "eqtl_marker_id", "chrom", "pos",
    "eqtl_ref", "eqtl_alt", "eqtl_ambiguous_at_cg", "eqtl_block_id", "eqtl_tissue",
    "eqtl_feature_id", "eqtl_feature_primary_id", "eqtl_feature_label",
    "eqtl_feature_class", "eqtl_feature_origin", "eqtl_expression_gene_id",
    "eqtl_expression_source", "eqtl_expression_matrix_scope",
    "eqtl_catalogue_layer", "eqtl_catalogue_role", "eqtl_is_top_variant",
    "eqtl_is_credible_set_member",
    "eqtl_priority_ppa_ge_0.1", "eqtl_ppa", "eqtl_ppa_tier", "eqtl_ppa_status",
    "eqtl_ppa_boundary_adjusted_count", "eqtl_marker_rank_by_ppa",
    "eqtl_block_max_ppa", "eqtl_block_top_ppa_tie_count",
    "eqtl_block_top_marker_is_unique", "eqtl_block_credible_set_size",
    "eqtl_block_credible_set_ppa_sum", "eqtl_input_interpretation_note",
    "eqtl_selected_sample_count", "eqtl_common_marker_count", "eqtl_total_reads",
    "eqtl_phase_scope", "eqtl_local_ps_by_sample", "gwas_analysis_scope",
    "gwas_branch", "gwas_candidate_type", "gwas_locus_id", "gwas_trait",
    "gwas_ld_variant_id", "gwas_ld_a1", "gwas_ld_a2", "gwas_a1_eqtl_relation",
    "gwas_a2_eqtl_relation", "gwas_allele_order_status", "gwas_beta_for_ld_a1",
    "gwas_se", "gwas_z_for_ld_a1", "gwas_p_lrt", "gwas_pip", "gwas_pip_rank",
    "gwas_credible_sets", "gwas_in_credible_set", "gwas_ambiguous_at_cg",
    "gwas_harmonization_status", "gwas_locus_status", "gwas_model_status",
    "gwas_variant_count", "gwas_n_locus_median", "gwas_lead_variant_id",
    "gwas_lead_p_lrt", "allele_match_status", "interpretation_note",
)
MISMATCH_COLUMNS = (
    "eqtl_block_id", "eqtl_marker_id", "eqtl_chrom", "eqtl_pos", "eqtl_ref",
    "eqtl_alt", "gwas_analysis_scope", "gwas_branch", "gwas_candidate_type",
    "gwas_locus_id", "gwas_trait", "gwas_ld_variant_id", "gwas_ld_a1",
    "gwas_ld_a2", "status",
)
AMBIGUOUS_COLUMNS = (
    "eqtl_block_id", "eqtl_marker_id", "chrom", "pos", "ref", "alt",
    "eqtl_feature_id", "eqtl_tissue", "gwas_branch", "status", "reason",
)
DUPLICATE_COLUMNS = (
    "result_layer", "duplicate_key", "eqtl_block_id", "eqtl_marker_id",
    "gwas_analysis_scope", "gwas_branch", "gwas_trait", "gwas_locus_id",
    "gwas_ld_variant_id", "action",
)
TRAIT_SUMMARY_COLUMNS = (
    "result_layer", "gwas_analysis_scope", "gwas_branch", "gwas_trait",
    "overlap_rows", "unique_variants", "unique_eqtl_blocks", "unique_eqtl_features",
    "unique_eqtl_tissues", "unique_gwas_loci", "min_gwas_p_lrt", "max_eqtl_ppa",
    "max_gwas_pip",
)
BRANCH_COMPARISON_COLUMNS = (
    "run_dir", "branch_output_dir", "gwas_branch", "gwas_manifest_rows",
    "gwas_locus_count", "gwas_primary_loci", "gwas_exploratory_loci",
    "gwas_variant_rows", "gwas_credible_rows", "priority_both_overlap_rows",
    "priority_both_unique_variants", "full_both_overlap_rows",
    "full_both_unique_variants", "priority_noncredible_diagnostic_rows",
    "allele_mismatch_qc_rows", "ambiguous_eqtl_priority_qc_rows",
    "duplicate_output_qc_rows",
)

OUTPUT_NAMES = {
    "priority_both": "dog5_expression_ase_ppa_ge_0.1_gwas_both_credible_sets_exact_overlap.tsv",
    "full_both": "dog5_expression_ase_gwas_both_95pct_credible_sets_exact_overlap.tsv",
    "noncredible": "dog5_expression_ase_ppa_ge_0.1_gwas_noncredible_locus_variant_overlap.tsv",
    "mismatch": "dog5_expression_ase_gwas_finemap_allele_mismatch_qc.tsv",
    "ambiguous": "dog5_expression_ase_gwas_finemap_ambiguous_variant_qc.tsv",
    "duplicates": "dog5_expression_ase_gwas_finemap_duplicate_key_qc.tsv",
    "trait_summary": "dog5_expression_ase_gwas_finemap_trait_summary.tsv",
    "summary": "dog5_expression_ase_gwas_finemap_overlap_summary.tsv",
}
BRANCH_COMPARISON_NAME = "dog5_expression_ase_gwas_finemap_branch_comparison.tsv"
BRANCHES = ("main", "strict_no_ambiguous")

AMBIGUOUS_ALLELE_SETS = {("A", "T"), ("C", "G")}
ALLOWED_EQTL_PPA_STATUS = {"finite", "finite_boundary_adjusted"}
ALLOWED_GWAS_BRANCHES = {"main", "strict_no_ambiguous"}
ALLOWED_ANALYSIS_SCOPES = {"primary", "exploratory"}
ALLOWED_HARMONIZATION = {
    "HARMONIZED_NONAMBIGUOUS",
    "HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED",
}
ALLOWED_FEATURE_ORIGINS = {
    "reference_annotation_gene",
    "dog5_manual_new_gene",
}
EXPECTED_EXPRESSION_MATRIX_SCOPE = "reference_annotation_plus_dog5_manual_new_gene"
EXPECTED_FEATURE_ORIGIN_RULE = (
    "gene-=>reference_annotation_gene;DOG5_NEWGENE_=>dog5_manual_new_gene"
)
EXPECTED_CATALOGUE_ROLE_BY_FLAGS = {
    (False, False): "finite_background",
    (True, False): "top_variant",
    (False, True): "credible_set",
    (True, True): "top_and_credible_set",
}


@dataclass(frozen=True)
class InputPaths:
    eqtl_export_dir: Path
    eqtl_priority: Path
    eqtl_credible_set: Path
    eqtl_export_summary: Path
    eqtl_export_completion: Path
    eqtl_final_qc: Path
    eqtl_final_qc_completion: Path
    gwas_run_dir: Path
    gwas_manifest: Path
    gwas_overall_summary: Path
    gwas_branch_summary: Path
    gwas_model_parameters: Path
    gwas_completion: Path


@dataclass(frozen=True)
class ManifestRecord:
    locus_id: str
    trait: str
    analysis_scope: str
    branch: str
    status: str
    model_status: str
    output_dir: Path
    variant_pip: Path
    credible_sets: Path
    credible_set_members: Path
    lead_variant_id: str
    lead_p_lrt: float
    variant_count: int
    n_locus_median: float
    reference_sample_count: int
    credible_set_count: int


@dataclass(frozen=True)
class GwasRecord:
    analysis_scope: str
    branch: str
    candidate_type: str
    locus_id: str
    trait: str
    chrom: str
    pos: int
    marker_id: str
    allele1: str
    allele0: str
    beta_for_ld_a1: float
    se: float
    z_for_ld_a1: float
    p_lrt: float
    pip: float
    pip_rank: int
    credible_sets: str
    in_credible_set: bool
    ambiguous_at_cg: str
    harmonization_status: str
    locus_status: str
    model_status: str
    variant_count: int
    n_locus_median: float
    lead_variant_id: str
    lead_p_lrt: float

    @property
    def position(self) -> Tuple[str, int]:
        return self.chrom, self.pos

    @property
    def variant(self) -> Tuple[str, int, Tuple[str, str]]:
        return self.chrom, self.pos, allele_set(self.allele1, self.allele0)


@dataclass
class LayerStats:
    rows: int = 0
    variants: Set[str] = field(default_factory=set)
    blocks: Set[str] = field(default_factory=set)
    features: Set[str] = field(default_factory=set)
    tissues: Set[str] = field(default_factory=set)
    loci: Set[str] = field(default_factory=set)
    min_p: float | None = None
    max_ppa: float | None = None
    max_pip: float | None = None

    def add(self, row: Mapping[str, object]) -> None:
        self.rows += 1
        self.variants.add(clean(row.get("normalized_variant_key")))
        self.blocks.add(clean(row.get("eqtl_block_id")))
        self.features.add(clean(row.get("eqtl_feature_id")))
        self.tissues.add(clean(row.get("eqtl_tissue")))
        self.loci.add(clean(row.get("gwas_locus_id")))
        p_value = optional_float(row.get("gwas_p_lrt"))
        ppa = optional_float(row.get("eqtl_ppa"))
        pip = optional_float(row.get("gwas_pip"))
        if p_value is not None:
            self.min_p = p_value if self.min_p is None else min(self.min_p, p_value)
        if ppa is not None:
            self.max_ppa = ppa if self.max_ppa is None else max(self.max_ppa, ppa)
        if pip is not None:
            self.max_pip = pip if self.max_pip is None else max(self.max_pip, pip)


def log(message: str) -> None:
    print("[{0} ] {1}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message), flush=True)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def normalize_chrom(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    return text if text.startswith("chr") else "chr" + text


def safe_int(value: object, field_name: str) -> int:
    try:
        return int(float(clean(value)))
    except (TypeError, ValueError):
        raise ValueError("field {0} is not an integer: {1!r}".format(field_name, value))


def safe_float(value: object, field_name: str) -> float:
    try:
        parsed = float(clean(value))
    except (TypeError, ValueError):
        raise ValueError("field {0} is not floating-point number: {1!r}".format(field_name, value))
    if not math.isfinite(parsed):
        raise ValueError("field {0} is not finite value: {1!r}".format(field_name, value))
    return parsed


def optional_float(value: object) -> float | None:
    text = clean(value)
    if not text or text.upper() in {"NA", "N/A", "."}:
        return None
    parsed = float(text)
    return parsed if math.isfinite(parsed) else None


def is_yes(value: object) -> bool:
    return clean(value).lower() in {"1", "yes", "true", "y"}


def allele_set(first: object, second: object) -> Tuple[str, str]:
    alleles = tuple(sorted((clean(first).upper(), clean(second).upper())))
    if not alleles[0] or not alleles[1] or alleles[0] == alleles[1]:
        raise ValueError("invalid of biallelic gene: {0!r}, {1!r}".format(first, second))
    return alleles


def variant_key(chrom: object, pos: object, first: object, second: object) -> Tuple[str, int, Tuple[str, str]]:
    normalized_chrom = normalize_chrom(chrom)
    normalized_pos = safe_int(pos, "pos")
    if not normalized_chrom or normalized_pos <= 0:
        raise ValueError("invalid coordinates: {0!r}:{1!r}".format(chrom, pos))
    return normalized_chrom, normalized_pos, allele_set(first, second)


def format_variant_key(key: Tuple[str, int, Tuple[str, str]]) -> str:
    return "{0}:{1}:{2}/{3}".format(key[0], key[1], key[2][0], key[2][1])


def open_text(path: Path) -> TextIO:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError("input file does not exist or is empty: {0}".format(path))
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def iter_tsv(path: Path, required: Sequence[str]) -> Iterator[Dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], required)
        for row in reader:
            yield dict(row)


def require_columns(path: Path, observed: Sequence[str], required: Sequence[str]) -> None:
    missing = sorted(set(required) - set(observed))
    if missing:
        raise ValueError("{0} missing required fields: {1}".format(path, ", ".join(missing)))


def load_key_values(path: Path) -> Dict[str, str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        raw_rows = [row for row in reader if row]
    if not raw_rows:
        raise ValueError("{0} is empty".format(path))
    first = raw_rows[0]
    has_named_header = len(first) >= 2 and first[1] == "value" and first[0] in {"key", "metric", "parameter"}
    if has_named_header:
        rows = raw_rows[1:]
    elif len(first) >= 2 and first[0] == "overall_status":
        rows = raw_rows
    else:
        raise ValueError("{0} does not use the supported two-column key-value format".format(path))
    result: Dict[str, str] = {}
    for row in rows:
        if len(row) < 2:
            raise ValueError("{0} exists insufficient two columns of key value row: {1!r}".format(path, row))
        key = clean(row[0])
        if not key or key in result:
            raise ValueError("{0} contains empty or duplicate key: {1!r}".format(path, key))
        result[key] = clean(row[1])
    return result


def load_equals_key_values(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if "=" not in line:
                raise ValueError("{0} is not key=value format: {1!r}".format(path, line))
            key, value = line.split("=", 1)
            key = key.strip()
            if not key or key in result:
                raise ValueError("{0} contains empty or duplicate key: {1!r}".format(path, key))
            result[key] = value.strip()
    if not result:
        raise ValueError("{0} is empty".format(path))
    return result


def expected_feature_origin(feature_id: str) -> str:
    if feature_id.startswith("gene-"):
        return "reference_annotation_gene"
    if feature_id.startswith("DOG5_NEWGENE_"):
        return "dog5_manual_new_gene"
    raise ValueError("cannot based on feature_id determine feature_origin: {0}".format(feature_id))


def make_input_paths(args: argparse.Namespace) -> InputPaths:
    eqtl_dir = Path(args.eqtl_export_dir).expanduser()
    final_qc = Path(args.eqtl_final_qc).expanduser()
    gwas_dir = Path(args.gwas_run_dir).expanduser()
    for label, path in (("--eqtl-export-dir", eqtl_dir), ("--eqtl-final-qc", final_qc), ("--gwas-run-dir", gwas_dir)):
        if not path.is_absolute():
            raise ValueError("{0} must be an absolute path".format(label))
    return InputPaths(
        eqtl_export_dir=eqtl_dir.resolve(),
        eqtl_priority=eqtl_dir / "no_gwas_expression_ase_priority_variants_ppa_ge_0.1.tsv.gz",
        eqtl_credible_set=eqtl_dir / "no_gwas_expression_ase_credible_set_variants.tsv.gz",
        eqtl_export_summary=eqtl_dir / "no_gwas_expression_ase_export_summary.tsv",
        eqtl_export_completion=eqtl_dir / "execution_completion.marker",
        eqtl_final_qc=final_qc.resolve(),
        eqtl_final_qc_completion=final_qc.parent.resolve() / "execution_completion.marker",
        gwas_run_dir=gwas_dir.resolve(),
        gwas_manifest=gwas_dir / "locus_run_manifest.tsv",
        gwas_overall_summary=gwas_dir / "overall_summary.tsv",
        gwas_branch_summary=gwas_dir / "branch_summary.tsv",
        gwas_model_parameters=gwas_dir / "model_parameters.tsv",
        gwas_completion=gwas_dir / "COMPLETED.ok",
    )


def validate_eqtl_row(
    row: Mapping[str, str],
    require_priority: bool,
    require_credible: bool,
    expression_matrix_scope: str,
) -> None:
    if clean(row.get("result_status")) != "valid":
        raise ValueError("eQTL catalogue contains non-valid result")
    if clean(row.get("result_interpretable")) != "1":
        raise ValueError("eQTL catalogue exists cannot interpretation result")
    ppa_status = clean(row.get("ppa_status"))
    if ppa_status not in ALLOWED_EQTL_PPA_STATUS:
        raise ValueError("eQTL catalogue contains an unsupported PPA status: {0}".format(ppa_status))
    adjusted = safe_int(row.get("ppa_boundary_adjusted_count"), "ppa_boundary_adjusted_count")
    if ppa_status == "finite" and adjusted != 0:
        raise ValueError("finite PPA recorded boundary adjustment count must be 0")
    if ppa_status == "finite_boundary_adjusted" and adjusted <= 0:
        raise ValueError("boundary adjustment PPA record missing positive of adjust count")
    ppa = safe_float(row.get("ppa"), "ppa")
    if not 0.0 <= ppa <= 1.0:
        raise ValueError("eQTL PPA not in[0,1]within range")
    if require_priority and ppa < 0.1:
        raise ValueError("priority input exists PPA<0.1 rows")
    if require_credible and not is_yes(row.get("is_credible_set_member")):
        raise ValueError("credible-set input contains non-credible set member")
    if is_yes(row.get("is_credible_set_member")) != is_yes(row.get("credible_set_member")):
        raise ValueError("eQTL catalogue of two credible set member field mismatch")
    if not clean(row.get("block_id")) or not clean(row.get("feature_id")):
        raise ValueError("eQTL catalogue contains empty block_id or feature_id")
    if clean(row.get("catalogue_layer")) != "no_gwas_expression_ase_regulatory_variant":
        raise ValueError("eQTL catalogue_layer is inconsistent with the Step 8.6 frozen output definitions")
    top_flag = is_yes(row.get("is_top_variant"))
    credible_flag = is_yes(row.get("is_credible_set_member"))
    expected_role = EXPECTED_CATALOGUE_ROLE_BY_FLAGS[(top_flag, credible_flag)]
    if clean(row.get("catalogue_role")) != expected_role:
        raise ValueError("eQTL catalogue_role is inconsistent with the top/credible marker")
    feature_id = clean(row.get("feature_id"))
    feature_origin = clean(row.get("feature_origin"))
    if feature_origin not in ALLOWED_FEATURE_ORIGINS:
        raise ValueError("eQTL catalogue contains an unsupported feature_origin: {0}".format(feature_origin))
    if feature_origin != expected_feature_origin(feature_id):
        raise ValueError("eQTL feature_origin is inconsistent with the feature_id prefix: {0}".format(feature_id))
    if clean(row.get("expression_matrix_scope")) != expression_matrix_scope:
        raise ValueError("eQTL expression_matrix_scope differs from the Step 8.6 summary")
    common_marker_count = safe_int(row.get("common_marker_count"), "common_marker_count")
    marker_rank = safe_int(row.get("marker_rank_by_ppa"), "marker_rank_by_ppa")
    if common_marker_count <= 0 or not 1 <= marker_rank <= common_marker_count:
        raise ValueError("eQTL marker_rank_by_ppa exceeds block marker count span")
    block_max_ppa = safe_float(row.get("block_max_ppa"), "block_max_ppa")
    if not 0.0 <= block_max_ppa <= 1.0 or ppa > block_max_ppa + 1e-12:
        raise ValueError("eQTL PPA differs from block_max_ppa")
    top_tie_count = safe_int(row.get("block_top_ppa_tie_count"), "block_top_ppa_tie_count")
    top_is_unique = safe_int(row.get("block_top_marker_is_unique"), "block_top_marker_is_unique")
    if top_tie_count <= 0 or top_is_unique not in {0, 1}:
        raise ValueError("eQTL block top tied metadata invalid")
    if top_is_unique != int(top_tie_count == 1):
        raise ValueError("eQTL block_top_marker_is_unique is inconsistent with the tied count")
    credible_set_size = safe_int(row.get("block_credible_set_size"), "block_credible_set_size")
    credible_set_ppa_sum = safe_float(row.get("block_credible_set_ppa_sum"), "block_credible_set_ppa_sum")
    if not 1 <= credible_set_size <= common_marker_count:
        raise ValueError("eQTL block_credible_set_size exceeds block marker count span")
    if not 0.0 < credible_set_ppa_sum <= 1.0 + 1e-9:
        raise ValueError("eQTL block_credible_set_ppa_sum not in valid span")
    variant_key(row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt"))


def validate_frozen_eqtl_inputs(paths: InputPaths) -> Tuple[Dict[str, str], Dict[str, str]]:
    for path in (
        paths.eqtl_priority,
        paths.eqtl_credible_set,
        paths.eqtl_export_summary,
        paths.eqtl_export_completion,
        paths.eqtl_final_qc,
        paths.eqtl_final_qc_completion,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("eQTL input file does not exist or is empty: {0}".format(path))
    export_summary = load_key_values(paths.eqtl_export_summary)
    export_completion = load_equals_key_values(paths.eqtl_export_completion)
    final_qc = load_key_values(paths.eqtl_final_qc)
    final_qc_completion = load_equals_key_values(paths.eqtl_final_qc_completion)
    if final_qc.get("final_qc_status") != "PASS" or final_qc.get("issue_count") != "0":
        raise ValueError("eQTL 8.7 final QC failed")
    if final_qc.get("run_id") != paths.eqtl_final_qc.parent.name:
        raise ValueError("eQTL Step 8.7 directory name differs from final_qc.run_id")
    if final_qc_completion.get("run_id") != final_qc.get("run_id"):
        raise ValueError("eQTL Step 8.7 completion marker differs from final_qc.run_id")
    if final_qc_completion.get("execution_status") != "success":
        raise ValueError("eQTL 8.7 completion marker is not execution_status=success")
    if final_qc_completion.get("final_qc_status") != final_qc.get("final_qc_status"):
        raise ValueError("eQTL Step 8.7 completion marker differs from final_qc_status")
    if export_summary.get("run_id") != paths.eqtl_export_dir.name:
        raise ValueError("eQTL export directory name differs from export_summary.run_id")
    if export_completion.get("run_id") != export_summary.get("run_id"):
        raise ValueError("eQTL Step 8.6 completion marker differs from export_summary.run_id")
    if export_completion.get("execution_status") != "success":
        raise ValueError("eQTL 8.6 completion marker is not execution_status=success")
    if export_summary.get("expression_matrix_scope") != EXPECTED_EXPRESSION_MATRIX_SCOPE:
        raise ValueError("eQTL 8.6 does not use the expression matrix containing reference annotation plus five-dog candidate genes")
    if export_summary.get("feature_origin_rule") != EXPECTED_FEATURE_ORIGIN_RULE:
        raise ValueError("eQTL Step 8.6 feature_origin_rule is inconsistent with the five-dog-specific annotation")
    if Path(clean(final_qc.get("source_export_dir"))).resolve() != paths.eqtl_export_dir:
        raise ValueError("eQTL 8.7 final QC not for current 8.6 export directory")
    if final_qc.get("source_8.6_run_id") != export_summary.get("run_id"):
        raise ValueError("eQTL Step 8.7 final QC does not match the Step 8.6 run_id")
    if final_qc.get("source_8.6_execution_status") != "success":
        raise ValueError("eQTL 8.7 unconfirmed 8.6 execution success")
    if Path(clean(final_qc.get("source_8.6_execution_completion_marker"))).resolve() != paths.eqtl_export_completion:
        raise ValueError("eQTL 8.7 referenced 8.6 completion marker mismatch")
    if final_qc.get("source_8.6_expression_matrix_scope") != EXPECTED_EXPRESSION_MATRIX_SCOPE:
        raise ValueError("eQTL 8.7 recorded expression matrix span mismatch")
    if final_qc.get("source_8.6_feature_origin_rule") != EXPECTED_FEATURE_ORIGIN_RULE:
        raise ValueError("eQTL 8.7 recorded feature_origin_rule mismatch")
    return export_summary, final_qc


def manifest_record(row: Mapping[str, str], run_dir: Path) -> ManifestRecord:
    locus_id = clean(row.get("locus_id"))
    branch = clean(row.get("branch"))
    scope = clean(row.get("analysis_scope"))
    if branch not in ALLOWED_GWAS_BRANCHES or scope not in ALLOWED_ANALYSIS_SCOPES:
        raise ValueError("Manifest contains an unsupported branch or analysis_scope: {0}, {1}".format(branch, scope))
    if clean(row.get("status")) != "PASS" or clean(row.get("model_status")) != "CONVERGED":
        raise ValueError("GWAS manifest contains a locus without PASS status or CONVERGED model status: {0}".format(locus_id))
    output_dir = Path(clean(row.get("output_dir"))).resolve()
    paths = {name: Path(clean(row.get(name))).resolve() for name in ("variant_pip", "credible_sets", "credible_set_members")}
    run_root = run_dir.resolve()
    for label, path in (("output_dir", output_dir), *paths.items()):
        try:
            path.relative_to(run_root)
        except ValueError:
            raise ValueError("manifest of {0} not in frozen GWAS run within directory: {1}".format(label, path))
        if label == "output_dir":
            if not path.is_dir():
                raise FileNotFoundError("manifest referenced output_dir directory does not exist: {0}".format(path))
        elif not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("manifest referenced {0} does not exist or is empty: {1}".format(label, path))
    completion = output_dir / "COMPLETED.ok"
    if not completion.is_file() or completion.stat().st_size == 0:
        raise FileNotFoundError("locus missing COMPLETED.ok: {0}".format(output_dir))
    run_status_path = output_dir / "run_status.tsv"
    if not run_status_path.is_file() or run_status_path.stat().st_size == 0:
        raise FileNotFoundError("locus missing run_status.tsv: {0}".format(output_dir))
    status_fields, status_rows = read_tsv(run_status_path)
    require_columns(run_status_path, status_fields, ("locus_id", "branch", "status", "overall_status", "model_status"))
    if len(status_rows) != 1:
        raise ValueError("locus run_status.tsv must be exactly has one row: {0}".format(output_dir))
    status_row = status_rows[0]
    if (clean(status_row.get("locus_id")), clean(status_row.get("branch"))) != (locus_id, branch):
        raise ValueError("run_status.tsv identity differs from the manifest: {0}".format(locus_id))
    if clean(status_row.get("status")) != "SUCCESS" or clean(status_row.get("model_status")) != "CONVERGED":
        raise ValueError("run_status.tsv is not SUCCESS/CONVERGED: {0}".format(locus_id))
    return ManifestRecord(
        locus_id=locus_id,
        trait=clean(row.get("trait")),
        analysis_scope=scope,
        branch=branch,
        status=clean(row.get("status")),
        model_status=clean(row.get("model_status")),
        output_dir=output_dir,
        variant_pip=paths["variant_pip"],
        credible_sets=paths["credible_sets"],
        credible_set_members=paths["credible_set_members"],
        lead_variant_id=clean(row.get("gwas_lead_variant_id")),
        lead_p_lrt=safe_float(row.get("gwas_lead_p_lrt"), "gwas_lead_p_lrt"),
        variant_count=safe_int(row.get("variant_count"), "variant_count"),
        n_locus_median=safe_float(row.get("n_locus_median"), "n_locus_median"),
        reference_sample_count=safe_int(row.get("reference_sample_count"), "reference_sample_count"),
        credible_set_count=safe_int(row.get("credible_set_count"), "credible_set_count"),
    )


def validate_top_level_gwas(paths: InputPaths, manifest_records: Sequence[ManifestRecord]) -> None:
    for path in (paths.gwas_completion, paths.gwas_manifest, paths.gwas_overall_summary,
                 paths.gwas_branch_summary, paths.gwas_model_parameters):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("GWAS frozen run file does not exist or is empty: {0}".format(path))
    completion = load_key_values(paths.gwas_completion)
    if completion.get("overall_status") != "PASS":
        raise ValueError("GWAS frozen run COMPLETED.ok is not PASS")
    overall = load_key_values(paths.gwas_overall_summary)
    main_loci = {record.locus_id for record in manifest_records if record.branch == "main"}
    strict_loci = {record.locus_id for record in manifest_records if record.branch == "strict_no_ambiguous"}
    if main_loci != strict_loci:
        raise ValueError("Locus sets differ between GWAS main and strict_no_ambiguous")
    expected_task_count = len(manifest_records)
    if safe_int(overall.get("selected_task_count"), "selected_task_count") != expected_task_count:
        raise ValueError("selected_task_count in GWAS overall_summary differs from the manifest")
    if safe_int(overall.get("success_task_count"), "success_task_count") != expected_task_count:
        raise ValueError("GWAS not all manifest task success")
    if safe_int(overall.get("failed_task_count"), "failed_task_count") != 0:
        raise ValueError("GWAS exists failed task")
    if safe_int(overall.get("skipped_task_count"), "skipped_task_count") != 0:
        raise ValueError("GWAS exists skip task, refuse unexpected inclusion of compare")
    params = load_key_values(paths.gwas_model_parameters)
    if params.get("z_input") != "z_for_ld_a1" or params.get("n_input") != "n_locus_median":
        raise ValueError("GWAS model_parameters is inconsistent with the variant_pip field definitions")
    if params.get("coverage") != "0.95":
        raise ValueError("GWAS SuSiE coverage is not 0.95")
    branches = {value.strip() for value in clean(params.get("branches")).split(",") if value.strip()}
    if branches != ALLOWED_GWAS_BRANCHES:
        raise ValueError("GWAS branch is not main and strict_no_ambiguous")
    branch_fields, branch_rows = read_tsv(paths.gwas_branch_summary)
    require_columns(paths.gwas_branch_summary, branch_fields, (
        "branch", "selected_task_count", "success_task_count", "failed_task_count",
        "skipped_task_count", "modeled_variant_count", "credible_set_count",
    ))
    observed_branches: Set[str] = set()
    for row in branch_rows:
        branch = clean(row.get("branch"))
        if branch not in ALLOWED_GWAS_BRANCHES or branch in observed_branches:
            raise ValueError("GWAS branch_summary contains empty, duplicate or unsupported branch: {0}".format(branch))
        observed_branches.add(branch)
        branch_count = sum(record.branch == branch for record in manifest_records)
        if safe_int(row.get("selected_task_count"), "selected_task_count") != branch_count:
            raise ValueError("branch_summary task count differs from the manifest: {0}".format(branch))
        if safe_int(row.get("success_task_count"), "success_task_count") != branch_count:
            raise ValueError("branch_summary success count differs from the manifest: {0}".format(branch))
        if safe_int(row.get("failed_task_count"), "failed_task_count") != 0 or safe_int(row.get("skipped_task_count"), "skipped_task_count") != 0:
            raise ValueError("branch_summary exists failed or skip task: {0}".format(branch))
    if observed_branches != ALLOWED_GWAS_BRANCHES:
        raise ValueError("branch_summary no at the same time contains two GWAS branch")


def split_cs_labels(value: object) -> Set[str]:
    return {label for label in re.split(r"[,;]", clean(value)) if label}


def load_gwas_records(paths: InputPaths) -> Tuple[List[GwasRecord], Dict[str, int], Dict[Tuple[str, str], int]]:
    fields, rows = read_tsv(paths.gwas_manifest)
    require_columns(paths.gwas_manifest, fields, GWAS_MANIFEST_REQUIRED_COLUMNS)
    records_manifest: List[ManifestRecord] = []
    seen_locus_branch: Set[Tuple[str, str]] = set()
    for row in rows:
        record = manifest_record(row, paths.gwas_run_dir)
        identity = (record.locus_id, record.branch)
        if identity in seen_locus_branch:
            raise ValueError("GWAS manifest contains duplicate locus_id+branch: {0}".format(identity))
        seen_locus_branch.add(identity)
        records_manifest.append(record)
    if not records_manifest:
        raise ValueError("GWAS locus_run_manifest is empty")
    validate_top_level_gwas(paths, records_manifest)

    branch_counts: Dict[str, int] = defaultdict(int)
    scope_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    for record in records_manifest:
        branch_counts[record.branch] += 1
        scope_counts[(record.analysis_scope, record.branch)] += 1

    all_records: List[GwasRecord] = []
    expected_member_pairs: Set[Tuple[str, str, str]] = set()
    observed_member_pairs: Set[Tuple[str, str, str]] = set()
    seen_member_pairs: Set[Tuple[str, str, str]] = set()
    for index, manifest in enumerate(records_manifest, 1):
        row_count = 0
        for row in iter_tsv(manifest.variant_pip, GWAS_VARIANT_REQUIRED_COLUMNS):
            row_count += 1
            if clean(row.get("locus_id")) != manifest.locus_id or clean(row.get("trait")) != manifest.trait:
                raise ValueError("variant_pip identity differs from the manifest: {0}".format(manifest.locus_id))
            if clean(row.get("analysis_scope")) != manifest.analysis_scope:
                raise ValueError("variant_pip analysis_scope differs from the manifest: {0}".format(manifest.locus_id))
            key = variant_key(row.get("chrom"), row.get("pos"), row.get("ld_a1"), row.get("ld_a2"))
            ambiguous = clean(row.get("ambiguous_at_cg")).lower()
            if ambiguous not in {"yes", "no"}:
                raise ValueError("ambiguous_at_cg is not yes/no: {0}".format(manifest.locus_id))
            harmonization = clean(row.get("harmonization_status"))
            if harmonization not in ALLOWED_HARMONIZATION:
                raise ValueError("GWAS contains an unsupported harmonization_status: {0}".format(harmonization))
            if manifest.branch == "strict_no_ambiguous" and ambiguous != "no":
                raise ValueError("strict_no_ambiguous still contains ambiguous site: {0}".format(manifest.locus_id))
            if ambiguous == "yes" and harmonization != "HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED":
                raise ValueError("ambiguous site no frequency parse status: {0}".format(manifest.locus_id))
            if ambiguous == "no" and harmonization != "HARMONIZED_NONAMBIGUOUS":
                raise ValueError("unambiguous site status mismatch: {0}".format(manifest.locus_id))
            in_cs = is_yes(row.get("in_credible_set"))
            cs_labels = split_cs_labels(row.get("credible_sets"))
            if in_cs and not cs_labels:
                raise ValueError("in_credible_set=yes but credible_sets is empty: {0}".format(manifest.locus_id))
            for cs_id in cs_labels:
                expected_member_pairs.add((manifest.locus_id + "|" + manifest.branch, clean(row.get("ld_variant_id")), cs_id))
            all_records.append(GwasRecord(
                analysis_scope=manifest.analysis_scope,
                branch=manifest.branch,
                candidate_type="susie_credible_set" if in_cs else "susie_noncredible_locus_variant",
                locus_id=manifest.locus_id,
                trait=manifest.trait,
                chrom=key[0],
                pos=key[1],
                marker_id=clean(row.get("ld_variant_id")),
                allele1=clean(row.get("ld_a1")).upper(),
                allele0=clean(row.get("ld_a2")).upper(),
                beta_for_ld_a1=safe_float(row.get("beta_for_ld_a1"), "beta_for_ld_a1"),
                se=safe_float(row.get("se"), "se"),
                z_for_ld_a1=safe_float(row.get("z_for_ld_a1"), "z_for_ld_a1"),
                p_lrt=safe_float(row.get("p_lrt"), "p_lrt"),
                pip=safe_float(row.get("pip"), "pip"),
                pip_rank=safe_int(row.get("pip_rank"), "pip_rank"),
                credible_sets=clean(row.get("credible_sets")),
                in_credible_set=in_cs,
                ambiguous_at_cg=ambiguous,
                harmonization_status=harmonization,
                locus_status=manifest.status,
                model_status=manifest.model_status,
                variant_count=manifest.variant_count,
                n_locus_median=manifest.n_locus_median,
                lead_variant_id=manifest.lead_variant_id,
                lead_p_lrt=manifest.lead_p_lrt,
            ))
        if row_count != manifest.variant_count:
            raise ValueError("variant_pip row count differs from manifest variant_count: {0}".format(manifest.locus_id))

        member_fields, member_rows = read_tsv(manifest.credible_set_members)
        require_columns(manifest.credible_set_members, member_fields, GWAS_MEMBER_REQUIRED_COLUMNS)
        for row in member_rows:
            if clean(row.get("locus_id")) != manifest.locus_id or clean(row.get("branch")) != manifest.branch:
                raise ValueError("credible_set_members identity differs from the manifest: {0}".format(manifest.locus_id))
            pair = (manifest.locus_id + "|" + manifest.branch, clean(row.get("ld_variant_id")), clean(row.get("cs_id")))
            if pair in seen_member_pairs:
                raise ValueError("credible_set_members contains duplicate member: {0}".format(pair))
            seen_member_pairs.add(pair)
            observed_member_pairs.add(pair)
        if {pair for pair in expected_member_pairs if pair[0] == manifest.locus_id + "|" + manifest.branch} != {
            pair for pair in observed_member_pairs if pair[0] == manifest.locus_id + "|" + manifest.branch
        }:
            raise ValueError("Credible-set members differ between variant_pip and credible_set_members: {0}".format(manifest.locus_id))

        cs_fields, cs_rows = read_tsv(manifest.credible_sets)
        require_columns(manifest.credible_sets, cs_fields, GWAS_CS_REQUIRED_COLUMNS)
        if len(cs_rows) != manifest.credible_set_count:
            raise ValueError("credible_sets row count differs from manifest credible_set_count: {0}".format(manifest.locus_id))
        if index % 25 == 0 or index == len(records_manifest):
            log("[2/6] checked GWAS locus={0}/{1}".format(index, len(records_manifest)))

    if expected_member_pairs != observed_member_pairs:
        raise ValueError("GWAS global credible set member set mismatch")
    metrics = {
        "gwas_locus_count": len({record.locus_id for record in records_manifest}),
        "gwas_manifest_rows": len(records_manifest),
        "gwas_variant_rows": len(all_records),
        "gwas_credible_rows": sum(record.in_credible_set for record in all_records),
        "gwas_primary_loci": len({record.locus_id for record in records_manifest if record.analysis_scope == "primary"}),
        "gwas_exploratory_loci": len({record.locus_id for record in records_manifest if record.analysis_scope == "exploratory"}),
        "gwas_main_variant_rows": sum(1 for record in all_records if record.branch == "main"),
        "gwas_strict_variant_rows": sum(1 for record in all_records if record.branch == "strict_no_ambiguous"),
    }
    for branch in BRANCHES:
        branch_records = [record for record in all_records if record.branch == branch]
        branch_manifest = [record for record in records_manifest if record.branch == branch]
        metrics["gwas_{0}_manifest_rows".format(branch)] = len(branch_manifest)
        metrics["gwas_{0}_locus_count".format(branch)] = len({record.locus_id for record in branch_manifest})
        metrics["gwas_{0}_primary_loci".format(branch)] = len({
            record.locus_id for record in branch_manifest if record.analysis_scope == "primary"
        })
        metrics["gwas_{0}_exploratory_loci".format(branch)] = len({
            record.locus_id for record in branch_manifest if record.analysis_scope == "exploratory"
        })
        metrics["gwas_{0}_variant_rows".format(branch)] = len(branch_records)
        metrics["gwas_{0}_credible_rows".format(branch)] = sum(
            record.in_credible_set for record in branch_records
        )
    return all_records, metrics, scope_counts


def gwas_allele_audit(eqtl_row: Mapping[str, str], gwas: GwasRecord) -> Tuple[str, str, str]:
    ref = clean(eqtl_row.get("ref")).upper()
    alt = clean(eqtl_row.get("alt")).upper()
    if gwas.allele1 == alt and gwas.allele0 == ref:
        return "A1_ALT", "A2_REF", "ld_a1_ALT_ld_a2_REF"
    if gwas.allele1 == ref and gwas.allele0 == alt:
        return "A1_REF", "A2_ALT", "ld_a1_REF_ld_a2_ALT"
    raise ValueError("unordered allele set matching after cannot parse LD A1/A2 orientation")


def build_overlap_row(layer: str, eqtl: Mapping[str, str], gwas: GwasRecord, priority_keys: Set[Tuple[str, str]]) -> Dict[str, object]:
    key = variant_key(eqtl.get("chrom"), eqtl.get("pos"), eqtl.get("ref"), eqtl.get("alt"))
    relation1, relation0, order_status = gwas_allele_audit(eqtl, gwas)
    eqtl_identity = (clean(eqtl.get("block_id")), clean(eqtl.get("marker_id")))
    notes = {
        "priority_both": "expression+ASE PPA>=0.1 and membership in both GWAS and expression+ASE 95% credible sets; not proof of colocalization or causality.",
        "full_both": "Exact overlap of the complete expression+ASE and GWAS 95% credible sets; not proof of colocalization or causality.",
        "noncredible": "Overlap of expression+ASE variants with PPA>=0.1 and non-credible-set variants within a GWAS locus; diagnostic only.",
    }
    if gwas.branch == "strict_no_ambiguous":
        notes["priority_both"] += "; strict_no_ambiguous sensitivity branch"
        notes["full_both"] += "; strict_no_ambiguous sensitivity branch"
        notes["noncredible"] += "; strict_no_ambiguous sensitivity branch"
    elif gwas.ambiguous_at_cg == "yes":
        notes["priority_both"] += "; this A/T or C/G site is frequency-resolved in the main branch"
        notes["full_both"] += "; this A/T or C/G site is frequency-resolved in the main branch"
        notes["noncredible"] += "; this A/T or C/G site is frequency-resolved in the main branch"
    return {
        "result_layer": layer,
        "normalized_variant_key": format_variant_key(key),
        "eqtl_marker_id": clean(eqtl.get("marker_id")),
        "chrom": key[0], "pos": key[1],
        "eqtl_ref": clean(eqtl.get("ref")).upper(), "eqtl_alt": clean(eqtl.get("alt")).upper(),
        "eqtl_ambiguous_at_cg": "yes" if key[2] in AMBIGUOUS_ALLELE_SETS else "no",
        "eqtl_block_id": clean(eqtl.get("block_id")), "eqtl_tissue": clean(eqtl.get("tissue")),
        "eqtl_feature_id": clean(eqtl.get("feature_id")), "eqtl_feature_primary_id": clean(eqtl.get("feature_primary_id")),
        "eqtl_feature_label": clean(eqtl.get("feature_label")), "eqtl_feature_class": clean(eqtl.get("feature_class")),
        "eqtl_feature_origin": clean(eqtl.get("feature_origin")),
        "eqtl_expression_gene_id": clean(eqtl.get("expression_gene_id")),
        "eqtl_expression_source": clean(eqtl.get("expression_source")),
        "eqtl_expression_matrix_scope": clean(eqtl.get("expression_matrix_scope")),
        "eqtl_catalogue_layer": clean(eqtl.get("catalogue_layer")),
        "eqtl_catalogue_role": clean(eqtl.get("catalogue_role")),
        "eqtl_is_top_variant": 1 if is_yes(eqtl.get("is_top_variant")) else 0,
        "eqtl_is_credible_set_member": 1 if is_yes(eqtl.get("is_credible_set_member")) else 0,
        "eqtl_priority_ppa_ge_0.1": 1 if eqtl_identity in priority_keys else 0,
        "eqtl_ppa": safe_float(eqtl.get("ppa"), "ppa"), "eqtl_ppa_tier": clean(eqtl.get("ppa_tier")),
        "eqtl_ppa_status": clean(eqtl.get("ppa_status")),
        "eqtl_ppa_boundary_adjusted_count": safe_int(eqtl.get("ppa_boundary_adjusted_count"), "ppa_boundary_adjusted_count"),
        "eqtl_marker_rank_by_ppa": safe_int(eqtl.get("marker_rank_by_ppa"), "marker_rank_by_ppa"),
        "eqtl_block_max_ppa": safe_float(eqtl.get("block_max_ppa"), "block_max_ppa"),
        "eqtl_block_top_ppa_tie_count": safe_int(
            eqtl.get("block_top_ppa_tie_count"), "block_top_ppa_tie_count"
        ),
        "eqtl_block_top_marker_is_unique": safe_int(
            eqtl.get("block_top_marker_is_unique"), "block_top_marker_is_unique"
        ),
        "eqtl_block_credible_set_size": safe_int(
            eqtl.get("block_credible_set_size"), "block_credible_set_size"
        ),
        "eqtl_block_credible_set_ppa_sum": safe_float(
            eqtl.get("block_credible_set_ppa_sum"), "block_credible_set_ppa_sum"
        ),
        "eqtl_input_interpretation_note": clean(eqtl.get("interpretation_note")),
        "eqtl_selected_sample_count": safe_int(eqtl.get("selected_sample_count"), "selected_sample_count"),
        "eqtl_common_marker_count": safe_int(eqtl.get("common_marker_count"), "common_marker_count"),
        "eqtl_total_reads": safe_int(eqtl.get("total_reads"), "total_reads"),
        "eqtl_phase_scope": clean(eqtl.get("phase_scope")), "eqtl_local_ps_by_sample": clean(eqtl.get("local_ps_by_sample")),
        "gwas_analysis_scope": gwas.analysis_scope, "gwas_branch": gwas.branch,
        "gwas_candidate_type": gwas.candidate_type, "gwas_locus_id": gwas.locus_id, "gwas_trait": gwas.trait,
        "gwas_ld_variant_id": gwas.marker_id, "gwas_ld_a1": gwas.allele1, "gwas_ld_a2": gwas.allele0,
        "gwas_a1_eqtl_relation": relation1, "gwas_a2_eqtl_relation": relation0,
        "gwas_allele_order_status": order_status, "gwas_beta_for_ld_a1": gwas.beta_for_ld_a1,
        "gwas_se": gwas.se, "gwas_z_for_ld_a1": gwas.z_for_ld_a1, "gwas_p_lrt": gwas.p_lrt,
        "gwas_pip": gwas.pip, "gwas_pip_rank": gwas.pip_rank, "gwas_credible_sets": gwas.credible_sets,
        "gwas_in_credible_set": 1 if gwas.in_credible_set else 0, "gwas_ambiguous_at_cg": gwas.ambiguous_at_cg,
        "gwas_harmonization_status": gwas.harmonization_status, "gwas_locus_status": gwas.locus_status,
        "gwas_model_status": gwas.model_status, "gwas_variant_count": gwas.variant_count,
        "gwas_n_locus_median": gwas.n_locus_median, "gwas_lead_variant_id": gwas.lead_variant_id,
        "gwas_lead_p_lrt": gwas.lead_p_lrt, "allele_match_status": "allele_set_match",
        "interpretation_note": notes[layer],
    }


def output_duplicate_key(layer: str, row: Mapping[str, object]) -> Tuple[str, ...]:
    return (
        layer, clean(row.get("eqtl_block_id")), clean(row.get("eqtl_marker_id")),
        clean(row.get("gwas_analysis_scope")), clean(row.get("gwas_branch")),
        clean(row.get("gwas_trait")), clean(row.get("gwas_locus_id")),
        clean(row.get("gwas_ld_variant_id")),
    )


def write_dict(writer: csv.DictWriter, columns: Sequence[str], row: Mapping[str, object]) -> None:
    writer.writerow({column: row.get(column, "") for column in columns})


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare exact shared variants between five-dog expression+ASE results and the current GWAS fine-mapping results.")
    parser.add_argument("--run-dir", type=Path, required=True, help="Absolute output directory for this run.")
    parser.add_argument("--eqtl-export-dir", type=Path, default=DEFAULT_EQTL_EXPORT_DIR, help="Frozen export directory from step 8.6.")
    parser.add_argument("--eqtl-final-qc", type=Path, default=DEFAULT_EQTL_FINAL_QC, help="Final QC summary from step 8.7.")
    parser.add_argument("--gwas-run-dir", type=Path, default=DEFAULT_GWAS_RUN_DIR, help="Completed run directory for 2.gwas-finemap/1.10.")
    args = parser.parse_args()
    if not args.run_dir.is_absolute():
        raise ValueError("--run-dir must be an absolute path")
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    final_paths = {
        branch: {key: run_dir / branch / name for key, name in OUTPUT_NAMES.items()}
        for branch in BRANCHES
    }
    comparison_path = run_dir / BRANCH_COMPARISON_NAME
    existing = [
        str(path)
        for branch_paths in final_paths.values()
        for path in branch_paths.values()
        if path.exists()
    ]
    if comparison_path.exists():
        existing.append(str(comparison_path))
    if existing:
        raise FileExistsError("main output already exists, refuse to overwrite: {0}".format(", ".join(existing)))

    paths = make_input_paths(args)
    log("[1/6] validate 8.6/8.7 expression+ASE frozen input")
    export_summary, final_qc = validate_frozen_eqtl_inputs(paths)
    log("[1/6] eQTL/ASE QC=PASS; export_run_id={0}".format(export_summary.get("run_id")))

    log("[2/6] read GWAS manifest and check each branch/locus of PIP and credible set")
    gwas_records, gwas_metrics, gwas_scope_counts = load_gwas_records(paths)
    gwas_by_position: MutableMapping[Tuple[str, int], List[GwasRecord]] = defaultdict(list)
    gwas_by_variant: MutableMapping[Tuple[str, int, Tuple[str, str]], List[GwasRecord]] = defaultdict(list)
    gwas_cs_by_variant: MutableMapping[Tuple[str, int, Tuple[str, str]], List[GwasRecord]] = defaultdict(list)
    for record in gwas_records:
        gwas_by_position[record.position].append(record)
        gwas_by_variant[record.variant].append(record)
        if record.in_credible_set:
            gwas_cs_by_variant[record.variant].append(record)
    log("[2/6] GWAS rows={0}; credible_rows={1}; primary_loci={2}; exploratory_loci={3}".format(
        gwas_metrics["gwas_variant_rows"], gwas_metrics["gwas_credible_rows"],
        gwas_metrics["gwas_primary_loci"], gwas_metrics["gwas_exploratory_loci"]))

    staging_dir = Path(tempfile.mkdtemp(prefix=".1.1_compare_staging.", dir=run_dir))
    staged_paths = {
        branch: {key: staging_dir / branch / name for key, name in OUTPUT_NAMES.items()}
        for branch in BRANCHES
    }
    for branch in BRANCHES:
        staged_paths[branch]["priority_both"].parent.mkdir(parents=True, exist_ok=True)
    layer_stats: MutableMapping[Tuple[str, str, str, str], LayerStats] = defaultdict(LayerStats)
    overall_stats: MutableMapping[Tuple[str, str], LayerStats] = defaultdict(LayerStats)
    seen_output_keys: MutableMapping[Tuple[str, str], Set[Tuple[str, ...]]] = defaultdict(set)
    priority_keys: Set[Tuple[str, str]] = set()
    priority_row_count = 0
    priority_noncredible_count = 0
    priority_boundary_adjusted_count = 0
    priority_origin_counts: MutableMapping[str, int] = defaultdict(int)
    full_cs_row_count = 0
    full_cs_boundary_adjusted_count = 0
    full_cs_origin_counts: MutableMapping[str, int] = defaultdict(int)
    mismatch_counts: MutableMapping[str, int] = defaultdict(int)
    ambiguous_count = 0
    duplicate_counts: MutableMapping[str, int] = defaultdict(int)

    try:
        with ExitStack() as stack:
            writers: Dict[Tuple[str, str], csv.DictWriter] = {}
            output_columns = {
                "priority_both": OVERLAP_COLUMNS, "full_both": OVERLAP_COLUMNS,
                "noncredible": OVERLAP_COLUMNS, "mismatch": MISMATCH_COLUMNS,
                "ambiguous": AMBIGUOUS_COLUMNS, "duplicates": DUPLICATE_COLUMNS,
            }
            for branch in BRANCHES:
                for key, columns in output_columns.items():
                    handle = stack.enter_context(
                        staged_paths[branch][key].open("w", encoding="utf-8", newline="")
                    )
                    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writers[(branch, key)] = writer

            log("[3/6] streaming process PPA>=0.1 expression+ASE candidate")
            seen_priority_input: Set[Tuple[str, str]] = set()
            seen_mismatch: Set[Tuple[str, ...]] = set()
            seen_ambiguous: Set[Tuple[str, str]] = set()
            for eqtl in iter_tsv(paths.eqtl_priority, EQTL_REQUIRED_COLUMNS):
                priority_row_count += 1
                validate_eqtl_row(
                    eqtl,
                    True,
                    False,
                    export_summary["expression_matrix_scope"],
                )
                priority_origin_counts[clean(eqtl.get("feature_origin"))] += 1
                if clean(eqtl.get("ppa_status")) == "finite_boundary_adjusted":
                    priority_boundary_adjusted_count += 1
                identity = (clean(eqtl.get("block_id")), clean(eqtl.get("marker_id")))
                if identity in seen_priority_input:
                    raise ValueError("eQTL priority contains duplicate block_id+marker_id: {0}".format(identity))
                seen_priority_input.add(identity)
                priority_keys.add(identity)
                if not is_yes(eqtl.get("is_credible_set_member")):
                    priority_noncredible_count += 1
                eqtl_variant = variant_key(eqtl.get("chrom"), eqtl.get("pos"), eqtl.get("ref"), eqtl.get("alt"))
                eqtl_position = (eqtl_variant[0], eqtl_variant[1])
                if eqtl_variant[2] in AMBIGUOUS_ALLELE_SETS and identity not in seen_ambiguous:
                    for branch in BRANCHES:
                        write_dict(writers[(branch, "ambiguous")], AMBIGUOUS_COLUMNS, {
                            "eqtl_block_id": identity[0], "eqtl_marker_id": identity[1],
                            "chrom": eqtl_variant[0], "pos": eqtl_variant[1],
                            "ref": clean(eqtl.get("ref")).upper(), "alt": clean(eqtl.get("alt")).upper(),
                            "eqtl_feature_id": clean(eqtl.get("feature_id")),
                            "eqtl_tissue": clean(eqtl.get("tissue")), "gwas_branch": branch,
                            "status": (
                                "main_branch_retained_for_frequency_resolution"
                                if branch == "main" else "strict_no_ambiguous_excluded"
                            ),
                            "reason": (
                                "The main branch allows frequency resolution of A/T or C/G sites matched to GWAS; unmatched sites do not constitute overlaps."
                                if branch == "main" else
                                "The strict_no_ambiguous branch excludes all ambiguous A/T or C/G sites."
                            ),
                        })
                    seen_ambiguous.add(identity)
                    ambiguous_count += 1
                for gwas in gwas_by_variant.get(eqtl_variant, []):
                    if gwas.in_credible_set and is_yes(eqtl.get("is_credible_set_member")):
                        layer = "priority_both"
                    elif not gwas.in_credible_set:
                        layer = "noncredible"
                    else:
                        continue
                    output_row = build_overlap_row(layer, eqtl, gwas, priority_keys)
                    dedup_key = output_duplicate_key(layer, output_row)
                    branch_layer = (gwas.branch, layer)
                    if dedup_key in seen_output_keys[branch_layer]:
                        write_dict(writers[(gwas.branch, "duplicates")], DUPLICATE_COLUMNS, {
                            "result_layer": layer, "duplicate_key": "|".join(dedup_key),
                            "eqtl_block_id": output_row["eqtl_block_id"], "eqtl_marker_id": output_row["eqtl_marker_id"],
                            "gwas_analysis_scope": output_row["gwas_analysis_scope"], "gwas_branch": output_row["gwas_branch"],
                            "gwas_trait": output_row["gwas_trait"], "gwas_locus_id": output_row["gwas_locus_id"],
                            "gwas_ld_variant_id": output_row["gwas_ld_variant_id"], "action": "duplicate_skipped",
                        })
                        duplicate_counts[gwas.branch] += 1
                        continue
                    seen_output_keys[branch_layer].add(dedup_key)
                    write_dict(writers[(gwas.branch, layer)], OVERLAP_COLUMNS, output_row)
                    layer_stats[(layer, gwas.analysis_scope, gwas.branch, gwas.trait)].add(output_row)
                    overall_stats[(gwas.branch, layer)].add(output_row)
                for gwas in gwas_by_position.get(eqtl_position, []):
                    if gwas.variant == eqtl_variant:
                        continue
                    mismatch_key = (identity[0], identity[1], gwas.analysis_scope, gwas.branch, gwas.locus_id, gwas.marker_id)
                    if mismatch_key in seen_mismatch:
                        continue
                    seen_mismatch.add(mismatch_key)
                    write_dict(writers[(gwas.branch, "mismatch")], MISMATCH_COLUMNS, {
                        "eqtl_block_id": identity[0], "eqtl_marker_id": identity[1],
                        "eqtl_chrom": eqtl_variant[0], "eqtl_pos": eqtl_variant[1],
                        "eqtl_ref": clean(eqtl.get("ref")).upper(), "eqtl_alt": clean(eqtl.get("alt")).upper(),
                        "gwas_analysis_scope": gwas.analysis_scope, "gwas_branch": gwas.branch,
                        "gwas_candidate_type": gwas.candidate_type, "gwas_locus_id": gwas.locus_id,
                        "gwas_trait": gwas.trait, "gwas_ld_variant_id": gwas.marker_id,
                        "gwas_ld_a1": gwas.allele1, "gwas_ld_a2": gwas.allele0,
                        "status": "same_position_allele_set_mismatch",
                    })
                    mismatch_counts[gwas.branch] += 1
            expected_priority_rows = safe_int(export_summary.get("priority_variant_rows_ppa_ge_0.1"), "priority_variant_rows_ppa_ge_0.1")
            if priority_row_count != expected_priority_rows:
                raise ValueError("eQTL priority row count differs from the Step 8.6 summary")
            log("[3/6] priority_rows={0}; main.priority_both={1}; strict.priority_both={2}; "
                "main.noncredible={3}; strict.noncredible={4}".format(
                    priority_row_count,
                    overall_stats[("main", "priority_both")].rows,
                    overall_stats[("strict_no_ambiguous", "priority_both")].rows,
                    overall_stats[("main", "noncredible")].rows,
                    overall_stats[("strict_no_ambiguous", "noncredible")].rows,
                ))

            log("[4/6] streaming process complete expression+ASE credible set")
            for eqtl in iter_tsv(paths.eqtl_credible_set, EQTL_REQUIRED_COLUMNS):
                full_cs_row_count += 1
                validate_eqtl_row(
                    eqtl,
                    False,
                    True,
                    export_summary["expression_matrix_scope"],
                )
                full_cs_origin_counts[clean(eqtl.get("feature_origin"))] += 1
                if clean(eqtl.get("ppa_status")) == "finite_boundary_adjusted":
                    full_cs_boundary_adjusted_count += 1
                eqtl_variant = variant_key(eqtl.get("chrom"), eqtl.get("pos"), eqtl.get("ref"), eqtl.get("alt"))
                for gwas in gwas_cs_by_variant.get(eqtl_variant, []):
                    output_row = build_overlap_row("full_both", eqtl, gwas, priority_keys)
                    dedup_key = output_duplicate_key("full_both", output_row)
                    branch_layer = (gwas.branch, "full_both")
                    if dedup_key in seen_output_keys[branch_layer]:
                        write_dict(writers[(gwas.branch, "duplicates")], DUPLICATE_COLUMNS, {
                            "result_layer": "full_both", "duplicate_key": "|".join(dedup_key),
                            "eqtl_block_id": output_row["eqtl_block_id"], "eqtl_marker_id": output_row["eqtl_marker_id"],
                            "gwas_analysis_scope": output_row["gwas_analysis_scope"], "gwas_branch": output_row["gwas_branch"],
                            "gwas_trait": output_row["gwas_trait"], "gwas_locus_id": output_row["gwas_locus_id"],
                            "gwas_ld_variant_id": output_row["gwas_ld_variant_id"], "action": "duplicate_skipped",
                        })
                        duplicate_counts[gwas.branch] += 1
                        continue
                    seen_output_keys[branch_layer].add(dedup_key)
                    write_dict(writers[(gwas.branch, "full_both")], OVERLAP_COLUMNS, output_row)
                    layer_stats[("full_both", gwas.analysis_scope, gwas.branch, gwas.trait)].add(output_row)
                    overall_stats[(gwas.branch, "full_both")].add(output_row)
                if full_cs_row_count % 1_000_000 == 0:
                    log("[4/6] credible_set_rows_processed={0}; main_exact_rows={1}; strict_exact_rows={2}".format(
                        full_cs_row_count,
                        overall_stats[("main", "full_both")].rows,
                        overall_stats[("strict_no_ambiguous", "full_both")].rows,
                    ))
            expected_cs_rows = safe_int(export_summary.get("credible_set_variant_rows"), "credible_set_variant_rows")
            if full_cs_row_count != expected_cs_rows:
                raise ValueError("eQTL credible-set row count differs from the Step 8.6 summary")
            log("[4/6] credible_set_rows={0}; main_full_both_rows={1}; strict_full_both_rows={2}".format(
                full_cs_row_count,
                overall_stats[("main", "full_both")].rows,
                overall_stats[("strict_no_ambiguous", "full_both")].rows,
            ))

        log("[5/6] write main with strict_no_ambiguous of branch summary and root directory compare summary")
        for branch in BRANCHES:
            with staged_paths[branch]["trait_summary"].open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=TRAIT_SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for (layer, scope, row_branch, trait), stats in sorted(layer_stats.items()):
                    if row_branch != branch:
                        continue
                    write_dict(writer, TRAIT_SUMMARY_COLUMNS, {
                        "result_layer": layer, "gwas_analysis_scope": scope, "gwas_branch": row_branch,
                        "gwas_trait": trait, "overlap_rows": stats.rows, "unique_variants": len(stats.variants),
                        "unique_eqtl_blocks": len(stats.blocks), "unique_eqtl_features": len(stats.features),
                        "unique_eqtl_tissues": len(stats.tissues), "unique_gwas_loci": len(stats.loci),
                        "min_gwas_p_lrt": "" if stats.min_p is None else "{0:.12g}".format(stats.min_p),
                        "max_eqtl_ppa": "" if stats.max_ppa is None else "{0:.12g}".format(stats.max_ppa),
                        "max_gwas_pip": "" if stats.max_pip is None else "{0:.12g}".format(stats.max_pip),
                    })

            summary_rows = [
                ("run_dir", str(run_dir)),
                ("branch_output_dir", str(final_paths[branch]["summary"].parent)),
                ("gwas_branch", branch),
                ("eqtl_priority_input", str(paths.eqtl_priority)),
                ("eqtl_credible_set_input", str(paths.eqtl_credible_set)),
                ("eqtl_export_summary_input", str(paths.eqtl_export_summary)),
                ("eqtl_export_completion_input", str(paths.eqtl_export_completion)),
                ("eqtl_final_qc_input", str(paths.eqtl_final_qc)),
                ("eqtl_final_qc_completion_input", str(paths.eqtl_final_qc_completion)),
                ("gwas_run_dir", str(paths.gwas_run_dir)),
                ("gwas_manifest_input", str(paths.gwas_manifest)),
                ("gwas_overall_summary_input", str(paths.gwas_overall_summary)),
                ("gwas_branch_summary_input", str(paths.gwas_branch_summary)),
                ("gwas_model_parameters_input", str(paths.gwas_model_parameters)),
                ("eqtl_final_qc_status", final_qc.get("final_qc_status", "")),
                ("eqtl_expression_matrix_scope", export_summary.get("expression_matrix_scope", "")),
                ("eqtl_feature_origin_rule", export_summary.get("feature_origin_rule", "")),
                ("eqtl_priority_rows", priority_row_count),
                ("eqtl_priority_reference_annotation_rows", priority_origin_counts["reference_annotation_gene"]),
                ("eqtl_priority_dog5_manual_new_gene_rows", priority_origin_counts["dog5_manual_new_gene"]),
                ("eqtl_priority_noncredible_rows", priority_noncredible_count),
                ("eqtl_priority_boundary_adjusted_rows", priority_boundary_adjusted_count),
                ("eqtl_full_credible_set_rows", full_cs_row_count),
                ("eqtl_full_credible_set_reference_annotation_rows", full_cs_origin_counts["reference_annotation_gene"]),
                ("eqtl_full_credible_set_dog5_manual_new_gene_rows", full_cs_origin_counts["dog5_manual_new_gene"]),
                ("eqtl_full_credible_set_boundary_adjusted_rows", full_cs_boundary_adjusted_count),
                ("gwas_manifest_rows", gwas_metrics["gwas_{0}_manifest_rows".format(branch)]),
                ("gwas_locus_count", gwas_metrics["gwas_{0}_locus_count".format(branch)]),
                ("gwas_primary_loci", gwas_metrics["gwas_{0}_primary_loci".format(branch)]),
                ("gwas_exploratory_loci", gwas_metrics["gwas_{0}_exploratory_loci".format(branch)]),
                ("gwas_variant_rows", gwas_metrics["gwas_{0}_variant_rows".format(branch)]),
                ("gwas_credible_rows", gwas_metrics["gwas_{0}_credible_rows".format(branch)]),
                ("priority_both_credible_overlap_rows", overall_stats[(branch, "priority_both")].rows),
                ("priority_both_credible_unique_variants", len(overall_stats[(branch, "priority_both")].variants)),
                ("full_both_credible_overlap_rows", overall_stats[(branch, "full_both")].rows),
                ("full_both_credible_unique_variants", len(overall_stats[(branch, "full_both")].variants)),
                ("priority_noncredible_diagnostic_rows", overall_stats[(branch, "noncredible")].rows),
                ("allele_mismatch_qc_rows", mismatch_counts[branch]),
                ("ambiguous_eqtl_priority_qc_rows", ambiguous_count),
                ("duplicate_output_qc_rows", duplicate_counts[branch]),
                ("matching_rule", "same_CHR_POS_unordered_allele_set_no_complement_guess"),
                ("branch_policy", "Main results from main; sensitivity results from strict_no_ambiguous; the branches are stored in separate directories and their counts are not combined."),
                ("interpretation_boundary", "Shared fine-map candidates; not formal colocalization, genome-wide significance, or proof of a causal variant or gene."),
            ]


            for scope in ("primary", "exploratory"):
                metric_name = "gwas_{0}_loci".format(scope)
                expected_count = gwas_metrics[
                    "gwas_{0}_{1}_loci".format(branch, scope)
                ]
                observed_count = gwas_scope_counts[(scope, branch)]
                if expected_count != observed_count:
                    raise ValueError(
                        "GWAS stratify locus count mismatch: branch={0}; scope={1}; "
                        "metrics={2}; scope_counts={3}".format(
                            branch,
                            scope,
                            expected_count,
                            observed_count,
                        )
                    )
                if sum(name == metric_name for name, _ in summary_rows) != 1:
                    raise ValueError(
                        "branch summary in metric is not unique value: {0}".format(
                            metric_name
                        )
                    )
            metric_names = [name for name, _ in summary_rows]
            duplicate_metric_names = sorted({
                name for name in metric_names if metric_names.count(name) > 1
            })
            if duplicate_metric_names:
                raise ValueError(
                    "branch summary contains duplicate metric: {0}".format(
                        ",".join(duplicate_metric_names)
                    )
                )
            with staged_paths[branch]["summary"].open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["metric", "value"])
                writer.writerows(summary_rows)

        comparison_staged_path = staging_dir / BRANCH_COMPARISON_NAME
        with comparison_staged_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=BRANCH_COMPARISON_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for branch in BRANCHES:
                write_dict(writer, BRANCH_COMPARISON_COLUMNS, {
                    "run_dir": str(run_dir),
                    "branch_output_dir": str(final_paths[branch]["summary"].parent),
                    "gwas_branch": branch,
                    "gwas_manifest_rows": gwas_metrics["gwas_{0}_manifest_rows".format(branch)],
                    "gwas_locus_count": gwas_metrics["gwas_{0}_locus_count".format(branch)],
                    "gwas_primary_loci": gwas_metrics["gwas_{0}_primary_loci".format(branch)],
                    "gwas_exploratory_loci": gwas_metrics["gwas_{0}_exploratory_loci".format(branch)],
                    "gwas_variant_rows": gwas_metrics["gwas_{0}_variant_rows".format(branch)],
                    "gwas_credible_rows": gwas_metrics["gwas_{0}_credible_rows".format(branch)],
                    "priority_both_overlap_rows": overall_stats[(branch, "priority_both")].rows,
                    "priority_both_unique_variants": len(overall_stats[(branch, "priority_both")].variants),
                    "full_both_overlap_rows": overall_stats[(branch, "full_both")].rows,
                    "full_both_unique_variants": len(overall_stats[(branch, "full_both")].variants),
                    "priority_noncredible_diagnostic_rows": overall_stats[(branch, "noncredible")].rows,
                    "allele_mismatch_qc_rows": mismatch_counts[branch],
                    "ambiguous_eqtl_priority_qc_rows": ambiguous_count,
                    "duplicate_output_qc_rows": duplicate_counts[branch],
                })

        for branch in BRANCHES:
            final_paths[branch]["summary"].parent.mkdir(parents=True, exist_ok=True)
            for key, final_path in final_paths[branch].items():
                os.replace(staged_paths[branch][key], final_path)
        os.replace(comparison_staged_path, comparison_path)
        shutil.rmtree(staging_dir)
    except Exception:
        log("[error] task failed; temporary directory retain used for troubleshoot: {0}".format(staging_dir))
        raise

    log("[6/6] completed / done")
    for branch in BRANCHES:
        log("[output] branch={0}; output_dir={1}".format(branch, final_paths[branch]["summary"].parent))
        for key, path in final_paths[branch].items():
            log("[output] branch={0}; {1}={2}".format(branch, key, path))
    log("[output] branch_comparison={0}".format(comparison_path))
    log("[done] five dogs expression+ASE with current GWAS fine-map exact overlap completed; "
        "main_priority_both_rows={0}; strict_priority_both_rows={1}; "
        "main_full_both_rows={2}; strict_full_both_rows={3}".format(
            overall_stats[("main", "priority_both")].rows,
            overall_stats[("strict_no_ambiguous", "priority_both")].rows,
            overall_stats[("main", "full_both")].rows,
            overall_stats[("strict_no_ambiguous", "full_both")].rows,
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
