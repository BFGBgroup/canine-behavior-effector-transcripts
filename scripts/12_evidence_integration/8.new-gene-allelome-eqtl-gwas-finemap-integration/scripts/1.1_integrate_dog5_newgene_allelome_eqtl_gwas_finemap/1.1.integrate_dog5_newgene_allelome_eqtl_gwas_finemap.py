#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
import traceback
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


DEFAULT_STEP6_LINKS = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740/dog5_allelome_eqtl_finemap_expression_gene_link_overlap.tsv.gz'
)
DEFAULT_STEP6_SUMMARY = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740/dog5_allelome_eqtl_finemap_run_summary.tsv'
)
DEFAULT_STEP3_BOTH_CREDIBLE = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/1.1_compare_dog5_eqtl_gwas_finemap/20260910_131929/main/dog5_expression_ase_ppa_ge_0.1_gwas_both_credible_sets_exact_overlap.tsv'
)
DEFAULT_STEP3_NONCREDIBLE = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/1.1_compare_dog5_eqtl_gwas_finemap/20260910_131929/main/dog5_expression_ase_ppa_ge_0.1_gwas_noncredible_locus_variant_overlap.tsv'
)
DEFAULT_STEP3_SUMMARY = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/1.1_compare_dog5_eqtl_gwas_finemap/20260910_131929/main/dog5_expression_ase_gwas_finemap_overlap_summary.tsv'
)
DEFAULT_GWAS_RUN_DIR = Path(
    '/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap'
)


STEP3_REQUIRED_COLUMNS = {
    "result_layer",
    "normalized_variant_key",
    "eqtl_marker_id",
    "chrom",
    "pos",
    "eqtl_ref",
    "eqtl_alt",
    "eqtl_block_id",
    "eqtl_tissue",
    "eqtl_feature_id",
    "eqtl_feature_primary_id",
    "eqtl_feature_label",
    "eqtl_feature_class",
    "eqtl_feature_origin",
    "eqtl_expression_gene_id",
    "eqtl_expression_source",
    "eqtl_expression_matrix_scope",
    "eqtl_is_credible_set_member",
    "eqtl_priority_ppa_ge_0.1",
    "eqtl_ppa",
    "eqtl_selected_sample_count",
    "eqtl_common_marker_count",
    "eqtl_total_reads",
    "eqtl_phase_scope",
    "gwas_analysis_scope",
    "gwas_branch",
    "gwas_candidate_type",
    "gwas_locus_id",
    "gwas_trait",
    "gwas_ld_variant_id",
    "gwas_ld_a1",
    "gwas_ld_a2",
    "gwas_p_lrt",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_credible_sets",
    "gwas_in_credible_set",
    "gwas_locus_status",
    "gwas_model_status",
    "gwas_lead_variant_id",
    "gwas_lead_p_lrt",
    "allele_match_status",
    "interpretation_note",
}

STEP6_REQUIRED_COLUMNS = {
    "block_id",
    "eqtl_tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "feature_origin",
    "expression_gene_id",
    "expression_source",
    "expression_matrix_scope",
    "eqtl_chrom",
    "selected_sample_count",
    "selected_samples",
    "common_marker_count",
    "total_reads",
    "animal_id",
    "link_tissue",
    "sample_tissue",
    "link_chromosome",
    "canonical_pair_id",
    "sample_pair_id",
    "evidence_stratum",
    "overlap_bp",
    "either_endpoint_multigene",
    "shared_variant_identity_checked",
    "source_link_file",
    "source_link_line",
    "pair_type",
    "interval_gap_bp",
    "mechanism",
    "linkage_score",
    "base_record_type",
    "base_gene_id",
    "base_gene_symbol",
    "base_gene_strand",
    "base_chr",
    "base_start",
    "base_end",
    "base_name",
    "base_gene_biotype",
    "base_endpoint_type",
    "base_manual_new_gene_hit",
    "base_manual_new_gene_ids",
    "base_manual_new_gene_strands",
    "base_manual_new_gene_class_codes",
    "base_manual_new_gene_novelty_categories",
    "base_manual_new_gene_candidate_categories",
    "base_manual_new_gene_dataset_support",
    "base_manual_new_gene_merge_statuses",
    "target_record_type",
    "target_gene_id",
    "target_gene_symbol",
    "target_gene_strand",
    "target_chr",
    "target_start",
    "target_end",
    "target_name",
    "target_gene_biotype",
    "target_endpoint_type",
    "target_manual_new_gene_hit",
    "target_manual_new_gene_ids",
    "target_manual_new_gene_strands",
    "target_manual_new_gene_class_codes",
    "target_manual_new_gene_novelty_categories",
    "target_manual_new_gene_candidate_categories",
    "target_manual_new_gene_dataset_support",
    "target_manual_new_gene_merge_statuses",
    "endpoint_hit",
    "endpoint_mapping_rules",
    "evidence_class",
    "evidence_level",
    "hit_endpoint_has_manual_new_gene",
    "link_contains_manual_new_gene",
    "new_gene_evidence_scope",
}

GWAS_MANIFEST_REQUIRED_COLUMNS = {
    "locus_id", "trait", "analysis_scope", "branch", "status",
    "model_status", "gwas_lead_variant_id", "gwas_lead_pos",
    "gwas_lead_p_lrt", "variant_count", "credible_set_count", "max_pip",
    "top_variant_id", "top_variant_pos", "variant_pip",
}

GWAS_VARIANT_REQUIRED_COLUMNS = {
    "locus_id", "trait", "analysis_scope", "chrom", "pos",
    "ld_variant_id", "ld_a1", "ld_a2", "p_lrt", "pip", "pip_rank",
    "credible_sets", "in_credible_set",
}

GWAS_CONTEXT_COLUMNS = [
    "gwas_overlap_snp_id",
    "gwas_overlap_snp_chrom",
    "gwas_overlap_snp_pos",
    "gwas_overlap_snp_p_lrt",
    "gwas_overlap_snp_pip",
    "gwas_overlap_snp_pip_rank",
    "gwas_overlap_snp_in_credible_set",
    "gwas_locus_lead_snp_id",
    "gwas_locus_lead_chrom",
    "gwas_locus_lead_pos",
    "gwas_locus_lead_allele1",
    "gwas_locus_lead_allele0",
    "gwas_locus_lead_p_lrt",
    "gwas_locus_lead_pip",
    "gwas_locus_lead_pip_rank",
    "gwas_locus_lead_credible_sets",
    "gwas_locus_lead_in_credible_set",
    "overlap_snp_equals_lead",
    "gwas_locus_top_pip_snp_id",
    "gwas_locus_top_pip_chrom",
    "gwas_locus_top_pip_pos",
    "gwas_locus_top_pip_allele1",
    "gwas_locus_top_pip_allele0",
    "gwas_locus_top_pip_p_lrt",
    "gwas_locus_top_pip",
    "gwas_locus_top_pip_rank",
    "gwas_locus_top_pip_credible_sets",
    "gwas_locus_top_pip_in_credible_set",
    "overlap_snp_equals_top_pip",
]

CANDIDATE_DERIVED_COLUMNS = [
    "gwas_finemap_support_class",
    "gwas_locus_lead_variant_id",
    "shared_snp_equals_gwas_locus_lead",
    *GWAS_CONTEXT_COLUMNS,
    "allelome_block_match",
    "multilayer_integration_scope",
    "exact_shared_snp_allelome_endpoint_rows",
    "shared_snp_equals_any_allelome_endpoint",
    "shared_snp_min_distance_to_allelome_endpoint_bp",
    "expanded_allelome_link_rows",
    "allelome_link_count",
    "supporting_dog_count",
    "supporting_dogs",
    "pcgene_link_count",
    "partner_gene_count",
    "partner_gene_ids",
    "partner_gene_symbols",
    "mechanisms",
    "endpoint_manual_new_gene_ids",
    "endpoint_manual_new_gene_strands",
    "endpoint_manual_new_gene_class_codes",
    "endpoint_manual_new_gene_novelty_categories",
    "endpoint_manual_new_gene_candidate_categories",
    "endpoint_manual_new_gene_dataset_support",
    "endpoint_manual_new_gene_merge_statuses",
    "evidence_tier",
    "evidence_tier_order",
    "main_result",
    "evidence_interpretation",
    "allelome_evidence_strata",
    "allelome_multigene_link_rows",
    "allelome_overlap_link_rows",
]

LINK_PREFIX_COLUMNS = [
    "candidate_id",
    "new_gene_id",
    "evidence_tier",
    "evidence_tier_order",
    "main_result",
    "gwas_finemap_support_class",
    "gwas_locus_lead_variant_id",
    "shared_snp_equals_gwas_locus_lead",
    *GWAS_CONTEXT_COLUMNS,
    "multilayer_integration_scope",
    "evidence_interpretation",
]

LINK_DERIVED_COLUMNS = [
    "new_gene_endpoint_role",
    "new_gene_endpoint_record_types",
    "new_gene_endpoint_names",
    "new_gene_endpoint_chromosomes",
    "new_gene_endpoint_starts0",
    "new_gene_endpoint_ends0",
    "endpoint_coordinate_system",
    "new_gene_endpoint_manual_ids",
    "new_gene_endpoint_manual_strands",
    "new_gene_endpoint_manual_class_codes",
    "new_gene_endpoint_manual_novelty_categories",
    "new_gene_endpoint_manual_candidate_categories",
    "new_gene_endpoint_manual_dataset_support",
    "new_gene_endpoint_manual_merge_statuses",
    "allelome_endpoint_variant_keys",
    "shared_finemap_snp_equals_allelome_endpoint",
    "shared_finemap_snp_min_distance_to_allelome_endpoint_bp",
    "partner_role",
    "partner_record_type",
    "partner_gene_id",
    "partner_gene_symbol",
    "partner_gene_strand",
    "partner_chrom",
    "partner_start",
    "partner_end",
    "partner_name",
    "partner_biotype",
    "partner_endpoint_type",
    "candidate_orientation",
    "direction_basis",
]

PAIR_COLUMNS = [
    "evidence_tier",
    "evidence_tier_order",
    "main_result",
    "candidate_id",
    "new_gene_id",
    "eqtl_block_id",
    "eqtl_tissue",
    "chrom",
    "normalized_variant_key",
    "gwas_trait",
    "gwas_locus_id",
    "gwas_ld_variant_id",
    "gwas_p_lrt",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_locus_lead_variant_id",
    "shared_snp_equals_gwas_locus_lead",
    *GWAS_CONTEXT_COLUMNS,
    "multilayer_integration_scope",
    "shared_snp_equals_any_allelome_endpoint",
    "shared_snp_min_distance_to_allelome_endpoint_bp",
    "partner_gene_id",
    "partner_gene_symbol",
    "mechanism",
    "expanded_link_rows",
    "supporting_dog_count",
    "supporting_dogs",
    "allelome_endpoint_count",
    "allelome_endpoints",
    "linkage_score_min",
    "linkage_score_max",
    "interval_gap_bp_min",
    "interval_gap_bp_max",
    "candidate_orientation",
    "direction_interpretation",
    "allelome_evidence_strata",
    "allelome_multigene_link_rows",
    "allelome_overlap_link_rows",
]

TIER_INFO = {
    ("priority_both", "yes"): (
        "A_BOTH_95PCT_CREDIBLE_PLUS_ALLELOME",
        1,
        "yes",
        "Exact sharing of eQTL and GWAS 95% credible-set variants, with an Allelome relationship in the same new-gene/eQTL block; primary evidence-tier candidates. Allelome endpoints are gene intervals, so three-layer overlap at the same SNP cannot be established; review flags must still inform interpretation.",
    ),
    ("noncredible", "yes"): (
        "B_GWAS_NONCREDIBLE_PLUS_ALLELOME",
        2,
        "no",
        "An eQTL high-PPA SNP matches a non-credible-set variant within a GWAS locus, with an Allelome relationship in the same new-gene/eQTL block; diagnostic only.",
    ),
    ("priority_both", "no"): (
        "C_BOTH_95PCT_CREDIBLE_NO_ALLELOME",
        3,
        "no",
        "Exact sharing of eQTL and GWAS 95% credible-set variants, without an Allelome relationship in the same current new-gene/eQTL block.",
    ),
    ("noncredible", "no"): (
        "D_GWAS_NONCREDIBLE_NO_ALLELOME",
        4,
        "no",
        "Diagnostic overlap with a GWAS non-credible-set variant, without an Allelome relationship in the same current new-gene/eQTL block.",
    ),
}

NEW_GENE_RE = re.compile(r"DOG5_NEWGENE_\d+")
LOG_HANDLE = None


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def log(message: str) -> None:
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S} ] {message}"
    print(line, flush=True)
    if LOG_HANDLE is not None:
        LOG_HANDLE.write(line + "\n")
        LOG_HANDLE.flush()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="", encoding="utf-8")
    return path.open(mode, newline="", encoding="utf-8")


def require_columns(
    path: Path,
    observed: Sequence[str],
    required: Iterable[str],
) -> None:
    missing = sorted(set(required) - set(observed))
    if missing:
        raise ValueError(f"{path} missing required fields: {','.join(missing)}")


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, [dict(row) for row in reader]


def read_key_value_contract(path: Path, key_field: str) -> Dict[str, str]:

    fields, rows = read_tsv(path)
    require_columns(path, fields, {key_field, "value"})
    values: Dict[str, str] = {}
    for row in rows:
        key = clean(row.get(key_field))
        value = clean(row.get("value"))
        if not key:
            raise ValueError(f"{path} contains empty {key_field}")
        if key in values:
            raise ValueError(
                f"{path} contains duplicate {key_field} or conflict value: {key!r}; "
                f"first={values[key]!r}; current={value!r}"
            )
        values[key] = value
    return values


def as_int(value: object, label: str) -> int:
    text = clean(value)
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"{label} is not an integer: {text!r}") from exc


def as_float(value: object, label: str) -> float:
    text = clean(value)
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(f"{label} is not numeric: {text!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} is not finite numeric value: {text!r}")
    return number


def split_values(*values: object) -> List[str]:
    result: List[str] = []
    seen: Set[str] = set()
    for value in values:
        for item in clean(value).split("|"):
            item = item.strip()
            if not item or item in {"NA", "."} or item in seen:
                continue
            seen.add(item)
            result.append(item)
    return result


def join_values(values: Iterable[object]) -> str:
    cleaned = sorted(
        {
            clean(value)
            for value in values
            if clean(value) and clean(value) not in {"NA", "."}
        }
    )
    return "|".join(cleaned) if cleaned else "NA"


def join_split_values(rows: Sequence[Mapping[str, object]], field: str) -> str:
    return join_values(
        item
        for row in rows
        for item in split_values(row.get(field))
    )


def canonical_variant(
    chrom: object,
    pos: object,
    allele_a: object,
    allele_b: object,
) -> str:
    contig = clean(chrom)
    position = as_int(pos, "variant position")
    alleles = sorted([clean(allele_a).upper(), clean(allele_b).upper()])
    if not contig or any(not allele for allele in alleles):
        raise ValueError("cannot construct standard variant key: chromosome or allele is empty")
    return f"{contig}:{position}:{alleles[0]}/{alleles[1]}"


def normalize_chrom(value: object) -> str:

    text = clean(value)
    if not text:
        return ""
    return text if text.startswith("chr") else f"chr{text}"


def build_gwas_main_indexes(
    rows: Sequence[Mapping[str, str]],
) -> Tuple[
    Dict[str, Mapping[str, str]],
    Dict[Tuple[str, str], Mapping[str, str]],
]:

    manifest_index: Dict[str, Mapping[str, str]] = {}
    variant_index: Dict[Tuple[str, str], Mapping[str, str]] = {}
    for row in rows:
        if clean(row.get("branch")) != "main":
            continue
        locus_id = clean(row.get("locus_id"))
        if not locus_id or locus_id in manifest_index:
            raise ValueError(f"GWAS main manifest contains empty or duplicate locus: {locus_id!r}")
        if clean(row.get("status")) != "PASS":
            raise ValueError(f"GWAS main locus is not PASS: {locus_id}")
        if clean(row.get("model_status")) != "CONVERGED":
            raise ValueError(f"GWAS main locus not converged: {locus_id}")
        variant_path = Path(clean(row.get("variant_pip")))
        if not variant_path.is_absolute() or not variant_path.is_file():
            raise FileNotFoundError(
                f"GWAS main locus of variant_pip does not exist: {variant_path}"
            )
        variant_fields, variant_rows = read_tsv(variant_path)
        require_columns(
            variant_path,
            variant_fields,
            GWAS_VARIANT_REQUIRED_COLUMNS,
        )
        if len(variant_rows) != as_int(row.get("variant_count"), "variant_count"):
            raise ValueError(f"GWAS variant_pip row count and manifest mismatch: {locus_id}")
        manifest_index[locus_id] = row
        for variant in variant_rows:
            if clean(variant.get("locus_id")) != locus_id:
                raise ValueError(f"GWAS variant_pip locus identity mismatch: {locus_id}")
            if clean(variant.get("trait")) != clean(row.get("trait")):
                raise ValueError(f"GWAS variant_pip trait mismatch: {locus_id}")
            if clean(variant.get("analysis_scope")) != clean(
                row.get("analysis_scope")
            ):
                raise ValueError(f"GWAS variant_pip analysis_scope mismatch: {locus_id}")
            key = (locus_id, clean(variant.get("ld_variant_id")))
            if not key[1] or key in variant_index:
                raise ValueError(f"GWAS main contains empty or duplicate locus/variant key: {key}")
            variant_index[key] = variant
    if not manifest_index:
        raise ValueError("GWAS manifest contains no main branch locus")
    return manifest_index, variant_index


def build_gwas_context(
    candidate: Mapping[str, object],
    manifest_index: Mapping[str, Mapping[str, str]],
    variant_index: Mapping[Tuple[str, str], Mapping[str, str]],
) -> Dict[str, object]:
    locus_id = clean(candidate.get("gwas_locus_id"))
    overlap_id = clean(candidate.get("gwas_ld_variant_id"))
    manifest = manifest_index.get(locus_id)
    if manifest is None:
        raise ValueError(f"GWAS main manifest missing Step 3 locus: {locus_id}")
    lead_id = clean(manifest.get("gwas_lead_variant_id"))
    top_id = clean(manifest.get("top_variant_id"))
    overlap = variant_index.get((locus_id, overlap_id))
    lead = variant_index.get((locus_id, lead_id))
    top = variant_index.get((locus_id, top_id))
    if overlap is None:
        raise ValueError(
            f"GWAS main variant_pip missing overlap SNP: locus={locus_id}, snp={overlap_id}"
        )
    if lead is None:
        raise ValueError(
            f"GWAS main variant_pip missing association lead: locus={locus_id}, snp={lead_id}"
        )
    if top is None:
        raise ValueError(
            f"GWAS main variant_pip missing top-PIP SNP: locus={locus_id}, snp={top_id}"
        )

    for label, observed, expected in [
        ("trait", overlap.get("trait"), candidate.get("gwas_trait")),
        ("analysis_scope", overlap.get("analysis_scope"), candidate.get("gwas_analysis_scope")),
        ("ld_a1", overlap.get("ld_a1"), candidate.get("gwas_ld_a1")),
        ("ld_a2", overlap.get("ld_a2"), candidate.get("gwas_ld_a2")),
    ]:
        if clean(observed) != clean(expected):
            raise ValueError(
                f"Step 3 overlap SNP with GWAS main variant_pip of {label} mismatch: "
                f"{clean(observed)!r} != {clean(expected)!r}"
            )
    if normalize_chrom(overlap.get("chrom")) != clean(candidate.get("chrom")):
        raise ValueError("Step 3 overlap SNP with GWAS main variant_pip of chromosome mismatch")
    if as_int(overlap.get("pos"), "GWAS overlap pos") != as_int(
        candidate.get("pos"), "Step 3 overlap pos"
    ):
        raise ValueError("Step 3 overlap SNP with GWAS main variant_pip of position mismatch")
    for label, observed, expected in [
        ("p_lrt", overlap.get("p_lrt"), candidate.get("gwas_p_lrt")),
        ("pip", overlap.get("pip"), candidate.get("gwas_pip")),
    ]:
        if not math.isclose(
            as_float(observed, f"GWAS overlap {label}"),
            as_float(expected, f"Step 3 overlap {label}"),
            rel_tol=1e-12,
            abs_tol=1e-15,
        ):
            raise ValueError(
                f"Step 3 overlap SNP with GWAS main variant_pip of {label} mismatch"
            )
    if as_int(overlap.get("pip_rank"), "GWAS overlap pip_rank") != as_int(
        candidate.get("gwas_pip_rank"), "Step 3 overlap pip_rank"
    ):
        raise ValueError("Step 3 overlap SNP with GWAS main variant_pip of pip_rank mismatch")
    expected_credible = (
        "yes" if clean(candidate.get("gwas_in_credible_set")) == "1" else "no"
    )
    if clean(overlap.get("in_credible_set")) != expected_credible:
        raise ValueError(
            "Step 3 overlap SNP with GWAS main variant_pip of credible set status mismatch"
        )
    if clean(overlap.get("credible_sets")) != clean(
        candidate.get("gwas_credible_sets")
    ):
        raise ValueError("Step 3 overlap SNP with GWAS main of credible_sets mismatch")
    if clean(candidate.get("gwas_branch")) != "main":
        raise ValueError("Step 3 candidate is not main branch")
    if clean(manifest.get("trait")) != clean(candidate.get("gwas_trait")):
        raise ValueError("GWAS manifest with Step 3 candidate of trait mismatch")
    if clean(manifest.get("analysis_scope")) != clean(
        candidate.get("gwas_analysis_scope")
    ):
        raise ValueError("GWAS manifest with Step 3 candidate of analysis_scope mismatch")
    if clean(manifest.get("status")) != clean(candidate.get("gwas_locus_status")):
        raise ValueError("GWAS manifest with Step 3 candidate of locus status mismatch")
    if clean(manifest.get("model_status")) != clean(
        candidate.get("gwas_model_status")
    ):
        raise ValueError("GWAS manifest with Step 3 candidate of model status mismatch")
    if clean(candidate.get("gwas_lead_variant_id")) != lead_id:
        raise ValueError("Step 3 explicit lead SNP with GWAS manifest mismatch")
    if as_int(lead.get("pos"), "GWAS lead pos") != as_int(
        manifest.get("gwas_lead_pos"), "manifest lead pos"
    ):
        raise ValueError("GWAS association lead position and manifest mismatch")
    if not math.isclose(
        as_float(candidate.get("gwas_lead_p_lrt"), "Step 3 lead p_lrt"),
        as_float(manifest.get("gwas_lead_p_lrt"), "manifest lead p_lrt"),
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        raise ValueError("Step 3 explicit lead P value and GWAS manifest mismatch")
    if not math.isclose(
        as_float(lead.get("p_lrt"), "GWAS lead variant p_lrt"),
        as_float(manifest.get("gwas_lead_p_lrt"), "manifest lead p_lrt"),
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        raise ValueError("GWAS association lead P value and manifest mismatch")
    if as_int(top.get("pos"), "GWAS top-PIP pos") != as_int(
        manifest.get("top_variant_pos"), "manifest top-PIP pos"
    ):
        raise ValueError("GWAS top-PIP variant position and manifest mismatch")
    if not math.isclose(
        as_float(top.get("pip"), "GWAS top-PIP"),
        as_float(manifest.get("max_pip"), "manifest max_pip"),
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        raise ValueError("GWAS top-PIP variant and manifest max_pip mismatch")

    return {
        "gwas_overlap_snp_id": overlap_id,
        "gwas_overlap_snp_chrom": normalize_chrom(overlap.get("chrom")),
        "gwas_overlap_snp_pos": as_int(overlap.get("pos"), "GWAS overlap pos"),
        "gwas_overlap_snp_p_lrt": clean(overlap.get("p_lrt")),
        "gwas_overlap_snp_pip": clean(overlap.get("pip")),
        "gwas_overlap_snp_pip_rank": as_int(
            overlap.get("pip_rank"), "GWAS overlap pip_rank"
        ),
        "gwas_overlap_snp_in_credible_set": clean(
            overlap.get("in_credible_set")
        ),
        "gwas_locus_lead_snp_id": lead_id,
        "gwas_locus_lead_chrom": normalize_chrom(lead.get("chrom")),
        "gwas_locus_lead_pos": as_int(lead.get("pos"), "GWAS lead pos"),
        "gwas_locus_lead_allele1": clean(lead.get("ld_a1")),
        "gwas_locus_lead_allele0": clean(lead.get("ld_a2")),
        "gwas_locus_lead_p_lrt": clean(lead.get("p_lrt")),
        "gwas_locus_lead_pip": clean(lead.get("pip")),
        "gwas_locus_lead_pip_rank": as_int(
            lead.get("pip_rank"), "GWAS lead pip_rank"
        ),
        "gwas_locus_lead_credible_sets": clean(lead.get("credible_sets"))
        or "NA",
        "gwas_locus_lead_in_credible_set": clean(
            lead.get("in_credible_set")
        ),
        "overlap_snp_equals_lead": "yes" if overlap_id == lead_id else "no",
        "gwas_locus_top_pip_snp_id": top_id,
        "gwas_locus_top_pip_chrom": normalize_chrom(top.get("chrom")),
        "gwas_locus_top_pip_pos": as_int(top.get("pos"), "GWAS top-PIP pos"),
        "gwas_locus_top_pip_allele1": clean(top.get("ld_a1")),
        "gwas_locus_top_pip_allele0": clean(top.get("ld_a2")),
        "gwas_locus_top_pip_p_lrt": clean(top.get("p_lrt")),
        "gwas_locus_top_pip": clean(top.get("pip")),
        "gwas_locus_top_pip_rank": as_int(top.get("pip_rank"), "GWAS top-PIP rank"),
        "gwas_locus_top_pip_credible_sets": clean(top.get("credible_sets")) or "NA",
        "gwas_locus_top_pip_in_credible_set": clean(top.get("in_credible_set")),
        "overlap_snp_equals_top_pip": "yes" if overlap_id == top_id else "no",
    }


def candidate_id(row: Mapping[str, str]) -> str:
    return "::".join(
        [
            clean(row.get("result_layer")),
            clean(row.get("gwas_analysis_scope")),
            clean(row.get("gwas_branch")),
            clean(row.get("eqtl_block_id")),
            clean(row.get("normalized_variant_key")),
            clean(row.get("gwas_locus_id")),
            clean(row.get("gwas_ld_variant_id")),
        ]
    )


def block_key(row: Mapping[str, str]) -> Tuple[str, str, str, str]:
    return (
        clean(row.get("eqtl_block_id")),
        clean(row.get("eqtl_expression_gene_id")),
        clean(row.get("eqtl_tissue")),
        clean(row.get("chrom")),
    )


def step6_block_key(row: Mapping[str, str]) -> Tuple[str, str, str, str]:
    return (
        clean(row.get("block_id")),
        clean(row.get("expression_gene_id")),
        clean(row.get("eqtl_tissue")),
        clean(row.get("eqtl_chrom")),
    )


def validate_step3_row(
    row: Mapping[str, str],
    expected_layer: str,
) -> None:
    if clean(row.get("result_layer")) != expected_layer:
        raise ValueError("Step 3 result level label with input file mismatch")
    if clean(row.get("eqtl_priority_ppa_ge_0.1")) != "1":
        raise ValueError("Step 3 prioritize result in encountered eQTL PPA<0.1 record")
    if clean(row.get("eqtl_is_credible_set_member")) != "1":
        raise ValueError("Step 3 results contain variants outside the eQTL credible sets")
    if clean(row.get("allele_match_status")) != "allele_set_match":
        raise ValueError("Step 3 result encountered allele set does not match")
    if clean(row.get("gwas_branch")) != "main":
        raise ValueError("Step 3 results contain a GWAS branch other than main")
    if clean(row.get("gwas_analysis_scope")) not in {"primary", "exploratory"}:
        raise ValueError("Step 3 GWAS analysis_scope is not primary/exploratory")
    if clean(row.get("gwas_locus_status")) != "PASS":
        raise ValueError("Step 3 results contain a non-PASS GWAS locus")
    if clean(row.get("gwas_model_status")) != "CONVERGED":
        raise ValueError("Step 3 result encountered not converged GWAS model")
    if clean(row.get("eqtl_phase_scope")) != clean(row.get("chrom")):
        raise ValueError("Step 3 results do not cover the whole chromosome phase_scope")

    ppa = as_float(row.get("eqtl_ppa"), "Step 3 eqtl_ppa")
    pip = as_float(row.get("gwas_pip"), "Step 3 gwas_pip")
    if not 0.1 <= ppa <= 1.0:
        raise ValueError(f"Step 3 eqtl_ppa not in[0.1,1]span: {ppa}")
    if not 0.0 <= pip <= 1.0:
        raise ValueError(f"Step 3 gwas_pip not in[0,1]span: {pip}")
    if as_int(row.get("gwas_pip_rank"), "Step 3 gwas_pip_rank") < 1:
        raise ValueError("Step 3 gwas_pip_rank less than 1")

    observed_variant = canonical_variant(
        row.get("chrom"),
        row.get("pos"),
        row.get("eqtl_ref"),
        row.get("eqtl_alt"),
    )
    if clean(row.get("normalized_variant_key")) != observed_variant:
        raise ValueError(
            "Step 3 normalized_variant_key with CHROM/POS/REF/ALT mismatch: "
            f"{clean(row.get('normalized_variant_key'))} != {observed_variant}"
        )

    if expected_layer == "priority_both":
        if clean(row.get("gwas_candidate_type")) != "susie_credible_set":
            raise ValueError("priority_both record is not GWAS credible set member")
        if clean(row.get("gwas_in_credible_set")) != "1":
            raise ValueError("priority_both recorded gwas_in_credible_set is not 1")
        if not clean(row.get("gwas_credible_sets")):
            raise ValueError("priority_both record missing gwas_credible_sets")
    else:
        if clean(row.get("gwas_candidate_type")) != (
            "susie_noncredible_locus_variant"
        ):
            raise ValueError("noncredible recorded GWAS incorrect candidate type")
        if clean(row.get("gwas_in_credible_set")) != "0":
            raise ValueError("noncredible recorded gwas_in_credible_set is not 0")

    if not clean(row.get("gwas_lead_variant_id")):
        raise ValueError("Step 3 result missing explicit gwas_lead_variant_id")


def is_new_gene_candidate(row: Mapping[str, str]) -> bool:
    gene_id = clean(row.get("eqtl_expression_gene_id"))
    if NEW_GENE_RE.fullmatch(gene_id) is None:
        return False
    if clean(row.get("eqtl_expression_source")) != (
        "manual_new_gene_log2_sizefactor_norm_counts"
    ):
        raise ValueError(f"{gene_id} of expression_source is not Manual new gene source")
    if clean(row.get("eqtl_feature_origin")) != "dog5_manual_new_gene":
        raise ValueError(f"{gene_id} of feature_origin is not five-dog-specific new gene")
    if clean(row.get("eqtl_expression_matrix_scope")) != (
        "reference_annotation_plus_dog5_manual_new_gene"
    ):
        raise ValueError(f"{gene_id} of expression_matrix_scope is not five-dog-specific matrix")
    if clean(row.get("eqtl_feature_id")) != gene_id:
        raise ValueError(f"{gene_id} of eqtl_feature_id with expression_gene_id mismatch")
    if clean(row.get("eqtl_feature_primary_id")) != gene_id:
        raise ValueError(
            f"{gene_id} of eqtl_feature_primary_id with expression_gene_id mismatch"
        )
    return True


def validate_step6_link(
    row: Mapping[str, str],
    candidate: Mapping[str, str],
) -> None:
    string_pairs = [
        ("block_id", "eqtl_block_id"),
        ("eqtl_tissue", "eqtl_tissue"),
        ("expression_gene_id", "eqtl_expression_gene_id"),
        ("expression_source", "eqtl_expression_source"),
        ("eqtl_chrom", "chrom"),
        ("feature_id", "eqtl_feature_id"),
        ("feature_primary_id", "eqtl_feature_primary_id"),
        ("feature_label", "eqtl_feature_label"),
        ("feature_class", "eqtl_feature_class"),
        ("feature_origin", "eqtl_feature_origin"),
        ("expression_matrix_scope", "eqtl_expression_matrix_scope"),
    ]
    for step6_field, step3_field in string_pairs:
        if clean(row.get(step6_field)) != clean(candidate.get(step3_field)):
            raise ValueError(
                f"Step 6 and Step 3 block field mismatch: {step6_field}/{step3_field}"
            )

    integer_pairs = [
        ("selected_sample_count", "eqtl_selected_sample_count"),
        ("common_marker_count", "eqtl_common_marker_count"),
        ("total_reads", "eqtl_total_reads"),
    ]
    for step6_field, step3_field in integer_pairs:
        if as_int(row.get(step6_field), step6_field) != as_int(
            candidate.get(step3_field),
            step3_field,
        ):
            raise ValueError(
                f"Step 6 and Step 3 block count field mismatch: {step6_field}/{step3_field}"
            )

    if clean(row.get("link_tissue")) != clean(candidate.get("eqtl_tissue")):
        raise ValueError("Step 6 link_tissue with Step 3 tissue mismatch")
    if clean(row.get("link_chromosome")) != clean(candidate.get("chrom")):
        raise ValueError("Step 6 link_chromosome with Step 3 chromosome mismatch")
    if clean(row.get("evidence_class")) != (
        "expression_gene_matches_allelome_endpoint"
    ):
        raise ValueError("Step 6 record is not expression gene endpoint evidence")
    if clean(row.get("evidence_level")) != "direct_feature_endpoint_evidence":
        raise ValueError("Step 6 record is not direct feature endpoint evidence")
    if clean(row.get("hit_endpoint_has_manual_new_gene")) != "yes":
        raise ValueError("Step 6 match endpoint no Manual new gene")
    if clean(row.get("link_contains_manual_new_gene")) != "yes":
        raise ValueError("Step 6 Allelome relationship no Manual new gene")
    if clean(row.get("new_gene_evidence_scope")) != (
        "direct_hit_endpoint_manual_new_gene"
    ):
        raise ValueError("Step 6 new gene evidence is not direct match endpoint")
    if (
        "expression_gene_id_equals_allelome_gene_endpoint"
        not in clean(row.get("endpoint_mapping_rules"))
    ):
        raise ValueError("Step 6 endpoint mapping rule is not Manual new gene direct mapping")


def endpoint_prefixes(hit: str) -> Tuple[List[str], str | None]:
    if hit == "base":
        return ["base"], "target"
    if hit == "target":
        return ["target"], "base"
    if hit == "both":
        return ["base", "target"], None
    raise ValueError(f"Step 6 endpoint_hit unsupported: {hit!r}")


def build_expanded_link(
    candidate_key_id: str,
    candidate: Mapping[str, str],
    row: Mapping[str, str],
) -> Dict[str, object]:
    hit = clean(row.get("endpoint_hit"))
    new_prefixes, partner_prefix = endpoint_prefixes(hit)
    new_gene_id = clean(candidate.get("eqtl_expression_gene_id"))

    for prefix in new_prefixes:
        if clean(row.get(f"{prefix}_record_type")) != "manual_new_gene":
            raise ValueError(
                f"{new_gene_id} of direct new gene endpoint is not manual_new_gene"
            )
        endpoint_ids = set(split_values(row.get(f"{prefix}_manual_new_gene_ids")))
        if (
            endpoint_ids != {new_gene_id}
            or clean(row.get(f"{prefix}_gene_id")) != new_gene_id
        ):
            raise ValueError(
                f"{new_gene_id} not in {prefix} endpoint Manual new gene ID in"
            )


    if partner_prefix is None:
        partner = {
            "partner_role": "none_both_endpoints_map_to_same_new_gene",
            "partner_record_type": "NA",
            "partner_gene_id": "NA",
            "partner_gene_symbol": "NA",
            "partner_gene_strand": "NA",
            "partner_chrom": "NA",
            "partner_start": "NA",
            "partner_end": "NA",
            "partner_name": "NA",
            "partner_biotype": "NA",
            "partner_endpoint_type": "same_new_gene_endpoint_pair",
            "candidate_orientation": "within_same_new_gene",
            "direction_basis": (
                "both_allelome_endpoints_map_to_same_manual_new_gene"
            ),
        }
    else:
        partner_family = clean(row.get(f"{partner_prefix}_endpoint_type"))
        if partner_family == "known_protein_coding":
            orientation = "new_gene_to_pcGene_display_orientation"
            direction_basis = (
                "posthoc_endpoint_class_orientation_not_causal_direction"
            )
        else:
            orientation = "undirected_new_gene_partner_pair"
            direction_basis = (
                "base_target_are_canonical_endpoints_not_causal_direction"
            )
        partner = {
            "partner_role": partner_prefix,
            "partner_record_type": clean(
                row.get(f"{partner_prefix}_record_type")
            )
            or "NA",
            "partner_gene_id": clean(row.get(f"{partner_prefix}_gene_id"))
            or "NA",
            "partner_gene_symbol": clean(
                row.get(f"{partner_prefix}_gene_symbol")
            )
            or "NA",
            "partner_gene_strand": clean(
                row.get(f"{partner_prefix}_gene_strand")
            )
            or "NA",
            "partner_chrom": clean(row.get(f"{partner_prefix}_chr")) or "NA",
            "partner_start": clean(row.get(f"{partner_prefix}_start")) or "NA",
            "partner_end": clean(row.get(f"{partner_prefix}_end")) or "NA",
            "partner_name": clean(row.get(f"{partner_prefix}_name")) or "NA",
            "partner_biotype": clean(
                row.get(f"{partner_prefix}_gene_biotype")
            )
            or "NA",
            "partner_endpoint_type": partner_family or "NA",
            "candidate_orientation": orientation,
            "direction_basis": direction_basis,
        }

    tier, tier_order, main_result, interpretation = TIER_INFO[
        (clean(candidate.get("result_layer")), "yes")
    ]
    lead_id = clean(candidate.get("gwas_lead_variant_id"))
    expanded: Dict[str, object] = {
        "candidate_id": candidate_key_id,
        "new_gene_id": new_gene_id,
        "evidence_tier": tier,
        "evidence_tier_order": tier_order,
        "main_result": main_result,
        "gwas_finemap_support_class": (
            "both_95pct_credible_sets"
            if clean(candidate.get("result_layer")) == "priority_both"
            else "gwas_noncredible_diagnostic"
        ),
        "gwas_locus_lead_variant_id": lead_id,
        "shared_snp_equals_gwas_locus_lead": (
            "yes"
            if clean(candidate.get("gwas_ld_variant_id")) == lead_id
            else "no"
        ),
        "multilayer_integration_scope": "same_newgene_eqtl_block_allelome_only",
        "evidence_interpretation": interpretation,
    }
    expanded.update(candidate)
    expanded.update(row)
    expanded.update(
        {
            "new_gene_endpoint_role": hit,
            "new_gene_endpoint_record_types": join_values(
                row.get(f"{prefix}_record_type") for prefix in new_prefixes
            ),
            "new_gene_endpoint_names": join_values(
                row.get(f"{prefix}_name") for prefix in new_prefixes
            ),
            "new_gene_endpoint_chromosomes": join_values(
                row.get(f"{prefix}_chr") for prefix in new_prefixes
            ),
            "new_gene_endpoint_starts0": join_values(
                row.get(f"{prefix}_start") for prefix in new_prefixes
            ),
            "new_gene_endpoint_ends0": join_values(
                row.get(f"{prefix}_end") for prefix in new_prefixes
            ),
            "endpoint_coordinate_system": "0-based_half-open",
            "new_gene_endpoint_manual_ids": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(f"{prefix}_manual_new_gene_ids")
                )
            ),
            "new_gene_endpoint_manual_strands": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(f"{prefix}_manual_new_gene_strands")
                )
            ),
            "new_gene_endpoint_manual_class_codes": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(f"{prefix}_manual_new_gene_class_codes")
                )
            ),
            "new_gene_endpoint_manual_novelty_categories": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(
                        f"{prefix}_manual_new_gene_novelty_categories"
                    )
                )
            ),
            "new_gene_endpoint_manual_candidate_categories": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(
                        f"{prefix}_manual_new_gene_candidate_categories"
                    )
                )
            ),
            "new_gene_endpoint_manual_dataset_support": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(f"{prefix}_manual_new_gene_dataset_support")
                )
            ),
            "new_gene_endpoint_manual_merge_statuses": join_values(
                item
                for prefix in new_prefixes
                for item in split_values(
                    row.get(f"{prefix}_manual_new_gene_merge_statuses")
                )
            ),
            "allelome_endpoint_variant_keys": "NA",
            "shared_finemap_snp_equals_allelome_endpoint": "not_evaluable_gene_endpoint",
            "shared_finemap_snp_min_distance_to_allelome_endpoint_bp": "NA",
        }
    )
    expanded.update(partner)
    return expanded


def format_number(value: float) -> str:
    return f"{value:.12g}"


def summarize_allelome_connection(
    rows: Sequence[Mapping[str, object]],
) -> Dict[str, object]:

    return {
        "multilayer_integration_scope": (
            "same_newgene_eqtl_block_allelome_only" if rows else "no_allelome_link"
        ),
        "exact_shared_snp_allelome_endpoint_rows": "NA",
        "shared_snp_equals_any_allelome_endpoint": (
            "not_evaluable_gene_endpoint" if rows else "no_allelome_link"
        ),
        "shared_snp_min_distance_to_allelome_endpoint_bp": "NA",
    }


def summarize_link_review(rows: Sequence[Mapping[str, object]]) -> Dict[str, object]:

    return {
        "allelome_evidence_strata": join_values(row.get("evidence_stratum") for row in rows),
        "allelome_multigene_link_rows": sum(
            clean(row.get("either_endpoint_multigene")) == "yes" for row in rows
        ),
        "allelome_overlap_link_rows": sum(
            as_int(row.get("overlap_bp"), "overlap_bp") > 0 for row in rows
        ),
    }


def build_pair_summary(
    rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    groups: MutableMapping[Tuple[str, ...], Dict[str, object]] = {}
    for row in rows:
        key = (
            clean(row.get("evidence_tier")),
            clean(row.get("candidate_id")),
            clean(row.get("partner_gene_id")),
            clean(row.get("partner_gene_symbol")),
            clean(row.get("mechanism")),
        )
        if any(not value or value == "NA" for value in key):
            raise ValueError(f"newgene-pcGene relationship contains empty aggregation key: {key}")
        group = groups.setdefault(
            key,
            {
                "first": row,
                "rows": 0,
                "dogs": set(),
                "endpoint_gene_ids": set(),
                "scores": [],
                "distances": [],
                "connection_rows": [],
            },
        )
        group["rows"] = int(group["rows"]) + 1
        group["dogs"].add(clean(row.get("animal_id")))

        group["endpoint_gene_ids"].update(
            split_values(row.get("new_gene_endpoint_manual_ids"))
        )
        group["scores"].append(
            as_float(row.get("linkage_score"), "linkage_score")
        )
        group["distances"].append(
            as_int(row.get("interval_gap_bp"), "interval_gap_bp")
        )
        group["connection_rows"].append(row)

    output: List[Dict[str, object]] = []
    for group in groups.values():
        first = group["first"]
        scores = group["scores"]
        distances = group["distances"]
        connection = summarize_allelome_connection(group["connection_rows"])
        output.append(
            {
                "evidence_tier": first["evidence_tier"],
                "evidence_tier_order": first["evidence_tier_order"],
                "main_result": first["main_result"],
                "candidate_id": first["candidate_id"],
                "new_gene_id": first["new_gene_id"],
                "eqtl_block_id": first["eqtl_block_id"],
                "eqtl_tissue": first["eqtl_tissue"],
                "chrom": first["chrom"],
                "normalized_variant_key": first["normalized_variant_key"],
                "gwas_trait": first["gwas_trait"],
                "gwas_locus_id": first["gwas_locus_id"],
                "gwas_ld_variant_id": first["gwas_ld_variant_id"],
                "gwas_p_lrt": first["gwas_p_lrt"],
                "gwas_pip": first["gwas_pip"],
                "gwas_pip_rank": first["gwas_pip_rank"],
                "gwas_locus_lead_variant_id": first[
                    "gwas_locus_lead_variant_id"
                ],
                "shared_snp_equals_gwas_locus_lead": first[
                    "shared_snp_equals_gwas_locus_lead"
                ],
                **{
                    field: first[field]
                    for field in GWAS_CONTEXT_COLUMNS
                },
                "multilayer_integration_scope": connection[
                    "multilayer_integration_scope"
                ],
                "shared_snp_equals_any_allelome_endpoint": connection[
                    "shared_snp_equals_any_allelome_endpoint"
                ],
                "shared_snp_min_distance_to_allelome_endpoint_bp": connection[
                    "shared_snp_min_distance_to_allelome_endpoint_bp"
                ],
                **summarize_link_review(group["connection_rows"]),
                "partner_gene_id": first["partner_gene_id"],
                "partner_gene_symbol": first["partner_gene_symbol"],
                "mechanism": first["mechanism"],
                "expanded_link_rows": group["rows"],
                "supporting_dog_count": len(group["dogs"]),
                "supporting_dogs": join_values(group["dogs"]),
                "allelome_endpoint_count": len(
                    group["endpoint_gene_ids"]
                ),
                "allelome_endpoints": join_values(
                    group["endpoint_gene_ids"]
                ),
                "linkage_score_min": format_number(min(scores)),
                "linkage_score_max": format_number(max(scores)),
                "interval_gap_bp_min": min(distances),
                "interval_gap_bp_max": max(distances),
                "candidate_orientation": first["candidate_orientation"],
                "direction_interpretation": (
                    "new_gene_to_pcGene specifies only the display orientation of endpoint categories; "
                    "Allelome enhancing or repressing relationships do not establish causal regulatory direction."
                ),
            }
        )
    return sorted(
        output,
        key=lambda row: (
            int(row["evidence_tier_order"]),
            clean(row["new_gene_id"]),
            clean(row["partner_gene_symbol"]),
            clean(row["mechanism"]),
        ),
    )


def chromosome_sort_key(value: object) -> Tuple[int, str]:
    chrom = clean(value)
    match = re.fullmatch(r"chr(\d+)", chrom)
    if match:
        return int(match.group(1)), chrom
    return 10_000, chrom


def candidate_sort_key(
    row: Mapping[str, object],
) -> Tuple[int, int, str, int, str, str]:
    chrom_order, chrom = chromosome_sort_key(row.get("chrom"))
    return (
        as_int(row.get("evidence_tier_order"), "evidence_tier_order"),
        chrom_order,
        chrom,
        as_int(row.get("pos"), "candidate pos"),
        clean(row.get("new_gene_id")),
        clean(row.get("eqtl_tissue")),
    )


def ordered_unique(fields: Iterable[str]) -> List[str]:
    result: List[str] = []
    seen: Set[str] = set()
    for field in fields:
        if field not in seen:
            seen.add(field)
            result.append(field)
    return result


def write_tsv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with gzip.open(path, "xt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(
    path: Path,
    inputs: Sequence[Tuple[str, Path, str]],
    outputs: Sequence[Tuple[str, Path, str]],
) -> None:
    rows: List[Dict[str, object]] = []
    for kind, items in (("input", inputs), ("output", outputs)):
        for label, item_path, purpose in items:
            rows.append(
                {
                    "kind": kind,
                    "label": label,
                    "path": str(item_path.resolve()),
                    "bytes": item_path.stat().st_size,
                    "purpose": purpose,
                }
            )
    write_tsv(
        path,
        ["kind", "label", "path", "bytes", "purpose"],
        rows,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Integrate five-dog new-gene, Allelome, eQTL, and GWAS fine-mapping evidence."
    )
    parser.add_argument(
        "--step6-links",
        type=Path,
        default=DEFAULT_STEP6_LINKS,
    )
    parser.add_argument(
        "--step6-summary",
        type=Path,
        default=DEFAULT_STEP6_SUMMARY,
    )
    parser.add_argument(
        "--step3-both-credible",
        type=Path,
        default=DEFAULT_STEP3_BOTH_CREDIBLE,
    )
    parser.add_argument(
        "--step3-noncredible",
        type=Path,
        default=DEFAULT_STEP3_NONCREDIBLE,
    )
    parser.add_argument(
        "--step3-summary",
        type=Path,
        default=DEFAULT_STEP3_SUMMARY,
    )
    parser.add_argument(
        "--gwas-run-dir",
        type=Path,
        default=DEFAULT_GWAS_RUN_DIR,
        help="Current GWAS fine-mapping run directory; only the main branch in locus_run_manifest is loaded.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
    )
    parser.add_argument(
        "--run-id",
        default=datetime.now().strftime("%Y%m%d_%H%M%S"),
    )
    return parser.parse_args()


def run(args: argparse.Namespace) -> int:
    global LOG_HANDLE

    gwas_manifest_path = args.gwas_run_dir / "locus_run_manifest.tsv"
    gwas_branch_summary_path = args.gwas_run_dir / "branch_summary.tsv"
    gwas_overall_summary_path = args.gwas_run_dir / "overall_summary.tsv"
    gwas_model_parameters_path = args.gwas_run_dir / "model_parameters.tsv"
    inputs = [
        (
            "step6_links",
            args.step6_links,
            "Per-dog direct relationships between expression genes and Allelome endpoints within the same eQTL block in step 6",
        ),
        (
            "step6_summary",
            args.step6_summary,
            "Step 6 PASS status, chromosome-wide phasing rules, and gene-link count requirements",
        ),
        (
            "step3_both_credible",
            args.step3_both_credible,
            "Step 3 exact overlaps with eQTL PPA>=0.1 and membership in both eQTL/GWAS 95% credible sets",
        ),
        (
            "step3_noncredible",
            args.step3_noncredible,
            "Step 3 diagnostic overlaps between high-PPA eQTL variants and non-credible-set variants within GWAS loci",
        ),
        (
            "step3_summary",
            args.step3_summary,
            "Step 3 fine-map status, thresholds, matching rules, and record-count requirements",
        ),
        ("gwas_manifest", gwas_manifest_path, "Locus identities, statuses, and per-locus file paths for each GWAS branch"),
        ("gwas_branch_summary", gwas_branch_summary_path, "GWAS main-branch task, variant, and credible-set count requirements"),
        ("gwas_overall_summary", gwas_overall_summary_path, "Completion status of the entire GWAS run"),
        ("gwas_model_parameters", gwas_model_parameters_path, "SuSiE-RSS version, coverage, and input-semantics requirements"),
    ]

    for label, path, _ in inputs:
        if not path.is_absolute():
            raise ValueError(f"{label} is not an absolute path: {path}")
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist or is not file: {path}")
    if not args.output_root.is_absolute():
        raise ValueError(f"output_root is not an absolute path: {args.output_root}")
    if not args.gwas_run_dir.is_absolute():
        raise ValueError(f"gwas_run_dir is not an absolute path: {args.gwas_run_dir}")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.run_id):
        raise ValueError(f"run_id contains not allowed character: {args.run_id!r}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    run_dir = args.output_root / args.run_id
    run_dir.mkdir(exist_ok=False)
    log_path = (
        run_dir
        / "1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap.log"
    )
    LOG_HANDLE = log_path.open("x", encoding="utf-8")

    all_path = run_dir / "dog5_newgene_eqtl_gwas_finemap_all_candidates.tsv"
    primary_path = (
        run_dir / "dog5_primary_both_credible_allelome_candidates.tsv"
    )
    secondary_path = (
        run_dir / "dog5_noncredible_diagnostic_allelome_candidates.tsv"
    )
    no_allelome_path = (
        run_dir / "dog5_genetic_candidates_without_allelome.tsv"
    )
    links_path = run_dir / "dog5_integrated_newgene_allelome_links.tsv.gz"
    pcgene_path = (
        run_dir / "dog5_integrated_newgene_pcGene_links.tsv.gz"
    )

    pair_path = run_dir / "dog5_candidate_newgene_pcGene_mechanism_pairs.tsv"
    summary_path = (
        run_dir / "dog5_newgene_allelome_eqtl_gwas_finemap_summary.tsv"
    )
    manifest_path = (
        run_dir / "dog5_newgene_allelome_eqtl_gwas_finemap_manifest.tsv"
    )

    try:
        log("[stage 1/6] check Step 6, Step 3 and GWAS main contract")
        step6_metrics = read_key_value_contract(
            args.step6_summary,
            "metric",
        )
        step3_metrics = read_key_value_contract(
            args.step3_summary,
            "metric",
        )
        gwas_model_metrics = read_key_value_contract(
            gwas_model_parameters_path,
            "parameter",
        )
        gwas_overall_metrics = read_key_value_contract(
            gwas_overall_summary_path,
            "metric",
        )
        branch_fields, branch_rows = read_tsv(gwas_branch_summary_path)
        require_columns(
            gwas_branch_summary_path,
            branch_fields,
            {"branch", "selected_task_count", "success_task_count", "failed_task_count", "modeled_variant_count", "credible_set_count"},
        )
        main_branch_rows = [
            row for row in branch_rows if clean(row.get("branch")) == "main"
        ]
        if len(main_branch_rows) != 1:
            raise ValueError("GWAS branch_summary must be and can only has one records main record")
        main_branch = main_branch_rows[0]
        if step6_metrics.get("overall_status") != "PASS":
            raise ValueError("Step 6 summary is not PASS")
        if step6_metrics.get("phase_scope_rule") != "phase_scope_must_equal_chrN":
            raise ValueError("Step 6 does not cover the whole chromosome phase_scope")
        if step6_metrics.get("local_ps_rule") != (
            "provenance_only_not_used_for_filtering"
        ):
            raise ValueError("Step 6 local PS rule is not only retain source")
        if step6_metrics.get("expression_matrix_scope") != (
            "reference_annotation_plus_dog5_manual_new_gene"
        ):
            raise ValueError("Step 6 is not reference annotation plus five-dog-specific new gene matrix")
        if step3_metrics.get("eqtl_final_qc_status") != "PASS":
            raise ValueError("Step 3 eQTL final QC is not PASS")
        if step3_metrics.get("eqtl_expression_matrix_scope") != (
            "reference_annotation_plus_dog5_manual_new_gene"
        ):
            raise ValueError("Step 3 is not reference annotation plus five-dog-specific new gene matrix")
        if step3_metrics.get("eqtl_feature_origin_rule") != (
            "gene-=>reference_annotation_gene;"
            "DOG5_NEWGENE_=>dog5_manual_new_gene"
        ):
            raise ValueError("Step 3 five-dog-specific new gene source rule mismatch")
        if step3_metrics.get("matching_rule") != (
            "same_CHR_POS_unordered_allele_set_no_complement_guess"
        ):
            raise ValueError("Step 3 variant matching does not use unordered allele sets")
        if step3_metrics.get("gwas_branch") != "main":
            raise ValueError("Step 3 summary is not main branch")
        if step3_metrics.get("allele_mismatch_qc_rows") != "0":
            raise ValueError("Step 3 exists allele does not match record")
        if step3_metrics.get("duplicate_output_qc_rows") != "0":
            raise ValueError("Step 3 contains duplicate output record")
        if gwas_model_metrics.get("coverage") != "0.95":
            raise ValueError("GWAS SuSiE-RSS coverage is not 0.95")
        if gwas_model_metrics.get("z_input") != "z_for_ld_a1":
            raise ValueError("GWAS SuSiE-RSS z input is not z_for_ld_a1")
        if gwas_model_metrics.get("n_input") != "n_locus_median":
            raise ValueError("GWAS SuSiE-RSS n input is not n_locus_median")
        if as_int(gwas_overall_metrics.get("failed_task_count"), "GWAS failed_task_count") != 0:
            raise ValueError("GWAS fine-map exists failed task")
        if as_int(main_branch.get("failed_task_count"), "GWAS main failed_task_count") != 0:
            raise ValueError("GWAS main branch exists failed task")
        if as_int(main_branch.get("selected_task_count"), "GWAS main selected_task_count") != as_int(
            main_branch.get("success_task_count"), "GWAS main success_task_count"
        ):
            raise ValueError("GWAS main branch select task count and success task count mismatch")

        log(
            "[stage 2/6] read Step 3 two level result and GWAS main per locus result, "
            "and filter Manual new gene candidate"
        )
        both_fields, both_rows = read_tsv(args.step3_both_credible)
        non_fields, non_rows = read_tsv(args.step3_noncredible)
        manifest_fields, manifest_rows = read_tsv(gwas_manifest_path)
        require_columns(
            args.step3_both_credible,
            both_fields,
            STEP3_REQUIRED_COLUMNS,
        )
        require_columns(
            args.step3_noncredible,
            non_fields,
            STEP3_REQUIRED_COLUMNS,
        )
        require_columns(
            gwas_manifest_path,
            manifest_fields,
            GWAS_MANIFEST_REQUIRED_COLUMNS,
        )
        manifest_index, gwas_variant_index = build_gwas_main_indexes(
            manifest_rows
        )
        if len(manifest_index) != as_int(
            main_branch.get("selected_task_count"),
            "GWAS main selected_task_count",
        ):
            raise ValueError("GWAS main manifest locus count and branch_summary mismatch")
        if len(gwas_variant_index) != as_int(
            main_branch.get("modeled_variant_count"),
            "GWAS main modeled_variant_count",
        ):
            raise ValueError("GWAS main per locus variant count with branch_summary mismatch")
        if sum(
            as_int(row.get("credible_set_count"), "credible_set_count")
            for row in manifest_index.values()
        ) != as_int(
            main_branch.get("credible_set_count"),
            "GWAS main credible_set_count",
        ):
            raise ValueError("GWAS main credible set count and branch_summary mismatch")
        if len(manifest_index) != as_int(
            step3_metrics.get("gwas_manifest_rows"),
            "Step 3 GWAS manifest row count",
        ):
            raise ValueError("GWAS main manifest with Step 3 summary row counts differ")
        if len(gwas_variant_index) != as_int(
            step3_metrics.get("gwas_variant_rows"),
            "Step 3 GWAS variant row count",
        ):
            raise ValueError("GWAS main variant count with Step 3 summary mismatch")
        if both_fields != non_fields:
            raise ValueError("Step 3 field order differs between the shared-credible-set and non-credible-set layers")
        if len(both_rows) != as_int(
            step3_metrics.get("priority_both_credible_overlap_rows"),
            "Step 3 both credible set row count",
        ):
            raise ValueError("Step 3 both credible set details and summary row counts differ")
        if len(non_rows) != as_int(
            step3_metrics.get("priority_noncredible_diagnostic_rows"),
            "Step 3 non-credible-set row count",
        ):
            raise ValueError("Step 3 row counts differ between non-credible-set details and summary")

        for row in both_rows:
            validate_step3_row(row, "priority_both")
        for row in non_rows:
            validate_step3_row(row, "noncredible")

        candidates = [
            row
            for row in both_rows + non_rows
            if is_new_gene_candidate(row)
        ]
        candidate_by_id: Dict[str, Dict[str, object]] = {}
        candidates_by_block: MutableMapping[
            Tuple[str, str, str, str],
            List[str],
        ] = defaultdict(list)
        block_signatures: Dict[Tuple[str, str, str, str], Tuple[str, ...]] = {}
        for row in candidates:
            candidate = dict(row)
            candidate.update(
                build_gwas_context(
                    candidate,
                    manifest_index,
                    gwas_variant_index,
                )
            )
            current_id = candidate_id(candidate)
            if not current_id or current_id in candidate_by_id:
                raise ValueError(f"Step 3 new gene candidate_id is empty or duplicate: {current_id}")
            current_block = block_key(candidate)
            if any(not value for value in current_block):
                raise ValueError(f"Step 3 new gene candidate contains empty block join key: {current_block}")
            signature = (
                clean(candidate.get("eqtl_feature_id")),
                clean(candidate.get("eqtl_feature_primary_id")),
                clean(candidate.get("eqtl_feature_label")),
                clean(candidate.get("eqtl_feature_class")),
                clean(candidate.get("eqtl_expression_source")),
                clean(candidate.get("eqtl_selected_sample_count")),
                clean(candidate.get("eqtl_common_marker_count")),
                clean(candidate.get("eqtl_total_reads")),
            )
            previous = block_signatures.setdefault(current_block, signature)
            if previous != signature:
                raise ValueError(
                    f"Step 3 same block expression gene metadata mismatch: {current_block}"
                )
            candidate_by_id[current_id] = candidate
            candidates_by_block[current_block].append(current_id)

        log(
            f"[stage 2/6] Step 3 input={len(both_rows) + len(non_rows)}; "
            f"GWAS main loci={len(manifest_index)}; "
            f"GWAS main variants={len(gwas_variant_index)}; "
            f"Manual new gene candidate={len(candidates)}; "
            "unique novel gene="
            f"{len({clean(row.get('eqtl_expression_gene_id')) for row in candidates})}; "
            f"unique block={len(candidates_by_block)}"
        )

        log("[stage 3/6] streaming scan Step 6 per-dog Allelome relationship and strict join")
        links_by_candidate: MutableMapping[
            str,
            List[Dict[str, object]],
        ] = defaultdict(list)
        expanded_rows: List[Dict[str, object]] = []
        seen_link_keys: Set[Tuple[str, str, str, str]] = set()
        observed_step6_rows = 0
        with open_text(args.step6_links, "r") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            step6_fields = list(reader.fieldnames or [])
            require_columns(
                args.step6_links,
                step6_fields,
                STEP6_REQUIRED_COLUMNS,
            )
            for line_number, row in enumerate(reader, start=2):
                observed_step6_rows += 1
                if observed_step6_rows % 100_000 == 0:
                    log(
                        "[stage 3/6] scanned Step 6"
                        f"{observed_step6_rows} row; current match expand"
                        f"{len(expanded_rows)} row"
                    )
                current_candidates = candidates_by_block.get(
                    step6_block_key(row)
                )
                if not current_candidates:
                    continue
                representative = candidate_by_id[current_candidates[0]]
                validate_step6_link(row, representative)
                for current_id in current_candidates:
                    candidate = candidate_by_id[current_id]
                    validate_step6_link(row, candidate)
                    unique_key = (
                        current_id,
                        clean(row.get("animal_id")),
                        clean(row.get("sample_pair_id")),
                        clean(row.get("endpoint_hit")),
                    )
                    if any(not value for value in unique_key):
                        raise ValueError(
                            f"{args.step6_links}:{line_number} contains empty expand key"
                        )
                    if unique_key in seen_link_keys:
                        raise ValueError(
                            f"{args.step6_links}:{line_number} expand relationship duplicate: "
                            f"{unique_key}"
                        )
                    seen_link_keys.add(unique_key)
                    expanded = build_expanded_link(
                        current_id,
                        candidate,
                        row,
                    )
                    expanded_rows.append(expanded)
                    links_by_candidate[current_id].append(expanded)

        expected_step6_rows = as_int(
            step6_metrics.get("expression_gene_link_overlap_rows"),
            "Step 6 expression_gene_link_overlap_rows",
        )
        if observed_step6_rows != expected_step6_rows:
            raise ValueError(
                "Step 6 gene-link details and summary row counts differ: "
                f"{observed_step6_rows} != {expected_step6_rows}"
            )
        log(
            f"[stage 3/6] Step 6 input={observed_step6_rows}; "
            f"new gene candidate expand relationship={len(expanded_rows)}"
        )

        log("[stage 4/6] summary candidate and generate A/B/C/D evidence layer")
        candidate_rows: List[Dict[str, object]] = []
        for current_id, candidate in candidate_by_id.items():
            rows = links_by_candidate.get(current_id, [])
            has_allelome = "yes" if rows else "no"
            connection = summarize_allelome_connection(rows)
            tier, tier_order, main_result, interpretation = TIER_INFO[
                (clean(candidate.get("result_layer")), has_allelome)
            ]
            lead_id = clean(candidate.get("gwas_lead_variant_id"))
            partner_gene_ids = {
                clean(row.get("partner_gene_id"))
                for row in rows
                if clean(row.get("partner_gene_id")) not in {"", "NA"}
                and clean(row.get("partner_record_type")) in {"known_gene", "manual_new_gene"}
            }
            partner_gene_symbols = {
                clean(row.get("partner_gene_symbol"))
                for row in rows
                if clean(row.get("partner_gene_symbol")) not in {"", "NA"}
                and clean(row.get("partner_record_type")) in {"known_gene", "manual_new_gene"}
            }
            result: Dict[str, object] = {
                "candidate_id": current_id,
                "new_gene_id": clean(
                    candidate.get("eqtl_expression_gene_id")
                ),
            }
            result.update(candidate)
            result.update(
                {
                    "gwas_finemap_support_class": (
                        "both_95pct_credible_sets"
                        if clean(candidate.get("result_layer"))
                        == "priority_both"
                        else "gwas_noncredible_diagnostic"
                    ),
                    "gwas_locus_lead_variant_id": lead_id,
                    "shared_snp_equals_gwas_locus_lead": (
                        "yes"
                        if clean(candidate.get("gwas_ld_variant_id")) == lead_id
                        else "no"
                    ),
                    "allelome_block_match": has_allelome,
                    "multilayer_integration_scope": connection[
                        "multilayer_integration_scope"
                    ],
                    "exact_shared_snp_allelome_endpoint_rows": connection[
                        "exact_shared_snp_allelome_endpoint_rows"
                    ],
                    "shared_snp_equals_any_allelome_endpoint": connection[
                        "shared_snp_equals_any_allelome_endpoint"
                    ],
                    "shared_snp_min_distance_to_allelome_endpoint_bp": (
                        connection[
                            "shared_snp_min_distance_to_allelome_endpoint_bp"
                        ]
                    ),
                    "expanded_allelome_link_rows": len(rows),
                    "allelome_link_count": len(
                        {
                            clean(row.get("sample_pair_id"))
                            for row in rows
                        }
                    ),
                    "supporting_dog_count": len(
                        {clean(row.get("animal_id")) for row in rows}
                    ),
                    "supporting_dogs": join_values(
                        row.get("animal_id") for row in rows
                    ),
                    "pcgene_link_count": sum(
                        clean(row.get("partner_endpoint_type")) == "known_protein_coding"
                        for row in rows
                    ),
                    "partner_gene_count": len(partner_gene_ids),
                    "partner_gene_ids": join_values(partner_gene_ids),
                    "partner_gene_symbols": join_values(
                        partner_gene_symbols
                    ),
                    "mechanisms": join_values(
                        row.get("mechanism") for row in rows
                    ),
                    "endpoint_manual_new_gene_ids": join_split_values(
                        rows,
                        "new_gene_endpoint_manual_ids",
                    ),
                    "endpoint_manual_new_gene_strands": join_split_values(
                        rows,
                        "new_gene_endpoint_manual_strands",
                    ),
                    "endpoint_manual_new_gene_class_codes": join_split_values(
                        rows,
                        "new_gene_endpoint_manual_class_codes",
                    ),
                    "endpoint_manual_new_gene_novelty_categories": (
                        join_split_values(
                            rows,
                            "new_gene_endpoint_manual_novelty_categories",
                        )
                    ),
                    "endpoint_manual_new_gene_candidate_categories": (
                        join_split_values(
                            rows,
                            "new_gene_endpoint_manual_candidate_categories",
                        )
                    ),
                    "endpoint_manual_new_gene_dataset_support": (
                        join_split_values(
                            rows,
                            "new_gene_endpoint_manual_dataset_support",
                        )
                    ),
                    "endpoint_manual_new_gene_merge_statuses": (
                        join_split_values(
                            rows,
                            "new_gene_endpoint_manual_merge_statuses",
                        )
                    ),
                    "evidence_tier": tier,
                    "evidence_tier_order": tier_order,
                    "main_result": main_result,
                    "evidence_interpretation": interpretation,
                    **summarize_link_review(rows),
                }
            )
            candidate_rows.append(result)

        candidate_rows.sort(key=candidate_sort_key)
        expanded_rows.sort(
            key=lambda row: (
                as_int(row.get("evidence_tier_order"), "evidence_tier_order"),
                chromosome_sort_key(row.get("chrom")),
                as_int(row.get("pos"), "expanded pos"),
                clean(row.get("new_gene_id")),
                clean(row.get("animal_id")),
                clean(row.get("canonical_pair_id")),
            )
        )
        pcgene_rows = [
            row
            for row in expanded_rows
            if clean(row.get("partner_endpoint_type")) == "known_protein_coding"
        ]
        pair_rows = build_pair_summary(pcgene_rows)
        tier_counts = Counter(
            clean(row.get("evidence_tier")) for row in candidate_rows
        )
        log(
            "[stage 4/6] "
            + "; ".join(
                f"{tier}={count}"
                for tier, count in sorted(tier_counts.items())
            )
        )

        log("[stage 5/6] stratify write candidate, per dog relationship and newgene-pcGene summary")
        candidate_fields = ordered_unique(
            ["candidate_id", "new_gene_id"]
            + both_fields
            + CANDIDATE_DERIVED_COLUMNS
        )
        link_fields = ordered_unique(
            LINK_PREFIX_COLUMNS
            + both_fields
            + step6_fields
            + LINK_DERIVED_COLUMNS
        )
        primary = [
            row
            for row in candidate_rows
            if row["evidence_tier"]
            == "A_BOTH_95PCT_CREDIBLE_PLUS_ALLELOME"
        ]
        secondary = [
            row
            for row in candidate_rows
            if row["evidence_tier"]
            == "B_GWAS_NONCREDIBLE_PLUS_ALLELOME"
        ]
        no_allelome = [
            row
            for row in candidate_rows
            if row["allelome_block_match"] == "no"
        ]
        write_tsv(all_path, candidate_fields, candidate_rows)
        write_tsv(primary_path, candidate_fields, primary)
        write_tsv(secondary_path, candidate_fields, secondary)
        write_tsv(no_allelome_path, candidate_fields, no_allelome)
        write_tsv_gz(links_path, link_fields, expanded_rows)
        write_tsv_gz(pcgene_path, link_fields, pcgene_rows)
        write_tsv(pair_path, PAIR_COLUMNS, pair_rows)

        primary_links = sum(
            row["evidence_tier"]
            == "A_BOTH_95PCT_CREDIBLE_PLUS_ALLELOME"
            for row in expanded_rows
        )
        secondary_links = sum(
            row["evidence_tier"]
            == "B_GWAS_NONCREDIBLE_PLUS_ALLELOME"
            for row in expanded_rows
        )
        exact_endpoint_rows = "NA"
        summary_rows = [
            {"metric": "run_id", "value": args.run_id},
            {"metric": "overall_status", "value": "PASS"},
            {"metric": "endpoint_coordinate_system", "value": "0-based_half-open"},
            {"metric": "allelome_link_count_unit", "value": "unique_sample_pair_id_within_candidate"},
            {"metric": "unique_newgene_pcGene_pairs", "value": len({
                (clean(row.get("new_gene_id")), clean(row.get("partner_gene_id")))
                for row in pcgene_rows
            })},
            *({"metric": key, "value": value} for key, value in summarize_link_review(expanded_rows).items()),
            {
                "metric": "step3_both_credible_input_rows",
                "value": len(both_rows),
            },
            {
                "metric": "step3_noncredible_input_rows",
                "value": len(non_rows),
            },
            {
                "metric": "gwas_main_locus_rows",
                "value": len(manifest_index),
            },
            {
                "metric": "gwas_main_variant_rows",
                "value": len(gwas_variant_index),
            },
            {
                "metric": "step3_new_gene_candidate_rows",
                "value": len(candidate_rows),
            },
            {
                "metric": "step3_new_gene_candidate_unique_genes",
                "value": len(
                    {clean(row.get("new_gene_id")) for row in candidate_rows}
                ),
            },
            {
                "metric": "step3_new_gene_candidate_unique_blocks",
                "value": len(candidates_by_block),
            },
            {
                "metric": "step6_expression_gene_link_input_rows",
                "value": observed_step6_rows,
            },
            {
                "metric": "tier_A_both_credible_plus_allelome_rows",
                "value": len(primary),
            },
            {
                "metric": "tier_B_noncredible_plus_allelome_rows",
                "value": len(secondary),
            },
            {
                "metric": "tier_C_both_credible_no_allelome_rows",
                "value": tier_counts.get(
                    "C_BOTH_95PCT_CREDIBLE_NO_ALLELOME",
                    0,
                ),
            },
            {
                "metric": "tier_D_noncredible_no_allelome_rows",
                "value": tier_counts.get(
                    "D_GWAS_NONCREDIBLE_NO_ALLELOME",
                    0,
                ),
            },
            {
                "metric": "integrated_allelome_link_rows",
                "value": len(expanded_rows),
            },
            {
                "metric": "primary_allelome_link_rows",
                "value": primary_links,
            },
            {
                "metric": "secondary_allelome_link_rows",
                "value": secondary_links,
            },
            {
                "metric": "integrated_newgene_pcGene_link_rows",
                "value": len(pcgene_rows),
            },
            {
                "metric": "candidate_pcGene_mechanism_pair_rows",
                "value": len(pair_rows),
            },
            {
                "metric": "shared_snp_equals_gwas_locus_lead_rows",
                "value": sum(
                    row["shared_snp_equals_gwas_locus_lead"] == "yes"
                    for row in candidate_rows
                ),
            },
            {
                "metric": "shared_snp_equals_gwas_top_pip_rows",
                "value": sum(
                    row["overlap_snp_equals_top_pip"] == "yes"
                    for row in candidate_rows
                ),
            },
            {
                "metric": "shared_snp_equals_allelome_endpoint_rows",
                "value": exact_endpoint_rows,
            },
            {
                "metric": "candidate_exact_three_layer_variant_rows",
                "value": "NA",
            },
            {
                "metric": "candidate_same_block_multilayer_rows",
                "value": sum(
                    row["multilayer_integration_scope"]
                    == "same_newgene_eqtl_block_allelome_only"
                    for row in candidate_rows
                ),
            },
            {
                "metric": "primary_join_rule",
                "value": "block_id+expression_gene_id+tissue+chrom",
            },
            {
                "metric": "candidate_identity_rule",
                "value": (
                    "result_layer+gwas_analysis_scope+gwas_branch+eqtl_block_id+"
                    "normalized_variant_key+gwas_locus_id+gwas_ld_variant_id"
                ),
            },
            {
                "metric": "gwas_lead_rule",
                "value": (
                    "read_explicit_gwas_lead_variant_id_from_main_manifest_then_"
                    "lookup_variant_pip_by_locus_id+ld_variant_id"
                ),
            },
            {
                "metric": "gwas_top_pip_rule",
                "value": (
                    "read_explicit_top_variant_id_and_max_pip_from_main_manifest_"
                    "then_lookup_variant_pip_by_locus_id+ld_variant_id"
                ),
            },
            {
                "metric": "eqtl_association_fdr_screening",
                "value": step6_metrics.get(
                    "association_fdr_screening",
                    "NA",
                ),
            },
            {
                "metric": "eqtl_null_causal_model",
                "value": step6_metrics.get("null_causal_model", "NA"),
            },
            {
                "metric": "phase_scope_rule",
                "value": step6_metrics.get("phase_scope_rule", "NA"),
            },
            {
                "metric": "local_ps_rule",
                "value": step6_metrics.get("local_ps_rule", "NA"),
            },
            {
                "metric": "direction_interpretation",
                "value": (
                    "new_gene_to_pcGene_is_display_orientation_only;"
                    "base_target_are_not_causal_direction"
                ),
            },
            {
                "metric": "integration_rule",
                "value": (
                    "Step3_exact_eQTL_GWAS_variant_overlap_joined_to_"
                    "Step6_Allelome_by_block_id+expression_gene_id+tissue+chrom"
                ),
            },
            {
                "metric": "variant_interpretation",
                "value": (
                    "Step6_newgene_endpoint_is_gene_interval;exact_variant_identity_"
                    "and_SNP_to_endpoint_distance_not_evaluable"
                ),
            },
            {
                "metric": "evidence_boundary",
                "value": (
                    "shared_finemap_candidate_plus_Allelome_convergence;"
                    "not_formal_colocalization_or_causal_proof"
                ),
            },
        ]
        write_tsv(summary_path, ["metric", "value"], summary_rows)

        log("[stage 6/6] write manifest and complete")
        outputs = [
            (
                "all_new_gene_candidates",
                all_path,
                "All Manual new-gene candidates from step 3 and their step 6 Allelome match status",
            ),
            (
                "primary_candidates",
                primary_path,
                "Primary evidence-tier candidates with exact sharing of both 95% credible sets and an Allelome relationship in the same new-gene/eQTL block",
            ),
            (
                "noncredible_diagnostic_candidates",
                secondary_path,
                "Candidates with diagnostic GWAS non-credible-set overlaps and Allelome relationships; a block may contain multiple candidate variants",
            ),
            (
                "without_allelome",
                no_allelome_path,
                "Candidates present in step 3 but without a direct same-block Allelome relationship in step 6",
            ),
            (
                "integrated_newgene_allelome_links",
                links_path,
                "Expanded Allelome relationships per candidate, dog, and canonical pair",
            ),
            (
                "integrated_newgene_pcGene_links",
                pcgene_path,
                "Candidate relationships with a new-gene endpoint and a pcGene partner",
            ),
            (
                "candidate_newgene_pcGene_mechanism_pairs",
                pair_path,
                "Within each fine-map candidate, aggregate across dogs by pcGene and mechanism, retaining the same-block integration scope and review flags.",
            ),
            (
                "integration_summary",
                summary_path,
                "Input reconciliation, A/B/C/D counts, joining rules, and interpretation limits",
            ),
            (
                "run_log",
                log_path,
                "Stage-specific run logs",
            ),
        ]
        write_manifest(manifest_path, inputs, outputs)
        log(f"result directory: {run_dir.resolve()}")
        log(
            f"both credible set plus Allelome candidate={len(primary)}; non-credible set plus Allelome diagnostic candidate={len(secondary)}; "
            f"candidate-newgene-pcGene-mechanism summary row={len(pair_rows)}"
        )
        log("task completed / done")
        log("execution_status=success")
        log(
            "completion only indicates that the main program exited normally; output completeness, format, "
            "counts and biological quality require separate review"
        )
        return 0
    except Exception as exc:
        log(f"task failed: {exc}")
        log("execution_status=failed")
        traceback.print_exc(file=LOG_HANDLE)
        traceback.print_exc(file=sys.stderr)
        return 1
    finally:
        if LOG_HANDLE is not None:
            LOG_HANDLE.close()
            LOG_HANDLE = None


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
