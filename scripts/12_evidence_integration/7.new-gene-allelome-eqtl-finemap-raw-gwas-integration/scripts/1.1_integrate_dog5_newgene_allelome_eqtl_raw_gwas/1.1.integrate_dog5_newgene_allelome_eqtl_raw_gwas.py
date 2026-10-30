#!/usr/bin/python3.10
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


DEFAULT_EQTL_GWAS_DETAIL = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/2.eqtl-finemap-compare-raw-gwas/1.1_export_dog5_eqtl_finemap_gwas_exact_overlap/20260917_195542/five_dogs_eqtl_finemap_gwas_exact_overlaps.csv'
)
DEFAULT_EQTL_GWAS_SUMMARY = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/2.eqtl-finemap-compare-raw-gwas/1.1_export_dog5_eqtl_finemap_gwas_exact_overlap/20260917_195542/five_dogs_eqtl_finemap_gwas_exact_overlap_overall_summary.csv'
)
DEFAULT_ALLELOME_EQTL_LINKS = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740/dog5_allelome_eqtl_finemap_expression_gene_link_overlap.tsv.gz'
)
DEFAULT_ALLELOME_EQTL_SUMMARY = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740/dog5_allelome_eqtl_finemap_run_summary.tsv'
)
DEFAULT_NEW_GENE_TRANSCRIPT_MAP = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv'
)
DEFAULT_NEW_GENE_BUILD_SUMMARY = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.summary.tsv'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/7.new-gene-allelome-eqtl-finemap-raw-gwas-integration/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas'
)

EXPECTED_EXPRESSION_MATRIX_SCOPE = (
    "reference_annotation_plus_dog5_manual_new_gene"
)
EXPECTED_FEATURE_ORIGIN_RULE = (
    "gene-=>reference_annotation_gene;DOG5_NEWGENE_=>dog5_manual_new_gene"
)
EXPECTED_NEW_GENE_PREFIX = "DOG5_NEWGENE_"
EXPECTED_NEW_GENE_SOURCE = "five-dogs"
EXPECTED_NEW_GENE_SELECTION_STAGE = (
    "2.1_transcript_first_preASE_candidate_assessment"
)

EQTL_GWAS_REQUIRED_COLUMNS = {
    "eqtl_finemap_snp",
    "normalized_variant_key",
    "chrom",
    "pos",
    "ref",
    "alt",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "feature_origin",
    "expression_matrix_scope",
    "expression_gene_id",
    "tissue",
    "block_id",
    "finemap_role",
    "is_top_variant",
    "is_credible_set_member",
    "finemap_ppa",
    "ppa_tier",
    "marker_rank_by_ppa",
    "block_max_ppa",
    "block_top_ppa_tie_count",
    "block_top_marker_is_unique",
    "block_credible_set_size",
    "block_credible_set_ppa_sum",
    "selected_sample_count",
    "common_marker_count",
    "total_reads",
    "phase_scope",
    "local_ps_by_sample",
    "gwas_trait",
    "gwas_marker_id",
    "gwas_p_lrt",
    "gwas_variant_p_tier",
    "interpretation_note",
    "gwas_p_value_status",
    "gwas_modelled_locus_count",
    "gwas_region_match_status",
    "is_p_lt_1e-5",
    "gwas_af",
    "gwas_beta",
    "gwas_se",
    "gwas_allele1",
    "gwas_allele0",
    "allele_order_status",
    "allele_match_status",
    "beta_alignment_status",
    "overlap_type",
}

ALLELOME_REQUIRED_COLUMNS = {
    "block_id",
    "eqtl_tissue",
    "feature_id",
    "feature_primary_id",
    "expression_gene_id",
    "eqtl_chrom",
    "selected_sample_count",
    "selected_samples",
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

NEW_GENE_TRANSCRIPT_MAP_REQUIRED_COLUMNS = {
    "final_gene_id",
    "final_transcript_id",
    "source_dataset",
    "selection_stage",
    "original_gene_id",
    "chr",
    "strand",
    "gene_start_0based",
    "gene_end_0based",
    "class_code",
    "novelty_category",
    "FEELnc_biotype",
}

CANDIDATE_COLUMNS = [
    "candidate_id",
    "block_id",
    "new_gene_id",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "tissue",
    "chrom",
    "eqtl_gwas_pos",
    "eqtl_gwas_ref",
    "eqtl_gwas_alt",
    "normalized_variant_key",
    "eqtl_finemap_snp",
    "finemap_role",
    "is_top_variant",
    "is_credible_set_member",
    "finemap_ppa",
    "ppa_tier",
    "marker_rank_by_ppa",
    "block_max_ppa",
    "block_top_ppa_tie_count",
    "block_top_marker_is_unique",
    "block_credible_set_size",
    "block_credible_set_ppa_sum",
    "selected_sample_count",
    "common_marker_count",
    "total_reads",
    "phase_scope",
    "local_ps_by_sample",
    "gwas_trait",
    "gwas_marker_id",
    "gwas_p_lrt",
    "gwas_variant_p_tier",
    "interpretation_note",
    "gwas_p_value_status",
    "gwas_modelled_locus_count",
    "gwas_region_match_status",
    "is_p_lt_1e-5",
    "gwas_af",
    "gwas_beta",
    "gwas_se",
    "gwas_allele1",
    "gwas_allele0",
    "allele_order_status",
    "allele_match_status",
    "beta_alignment_status",
    "source_overlap_type",
    "new_gene_chr",
    "new_gene_start0",
    "new_gene_end0",
    "new_gene_span_bp",
    "new_gene_strand",
    "new_gene_merge_status",
    "new_gene_dataset_support",
    "new_gene_class_codes",
    "new_gene_novelty_categories",
    "new_gene_feelnc_biotypes",
    "eqtl_gwas_snp_inside_new_gene_span",
    "eqtl_gwas_snp_distance_to_new_gene_span_bp",
    "allelome_block_match",
    "allelome_link_count",
    "supporting_dog_count",
    "supporting_dogs",
    "pcgene_link_count",
    "partner_gene_count",
    "partner_gene_ids",
    "partner_gene_symbols",
    "mechanisms",
    "pcgene_expanded_link_count",
    "unique_pcgene_count",
    "exact_endpoint_link_count",
]

CANDIDATE_SUMMARY_COLUMN_COUNT = 12

EXPANDED_LINK_COLUMNS = CANDIDATE_COLUMNS[
    :-CANDIDATE_SUMMARY_COLUMN_COUNT
] + [
    "animal_id",
    "sample_tissue",
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
    "new_gene_endpoint_role",
    "allelome_endpoint_gene_id",
    "allelome_endpoint_start0",
    "allelome_endpoint_end0",
    "allelome_endpoint_strand",
    "allelome_endpoint_chrom",
    "allelome_endpoint_pos",
    "allelome_endpoint_ref",
    "allelome_endpoint_alt",
    "allelome_endpoint_name",
    "eqtl_gwas_snp_equals_allelome_endpoint",
    "eqtl_gwas_snp_distance_to_allelome_endpoint_bp",
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
    "candidate_direction",
    "pair_orientation",
    "causal_direction_supported",
    "evidence_tier",
    "direction_basis",
    "integration_level",
]


class RunLogger:
    def __init__(self, path: Path) -> None:
        self.handle = path.open("w", encoding="utf-8")

    def log(self, message: str) -> None:
        line = f"[{datetime.now():%Y-%m-%d %H:%M:%S} ] {message}"
        print(line, flush=True)
        print(line, file=self.handle, flush=True)

    def close(self) -> None:
        self.handle.close()


LOGGER: RunLogger | None = None


def log(message: str) -> None:
    if LOGGER is None:
        print(f"[{datetime.now():%Y-%m-%d %H:%M:%S} ] {message}", flush=True)
    else:
        LOGGER.log(message)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def output_value(value: object) -> str:
    text = clean(value)
    return text if text and text not in {".", "NA"} else "NA"


def parse_int(value: object, label: str) -> int:
    try:
        return int(float(clean(value)))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not an integer: {value!r}") from exc


def parse_float(value: object, label: str) -> float:
    try:
        result = float(clean(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite numeric value: {value!r}")
    return result


def normalize_chrom(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    if text.startswith("chr") or text.startswith("NW_"):
        return text
    if text.isdigit() or text in {"X", "Y", "M", "MT"}:
        return f"chr{text}"
    return text


def split_values(*values: object) -> List[str]:
    result: List[str] = []
    seen: Set[str] = set()
    for value in values:
        for item in re.split(r"[|;,]", clean(value)):
            item = item.strip()
            if item and item not in {"NA", ".", "-"} and item not in seen:
                seen.add(item)
                result.append(item)
    return result


def join_values(values: Iterable[object]) -> str:
    result: List[str] = []
    seen: Set[str] = set()
    for value in values:
        text = clean(value)
        if text and text not in {"NA", ".", "-"} and text not in seen:
            seen.add(text)
            result.append(text)
    return "|".join(result) if result else "NA"


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} does not exist or is empty: {path}")


def require_columns(
    path: Path, observed: Sequence[str], required: Set[str]
) -> None:
    missing = sorted(required - set(observed))
    if missing:
        raise ValueError(f"{path} missing required fields: {','.join(missing)}")


def read_metric_table(
    path: Path, delimiter: str, key_column: str
) -> Dict[str, str]:
    encoding = "utf-8-sig" if delimiter == "," else "utf-8"
    with path.open("r", encoding=encoding, newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        require_columns(path, reader.fieldnames or [], {key_column, "value"})
        result: Dict[str, str] = {}
        for row in reader:
            key = clean(row.get(key_column))
            if not key or key in result:
                raise ValueError(f"{path} contains empty or duplicate metric: {key!r}")
            result[key] = clean(row.get("value"))
    return result


def allele_set(ref: object, alt: object) -> Tuple[str, str]:
    alleles = sorted([clean(ref).upper(), clean(alt).upper()])
    if not alleles[0] or not alleles[1] or alleles[0] == alleles[1]:
        raise ValueError(f"REF/ALT invalid: {ref!r}, {alt!r}")
    return alleles[0], alleles[1]


def point_interval_distance(pos1: int, start0: int, end0: int) -> int:
    pos0 = pos1 - 1
    if start0 <= pos0 < end0:
        return 0
    if pos0 < start0:
        return start0 - pos0
    return pos0 - end0 + 1


def load_new_gene_metadata(
    transcript_map_path: Path,
    summary_path: Path,
) -> Dict[str, Dict[str, object]]:
    summary = read_metric_table(summary_path, "\t", "metric")
    required_metrics = {
        "source_dataset",
        "selection_stage",
        "candidate_gene_count",
        "candidate_transcript_count",
        "final_gene_id_prefix",
        "transcript_map",
    }
    missing_metrics = sorted(required_metrics - set(summary))
    if missing_metrics:
        raise ValueError(
            f"{summary_path} missing metric: {','.join(missing_metrics)}"
        )
    expected_summary_values = {
        "source_dataset": EXPECTED_NEW_GENE_SOURCE,
        "selection_stage": EXPECTED_NEW_GENE_SELECTION_STAGE,
        "final_gene_id_prefix": EXPECTED_NEW_GENE_PREFIX,
    }
    for key, expected in expected_summary_values.items():
        if summary[key] != expected:
            raise ValueError(
                f"five-only new-gene summary {key} unexpected: "
                f"observed={summary[key]!r}, expected={expected!r}"
            )
    recorded_map = Path(summary["transcript_map"]).resolve()
    if recorded_map != transcript_map_path.resolve():
        raise ValueError(
            "five-only new-gene summary recorded transcript map does not match the input: "
            f"summary={recorded_map}, input={transcript_map_path.resolve()}"
        )
    expected_gene_count = parse_int(
        summary["candidate_gene_count"],
        "new-gene summary candidate_gene_count",
    )
    expected_transcript_count = parse_int(
        summary["candidate_transcript_count"],
        "new-gene summary candidate_transcript_count",
    )


    result: Dict[str, Dict[str, object]] = {}
    identities: Dict[str, Tuple[str, str, int, int, str]] = {}
    transcript_ids: Set[str] = set()
    class_codes: MutableMapping[str, List[str]] = defaultdict(list)
    novelty_categories: MutableMapping[str, List[str]] = defaultdict(list)
    feelnc_biotypes: MutableMapping[str, List[str]] = defaultdict(list)
    observed_rows = 0
    with transcript_map_path.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            transcript_map_path,
            reader.fieldnames or [],
            NEW_GENE_TRANSCRIPT_MAP_REQUIRED_COLUMNS,
        )
        for line_number, row in enumerate(reader, start=2):
            observed_rows += 1
            gene_id = clean(row.get("final_gene_id"))
            transcript_id = clean(row.get("final_transcript_id"))
            if not gene_id.startswith(EXPECTED_NEW_GENE_PREFIX):
                raise ValueError(
                    f"{transcript_map_path}:{line_number} final_gene_id prefix unexpected"
                )
            if not transcript_id or transcript_id in transcript_ids:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} "
                    "final_transcript_id missing or duplicate"
                )
            transcript_ids.add(transcript_id)
            if (
                clean(row.get("source_dataset")) != EXPECTED_NEW_GENE_SOURCE
                or clean(row.get("selection_stage"))
                != EXPECTED_NEW_GENE_SELECTION_STAGE
            ):
                raise ValueError(
                    f"{transcript_map_path}:{line_number} "
                    "data source or filter stage unexpected"
                )
            chrom = normalize_chrom(row.get("chr"))
            strand = clean(row.get("strand"))
            start0 = parse_int(
                row.get("gene_start_0based"),
                f"{transcript_map_path}:{line_number} gene_start_0based",
            )
            end0 = parse_int(
                row.get("gene_end_0based"),
                f"{transcript_map_path}:{line_number} gene_end_0based",
            )
            original_gene_id = clean(row.get("original_gene_id"))
            if (
                not chrom
                or strand not in {"+", "-"}
                or start0 < 0
                or end0 <= start0
                or not original_gene_id
            ):
                raise ValueError(
                    f"{transcript_map_path}:{line_number} new gene metadata invalid"
                )
            identity = (chrom, strand, start0, end0, original_gene_id)
            if gene_id in identities and identities[gene_id] != identity:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} "
                    "the same final_gene_id metadata mismatch"
                )
            identities[gene_id] = identity
            if gene_id not in result:
                result[gene_id] = {
                    "new_gene_id": gene_id,
                    "new_gene_chr": chrom,
                    "new_gene_start0": start0,
                    "new_gene_end0": end0,
                    "new_gene_span_bp": end0 - start0,
                    "new_gene_strand": strand,
                    "new_gene_merge_status": (
                        "five_only_not_cross_resource_merged"
                    ),
                    "new_gene_dataset_support": EXPECTED_NEW_GENE_SOURCE,
                }
            for values, value, label in (
                (class_codes[gene_id], clean(row.get("class_code")), "class_code"),
                (
                    novelty_categories[gene_id],
                    clean(row.get("novelty_category")),
                    "novelty_category",
                ),
                (
                    feelnc_biotypes[gene_id],
                    clean(row.get("FEELnc_biotype")),
                    "FEELnc_biotype",
                ),
            ):
                if not value:
                    raise ValueError(
                        f"{transcript_map_path}:{line_number} {label} missing"
                    )
                if value not in values:
                    values.append(value)

    if (
        observed_rows != expected_transcript_count
        or len(result) != expected_gene_count
    ):
        raise ValueError(
            "five-only transcript map count with build summary mismatch: "
            f"rows={observed_rows}/{expected_transcript_count}, "
            f"genes={len(result)}/{expected_gene_count}"
        )
    for gene_id, meta in result.items():
        meta["new_gene_class_codes"] = join_values(class_codes[gene_id])
        meta["new_gene_novelty_categories"] = join_values(
            novelty_categories[gene_id]
        )
        meta["new_gene_feelnc_biotypes"] = join_values(
            feelnc_biotypes[gene_id]
        )
    return result


def load_new_gene_candidates(
    detail_path: Path,
    summary_path: Path,
    metadata: Mapping[str, Mapping[str, object]],
) -> Tuple[List[Dict[str, object]], int]:
    summary = read_metric_table(summary_path, ",", "metric")
    required_metrics = {
        "exact_overlap_rows",
        "gwas_candidate_p_threshold",
        "allele_matching_rule",
        "proximity_analysis_applied",
        "expression_matrix_scope",
        "feature_origin_rule",
        "suggestive_p_threshold",
    }
    missing_metrics = sorted(required_metrics - set(summary))
    if missing_metrics:
        raise ValueError(
            f"{summary_path} missing metric: {','.join(missing_metrics)}"
        )
    if summary["allele_matching_rule"] != "CHROM_POS_unordered_allele_set":
        raise ValueError("2 results use an unexpected allele-matching rule")
    if summary["proximity_analysis_applied"] != "0":
        raise ValueError("2 result unexpected perform position nearby analysis")
    if summary["expression_matrix_scope"] != EXPECTED_EXPRESSION_MATRIX_SCOPE:
        raise ValueError("2 results do not use the current five-dog expression-matrix scope")
    if summary["feature_origin_rule"] != EXPECTED_FEATURE_ORIGIN_RULE:
        raise ValueError("2 result feature source rule is not current five-dog-specific contract")
    expected_rows = parse_int(
        summary["exact_overlap_rows"], "2 result exact_overlap_rows"
    )
    gwas_threshold = parse_float(
        summary["gwas_candidate_p_threshold"],
        "2 result gwas_candidate_p_threshold",
    )
    cumulative_p_threshold = parse_float(
        summary["suggestive_p_threshold"],
        "2 result suggestive_p_threshold",
    )
    if not math.isclose(
        cumulative_p_threshold,
        1e-5,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise ValueError("2 result cumulative P<1e-5 threshold is not 1e-5")


    if not math.isclose(gwas_threshold, 1e-4, rel_tol=0, abs_tol=1e-15):
        raise ValueError("2 result input pool is not P<=1e-4")
    for key, expected in (("genomewide_p_threshold", 2.97e-6), ("exploratory_point_p_threshold", 5e-5)):
        if not math.isclose(parse_float(summary.get(key), key), expected, rel_tol=0, abs_tol=1e-15):
            raise ValueError("2 result single site threshold mismatch: " + key)
    candidates: List[Dict[str, object]] = []
    seen_candidate_ids: Set[str] = set()
    observed_rows = 0
    with detail_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(
            detail_path, reader.fieldnames or [], EQTL_GWAS_REQUIRED_COLUMNS
        )
        for line_number, row in enumerate(reader, start=2):
            observed_rows += 1
            if clean(row.get("allele_match_status")) != "allele_set_match":
                raise ValueError(
                    f"{detail_path}:{line_number} is not allele set exact matching"
                )
            if (
                clean(row.get("overlap_type"))
                != "exact_variant_p_le_1e-4_allele_set_match"
            ):
                raise ValueError(
                    f"{detail_path}:{line_number} overlap_type is not expected value"
                )
            if clean(row.get("is_credible_set_member")) != "1":
                raise ValueError(
                    f"{detail_path}:{line_number} eQTL variant not in credible set"
                )
            ppa = parse_float(
                row.get("finemap_ppa"),
                f"{detail_path}:{line_number} finemap_ppa",
            )
            gwas_p = parse_float(
                row.get("gwas_p_lrt"),
                f"{detail_path}:{line_number} gwas_p_lrt",
            )
            if not 0.1 <= ppa <= 1.0 + 1e-8:
                raise ValueError(
                    f"{detail_path}:{line_number} PPA not in[0.1,1]span"
                )
            if not 0 <= gwas_p <= gwas_threshold:
                raise ValueError(
                    f"{detail_path}:{line_number} GWAS p value exceeds input threshold"
                )
            expected_tier = ("genomewide" if gwas_p < 2.97e-6 else
                             "suggestive" if gwas_p < 1e-5 else
                             "exploratory_point" if gwas_p < 5e-5 else
                             "extended_exploratory_point")
            if row.get("gwas_p_value_status") != "valid" or row.get("gwas_variant_p_tier") != expected_tier:
                raise ValueError(f"{detail_path}:{line_number} single site category with P value mismatch")
            observed_strong_flag = clean(
                row.get("is_p_lt_1e-5")
            )
            expected_strong_flag = (
                "1" if gwas_p < cumulative_p_threshold else "0"
            )
            if observed_strong_flag != expected_strong_flag:
                raise ValueError(
                    f"{detail_path}:{line_number} cumulative P<1e-5p value flag with p_lrt mismatch"
                )
            if (
                clean(row.get("expression_matrix_scope"))
                != EXPECTED_EXPRESSION_MATRIX_SCOPE
            ):
                raise ValueError(
                    f"{detail_path}:{line_number} expression_matrix_scope unexpected"
                )
            chrom = normalize_chrom(row.get("chrom"))
            if normalize_chrom(row.get("phase_scope")) != chrom:
                raise ValueError(
                    f"{detail_path}:{line_number} phase_scope is not all chromosome"
                )
            gene_id = clean(row.get("expression_gene_id"))
            if gene_id not in metadata:
                if gene_id.startswith(EXPECTED_NEW_GENE_PREFIX):
                    raise ValueError(f"novel candidate gene missing five dogs annotation mapping: {gene_id}")
                continue
            if (
                clean(row.get("feature_origin"))
                != "dog5_manual_new_gene"
                or not gene_id.startswith(EXPECTED_NEW_GENE_PREFIX)
            ):
                raise ValueError(
                    f"{detail_path}:{line_number} new gene source field unexpected"
                )
            meta = metadata[gene_id]
            if normalize_chrom(meta["new_gene_chr"]) != chrom:
                raise ValueError(
                    f"{detail_path}:{line_number} eQTL and novel gene chromosome mismatch"
                )
            pos = parse_int(row.get("pos"), f"{detail_path}:{line_number} pos")
            ref = clean(row.get("ref")).upper()
            alt = clean(row.get("alt")).upper()
            allele_set(ref, alt)
            distance = point_interval_distance(
                pos,
                int(meta["new_gene_start0"]),
                int(meta["new_gene_end0"]),
            )
            candidate_id = "::".join(
                [
                    clean(row.get("block_id")),
                    clean(row.get("normalized_variant_key")),
                    clean(row.get("gwas_trait")),
                    clean(row.get("gwas_marker_id")),
                ]
            )
            if not all(candidate_id.split("::")) or candidate_id in seen_candidate_ids:
                raise ValueError(
                    f"{detail_path}:{line_number} candidate_id missing or duplicate"
                )
            seen_candidate_ids.add(candidate_id)
            candidate: Dict[str, object] = {
                "candidate_id": candidate_id,
                "block_id": clean(row.get("block_id")),
                "new_gene_id": gene_id,
                "feature_id": clean(row.get("feature_id")),
                "feature_primary_id": clean(
                    row.get("feature_primary_id")
                ),
                "feature_label": clean(row.get("feature_label")),
                "feature_class": clean(row.get("feature_class")),
                "tissue": clean(row.get("tissue")),
                "chrom": chrom,
                "eqtl_gwas_pos": pos,
                "eqtl_gwas_ref": ref,
                "eqtl_gwas_alt": alt,
                "normalized_variant_key": clean(
                    row.get("normalized_variant_key")
                ),
                "eqtl_finemap_snp": clean(row.get("eqtl_finemap_snp")),
                "finemap_role": clean(row.get("finemap_role")),
                "is_top_variant": clean(row.get("is_top_variant")),
                "is_credible_set_member": clean(
                    row.get("is_credible_set_member")
                ),
                "finemap_ppa": f"{ppa:.15g}",
                "ppa_tier": clean(row.get("ppa_tier")),
                "marker_rank_by_ppa": clean(row.get("marker_rank_by_ppa")),
                "block_max_ppa": clean(row.get("block_max_ppa")),
                "block_top_ppa_tie_count": clean(row.get("block_top_ppa_tie_count")),
                "block_top_marker_is_unique": clean(row.get("block_top_marker_is_unique")),
                "block_credible_set_size": clean(row.get("block_credible_set_size")),
                "block_credible_set_ppa_sum": clean(row.get("block_credible_set_ppa_sum")),
                "selected_sample_count": clean(
                    row.get("selected_sample_count")
                ),
                "common_marker_count": clean(row.get("common_marker_count")),
                "total_reads": clean(row.get("total_reads")),
                "phase_scope": clean(row.get("phase_scope")),
                "local_ps_by_sample": output_value(
                    row.get("local_ps_by_sample")
                ),
                "gwas_trait": clean(row.get("gwas_trait")),
                "gwas_marker_id": clean(row.get("gwas_marker_id")),
                "gwas_p_lrt": f"{gwas_p:.15g}",
                "gwas_variant_p_tier": clean(row.get("gwas_variant_p_tier")),
                "interpretation_note": clean(row.get("interpretation_note")),
                "gwas_p_value_status": clean(row.get("gwas_p_value_status")),
                "gwas_modelled_locus_count": clean(row.get("gwas_modelled_locus_count")),
                "gwas_region_match_status": clean(row.get("gwas_region_match_status")),
                "is_p_lt_1e-5": observed_strong_flag,
                "gwas_af": clean(row.get("gwas_af")),
                "gwas_beta": clean(row.get("gwas_beta")),
                "gwas_se": clean(row.get("gwas_se")),
                "gwas_allele1": clean(row.get("gwas_allele1")),
                "gwas_allele0": clean(row.get("gwas_allele0")),
                "allele_order_status": clean(
                    row.get("allele_order_status")
                ),
                "allele_match_status": clean(
                    row.get("allele_match_status")
                ),
                "beta_alignment_status": clean(
                    row.get("beta_alignment_status")
                ),
                "source_overlap_type": clean(row.get("overlap_type")),
                "eqtl_gwas_snp_inside_new_gene_span": (
                    "yes" if distance == 0 else "no"
                ),
                "eqtl_gwas_snp_distance_to_new_gene_span_bp": distance,
            }
            candidate.update(meta)
            candidates.append(candidate)

    if observed_rows != expected_rows:
        raise ValueError(
            "2 number details row count and summary mismatch: "
            f"observed={observed_rows}, expected={expected_rows}"
        )
    return candidates, observed_rows


def endpoint_contains_gene(
    row: Mapping[str, str], role: str, gene_id: str
) -> bool:
    return (clean(row.get(f"{role}_record_type")) == "manual_new_gene"
            and clean(row.get(f"{role}_gene_id")) == gene_id)


def endpoint_fields(
    row: Mapping[str, str], role: str
) -> Dict[str, str]:
    return {
        "record_type": clean(row.get(f"{role}_record_type")),
        "gene_id": clean(row.get(f"{role}_gene_id")),
        "gene_symbol": clean(row.get(f"{role}_gene_symbol")),
        "gene_strand": clean(row.get(f"{role}_gene_strand")),
        "chrom": normalize_chrom(row.get(f"{role}_chr")),
        "start": clean(row.get(f"{role}_start")),
        "end": clean(row.get(f"{role}_end")),
        "name": clean(row.get(f"{role}_name")),
        "biotype": clean(row.get(f"{role}_gene_biotype")),
        "endpoint_type": clean(row.get(f"{role}_endpoint_type")),
    }


def build_expanded_link(
    candidate: Mapping[str, object],
    row: Mapping[str, str],
) -> Dict[str, object]:
    hit = clean(row.get("endpoint_hit"))
    if hit not in {"base", "target", "both"}:
        raise ValueError(f"unsupported of endpoint_hit: {hit!r}")


    if hit not in {"base", "target"}:
        raise ValueError("current direct gene endpoint relationship must be unique match base or target")
    if not endpoint_contains_gene(row, hit, str(candidate["new_gene_id"])):
        raise ValueError("match endpoint is not current expression_gene_id corresponding manual_new_gene")
    partner_role = "target" if hit == "base" else "base"
    new_endpoint = endpoint_fields(row, hit)
    partner = endpoint_fields(row, partner_role)
    endpoint_chrom = new_endpoint["chrom"]
    if (endpoint_chrom != normalize_chrom(candidate["chrom"])
            or parse_int(new_endpoint["start"], "gene start")
            != int(candidate["new_gene_start0"])
            or parse_int(new_endpoint["end"], "gene end")
            != int(candidate["new_gene_end0"])
            or new_endpoint["gene_strand"] != candidate["new_gene_strand"]):
        raise ValueError("Allelome new gene endpoint and five dogs annotation of chromosome, interval or strand mismatch")
    feelnc = set(split_values(candidate["new_gene_feelnc_biotypes"]))
    if "lncRNA" in feelnc and partner["endpoint_type"] == "known_protein_coding":
        direction = "new_lncRNA_to_pcGene"
        pair_orientation = "lncRNA_pcGene_class_orientation"
        direction_basis = (
            "posthoc_class_orientation_lncRNA_candidate_regulator_"
            "pcGene_putative_target"
        )
    else:
        direction = "undirected_pair"
        pair_orientation = "undirected_pair"
        direction_basis = (
            "base_target_are_canonical_endpoints_not_causal_direction"
        )

    expanded = {
        key: candidate[key]
        for key in CANDIDATE_COLUMNS[:-CANDIDATE_SUMMARY_COLUMN_COUNT]
    }
    expanded.update(
        {
            "animal_id": clean(row.get("animal_id")),
            "sample_tissue": clean(row.get("sample_tissue")),
            "canonical_pair_id": clean(row.get("canonical_pair_id")),
            "sample_pair_id": clean(row.get("sample_pair_id")),
            "evidence_stratum": clean(row.get("evidence_stratum")),
            "overlap_bp": clean(row.get("overlap_bp")),
            "either_endpoint_multigene": clean(row.get("either_endpoint_multigene")),
            "shared_variant_identity_checked": clean(row.get("shared_variant_identity_checked")),
            "source_link_file": clean(row.get("source_link_file")),
            "source_link_line": clean(row.get("source_link_line")),
            "pair_type": clean(row.get("pair_type")),
            "interval_gap_bp": clean(row.get("interval_gap_bp")),
            "mechanism": clean(row.get("mechanism")),
            "linkage_score": clean(row.get("linkage_score")),
            "new_gene_endpoint_role": hit,
            "allelome_endpoint_gene_id": new_endpoint["gene_id"],
            "allelome_endpoint_start0": new_endpoint["start"],
            "allelome_endpoint_end0": new_endpoint["end"],
            "allelome_endpoint_strand": new_endpoint["gene_strand"],
            "allelome_endpoint_chrom": endpoint_chrom,
            "allelome_endpoint_pos": "NA",
            "allelome_endpoint_ref": "NA",
            "allelome_endpoint_alt": "NA",
            "allelome_endpoint_name": new_endpoint["name"],
            "eqtl_gwas_snp_equals_allelome_endpoint": "not_evaluable_gene_endpoint",
            "eqtl_gwas_snp_distance_to_allelome_endpoint_bp": "NA",
            "partner_role": partner_role,
            "partner_record_type": partner["record_type"],
            "partner_gene_id": output_value(partner["gene_id"]),
            "partner_gene_symbol": output_value(partner["gene_symbol"]),
            "partner_gene_strand": output_value(partner["gene_strand"]),
            "partner_chrom": output_value(partner["chrom"]),
            "partner_start": output_value(partner["start"]),
            "partner_end": output_value(partner["end"]),
            "partner_name": output_value(partner["name"]),
            "partner_biotype": output_value(partner["biotype"]),
            "partner_endpoint_type": output_value(partner["endpoint_type"]),
            "candidate_direction": direction,
            "pair_orientation": pair_orientation,
            "causal_direction_supported": "no",
            "evidence_tier": "DIRECT_NEW_GENE_ENDPOINT_BLOCK_SUPPORT",
            "direction_basis": direction_basis,
            "integration_level": "block_level_new_gene_allelome_relation",
        }
    )
    return expanded


def scan_allelome_links(
    links_path: Path,
    summary_path: Path,
    candidates: Sequence[Mapping[str, object]],
) -> Tuple[List[Dict[str, object]], int]:
    summary = read_metric_table(summary_path, "\t", "metric")
    required_metrics = {
        "overall_status",
        "expression_gene_link_overlap_rows",
        "phase_scope_rule",
        "local_ps_rule",
    }
    missing_metrics = sorted(required_metrics - set(summary))
    if missing_metrics:
        raise ValueError(
            f"{summary_path} missing metric: {','.join(missing_metrics)}"
        )
    if summary["overall_status"] != "PASS":
        raise ValueError("6 result is not PASS")
    if summary["phase_scope_rule"] != "phase_scope_must_equal_chrN":
        raise ValueError("6 result is not all chromosome phase_scope")
    if summary["local_ps_rule"] != "provenance_only_not_used_for_filtering":
        raise ValueError("6 result local PS rule is not expected value")
    expected_rows = parse_int(
        summary["expression_gene_link_overlap_rows"],
        "6 result expression_gene_link_overlap_rows",
    )

    candidates_by_block: MutableMapping[
        str, List[Mapping[str, object]]
    ] = defaultdict(list)
    for candidate in candidates:
        candidates_by_block[str(candidate["block_id"])].append(candidate)

    expanded_rows: List[Dict[str, object]] = []
    seen_expanded_keys: Set[Tuple[str, str, str]] = set()
    observed_rows = 0
    with gzip.open(links_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            links_path, reader.fieldnames or [], ALLELOME_REQUIRED_COLUMNS
        )
        for line_number, row in enumerate(reader, start=2):
            observed_rows += 1
            block_id = clean(row.get("block_id"))
            block_candidates = candidates_by_block.get(block_id)
            if not block_candidates:
                continue
            if (
                clean(row.get("evidence_class"))
                != "expression_gene_matches_allelome_endpoint"
                or clean(row.get("evidence_level"))
                != "direct_feature_endpoint_evidence"
            ):
                raise ValueError(
                    f"{links_path}:{line_number} is not direct gene endpoint evidence"
                )
            if (
                clean(row.get("hit_endpoint_has_manual_new_gene")) != "yes"
                or clean(row.get("new_gene_evidence_scope"))
                != "direct_hit_endpoint_manual_new_gene"
            ):
                raise ValueError(
                    f"{links_path}:{line_number} no direct match Manual new gene"
                )
            if (
                "expression_gene_id_equals_allelome_gene_endpoint"
                not in clean(row.get("endpoint_mapping_rules"))
            ):
                raise ValueError(
                    f"{links_path}:{line_number} new gene mapping rule is not expected value"
                )
            for candidate in block_candidates:
                comparisons = {
                    "tissue": clean(row.get("eqtl_tissue")),
                    "expression_gene_id": clean(
                        row.get("expression_gene_id")
                    ),
                    "chrom": normalize_chrom(row.get("eqtl_chrom")),
                }
                expected = {
                    "tissue": str(candidate["tissue"]),
                    "expression_gene_id": str(candidate["new_gene_id"]),
                    "chrom": str(candidate["chrom"]),
                }
                if comparisons != expected:
                    raise ValueError(
                        f"{links_path}:{line_number} block join field mismatch: "
                        f"{comparisons!r} != {expected!r}"
                    )
                if clean(row.get("link_tissue")) != str(candidate["tissue"]):
                    raise ValueError(
                        f"{links_path}:{line_number} link_tissue mismatch"
                    )
                if normalize_chrom(row.get("link_chromosome")) != str(
                    candidate["chrom"]
                ):
                    raise ValueError(
                        f"{links_path}:{line_number} link chromosome mismatch"
                    )
                if clean(row.get("animal_id")) not in split_values(row.get("selected_samples")):
                    raise ValueError("Allelome animal not in this Expression–ASE block of selected_samples in")
                expanded = build_expanded_link(candidate, row)
                expanded_key = (
                    str(candidate["candidate_id"]),
                    clean(row.get("sample_pair_id")),
                    str(expanded["new_gene_endpoint_role"]),
                )
                if expanded_key in seen_expanded_keys:
                    raise ValueError(
                        f"{links_path}:{line_number} expand relationship duplicate key"
                    )
                seen_expanded_keys.add(expanded_key)
                expanded_rows.append(expanded)

    if observed_rows != expected_rows:
        raise ValueError(
            "6 number gene-link details row count and summary mismatch: "
            f"observed={observed_rows}, expected={expected_rows}"
        )
    return expanded_rows, observed_rows


def add_link_summaries(
    candidates: Sequence[Mapping[str, object]],
    expanded_rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    links_by_candidate: MutableMapping[
        str, List[Mapping[str, object]]
    ] = defaultdict(list)
    for row in expanded_rows:
        links_by_candidate[str(row["candidate_id"])].append(row)

    output: List[Dict[str, object]] = []
    for candidate in candidates:
        rows = links_by_candidate.get(str(candidate["candidate_id"]), [])
        partner_ids = [
            row["partner_gene_id"]
            for row in rows
            if row["partner_record_type"] in {"known_gene", "manual_new_gene"}
            and row["partner_gene_id"] != "NA"
        ]
        partner_symbols = [
            row["partner_gene_symbol"]
            for row in rows
            if row["partner_record_type"] in {"known_gene", "manual_new_gene"}
            and row["partner_gene_symbol"] != "NA"
        ]
        pcgene_rows = [
            row
            for row in rows
            if row["pair_orientation"] == "lncRNA_pcGene_class_orientation"
        ]
        result = dict(candidate)
        result.update(
            {
                "allelome_block_match": "yes" if rows else "no",
                "allelome_link_count": len(
                    {str(row["sample_pair_id"]) for row in rows}
                ),
                "supporting_dog_count": len(
                    {str(row["animal_id"]) for row in rows}
                ),
                "supporting_dogs": join_values(
                    row["animal_id"] for row in rows
                ),
                "pcgene_link_count": sum(
                    row["candidate_direction"] == "new_lncRNA_to_pcGene"
                    for row in rows
                ),
                "partner_gene_count": len(set(partner_ids)),
                "partner_gene_ids": join_values(partner_ids),
                "partner_gene_symbols": join_values(partner_symbols),
                "mechanisms": join_values(
                    row["mechanism"] for row in rows
                ),
                "pcgene_expanded_link_count": len(pcgene_rows),
                "unique_pcgene_count": len(
                    {
                        str(row["partner_gene_id"])
                        for row in pcgene_rows
                        if row["partner_gene_id"] != "NA"
                    }
                ),

                "exact_endpoint_link_count": "NA",
            }
        )
        output.append(result)
    return output


def chromosome_sort_key(value: object) -> Tuple[int, str]:
    chrom = normalize_chrom(value)
    match = re.fullmatch(r"chr(\d+)", chrom)
    if match:
        return int(match.group(1)), chrom
    return 10_000, chrom


def candidate_sort_key(
    row: Mapping[str, object]
) -> Tuple[int, str, int, str, str]:
    chrom_order, chrom = chromosome_sort_key(row["chrom"])
    return (
        chrom_order,
        chrom,
        int(row["eqtl_gwas_pos"]),
        str(row["tissue"]),
        str(row["new_gene_id"]),
    )


def expanded_sort_key(
    row: Mapping[str, object]
) -> Tuple[int, str, int, str, str, str, str]:
    base = candidate_sort_key(row)
    return base + (
        str(row["animal_id"]),
        str(row["canonical_pair_id"]),
    )


def write_tsv(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=columns,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=columns,
            extrasaction="ignore",
            lineterminator="\n",
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
        description=(
            "Integrate Allelome, eQTL fine-mapping, and raw GWAS evidence for five-dog new genes."
        )
    )
    parser.add_argument(
        "--eqtl-gwas-detail", type=Path, default=DEFAULT_EQTL_GWAS_DETAIL
    )
    parser.add_argument(
        "--eqtl-gwas-summary", type=Path, default=DEFAULT_EQTL_GWAS_SUMMARY
    )
    parser.add_argument(
        "--allelome-eqtl-links",
        type=Path,
        default=DEFAULT_ALLELOME_EQTL_LINKS,
    )
    parser.add_argument(
        "--allelome-eqtl-summary",
        type=Path,
        default=DEFAULT_ALLELOME_EQTL_SUMMARY,
    )
    parser.add_argument(
        "--new-gene-transcript-map",
        type=Path,
        default=DEFAULT_NEW_GENE_TRANSCRIPT_MAP,
    )
    parser.add_argument(
        "--new-gene-build-summary",
        type=Path,
        default=DEFAULT_NEW_GENE_BUILD_SUMMARY,
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--run-id", required=True)
    return parser.parse_args()


def main() -> int:
    global LOGGER
    args = parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.run_id):
        raise ValueError("--run-id can only contains letter, number, point, underscore and hyphen")

    args.output_root.mkdir(parents=True, exist_ok=True)
    run_dir = args.output_root / args.run_id
    if run_dir.exists():
        raise FileExistsError(f"run directory already exists, refuse to overwrite: {run_dir}")
    run_dir.mkdir()
    log_path = (
        run_dir
        / "1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas.log"
    )
    LOGGER = RunLogger(log_path)

    all_path = run_dir / "dog5_newgene_eqtl_raw_gwas_all_candidates.tsv"
    matched_path = (
        run_dir
        / "dog5_newgene_eqtl_raw_gwas_allelome_matched_blocks.tsv"
    )
    unmatched_path = (
        run_dir
        / "dog5_newgene_eqtl_raw_gwas_allelome_unmatched_blocks.tsv"
    )
    expanded_path = (
        run_dir
        / "dog5_newgene_eqtl_raw_gwas_allelome_expanded_links.tsv.gz"
    )
    pcgene_path = run_dir / "dog5_newgene_to_pcGene_candidate_links.tsv.gz"
    summary_path = (
        run_dir
        / "dog5_newgene_allelome_eqtl_raw_gwas_integration_summary.tsv"
    )
    manifest_path = (
        run_dir
        / "dog5_newgene_allelome_eqtl_raw_gwas_integration_manifest.tsv"
    )

    inputs = [
        (
            "eqtl_gwas_detail",
            args.eqtl_gwas_detail,
            "Exact overlap between result 2 eQTL fine-map SNPs with PPA>=0.1 and raw GWAS variants with p_lrt<=1e-4",
        ),
        (
            "eqtl_gwas_summary",
            args.eqtl_gwas_summary,
            "Input requirements, thresholds, and count summary for result 2",
        ),
        (
            "allelome_eqtl_links",
            args.allelome_eqtl_links,
            "Block-level direct relationships between expression genes and Allelome endpoints in result 6",
        ),
        (
            "allelome_eqtl_summary",
            args.allelome_eqtl_summary,
            "PASS status, chromosome-wide phasing rules, and count summary for result 6",
        ),
        (
            "new_gene_transcript_map",
            args.new_gene_transcript_map,
            "Five-dog DOG5_NEWGENE transcripts, gene spans, strands, and FEELnc types",
        ),
        (
            "new_gene_build_summary",
            args.new_gene_build_summary,
            "Five-dog new-gene construction provenance, ID prefix, and count requirements",
        ),
    ]

    try:
        log("[stage 1/6] check all absolute input path")
        for label, path, _ in inputs:
            require_file(path, label)
            log(f"[stage 1/6] {label}: {path.resolve()}")

        log("[stage 2/6] read five-dog-specific DOG5_NEWGENE metadata")
        metadata = load_new_gene_metadata(
            args.new_gene_transcript_map,
            args.new_gene_build_summary,
        )
        log(f"[stage 2/6] DOG5_NEWGENE count={len(metadata)}")

        log("[stage 3/6] read and filter 2 number eQTL fine-map with raw GWAS exact intersection")
        candidates, eqtl_gwas_input_rows = load_new_gene_candidates(
            args.eqtl_gwas_detail,
            args.eqtl_gwas_summary,
            metadata,
        )
        log(
            "[stage 3/6] 2 number input row count="
            f"{eqtl_gwas_input_rows}; new gene candidate row count={len(candidates)}; "
            "unique novel gene="
            f"{len({str(row['new_gene_id']) for row in candidates})}"
        )

        log("[stage 4/6] scan 6 number block level Allelome relationship and perform strict join")
        expanded_rows, allelome_input_rows = scan_allelome_links(
            args.allelome_eqtl_links,
            args.allelome_eqtl_summary,
            candidates,
        )
        log(
            "[stage 4/6] 6 number input row count="
            f"{allelome_input_rows}; expand match relationship={len(expanded_rows)}"
        )

        log("[stage 5/6] aggregation block result and stratify write")
        candidate_rows = add_link_summaries(candidates, expanded_rows)
        candidate_rows.sort(key=candidate_sort_key)
        matched_rows = [
            row
            for row in candidate_rows
            if row["allelome_block_match"] == "yes"
        ]
        unmatched_rows = [
            row
            for row in candidate_rows
            if row["allelome_block_match"] == "no"
        ]
        expanded_rows.sort(key=expanded_sort_key)
        pcgene_rows = [
            row
            for row in expanded_rows
            if row["pair_orientation"] == "lncRNA_pcGene_class_orientation"
        ]

        write_tsv(all_path, CANDIDATE_COLUMNS, candidate_rows)
        write_tsv(matched_path, CANDIDATE_COLUMNS, matched_rows)
        write_tsv(unmatched_path, CANDIDATE_COLUMNS, unmatched_rows)
        write_tsv_gz(
            expanded_path, EXPANDED_LINK_COLUMNS, expanded_rows
        )
        write_tsv_gz(pcgene_path, EXPANDED_LINK_COLUMNS, pcgene_rows)

        summary_rows = [
            {"metric": "run_id", "value": args.run_id},
            {"metric": "overall_status", "value": "PASS"},
            {
                "metric": "eqtl_raw_gwas_input_rows",
                "value": eqtl_gwas_input_rows,
            },
            {
                "metric": "new_gene_candidate_rows",
                "value": len(candidate_rows),
            },
            {
                "metric": "new_gene_candidate_unique_genes",
                "value": len(
                    {str(row["new_gene_id"]) for row in candidate_rows}
                ),
            },
            {
                "metric": "allelome_gene_link_input_rows",
                "value": allelome_input_rows,
            },
            {
                "metric": "matched_candidate_rows",
                "value": len(matched_rows),
            },
            {
                "metric": "matched_unique_new_genes",
                "value": len(
                    {str(row["new_gene_id"]) for row in matched_rows}
                ),
            },
            {
                "metric": "unmatched_candidate_rows",
                "value": len(unmatched_rows),
            },
            {
                "metric": "expanded_allelome_link_rows",
                "value": len(expanded_rows),
            },
            {
                "metric": "new_lncRNA_to_pcGene_rows",
                "value": len(pcgene_rows),
            },
            {
                "metric": "lncRNA_pcGene_class_oriented_rows",
                "value": len(pcgene_rows),
            },
            {
                "metric": "unique_lncRNA_pcGene_class_oriented_pairs",
                "value": len(
                    {
                        (
                            str(row["new_gene_id"]),
                            str(row["partner_gene_id"]),
                        )
                        for row in pcgene_rows
                    }
                ),
            },
            {
                "metric": "exact_lncRNA_pcGene_endpoint_rows",
                "value": "NA_gene_endpoints_not_single_variants",
            },
            {
                "metric": "exact_same_snp_endpoint_rows",
                "value": "NA_gene_endpoints_not_single_variants",
            },
            {
                "metric": "primary_join_rule",
                "value": (
                    "block_id+tissue+expression_gene_id+chrom"
                ),
            },
            {
                "metric": "link_count_unit",
                "value": "candidate_id+sample_pair_id+endpoint_role; canonical_pair_id_excludes_animal_and_tissue",
            },
            {
                "metric": "review_strata_policy",
                "value": "retain_all_links_and_copy_evidence_stratum; overlap_and_multigene_flags_are_not_causal_evidence",
            },
            {
                "metric": "gwas_filter_scope",
                "value": "raw_GWAS_p_lrt<=1e-4_pool; cumulative_flag_p_lrt<1e-5;four_variant_p_tiers_retained; credible_set_is_Expression_ASE_only",
            },
            {
                "metric": "variant_interpretation",
                "value": (
                    "Allelome_endpoint_is_gene_interval_not_single_SNP; "
                    "exact_SNP_identity_not_evaluable; inherited_model_membership_annotation_only; no_GWAS_PIP_or_CS_filter"
                ),
            },
            {
                "metric": "direction_interpretation",
                "value": (
                    "base_target_are_undirected; pair_orientation_is_"
                    "posthoc_class_orientation_only; causal_direction_"
                    "supported_is_no"
                ),
            },
            {
                "metric": "evidence_boundary",
                "value": (
                    "multilayer_candidate_convergence_not_colocalization_"
                    "or_causal_proof"
                ),
            },
        ]
        write_tsv(summary_path, ["metric", "value"], summary_rows)

        outputs = [
            (
                "all_new_gene_candidates",
                all_path,
                "All Manual new-gene candidates from result 2 and their Allelome match status",
            ),
            (
                "matched_blocks",
                matched_path,
                "New-gene candidates with block-level Allelome relationships",
            ),
            (
                "unmatched_blocks",
                unmatched_path,
                "New-gene candidates without block-level Allelome relationships, retained for provenance",
            ),
            (
                "expanded_links",
                expanded_path,
                "Expanded relationships per candidate, dog, and canonical Allelome link",
            ),
            (
                "new_lncRNA_to_pcGene_links",
                pcgene_path,
                "Subset oriented by post hoc endpoint categories, with a lncRNA new gene and a pcGene partner; this does not indicate causal direction.",
            ),
            (
                "integration_summary",
                summary_path,
                "Input reconciliation, match counts, and interpretation limits",
            ),
            (
                "run_log",
                log_path,
                "Stage-specific run logs",
            ),
        ]
        write_manifest(manifest_path, inputs, outputs)

        log(
            "[stage 5/6] all novel gene candidate="
            f"{len(candidate_rows)}; match block={len(matched_rows)}; "
            f"no match block={len(unmatched_rows)}; "
            f"expand links={len(expanded_rows)}; "
            f"lncRNA-pcGene posterior category orient={len(pcgene_rows)}"
        )
        log("[stage 6/6] all output completed")
        log(f"result directory: {run_dir.resolve()}")
        log("task completed / done")
        return 0
    except Exception as exc:
        log(f"task failed: {exc}")
        for line in traceback.format_exc().rstrip().splitlines():
            log(line)
        return 1
    finally:
        if LOGGER is not None:
            LOGGER.close()
            LOGGER = None


if __name__ == "__main__":
    raise SystemExit(main())
