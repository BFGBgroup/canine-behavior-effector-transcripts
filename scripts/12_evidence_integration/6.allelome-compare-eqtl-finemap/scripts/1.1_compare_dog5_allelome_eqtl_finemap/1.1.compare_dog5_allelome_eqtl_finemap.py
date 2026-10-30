#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import json
import csv
import gzip
import math
import re
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


DEFAULT_ALLELOME_LINKS = Path(
    '/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_links.annotated.tsv'
)
DEFAULT_ALLELOME_SUMMARY = Path(
    '/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_annotation.summary.tsv'
)
DEFAULT_CANDIDATE_BLOCKS = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs/20260821_181104/no_gwas_plasma_candidate_blocks.tsv'
)
DEFAULT_VALID_BLOCKS = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.5_no_gwas_build_expression_ase_regulatory_catalogue/20260827_201055/no_gwas_expression_ase_valid_block_catalogue.tsv'
)
DEFAULT_CREDIBLE_VARIANTS = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.6_no_gwas_export_expression_ase_catalogue_deliverables/20260827_224912/no_gwas_expression_ase_credible_set_variants.tsv.gz'
)
DEFAULT_FINAL_QC_SUMMARY = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.7_no_gwas_final_catalogue_qc_report/20260828_122551/no_gwas_expression_ase_final_qc_summary.tsv'
)
DEFAULT_NEW_GENE_GTF = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf'
)
DEFAULT_NEW_GENE_TRANSCRIPT_MAP = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv'
)
DEFAULT_NEW_GENE_BUILD_SUMMARY = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.summary.tsv'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap'
)

DEFAULT_PRIORITY_PPA = 0.1
INTERVAL_BIN_SIZE = 1_000_000
EXPECTED_EXPRESSION_MATRIX_SCOPE = (
    "reference_annotation_plus_dog5_manual_new_gene"
)
EXPECTED_FEATURE_ORIGIN_RULE = (
    "gene-=>reference_annotation_gene;"
    "DOG5_NEWGENE_=>dog5_manual_new_gene"
)
EXPECTED_NEW_GENE_PREFIX = "DOG5_NEWGENE_"
EXPECTED_NEW_GENE_SOURCE = "five-dogs"
EXPECTED_NEW_GENE_SELECTION_STAGE = (
    "2.1_transcript_first_preASE_candidate_assessment"
)
AUTOSOMES = {f"chr{number}" for number in range(1, 39)}
KNOWN_DOGS = {
    "Ace_396226",
    "CodRou_346627",
    "Crak_358336",
    "Dozer_386519",
    "SanMcg_349891",
}
KNOWN_TISSUES = {
    "cerebellum",
    "cortex_frontal",
    "cortex_occipital",
    "pituitary_gland",
}
VARIANT_RE = re.compile(
    r"(?:\(NA\))?((?:chr|NW_)[^:|\s]+):(\d+)_([A-Za-z]+)>([A-Za-z]+)"
)
GTF_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')

LINK_REQUIRED_COLUMNS = {
    "animal_id",
    "tissue",
    "source_link_file",
    "source_link_line",
    "evidence_stratum",
    "overlap_bp",
    "either_endpoint_multigene",
    "shared_variant_identity_checked",
    "interval_gap_bp",
    "pair_type",
    "mechanism",
    "linkage_score",
    "chr_base",
    "start_base",
    "end_base",
    "name_base",
    "base_record_type",
    "base_gene_id",
    "base_gene_symbol",
    "base_gene_strand",
    "base_annotation_source",
    "base_gene_biotype",
    "base_endpoint_type",
    "chr_target",
    "start_target",
    "end_target",
    "name_target",
    "target_record_type",
    "target_gene_id",
    "target_gene_symbol",
    "target_gene_strand",
    "target_annotation_source",
    "target_gene_biotype",
    "target_endpoint_type",
}
CANDIDATE_BLOCK_REQUIRED_COLUMNS = {
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
    "selected_samples",
    "selected_phase_scopes",
    "plasma_input_allowed",
}
VALID_BLOCK_REQUIRED_COLUMNS = {
    "block_id",
    "feature_origin",
    "expression_matrix_scope",
    "interpretation_note",
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
    "result_status",
    "result_interpretable",
    "ppa_status",
    "ppa_sum",
    "credible_set_size",
    "credible_set_ppa_sum",
    "max_ppa",
    "top_ppa_tie_count",
    "top_marker_is_unique",
    "top_marker_id",
    "top_marker_pos",
}
CREDIBLE_REQUIRED_COLUMNS = {
    "catalogue_role",
    "feature_origin",
    "expression_matrix_scope",
    "is_top_variant",
    "is_credible_set_member",
    "ppa_tier",
    "ppa_status",
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
}
NEW_GENE_TRANSCRIPT_MAP_REQUIRED_COLUMNS = {
    "final_gene_id",
    "final_transcript_id",
    "source_dataset",
    "selection_stage",
    "original_gene_id",
    "original_transcript_id",
    "chr",
    "strand",
    "gene_start_0based",
    "gene_end_0based",
    "transcript_start_0based",
    "transcript_end_0based",
    "class_code",
    "novelty_category",
    "gene_locus_filter_reason",
}

BLOCK_COLUMNS = [
    "block_id",
    "feature_origin",
    "expression_matrix_scope",
    "block_interpretation_note",
    "eqtl_tissue",
    "feature_id",
    "feature_primary_id",
    "feature_label",
    "feature_class",
    "expression_gene_id",
    "expression_source",
    "eqtl_chrom",
    "selected_sample_count",
    "selected_samples",
    "common_marker_count",
    "total_reads",
    "credible_set_size",
    "credible_set_ppa_sum",
    "max_ppa",
    "top_marker_id",
    "top_marker_pos",
    "top_ppa_tie_count",
    "top_marker_is_unique",
    "block_ppa_status",
]
VARIANT_COLUMNS = [
    "eqtl_pos",
    "marker_id",
    "ref",
    "alt",
    "unordered_ref_alt",
    "ppa",
    "ppa_tier",
    "variant_ppa_status",
    "marker_rank_by_ppa",
    "is_top_variant",
    "catalogue_role",
    "phase_scope",
    "local_ps_by_sample",
    "variant_interpretation_note",
]
LINK_COLUMNS = [
    "animal_id",
    "link_tissue",
    "sample_tissue",
    "link_chromosome",
    "canonical_pair_id",
    "sample_pair_id", "evidence_stratum", "overlap_bp", "either_endpoint_multigene",
    "shared_variant_identity_checked", "source_link_file", "source_link_line",
    "pair_type",
    "interval_gap_bp",
    "mechanism",
    "linkage_score",
]
ENDPOINT_SUFFIXES = [
    "record_type",
    "gene_id",
    "gene_symbol",
    "gene_strand",
    "chr",
    "start",
    "end",
    "name",
    "annotation_source",
    "gene_biotype",
    "endpoint_type",
    "manual_new_gene_hit",
    "manual_new_gene_hit_n",
    "manual_new_gene_ids",
    "manual_new_gene_strands",
    "manual_new_gene_class_codes",
    "manual_new_gene_novelty_categories",
    "manual_new_gene_candidate_categories",
    "manual_new_gene_dataset_support",
    "manual_new_gene_merge_statuses",
]
ENDPOINT_COLUMNS = [
    f"{prefix}_{suffix}"
    for prefix in ("base", "target")
    for suffix in ENDPOINT_SUFFIXES
]
EVIDENCE_COLUMNS = [
    "endpoint_hit",
    "endpoint_mapping_rules",
    "evidence_class",
    "evidence_level",
    "hit_endpoint_gene_id",
    "hit_endpoint_gene_symbol",
    "hit_endpoint_has_manual_new_gene",
    "link_contains_manual_new_gene",
    "new_gene_evidence_scope",
]
GENE_OVERLAP_COLUMNS = BLOCK_COLUMNS + LINK_COLUMNS + ENDPOINT_COLUMNS + EVIDENCE_COLUMNS
VARIANT_OVERLAP_COLUMNS = (
    BLOCK_COLUMNS
    + VARIANT_COLUMNS
    + LINK_COLUMNS
    + ENDPOINT_COLUMNS
    + EVIDENCE_COLUMNS
)
BLOCK_SUMMARY_COLUMNS = BLOCK_COLUMNS + [
    "credible_variant_rows_observed",
    "priority_variant_rows_ppa_ge_threshold",
    "expression_gene_link_overlap_rows",
    "exact_variant_endpoint_overlap_rows",
    "priority_gene_span_context_rows",
    "priority_ppa_threshold",
]


def set_large_csv_fields() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


@dataclass(frozen=True, slots=True)
class Endpoint:
    role: str
    record_type: str
    gene_id: str
    gene_symbol: str
    gene_strand: str
    chrom: str
    start: int
    end: int
    name: str
    annotation_source: str
    gene_biotype: str
    endpoint_type: str
    variant_pos: int | None
    variant_ref: str
    variant_alt: str
    manual_new_gene_hit: str
    manual_new_gene_hit_n: int
    manual_new_gene_ids: str
    manual_new_gene_strands: str
    manual_new_gene_class_codes: str
    manual_new_gene_novelty_categories: str
    manual_new_gene_candidate_categories: str
    manual_new_gene_dataset_support: str
    manual_new_gene_merge_statuses: str


@dataclass(frozen=True, slots=True)
class Link:
    animal_id: str
    tissue: str
    sample_tissue: str
    chromosome: str
    canonical_pair_id: str
    sample_pair_id: str
    evidence_stratum: str
    overlap_bp: str
    either_endpoint_multigene: str
    shared_variant_identity_checked: str
    source_link_file: str
    source_link_line: str
    pair_type: str
    interval_gap_bp: str
    mechanism: str
    linkage_score: str
    base: Endpoint
    target: Endpoint


@dataclass(frozen=True, slots=True)
class EndpointRef:
    link: Link
    role: str

    @property
    def endpoint(self) -> Endpoint:
        return self.link.base if self.role == "base" else self.link.target


def log(message: str) -> None:

    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S} ] {message}", flush=True)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def normalize_chrom(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    if text.startswith("chr") or text.startswith("NW_"):
        return text
    if text.isdigit() or text in {"X", "Y", "M", "MT"}:
        return f"chr{text}"
    return text


def parse_int(value: object, label: str) -> int:
    try:
        return int(float(clean(value)))
    except ValueError as exc:
        raise ValueError(f"{label} is not an integer: {value!r}") from exc


def parse_float(value: object, label: str) -> float:
    try:
        result = float(clean(value))
    except ValueError as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite numeric value: {value!r}")
    return result


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


def allele_key(ref: object, alt: object) -> Tuple[str, str]:
    values = sorted([clean(ref).upper(), clean(alt).upper()])
    if len(values) != 2 or not values[0] or not values[1] or values[0] == values[1]:
        raise ValueError(f"REF/ALT invalid: {ref!r}, {alt!r}")
    return values[0], values[1]


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
    return "|".join(result)


def join_strands(values: Iterable[object]) -> str:
    result: List[str] = []
    seen: Set[str] = set()
    for value in values:
        strand = clean(value)
        if strand in {"+", "-"} and strand not in seen:
            seen.add(strand)
            result.append(strand)
    return "|".join(result)


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} does not exist or is empty: {path}")


def require_columns(path: Path, observed: Sequence[str], required: Set[str]) -> None:
    missing = sorted(required - set(observed))
    if missing:
        raise ValueError(f"{path} missing required fields: {','.join(missing)}")


def read_key_value(path: Path, key_column: str) -> Dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], {key_column, "value"})
        result: Dict[str, str] = {}
        for row in reader:
            key = clean(row.get(key_column))
            if not key or key in result:
                raise ValueError(f"{path} contains empty or duplicate {key_column}: {key!r}")
            result[key] = clean(row.get("value"))
    return result


def parse_gtf_attributes(text: str) -> Dict[str, str]:
    return {
        match.group(1): match.group(2)
        for match in GTF_ATTRIBUTE_RE.finditer(text)
    }


def candidate_category(class_codes: Sequence[str], novelty: Sequence[str]) -> str:
    code_set = set(class_codes)
    novelty_set = set(novelty)
    categories: List[str] = []
    if "u" in code_set or "intergenic_novel_transcript" in novelty_set:
        categories.append("intergenic_new_gene_candidate")
    if "x" in code_set or any("antisense" in item for item in novelty_set):
        categories.append("antisense_candidate")
    if "i" in code_set or any("intronic" in item for item in novelty_set):
        categories.append("intronic_candidate")
    if "j" in code_set or "known_gene_novel_isoform" in novelty_set:
        categories.append("known_gene_novel_isoform")
    return "|".join(categories or ["other_candidate_transcript_locus"])


def validate_allelome_summary(path: Path) -> Dict[str, str]:

    summary = read_key_value(path, "metric")
    parse_int(summary["annotated_rows"], "annotated_rows")
    return summary


def validate_final_qc(
    path: Path,
    credible_path: Path,
    valid_blocks_path: Path,
) -> Dict[str, str]:
    summary = read_key_value(path, "key")
    expected = {
        "final_qc_status": "PASS",
        "issue_count": "0",
        "source_8.6_expression_matrix_scope": EXPECTED_EXPRESSION_MATRIX_SCOPE,
        "source_8.6_feature_origin_rule": EXPECTED_FEATURE_ORIGIN_RULE,
        "source_8.6_association_fdr_screening": "not_performed",
        "source_8.6_null_causal_model": "not_included",
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise ValueError(
                f"eQTL final QC {key} unexpected: "
                f"observed={summary.get(key)!r}, expected={value!r}"
            )
    valid_blocks = parse_int(
        summary.get("source_8.6_valid_blocks"),
        "eQTL final QC source_8.6_valid_blocks",
    )
    credible_rows = parse_int(
        summary.get("source_8.6_credible_set_variant_rows"),
        "eQTL final QC source_8.6_credible_set_variant_rows",
    )
    observed_credible_rows = parse_int(
        summary.get("observed_rows_credible_set_variants"),
        "eQTL final QC observed_rows_credible_set_variants",
    )
    if valid_blocks < 1 or credible_rows < 1:
        raise ValueError("eQTL final QC recorded valid block or credible set row count is not positive number")
    if observed_credible_rows != credible_rows:
        raise ValueError(
            "eQTL final QC of credible set input row count and actual audit row counts differ: "
            f"source={credible_rows}, observed={observed_credible_rows}"
        )
    export_dir = Path(summary["source_export_dir"]).resolve()
    if credible_path.resolve().parent != export_dir:
        raise ValueError(
            "credible set file does not belong to final QC recorded 8.6 output directory: "
            f"credible={credible_path.parent}, qc={export_dir}"
        )
    recorded_valid = Path(summary["source_8.6_source_valid_block_catalogue"]).resolve()
    if recorded_valid != valid_blocks_path.resolve():
        raise ValueError(
            "valid block file and final QC record mismatch: "
            f"input={valid_blocks_path}, qc={recorded_valid}"
        )
    return summary


def load_new_gene_loci(
    gtf_path: Path,
    transcript_map_path: Path,
    summary_path: Path,
) -> Tuple[Dict[str, int], Mapping[Tuple[str, int], Sequence[Dict[str, object]]]]:
    summary = read_key_value(summary_path, "metric")
    expected_loci = parse_int(
        summary.get("candidate_gene_count"),
        "new-gene summary candidate_gene_count",
    )
    expected_transcripts = parse_int(
        summary.get("candidate_transcript_count"),
        "new-gene summary candidate_transcript_count",
    )
    expected_rows = parse_int(
        summary.get("appended_gtf_row_count"),
        "new-gene summary appended_gtf_row_count",
    )
    expected_summary_values = {
        "source_dataset": EXPECTED_NEW_GENE_SOURCE,
        "selection_stage": EXPECTED_NEW_GENE_SELECTION_STAGE,
        "final_gene_id_prefix": EXPECTED_NEW_GENE_PREFIX,
    }
    for key, expected in expected_summary_values.items():
        if clean(summary.get(key)) != expected:
            raise ValueError(
                f"five-only new-gene summary {key} unexpected: "
                f"observed={summary.get(key)!r}, expected={expected!r}"
            )
    for key, observed in (
        ("new_only_gtf", gtf_path),
        ("transcript_map", transcript_map_path),
    ):
        recorded = Path(clean(summary.get(key))).resolve()
        if recorded != observed.resolve():
            raise ValueError(
                f"five-only new-gene summary {key} does not match the input: "
                f"summary={recorded}, input={observed}"
            )

    metadata: Dict[str, Dict[str, object]] = {}
    transcript_ids: Set[str] = set()
    observed_map_rows = 0
    with transcript_map_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            transcript_map_path,
            reader.fieldnames or [],
            NEW_GENE_TRANSCRIPT_MAP_REQUIRED_COLUMNS,
        )
        for line_number, row in enumerate(reader, start=2):
            observed_map_rows += 1
            locus_id = clean(row.get("final_gene_id"))
            transcript_id = clean(row.get("final_transcript_id"))
            if not locus_id.startswith(EXPECTED_NEW_GENE_PREFIX):
                raise ValueError(
                    f"{transcript_map_path}:{line_number} final_gene_id prefix unexpected"
                )
            if not transcript_id or transcript_id in transcript_ids:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} final_transcript_id missing or duplicate"
                )
            transcript_ids.add(transcript_id)
            if (
                clean(row.get("source_dataset")) != EXPECTED_NEW_GENE_SOURCE
                or clean(row.get("selection_stage"))
                != EXPECTED_NEW_GENE_SELECTION_STAGE
            ):
                raise ValueError(
                    f"{transcript_map_path}:{line_number} data source or filter stage unexpected"
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
            transcript_start0 = parse_int(
                row.get("transcript_start_0based"),
                f"{transcript_map_path}:{line_number} transcript_start_0based",
            )
            transcript_end0 = parse_int(
                row.get("transcript_end_0based"),
                f"{transcript_map_path}:{line_number} transcript_end_0based",
            )
            if strand not in {"+", "-"} or start0 < 0 or end0 <= start0:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} new gene interval or strand invalid"
                )
            if not start0 <= transcript_start0 < transcript_end0 <= end0:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} transcript interval not in gene span within"
                )
            original_gene_id = clean(row.get("original_gene_id"))
            class_code = clean(row.get("class_code"))
            novelty = clean(row.get("novelty_category"))
            filter_reason = clean(row.get("gene_locus_filter_reason"))
            if not original_gene_id or not class_code or not novelty or not filter_reason:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} five dogs candidate provenance field missing"
                )
            if locus_id not in metadata:
                metadata[locus_id] = {
                    "unified_locus_id": locus_id,
                    "chrom": chrom,
                    "strand": strand,
                    "start0": start0,
                    "end0": end0,
                    "original_gene_id": original_gene_id,
                    "merge_status": "five_only_not_cross_resource_merged",
                    "dataset_support": EXPECTED_NEW_GENE_SOURCE,
                    "class_codes": [],
                    "novelty": [],
                }
            meta = metadata[locus_id]
            expected_identity = (
                str(meta["chrom"]),
                str(meta["strand"]),
                int(meta["start0"]),
                int(meta["end0"]),
                str(meta["original_gene_id"]),
            )
            observed_identity = (
                chrom,
                strand,
                start0,
                end0,
                original_gene_id,
            )
            if observed_identity != expected_identity:
                raise ValueError(
                    f"{transcript_map_path}:{line_number} the same final_gene_id metadata mismatch"
                )
            for key, value in (
                ("class_codes", class_code),
                ("novelty", novelty),
            ):
                values = meta[key]
                if not isinstance(values, list):
                    raise TypeError(f"{locus_id} of {key} aggregation container is not list")
                if value not in values:
                    values.append(value)

    if observed_map_rows != expected_transcripts or len(metadata) != expected_loci:
        raise ValueError(
            "five-only transcript map count with build summary mismatch: "
            f"rows={observed_map_rows}/{expected_transcripts}, "
            f"genes={len(metadata)}/{expected_loci}"
        )
    for meta in metadata.values():
        class_codes = list(meta["class_codes"])
        novelty = list(meta["novelty"])
        meta["candidate_category"] = candidate_category(class_codes, novelty)

    transcript_spans: Dict[str, List[object]] = {}
    observed_gtf_transcripts: Set[str] = set()
    observed_rows = 0
    observed_transcripts = 0
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gtf_path}:{line_number} GTF column count is not 9: {len(fields)}"
                )
            observed_rows += 1
            chrom, _, feature, start_text, end_text, _, strand, _, attributes = fields
            if feature not in {"transcript", "exon"}:
                raise ValueError(
                    f"{gtf_path}:{line_number} encountered unsupported of feature: {feature}"
                )
            attrs = parse_gtf_attributes(attributes)
            locus_id = clean(attrs.get("gene_id"))
            if not locus_id or locus_id not in metadata:
                raise ValueError(
                    f"{gtf_path}:{line_number} gene_id not in five-only transcript map"
                )
            if (
                clean(attrs.get("source_dataset")) != EXPECTED_NEW_GENE_SOURCE
                or clean(attrs.get("selection_stage"))
                != EXPECTED_NEW_GENE_SELECTION_STAGE
            ):
                raise ValueError(
                    f"{gtf_path}:{line_number} data source or filter stage unexpected"
                )
            normalized_chrom = normalize_chrom(chrom)
            meta = metadata[locus_id]
            if normalized_chrom != meta["chrom"] or strand != meta["strand"]:
                raise ValueError(
                    f"{gtf_path}:{line_number} with transcript map chromosome or strand mismatch"
                )
            if feature != "transcript":
                continue
            observed_transcripts += 1
            transcript_id = clean(attrs.get("transcript_id"))
            if (
                not transcript_id
                or transcript_id not in transcript_ids
                or transcript_id in observed_gtf_transcripts
            ):
                raise ValueError(
                    f"{gtf_path}:{line_number} transcript_id missing, duplicate or not in map"
                )
            observed_gtf_transcripts.add(transcript_id)
            start0 = parse_int(start_text, "GTF start") - 1
            end0 = parse_int(end_text, "GTF end")
            if locus_id not in transcript_spans:
                transcript_spans[locus_id] = [normalized_chrom, strand, start0, end0]
            else:
                transcript_spans[locus_id][2] = min(
                    int(transcript_spans[locus_id][2]), start0
                )
                transcript_spans[locus_id][3] = max(
                    int(transcript_spans[locus_id][3]), end0
                )

    if (
        observed_rows != expected_rows
        or observed_transcripts != expected_transcripts
        or observed_gtf_transcripts != transcript_ids
    ):
        raise ValueError(
            "five-only new-only GTF count or transcript set and build summary mismatch: "
            f"rows={observed_rows}/{expected_rows}, "
            f"transcripts={observed_transcripts}/{expected_transcripts}"
        )
    if set(transcript_spans) != set(metadata):
        raise ValueError("five-only new-only GTF with transcript map of gene ID set mismatch")

    bins: MutableMapping[Tuple[str, int], List[Dict[str, object]]] = defaultdict(list)
    autosomal_loci = 0
    for locus_id, meta in metadata.items():
        _, _, start0, end0 = transcript_spans[locus_id]
        if int(start0) != meta["start0"] or int(end0) != meta["end0"]:
            raise ValueError(
                f"{locus_id} of GTF transcript union with five-only map boundary mismatch"
            )
        chrom = str(meta["chrom"])
        if chrom not in AUTOSOMES:
            continue
        autosomal_loci += 1
        for bin_id in range(
            int(meta["start0"]) // INTERVAL_BIN_SIZE,
            (int(meta["end0"]) - 1) // INTERVAL_BIN_SIZE + 1,
        ):
            bins[(chrom, bin_id)].append(meta)
    return {
        "total_loci": expected_loci,
        "autosomal_loci": autosomal_loci,
        "non_autosomal_loci": expected_loci - autosomal_loci,
        "transcripts": expected_transcripts,
        "gtf_rows": expected_rows,
    }, bins


def find_new_gene_hits(
    bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    chrom: str,
    pos1: int,
) -> List[Dict[str, object]]:
    point0 = pos1 - 1
    return sorted(
        [
            row
            for row in bins.get((chrom, point0 // INTERVAL_BIN_SIZE), [])
            if int(row["start0"]) <= point0 < int(row["end0"])
        ],
        key=lambda row: clean(row["unified_locus_id"]),
    )


def summarize_new_gene_hits(
    hits: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    return {
        "manual_new_gene_hit": "yes" if hits else "no",
        "manual_new_gene_hit_n": len(hits),
        "manual_new_gene_ids": join_values(
            row["unified_locus_id"] for row in hits
        ),
        "manual_new_gene_strands": join_strands(row["strand"] for row in hits),
        "manual_new_gene_class_codes": join_values(
            code for row in hits for code in row["class_codes"]
        ),
        "manual_new_gene_novelty_categories": join_values(
            item for row in hits for item in row["novelty"]
        ),
        "manual_new_gene_candidate_categories": join_values(
            row["candidate_category"] for row in hits
        ),
        "manual_new_gene_dataset_support": join_values(
            row["dataset_support"] for row in hits
        ),
        "manual_new_gene_merge_statuses": join_values(
            row["merge_status"] for row in hits
        ),
    }


def load_candidate_blocks(path: Path) -> Dict[str, Dict[str, object]]:
    result: Dict[str, Dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(
            path,
            reader.fieldnames or [],
            CANDIDATE_BLOCK_REQUIRED_COLUMNS,
        )
        for line_number, row in enumerate(reader, start=2):
            block_id = clean(row.get("block_id"))
            if not block_id or block_id in result:
                raise ValueError(f"{path}:{line_number} block_id missing or duplicate")
            samples = tuple(
                item for item in clean(row.get("selected_samples")).split(",") if item
            )
            selected_count = parse_int(
                row.get("selected_sample_count"),
                f"{path}:{line_number} selected_sample_count",
            )
            if len(samples) != selected_count or len(set(samples)) != len(samples):
                raise ValueError(
                    f"{path}:{line_number} selected_samples unexpected count or uniqueness"
                )
            if not set(samples).issubset(KNOWN_DOGS):
                raise ValueError(
                    f"{path}:{line_number} encountered unknown five dogs sample: "
                    f"{sorted(set(samples) - KNOWN_DOGS)}"
                )
            chrom = normalize_chrom(row.get("chrom"))
            tissue = clean(row.get("tissue"))
            if chrom not in AUTOSOMES or tissue not in KNOWN_TISSUES:
                raise ValueError(
                    f"{path}:{line_number} chromosome or tissue invalid: {chrom}, {tissue}"
                )
            scopes: Dict[str, str] = {}
            for token in clean(row.get("selected_phase_scopes")).split(","):
                sample, separator, scope = token.partition(":")
                if not separator or sample in scopes:
                    raise ValueError(
                        f"{path}:{line_number} selected_phase_scopes format unexpected"
                    )
                scopes[sample] = normalize_chrom(scope)
            if set(scopes) != set(samples) or any(
                scope != chrom for scope in scopes.values()
            ):
                raise ValueError(
                    f"{path}:{line_number} is not per sample all chromosome phase_scope"
                )
            result[block_id] = {
                "block_id": block_id,
                "tissue": tissue,
                "feature_id": clean(row.get("feature_id")),
                "feature_primary_id": clean(row.get("feature_primary_id")),
                "feature_label": clean(row.get("feature_label")),
                "feature_class": clean(row.get("feature_class")),
                "expression_gene_id": clean(row.get("expression_gene_id")),
                "expression_source": clean(row.get("expression_source")),
                "chrom": chrom,
                "selected_sample_count": selected_count,
                "selected_samples_tuple": samples,
                "selected_samples": ",".join(samples),
                "plasma_input_allowed": clean(row.get("plasma_input_allowed")),
            }
    return result


def load_valid_blocks(
    path: Path,
    candidate_blocks: Mapping[str, Dict[str, object]],
    expected_valid_blocks: int,
) -> Dict[str, Dict[str, object]]:
    result: Dict[str, Dict[str, object]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], VALID_BLOCK_REQUIRED_COLUMNS)
        for line_number, row in enumerate(reader, start=2):
            block_id = clean(row.get("block_id"))
            if not block_id or block_id in result:
                raise ValueError(f"{path}:{line_number} block_id missing or duplicate")
            if block_id not in candidate_blocks:
                raise ValueError(
                    f"{path}:{line_number} block_id not in 8.2 candidate blocks in"
                )
            candidate = candidate_blocks[block_id]
            if candidate["plasma_input_allowed"] != "1":
                raise ValueError(
                    f"{path}:{line_number} valid block in 8.2 not permitted to enter PLASMA"
                )
            if (
                clean(row.get("result_status")) != "valid"
                or clean(row.get("result_interpretable")) != "1"
            ):
                raise ValueError(f"{path}:{line_number} is not interpretable valid block")
            expression_gene_id = clean(row.get("expression_gene_id"))
            if expression_gene_id.startswith(EXPECTED_NEW_GENE_PREFIX):
                expected_feature_origin = "dog5_manual_new_gene"
            elif expression_gene_id.startswith("gene-"):
                expected_feature_origin = "reference_annotation_gene"
            else:
                raise ValueError(
                    f"{path}:{line_number} expression_gene_id prefix does not belong to current annotation"
                )
            feature_origin = clean(row.get("feature_origin"))
            expression_matrix_scope = clean(row.get("expression_matrix_scope"))
            interpretation_note = clean(row.get("interpretation_note"))
            if feature_origin != expected_feature_origin:
                raise ValueError(
                    f"{path}:{line_number} feature_origin with expression_gene_id mismatch"
                )
            if expression_matrix_scope != EXPECTED_EXPRESSION_MATRIX_SCOPE:
                raise ValueError(
                    f"{path}:{line_number} expression_matrix_scope is not current five-dog-specific span"
                )
            if not interpretation_note:
                raise ValueError(f"{path}:{line_number} interpretation_note missing")
            comparisons = {
                "tissue": clean(row.get("tissue")),
                "feature_id": clean(row.get("feature_id")),
                "feature_primary_id": clean(row.get("feature_primary_id")),
                "feature_label": clean(row.get("feature_label")),
                "feature_class": clean(row.get("feature_class")),
                "expression_gene_id": clean(row.get("expression_gene_id")),
                "expression_source": clean(row.get("expression_source")),
                "chrom": normalize_chrom(row.get("chrom")),
                "selected_sample_count": parse_int(
                    row.get("selected_sample_count"),
                    f"{path}:{line_number} selected_sample_count",
                ),
            }
            for key, observed in comparisons.items():
                if observed != candidate[key]:
                    raise ValueError(
                        f"{path}:{line_number} with 8.2 of {key} mismatch: "
                        f"{observed!r} != {candidate[key]!r}"
                    )
            block = dict(candidate)
            ppa_status = clean(row.get("ppa_status"))
            top_marker_is_unique = clean(row.get("top_marker_is_unique"))
            if ppa_status not in {"finite", "finite_boundary_adjusted"}:
                raise ValueError(
                    f"{path}:{line_number} ppa_status invalid: {ppa_status!r}"
                )
            if top_marker_is_unique not in {"0", "1"}:
                raise ValueError(
                    f"{path}:{line_number} top_marker_is_unique invalid"
                )
            block.update(
                {
                    "feature_origin": feature_origin,
                    "expression_matrix_scope": expression_matrix_scope,
                    "interpretation_note": interpretation_note,
                    "common_marker_count": parse_int(
                        row.get("common_marker_count"),
                        f"{path}:{line_number} common_marker_count",
                    ),
                    "total_reads": parse_int(
                        row.get("total_reads"),
                        f"{path}:{line_number} total_reads",
                    ),
                    "ppa_status": ppa_status,
                    "ppa_sum": parse_float(
                        row.get("ppa_sum"),
                        f"{path}:{line_number} ppa_sum",
                    ),
                    "credible_set_size": parse_int(
                        row.get("credible_set_size"),
                        f"{path}:{line_number} credible_set_size",
                    ),
                    "credible_set_ppa_sum": parse_float(
                        row.get("credible_set_ppa_sum"),
                        f"{path}:{line_number} credible_set_ppa_sum",
                    ),
                    "max_ppa": parse_float(
                        row.get("max_ppa"),
                        f"{path}:{line_number} max_ppa",
                    ),
                    "top_ppa_tie_count": parse_int(
                        row.get("top_ppa_tie_count"),
                        f"{path}:{line_number} top_ppa_tie_count",
                    ),
                    "top_marker_is_unique": top_marker_is_unique,
                    "top_marker_id": clean(row.get("top_marker_id")),
                    "top_marker_pos": parse_int(
                        row.get("top_marker_pos"),
                        f"{path}:{line_number} top_marker_pos",
                    ),
                }
            )
            result[block_id] = block
    if len(result) != expected_valid_blocks:
        raise ValueError(
            "valid block unexpected count: "
            f"observed={len(result)}, expected={expected_valid_blocks}"
        )
    return result


def build_gene_block_index(
    blocks: Mapping[str, Dict[str, object]],
) -> Mapping[Tuple[str, str, str, str], Sequence[str]]:
    index: MutableMapping[Tuple[str, str, str, str], List[str]] = defaultdict(list)
    for block_id, block in blocks.items():
        for animal_id in block["selected_samples_tuple"]:
            key = (
                str(animal_id),
                str(block["tissue"]),
                str(block["chrom"]),
                str(block["expression_gene_id"]),
            )
            index[key].append(block_id)
    return index


def parse_unannotated_variant(
    role: str,
    row: Mapping[str, str],
) -> Tuple[int, str, str]:
    text = clean(row.get(f"name_{role}")) or clean(
        row.get(f"{role}_gene_id")
    )
    match = VARIANT_RE.search(text)
    if match is None:
        raise ValueError(f"cannot parse the unannotated endpoint variant: {text!r}")
    chrom = normalize_chrom(match.group(1))
    pos = int(match.group(2))
    ref = match.group(3).upper()
    alt = match.group(4).upper()
    if chrom != normalize_chrom(row.get(f"chr_{role}")):
        raise ValueError(f"unannotated endpoint chromosome does not match its name: {text}")
    if (
        pos - 1 != parse_int(row.get(f"start_{role}"), f"{role} start")
        or pos != parse_int(row.get(f"end_{role}"), f"{role} end")
    ):
        raise ValueError(f"unannotated endpoint position does not match its name: {text}")
    allele_key(ref, alt)
    return pos, ref, alt


def endpoint_from_row(
    role: str,
    row: Mapping[str, str],
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
) -> Endpoint:
    record_type = clean(row.get(f"{role}_record_type"))
    chrom = normalize_chrom(row.get(f"chr_{role}"))
    start = parse_int(row.get(f"start_{role}"), f"{role} start")
    end = parse_int(row.get(f"end_{role}"), f"{role} end")
    if chrom not in AUTOSOMES:
        raise ValueError(f"Allelome endpoint is not chr1-chr38: {chrom}")

    variant_pos: int | None = None
    variant_ref = ""
    variant_alt = ""
    new_gene_hits: List[Dict[str, object]] = []
    if record_type in {"known_gene", "manual_new_gene"}:
        if start < 0 or end <= start:
            raise ValueError(f"annotated_gene endpoint interval invalid: {chrom}:{start}-{end}")
        if not clean(row.get(f"{role}_gene_id")):
            raise ValueError("annotated_gene endpoint missing gene_id")
    elif record_type == "unannotated_variant":
        variant_pos, variant_ref, variant_alt = parse_unannotated_variant(role, row)
        new_gene_hits = find_new_gene_hits(new_gene_bins, chrom, variant_pos)
    else:
        raise ValueError(f"unsupported of Allelome record_type: {record_type!r}")
    if record_type == "manual_new_gene":

        new_gene_hits = [hit for hit in find_new_gene_hits(new_gene_bins, chrom, start + 1)
                         if hit["unified_locus_id"] == clean(row.get(f"{role}_gene_id"))]
        if not new_gene_hits:
            raise ValueError("manual_new_gene endpoint ID not in five dogs GTF in")
    annotation = summarize_new_gene_hits(new_gene_hits)
    return Endpoint(
        role=role,
        record_type=record_type,
        gene_id=clean(row.get(f"{role}_gene_id")),
        gene_symbol=clean(row.get(f"{role}_gene_symbol")),
        gene_strand=clean(row.get(f"{role}_gene_strand")),
        chrom=chrom,
        start=start,
        end=end,
        name=clean(row.get(f"name_{role}")),
        annotation_source=clean(row.get(f"{role}_annotation_source")),
        gene_biotype=clean(row.get(f"{role}_gene_biotype")),
        endpoint_type=clean(row.get(f"{role}_endpoint_type")),
        variant_pos=variant_pos,
        variant_ref=variant_ref,
        variant_alt=variant_alt,
        manual_new_gene_hit=str(annotation["manual_new_gene_hit"]),
        manual_new_gene_hit_n=int(annotation["manual_new_gene_hit_n"]),
        manual_new_gene_ids=str(annotation["manual_new_gene_ids"]),
        manual_new_gene_strands=str(annotation["manual_new_gene_strands"]),
        manual_new_gene_class_codes=str(
            annotation["manual_new_gene_class_codes"]
        ),
        manual_new_gene_novelty_categories=str(
            annotation["manual_new_gene_novelty_categories"]
        ),
        manual_new_gene_candidate_categories=str(
            annotation["manual_new_gene_candidate_categories"]
        ),
        manual_new_gene_dataset_support=str(
            annotation["manual_new_gene_dataset_support"]
        ),
        manual_new_gene_merge_statuses=str(
            annotation["manual_new_gene_merge_statuses"]
        ),
    )


def block_output(block: Mapping[str, object]) -> Dict[str, object]:
    return {
        "block_id": block["block_id"],
        "feature_origin": block["feature_origin"],
        "expression_matrix_scope": block["expression_matrix_scope"],
        "block_interpretation_note": block["interpretation_note"],
        "eqtl_tissue": block["tissue"],
        "feature_id": block["feature_id"],
        "feature_primary_id": block["feature_primary_id"],
        "feature_label": block["feature_label"],
        "feature_class": block["feature_class"],
        "expression_gene_id": block["expression_gene_id"],
        "expression_source": block["expression_source"],
        "eqtl_chrom": block["chrom"],
        "selected_sample_count": block["selected_sample_count"],
        "selected_samples": block["selected_samples"],
        "common_marker_count": block["common_marker_count"],
        "total_reads": block["total_reads"],
        "credible_set_size": block["credible_set_size"],
        "credible_set_ppa_sum": block["credible_set_ppa_sum"],
        "max_ppa": block["max_ppa"],
        "top_marker_id": block["top_marker_id"],
        "top_marker_pos": block["top_marker_pos"],
        "top_ppa_tie_count": block["top_ppa_tie_count"],
        "top_marker_is_unique": block["top_marker_is_unique"],
        "block_ppa_status": block["ppa_status"],
    }


def variant_output(row: Mapping[str, str], ppa: float) -> Dict[str, object]:
    ref = clean(row.get("ref")).upper()
    alt = clean(row.get("alt")).upper()
    alleles = allele_key(ref, alt)
    return {
        "eqtl_pos": clean(row.get("pos")),
        "marker_id": clean(row.get("marker_id")),
        "ref": ref,
        "alt": alt,
        "unordered_ref_alt": "|".join(alleles),
        "ppa": f"{ppa:.15g}",
        "ppa_tier": clean(row.get("ppa_tier")),
        "variant_ppa_status": clean(row.get("ppa_status")),
        "marker_rank_by_ppa": clean(row.get("marker_rank_by_ppa")),
        "is_top_variant": clean(row.get("is_top_variant")),
        "catalogue_role": clean(row.get("catalogue_role")),
        "phase_scope": clean(row.get("phase_scope")),
        "local_ps_by_sample": clean(row.get("local_ps_by_sample")),
        "variant_interpretation_note": clean(row.get("interpretation_note")),
    }


def link_output(link: Link) -> Dict[str, object]:
    return {
        "animal_id": link.animal_id,
        "link_tissue": link.tissue,
        "sample_tissue": link.sample_tissue,
        "link_chromosome": link.chromosome,
        "canonical_pair_id": link.canonical_pair_id,
        "sample_pair_id": link.sample_pair_id,
        "evidence_stratum": link.evidence_stratum,
        "overlap_bp": link.overlap_bp,
        "either_endpoint_multigene": link.either_endpoint_multigene,
        "shared_variant_identity_checked": link.shared_variant_identity_checked,
        "source_link_file": link.source_link_file,
        "source_link_line": link.source_link_line,
        "pair_type": link.pair_type,
        "interval_gap_bp": link.interval_gap_bp,
        "mechanism": link.mechanism,
        "linkage_score": link.linkage_score,
    }


def endpoint_output(prefix: str, endpoint: Endpoint) -> Dict[str, object]:
    return {
        f"{prefix}_record_type": endpoint.record_type,
        f"{prefix}_gene_id": endpoint.gene_id,
        f"{prefix}_gene_symbol": endpoint.gene_symbol,
        f"{prefix}_gene_strand": endpoint.gene_strand,
        f"{prefix}_chr": endpoint.chrom,
        f"{prefix}_start": endpoint.start,
        f"{prefix}_end": endpoint.end,
        f"{prefix}_name": endpoint.name,
        f"{prefix}_annotation_source": endpoint.annotation_source,
        f"{prefix}_gene_biotype": endpoint.gene_biotype,
        f"{prefix}_endpoint_type": endpoint.endpoint_type,
        f"{prefix}_manual_new_gene_hit": endpoint.manual_new_gene_hit,
        f"{prefix}_manual_new_gene_hit_n": endpoint.manual_new_gene_hit_n,
        f"{prefix}_manual_new_gene_ids": endpoint.manual_new_gene_ids,
        f"{prefix}_manual_new_gene_strands": endpoint.manual_new_gene_strands,
        f"{prefix}_manual_new_gene_class_codes": (
            endpoint.manual_new_gene_class_codes
        ),
        f"{prefix}_manual_new_gene_novelty_categories": (
            endpoint.manual_new_gene_novelty_categories
        ),
        f"{prefix}_manual_new_gene_candidate_categories": (
            endpoint.manual_new_gene_candidate_categories
        ),
        f"{prefix}_manual_new_gene_dataset_support": (
            endpoint.manual_new_gene_dataset_support
        ),
        f"{prefix}_manual_new_gene_merge_statuses": (
            endpoint.manual_new_gene_merge_statuses
        ),
    }


def endpoint_hit_output(
    link: Link,
    roles: Sequence[str],
    mapping_rules: Mapping[str, str],
    evidence_class: str,
    evidence_level: str,
) -> Dict[str, object]:
    role_set = set(roles)
    endpoint_hit = "both" if role_set == {"base", "target"} else roles[0]
    hit_endpoints = [
        link.base if role == "base" else link.target for role in roles
    ]
    hit_has_new = any(
        endpoint.manual_new_gene_hit == "yes" for endpoint in hit_endpoints
    )
    link_has_new = any(
        endpoint.manual_new_gene_hit == "yes"
        for endpoint in (link.base, link.target)
    )
    if hit_has_new:
        new_gene_scope = "direct_hit_endpoint_manual_new_gene"
    elif link_has_new:
        new_gene_scope = "other_endpoint_manual_new_gene_context"
    else:
        new_gene_scope = "no_manual_new_gene_context"
    return {
        "endpoint_hit": endpoint_hit,
        "endpoint_mapping_rules": "|".join(
            f"{role}:{mapping_rules[role]}" for role in sorted(role_set)
        ),
        "evidence_class": evidence_class,
        "evidence_level": evidence_level,
        "hit_endpoint_gene_id": join_values(
            endpoint.gene_id for endpoint in hit_endpoints
        ),
        "hit_endpoint_gene_symbol": join_values(
            endpoint.gene_symbol for endpoint in hit_endpoints
        ),
        "hit_endpoint_has_manual_new_gene": "yes" if hit_has_new else "no",
        "link_contains_manual_new_gene": "yes" if link_has_new else "no",
        "new_gene_evidence_scope": new_gene_scope,
    }


def complete_overlap_output(
    block: Mapping[str, object],
    link: Link,
    evidence: Mapping[str, object],
    variant: Mapping[str, object] | None = None,
) -> Dict[str, object]:
    output: Dict[str, object] = {}
    output.update(block_output(block))
    if variant is not None:
        output.update(variant)
    output.update(link_output(link))
    output.update(endpoint_output("base", link.base))
    output.update(endpoint_output("target", link.target))
    output.update(evidence)
    return output


def open_gzip_writer(path: Path, fieldnames: Sequence[str]):
    handle = gzip.open(path, "wt", encoding="utf-8", newline="")
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fieldnames,
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    return handle, writer


def write_plain_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def matching_gene_ids(endpoint: Endpoint) -> Mapping[str, str]:
    if endpoint.record_type in {"known_gene", "manual_new_gene"}:
        return {
            endpoint.gene_id: "expression_gene_id_equals_allelome_gene_endpoint"
        }
    return {
        gene_id: (
            "expression_gene_id_equals_manual_new_gene_annotation_"
            "of_unannotated_endpoint"
        )
        for gene_id in split_values(endpoint.manual_new_gene_ids)
    }


def load_links_and_write_gene_overlaps(
    links_path: Path,
    blocks: Mapping[str, Dict[str, object]],
    gene_block_index: Mapping[Tuple[str, str, str, str], Sequence[str]],
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    gene_writer: csv.DictWriter,
    block_metrics: MutableMapping[str, Counter],
) -> Tuple[
    Mapping[Tuple[str, str, int, Tuple[str, str]], Sequence[EndpointRef]],
    Mapping[Tuple[str, str, int], Sequence[EndpointRef]],
    Counter,
]:
    exact_index: MutableMapping[
        Tuple[str, str, int, Tuple[str, str]], List[EndpointRef]
    ] = defaultdict(list)
    interval_index: MutableMapping[
        Tuple[str, str, int], List[EndpointRef]
    ] = defaultdict(list)
    metrics = Counter()
    seen_links: Set[str] = set()
    hit_links: Set[str] = set()
    hit_blocks: Set[str] = set()

    with links_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(links_path, reader.fieldnames or [], LINK_REQUIRED_COLUMNS)
        for line_number, row in enumerate(reader, start=2):
            metrics["allelome_link_rows"] += 1
            chrom_key = normalize_chrom(row["chr_base"])
            endpoint_ids = sorted([row["base_gene_id"], row["target_gene_id"]])
            canonical_pair_id = json.dumps([chrom_key, *endpoint_ids], separators=(",", ":"))
            sample_pair_id = json.dumps([row["animal_id"], row["tissue"], chrom_key, *endpoint_ids], separators=(",", ":"))
            if not canonical_pair_id or sample_pair_id in seen_links:
                raise ValueError(
                    f"{links_path}:{line_number} canonical_pair_id missing or duplicate"
                )
            seen_links.add(sample_pair_id)
            animal_id = clean(row.get("animal_id"))
            tissue = clean(row.get("tissue"))
            chromosome = normalize_chrom(row["chr_base"])
            if (
                animal_id not in KNOWN_DOGS
                or tissue not in KNOWN_TISSUES
                or chromosome not in AUTOSOMES
            ):
                raise ValueError(
                    f"{links_path}:{line_number} dog, tissue or chromosome invalid"
                )
            mechanism = clean(row.get("mechanism"))
            if mechanism not in {"enhancing", "repressing"}:
                raise ValueError(
                    f"{links_path}:{line_number} mechanism invalid: {mechanism}"
                )
            parse_float(row.get("linkage_score"), "linkage_score")
            base = endpoint_from_row("base", row, new_gene_bins)
            target = endpoint_from_row("target", row, new_gene_bins)
            if base.chrom != chromosome or target.chrom != chromosome:
                raise ValueError(
                    f"{links_path}:{line_number} link with endpoint chromosome mismatch"
                )
            link = Link(
                animal_id=animal_id,
                tissue=tissue,
                sample_tissue=f"{animal_id}|{tissue}",
                chromosome=chromosome,
                canonical_pair_id=canonical_pair_id,
                sample_pair_id=sample_pair_id,
                evidence_stratum=row["evidence_stratum"],
                overlap_bp=row["overlap_bp"],
                either_endpoint_multigene=row["either_endpoint_multigene"],
                shared_variant_identity_checked=row["shared_variant_identity_checked"],
                source_link_file=row["source_link_file"],
                source_link_line=row["source_link_line"],
                pair_type=clean(row.get("pair_type")),
                interval_gap_bp=clean(row.get("interval_gap_bp")),
                mechanism=mechanism,
                linkage_score=clean(row.get("linkage_score")),
                base=base,
                target=target,
            )

            grouped_blocks: MutableMapping[str, Dict[str, str]] = defaultdict(dict)
            for role, endpoint in (("base", base), ("target", target)):
                endpoint_ref = EndpointRef(link=link, role=role)
                if endpoint.record_type == "unannotated_variant":
                    assert endpoint.variant_pos is not None
                    key = (
                        tissue,
                        chromosome,
                        endpoint.variant_pos,
                        allele_key(endpoint.variant_ref, endpoint.variant_alt),
                    )
                    exact_index[key].append(endpoint_ref)
                else:
                    for bin_id in range(
                        endpoint.start // INTERVAL_BIN_SIZE,
                        (endpoint.end - 1) // INTERVAL_BIN_SIZE + 1,
                    ):
                        interval_index[(tissue, chromosome, bin_id)].append(
                            endpoint_ref
                        )

                for gene_id, rule in matching_gene_ids(endpoint).items():
                    block_key = (
                        animal_id,
                        tissue,
                        chromosome,
                        gene_id,
                    )
                    for block_id in gene_block_index.get(block_key, []):
                        grouped_blocks[block_id][role] = rule

            for block_id, role_rules in grouped_blocks.items():
                roles = sorted(role_rules)
                link_evidence = endpoint_hit_output(
                    link,
                    roles,
                    role_rules,
                    "expression_gene_matches_allelome_endpoint",
                    "direct_feature_endpoint_evidence",
                )
                gene_writer.writerow(
                    complete_overlap_output(
                        blocks[block_id],
                        link,
                        link_evidence,
                    )
                )
                metrics["expression_gene_link_overlap_rows"] += 1
                block_metrics[block_id]["expression_gene_link_overlap_rows"] += 1
                hit_links.add(sample_pair_id)
                hit_blocks.add(block_id)

            if metrics["allelome_link_rows"] % 50_000 == 0:
                log(
                    "[stage 4/7] Allelome links="
                    f"{metrics['allelome_link_rows']}; "
                    "expression-gene overlaps="
                    f"{metrics['expression_gene_link_overlap_rows']}"
                )

    metrics["expression_gene_hit_links"] = len(hit_links)
    metrics["expression_gene_hit_blocks"] = len(hit_blocks)
    return exact_index, interval_index, metrics


def validate_credible_row(
    path: Path,
    line_number: int,
    row: Mapping[str, str],
    block: Mapping[str, object],
) -> Tuple[int, float, Tuple[str, str]]:
    if clean(row.get("is_credible_set_member")) != "1":
        raise ValueError(f"{path}:{line_number} is not 95% credible-set member")
    if clean(row.get("credible_set_member")) != "1":
        raise ValueError(f"{path}:{line_number} credible_set_member is not 1")
    if (
        clean(row.get("result_status")) != "valid"
        or clean(row.get("result_interpretable")) != "1"
    ):
        raise ValueError(f"{path}:{line_number} is not interpretable valid result")
    comparisons = {
        "feature_origin": clean(row.get("feature_origin")),
        "expression_matrix_scope": clean(row.get("expression_matrix_scope")),
        "tissue": clean(row.get("tissue")),
        "feature_id": clean(row.get("feature_id")),
        "feature_primary_id": clean(row.get("feature_primary_id")),
        "feature_label": clean(row.get("feature_label")),
        "feature_class": clean(row.get("feature_class")),
        "expression_gene_id": clean(row.get("expression_gene_id")),
        "expression_source": clean(row.get("expression_source")),
        "chrom": normalize_chrom(row.get("chrom")),
        "selected_sample_count": parse_int(
            row.get("selected_sample_count"),
            f"{path}:{line_number} selected_sample_count",
        ),
        "common_marker_count": parse_int(
            row.get("common_marker_count"),
            f"{path}:{line_number} common_marker_count",
        ),
        "total_reads": parse_int(
            row.get("total_reads"),
            f"{path}:{line_number} total_reads",
        ),
    }
    for key, observed in comparisons.items():
        if observed != block[key]:
            raise ValueError(
                f"{path}:{line_number} with valid block of {key} mismatch: "
                f"{observed!r} != {block[key]!r}"
            )
    if not clean(row.get("interpretation_note")):
        raise ValueError(f"{path}:{line_number} variant interpretation_note missing")
    chrom = str(block["chrom"])
    if normalize_chrom(row.get("phase_scope")) != chrom:
        raise ValueError(
            f"{path}:{line_number} phase_scope is not block all chromosome span"
        )
    pos = parse_int(row.get("pos"), f"{path}:{line_number} pos")
    if pos < 1:
        raise ValueError(f"{path}:{line_number} pos invalid")
    ppa = parse_float(row.get("ppa"), f"{path}:{line_number} ppa")
    if not 0.0 <= ppa <= 1.0:
        raise ValueError(f"{path}:{line_number} PPA not in 0-1")
    if clean(row.get("ppa_status")) != block["ppa_status"]:
        raise ValueError(
            f"{path}:{line_number} variant with block of ppa_status mismatch"
        )
    if clean(row.get("ppa_tier")) != expected_ppa_tier(ppa):
        raise ValueError(
            f"{path}:{line_number} ppa_tier with PPA numeric value mismatch"
        )
    marker_rank = parse_int(
        row.get("marker_rank_by_ppa"),
        f"{path}:{line_number} marker_rank_by_ppa",
    )
    is_top = clean(row.get("is_top_variant"))
    if marker_rank < 1 or is_top not in {"0", "1"}:
        raise ValueError(
            f"{path}:{line_number} marker rank or is_top_variant invalid"
        )
    if (marker_rank == 1) != (is_top == "1"):
        raise ValueError(
            f"{path}:{line_number} marker rank with is_top_variant mismatch"
        )
    expected_role = "top_and_credible_set" if is_top == "1" else "credible_set"
    if clean(row.get("catalogue_role")) != expected_role:
        raise ValueError(
            f"{path}:{line_number} catalogue_role with top/credible flag mismatch"
        )
    alleles = allele_key(row.get("ref"), row.get("alt"))
    expected_marker = (
        f"{chrom}:{pos}_{clean(row.get('ref')).upper()}>"
        f"{clean(row.get('alt')).upper()}"
    )
    if clean(row.get("marker_id")) != expected_marker:
        raise ValueError(
            f"{path}:{line_number} marker_id with CHROM/POS/REF/ALT mismatch"
        )
    return pos, ppa, alleles


def process_credible_variants(
    path: Path,
    blocks: Mapping[str, Dict[str, object]],
    exact_index: Mapping[
        Tuple[str, str, int, Tuple[str, str]], Sequence[EndpointRef]
    ],
    interval_index: Mapping[Tuple[str, str, int], Sequence[EndpointRef]],
    expected_credible_rows: int,
    priority_ppa: float,
    exact_writer: csv.DictWriter,
    context_writer: csv.DictWriter,
    block_metrics: MutableMapping[str, Counter],
) -> Counter:
    metrics = Counter()
    exact_hit_links: Set[str] = set()
    exact_hit_markers: Set[str] = set()
    context_hit_links: Set[str] = set()
    context_hit_markers: Set[str] = set()

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], CREDIBLE_REQUIRED_COLUMNS)
        for line_number, row in enumerate(reader, start=2):
            metrics["credible_variant_rows"] += 1
            block_id = clean(row.get("block_id"))
            if block_id not in blocks:
                raise ValueError(
                    f"{path}:{line_number} block_id not in valid blocks: {block_id}"
                )
            block = blocks[block_id]
            pos, ppa, alleles = validate_credible_row(
                path, line_number, row, block
            )
            block_metrics[block_id]["credible_variant_rows_observed"] += 1
            selected_samples = set(block["selected_samples_tuple"])
            tissue = str(block["tissue"])
            chrom = str(block["chrom"])
            variant_fields = variant_output(row, ppa)

            exact_key = (tissue, chrom, pos, alleles)
            for endpoint_ref in exact_index.get(exact_key, []):
                if endpoint_ref.link.animal_id not in selected_samples:
                    continue
                endpoint = endpoint_ref.endpoint
                evidence = endpoint_hit_output(
                    endpoint_ref.link,
                    [endpoint_ref.role],
                    {
                        endpoint_ref.role: (
                            "credible_variant_chrom_pos_unordered_ref_alt_"
                            "equals_unannotated_ase_endpoint"
                        )
                    },
                    "credible_variant_exact_unannotated_endpoint",
                    "direct_variant_endpoint_evidence",
                )
                exact_writer.writerow(
                    complete_overlap_output(
                        block,
                        endpoint_ref.link,
                        evidence,
                        variant_fields,
                    )
                )
                metrics["exact_variant_endpoint_overlap_rows"] += 1
                block_metrics[block_id]["exact_variant_endpoint_overlap_rows"] += 1
                exact_hit_links.add(endpoint_ref.link.sample_pair_id)
                exact_hit_markers.add(clean(row.get("marker_id")))
                if endpoint.manual_new_gene_hit == "yes":
                    metrics["exact_endpoint_manual_new_gene_rows"] += 1

            if ppa >= priority_ppa:
                metrics["priority_variant_rows"] += 1
                block_metrics[block_id][
                    "priority_variant_rows_ppa_ge_threshold"
                ] += 1
                point0 = pos - 1
                interval_key = (
                    tissue,
                    chrom,
                    point0 // INTERVAL_BIN_SIZE,
                )
                for endpoint_ref in interval_index.get(interval_key, []):
                    if endpoint_ref.link.animal_id not in selected_samples:
                        continue
                    endpoint = endpoint_ref.endpoint
                    if not endpoint.start <= point0 < endpoint.end:
                        continue
                    same_gene = endpoint.gene_id == block["expression_gene_id"]
                    evidence_class = (
                        "priority_variant_inside_expression_gene_endpoint"
                        if same_gene
                        else "priority_variant_inside_other_gene_endpoint"
                    )
                    evidence = endpoint_hit_output(
                        endpoint_ref.link,
                        [endpoint_ref.role],
                        {
                            endpoint_ref.role: (
                                "ppa_ge_threshold_variant_pos0_inside_"
                                "gene_span_start0_end0"
                            )
                        },
                        evidence_class,
                        "positional_context_not_direct_endpoint_evidence",
                    )
                    context_writer.writerow(
                        complete_overlap_output(
                            block,
                            endpoint_ref.link,
                            evidence,
                            variant_fields,
                        )
                    )
                    metrics["priority_gene_span_context_rows"] += 1
                    block_metrics[block_id][
                        "priority_gene_span_context_rows"
                    ] += 1
                    context_hit_links.add(endpoint_ref.link.sample_pair_id)
                    context_hit_markers.add(clean(row.get("marker_id")))

            if metrics["credible_variant_rows"] % 1_000_000 == 0:
                log(
                    "[stage 5/7] credible variants="
                    f"{metrics['credible_variant_rows']}/"
                    f"{expected_credible_rows}; "
                    f"exact rows={metrics['exact_variant_endpoint_overlap_rows']}; "
                    "priority context rows="
                    f"{metrics['priority_gene_span_context_rows']}"
                )

    if metrics["credible_variant_rows"] != expected_credible_rows:
        raise ValueError(
            "credible set row count and final QC mismatch: "
            f"observed={metrics['credible_variant_rows']}, "
            f"expected={expected_credible_rows}"
        )
    missing_or_mismatched = [
        block_id
        for block_id, block in blocks.items()
        if block_metrics[block_id]["credible_variant_rows_observed"]
        != block["credible_set_size"]
    ]
    if missing_or_mismatched:
        raise ValueError(
            "by block statistics of credible set row count and valid block credible_set_size mismatch: "
            f"{missing_or_mismatched[:5]}"
        )
    metrics["exact_hit_links"] = len(exact_hit_links)
    metrics["exact_hit_markers"] = len(exact_hit_markers)
    metrics["priority_context_hit_links"] = len(context_hit_links)
    metrics["priority_context_hit_markers"] = len(context_hit_markers)
    return metrics


def write_block_summary(
    path: Path,
    blocks: Mapping[str, Dict[str, object]],
    block_metrics: Mapping[str, Counter],
    priority_ppa: float,
) -> None:
    handle, writer = open_gzip_writer(path, BLOCK_SUMMARY_COLUMNS)
    try:
        for block_id in sorted(blocks):
            output = block_output(blocks[block_id])
            counts = block_metrics[block_id]
            output.update(
                {
                    "credible_variant_rows_observed": counts[
                        "credible_variant_rows_observed"
                    ],
                    "priority_variant_rows_ppa_ge_threshold": counts[
                        "priority_variant_rows_ppa_ge_threshold"
                    ],
                    "expression_gene_link_overlap_rows": counts[
                        "expression_gene_link_overlap_rows"
                    ],
                    "exact_variant_endpoint_overlap_rows": counts[
                        "exact_variant_endpoint_overlap_rows"
                    ],
                    "priority_gene_span_context_rows": counts[
                        "priority_gene_span_context_rows"
                    ],
                    "priority_ppa_threshold": priority_ppa,
                }
            )
            writer.writerow(output)
    finally:
        handle.close()


def write_manifest(
    path: Path,
    inputs: Sequence[Tuple[str, Path, str]],
    outputs: Sequence[Tuple[str, Path, str]],
) -> None:
    rows = []
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
    write_plain_tsv(
        path,
        ["kind", "label", "path", "bytes", "purpose"],
        rows,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Integrate five-dog Allelome links across all chromosomes with joint expression+ASE eQTL fine-mapping results."
        )
    )
    parser.add_argument("--allelome-links", type=Path, default=DEFAULT_ALLELOME_LINKS)
    parser.add_argument(
        "--allelome-summary", type=Path, default=DEFAULT_ALLELOME_SUMMARY
    )
    parser.add_argument(
        "--candidate-blocks", type=Path, default=DEFAULT_CANDIDATE_BLOCKS
    )
    parser.add_argument("--valid-blocks", type=Path, default=DEFAULT_VALID_BLOCKS)
    parser.add_argument(
        "--credible-variants", type=Path, default=DEFAULT_CREDIBLE_VARIANTS
    )
    parser.add_argument(
        "--final-qc-summary", type=Path, default=DEFAULT_FINAL_QC_SUMMARY
    )
    parser.add_argument("--new-gene-gtf", type=Path, default=DEFAULT_NEW_GENE_GTF)
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
    parser.add_argument(
        "--priority-ppa",
        type=float,
        default=DEFAULT_PRIORITY_PPA,
        help="Minimum PPA for gene-span positional-context analysis; default 0.1.",
    )
    return parser.parse_args()


def main() -> int:
    set_large_csv_fields()
    args = parse_args()
    if not 0.0 <= args.priority_ppa <= 1.0:
        raise ValueError("--priority-ppa must be located in 0-1")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.run_id):
        raise ValueError("--run-id can only contains letter, number, point, underscore and hyphen")

    args.output_root.mkdir(parents=True, exist_ok=True)
    run_dir = args.output_root / args.run_id
    if run_dir.exists():
        raise FileExistsError(f"run directory already exists, refuse to overwrite: {run_dir}")
    run_dir.mkdir()

    inputs = [
        ("allelome_links", args.allelome_links, "Canonical Allelome links for five dogs"),
        ("allelome_summary", args.allelome_summary, "Allelome 1.4 PASS summary"),
        ("candidate_blocks", args.candidate_blocks, "Dogs actually included in each 8.2 block and chromosome-wide phase_scope"),
        ("valid_blocks", args.valid_blocks, "8.5 interpretable PLASMA block summary"),
        ("credible_variants", args.credible_variants, "8.6 final 95% credible-set variants"),
        ("final_qc_summary", args.final_qc_summary, "8.7 final PASS and count-reconciliation summary"),
        ("new_gene_gtf", args.new_gene_gtf, "Five-dog DOG5_NEWGENE structure GTF"),
        (
            "new_gene_transcript_map",
            args.new_gene_transcript_map,
            "Five-dog DOG5_NEWGENE transcript and candidate provenance map",
        ),
        ("new_gene_build_summary", args.new_gene_build_summary, "Five-dog GTF construction summary"),
    ]
    gene_path = (
        run_dir
        / "dog5_allelome_eqtl_finemap_expression_gene_link_overlap.tsv.gz"
    )
    exact_path = (
        run_dir
        / "dog5_allelome_eqtl_finemap_credible_variant_exact_endpoint_overlap.tsv.gz"
    )
    context_path = (
        run_dir
        / "dog5_allelome_eqtl_finemap_priority_variant_gene_span_context.tsv.gz"
    )
    block_summary_path = (
        run_dir / "dog5_allelome_eqtl_finemap_block_summary.tsv.gz"
    )
    run_summary_path = run_dir / "dog5_allelome_eqtl_finemap_run_summary.tsv"
    manifest_path = run_dir / "dog5_allelome_eqtl_finemap_manifest.tsv"
    completion_marker_path = run_dir / "execution_completion.marker"

    try:
        log("[stage 1/7] start and check all input path")
        for label, path, _ in inputs:
            require_file(path, label)
            log(f"[stage 1/7] {label}: {path.resolve()}")
        validate_allelome_summary(args.allelome_summary)
        final_qc = validate_final_qc(
            args.final_qc_summary,
            args.credible_variants,
            args.valid_blocks,
        )

        expected_valid_blocks = parse_int(
            final_qc["source_8.6_valid_blocks"],
            "final QC source_8.6_valid_blocks",
        )
        expected_credible_rows = parse_int(
            final_qc["source_8.6_credible_set_variant_rows"],
            "final QC source_8.6_credible_set_variant_rows",
        )

        log("[stage 2/7] read and closure five-dog-specific DOG5_NEWGENE annotation")
        new_gene_metrics, new_gene_bins = load_new_gene_loci(
            args.new_gene_gtf,
            args.new_gene_transcript_map,
            args.new_gene_build_summary,
        )
        log(
            "[stage 2/7] DOG5_NEWGENE loci="
            f"{new_gene_metrics['total_loci']}; autosomal="
            f"{new_gene_metrics['autosomal_loci']}; non-autosomal="
            f"{new_gene_metrics['non_autosomal_loci']}"
        )

        log("[stage 3/7] read 8.2 candidate blocks and 8.5 valid blocks")
        candidate_blocks = load_candidate_blocks(args.candidate_blocks)
        blocks = load_valid_blocks(
            args.valid_blocks,
            candidate_blocks,
            expected_valid_blocks,
        )
        gene_block_index = build_gene_block_index(blocks)
        feature_origin_counts = Counter(
            str(block["feature_origin"]) for block in blocks.values()
        )
        log(
            f"[stage 3/7] candidate blocks={len(candidate_blocks)}; "
            f"valid blocks={len(blocks)}; "
            f"selected-sample gene keys={len(gene_block_index)}"
        )

        block_metrics: MutableMapping[str, Counter] = defaultdict(Counter)
        gene_handle, gene_writer = open_gzip_writer(
            gene_path, GENE_OVERLAP_COLUMNS
        )
        exact_handle, exact_writer = open_gzip_writer(
            exact_path, VARIANT_OVERLAP_COLUMNS
        )
        context_handle, context_writer = open_gzip_writer(
            context_path, VARIANT_OVERLAP_COLUMNS
        )
        try:
            log("[stage 4/7] scan Allelome links, output block level expression gene relationship and construct endpoint index")
            exact_index, interval_index, link_metrics = (
                load_links_and_write_gene_overlaps(
                    args.allelome_links,
                    blocks,
                    gene_block_index,
                    new_gene_bins,
                    gene_writer,
                    block_metrics,
                )
            )
            gene_handle.close()
            log(
                "[stage 4/7] Allelome completed: "
                f"exact keys={len(exact_index)}; interval bins={len(interval_index)}"
            )

            log(
                "[stage 5/7] scanning all credible sets (95%) in streaming mode; "
                "local PS is retained only for provenance and is not used for filtering"
            )
            variant_metrics = process_credible_variants(
                args.credible_variants,
                blocks,
                exact_index,
                interval_index,
                expected_credible_rows,
                args.priority_ppa,
                exact_writer,
                context_writer,
                block_metrics,
            )
        finally:
            if not gene_handle.closed:
                gene_handle.close()
            exact_handle.close()
            context_handle.close()

        log("[stage 6/7] write per block closure table, run summary and list")
        write_block_summary(
            block_summary_path,
            blocks,
            block_metrics,
            args.priority_ppa,
        )
        summary_rows = [
            {"metric": "run_id", "value": args.run_id},
            {"metric": "overall_status", "value": "PASS"},
            {"metric": "allelome_link_rows", "value": link_metrics["allelome_link_rows"]},
            {
                "metric": "valid_blocks",
                "value": len(blocks),
            },
            {
                "metric": "credible_variant_rows",
                "value": variant_metrics["credible_variant_rows"],
            },
            {
                "metric": "expression_matrix_scope",
                "value": EXPECTED_EXPRESSION_MATRIX_SCOPE,
            },
            {
                "metric": "reference_annotation_valid_blocks",
                "value": feature_origin_counts["reference_annotation_gene"],
            },
            {
                "metric": "dog5_manual_new_gene_valid_blocks",
                "value": feature_origin_counts["dog5_manual_new_gene"],
            },
            {
                "metric": "dog5_manual_new_gene_loci_total",
                "value": new_gene_metrics["total_loci"],
            },
            {
                "metric": "dog5_manual_new_gene_loci_autosomal",
                "value": new_gene_metrics["autosomal_loci"],
            },
            {
                "metric": "dog5_manual_new_gene_loci_non_autosomal_excluded",
                "value": new_gene_metrics["non_autosomal_loci"],
            },
            {
                "metric": "priority_ppa_threshold",
                "value": args.priority_ppa,
            },
            {
                "metric": "priority_variant_rows",
                "value": variant_metrics["priority_variant_rows"],
            },
            {
                "metric": "expression_gene_link_overlap_rows",
                "value": link_metrics["expression_gene_link_overlap_rows"],
            },
            {
                "metric": "expression_gene_hit_links",
                "value": link_metrics["expression_gene_hit_links"],
            },
            {
                "metric": "expression_gene_hit_blocks",
                "value": link_metrics["expression_gene_hit_blocks"],
            },
            {
                "metric": "exact_variant_endpoint_overlap_rows",
                "value": variant_metrics["exact_variant_endpoint_overlap_rows"],
            },
            {
                "metric": "exact_hit_links",
                "value": variant_metrics["exact_hit_links"],
            },
            {
                "metric": "exact_hit_markers",
                "value": variant_metrics["exact_hit_markers"],
            },
            {
                "metric": "exact_endpoint_manual_new_gene_rows",
                "value": variant_metrics["exact_endpoint_manual_new_gene_rows"],
            },
            {
                "metric": "priority_gene_span_context_rows",
                "value": variant_metrics["priority_gene_span_context_rows"],
            },
            {
                "metric": "priority_context_hit_links",
                "value": variant_metrics["priority_context_hit_links"],
            },
            {
                "metric": "priority_context_hit_markers",
                "value": variant_metrics["priority_context_hit_markers"],
            },
            {
                "metric": "phase_scope_rule",
                "value": "phase_scope_must_equal_chrN",
            },
            {
                "metric": "local_ps_rule",
                "value": "provenance_only_not_used_for_filtering",
            },
            {
                "metric": "association_fdr_screening",
                "value": final_qc["source_8.6_association_fdr_screening"],
            },
            {
                "metric": "null_causal_model",
                "value": final_qc["source_8.6_null_causal_model"],
            },
            {
                "metric": "interpretation_scope",
                "value": final_qc["source_8.6_interpretation_scope"],
            },
            {
                "metric": "not_interpretable_as",
                "value": final_qc["source_8.6_not_interpretable_as"],
            },
        ]
        write_plain_tsv(run_summary_path, ["metric", "value"], summary_rows)

        outputs = [
            ("expression_gene_link_overlap", gene_path, "Block-level direct relationships between expression genes and Allelome endpoints"),
            ("credible_variant_exact_endpoint_overlap", exact_path, "Exact matches between all credible-set SNPs and unannotated ASE endpoints"),
            ("priority_variant_gene_span_context", context_path, "Gene-span positional context for SNPs above the PPA threshold"),
            ("block_summary", block_summary_path, "Input and match-count reconciliation for each valid block"),
            ("run_summary", run_summary_path, "Run counts and interpretation limits"),
        ]
        write_manifest(manifest_path, inputs, outputs)

        log("[stage 7/7] all output completed")
        log(f"result directory: {run_dir.resolve()}")
        completion_marker_path.write_text(
            "execution_status=success\n"
            f"run_id={args.run_id}\n"
            f"completed_at={datetime.now():%Y-%m-%d %H:%M:%S}\n",
            encoding="utf-8",
        )
        log("task completed / done")
        log("execution_status=success")
        log("completion only indicates that the main program exited normally; output completeness, format, counts and biological quality require separate review")
        return 0
    except Exception as exc:
        log(f"task failed: {exc}")
        for line in traceback.format_exc().rstrip().splitlines():
            log(line)
        log("execution_status=failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
