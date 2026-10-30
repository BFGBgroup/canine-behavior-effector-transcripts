#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import bisect
import json
import csv
import gzip
import math
import re
import sys
import traceback
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Sequence, Set, TextIO, Tuple


DEFAULT_ALLELOME_LINKS = Path("/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_links.annotated.tsv")
DEFAULT_ALLELOME_SUMMARY = Path("/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_annotation.summary.tsv")
DEFAULT_NEW_GENE_GTF = Path("/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf")
DEFAULT_NEW_GENE_METADATA = Path("/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv")
DEFAULT_NEW_GENE_BUILD_SUMMARY = Path("/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.summary.tsv")
DEFAULT_GWAS_RUN_DIR = Path(
    '/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352'
)
DEFAULT_OUTPUT_ROOT = Path("/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap")

NEW_GENE_BIN_SIZE = 1_000_000
AUTOSOMES = {f"chr{number}" for number in range(1, 39)}
BRANCHES = ("main", "strict_no_ambiguous")
ALLOWED_ANALYSIS_SCOPES = {"primary", "exploratory"}
ALLOWED_HARMONIZATION = {
    "HARMONIZED_NONAMBIGUOUS",
    "HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED",
}
VARIANT_RE = re.compile(
    r"(?:\(NA\))?(chr[^:|\s]+):(\d+)_([A-Za-z]+)>([A-Za-z]+)"
)
GTF_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')

LINK_REQUIRED_COLUMNS = {
    "animal_id",
    "tissue",
    "source_link_file", "source_link_line", "source_metadata_file",
    "overlap_bp", "strand_relation", "either_endpoint_multigene",
    "both_endpoints_unambiguous_count_positive", "identical_haplotype_counts",
    "evidence_stratum", "shared_variant_identity_checked",
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
NEW_GENE_METADATA_REQUIRED_COLUMNS = {
    "final_gene_id", "source_dataset", "chr", "strand", "gene_start_0based",
    "gene_end_0based", "class_code", "novelty_category", "gene_locus_filter_reason",
}
GWAS_MANIFEST_REQUIRED_COLUMNS = {
    "locus_id",
    "trait",
    "analysis_scope",
    "branch",
    "status",
    "model_status",
    "gwas_lead_variant_id",
    "gwas_lead_pos",
    "gwas_lead_p_lrt",
    "variant_count",
    "reference_sample_count",
    "n_locus_median",
    "coverage",
    "credible_set_count",
    "max_pip",
    "top_variant_id",
    "top_variant_pos",
    "output_dir",
    "variant_pip",
    "credible_sets",
    "credible_set_members",
}
GWAS_VARIANT_REQUIRED_COLUMNS = {
    "matrix_variant_index", "locus_id", "trait", "analysis_scope",
    "chrom", "pos", "ld_variant_id", "ld_a1", "ld_a2",
    "beta_for_gwas_allele1", "beta_for_ld_a1", "se", "z_for_ld_a1",
    "p_lrt", "n_obs", "ambiguous_at_cg", "harmonization_status",
    "pip", "pip_rank", "credible_sets", "in_credible_set",
}
GWAS_MEMBER_REQUIRED_COLUMNS = {
    "locus_id", "trait", "analysis_scope", "branch", "cs_index", "cs_id",
    "matrix_variant_index", "chrom", "pos", "ld_variant_id", "ld_a1",
    "ld_a2", "pip", "z_for_ld_a1", "p_lrt",
}
GWAS_CS_REQUIRED_COLUMNS = {
    "locus_id", "trait", "analysis_scope", "branch", "cs_index", "cs_id",
    "cs_size", "coverage",
}

GWAS_COLUMNS = [
    "candidate_type",
    "gwas_analysis_scope",
    "gwas_branch",
    "trait",
    "locus_id",
    "gwas_ld_variant_id",
    "gwas_chrom",
    "gwas_pos",
    "gwas_ld_a1",
    "gwas_ld_a2",
    "gwas_unordered_alleles",
    "gwas_beta_for_ld_a1",
    "gwas_se",
    "gwas_z_for_ld_a1",
    "gwas_variant_p_lrt",
    "gwas_n_obs",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_credible_sets",
    "gwas_in_credible_set",
    "gwas_ambiguous_at_cg",
    "gwas_harmonization_status",
    "gwas_locus_status",
    "gwas_model_status",
    "gwas_variant_count",
    "gwas_reference_sample_count",
    "gwas_n_locus_median",
    "gwas_locus_lead_variant_id",
    "gwas_locus_lead_pos",
    "gwas_locus_lead_p_lrt",
    "gwas_locus_top_variant_id",
    "gwas_locus_top_pos",
    "gwas_locus_max_pip",
    "candidate_key",
]
LINK_COLUMNS = [
    "canonical_pair_id",
    "sample_pair_id",
    "source_link_file", "source_link_line", "source_metadata_file",
    "overlap_bp", "strand_relation", "either_endpoint_multigene",
    "both_endpoints_unambiguous_count_positive", "identical_haplotype_counts",
    "evidence_stratum", "shared_variant_identity_checked",
    "animal_id",
    "tissue",
    "sample_tissue",
    "link_chromosome",
    "pair_type",
    "interval_gap_bp",
    "mechanism",
    "linkage_score",
    "endpoint_hit",
    "base_mapping_rule",
    "target_mapping_rule",
]
ENDPOINT_SUFFIXES = [
    "record_type",
    "gene_id",
    "gene_symbol",
    "gene_strand",
    "chr",
    "start",
    "end",
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
    "manual_new_gene_filter_reasons",
]
ENDPOINT_COLUMNS = [
    f"{prefix}_{suffix}"
    for prefix in ("base", "target")
    for suffix in ENDPOINT_SUFFIXES
]
OVERLAP_EXTRA_COLUMNS = [
    "link_contains_manual_new_gene",
    "gwas_hit_endpoint_has_manual_new_gene",
    "new_gene_evidence_scope",
]
OVERLAP_COLUMNS = GWAS_COLUMNS + LINK_COLUMNS + ENDPOINT_COLUMNS + OVERLAP_EXTRA_COLUMNS
TARGET_STRUCTURE_EXTRA_COLUMNS = [
    "evidence_tier",
    "ase_gwas_mapping_class",
    "unannotated_endpoint_role",
    "paired_pcgene_endpoint_role",
    "shared_manual_new_gene_id",
    "shared_new_gene_chrom",
    "shared_new_gene_start0",
    "shared_new_gene_end0",
    "shared_new_gene_strand",
    "shared_new_gene_length_bp",
    "shared_new_gene_class_codes",
    "shared_new_gene_novelty_categories",
    "shared_new_gene_candidate_category",
    "shared_new_gene_dataset_support",
    "shared_new_gene_filter_reason",
    "unannotated_ase_chrom",
    "unannotated_ase_pos",
    "unannotated_ase_ref",
    "unannotated_ase_alt",
    "ase_gwas_distance_bp",
    "gwas_in_paired_pcgene",
    "gwas_is_locus_top_variant",
    "target_structure_pass",
]
TARGET_STRUCTURE_COLUMNS = OVERLAP_COLUMNS + TARGET_STRUCTURE_EXTRA_COLUMNS
HIT_LINK_SUMMARY_COLUMNS = [
    "canonical_pair_id",
    "sample_pair_id",
    "source_link_file", "source_link_line", "source_metadata_file",
    "overlap_bp", "strand_relation", "either_endpoint_multigene",
    "both_endpoints_unambiguous_count_positive", "identical_haplotype_counts",
    "evidence_stratum", "shared_variant_identity_checked",
    "animal_id",
    "tissue",
    "sample_tissue",
    "chromosome",
    "pair_type",
    "mechanism",
    "linkage_score",
    "credible_overlap_row_n",
    "credible_candidate_n",
    "analysis_scopes",
    "traits",
    "endpoint_hits",
    "minimum_gwas_variant_p_lrt",
    "link_contains_manual_new_gene",
    "direct_new_gene_endpoint_hit_n",
    "manual_new_gene_ids",
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


def normalize_chrom(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    if text.startswith("chr"):
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
        return float(clean(value))
    except ValueError as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc


def is_yes(value: object) -> bool:
    return clean(value).lower() in {"yes", "true", "1", "y"}


def allele_tuple(*values: object) -> Tuple[str, ...]:
    return tuple(sorted({clean(value).upper() for value in values if clean(value)}))


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


def open_text(path: Path) -> TextIO:

    require_file(path, "TSV input")
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def iter_tsv(path: Path, required: Set[str]) -> Iterator[Dict[str, str]]:

    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], required)
        for row in reader:
            yield dict(row)


def load_key_values(path: Path) -> Dict[str, str]:

    require_file(path, "key value input")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        raw_rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    if not raw_rows:
        raise ValueError(f"{path} is empty")
    first = raw_rows[0]
    has_header = (
        len(first) >= 2
        and first[1] == "value"
        and first[0] in {"key", "metric", "parameter"}
    )
    rows = raw_rows[1:] if has_header else raw_rows
    if not has_header and (len(first) < 2 or first[0] != "overall_status"):
        raise ValueError(f"{path} does not use the supported two-column key-value format")
    result: Dict[str, str] = {}
    for row in rows:
        if len(row) < 2:
            raise ValueError(f"{path} exists insufficient two columns of key value row: {row!r}")
        key = clean(row[0])
        if not key or key in result:
            raise ValueError(f"{path} contains empty or duplicate key: {key!r}")
        result[key] = clean(row[1])
    return result


def read_two_column_table(path: Path, key_name: str) -> Dict[str, str]:
    require_file(path, "two-column table")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], {key_name, "value"})
        result: Dict[str, str] = {}
        for row in reader:
            key = clean(row.get(key_name))
            if not key or key in result:
                raise ValueError(f"{path} contains empty or duplicate {key_name}: {key!r}")
            result[key] = clean(row.get("value"))
    return result


def candidate_category(class_codes: Sequence[str], novelty: Sequence[str]) -> str:
    code_set = set(class_codes)
    novelty_set = set(novelty)
    categories: List[str] = []
    if "u" in code_set or "intergenic_novel_transcript" in novelty_set:
        categories.append("intergenic_new_gene_candidate")
    if "x" in code_set or any("antisense" in value for value in novelty_set):
        categories.append("antisense_candidate")
    if "i" in code_set or any("intronic" in value for value in novelty_set):
        categories.append("intronic_candidate")
    if "j" in code_set or "known_gene_novel_isoform" in novelty_set:
        categories.append("known_gene_novel_isoform")
    return "|".join(categories or ["other_candidate_transcript_locus"])


def parse_gtf_attributes(text: str) -> Dict[str, str]:
    return {
        match.group(1): match.group(2)
        for match in GTF_ATTRIBUTE_RE.finditer(text)
    }


def validate_new_gene_build_summary(summary_path, gtf_path, metadata_path):
    summary = read_two_column_table(summary_path, "metric")
    if summary.get("source_dataset") != "five-dogs":
        raise ValueError("requires five-dog-specific GTF build result")
    for key, path in (("new_only_gtf", gtf_path), ("transcript_map", metadata_path)):
        if Path(summary[key]).resolve() != path.resolve():
            raise ValueError(f"GTF build input path mismatch: {key}")
    return (int(summary["candidate_gene_count"]), int(summary["candidate_transcript_count"]),
            int(summary["appended_gtf_row_count"]))


def load_new_gene_metadata(metadata_path):

    metadata = {}
    for row in iter_tsv(metadata_path, NEW_GENE_METADATA_REQUIRED_COLUMNS):
        gid = row["final_gene_id"]
        if row["source_dataset"] != "five-dogs" or not gid.startswith("DOG5_NEWGENE_"):
            raise ValueError(f"non-five-dog novel gene: {gid}")
        record = metadata.setdefault(gid, {
            "gene_id": gid, "chrom": normalize_chrom(row["chr"]), "strand": row["strand"],
            "metadata_start0": int(row["gene_start_0based"]), "metadata_end0": int(row["gene_end_0based"]),
            "dataset_support": row["source_dataset"], "class_codes": [], "novelty": [],
            "filter_reason": row["gene_locus_filter_reason"],
        })
        if (record["chrom"], record["strand"], record["metadata_start0"], record["metadata_end0"]) != (
                normalize_chrom(row["chr"]), row["strand"], int(row["gene_start_0based"]), int(row["gene_end_0based"])):
            raise ValueError(f"same gene transcript of gene span mismatch: {gid}")
        record["class_codes"] = sorted(set(record["class_codes"]) | set(split_values(row["class_code"])))
        record["novelty"] = sorted(set(record["novelty"]) | set(split_values(row["novelty_category"])))
        record["candidate_category"] = candidate_category(record["class_codes"], record["novelty"])
    return metadata


def load_new_gene_loci(
    gtf_path: Path,
    metadata_path: Path,
    build_summary_path: Path,
) -> Tuple[List[Dict[str, object]], Dict[Tuple[str, int], List[Dict[str, object]]]]:
    expected_loci, expected_transcripts, expected_gtf_rows = (
        validate_new_gene_build_summary(
            build_summary_path,
            gtf_path,
            metadata_path,
        )
    )
    metadata = load_new_gene_metadata(metadata_path)
    if len(metadata) != expected_loci:
        raise ValueError(
            "new gene metadata row count and Manual GTF summary mismatch: "
            f"observed={len(metadata)}, expected={expected_loci}"
        )

    gene_spans: Dict[str, Dict[str, object]] = {}
    seen_transcripts: Set[str] = set()
    observed_gtf_rows = 0
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
            observed_gtf_rows += 1
            chrom, _, feature, start_text, end_text, _, strand, _, attr_text = fields
            if feature not in {"transcript", "exon"}:
                raise ValueError(
                    f"{gtf_path}:{line_number} contains non-transcript/exon record: "
                    f"{feature!r}"
                )
            attrs = parse_gtf_attributes(attr_text)
            locus_id = clean(attrs.get("gene_id"))
            if not locus_id.startswith("DOG5_NEWGENE_"):
                raise ValueError(f"non-five dogs candidate GTF gene_id: {locus_id}")
            if locus_id not in metadata:
                raise ValueError(
                    f"{gtf_path}:{line_number} GTF gene_id not in five dogs transcript_map in: "
                    f"{locus_id}"
                )
            normalized_chrom = normalize_chrom(chrom)
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"{gtf_path}:{line_number} strand invalid: {strand!r}"
                )
            start1 = parse_int(start_text, f"{gtf_path}:{line_number} start")
            end1 = parse_int(end_text, f"{gtf_path}:{line_number} end")
            if start1 < 1 or end1 < start1:
                raise ValueError(
                    f"{gtf_path}:{line_number} GTF interval invalid: {start1}-{end1}"
                )
            metadata_row = metadata[locus_id]
            if (
                normalized_chrom != metadata_row["chrom"]
                or strand != metadata_row["strand"]
            ):
                raise ValueError(
                    f"{gtf_path}:{line_number} GTF with five dogs transcript_map of chromosome or strand"
                    f"mismatch: {locus_id}"
                )
            if feature != "transcript":
                continue

            observed_transcripts += 1
            transcript_id = clean(attrs.get("transcript_id"))
            if not transcript_id or transcript_id in seen_transcripts:
                raise ValueError(
                    f"{gtf_path}:{line_number} transcript_id missing or duplicate: "
                    f"{transcript_id!r}"
                )
            seen_transcripts.add(transcript_id)
            start0 = start1 - 1
            span = gene_spans.setdefault(
                locus_id,
                {
                    "chrom": normalized_chrom,
                    "strand": strand,
                    "start0": start0,
                    "end0": end1,
                },
            )
            span["start0"] = min(int(span["start0"]), start0)
            span["end0"] = max(int(span["end0"]), end1)

    if observed_gtf_rows != expected_gtf_rows:
        raise ValueError(
            "Manual new-only GTF record count and build summary mismatch: "
            f"observed={observed_gtf_rows}, expected={expected_gtf_rows}"
        )
    if observed_transcripts != expected_transcripts:
        raise ValueError(
            "Manual new-only GTF transcript count and build summary mismatch: "
            f"observed={observed_transcripts}, expected={expected_transcripts}"
        )
    if set(gene_spans) != set(metadata):
        missing = sorted(set(metadata) - set(gene_spans))
        extra = sorted(set(gene_spans) - set(metadata))
        raise ValueError(
            "Manual GTF with five dogs transcript_map of gene ID set mismatch: "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )

    rows: List[Dict[str, object]] = []
    bins: Dict[Tuple[str, int], List[Dict[str, object]]] = defaultdict(list)
    for locus_id in sorted(metadata):
        record = dict(metadata[locus_id])
        span = gene_spans[locus_id]
        start0 = int(span["start0"])
        end0 = int(span["end0"])
        metadata_start0 = int(record.pop("metadata_start0"))
        metadata_end0 = int(record.pop("metadata_end0"))
        if start0 != metadata_start0 or end0 != metadata_end0:
            raise ValueError(
                f"{locus_id} of GTF gene span with 2.8 locus boundary mismatch: "
                f"gtf={start0}-{end0}, metadata={metadata_start0}-{metadata_end0}"
            )
        record["start0"] = start0
        record["end0"] = end0
        rows.append(record)
        chrom = str(record["chrom"])
        if chrom not in AUTOSOMES:
            continue
        for bin_id in range(
            start0 // NEW_GENE_BIN_SIZE,
            (end0 - 1) // NEW_GENE_BIN_SIZE + 1,
        ):
            bins[(chrom, bin_id)].append(record)
    return rows, bins


def find_new_gene_loci(
    bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    chrom: str,
    pos1: int,
) -> List[Dict[str, object]]:
    point0 = pos1 - 1
    return sorted(
        [
            row
            for row in bins.get((chrom, point0 // NEW_GENE_BIN_SIZE), [])
            if int(row["start0"]) <= point0 < int(row["end0"])
        ],
        key=lambda row: clean(row["gene_id"]),
    )


def summarize_new_gene_hits(hits: Sequence[Dict[str, object]]) -> Dict[str, str]:
    return {
        "manual_new_gene_hit": "yes" if hits else "no",
        "manual_new_gene_hit_n": str(len(hits)),
        "manual_new_gene_ids": join_values(row["gene_id"] for row in hits),
        "manual_new_gene_strands": join_strands(row["strand"] for row in hits),
        "manual_new_gene_class_codes": join_values(
            code for row in hits for code in row["class_codes"]
        ),
        "manual_new_gene_novelty_categories": join_values(
            value for row in hits for value in row["novelty"]
        ),
        "manual_new_gene_candidate_categories": join_values(
            value
            for row in hits
            for value in split_values(row["candidate_category"])
        ),
        "manual_new_gene_dataset_support": join_values(
            row["dataset_support"] for row in hits
        ),
        "manual_new_gene_filter_reasons": join_values(
            row["filter_reason"] for row in hits
        ),
    }


def split_cs_labels(value: object) -> Set[str]:

    return {label for label in re.split(r"[,;]", clean(value)) if label}


def finite_float(value: object, label: str) -> float:

    parsed = parse_float(value, label)
    if not math.isfinite(parsed):
        raise ValueError(f"{label} is not finite value: {value!r}")
    return parsed


def validate_gwas_top_level(
    run_dir: Path,
    manifest_rows: Sequence[Mapping[str, str]],
) -> Dict[str, Dict[str, int]]:

    completion_path = run_dir / "COMPLETED.ok"
    overall_path = run_dir / "overall_summary.tsv"
    branch_path = run_dir / "branch_summary.tsv"
    parameters_path = run_dir / "model_parameters.tsv"
    manifest_path = run_dir / "locus_run_manifest.tsv"
    for path in (
        completion_path,
        overall_path,
        branch_path,
        parameters_path,
        manifest_path,
    ):
        require_file(path, "GWAS frozen run input")

    completion = load_key_values(completion_path)
    if completion.get("overall_status") != "PASS":
        raise ValueError("GWAS COMPLETED.ok is not PASS")
    expected_branches = set(BRANCHES)
    completion_branches = {
        value.strip()
        for value in completion.get("branches", "").split(",")
        if value.strip()
    }
    if completion_branches != expected_branches:
        raise ValueError(
            f"GWAS COMPLETED.ok branch unexpected: {completion_branches}"
        )

    overall = load_key_values(overall_path)
    task_count = len(manifest_rows)
    for key, expected in {
        "selected_task_count": task_count,
        "success_task_count": task_count,
        "failed_task_count": 0,
        "skipped_task_count": 0,
    }.items():
        observed = parse_int(overall.get(key), f"overall_summary {key}")
        if observed != expected:
            raise ValueError(
                f"GWAS overall_summary {key} unexpected: "
                f"observed={observed}, expected={expected}"
            )

    parameters = load_key_values(parameters_path)
    if parameters.get("coverage") != "0.95":
        raise ValueError("GWAS SuSiE coverage is not 0.95")
    if parameters.get("z_input") != "z_for_ld_a1":
        raise ValueError("GWAS z_input is not z_for_ld_a1")
    if parameters.get("n_input") != "n_locus_median":
        raise ValueError("GWAS n_input is not n_locus_median")
    parameter_branches = {
        value.strip()
        for value in parameters.get("branches", "").split(",")
        if value.strip()
    }
    if parameter_branches != expected_branches:
        raise ValueError(f"GWAS model_parameters branch unexpected: {parameter_branches}")

    branch_fields, branch_rows = read_tsv(branch_path)
    require_columns(
        branch_path,
        branch_fields,
        {
            "branch", "selected_task_count", "success_task_count",
            "failed_task_count", "skipped_task_count", "modeled_variant_count",
            "credible_set_count",
        },
    )
    result: Dict[str, Dict[str, int]] = {}
    for row in branch_rows:
        branch = clean(row.get("branch"))
        if branch not in expected_branches or branch in result:
            raise ValueError(f"GWAS branch_summary contains empty, duplicate or invalid branch: {branch!r}")
        manifest_count = sum(clean(item.get("branch")) == branch for item in manifest_rows)
        branch_values = {
            "selected_task_count": parse_int(row.get("selected_task_count"), "selected_task_count"),
            "success_task_count": parse_int(row.get("success_task_count"), "success_task_count"),
            "failed_task_count": parse_int(row.get("failed_task_count"), "failed_task_count"),
            "skipped_task_count": parse_int(row.get("skipped_task_count"), "skipped_task_count"),
            "modeled_variant_count": parse_int(row.get("modeled_variant_count"), "modeled_variant_count"),
            "credible_set_count": parse_int(row.get("credible_set_count"), "credible_set_count"),
        }
        if (
            branch_values["selected_task_count"] != manifest_count
            or branch_values["success_task_count"] != manifest_count
            or branch_values["failed_task_count"] != 0
            or branch_values["skipped_task_count"] != 0
        ):
            raise ValueError(f"GWAS branch_summary task statistics with manifest mismatch: {branch}")
        result[branch] = branch_values
    if set(result) != expected_branches:
        raise ValueError("GWAS branch_summary no at the same time contains main and strict_no_ambiguous")
    return result


def validate_locus_artifacts(
    run_dir: Path,
    manifest: Mapping[str, str],
) -> Dict[str, Path]:

    locus_id = clean(manifest.get("locus_id"))
    branch = clean(manifest.get("branch"))
    output_dir = Path(clean(manifest.get("output_dir"))).resolve()
    paths = {
        "output_dir": output_dir,
        "variant_pip": Path(clean(manifest.get("variant_pip"))).resolve(),
        "credible_sets": Path(clean(manifest.get("credible_sets"))).resolve(),
        "credible_set_members": Path(clean(manifest.get("credible_set_members"))).resolve(),
    }
    run_root = run_dir.resolve()
    for label, path in paths.items():
        try:
            path.relative_to(run_root)
        except ValueError as exc:
            raise ValueError(f"manifest {label} not in frozen GWAS run in: {path}") from exc
        if label == "output_dir":
            if not path.is_dir():
                raise FileNotFoundError(f"manifest output_dir does not exist: {path}")
        else:
            require_file(path, f"GWAS locus {label}")

    completion_path = output_dir / "COMPLETED.ok"
    status_path = output_dir / "run_status.tsv"
    completion_fields, completion_rows = read_tsv(completion_path)
    require_columns(
        completion_path,
        completion_fields,
        {"overall_status", "locus_id", "branch"},
    )
    if len(completion_rows) != 1 or (
        clean(completion_rows[0].get("overall_status")),
        clean(completion_rows[0].get("locus_id")),
        clean(completion_rows[0].get("branch")),
    ) != ("PASS", locus_id, branch):
        raise ValueError(f"locus COMPLETED.ok with manifest mismatch: {locus_id}|{branch}")
    status_fields, status_rows = read_tsv(status_path)
    require_columns(
        status_path,
        status_fields,
        {"locus_id", "branch", "status", "overall_status", "model_status"},
    )
    if len(status_rows) != 1 or (
        clean(status_rows[0].get("locus_id")),
        clean(status_rows[0].get("branch")),
        clean(status_rows[0].get("status")),
        clean(status_rows[0].get("overall_status")),
        clean(status_rows[0].get("model_status")),
    ) != (locus_id, branch, "SUCCESS", "PASS", "CONVERGED"):
        raise ValueError(f"locus run_status.tsv with manifest mismatch: {locus_id}|{branch}")
    return paths


def build_gwas_record(
    row: Mapping[str, str],
    candidate_type: str,
    manifest: Mapping[str, str],
) -> Dict[str, object]:

    locus_id = clean(manifest.get("locus_id"))
    trait = clean(manifest.get("trait"))
    branch = clean(manifest.get("branch"))
    analysis_scope = clean(manifest.get("analysis_scope"))
    chrom = normalize_chrom(row.get("chrom"))
    pos = parse_int(row.get("pos"), f"{locus_id} pos")
    ld_a1 = clean(row.get("ld_a1")).upper()
    ld_a2 = clean(row.get("ld_a2")).upper()
    alleles = allele_tuple(ld_a1, ld_a2)
    marker_id = clean(row.get("ld_variant_id"))
    if chrom not in AUTOSOMES or pos <= 0 or len(alleles) != 2:
        raise ValueError(
            f"GWAS fine-map coordinates or allele invalid: {locus_id}, {chrom}:{pos}, {alleles}"
        )
    if not marker_id or not trait:
        raise ValueError(f"GWAS fine-map missing key identifiers: {locus_id}")
    p_lrt = finite_float(row.get("p_lrt"), f"{locus_id} p_lrt")
    pip = finite_float(row.get("pip"), f"{locus_id} pip")
    if not 0 <= p_lrt <= 1 or not 0 <= pip <= 1:
        raise ValueError(f"GWAS p_lrt or PIP exceeds [0,1]: {locus_id}|{marker_id}")
    variant_count = parse_int(manifest.get("variant_count"), f"{locus_id} variant_count")
    pip_rank = parse_int(row.get("pip_rank"), f"{locus_id} pip_rank")
    if pip_rank < 1 or pip_rank > variant_count:
        raise ValueError(f"GWAS pip_rank exceeds locus variant count: {locus_id}|{marker_id}")
    candidate_key = "|".join(
        (branch, analysis_scope, candidate_type, locus_id, trait, marker_id)
    )
    return {
        "candidate_type": candidate_type,
        "gwas_analysis_scope": analysis_scope,
        "gwas_branch": branch,
        "trait": trait,
        "locus_id": locus_id,
        "gwas_ld_variant_id": marker_id,
        "gwas_chrom": chrom,
        "gwas_pos": pos,
        "gwas_ld_a1": ld_a1,
        "gwas_ld_a2": ld_a2,
        "gwas_unordered_alleles": "/".join(alleles),
        "gwas_beta_for_ld_a1": clean(row.get("beta_for_ld_a1")),
        "gwas_se": clean(row.get("se")),
        "gwas_z_for_ld_a1": clean(row.get("z_for_ld_a1")),
        "gwas_variant_p_lrt": clean(row.get("p_lrt")),
        "gwas_n_obs": clean(row.get("n_obs")),
        "gwas_pip": clean(row.get("pip")),
        "gwas_pip_rank": clean(row.get("pip_rank")),
        "gwas_credible_sets": clean(row.get("credible_sets")),
        "gwas_in_credible_set": "yes" if candidate_type == "susie_credible_set" else "no",
        "gwas_ambiguous_at_cg": clean(row.get("ambiguous_at_cg")).lower(),
        "gwas_harmonization_status": clean(row.get("harmonization_status")),
        "gwas_locus_status": clean(manifest.get("status")),
        "gwas_model_status": clean(manifest.get("model_status")),
        "gwas_variant_count": str(variant_count),
        "gwas_reference_sample_count": clean(manifest.get("reference_sample_count")),
        "gwas_n_locus_median": clean(manifest.get("n_locus_median")),
        "gwas_locus_lead_variant_id": clean(manifest.get("gwas_lead_variant_id")),
        "gwas_locus_lead_pos": clean(manifest.get("gwas_lead_pos")),
        "gwas_locus_lead_p_lrt": clean(manifest.get("gwas_lead_p_lrt")),
        "gwas_locus_top_variant_id": clean(manifest.get("top_variant_id")),
        "gwas_locus_top_pos": clean(manifest.get("top_variant_pos")),
        "gwas_locus_max_pip": clean(manifest.get("max_pip")),
        "candidate_key": candidate_key,
        "_alleles": alleles,
    }


def load_gwas_candidates(
    run_dir: Path,
) -> Tuple[Dict[str, Dict[str, object]], Dict[str, Dict[str, int]]]:

    manifest_path = run_dir / "locus_run_manifest.tsv"
    manifest_fields, manifest_rows = read_tsv(manifest_path)
    require_columns(
        manifest_path,
        manifest_fields,
        GWAS_MANIFEST_REQUIRED_COLUMNS,
    )
    if not manifest_rows:
        raise ValueError("GWAS locus_run_manifest.tsv is empty")
    branch_expectations = validate_gwas_top_level(run_dir, manifest_rows)
    branch_candidates: Dict[str, List[Dict[str, object]]] = {
        branch: [] for branch in BRANCHES
    }
    seen_manifest: Set[Tuple[str, str]] = set()
    seen_candidate_keys: Set[str] = set()
    branch_scope_loci: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    branch_credible_membership_rows = Counter()

    for manifest_index, manifest in enumerate(manifest_rows, start=1):
        locus_id = clean(manifest.get("locus_id"))
        trait = clean(manifest.get("trait"))
        branch = clean(manifest.get("branch"))
        analysis_scope = clean(manifest.get("analysis_scope"))
        identity = (locus_id, branch)
        if not locus_id or not trait or identity in seen_manifest:
            raise ValueError(f"GWAS manifest contains empty or duplicate locus_id+branch: {identity}")
        seen_manifest.add(identity)
        if branch not in BRANCHES or analysis_scope not in ALLOWED_ANALYSIS_SCOPES:
            raise ValueError(f"GWAS manifest branch or analysis_scope invalid: {identity}|{analysis_scope}")
        if clean(manifest.get("status")) != "PASS" or clean(manifest.get("model_status")) != "CONVERGED":
            raise ValueError(f"GWAS manifest locus is not PASS/CONVERGED: {identity}")
        if clean(manifest.get("coverage")) != "0.95":
            raise ValueError(f"GWAS manifest locus coverage is not 0.95: {identity}")
        paths = validate_locus_artifacts(run_dir, manifest)
        variant_count = parse_int(manifest.get("variant_count"), f"{identity} variant_count")
        credible_set_count = parse_int(
            manifest.get("credible_set_count"), f"{identity} credible_set_count"
        )
        branch_scope_loci[(branch, analysis_scope)].add(locus_id)

        expected_members: Set[Tuple[str, str, str]] = set()
        observed_variant_count = 0
        for row in iter_tsv(paths["variant_pip"], GWAS_VARIANT_REQUIRED_COLUMNS):
            observed_variant_count += 1
            if (
                clean(row.get("locus_id")) != locus_id
                or clean(row.get("trait")) != trait
                or clean(row.get("analysis_scope")) != analysis_scope
            ):
                raise ValueError(f"variant_pip identity differs from the manifest: {identity}")
            ambiguous = clean(row.get("ambiguous_at_cg")).lower()
            harmonization = clean(row.get("harmonization_status"))
            if ambiguous not in {"yes", "no"} or harmonization not in ALLOWED_HARMONIZATION:
                raise ValueError(f"GWAS ambiguous/harmonize status invalid: {identity}")
            if branch == "strict_no_ambiguous" and ambiguous != "no":
                raise ValueError(f"strict_no_ambiguous still contains ambiguous site: {identity}")
            if ambiguous == "yes" and harmonization != "HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED":
                raise ValueError(f"ambiguous sites are not marked as resolved by frequency: {identity}")
            if ambiguous == "no" and harmonization != "HARMONIZED_NONAMBIGUOUS":
                raise ValueError(f"unambiguous site harmonize status mismatch: {identity}")
            in_cs = is_yes(row.get("in_credible_set"))
            cs_labels = split_cs_labels(row.get("credible_sets"))
            if in_cs != bool(cs_labels):
                raise ValueError(f"variant_pip in_credible_set with credible_sets mismatch: {identity}")
            for cs_id in cs_labels:
                expected_members.add(
                    (
                        clean(row.get("matrix_variant_index")),
                        clean(row.get("ld_variant_id")),
                        cs_id,
                    )
                )
            candidate_type = (
                "susie_credible_set" if in_cs else "susie_noncredible_locus_variant"
            )
            record = build_gwas_record(row, candidate_type, manifest)
            key = str(record["candidate_key"])
            if key in seen_candidate_keys:
                raise ValueError(f"duplicate GWAS candidate_key: {key}")
            seen_candidate_keys.add(key)
            branch_candidates[branch].append(record)
        if observed_variant_count != variant_count:
            raise ValueError(
                f"variant_pip row count differs from manifest variant_count: "
                f"{identity}; observed={observed_variant_count}; expected={variant_count}"
            )

        member_fields, member_rows = read_tsv(paths["credible_set_members"])
        require_columns(
            paths["credible_set_members"],
            member_fields,
            GWAS_MEMBER_REQUIRED_COLUMNS,
        )
        observed_members: Set[Tuple[str, str, str]] = set()
        for row in member_rows:
            if (
                clean(row.get("locus_id")) != locus_id
                or clean(row.get("trait")) != trait
                or clean(row.get("analysis_scope")) != analysis_scope
                or clean(row.get("branch")) != branch
            ):
                raise ValueError(f"credible_set_members identity differs from the manifest: {identity}")
            member = (
                clean(row.get("matrix_variant_index")),
                clean(row.get("ld_variant_id")),
                clean(row.get("cs_id")),
            )
            if not all(member) or member in observed_members:
                raise ValueError(f"credible_set_members contains empty or duplicate member: {identity}|{member}")
            observed_members.add(member)
        if observed_members != expected_members:
            raise ValueError(
                f"variant_pip with credible_set_members member mismatch: {identity}; "
                f"missing={len(expected_members - observed_members)}; "
                f"extra={len(observed_members - expected_members)}"
            )
        branch_credible_membership_rows[branch] += len(observed_members)

        cs_fields, cs_rows = read_tsv(paths["credible_sets"])
        require_columns(paths["credible_sets"], cs_fields, GWAS_CS_REQUIRED_COLUMNS)
        if len(cs_rows) != credible_set_count:
            raise ValueError(f"credible_sets row count and manifest mismatch: {identity}")
        for row in cs_rows:
            if (
                clean(row.get("locus_id")) != locus_id
                or clean(row.get("trait")) != trait
                or clean(row.get("analysis_scope")) != analysis_scope
                or clean(row.get("branch")) != branch
                or clean(row.get("coverage")) != "0.95"
            ):
                raise ValueError(f"credible_sets with manifest identity or coverage mismatch: {identity}")
        if manifest_index % 25 == 0 or manifest_index == len(manifest_rows):
            log(f"[stage 3/6] checked GWAS locus={manifest_index}/{len(manifest_rows)}")

    branch_results: Dict[str, Dict[str, object]] = {}
    branch_metrics: Dict[str, Dict[str, int]] = {}
    for branch in BRANCHES:
        candidates = branch_candidates[branch]
        by_chrom: Dict[str, List[Dict[str, object]]] = defaultdict(list)
        by_position: Dict[Tuple[str, int], List[Dict[str, object]]] = defaultdict(list)
        for row in candidates:
            chrom = str(row["gwas_chrom"])
            pos = int(row["gwas_pos"])
            by_chrom[chrom].append(row)
            by_position[(chrom, pos)].append(row)
        positions: Dict[str, List[int]] = {}
        for chrom, rows in by_chrom.items():
            rows.sort(key=lambda row: (int(row["gwas_pos"]), str(row["candidate_key"])))
            positions[chrom] = [int(row["gwas_pos"]) for row in rows]
        type_counts = Counter(str(row["candidate_type"]) for row in candidates)
        observed_cs_count = sum(
            parse_int(row.get("credible_set_count"), "credible_set_count")
            for row in manifest_rows
            if clean(row.get("branch")) == branch
        )
        if len(candidates) != branch_expectations[branch]["modeled_variant_count"]:
            raise ValueError(f"{branch} modeling variant count with branch_summary mismatch")
        if observed_cs_count != branch_expectations[branch]["credible_set_count"]:
            raise ValueError(f"{branch} credible set count and branch_summary mismatch")
        branch_results[branch] = {
            "candidates": candidates,
            "by_chrom": dict(by_chrom),
            "positions": positions,
            "by_position": dict(by_position),
            "candidate_counts": dict(type_counts),
        }
        branch_metrics[branch] = {
            "locus_count": sum(clean(row.get("branch")) == branch for row in manifest_rows),
            "primary_locus_count": len(branch_scope_loci[(branch, "primary")]),
            "exploratory_locus_count": len(branch_scope_loci[(branch, "exploratory")]),
            "variant_count": len(candidates),
            "credible_candidate_count": type_counts.get("susie_credible_set", 0),
            "noncredible_candidate_count": type_counts.get("susie_noncredible_locus_variant", 0),
            "credible_set_count": observed_cs_count,
            "credible_membership_rows": branch_credible_membership_rows[branch],
        }
    return branch_results, branch_metrics


def parse_unannotated_variant(endpoint: Mapping[str, object]) -> Dict[str, object]:
    match = None
    for value in (
        endpoint.get("gene_id"),
        endpoint.get("gene_symbol"),
        endpoint.get("name"),
    ):
        match = VARIANT_RE.search(clean(value))
        if match:
            break
    if match is None:
        raise ValueError(f"cannot parse the unannotated endpoint CHROM:POS_REF>ALT: {endpoint}")
    chrom = normalize_chrom(match.group(1))
    pos = int(match.group(2))
    ref = match.group(3).upper()
    alt = match.group(4).upper()
    endpoint_chrom = normalize_chrom(endpoint.get("chrom"))
    start = parse_int(endpoint.get("start"), "unannotated endpoint start")
    end = parse_int(endpoint.get("end"), "unannotated endpoint end")
    if endpoint_chrom != chrom or start != pos - 1 or end != pos:
        raise ValueError(
            f"unannotated endpoint fields differ: parsed={chrom}:{pos}_{ref}>{alt}, "
            f"columns={endpoint_chrom}:{start}-{end}"
        )
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "alleles": allele_tuple(ref, alt),
    }


def endpoint_from_row(row: Mapping[str, str], prefix: str) -> Dict[str, object]:
    return {
        "record_type": clean(row.get(f"{prefix}_record_type")),
        "gene_id": clean(row.get(f"{prefix}_gene_id")),
        "gene_symbol": clean(row.get(f"{prefix}_gene_symbol")),
        "gene_strand": clean(row.get(f"{prefix}_gene_strand")),
        "chrom": normalize_chrom(row.get(f"chr_{prefix}")),
        "start": clean(row.get(f"start_{prefix}")),
        "end": clean(row.get(f"end_{prefix}")),
        "name": clean(row.get(f"name_{prefix}")),
        "annotation_source": clean(row.get(f"{prefix}_annotation_source")),
        "gene_biotype": clean(row.get(f"{prefix}_gene_biotype")),
        "endpoint_type": clean(row.get(f"{prefix}_endpoint_type")),
    }


def annotate_endpoint(
    endpoint: Mapping[str, object],
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
) -> Tuple[Dict[str, str], Dict[str, object] | None]:
    record_type = clean(endpoint.get("record_type"))
    if record_type == "manual_new_gene":

        hits = [r for r in find_new_gene_loci(new_gene_bins, str(endpoint["chrom"]), int(endpoint["start"]) + 1)
                if r["gene_id"] == endpoint["gene_id"]]
        if len(hits) != 1 or (hits[0]["start0"], hits[0]["end0"], hits[0]["strand"]) != (
                int(endpoint["start"]), int(endpoint["end"]), endpoint["gene_strand"]):
            raise ValueError(f"New-gene endpoint does not match the GTF: {endpoint['gene_id']}")
        return summarize_new_gene_hits(hits), None
    if record_type == "known_gene":
        return summarize_new_gene_hits([]), None
    if record_type != "unannotated_variant":
        raise ValueError(f"unsupported endpoint record_type: {record_type}")
    variant = parse_unannotated_variant(endpoint)
    hits = find_new_gene_loci(
        new_gene_bins, str(variant["chrom"]), int(variant["pos"])
    )
    return summarize_new_gene_hits(hits), variant


def match_endpoint(
    endpoint: Mapping[str, object],
    by_chrom: Mapping[str, Sequence[Dict[str, object]]],
    positions: Mapping[str, Sequence[int]],
    by_position: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
) -> Tuple[List[Dict[str, object]], str]:
    record_type = clean(endpoint.get("record_type"))
    chrom = normalize_chrom(endpoint.get("chrom"))
    if record_type == "unannotated_variant":
        variant = parse_unannotated_variant(endpoint)
        hits = [
            row
            for row in by_position.get(
                (str(variant["chrom"]), int(variant["pos"])), []
            )
            if tuple(row["_alleles"]) == tuple(variant["alleles"])
        ]
        return hits, "point_exact_chrom_pos_unordered_allele_set"
    if record_type in {"known_gene", "manual_new_gene"}:
        start0 = parse_int(endpoint.get("start"), "annotated endpoint START0")
        end0 = parse_int(endpoint.get("end"), "annotated endpoint END0")
        if chrom not in AUTOSOMES or start0 < 0 or end0 <= start0:
            raise ValueError(
                f"annotated gene endpoint interval invalid: {chrom}:{start0}-{end0}"
            )
        chrom_rows = by_chrom.get(chrom, [])
        chrom_positions = positions.get(chrom, [])
        left = bisect.bisect_right(chrom_positions, start0)
        right = bisect.bisect_right(chrom_positions, end0)
        return (
            list(chrom_rows[left:right]),
            "gene_span_start0_end0_contains_1based_gwas_pos",
        )
    raise ValueError(f"unsupported endpoint record_type: {record_type}")


def endpoint_output(
    prefix: str,
    endpoint: Mapping[str, object],
    annotation: Mapping[str, str],
) -> Dict[str, str]:
    result = {
        f"{prefix}_record_type": clean(endpoint.get("record_type")),
        f"{prefix}_gene_id": clean(endpoint.get("gene_id")),
        f"{prefix}_gene_symbol": clean(endpoint.get("gene_symbol")),
        f"{prefix}_gene_strand": clean(endpoint.get("gene_strand")),
        f"{prefix}_chr": clean(endpoint.get("chrom")),
        f"{prefix}_start": clean(endpoint.get("start")),
        f"{prefix}_end": clean(endpoint.get("end")),
        f"{prefix}_annotation_source": clean(endpoint.get("annotation_source")),
        f"{prefix}_gene_biotype": clean(endpoint.get("gene_biotype")),
        f"{prefix}_endpoint_type": clean(endpoint.get("endpoint_type")),
    }
    result.update(
        {f"{prefix}_{key}": value for key, value in annotation.items()}
    )
    return result


def gwas_output(row: Mapping[str, object]) -> Dict[str, object]:
    return {column: row.get(column, "") for column in GWAS_COLUMNS}


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


def new_gene_scope(
    endpoint_hit: str,
    base_annotation: Mapping[str, str],
    target_annotation: Mapping[str, str],
) -> Tuple[str, str, str]:
    base_has = base_annotation.get("manual_new_gene_hit") == "yes"
    target_has = target_annotation.get("manual_new_gene_hit") == "yes"
    link_has = base_has or target_has
    hit_has = (
        (endpoint_hit in {"base", "both"} and base_has)
        or (endpoint_hit in {"target", "both"} and target_has)
    )
    if hit_has:
        scope = "gwas_hit_endpoint_with_new_gene_annotation"
    elif link_has:
        scope = "link_context_only_gwas_hits_other_endpoint"
    else:
        scope = "no_manual_new_gene_context"
    return ("yes" if link_has else "no", "yes" if hit_has else "no", scope)


def candidates_in_interval(
    chrom: str,
    start0: int,
    end0: int,
    by_chrom: Mapping[str, Sequence[Dict[str, object]]],
    positions: Mapping[str, Sequence[int]],
) -> List[Dict[str, object]]:

    chrom_rows = by_chrom.get(chrom, [])
    chrom_positions = positions.get(chrom, [])
    left = bisect.bisect_right(chrom_positions, start0)
    right = bisect.bisect_right(chrom_positions, end0)
    return list(chrom_rows[left:right])


def unannotated_pcgene_pair(
    base: Mapping[str, object],
    target: Mapping[str, object],
) -> Tuple[str, Mapping[str, object], str, Mapping[str, object]] | None:

    for unannotated_role, unannotated, pcgene_role, pcgene in (
        ("base", base, "target", target),
        ("target", target, "base", base),
    ):
        if (
            clean(unannotated.get("record_type")) == "unannotated_variant"
            and clean(pcgene.get("record_type")) == "known_gene"
            and clean(pcgene.get("gene_biotype")) == "protein_coding"
        ):
            return unannotated_role, unannotated, pcgene_role, pcgene
    return None


def point_in_endpoint_gene(pos1: int, endpoint: Mapping[str, object]) -> bool:
    start0 = parse_int(endpoint.get("start"), "paired pcGene START0")
    end0 = parse_int(endpoint.get("end"), "paired pcGene END0")
    if start0 < 0 or end0 <= start0:
        raise ValueError(f"paired pcGene interval invalid: {start0}-{end0}")
    return start0 < pos1 <= end0


def target_structure_output(
    candidate: Mapping[str, object],
    row: Mapping[str, str],
    base: Mapping[str, object],
    target: Mapping[str, object],
    base_annotation: Mapping[str, str],
    target_annotation: Mapping[str, str],
    unannotated_role: str,
    unannotated: Mapping[str, object],
    pcgene_role: str,
    pcgene: Mapping[str, object],
    ase_variant: Mapping[str, object],
    shared_new_gene: Mapping[str, object],
    mapping_class: str,
) -> Dict[str, object]:
    exact = mapping_class == "exact_variant_match"
    candidate_type = clean(candidate.get("candidate_type"))
    gwas_pos = int(candidate["gwas_pos"])
    gwas_in_pcgene = point_in_endpoint_gene(gwas_pos, pcgene)
    if candidate_type == "susie_noncredible_locus_variant":
        evidence_tier = (
            "Diagnostic_noncredible_exact"
            if exact
            else "Diagnostic_noncredible_same_new_gene"
        )
    else:
        evidence_tier = "Tier1_exact" if exact else "Tier2_same_new_gene"
    output: Dict[str, object] = {}
    output.update(gwas_output(candidate))
    output.update(
        {
            "canonical_pair_id": clean(row.get("canonical_pair_id")),
            **{k: row[k] for k in ("sample_pair_id", "source_link_file", "source_link_line", "source_metadata_file",
                                "overlap_bp", "strand_relation", "either_endpoint_multigene",
                                "both_endpoints_unambiguous_count_positive", "identical_haplotype_counts",
                                "evidence_stratum", "shared_variant_identity_checked")},
            "animal_id": clean(row.get("animal_id")),
            "tissue": clean(row.get("tissue")),
            "sample_tissue": clean(row.get("sample_tissue")),
            "link_chromosome": normalize_chrom(row.get("chromosome")),
            "pair_type": clean(row.get("pair_type")),
            "interval_gap_bp": clean(row.get("interval_gap_bp")),
            "mechanism": clean(row.get("mechanism")),
            "linkage_score": clean(row.get("linkage_score")),
            "endpoint_hit": unannotated_role if exact else "",
            "base_mapping_rule": (
                "point_exact_chrom_pos_unordered_allele_set"
                if exact and unannotated_role == "base"
                else ""
            ),
            "target_mapping_rule": (
                "point_exact_chrom_pos_unordered_allele_set"
                if exact and unannotated_role == "target"
                else ""
            ),
        }
    )
    output.update(endpoint_output("base", base, base_annotation))
    output.update(endpoint_output("target", target, target_annotation))
    output.update(
        {
            "link_contains_manual_new_gene": "yes",
            "gwas_hit_endpoint_has_manual_new_gene": "yes" if exact else "no",
            "new_gene_evidence_scope": (
                "exact_gwas_ase_variant_in_same_manual_new_gene"
                if exact
                else "nonexact_gwas_ase_in_same_manual_new_gene"
            ),
            "evidence_tier": evidence_tier,
            "ase_gwas_mapping_class": mapping_class,
            "unannotated_endpoint_role": unannotated_role,
            "paired_pcgene_endpoint_role": pcgene_role,
            "shared_manual_new_gene_id": clean(
                shared_new_gene.get("gene_id")
            ),
            "shared_new_gene_chrom": clean(shared_new_gene.get("chrom")),
            "shared_new_gene_start0": shared_new_gene.get("start0", ""),
            "shared_new_gene_end0": shared_new_gene.get("end0", ""),
            "shared_new_gene_strand": clean(shared_new_gene.get("strand")),
            "shared_new_gene_length_bp": (
                int(shared_new_gene["end0"]) - int(shared_new_gene["start0"])
            ),
            "shared_new_gene_class_codes": join_values(
                shared_new_gene.get("class_codes", [])
            ),
            "shared_new_gene_novelty_categories": join_values(
                shared_new_gene.get("novelty", [])
            ),
            "shared_new_gene_candidate_category": clean(
                shared_new_gene.get("candidate_category")
            ),
            "shared_new_gene_dataset_support": clean(
                shared_new_gene.get("dataset_support")
            ),
            "shared_new_gene_filter_reason": clean(
                shared_new_gene.get("filter_reason")
            ),
            "unannotated_ase_chrom": clean(ase_variant.get("chrom")),
            "unannotated_ase_pos": ase_variant.get("pos", ""),
            "unannotated_ase_ref": clean(ase_variant.get("ref")),
            "unannotated_ase_alt": clean(ase_variant.get("alt")),
            "ase_gwas_distance_bp": abs(gwas_pos - int(ase_variant["pos"])),
            "gwas_in_paired_pcgene": "yes" if gwas_in_pcgene else "no",
            "gwas_is_locus_top_variant": (
                "yes"
                if clean(candidate.get("gwas_ld_variant_id"))
                == clean(candidate.get("gwas_locus_top_variant_id"))
                else "no"
            ),
            "target_structure_pass": "no" if gwas_in_pcgene else "yes",
        }
    )
    return output


def validate_allelome_summary(path: Path) -> None:

    read_two_column_table(path, "metric")


def process_links(
    links_path: Path,
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    by_chrom: Mapping[str, Sequence[Dict[str, object]]],
    positions: Mapping[str, Sequence[int]],
    by_position: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    run_dir: Path,
    candidates: Sequence[Mapping[str, object]],
) -> Dict[str, int]:
    output_paths = {
        "credible": run_dir / "dog5_allelome_gwas_credible_set_overlap.tsv.gz",
        "noncredible": (
            run_dir / "dog5_allelome_gwas_noncredible_diagnostic.tsv.gz"
        ),
        "direct_newgene": run_dir / "dog5_allelome_gwas_finemap_direct_new_gene_pcgene.tsv.gz",
        "newgene": (
            run_dir / "dog5_allelome_gwas_finemap_new_gene_context.tsv.gz"
        ),
        "hit_links": (
            run_dir / "dog5_allelome_gwas_finemap_hit_link_summary.tsv.gz"
        ),
        "target_exact": (
            run_dir / "dog5_allelome_gwas_finemap_exact_variant.tsv.gz"
        ),
        "target_nonexact": (
            run_dir
            / "dog5_allelome_gwas_finemap_same_new_gene_nonexact.tsv.gz"
        ),
        "target_union": (
            run_dir / "dog5_allelome_gwas_finemap_target_union.tsv.gz"
        ),
        "target_excluded": (
            run_dir
            / "dog5_allelome_gwas_finemap_excluded_gwas_in_paired_pcgene.tsv.gz"
        ),
        "target_noncredible": (
            run_dir
            / "dog5_allelome_gwas_finemap_same_new_gene_noncredible_diagnostic.tsv.gz"
        ),
    }
    handles: List[object] = []
    writers: Dict[str, csv.DictWriter] = {}
    for key in ("credible", "noncredible", "newgene", "direct_newgene"):
        handle, writer = open_gzip_writer(output_paths[key], OVERLAP_COLUMNS)
        handles.append(handle)
        writers[key] = writer
    hit_handle, hit_writer = open_gzip_writer(
        output_paths["hit_links"], HIT_LINK_SUMMARY_COLUMNS
    )
    handles.append(hit_handle)
    for key in (
        "target_exact",
        "target_nonexact",
        "target_union",
        "target_excluded",
        "target_noncredible",
    ):
        handle, writer = open_gzip_writer(
            output_paths[key], TARGET_STRUCTURE_COLUMNS
        )
        handles.append(handle)
        writers[key] = writer

    metrics = Counter()
    hit_candidate_keys: Set[str] = set()
    credible_candidate_types = {"susie_credible_set"}
    seen_links: Set[str] = set()

    try:
        with links_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(
                links_path, reader.fieldnames or [], LINK_REQUIRED_COLUMNS
            )
            for line_number, row in enumerate(reader, start=2):
                metrics["allelome_link_rows"] += 1
                ids = sorted((row["base_gene_id"], row["target_gene_id"]))
                row["canonical_pair_id"] = json.dumps([row["chr_base"], *ids], separators=(",", ":"))
                row["sample_pair_id"] = json.dumps([row["animal_id"], row["tissue"], row["chr_base"], *ids], separators=(",", ":"))
                row["sample_tissue"] = row["animal_id"] + "|" + row["tissue"]
                row["chromosome"] = row["chr_base"]
                if row["chr_base"] != row["chr_target"]:
                    raise ValueError("Allelome both endpoints across chromosome")
                canonical_pair_id = row["canonical_pair_id"]
                sample_pair_id = row["sample_pair_id"]
                if sample_pair_id in seen_links:
                    raise ValueError(
                        f"{links_path}:{line_number} sample_pair_id duplicate"
                    )
                seen_links.add(sample_pair_id)
                base = endpoint_from_row(row, "base")
                target = endpoint_from_row(row, "target")
                base_annotation, _ = annotate_endpoint(base, new_gene_bins)
                target_annotation, _ = annotate_endpoint(target, new_gene_bins)


                pair = unannotated_pcgene_pair(base, target)
                if pair is not None:
                    (
                        unannotated_role,
                        unannotated,
                        pcgene_role,
                        pcgene,
                    ) = pair
                    ase_variant = parse_unannotated_variant(unannotated)
                    shared_new_genes = find_new_gene_loci(
                        new_gene_bins,
                        str(ase_variant["chrom"]),
                        int(ase_variant["pos"]),
                    )
                    for shared_new_gene in shared_new_genes:
                        interval_candidates = candidates_in_interval(
                            str(shared_new_gene["chrom"]),
                            int(shared_new_gene["start0"]),
                            int(shared_new_gene["end0"]),
                            by_chrom,
                            positions,
                        )
                        for candidate in interval_candidates:
                            candidate_type = str(candidate["candidate_type"])
                            analysis_scope = str(candidate["gwas_analysis_scope"])
                            same_position = (
                                int(candidate["gwas_pos"])
                                == int(ase_variant["pos"])
                            )
                            same_alleles = (
                                tuple(candidate["_alleles"])
                                == tuple(ase_variant["alleles"])
                            )
                            if same_position and not same_alleles:
                                metrics[
                                    "target_same_position_allele_mismatch_rows"
                                ] += 1
                                continue
                            mapping_class = (
                                "exact_variant_match"
                                if same_position and same_alleles
                                else "same_new_gene_span_nonexact"
                            )
                            target_output = target_structure_output(
                                candidate,
                                row,
                                base,
                                target,
                                base_annotation,
                                target_annotation,
                                unannotated_role,
                                unannotated,
                                pcgene_role,
                                pcgene,
                                ase_variant,
                                shared_new_gene,
                                mapping_class,
                            )
                            if candidate_type not in credible_candidate_types:
                                writers["target_noncredible"].writerow(
                                    target_output
                                )
                                metrics[
                                    "target_noncredible_diagnostic_rows"
                                ] += 1
                                metrics[
                                    f"target_noncredible_diagnostic_rows_{analysis_scope}"
                                ] += 1
                                continue
                            if (
                                target_output["target_structure_pass"]
                                == "no"
                            ):
                                writers["target_excluded"].writerow(
                                    target_output
                                )
                                metrics[
                                    "target_excluded_gwas_in_pcgene_rows"
                                ] += 1
                                metrics[
                                    f"target_excluded_gwas_in_pcgene_rows_{analysis_scope}"
                                ] += 1
                            elif mapping_class == "exact_variant_match":
                                writers["target_exact"].writerow(target_output)
                                writers["target_union"].writerow(target_output)
                                metrics["target_exact_rows"] += 1
                                metrics["target_union_rows"] += 1
                                metrics[f"target_exact_rows_{analysis_scope}"] += 1
                            else:
                                writers["target_nonexact"].writerow(
                                    target_output
                                )
                                writers["target_union"].writerow(target_output)
                                metrics["target_nonexact_rows"] += 1
                                metrics["target_union_rows"] += 1
                                metrics[f"target_nonexact_rows_{analysis_scope}"] += 1

                base_hits, base_rule = match_endpoint(
                    base, by_chrom, positions, by_position
                )
                target_hits, target_rule = match_endpoint(
                    target, by_chrom, positions, by_position
                )
                grouped: MutableMapping[
                    str, Dict[str, object]
                ] = {}
                for endpoint_role, hits in (
                    ("base", base_hits),
                    ("target", target_hits),
                ):
                    for candidate in hits:
                        key = str(candidate["candidate_key"])
                        if key not in grouped:
                            grouped[key] = {
                                "candidate": candidate,
                                "roles": set(),
                            }
                        grouped[key]["roles"].add(endpoint_role)

                credible_rows: List[Dict[str, object]] = []
                for key in sorted(grouped):
                    candidate = grouped[key]["candidate"]
                    roles = grouped[key]["roles"]
                    endpoint_hit = (
                        "both"
                        if roles == {"base", "target"}
                        else next(iter(roles))
                    )
                    link_has, hit_has, scope = new_gene_scope(
                        endpoint_hit, base_annotation, target_annotation
                    )
                    output: Dict[str, object] = {}
                    output.update(gwas_output(candidate))
                    output.update(
                        {
                            "canonical_pair_id": canonical_pair_id,
                            **{k: row[k] for k in ("sample_pair_id", "source_link_file", "source_link_line", "source_metadata_file",
                                "overlap_bp", "strand_relation", "either_endpoint_multigene",
                                "both_endpoints_unambiguous_count_positive", "identical_haplotype_counts",
                                "evidence_stratum", "shared_variant_identity_checked")},
                            "animal_id": clean(row.get("animal_id")),
                            "tissue": clean(row.get("tissue")),
                            "sample_tissue": clean(row.get("sample_tissue")),
                            "link_chromosome": normalize_chrom(
                                row.get("chromosome")
                            ),
                            "pair_type": clean(
                                row.get("pair_type")
                            ),
                            "interval_gap_bp": clean(
                                row.get("interval_gap_bp")
                            ),
                            "mechanism": clean(row.get("mechanism")),
                            "linkage_score": clean(row.get("linkage_score")),
                            "endpoint_hit": endpoint_hit,
                            "base_mapping_rule": (
                                base_rule
                                if endpoint_hit in {"base", "both"}
                                else ""
                            ),
                            "target_mapping_rule": (
                                target_rule
                                if endpoint_hit in {"target", "both"}
                                else ""
                            ),
                        }
                    )
                    output.update(endpoint_output("base", base, base_annotation))
                    output.update(
                        endpoint_output("target", target, target_annotation)
                    )
                    output.update(
                        {
                            "link_contains_manual_new_gene": link_has,
                            "gwas_hit_endpoint_has_manual_new_gene": hit_has,
                            "new_gene_evidence_scope": scope,
                        }
                    )
                    if hit_has == "yes" and any(
                        e["record_type"] == "manual_new_gene" and role in roles
                        for role, e in (("base", base), ("target", target))
                    ):
                        output["new_gene_evidence_scope"] = "gwas_hit_candidate_new_gene_span"
                    candidate_type = str(candidate["candidate_type"])
                    analysis_scope = str(candidate["gwas_analysis_scope"])
                    if candidate_type == "susie_credible_set":
                        writers["credible"].writerow(output)

                        for role, e, other in (("base", base, target), ("target", target, base)):
                            if (e["record_type"] == "manual_new_gene" and endpoint_hit == role
                                    and other["record_type"] == "known_gene"
                                    and other["gene_biotype"] == "protein_coding"):
                                writers["direct_newgene"].writerow(output)
                                metrics["credible_direct_new_gene_pcgene_rows"] += 1
                        metrics["credible_overlap_rows"] += 1
                        metrics[f"credible_overlap_rows_{analysis_scope}"] += 1
                    else:
                        writers["noncredible"].writerow(output)
                        metrics["noncredible_overlap_rows"] += 1
                        metrics[f"noncredible_overlap_rows_{analysis_scope}"] += 1
                    hit_candidate_keys.add(str(candidate["candidate_key"]))
                    if candidate_type in credible_candidate_types:
                        credible_rows.append(output)
                        if scope != "no_manual_new_gene_context":
                            writers["newgene"].writerow(output)
                            metrics["credible_new_gene_context_rows"] += 1
                        if hit_has == "yes":
                            metrics["credible_direct_new_gene_rows"] += 1

                if credible_rows:
                    metrics["credible_hit_links"] += 1
                    credible_candidates = {
                        str(item["candidate_key"]): item for item in credible_rows
                    }
                    direct_rows = [
                        item
                        for item in credible_rows
                        if item["gwas_hit_endpoint_has_manual_new_gene"] == "yes"
                    ]
                    gene_ids = split_values(
                        *[
                            item.get("base_manual_new_gene_ids", "")
                            for item in credible_rows
                        ],
                        *[
                            item.get("target_manual_new_gene_ids", "")
                            for item in credible_rows
                        ],
                    )
                    hit_writer.writerow(
                        {
                            "canonical_pair_id": canonical_pair_id,
                            **{k: row[k] for k in ("sample_pair_id", "source_link_file", "source_link_line", "source_metadata_file",
                                "overlap_bp", "strand_relation", "either_endpoint_multigene",
                                "both_endpoints_unambiguous_count_positive", "identical_haplotype_counts",
                                "evidence_stratum", "shared_variant_identity_checked")},
                            "animal_id": clean(row.get("animal_id")),
                            "tissue": clean(row.get("tissue")),
                            "sample_tissue": clean(row.get("sample_tissue")),
                            "chromosome": normalize_chrom(row.get("chromosome")),
                            "pair_type": clean(
                                row.get("pair_type")
                            ),
                            "mechanism": clean(row.get("mechanism")),
                            "linkage_score": clean(row.get("linkage_score")),
                            "credible_overlap_row_n": len(credible_rows),
                            "credible_candidate_n": sum(
                                item["candidate_type"] == "susie_credible_set"
                                for item in credible_candidates.values()
                            ),
                            "analysis_scopes": join_values(
                                item["gwas_analysis_scope"]
                                for item in credible_candidates.values()
                            ),
                            "traits": join_values(
                                item["trait"]
                                for item in credible_candidates.values()
                            ),
                            "endpoint_hits": join_values(
                                item["endpoint_hit"] for item in credible_rows
                            ),
                            "minimum_gwas_variant_p_lrt": min(
                                parse_float(
                                    item["gwas_variant_p_lrt"],
                                    "gwas_variant_p_lrt",
                                )
                                for item in credible_rows
                            ),
                            "link_contains_manual_new_gene": (
                                "yes"
                                if any(
                                    item["link_contains_manual_new_gene"] == "yes"
                                    for item in credible_rows
                                )
                                else "no"
                            ),
                            "direct_new_gene_endpoint_hit_n": len(direct_rows),
                            "manual_new_gene_ids": "|".join(gene_ids),
                        }
                    )

                if metrics["allelome_link_rows"] % 50_000 == 0:
                    log(
                        "[stage 4/6] links="
                        f"{metrics['allelome_link_rows']}; "
                        f"credible_rows={metrics['credible_overlap_rows']}; "
                        f"diagnostic_rows={metrics['noncredible_overlap_rows']}"
                    )
    finally:
        for handle in handles:
            handle.close()

    credible_candidates = [
        row
        for row in candidates
        if row["candidate_type"] == "susie_credible_set"
    ]
    unmatched = [
        gwas_output(row)
        for row in credible_candidates
        if str(row["candidate_key"]) not in hit_candidate_keys
    ]
    unmatched_path = (
        run_dir / "dog5_allelome_gwas_finemap_unmatched_credible_set_candidates.tsv.gz"
    )
    unmatched_handle, unmatched_writer = open_gzip_writer(
        unmatched_path, GWAS_COLUMNS
    )
    try:
        unmatched_writer.writerows(unmatched)
    finally:
        unmatched_handle.close()
    metrics["unmatched_credible_set_candidates"] = len(unmatched)
    metrics["unique_hit_candidates_all_layers"] = len(hit_candidate_keys)
    return dict(metrics)


def write_candidate_inventory(
    path: Path, candidates: Sequence[Mapping[str, object]]
) -> None:
    handle, writer = open_gzip_writer(path, GWAS_COLUMNS)
    try:
        for row in candidates:
            writer.writerow(gwas_output(row))
    finally:
        handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare five-dog Allelome links across all chromosomes with stratified GWAS fine-mapping candidates."
    )
    parser.add_argument(
        "--allelome-links", type=Path, default=DEFAULT_ALLELOME_LINKS
    )
    parser.add_argument(
        "--allelome-summary", type=Path, default=DEFAULT_ALLELOME_SUMMARY
    )
    parser.add_argument(
        "--new-gene-gtf", type=Path, default=DEFAULT_NEW_GENE_GTF
    )
    parser.add_argument(
        "--new-gene-metadata", type=Path, default=DEFAULT_NEW_GENE_METADATA
    )
    parser.add_argument(
        "--new-gene-build-summary",
        type=Path,
        default=DEFAULT_NEW_GENE_BUILD_SUMMARY,
    )
    parser.add_argument(
        "--gwas-run-dir", type=Path, default=DEFAULT_GWAS_RUN_DIR
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    return parser.parse_args()


def main() -> int:
    global LOGGER
    args = parse_args()
    run_dir = args.output_root / args.run_id
    if run_dir.exists():
        raise FileExistsError(f"refuse to overwrite existing run directory: {run_dir}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir()
    LOGGER = RunLogger(run_dir / "1.1.compare_dog5_allelome_gwas_finemap.log")

    try:
        log("[stage 1/6] validate fixed input with upstream status")
        for path, label in (
            (args.allelome_links, "Allelome links"),
            (args.allelome_summary, "Allelome summary"),
            (args.new_gene_gtf, "Manual new-only GTF"),
            (args.new_gene_metadata, "Five-dog new-gene metadata table"),
            (args.new_gene_build_summary, "Manual GTF construction summary table"),
        ):
            require_file(path, label)
        if not args.gwas_run_dir.is_dir():
            raise FileNotFoundError(
                f"GWAS fine-map frozen run directory does not exist: {args.gwas_run_dir}"
            )
        validate_allelome_summary(args.allelome_summary)

        log("[stage 2/6] read Manual GTF gene span with five dogs transcript_map")
        new_gene_rows, new_gene_bins = load_new_gene_loci(
            args.new_gene_gtf,
            args.new_gene_metadata,
            args.new_gene_build_summary,
        )
        log(
            f"[stage 2/6] five_only_new_gene_rows={len(new_gene_rows)}; "
            f"autosomal_rows={sum(row['chrom'] in AUTOSOMES for row in new_gene_rows)}"
        )

        log("[stage 3/6] read and stratify validate GWAS fine-map candidate")
        branch_inputs, branch_gwas_metrics = load_gwas_candidates(args.gwas_run_dir)
        log(
            "[stage 3/6] GWAS branch validate completed: "
            + "; ".join(
                f"{branch}: loci={branch_gwas_metrics[branch]['locus_count']}, "
                f"variants={branch_gwas_metrics[branch]['variant_count']}, "
                f"credible={branch_gwas_metrics[branch]['credible_candidate_count']}"
                for branch in BRANCHES
            )
        )

        branch_process_metrics: Dict[str, Dict[str, int]] = {}
        for branch_index, branch in enumerate(BRANCHES, start=1):
            branch_dir = run_dir / branch
            branch_dir.mkdir()
            branch_data = branch_inputs[branch]
            candidates = branch_data["candidates"]
            write_candidate_inventory(
                branch_dir / "dog5_gwas_finemap_candidate_inventory.tsv.gz",
                candidates,
            )
            log(
                f"[stage 4/6] branch {branch_index}/{len(BRANCHES)} {branch}: "
                "streaming matching Allelome both endpoints and GWAS fine-map"
            )
            metrics = process_links(
                args.allelome_links,
                new_gene_bins,
                branch_data["by_chrom"],
                branch_data["positions"],
                branch_data["by_position"],
                branch_dir,
                candidates,
            )
            branch_process_metrics[branch] = metrics

            gwas_metrics = branch_gwas_metrics[branch]
            summary_rows = [
                ("overall_status", "PASS"),
                ("gwas_branch", branch),
                (
                    "branch_interpretation_role",
                    "primary_result" if branch == "main" else "sensitivity_analysis",
                ),
                ("allelome_link_rows", metrics["allelome_link_rows"]),
                ("unified_new_gene_rows", len(new_gene_rows)),
                ("gwas_locus_count", gwas_metrics["locus_count"]),
                ("gwas_primary_locus_count", gwas_metrics["primary_locus_count"]),
                ("gwas_exploratory_locus_count", gwas_metrics["exploratory_locus_count"]),
                ("gwas_variant_count", gwas_metrics["variant_count"]),
                ("gwas_credible_set_count", gwas_metrics["credible_set_count"]),
                ("gwas_credible_set_candidates", gwas_metrics["credible_candidate_count"]),
                ("gwas_credible_membership_rows", gwas_metrics["credible_membership_rows"]),
                (
                    "gwas_noncredible_locus_diagnostic_candidates",
                    gwas_metrics["noncredible_candidate_count"],
                ),
                ("credible_overlap_rows", metrics.get("credible_overlap_rows", 0)),
                ("credible_overlap_rows_primary_scope", metrics.get("credible_overlap_rows_primary", 0)),
                ("credible_overlap_rows_exploratory_scope", metrics.get("credible_overlap_rows_exploratory", 0)),
                (
                    "noncredible_diagnostic_overlap_rows",
                    metrics.get("noncredible_overlap_rows", 0),
                ),
                ("noncredible_diagnostic_overlap_rows_primary_scope", metrics.get("noncredible_overlap_rows_primary", 0)),
                ("noncredible_diagnostic_overlap_rows_exploratory_scope", metrics.get("noncredible_overlap_rows_exploratory", 0)),
                ("credible_hit_links", metrics.get("credible_hit_links", 0)),
                ("credible_new_gene_context_rows", metrics.get("credible_new_gene_context_rows", 0)),
                ("credible_direct_new_gene_rows", metrics.get("credible_direct_new_gene_rows", 0)),
                ("unmatched_credible_set_candidates", metrics.get("unmatched_credible_set_candidates", 0)),
                ("unique_hit_candidates_all_layers", metrics.get("unique_hit_candidates_all_layers", 0)),
                ("target_exact_rows", metrics.get("target_exact_rows", 0)),
                ("target_exact_rows_primary_scope", metrics.get("target_exact_rows_primary", 0)),
                ("target_exact_rows_exploratory_scope", metrics.get("target_exact_rows_exploratory", 0)),
                ("target_same_new_gene_nonexact_rows", metrics.get("target_nonexact_rows", 0)),
                ("target_same_new_gene_nonexact_rows_primary_scope", metrics.get("target_nonexact_rows_primary", 0)),
                ("target_same_new_gene_nonexact_rows_exploratory_scope", metrics.get("target_nonexact_rows_exploratory", 0)),
                ("target_union_rows", metrics.get("target_union_rows", 0)),
                (
                    "target_excluded_gwas_in_paired_pcgene_rows",
                    metrics.get("target_excluded_gwas_in_pcgene_rows", 0),
                ),
                (
                    "target_same_position_allele_mismatch_rows",
                    metrics.get("target_same_position_allele_mismatch_rows", 0),
                ),
                (
                    "target_noncredible_diagnostic_rows",
                    metrics.get("target_noncredible_diagnostic_rows", 0),
                ),
                ("credible_direct_new_gene_pcgene_rows", metrics.get("credible_direct_new_gene_pcgene_rows", 0)),
                ("gwas_susie_coverage", "0.95"),
                (
                    "interpretation_boundary",
                    "positional_prioritization_not_colocalization_or_causal_proof",
                ),
            ]
            write_plain_tsv(
                branch_dir / "dog5_allelome_gwas_finemap_summary.tsv",
                ["metric", "value"],
                ({"metric": key, "value": value} for key, value in summary_rows),
            )
            branch_manifest_rows = [
                ("run_id", args.run_id),
                ("run_directory", run_dir),
                ("branch_directory", branch_dir),
                ("gwas_branch", branch),
                ("branch_interpretation_role", "primary_result" if branch == "main" else "sensitivity_analysis"),
                ("allelome_links", args.allelome_links),
                ("allelome_summary", args.allelome_summary),
                ("new_gene_gtf", args.new_gene_gtf),
                ("new_gene_metadata", args.new_gene_metadata),
                ("new_gene_build_summary", args.new_gene_build_summary),
                ("gwas_run_dir", args.gwas_run_dir),
                ("gwas_locus_run_manifest", args.gwas_run_dir / "locus_run_manifest.tsv"),
                ("gwas_branch_summary", args.gwas_run_dir / "branch_summary.tsv"),
                ("gwas_model_parameters", args.gwas_run_dir / "model_parameters.tsv"),
                ("annotated_gene_mapping_rule", "START0<POS1<=END0"),
                ("unannotated_endpoint_mapping_rule", "exact_CHROM_POS_unordered_allele_set"),
                (
                    "new_gene_annotation_rule",
                    "unannotated_variant_point0=POS1-1; "
                    "GTF_transcript_union_gene_start0<=point0<"
                    "GTF_transcript_union_gene_end0",
                ),
                (
                    "new_gene_annotation_source",
                    "Manual new-only GTF transcript records grouped by gene_id; "
                    "five-only transcript_map retained as metadata and boundary cross-check",
                ),
                (
                    "new_gene_interpretation",
                    "manual_new_gene endpoints already aggregated; unannotated_variant post_hoc annotation remains separate",
                ),
                ("primary_candidate_layers", "susie_95pct_credible_set"),
                (
                    "primary_overlap_output",
                    branch_dir / "dog5_allelome_gwas_credible_set_overlap.tsv.gz",
                ),
                ("direct_new_gene_pcgene_output", branch_dir / "dog5_allelome_gwas_finemap_direct_new_gene_pcgene.tsv.gz"),
                ("direct_new_gene_pcgene_rule", "credible SNP hits manual_new_gene span only; paired known protein_coding gene not hit; not exact ASE SNP evidence"),
                ("primary_overlap_definition", "susie_95pct_credible_set_overlap"),
                (
                    "target_exact_output",
                    branch_dir / "dog5_allelome_gwas_finemap_exact_variant.tsv.gz",
                ),
                (
                    "target_same_new_gene_nonexact_output",
                    branch_dir / "dog5_allelome_gwas_finemap_same_new_gene_nonexact.tsv.gz",
                ),
                (
                    "target_union_output",
                    branch_dir / "dog5_allelome_gwas_finemap_target_union.tsv.gz",
                ),
                (
                    "target_excluded_output",
                    branch_dir / "dog5_allelome_gwas_finemap_excluded_gwas_in_paired_pcgene.tsv.gz",
                ),
                (
                    "target_noncredible_diagnostic_output",
                    branch_dir / "dog5_allelome_gwas_finemap_same_new_gene_noncredible_diagnostic.tsv.gz",
                ),
                (
                    "target_structure_definition",
                    "unannotated_ASE--pcGene_link; ASE_and_GWAS_share_Manual_new_gene; "
                    "GWAS_outside_paired_pcGene",
                ),
                (
                    "target_evidence_hierarchy",
                    "Tier1_exact_CHROM_POS_unordered_allele_set;"
                    "Tier2_same_Manual_new_gene_nonexact",
                ),
                ("diagnostic_candidate_layer", "susie_noncredible_locus_variant"),
                (
                    "allelome_direction_interpretation",
                    "base_target_are_canonical_pair_endpoints_not_causal_direction",
                ),
            ]
            write_plain_tsv(
                branch_dir / "dog5_allelome_gwas_finemap_manifest.tsv",
                ["item", "value"],
                ({"item": key, "value": value} for key, value in branch_manifest_rows),
            )
            log(
                f"[stage 5/6] {branch} summary: "
                f"credible_rows={metrics.get('credible_overlap_rows', 0)}; "
                f"credible_hit_links={metrics.get('credible_hit_links', 0)}; "
                f"target_exact={metrics.get('target_exact_rows', 0)}; "
                f"target_nonexact={metrics.get('target_nonexact_rows', 0)}"
            )

        comparison_columns = [
            "branch", "interpretation_role", "locus_count", "primary_locus_count",
            "exploratory_locus_count", "variant_count", "credible_set_count",
            "credible_candidate_count", "credible_membership_rows",
            "noncredible_candidate_count", "credible_overlap_rows",
            "credible_overlap_rows_primary_scope", "credible_overlap_rows_exploratory_scope",
            "noncredible_diagnostic_overlap_rows", "credible_hit_links",
            "target_exact_rows", "target_exact_rows_primary_scope",
            "target_exact_rows_exploratory_scope", "target_same_new_gene_nonexact_rows",
            "target_same_new_gene_nonexact_rows_primary_scope",
            "target_same_new_gene_nonexact_rows_exploratory_scope",
        ]
        comparison_rows = []
        for branch in BRANCHES:
            gwas_metrics = branch_gwas_metrics[branch]
            metrics = branch_process_metrics[branch]
            comparison_rows.append(
                {
                    "branch": branch,
                    "interpretation_role": "primary_result" if branch == "main" else "sensitivity_analysis",
                    "locus_count": gwas_metrics["locus_count"],
                    "primary_locus_count": gwas_metrics["primary_locus_count"],
                    "exploratory_locus_count": gwas_metrics["exploratory_locus_count"],
                    "variant_count": gwas_metrics["variant_count"],
                    "credible_set_count": gwas_metrics["credible_set_count"],
                    "credible_candidate_count": gwas_metrics["credible_candidate_count"],
                    "credible_membership_rows": gwas_metrics["credible_membership_rows"],
                    "noncredible_candidate_count": gwas_metrics["noncredible_candidate_count"],
                    "credible_overlap_rows": metrics.get("credible_overlap_rows", 0),
                    "credible_overlap_rows_primary_scope": metrics.get("credible_overlap_rows_primary", 0),
                    "credible_overlap_rows_exploratory_scope": metrics.get("credible_overlap_rows_exploratory", 0),
                    "noncredible_diagnostic_overlap_rows": metrics.get("noncredible_overlap_rows", 0),
                    "credible_hit_links": metrics.get("credible_hit_links", 0),
                    "target_exact_rows": metrics.get("target_exact_rows", 0),
                    "target_exact_rows_primary_scope": metrics.get("target_exact_rows_primary", 0),
                    "target_exact_rows_exploratory_scope": metrics.get("target_exact_rows_exploratory", 0),
                    "target_same_new_gene_nonexact_rows": metrics.get("target_nonexact_rows", 0),
                    "target_same_new_gene_nonexact_rows_primary_scope": metrics.get("target_nonexact_rows_primary", 0),
                    "target_same_new_gene_nonexact_rows_exploratory_scope": metrics.get("target_nonexact_rows_exploratory", 0),
                }
            )
        write_plain_tsv(
            run_dir / "dog5_allelome_gwas_finemap_branch_comparison.tsv",
            comparison_columns,
            comparison_rows,
        )
        write_plain_tsv(
            run_dir / "dog5_allelome_gwas_finemap_run_manifest.tsv",
            ["item", "value"],
            (
                {"item": key, "value": value}
                for key, value in [
                    ("run_id", args.run_id),
                    ("run_directory", run_dir),
                    ("gwas_run_dir", args.gwas_run_dir),
                    ("branch_policy", "main_primary_result;strict_no_ambiguous_sensitivity;never_sum_branches"),
                    ("main_directory", run_dir / "main"),
                    ("strict_no_ambiguous_directory", run_dir / "strict_no_ambiguous"),
                    ("branch_comparison", run_dir / "dog5_allelome_gwas_finemap_branch_comparison.tsv"),
                ]
            ),
        )

        log(
            "[stage 6/6] two GWAS branch all completed independent output; "
            "main is main results, strict_no_ambiguous for sensitivity analysis only"
        )
        log(f"task completed / done: {run_dir}")
        return 0
    except Exception as exc:
        log(f"ERROR: {exc}")
        traceback.print_exc()
        return 1
    finally:
        if LOGGER is not None:
            LOGGER.close()


if __name__ == "__main__":
    sys.exit(main())
