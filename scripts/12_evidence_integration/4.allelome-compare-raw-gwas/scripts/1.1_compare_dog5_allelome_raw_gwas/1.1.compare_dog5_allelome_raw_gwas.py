#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import json
import os
import signal
import traceback
from contextlib import ExitStack, redirect_stdout, redirect_stderr
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


SCRIPT_DIR = Path("/path/to/repository/scripts/12_evidence_integration/4.allelome-compare-raw-gwas/scripts/1.1_compare_dog5_allelome_raw_gwas")
INPUT_RUN = Path("/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429")
DEFAULT_ALLELOME_LINKS = INPUT_RUN / "dog5_allelome_links.annotated.tsv"
DEFAULT_ALLELOME_SUMMARY = INPUT_RUN / "dog5_allelome_annotation.summary.tsv"
GTF_RUN = Path("/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722")
DEFAULT_NEW_GENE_GTF = GTF_RUN / "five_only_new_gene_loci.new_only.gtf"
DEFAULT_NEW_GENE_BUILD_SUMMARY = GTF_RUN / "five_only_new_gene_loci.summary.tsv"
DEFAULT_GWAS_DIR = Path("/path/to/project/legacy/reference/gwas/summary_statistics")
DEFAULT_OUTPUT_ROOT = Path("/path/to/project/analysis_20260817/results/5-dogs/4.allelome-compare-raw-gwas/1.1_compare_dog5_allelome_raw_gwas")

GWAS_PATTERN = "score_*_GWAS_summary_statistics.tsv.gz"
GWAS_SUFFIX = "_GWAS_summary_statistics.tsv.gz"
EXPECTED_GWAS_FILE_COUNT = 14
GWAS_STRICT_THRESHOLD = 5e-8
GWAS_STRONG_THRESHOLD = 1e-5
GWAS_OLD_COMPARABLE_THRESHOLD = 5e-5
GWAS_CANDIDATE_THRESHOLD = 1e-4
NEW_GENE_BIN_SIZE = 1_000_000
AUTOSOMES = {f"chr{number}" for number in range(1, 39)}
VARIANT_RE = re.compile(
    r"(?:\(NA\))?(chr[^:|\s]+):(\d+)_([A-Za-z]+)>([A-Za-z]+)"
)
GTF_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')


TRACE_COLUMNS = [
    "canonical_pair_id", "task_id", "source_link_file", "source_link_line",
    "source_metadata_file", "overlap_bp", "strand_relation",
    "either_endpoint_multigene", "both_endpoints_unambiguous_count_positive",
    "identical_haplotype_counts", "evidence_stratum", "shared_variant_identity_checked",
]
LINK_REQUIRED_COLUMNS = {
    "animal_id", "tissue", "pair_type", "interval_gap_bp", "mechanism", "linkage_score",
    *[x for x in TRACE_COLUMNS if x != "canonical_pair_id"],
    *[f"{name}_{side}" for side in ("base", "target") for name in ("chr", "start", "end", "name")],
    *[f"{side}_{name}" for side in ("base", "target") for name in
      ("record_type", "gene_id", "gene_symbol", "gene_strand", "annotation_source",
       "gene_biotype", "endpoint_type")],
}

GWAS_REQUIRED_COLUMNS = {
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
}

GWAS_OUTPUT_COLUMNS = [
    "trait",
    "source_file",
    "gwas_marker_id",
    "gwas_chrom",
    "gwas_pos",
    "gwas_allele1",
    "gwas_allele0",
    "gwas_unordered_alleles",
    "gwas_af",
    "gwas_beta",
    "gwas_se",
    "gwas_p_wald",
    "gwas_p_lrt",
    "gwas_p_score",
    "gwas_or",
    "gwas_or_95ci_lower",
    "gwas_or_95ci_upper",
    "gwas_p_tier",
    "passes_old_p_lt_5e-5",
    "trait_variant_key",
    "gwas_manual_new_gene_hit",
    "gwas_manual_new_gene_hit_n",
    "gwas_manual_new_gene_ids",
    "gwas_manual_new_gene_strands",
    "gwas_manual_new_gene_class_codes",
    "gwas_manual_new_gene_novelty_categories",
    "gwas_manual_new_gene_candidate_categories",
    "gwas_manual_new_gene_dataset_support",
    "gwas_manual_new_gene_filter_reasons",
]

OVERLAP_COLUMNS = GWAS_OUTPUT_COLUMNS + [
    "sample_pair_id",
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
    "base_record_type",
    "base_gene_id",
    "base_gene_symbol",
    "base_gene_strand",
    "base_chr",
    "base_start",
    "base_end",
    "base_annotation_source",
    "base_gene_biotype",
    "base_endpoint_type",
    "base_manual_new_gene_hit",
    "base_manual_new_gene_hit_n",
    "base_manual_new_gene_ids",
    "base_manual_new_gene_strands",
    "base_manual_new_gene_class_codes",
    "base_manual_new_gene_novelty_categories",
    "base_manual_new_gene_candidate_categories",
    "base_manual_new_gene_dataset_support",
    "base_manual_new_gene_filter_reasons",
    "target_record_type",
    "target_gene_id",
    "target_gene_symbol",
    "target_gene_strand",
    "target_chr",
    "target_start",
    "target_end",
    "target_annotation_source",
    "target_gene_biotype",
    "target_endpoint_type",
    "target_manual_new_gene_hit",
    "target_manual_new_gene_hit_n",
    "target_manual_new_gene_ids",
    "target_manual_new_gene_strands",
    "target_manual_new_gene_class_codes",
    "target_manual_new_gene_novelty_categories",
    "target_manual_new_gene_candidate_categories",
    "target_manual_new_gene_dataset_support",
    "target_manual_new_gene_filter_reasons",
    "link_contains_manual_new_gene",
    "gwas_hit_endpoint_has_manual_new_gene",
    "gwas_ase_shared_manual_new_gene_ids",
    "new_gene_evidence_scope",
]

ENDPOINT_ANNOTATION_COLUMNS = [
    "endpoint_key",
    "animal_id",
    "tissue",
    "sample_tissue",
    "endpoint_roles",
    "chrom",
    "pos_1based",
    "ref",
    "alt",
    "unordered_alleles",
    "source_gene_id",
    "source_gene_symbol",
    "manual_new_gene_hit",
    "manual_new_gene_hit_n",
    "manual_new_gene_ids",
    "manual_new_gene_strands",
    "manual_new_gene_class_codes",
    "manual_new_gene_novelty_categories",
    "manual_new_gene_candidate_categories",
    "manual_new_gene_dataset_support",
    "manual_new_gene_filter_reasons",
    "annotation_rule",
]

NEW_GENE_SPAN_OVERLAP_COLUMNS = GWAS_OUTPUT_COLUMNS + [
    "sample_pair_id",
    "animal_id",
    "tissue",
    "sample_tissue",
    "link_chromosome",
    "pair_type",
    "interval_gap_bp",
    "mechanism",
    "linkage_score",
    "unannotated_endpoint_role",
    "unannotated_endpoint_key",
    "ase_chrom",
    "ase_pos_1based",
    "ase_ref",
    "ase_alt",
    "ase_unordered_alleles",
    "ase_manual_new_gene_ids",
    "shared_manual_new_gene_n",
    "shared_manual_new_gene_ids",
    "new_gene_gwas_match_class",
    "same_position",
    "unordered_allele_set_match",
    "ase_gwas_distance_bp",
    "partner_endpoint_role",
    "partner_record_type",
    "partner_gene_id",
    "partner_gene_symbol",
    "partner_gene_strand",
    "partner_chr",
    "partner_start",
    "partner_end",
    "partner_gene_biotype",
    "partner_endpoint_type",
    "partner_is_protein_coding_gene",
    "gwas_hits_partner_endpoint",
]

HIT_LINK_COLUMNS = [
    "sample_pair_id",
    "animal_id",
    "tissue",
    "sample_tissue",
    "chromosome",
    "pair_type",
    "mechanism",
    "linkage_score",
    "overlap_row_n",
    "unique_trait_variant_n",
    "traits",
    "endpoint_hits",
    "minimum_gwas_p_lrt",
    "has_old_p_lt_5e-5_hit",
    "link_contains_manual_new_gene",
    "direct_new_gene_endpoint_hit_n",
    "manual_new_gene_ids",
]

OVERLAP_COLUMNS += TRACE_COLUMNS
NEW_GENE_SPAN_OVERLAP_COLUMNS += TRACE_COLUMNS
HIT_LINK_COLUMNS += TRACE_COLUMNS


def log(message: str) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ] {message}"
    print(line, flush=True)


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def normalize_chrom(value: object) -> str:
    text = clean(value)
    if not text:
        return ""
    if text.startswith("chr"):
        return text
    if text.isdigit():
        return f"chr{text}"
    return text


def parse_int(value: object, label: str) -> int:
    text = clean(value)
    try:
        number = int(text)
    except ValueError as exc:
        raise ValueError(f"{label} is not an integer: {text!r}") from exc
    return number


def parse_float(value: object, label: str) -> float:
    text = clean(value)
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(f"{label} is not numeric: {text!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} is not finite numeric value: {text!r}")
    return number


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} does not exist or is empty: {path}")


def require_columns(
    path: Path,
    observed: Sequence[str],
    required: Set[str],
) -> None:
    missing = sorted(required - set(observed))
    if missing:
        raise ValueError(f"{path} missing required fields: {','.join(missing)}")


def read_metric_table(path: Path) -> Dict[str, str]:
    require_file(path, "summary table")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames or [], {"metric", "value"})
        return {
            clean(row["metric"]): clean(row["value"])
            for row in reader
            if clean(row.get("metric"))
        }


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


def allele_tuple(*values: object) -> Tuple[str, ...]:
    alleles = sorted({clean(value).upper() for value in values if clean(value)})
    return tuple(alleles)


def p_tier(p_value: float) -> str:
    if p_value <= GWAS_STRICT_THRESHOLD:
        return "strict_p_le_5e-8"
    if p_value <= GWAS_STRONG_THRESHOLD:
        return "strong_suggestive_p_le_1e-5"
    if p_value <= GWAS_OLD_COMPARABLE_THRESHOLD:
        return "loose_suggestive_p_le_5e-5"
    return "very_loose_p_le_1e-4"


def trait_from_path(path: Path) -> str:
    if not path.name.endswith(GWAS_SUFFIX):
        raise ValueError(f"cannot from filename parse GWAS trait: {path.name}")
    return path.name[: -len(GWAS_SUFFIX)]


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
    if not categories:
        categories.append("other_candidate_transcript_locus")
    return "|".join(categories)


def parse_gtf_attributes(text: str) -> Dict[str, str]:
    return {
        match.group(1): match.group(2)
        for match in GTF_ATTRIBUTE_RE.finditer(text)
    }


def load_new_gene_loci(gtf_path: Path, build_summary_path: Path):

    summary = read_metric_table(build_summary_path)
    if summary.get("source_dataset") != "five-dogs":
        raise ValueError("new gene build summary is not five-dogs")
    if Path(summary.get("new_only_gtf", "")).resolve() != gtf_path.resolve():
        raise ValueError("Build summary does not match the specified new-only GTF path")
    genes = {}
    with gtf_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF non-9 column: {gtf_path}:{line_number}")
            if fields[2] not in {"transcript", "exon"}:
                raise ValueError(f"non-transcript/exon record: {line_number}")
            if fields[2] != "transcript":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            if not gene_id.startswith("DOG5_NEWGENE_") or attrs.get("source_dataset") != "five-dogs":
                raise ValueError(f"is not five-dog novel gene: {gene_id}")
            chrom, strand = normalize_chrom(fields[0]), fields[6]
            start0, end0 = int(fields[3]) - 1, int(fields[4])
            if strand not in {"+", "-"} or start0 < 0 or end0 <= start0:
                raise ValueError(f"GTF interval or strand invalid: {gene_id}")
            row = genes.setdefault(gene_id, dict(gene_id=gene_id, chrom=chrom,
                strand=strand, start0=start0, end0=end0, dataset_support=attrs["source_dataset"],
                class_codes=set(), novelty=set(), filter_reasons=set()))
            if (row["chrom"], row["strand"]) != (chrom, strand):
                raise ValueError(f"the same gene across chromosome or across strand: {gene_id}")
            row["start0"] = min(row["start0"], start0)
            row["end0"] = max(row["end0"], end0)
            row["class_codes"].update(split_values(attrs["class_code"]))
            row["novelty"].update(split_values(attrs["novelty_category"]))
            row["filter_reasons"].update(split_values(attrs["gene_locus_filter_reason"]))
    bins = defaultdict(list)
    rows = []
    for gene_id in sorted(genes):
        row = genes[gene_id]
        row["class_codes"] = sorted(row["class_codes"])
        row["novelty"] = sorted(row["novelty"])
        row["filter_reason"] = "|".join(sorted(row.pop("filter_reasons")))
        row["candidate_category"] = candidate_category(row["class_codes"], row["novelty"])
        rows.append(row)
        if row["chrom"] in AUTOSOMES:
            for bin_id in range(row["start0"] // NEW_GENE_BIN_SIZE,
                                (row["end0"] - 1) // NEW_GENE_BIN_SIZE + 1):
                bins[(row["chrom"], bin_id)].append(row)
    return rows, bins


def find_new_gene_loci(
    bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    chrom: str,
    pos1: int,
) -> List[Dict[str, object]]:
    point0 = pos1 - 1
    candidates = bins.get((chrom, point0 // NEW_GENE_BIN_SIZE), [])
    hits = [
        row
        for row in candidates
        if int(row["start0"]) <= point0 < int(row["end0"])
    ]
    return sorted(hits, key=lambda row: clean(row["gene_id"]))


def summarize_new_gene_hits(hits: Sequence[Dict[str, object]]) -> Dict[str, str]:
    return {
        "manual_new_gene_hit": "yes" if hits else "no",
        "manual_new_gene_hit_n": str(len(hits)),
        "manual_new_gene_ids": join_values(
            row["gene_id"] for row in hits
        ),
        "manual_new_gene_strands": join_strands(row["strand"] for row in hits),
        "manual_new_gene_class_codes": join_values(
            code for row in hits for code in row["class_codes"]
        ),
        "manual_new_gene_novelty_categories": join_values(
            value for row in hits for value in row["novelty"]
        ),
        "manual_new_gene_candidate_categories": join_values(
            category
            for row in hits
            for category in split_values(row["candidate_category"])
        ),
        "manual_new_gene_dataset_support": join_values(
            row["dataset_support"] for row in hits
        ),
        "manual_new_gene_filter_reasons": join_values(
            row["filter_reason"] for row in hits
        ),
    }


def annotate_gwas_new_gene_hits(
    candidates: Sequence[Dict[str, object]],
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
) -> Tuple[Dict[str, List[Dict[str, object]]], int]:

    by_new_gene: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    annotated_candidate_count = 0
    for candidate in candidates:
        hits = find_new_gene_loci(
            new_gene_bins,
            str(candidate["gwas_chrom"]),
            int(candidate["gwas_pos"]),
        )
        annotation = summarize_new_gene_hits(hits)
        candidate.update(prefixed_annotation("gwas", annotation))
        if hits:
            annotated_candidate_count += 1
        for hit in hits:
            by_new_gene[str(hit["gene_id"])].append(candidate)
    return dict(by_new_gene), annotated_candidate_count


def discover_gwas_files(gwas_dir: Path) -> List[Path]:
    files = sorted(gwas_dir.glob(GWAS_PATTERN))
    if len(files) != EXPECTED_GWAS_FILE_COUNT:
        raise ValueError(
            f"GWAS file count unexpected: observed={len(files)}, "
            f"expected={EXPECTED_GWAS_FILE_COUNT}, pattern={GWAS_PATTERN}"
        )
    return files


def load_gwas_candidates(
    gwas_dir: Path,
) -> Tuple[
    List[Dict[str, object]],
    Dict[str, List[Dict[str, object]]],
    Dict[str, List[int]],
    Dict[Tuple[str, int], List[Dict[str, object]]],
]:
    candidates: List[Dict[str, object]] = []
    seen_keys: Set[str] = set()
    files = discover_gwas_files(gwas_dir)

    for file_index, path in enumerate(files, start=1):
        trait = trait_from_path(path)
        log(
            f"[stage 2/5] scan GWAS {file_index}/{len(files)}: {path.name}"
        )
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(path, reader.fieldnames or [], GWAS_REQUIRED_COLUMNS)
            for line_number, row in enumerate(reader, start=2):
                p_text = clean(row.get("p_lrt"))
                if not p_text:
                    continue
                p_value = parse_float(p_text, f"{path}:{line_number} p_lrt")
                if p_value > GWAS_CANDIDATE_THRESHOLD:
                    continue
                chrom = normalize_chrom(row.get("chr"))
                pos = parse_int(row.get("ps"), f"{path}:{line_number} ps")
                alleles = allele_tuple(row.get("allele1"), row.get("allele0"))
                if chrom not in AUTOSOMES or pos <= 0 or len(alleles) != 2:
                    raise ValueError(
                        f"{path}:{line_number} GWAS candidate coordinates or allele invalid: "
                        f"chrom={chrom}, pos={pos}, alleles={alleles}"
                    )
                trait_variant_key = (
                    f"{trait}|{chrom}:{pos}:{'/'.join(alleles)}"
                )
                if trait_variant_key in seen_keys:
                    raise ValueError(
                        f"duplicate GWAS trait-variant candidate: {trait_variant_key}"
                    )
                seen_keys.add(trait_variant_key)
                candidates.append(
                    {
                        "trait": trait,
                        "source_file": path.name,
                        "gwas_marker_id": clean(row.get("rs")),
                        "gwas_chrom": chrom,
                        "gwas_pos": pos,
                        "gwas_allele1": clean(row.get("allele1")).upper(),
                        "gwas_allele0": clean(row.get("allele0")).upper(),
                        "gwas_unordered_alleles": "/".join(alleles),
                        "gwas_af": clean(row.get("af")),
                        "gwas_beta": clean(row.get("beta")),
                        "gwas_se": clean(row.get("se")),
                        "gwas_p_wald": clean(row.get("p_wald")),
                        "gwas_p_lrt": p_value,
                        "gwas_p_score": clean(row.get("p_score")),
                        "gwas_or": clean(row.get("OR")),
                        "gwas_or_95ci_lower": clean(row.get("OR_95CI_LOWER")),
                        "gwas_or_95ci_upper": clean(row.get("OR_95CI_UPPER")),
                        "gwas_p_tier": p_tier(p_value),
                        "passes_old_p_lt_5e-5": (
                            "yes"
                            if p_value < GWAS_OLD_COMPARABLE_THRESHOLD
                            else "no"
                        ),
                        "trait_variant_key": trait_variant_key,
                        "_alleles": alleles,
                    }
                )

    by_chrom: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    by_position: Dict[Tuple[str, int], List[Dict[str, object]]] = defaultdict(list)
    for row in candidates:
        chrom = str(row["gwas_chrom"])
        pos = int(row["gwas_pos"])
        by_chrom[chrom].append(row)
        by_position[(chrom, pos)].append(row)
    positions: Dict[str, List[int]] = {}
    for chrom, rows in by_chrom.items():
        rows.sort(
            key=lambda row: (
                int(row["gwas_pos"]),
                str(row["trait"]),
                str(row["gwas_marker_id"]),
            )
        )
        positions[chrom] = [int(row["gwas_pos"]) for row in rows]
    return candidates, dict(by_chrom), positions, dict(by_position)


def parse_unannotated_variant(endpoint: Mapping[str, object]) -> Dict[str, object]:
    texts = [
        clean(endpoint.get("gene_id")),
        clean(endpoint.get("gene_symbol")),
        clean(endpoint.get("name")),
    ]
    match = None
    for text in texts:
        match = VARIANT_RE.search(text)
        if match:
            break
    if match is None:
        raise ValueError(
            f"cannot parse the unannotated endpoint CHROM:POS_REF>ALT: {texts}"
        )
    chrom = normalize_chrom(match.group(1))
    pos = int(match.group(2))
    ref = match.group(3).upper()
    alt = match.group(4).upper()
    endpoint_chr = normalize_chrom(endpoint.get("chrom"))
    start = parse_int(endpoint.get("start"), "unannotated endpoint start")
    end = parse_int(endpoint.get("end"), "unannotated endpoint end")
    if endpoint_chr != chrom or start != pos - 1 or end != pos:
        raise ValueError(
            f"unannotated endpoint fields differ: parsed={chrom}:{pos}_{ref}>{alt}, "
            f"columns={endpoint_chr}:{start}-{end}"
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
        "role": prefix,
        "record_type": clean(row[f"{prefix}_record_type"]),
        "gene_id": clean(row[f"{prefix}_gene_id"]),
        "gene_symbol": clean(row[f"{prefix}_gene_symbol"]),
        "gene_strand": clean(row[f"{prefix}_gene_strand"]),
        "chrom": normalize_chrom(row[f"chr_{prefix}"]),
        "start": clean(row[f"start_{prefix}"]),
        "end": clean(row[f"end_{prefix}"]),
        "name": clean(row[f"name_{prefix}"]),
        "annotation_source": clean(row[f"{prefix}_annotation_source"]),
        "gene_biotype": clean(row[f"{prefix}_gene_biotype"]),
        "endpoint_type": clean(row[f"{prefix}_endpoint_type"]),
    }


def annotate_endpoint(
    endpoint: Mapping[str, object],
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
) -> Tuple[
    Dict[str, str],
    Dict[str, object] | None,
    List[Dict[str, object]],
]:
    if endpoint["record_type"] == "manual_new_gene":

        hits = [r for r in find_new_gene_loci(new_gene_bins, str(endpoint["chrom"]),
                                             int(endpoint["start"]) + 1)
                if r["gene_id"] == endpoint["gene_id"]]
        if len(hits) != 1 or (hits[0]["start0"], hits[0]["end0"], hits[0]["strand"]) != (
                int(endpoint["start"]), int(endpoint["end"]), endpoint["gene_strand"]):
            raise ValueError(f"New-gene endpoint does not match the current GTF: {endpoint['gene_id']}")
        return summarize_new_gene_hits(hits), None, hits
    if endpoint["record_type"] == "known_gene":
        return summarize_new_gene_hits([]), None, []
    if endpoint["record_type"] != "unannotated_variant":
        raise ValueError(f"unsupported endpoint record_type: {endpoint['record_type']}")
    variant = parse_unannotated_variant(endpoint)
    hits = find_new_gene_loci(
        new_gene_bins,
        str(variant["chrom"]),
        int(variant["pos"]),
    )
    return summarize_new_gene_hits(hits), variant, hits


def shared_new_gene_gwas_candidates(
    endpoint_new_gene_hits: Sequence[Dict[str, object]],
    gwas_by_new_gene: Mapping[str, Sequence[Dict[str, object]]],
) -> List[Tuple[Dict[str, object], List[str]]]:

    matches: Dict[str, Dict[str, object]] = {}
    for locus in endpoint_new_gene_hits:
        locus_id = str(locus["gene_id"])
        for candidate in gwas_by_new_gene.get(locus_id, []):
            key = str(candidate["trait_variant_key"])
            entry = matches.setdefault(
                key,
                {"candidate": candidate, "shared_ids": set()},
            )
            entry["shared_ids"].add(locus_id)
    return [
        (
            matches[key]["candidate"],
            sorted(matches[key]["shared_ids"]),
        )
        for key in sorted(matches)
    ]


def classify_new_gene_gwas_match(
    variant: Mapping[str, object],
    candidate: Mapping[str, object],
) -> Tuple[str, str, str, int]:

    same_position = (
        str(variant["chrom"]) == str(candidate["gwas_chrom"])
        and int(variant["pos"]) == int(candidate["gwas_pos"])
    )
    allele_match = (
        tuple(variant["alleles"]) == tuple(candidate["_alleles"])
        if same_position
        else None
    )
    if same_position and allele_match:
        match_class = "exact_unannotated_ase_variant"
    elif same_position:
        match_class = "same_position_allele_mismatch"
    else:
        match_class = "same_manual_new_gene_span_nonexact"
    return (
        match_class,
        "yes" if same_position else "no",
        "yes" if allele_match else "no" if same_position else "not_applicable",
        abs(int(candidate["gwas_pos"]) - int(variant["pos"])),
    )


def build_new_gene_span_overlap_row(
    link_row: Mapping[str, str],
    endpoint_role: str,
    variant: Mapping[str, object],
    endpoint_annotation: Mapping[str, str],
    partner_role: str,
    partner: Mapping[str, object],
    partner_hit_keys: Set[str],
    candidate: Mapping[str, object],
    shared_new_gene_ids: Sequence[str],
) -> Dict[str, object]:
    (
        match_class,
        same_position,
        allele_set_match,
        distance_bp,
    ) = classify_new_gene_gwas_match(variant, candidate)
    sample_tissue = clean(link_row["sample_tissue"])
    return {
        **output_gwas_row(candidate),
        **{k: link_row[k] for k in TRACE_COLUMNS},
        "sample_pair_id": clean(link_row["sample_pair_id"]),
        "animal_id": clean(link_row["animal_id"]),
        "tissue": clean(link_row["tissue"]),
        "sample_tissue": sample_tissue,
        "link_chromosome": clean(link_row["chromosome"]),
        "pair_type": clean(link_row["pair_type"]),
        "interval_gap_bp": clean(link_row["interval_gap_bp"]),
        "mechanism": clean(link_row["mechanism"]),
        "linkage_score": clean(link_row["linkage_score"]),
        "unannotated_endpoint_role": endpoint_role,
        "unannotated_endpoint_key": (
            f"{sample_tissue}|{variant['chrom']}:{variant['pos']}_"
            f"{variant['ref']}>{variant['alt']}"
        ),
        "ase_chrom": variant["chrom"],
        "ase_pos_1based": variant["pos"],
        "ase_ref": variant["ref"],
        "ase_alt": variant["alt"],
        "ase_unordered_alleles": "/".join(variant["alleles"]),
        "ase_manual_new_gene_ids": endpoint_annotation[
            "manual_new_gene_ids"
        ],
        "shared_manual_new_gene_n": len(shared_new_gene_ids),
        "shared_manual_new_gene_ids": "|".join(shared_new_gene_ids),
        "new_gene_gwas_match_class": match_class,
        "same_position": same_position,
        "unordered_allele_set_match": allele_set_match,
        "ase_gwas_distance_bp": distance_bp,
        "partner_endpoint_role": partner_role,
        "partner_record_type": clean(partner["record_type"]),
        "partner_gene_id": clean(partner["gene_id"]),
        "partner_gene_symbol": clean(partner["gene_symbol"]),
        "partner_gene_strand": clean(partner["gene_strand"]),
        "partner_chr": clean(partner["chrom"]),
        "partner_start": clean(partner["start"]),
        "partner_end": clean(partner["end"]),
        "partner_gene_biotype": clean(partner["gene_biotype"]),
        "partner_endpoint_type": clean(partner["endpoint_type"]),
        "partner_is_protein_coding_gene": (
            "yes"
            if clean(partner["record_type"]) == "known_gene"
            and (
                clean(partner["gene_biotype"]) == "protein_coding"
                or clean(partner["endpoint_type"]) == "known_protein_coding"
            )
            else "no"
        ),
        "gwas_hits_partner_endpoint": (
            "yes"
            if str(candidate["trait_variant_key"]) in partner_hit_keys
            else "no"
        ),
    }


def match_endpoint_to_gwas(
    endpoint: Mapping[str, object],
    gwas_by_chrom: Mapping[str, Sequence[Dict[str, object]]],
    gwas_positions: Mapping[str, Sequence[int]],
    gwas_by_position: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
) -> Tuple[List[Dict[str, object]], str]:
    record_type = str(endpoint["record_type"])
    chrom = str(endpoint["chrom"])
    if record_type == "unannotated_variant":
        variant = parse_unannotated_variant(endpoint)
        candidates = gwas_by_position.get(
            (str(variant["chrom"]), int(variant["pos"])),
            [],
        )
        hits = [
            row
            for row in candidates
            if tuple(row["_alleles"]) == tuple(variant["alleles"])
        ]
        return hits, "point_exact_chrom_pos_unordered_allele_set"
    if record_type in {"known_gene", "manual_new_gene"}:
        start0 = parse_int(endpoint["start"], "annotated endpoint START0")
        end0 = parse_int(endpoint["end"], "annotated endpoint END0")
        if chrom not in AUTOSOMES or start0 < 0 or end0 <= start0:
            raise ValueError(
                f"annotated gene endpoint interval invalid: {chrom}:{start0}-{end0}"
            )
        rows = gwas_by_chrom.get(chrom, [])
        positions = gwas_positions.get(chrom, [])
        left = bisect.bisect_right(positions, start0)
        right = bisect.bisect_right(positions, end0)
        return (
            list(rows[left:right]),
            "gene_span_start0_end0_contains_1based_gwas_pos",
        )
    raise ValueError(f"unsupported endpoint record_type: {record_type}")


def prefixed_annotation(prefix: str, annotation: Mapping[str, str]) -> Dict[str, str]:
    return {f"{prefix}_{key}": value for key, value in annotation.items()}


def endpoint_output(prefix: str, endpoint: Mapping[str, object]) -> Dict[str, str]:
    return {
        f"{prefix}_record_type": clean(endpoint["record_type"]),
        f"{prefix}_gene_id": clean(endpoint["gene_id"]),
        f"{prefix}_gene_symbol": clean(endpoint["gene_symbol"]),
        f"{prefix}_gene_strand": clean(endpoint["gene_strand"]),
        f"{prefix}_chr": clean(endpoint["chrom"]),
        f"{prefix}_start": clean(endpoint["start"]),
        f"{prefix}_end": clean(endpoint["end"]),
        f"{prefix}_annotation_source": clean(endpoint["annotation_source"]),
        f"{prefix}_gene_biotype": clean(endpoint["gene_biotype"]),
        f"{prefix}_endpoint_type": clean(endpoint["endpoint_type"]),
    }


def write_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def open_gzip_writer(
    path: Path,
    fieldnames: Sequence[str],
):
    handle = gzip.open(path, "xt", encoding="utf-8", newline="")
    writer = csv.DictWriter(
        handle,
        fieldnames=list(fieldnames),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    return handle, writer


def output_gwas_row(row: Mapping[str, object]) -> Dict[str, object]:
    return {column: row.get(column, "") for column in GWAS_OUTPUT_COLUMNS}


def update_endpoint_annotation_inventory(
    inventory: MutableMapping[Tuple[str, str, int, str, str], Dict[str, object]],
    row: Mapping[str, str],
    endpoint: Mapping[str, object],
    role: str,
    variant: Mapping[str, object] | None,
    annotation: Mapping[str, str],
) -> None:
    if variant is None:
        return
    key = (
        clean(row["sample_tissue"]),
        str(variant["chrom"]),
        int(variant["pos"]),
        str(variant["ref"]),
        str(variant["alt"]),
    )
    if key not in inventory:
        inventory[key] = {
            "endpoint_key": (
                f"{key[0]}|{key[1]}:{key[2]}_{key[3]}>{key[4]}"
            ),
            "animal_id": clean(row["animal_id"]),
            "tissue": clean(row["tissue"]),
            "sample_tissue": clean(row["sample_tissue"]),
            "_roles": set(),
            "chrom": key[1],
            "pos_1based": key[2],
            "ref": key[3],
            "alt": key[4],
            "unordered_alleles": "/".join(allele_tuple(key[3], key[4])),
            "source_gene_id": clean(endpoint["gene_id"]),
            "source_gene_symbol": clean(endpoint["gene_symbol"]),
            **annotation,
            "annotation_rule": (
                "unannotated_variant_pos1_to_manual_gtf_transcript_union_gene_span"
            ),
        }
    inventory[key]["_roles"].add(role)


def flatten_endpoint_inventory(
    inventory: Mapping[Tuple[str, str, int, str, str], Dict[str, object]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for key in sorted(inventory):
        row = dict(inventory[key])
        row["endpoint_roles"] = "|".join(sorted(row.pop("_roles")))
        rows.append(row)
    return rows


def update_hit_link_summary(
    summaries: MutableMapping[str, Dict[str, object]],
    overlap: Mapping[str, object],
) -> None:
    pair_id = str(overlap["sample_pair_id"])
    if pair_id not in summaries:
        summaries[pair_id] = {
            "sample_pair_id": pair_id,
            **{k: overlap[k] for k in TRACE_COLUMNS},
            "animal_id": overlap["animal_id"],
            "tissue": overlap["tissue"],
            "sample_tissue": overlap["sample_tissue"],
            "chromosome": overlap["link_chromosome"],
            "pair_type": overlap["pair_type"],
            "mechanism": overlap["mechanism"],
            "linkage_score": overlap["linkage_score"],
            "overlap_row_n": 0,
            "_trait_variants": set(),
            "_traits": set(),
            "_endpoint_hits": set(),
            "_minimum_p": float("inf"),
            "_old_hit": False,
            "_contains_new_gene": False,
            "_direct_new_gene_n": 0,
            "_new_gene_ids": set(),
        }
    summary = summaries[pair_id]
    summary["overlap_row_n"] = int(summary["overlap_row_n"]) + 1
    summary["_trait_variants"].add(overlap["trait_variant_key"])
    summary["_traits"].add(overlap["trait"])
    summary["_endpoint_hits"].add(overlap["endpoint_hit"])
    summary["_minimum_p"] = min(
        float(summary["_minimum_p"]),
        float(overlap["gwas_p_lrt"]),
    )
    if overlap["passes_old_p_lt_5e-5"] == "yes":
        summary["_old_hit"] = True
    if overlap["link_contains_manual_new_gene"] == "yes":
        summary["_contains_new_gene"] = True
    if overlap["gwas_hit_endpoint_has_manual_new_gene"] == "yes":
        summary["_direct_new_gene_n"] = int(summary["_direct_new_gene_n"]) + 1
    for value in (
        clean(overlap["base_manual_new_gene_ids"]),
        clean(overlap["target_manual_new_gene_ids"]),
    ):
        summary["_new_gene_ids"].update(split_values(value))


def flatten_hit_link_summaries(
    summaries: Mapping[str, Dict[str, object]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for pair_id in sorted(summaries):
        source = summaries[pair_id]
        rows.append(
            {
                **{
                    key: source[key]
                    for key in [
                        "sample_pair_id",
                        "animal_id",
                        "tissue",
                        "sample_tissue",
                        "chromosome",
                        "pair_type",
                        "mechanism",
                        "linkage_score",
                        "overlap_row_n",
                    ]
                },
                **{k: source[k] for k in TRACE_COLUMNS},
                "unique_trait_variant_n": len(source["_trait_variants"]),
                "traits": "|".join(sorted(source["_traits"])),
                "endpoint_hits": "|".join(sorted(source["_endpoint_hits"])),
                "minimum_gwas_p_lrt": source["_minimum_p"],
                "has_old_p_lt_5e-5_hit": (
                    "yes" if source["_old_hit"] else "no"
                ),
                "link_contains_manual_new_gene": (
                    "yes" if source["_contains_new_gene"] else "no"
                ),
                "direct_new_gene_endpoint_hit_n": source["_direct_new_gene_n"],
                "manual_new_gene_ids": "|".join(
                    sorted(source["_new_gene_ids"])
                ),
            }
        )
    return rows


def process_links(
    links_path: Path,
    expected_link_rows: int,
    new_gene_bins: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    gwas_candidates: Sequence[Dict[str, object]],
    gwas_by_chrom: Mapping[str, Sequence[Dict[str, object]]],
    gwas_positions: Mapping[str, Sequence[int]],
    gwas_by_position: Mapping[Tuple[str, int], Sequence[Dict[str, object]]],
    gwas_by_new_gene: Mapping[str, Sequence[Dict[str, object]]],
    run_dir: Path,
) -> Dict[str, object]:
    all_path = run_dir / "dog5_allelome_raw_gwas_overlap_all_p_le_1e-4.tsv.gz"
    old_path = run_dir / "dog5_allelome_raw_gwas_overlap_old_p_lt_5e-5.tsv.gz"
    new_gene_path = (
        run_dir / "dog5_allelome_raw_gwas_overlap_new_gene_context.tsv.gz"
    )
    span_all_path = (
        run_dir
        / "dog5_unannotated_ase_gwas_new_gene_span_overlap_all_p_le_1e-4.tsv.gz"
    )
    span_strong_path = (
        run_dir
        / "dog5_unannotated_ase_gwas_new_gene_span_overlap_strong_p_le_1e-5.tsv.gz"
    )
    span_strong_pcgene_path = (
        run_dir
        / "dog5_unannotated_ase_to_pcgene_gwas_new_gene_span_strong_"
        "p_le_1e-5.tsv.gz"
    )

    all_handle, all_writer = open_gzip_writer(all_path, OVERLAP_COLUMNS)
    old_handle, old_writer = open_gzip_writer(old_path, OVERLAP_COLUMNS)
    new_handle, new_writer = open_gzip_writer(new_gene_path, OVERLAP_COLUMNS)
    span_all_handle, span_all_writer = open_gzip_writer(
        span_all_path, NEW_GENE_SPAN_OVERLAP_COLUMNS
    )
    span_strong_handle, span_strong_writer = open_gzip_writer(
        span_strong_path, NEW_GENE_SPAN_OVERLAP_COLUMNS
    )
    span_strong_pcgene_handle, span_strong_pcgene_writer = open_gzip_writer(
        span_strong_pcgene_path, NEW_GENE_SPAN_OVERLAP_COLUMNS
    )

    link_count = 0
    overlap_count = 0
    old_overlap_count = 0
    new_gene_context_count = 0
    direct_new_gene_count = 0
    seen_pair_ids: Set[str] = set()
    hit_marker_keys: Set[str] = set()
    hit_link_summaries: Dict[str, Dict[str, object]] = {}
    endpoint_inventory: Dict[
        Tuple[str, str, int, str, str],
        Dict[str, object],
    ] = {}
    overlap_tier_counts: Counter[str] = Counter()
    span_overlap_count = 0
    span_strong_overlap_count = 0
    span_exact_count = 0
    span_nonexact_count = 0
    span_allele_mismatch_count = 0
    span_hit_link_ids: Set[str] = set()
    span_strong_hit_link_ids: Set[str] = set()
    span_hit_trait_variants: Set[str] = set()
    span_strong_hit_trait_variants: Set[str] = set()
    span_strong_match_class_counts: Counter[str] = Counter()
    span_strong_pcgene_count = 0
    span_strong_pcgene_hit_link_ids: Set[str] = set()
    span_strong_pcgene_hit_trait_variants: Set[str] = set()

    try:
        with links_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(links_path, reader.fieldnames or [], LINK_REQUIRED_COLUMNS)
            for line_number, row in enumerate(reader, start=2):
                link_count += 1
                endpoint_ids = sorted((row["base_gene_id"], row["target_gene_id"]))
                row["canonical_pair_id"] = json.dumps(
                    [row["chr_base"], *endpoint_ids], separators=(",", ":"))
                row["sample_pair_id"] = json.dumps(
                    [row["animal_id"], row["tissue"], row["chr_base"], *endpoint_ids],
                    separators=(",", ":"))
                row["sample_tissue"] = row["animal_id"] + "|" + row["tissue"]
                row["chromosome"] = row["chr_base"]
                if row["chr_base"] != row["chr_target"]:
                    raise ValueError("Allelome endpoint across chromosome")
                pair_id = clean(row["sample_pair_id"])
                if not pair_id or pair_id in seen_pair_ids:
                    raise ValueError(
                        f"{links_path}:{line_number} sample_pair_id "
                        f"missing or duplicate: {pair_id!r}"
                    )
                seen_pair_ids.add(pair_id)
                if clean(row["mechanism"]) not in {"enhancing", "repressing"}:
                    raise ValueError(
                        f"{links_path}:{line_number} mechanism invalid: "
                        f"{row['mechanism']!r}"
                    )

                base = endpoint_from_row(row, "base")
                target = endpoint_from_row(row, "target")
                (
                    base_annotation,
                    base_variant,
                    base_new_gene_hits,
                ) = annotate_endpoint(base, new_gene_bins)
                (
                    target_annotation,
                    target_variant,
                    target_new_gene_hits,
                ) = annotate_endpoint(target, new_gene_bins)
                update_endpoint_annotation_inventory(
                    endpoint_inventory,
                    row,
                    base,
                    "base",
                    base_variant,
                    base_annotation,
                )
                update_endpoint_annotation_inventory(
                    endpoint_inventory,
                    row,
                    target,
                    "target",
                    target_variant,
                    target_annotation,
                )

                base_hits, base_rule = match_endpoint_to_gwas(
                    base,
                    gwas_by_chrom,
                    gwas_positions,
                    gwas_by_position,
                )
                target_hits, target_rule = match_endpoint_to_gwas(
                    target,
                    gwas_by_chrom,
                    gwas_positions,
                    gwas_by_position,
                )
                base_hit_keys = {
                    str(candidate["trait_variant_key"])
                    for candidate in base_hits
                }
                target_hit_keys = {
                    str(candidate["trait_variant_key"])
                    for candidate in target_hits
                }

                endpoint_span_specs = [
                    (
                        "base",
                        base_variant,
                        base_annotation,
                        base_new_gene_hits,
                        "target",
                        target,
                        target_hit_keys,
                    ),
                    (
                        "target",
                        target_variant,
                        target_annotation,
                        target_new_gene_hits,
                        "base",
                        base,
                        base_hit_keys,
                    ),
                ]
                for (
                    endpoint_role,
                    variant,
                    endpoint_annotation,
                    endpoint_new_gene_hits,
                    partner_role,
                    partner,
                    partner_hit_keys,
                ) in endpoint_span_specs:
                    if variant is None or not endpoint_new_gene_hits:
                        continue
                    for candidate, shared_ids in shared_new_gene_gwas_candidates(
                        endpoint_new_gene_hits,
                        gwas_by_new_gene,
                    ):
                        span_overlap = build_new_gene_span_overlap_row(
                            row,
                            endpoint_role,
                            variant,
                            endpoint_annotation,
                            partner_role,
                            partner,
                            partner_hit_keys,
                            candidate,
                            shared_ids,
                        )
                        span_all_writer.writerow(
                            {
                                column: span_overlap.get(column, "")
                                for column in NEW_GENE_SPAN_OVERLAP_COLUMNS
                            }
                        )
                        span_overlap_count += 1
                        pair_id_for_span = clean(row["sample_pair_id"])
                        candidate_key_for_span = str(
                            candidate["trait_variant_key"]
                        )
                        span_hit_link_ids.add(pair_id_for_span)
                        span_hit_trait_variants.add(candidate_key_for_span)
                        match_class = str(
                            span_overlap["new_gene_gwas_match_class"]
                        )
                        if match_class == "exact_unannotated_ase_variant":
                            span_exact_count += 1
                        elif match_class == "same_position_allele_mismatch":
                            span_allele_mismatch_count += 1
                        else:
                            span_nonexact_count += 1
                        if float(candidate["gwas_p_lrt"]) <= GWAS_STRONG_THRESHOLD:
                            span_strong_writer.writerow(
                                {
                                    column: span_overlap.get(column, "")
                                    for column in NEW_GENE_SPAN_OVERLAP_COLUMNS
                                }
                            )
                            span_strong_overlap_count += 1
                            span_strong_hit_link_ids.add(pair_id_for_span)
                            span_strong_hit_trait_variants.add(
                                candidate_key_for_span
                            )
                            span_strong_match_class_counts[match_class] += 1
                            if (
                                span_overlap[
                                    "partner_is_protein_coding_gene"
                                ]
                                == "yes"
                            ):
                                span_strong_pcgene_writer.writerow(
                                    {
                                        column: span_overlap.get(column, "")
                                        for column in NEW_GENE_SPAN_OVERLAP_COLUMNS
                                    }
                                )
                                span_strong_pcgene_count += 1
                                span_strong_pcgene_hit_link_ids.add(
                                    pair_id_for_span
                                )
                                span_strong_pcgene_hit_trait_variants.add(
                                    candidate_key_for_span
                                )

                hits: Dict[str, Dict[str, object]] = {}
                for candidate in base_hits:
                    key = str(candidate["trait_variant_key"])
                    hits.setdefault(
                        key,
                        {"candidate": candidate, "base": "", "target": ""},
                    )["base"] = base_rule
                for candidate in target_hits:
                    key = str(candidate["trait_variant_key"])
                    hits.setdefault(
                        key,
                        {"candidate": candidate, "base": "", "target": ""},
                    )["target"] = target_rule

                for key in sorted(hits):
                    match = hits[key]
                    candidate = match["candidate"]
                    base_hit = bool(match["base"])
                    target_hit = bool(match["target"])
                    endpoint_hit = (
                        "both"
                        if base_hit and target_hit
                        else "base"
                        if base_hit
                        else "target"
                    )
                    base_new = base_annotation["manual_new_gene_hit"] == "yes"
                    target_new = (
                        target_annotation["manual_new_gene_hit"] == "yes"
                    )
                    link_contains_new = base_new or target_new
                    hit_endpoint_has_new = (
                        (base_hit and base_new) or (target_hit and target_new)
                    )
                    endpoint_new_gene_ids = set(
                        split_values(
                            base_annotation["manual_new_gene_ids"],
                            target_annotation["manual_new_gene_ids"],
                        )
                    )
                    candidate_new_gene_ids = set(
                        split_values(candidate["gwas_manual_new_gene_ids"])
                    )
                    shared_link_new_gene_ids = sorted(
                        endpoint_new_gene_ids & candidate_new_gene_ids
                    )
                    if hit_endpoint_has_new:
                        evidence_scope = (
                            "gwas_hit_candidate_new_gene_span"
                            if ((base_hit and base["record_type"] == "manual_new_gene")
                                or (target_hit and target["record_type"] == "manual_new_gene"))
                            else "exact_gwas_hit_unannotated_variant_with_new_gene_annotation"
                        )
                    elif shared_link_new_gene_ids:
                        evidence_scope = (
                            "gwas_and_unannotated_ase_share_manual_new_gene_span"
                        )
                    elif link_contains_new:
                        evidence_scope = (
                            "link_context_only_gwas_hits_other_endpoint"
                        )
                    else:
                        evidence_scope = "no_manual_new_gene_context"

                    overlap: Dict[str, object] = {
                        **output_gwas_row(candidate),
                        **{k: row[k] for k in TRACE_COLUMNS},
                        "sample_pair_id": pair_id,
                        "animal_id": clean(row["animal_id"]),
                        "tissue": clean(row["tissue"]),
                        "sample_tissue": clean(row["sample_tissue"]),
                        "link_chromosome": clean(row["chromosome"]),
                        "pair_type": clean(
                            row["pair_type"]
                        ),
                        "interval_gap_bp": clean(
                            row["interval_gap_bp"]
                        ),
                        "mechanism": clean(row["mechanism"]),
                        "linkage_score": clean(row["linkage_score"]),
                        "endpoint_hit": endpoint_hit,
                        "base_mapping_rule": match["base"],
                        "target_mapping_rule": match["target"],
                        **endpoint_output("base", base),
                        **prefixed_annotation("base", base_annotation),
                        **endpoint_output("target", target),
                        **prefixed_annotation("target", target_annotation),
                        "link_contains_manual_new_gene": (
                            "yes" if link_contains_new else "no"
                        ),
                        "gwas_hit_endpoint_has_manual_new_gene": (
                            "yes" if hit_endpoint_has_new else "no"
                        ),
                        "gwas_ase_shared_manual_new_gene_ids": "|".join(
                            shared_link_new_gene_ids
                        ),
                        "new_gene_evidence_scope": evidence_scope,
                    }
                    all_writer.writerow(
                        {column: overlap.get(column, "") for column in OVERLAP_COLUMNS}
                    )
                    overlap_count += 1
                    overlap_tier_counts[str(candidate["gwas_p_tier"])] += 1
                    hit_marker_keys.add(str(candidate["trait_variant_key"]))
                    update_hit_link_summary(hit_link_summaries, overlap)

                    if candidate["passes_old_p_lt_5e-5"] == "yes":
                        old_writer.writerow(
                            {
                                column: overlap.get(column, "")
                                for column in OVERLAP_COLUMNS
                            }
                        )
                        old_overlap_count += 1
                    if link_contains_new:
                        new_writer.writerow(
                            {
                                column: overlap.get(column, "")
                                for column in OVERLAP_COLUMNS
                            }
                        )
                        new_gene_context_count += 1
                    if hit_endpoint_has_new:
                        direct_new_gene_count += 1

                if link_count % 50_000 == 0:
                    log(
                        f"[stage 3/5] links={link_count}/{expected_link_rows}, "
                        f"overlap_rows={overlap_count}, "
                        f"same_new_gene_rows={span_overlap_count}"
                    )
    finally:
        all_handle.close()
        old_handle.close()
        new_handle.close()
        span_all_handle.close()
        span_strong_handle.close()
        span_strong_pcgene_handle.close()

    endpoint_rows = flatten_endpoint_inventory(endpoint_inventory)
    endpoint_path = (
        run_dir / "dog5_unannotated_endpoint_new_gene_annotation.tsv.gz"
    )
    endpoint_handle, endpoint_writer = open_gzip_writer(
        endpoint_path, ENDPOINT_ANNOTATION_COLUMNS
    )
    try:
        for row in endpoint_rows:
            endpoint_writer.writerow(
                {
                    column: row.get(column, "")
                    for column in ENDPOINT_ANNOTATION_COLUMNS
                }
            )
    finally:
        endpoint_handle.close()

    hit_link_rows = flatten_hit_link_summaries(hit_link_summaries)
    hit_link_path = run_dir / "dog5_allelome_raw_gwas_hit_link_summary.tsv.gz"
    hit_link_handle, hit_link_writer = open_gzip_writer(
        hit_link_path, HIT_LINK_COLUMNS
    )
    try:
        for row in hit_link_rows:
            hit_link_writer.writerow(
                {column: row.get(column, "") for column in HIT_LINK_COLUMNS}
            )
    finally:
        hit_link_handle.close()

    unmatched_rows = [
        output_gwas_row(row)
        for row in gwas_candidates
        if row["trait_variant_key"] not in hit_marker_keys
    ]
    unmatched_path = (
        run_dir / "dog5_raw_gwas_candidate_markers_without_allelome_hit.tsv.gz"
    )
    unmatched_handle, unmatched_writer = open_gzip_writer(
        unmatched_path, GWAS_OUTPUT_COLUMNS
    )
    try:
        for row in unmatched_rows:
            unmatched_writer.writerow(row)
    finally:
        unmatched_handle.close()

    endpoint_new_gene_count = sum(
        row["manual_new_gene_hit"] == "yes" for row in endpoint_rows
    )
    return {
        "link_count": link_count,
        "overlap_count": overlap_count,
        "old_overlap_count": old_overlap_count,
        "new_gene_context_count": new_gene_context_count,
        "direct_new_gene_count": direct_new_gene_count,
        "hit_link_count": len(hit_link_rows),
        "hit_canonical_pair_count": len({r["canonical_pair_id"] for r in hit_link_rows}),
        "hit_marker_count": len(hit_marker_keys),
        "unmatched_marker_count": len(unmatched_rows),
        "unique_unannotated_endpoint_count": len(endpoint_rows),
        "new_gene_annotated_endpoint_count": endpoint_new_gene_count,
        "overlap_tier_counts": overlap_tier_counts,
        "span_overlap_count": span_overlap_count,
        "span_strong_overlap_count": span_strong_overlap_count,
        "span_exact_count": span_exact_count,
        "span_nonexact_count": span_nonexact_count,
        "span_allele_mismatch_count": span_allele_mismatch_count,
        "span_hit_link_count": len(span_hit_link_ids),
        "span_strong_hit_link_count": len(span_strong_hit_link_ids),
        "span_hit_trait_variant_count": len(span_hit_trait_variants),
        "span_strong_hit_trait_variant_count": len(
            span_strong_hit_trait_variants
        ),
        "span_strong_match_class_counts": span_strong_match_class_counts,
        "span_strong_pcgene_count": span_strong_pcgene_count,
        "span_strong_pcgene_hit_link_count": len(
            span_strong_pcgene_hit_link_ids
        ),
        "span_strong_pcgene_hit_trait_variant_count": len(
            span_strong_pcgene_hit_trait_variants
        ),
    }


def write_candidate_inventory(
    run_dir: Path,
    candidates: Sequence[Dict[str, object]],
) -> Path:
    path = run_dir / "dog5_raw_gwas_candidate_markers_p_le_1e-4.tsv.gz"
    handle, writer = open_gzip_writer(path, GWAS_OUTPUT_COLUMNS)
    try:
        for row in candidates:
            writer.writerow(output_gwas_row(row))
    finally:
        handle.close()
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare current five-dog known_gene, manual_new_gene, and unannotated SNP records with raw GWAS results, "
            "distinguishing gene-span overlap from exact SNP matches while retaining risk flags."
        )
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
        "--new-gene-build-summary",
        type=Path,
        default=DEFAULT_NEW_GENE_BUILD_SUMMARY,
    )
    parser.add_argument("--gwas-dir", type=Path, default=DEFAULT_GWAS_DIR)
    parser.add_argument(
        "--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT
    )
    parser.add_argument(
        "--run-id",
        default=datetime.now().strftime("%Y%m%d_%H%M%S"),
        help="Run directory name; defaults to the current time in YYYYMMDD_HHMMSS format.",
    )
    return parser.parse_args()


def run(args: argparse.Namespace) -> int:
    for path, label in [
        (args.allelome_links, "Allelome link table"),
        (args.allelome_summary, "Allelome summary table"),
        (args.new_gene_gtf, "Manual new-only GTF"),
        (args.new_gene_build_summary, "Manual GTF construction summary table"),
    ]:
        require_file(path, label)
    if not args.gwas_dir.is_dir():
        raise FileNotFoundError(f"GWAS directory does not exist: {args.gwas_dir}")
    if not re.fullmatch(r"\d{8}_\d{6}", args.run_id):
        raise ValueError(f"run-id invalid: {args.run_id!r}")

    run_dir = args.output_root / args.run_id
    if run_dir.exists():
        raise FileExistsError(f"run directory already exists, refuse to overwrite: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=False)
    try:
        log("[stage 1/5] read current Allelome summary, all endpoint all are 0-based half-open")
        link_summary = read_metric_table(args.allelome_summary)
        expected_link_rows = parse_int(link_summary.get("annotated_rows", ""), "annotated_rows")
        marker = args.allelome_links.parent / "execution_completion.marker"
        marker_values = dict(line.split("=", 1) for line in marker.read_text().splitlines() if "=" in line)
        if marker_values.get("execution_status") != "success":
            raise ValueError("upstream Allelome run not completed")
        log("[stage 1/5] read Manual GTF gene span with five-dog novel gene attribute")
        new_gene_rows, new_gene_bins = load_new_gene_loci(
            args.new_gene_gtf,
            args.new_gene_build_summary,
        )
        autosomal_new_gene_rows = sum(
            row["chrom"] in AUTOSOMES for row in new_gene_rows
        )
        log(
            f"[stage 1/5] five_only_new_gene_rows={len(new_gene_rows)}, "
            f"autosomal_rows={autosomal_new_gene_rows}"
        )

        log("[stage 2/5] read p_lrt<=1e-4 raw GWAS candidate")
        (
            gwas_candidates,
            gwas_by_chrom,
            gwas_positions,
            gwas_by_position,
        ) = load_gwas_candidates(args.gwas_dir)
        (
            gwas_by_new_gene,
            gwas_new_gene_candidate_count,
        ) = annotate_gwas_new_gene_hits(gwas_candidates, new_gene_bins)
        candidate_inventory_path = write_candidate_inventory(
            run_dir, gwas_candidates
        )
        candidate_tier_counts = Counter(
            str(row["gwas_p_tier"]) for row in gwas_candidates
        )
        log(
            f"[stage 2/5] gwas_candidate_rows={len(gwas_candidates)}, "
            f"gwas_in_new_gene_rows={gwas_new_gene_candidate_count}, "
            f"candidate_inventory={candidate_inventory_path}"
        )

        log("[stage 3/5] streaming matching Allelome both endpoints and GWAS")
        result = process_links(
            args.allelome_links,
            expected_link_rows,
            new_gene_bins,
            gwas_candidates,
            gwas_by_chrom,
            gwas_positions,
            gwas_by_position,
            gwas_by_new_gene,
            run_dir,
        )

        log("[stage 4/5] write summary with run list")
        summary_rows: List[Dict[str, object]] = [
            {"metric": "summary_scope", "value": "descriptive_counts_not_postflight_QC"},
            {"metric": "allelome_link_rows", "value": result["link_count"]},
            {
                "metric": "five_only_new_gene_rows_all_contigs",
                "value": len(new_gene_rows),
            },
            {
                "metric": "five_only_new_gene_rows_autosomes",
                "value": autosomal_new_gene_rows,
            },
            {
                "metric": "gwas_file_count",
                "value": EXPECTED_GWAS_FILE_COUNT,
            },
            {
                "metric": "gwas_candidate_rows_p_le_1e-4",
                "value": len(gwas_candidates),
            },
            {
                "metric": "gwas_candidate_rows_in_manual_new_gene_spans",
                "value": gwas_new_gene_candidate_count,
            },
            {
                "metric": "overlap_rows_p_le_1e-4",
                "value": result["overlap_count"],
            },
            {
                "metric": "overlap_rows_old_p_lt_5e-5",
                "value": result["old_overlap_count"],
            },
            {
                "metric": "unique_hit_trait_variants",
                "value": result["hit_marker_count"],
            },
            {
                "metric": "unique_hit_sample_links",
                "value": result["hit_link_count"],
            },
            {"metric": "unique_hit_endpoint_pairs", "value": result["hit_canonical_pair_count"]},
            {
                "metric": "unmatched_gwas_candidate_rows",
                "value": result["unmatched_marker_count"],
            },
            {
                "metric": "unique_unannotated_endpoints",
                "value": result["unique_unannotated_endpoint_count"],
            },
            {
                "metric": "new_gene_annotated_unannotated_endpoints",
                "value": result["new_gene_annotated_endpoint_count"],
            },
            {
                "metric": "overlap_rows_with_new_gene_link_context",
                "value": result["new_gene_context_count"],
            },
            {
                "metric": "gwas_hit_new_gene_endpoint_rows",
                "value": result["direct_new_gene_count"],
            },
            {
                "metric": "same_new_gene_span_rows_p_le_1e-4",
                "value": result["span_overlap_count"],
            },
            {
                "metric": "same_new_gene_span_unique_trait_variants_p_le_1e-4",
                "value": result["span_hit_trait_variant_count"],
            },
            {
                "metric": "same_new_gene_span_unique_allelome_links_p_le_1e-4",
                "value": result["span_hit_link_count"],
            },
            {
                "metric": "same_new_gene_span_exact_rows_p_le_1e-4",
                "value": result["span_exact_count"],
            },
            {
                "metric": "same_new_gene_span_nonexact_rows_p_le_1e-4",
                "value": result["span_nonexact_count"],
            },
            {
                "metric": "same_new_gene_span_allele_mismatch_rows_p_le_1e-4",
                "value": result["span_allele_mismatch_count"],
            },
            {
                "metric": "same_new_gene_span_rows_p_le_1e-5",
                "value": result["span_strong_overlap_count"],
            },
            {
                "metric": "same_new_gene_span_unique_trait_variants_p_le_1e-5",
                "value": result["span_strong_hit_trait_variant_count"],
            },
            {
                "metric": "same_new_gene_span_unique_allelome_links_p_le_1e-5",
                "value": result["span_strong_hit_link_count"],
            },
            {
                "metric": "same_new_gene_span_exact_rows_p_le_1e-5",
                "value": result["span_strong_match_class_counts"][
                    "exact_unannotated_ase_variant"
                ],
            },
            {
                "metric": "same_new_gene_span_nonexact_rows_p_le_1e-5",
                "value": result["span_strong_match_class_counts"][
                    "same_manual_new_gene_span_nonexact"
                ],
            },
            {
                "metric": "same_new_gene_span_allele_mismatch_rows_p_le_1e-5",
                "value": result["span_strong_match_class_counts"][
                    "same_position_allele_mismatch"
                ],
            },
            {
                "metric": "same_new_gene_span_unannotated_to_pcgene_rows_p_le_1e-5",
                "value": result["span_strong_pcgene_count"],
            },
            {
                "metric": (
                    "same_new_gene_span_unannotated_to_pcgene_"
                    "unique_trait_variants_p_le_1e-5"
                ),
                "value": result[
                    "span_strong_pcgene_hit_trait_variant_count"
                ],
            },
            {
                "metric": (
                    "same_new_gene_span_unannotated_to_pcgene_"
                    "unique_allelome_links_p_le_1e-5"
                ),
                "value": result["span_strong_pcgene_hit_link_count"],
            },
        ]
        for tier in [
            "strict_p_le_5e-8",
            "strong_suggestive_p_le_1e-5",
            "loose_suggestive_p_le_5e-5",
            "very_loose_p_le_1e-4",
        ]:
            summary_rows.append(
                {
                    "metric": f"gwas_candidate_rows_{tier}",
                    "value": candidate_tier_counts[tier],
                }
            )
            summary_rows.append(
                {
                    "metric": f"overlap_rows_{tier}",
                    "value": result["overlap_tier_counts"][tier],
                }
            )
        write_tsv(
            run_dir / "dog5_allelome_raw_gwas_summary.tsv",
            ["metric", "value"],
            summary_rows,
        )

        manifest_rows = [
            {"item": "run_id", "value": args.run_id},
            {"item": "run_directory", "value": run_dir},
            {"item": "allelome_links", "value": args.allelome_links},
            {"item": "allelome_summary", "value": args.allelome_summary},
            {"item": "new_gene_gtf", "value": args.new_gene_gtf},
            {
                "item": "new_gene_build_summary",
                "value": args.new_gene_build_summary,
            },
            {"item": "gwas_directory", "value": args.gwas_dir},
            {"item": "gwas_pattern", "value": GWAS_PATTERN},
            {"item": "coordinate_convention", "value": "Allelome START0/END0 0-based half-open; GWAS POS 1-based"},
            {"item": "canonical_pair_definition", "value": "chromosome + unordered GeneID pair"},
            {"item": "sample_pair_definition", "value": "animal + tissue + canonical pair; unique link counts are sample-level"},
            {"item": "risk_policy", "value": "retain all links and upstream evidence_stratum; no causal claim"},
            {"item": "both_policy", "value": "one row per trait-variant and sample pair; both endpoint hits are not duplicated"},

            {
                "item": "gwas_candidate_rule",
                "value": "p_lrt<=1e-4",
            },
            {
                "item": "old_comparable_rule",
                "value": "p_lrt<5e-5",
            },
            {
                "item": "unannotated_endpoint_gwas_rule",
                "value": "exact_CHROM_POS_unordered_allele_set",
            },
            {
                "item": "annotated_gene_gwas_rule",
                "value": "START0<POS1<=END0",
            },
            {
                "item": "new_gene_annotation_rule",
                "value": (
                    "unannotated_variant_point0=POS1-1; "
                    "GTF_transcript_union_gene_start0<=point0<"
                    "GTF_transcript_union_gene_end0"
                ),
            },
            {
                "item": "new_gene_annotation_source",
                "value": (
                    "Manual new-only GTF transcript records grouped by gene_id; "
                    "five-dogs class_code/novelty_category/filter_reason from same GTF"
                ),
            },
            {
                "item": "new_gene_interpretation",
                "value": (
                    "known_gene and manual_new_gene use gene spans; "
                    "GWAS and unannotated ASE points are additionally compared "
                    "within shared Manual new-gene spans"
                ),
            },
            {
                "item": "gwas_new_gene_annotation_rule",
                "value": (
                    "gwas_point0=POS1-1; GTF_transcript_union_gene_start0<="
                    "point0<GTF_transcript_union_gene_end0"
                ),
            },
            {
                "item": "same_new_gene_span_rule",
                "value": (
                    "intersection(gwas_manual_new_gene_ids,"
                    "ase_manual_new_gene_ids) is non-empty"
                ),
            },
            {
                "item": "same_position_exact_rule",
                "value": "exact_CHROM_POS_unordered_allele_set",
            },
            {
                "item": "same_new_gene_span_primary_rule",
                "value": "gwas_p_lrt<=1e-5",
            },
            {
                "item": "same_new_gene_span_all_output",
                "value": (
                    run_dir
                    / "dog5_unannotated_ase_gwas_new_gene_span_overlap_all_"
                    "p_le_1e-4.tsv.gz"
                ),
            },
            {
                "item": "same_new_gene_span_primary_output",
                "value": (
                    run_dir
                    / "dog5_unannotated_ase_gwas_new_gene_span_overlap_strong_"
                    "p_le_1e-5.tsv.gz"
                ),
            },
            {
                "item": "same_new_gene_span_unannotated_to_pcgene_output",
                "value": (
                    run_dir
                    / "dog5_unannotated_ase_to_pcgene_gwas_new_gene_span_"
                    "strong_p_le_1e-5.tsv.gz"
                ),
            },
        ]
        write_tsv(
            run_dir / "dog5_allelome_raw_gwas_manifest.tsv",
            ["item", "value"],
            manifest_rows,
        )

        log(
            f"[stage 5/5] summary: overlap_rows={result['overlap_count']}, "
            f"hit_links={result['hit_link_count']}, "
            f"direct_new_gene_rows={result['direct_new_gene_count']}, "
            f"same_new_gene_rows_p_le_1e-5="
            f"{result['span_strong_overlap_count']}, "
            f"unannotated_to_pcgene_rows_p_le_1e-5="
            f"{result['span_strong_pcgene_count']}"
        )
        return 0
    except Exception as exc:
        log(f"[ERROR] {type(exc).__name__}: {exc}")
        raise

def stop(signum, frame):
    raise InterruptedError(f"received_signal={signum}")


def main() -> int:
    args = parse_args()
    if not re.fullmatch(r"\d{8}_\d{6}", args.run_id):
        raise ValueError("run-id must be YYYYMMDD_HHMMSS")
    log_dir = SCRIPT_DIR / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_root / args.run_id
    prefix = f"1.1.compare_dog5_allelome_raw_gwas.{args.run_id}"
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    with ExitStack() as stack:
        master = stack.enter_context((log_dir / (prefix + ".log")).open("x", encoding="utf-8", buffering=1))
        stack.enter_context(redirect_stdout(master))
        stack.enter_context(redirect_stderr(master))
        handles, writers = {}, {}
        for state in ("success", "failed"):
            handles[state] = stack.enter_context((log_dir / (prefix + f".{state}.tsv")).open(
                "x", encoding="utf-8", newline="", buffering=1))
            writers[state] = csv.DictWriter(handles[state],
                fieldnames=["run_id", "execution_status", "output_path", "reason"], delimiter="\t")
            writers[state].writeheader()
        success_header_end = handles["success"].tell()
        marker_owned = False
        signal.signal(signal.SIGTERM, stop)
        signal.signal(signal.SIGINT, stop)
        try:
            log(f"start; serial=1 threads=1(retain the original serial settings); output={output}")
            for name in ("allelome_links", "allelome_summary", "new_gene_gtf", "new_gene_build_summary", "gwas_dir"):
                log(f"{name}={getattr(args, name)}")
            run(args)
            with (output / "execution_completion.marker").open("x", encoding="utf-8") as handle:
                marker_owned = True
                handle.write(f"execution_status=success\nrun_id={args.run_id}\n")
            writers["success"].writerow(dict(run_id=args.run_id, execution_status="success", output_path=output, reason="normal_end"))
            log("task completed / done; execution_status=success")
            log("please separately check output format, count and biological interpretation; execution completed does not imply result QC passed.")
            return 0
        except BaseException as exc:
            if marker_owned:
                (output / "execution_completion.marker").unlink(missing_ok=True)
            handles["success"].seek(success_header_end)
            handles["success"].truncate()
            writers["failed"].writerow(dict(run_id=args.run_id, execution_status="failed", output_path=output, reason=str(exc)))
            traceback.print_exc()
            log("execution_status=failed")
            raise


if __name__ == "__main__":
    sys.exit(main())
