#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import gzip
import html
import re
import sys
import traceback
from contextlib import redirect_stdout, redirect_stderr
import zipfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, Mapping, Sequence
from xml.sax.saxutils import escape as xml_escape


DEFAULT_OVERLAP = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap/20260910_180641/main/dog5_allelome_gwas_credible_set_overlap.tsv.gz'
)
DEFAULT_INVENTORY = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap/20260910_180641/main/dog5_gwas_finemap_candidate_inventory.tsv.gz'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/ncrna-mrna'
)
BRANCH = "main"
PAIR_CLASS = "known_lncRNA--known_protein_coding"
RUN_LOG_NAME = "dog5_annotated_ncrna_pcgene_gwas_finemap.log"
SVG_ROW_LIMIT = 50

REQUIRED_OVERLAP_FIELDS = {
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
    "gwas_variant_p_lrt",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_credible_sets",
    "gwas_in_credible_set",
    "gwas_variant_count",
    "canonical_pair_id",
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
    "target_record_type",
    "base_gene_id",
    "base_gene_symbol",
    "base_gene_strand",
    "base_chr",
    "base_start",
    "base_end",
    "base_annotation_source",
    "base_gene_biotype",
    "base_endpoint_type",
    "target_gene_id",
    "target_gene_symbol",
    "target_gene_strand",
    "target_chr",
    "target_start",
    "target_end",
    "target_annotation_source",
    "target_gene_biotype",
    "target_endpoint_type",
}

REQUIRED_INVENTORY_FIELDS = {
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
    "gwas_variant_p_lrt",
    "gwas_pip",
    "gwas_in_credible_set",
    "gwas_credible_sets",
}

DERIVED_FIELDS = [
    "evidence_scope",
    "relation_direction",
    "gwas_variant_coordinate",
    "gwas_variant_display",
    "gwas_pip_rank_display",
    "gwas_credible_set_pip_rank",
    "gwas_credible_set_size",
    "gwas_credible_set_pip_rank_display",
    "gwas_credible_set_min_p_snp",
    "gwas_credible_set_min_p_snp_coordinate",
    "gwas_credible_set_min_p_snp_display",
    "gwas_credible_set_min_p_lrt",
    "ncrna_side",
    "ncrna_gene_id",
    "ncrna_gene_symbol",
    "ncrna_strand",
    "ncrna_chr",
    "ncrna_start",
    "ncrna_end",
    "ncrna_interval",
    "ncrna_biotype",
    "pcgene_side",
    "pcgene_gene_id",
    "pcgene_gene_symbol",
    "pcgene_strand",
    "pcgene_chr",
    "pcgene_start",
    "pcgene_end",
    "pcgene_interval",
    "pcgene_biotype",
    "gwas_hit_on_ncrna",
    "gwas_hit_on_pcgene",
    "gwas_hit_relation",
    "ncrna_mapping_rule",
    "pcgene_mapping_rule",
    "ncrna_pcgene_pair_key",
    "ncrna_pcgene_pair_label",
]

FOCUS_FIELDS = [
    "trait",
    "gwas_analysis_scope",
    "locus_id",
    "gwas_variant_display",
    "gwas_variant_p_lrt",
    "gwas_pip",
    "gwas_pip_rank_display",
    "gwas_credible_sets",
    "gwas_credible_set_pip_rank_display",
    "gwas_credible_set_min_p_snp_display",
    "gwas_credible_set_min_p_lrt",
    "ncrna_gene_id",
    "ncrna_gene_symbol",
    "ncrna_interval",
    "pcgene_gene_id",
    "pcgene_gene_symbol",
    "mechanism",
    "linkage_score",
    "sample_tissue",
    "interval_gap_bp",
    "gwas_hit_relation", "evidence_stratum", "overlap_bp", "either_endpoint_multigene",
    "shared_variant_identity_checked",
]

SUMMARY_FIELDS = ["metric", "value", "definition"]
PAIR_SUMMARY_FIELDS = [
    "ncrna_gene_id",
    "ncrna_gene_symbol",
    "pcgene_gene_id",
    "pcgene_gene_symbol",
    "pair_key",
    "credible_overlap_row_count",
    "unique_gwas_variant_count",
    "traits",
    "analysis_scopes",
    "animals",
    "tissues",
    "gwas_hit_on_ncrna_rows",
    "enhancing_rows",
    "repressing_rows",
    "minimum_gwas_variant_p_lrt",
    "minimum_credible_set_p_lrt",
    "maximum_linkage_score",
]

SVG_LABELS = {
    "trait": "GWAS Trait",
    "gwas_analysis_scope": "GWAS scope",
    "locus_id": "Locus",
    "gwas_variant_display": "Fine-mapped SNP (ID | chr:pos:alleles)",
    "gwas_variant_p_lrt": "Matched SNP P LRT",
    "gwas_pip": "Matched SNP PIP",
    "gwas_pip_rank_display": "Locus PIP Rank (rank/total)",
    "gwas_credible_sets": "Credible Set",
    "gwas_credible_set_pip_rank_display": "Credible-set PIP Rank (rank/size)",
    "gwas_credible_set_min_p_snp_display": "Credible-set Min-p SNP (ID | chr:pos:alleles)",
    "gwas_credible_set_min_p_lrt": "Credible-set Min P LRT",
    "ncrna_gene_id": "ncRNA ID",
    "ncrna_gene_symbol": "ncRNA",
    "ncrna_interval": "ncRNA interval",
    "pcgene_gene_id": "pcGene ID",
    "pcgene_gene_symbol": "pcGene",
    "mechanism": "Allelome Pattern",
    "linkage_score": "Allelome Score",
    "sample_tissue": "Sample Tissue",
    "interval_gap_bp": "Interval Distance (bp)",
}
SVG_WIDTHS = {
    "trait": 185,
    "gwas_analysis_scope": 95,
    "locus_id": 175,
    "gwas_variant_display": 245,
    "gwas_variant_p_lrt": 115,
    "gwas_pip": 100,
    "gwas_pip_rank_display": 100,
    "gwas_credible_sets": 90,
    "gwas_credible_set_pip_rank_display": 145,
    "gwas_credible_set_min_p_snp_display": 245,
    "gwas_credible_set_min_p_lrt": 115,
    "ncrna_gene_id": 135,
    "ncrna_gene_symbol": 125,
    "ncrna_interval": 190,
    "pcgene_gene_id": 135,
    "pcgene_gene_symbol": 125,
    "mechanism": 105,
    "linkage_score": 100,
    "sample_tissue": 190,
    "interval_gap_bp": 120,
}
SVG_GROUPS = (
    (
        "GWAS fine-mapping evidence",
        (
            "trait",
            "gwas_analysis_scope",
            "locus_id",
            "gwas_variant_display",
            "gwas_variant_p_lrt",
            "gwas_pip",
            "gwas_pip_rank_display",
            "gwas_credible_sets",
            "gwas_credible_set_pip_rank_display",
            "gwas_credible_set_min_p_snp_display",
            "gwas_credible_set_min_p_lrt",
        ),
    ),
    (
        "Conventional ncRNA–pcGene annotation",
        (
            "ncrna_gene_id",
            "ncrna_gene_symbol",
            "ncrna_interval",
            "pcgene_gene_id",
            "pcgene_gene_symbol",
        ),
    ),
    (
        "Allelome association evidence",
        ("mechanism", "linkage_score", "sample_tissue", "interval_gap_bp"),
    ),
)
SVG_GROUP_FILLS = {
    "GWAS fine-mapping evidence": "#e8f0ea",
    "Conventional ncRNA–pcGene annotation": "#eef2f7",
    "Allelome association evidence": "#f7efe7",
}


def now_string() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class RunLogger:


    def __init__(self, path: Path) -> None:
        self.handle = sys.stdout

    def log(self, message: str) -> None:
        line = f"[{now_string()} ] {message}"
        print(line, flush=True)

    def __enter__(self) -> "RunLogger":
        return self

    def __exit__(self, *_: object) -> None:
        self.handle.flush()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"input file no column header: {path}")
        fields = list(reader.fieldnames)
        rows = [{key: (value or "") for key, value in row.items()} for row in reader]
    return fields, rows


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def as_float(value: object) -> float | None:
    try:
        text = str(value).strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def as_int(value: object) -> int | None:
    try:
        text = str(value).strip()
        return int(float(text)) if text else None
    except (TypeError, ValueError):
        return None


def unique(values: Iterable[object]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            result.append(text)
    return result


def is_newgene_auxiliary_field(field: str) -> bool:
    lowered = field.strip().lower()
    return "new_gene" in lowered or "manual_new_gene" in lowered


def retained_input_fields(fields: Sequence[str]) -> list[str]:
    return [field for field in fields if not is_newgene_auxiliary_field(field)]


def normalize_chr(value: object) -> str:
    return str(value or "").strip().lower().removeprefix("chr")


def side_value(row: Mapping[str, str], side: str, name: str) -> str:
    return row.get(f"{side}_{name}", "").strip()


def side_is_ncrna(row: Mapping[str, str], side: str) -> bool:
    return side_value(row, side, "record_type") == "known_gene" and side_value(row, side, "endpoint_type") == "known_lncRNA"


def side_is_pcgene(row: Mapping[str, str], side: str) -> bool:
    return side_value(row, side, "record_type") == "known_gene" and side_value(row, side, "gene_biotype") == "protein_coding"


def standard_pair(row: Mapping[str, str]) -> bool:
    if row.get("pair_type", "").strip() != PAIR_CLASS:
        return False
    return (side_is_ncrna(row, "base") and side_is_pcgene(row, "target")) or (
        side_is_ncrna(row, "target") and side_is_pcgene(row, "base")
    )


def coordinate_contains(row: Mapping[str, str], side: str) -> bool:
    gwas_chr = normalize_chr(row.get("gwas_chrom", ""))
    gene_chr = normalize_chr(side_value(row, side, "chr"))
    position = as_int(row.get("gwas_pos", ""))
    start = as_int(side_value(row, side, "start"))
    end = as_int(side_value(row, side, "end"))

    return bool(
        gwas_chr
        and gene_chr
        and gwas_chr == gene_chr
        and position is not None
        and start is not None
        and end is not None
        and start < position <= end
    )


def gwas_hits_side(row: Mapping[str, str], side: str) -> bool:

    return coordinate_contains(row, side)


def split_set_ids(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"[,;]", value or "") if part.strip()]


def format_coordinate(row: Mapping[str, str]) -> str:
    chrom = str(row.get("gwas_chrom", "")).strip()
    chrom = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    pos = str(row.get("gwas_pos", "")).strip()
    a1 = str(row.get("gwas_ld_a1", "")).strip()
    a2 = str(row.get("gwas_ld_a2", "")).strip()
    return f"{chrom}:{pos}:{a1}/{a2}"


def build_inventory_index(
    rows: Sequence[Mapping[str, str]],
) -> tuple[
    dict[tuple[str, str, str, str], Mapping[str, str]],
    dict[tuple[str, str, str], Mapping[str, str]],
    dict[tuple[str, str, str, str], tuple[int, int]],
]:
    variant_index: dict[tuple[str, str, str, str], Mapping[str, str]] = {}
    minimum_index: dict[tuple[str, str, str], Mapping[str, str]] = {}
    credible_set_members: dict[tuple[str, str, str], list[tuple[str, float]]] = defaultdict(list)
    for row in rows:
        if row.get("gwas_branch") != BRANCH:
            continue
        if row.get("candidate_type") != "susie_credible_set":
            continue
        if row.get("gwas_in_credible_set") != "yes":
            continue
        scope = row.get("gwas_analysis_scope", "")
        locus = row.get("locus_id", "")
        variant = row.get("gwas_ld_variant_id", "")
        variant_key = (scope, BRANCH, locus, variant)
        if variant_key in variant_index:
            raise ValueError(f"candidate inventory contains duplicate credible set SNP: {scope}/{locus}/{variant}")
        variant_index[variant_key] = row
        pip = as_float(row.get("gwas_pip", ""))
        if pip is None:
            raise ValueError(f"candidate inventory credible set SNP of PIP is not numeric: {scope}/{locus}/{variant}")
        for cs_id in split_set_ids(row.get("gwas_credible_sets", "")):
            key = (scope, locus, cs_id)
            credible_set_members[key].append((variant, pip))
            old = minimum_index.get(key)
            if old is None or (as_float(row.get("gwas_variant_p_lrt", "")) or float("inf")) < (as_float(old.get("gwas_variant_p_lrt", "")) or float("inf")):
                minimum_index[key] = row


    credible_set_rank_index: dict[tuple[str, str, str, str], tuple[int, int]] = {}
    for (scope, locus, cs_id), members in credible_set_members.items():
        variants = [variant for variant, _pip in members]
        if len(variants) != len(set(variants)):
            raise ValueError(f"candidate inventory of credible set member duplicate: {scope}/{locus}/{cs_id}")
        set_size = len(members)
        pip_values = [pip for _variant, pip in members]
        for variant, pip in members:
            rank = 1 + sum(other_pip > pip for other_pip in pip_values)
            credible_set_rank_index[(scope, locus, cs_id, variant)] = (rank, set_size)
    return variant_index, minimum_index, credible_set_rank_index


def oriented_enriched_row(
    row: Mapping[str, str],
    fields: Sequence[str],
    variant_index: Mapping[tuple[str, str, str, str], Mapping[str, str]],
    minimum_index: Mapping[tuple[str, str, str], Mapping[str, str]],
    credible_set_rank_index: Mapping[tuple[str, str, str, str], tuple[int, int]],
) -> dict[str, str]:
    if side_is_ncrna(row, "base") and side_is_pcgene(row, "target"):
        ncrna_side, pc_side = "base", "target"
    elif side_is_ncrna(row, "target") and side_is_pcgene(row, "base"):
        ncrna_side, pc_side = "target", "base"
    else:
        raise ValueError(f"cannot orient annotated_ncRNA--pcGene: {row.get('canonical_pair_id', '')}")

    scope = row.get("gwas_analysis_scope", "")
    locus = row.get("locus_id", "")
    variant_id = row.get("gwas_ld_variant_id", "")
    candidate = variant_index.get((scope, BRANCH, locus, variant_id))
    if candidate is None:
        raise ValueError(f"candidate inventory missing overlap SNP: {scope}/{locus}/{variant_id}")
    for field in ("trait", "gwas_chrom", "gwas_pos", "gwas_ld_a1", "gwas_ld_a2", "gwas_variant_p_lrt", "gwas_pip"):
        if candidate.get(field, "") != row.get(field, ""):
            raise ValueError(f"overlap with candidate inventory of {field} mismatch: {variant_id}")

    cs_ids = split_set_ids(row.get("gwas_credible_sets", ""))
    minima = [minimum_index[(scope, locus, cs_id)] for cs_id in cs_ids if (scope, locus, cs_id) in minimum_index]
    if not minima:
        raise ValueError(f"cannot find credible set minimum p_lrt SNP: {scope}/{locus}/{cs_ids}")
    minimum = min(minima, key=lambda item: as_float(item.get("gwas_variant_p_lrt", "")) or float("inf"))
    credible_set_ranks: list[str] = []
    credible_set_sizes: list[str] = []
    credible_set_rank_displays: list[str] = []
    for cs_id in cs_ids:
        rank_and_size = credible_set_rank_index.get((scope, locus, cs_id, variant_id))
        if rank_and_size is None:
            raise ValueError(f"cannot find the within-credible-set PIP ranking: {scope}/{locus}/{cs_id}/{variant_id}")
        rank, set_size = rank_and_size
        credible_set_ranks.append(f"{cs_id}:{rank}")
        credible_set_sizes.append(f"{cs_id}:{set_size}")
        credible_set_rank_displays.append(f"{cs_id}:{rank}/{set_size}")
    ncrna_hit = gwas_hits_side(row, ncrna_side)
    pc_hit = gwas_hits_side(row, pc_side)
    hit_relation = "both" if ncrna_hit and pc_hit else "ncrna_only" if ncrna_hit else "pcgene_only" if pc_hit else "neither"
    ncrna_id = side_value(row, ncrna_side, "gene_id") or side_value(row, ncrna_side, "gene_symbol")
    pc_id = side_value(row, pc_side, "gene_id") or side_value(row, pc_side, "gene_symbol")
    ncrna_symbol = side_value(row, ncrna_side, "gene_symbol") or ncrna_id
    pc_symbol = side_value(row, pc_side, "gene_symbol") or pc_id

    result = {field: str(row.get(field, "")) for field in fields}
    result.update(
        {
            "evidence_scope": "main_susie_credible_set",
            "relation_direction": "undirected_Allelome_pair",
            "gwas_variant_coordinate": format_coordinate(row),
            "gwas_variant_display": f"{variant_id} | {format_coordinate(row)}",
            "gwas_pip_rank_display": f"{row.get('gwas_pip_rank', '')}/{row.get('gwas_variant_count', '')}",
            "gwas_credible_set_pip_rank": ";".join(credible_set_ranks),
            "gwas_credible_set_size": ";".join(credible_set_sizes),
            "gwas_credible_set_pip_rank_display": ";".join(credible_set_rank_displays),
            "gwas_credible_set_min_p_snp": minimum.get("gwas_ld_variant_id", "NA"),
            "gwas_credible_set_min_p_snp_coordinate": format_coordinate({
                "gwas_chrom": minimum.get("gwas_chrom", ""),
                "gwas_pos": minimum.get("gwas_pos", ""),
                "gwas_ld_a1": minimum.get("gwas_ld_a1", ""),
                "gwas_ld_a2": minimum.get("gwas_ld_a2", ""),
            }),
            "gwas_credible_set_min_p_snp_display": f"{minimum.get('gwas_ld_variant_id', 'NA')} | {format_coordinate({'gwas_chrom': minimum.get('gwas_chrom', ''), 'gwas_pos': minimum.get('gwas_pos', ''), 'gwas_ld_a1': minimum.get('gwas_ld_a1', ''), 'gwas_ld_a2': minimum.get('gwas_ld_a2', '')})}",
            "gwas_credible_set_min_p_lrt": minimum.get("gwas_variant_p_lrt", "NA"),
            "ncrna_side": ncrna_side,
            "ncrna_gene_id": ncrna_id,
            "ncrna_gene_symbol": ncrna_symbol,
            "ncrna_strand": side_value(row, ncrna_side, "gene_strand"),
            "ncrna_chr": side_value(row, ncrna_side, "chr"),
            "ncrna_start": side_value(row, ncrna_side, "start"),
            "ncrna_end": side_value(row, ncrna_side, "end"),
            "ncrna_interval": f"{side_value(row, ncrna_side, 'chr')}:{int(side_value(row, ncrna_side, 'start')) + 1}-{side_value(row, ncrna_side, 'end')}",
            "ncrna_biotype": side_value(row, ncrna_side, "gene_biotype"),
            "pcgene_side": pc_side,
            "pcgene_gene_id": pc_id,
            "pcgene_gene_symbol": pc_symbol,
            "pcgene_strand": side_value(row, pc_side, "gene_strand"),
            "pcgene_chr": side_value(row, pc_side, "chr"),
            "pcgene_start": side_value(row, pc_side, "start"),
            "pcgene_end": side_value(row, pc_side, "end"),
            "pcgene_interval": f"{side_value(row, pc_side, 'chr')}:{int(side_value(row, pc_side, 'start')) + 1}-{side_value(row, pc_side, 'end')}",
            "pcgene_biotype": side_value(row, pc_side, "gene_biotype"),
            "gwas_hit_on_ncrna": "yes" if ncrna_hit else "no",
            "gwas_hit_on_pcgene": "yes" if pc_hit else "no",
            "gwas_hit_relation": hit_relation,
            "ncrna_mapping_rule": row.get(f"{ncrna_side}_mapping_rule", ""),
            "pcgene_mapping_rule": row.get(f"{pc_side}_mapping_rule", ""),
            "ncrna_pcgene_pair_key": f"{ncrna_id}|{pc_id}",
            "ncrna_pcgene_pair_label": f"{ncrna_symbol} ↔ {pc_symbol}",
        }
    )
    return result


def build_pair_summary(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row.get("ncrna_pcgene_pair_key", "")].append(row)
    output: list[dict[str, str]] = []
    for key, members in groups.items():
        first = members[0]
        p_values = [value for value in (as_float(row.get("gwas_variant_p_lrt", "")) for row in members) if value is not None]
        min_cs = [value for value in (as_float(row.get("gwas_credible_set_min_p_lrt", "")) for row in members) if value is not None]
        links = [value for value in (as_float(row.get("linkage_score", "")) for row in members) if value is not None]
        output.append(
            {
                "ncrna_gene_id": first.get("ncrna_gene_id", ""),
                "ncrna_gene_symbol": first.get("ncrna_gene_symbol", ""),
                "pcgene_gene_id": first.get("pcgene_gene_id", ""),
                "pcgene_gene_symbol": first.get("pcgene_gene_symbol", ""),
                "pair_key": key,
                "credible_overlap_row_count": str(len(members)),
                "unique_gwas_variant_count": str(len(unique(f"{row.get('trait','')}|{row.get('gwas_ld_variant_id','')}" for row in members))),
                "traits": ";".join(unique(row.get("trait", "") for row in members)),
                "analysis_scopes": ";".join(unique(row.get("gwas_analysis_scope", "") for row in members)),
                "animals": ";".join(unique(row.get("animal_id", "") for row in members)),
                "tissues": ";".join(unique(row.get("tissue", "") for row in members)),
                "gwas_hit_on_ncrna_rows": str(sum(row.get("gwas_hit_on_ncrna") == "yes" for row in members)),
                "enhancing_rows": str(sum(row.get("mechanism", "").lower() == "enhancing" for row in members)),
                "repressing_rows": str(sum(row.get("mechanism", "").lower() == "repressing" for row in members)),
                "minimum_gwas_variant_p_lrt": f"{min(p_values):.12g}" if p_values else "",
                "minimum_credible_set_p_lrt": f"{min(min_cs):.12g}" if min_cs else "",
                "maximum_linkage_score": f"{max(links):.6g}" if links else "",
            }
        )
    output.sort(key=lambda row: (-int(row["credible_overlap_row_count"]), as_float(row["minimum_gwas_variant_p_lrt"]) or float("inf")))
    return output


def inventory_summary(rows: Sequence[Mapping[str, str]]) -> dict[str, str]:

    credible = [
        row
        for row in rows
        if row.get("gwas_branch") == BRANCH
        and row.get("candidate_type") == "susie_credible_set"
        and row.get("gwas_in_credible_set") == "yes"
    ]
    unique_keys = {
        (
            row.get("gwas_chrom", ""),
            row.get("gwas_pos", ""),
            row.get("gwas_ld_variant_id", ""),
        )
        for row in credible
    }
    credible_set_keys = {
        (row.get("gwas_analysis_scope", ""), row.get("locus_id", ""), cs_id)
        for row in credible
        for cs_id in split_set_ids(row.get("gwas_credible_sets", ""))
    }
    low_p = [
        row
        for row in credible
        if (value := as_float(row.get("gwas_variant_p_lrt", ""))) is not None
        and value < 1e-5
    ]
    low_p_keys = {
        (
            row.get("gwas_chrom", ""),
            row.get("gwas_pos", ""),
            row.get("gwas_ld_variant_id", ""),
        )
        for row in low_p
    }
    minimum = min(
        credible,
        key=lambda row: as_float(row.get("gwas_variant_p_lrt", ""))
        if as_float(row.get("gwas_variant_p_lrt", "")) is not None
        else float("inf"),
        default=None,
    )
    if minimum is None:
        return {
            "credible_set_count": "0",
            "credible_set_membership_rows": "0",
            "unique_credible_set_snps": "0",
            "p_lrt_lt_1e-5_membership_rows": "0",
            "p_lrt_lt_1e-5_unique_snps": "0",
            "min_gwas_p_lrt": "NA",
            "min_gwas_p_lrt_snp": "NA",
            "min_gwas_p_lrt_snp_display": "NA",
        }
    return {
        "credible_set_count": str(len(credible_set_keys)),
        "credible_set_membership_rows": str(len(credible)),
        "unique_credible_set_snps": str(len(unique_keys)),
        "p_lrt_lt_1e-5_membership_rows": str(len(low_p)),
        "p_lrt_lt_1e-5_unique_snps": str(len(low_p_keys)),
        "min_gwas_p_lrt": minimum.get("gwas_variant_p_lrt", "NA"),
        "min_gwas_p_lrt_snp": minimum.get("gwas_ld_variant_id", "NA"),
        "min_gwas_p_lrt_snp_display": (
            f"{minimum.get('gwas_ld_variant_id', 'NA')} | {format_coordinate(minimum)}"
        ),
    }


def build_summary(rows: Sequence[Mapping[str, str]], focus_rows: Sequence[Mapping[str, str]], inventory_rows: Sequence[Mapping[str, str]], overlap_rows: int, run_dir: Path) -> list[dict[str, str]]:
    gwas_summary = inventory_summary(inventory_rows)
    metrics = [
        ("overlap_input_rows", overlap_rows, "Total input rows for main-branch credible-set overlaps."),
        ("candidate_inventory_rows", len(inventory_rows), "Total candidate inventory input rows."),
        ("main_credible_set_count", gwas_summary["credible_set_count"], "Number of SuSiE credible sets in the main branch."),
        ("main_credible_set_membership_rows", gwas_summary["credible_set_membership_rows"], "Expanded SuSiE credible-set member rows in the main branch."),
        ("main_unique_credible_set_snps", gwas_summary["unique_credible_set_snps"], "Unique SNPs in main-branch credible sets, deduplicated by chromosome, position, and marker ID."),
        ("main_p_lrt_lt_1e-5_membership_rows", gwas_summary["p_lrt_lt_1e-5_membership_rows"], "Member rows with p_lrt<1e-5 in main credible sets."),
        ("main_p_lrt_lt_1e-5_unique_snps", gwas_summary["p_lrt_lt_1e-5_unique_snps"], "Unique SNPs with p_lrt<1e-5 in main credible sets."),
        ("main_min_gwas_p_lrt", gwas_summary["min_gwas_p_lrt"], "Lowest p_lrt among all main credible-set members."),
        ("main_min_gwas_p_lrt_snp", gwas_summary["min_gwas_p_lrt_snp_display"], "SNP with the lowest p_lrt across all main credible sets."),
        ("annotated_ncrna_pcgene_credible_rows", len(rows), "Expanded annotated_ncRNA--pcGene rows in main SuSiE credible sets."),
        ("focus_gwas_hit_on_ncrna_rows", len(focus_rows), "Rows above in which the GWAS fine-map SNP falls within the ncRNA gene span."),
        ("unique_annotated_ncrna_pcgene_pairs", len(unique(row.get("ncrna_pcgene_pair_key", "") for row in rows)), "Number of relationships deduplicated by ncRNA gene ID and pcGene ID."),
        ("unique_focus_pairs", len(unique(row.get("ncrna_pcgene_pair_key", "") for row in focus_rows)), "Number of deduplicated relationships with at least one SNP overlapping the ncRNA."),
        ("unique_ncrna_count", len(unique(row.get("ncrna_gene_id", "") for row in rows)), "Number of annotated ncRNAs involved."),
        ("unique_pcgene_count", len(unique(row.get("pcgene_gene_id", "") for row in rows)), "Number of pcGenes involved."),
        ("primary_rows", sum(row.get("gwas_analysis_scope") == "primary" for row in rows), "Expanded primary-scope rows."),
        ("exploratory_rows", sum(row.get("gwas_analysis_scope") == "exploratory" for row in rows), "Expanded exploratory-scope rows."),
        ("minimum_matched_gwas_p_lrt", min((as_float(row.get("gwas_variant_p_lrt", "")) for row in rows if as_float(row.get("gwas_variant_p_lrt", "")) is not None), default="NA"), "Lowest p_lrt among matched fine-map SNPs."),
        ("output_run_dir", str(run_dir), "Timestamped output directory for this run."),
    ]
    return [{"metric": str(key), "value": str(value), "definition": definition} for key, value, definition in metrics]


def write_html(path: Path, title: str, fields: Sequence[str], rows: Sequence[Mapping[str, str]], note: str) -> None:
    headers = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields) + "</tr>"
        for row in rows
    ) or '<tr><td colspan="2">No eligible records</td></tr>'
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font-family:Arial,'Microsoft YaHei',sans-serif;margin:24px;color:#1f2933}}h1{{font-size:20px}}p{{color:#52606d;font-size:13px}}table{{border-collapse:collapse;min-width:1800px;font-size:12px}}th,td{{border:1px solid #cbd2d9;padding:5px 7px;white-space:nowrap}}th{{background:#e8f0ec;position:sticky;top:0}}tr:nth-child(even){{background:#f8faf9}}</style></head><body><h1>{html.escape(title)}</h1><p>{html.escape(note)}</p><table><thead><tr>{headers}</tr></thead><tbody>{body}</tbody></table></body></html>"""
    path.write_text(document, encoding="utf-8")


def svg_escape(value: object) -> str:
    return xml_escape(str(value or ""), {"\"": "&quot;", "'": "&apos;"})


def svg_label_lines(label: str, max_chars: int = 18) -> list[str]:
    if len(label) <= max_chars:
        return [label]
    words = label.split(" ")
    if len(words) > 1:
        first: list[str] = []
        second: list[str] = []
        for word in words:
            if not second and sum(len(item) for item in first) + len(first) + len(word) <= max_chars:
                first.append(word)
            else:
                second.append(word)
        if first and second:
            return [" ".join(first), " ".join(second)]
    return [label[: max_chars - 1] + "…", label[max_chars - 1 : max_chars * 2 - 1]]


def svg_value(field: str, value: object) -> str:
    text = str(value or "")
    if field in {"gwas_variant_p_lrt", "gwas_credible_set_min_p_lrt"}:
        number = as_float(text)
        return f"{number:.3e}" if number is not None else text
    if field in {"gwas_pip", "linkage_score"}:
        number = as_float(text)
        return f"{number:.3f}" if number is not None else text
    if field in {"interval_gap_bp"}:
        number = as_int(text)
        return f"{number:,}" if number is not None else text
    return text


def svg_truncate(value: str, width: int) -> str:
    limit = max(8, int(width / 7.0))
    return value if len(value) <= limit else value[: limit - 1] + "…"


def write_svg_table(path: Path, title: str, rows: Sequence[Mapping[str, str]], total_rows: int, summary_text: str) -> None:
    fields = [field for field in FOCUS_FIELDS if field in SVG_LABELS]
    display_rows = list(rows[:SVG_ROW_LIMIT])
    widths = [SVG_WIDTHS[field] for field in fields]
    left = 20
    title_height = 84
    group_height = 24
    header_height = 48
    row_height = 31
    footer_height = 32
    width = left * 2 + sum(widths)
    height = title_height + group_height + header_height + row_height * max(1, len(display_rows)) + footer_height
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    parts.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    parts.append(f'<text x="{width / 2:.1f}" y="30" font-size="16" font-family="Helvetica,Arial,sans-serif" font-weight="400" fill="#111111" text-anchor="middle">{svg_escape(title)}</text>')
    parts.append(f'<text x="{width / 2:.1f}" y="50" font-size="11" font-family="Helvetica,Arial,sans-serif" fill="#111111" text-anchor="middle">{svg_escape(summary_text)}</text>')
    parts.append(f'<text x="{width / 2:.1f}" y="68" font-size="11" font-family="Helvetica,Arial,sans-serif" fill="#111111" text-anchor="middle">main branch; focus rows are GWAS SNPs inside conventional ncRNA spans; displaying {len(display_rows)}/{total_rows}</text>')
    positions: dict[str, tuple[int, int]] = {}
    x = left
    for field, field_width in zip(fields, widths):
        positions[field] = (x, field_width)
        x += field_width
    group_y = title_height
    header_y = group_y + group_height
    data_y = header_y + header_height
    used: set[str] = set()
    for group_label, group_fields in SVG_GROUPS:
        present = [field for field in group_fields if field in positions]
        if not present:
            continue
        start_x = positions[present[0]][0]
        end_x = positions[present[-1]][0] + positions[present[-1]][1]
        fill = SVG_GROUP_FILLS.get(group_label, "#f0ebe5")
        parts.append(f'<rect x="{start_x}" y="{group_y}" width="{end_x - start_x}" height="{group_height}" fill="{fill}" stroke="#9e9e9e"/>')
        parts.append(f'<text x="{(start_x + end_x) / 2:.1f}" y="{group_y + 16}" font-size="10.5" font-family="Helvetica,Arial,sans-serif" font-weight="700" fill="#111111" text-anchor="middle">{svg_escape(group_label)}</text>')
        used.update(present)
    for field in fields:
        if field not in used:
            x_i, width_i = positions[field]
            parts.append(f'<rect x="{x_i}" y="{group_y}" width="{width_i}" height="{group_height}" fill="#f0ebe5" stroke="#9e9e9e"/>')
    for field in fields:
        x_i, width_i = positions[field]
        parts.append(f'<rect x="{x_i}" y="{header_y}" width="{width_i}" height="{header_height}" fill="#f7f7f7" stroke="#9e9e9e"/>')
        lines = svg_label_lines(SVG_LABELS[field])
        first_y = header_y + (19 if len(lines) == 1 else 14)
        for line_index, line in enumerate(lines[:2]):
            parts.append(f'<text x="{x_i + width_i / 2:.1f}" y="{first_y + line_index * 12}" font-size="9.5" font-family="Helvetica,Arial,sans-serif" font-weight="700" fill="#111111" text-anchor="middle">{svg_escape(line)}</text>')
    previous_group: tuple[str, str] | None = None
    for index, row in enumerate(display_rows):
        y = data_y + index * row_height
        fill = "#ffffff" if index % 2 == 0 else "#fafafa"
        current_group = (row.get("trait", ""), row.get("gwas_variant_display", ""))
        if previous_group is not None and current_group != previous_group:
            parts.append(f'<line x1="{left}" y1="{y}" x2="{left + sum(widths)}" y2="{y}" stroke="#777777" stroke-width="1.1"/>')
        previous_group = current_group
        for field in fields:
            x_i, width_i = positions[field]
            value = svg_value(field, row.get(field, ""))
            shown = svg_truncate(value, width_i)
            parts.append(f'<g><title>{svg_escape(value)}</title><rect x="{x_i}" y="{y}" width="{width_i}" height="{row_height}" fill="{fill}" stroke="#bdbdbd"/>')
            if field == "mechanism":
                chip_fill = "#fff0e0" if value.lower() == "enhancing" else "#fbe7e8" if value.lower() == "repressing" else fill
                parts.append(f'<rect x="{x_i + 6}" y="{y + 7}" width="{max(12, width_i - 12)}" height="{row_height - 14}" rx="8" fill="{chip_fill}"/>')
            text_color = "#9a4d00" if field == "mechanism" and value.lower() == "enhancing" else "#8b1e2d" if field == "mechanism" and value.lower() == "repressing" else "#111111"
            parts.append(f'<text x="{x_i + width_i / 2:.1f}" y="{y + 20}" font-size="9.2" font-family="Helvetica,Arial,sans-serif" fill="{text_color}" text-anchor="middle">{svg_escape(shown)}</text></g>')
    footer_y = data_y + row_height * max(1, len(display_rows)) + 21
    parts.append(f'<text x="{left}" y="{footer_y}" font-size="9.5" font-family="Helvetica,Arial,sans-serif" fill="#52606d">Full evidence tables: TSV/CSV/XLSX; SVG shows only the first {SVG_ROW_LIMIT} rows ranked by matched SNP p_lrt.</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def write_svg_pair_barplot(path: Path, pair_rows: Sequence[Mapping[str, str]]) -> None:
    selected = list(pair_rows[:20])
    left = 360
    plot_width = 820
    row_height = 30
    title_height = 70
    width = left + plot_width + 120
    height = title_height + max(1, len(selected)) * row_height + 35
    max_count = max((int(row["credible_overlap_row_count"]) for row in selected), default=1)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="#ffffff"/>']
    parts.append(f'<text x="{width / 2:.1f}" y="28" font-size="16" font-family="Helvetica,Arial,sans-serif" text-anchor="middle">Top ncRNA–pcGene pairs by credible-set overlap rows</text>')
    parts.append('<text x="20" y="51" font-size="10.5" font-family="Helvetica,Arial,sans-serif" fill="#52606d">Bars count expanded GWAS credible-set × Allelome rows; orange marks rows whose SNP lies inside the ncRNA span.</text>')
    for index, row in enumerate(selected):
        y = title_height + index * row_height
        label = f"{row.get('ncrna_gene_symbol') or row.get('ncrna_gene_id')} ↔ {row.get('pcgene_gene_symbol') or row.get('pcgene_gene_id')}"
        if len(label) > 48:
            label = label[:47] + "…"
        total = int(row["credible_overlap_row_count"])
        hit = int(row["gwas_hit_on_ncrna_rows"])
        total_width = max(1, int(plot_width * total / max_count))
        hit_width = int(plot_width * hit / max_count)
        parts.append(f'<text x="20" y="{y + 19}" font-size="10.5" font-family="Helvetica,Arial,sans-serif">{svg_escape(label)}</text>')
        parts.append(f'<rect x="{left}" y="{y + 5}" width="{total_width}" height="17" fill="#b0b0b0"/>')
        if hit_width:
            parts.append(f'<rect x="{left}" y="{y + 5}" width="{hit_width}" height="17" fill="#d97706"/>')
        parts.append(f'<text x="{left + total_width + 6}" y="{y + 18}" font-size="9.5" font-family="Helvetica,Arial,sans-serif">{total} / {hit}</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def xlsx_cell(value: object) -> str:
    return f'<c t="inlineStr"><is><t>{xml_escape(str(value if value is not None else ""))}</t></is></c>'


def write_xlsx(path: Path, sheets: Sequence[tuple[str, Sequence[str], Sequence[Mapping[str, object]]]]) -> None:
    content_types = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>', '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>']
    sheet_nodes: list[str] = []
    rel_nodes: list[str] = []
    for index, (name, fields, _rows) in enumerate(sheets, 1):
        sheet_nodes.append(f'<sheet name="{xml_escape(name)}" sheetId="{index}" r:id="rId{index}"/>')
        rel_nodes.append(f'<Relationship Id="rId{index}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{index}.xml"/>')
        content_types.append(f'<Override PartName="/xl/worksheets/sheet{index}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
    content_types.append('</Types>')
    workbook = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' + "".join(sheet_nodes) + '</sheets></workbook>'
    rels = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' + "".join(rel_nodes) + '</Relationships>'
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", "".join(content_types))
        archive.writestr("_rels/.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", rels)
        for index, (_name, fields, rows) in enumerate(sheets, 1):
            xml_rows = ["<row>" + "".join(xlsx_cell(field) for field in fields) + "</row>"]
            xml_rows.extend("<row>" + "".join(xlsx_cell(row.get(field, "")) for field in fields) + "</row>" for row in rows)
            sheet_xml = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + "".join(xml_rows) + '</sheetData></worksheet>'
            archive.writestr(f"xl/worksheets/sheet{index}.xml", sheet_xml)


def xls_cell(value: object) -> str:
    return f'<Cell><Data ss:Type="String">{xml_escape(str(value if value is not None else ""))}</Data></Cell>'


def write_xls(path: Path, sheets: Sequence[tuple[str, Sequence[str], Sequence[Mapping[str, object]]]]) -> None:
    parts = ['<?xml version="1.0"?>', '<?mso-application progid="Excel.Sheet"?>', '<Workbook xmlns="urn:schemas-microsoft-com:office:spreadsheet" xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet">']
    for name, fields, rows in sheets:
        parts.append(f'<Worksheet ss:Name="{xml_escape(name)}"><Table><Row>{"".join(xls_cell(field) for field in fields)}</Row>')
        parts.extend(f'<Row>{"".join(xls_cell(row.get(field, "")) for field in fields)}</Row>' for row in rows)
        parts.append('</Table></Worksheet>')
    parts.append('</Workbook>')
    path.write_text("".join(parts), encoding="utf-8")


def write_manifest(path: Path, overlap_path: Path, inventory_path: Path, run_dir: Path, raw_fields: Sequence[str], output_fields: Sequence[str], counts: Mapping[str, object]) -> None:
    rows = [
        {"key": "script", "value": "build_annotated_ncrna_pcgene_gwas_finemap_allelome.py"},
        {"key": "overlap_input", "value": str(overlap_path)},
        {"key": "candidate_inventory_input", "value": str(inventory_path)},
        {"key": "output_run_dir", "value": str(run_dir)},
        {"key": "gwas_branch_rule", "value": "gwas_branch == main"},
        {"key": "credible_set_rule", "value": "candidate_type == susie_credible_set and gwas_in_credible_set == yes"},
        {"key": "credible_set_pip_rank_rule", "value": "within gwas_analysis_scope+locus_id+credible_set_id; descending PIP; ties use minimum rank"},
        {"key": "pair_class_rule", "value": PAIR_CLASS},
        {"key": "coordinate_rule", "value": "START0 < 1-based POS1 <= END0"},
        {"key": "raw_input_columns", "value": str(len(raw_fields))},
        {"key": "retained_core_columns", "value": str(len(output_fields))},
        {"key": "excluded_newgene_columns", "value": ",".join(field for field in raw_fields if is_newgene_auxiliary_field(field))},
    ]
    rows.extend({"key": key, "value": str(value)} for key, value in counts.items())
    write_tsv(path, ["key", "value"], rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract annotated ncRNA-pcGene Allelome links from the main GWAS fine-mapping results.")
    parser.add_argument("--overlap", type=Path, default=DEFAULT_OVERLAP, help="Absolute path to the main credible-set overlap TSV/TSV.GZ file.")
    parser.add_argument("--candidate-inventory", type=Path, default=DEFAULT_INVENTORY, help="Absolute path to the main candidate inventory TSV/TSV.GZ file.")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT, help="Absolute path to the timestamped output root directory.")
    parser.add_argument("--run-id", default=None, help="Optional timestamped directory name in YYYYMMDD_HHMMSS format.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    overlap_path = args.overlap.expanduser()
    inventory_path = args.candidate_inventory.expanduser()
    output_root = args.output_root.expanduser()
    for path in (overlap_path, inventory_path, output_root):
        if not path.is_absolute():
            print(f"error: path must be absolute: {path}", file=sys.stderr)
            return 2
    for path in (overlap_path, inventory_path):
        if not path.is_file():
            print(f"error: input file does not exist: {path}", file=sys.stderr)
            return 2
    run_id = args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    if not re.fullmatch(r"\d{8}_\d{6}", run_id):
        print(f"error: --run-id must be such as YYYYMMDD_HHMMSS: {run_id}", file=sys.stderr)
        return 2
    run_dir = output_root / run_id
    if run_dir.exists():
        print(f"error: output directory already exists; exiting to avoid overwriting it: {run_dir}", file=sys.stderr)
        return 2
    run_dir.mkdir(parents=True, exist_ok=False)

    with RunLogger(run_dir / RUN_LOG_NAME) as log:
        log.log("task started / start; concurrent_tasks=1; threads_per_task=1(retain serial)")
        log.log(f"[stage 1/5] read fine-map overlap: {overlap_path}")
        overlap_fields, overlap_rows = read_tsv(overlap_path)
        missing_overlap = sorted(REQUIRED_OVERLAP_FIELDS.difference(overlap_fields))
        if missing_overlap:
            log.log("overlap input field missing: " + ", ".join(missing_overlap))
            return 2
        log.log(f"overlap column header check passed: {len(overlap_fields)} column; read {len(overlap_rows)} row")

        log.log(f"[stage 2/5] read candidate inventory: {inventory_path}")
        inventory_fields, inventory_rows = read_tsv(inventory_path)
        missing_inventory = sorted(REQUIRED_INVENTORY_FIELDS.difference(inventory_fields))
        if missing_inventory:
            log.log("candidate inventory field missing: " + ", ".join(missing_inventory))
            return 2
        try:
            variant_index, minimum_index, credible_set_rank_index = build_inventory_index(inventory_rows)
        except ValueError as exc:
            log.log(f"candidate inventory index failed: {exc}")
            return 2
        log.log(f"candidate inventory core credible set index: {len(variant_index)} items SNP; {len(minimum_index)} items credible set")

        log.log("[stage 3/5] filter main credible set of annotated_ncRNA--pcGene relationship and pad credible set minimum p_lrt and internal PIP ranking")
        output_input = retained_input_fields(overlap_fields)
        rows: list[dict[str, str]] = []
        skipped = defaultdict(int)
        for row in overlap_rows:
            if row.get("gwas_branch") != BRANCH:
                skipped["non_main_branch"] += 1
                continue
            if row.get("candidate_type") != "susie_credible_set":
                skipped["non_credible_candidate"] += 1
                continue
            if row.get("gwas_in_credible_set") != "yes":
                skipped["not_in_credible_set"] += 1
                continue
            if not standard_pair(row):
                skipped["not_annotated_ncrna_pcgene"] += 1
                continue
            try:
                rows.append(
                    oriented_enriched_row(
                        row,
                        output_input,
                        variant_index,
                        minimum_index,
                        credible_set_rank_index,
                    )
                )
            except ValueError as exc:
                raise ValueError(f"input association or orient error, do not silently discard: {exc}") from exc
        focus_rows = [row for row in rows if row.get("gwas_hit_on_ncrna") == "yes"]
        rows.sort(key=lambda row: (as_float(row.get("gwas_variant_p_lrt", "")) or float("inf"), -(as_float(row.get("gwas_pip", "")) or 0.0), row.get("gwas_variant_display", "")))
        focus_rows.sort(key=lambda row: (as_float(row.get("gwas_variant_p_lrt", "")) or float("inf"), -(as_float(row.get("gwas_pip", "")) or 0.0), row.get("gwas_variant_display", "")))
        log.log(f"retain {len(rows)} regular records ncRNA–pcGene credible overlap expanded rows; focus={len(focus_rows)}; skip={dict(skipped)}")

        log.log("[stage 4/5] write complete evidence, relationship summary and SVG/Excel")
        pair_rows = build_pair_summary(rows)
        summary_rows = build_summary(rows, focus_rows, inventory_rows, len(overlap_rows), run_dir)
        all_fields = output_input + DERIVED_FIELDS
        focus_fields = [field for field in FOCUS_FIELDS if field in all_fields]
        counts = {"overlap_rows": len(overlap_rows), "inventory_rows": len(inventory_rows), "all_rows": len(rows), "focus_rows": len(focus_rows), "unique_pairs": len(pair_rows), "skipped_rows": sum(skipped.values())}
        write_manifest(run_dir / "analysis_manifest.tsv", overlap_path, inventory_path, run_dir, overlap_fields, output_input, counts)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_all.tsv", all_fields, rows)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_focus.tsv", all_fields, focus_rows)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_pair_summary.tsv", PAIR_SUMMARY_FIELDS, pair_rows)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_summary.tsv", SUMMARY_FIELDS, summary_rows)
        write_csv(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_focus_table.csv", focus_fields, focus_rows)
        write_html(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_focus_table.html", "Dog5 main GWAS fine-mapping × conventional ncRNA–pcGene Allelome", focus_fields, focus_rows, "main branch; SuSiE 95% credible-set overlap; focus contains records with a GWAS SNP within the ncRNA gene span.")
        write_xlsx(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_tables.xlsx", [("All", all_fields, rows), ("Focus", focus_fields, focus_rows), ("Pair_Summary", PAIR_SUMMARY_FIELDS, pair_rows), ("Summary", SUMMARY_FIELDS, summary_rows)])
        write_xls(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_tables.xls", [("All", all_fields, rows), ("Focus", focus_fields, focus_rows), ("Pair_Summary", PAIR_SUMMARY_FIELDS, pair_rows), ("Summary", SUMMARY_FIELDS, summary_rows)])
        gwas_summary = inventory_summary(inventory_rows)
        minimum_text = min((as_float(row.get("gwas_variant_p_lrt", "")) for row in rows if as_float(row.get("gwas_variant_p_lrt", "")) is not None), default=None)
        matched_text = f"{minimum_text:.6g}" if minimum_text is not None else "NA"
        summary_text = (
            f"main expanded annotated ncRNA–pcGene rows={len(rows)}; focus rows={len(focus_rows)}; "
            f"unique pairs={len(pair_rows)}; all-main credible sets={gwas_summary['credible_set_count']}; "
            f"credible-set SNPs={gwas_summary['unique_credible_set_snps']}; "
            f"p_lrt<1e-5 SNPs={gwas_summary['p_lrt_lt_1e-5_unique_snps']}; "
            f"minimum matched SNP p_lrt={matched_text}"
        )
        write_svg_table(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_focus_table.svg", "Five-dog main Allelome × GWAS fine-mapping: conventional ncRNA–pcGene", focus_rows, len(focus_rows), summary_text)
        write_svg_pair_barplot(run_dir / "dog5_annotated_ncrna_pcgene_gwas_finemap_pair_barplot.svg", pair_rows)

        log.log("[stage 5/5] completed write and record execution status(not a result QC)")
        log.log("GWAS fine-map core field retain; add within-credible-set PIP ranking; newgene auxiliary columns are absent from the output; SVG use summary style of group column header.")
        log.log(f"output directory: {run_dir}")
        (run_dir / "execution_completion.marker").write_text(
            "execution_status=success\n",
            encoding="utf-8",
        )
        log.log("task completed / done")
        log.log("execution_status=success")
        log.log("reminder: Successful execution does not imply that output completeness, format, counts, or biological quality meet requirements; review the outputs separately.")
    return 0


if __name__ == "__main__":

    log_dir = Path("/path/to/repository/scripts/12_evidence_integration/5.allelome-compare-gwas-finemap/ncrna-mrna/scripts/log")
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    with (log_dir / f"ncrna_finemap.{stamp}.log").open("x", encoding="utf-8") as handle:
        with redirect_stdout(handle), redirect_stderr(handle):
            success = log_dir / f"ncrna_finemap.{stamp}.success.tsv"
            failed = log_dir / f"ncrna_finemap.{stamp}.failed.tsv"
            write_tsv(success, ["task", "exit_code"], [])
            write_tsv(failed, ["task", "exit_code"], [])
            try:
                code = main()
            except BaseException:
                traceback.print_exc()
                write_tsv(failed, ["task", "exit_code"], [{"task": "summary", "exit_code": 1}])
                raise
            write_tsv(success if code == 0 else failed, ["task", "exit_code"], [{"task": "summary", "exit_code": code}])
    raise SystemExit(code)
