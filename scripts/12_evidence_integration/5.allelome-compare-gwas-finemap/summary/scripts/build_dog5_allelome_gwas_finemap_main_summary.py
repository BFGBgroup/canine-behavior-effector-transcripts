#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import sys
import traceback
from contextlib import redirect_stdout, redirect_stderr
import csv
import gzip
import html
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


DEFAULT_DIRECT_LINKS = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap/20260910_180641/main/dog5_allelome_gwas_finemap_direct_new_gene_pcgene.tsv.gz'
)
DEFAULT_CANDIDATE_INVENTORY = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap/20260910_180641/main/dog5_gwas_finemap_candidate_inventory.tsv.gz'
)
DEFAULT_COMPARISON_MANIFEST = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap/20260910_180641/main/dog5_allelome_gwas_finemap_manifest.tsv'
)
DEFAULT_NEW_GENE_ANNOTATION = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/summary'
)
DEFAULT_WORKBOOK_NAME = "5dogs_gwas_finemap_allelome_main_tables.xls"
DEFAULT_FULL_TSV_NAME = "5dogs_gwas_finemap_allelome_main_all.tsv"
DEFAULT_FOCUS_TSV_NAME = "5dogs_gwas_finemap_allelome_main_focus.tsv"
DEFAULT_DISPLAY_CSV_NAME = "5dogs_gwas_finemap_allelome_main_focus_table.csv"
DEFAULT_DISPLAY_HTML_NAME = "5dogs_gwas_finemap_allelome_main_focus_table.html"
DEFAULT_DISPLAY_XLSX_NAME = "5dogs_gwas_finemap_allelome_main_focus_table.xlsx"
DEFAULT_DISPLAY_SVG_NAME = "5dogs_gwas_finemap_allelome_main_focus_table.svg"
DEFAULT_GWAS_MEMBERS_NAME = "dog5_allelome_gwas_finemap_main_gwas_credible_set_members.tsv"
DEFAULT_GWAS_SUMMARY_NAME = "dog5_allelome_gwas_finemap_main_gwas_credible_set_summary.tsv"
DEFAULT_SUMMARY_NAME = "dog5_allelome_gwas_finemap_main_summary.tsv"
DEFAULT_MANIFEST_NAME = "analysis_manifest.tsv"
DEFAULT_LOG_NAME = "dog5_allelome_gwas_finemap_main_summary.log"

BRANCH = "main"
GWAS_MEMBER_REQUIRED_FIELDS = {
    "locus_id",
    "trait",
    "analysis_scope",
    "branch",
    "cs_id",
    "chrom",
    "pos",
    "ld_variant_id",
    "ld_a1",
    "ld_a2",
    "pip",
    "z_for_ld_a1",
    "p_lrt",
}

GWAS_CREDIBLE_SET_MEMBER_COLUMNS = [
    "branch", "analysis_scope", "trait", "locus_id", "cs_index", "cs_id",
    "matrix_variant_index", "chrom", "pos", "ld_variant_id", "ld_a1", "ld_a2",
    "p_lrt", "pip", "z_for_ld_a1", "p_lrt_lt_1e-5", "credible_snp_key",
]

GWAS_CREDIBLE_SUMMARY_COLUMNS = [
    "scope", "gwas_locus_count", "credible_set_count",
    "credible_set_membership_rows", "unique_credible_set_snps",
    "duplicate_membership_rows", "p_lrt_lt_1e-5_membership_rows",
    "p_lrt_lt_1e-5_unique_snps", "min_gwas_p_lrt", "min_gwas_p_lrt_snp",
    "min_gwas_p_lrt_snp_coordinate", "min_gwas_p_lrt_chrom",
    "min_gwas_p_lrt_pos", "min_gwas_p_lrt_trait", "min_gwas_p_lrt_locus",
    "min_gwas_p_lrt_pip", "min_gwas_p_lrt_z_for_ld_a1",
]

SUMMARY_COLUMNS = ["metric", "value"]
MANIFEST_COLUMNS = ["role", "path", "exists", "size_bytes", "rows", "schema_status"]

STANDARD_GENE_COORDINATE_SOURCES = {
    "gtf_unique_gene_id",
    "gtf_unique_gene_id_non_primary_contig",
}


def log(message: str, log_handle=None) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ] {message}"
    print(line, flush=True)
    if log_handle is not None and log_handle is not sys.stdout:
        log_handle.write(line + "\n")
        log_handle.flush()


def read_tsv_rows(path: Path) -> List[Dict[str, str]]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_item_value_table(path: Path) -> Dict[str, str]:
    rows = read_tsv_rows(path)
    if not rows:
        raise ValueError(f"item/value table is empty: {path}")
    if list(rows[0])[:2] != ["item", "value"]:
        raise ValueError(f"column header is not item/value: {path}")
    return {
        row.get("item", ""): row.get("value", "")
        for row in rows
        if row.get("item", "")
    }


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def pick_first(row: Dict[str, str], candidates: Sequence[str]) -> str:
    for key in candidates:
        value = row.get(key, "")
        if value not in ("", None):
            return str(value)
    return ""


def to_int(text: str) -> int:
    return int(float(str(text).strip()))


def normalize_chr(value: str) -> str:
    text = str(value).strip()
    if not text:
        return ""
    if text.startswith("chr"):
        return text
    return f"chr{text}"


def split_annotation_values(value: str) -> List[str]:
    values: List[str] = []
    seen = set()
    for part in re.split(r"[|;]", str(value)):
        text = part.strip()
        if not text or text in {"NA", "."} or text in seen:
            continue
        values.append(text)
        seen.add(text)
    return values


def build_new_gene_feelnc_index(rows, path):

    index = defaultdict(set)
    for row in rows:
        index[row["final_gene_id"]].add(row["FEELnc_biotype"] or "NA")
    return {key: "|".join(sorted(values)) for key, values in index.items()}


def lookup_new_gene_feelnc_biotypes(
    unified_locus_ids: str,
    feelnc_index: Dict[str, str],
) -> str:
    locus_ids = split_annotation_values(unified_locus_ids)
    missing = [locus_id for locus_id in locus_ids if locus_id not in feelnc_index]
    if missing:
        raise ValueError(
            "unified novel gene annotation table missing locus: "
            f"{', '.join(missing)}"
        )
    biotypes: List[str] = []
    for locus_id in locus_ids:
        for biotype in split_annotation_values(feelnc_index[locus_id]):
            if biotype not in biotypes:
                biotypes.append(biotype)
    return "|".join(biotypes) if biotypes else "NA"


INVENTORY_REQUIRED_FIELDS = {
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
    "gwas_variant_p_lrt",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_variant_count",
    "gwas_credible_sets",
    "gwas_in_credible_set",
    "gwas_harmonization_status",
    "gwas_locus_lead_variant_id",
    "gwas_locus_lead_pos",
    "gwas_locus_lead_p_lrt",
    "gwas_locus_top_variant_id",
    "gwas_locus_top_pos",
    "gwas_locus_max_pip",
}


def build_candidate_inventory_index(
    rows: Sequence[Dict[str, str]],
    path: Path,
) -> Dict[Tuple[str, str, str, str], Dict[str, str]]:
    if not rows:
        raise ValueError(f"GWAS fine-map candidate inventory is empty: {path}")
    missing = sorted(INVENTORY_REQUIRED_FIELDS - set(rows[0]))
    if missing:
        raise ValueError(
            "GWAS fine-map candidate inventory missing required fields: "
            f"{', '.join(missing)}"
        )

    index: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
    for row in rows:
        key = (
            row.get("gwas_analysis_scope", ""),
            row.get("gwas_branch", ""),
            row.get("locus_id", ""),
            row.get("gwas_ld_variant_id", ""),
        )
        if not all(key):
            raise ValueError(
                "candidate inventory contains empty scope, branch, locus_id or variant_id"
            )
        if key in index:
            raise ValueError(
                "candidate inventory contains duplicate scope+branch+locus_id+variant_id: "
                f"{' / '.join(key)}"
            )
        index[key] = row
    return index


def read_tsv_schema_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:

    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open("r", encoding="utf-8", newline="")
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def format_gwas_coordinate(row: Mapping[str, str]) -> str:

    chrom = str(row.get("chrom", row.get("gwas_chrom", ""))).strip()
    pos = str(row.get("pos", row.get("gwas_pos", ""))).strip()
    a1 = str(row.get("ld_a1", row.get("gwas_allele1", ""))).strip().upper()
    a2 = str(row.get("ld_a2", row.get("gwas_allele0", ""))).strip().upper()
    if not chrom or not pos or not a1 or not a2:
        return "NA"
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    alleles = sorted((a1, a2))
    return f"{chrom}:{pos}:{alleles[0]}/{alleles[1]}"


def read_gwas_credible_members(
    gwas_run_dir: Path,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:

    manifest_path = gwas_run_dir / "locus_run_manifest.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"GWAS locus manifest does not exist: {manifest_path}")
    manifest_fields, manifest_rows = read_tsv_schema_rows(manifest_path)
    required_manifest = {
        "locus_id", "trait", "analysis_scope", "branch", "status",
        "model_status", "credible_set_count", "credible_set_members",
    }
    missing_manifest = sorted(required_manifest - set(manifest_fields))
    if missing_manifest:
        raise ValueError(
            f"{manifest_path} missing credible set field: {', '.join(missing_manifest)}"
        )
    main_manifest = [row for row in manifest_rows if row.get("branch") == BRANCH]
    if not main_manifest:
        raise ValueError(f"{manifest_path} contains no main branch locus")
    invalid_main = [
        row.get("locus_id", "")
        for row in main_manifest
        if row.get("status") != "PASS" or row.get("model_status") != "CONVERGED"
    ]
    if invalid_main:
        raise ValueError(
            "GWAS main manifest contains non-PASS/CONVERGED locus: "
            + ", ".join(invalid_main[:5])
        )

    members: List[Dict[str, str]] = []
    for manifest_row in main_manifest:
        member_path = Path(manifest_row.get("credible_set_members", ""))
        if not member_path.is_absolute() or not member_path.is_file():
            raise FileNotFoundError(
                f"GWAS credible_set_members does not exist or is not an absolute path: {member_path}"
            )
        member_fields, member_rows = read_tsv_schema_rows(member_path)
        missing_member = sorted(GWAS_MEMBER_REQUIRED_FIELDS - set(member_fields))
        if missing_member:
            raise ValueError(
                f"{member_path} missing credible set member field: {', '.join(missing_member)}"
            )
        for member_row in member_rows:
            for key in ("locus_id", "trait", "analysis_scope", "branch"):
                if member_row.get(key, "") != manifest_row.get(key, ""):
                    raise ValueError(
                        f"{member_path} of {key} with manifest mismatch"
                    )
            try:
                p_lrt = float(member_row.get("p_lrt", ""))
            except ValueError as exc:
                raise ValueError(f"{member_path} contains non-numeric value p_lrt") from exc
            if not math.isfinite(p_lrt):
                raise ValueError(f"{member_path} contains non-finite p_lrt")
            chrom = member_row.get("chrom", "")
            pos = member_row.get("pos", "")
            variant_id = member_row.get("ld_variant_id", "")
            a1 = member_row.get("ld_a1", "").strip().upper()
            a2 = member_row.get("ld_a2", "").strip().upper()
            if not chrom or not pos or not variant_id or not a1 or not a2:
                raise ValueError(
                    f"{member_path} contains empty chrom/pos/ld_variant_id/ld_a1/ld_a2"
                )
            members.append(
                {
                    "branch": BRANCH,
                    "analysis_scope": manifest_row.get("analysis_scope", ""),
                    "trait": manifest_row.get("trait", ""),
                    "locus_id": manifest_row.get("locus_id", ""),
                    "cs_index": member_row.get("cs_index", ""),
                    "cs_id": member_row.get("cs_id", ""),
                    "matrix_variant_index": member_row.get("matrix_variant_index", ""),
                    "chrom": chrom,
                    "pos": pos,
                    "ld_variant_id": variant_id,
                    "ld_a1": a1,
                    "ld_a2": a2,
                    "p_lrt": member_row.get("p_lrt", ""),
                    "pip": member_row.get("pip", ""),
                    "z_for_ld_a1": member_row.get("z_for_ld_a1", ""),
                    "p_lrt_lt_1e-5": "yes" if p_lrt < 1e-5 else "no",
                    "credible_snp_key": f"chr{str(chrom).removeprefix('chr')}:{pos}:{variant_id}",
                }
            )
    return main_manifest, members


def summarize_gwas_credible_members(
    manifest_rows: Sequence[Mapping[str, str]],
    members: Sequence[Mapping[str, str]],
) -> List[Dict[str, object]]:

    output: List[Dict[str, object]] = []
    for scope in ("all_main", "primary", "exploratory"):
        if scope == "all_main":
            selected_manifest = list(manifest_rows)
            selected_members = list(members)
        else:
            selected_manifest = [
                row for row in manifest_rows if row.get("analysis_scope") == scope
            ]
            selected_members = [
                row for row in members if row.get("analysis_scope") == scope
            ]
        unique_keys = {
            row.get("credible_snp_key", "")
            for row in selected_members
            if row.get("credible_snp_key", "")
        }
        low_p_members = [
            row for row in selected_members if row.get("p_lrt_lt_1e-5") == "yes"
        ]
        low_p_keys = {
            row.get("credible_snp_key", "")
            for row in low_p_members
            if row.get("credible_snp_key", "")
        }
        min_row = min(
            selected_members,
            key=lambda row: float(row.get("p_lrt", "inf")),
            default=None,
        )
        min_values: Dict[str, object]
        if min_row is None:
            min_values = {
                "min_gwas_p_lrt": "NA", "min_gwas_p_lrt_snp": "NA",
                "min_gwas_p_lrt_snp_coordinate": "NA", "min_gwas_p_lrt_chrom": "NA",
                "min_gwas_p_lrt_pos": "NA", "min_gwas_p_lrt_trait": "NA",
                "min_gwas_p_lrt_locus": "NA", "min_gwas_p_lrt_pip": "NA",
                "min_gwas_p_lrt_z_for_ld_a1": "NA",
            }
        else:
            min_values = {
                "min_gwas_p_lrt": min_row.get("p_lrt", "NA"),
                "min_gwas_p_lrt_snp": min_row.get("ld_variant_id", "NA"),
                "min_gwas_p_lrt_snp_coordinate": format_gwas_coordinate(min_row),
                "min_gwas_p_lrt_chrom": min_row.get("chrom", "NA"),
                "min_gwas_p_lrt_pos": min_row.get("pos", "NA"),
                "min_gwas_p_lrt_trait": min_row.get("trait", "NA"),
                "min_gwas_p_lrt_locus": min_row.get("locus_id", "NA"),
                "min_gwas_p_lrt_pip": min_row.get("pip", "NA"),
                "min_gwas_p_lrt_z_for_ld_a1": min_row.get("z_for_ld_a1", "NA"),
            }
        output.append(
            {
                "scope": scope,
                "gwas_locus_count": len(selected_manifest),
                "credible_set_count": sum(
                    int(row.get("credible_set_count", "0"))
                    for row in selected_manifest
                ),
                "credible_set_membership_rows": len(selected_members),
                "unique_credible_set_snps": len(unique_keys),
                "duplicate_membership_rows": len(selected_members) - len(unique_keys),
                "p_lrt_lt_1e-5_membership_rows": len(low_p_members),
                "p_lrt_lt_1e-5_unique_snps": len(low_p_keys),
                **min_values,
            }
        )
    return output


def enrich_credible_set_context(rows, members):

    sets = defaultdict(dict)
    for member in members:
        sets[(member["locus_id"], member["cs_id"])][member["ld_variant_id"]] = member
    result = []
    for row in rows:
        for cs in row["gwas_credible_sets"].replace(",", ";").split(";"):
            if not cs:
                continue
            pool = sets[(row["gwas_locus_id"], cs)]
            match = pool[row["gwas_snp_id"]]
            minimum = min(pool.values(), key=lambda m: (float(m["p_lrt"]), m["ld_variant_id"]))
            rank = 1 + sum(float(m["pip"]) > float(match["pip"]) for m in pool.values())
            record = dict(row)
            record.update({
                "gwas_credible_sets": cs,
                "gwas_variant_coordinate": format_gwas_coordinate(match),
                "gwas_variant_display": f'{match["ld_variant_id"]} | {format_gwas_coordinate(match)}',
                "gwas_pip_rank_display": f'{row["gwas_pip_rank"]}/{row["gwas_variant_count"]}',
                "gwas_cs_pip_rank": str(rank),
                "gwas_cs_size": str(len(pool)),
                "gwas_cs_pip_rank_display": f'{cs}: {rank}/{len(pool)}',
                "gwas_credible_set_min_p_snp": minimum["ld_variant_id"],
                "gwas_credible_set_min_p_snp_coordinate": format_gwas_coordinate(minimum),
                "gwas_credible_set_min_p_snp_display": f'{minimum["ld_variant_id"]} | {format_gwas_coordinate(minimum)}',
                "gwas_credible_set_min_p_lrt": minimum["p_lrt"],
            })
            result.append(record)
    return result


def build_current_enriched_rows(input_rows, feelnc_index, candidate_index):

    output = []
    for row in input_rows:
        side = row["endpoint_hit"]
        other = "target" if side == "base" else "base"
        if (side not in {"base", "target"} or row["gwas_branch"] != "main"
                or row["gwas_in_credible_set"] != "yes"
                or row[side + "_record_type"] != "manual_new_gene"
                or row[other + "_gene_biotype"] != "protein_coding"):
            raise ValueError("input is not main credible set direct new gene—pcGene link")
        gene = row[side + "_gene_id"]
        key = (row["gwas_analysis_scope"], row["gwas_branch"], row["locus_id"], row["gwas_ld_variant_id"])
        candidate = candidate_index[key]
        record = dict(row)

        record["gwas_pip_rank"] = candidate["gwas_pip_rank"]
        record["gwas_variant_count"] = candidate["gwas_variant_count"]
        record.update({
            "gwas_trait": row["trait"], "gwas_locus_id": row["locus_id"],
            "gwas_snp_id": row["gwas_ld_variant_id"], "gwas_variant": row["gwas_ld_variant_id"],
            "gwas_chr": row["gwas_chrom"], "gwas_allele1": row["gwas_ld_a1"],
            "gwas_allele0": row["gwas_ld_a2"],
            "manual_new_gene": gene,
            "new_gene_span": f'{row[side+"_chr"]}:{int(row[side+"_start"])+1}–{row[side+"_end"]}',
            "new_gene_strand": row[side + "_gene_strand"],
            "protein_coding_gene": row[other + "_gene_symbol"] or row[other + "_gene_id"],
            "protein_coding_gene_id": row[other + "_gene_id"],
            "new_gene_class_code": row[side + "_manual_new_gene_class_codes"],
            "new_gene_feelnc_biotype": feelnc_index[gene],
            "sample": row["animal_id"],
        })
        output.append(record)
    return output


def chrom_sort_key(chr_name: str) -> Tuple[int, str]:
    text = chr_name.replace("chr", "")
    if text.isdigit():
        return (0, f"{int(text):04d}")
    return (1, text)


def row_sort_key(row: Dict[str, str]) -> Tuple[Tuple[int, str], int, str, str]:
    pos = to_int(row.get("gwas_pos", "0") or "0")
    return (
        chrom_sort_key(row.get("gwas_chr", "")),
        pos,
        row.get("gwas_trait", ""),
        row.get("sample_tissue", ""),
    )


FOCUS_COLUMNS = [
    "gwas_trait", "gwas_analysis_scope", "gwas_branch", "gwas_locus_id",
    "gwas_variant_display", "gwas_variant_p_lrt", "gwas_pip",
    "gwas_pip_rank_display", "gwas_credible_sets", "gwas_cs_pip_rank",
    "gwas_cs_size", "gwas_cs_pip_rank_display",
    "gwas_credible_set_min_p_snp_display", "gwas_credible_set_min_p_lrt",
    "manual_new_gene", "new_gene_span", "new_gene_strand",
    "new_gene_class_code", "new_gene_feelnc_biotype",
    "protein_coding_gene", "protein_coding_gene_id", "sample_tissue",
    "mechanism", "linkage_score", "evidence_stratum",
    "overlap_bp", "either_endpoint_multigene", "shared_variant_identity_checked",
    "candidate_key", "sample_pair_id", "canonical_pair_id",
]
FULL_COLUMNS = FOCUS_COLUMNS


DISPLAY_COLUMNS = [
    "gwas_trait",
    "gwas_variant_display",
    "gwas_variant_p_lrt",
    "gwas_pip",
    "gwas_pip_rank_display",
    "gwas_credible_sets",
    "gwas_credible_set_min_p_snp_display",
    "gwas_credible_set_min_p_lrt",
    "evidence_stratum",
    "new_gene_strand",
    "new_gene_span",
    "gwas_cs_pip_rank_display",
    "protein_coding_gene",
    "mechanism",
    "linkage_score",
    "sample_tissue",
    "manual_new_gene",
    "new_gene_class_code",
    "new_gene_feelnc_biotype",
]


DISPLAY_LABELS = {
    "gwas_trait": "GWAS Trait",
    "gwas_variant_display": "Matched GWAS SNP (ID | chr:pos:alleles)",
    "gwas_variant_p_lrt": "Matched SNP P LRT",
    "gwas_pip": "Matched SNP PIP",
    "gwas_pip_rank_display": "Locus PIP Rank (rank/total)",
    "gwas_credible_sets": "Credible Set",
    "gwas_credible_set_min_p_snp_display": "Credible-set Min-p SNP (ID | chr:pos:alleles)",
    "gwas_credible_set_min_p_lrt": "Credible-set Min P LRT",
    "evidence_stratum": "Count Evidence Stratum",
    "new_gene_strand": "New Gene Strand",
    "new_gene_span": "New Gene Span (1-based)",
    "gwas_cs_pip_rank_display": "CS PIP Rank (rank/size)",
    "protein_coding_gene": "Paired Protein-Coding Gene",
    "mechanism": "Allelome Pattern",
    "linkage_score": "Allelome Score",
    "sample_tissue": "Sample Tissue",
    "manual_new_gene": "Matched Manual New Gene",
    "new_gene_class_code": "New Gene Class Code",
    "new_gene_feelnc_biotype": "New Gene FEELnc Biotype",
}

DISPLAY_COLUMN_WIDTHS = {
    "gwas_trait": 24,
    "gwas_variant_display": 36,
    "gwas_variant_p_lrt": 20,
    "gwas_pip": 20,
    "gwas_pip_rank_display": 18,
    "gwas_credible_sets": 14,
    "gwas_credible_set_min_p_snp_display": 38,
    "gwas_credible_set_min_p_lrt": 20,
    "evidence_stratum": 22,
    "new_gene_strand": 22,
    "new_gene_span": 28,
    "gwas_cs_pip_rank_display": 22,
    "protein_coding_gene": 28,
    "mechanism": 15,
    "linkage_score": 16,
    "sample_tissue": 28,
    "manual_new_gene": 30,
    "new_gene_class_code": 20,
    "new_gene_feelnc_biotype": 24,
}


def sanitize_sheet_name(name: str) -> str:
    invalid = '[]:*?/\\'
    cleaned = "".join("_" if ch in invalid else ch for ch in name)
    return cleaned[:31] or "Sheet1"


def xml_cell(value: str, style_id: str = "Body") -> str:
    text = escape("" if value is None else str(value))
    return (
        f'<Cell ss:StyleID="{style_id}">'
        f'<Data ss:Type="String">{text}</Data>'
        f"</Cell>"
    )


def xml_row(values: Sequence[str], style_id: str = "Body") -> str:
    return "<Row>" + "".join(xml_cell(value, style_id=style_id) for value in values) + "</Row>"


def write_excel_workbook(
    path: Path,
    *,
    all_rows: Sequence[Dict[str, str]],
    focus_rows: Sequence[Dict[str, str]],
) -> None:
    note_lines = [
        "Full_Results: all expanded direct-new-gene–pcGene rows from upstream direct-newgene results.",
        "Nature_Focus: the same records as Full_Results, without aggregation or removal of repeated biological relationships.",
        "This workbook summarizes only the main branch: GWAS credible-set SNPs overlap the manual_new_gene endpoint but not the paired protein_coding endpoint.",
        "No secondary filtering by member p_lrt; primary and exploratory are defined upstream, and main does not mean primary.",
        "The denominator for Locus PIP Rank is the full locus candidate count; the denominator for CS PIP Rank is the current credible-set size, using competitive PIP ranks.",
        "The lowest credible-set p value is selected from members of the row's locus and CS, not from the global minimum across all loci.",
        "PIP is the posterior inclusion probability of a candidate variant conditional on the current fine-mapping locus and model, not a genome-wide significance P value.",
        "The current evidence does not identify a causal regulator/target; base/target specifies only endpoint order.",
        "mechanism retains the enhancing / repressing patterns from Allelome, without interpreting them as established causal regulatory directions.",
        "A GWAS credible-set variant falls within an aggregated new-gene endpoint interval; this does not establish identity with a particular ASE site.",
        "Class Code retains upstream annotations; FEELnc labels come from the five-dog transcript_map, taking the union across transcripts of the same gene.",
        "New Gene Span displays a 1-based inclusive interval converted from the upstream 0-based half-open interval using start+1 and end.",
        "Multi-gene counting risk flags are retained; shared_variant_identity_checked=no must not be interpreted as absence of shared ASE sites.",
        "The workbook uses Excel 2003 XML .xls format, written with the standard library and directly readable in Excel.",
    ]

    all_headers = list(dict.fromkeys(key for row in all_rows for key in row)) or FULL_COLUMNS
    focus_headers = FOCUS_COLUMNS

    xml_parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<?mso-application progid="Excel.Sheet"?>',
        '<Workbook xmlns="urn:schemas-microsoft-com:office:spreadsheet"'
        ' xmlns:o="urn:schemas-microsoft-com:office:office"'
        ' xmlns:x="urn:schemas-microsoft-com:office:excel"'
        ' xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet"'
        ' xmlns:html="http://www.w3.org/TR/REC-html40">',
        "<Styles>",
        '<Style ss:ID="Default" ss:Name="Normal"><Alignment ss:Vertical="Center" ss:WrapText="1"/>'
        '<Borders/>'
        '<Font ss:FontName="Arial" ss:Size="10"/>'
        '<Interior/>'
        '<NumberFormat/>'
        '<Protection/>'
        "</Style>",
        '<Style ss:ID="Title"><Alignment ss:Horizontal="Center" ss:Vertical="Center" ss:WrapText="1"/>'
        '<Font ss:FontName="Arial" ss:Size="12" ss:Bold="1"/>'
        '<Interior ss:Color="#FFFFFF" ss:Pattern="Solid"/>'
        '<Borders>'
        '<Border ss:Position="Bottom" ss:LineStyle="Continuous" ss:Weight="2"/>'
        "</Borders></Style>",
        '<Style ss:ID="Header"><Alignment ss:Horizontal="Center" ss:Vertical="Center" ss:WrapText="1"/>'
        '<Font ss:FontName="Arial" ss:Size="10" ss:Bold="1"/>'
        '<Interior ss:Color="#E6E6E6" ss:Pattern="Solid"/>'
        '<Borders>'
        '<Border ss:Position="Bottom" ss:LineStyle="Continuous" ss:Weight="1"/>'
        '<Border ss:Position="Left" ss:LineStyle="Continuous" ss:Weight="1"/>'
        '<Border ss:Position="Right" ss:LineStyle="Continuous" ss:Weight="1"/>'
        '<Border ss:Position="Top" ss:LineStyle="Continuous" ss:Weight="1"/>'
        "</Borders></Style>",
        '<Style ss:ID="Body"><Alignment ss:Vertical="Center" ss:WrapText="1"/>'
        '<Font ss:FontName="Arial" ss:Size="10"/>'
        '<Borders>'
        '<Border ss:Position="Bottom" ss:LineStyle="Continuous" ss:Weight="1"/>'
        '<Border ss:Position="Left" ss:LineStyle="Continuous" ss:Weight="1"/>'
        '<Border ss:Position="Right" ss:LineStyle="Continuous" ss:Weight="1"/>'
        '<Border ss:Position="Top" ss:LineStyle="Continuous" ss:Weight="1"/>'
        "</Borders></Style>",
        '<Style ss:ID="Note"><Alignment ss:Vertical="Center" ss:WrapText="1"/>'
        '<Font ss:FontName="Arial" ss:Size="10" ss:Italic="1"/>'
        '<Borders>'
        '<Border ss:Position="Bottom" ss:LineStyle="Continuous" ss:Weight="1"/>'
        "</Borders></Style>",
        "</Styles>",
    ]

    def add_sheet(sheet_name: str, title: str, headers: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
        xml_parts.append(f'<Worksheet ss:Name="{escape(sanitize_sheet_name(sheet_name))}">')
        xml_parts.append("<Table>")
        xml_parts.append(
            f'<Row><Cell ss:StyleID="Title" ss:MergeAcross="{max(len(headers) - 1, 0)}">'
            f'<Data ss:Type="String">{escape(title)}</Data></Cell></Row>'
        )
        xml_parts.append(xml_row(headers, style_id="Header"))
        for row in rows:
            xml_parts.append(xml_row([row.get(header, "") for header in headers], style_id="Body"))
        xml_parts.append("</Table></Worksheet>")

    add_sheet(
        "Full_Results",
        "5-dogs GWAS fine-mapping × Allelome association: full table",
        all_headers,
        all_rows,
    )
    add_sheet(
        "Nature_Focus",
        "5-dogs GWAS fine-mapping × Allelome association: focus table",
        focus_headers,
        focus_rows,
    )

    xml_parts.append(f'<Worksheet ss:Name="{escape(sanitize_sheet_name("Notes"))}"><Table>')
    xml_parts.append(
        '<Row><Cell ss:StyleID="Title" ss:MergeAcross="0">'
        '<Data ss:Type="String">Workbook notes</Data></Cell></Row>'
    )
    for line in note_lines:
        xml_parts.append(xml_row([line], style_id="Note"))
    xml_parts.append("</Table></Worksheet>")
    xml_parts.append("</Workbook>")

    path.write_text("\n".join(xml_parts), encoding="utf-8")


def yes_no_label(value: str) -> str:
    if str(value).strip().lower() == "yes":
        return "Yes"
    if str(value).strip().lower() == "no":
        return "No"
    return value or "NA"


def mechanism_label(value: str) -> str:
    text = str(value).strip()
    if not text:
        return "NA"
    return text.capitalize()


def build_display_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    display_rows: List[Dict[str, str]] = []
    for row in rows:
        display_rows.append(
            {
                "gwas_trait": row.get("gwas_trait", "") or "NA",
                "gwas_variant_display": row.get("gwas_variant_display", "") or "NA",
                "gwas_variant_p_lrt": row.get(
                    "gwas_variant_p_lrt", ""
                )
                or "NA",
                "gwas_pip": row.get("gwas_pip", "") or "NA",
                "gwas_pip_rank_display": row.get("gwas_pip_rank_display", "") or "NA",
                "gwas_credible_sets": row.get(
                    "gwas_credible_sets", ""
                )
                or "NA",
                "gwas_credible_set_min_p_snp_display": row.get(
                    "gwas_credible_set_min_p_snp_display", ""
                ) or "NA",
                "gwas_credible_set_min_p_lrt": row.get(
                    "gwas_credible_set_min_p_lrt", ""
                ) or "NA",
                "evidence_stratum": row.get("evidence_stratum", "") or "NA",
                "new_gene_strand": row.get("new_gene_strand", ""),
                "new_gene_span": row.get(
                    "new_gene_span", ""
                )
                or "NA",
                "gwas_cs_pip_rank_display": row.get(
                    "gwas_cs_pip_rank_display", ""
                )
                or "NA",
                "protein_coding_gene": row.get(
                    "protein_coding_gene", ""
                )
                or "NA",
                "mechanism": mechanism_label(row.get("mechanism", "")),
                "linkage_score": row.get("linkage_score", "") or "NA",
                "sample_tissue": (row.get("sample_tissue", "") or "NA").replace("_", " "),
                "manual_new_gene": row.get("manual_new_gene", ""),
                "new_gene_class_code": row.get(
                    "new_gene_class_code", ""
                )
                or "NA",
                "new_gene_feelnc_biotype": row.get(
                    "new_gene_feelnc_biotype", ""
                )
                or "NA",
            }
        )
    return display_rows


def build_display_summary(all_rows, focus_rows, gwas_summary):
    return {
        "main_overlap_rows": str(len(all_rows)),
        "unique_gwas_variants": str(len({(r["gwas_trait"], r["gwas_snp_id"]) for r in all_rows})),
        "unique_new_genes": str(len({r["manual_new_gene"] for r in all_rows})),
        "gwas_unique_credible_set_snps": str(gwas_summary["unique_credible_set_snps"]),
        "gwas_p_lt_1e-5_unique_snps": str(gwas_summary["p_lrt_lt_1e-5_unique_snps"]),
    }


def write_csv_table(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_html_table(
    path: Path,
    *,
    title: str,
    fieldnames: Sequence[str],
    rows: Sequence[Dict[str, str]],
    summary: Dict[str, str],
) -> None:
    css = """
body{font-family:Arial,sans-serif;margin:24px;background:#f7f7f4;color:#222}
h1{font-size:22px;margin:0 0 8px}
p.meta{margin:0 0 16px;color:#555;font-size:13px}
.table-wrap{max-width:100%;overflow-x:auto;border:1px solid #cfcfcf;background:#fff}
table{border-collapse:collapse;min-width:2350px;width:100%;font-size:12px;background:#fff}
th,td{border:1px solid #cfcfcf;padding:6px 8px;text-align:left;vertical-align:top}
th{background:#e6efe8;font-weight:700;position:sticky;top:0;z-index:1}
tr.group-head th{background:#dce7df;text-align:center;top:0}
thead tr:nth-child(2) th{top:29px}
tr:nth-child(even){background:#fafafa}
"""
    lines = [
        "<!DOCTYPE html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='utf-8' />",
        f"<title>{html.escape(title)}</title>",
        f"<style>{css}</style>",
        "</head>",
        "<body>",
        f"<h1>{html.escape(title)}</h1>",
        (
            "<p class='meta'>"
            f"main expanded rows = {html.escape(summary.get('main_overlap_rows', 'NA'))}; "
            f"unique GWAS variants = {html.escape(summary.get('unique_gwas_variants', 'NA'))}; "
            f"all-main credible-set SNPs = {html.escape(summary.get('gwas_unique_credible_set_snps', 'NA'))}; "
            f"all-main p_lrt&lt;1e-5 SNPs = {html.escape(summary.get('gwas_p_lt_1e-5_unique_snps', 'NA'))}; "
            "Rows are retained without biological-relation aggregation."
            "</p>"
        ),
        "<div class='table-wrap'><table>",
        (
            "<thead><tr class='group-head'>"
            "<th colspan='8'>GWAS fine-mapping evidence</th>"
            "<th colspan='8'>Allelome association evidence</th>"
            "<th colspan='3'>Manual new-gene evidence</th>"
            "</tr><tr>"
        ),
    ]
    for header in fieldnames:
        lines.append(f"<th>{html.escape(DISPLAY_LABELS.get(header, header))}</th>")
    lines.append("</tr></thead><tbody>")
    for row in rows:
        lines.append("<tr>")
        for header in fieldnames:
            lines.append(f"<td>{html.escape(row.get(header, ''))}</td>")
        lines.append("</tr>")
    lines.extend(["</tbody></table></div>", "</body>", "</html>"])
    path.write_text("\n".join(lines), encoding="utf-8")


def excel_col_name(index: int) -> str:
    name = ""
    current = index
    while current >= 0:
        current, remainder = divmod(current, 26)
        name = chr(ord("A") + remainder) + name
        current -= 1
    return name


def xlsx_inline_cell(ref: str, value: str) -> str:
    return f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'


def write_xlsx_table(
    path: Path,
    *,
    title: str,
    fieldnames: Sequence[str],
    rows: Sequence[Dict[str, str]],
) -> None:
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    header_row = [DISPLAY_LABELS.get(field, field) for field in fieldnames]
    data_rows = [header_row] + [[row.get(field, "") for field in fieldnames] for row in rows]
    sheet_rows = []
    for row_idx, row in enumerate(data_rows, start=1):
        cells = []
        for col_idx, value in enumerate(row):
            ref = f"{excel_col_name(col_idx)}{row_idx}"
            cells.append(xlsx_inline_cell(ref, str(value)))
        sheet_rows.append(f'<row r="{row_idx}">{"".join(cells)}</row>')

    max_col = excel_col_name(len(fieldnames) - 1)
    column_xml = "".join(
        (
            f'<col min="{idx}" max="{idx}" '
            f'width="{DISPLAY_COLUMN_WIDTHS.get(field, 24)}" customWidth="1"/>'
        )
        for idx, field in enumerate(fieldnames, start=1)
    )
    worksheet = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <dimension ref="A1:{max_col}{len(data_rows)}"/>
  <sheetViews><sheetView workbookViewId="0"/></sheetViews>
  <sheetFormatPr defaultRowHeight="15"/>
  <cols>
    {column_xml}
  </cols>
  <sheetData>
    {''.join(sheet_rows)}
  </sheetData>
</worksheet>
"""
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"""
    workbook = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="focus_overlap" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>
"""
    workbook_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>
"""
    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="1"><fill><patternFill patternType="none"/></fill></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf/></cellStyleXfs>
  <cellXfs count="1"><xf xfId="0"/></cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>
"""
    core = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>{escape(title)}</dc:title>
  <dc:creator>Codex</dc:creator>
  <cp:lastModifiedBy>Codex</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>
</cp:coreProperties>
"""
    app = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Codex</Application>
</Properties>
"""
    with ZipFile(path, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("docProps/core.xml", core)
        zf.writestr("docProps/app.xml", app)
        zf.writestr("xl/workbook.xml", workbook)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        zf.writestr("xl/styles.xml", styles)
        zf.writestr("xl/worksheets/sheet1.xml", worksheet)


def wrap_text(value: str, width: int) -> List[str]:
    pieces = re.split(r"([,; ])", value)
    lines: List[str] = []
    current = ""
    for piece in pieces:
        if len(current) + len(piece) <= width:
            current += piece
        else:
            if current:
                lines.append(current.strip())
            current = piece.strip()
    if current:
        lines.append(current.strip())
    return lines or [value]


def wrap_svg_value(column: str, value: str, width: int) -> List[str]:
    text = str(value).strip() or "NA"
    if column == "manual_new_gene":
        parts = [part.strip() for part in text.split(";") if part.strip()]
        return parts or [text]
    return wrap_text(text, width)


def svg_badge_style(column: str, value: str) -> Tuple[str, str]:
    if column in {
        "new_gene_strand",
        "gwas_is_locus_top_variant",
    }:
        if value == "Yes":
            return "#e7f3ec", "#1f5a3a"
        return "#f1f3f5", "#5c6773"
    if column == "evidence_stratum":
        if value == "nonoverlap_no_multigene_flag":
            return "#e7f3ec", "#1f5a3a"
        return "#eef2f7", "#34495e"
    if column == "mechanism":
        if value == "Enhancing":
            return "#e7f3ec", "#1f5a3a"
        if value == "Repressing":
            return "#fff1e6", "#9a4d00"
    return "#f1f3f5", "#5c6773"


def write_svg_table(
    path: Path,
    *,
    title: str,
    rows: Sequence[Dict[str, str]],
    summary: Dict[str, str],
) -> None:
    column_specs = [
        ("gwas_trait", 170, 20, "text_left"),
        ("gwas_variant_display", 285, 44, "text_left"),
        ("gwas_variant_p_lrt", 110, 14, "numeric"),
        ("gwas_pip", 110, 14, "numeric"),
        ("gwas_pip_rank_display", 100, 14, "numeric"),
        ("gwas_credible_sets", 80, 10, "text_left"),
        ("gwas_credible_set_min_p_snp_display", 285, 44, "text_left"),
        ("gwas_credible_set_min_p_lrt", 110, 14, "numeric"),
        ("evidence_stratum", 205, 25, "text_left"),
        ("new_gene_strand", 65, 8, "text_left"),
        ("new_gene_span", 175, 21, "text_left"),
        ("gwas_cs_pip_rank_display", 110, 12, "numeric"),
        ("protein_coding_gene", 155, 19, "text_left"),
        ("mechanism", 110, 15, "badge"),
        ("linkage_score", 90, 11, "numeric"),
        ("sample_tissue", 185, 22, "text_left"),
        ("manual_new_gene", 220, 25, "text_left"),
        ("new_gene_class_code", 120, 16, "text_left"),
        ("new_gene_feelnc_biotype", 155, 20, "text_left"),
    ]
    group_specs = [
        ("GWAS fine-mapping evidence", 0, 8, "#e8f0ea"),
        ("Allelome association evidence", 8, 16, "#eef2f7"),
        ("Manual new-gene evidence", 16, 19, "#f7efe7"),
    ]
    left_margin = 44
    top_margin = 36
    title_height = 70
    group_header_height = 24
    header_height = 38
    line_height = 11
    row_padding = 10
    width = left_margin * 2 + sum(spec[1] for spec in column_specs)
    row_line_map: List[List[List[str]]] = []
    row_heights: List[int] = []
    for row in rows:
        wrapped_cells: List[List[str]] = []
        max_lines = 1
        for key, _col_width, wrap_width, align in column_specs:
            if align == "badge":
                lines = [row.get(key, "NA")]
            else:
                lines = wrap_svg_value(key, row.get(key, "NA"), wrap_width)
            wrapped_cells.append(lines)
            max_lines = max(max_lines, len(lines))
        row_line_map.append(wrapped_cells)
        row_heights.append(max(header_height, row_padding + line_height * max_lines))

    total_row_height = group_header_height + header_height + sum(row_heights)
    height = top_margin + title_height + total_row_height + 22

    x_positions = []
    current_x = left_margin
    for _key, col_width, _wrap_width, _align in column_specs:
        x_positions.append(current_x)
        current_x += col_width

    subtitle = escape(
        "Main 95% credible-set SNPs in candidate new-gene spans; "
        "paired pcGene not hit; positional association, not exact ASE SNP evidence."
    )
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{width / 2:.1f}" y="{top_margin}" font-size="16" font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">{escape(title)}</text>',
        (
            f'<text x="{width / 2:.1f}" y="{top_margin + 18}" font-size="11.5" font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">'
            f"{subtitle}"
            "</text>"
        ),
        (
            f'<text x="{width / 2:.1f}" y="{top_margin + 36}" font-size="11.5" font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">'
            f"main expanded rows={escape(summary.get('main_overlap_rows', 'NA'))}; "
            f"unique GWAS variants={escape(summary.get('unique_gwas_variants', 'NA'))}; "
            f"all-main credible-set SNPs={escape(summary.get('gwas_unique_credible_set_snps', 'NA'))}; "
            f"p_lrt&lt;1e-5 SNPs={escape(summary.get('gwas_p_lt_1e-5_unique_snps', 'NA'))}; "
            "</text>"
        ),
    ]

    header_y = top_margin + title_height
    table_width = sum(spec[1] for spec in column_specs)
    table_bottom = height - 22
    svg.append(f'<line x1="{left_margin}" y1="{header_y}" x2="{left_margin + table_width}" y2="{header_y}" stroke="#111111" stroke-width="1.1"/>')

    for group_label, start_idx, end_idx, fill in group_specs:
        group_x = x_positions[start_idx]
        group_width = sum(column_specs[idx][1] for idx in range(start_idx, end_idx))
        svg.append(
            f'<rect x="{group_x}" y="{header_y}" width="{group_width}" '
            f'height="{group_header_height}" fill="{fill}"/>'
        )
        svg.append(
            f'<text x="{group_x + group_width / 2:.1f}" y="{header_y + 16}" '
            'font-size="10.5" font-family="Helvetica, Arial, sans-serif" '
            'font-weight="700" fill="#111111" text-anchor="middle">'
            f"{escape(group_label)}</text>"
        )

    column_header_y = header_y + group_header_height
    svg.append(
        f'<line x1="{left_margin}" y1="{column_header_y}" '
        f'x2="{left_margin + table_width}" y2="{column_header_y}" '
        'stroke="#aaaaaa" stroke-width="0.7"/>'
    )
    svg.append(
        f'<line x1="{left_margin}" y1="{column_header_y + header_height}" '
        f'x2="{left_margin + table_width}" y2="{column_header_y + header_height}" '
        'stroke="#111111" stroke-width="0.8"/>'
    )

    for idx, (key, col_width, wrap_width, align) in enumerate(column_specs):
        x = x_positions[idx]
        label = DISPLAY_LABELS.get(key, key)
        lines = wrap_text(label, wrap_width)
        for line_idx, line in enumerate(lines[:2]):
            svg.append(
                f'<text x="{x + col_width / 2:.1f}" y="{column_header_y + 14 + line_idx * line_height}" '
                'font-size="10" font-family="Helvetica, Arial, sans-serif" font-weight="700" fill="#111111" '
                'text-anchor="middle">'
                f"{escape(line)}</text>"
            )

    group_boundaries = [x_positions[idx] for idx in (8, 16)]
    current_y = column_header_y + header_height
    previous_group_key: Tuple[str, str] | None = None
    visual_group_index = -1
    for row_idx, (row, wrapped_cells, row_height) in enumerate(zip(rows, row_line_map, row_heights), start=1):
        y = current_y
        group_key = (row.get("gwas_variant_display", ""), row.get("gwas_trait", ""))
        if group_key != previous_group_key:
            visual_group_index += 1
            previous_group_key = group_key
            if row_idx > 1:
                svg.append(
                    f'<line x1="{left_margin}" y1="{y}" '
                    f'x2="{left_margin + table_width}" y2="{y}" '
                    'stroke="#777777" stroke-width="1.1"/>'
                )
        if visual_group_index % 2 == 1:
            svg.append(
                f'<rect x="{left_margin}" y="{y}" width="{table_width}" '
                f'height="{row_height}" fill="#fafafa"/>'
            )
        for col_idx, (key, col_width, _wrap_width, align) in enumerate(column_specs):
            x = x_positions[col_idx]
            lines = wrapped_cells[col_idx]
            if align == "badge":
                badge_fill, badge_text = svg_badge_style(key, row.get(key, "NA"))
                badge_width = min(col_width - 16, max(74, len(row.get(key, "NA")) * 5.8 + 18))
                badge_x = x + 8
                badge_y = y + 8
                svg.append(f'<rect x="{badge_x}" y="{badge_y}" width="{badge_width}" height="18" rx="9" fill="{badge_fill}"/>')
                svg.append(
                    f'<text x="{badge_x + badge_width / 2:.1f}" y="{badge_y + 12}" font-size="9" font-family="Arial" '
                    f'font-weight="700" fill="{badge_text}" text-anchor="middle">{escape(row.get(key, "NA"))}</text>'
                )
                continue
            text_x = x + (8 if align == "text_left" else col_width / 2)
            text_start_y = y + (row_height - line_height * len(lines)) / 2 + line_height * 0.75
            for line_idx, line in enumerate(lines):
                text_y = text_start_y + line_idx * line_height
                if align == "numeric":
                    svg.append(
                        f'<text x="{text_x:.1f}" y="{text_y:.1f}" font-size="9.4" '
                        'font-family="Menlo,Consolas,monospace" fill="#111111" text-anchor="middle" dominant-baseline="middle">'
                        f"{escape(line)}</text>"
                    )
                else:
                    text_anchor = "start" if align == "text_left" else "middle"
                    svg.append(
                        f'<text x="{text_x:.1f}" y="{text_y:.1f}" font-size="9.5" '
                        'font-family="Helvetica, Arial, sans-serif" fill="#111111" '
                        f'text-anchor="{text_anchor}" dominant-baseline="middle">'
                        f"{escape(line)}</text>"
                    )
        current_y += row_height
        svg.append(f'<line x1="{left_margin}" y1="{current_y}" x2="{left_margin + table_width}" y2="{current_y}" stroke="#dddddd" stroke-width="0.7"/>')

    for boundary_x in group_boundaries:
        svg.append(
            f'<line x1="{boundary_x}" y1="{header_y}" x2="{boundary_x}" '
            f'y2="{table_bottom}" stroke="#888888" stroke-width="0.8"/>'
        )
    svg.append(f'<line x1="{left_margin}" y1="{table_bottom}" x2="{left_margin + table_width}" y2="{table_bottom}" stroke="#111111" stroke-width="1.0"/>')
    svg.append("</svg>")
    path.write_text("\n".join(svg), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize five-dog main Allelome x GWAS fine-mapping results and complete GWAS credible sets."
    )
    parser.add_argument(
        "--direct-links",
        type=Path,
        default=DEFAULT_DIRECT_LINKS,
        help="Direct new-gene-pcGene TSV.GZ generated by the current step 1.1.",
    )
    parser.add_argument(
        "--candidate-inventory",
        type=Path,
        default=DEFAULT_CANDIDATE_INVENTORY,
        help="GWAS fine-mapping candidate inventory TSV.GZ generated by the current step 1.1.",
    )
    parser.add_argument(
        "--comparison-manifest",
        type=Path,
        default=DEFAULT_COMPARISON_MANIFEST,
        help="The item/value manifest from the current step 1.1 main output.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Root directory for timestamped results.",
    )
    parser.add_argument(
        "--new-gene-annotation",
        type=Path,
        default=DEFAULT_NEW_GENE_ANNOTATION,
        help="Five-dog new-gene transcript_map TSV file.",
    )
    parser.add_argument(
        "--run-id",
        default=datetime.now().strftime("%Y%m%d_%H%M%S"),
        help="Directory name for this result set; defaults to the current timestamp.",
    )
    parser.add_argument("--workbook-name", default=DEFAULT_WORKBOOK_NAME, help="Excel workbook filename.")
    parser.add_argument("--full-tsv-name", default=DEFAULT_FULL_TSV_NAME, help="Filename for the complete TSV table.")
    parser.add_argument("--focus-tsv-name", default=DEFAULT_FOCUS_TSV_NAME, help="Filename for the focus TSV table.")
    parser.add_argument("--display-csv-name", default=DEFAULT_DISPLAY_CSV_NAME, help="Filename for the display CSV.")
    parser.add_argument("--display-html-name", default=DEFAULT_DISPLAY_HTML_NAME, help="Filename for the display HTML.")
    parser.add_argument("--display-xlsx-name", default=DEFAULT_DISPLAY_XLSX_NAME, help="Filename for the display XLSX.")
    parser.add_argument("--display-svg-name", default=DEFAULT_DISPLAY_SVG_NAME, help="Filename for the display SVG.")
    args = parser.parse_args()

    for path in (
        args.direct_links,
        args.candidate_inventory,
        args.comparison_manifest,
        args.new_gene_annotation,
    ):
        if not path.is_file():
            raise FileNotFoundError(f"input file does not exist: {path}")

    comparison_manifest = read_item_value_table(args.comparison_manifest)
    if comparison_manifest.get("gwas_branch") != BRANCH:
        raise ValueError("comparison manifest is not main branch")
    manifest_target = Path(comparison_manifest.get("direct_new_gene_pcgene_output", ""))
    if manifest_target != args.direct_links:
        raise ValueError(
            "comparison manifest of direct_new_gene_pcgene_output with --direct-links mismatch"
        )
    gwas_run_dir = Path(comparison_manifest.get("gwas_run_dir", ""))
    if not gwas_run_dir.is_absolute() or not gwas_run_dir.is_dir():
        raise FileNotFoundError(
            f"comparison manifest in gwas_run_dir invalid: {gwas_run_dir}"
        )

    run_dir = args.output_root / args.run_id
    if run_dir.exists():
        raise FileExistsError(f"refuse to overwrite existing result directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=False)
    output_paths = {
        "workbook": run_dir / args.workbook_name,
        "full": run_dir / args.full_tsv_name,
        "focus": run_dir / args.focus_tsv_name,
        "display_csv": run_dir / args.display_csv_name,
        "display_html": run_dir / args.display_html_name,
        "display_xlsx": run_dir / args.display_xlsx_name,
        "display_svg": run_dir / args.display_svg_name,
        "gwas_members": run_dir / DEFAULT_GWAS_MEMBERS_NAME,
        "gwas_summary": run_dir / DEFAULT_GWAS_SUMMARY_NAME,
        "summary": run_dir / DEFAULT_SUMMARY_NAME,
        "manifest": run_dir / DEFAULT_MANIFEST_NAME,
    }
    with redirect_stdout(sys.stdout) as log_handle:
        log(f"[start] branch=main; concurrent_tasks=1; threads_per_task=1; direct_links={args.direct_links}", log_handle)
        log("[1/4] read and validate Allelome main input and unified novel gene annotation", log_handle)
        target_rows = read_tsv_rows(args.direct_links)
        inventory_rows = read_tsv_rows(args.candidate_inventory)
        candidate_index = build_candidate_inventory_index(
            inventory_rows, args.candidate_inventory
        )
        annotation_rows = read_tsv_rows(args.new_gene_annotation)
        feelnc_index = build_new_gene_feelnc_index(
            annotation_rows, args.new_gene_annotation
        )
        log(
            f"[input] target_rows={len(target_rows)}; "
            f"candidate_inventory_rows={len(inventory_rows)}; "
            f"new_gene_rows={len(annotation_rows)}",
            log_handle,
        )

        log("[2/4] read complete GWAS main credible set and supplement matching SNP context", log_handle)
        gwas_manifest_rows, gwas_members = read_gwas_credible_members(gwas_run_dir)
        gwas_summary_rows = summarize_gwas_credible_members(
            gwas_manifest_rows, gwas_members
        )
        gwas_all_summary = next(
            row for row in gwas_summary_rows if row.get("scope") == "all_main"
        )
        enriched_rows = build_current_enriched_rows(
            target_rows, feelnc_index, candidate_index
        )
        enriched_rows = sorted(enriched_rows, key=row_sort_key)
        evidence_rows = enrich_credible_set_context(enriched_rows, gwas_members)
        log(
            f"[validated] main_loci={len(gwas_manifest_rows)}; "
            f"credible_memberships={len(gwas_members)}; "
            f"unique_credible_snps={gwas_all_summary['unique_credible_set_snps']}; "
            f"direct_new_gene_pcgene_rows={len(evidence_rows)}",
            log_handle,
        )

        focus_rows = list(evidence_rows)
        display_rows = build_display_rows(focus_rows)
        display_summary = build_display_summary(
            evidence_rows, focus_rows, gwas_all_summary
        )

        log("[3/4] write credible set member, summary table and Allelome evidence table", log_handle)
        write_tsv(output_paths["full"], list(dict.fromkeys(k for r in evidence_rows for k in r)), evidence_rows)
        write_tsv(output_paths["focus"], FOCUS_COLUMNS, focus_rows)
        write_tsv(
            output_paths["gwas_members"],
            GWAS_CREDIBLE_SET_MEMBER_COLUMNS,
            gwas_members,
        )
        write_tsv(
            output_paths["gwas_summary"],
            GWAS_CREDIBLE_SUMMARY_COLUMNS,
            gwas_summary_rows,
        )
        summary_metrics = [
            ("run_id", args.run_id),
            ("gwas_branch", BRANCH),
            ("direct_new_gene_pcgene_input_rows", len(target_rows)),
            ("credible_set_expanded_rows", len(evidence_rows)),
            ("target_unique_gwas_snps", len({row.get("gwas_snp_id", "") for row in evidence_rows})),
            ("gwas_run_dir", str(gwas_run_dir)),
            ("gwas_locus_count", gwas_all_summary["gwas_locus_count"]),
            ("gwas_credible_set_count", gwas_all_summary["credible_set_count"]),
            ("gwas_credible_set_membership_rows", gwas_all_summary["credible_set_membership_rows"]),
            ("gwas_unique_credible_set_snps", gwas_all_summary["unique_credible_set_snps"]),
            ("gwas_p_lrt_lt_1e-5_unique_snps", gwas_all_summary["p_lrt_lt_1e-5_unique_snps"]),
            ("gwas_min_p_lrt", gwas_all_summary["min_gwas_p_lrt"]),
            ("gwas_min_p_lrt_snp", gwas_all_summary["min_gwas_p_lrt_snp"]),
            ("gwas_min_p_lrt_snp_coordinate", gwas_all_summary["min_gwas_p_lrt_snp_coordinate"]),
            ("interpretation_boundary", "Positional association between Allelome and GWAS fine-map candidates; not proof of colocalization or causality."),
        ]
        write_tsv(
            output_paths["summary"],
            SUMMARY_COLUMNS,
            [{"metric": key, "value": value} for key, value in summary_metrics],
        )

        log("[4/4] write Excel, CSV, HTML, XLSX, SVG and manifest", log_handle)
        write_excel_workbook(
            output_paths["workbook"], all_rows=evidence_rows, focus_rows=focus_rows
        )
        write_csv_table(output_paths["display_csv"], DISPLAY_COLUMNS, display_rows)
        title = "Five-dog main Allelome x GWAS fine-mapping summary"
        write_html_table(
            output_paths["display_html"], title=title,
            fieldnames=DISPLAY_COLUMNS, rows=display_rows, summary=display_summary,
        )
        write_xlsx_table(
            output_paths["display_xlsx"], title=title,
            fieldnames=DISPLAY_COLUMNS, rows=display_rows,
        )
        write_svg_table(
            output_paths["display_svg"], title=title,
            rows=display_rows, summary=display_summary,
        )

        manifest_rows = []
        for role, path, row_count in (
            ("input_direct_links", args.direct_links, len(target_rows)),
            ("input_candidate_inventory", args.candidate_inventory, len(inventory_rows)),
            ("input_comparison_manifest", args.comparison_manifest, len(comparison_manifest)),
            ("input_new_gene_annotation", args.new_gene_annotation, len(annotation_rows)),
            ("input_gwas_locus_manifest", gwas_run_dir / "locus_run_manifest.tsv", len(gwas_manifest_rows)),
        ):
            manifest_rows.append(
                {"role": role, "path": str(path), "exists": "yes",
                 "size_bytes": path.stat().st_size, "rows": row_count,
                 "schema_status": "PASS"}
            )
        for role, path in output_paths.items():
            if role == "manifest":
                continue
            manifest_rows.append(
                {"role": f"output_{role}", "path": str(path), "exists": "yes",
                 "size_bytes": path.stat().st_size, "rows": "NA",
                 "schema_status": "written"}
            )
        write_tsv(output_paths["manifest"], MANIFEST_COLUMNS, manifest_rows)
        (run_dir / "execution_completion.marker").write_text("execution_status=success\n", encoding="utf-8")
        log(f"[done] outputs={len(output_paths)}; output_dir={run_dir}; task completed / done; execution_status=success", log_handle)
        log("reminder: execution completed does not imply result QC passed, please check separately format, count and biological interpretation.", log_handle)
    return 0


if __name__ == "__main__":

    log_dir = Path("/path/to/repository/scripts/12_evidence_integration/5.allelome-compare-gwas-finemap/summary/scripts/log")
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    with (log_dir / f"main_summary.{stamp}.log").open("x", encoding="utf-8") as handle:
        with redirect_stdout(handle), redirect_stderr(handle):
            success_path = log_dir / f"main_summary.{stamp}.success.tsv"
            failed_path = log_dir / f"main_summary.{stamp}.failed.tsv"
            write_tsv(success_path, ["task", "exit_code"], [])
            write_tsv(failed_path, ["task", "exit_code"], [])
            try:
                code = main()
            except BaseException:
                traceback.print_exc()
                write_tsv(failed_path, ["task", "exit_code"], [{"task": "main_summary", "exit_code": "1"}])
                raise
            write_tsv(success_path, ["task", "exit_code"], [{"task": "main_summary", "exit_code": str(code)}])
    raise SystemExit(code)
