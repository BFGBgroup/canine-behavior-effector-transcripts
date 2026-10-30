#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import gzip
import html
import re
import sys
import signal
import traceback
from contextlib import ExitStack, redirect_stdout, redirect_stderr
import zipfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, Mapping, Sequence
from xml.sax.saxutils import escape as xml_escape


DEFAULT_INPUT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/4.allelome-compare-raw-gwas/1.1_compare_dog5_allelome_raw_gwas/20260910_171716/dog5_allelome_raw_gwas_overlap_all_p_le_1e-4.tsv.gz'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/4.allelome-compare-raw-gwas/ncrna-mrna'
)
GWAS_P_LRT_THRESHOLD = 5e-5
PAIR_CLASS = "annotated_ncRNA--pcGene"
SCRIPT_DIR = Path("/path/to/repository/scripts/12_evidence_integration/4.allelome-compare-raw-gwas/ncrna-mrna/scripts")
SVG_ROW_LIMIT = 100
BAR_LIMIT = 20


SVG_FIELD_LABELS = {
    "trait": "Trait",
    "gwas_marker_id": "GWAS SNP",
    "gwas_p_lrt": "GWAS p_lrt",
    "ncrna_gene_symbol": "ncRNA",
    "ncrna_interval": "ncRNA interval",
    "pcgene_gene_symbol": "pcGene",
    "mechanism": "Allelome mode",
    "linkage_score": "Linkage",
    "sample_tissue": "Sample / tissue",
    "interval_gap_bp": "Gene gap (bp)",
    "evidence_stratum": "Evidence stratum",
    "gwas_hit_relation": "GWAS hit side",
}
SVG_FIELD_WIDTHS = {
    "trait": 145,
    "gwas_marker_id": 145,
    "gwas_p_lrt": 105,
    "ncrna_gene_symbol": 125,
    "ncrna_interval": 180,
    "pcgene_gene_symbol": 130,
    "mechanism": 95,
    "linkage_score": 75,
    "sample_tissue": 155,
    "interval_gap_bp": 95,
    "evidence_stratum": 230,
    "gwas_hit_relation": 100,
}
SVG_FIELD_GROUPS = (
    ("GWAS evidence", ("trait", "gwas_marker_id", "gwas_p_lrt")),
    ("ncRNA–pcGene relation", ("ncrna_gene_symbol", "ncrna_interval", "pcgene_gene_symbol")),
    ("Allelome evidence", ("mechanism", "linkage_score")),
    ("Context", ("sample_tissue", "interval_gap_bp")),
    ("Review", ("gwas_hit_relation", "evidence_stratum")),
)

REQUIRED_FIELDS = {
    "trait_variant_key", "sample_pair_id", "evidence_stratum", "source_link_file", "source_link_line",
    "trait",
    "gwas_marker_id",
    "gwas_chrom",
    "gwas_pos",
    "gwas_p_lrt",
    "gwas_p_tier",
    "gwas_allele1",
    "gwas_allele0",
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
    "base_gene_id",
    "base_gene_symbol",
    "base_gene_strand",
    "base_chr",
    "base_start",
    "base_end",
    "base_annotation_source",
    "base_gene_biotype",
    "base_record_type",
    "target_gene_id",
    "target_gene_symbol",
    "target_gene_strand",
    "target_chr",
    "target_start",
    "target_end",
    "target_annotation_source",
    "target_gene_biotype",
    "target_record_type",
}

DERIVED_FIELDS = [
    "evidence_scope",
    "relation_direction",
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
    "gwas_marker_id",
    "gwas_chrom",
    "gwas_pos",
    "gwas_allele1",
    "gwas_allele0",
    "gwas_unordered_alleles",
    "gwas_p_lrt",
    "gwas_p_wald",
    "gwas_p_score",
    "gwas_af",
    "gwas_beta",
    "gwas_se",
    "gwas_p_tier",
    "gwas_hit_on_ncrna",
    "gwas_hit_on_pcgene",
    "gwas_hit_relation",
    "ncrna_gene_id",
    "ncrna_gene_symbol",
    "ncrna_strand",
    "ncrna_chr",
    "ncrna_start",
    "ncrna_end",
    "ncrna_interval",
    "pcgene_gene_id",
    "pcgene_gene_symbol",
    "pcgene_strand",
    "pcgene_chr",
    "pcgene_start",
    "pcgene_end",
    "pcgene_interval",
    "mechanism",
    "linkage_score",
    "interval_gap_bp",
    "animal_id",
    "tissue",
    "sample_tissue",
    "canonical_pair_id",
    "sample_pair_id", "trait_variant_key", "source_link_file", "source_link_line",
    "evidence_stratum", "either_endpoint_multigene", "overlap_bp", "strand_relation",
    "ncrna_pcgene_pair_key",
    "ncrna_pcgene_pair_label",
]

SUMMARY_FIELDS = [
    "metric",
    "value",
    "definition",
]

PAIR_SUMMARY_FIELDS = [
    "ncrna_gene_id",
    "ncrna_gene_symbol",
    "pcgene_gene_id",
    "pcgene_gene_symbol",
    "pair_key",
    "row_count",
    "unique_gwas_variant_count",
    "traits",
    "animals",
    "tissues",
    "gwas_hit_on_ncrna_rows",
    "gwas_hit_on_pcgene_rows",
    "enhancing_rows",
    "repressing_rows",
    "minimum_gwas_p_lrt",
    "minimum_interval_gap_bp",
    "maximum_linkage_score",
]


def stamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class RunLogger:


    def __init__(self, path: Path, run_dir: Path, run_id: str) -> None:
        self.path, self.run_dir, self.run_id = path, run_dir, run_id
        self.stack = ExitStack()

    def log(self, message: str) -> None:
        print(f"[{stamp()} ] {message}", flush=True)

    def __enter__(self) -> "RunLogger":
        try:
            self.handle = self.stack.enter_context(self.path.open("x", encoding="utf-8"))
            prefix = str(self.path)[:-4]
            self.success = self.stack.enter_context(Path(prefix + ".success.tsv").open("x", encoding="utf-8"))
            self.failed = self.stack.enter_context(Path(prefix + ".failed.tsv").open("x", encoding="utf-8"))
            self.success.write("run_id\toutput_dir\n")
            self.failed.write("run_id\terror\n")
            self.success.flush()
            self.failed.flush()
            self.stack.enter_context(redirect_stdout(self.handle))
            self.stack.enter_context(redirect_stderr(self.handle))
            return self
        except BaseException:
            self.stack.close()
            raise

    def __exit__(self, exc_type, exc, tb) -> None:
        try:
            if exc_type is None:
                (self.run_dir / "execution_completion.marker").write_text("execution_status=success\n", encoding="utf-8")
                self.success.write(f"{self.run_id}\t{self.run_dir}\n")
                self.success.flush()
                self.log("task completed / done; execution_status=success")
                self.log("execution completed does not imply result QC passed; please separately check format, count and biological interpretation.")
            else:
                self.failed.write(f"{self.run_id}\t{str(exc).replace(chr(9), ' ').replace(chr(10), ' ')}\n")
                self.failed.flush()
                traceback.print_exception(exc_type, exc, tb)
        finally:
            self.stack.close()


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


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
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        text = str(value).strip()
        if text and text not in seen:
            seen.add(text)
            result.append(text)
    return result


def is_newgene_auxiliary_field(field: str) -> bool:

    lowered = field.strip().lower()
    return "new_gene" in lowered or "manual_new_gene" in lowered


def output_input_fields(fields: Sequence[str]) -> list[str]:
    return [field for field in fields if not is_newgene_auxiliary_field(field)]


def norm_chr(value: object) -> str:
    return str(value or "").strip().lower().removeprefix("chr")


def side_is_ncrna(row: Mapping[str, str], side: str) -> bool:

    return row[f"{side}_record_type"] == "known_gene" and row[f"{side}_gene_biotype"] in {"lncRNA", "misc_RNA"}


def side_is_pcgene(row: Mapping[str, str], side: str) -> bool:
    return row[f"{side}_record_type"] == "known_gene" and row[f"{side}_gene_biotype"] == "protein_coding"


def standard_pair(row: Mapping[str, str]) -> bool:

    return (side_is_ncrna(row, "base") and side_is_pcgene(row, "target")) or (
        side_is_ncrna(row, "target") and side_is_pcgene(row, "base")
    )


def side_value(row: Mapping[str, str], side: str, name: str) -> str:
    return row.get(f"{side}_{name}", "").strip()


def interval_for(row: Mapping[str, str], side: str) -> str:
    chromosome = side_value(row, side, "chr")
    start = side_value(row, side, "start")
    end = side_value(row, side, "end")
    return f"{chromosome}:{int(start) + 1}-{end}" if chromosome and start and end else ""


def coordinate_contains(row: Mapping[str, str], side: str) -> bool:
    gwas_chr = norm_chr(row.get("gwas_chrom", ""))
    side_chr = norm_chr(side_value(row, side, "chr"))
    position = as_int(row.get("gwas_pos", ""))
    start = as_int(side_value(row, side, "start"))
    end = as_int(side_value(row, side, "end"))


    return bool(gwas_chr and side_chr and gwas_chr == side_chr and position is not None and start is not None and end is not None and start < position <= end)


def gwas_hits_side(row: Mapping[str, str], side: str) -> bool:

    return coordinate_contains(row, side)


def oriented_row(row: Mapping[str, str], fields: Sequence[str]) -> dict[str, str]:
    if side_is_ncrna(row, "base") and side_is_pcgene(row, "target"):
        ncrna_side, pc_side = "base", "target"
    elif side_is_ncrna(row, "target") and side_is_pcgene(row, "base"):
        ncrna_side, pc_side = "target", "base"
    else:
        raise ValueError("Cannot assign annotated_ncRNA and pcGene roles to this pair: " + row.get("canonical_pair_id", ""))

    ncrna_hit = gwas_hits_side(row, ncrna_side)
    pc_hit = gwas_hits_side(row, pc_side)
    if ncrna_hit and pc_hit:
        hit_relation = "both"
    elif ncrna_hit:
        hit_relation = "ncrna_only"
    elif pc_hit:
        hit_relation = "pcgene_only"
    else:
        hit_relation = "neither"

    result = {field: str(row.get(field, "")) for field in fields}
    result.update(
        {
            "evidence_scope": PAIR_CLASS,
            "relation_direction": "undirected_Allelome_pair",
            "ncrna_side": ncrna_side,
            "ncrna_gene_id": side_value(row, ncrna_side, "gene_id"),
            "ncrna_gene_symbol": side_value(row, ncrna_side, "gene_symbol"),
            "ncrna_strand": side_value(row, ncrna_side, "gene_strand"),
            "ncrna_chr": side_value(row, ncrna_side, "chr"),
            "ncrna_start": side_value(row, ncrna_side, "start"),
            "ncrna_end": side_value(row, ncrna_side, "end"),
            "ncrna_interval": interval_for(row, ncrna_side),
            "ncrna_biotype": side_value(row, ncrna_side, "gene_biotype"),
            "pcgene_side": pc_side,
            "pcgene_gene_id": side_value(row, pc_side, "gene_id"),
            "pcgene_gene_symbol": side_value(row, pc_side, "gene_symbol"),
            "pcgene_strand": side_value(row, pc_side, "gene_strand"),
            "pcgene_chr": side_value(row, pc_side, "chr"),
            "pcgene_start": side_value(row, pc_side, "start"),
            "pcgene_end": side_value(row, pc_side, "end"),
            "pcgene_interval": interval_for(row, pc_side),
            "pcgene_biotype": side_value(row, pc_side, "gene_biotype"),
            "gwas_hit_on_ncrna": "yes" if ncrna_hit else "no",
            "gwas_hit_on_pcgene": "yes" if pc_hit else "no",
            "gwas_hit_relation": hit_relation,
            "ncrna_mapping_rule": row.get(f"{ncrna_side}_mapping_rule", ""),
            "pcgene_mapping_rule": row.get(f"{pc_side}_mapping_rule", ""),
            "ncrna_pcgene_pair_key": "|".join(
                [
                    side_value(row, ncrna_side, "gene_id") or side_value(row, ncrna_side, "gene_symbol"),
                    side_value(row, pc_side, "gene_id") or side_value(row, pc_side, "gene_symbol"),
                ]
            ),
            "ncrna_pcgene_pair_label": " ↔ ".join(
                [
                    side_value(row, ncrna_side, "gene_symbol") or side_value(row, ncrna_side, "gene_id"),
                    side_value(row, pc_side, "gene_symbol") or side_value(row, pc_side, "gene_id"),
                ]
            ),
        }
    )
    return result


def build_pair_summary(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row.get("ncrna_pcgene_pair_key", "")].append(row)

    result: list[dict[str, str]] = []
    for key, members in groups.items():
        p_values = [p for p in (as_float(row.get("gwas_p_lrt", "")) for row in members) if p is not None]
        distances = [d for d in (as_int(row.get("interval_gap_bp", "")) for row in members) if d is not None]
        linkages = [x for x in (as_float(row.get("linkage_score", "")) for row in members) if x is not None]
        first = members[0]
        result.append(
            {
                "ncrna_gene_id": first.get("ncrna_gene_id", ""),
                "ncrna_gene_symbol": first.get("ncrna_gene_symbol", ""),
                "pcgene_gene_id": first.get("pcgene_gene_id", ""),
                "pcgene_gene_symbol": first.get("pcgene_gene_symbol", ""),
                "pair_key": key,
                "row_count": str(len(members)),
                "unique_gwas_variant_count": str(len(unique(x["trait_variant_key"] for x in members))),
                "traits": ";".join(unique(x.get("trait", "") for x in members)),
                "animals": ";".join(unique(x.get("animal_id", "") for x in members)),
                "tissues": ";".join(unique(x.get("tissue", "") for x in members)),
                "gwas_hit_on_ncrna_rows": str(sum(x.get("gwas_hit_on_ncrna") == "yes" for x in members)),
                "gwas_hit_on_pcgene_rows": str(sum(x.get("gwas_hit_on_pcgene") == "yes" for x in members)),
                "enhancing_rows": str(sum(x.get("mechanism", "").lower() == "enhancing" for x in members)),
                "repressing_rows": str(sum(x.get("mechanism", "").lower() == "repressing" for x in members)),
                "minimum_gwas_p_lrt": f"{min(p_values):.12g}" if p_values else "",
                "minimum_interval_gap_bp": str(min(distances)) if distances else "",
                "maximum_linkage_score": f"{max(linkages):.6g}" if linkages else "",
            }
        )
    result.sort(key=lambda x: (-int(x["row_count"]), as_float(x["minimum_gwas_p_lrt"]) or float("inf"), x["pair_key"]))
    return result


def build_summary(rows_all: Sequence[Mapping[str, str]], rows_focus: Sequence[Mapping[str, str]], input_rows: int, run_dir: Path, input_path: Path) -> list[dict[str, str]]:
    pairs_all = unique(row.get("ncrna_pcgene_pair_key", "") for row in rows_all)
    pairs_focus = unique(row.get("ncrna_pcgene_pair_key", "") for row in rows_focus)
    metrics = [
        ("input_rows", input_rows, "All expanded GWAS×Allelome rows in the original input."),
        ("standard_annotated_ncrna_pcgene_rows", len(rows_all), "Expanded annotated_ncRNA--pcGene rows after filtering for p_lrt < 5e-5."),
        ("focus_p_lt_1e5_rows", sum(float(r["gwas_p_lrt"]) < 1e-5 for r in rows_focus), "Expanded focus rows with p_lrt<1e-5."),
        ("focus_both_rows", sum(r["gwas_hit_relation"] == "both" for r in rows_focus), "The same marker overlaps both endpoint gene spans and is retained as one record."),
        ("focus_multigene_review_rows", sum(r["evidence_stratum"] == "nonoverlap_multigene_review" for r in rows_focus), "Non-overlapping multi-gene assignment records requiring review in focus."),
        ("focus_overlap_review_rows", sum(r["evidence_stratum"] == "overlap_review" for r in rows_focus), "Endpoint-interval overlap records requiring review in focus."),
        ("focus_gwas_hit_on_ncrna_rows", len(rows_focus), "Rows above in which the GWAS site falls within an annotated ncRNA gene span."),
        ("unique_annotated_ncrna_pcgene_pairs", len(pairs_all), "Number of relationships deduplicated by ncRNA gene ID and pcGene ID."),
        ("unique_focus_pairs", len(pairs_focus), "Number of deduplicated relationships with at least one GWAS site overlapping the ncRNA."),
        ("unique_ncrna_count", len(unique(row.get("ncrna_gene_id", "") for row in rows_all)), "Number of annotated noncoding RNAs involved."),
        ("unique_pcgene_count", len(unique(row.get("pcgene_gene_id", "") for row in rows_all)), "Number of protein-coding genes involved."),
        ("ncrna_hit_rows", sum(row.get("gwas_hit_on_ncrna") == "yes" for row in rows_all), "Expanded rows with a GWAS site overlapping the ncRNA gene span."),
        ("pcgene_hit_rows", sum(row.get("gwas_hit_on_pcgene") == "yes" for row in rows_all), "Expanded rows with a GWAS site overlapping the pcGene gene span."),
        ("enhancing_rows", sum(row.get("mechanism", "").lower() == "enhancing" for row in rows_all), "Expanded rows with an enhancing Allelome mechanism."),
        ("repressing_rows", sum(row.get("mechanism", "").lower() == "repressing" for row in rows_all), "Expanded rows with a repressing Allelome mechanism."),
        ("gwas_p_lrt_threshold", f"<{GWAS_P_LRT_THRESHOLD:g}", "Raw GWAS filtering threshold retained from the historical script."),
        ("input_file", str(input_path), "Absolute path actually read."),
        ("run_output_directory", str(run_dir), "Timestamped output directory for this run."),
    ]
    return [{"metric": str(metric), "value": str(value), "definition": definition} for metric, value, definition in metrics]


def write_html(path: Path, title: str, fields: Sequence[str], rows: Sequence[Mapping[str, str]], note: str) -> None:
    headers = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body_parts: list[str] = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields)
        body_parts.append(f"<tr>{cells}</tr>")
    body = "\n".join(body_parts) if body_parts else '<tr><td colspan="2">No eligible records</td></tr>'
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font-family:Arial,'Microsoft YaHei',sans-serif;margin:24px;color:#1f2933}}
h1{{font-size:20px}}p.note{{color:#52606d;font-size:13px}}
table{{border-collapse:collapse;min-width:1250px;font-size:12px}}th,td{{border:1px solid #cbd2d9;padding:5px 7px;white-space:nowrap}}
th{{background:#e8f0ec;position:sticky;top:0}}tr:nth-child(even){{background:#f8faf9}}</style></head>
<body><h1>{html.escape(title)}</h1><p class="note">{html.escape(note)}</p>
<table><thead><tr>{headers}</tr></thead><tbody>{body}</tbody></table></body></html>"""
    path.write_text(document, encoding="utf-8")


def svg_text(text: object) -> str:
    return xml_escape(str(text), {"\"": "&quot;", "'": "&apos;"})


def svg_label_lines(label: str, max_chars: int = 16) -> list[str]:

    if len(label) <= max_chars:
        return [label]
    words = label.split(" ")
    if len(words) > 1:
        first: list[str] = []
        second: list[str] = []
        for word in words:
            target = first if not second and sum(len(x) for x in first) + len(first) + len(word) <= max_chars else second
            target.append(word)
        if first and second:
            return [" ".join(first), " ".join(second)]
    return [label[: max_chars - 1] + "…", label[max_chars - 1 : max_chars * 2 - 1]]


def svg_display_value(field: str, value: object) -> str:

    text = str(value or "")
    if field == "gwas_p_lrt":
        number = as_float(text)
        return f"{number:.3e}" if number is not None else text
    if field == "interval_gap_bp":
        number = as_int(text)
        return f"{number:,}" if number is not None else text
    return text


def svg_truncate(text: str, max_chars: int) -> str:
    if len(text) <= max_chars:
        return text
    return text[: max(1, max_chars - 1)] + "…"


def write_svg_table(path: Path, title: str, fields: Sequence[str], rows: Sequence[Mapping[str, str]], note: str) -> None:
    display_rows = list(rows[:SVG_ROW_LIMIT])
    widths = [SVG_FIELD_WIDTHS.get(field, 120) for field in fields]
    row_height = 25
    group_height = 24
    header_height = 46
    title_height = 72
    footer_height = 32
    left_margin = 15
    width = max(1000, sum(widths) + left_margin * 2)
    height = title_height + group_height + header_height + row_height * max(1, len(display_rows)) + footer_height
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    parts.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    parts.append(f'<text x="{left_margin}" y="27" font-family="Arial,Microsoft YaHei" font-size="17" font-weight="bold">{svg_text(title)}</text>')
    parts.append(f'<text x="{left_margin}" y="51" font-family="Arial,Microsoft YaHei" font-size="11" fill="#52606d">{svg_text(note)}</text>')

    table_top = title_height
    group_y = table_top
    header_y = table_top + group_height
    data_y = header_y + header_height


    group_x = left_margin
    field_positions: dict[str, tuple[int, int]] = {}
    for field, width_i in zip(fields, widths):
        field_positions[field] = (group_x, width_i)
        group_x += width_i
    used_fields: set[str] = set()
    for group_label, group_fields in SVG_FIELD_GROUPS:
        present = [field for field in group_fields if field in field_positions]
        if not present:
            continue
        start_x = field_positions[present[0]][0]
        end_x = field_positions[present[-1]][0] + field_positions[present[-1]][1]
        group_width = end_x - start_x
        parts.append(f'<rect x="{start_x}" y="{group_y}" width="{group_width}" height="{group_height}" fill="#dfeae4" stroke="#9fb3aa"/>')
        parts.append(f'<text x="{start_x + group_width / 2:.1f}" y="{group_y + 16}" text-anchor="middle" font-family="Arial,Microsoft YaHei" font-size="10" font-weight="bold">{svg_text(group_label)}</text>')
        used_fields.update(present)

    for field in fields:
        if field in used_fields:
            continue
        x_i, width_i = field_positions[field]
        parts.append(f'<rect x="{x_i}" y="{group_y}" width="{width_i}" height="{group_height}" fill="#dfeae4" stroke="#9fb3aa"/>')

    for field, width_i in zip(fields, widths):
        x, _ = field_positions[field]
        parts.append(f'<rect x="{x}" y="{header_y}" width="{width_i}" height="{header_height}" fill="#eef5f1" stroke="#b7c7bf"/>')
        label = SVG_FIELD_LABELS.get(field, field)
        lines = svg_label_lines(label)
        first_y = header_y + (20 if len(lines) == 1 else 15)
        for line_index, line in enumerate(lines[:2]):
            parts.append(f'<text x="{x + width_i / 2:.1f}" y="{first_y + line_index * 15}" text-anchor="middle" font-family="Arial,Microsoft YaHei" font-size="10" font-weight="bold">{svg_text(line)}</text>')

    for index, row in enumerate(display_rows):
        y = data_y + index * row_height
        fill = "#ffffff" if index % 2 == 0 else "#f7faf8"
        for field, width_i in zip(fields, widths):
            x, _ = field_positions[field]
            value = svg_display_value(field, row.get(field, ""))
            max_chars = max(8, int(width_i / 7.0))
            shown = svg_truncate(value, max_chars)
            parts.append(f'<g><title>{svg_text(value)}</title><rect x="{x}" y="{y}" width="{width_i}" height="{row_height}" fill="{fill}" stroke="#d9e2dc"/>')
            if field == "mechanism":
                mechanism_fill = "#fff1e6" if value == "enhancing" else "#fbe9ea" if value == "repressing" else fill
                parts.append(f'<rect x="{x + 5}" y="{y + 5}" width="{max(10, width_i - 10)}" height="{row_height - 10}" rx="8" fill="{mechanism_fill}"/>')
            text_fill = "#8a3b12" if field == "mechanism" and value == "enhancing" else "#8b1e2d" if field == "mechanism" and value == "repressing" else "#1f2933"
            parts.append(f'<text x="{x + width_i / 2:.1f}" y="{y + 17}" text-anchor="middle" font-family="Arial,Microsoft YaHei" font-size="9.5" fill="{text_fill}">{svg_text(shown)}</text></g>')

    footer_y = data_y + row_height * max(1, len(display_rows)) + 21
    parts.append(f'<text x="{left_margin}" y="{footer_y}" font-family="Arial,Microsoft YaHei" font-size="10" fill="#52606d">Showing {len(display_rows)} / {len(rows)} rows; full fields are available in TSV, CSV, and XLSX.</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def write_svg_pair_barplot(path: Path, pair_rows: Sequence[Mapping[str, str]], note: str) -> None:
    selected = list(pair_rows[:BAR_LIMIT])
    left = 310
    plot_width = 830
    row_height = 30
    title_height = 86
    height = title_height + max(1, len(selected)) * row_height + 40
    width = left + plot_width + 100
    max_count = max((int(row["row_count"]) for row in selected), default=1)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    parts.append('<rect width="100%" height="100%" fill="#ffffff"/>')
    parts.append('<text x="20" y="29" font-family="Arial,Microsoft YaHei" font-size="18" font-weight="bold">Expanded GWAS rows for annotated ncRNA–pcGene relationships</text>')
    parts.append(f'<text x="20" y="52" font-family="Arial,Microsoft YaHei" font-size="12" fill="#52606d">{svg_text(note)}</text>')
    parts.append('<rect x="20" y="66" width="14" height="14" fill="#9aa5b1"/><text x="40" y="78" font-family="Arial" font-size="11">All expanded GWAS×Allelome rows</text>')
    parts.append('<rect x="220" y="66" width="14" height="14" fill="#d97706"/><text x="240" y="78" font-family="Arial" font-size="11">GWAS within ncRNA gene span</text>')
    for index, row in enumerate(selected):
        y = title_height + index * row_height
        label = f"{row.get('ncrna_gene_symbol') or row.get('ncrna_gene_id')} ↔ {row.get('pcgene_gene_symbol') or row.get('pcgene_gene_id')}"
        if len(label) > 43:
            label = label[:40] + "..."
        total = int(row["row_count"])
        hit = int(row["gwas_hit_on_ncrna_rows"])
        total_w = max(1, int(plot_width * total / max_count))
        hit_w = int(plot_width * hit / max_count)
        parts.append(f'<text x="20" y="{y + 19}" font-family="Arial,Microsoft YaHei" font-size="11">{svg_text(label)}</text>')
        parts.append(f'<rect x="{left}" y="{y + 4}" width="{total_w}" height="18" fill="#9aa5b1"/>')
        if hit_w:
            parts.append(f'<rect x="{left}" y="{y + 4}" width="{hit_w}" height="18" fill="#d97706"/>')
        parts.append(f'<text x="{left + total_w + 6}" y="{y + 18}" font-family="Arial" font-size="10">{total} / {hit}</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def xlsx_cell(value: object) -> str:
    text = xml_escape(str(value if value is not None else ""))
    return f'<c t="inlineStr"><is><t>{text}</t></is></c>'


def write_xlsx(path: Path, sheets: Sequence[tuple[str, Sequence[str], Sequence[Mapping[str, object]]]]) -> None:
    workbook = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>{}</sheets></workbook>'''
    rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{}</Relationships>'''
    sheet_nodes: list[str] = []
    rel_nodes: list[str] = []
    content_types = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>', '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">', '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>', '<Default Extension="xml" ContentType="application/xml"/>', '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>']
    for index, (name, fields, rows) in enumerate(sheets, 1):
        sheet_nodes.append(f'<sheet name="{xml_escape(name)}" sheetId="{index}" r:id="rId{index}"/>')
        rel_nodes.append(f'<Relationship Id="rId{index}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{index}.xml"/>')
        content_types.append(f'<Override PartName="/xl/worksheets/sheet{index}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
    content_types.append('</Types>')
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", "".join(content_types))
        archive.writestr("_rels/.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        archive.writestr("xl/workbook.xml", workbook.format("".join(sheet_nodes)))
        archive.writestr("xl/_rels/workbook.xml.rels", rels.format("".join(rel_nodes)))
        for index, (_name, fields, rows) in enumerate(sheets, 1):
            row_xml: list[str] = []
            row_xml.append("<row>" + "".join(xlsx_cell(field) for field in fields) + "</row>")
            for row in rows:
                row_xml.append("<row>" + "".join(xlsx_cell(row.get(field, "")) for field in fields) + "</row>")
            sheet_xml = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + "".join(row_xml) + "</sheetData></worksheet>"
            archive.writestr(f"xl/worksheets/sheet{index}.xml", sheet_xml)


def xls_cell(value: object) -> str:
    return f'<Cell><Data ss:Type="String">{xml_escape(str(value if value is not None else ""))}</Data></Cell>'


def write_xls(path: Path, sheets: Sequence[tuple[str, Sequence[str], Sequence[Mapping[str, object]]]]) -> None:
    parts = ['<?xml version="1.0"?>', '<?mso-application progid="Excel.Sheet"?>', '<Workbook xmlns="urn:schemas-microsoft-com:office:spreadsheet" xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet">']
    for name, fields, rows in sheets:
        parts.append(f'<Worksheet ss:Name="{xml_escape(name)}"><Table>')
        parts.append('<Row>' + "".join(xls_cell(field) for field in fields) + '</Row>')
        for row in rows:
            parts.append('<Row>' + "".join(xls_cell(row.get(field, "")) for field in fields) + '</Row>')
        parts.append('</Table></Worksheet>')
    parts.append('</Workbook>')
    path.write_text("".join(parts), encoding="utf-8")


def write_manifest(path: Path, run_dir: Path, input_path: Path, fields: Sequence[str], output_fields: Sequence[str], counts: Mapping[str, object]) -> None:
    rows = [
        {"key": "script", "value": "build_annotated_ncrna_pcgene_raw_gwas_allelome.py"},
        {"key": "input_path", "value": str(input_path)},
        {"key": "output_run_dir", "value": str(run_dir)},
        {"key": "gwas_p_lrt_rule", "value": f"gwas_p_lrt < {GWAS_P_LRT_THRESHOLD:g}"},
        {"key": "pair_class_rule", "value": PAIR_CLASS},
        {"key": "coordinate_rule", "value": "Same chromosome and START0 < 1-based GWAS position <= END0 (following the upstream manifest)."},
        {"key": "relation_direction", "value": "Allelome pair is undirected; base/target orientation is normalized"},
        {"key": "raw_input_columns", "value": str(len(fields))},
        {"key": "retained_output_input_columns", "value": str(len(output_fields))},
        {"key": "excluded_newgene_auxiliary_columns", "value": ",".join(field for field in fields if is_newgene_auxiliary_field(field))},
    ]
    rows.extend({"key": key, "value": str(value)} for key, value in counts.items())
    write_tsv(path, ["key", "value"], rows)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description="Extract raw-GWAS-Allelome links for annotated_ncRNA--pcGene pairs and generate tables and figures.")
    result.add_argument("--overlap-detail", type=Path, default=DEFAULT_INPUT, help="Absolute path to the raw GWAS-Allelome overlap TSV/TSV.GZ file.")
    result.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT, help="Absolute path to the timestamped output root directory.")
    result.add_argument("--run-id", default=None, help="Optional timestamped directory name, for example 20260816_120000.")
    return result


def main() -> int:
    args = parser().parse_args()
    def interrupted(signum, frame):
        raise InterruptedError(f"received termination signal {signum}")
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    input_path = args.overlap_detail.expanduser()
    output_root = args.output_root.expanduser()
    if not input_path.is_absolute() or not output_root.is_absolute():
        print("error: --overlap-detail and --output-root must use absolute path.", file=sys.stderr)
        return 2
    if not input_path.is_file():
        print(f"error: input file does not exist: {input_path}", file=sys.stderr)
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

    log_dir = SCRIPT_DIR / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{Path(__file__).stem}.{run_id}.log"
    with RunLogger(log_path, run_dir, run_id) as log:
        log.log("task started / start; retain serial 1 task, per-task 1 thread(only standard library, no numerical-library thread pool)")
        log.log(f"[stage 1/5] read input: {input_path}")
        fields, input_rows = read_tsv(input_path)
        missing = sorted(REQUIRED_FIELDS.difference(fields))
        if missing:
            log.log("input field missing: " + ", ".join(missing))
            raise ValueError("input field missing: " + ", ".join(missing))
        log.log(f"input table header check passed: {len(fields)} column; read {len(input_rows)} row")

        retained_input_fields = output_input_fields(fields)
        log.log(
            "[stage 2/5] filter annotated_ncRNA--pcGene and standardize the orientation of both endpoints; "
            f"retain {len(retained_input_fields)} items GWAS/Allelome core input column, "
            "exclude upstream newgene auxiliary column"
        )
        standard_rows: list[dict[str, str]] = []
        malformed_pair_rows = 0
        for row in input_rows:
            if not standard_pair(row):
                continue
            p_lrt = as_float(row.get("gwas_p_lrt", ""))
            if p_lrt is None or not p_lrt < GWAS_P_LRT_THRESHOLD:
                continue
            try:
                standard_rows.append(oriented_row(row, retained_input_fields))
            except ValueError:
                malformed_pair_rows += 1
        focus_rows = [row for row in standard_rows if row.get("gwas_hit_on_ncrna") == "yes"]
        log.log(f"standard pair after filtering: {len(standard_rows)} row; GWAS match regular ncRNA: {len(focus_rows)} row; cannot orient: {malformed_pair_rows} row")

        log.log("[stage 3/5] build deduplicate relationship summary and run list")
        pair_rows = build_pair_summary(standard_rows)
        summary_rows = build_summary(standard_rows, focus_rows, len(input_rows), run_dir, input_path)
        counts = {
            "input_rows": len(input_rows),
            "standard_rows": len(standard_rows),
            "focus_rows": len(focus_rows),
            "unique_pairs": len(pair_rows),
            "malformed_pair_rows": malformed_pair_rows,
        }
        write_manifest(run_dir / "analysis_manifest.tsv", run_dir, input_path, fields, retained_input_fields, counts)

        log.log("[stage 4/5] write TSV, CSV, HTML, XLSX, XLS and SVG")
        all_fields = list(retained_input_fields) + DERIVED_FIELDS
        focus_fields = list(dict.fromkeys(field for field in FOCUS_FIELDS if field in all_fields))
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_raw_gwas_allelome_all.tsv", all_fields, standard_rows)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_raw_gwas_allelome_focus.tsv", all_fields, focus_rows)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_pair_summary.tsv", PAIR_SUMMARY_FIELDS, pair_rows)
        write_tsv(run_dir / "dog5_annotated_ncrna_pcgene_summary.tsv", SUMMARY_FIELDS, summary_rows)
        write_csv(run_dir / "dog5_annotated_ncrna_pcgene_focus_table.csv", focus_fields, focus_rows)
        write_html(
            run_dir / "dog5_annotated_ncrna_pcgene_focus_table.html",
            "Dog5 annotated ncRNA–pcGene raw GWAS–Allelome focus",
            focus_fields,
            focus_rows,
            f"Retain core GWAS/Allelome columns; p_lrt < {GWAS_P_LRT_THRESHOLD:g}; show only rows with GWAS sites within annotated ncRNA gene spans.",
        )
        write_xlsx(
            run_dir / "dog5_annotated_ncrna_pcgene_tables.xlsx",
            [
                ("focus", focus_fields, focus_rows),
                ("all", all_fields, standard_rows),
                ("pair_summary", PAIR_SUMMARY_FIELDS, pair_rows),
                ("summary", SUMMARY_FIELDS, summary_rows),
            ],
        )
        write_xls(
            run_dir / "dog5_annotated_ncrna_pcgene_tables.xls",
            [
                ("focus", focus_fields, focus_rows),
                ("pair_summary", PAIR_SUMMARY_FIELDS, pair_rows),
                ("summary", SUMMARY_FIELDS, summary_rows),
            ],
        )
        figure_fields = [field for field in [
            "trait", "gwas_marker_id", "gwas_p_lrt", "ncrna_gene_symbol", "ncrna_interval", "pcgene_gene_symbol", "mechanism", "linkage_score", "sample_tissue", "interval_gap_bp", "gwas_hit_relation", "evidence_stratum"
        ] if field in focus_fields]
        write_svg_table(
            run_dir / "dog5_annotated_ncrna_pcgene_focus_table.svg",
            "Dog5 conventional annotated ncRNA–pcGene GWAS overlap",
            figure_fields,
            sorted(focus_rows, key=lambda row: as_float(row.get("gwas_p_lrt", "")) or float("inf")),
            f"p_lrt < {GWAS_P_LRT_THRESHOLD:g}; focus rows are GWAS hits inside ncRNA spans; SVG displays up to {SVG_ROW_LIMIT} rows.",
        )
        write_svg_pair_barplot(
            run_dir / "dog5_annotated_ncrna_pcgene_pair_barplot.svg",
            pair_rows,
            "Gray: all expanded raw-GWAS×Allelome rows for this ncRNA–pcGene pair; orange: rows with a GWAS site within the ncRNA gene span; pairs are deduplicated by gene ID.",
        )

        log.log("[stage 5/5] file writing completed")
        log.log("GWAS/Allelome core field in all/focus TSV retained in; upstream newgene auxiliary column already exclude; regular ncRNA–pcGene relationship by unordered pair normalization.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
