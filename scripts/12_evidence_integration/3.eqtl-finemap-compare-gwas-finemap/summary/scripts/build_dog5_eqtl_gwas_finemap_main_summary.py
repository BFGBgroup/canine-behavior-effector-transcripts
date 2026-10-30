#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import signal
import traceback
from collections import defaultdict
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple
from xml.sax.saxutils import escape as xml_escape


DEFAULT_INPUT_RUN_DIR = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/1.1_compare_dog5_eqtl_gwas_finemap/20260910_131929'
)
DEFAULT_OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/summary'
)
SCRIPT_DIR = Path(
    '/path/to/repository/scripts/12_evidence_integration/3.eqtl-finemap-compare-gwas-finemap/summary/scripts'
)
LOG_PREFIX = "dog5_eqtl_gwas_finemap_main_summary"
BRANCH = "main"
INPUT_COMPLETION_NAME = "execution_completion.marker"
BRANCH_COMPARISON_NAME = "dog5_expression_ase_gwas_finemap_branch_comparison.tsv"
INPUT_NAMES = {
    "priority": "dog5_expression_ase_ppa_ge_0.1_gwas_both_credible_sets_exact_overlap.tsv",
    "full": "dog5_expression_ase_gwas_both_95pct_credible_sets_exact_overlap.tsv",
    "noncredible": "dog5_expression_ase_ppa_ge_0.1_gwas_noncredible_locus_variant_overlap.tsv",
    "overlap_summary": "dog5_expression_ase_gwas_finemap_overlap_summary.tsv",
    "trait_summary": "dog5_expression_ase_gwas_finemap_trait_summary.tsv",
}

COMMON_REQUIRED_FIELDS = {
    "result_layer",
    "normalized_variant_key",
    "eqtl_marker_id",
    "chrom",
    "pos",
    "eqtl_ref",
    "eqtl_alt",
    "eqtl_ambiguous_at_cg",
    "eqtl_block_id",
    "eqtl_tissue",
    "eqtl_feature_id",
    "eqtl_feature_label",
    "eqtl_feature_origin",
    "eqtl_expression_gene_id",
    "eqtl_expression_matrix_scope",
    "eqtl_ppa",
    "eqtl_ppa_tier",
    "gwas_analysis_scope",
    "gwas_branch",
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
    "gwas_variant_count",
    "gwas_ambiguous_at_cg",
    "gwas_harmonization_status",
    "gwas_locus_status",
    "gwas_model_status",
    "allele_match_status",
    "interpretation_note",
}

PRIORITY_DERIVED_COLUMNS = [
    "gwas_credible_set_min_p_snp",
    "gwas_credible_set_min_p_snp_coordinate",
    "gwas_credible_set_min_p_lrt",
]


DISPLAY_COLUMNS = [
    "gwas_trait",
    "gwas_ld_variant_id",
    "gwas_p_lrt",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_credible_sets",
    "gwas_credible_set_min_p_snp",
    "gwas_credible_set_min_p_lrt",
    "eqtl_tissue",
    "eqtl_feature_label",
    "eqtl_marker_id",
    "eqtl_ppa",
]

DISPLAY_LABELS = {
    "gwas_trait": "GWAS Trait",
    "gwas_ld_variant_id": "Matched GWAS SNP",
    "gwas_p_lrt": "Matched SNP P LRT",
    "gwas_pip": "Matched SNP PIP",
    "gwas_pip_rank": "Locus PIP Rank (rank/total)",
    "gwas_credible_sets": "GWAS CS",
    "gwas_in_credible_set": "In Credible Set",
    "gwas_credible_set_min_p_snp": "Credible-set Min-p SNP",
    "gwas_credible_set_min_p_lrt": "Credible-set Min P LRT",
    "normalized_variant_key": "Variant",
    "eqtl_tissue": "Expression+ASE Tissue",
    "eqtl_feature_label": "Expression/ASE Feature",
    "eqtl_marker_id": "Expression+ASE Fine-map SNP",
    "eqtl_ppa": "Expression+ASE PPA",
    "eqtl_ppa_tier": "Expression+ASE PPA Tier",
    "interpretation_note": "Interpretation",
}


DISPLAY_COLUMN_SPECS = [
    ("gwas_trait", 150, 22, "text_left"),
    ("gwas_ld_variant_id", 285, 44, "text_left"),
    ("gwas_p_lrt", 100, 14, "numeric"),
    ("gwas_pip", 105, 14, "numeric"),
    ("gwas_pip_rank", 100, 14, "numeric"),
    ("gwas_credible_sets", 90, 12, "text_left"),
    ("gwas_credible_set_min_p_snp", 285, 44, "text_left"),
    ("gwas_credible_set_min_p_lrt", 105, 14, "numeric"),
    ("eqtl_tissue", 125, 18, "text_left"),
    ("eqtl_feature_label", 160, 20, "text_left"),
    ("eqtl_marker_id", 190, 24, "text_left"),
    ("eqtl_ppa", 105, 14, "numeric"),
]

DISPLAY_GROUPS = [
    ("GWAS fine-mapping evidence", 0, 8, "#e8f0ea"),
    ("expression+ASE fine-mapping evidence", 8, 12, "#eef2f7"),
]

PRIORITY_EVIDENCE_COLUMNS = [
    "result_layer",
    "normalized_variant_key",
    "eqtl_marker_id",
    "eqtl_block_id",
    "eqtl_tissue",
    "eqtl_feature_id",
    "eqtl_feature_label",
    "eqtl_feature_origin",
    "eqtl_expression_matrix_scope",
    "eqtl_is_top_variant",
    "eqtl_is_credible_set_member",
    "eqtl_ppa",
    "eqtl_ppa_tier",
    "eqtl_ppa_status",
    "eqtl_marker_rank_by_ppa",
    "eqtl_block_max_ppa",
    "eqtl_block_top_ppa_tie_count",
    "eqtl_block_top_marker_is_unique",
    "eqtl_block_credible_set_size",
    "eqtl_block_credible_set_ppa_sum",
    "eqtl_selected_sample_count",
    "eqtl_common_marker_count",
    "eqtl_total_reads",
    "eqtl_input_interpretation_note",
    "gwas_analysis_scope",
    "gwas_trait",
    "gwas_locus_id",
    "gwas_ld_variant_id",
    "gwas_ld_a1",
    "gwas_ld_a2",
    "gwas_p_lrt",
    "gwas_pip",
    "gwas_pip_rank",
    "gwas_variant_count",
    "gwas_credible_sets",
    "gwas_in_credible_set",
    "gwas_credible_set_min_p_snp",
    "gwas_credible_set_min_p_snp_coordinate",
    "gwas_credible_set_min_p_lrt",
    "gwas_lead_variant_id",
    "gwas_lead_p_lrt",
    "gwas_beta_for_ld_a1",
    "gwas_se",
    "gwas_z_for_ld_a1",
    "gwas_harmonization_status",
    "allele_match_status",
    "interpretation_note",
]

GWAS_MEMBER_REQUIRED_FIELDS = {
    "locus_id",
    "trait",
    "analysis_scope",
    "branch",
    "cs_index",
    "cs_id",
    "matrix_variant_index",
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
    "branch",
    "analysis_scope",
    "trait",
    "locus_id",
    "cs_index",
    "cs_id",
    "matrix_variant_index",
    "chrom",
    "pos",
    "ld_variant_id",
    "ld_a1",
    "ld_a2",
    "p_lrt",
    "pip",
    "z_for_ld_a1",
    "p_lrt_lt_1e-5",
    "credible_snp_key",
]

GWAS_CREDIBLE_SUMMARY_COLUMNS = [
    "scope",
    "gwas_locus_count",
    "credible_set_count",
    "credible_set_membership_rows",
    "unique_credible_set_snps",
    "duplicate_membership_rows",
    "p_lrt_lt_1e-5_membership_rows",
    "p_lrt_lt_1e-5_unique_snps",
    "min_gwas_p_lrt",
    "min_gwas_p_lrt_snp",
    "min_gwas_p_lrt_snp_coordinate",
    "min_gwas_p_lrt_chrom",
    "min_gwas_p_lrt_pos",
    "min_gwas_p_lrt_trait",
    "min_gwas_p_lrt_locus",
    "min_gwas_p_lrt_pip",
    "min_gwas_p_lrt_z_for_ld_a1",
]

TRAIT_SUMMARY_COLUMNS = [
    "result_layer",
    "gwas_analysis_scope",
    "gwas_trait",
    "overlap_rows",
    "unique_variants",
    "unique_eqtl_blocks",
    "unique_eqtl_features",
    "unique_eqtl_tissues",
    "unique_gwas_loci",
    "min_gwas_p_lrt",
    "max_eqtl_ppa",
    "max_gwas_pip",
]

TISSUE_SUMMARY_COLUMNS = [
    "result_layer",
    "eqtl_tissue",
    "overlap_rows",
    "unique_variants",
    "unique_eqtl_blocks",
    "unique_eqtl_features",
    "unique_gwas_traits",
    "unique_gwas_loci",
    "primary_rows",
    "exploratory_rows",
]

SUMMARY_COLUMNS = ["metric", "value"]
MANIFEST_COLUMNS = ["role", "path", "exists", "size_bytes", "rows", "schema_status"]
ALLOWED_FEATURE_ORIGINS = {"reference_annotation_gene", "dog5_manual_new_gene"}
EXPECTED_EXPRESSION_MATRIX_SCOPE = "reference_annotation_plus_dog5_manual_new_gene"


def log(message: str, log_handle=None) -> None:
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ] {message}"
    print(line, flush=True)
    if log_handle is not None:
        log_handle.write(line + "\n")
        log_handle.flush()


@contextmanager
def execution_log(run_id: str, run_dir: Path):

    log_dir = SCRIPT_DIR / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{LOG_PREFIX}.{run_id}.log"
    success_path = log_dir / f"{LOG_PREFIX}.{run_id}.success.tsv"
    failed_path = log_dir / f"{LOG_PREFIX}.{run_id}.failed.tsv"
    for path in (log_path, success_path, failed_path):
        if path.exists():
            raise FileExistsError(f"refuse to overwrite existing execution record: {path}")
    fields = [
        "run_id", "task", "status", "exit_code", "started_at", "finished_at",
        "output_path", "task_log", "reason",
    ]
    started_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def stop_on_signal(signum, _frame):
        raise SystemExit(128 + signum)


    with log_path.open("x", encoding="utf-8", buffering=1) as handle:
        with redirect_stdout(handle), redirect_stderr(handle):
            write_tsv(success_path, fields, [])
            write_tsv(failed_path, fields, [])
            previous_handlers = {
                signum: signal.signal(signum, stop_on_signal)
                for signum in (signal.SIGINT, signal.SIGTERM)
            }
            status_row = {
                "run_id": run_id, "task": LOG_PREFIX, "started_at": started_at,
                "output_path": str(run_dir), "task_log": str(log_path),
            }
            completion_path = run_dir / "execution_completion.marker"
            completion_started = False
            log(f"[start] run_id={run_id}; output_dir={run_dir}; master_log={log_path}")
            try:

                yield None

                finished_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                write_tsv(success_path, fields, [{
                    **status_row, "status": "success", "exit_code": 0,
                    "finished_at": finished_at, "reason": "normal_end",
                }])
                completion_started = True
                completion_path.write_text(
                    f"run_id={run_id}\nexecution_status=success\nfinished_at={finished_at}\n",
                    encoding="utf-8",
                )
                log("[finish] total=1; success=1; failed=0")
                log("task completed / done; execution_status=success")
                log("reminder: execution completion does not establish output completeness, format, count and biological quality already passed, still requires separately review.")
            except BaseException as exc:
                code = exc.code if isinstance(exc, SystemExit) and isinstance(exc.code, int) else 1
                code = code or 1
                reason = str(exc).replace("\t", " ").replace("\n", " ")
                traceback.print_exc()


                if completion_started:
                    completion_path.unlink(missing_ok=True)
                write_tsv(success_path, fields, [])
                write_tsv(failed_path, fields, [{
                    **status_row, "status": "failed", "exit_code": code,
                    "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "reason": reason,
                }])
                log(f"execution_status=failed; total=1; success=0; failed=1; exit_code={code}")
                if isinstance(exc, SystemExit):
                    raise SystemExit(code) from exc
                raise
            finally:
                for signum, handler in previous_handlers.items():
                    signal.signal(signum, handler)


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        return fieldnames, list(reader)


def read_gzip_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        return fieldnames, list(reader)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def display_value(
    field: str,
    value: object,
    row: Mapping[str, object] | None = None,
) -> str:

    text = str(value if value not in (None, "") else "NA")
    if field == "gwas_trait" and row is not None:
        return f"{text}; [{row.get('gwas_analysis_scope', 'NA')}]"
    if field == "gwas_ld_variant_id" and row is not None:

        coordinate = str(row.get("normalized_variant_key", "")).strip()
        if coordinate:
            return f"{text} | {coordinate}"
    if field == "gwas_credible_set_min_p_snp" and row is not None:
        coordinate = str(
            row.get("gwas_credible_set_min_p_snp_coordinate", "")
        ).strip()
        if coordinate:
            return f"{text} | {coordinate}"
    if field == "gwas_pip_rank" and row is not None:

        total = str(row.get("gwas_variant_count", "")).strip()
        if total:
            return f"{text}/{total}"
    if field == "gwas_in_credible_set":
        return {"1": "Yes", "0": "No"}.get(text, text)
    if field == "allele_match_status" and text == "allele_set_match":
        return "Match"
    return text


def build_display_summary(
    rows: Sequence[Mapping[str, str]],
) -> Dict[str, str]:

    summary = {
        "rows": str(len(rows)),
        "unique_variants": str(len({row.get("normalized_variant_key", "") for row in rows if row.get("normalized_variant_key")})),
        "unique_traits": str(len({row.get("gwas_trait", "") for row in rows if row.get("gwas_trait")})),
        "unique_tissues": str(len({row.get("eqtl_tissue", "") for row in rows if row.get("eqtl_tissue")})),
        "ambiguous_rows": str(
            sum(
                row.get("eqtl_ambiguous_at_cg") == "yes"
                or row.get("gwas_ambiguous_at_cg") == "yes"
                for row in rows
            )
        ),
    }
    return summary


def write_html(
    path: Path,
    title: str,
    fieldnames: Sequence[str],
    rows: Sequence[Mapping[str, object]],
    summary: Mapping[str, str],
) -> None:

    group_header = []
    for label, start, end, _fill in DISPLAY_GROUPS:
        group_header.append(
            f"<th colspan='{end - start}'>{html.escape(label)}</th>"
        )
    head = "".join(
        f"<th>{html.escape(DISPLAY_LABELS.get(field, field))}</th>"
        for field in fieldnames
    )
    body = []
    previous_group: Tuple[str, str] | None = None
    for row in rows:
        group = (str(row.get("gwas_locus_id", "")), str(row.get("gwas_trait", "")))
        row_class = " group-break" if previous_group is not None and group != previous_group else ""
        cells = "".join(
            f"<td>{html.escape(display_value(field, row.get(field, ''), row))}</td>"
            for field in fieldnames
        )
        body.append(f"<tr class='{row_class.strip()}'>{cells}</tr>")
        previous_group = group
    css = """
body{font-family:Arial,sans-serif;margin:24px;background:#f7f7f4;color:#222}
h1{font-size:22px;margin:0 0 8px}
p.meta{margin:0 0 16px;color:#555;font-size:13px}
.table-wrap{max-width:100%;overflow-x:auto;border:1px solid #cfcfcf;background:#fff}
table{border-collapse:collapse;min-width:1800px;width:100%;font-size:12px;background:#fff}
th,td{border:1px solid #cfcfcf;padding:6px 8px;text-align:left;vertical-align:top}
th{background:#e6efe8;font-weight:700;position:sticky;top:0;z-index:1}
tr.group-head th{background:#dce7df;text-align:center;top:0}
thead tr:nth-child(2) th{top:29px}
tr:nth-child(even){background:#fafafa}
tr.group-break td{border-top:2px solid #777}
"""
    content = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        f"<title>{html.escape(title)}</title>",
        f"<style>{css}</style></head><body>",
        f"<h1>{html.escape(title)}</h1>",
        (
            "<p class='meta'>"
            f"selected rows = {html.escape(summary.get('rows', 'NA'))}; "
            f"unique variants = {html.escape(summary.get('unique_variants', 'NA'))}; "
            f"traits = {html.escape(summary.get('unique_traits', 'NA'))}; "
            f"tissues = {html.escape(summary.get('unique_tissues', 'NA'))}; "
            f"ambiguous rows = {html.escape(summary.get('ambiguous_rows', 'NA'))}; "
            f"branch = {html.escape(BRANCH)}."
            " Rows are retained without biological-relation aggregation."
            " Locus PIP rank is not within-CS rank. Min-p SNP is from the matched GWAS CS."
            " No additional p-value filter; shared candidates, not proof of colocalization or causality."
            " Frequency-resolved A/T or C/G variants remain in main; strict results are not pooled."
            "</p>"
        ),
        "<div class='table-wrap'><table>",
        "<thead><tr class='group-head'>",
        "".join(group_header),
        f"</tr><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody>",
        "</table></div>",
        "</body></html>",
    ]
    path.write_text("\n".join(content), encoding="utf-8")


def wrap_text(value: object, width: int) -> List[str]:

    text = str(value).strip() or "NA"
    chunks: List[str] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        while len(part) > width:
            chunks.append(part[:width])
            part = part[width:]
        if part:
            chunks.append(part)
    return chunks or ["NA"]


def svg_badge_style(field: str, value: str) -> Tuple[str, str]:

    if field == "gwas_in_credible_set":
        if value == "Yes":
            return "#e7f3ec", "#1f5a3a"
        return "#f1f3f5", "#5c6773"
    if field == "allele_match_status":
        if value == "Match":
            return "#e7f3ec", "#1f5a3a"
        return "#fde8e7", "#9b2c2c"
    if field == "gwas_harmonization_status":
        if "AMBIGUOUS" in value:
            return "#fff1e6", "#9a4d00"
        if "HARMONIZED" in value:
            return "#e7f3ec", "#1f5a3a"
    return "#f1f3f5", "#5c6773"


def write_svg_table(
    path: Path,
    *,
    title: str,
    rows: Sequence[Mapping[str, str]],
    summary: Mapping[str, str],
) -> None:

    left_margin = 30
    right_margin = 30
    top_margin = 36
    title_height = 88
    group_header_height = 24
    header_height = 38
    line_height = 12
    row_padding = 16
    column_specs = DISPLAY_COLUMN_SPECS
    table_width = sum(spec[1] for spec in column_specs)
    width = left_margin + table_width + right_margin

    x_positions: List[int] = []
    current_x = left_margin
    for _key, column_width, _wrap_width, _align in column_specs:
        x_positions.append(current_x)
        current_x += column_width

    row_line_map: List[List[List[str]]] = []
    row_heights: List[int] = []
    for row in rows:
        wrapped_cells: List[List[str]] = []
        max_lines = 1
        for key, _column_width, wrap_width, align in column_specs:
            if align == "badge":
                lines = [display_value(key, row.get(key, ""), row)]
            else:
                lines = wrap_text(display_value(key, row.get(key, ""), row), wrap_width)
            wrapped_cells.append(lines)
            max_lines = max(max_lines, len(lines))
        row_line_map.append(wrapped_cells)
        row_heights.append(row_padding + line_height * max_lines)

    height = (
        top_margin
        + title_height
        + group_header_height
        + header_height
        + sum(row_heights)
        + 22
    )
    subtitle_one = "main branch; PPA >= 0.1 and expression+ASE/GWAS credible-set overlap"
    subtitle_two = (
        f"selected rows={summary.get('rows', 'NA')}; "
        f"unique variants={summary.get('unique_variants', 'NA')}; "
        f"ambiguous rows={summary.get('ambiguous_rows', 'NA')}"
    )
    subtitle_three = (
        "Rank: all SNPs in the same locus, not within-CS rank. "
        "Min-p SNP: matched GWAS CS. No additional p-value filter; not colocalization evidence."
    )
    svg = [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="100%" '
            f'height="{height}" viewBox="0 0 {width} {height}" '
            'preserveAspectRatio="xMidYMin meet" '
            'style="display:block;max-width:100%;height:auto">'
        ),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        (
            f'<text x="{width / 2:.1f}" y="{top_margin}" font-size="16" '
            'font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" '
            'text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">'
            f"{xml_escape(title)}</text>"
        ),
        (
            f'<text x="{width / 2:.1f}" y="{top_margin + 18}" font-size="11.5" '
            'font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" '
            'text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">'
            f"{xml_escape(subtitle_one)}</text>"
        ),
        (
            f'<text x="{width / 2:.1f}" y="{top_margin + 36}" font-size="11.5" '
            'font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" '
            'text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">'
            f"{xml_escape(subtitle_two)}</text>"
        ),
        (
            f'<text x="{width / 2:.1f}" y="{top_margin + 54}" font-size="11.5" '
            'font-family="Helvetica, Arial, sans-serif" font-weight="400" fill="#111111" '
            'text-anchor="middle" dominant-baseline="middle" font-variant-numeric="tabular-nums">'
            f"{xml_escape(subtitle_three)}</text>"
        ),
    ]

    header_y = top_margin + title_height
    svg.append(
        f'<line x1="{left_margin}" y1="{header_y}" '
        f'x2="{left_margin + table_width}" y2="{header_y}" '
        'stroke="#111111" stroke-width="1.1"/>'
    )
    for group_label, start, end, fill in DISPLAY_GROUPS:
        group_x = x_positions[start]
        group_width = sum(column_specs[idx][1] for idx in range(start, end))
        svg.append(
            f'<rect x="{group_x}" y="{header_y}" width="{group_width}" '
            f'height="{group_header_height}" fill="{fill}"/>'
        )
        svg.append(
            f'<text x="{group_x + group_width / 2:.1f}" y="{header_y + 16}" '
            'font-size="10.5" font-family="Helvetica, Arial, sans-serif" '
            'font-weight="700" fill="#111111" text-anchor="middle">'
            f"{xml_escape(group_label)}</text>"
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
    for key, column_width, wrap_width, _align in column_specs:
        column_index = [spec[0] for spec in column_specs].index(key)
        x = x_positions[column_index]
        lines = wrap_text(DISPLAY_LABELS.get(key, key), wrap_width)
        for line_index, line in enumerate(lines[:2]):
            svg.append(
                f'<text x="{x + column_width / 2:.1f}" '
                f'y="{column_header_y + 14 + line_index * line_height}" '
                'font-size="10" font-family="Helvetica, Arial, sans-serif" '
                'font-weight="700" fill="#111111" text-anchor="middle">'
                f"{xml_escape(line)}</text>"
            )

    current_y = column_header_y + header_height
    previous_group: Tuple[str, str] | None = None
    visual_group_index = -1
    for row, wrapped_cells, row_height in zip(rows, row_line_map, row_heights):
        group = (row.get("gwas_locus_id", ""), row.get("gwas_trait", ""))
        if group != previous_group:
            visual_group_index += 1
            if previous_group is not None:
                svg.append(
                    f'<line x1="{left_margin}" y1="{current_y}" '
                    f'x2="{left_margin + table_width}" y2="{current_y}" '
                    'stroke="#777777" stroke-width="1.1"/>'
                )
            previous_group = group
        if visual_group_index % 2 == 1:
            svg.append(
                f'<rect x="{left_margin}" y="{current_y}" width="{table_width}" '
                f'height="{row_height}" fill="#fafafa"/>'
            )
        for column_index, (key, column_width, _wrap_width, align) in enumerate(column_specs):
            x = x_positions[column_index]
            value = display_value(key, row.get(key, ""), row)
            lines = wrapped_cells[column_index]
            if align == "badge":
                badge_fill, badge_text = svg_badge_style(key, value)
                badge_width = min(column_width - 16, max(74, len(value) * 5.8 + 18))
                badge_x = x + 8
                badge_y = current_y + 8
                svg.append(
                    f'<rect x="{badge_x}" y="{badge_y}" width="{badge_width}" '
                    f'height="18" rx="9" fill="{badge_fill}"/>'
                )
                svg.append(
                    f'<text x="{badge_x + badge_width / 2:.1f}" y="{badge_y + 12}" '
                    'font-size="9" font-family="Arial" font-weight="700" '
                    f'fill="{badge_text}" text-anchor="middle">{xml_escape(value)}</text>'
                )
                continue
            text_x = x + (8 if align == "text_left" else column_width / 2)
            text_start_y = current_y + (row_height - line_height * len(lines)) / 2 + line_height * 0.75
            for line_index, line in enumerate(lines):
                text_y = text_start_y + line_index * line_height
                if align == "numeric":
                    svg.append(
                        f'<text x="{text_x:.1f}" y="{text_y:.1f}" font-size="9.4" '
                        'font-family="Menlo,Consolas,monospace" fill="#111111" '
                        'text-anchor="middle" dominant-baseline="middle">'
                        f"{xml_escape(line)}</text>"
                    )
                else:
                    svg.append(
                        f'<text x="{text_x:.1f}" y="{text_y:.1f}" font-size="9.5" '
                        'font-family="Helvetica, Arial, sans-serif" fill="#111111" '
                        'text-anchor="start" dominant-baseline="middle">'
                        f"{xml_escape(line)}</text>"
                    )
        current_y += row_height
        svg.append(
            f'<line x1="{left_margin}" y1="{current_y}" '
            f'x2="{left_margin + table_width}" y2="{current_y}" '
            'stroke="#dddddd" stroke-width="0.7"/>'
        )

    for _label, start, _end, _fill in DISPLAY_GROUPS[1:]:
        boundary_x = x_positions[start]
        svg.append(
            f'<line x1="{boundary_x}" y1="{header_y}" '
            f'x2="{boundary_x}" y2="{height - 22}" '
            'stroke="#888888" stroke-width="0.8"/>'
        )
    svg.append(
        f'<line x1="{left_margin}" y1="{height - 22}" '
        f'x2="{left_margin + table_width}" y2="{height - 22}" '
        'stroke="#111111" stroke-width="1.0"/>'
    )
    svg.append("</svg>")
    path.write_text("\n".join(svg), encoding="utf-8")


def validate_rows(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, str]], expected_layer: str) -> None:
    if not rows:
        raise ValueError(f"input table is empty: {path}")
    missing = sorted(COMMON_REQUIRED_FIELDS - set(fieldnames))
    if missing:
        raise ValueError(f"{path} missing fields: {', '.join(missing)}")
    wrong_layer = sorted({row.get("result_layer", "") for row in rows if row.get("result_layer") != expected_layer})
    if wrong_layer:
        raise ValueError(f"{path} contains non-expected result_layer: {wrong_layer}")
    wrong_branch = sorted({row.get("gwas_branch", "") for row in rows if row.get("gwas_branch") != BRANCH})
    if wrong_branch:
        raise ValueError(f"{path} is not main branch: {wrong_branch}")
    if any(not row.get("normalized_variant_key") for row in rows):
        raise ValueError(f"{path} contains empty normalized_variant_key")
    invalid_origins = sorted(
        {
            row.get("eqtl_feature_origin", "")
            for row in rows
            if row.get("eqtl_feature_origin", "") not in ALLOWED_FEATURE_ORIGINS
        }
    )
    if invalid_origins:
        raise ValueError(f"{path} contains non-expected eqtl_feature_origin: {invalid_origins}")
    invalid_scopes = sorted(
        {
            row.get("eqtl_expression_matrix_scope", "")
            for row in rows
            if row.get("eqtl_expression_matrix_scope", "")
            != EXPECTED_EXPRESSION_MATRIX_SCOPE
        }
    )
    if invalid_scopes:
        raise ValueError(
            f"{path} contains non-expected eqtl_expression_matrix_scope: {invalid_scopes}"
        )


def read_metric_table(path: Path) -> Tuple[Dict[str, str], int]:
    fieldnames, rows = read_tsv(path)
    if fieldnames[:2] != ["metric", "value"]:
        raise ValueError(f"summary table field is not metric/value: {path}")
    return ({row.get("metric", ""): row.get("value", "") for row in rows if row.get("metric")}, len(rows))


def read_completion_marker(path: Path) -> Dict[str, str]:

    values: Dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        if "=" not in line:
            raise ValueError(f"completion marker contains non-key=value row: {path}: {line}")
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    if values.get("execution_status") != "success":
        raise ValueError(
            f"upstream completion marker is not execution_status=success: {path}"
        )
    return values


def read_trait_summary(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    fieldnames, rows = read_tsv(path)
    required = set(TRAIT_SUMMARY_COLUMNS)
    missing = sorted(required - set(fieldnames))
    if missing:
        raise ValueError(f"upstream trait summary missing fields: {', '.join(missing)}")
    wrong_branch = sorted({row.get("gwas_branch", "") for row in rows if row.get("gwas_branch") != BRANCH})
    if wrong_branch:
        raise ValueError(f"upstream trait summary contains non-main branch: {wrong_branch}")
    return fieldnames, rows


def read_branch_comparison(path: Path) -> Tuple[Dict[str, str], int]:
    fieldnames, rows = read_tsv(path)
    required = {
        "gwas_branch",
        "gwas_locus_count",
        "gwas_credible_rows",
        "priority_both_overlap_rows",
        "full_both_overlap_rows",
    }
    missing = sorted(required - set(fieldnames))
    if missing:
        raise ValueError(f"branch comparison missing fields: {', '.join(missing)}")
    matches = [row for row in rows if row.get("gwas_branch") == BRANCH]
    if len(matches) != 1:
        raise ValueError(f"branch comparison in main row count is not 1: {len(matches)}")
    return matches[0], len(rows)


def read_gwas_credible_members(
    gwas_run_dir: Path,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:

    manifest_path = gwas_run_dir / "locus_run_manifest.tsv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"GWAS locus manifest does not exist: {manifest_path}")
    manifest_fields, manifest_rows = read_tsv(manifest_path)
    required_manifest = {
        "locus_id",
        "trait",
        "analysis_scope",
        "branch",
        "status",
        "model_status",
        "credible_set_count",
        "credible_set_members",
    }
    missing_manifest = sorted(required_manifest - set(manifest_fields))
    if missing_manifest:
        raise ValueError(
            f"{manifest_path} missing credible set field: {', '.join(missing_manifest)}"
        )
    main_manifest = [row for row in manifest_rows if row.get("branch") == BRANCH]
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
    if not main_manifest:
        raise ValueError(f"{manifest_path} contains no main branch locus")

    members: List[Dict[str, str]] = []
    for manifest_row in main_manifest:
        member_path = Path(manifest_row.get("credible_set_members", ""))
        if not member_path.is_absolute() or not member_path.is_file():
            raise FileNotFoundError(
                f"GWAS credible_set_members does not exist or is not an absolute path: {member_path}"
            )
        member_fields, member_rows = read_gzip_tsv(member_path)
        missing_member = sorted(GWAS_MEMBER_REQUIRED_FIELDS - set(member_fields))
        if missing_member:
            raise ValueError(
                f"{member_path} missing credible set member field: {', '.join(missing_member)}"
            )
        for member_row in member_rows:
            for key in ("locus_id", "trait", "analysis_scope", "branch"):
                if member_row.get(key, "") != manifest_row.get(key, ""):
                    raise ValueError(
                        f"{member_path} of {key} with manifest mismatch: "
                        f"{member_row.get(key, '')}!={manifest_row.get(key, '')}"
                    )
            try:
                p_lrt = float(member_row.get("p_lrt", ""))
            except ValueError as exc:
                raise ValueError(f"{member_path} contains non-numeric value p_lrt") from exc
            if not math.isfinite(p_lrt):
                raise ValueError(f"{member_path} contains non-finite p_lrt")
            chrom = member_row.get("chrom", "")
            pos = member_row.get("pos", "")
            ld_variant_id = member_row.get("ld_variant_id", "")
            ld_a1 = member_row.get("ld_a1", "").strip().upper()
            ld_a2 = member_row.get("ld_a2", "").strip().upper()
            if not chrom or not pos or not ld_variant_id or not ld_a1 or not ld_a2:
                raise ValueError(
                    f"{member_path} contains empty of chrom/pos/ld_variant_id/ld_a1/ld_a2"
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
                    "ld_variant_id": ld_variant_id,
                    "ld_a1": ld_a1,
                    "ld_a2": ld_a2,
                    "p_lrt": member_row.get("p_lrt", ""),
                    "pip": member_row.get("pip", ""),
                    "z_for_ld_a1": member_row.get("z_for_ld_a1", ""),
                    "p_lrt_lt_1e-5": "yes" if p_lrt < 1e-5 else "no",
                    "credible_snp_key": f"{chrom}:{pos}:{ld_variant_id}",
                }
            )
    return main_manifest, members


def format_gwas_credible_coordinate(row: Mapping[str, str]) -> str:

    chrom = str(row.get("chrom", "")).strip()
    pos = str(row.get("pos", "")).strip()
    alleles = sorted(
        (
            str(row.get("ld_a1", "")).strip().upper(),
            str(row.get("ld_a2", "")).strip().upper(),
        )
    )
    if not chrom or not pos or len(alleles) != 2 or not all(alleles):
        return "NA"
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return f"{chrom}:{pos}:{alleles[0]}/{alleles[1]}"


def enrich_priority_display_rows(
    rows: Sequence[Mapping[str, str]],
    members: Sequence[Mapping[str, str]],
) -> List[Dict[str, str]]:

    member_index: Dict[Tuple[str, str], List[Mapping[str, str]]] = defaultdict(list)
    for member in members:
        member_index[(member.get("locus_id", ""), member.get("cs_id", ""))].append(member)

    enriched: List[Dict[str, str]] = []
    for row in rows:
        credible_set_ids = [
            part.strip()
            for part in str(row.get("gwas_credible_sets", ""))
            .replace(",", ";")
            .split(";")
            if part.strip()
        ]
        candidates: List[Mapping[str, str]] = []
        for cs_id in credible_set_ids:
            candidates.extend(
                member_index.get((row.get("gwas_locus_id", ""), cs_id), [])
            )
        if not candidates:
            raise ValueError(
                "priority row cannot find corresponding credible set member: "
                f"locus={row.get('gwas_locus_id', '')}; "
                f"credible_sets={row.get('gwas_credible_sets', '')}"
            )
        minimum = min(candidates, key=lambda member: float(member.get("p_lrt", "inf")))
        output_row = dict(row)
        output_row["gwas_credible_set_min_p_snp"] = minimum.get(
            "ld_variant_id", "NA"
        )
        output_row["gwas_credible_set_min_p_snp_coordinate"] = (
            format_gwas_credible_coordinate(minimum)
        )
        output_row["gwas_credible_set_min_p_lrt"] = minimum.get("p_lrt", "NA")
        enriched.append(output_row)
    return enriched


def summarize_gwas_credible_members(
    manifest_rows: Sequence[Mapping[str, str]],
    members: Sequence[Mapping[str, str]],
) -> List[Dict[str, object]]:

    output: List[Dict[str, object]] = []
    scopes = ("all_main", "primary", "exploratory")
    for scope in scopes:
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
            if row.get("credible_snp_key")
        }
        low_p_members = [
            row for row in selected_members if row.get("p_lrt_lt_1e-5") == "yes"
        ]
        low_p_keys = {
            row.get("credible_snp_key", "")
            for row in low_p_members
            if row.get("credible_snp_key")
        }
        min_row = min(
            selected_members,
            key=lambda row: float(row.get("p_lrt", "inf")),
            default=None,
        )
        if min_row is None:
            min_values = {
                "min_gwas_p_lrt": "NA",
                "min_gwas_p_lrt_snp": "NA",
                "min_gwas_p_lrt_snp_coordinate": "NA",
                "min_gwas_p_lrt_chrom": "NA",
                "min_gwas_p_lrt_pos": "NA",
                "min_gwas_p_lrt_trait": "NA",
                "min_gwas_p_lrt_locus": "NA",
                "min_gwas_p_lrt_pip": "NA",
                "min_gwas_p_lrt_z_for_ld_a1": "NA",
            }
        else:
            min_values = {
                "min_gwas_p_lrt": min_row.get("p_lrt", "NA"),
                "min_gwas_p_lrt_snp": min_row.get("ld_variant_id", "NA"),
                "min_gwas_p_lrt_snp_coordinate": format_gwas_credible_coordinate(min_row),
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


def unique_count(rows: Sequence[Mapping[str, str]], field: str) -> int:
    return len({row.get(field, "") for row in rows if row.get(field, "")})


def float_values(rows: Sequence[Mapping[str, str]], field: str) -> List[float]:
    values: List[float] = []
    for row in rows:
        text = row.get(field, "")
        if text in {"", "NA", "."}:
            continue
        values.append(float(text))
    return values


def sort_key(row: Mapping[str, str]) -> Tuple[str, int, str, str, str]:
    text = row.get("chrom", "")
    chrom = text[3:] if text.startswith("chr") else text
    try:
        chrom_key = f"0{int(chrom):03d}"
    except ValueError:
        chrom_key = f"1{text}"
    try:
        position = int(float(row.get("pos", "0")))
    except ValueError:
        position = 0
    return (
        chrom_key,
        position,
        row.get("gwas_trait", ""),
        row.get("eqtl_tissue", ""),
        row.get("eqtl_feature_label", ""),
    )


def association_key(row: Mapping[str, str]) -> Tuple[str, ...]:
    return tuple(
        row.get(field, "")
        for field in (
            "result_layer",
            "normalized_variant_key",
            "eqtl_block_id",
            "eqtl_feature_id",
            "eqtl_tissue",
            "gwas_analysis_scope",
            "gwas_locus_id",
            "gwas_trait",
        )
    )


def aggregate_trait(rows: Sequence[Mapping[str, str]]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str, str], List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row.get("result_layer", ""), row.get("gwas_analysis_scope", ""), row.get("gwas_trait", ""))].append(row)
    output: List[Dict[str, object]] = []
    for (layer, scope, trait), group in sorted(groups.items()):
        p_values = float_values(group, "gwas_p_lrt")
        ppa_values = float_values(group, "eqtl_ppa")
        pip_values = float_values(group, "gwas_pip")
        output.append(
            {
                "result_layer": layer,
                "gwas_analysis_scope": scope,
                "gwas_trait": trait,
                "overlap_rows": len(group),
                "unique_variants": unique_count(group, "normalized_variant_key"),
                "unique_eqtl_blocks": unique_count(group, "eqtl_block_id"),
                "unique_eqtl_features": unique_count(group, "eqtl_feature_id"),
                "unique_eqtl_tissues": unique_count(group, "eqtl_tissue"),
                "unique_gwas_loci": unique_count(group, "gwas_locus_id"),
                "min_gwas_p_lrt": min(p_values) if p_values else "NA",
                "max_eqtl_ppa": max(ppa_values) if ppa_values else "NA",
                "max_gwas_pip": max(pip_values) if pip_values else "NA",
            }
        )
    return output


def aggregate_tissue(rows: Sequence[Mapping[str, str]]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str], List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row.get("result_layer", ""), row.get("eqtl_tissue", ""))].append(row)
    output: List[Dict[str, object]] = []
    for (layer, tissue), group in sorted(groups.items()):
        output.append(
            {
                "result_layer": layer,
                "eqtl_tissue": tissue,
                "overlap_rows": len(group),
                "unique_variants": unique_count(group, "normalized_variant_key"),
                "unique_eqtl_blocks": unique_count(group, "eqtl_block_id"),
                "unique_eqtl_features": unique_count(group, "eqtl_feature_id"),
                "unique_gwas_traits": unique_count(group, "gwas_trait"),
                "unique_gwas_loci": unique_count(group, "gwas_locus_id"),
                "primary_rows": sum(row.get("gwas_analysis_scope") == "primary" for row in group),
                "exploratory_rows": sum(row.get("gwas_analysis_scope") == "exploratory" for row in group),
            }
        )
    return output


def compare_trait_summary(
    upstream_rows: Sequence[Mapping[str, str]],
    derived_rows: Sequence[Mapping[str, object]],
) -> List[str]:

    key_fields = ("result_layer", "gwas_analysis_scope", "gwas_trait")
    count_fields = (
        "overlap_rows",
        "unique_variants",
        "unique_eqtl_blocks",
        "unique_eqtl_features",
        "unique_eqtl_tissues",
        "unique_gwas_loci",
    )
    upstream_index = {
        tuple(row.get(field, "") for field in key_fields): row
        for row in upstream_rows
    }
    derived_index = {
        tuple(str(row.get(field, "")) for field in key_fields): row
        for row in derived_rows
    }
    mismatches: List[str] = []
    if set(upstream_index) != set(derived_index):
        mismatches.append(
            f"trait_summary_key_set:{len(upstream_index)}!={len(derived_index)}"
        )
    for key in sorted(set(upstream_index) & set(derived_index)):
        upstream = upstream_index[key]
        derived = derived_index[key]
        for field in count_fields:
            if str(upstream.get(field, "")) != str(derived.get(field, "")):
                mismatches.append(
                    f"trait_summary:{'|'.join(key)}:{field}:"
                    f"{upstream.get(field, '')}!={derived.get(field, '')}"
                )
    return mismatches


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize overlap between five-dog expression+ASE results and GWAS fine-mapping results in the main branch."
    )
    parser.add_argument("--input-run-dir", type=Path, default=DEFAULT_INPUT_RUN_DIR, help="Timestamped directory containing step 1.1 main output.")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT, help="Root directory for summary results.")
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S"), help="Directory name for this summary run.")
    args = parser.parse_args()
    if not args.input_run_dir.is_absolute() or not args.output_root.is_absolute():
        parser.error("--input-run-dir and --output-root must use absolute path")
    if (len(args.run_id) != 15 or args.run_id[8] != "_"
            or not (args.run_id[:8] + args.run_id[9:]).isdigit()):
        parser.error("--run-id must be YYYYMMDD_HHMMSS")
    datetime.strptime(args.run_id, "%Y%m%d_%H%M%S")

    main_dir = args.input_run_dir / BRANCH
    input_paths = {
        key: main_dir / name for key, name in INPUT_NAMES.items()
    }
    input_completion_path = args.input_run_dir / INPUT_COMPLETION_NAME
    branch_comparison_path = args.input_run_dir / BRANCH_COMPARISON_NAME
    run_dir = args.output_root / args.run_id

    with execution_log(args.run_id, run_dir) as log_handle:
        for path in [*input_paths.values(), branch_comparison_path, input_completion_path]:
            if not path.is_file():
                raise FileNotFoundError(f"input file does not exist: {path}")
        input_completion = read_completion_marker(input_completion_path)
        if input_completion.get("run_id") != args.input_run_dir.name:
            raise ValueError(
                "upstream completion marker run_id does not match the input directory name: "
                f"{input_completion.get('run_id', '')}!={args.input_run_dir.name}"
            )
        if run_dir.exists():
            raise FileExistsError(f"refuse to overwrite existing result directory: {run_dir}")
        run_dir.mkdir(parents=True, exist_ok=False)
        log(
            f"[start] input_run_dir={args.input_run_dir}; branch={BRANCH}; "
            f"input_execution_status={input_completion['execution_status']}",
            log_handle,
        )
        log(
            "[resource] concurrent_tasks=1; threads_per_task=1; "
            "baseline=preserved_serial",
            log_handle,
        )
        log("[1/4] read and validate main input table", log_handle)
        loaded: Dict[str, Tuple[List[str], List[Dict[str, str]]]] = {}
        layers = {"priority": "priority_both", "full": "full_both", "noncredible": "noncredible"}
        for key in ("priority", "full", "noncredible"):
            fields, rows = read_tsv(input_paths[key])
            validate_rows(input_paths[key], fields, rows, layers[key])
            loaded[key] = (fields, rows)
            log(f"[input] {key}_rows={len(rows)}; path={input_paths[key]}", log_handle)

        overlap_metrics, overlap_summary_rows = read_metric_table(input_paths["overlap_summary"])
        _upstream_trait_fields, upstream_trait_rows = read_trait_summary(input_paths["trait_summary"])
        branch_row, branch_comparison_rows = read_branch_comparison(branch_comparison_path)
        log(f"[input] upstream_trait_summary_rows={len(upstream_trait_rows)}", log_handle)
        gwas_run_dir_text = overlap_metrics.get("gwas_run_dir", "")
        if not gwas_run_dir_text:
            raise ValueError(
                f"{input_paths['overlap_summary']} missing gwas_run_dir, cannot read complete GWAS credible set"
            )
        gwas_run_dir = Path(gwas_run_dir_text)
        if not gwas_run_dir.is_absolute():
            raise ValueError(f"overlap summary in gwas_run_dir is not an absolute path: {gwas_run_dir}")
        gwas_manifest_rows, gwas_credible_members = read_gwas_credible_members(gwas_run_dir)
        gwas_credible_summary_rows = summarize_gwas_credible_members(
            gwas_manifest_rows,
            gwas_credible_members,
        )
        gwas_all_summary = next(
            row for row in gwas_credible_summary_rows if row.get("scope") == "all_main"
        )
        log(
            f"[input] gwas_main_loci={len(gwas_manifest_rows)}; "
            f"gwas_credible_members={len(gwas_credible_members)}; "
            f"gwas_unique_credible_snps={gwas_all_summary['unique_credible_set_snps']}",
            log_handle,
        )

        all_rows = loaded["priority"][1] + loaded["full"][1] + loaded["noncredible"][1]
        duplicate_keys = len(all_rows) - len({association_key(row) for row in all_rows})
        allele_mismatch_rows = sum(row.get("allele_match_status") != "allele_set_match" for row in all_rows)
        if duplicate_keys:
            raise ValueError(f"main input contains duplicate association key: {duplicate_keys}")
        if allele_mismatch_rows:
            raise ValueError(f"main input contains non-allele_set_match record: {allele_mismatch_rows}")

        log("[2/4] check upstream overlap summary, trait summary and branch comparison", log_handle)
        derived_counts = {
            "priority_both_credible_overlap_rows": len(loaded["priority"][1]),
            "priority_both_credible_unique_variants": unique_count(loaded["priority"][1], "normalized_variant_key"),
            "full_both_credible_overlap_rows": len(loaded["full"][1]),
            "full_both_credible_unique_variants": unique_count(loaded["full"][1], "normalized_variant_key"),
            "priority_noncredible_diagnostic_rows": len(loaded["noncredible"][1]),
        }


        upstream_mismatches: List[str] = []
        for key, expected in derived_counts.items():
            if key not in overlap_metrics:
                upstream_mismatches.append(f"overlap_summary:{key}:missing")
            elif str(expected) != str(overlap_metrics[key]):
                upstream_mismatches.append(f"overlap_summary:{key}:{overlap_metrics[key]}!={expected}")
        branch_checks = {
            "gwas_locus_count": len(gwas_manifest_rows),
            "gwas_credible_rows": len(gwas_credible_members),
            "priority_both_overlap_rows": len(loaded["priority"][1]),
            "full_both_overlap_rows": len(loaded["full"][1]),
            "priority_noncredible_diagnostic_rows": len(loaded["noncredible"][1]),
        }
        for key, expected in branch_checks.items():
            if key not in branch_row:
                upstream_mismatches.append(f"branch_comparison:{key}:missing")
            elif str(branch_row.get(key, "")) != str(expected):
                upstream_mismatches.append(f"branch_comparison:{key}:{branch_row.get(key, '')}!={expected}")
        if upstream_mismatches:
            raise ValueError("upstream summary counts differ: " + "; ".join(upstream_mismatches))

        trait_rows = aggregate_trait(all_rows)
        upstream_mismatches.extend(compare_trait_summary(upstream_trait_rows, trait_rows))
        if upstream_mismatches:
            raise ValueError("upstream summary counts differ: " + "; ".join(upstream_mismatches))
        tissue_rows = aggregate_tissue(all_rows)
        log(
            f"[3/4] generate trait/tissue and GWAS credible set summary; "
            f"trait_rows={len(trait_rows)}; tissue_rows={len(tissue_rows)}; "
            f"credible_summary_rows={len(gwas_credible_summary_rows)}",
            log_handle,
        )

        output_paths = {
            "priority_focus": run_dir / "dog5_eqtl_gwas_finemap_main_priority_focus.tsv",
            "priority_csv": run_dir / "dog5_eqtl_gwas_finemap_main_priority_focus.csv",
            "priority_html": run_dir / "dog5_eqtl_gwas_finemap_main_priority_focus.html",
            "priority_svg": run_dir / "dog5_eqtl_gwas_finemap_main_priority_focus.svg",
            "priority_evidence": run_dir / "dog5_eqtl_gwas_finemap_main_priority_evidence.tsv",
            "gwas_credible_members": run_dir / "dog5_eqtl_gwas_finemap_main_gwas_credible_set_members.tsv",
            "gwas_credible_summary": run_dir / "dog5_eqtl_gwas_finemap_main_gwas_credible_set_summary.tsv",
            "full": run_dir / "dog5_eqtl_gwas_finemap_main_full_overlap.tsv",
            "noncredible": run_dir / "dog5_eqtl_gwas_finemap_main_noncredible_diagnostic.tsv",
            "trait": run_dir / "dog5_eqtl_gwas_finemap_main_trait_summary.tsv",
            "tissue": run_dir / "dog5_eqtl_gwas_finemap_main_tissue_summary.tsv",
            "summary": run_dir / "dog5_eqtl_gwas_finemap_main_summary.tsv",
            "manifest": run_dir / "analysis_manifest.tsv",
        }
        priority_rows = sorted(loaded["priority"][1], key=sort_key)
        full_rows = sorted(loaded["full"][1], key=sort_key)
        noncredible_rows = sorted(loaded["noncredible"][1], key=sort_key)
        priority_reference_rows = [
            row
            for row in priority_rows
            if row.get("eqtl_feature_origin") == "reference_annotation_gene"
        ]
        priority_new_gene_rows = [
            row
            for row in priority_rows
            if row.get("eqtl_feature_origin") == "dog5_manual_new_gene"
        ]
        full_reference_rows = [
            row
            for row in full_rows
            if row.get("eqtl_feature_origin") == "reference_annotation_gene"
        ]
        full_new_gene_rows = [
            row
            for row in full_rows
            if row.get("eqtl_feature_origin") == "dog5_manual_new_gene"
        ]
        priority_ppa_values = float_values(priority_rows, "eqtl_ppa")
        priority_p_values = float_values(priority_rows, "gwas_p_lrt")
        priority_pip_values = float_values(priority_rows, "gwas_pip")
        priority_low_p_rows = [
            row for row in priority_rows if float(row.get("gwas_p_lrt", "inf")) < 1e-5
        ]
        priority_display_rows = enrich_priority_display_rows(
            priority_rows,
            gwas_credible_members,
        )
        display_summary = build_display_summary(priority_rows)

        focus_fields = loaded["priority"][0] + PRIORITY_DERIVED_COLUMNS
        write_tsv(output_paths["priority_focus"], focus_fields, priority_display_rows)
        write_csv(output_paths["priority_csv"], focus_fields, priority_display_rows)
        write_tsv(
            output_paths["priority_evidence"],
            PRIORITY_EVIDENCE_COLUMNS,
            priority_display_rows,
        )
        write_html(
            output_paths["priority_html"],
            "Dog5 main expression+ASE × GWAS fine-mapping priority overlap",
            DISPLAY_COLUMNS,
            priority_display_rows,
            display_summary,
        )
        write_svg_table(
            output_paths["priority_svg"],
            title="Dog5 main expression+ASE x GWAS fine-mapping priority overlap",
            rows=priority_display_rows,
            summary=display_summary,
        )
        write_tsv(output_paths["full"], loaded["full"][0], full_rows)
        write_tsv(output_paths["noncredible"], loaded["noncredible"][0], noncredible_rows)
        write_tsv(output_paths["trait"], TRAIT_SUMMARY_COLUMNS, trait_rows)
        write_tsv(output_paths["tissue"], TISSUE_SUMMARY_COLUMNS, tissue_rows)
        write_tsv(
            output_paths["gwas_credible_members"],
            GWAS_CREDIBLE_SET_MEMBER_COLUMNS,
            gwas_credible_members,
        )
        write_tsv(
            output_paths["gwas_credible_summary"],
            GWAS_CREDIBLE_SUMMARY_COLUMNS,
            gwas_credible_summary_rows,
        )

        summary_metrics = [
            ("run_id", args.run_id),
            ("input_run_dir", str(args.input_run_dir)),
            ("input_execution_status", input_completion["execution_status"]),
            ("gwas_branch", BRANCH),
            ("branch_policy", "Summarizes main only; retains upstream frequency-resolved A/T and C/G variants; does not combine strict_no_ambiguous."),
            ("layer_policy", "priority_both is a subset of full_both; tiers are counted separately and must not be added together."),
            ("gwas_pip_rank_definition", "Descending competitive PIP rank among all candidate SNPs in the same GWAS locus; the denominator is gwas_variant_count, not credible-set size."),
            ("gwas_min_p_display_scope", "focus shows only the lowest p_lrt in the GWAS credible set containing the matched SNP; the global minimum is retained only in the overall GWAS statistics table."),
            ("priority_rows", len(priority_rows)),
            ("priority_unique_variants", unique_count(priority_rows, "normalized_variant_key")),
            ("priority_reference_annotation_rows", len(priority_reference_rows)),
            ("priority_reference_annotation_unique_variants", unique_count(priority_reference_rows, "normalized_variant_key")),
            ("priority_dog5_manual_new_gene_rows", len(priority_new_gene_rows)),
            ("priority_dog5_manual_new_gene_unique_variants", unique_count(priority_new_gene_rows, "normalized_variant_key")),
            ("full_rows", len(full_rows)),
            ("full_unique_variants", unique_count(full_rows, "normalized_variant_key")),
            ("full_reference_annotation_rows", len(full_reference_rows)),
            ("full_dog5_manual_new_gene_rows", len(full_new_gene_rows)),
            ("noncredible_diagnostic_rows", len(noncredible_rows)),
            ("noncredible_unique_variants", unique_count(noncredible_rows, "normalized_variant_key")),
            ("priority_eqtl_ppa_min", min(priority_ppa_values)),
            ("priority_eqtl_ppa_max", max(priority_ppa_values)),
            ("priority_gwas_p_lrt_min", min(priority_p_values)),
            ("priority_gwas_p_lrt_max", max(priority_p_values)),
            ("priority_gwas_pip_min", min(priority_pip_values)),
            ("priority_gwas_pip_max", max(priority_pip_values)),
            ("priority_gwas_unique_snps", unique_count(priority_rows, "gwas_ld_variant_id")),
            ("priority_gwas_p_lrt_lt_1e-5_rows", len(priority_low_p_rows)),
            ("priority_gwas_p_lrt_lt_1e-5_unique_snps", unique_count(priority_low_p_rows, "gwas_ld_variant_id")),
            ("gwas_run_dir", str(gwas_run_dir)),
            ("gwas_credible_set_count", gwas_all_summary["credible_set_count"]),
            ("gwas_credible_set_membership_rows", gwas_all_summary["credible_set_membership_rows"]),
            ("gwas_unique_credible_set_snps", gwas_all_summary["unique_credible_set_snps"]),
            ("gwas_duplicate_credible_membership_rows", gwas_all_summary["duplicate_membership_rows"]),
            ("gwas_p_lrt_lt_1e-5_membership_rows", gwas_all_summary["p_lrt_lt_1e-5_membership_rows"]),
            ("gwas_p_lrt_lt_1e-5_unique_snps", gwas_all_summary["p_lrt_lt_1e-5_unique_snps"]),
            ("gwas_min_p_lrt", gwas_all_summary["min_gwas_p_lrt"]),
            ("gwas_min_p_lrt_snp", gwas_all_summary["min_gwas_p_lrt_snp"]),
            ("gwas_min_p_lrt_snp_coordinate", gwas_all_summary["min_gwas_p_lrt_snp_coordinate"]),
            ("gwas_min_p_lrt_chrom", gwas_all_summary["min_gwas_p_lrt_chrom"]),
            ("gwas_min_p_lrt_pos", gwas_all_summary["min_gwas_p_lrt_pos"]),
            ("gwas_min_p_lrt_locus", gwas_all_summary["min_gwas_p_lrt_locus"]),
            ("gwas_min_p_lrt_trait", gwas_all_summary["min_gwas_p_lrt_trait"]),
            ("gwas_min_p_lrt_pip", gwas_all_summary["min_gwas_p_lrt_pip"]),
            ("gwas_min_p_lrt_z_for_ld_a1", gwas_all_summary["min_gwas_p_lrt_z_for_ld_a1"]),
            ("allele_mismatch_rows", allele_mismatch_rows),
            ("duplicate_association_rows", duplicate_keys),
            ("upstream_allele_mismatch_qc_rows", overlap_metrics["allele_mismatch_qc_rows"]),
            ("upstream_duplicate_output_qc_rows", overlap_metrics["duplicate_output_qc_rows"]),
            ("qc_count_definition", "upstream_* describes upstream unmatched or skipped records; allele_mismatch_rows and duplicate_association_rows describe only retained results, so these counts need not be equal."),
            ("upstream_summary_validation", "PASS"),
            ("interpretation_boundary", "Shared expression+ASE and GWAS fine-map candidates; not formal colocalization, genome-wide significance, or proof of causality."),
        ]
        write_tsv(output_paths["summary"], SUMMARY_COLUMNS, [{"metric": key, "value": value} for key, value in summary_metrics])

        manifest_rows = []
        for role, path, rows, schema_status, row_count in [
            ("input_priority", input_paths["priority"], priority_rows, "PASS", len(priority_rows)),
            ("input_full", input_paths["full"], full_rows, "PASS", len(full_rows)),
            ("input_noncredible", input_paths["noncredible"], noncredible_rows, "PASS", len(noncredible_rows)),
            ("input_overlap_summary", input_paths["overlap_summary"], [], "PASS", overlap_summary_rows),
            ("input_trait_summary", input_paths["trait_summary"], upstream_trait_rows, "PASS", len(upstream_trait_rows)),
            ("input_branch_comparison", branch_comparison_path, [], "PASS", branch_comparison_rows),
            ("input_execution_completion", input_completion_path, [], "PASS; execution_status=success", "NA"),
            ("input_gwas_locus_manifest", gwas_run_dir / "locus_run_manifest.tsv", gwas_manifest_rows, "PASS; member paths validated", len(gwas_manifest_rows)),
        ]:
            manifest_rows.append({"role": role, "path": str(path), "exists": "yes", "size_bytes": path.stat().st_size, "rows": row_count, "schema_status": schema_status})
        output_row_counts = {
            "priority_focus": len(priority_rows),
            "priority_csv": len(priority_rows),
            "priority_html": "NA",
            "priority_svg": "NA",
            "priority_evidence": len(priority_rows),
            "gwas_credible_members": len(gwas_credible_members),
            "gwas_credible_summary": len(gwas_credible_summary_rows),
            "full": len(full_rows),
            "noncredible": len(noncredible_rows),
            "trait": len(trait_rows),
            "tissue": len(tissue_rows),
            "summary": len(summary_metrics),
        }
        for role, path in output_paths.items():
            if role == "manifest":
                continue
            manifest_rows.append({"role": f"output_{role}", "path": str(path), "exists": "yes", "size_bytes": path.stat().st_size, "rows": output_row_counts[role], "schema_status": "written"})
        write_tsv(output_paths["manifest"], MANIFEST_COLUMNS, manifest_rows)
        log(
            f"[4/4] outputs={len(output_paths) + 1}; output_dir={run_dir}",
            log_handle,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
