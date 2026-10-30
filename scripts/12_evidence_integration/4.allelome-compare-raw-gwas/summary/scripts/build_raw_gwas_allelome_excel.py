#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import gzip
import html
import re
import os
import signal
import traceback
from contextlib import ExitStack, redirect_stdout, redirect_stderr
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


SCRIPT_DIR = Path("/path/to/repository/scripts/12_evidence_integration/4.allelome-compare-raw-gwas/summary/scripts")
DEFAULT_OVERLAP_DETAIL = Path("/path/to/project/analysis_20260817/results/5-dogs/4.allelome-compare-raw-gwas/1.1_compare_dog5_allelome_raw_gwas/20260910_171716/dog5_allelome_raw_gwas_overlap_all_p_le_1e-4.tsv.gz")
DEFAULT_NEW_GENE_ANNOTATION = Path("/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv")
DEFAULT_OUTPUT_ROOT = Path("/path/to/project/analysis_20260817/results/5-dogs/4.allelome-compare-raw-gwas/summary")
GWAS_P_LRT_THRESHOLD = 5e-5
DEFAULT_WORKBOOK_NAME = "5dogs_raw_gwas_allelome_tables.xls"
DEFAULT_FULL_TSV_NAME = "5dogs_raw_gwas_allelome_all.tsv"
DEFAULT_FOCUS_TSV_NAME = "5dogs_raw_gwas_allelome_nature_focus.tsv"
DEFAULT_DISPLAY_CSV_NAME = "5dogs_raw_gwas_allelome_focus_table.csv"
DEFAULT_DISPLAY_HTML_NAME = "5dogs_raw_gwas_allelome_focus_table.html"
DEFAULT_DISPLAY_XLSX_NAME = "5dogs_raw_gwas_allelome_focus_table.xlsx"
DEFAULT_DISPLAY_SVG_NAME = "5dogs_raw_gwas_allelome_focus_table.svg"

def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ] {message}", flush=True)


def read_tsv_rows(path: Path) -> List[Dict[str, str]]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


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

    required = {"final_gene_id", "source_dataset", "FEELnc_biotype"}
    if not rows or required - set(rows[0]):
        raise ValueError(f"transcript_map missing fields: {path}")
    grouped = defaultdict(set)
    for row in rows:
        gid = row["final_gene_id"]
        if not gid.startswith("DOG5_NEWGENE_") or row["source_dataset"] != "five-dogs":
            raise ValueError(f"non-five dogs candidate gene: {gid}")
        grouped[gid].update(split_annotation_values(row["FEELnc_biotype"]))
    return {gid: "|".join(sorted(values)) or "NA" for gid, values in grouped.items()}

CURRENT_REQUIRED_FIELDS = {
    "trait", "gwas_marker_id", "gwas_chrom", "gwas_pos", "gwas_unordered_alleles",
    "gwas_p_lrt", "trait_variant_key", "endpoint_hit", "interval_gap_bp",
    "sample_pair_id", "canonical_pair_id", "sample_tissue", "mechanism", "linkage_score",
    "evidence_stratum", "source_link_file", "source_link_line",
    *[f"{side}_{key}" for side in ("base", "target") for key in
      ("record_type", "gene_id", "gene_symbol", "gene_biotype", "chr", "start", "end",
       "gene_strand", "manual_new_gene_class_codes")],
}

def validate_current_schema(rows, path):

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        fields = csv.DictReader(handle, delimiter="\t").fieldnames or []
    if CURRENT_REQUIRED_FIELDS - set(fields):
        raise ValueError(f"input field missing: {sorted(CURRENT_REQUIRED_FIELDS - set(fields))}")
    return fields

def build_current_enriched_rows(input_rows, feelnc_index):
    enriched = []
    for row in input_rows:
        if not 0 <= float(row["gwas_p_lrt"]) < GWAS_P_LRT_THRESHOLD:
            continue
        for side, partner in (("base", "target"), ("target", "base")):
            if not (row[f"{side}_record_type"] == "manual_new_gene"
                    and row[f"{partner}_record_type"] == "known_gene"
                    and row[f"{partner}_gene_biotype"] == "protein_coding"
                    and row["endpoint_hit"] == side):
                continue
            gid = row[f"{side}_gene_id"]
            if gid not in feelnc_index:
                raise ValueError(f"transcript_map missing novel gene: {gid}")

            interval = f'{row[f"{side}_chr"]}:{int(row[f"{side}_start"])+1}–{row[f"{side}_end"]} ({row[f"{side}_gene_strand"]})'
            enriched.append({
                **row,
                "gwas_trait": row["trait"], "gwas_p_value": row["gwas_p_lrt"],
                "gwas_variant": row["gwas_marker_id"] + " | " + row["gwas_chrom"] + ":" + row["gwas_pos"] + ":" + row["gwas_unordered_alleles"],
                "mapping_rule": "Gene span", "new_gene_interval": interval,
                "distance": row["interval_gap_bp"],
                "protein_coding_gene": row[f"{partner}_gene_symbol"] or row[f"{partner}_gene_id"],
                "protein_coding_gene_id": row[f"{partner}_gene_id"],
                "manual_new_gene": gid,
                "new_gene_class_code": row[f"{side}_manual_new_gene_class_codes"],
                "new_gene_feelnc_biotype": feelnc_index[gid],
            })
    return enriched

def row_sort_key(row):
    return (int(row["gwas_chrom"].removeprefix("chr")), int(row["gwas_pos"]),
            row["gwas_trait"], row["sample_tissue"], row["sample_pair_id"])

DISPLAY_COLUMNS = [
    "gwas_trait", "gwas_p_value", "gwas_variant", "mapping_rule",
    "new_gene_interval", "distance", "protein_coding_gene", "mechanism",
    "linkage_score", "sample_tissue", "manual_new_gene", "new_gene_class_code",
    "new_gene_feelnc_biotype", "evidence_stratum",
]
DISPLAY_LABELS = {
    "gwas_trait": "GWAS Trait", "gwas_p_value": "GWAS P LRT",
    "gwas_variant": "GWAS SNP (ID | chr:pos:alleles)", "mapping_rule": "Mapping Level",
    "new_gene_interval": "New Gene Span (1-based)", "distance": "Gene Interval Gap (bp)",
    "protein_coding_gene": "Paired Protein-Coding Gene", "mechanism": "Allelome Pattern",
    "linkage_score": "Allelome Score", "sample_tissue": "Sample Tissue",
    "manual_new_gene": "Candidate New Gene", "new_gene_class_code": "Class Code",
    "new_gene_feelnc_biotype": "FEELnc Biotype", "evidence_stratum": "Evidence Stratum",
}
DISPLAY_COLUMN_WIDTHS = {key: 26 for key in DISPLAY_COLUMNS}

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
        "Full_Results retains all upstream matches; Nature_Focus is the new-gene–protein-coding-gene display subset.",
        "focus criteria: p_lrt<5e-5, with GWAS overlapping the manual_new_gene endpoint but not the paired protein_coding known_gene endpoint.",
        "A gene-span overlap is not an exact match to the same ASE SNP and cannot establish causal regulation.",
        "Original intervals are 0-based half-open; displayed intervals are 1-based inclusive; gap is taken from upstream interval_gap_bp.",
        "FEELnc annotations are aggregated by final_gene_id in the five-dog transcript_map; labels from different transcripts may coexist.",
        "Retain evidence_stratum; multi-gene assignment risks are not automatically excluded and require separate review.",
        "source_link_file/source_link_line and sample_pair_id provide provenance; canonical_pair_id excludes animal and tissue.",
        "enhancing/repressing are Allelome patterns, not established regulatory directions.",
        "This file uses Excel 2003 XML format; a standard XLSX display workbook is also provided.",
    ]

    all_headers = list(all_rows[0]) if all_rows else sorted(CURRENT_REQUIRED_FIELDS)
    focus_headers = list(focus_rows[0]) if focus_rows else list(dict.fromkeys(all_headers + DISPLAY_COLUMNS))

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
        "5-dogs raw GWAS × Allelome association: full table",
        all_headers,
        all_rows,
    )
    add_sheet(
        "Nature_Focus",
        "5-dogs raw GWAS × Allelome association: focus table",
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


def mechanism_label(value: str) -> str:
    text = str(value).strip()
    if not text:
        return "NA"
    return text.capitalize()


def build_display_rows(rows):
    result = []
    for row in rows:
        d = {key: str(row.get(key, "")) or "NA" for key in DISPLAY_COLUMNS}
        d["mechanism"] = mechanism_label(row["mechanism"])
        result.append(d)
    return result

def build_display_summary(all_rows, focus_rows):
    return {
        "all_overlap_rows": str(len(all_rows)), "focus_rows": str(len(focus_rows)),
        "unique_trait_variants": str(len({r["trait_variant_key"] for r in focus_rows})),
        "unique_new_genes": str(len({r["manual_new_gene"] for r in focus_rows})),
        "focus_p_lt_1e-5_rows": str(sum(float(r["gwas_p_lrt"]) < 1e-5 for r in focus_rows)),
        "nonoverlap_multigene_review_rows": str(sum(r["evidence_stratum"] == "nonoverlap_multigene_review" for r in focus_rows)),
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
table{border-collapse:collapse;min-width:1600px;width:100%;font-size:12px;background:#fff}
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
            f"selected expanded rows = {html.escape(summary.get('focus_rows', 'NA'))}; "
            f"unique trait-variants = {html.escape(summary.get('unique_trait_variants', 'NA'))}; "
            f"unique new genes = {html.escape(summary.get('unique_new_genes', 'NA'))}. "
            "Rows are retained without biological-relation aggregation."
            "</p>"
        ),
        "<div class='table-wrap'><table>",
        (
            "<thead><tr class='group-head'>"
            "<th colspan='4'>GWAS evidence</th>"
            "<th colspan='6'>Allelome association evidence</th>"
            "<th colspan='3'>Manual new-gene evidence</th><th>Review flags</th>"
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
    pieces = re.split(r"([,; _|:])", value)
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
    if column == "mapping_rule":
        if value == "Yes":
            return "#e7f3ec", "#1f5a3a"
        return "#f1f3f5", "#5c6773"
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
        ("gwas_p_value", 110, 14, "numeric"),
        ("gwas_variant", 155, 18, "text_left"),
        ("mapping_rule", 120, 15, "text_left"),
        ("new_gene_interval", 175, 21, "text_left"),
        ("distance", 110, 12, "numeric"),
        ("protein_coding_gene", 155, 19, "text_left"),
        ("mechanism", 110, 15, "badge"),
        ("linkage_score", 90, 11, "numeric"),
        ("sample_tissue", 185, 22, "text_left"),
        ("manual_new_gene", 220, 25, "text_left"),
        ("new_gene_class_code", 120, 16, "text_left"),
        ("new_gene_feelnc_biotype", 155, 20, "text_left"),
        ("evidence_stratum", 210, 26, "text_left"),
    ]
    group_specs = [
        ("GWAS evidence", 0, 4, "#e8f0ea"),
        ("Allelome association evidence", 4, 10, "#eef2f7"),
        ("Manual new-gene evidence", 10, 13, "#f7efe7"),
        ("Review flags", 13, 14, "#f2eeee"),
    ]
    left_margin = 44
    top_margin = 36
    title_height = 70
    group_header_height = 24
    line_height = 11
    header_height = max(38, 10 + line_height * max(
        len(wrap_text(DISPLAY_LABELS[key], wrap_width))
        for key, _, wrap_width, _ in column_specs
    ))
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
        "p_lrt < 5e-5; GWAS within candidate new-gene span and outside paired protein-coding gene; risk flags retained."
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
            f"selected expanded rows={escape(summary.get('focus_rows', 'NA'))}; "
            f"unique trait-variants={escape(summary.get('unique_trait_variants', 'NA'))}; "
            f"unique new genes={escape(summary.get('unique_new_genes', 'NA'))}"
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
        for line_idx, line in enumerate(lines):
            svg.append(
                f'<text x="{x + col_width / 2:.1f}" y="{column_header_y + 14 + line_idx * line_height}" '
                'font-size="10" font-family="Helvetica, Arial, sans-serif" font-weight="700" fill="#111111" '
                'text-anchor="middle">'
                f"{escape(line)}</text>"
            )

    group_boundaries = [x_positions[idx] for idx in (4, 10, 13)]
    current_y = column_header_y + header_height
    previous_group_key: Tuple[str, str] | None = None
    visual_group_index = -1
    for row_idx, (row, wrapped_cells, row_height) in enumerate(zip(rows, row_line_map, row_heights), start=1):
        y = current_y
        group_key = (row.get("gwas_variant", ""), row.get("gwas_trait", ""))
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
        description="Organize current five-dog raw GWAS x Allelome x Manual new gene records into display tables."
    )
    parser.add_argument(
        "--overlap-detail",
        type=Path,
        default=DEFAULT_OVERLAP_DETAIL,
        help="Current complete step 1.1 Allelome-GWAS match table.",
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
        help="Five-dog transcript_map annotation table.",
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


    if not re.fullmatch(r"\d{8}_\d{6}", args.run_id):
        parser.error("--run-id must be YYYYMMDD_HHMMSS")
    for key, value in vars(args).items():
        if key.endswith("_name") and (Path(value).name != value or value in {"", ".", ".."}):
            parser.error("output filenames must not contain directory components")
    names = [value for key, value in vars(args).items() if key.endswith("_name")]
    if len(names) != len(set(names)):
        parser.error("output filenames must be unique")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[key] = "1"
    log_dir = SCRIPT_DIR / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    prefix = log_dir / f"{Path(__file__).stem}.{args.run_id}"
    log_path = Path(str(prefix) + ".log")
    success_path = Path(str(prefix) + ".success.tsv")
    failed_path = Path(str(prefix) + ".failed.tsv")
    run_dir = args.output_root.resolve() / args.run_id
    def stop(signum, frame):
        raise InterruptedError(f"received signal {signum}")
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    with ExitStack() as stack:
        stream = stack.enter_context(log_path.open("x", encoding="utf-8"))
        success = stack.enter_context(success_path.open("x", encoding="utf-8"))
        failed = stack.enter_context(failed_path.open("x", encoding="utf-8"))
        success.write("run_id\toutput_dir\n")
        failed.write("run_id\terror\n")
        success.flush()
        failed.flush()
        stack.enter_context(redirect_stdout(stream))
        stack.enter_context(redirect_stderr(stream))
        created = False
        try:
            log(f"input: {args.overlap_detail}")
            log(f"annotation: {args.new_gene_annotation}")
            marker = args.overlap_detail.parent / "execution_completion.marker"
            if not marker.is_file() or not re.search(r"execution_status[=\t]success", marker.read_text()):
                raise ValueError(f"upstream missing success execution marker: {marker}")
            overlap_rows = read_tsv_rows(args.overlap_detail)
            fields = validate_current_schema(overlap_rows, args.overlap_detail)
            feelnc = build_new_gene_feelnc_index(read_tsv_rows(args.new_gene_annotation), args.new_gene_annotation)
            focus_rows = sorted(build_current_enriched_rows(overlap_rows, feelnc), key=row_sort_key)
            display_rows = build_display_rows(focus_rows)
            summary = build_display_summary(overlap_rows, focus_rows)
            log(f"statistics: {summary}")
            run_dir.mkdir(parents=True, exist_ok=False)
            created = True
            focus_fields = list(focus_rows[0]) if focus_rows else list(dict.fromkeys(fields + DISPLAY_COLUMNS + ["protein_coding_gene_id"]))
            write_tsv(run_dir / args.full_tsv_name, fields, overlap_rows)
            write_tsv(run_dir / args.focus_tsv_name, focus_fields, focus_rows)
            write_excel_workbook(run_dir / args.workbook_name, all_rows=overlap_rows, focus_rows=focus_rows)
            write_csv_table(run_dir / args.display_csv_name, DISPLAY_COLUMNS, display_rows)
            title = "BarkBase raw GWAS x candidate new-gene–protein-coding Allelome links"
            write_html_table(run_dir / args.display_html_name, title=title, fieldnames=DISPLAY_COLUMNS, rows=display_rows, summary=summary)
            write_xlsx_table(run_dir / args.display_xlsx_name, title=title, fieldnames=DISPLAY_COLUMNS, rows=display_rows)
            write_svg_table(run_dir / args.display_svg_name, title=title, rows=display_rows, summary=summary)
            write_tsv(run_dir / "summary_counts.tsv", ["metric", "value"],
                      [{"metric": key, "value": value} for key, value in summary.items()])
            write_tsv(run_dir / "input_manifest.tsv", ["key", "value"], [
                {"key": "overlap_detail", "value": str(args.overlap_detail.resolve())},
                {"key": "new_gene_annotation", "value": str(args.new_gene_annotation.resolve())},
                {"key": "entry_script", "value": str(Path(__file__).resolve())},
                {"key": "p_lrt_strict_upper_bound", "value": str(GWAS_P_LRT_THRESHOLD)},
                {"key": "focus_policy", "value": "manual_new_gene hit; paired known protein_coding gene not hit; no causal direction"},
                {"key": "threads", "value": "1"},
            ])
            log("execution completed; result content and figure layout still require independent review.")
            (run_dir / "execution_completion.marker").write_text("execution_status=success\n", encoding="utf-8")
            success.write(f"{args.run_id}\t{run_dir}\n")
            success.flush()
        except BaseException as error:
            if created:
                (run_dir / "execution_completion.marker").unlink(missing_ok=True)
            failed.write(f"{args.run_id}\t{str(error).replace(chr(9), ' ').replace(chr(10), ' ')}\n")
            failed.flush()
            traceback.print_exc()
            raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
