#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


PAPER_CLAIMED_LOCUS_COUNT = 21
TABLE1_LOCUS_COUNT = 19
GENOMEWIDE_THRESHOLD = 2.97e-6
SUGGESTIVE_THRESHOLD = 1e-5
EXPLORATORY_THRESHOLD = 5e-5


TRAIT_METADATA = (
    ("score_attachment_attention_seeking", "Attachment/attention seeking", "continuous", 1187, "", "no"),
    ("score_chasing", "Chasing", "continuous", 1187, "", "no"),
    ("score_dog_directed_aggression", "Dog-directed aggression", "binary", 600, "", "yes"),
    ("score_dog_directed_fear", "Dog-directed fear", "binary", 483, "", "yes"),
    ("score_dog_rivalry", "Dog rivalry", "binary", 1030, "", "no"),
    ("score_energy", "Energy", "binary", 1187, "", "yes"),
    ("score_excitability", "Excitability", "continuous", 1187, "", "no"),
    ("score_nonsocial_fear", "Nonsocial fear", "binary", 436, "", "yes"),
    ("score_owner_directed_aggression", "Owner-directed aggression", "binary", 610, "", "no"),
    ("score_separation_related_problems", "Separation-related problems", "binary", 640, "", "yes"),
    ("score_stranger_directed_aggression", "Stranger-directed aggression", "binary", 648, "", "no"),
    ("score_stranger_directed_fear", "Stranger-directed fear", "binary", 825, "795", "yes"),
    ("score_touch_sensitivity", "Touch sensitivity", "binary", 703, "", "yes"),
    ("score_trainability", "Trainability", "continuous", 1187, "", "yes"),
)

@dataclass(frozen=True)
class Table1Locus:
    trait: str
    chrom: int
    pos: int
    allele_a: str
    allele_b: str
    paper_p_displayed: float
    roi_start_raw: int | None
    roi_end_raw: int | None
    proximal_gene: str


TABLE1_LOCI = (
    Table1Locus("score_dog_directed_aggression", 24, 36839070, "G", "A", 2.41e-6, 36737986, 36933844, "PTPN1"),
    Table1Locus("score_dog_directed_fear", 5, 22838173, "G", "A", 1.10e-5, 22835142, 22838173, "ZC3H12C"),
    Table1Locus("score_dog_directed_fear", 15, 14837906, "A", "G", 5.61e-7, 13933241, 15420702, "PRDX1"),
    Table1Locus("score_dog_directed_fear", 31, 26842062, "C", "T", 1.35e-5, None, None, "HUNK"),
    Table1Locus("score_nonsocial_fear", 10, 20464746, "A", "G", 4.26e-6, 20382051, 20506103, "FBLN1"),
    Table1Locus("score_nonsocial_fear", 12, 58861652, "T", "C", 1.54e-5, None, None, "ASCC3"),
    Table1Locus("score_nonsocial_fear", 20, 24642572, "A", "G", 6.61e-6, None, None, "KBTBD8"),
    Table1Locus("score_stranger_directed_fear", 7, 67887313, "C", "A", 2.23e-5, 66388377, 69121669, "ADCYAP1"),
    Table1Locus("score_stranger_directed_fear", 10, 69081495, "C", "A", 1.44e-6, 68695939, 69188338, "ADD2"),
    Table1Locus("score_stranger_directed_fear", 35, 8272330, "C", "T", 1.54e-6, 8272330, 8596602, "SLC35B3"),
    Table1Locus("score_energy", 15, 41054337, "G", "A", 6.55e-7, 39826935, 41688715, "PMCH"),
    Table1Locus("score_separation_related_problems", 15, 36773135, "T", "C", 2.44e-6, 35227885, 37839330, "NEDD1"),
    Table1Locus("score_separation_related_problems", 33, 22432039, "G", "A", 1.10e-5, None, None, "IGSF11"),
    Table1Locus("score_touch_sensitivity", 17, 20787827, "T", "C", 1.14e-5, 20787827, 20810875, "SLC35F6"),
    Table1Locus("score_touch_sensitivity", 22, 9094733, "T", "C", 4.14e-7, 7712048, 11232729, "VWA8"),
    Table1Locus("score_touch_sensitivity", 25, 47544666, "G", "T", 1.45e-5, None, None, "ACKR3"),
    Table1Locus("score_trainability", 25, 49735848, "C", "T", 2.11e-6, 49713759, 49741131, "ROMO1"),
    Table1Locus("score_trainability", 6, 65510825, "C", "T", 5.32e-7, None, None, "ADGRL2"),
    Table1Locus("score_trainability", 27, 20876479, "G", "A", 1.68e-6, 20876479, 20435874, "ITPR2"),
)


ISSUE_COLUMNS = (
    "severity",
    "component",
    "path",
    "item",
    "reason",
    "observed",
    "expected",
)


def timestamp() -> str:
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S ]")


def log(message: str) -> None:
    print(f"{timestamp()} {message}", flush=True)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a summary-first dog behavior GWAS candidate manifest."
    )
    parser.add_argument("--validation-dir", required=True, type=Path)
    parser.add_argument("--gwas-dir", required=True, type=Path)
    parser.add_argument("--article-xlsx", required=True, type=Path)
    parser.add_argument("--article-pdf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def add_issue(
    issues: list[dict],
    severity: str,
    component: str,
    path: Path | str,
    item: str,
    reason: str,
    observed: object = "",
    expected: object = "",
) -> None:
    issues.append(
        {
            "severity": severity,
            "component": component,
            "path": str(path),
            "item": item,
            "reason": reason,
            "observed": observed,
            "expected": expected,
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite value: {value}")
    return result


def normalize_chromosome(value: str) -> int:
    cleaned = value.strip().removeprefix("chr")
    chrom = int(cleaned)
    if not 1 <= chrom <= 38:
        raise ValueError(f"non-autosomal chromosome: {value}")
    return chrom


def allele_set(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a.strip().upper(), b.strip().upper())))


def check_inputs(args: argparse.Namespace, issues: list[dict]) -> list[dict]:
    inputs = [
        ("step11_completion", args.validation_dir / "COMPLETED.ok", "core"),
        ("step11_summary", args.validation_dir / "validation_summary.tsv", "core"),
        ("step11_gwas_qc", args.validation_dir / "gwas_file_qc.tsv", "core"),
        ("article_dataset_s1", args.article_xlsx, "auxiliary_publication_audit"),
        ("article_pdf", args.article_pdf, "auxiliary_publication_audit"),
    ]
    rows = []
    for item, path, role in inputs:
        if not path.is_file():
            severity = "FAIL" if role == "core" else "INFO"
            reason = "required_file_missing" if role == "core" else "optional_audit_file_missing"
            status = "FAIL" if role == "core" else "OPTIONAL_MISSING"
            add_issue(issues, severity, "input", path, item, reason)
            rows.append(
                {
                    "item": item,
                    "role": role,
                    "path": path,
                    "exists": "no",
                    "status": status,
                }
            )
            continue
        stat = path.stat()
        rows.append(
            {
                "item": item,
                "role": role,
                "path": path,
                "exists": "yes",
                "size_bytes": stat.st_size,
                "modified_at": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                "sha256": sha256(path),
                "status": "PASS",
            }
        )

    summary_path = args.validation_dir / "validation_summary.tsv"
    if summary_path.is_file():
        with summary_path.open("r", encoding="utf-8", newline="") as handle:
            summary = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        observed = summary.get("overall_status", "")
        if observed not in {"PASS", "PASS_WITH_WARNINGS"}:
            add_issue(
                issues,
                "FAIL",
                "lineage",
                summary_path,
                "overall_status",
                "step11_validation_not_passed",
                observed,
                "PASS_or_PASS_WITH_WARNINGS",
            )
    return rows


def xlsx_rows(path: Path) -> dict[str, list[dict[str, str]]]:
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    ns = {"m": main_ns}
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = [
                "".join(node.text or "" for node in item.findall(".//m:t", ns))
                for item in root.findall("m:si", ns)
            ]

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target_by_id = {node.attrib["Id"]: node.attrib["Target"] for node in relationships}
        result: dict[str, list[dict[str, str]]] = {}

        for sheet in workbook.findall(".//m:sheet", ns):
            name = sheet.attrib["name"]
            relation_id = sheet.attrib[f"{{{rel_ns}}}id"]
            target = "xl/" + target_by_id[relation_id].lstrip("/")
            root = ET.fromstring(archive.read(target))
            matrix = []
            for row in root.findall(".//m:sheetData/m:row", ns):
                values: dict[str, str] = {}
                for cell in row.findall("m:c", ns):
                    column = "".join(char for char in cell.attrib["r"] if char.isalpha())
                    cell_type = cell.attrib.get("t")
                    value_node = cell.find("m:v", ns)
                    value = "" if value_node is None else value_node.text or ""
                    if cell_type == "s" and value:
                        value = shared[int(value)]
                    elif cell_type == "inlineStr":
                        value = "".join(node.text or "" for node in cell.findall(".//m:t", ns))
                    values[column] = value
                matrix.append(values)

            header_index = next(
                index for index, row in enumerate(matrix) if row.get("A", "").strip() == "CHR"
            )
            header = matrix[header_index]
            column_names = {column: value.strip() for column, value in header.items()}
            records = []
            for row in matrix[header_index + 1 :]:
                if not row.get("A", "").strip():
                    continue
                records.append(
                    {
                        column_names[column]: row.get(column, "").strip()
                        for column in column_names
                    }
                )
            result[name] = records
        return result


def publication_audit_positions(
    sheets: dict[str, list[dict[str, str]]],
    article_path: Path,
    issues: list[dict],
) -> set[tuple[int, int]]:

    positions = {(locus.chrom, locus.pos) for locus in TABLE1_LOCI}
    for sheet_name, records in sheets.items():
        for row_number, row in enumerate(records, start=1):
            try:
                positions.add((normalize_chromosome(row["CHR"]), int(float(row["BP"]))))
            except (KeyError, TypeError, ValueError) as exc:
                add_issue(
                    issues,
                    "INFO",
                    "dataset_s1",
                    article_path,
                    f"{sheet_name}:record{row_number}",
                    "invalid_optional_audit_coordinate",
                    repr(exc),
                    "valid_CHR_and_BP",
                )
    return positions


def load_gwas(
    args: argparse.Namespace,
    issues: list[dict],
    audit_positions: set[tuple[int, int]],
) -> tuple[dict[str, list[dict]], dict[str, list[dict]], list[dict]]:
    candidate_rows_by_trait: dict[str, list[dict]] = {}
    audit_rows_by_trait: dict[str, list[dict]] = {}
    file_rows = []

    for trait, _label, analysis_type, _n, _sensitivity, _primary in TRAIT_METADATA:
        path = args.gwas_dir / f"{trait}_GWAS_summary_statistics.tsv.gz"
        if not path.is_file():
            add_issue(issues, "FAIL", "gwas", path, trait, "required_file_missing")
            candidate_rows_by_trait[trait] = []
            audit_rows_by_trait[trait] = []
            continue

        candidates = []
        audit_rows = []
        total_records = 0
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "chr", "rs", "ps", "n_miss", "allele1", "allele0", "af",
                "beta", "se", "p_wald", "p_lrt", "p_score",
            }
            missing = sorted(required.difference(reader.fieldnames or []))
            if missing:
                add_issue(issues, "FAIL", "gwas", path, trait, "missing_columns", ",".join(missing), ",".join(sorted(required)))
                candidate_rows_by_trait[trait] = []
                audit_rows_by_trait[trait] = []
                continue
            for line_number, row in enumerate(reader, start=2):
                total_records += 1
                try:
                    chrom = normalize_chromosome(row["chr"])
                    pos = int(row["ps"])
                    p_wald = finite_float(row["p_wald"])
                    p_lrt = finite_float(row["p_lrt"])
                    p_score = finite_float(row["p_score"])
                    beta = finite_float(row["beta"])
                    se = finite_float(row["se"])
                    n_miss = int(row["n_miss"])
                    a1 = row["allele1"].upper()
                    a0 = row["allele0"].upper()
                    af = finite_float(row["af"])
                except (TypeError, ValueError, KeyError) as exc:
                    add_issue(issues, "FAIL", "gwas", path, f"{trait}:line{line_number}", "invalid_required_value", repr(exc), "valid_numeric_and_allele_fields")
                    continue
                parsed = {
                    "trait": trait,
                    "analysis_type": analysis_type,
                    "chrom": chrom,
                    "pos": pos,
                    "rs": row["rs"],
                    "allele1": a1,
                    "allele0": a0,
                    "allele_set": allele_set(a1, a0),
                    "af": af,
                    "beta": beta,
                    "se": se,
                    "n_miss": n_miss,
                    "p_wald": p_wald,
                    "p_lrt": p_lrt,
                    "p_score": p_score,
                    "input_file": str(path),
                }
                if p_lrt < EXPLORATORY_THRESHOLD:
                    candidates.append(parsed)
                if (chrom, pos) in audit_positions:
                    audit_rows.append(parsed)
        candidate_rows_by_trait[trait] = candidates
        audit_rows_by_trait[trait] = audit_rows
        file_rows.append(
            {
                "trait": trait,
                "input_file": path,
                "total_records": total_records,
                "candidate_variant_count_p_lrt_lt_5e5": len(candidates),
                "publication_audit_position_row_count": len(audit_rows),
                "status": "PASS",
            }
        )
    return candidate_rows_by_trait, audit_rows_by_trait, file_rows


def reconcile_dataset_s1(
    sheets: dict[str, list[dict[str, str]]],
    audit_rows_by_trait: dict[str, list[dict]],
    article_path: Path,
    issues: list[dict],
) -> list[dict]:
    if not sheets:
        return []

    expected = {"Quantitative_Traits": 25, "Binary_Traits": 127}
    for sheet, count in expected.items():
        observed = len(sheets.get(sheet, []))
        if observed != count:
            add_issue(issues, "INFO", "dataset_s1", article_path, sheet, "record_count_mismatch", observed, count)

    same_trait_index = {}
    cross_trait_index = defaultdict(list)
    for trait, rows in audit_rows_by_trait.items():
        for row in rows:
            key = (row["chrom"], row["pos"], row["allele_set"])
            same_trait_index[(trait, *key)] = row
            cross_trait_index[key].append(row)

    output = []
    for sheet_name, records in sheets.items():
        for row_number, row in enumerate(records, start=1):
            try:
                trait = row["Trait"]
                chrom = normalize_chromosome(row["CHR"])
                pos = int(float(row["BP"]))
                a1 = row["ALLEL1"].upper()
                a0 = row["ALLEL0"].upper()
                published_p = finite_float(row["Pvalue"])
            except (KeyError, TypeError, ValueError) as exc:
                add_issue(issues, "FAIL", "dataset_s1", article_path, f"{sheet_name}:record{row_number}", "invalid_required_value", repr(exc), "valid_dataset_s1_record")
                continue

            key = (trait, chrom, pos, allele_set(a1, a0))
            match = same_trait_index.get(key)
            cross_matches = cross_trait_index.get((chrom, pos, allele_set(a1, a0)), [])
            exact_p_traits = sorted(
                {
                    item["trait"]
                    for item in cross_matches
                    if math.isclose(published_p, item["p_lrt"], rel_tol=1e-9, abs_tol=1e-300)
                }
            )
            if match is None:
                status = "TRAIT_OR_VARIANT_MISMATCH"
            elif math.isclose(published_p, match["p_lrt"], rel_tol=1e-9, abs_tol=1e-300):
                status = "EXACT_P_LRT_MATCH"
            else:
                status = "P_LRT_VALUE_MISMATCH"

            output.append(
                {
                    "sheet": sheet_name,
                    "sheet_record_number": row_number,
                    "trait": trait,
                    "chrom": chrom,
                    "pos": pos,
                    "published_snp": row.get("SNP", ""),
                    "published_allele1": a1,
                    "published_allele0": a0,
                    "published_pvalue": published_p,
                    "raw_gwas_p_lrt": "" if match is None else match["p_lrt"],
                    "raw_gwas_p_score": "" if match is None else match["p_score"],
                    "exact_p_lrt_match_traits_at_variant": "|".join(exact_p_traits),
                    "match_status": status,
                }
            )

    mismatch_count = sum(row["match_status"] != "EXACT_P_LRT_MATCH" for row in output)
    if mismatch_count:
        add_issue(
            issues,
            "INFO",
            "dataset_s1",
            article_path,
            "full_crosscheck",
            "published_rows_not_exactly_matched_to_same_trait_raw_p_lrt",
            mismatch_count,
            0,
        )
    return output


def p_lrt_tier(p_lrt: float) -> str:
    if p_lrt < GENOMEWIDE_THRESHOLD:
        return "genomewide_p_lrt_lt_2.97e-6"
    if p_lrt < SUGGESTIVE_THRESHOLD:
        return "primary_suggestive_p_lrt_2.97e-6_to_1e-5"
    if p_lrt < EXPLORATORY_THRESHOLD:
        return "exploratory_p_lrt_1e-5_to_5e-5"
    return "outside_candidate_threshold"


def closest_raw_p_test(displayed_p: float, match: dict | None) -> tuple[str, float | str]:

    if match is None or displayed_p <= 0:
        return "", ""
    candidates = {
        "p_wald": match["p_wald"],
        "p_lrt": match["p_lrt"],
        "p_score": match["p_score"],
    }
    valid = {name: value for name, value in candidates.items() if value > 0}
    if not valid:
        return "", ""
    name = min(valid, key=lambda item: abs(math.log10(displayed_p) - math.log10(valid[item])))
    return name, valid[name]


def reconcile_table1(
    audit_rows_by_trait: dict[str, list[dict]],
    article_pdf: Path,
    issues: list[dict],
) -> list[dict]:
    output = []
    if len(TABLE1_LOCI) != TABLE1_LOCUS_COUNT:
        add_issue(issues, "INFO", "table1", article_pdf, "transcribed_locus_count", "embedded_table_count_changed", len(TABLE1_LOCI), TABLE1_LOCUS_COUNT)

    for index, locus in enumerate(TABLE1_LOCI, start=1):
        candidates = [
            row
            for row in audit_rows_by_trait.get(locus.trait, [])
            if row["chrom"] == locus.chrom and row["pos"] == locus.pos
        ]
        exact = [row for row in candidates if row["allele_set"] == allele_set(locus.allele_a, locus.allele_b)]
        match = exact[0] if len(exact) == 1 else None
        locus_id = f"PNAS_TABLE1_{index:02d}__{locus.trait}__chr{locus.chrom}_{locus.pos}"
        roi_start = ""
        roi_end = ""
        raw_reversed = "no"
        if locus.roi_start_raw is not None and locus.roi_end_raw is not None:
            roi_start = min(locus.roi_start_raw, locus.roi_end_raw)
            roi_end = max(locus.roi_start_raw, locus.roi_end_raw)
            raw_reversed = "yes" if locus.roi_start_raw > locus.roi_end_raw else "no"

        status = "PASS" if match is not None else "FAIL"
        if match is None:
            add_issue(
                issues,
                "INFO",
                "table1",
                article_pdf,
                locus_id,
                "optional_lead_not_exactly_matched_in_raw_gwas",
                len(exact),
                1,
            )
        closest_test, closest_value = closest_raw_p_test(locus.paper_p_displayed, match)

        output.append(
            {
                "table1_locus_id": locus_id,
                "table1_row_index": index,
                "trait": locus.trait,
                "chrom": locus.chrom,
                "lead_pos": locus.pos,
                "table1_allele_a": locus.allele_a,
                "table1_allele_b": locus.allele_b,
                "table1_p_displayed": locus.paper_p_displayed,
                "raw_lead_rs": "" if match is None else match["rs"],
                "raw_allele1": "" if match is None else match["allele1"],
                "raw_allele0": "" if match is None else match["allele0"],
                "raw_beta": "" if match is None else match["beta"],
                "raw_se": "" if match is None else match["se"],
                "raw_p_wald": "" if match is None else match["p_wald"],
                "raw_p_lrt": "" if match is None else match["p_lrt"],
                "raw_p_score": "" if match is None else match["p_score"],
                "table1_p_closest_raw_test": closest_test,
                "table1_p_closest_raw_value": closest_value,
                "p_lrt_tier": "" if match is None else p_lrt_tier(match["p_lrt"]),
                "official_roi_start_raw": "" if locus.roi_start_raw is None else locus.roi_start_raw,
                "official_roi_end_raw": "" if locus.roi_end_raw is None else locus.roi_end_raw,
                "official_roi_start": roi_start,
                "official_roi_end": roi_end,
                "official_roi_was_reversed_in_table": raw_reversed,
                "proximal_candidate_gene": locus.proximal_gene,
                "match_status": status,
            }
        )
    return output


def candidate_variants(rows_by_trait: dict[str, list[dict]]) -> list[dict]:

    output = []
    for trait in sorted(rows_by_trait):
        for row in sorted(
            rows_by_trait[trait],
            key=lambda item: (item["chrom"], item["pos"], item["p_lrt"], item["rs"]),
        ):
            tier = p_lrt_tier(row["p_lrt"])
            scope = "primary" if row["p_lrt"] < SUGGESTIVE_THRESHOLD else "exploratory"
            output.append(
                {
                    **row,
                    "analysis_scope": scope,
                    "p_lrt_tier": tier,
                }
            )
    return output


def table1_reconciliation_audit(table_rows: list[dict]) -> list[dict]:

    output = []
    for row in table_rows:
        output.append(
            {
                "audit_record_type": "TABLE1_LOCUS",
                "trait": row["trait"],
                "chrom": row["chrom"],
                "start": row["official_roi_start"],
                "end": row["official_roi_end"],
                "lead_pos": row["lead_pos"],
                "lead_p_lrt": row["raw_p_lrt"],
                "table1_locus_ids": row["table1_locus_id"],
                "reconciliation_status": (
                    "TABLE1_AND_RAW_GWAS_EXACT_MATCH"
                    if row["match_status"] == "PASS"
                    else "TABLE1_NOT_EXACTLY_MATCHED_TO_RAW_GWAS"
                ),
            }
        )
    return output


def main() -> int:
    args = arguments()
    issues: list[dict] = []

    if args.output_dir.exists():
        print(f"Output directory already exists: {args.output_dir}", file=sys.stderr)
        return 2
    args.output_dir.mkdir(parents=True)

    log("Checking inputs and validated step 1.1 lineage")
    input_rows = check_inputs(args, issues)
    if any(issue["severity"] == "FAIL" for issue in issues):
        write_tsv(args.output_dir / "manifest_issues.tsv", ISSUE_COLUMNS, issues)
        return 1

    sheets: dict[str, list[dict[str, str]]] = {}
    if args.article_xlsx.is_file():
        log("Parsing optional PNAS Dataset S1 audit input")
        try:
            sheets = xlsx_rows(args.article_xlsx)
        except (OSError, KeyError, ValueError, zipfile.BadZipFile, ET.ParseError) as exc:
            add_issue(
                issues,
                "INFO",
                "dataset_s1",
                args.article_xlsx,
                "workbook",
                "optional_xlsx_parse_error",
                repr(exc),
                "valid_xlsx",
            )
    else:
        log("Optional PNAS Dataset S1 audit skipped: file not found")

    audit_positions = publication_audit_positions(sheets, args.article_xlsx, issues)
    log("Reading complete GWAS summary files and retaining only p_lrt candidates below 5e-5")
    candidate_rows_by_trait, audit_rows_by_trait, gwas_rows = load_gwas(
        args,
        issues,
        audit_positions,
    )
    candidate_rows = candidate_variants(candidate_rows_by_trait)

    log("Reconciling optional publication records against complete GWAS summary rows")
    dataset_rows = reconcile_dataset_s1(sheets, audit_rows_by_trait, args.article_xlsx, issues)
    if args.article_pdf.is_file():
        table_rows = reconcile_table1(audit_rows_by_trait, args.article_pdf, issues)
    else:
        table_rows = []
        log("Optional PNAS Table 1 audit skipped: PDF not found")
    audit_rows = table1_reconciliation_audit(table_rows)

    if table_rows:
        add_issue(
            issues,
            "INFO",
            "paper_internal_consistency",
            args.article_pdf,
            "locus_count",
            "publication_locus_count_differs_from_machine_readable_table",
            f"paper_claim={PAPER_CLAIMED_LOCUS_COUNT};table1={len(table_rows)}",
            "publication_audit_only_not_main_candidate_selection",
        )

    trait_rows = []
    for trait, label, analysis_type, n_main, n_sensitivity, primary in TRAIT_METADATA:
        trait_rows.append(
            {
                "trait": trait,
                "trait_label": label,
                "analysis_type": analysis_type,
                "n_main": n_main,
                "n_sensitivity": n_sensitivity,
                "paper_primary_trait": primary,
                "locus_selection_p": "p_lrt",
                "fine_mapping_z": "beta/se",
                "genomewide_threshold": GENOMEWIDE_THRESHOLD,
                "suggestive_threshold": SUGGESTIVE_THRESHOLD,
                "exploratory_threshold": EXPLORATORY_THRESHOLD,
                "sample_size_source": "PNAS_main_text_Figure1_and_SI_Appendix_Table_S2",
                "sample_size_note": "SI_Table_S2_main_825;Figure1_sensitivity_795" if trait == "score_stranger_directed_fear" else "",
            }
        )

    input_columns = ("item", "role", "path", "exists", "size_bytes", "modified_at", "sha256", "status")
    gwas_columns = (
        "trait", "input_file", "total_records", "candidate_variant_count_p_lrt_lt_5e5",
        "publication_audit_position_row_count", "status",
    )
    trait_columns = (
        "trait", "trait_label", "analysis_type", "n_main", "n_sensitivity", "paper_primary_trait",
        "locus_selection_p", "fine_mapping_z", "genomewide_threshold", "suggestive_threshold",
        "exploratory_threshold", "sample_size_source", "sample_size_note",
    )
    dataset_columns = (
        "sheet", "sheet_record_number", "trait", "chrom", "pos", "published_snp",
        "published_allele1", "published_allele0", "published_pvalue", "raw_gwas_p_lrt",
        "raw_gwas_p_score", "exact_p_lrt_match_traits_at_variant", "match_status",
    )
    table_columns = (
        "table1_locus_id", "table1_row_index", "trait", "chrom", "lead_pos", "table1_allele_a",
        "table1_allele_b", "table1_p_displayed", "raw_lead_rs", "raw_allele1", "raw_allele0",
        "raw_beta", "raw_se", "raw_p_wald", "raw_p_lrt", "raw_p_score",
        "table1_p_closest_raw_test", "table1_p_closest_raw_value", "p_lrt_tier",
        "official_roi_start_raw", "official_roi_end_raw", "official_roi_start", "official_roi_end",
        "official_roi_was_reversed_in_table", "proximal_candidate_gene", "match_status",
    )
    candidate_columns = (
        "trait", "analysis_type", "chrom", "pos", "rs", "allele1", "allele0", "af", "beta", "se",
        "n_miss", "p_wald", "p_lrt", "p_score", "analysis_scope", "p_lrt_tier", "input_file",
    )
    audit_columns = (
        "audit_record_type", "trait", "chrom", "start", "end", "lead_pos", "lead_p_lrt",
        "table1_locus_ids", "reconciliation_status",
    )

    write_tsv(args.output_dir / "input_lineage.tsv", input_columns, input_rows)
    write_tsv(args.output_dir / "gwas_scan_summary.tsv", gwas_columns, gwas_rows)
    write_tsv(args.output_dir / "trait_metadata.tsv", trait_columns, trait_rows)
    write_tsv(args.output_dir / "pnas_reported_variants.tsv", dataset_columns, dataset_rows)
    write_tsv(args.output_dir / "pnas_table1_loci.tsv", table_columns, table_rows)
    write_tsv(args.output_dir / "gwas_p_lrt_candidate_variants.tsv", candidate_columns, candidate_rows)
    write_tsv(args.output_dir / "locus_reconciliation_audit.tsv", audit_columns, audit_rows)
    write_tsv(args.output_dir / "manifest_issues.tsv", ISSUE_COLUMNS, issues)

    severity = Counter(issue["severity"] for issue in issues)
    if severity["FAIL"]:
        overall_status = "FAIL"
    elif severity["WARN"]:
        overall_status = "PASS_WITH_WARNINGS"
    else:
        overall_status = "PASS"

    candidate_scope_counts = Counter(row["analysis_scope"] for row in candidate_rows)
    p_lrt_tier_counts = Counter(row["p_lrt_tier"] for row in candidate_rows)
    official_roi_count = sum(row["official_roi_start"] != "" for row in table_rows)
    if not dataset_rows and not table_rows:
        publication_audit_status = "NOT_RUN"
    elif severity["INFO"]:
        publication_audit_status = "COMPLETED_WITH_INFORMATIONAL_DIFFERENCES"
    else:
        publication_audit_status = "COMPLETED"
    summary_rows = [
        {"metric": "overall_status", "value": overall_status},
        {"metric": "main_analysis_source", "value": "complete_gwas_summary_statistics"},
        {"metric": "gwas_summary_p_lrt_candidate_variant_count", "value": len(candidate_rows)},
        {"metric": "primary_p_lrt_lt_1e5_candidate_variant_count", "value": candidate_scope_counts["primary"]},
        {"metric": "genomewide_p_lrt_lt_2_97e6_candidate_variant_count", "value": p_lrt_tier_counts["genomewide_p_lrt_lt_2.97e-6"]},
        {"metric": "primary_suggestive_p_lrt_2_97e6_to_1e5_candidate_variant_count", "value": p_lrt_tier_counts["primary_suggestive_p_lrt_2.97e-6_to_1e-5"]},
        {"metric": "exploratory_p_lrt_1e5_to_5e5_candidate_variant_count", "value": candidate_scope_counts["exploratory"]},
        {"metric": "warning_count", "value": severity["WARN"]},
        {"metric": "failure_count", "value": severity["FAIL"]},
        {"metric": "informational_audit_count", "value": severity["INFO"]},
        {"metric": "selection_p", "value": "p_lrt"},
        {"metric": "fine_mapping_z", "value": "beta/se"},
        {"metric": "region_window_bp", "value": "not_applied_in_step_1.2"},
        {"metric": "independent_locus_definition_status", "value": "not_built_in_step_1.2"},
        {"metric": "region_definition_status", "value": "deferred_to_ld_region_step"},
        {"metric": "publication_audit_status", "value": publication_audit_status},
        {"metric": "publication_paper_claimed_locus_count", "value": PAPER_CLAIMED_LOCUS_COUNT if table_rows else ""},
        {"metric": "publication_table1_locus_count", "value": len(table_rows)},
        {"metric": "publication_table1_exact_raw_gwas_matches", "value": sum(row["match_status"] == "PASS" for row in table_rows)},
        {"metric": "publication_dataset_s1_record_count", "value": len(dataset_rows)},
        {"metric": "publication_dataset_s1_exact_p_lrt_matches", "value": sum(row["match_status"] == "EXACT_P_LRT_MATCH" for row in dataset_rows)},
        {"metric": "publication_table1_official_roi_count", "value": official_roi_count},
        {"metric": "publication_table1_missing_official_roi_count", "value": len(table_rows) - official_roi_count},
    ]
    write_tsv(args.output_dir / "manifest_summary.tsv", ("metric", "value"), summary_rows)

    publication_audit_text = (
        "PUBLICATION AUDIT ONLY\n"
        "Primary candidates in this step are determined only by complete GWAS summary statistics and the p_lrt threshold.\n"
        "PNAS Dataset S1, Table 1, the published locus count, and ROI are used only for auxiliary checks, not candidate retention or main-workflow status decisions.\n"
        "This step does not construct fixed windows or independent LD loci; regions should subsequently be defined using the quality-controlled LD reference panel.\n"
    )
    (args.output_dir / "PUBLICATION_AUDIT.txt").write_text(publication_audit_text, encoding="utf-8")

    if overall_status == "FAIL":
        log(f"Manifest build failed with {severity['FAIL']} failure(s)")
        return 1

    (args.output_dir / "COMPLETED.ok").write_text(
        f"status={overall_status}\ncompleted_at={datetime.now().isoformat(timespec='seconds')}\n",
        encoding="utf-8",
    )
    log(
        "Summary-first candidate manifest completed: "
        f"p_lrt_candidate_variants={len(candidate_rows)}, "
        f"primary={candidate_scope_counts['primary']}, "
        f"exploratory={candidate_scope_counts['exploratory']}, "
        f"warnings={severity['WARN']}, publication_info={severity['INFO']}"
    )
    log(f"Output directory: {args.output_dir}")
    log("task completed / done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
