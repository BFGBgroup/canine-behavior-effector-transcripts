#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime
from pathlib import Path


TRAITS = (
    ("score_attachment_attention_seeking", "continuous", 468649),
    ("score_chasing", "continuous", 468649),
    ("score_dog_directed_aggression", "binary", 468647),
    ("score_dog_directed_fear", "binary", 468641),
    ("score_dog_rivalry", "binary", 468649),
    ("score_energy", "binary", 468649),
    ("score_excitability", "continuous", 468649),
    ("score_nonsocial_fear", "binary", 468641),
    ("score_owner_directed_aggression", "binary", 468649),
    ("score_separation_related_problems", "binary", 468648),
    ("score_stranger_directed_aggression", "binary", 468649),
    ("score_stranger_directed_fear", "binary", 468607),
    ("score_touch_sensitivity", "binary", 468648),
    ("score_trainability", "continuous", 468649),
)
BASE_COLUMNS = (
    "chr",
    "rs",
    "ps",
    "n_miss",
    "allele1",
    "allele0",
    "af",
    "beta",
    "se",
    "logl_H1",
    "l_remle",
    "l_mle",
    "p_wald",
    "p_lrt",
    "p_score",
)
BINARY_COLUMNS = ("OR", "OR_95CI_LOWER", "OR_95CI_UPPER")
GWAS_CHROMS = {str(index) for index in range(1, 39)}
LD_GROUPS = {"GRLS": 2846, "Oldies": 200, "PILOT": 178}
LD_SAMPLE_COUNT = 3224
LD_VARIANT_COUNT = 913984
LD_AUTOSOME_COUNT = 884031
XLSX_SHEETS = ("Quantitative_Traits", "Binary_Traits")

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
        description="Check the inputs and software environment required to rerun canine behavioral GWAS SuSiE-RSS analysis."
    )
    parser.add_argument("--gwas-dir", required=True, type=Path)
    parser.add_argument("--gwas-readme", required=True, type=Path)
    parser.add_argument("--article-xlsx", required=True, type=Path)
    parser.add_argument("--ld-prefix", required=True, type=Path)
    parser.add_argument("--ld-readme", required=True, type=Path)
    parser.add_argument("--plink-bin", required=True, type=Path)
    parser.add_argument("--rscript-bin", required=True, type=Path)
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


def build_manifest(inputs: list[tuple[str, str, Path]], issues: list[dict]) -> list[dict]:
    rows = []
    for kind, item, path in inputs:
        if not path.is_file():
            add_issue(
                issues,
                "FAIL",
                "input",
                path,
                item,
                "required_file_missing",
                "missing",
                "regular_file",
            )
            rows.append(
                {
                    "input_kind": kind,
                    "item": item,
                    "path": path,
                    "exists": "no",
                    "status": "FAIL",
                }
            )
            continue
        stat = path.stat()
        rows.append(
            {
                "input_kind": kind,
                "item": item,
                "path": path,
                "exists": "yes",
                "size_bytes": stat.st_size,
                "modified_at": datetime.fromtimestamp(stat.st_mtime).strftime(
                    "%Y-%m-%d %H:%M:%S"
                ),
                "sha256": sha256(path),
                "status": "PASS",
            }
        )
    return rows


def finite_float(value: str) -> float:
    if value in {"", "NA"}:
        raise ValueError
    result = float(value)
    if not math.isfinite(result):
        raise ValueError
    return result


def scan_gwas(
    trait: str,
    expected_type: str,
    expected_records: int,
    path: Path,
    issues: list[dict],
) -> dict:
    result = {
        "trait": trait,
        "input_file": path,
        "expected_analysis_type": expected_type,
        "expected_record_count": expected_records,
        "status": "FAIL",
    }
    if not path.is_file():
        return result

    record_count = 0
    invalid_rows = 0
    first_invalid_line = ""
    chromosomes = set()
    seen_rs = set()
    duplicate_rs = 0
    n_miss_min = None
    n_miss_max = None
    min_p_lrt = math.inf
    min_p_lrt_rs = ""
    threshold_counts = Counter()

    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = tuple(reader.fieldnames or ())
            missing = [column for column in BASE_COLUMNS if column not in header]
            binary_present = [column for column in BINARY_COLUMNS if column in header]
            if len(binary_present) == len(BINARY_COLUMNS):
                detected_type = "binary"
            elif not binary_present:
                detected_type = "continuous"
            else:
                detected_type = "partial_binary_columns"

            if missing:
                add_issue(
                    issues,
                    "FAIL",
                    "gwas",
                    path,
                    trait,
                    "missing_columns",
                    ",".join(missing),
                    ",".join(BASE_COLUMNS),
                )
            if detected_type != expected_type:
                add_issue(
                    issues,
                    "FAIL",
                    "gwas",
                    path,
                    trait,
                    "analysis_type_mismatch",
                    detected_type,
                    expected_type,
                )
            result.update(
                {
                    "detected_analysis_type": detected_type,
                    "header_columns": "|".join(header),
                }
            )
            if missing or detected_type != expected_type:
                return result

            for line_number, row in enumerate(reader, start=2):
                record_count += 1
                if None in row:
                    invalid_rows += 1
                    if not first_invalid_line:
                        first_invalid_line = str(line_number)
                    continue
                try:
                    chromosome = (row["chr"] or "").strip()
                    rs_id = (row["rs"] or "").strip()
                    position = int((row["ps"] or "").strip())
                    n_miss = int((row["n_miss"] or "").strip())
                    allele1 = (row["allele1"] or "").strip()
                    allele0 = (row["allele0"] or "").strip()
                    af = finite_float((row["af"] or "").strip())
                    finite_float((row["beta"] or "").strip())
                    se = finite_float((row["se"] or "").strip())
                    p_lrt = finite_float((row["p_lrt"] or "").strip())
                    if (
                        chromosome not in GWAS_CHROMS
                        or not rs_id
                        or position < 1
                        or n_miss < 0
                        or not allele1
                        or not allele0
                        or allele1 == allele0
                        or not 0.0 <= af <= 1.0
                        or se <= 0.0
                        or not 0.0 <= p_lrt <= 1.0
                    ):
                        raise ValueError
                except (KeyError, TypeError, ValueError):
                    invalid_rows += 1
                    if not first_invalid_line:
                        first_invalid_line = str(line_number)
                    continue

                chromosomes.add(chromosome)
                if rs_id in seen_rs:
                    duplicate_rs += 1
                else:
                    seen_rs.add(rs_id)
                n_miss_min = n_miss if n_miss_min is None else min(n_miss_min, n_miss)
                n_miss_max = n_miss if n_miss_max is None else max(n_miss_max, n_miss)
                if p_lrt < min_p_lrt:
                    min_p_lrt = p_lrt
                    min_p_lrt_rs = rs_id
                if p_lrt < 2.97e-6:
                    threshold_counts["genomewide"] += 1
                if p_lrt < 1e-5:
                    threshold_counts["suggestive"] += 1
                if p_lrt < 5e-5:
                    threshold_counts["exploratory"] += 1
    except (OSError, UnicodeError, csv.Error) as exc:
        add_issue(
            issues,
            "FAIL",
            "gwas",
            path,
            trait,
            "gzip_or_tsv_read_error",
            repr(exc),
            "readable_gzip_tsv",
        )
        return result

    checks = (
        ("record_count", record_count, expected_records),
        (
            "chromosome_set",
            ",".join(sorted(chromosomes, key=lambda value: int(value))),
            ",".join(str(i) for i in range(1, 39)),
        ),
        ("invalid_rows", invalid_rows, 0),
        ("duplicate_rs", duplicate_rs, 0),
    )
    failed = False
    for item, observed, expected in checks:
        if observed != expected:
            failed = True
            add_issue(
                issues,
                "FAIL",
                "gwas",
                path,
                trait,
                item,
                f"{observed};first_invalid_line={first_invalid_line}"
                if item == "invalid_rows"
                else observed,
                expected,
            )

    result.update(
        {
            "record_count": record_count,
            "chromosome_values": ",".join(
                sorted(chromosomes, key=lambda value: int(value) if value.isdigit() else 999)
            ),
            "duplicate_rs_count": duplicate_rs,
            "invalid_row_count": invalid_rows,
            "n_miss_min": "" if n_miss_min is None else n_miss_min,
            "n_miss_max": "" if n_miss_max is None else n_miss_max,
            "min_p_lrt": "" if not math.isfinite(min_p_lrt) else f"{min_p_lrt:.12g}",
            "min_p_lrt_rs": min_p_lrt_rs,
            "genomewide_p_lrt_lt_2_97e_6": threshold_counts["genomewide"],
            "suggestive_p_lrt_lt_1e_5": threshold_counts["suggestive"],
            "exploratory_p_lrt_lt_5e_5": threshold_counts["exploratory"],
            "status": "FAIL" if failed else "PASS",
        }
    )
    return result


def inspect_xlsx(path: Path, issues: list[dict]) -> dict:
    result = {
        "input_file": path,
        "required_sheet_names": "|".join(XLSX_SHEETS),
        "status": "FAIL",
    }
    if not path.is_file():
        return result
    try:
        with zipfile.ZipFile(path) as archive:
            corrupt = archive.testzip()
            if corrupt:
                raise zipfile.BadZipFile(corrupt)
            root = ET.fromstring(archive.read("xl/workbook.xml"))
        namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        sheets = tuple(
            node.attrib["name"] for node in root.findall("x:sheets/x:sheet", namespace)
        )
    except (KeyError, OSError, zipfile.BadZipFile, ET.ParseError) as exc:
        add_issue(
            issues,
            "FAIL",
            "article_xlsx",
            path,
            "Dataset_S1",
            "invalid_xlsx_structure",
            repr(exc),
            "valid_xlsx",
        )
        return result
    missing = [sheet for sheet in XLSX_SHEETS if sheet not in sheets]
    if missing:
        add_issue(
            issues,
            "FAIL",
            "article_xlsx",
            path,
            "Dataset_S1",
            "required_sheets_missing",
            ",".join(missing),
            "|".join(XLSX_SHEETS),
        )
    result.update(
        {
            "zip_integrity": "PASS",
            "sheet_names": "|".join(sheets),
            "status": "FAIL" if missing else "PASS",
        }
    )
    return result


def inspect_ld(prefix: Path, issues: list[dict]) -> tuple[dict, list[dict], list[dict]]:
    bed = Path(f"{prefix}.bed")
    bim = Path(f"{prefix}.bim")
    fam = Path(f"{prefix}.fam")
    if not all(path.is_file() for path in (bed, bim, fam)):
        return {"ld_prefix": prefix, "status": "FAIL"}, [], []

    group_counts = Counter()
    sample_ids = set()
    duplicate_samples = 0
    invalid_fam = 0
    sex_values = set()
    phenotype_values = set()
    with fam.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 6:
                invalid_fam += 1
                continue
            fid, iid, _father, _mother, sex, phenotype = fields
            group_counts[fid] += 1
            sex_values.add(sex)
            phenotype_values.add(phenotype)
            key = (fid, iid)
            duplicate_samples += key in sample_ids
            sample_ids.add(key)

    chromosome_counts = Counter()
    variant_ids = set()
    positions = set()
    duplicate_ids = 0
    duplicate_positions = 0
    invalid_bim = 0
    with bim.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 6:
                invalid_bim += 1
                continue
            chrom, variant_id, _cm, position, allele1, allele2 = fields
            chromosome_counts[chrom] += 1
            duplicate_ids += variant_id in variant_ids
            variant_ids.add(variant_id)
            position_key = (chrom, position)
            duplicate_positions += position_key in positions
            positions.add(position_key)
            if not allele1 or not allele2 or allele1 == allele2:
                invalid_bim += 1

    samples = sum(group_counts.values())
    variants = sum(chromosome_counts.values())
    autosomes = sum(
        count
        for chrom, count in chromosome_counts.items()
        if chrom.removeprefix("chr") in GWAS_CHROMS
    )
    with bed.open("rb") as handle:
        magic = handle.read(3).hex()
    expected_bed_size = 3 + ((samples + 3) // 4) * variants

    checks = (
        ("sample_count", samples, LD_SAMPLE_COUNT),
        ("variant_count", variants, LD_VARIANT_COUNT),
        ("autosome_variant_count", autosomes, LD_AUTOSOME_COUNT),
        ("duplicate_samples", duplicate_samples, 0),
        ("duplicate_variant_ids", duplicate_ids, 0),
        ("duplicate_positions", duplicate_positions, 0),
        ("invalid_fam_rows", invalid_fam, 0),
        ("invalid_bim_rows", invalid_bim, 0),
        ("bed_magic", magic, "6c1b01"),
        ("bed_size", bed.stat().st_size, expected_bed_size),
    )
    failed = False
    for item, observed, expected in checks:
        if observed != expected:
            failed = True
            add_issue(
                issues,
                "FAIL",
                "ld_reference",
                prefix,
                item,
                "ld_contract_mismatch",
                observed,
                expected,
            )

    group_rows = []
    for group in sorted(set(group_counts) | set(LD_GROUPS)):
        observed = group_counts.get(group, 0)
        expected = LD_GROUPS.get(group, "not_expected")
        status = "PASS" if observed == expected else "FAIL"
        group_rows.append(
            {
                "fid_group": group,
                "sample_count": observed,
                "expected_sample_count": expected,
                "status": status,
            }
        )
        if status == "FAIL":
            failed = True
            add_issue(
                issues,
                "FAIL",
                "ld_reference",
                fam,
                group,
                "sample_group_count_mismatch",
                observed,
                expected,
            )

    chromosome_rows = [
        {
            "chromosome": chrom,
            "variant_count": chromosome_counts[chrom],
            "is_autosome": "yes" if chrom.removeprefix("chr") in GWAS_CHROMS else "no",
        }
        for chrom in sorted(chromosome_counts)
    ]
    add_issue(
        issues,
        "WARN",
        "chromosome_naming",
        prefix,
        "gwas_ld_style",
        "GWAS=1..38;LD=chr1..chr38",
        "Match after explicit harmonization in Step1.4.",
    )

    qc = {
        "ld_prefix": prefix,
        "sample_count": samples,
        "expected_sample_count": LD_SAMPLE_COUNT,
        "variant_count": variants,
        "expected_variant_count": LD_VARIANT_COUNT,
        "autosome_variant_count": autosomes,
        "expected_autosome_variant_count": LD_AUTOSOME_COUNT,
        "duplicate_sample_count": duplicate_samples,
        "duplicate_variant_id_count": duplicate_ids,
        "duplicate_position_count": duplicate_positions,
        "bed_magic_hex": magic,
        "bed_size_bytes": bed.stat().st_size,
        "expected_bed_size_bytes": expected_bed_size,
        "sex_values": ",".join(sorted(sex_values)),
        "phenotype_values": ",".join(sorted(phenotype_values)),
        "status": "FAIL" if failed else "PASS",
    }
    return qc, group_rows, chromosome_rows


def version(command: tuple[str, ...]) -> tuple[bool, str]:
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return False, repr(exc)
    output = (completed.stdout + "\n" + completed.stderr).strip().replace("\n", " | ")
    return completed.returncode == 0, output


def software_rows(plink: Path, rscript: Path, issues: list[dict]) -> list[dict]:
    rows = [
        {
            "component": "python",
            "path": sys.executable,
            "version": sys.version.split()[0],
            "status": "PASS",
            "note": "The validation program uses only the Python standard library.",
        }
    ]
    checks = (
        ("plink", plink, (str(plink), "--version"), "Used for subsequent LD quality control."),
        ("Rscript", rscript, (str(rscript), "--version"), "Used for subsequent SuSiE-RSS analysis."),
        (
            "susieR",
            rscript,
            (str(rscript), "-e", 'cat(as.character(packageVersion("susieR")))'),
            "Pinning version 0.16.5 is recommended for a formal rerun.",
        ),
    )
    for component, path, command, note in checks:
        ok, observed_version = version(command)
        status = "PASS" if ok else "FAIL"
        if component == "susieR" and ok and observed_version != "0.16.5":
            status = "WARN"
        rows.append(
            {
                "component": component,
                "path": path,
                "version": observed_version,
                "status": status,
                "note": note,
            }
        )
        if status != "PASS":
            add_issue(
                issues,
                status,
                "software",
                path,
                component,
                "version_check" if ok else "version_command_failed",
                observed_version,
                "0.16.5" if component == "susieR" else "working_binary",
            )
    return rows


def run(args: argparse.Namespace) -> int:
    if args.output_dir.exists():
        raise FileExistsError(f"output directory already exists, refuse to overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    issues: list[dict] = []

    contracts = [
        {
            "trait": trait,
            "expected_analysis_type": analysis_type,
            "expected_record_count": records,
            "input_file": args.gwas_dir / f"{trait}_GWAS_summary_statistics.tsv.gz",
        }
        for trait, analysis_type, records in TRAITS
    ]
    expected_files = {Path(row["input_file"]).name for row in contracts}
    observed_files = (
        {path.name for path in args.gwas_dir.glob("*.tsv.gz")}
        if args.gwas_dir.is_dir()
        else set()
    )
    if expected_files != observed_files:
        add_issue(
            issues,
            "FAIL",
            "gwas",
            args.gwas_dir,
            "file_set",
            "|".join(sorted(observed_files)),
            "|".join(sorted(expected_files)),
        )

    input_paths = [
        ("gwas", trait, args.gwas_dir / f"{trait}_GWAS_summary_statistics.tsv.gz")
        for trait, _analysis_type, _records in TRAITS
    ]
    input_paths.extend(
        [
            ("metadata", "gwas_readme", args.gwas_readme),
            ("article", "pnas_dataset_s1", args.article_xlsx),
            ("ld", "plink_bed", Path(f"{args.ld_prefix}.bed")),
            ("ld", "plink_bim", Path(f"{args.ld_prefix}.bim")),
            ("ld", "plink_fam", Path(f"{args.ld_prefix}.fam")),
            ("metadata", "ld_readme", args.ld_readme),
        ]
    )

    log("stage 1/6: record input list and SHA256")
    manifest = build_manifest(input_paths, issues)
    write_tsv(
        args.output_dir / "input_file_manifest.tsv",
        ("input_kind", "item", "path", "exists", "size_bytes", "modified_at", "sha256", "status"),
        manifest,
    )
    write_tsv(
        args.output_dir / "trait_input_contract.tsv",
        ("trait", "expected_analysis_type", "expected_record_count", "input_file"),
        contracts,
    )

    log("stage 2/6: row by row check 14 items GWAS file")
    gwas_rows = []
    for index, (trait, analysis_type, records) in enumerate(TRAITS, start=1):
        log(f"GWAS progress {index}/{len(TRAITS)}: {trait}")
        gwas_rows.append(
            scan_gwas(
                trait,
                analysis_type,
                records,
                args.gwas_dir / f"{trait}_GWAS_summary_statistics.tsv.gz",
                issues,
            )
        )
    write_tsv(
        args.output_dir / "gwas_file_qc.tsv",
        (
            "trait",
            "input_file",
            "expected_analysis_type",
            "detected_analysis_type",
            "header_columns",
            "record_count",
            "expected_record_count",
            "chromosome_values",
            "duplicate_rs_count",
            "invalid_row_count",
            "n_miss_min",
            "n_miss_max",
            "min_p_lrt",
            "min_p_lrt_rs",
            "genomewide_p_lrt_lt_2_97e_6",
            "suggestive_p_lrt_lt_1e_5",
            "exploratory_p_lrt_lt_5e_5",
            "status",
        ),
        gwas_rows,
    )

    log("stage 3/6: check PNAS supplement table")
    article = inspect_xlsx(args.article_xlsx, issues)
    write_tsv(
        args.output_dir / "article_workbook_qc.tsv",
        ("input_file", "zip_integrity", "sheet_names", "required_sheet_names", "status"),
        [article],
    )

    log("stage 4/6: check PLINK reference panel")
    ld_qc, ld_groups, ld_chromosomes = inspect_ld(args.ld_prefix, issues)
    write_tsv(
        args.output_dir / "ld_reference_qc.tsv",
        (
            "ld_prefix",
            "sample_count",
            "expected_sample_count",
            "variant_count",
            "expected_variant_count",
            "autosome_variant_count",
            "expected_autosome_variant_count",
            "duplicate_sample_count",
            "duplicate_variant_id_count",
            "duplicate_position_count",
            "bed_magic_hex",
            "bed_size_bytes",
            "expected_bed_size_bytes",
            "sex_values",
            "phenotype_values",
            "status",
        ),
        [ld_qc],
    )
    write_tsv(
        args.output_dir / "ld_sample_group_counts.tsv",
        ("fid_group", "sample_count", "expected_sample_count", "status"),
        ld_groups,
    )
    write_tsv(
        args.output_dir / "ld_chromosome_counts.tsv",
        ("chromosome", "variant_count", "is_autosome"),
        ld_chromosomes,
    )

    log("stage 5/6: record software version")
    software = software_rows(args.plink_bin, args.rscript_bin, issues)
    write_tsv(
        args.output_dir / "software_versions.tsv",
        ("component", "path", "version", "status", "note"),
        software,
    )

    failures = sum(row["severity"] == "FAIL" for row in issues)
    warnings = sum(row["severity"] == "WARN" for row in issues)
    passed_gwas = sum(row.get("status") == "PASS" for row in gwas_rows)
    overall = "FAIL" if failures else ("PASS_WITH_WARNINGS" if warnings else "PASS")
    summary = [
        {
            "metric": "overall_status",
            "value": overall,
            "status": overall,
            "details": "FAIL blocks progression to the next step; WARN must be recorded but does not prevent completion.",
        },
        {
            "metric": "gwas_files_passed",
            "value": passed_gwas,
            "status": "PASS" if passed_gwas == len(TRAITS) else "FAIL",
            "details": f"expected={len(TRAITS)}",
        },
        {
            "metric": "ld_samples",
            "value": ld_qc.get("sample_count", ""),
            "status": ld_qc.get("status", "FAIL"),
            "details": f"expected={LD_SAMPLE_COUNT}",
        },
        {
            "metric": "ld_variants",
            "value": ld_qc.get("variant_count", ""),
            "status": ld_qc.get("status", "FAIL"),
            "details": f"expected={LD_VARIANT_COUNT}",
        },
        {
            "metric": "article_workbook",
            "value": article.get("sheet_names", ""),
            "status": article.get("status", "FAIL"),
            "details": "Dataset S1",
        },
        {
            "metric": "fail_issue_count",
            "value": failures,
            "status": "PASS" if failures == 0 else "FAIL",
            "details": "validation_issues.tsv",
        },
        {
            "metric": "warning_issue_count",
            "value": warnings,
            "status": "WARN" if warnings else "PASS",
            "details": "validation_issues.tsv",
        },
        {
            "metric": "output_directory",
            "value": args.output_dir,
            "status": "PASS",
            "details": "Separate directory for this run.",
        },
    ]

    log("stage 6/6: write the issue list and overview")
    write_tsv(args.output_dir / "validation_issues.tsv", ISSUE_COLUMNS, issues)
    write_tsv(
        args.output_dir / "validation_summary.tsv",
        ("metric", "value", "status", "details"),
        summary,
    )
    if failures:
        log(f"check failed: FAIL={failures}, WARN={warnings}")
        return 1

    (args.output_dir / "COMPLETED.ok").write_text(
        f"completed_at\t{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"overall_status\t{overall}\n",
        encoding="utf-8",
    )
    log(f"check passed: GWAS={passed_gwas}/{len(TRAITS)}, WARN={warnings}")
    log("task completed / done")
    return 0


def main() -> int:
    try:
        return run(arguments())
    except Exception as exc:
        log(f"unhandled exception: {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
