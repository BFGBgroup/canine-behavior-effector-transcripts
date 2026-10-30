#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


CANDIDATE_REQUIRED_COLUMNS = [
    "trait",
    "analysis_type",
    "chrom",
    "pos",
    "rs",
    "allele1",
    "allele0",
    "af",
    "beta",
    "se",
    "n_miss",
    "p_wald",
    "p_lrt",
    "p_score",
    "analysis_scope",
    "p_lrt_tier",
    "input_file",
]

TRAIT_REQUIRED_COLUMNS = [
    "trait",
    "analysis_type",
    "n_main",
    "locus_selection_p",
    "fine_mapping_z",
]

DERIVED_COLUMNS = [
    "candidate_row_id",
    "normalized_chrom",
    "normalized_pos1",
    "normalized_allele_pair",
    "candidate_match_key",
    "trait_metadata_present",
    "trait_metadata_analysis_type",
    "trait_metadata_n_main",
    "trait_metadata_match",
    "original_ld_variant_id",
    "original_ld_alleles",
    "final_ld_variant_id",
    "final_ld_alleles",
    "position_present_in_original_ld",
    "position_present_in_final_ld",
    "original_position_records",
    "final_position_records",
    "ld_compatibility_status",
    "removal_reason",
    "eligible_for_primary_ld",
]


VariantKey = tuple[str, int, str, str]
PositionKey = tuple[str, int]


@dataclass(frozen=True)
class BimRecord:
    chrom: str
    variant_id: str
    pos: int
    allele1: str
    allele2: str


@dataclass
class BimAudit:
    row_count: int
    exact_records: dict[VariantKey, BimRecord]
    position_records: dict[PositionKey, list[BimRecord]]
    duplicate_candidate_keys: list[str]


def timestamp() -> str:
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S ]")


def log(message: str) -> None:
    print(f"{timestamp()} {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit compatibility between GWAS candidate sites and the final PLINK LD reference panel."
    )
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--trait-metadata", type=Path, required=True)
    parser.add_argument("--step12-completion", type=Path, required=True)
    parser.add_argument("--candidate-producer", type=Path, required=True)
    parser.add_argument("--original-ld-prefix", type=Path, required=True)
    parser.add_argument("--final-ld-prefix", type=Path, required=True)
    parser.add_argument("--preliminary-maf-remove", type=Path, required=True)
    parser.add_argument("--preliminary-missingness-remove", type=Path, required=True)
    parser.add_argument("--final-qc-reasons", type=Path, required=True)
    parser.add_argument("--step15-completion", type=Path, required=True)
    parser.add_argument("--ld-producer", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def prefixed_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def require_file(path: Path, label: str, allow_empty: bool = False) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if not allow_empty and path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")


def read_tsv(
    path: Path,
    required_columns: Iterable[str],
    label: str,
) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        if columns is None:
            raise ValueError(f"{label} no column header: {path}")
        if len(columns) != len(set(columns)):
            raise ValueError(f"{label} contains duplicate column name: {columns}")
        missing = [column for column in required_columns if column not in columns]
        if missing:
            raise ValueError(f"{label} missing required columns: {','.join(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError(f"{label} no data row: {path}")
    return list(columns), rows


def normalize_chrom(value: str) -> str:
    text = value.strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if not text:
        raise ValueError("chromosome field is empty")
    if text.isdigit():
        return str(int(text))
    return text.upper()


def normalize_allele(value: str) -> str:
    allele = value.strip().upper()
    if allele in {"", ".", "NA", "NAN"}:
        raise ValueError(f"invalid allele: {value!r}")
    return allele


def normalize_position(value: str) -> int:
    try:
        position = int(value)
    except ValueError as exc:
        raise ValueError(f"invalid 1-based coordinates: {value!r}") from exc
    if position <= 0:
        raise ValueError(f"coordinates must be greater than 0: {value!r}")
    return position


def make_variant_key(
    chrom: str,
    pos: str | int,
    allele_a: str,
    allele_b: str,
) -> VariantKey:
    normalized_chrom = normalize_chrom(chrom)
    normalized_pos = normalize_position(str(pos))
    left, right = sorted((normalize_allele(allele_a), normalize_allele(allele_b)))
    return normalized_chrom, normalized_pos, left, right


def key_text(key: VariantKey) -> str:
    chrom, pos, allele_a, allele_b = key
    return f"{chrom}:{pos}:{allele_a}/{allele_b}"


def record_text(records: Iterable[BimRecord]) -> str:
    values = {
        f"{record.variant_id}:{record.allele1}/{record.allele2}"
        for record in records
    }
    return ";".join(sorted(values))


def scan_bim(
    path: Path,
    candidate_keys: set[VariantKey],
    candidate_positions: set[PositionKey],
    label: str,
) -> BimAudit:
    exact_records: dict[VariantKey, BimRecord] = {}
    position_records: dict[PositionKey, list[BimRecord]] = defaultdict(list)
    duplicate_candidate_keys: list[str] = []
    row_count = 0

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(
                    f"{label} row {line_number} does not have 6 PLINK BIM columns: {path}"
                )
            chrom, variant_id, _cm, pos_text, allele1, allele2 = fields
            key = make_variant_key(chrom, pos_text, allele1, allele2)
            position_key = key[0], key[1]
            record = BimRecord(
                chrom=chrom,
                variant_id=variant_id,
                pos=key[1],
                allele1=normalize_allele(allele1),
                allele2=normalize_allele(allele2),
            )
            row_count += 1

            if position_key in candidate_positions:
                position_records[position_key].append(record)
            if key in candidate_keys:
                if key in exact_records:
                    duplicate_candidate_keys.append(key_text(key))
                else:
                    exact_records[key] = record

    if row_count == 0:
        raise ValueError(f"{label} no valid BIM record: {path}")

    return BimAudit(
        row_count=row_count,
        exact_records=exact_records,
        position_records=dict(position_records),
        duplicate_candidate_keys=duplicate_candidate_keys,
    )


def read_variant_id_set(path: Path) -> set[str]:
    values: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip().split()[0] if line.strip() else ""
            if value and value != "variant_id":
                values.add(value)
    return values


def read_final_qc_reasons(path: Path) -> dict[str, str]:
    _columns, rows = read_tsv(path, ["variant_id", "reasons"], "final QC reason table")
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        variant_id = row["variant_id"].strip()
        reasons = row["reasons"].strip()
        if not variant_id or not reasons:
            raise ValueError(f"final QC reason table contains empty variant_id or reasons: {path}")
        grouped[variant_id].update(part for part in reasons.split(";") if part)
    return {
        variant_id: ";".join(sorted(reasons))
        for variant_id, reasons in grouped.items()
    }


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, object]]) -> None:
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


def classify_candidate(
    key: VariantKey,
    original_audit: BimAudit,
    final_audit: BimAudit,
    preliminary_maf_remove: set[str],
    preliminary_missingness_remove: set[str],
    final_qc_reasons: dict[str, str],
) -> tuple[str, str, str]:
    original_record = original_audit.exact_records.get(key)
    final_record = final_audit.exact_records.get(key)
    position_key = key[0], key[1]

    if final_record is not None:
        if original_record is None:
            return "FINAL_PRESENT_NOT_IN_ORIGINAL_LD", "unexpected_final_only", "yes"
        return "EXACT_IN_FINAL_LD", "", "yes"

    if original_record is not None:
        variant_id = original_record.variant_id
        reasons: list[str] = []
        if variant_id in preliminary_missingness_remove:
            reasons.append("PRELIMINARY_MISSINGNESS")
        if variant_id in preliminary_maf_remove:
            reasons.append("PRELIMINARY_MAF")
        if variant_id in final_qc_reasons:
            reasons.extend(
                f"FINAL_{reason}" for reason in final_qc_reasons[variant_id].split(";")
            )
        reason_text = ";".join(dict.fromkeys(reasons))

        if variant_id in final_qc_reasons:
            return "REMOVED_BY_FINAL_QC", reason_text, "no"
        if variant_id in preliminary_maf_remove:
            return "REMOVED_BY_PRELIM_MAF", reason_text, "no"
        if variant_id in preliminary_missingness_remove:
            return "REMOVED_BY_PRELIM_MISSINGNESS", reason_text, "no"
        return "REMOVED_AFTER_QC_UNCLASSIFIED", "not_listed_in_qc_removal_files", "no"

    original_position_present = position_key in original_audit.position_records
    final_position_present = position_key in final_audit.position_records
    if original_position_present or final_position_present:
        return "POSITION_PRESENT_ALLELE_MISMATCH", "coordinate_match_but_allele_set_differs", "no"
    return "NOT_IN_ORIGINAL_LD", "coordinate_absent_from_original_bim", "no"


def main() -> int:
    args = parse_args()
    log("start candidate site with LD panel compatibility audit")

    original_paths = {
        suffix: prefixed_path(args.original_ld_prefix, suffix)
        for suffix in (".bed", ".bim", ".fam")
    }
    final_paths = {
        suffix: prefixed_path(args.final_ld_prefix, suffix)
        for suffix in (".bed", ".bim", ".fam")
    }

    log("stage 1/5: check input file and upstream completion marker")
    nonempty_inputs = [
        (args.candidate_table, "Candidate variant table"),
        (args.trait_metadata, "Trait metadata table"),
        (args.step12_completion, "Step 1.2 completion marker"),
        (args.candidate_producer, "Candidate variant generation script"),
        (original_paths[".bed"], "Original LD BED"),
        (original_paths[".bim"], "Original LD BIM"),
        (original_paths[".fam"], "Original LD FAM"),
        (final_paths[".bed"], "Final LD BED"),
        (final_paths[".bim"], "Final LD BIM"),
        (final_paths[".fam"], "Final LD FAM"),
        (args.preliminary_maf_remove, "Initial MAF exclusion list"),
        (args.final_qc_reasons, "Final QC reason table"),
        (args.step15_completion, "Step 1.5 completion marker"),
        (args.ld_producer, "Final LD generation script"),
    ]
    for path, label in nonempty_inputs:
        require_file(path, label)
    require_file(
        args.preliminary_missingness_remove,
        "initial screening missingness remove list",
        allow_empty=True,
    )

    if args.output_dir.exists():
        raise FileExistsError(f"refuse to overwrite existing output directory: {args.output_dir}")

    log("stage 2/5: read candidate table and trait metadata")
    candidate_columns, candidate_rows = read_tsv(
        args.candidate_table,
        CANDIDATE_REQUIRED_COLUMNS,
        "GWAS candidate site table",
    )
    _trait_columns, trait_rows = read_tsv(
        args.trait_metadata,
        TRAIT_REQUIRED_COLUMNS,
        "trait metadata table",
    )

    trait_metadata: dict[str, dict[str, str]] = {}
    for row in trait_rows:
        trait = row["trait"].strip()
        if not trait:
            raise ValueError("trait metadata table contains empty trait")
        if trait in trait_metadata:
            raise ValueError(f"trait metadata table contains duplicate trait: {trait}")
        trait_metadata[trait] = row

    prepared_candidates: list[tuple[int, dict[str, str], VariantKey]] = []
    candidate_keys: set[VariantKey] = set()
    candidate_positions: set[PositionKey] = set()
    candidate_identity_counts: Counter[tuple[str, VariantKey]] = Counter()

    for row_number, row in enumerate(candidate_rows, start=1):
        key = make_variant_key(row["chrom"], row["pos"], row["allele1"], row["allele0"])
        prepared_candidates.append((row_number, row, key))
        candidate_keys.add(key)
        candidate_positions.add((key[0], key[1]))
        candidate_identity_counts[(row["trait"].strip(), key)] += 1

    duplicate_candidate_identities = [
        f"{trait}|{key_text(key)}"
        for (trait, key), count in candidate_identity_counts.items()
        if count > 1
    ]
    log(
        f"read candidate {len(prepared_candidates)} row, unique candidate key {len(candidate_keys)} items, "
        f"trait metadata {len(trait_metadata)} records"
    )

    log("stage 3/5: scan original and final BIM, only retain candidate related index")
    original_audit = scan_bim(
        original_paths[".bim"],
        candidate_keys,
        candidate_positions,
        "original LD BIM",
    )
    final_audit = scan_bim(
        final_paths[".bim"],
        candidate_keys,
        candidate_positions,
        "final LD BIM",
    )
    log(
        f"original BIM {original_audit.row_count} row; final BIM {final_audit.row_count} row"
    )

    preliminary_maf_remove = read_variant_id_set(args.preliminary_maf_remove)
    preliminary_missingness_remove = read_variant_id_set(
        args.preliminary_missingness_remove
    )
    final_qc_reasons = read_final_qc_reasons(args.final_qc_reasons)

    log("stage 4/5: per candidate classification and trace missing reason")
    output_rows: list[dict[str, object]] = []
    metadata_errors = 0
    unexpected_final_only = 0

    for row_number, source_row, key in prepared_candidates:
        trait = source_row["trait"].strip()
        metadata = trait_metadata.get(trait)
        metadata_present = metadata is not None
        metadata_analysis_type = metadata["analysis_type"].strip() if metadata else ""
        metadata_n_main = metadata["n_main"].strip() if metadata else ""
        metadata_match = (
            metadata_present
            and metadata_analysis_type == source_row["analysis_type"].strip()
        )
        if not metadata_match:
            metadata_errors += 1

        status, removal_reason, eligible = classify_candidate(
            key,
            original_audit,
            final_audit,
            preliminary_maf_remove,
            preliminary_missingness_remove,
            final_qc_reasons,
        )
        if status == "FINAL_PRESENT_NOT_IN_ORIGINAL_LD":
            unexpected_final_only += 1

        position_key = key[0], key[1]
        original_record = original_audit.exact_records.get(key)
        final_record = final_audit.exact_records.get(key)
        derived = {
            "candidate_row_id": row_number,
            "normalized_chrom": key[0],
            "normalized_pos1": key[1],
            "normalized_allele_pair": f"{key[2]}/{key[3]}",
            "candidate_match_key": key_text(key),
            "trait_metadata_present": "yes" if metadata_present else "no",
            "trait_metadata_analysis_type": metadata_analysis_type,
            "trait_metadata_n_main": metadata_n_main,
            "trait_metadata_match": "yes" if metadata_match else "no",
            "original_ld_variant_id": original_record.variant_id if original_record else "",
            "original_ld_alleles": (
                f"{original_record.allele1}/{original_record.allele2}"
                if original_record
                else ""
            ),
            "final_ld_variant_id": final_record.variant_id if final_record else "",
            "final_ld_alleles": (
                f"{final_record.allele1}/{final_record.allele2}"
                if final_record
                else ""
            ),
            "position_present_in_original_ld": (
                "yes" if position_key in original_audit.position_records else "no"
            ),
            "position_present_in_final_ld": (
                "yes" if position_key in final_audit.position_records else "no"
            ),
            "original_position_records": record_text(
                original_audit.position_records.get(position_key, [])
            ),
            "final_position_records": record_text(
                final_audit.position_records.get(position_key, [])
            ),
            "ld_compatibility_status": status,
            "removal_reason": removal_reason,
            "eligible_for_primary_ld": eligible,
        }
        output_rows.append({**source_row, **derived})

    status_counts = Counter(
        str(row["ld_compatibility_status"]) for row in output_rows
    )
    exact_count = status_counts["EXACT_IN_FINAL_LD"]
    missing_count = len(output_rows) - exact_count
    unclassified_count = status_counts["REMOVED_AFTER_QC_UNCLASSIFIED"]

    fatal_count = (
        len(duplicate_candidate_identities)
        + len(original_audit.duplicate_candidate_keys)
        + len(final_audit.duplicate_candidate_keys)
        + metadata_errors
        + unexpected_final_only
    )
    if fatal_count:
        overall_status = "FAIL"
    elif missing_count or unclassified_count:
        overall_status = "PASS_WITH_WARNINGS"
    else:
        overall_status = "PASS"

    log("stage 5/5: write audit result and summary")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    output_columns = candidate_columns + DERIVED_COLUMNS
    write_tsv(
        args.output_dir / "candidate_ld_compatibility.tsv",
        output_columns,
        output_rows,
    )

    missing_rows = [
        row for row in output_rows if row["ld_compatibility_status"] != "EXACT_IN_FINAL_LD"
    ]
    hwe_rows = [
        row
        for row in output_rows
        if "FINAL_HWE" in str(row["removal_reason"]).split(";")
    ]
    write_tsv(
        args.output_dir / "candidate_ld_missing.tsv",
        output_columns,
        missing_rows,
    )
    write_tsv(
        args.output_dir / "candidate_ld_hwe_removed.tsv",
        output_columns,
        hwe_rows,
    )

    by_scope_rows: list[dict[str, object]] = []
    by_scope_counts = Counter(
        (
            str(row["analysis_scope"]),
            str(row["ld_compatibility_status"]),
        )
        for row in output_rows
    )
    for (scope, status), count in sorted(by_scope_counts.items()):
        by_scope_rows.append(
            {"analysis_scope": scope, "ld_compatibility_status": status, "count": count}
        )
    write_tsv(
        args.output_dir / "candidate_ld_compatibility_by_scope.tsv",
        ["analysis_scope", "ld_compatibility_status", "count"],
        by_scope_rows,
    )

    summary_rows: list[dict[str, object]] = [
        {
            "metric": "overall_status",
            "value": overall_status,
            "status": "PASS" if overall_status != "FAIL" else "FAIL",
            "details": "Static compatibility audit of candidates against the final LD panel.",
        },
        {
            "metric": "candidate_rows",
            "value": len(output_rows),
            "status": "PASS",
            "details": "Step 1.2 candidate data rows.",
        },
        {
            "metric": "candidate_unique_variant_keys",
            "value": len(candidate_keys),
            "status": "INFO",
            "details": "Deduplicated across traits by chromosome/position/unordered allele set.",
        },
        {
            "metric": "original_bim_rows",
            "value": original_audit.row_count,
            "status": "PASS",
            "details": str(original_paths[".bim"]),
        },
        {
            "metric": "final_bim_rows",
            "value": final_audit.row_count,
            "status": "PASS",
            "details": str(final_paths[".bim"]),
        },
        {
            "metric": "exact_in_final_ld",
            "value": exact_count,
            "status": "PASS",
            "details": "Directly eligible as candidate variants in the main LD panel.",
        },
        {
            "metric": "not_exact_in_final_ld",
            "value": missing_count,
            "status": "WARN" if missing_count else "PASS",
            "details": "Retain reason labels; not directly eligible as main-LD leads.",
        },
        {
            "metric": "hwe_removed_candidate_rows",
            "value": len(hwe_rows),
            "status": "WARN" if hwe_rows else "PASS",
            "details": "removal_reason contains FINAL_HWE.",
        },
        {
            "metric": "unclassified_missing_rows",
            "value": unclassified_count,
            "status": "WARN" if unclassified_count else "PASS",
            "details": "Present in the original BIM but unexplained by the supplied exclusion lists.",
        },
        {
            "metric": "metadata_mismatch_rows",
            "value": metadata_errors,
            "status": "FAIL" if metadata_errors else "PASS",
            "details": "Candidate trait or analysis_type is inconsistent with trait_metadata.",
        },
        {
            "metric": "duplicate_candidate_trait_variant_keys",
            "value": len(duplicate_candidate_identities),
            "status": "FAIL" if duplicate_candidate_identities else "PASS",
            "details": ";".join(sorted(duplicate_candidate_identities)),
        },
        {
            "metric": "duplicate_original_bim_candidate_keys",
            "value": len(original_audit.duplicate_candidate_keys),
            "status": "FAIL" if original_audit.duplicate_candidate_keys else "PASS",
            "details": ";".join(sorted(set(original_audit.duplicate_candidate_keys))),
        },
        {
            "metric": "duplicate_final_bim_candidate_keys",
            "value": len(final_audit.duplicate_candidate_keys),
            "status": "FAIL" if final_audit.duplicate_candidate_keys else "PASS",
            "details": ";".join(sorted(set(final_audit.duplicate_candidate_keys))),
        },
        {
            "metric": "final_only_unexpected_rows",
            "value": unexpected_final_only,
            "status": "FAIL" if unexpected_final_only else "PASS",
            "details": "Exact match in the final BIM but no match in the original BIM.",
        },
    ]
    for status, count in sorted(status_counts.items()):
        summary_rows.append(
            {
                "metric": f"status_count::{status}",
                "value": count,
                "status": "INFO",
                "details": "Candidate row count.",
            }
        )
    write_tsv(
        args.output_dir / "candidate_ld_compatibility_summary.tsv",
        ["metric", "value", "status", "details"],
        summary_rows,
    )

    lineage_rows = [
        {
            "role": "candidate_table",
            "path": str(args.candidate_table),
            "producer": str(args.candidate_producer),
        },
        {
            "role": "trait_metadata",
            "path": str(args.trait_metadata),
            "producer": str(args.candidate_producer),
        },
        {
            "role": "step12_completion",
            "path": str(args.step12_completion),
            "producer": str(args.candidate_producer),
        },
        {
            "role": "original_ld_prefix",
            "path": str(args.original_ld_prefix),
            "producer": "official_GRLS_release_input",
        },
        {
            "role": "final_ld_prefix",
            "path": str(args.final_ld_prefix),
            "producer": str(args.ld_producer),
        },
        {
            "role": "preliminary_maf_remove",
            "path": str(args.preliminary_maf_remove),
            "producer": str(args.ld_producer),
        },
        {
            "role": "preliminary_missingness_remove",
            "path": str(args.preliminary_missingness_remove),
            "producer": str(args.ld_producer),
        },
        {
            "role": "final_qc_reasons",
            "path": str(args.final_qc_reasons),
            "producer": str(args.ld_producer),
        },
        {
            "role": "step15_completion",
            "path": str(args.step15_completion),
            "producer": str(args.ld_producer),
        },
    ]
    write_tsv(
        args.output_dir / "input_lineage.tsv",
        ["role", "path", "producer"],
        lineage_rows,
    )

    if overall_status == "FAIL":
        log(f"audit failed: fatal_count={fatal_count}, completion marker was not written")
        return 2

    completion_lines = [
        f"overall_status\t{overall_status}",
        f"candidate_rows\t{len(output_rows)}",
        f"exact_in_final_ld\t{exact_count}",
        f"not_exact_in_final_ld\t{missing_count}",
        f"output_dir\t{args.output_dir}",
    ]
    (args.output_dir / "COMPLETED.ok").write_text(
        "\n".join(completion_lines) + "\n",
        encoding="utf-8",
    )
    log(
        f"audit completed: status={overall_status}, exact match={exact_count}, "
        f"no exact match={missing_count}"
    )
    log("task completed / done")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        log(f"ERROR：{exc}")
        raise
