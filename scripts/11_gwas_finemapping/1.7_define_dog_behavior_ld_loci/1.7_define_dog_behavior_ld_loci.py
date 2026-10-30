#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import gzip
import math
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


EXACT_STATUS = "EXACT_IN_FINAL_LD"
PRIMARY_SCOPE = "primary"
EXPLORATORY_SCOPE = "exploratory"


@dataclass(frozen=True)
class BimVariant:
    chrom: str
    variant_id: str
    pos: int
    allele1: str
    allele2: str


def timestamp() -> str:
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S ]")


def log(message: str) -> None:
    print(f"{timestamp()} {message}", flush=True)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Define canine behavioral GWAS loci using the final primary LD panel and the r2>=0.7 rule."
    )
    parser.add_argument("--compatibility-table", required=True, type=Path)
    parser.add_argument("--trait-metadata", required=True, type=Path)
    parser.add_argument("--pnas-table1-loci", required=True, type=Path)
    parser.add_argument("--step16-completion", required=True, type=Path)
    parser.add_argument("--candidate-producer", required=True, type=Path)
    parser.add_argument("--audit-producer", required=True, type=Path)
    parser.add_argument("--ld-prefix", required=True, type=Path)
    parser.add_argument("--step15-completion", required=True, type=Path)
    parser.add_argument("--ld-producer", required=True, type=Path)
    parser.add_argument("--plink2-bin", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--ld-r2-threshold", type=float, default=0.7)
    parser.add_argument("--ld-window-kb", type=int, default=9_999_999)
    parser.add_argument("--ld-window-variant-count", type=int, default=99_999_999)
    parser.add_argument("--fallback-window-bp", type=int, default=1_000_000)
    return parser.parse_args()


def normalize_chrom(value: str) -> str:
    chrom = value.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    require(chrom.isdigit(), f"unexpected autosome encoding: {value}")
    number = int(chrom)
    require(1 <= number <= 38, f"chromosome not in dog autosome 1-38 span: {value}")
    return str(number)


def chrom_sort_key(chrom: str) -> int:
    return int(normalize_chrom(chrom))


def allele_pair(allele_a: str, allele_b: str) -> str:
    alleles = sorted([allele_a.strip().upper(), allele_b.strip().upper()])
    require(all(alleles), "allele must not be empty")
    return "/".join(alleles)


def parse_int(value: str, field: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise RuntimeError(f"field {field} is not an integer: {value}") from exc
    return parsed


def parse_float(value: str, field: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise RuntimeError(f"field {field} is not numeric: {value}") from exc
    require(math.isfinite(parsed), f"field {field} is not finite numeric value: {value}")
    return parsed


def read_tsv(path: Path, required_fields: Iterable[str]) -> tuple[list[dict], list[str]]:
    require(path.is_file() and path.stat().st_size > 0, f"input does not exist or is empty: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"TSV missing header: {path}")
        missing = sorted(set(required_fields) - set(reader.fieldnames))
        require(not missing, f"TSV missing required fields {missing}: {path}")
        rows = list(reader)
    require(rows, f"TSV no data row: {path}")
    return rows, list(reader.fieldnames)


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_tsv_gz(path: Path, fieldnames: list[str], rows: Iterable[dict]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_completion_values(path: Path) -> dict[str, str]:
    require(path.is_file() and path.stat().st_size > 0, f"completion marker does not exist or is empty: {path}")
    values = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t", 1)
            if len(fields) == 2:
                values[fields[0]] = fields[1]
    require(values, f"completion marker no can parse of key value: {path}")
    return values


def validate_step16_completion(path: Path) -> dict[str, str]:
    values = read_completion_values(path)
    status = values.get("overall_status", "")
    require(status in {"PASS", "PASS_WITH_WARNINGS"}, f"step 1.6 failed: {path} status={status}")
    require(parse_int(values.get("candidate_rows", ""), "candidate_rows") > 0, "step 1.6 candidate count invalid")
    return values


def validate_step15_completion(path: Path, expected_prefix: Path) -> dict[str, str]:
    values = read_completion_values(path)
    required_keys = {"completed_at", "primary_prefix", "primary_samples", "primary_variants"}
    missing = sorted(required_keys - set(values))
    require(not missing, f"step 1.5 completion marker missing fields {missing}: {path}")
    require(Path(values["primary_prefix"]) == expected_prefix, "step 1.5 completion marker in main LD prefix mismatch")
    require(parse_int(values["primary_samples"], "primary_samples") > 0, "step 1.5 sample count invalid")
    require(parse_int(values["primary_variants"], "primary_variants") > 0, "step 1.5 variant count invalid")
    return values


def read_bim(prefix: Path) -> tuple[dict[str, BimVariant], dict[str, tuple[int, int]], int]:
    path = Path(f"{prefix}.bim")
    require(path.is_file() and path.stat().st_size > 0, f"BIM does not exist or is empty: {path}")
    variants: dict[str, BimVariant] = {}
    chrom_bounds: dict[str, list[int]] = {}
    row_count = 0
    with path.open("r", encoding="utf-8") as handle:
        for row_count, raw_line in enumerate(handle, start=1):
            fields = raw_line.split()
            require(len(fields) >= 6, f"BIM row {row_count} has fewer than 6 columns: {path}")
            chrom = normalize_chrom(fields[0])
            variant_id = fields[1]
            pos = parse_int(fields[3], "BIM_POS")
            require(pos >= 1, f"BIM coordinates must be 1-based positive integer: {variant_id} {pos}")
            require(variant_id not in variants, f"final BIM contains duplicate variant ID: {variant_id}")
            variants[variant_id] = BimVariant(
                chrom=chrom,
                variant_id=variant_id,
                pos=pos,
                allele1=fields[4].upper(),
                allele2=fields[5].upper(),
            )
            if chrom not in chrom_bounds:
                chrom_bounds[chrom] = [pos, pos]
            else:
                chrom_bounds[chrom][0] = min(chrom_bounds[chrom][0], pos)
                chrom_bounds[chrom][1] = max(chrom_bounds[chrom][1], pos)
    require(row_count > 0, f"BIM no data row: {path}")
    return variants, {key: tuple(value) for key, value in chrom_bounds.items()}, row_count


def audit_fam(prefix: Path) -> tuple[int, int]:
    path = Path(f"{prefix}.fam")
    require(path.is_file() and path.stat().st_size > 0, f"FAM does not exist or is empty: {path}")
    sample_count = 0
    nonfounder_count = 0
    with path.open("r", encoding="utf-8") as handle:
        for sample_count, raw_line in enumerate(handle, start=1):
            fields = raw_line.split()
            require(len(fields) >= 6, f"FAM row {sample_count} has fewer than 6 columns: {path}")
            if fields[2] != "0" or fields[3] != "0":
                nonfounder_count += 1
    require(sample_count > 0, f"FAM no sample: {path}")
    require(
        nonfounder_count == 0,
        f"PLINK 2 LD uses only founders by default, but FAM contains {nonfounder_count} non-founder samples",
    )
    return sample_count, nonfounder_count


def validate_bed(prefix: Path, sample_count: int, variant_count: int) -> tuple[int, int, str]:
    path = Path(f"{prefix}.bed")
    require(path.is_file() and path.stat().st_size > 0, f"BED does not exist or is empty: {path}")
    with path.open("rb") as handle:
        magic = handle.read(3)
    require(magic == b"\x6c\x1b\x01", f"BED is not PLINK SNP-major binary format: {path}")
    bytes_per_variant = (sample_count + 3) // 4
    expected_size = 3 + bytes_per_variant * variant_count
    actual_size = path.stat().st_size
    require(
        actual_size == expected_size,
        f"BED size differs from the value expected from FAM/BIM: actual={actual_size} expected={expected_size}",
    )
    return actual_size, expected_size, magic.hex()


def plink2_version(plink2_bin: Path) -> str:
    completed = subprocess.run(
        [str(plink2_bin), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    require(completed.returncode == 0, f"cannot read PLINK2 version: {plink2_bin}")
    version = (completed.stdout or completed.stderr).strip()
    require(version.startswith("PLINK v2."), f"specified program is not PLINK2: {version}")
    return version


def run_command(command: list[object], label: str) -> None:
    command_text = " ".join(shlex.quote(str(item)) for item in command)
    log(f"RUN {label}: {command_text}")
    completed = subprocess.run([str(item) for item in command], check=False)
    require(completed.returncode == 0, f"{label} failed, exit code={completed.returncode}")


def interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def interval_union_length(intervals: Iterable[tuple[int, int]]) -> int:
    sorted_intervals = sorted(intervals)
    if not sorted_intervals:
        return 0
    for start, end in sorted_intervals:
        require(start <= end, f"interval start greater than end: {start}-{end}")
    current_start, current_end = sorted_intervals[0]
    total = 0
    for start, end in sorted_intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
            continue
        total += current_end - current_start + 1
        current_start, current_end = start, end
    return total + current_end - current_start + 1


def core_intervals_overlap(region_a: dict, region_b: dict) -> bool:
    return (
        region_a["chrom"] == region_b["chrom"]
        and interval_overlap(
            region_a["core_start"],
            region_a["core_end"],
            region_b["core_start"],
            region_b["core_end"],
        )
        > 0
    )


def pnas_roi_bounds(row: dict) -> tuple[int, int, bool]:
    start_text = row["official_roi_start"].strip()
    end_text = row["official_roi_end"].strip()
    require(bool(start_text) == bool(end_text), f"PNAS ROI only single side coordinates: {row['table1_locus_id']}")
    lead_pos = parse_int(row["lead_pos"], "lead_pos")
    if not start_text:
        return lead_pos, lead_pos, False
    start = parse_int(start_text, "official_roi_start")
    end = parse_int(end_text, "official_roi_end")
    require(start >= 1, f"PNAS standard ROI start must be positive integer: {row['table1_locus_id']}")
    require(start <= end, f"PNAS standard ROI start greater than end: {row['table1_locus_id']}")
    return start, end, True


def parse_vcor(
    path: Path,
    lead_ids: set[str],
    bim_variants: dict[str, BimVariant],
    r2_threshold: float,
) -> tuple[dict[str, dict[str, float]], int]:
    require(path.is_file() and path.stat().st_size > 0, f"PLINK2 not generated valid VCOR: {path}")
    ld_members: dict[str, dict[str, float]] = {
        lead_id: {lead_id: 1.0} for lead_id in sorted(lead_ids)
    }
    pair_count = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        header_line = handle.readline().rstrip("\n")
        require(header_line, f"VCOR missing header: {path}")
        delimiter = "\t" if "\t" in header_line else None
        header = header_line.split(delimiter)
        header[0] = header[0].lstrip("#")
        required = {
            "CHROM_A",
            "POS_A",
            "ID_A",
            "CHROM_B",
            "POS_B",
            "ID_B",
            "UNPHASED_R2",
        }
        require(required.issubset(set(header)), f"VCOR header and PLINK2 official format does not match: {header}")
        for line_number, raw_line in enumerate(handle, start=2):
            values = raw_line.rstrip("\n").split(delimiter)
            require(len(values) == len(header), f"VCOR row {line_number} has an unexpected column count")
            row = dict(zip(header, values))
            lead_id = row["ID_A"]
            member_id = row["ID_B"]
            require(lead_id in lead_ids, f"VCOR No. one variant is not specified lead: {lead_id}")
            require(lead_id in bim_variants and member_id in bim_variants, "VCOR contains BIM external variants")
            lead_variant = bim_variants[lead_id]
            member_variant = bim_variants[member_id]
            require(lead_variant.chrom == member_variant.chrom, "VCOR encountered across chromosome pair")
            require(normalize_chrom(row["CHROM_A"]) == lead_variant.chrom, "VCOR lead chromosome mismatch")
            require(normalize_chrom(row["CHROM_B"]) == member_variant.chrom, "VCOR member chromosome mismatch")
            require(parse_int(row["POS_A"], "POS_A") == lead_variant.pos, "VCOR lead coordinates mismatch")
            require(parse_int(row["POS_B"], "POS_B") == member_variant.pos, "VCOR member coordinates mismatch")
            r2_value = parse_float(row["UNPHASED_R2"], "UNPHASED_R2")
            require(r2_value + 1e-12 >= r2_threshold, f"VCOR contains below threshold of r2: {r2_value}")
            previous = ld_members[lead_id].get(member_id)
            if previous is None or r2_value > previous:
                ld_members[lead_id][member_id] = r2_value
            pair_count += 1
    return ld_members, pair_count


def make_seed_region(
    row: dict,
    ld_members: dict[str, dict[str, float]],
    bim_variants: dict[str, BimVariant],
    chrom_bounds: dict[str, tuple[int, int]],
    fallback_window_bp: int,
) -> dict:
    lead_id = row["final_ld_variant_id"]
    members = ld_members[lead_id]
    member_positions = [bim_variants[member_id].pos for member_id in members]
    core_start = min(member_positions)
    core_end = max(member_positions)
    proxy_count = len(members) - 1
    chrom = row["_chrom"]
    chrom_min, chrom_max = chrom_bounds[chrom]
    if proxy_count > 0:
        analysis_start = core_start
        analysis_end = core_end
        region_source = "R2_GE_THRESHOLD_BOUNDARIES"
    else:
        analysis_start = max(chrom_min, row["_pos"] - fallback_window_bp)
        analysis_end = min(chrom_max, row["_pos"] + fallback_window_bp)
        region_source = "FALLBACK_LEAD_PLUS_MINUS_WINDOW_NO_R2_PROXY"
    return {
        "chrom": chrom,
        "core_start": core_start,
        "core_end": core_end,
        "analysis_start": analysis_start,
        "analysis_end": analysis_end,
        "region_source": region_source,
        "component_region_sources": {region_source},
        "component_core_intervals": [(core_start, core_end, region_source)],
        "merge_event_count": 0,
        "ld_core_overlap_merge_count": 0,
        "no_proxy_core_overlap_merge_count": 0,
        "component_rows": [row],
        "member_ids": set(members),
        "fallback_component_count": 1 if proxy_count == 0 else 0,
    }


def merge_regions(seed_regions: list[dict]) -> list[dict]:
    if not seed_regions:
        return []
    ordered = sorted(
        seed_regions,
        key=lambda row: (
            row["component_rows"][0]["trait"],
            chrom_sort_key(row["chrom"]),
            row["core_start"],
            row["core_end"],
            row["analysis_start"],
            row["analysis_end"],
            min(item["_p_lrt"] for item in row["component_rows"]),
        ),
    )
    merged: list[dict] = []
    for region in ordered:
        region_trait = region["component_rows"][0]["trait"]
        current_trait = "" if not merged else merged[-1]["component_rows"][0]["trait"]
        if (
            not merged
            or current_trait != region_trait
            or not core_intervals_overlap(merged[-1], region)
        ):
            merged.append(
                {
                    **region,
                    "component_rows": list(region["component_rows"]),
                    "member_ids": set(region["member_ids"]),
                    "component_region_sources": set(region["component_region_sources"]),
                    "component_core_intervals": list(region["component_core_intervals"]),
                }
            )
            continue
        current = merged[-1]
        core_overlap = any(
            left_source == "R2_GE_THRESHOLD_BOUNDARIES"
            and right_source == "R2_GE_THRESHOLD_BOUNDARIES"
            and interval_overlap(left_start, left_end, right_start, right_end) > 0
            for left_start, left_end, left_source in current["component_core_intervals"]
            for right_start, right_end, right_source in region["component_core_intervals"]
        )
        current["core_start"] = min(current["core_start"], region["core_start"])
        current["core_end"] = max(current["core_end"], region["core_end"])
        current["analysis_start"] = min(current["analysis_start"], region["analysis_start"])
        current["analysis_end"] = max(current["analysis_end"], region["analysis_end"])
        current["component_rows"].extend(region["component_rows"])
        current["member_ids"].update(region["member_ids"])
        current["component_region_sources"].update(region["component_region_sources"])
        current["component_core_intervals"].extend(region["component_core_intervals"])
        current["merge_event_count"] += 1 + region["merge_event_count"]
        current["ld_core_overlap_merge_count"] += region["ld_core_overlap_merge_count"]
        current["no_proxy_core_overlap_merge_count"] += region[
            "no_proxy_core_overlap_merge_count"
        ]
        if core_overlap:
            current["ld_core_overlap_merge_count"] += 1
        else:
            current["no_proxy_core_overlap_merge_count"] += 1
        current["fallback_component_count"] += region["fallback_component_count"]
        current["region_source"] = "MERGED_COMPONENT_REGIONS"
    for region in merged:
        region["component_rows"].sort(
            key=lambda row: (row["_p_lrt"], row["_pos"], row["candidate_row_id"])
        )
        region["representative_lead"] = region["component_rows"][0]
        if region["merge_event_count"] == 0:
            region["merge_basis"] = "NONE"
        elif region["no_proxy_core_overlap_merge_count"] and region[
            "ld_core_overlap_merge_count"
        ]:
            region["merge_basis"] = "LD_CORE_AND_NO_PROXY_LEAD_POSITION_OVERLAP"
        elif region["no_proxy_core_overlap_merge_count"]:
            region["merge_basis"] = "NO_PROXY_LEAD_POSITION_OVERLAP"
        else:
            region["merge_basis"] = "LD_CORE_OVERLAP"
    return merged


def assign_locus_ids(regions: list[dict], scope: str) -> None:
    regions.sort(
        key=lambda row: (
            row["representative_lead"]["trait"],
            chrom_sort_key(row["chrom"]),
            row["analysis_start"],
            row["analysis_end"],
        )
    )
    counters: Counter[tuple[str, str]] = Counter()
    scope_code = "PRIMARY" if scope == PRIMARY_SCOPE else "EXPLORATORY"
    for region in regions:
        trait = region["representative_lead"]["trait"]
        counters[(trait, region["chrom"])] += 1
        region["analysis_scope"] = scope
        region["locus_id"] = (
            f"DOG_GWAS_{scope_code}__{trait}__chr{region['chrom']}_"
            f"{counters[(trait, region['chrom'])]:03d}"
        )


def choose_primary_overlap(region: dict, primary_regions: list[dict]) -> dict | None:
    candidates = []
    for primary in primary_regions:
        if primary["representative_lead"]["trait"] != region["representative_lead"]["trait"]:
            continue
        overlap = interval_overlap(
            region["analysis_start"],
            region["analysis_end"],
            primary["analysis_start"],
            primary["analysis_end"],
        )
        if overlap > 0 and region["chrom"] == primary["chrom"]:
            candidates.append((overlap, -primary["representative_lead"]["_p_lrt"], primary))
    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return candidates[0][2]


def find_containing_loci(row: dict, regions: list[dict]) -> list[dict]:
    matches = [
        region
        for region in regions
        if region["representative_lead"]["trait"] == row["trait"]
        and region["chrom"] == row["_chrom"]
        and region["analysis_start"] <= row["_pos"] <= region["analysis_end"]
    ]
    matches.sort(
        key=lambda region: (
            0 if region["analysis_scope"] == PRIMARY_SCOPE else 1,
            region["representative_lead"]["_p_lrt"],
            region["locus_id"],
        )
    )
    return matches


def format_number(value: float) -> str:
    return f"{value:.12g}"


def main() -> None:
    args = parse_args()
    require(args.threads >= 1, "--threads must be>=1")
    require(0.0 < args.ld_r2_threshold <= 1.0, "--ld-r2-threshold must be located in(0,1]")
    require(args.ld_window_kb >= 1, "--ld-window-kb must be>=1")
    require(args.ld_window_variant_count >= 2, "--ld-window-variant-count must be>=2")
    require(args.fallback_window_bp >= 0, "--fallback-window-bp must be>=0")
    require(not args.output_dir.exists(), f"refuse to overwrite existing output directory: {args.output_dir}")

    log("stage 1/7: check upstream completion marker, actual input and PLINK2 version")
    step16_completion = validate_step16_completion(args.step16_completion)
    step15_completion = validate_step15_completion(args.step15_completion, args.ld_prefix)
    for required_path in (
        args.candidate_producer,
        args.audit_producer,
        args.ld_producer,
        args.plink2_bin,
        Path(f"{args.ld_prefix}.bed"),
        Path(f"{args.ld_prefix}.bim"),
        Path(f"{args.ld_prefix}.fam"),
    ):
        require(required_path.is_file() and required_path.stat().st_size > 0, f"required file does not exist or is empty: {required_path}")
    version = plink2_version(args.plink2_bin)

    compatibility_fields = {
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
        "candidate_row_id",
        "normalized_chrom",
        "normalized_pos1",
        "normalized_allele_pair",
        "candidate_match_key",
        "final_ld_variant_id",
        "final_ld_alleles",
        "ld_compatibility_status",
        "removal_reason",
        "eligible_for_primary_ld",
    }
    trait_fields = {
        "trait",
        "trait_label",
        "analysis_type",
        "n_main",
        "genomewide_threshold",
        "suggestive_threshold",
        "exploratory_threshold",
    }
    pnas_fields = {
        "table1_locus_id",
        "trait",
        "chrom",
        "lead_pos",
        "table1_allele_a",
        "table1_allele_b",
        "table1_p_displayed",
        "raw_p_lrt",
        "official_roi_start",
        "official_roi_end",
        "proximal_candidate_gene",
        "match_status",
    }
    compatibility_rows, compatibility_columns = read_tsv(
        args.compatibility_table, compatibility_fields
    )
    trait_rows, _ = read_tsv(args.trait_metadata, trait_fields)
    pnas_rows, _ = read_tsv(args.pnas_table1_loci, pnas_fields)
    trait_by_id = {row["trait"]: row for row in trait_rows}
    require(
        len(compatibility_rows)
        == parse_int(step16_completion["candidate_rows"], "candidate_rows"),
        "Candidate row count in the Step 1.6 completion marker differs from the compatibility table",
    )
    require(len(trait_by_id) == len(trait_rows), "trait_metadata contains duplicate trait")
    require(len(trait_rows) == 14, f"expected 14 traits, actual {len(trait_rows)}")
    require(len(pnas_rows) == 19, f"expected 19 records PNAS Table 1 locus, actual {len(pnas_rows)}")
    pnas_by_candidate_key = {}
    for paper_row in pnas_rows:
        require(paper_row["trait"] in trait_by_id, f"PNAS locus of trait not in metadata in: {paper_row['trait']}")
        paper_chrom = normalize_chrom(paper_row["chrom"])
        paper_pos = parse_int(paper_row["lead_pos"], "lead_pos")
        pnas_roi_bounds(paper_row)
        paper_key = (
            paper_row["trait"],
            paper_chrom,
            paper_pos,
            allele_pair(paper_row["table1_allele_a"], paper_row["table1_allele_b"]),
        )
        require(paper_key not in pnas_by_candidate_key, f"PNAS Table 1 candidate duplicate key: {paper_key}")
        pnas_by_candidate_key[paper_key] = paper_row

    log("stage 2/7: read final main LD BIM/FAM and validate candidate key")
    bim_variants, chrom_bounds, bim_row_count = read_bim(args.ld_prefix)
    sample_count, nonfounder_count = audit_fam(args.ld_prefix)
    require(
        bim_row_count == parse_int(step15_completion["primary_variants"], "primary_variants"),
        "Variant count in the Step 1.5 completion marker differs from the BIM",
    )
    require(
        sample_count == parse_int(step15_completion["primary_samples"], "primary_samples"),
        "Sample count in the Step 1.5 completion marker differs from the FAM",
    )
    bed_actual_size, bed_expected_size, bed_magic_hex = validate_bed(
        args.ld_prefix, sample_count, bim_row_count
    )
    seen_trait_candidate_keys = set()
    seen_candidate_row_ids = set()
    exact_rows = []
    for row in compatibility_rows:
        require(row["trait"] in trait_by_id, f"candidate trait not in metadata in: {row['trait']}")
        require(row["analysis_scope"] in {PRIMARY_SCOPE, EXPLORATORY_SCOPE}, f"analysis_scope invalid: {row['analysis_scope']}")
        require(row["analysis_type"] == trait_by_id[row["trait"]]["analysis_type"], f"analysis_type mismatch: {row['trait']}")
        row["_chrom"] = normalize_chrom(row["normalized_chrom"])
        row["_pos"] = parse_int(row["normalized_pos1"], "normalized_pos1")
        row["_p_lrt"] = parse_float(row["p_lrt"], "p_lrt")
        row["_candidate_alleles"] = allele_pair(row["allele1"], row["allele0"])
        require(row["_chrom"] == normalize_chrom(row["chrom"]), "Candidate original and standardized chromosomes differ")
        require(row["_pos"] == parse_int(row["pos"], "pos"), "Candidate original and standardized coordinates differ")
        require(row["_candidate_alleles"] == row["normalized_allele_pair"], "candidate standard allele mismatch")
        duplicate_key = (row["trait"], row["candidate_match_key"])
        require(duplicate_key not in seen_trait_candidate_keys, f"the same trait candidate duplicate key: {duplicate_key}")
        seen_trait_candidate_keys.add(duplicate_key)
        require(row["candidate_row_id"] not in seen_candidate_row_ids, f"candidate_row_id duplicate: {row['candidate_row_id']}")
        seen_candidate_row_ids.add(row["candidate_row_id"])
        if row["ld_compatibility_status"] == EXACT_STATUS:
            require(row["eligible_for_primary_ld"] == "yes", "exactly matched candidates are not marked eligible=yes")
            variant_id = row["final_ld_variant_id"]
            require(variant_id in bim_variants, f"exact candidate ID not in final BIM: {variant_id}")
            variant = bim_variants[variant_id]
            require(variant.chrom == row["_chrom"] and variant.pos == row["_pos"], f"Candidate and BIM coordinates differ: {variant_id}")
            require(allele_pair(variant.allele1, variant.allele2) == row["_candidate_alleles"], f"Candidate and BIM alleles differ: {variant_id}")
            exact_rows.append(row)
        else:
            require(row["eligible_for_primary_ld"] == "no", "nonexact candidates are not marked eligible=no")

    lead_ids = sorted({row["final_ld_variant_id"] for row in exact_rows})
    require(exact_rows and lead_ids, "No candidates are available for the main LD panel")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    lead_list = args.output_dir / "unique_ld_seed_variants.txt"
    with lead_list.open("w", encoding="utf-8") as handle:
        for variant_id in lead_ids:
            handle.write(f"{variant_id}\n")

    log(
        f"read candidate {len(compatibility_rows)} row, exact match {len(exact_rows)} row, "
        f"unique LD seed {len(lead_ids)} items; final panel {sample_count} dogs/{bim_row_count} variants"
    )
    log("stage 3/7: use PLINK2 compute all seed with same chromosome marker of height LD pair")
    ld_output_prefix = args.output_dir / "candidate_seed_high_ld_r2"
    command = [
        args.plink2_bin,
        "--bfile",
        args.ld_prefix,
        "--dog",
        "--allow-extra-chr",
        "--ld-snp-list",
        lead_list,
        "--r2-unphased",
        "cols=chrom,pos,id",
        "--ld-window-kb",
        args.ld_window_kb,
        "--ld-window",
        args.ld_window_variant_count,
        "--ld-window-r2",
        args.ld_r2_threshold,
        "--threads",
        args.threads,
        "--out",
        ld_output_prefix,
    ]
    run_command(command, "PLINK2 unphased r2 locus boundaries")
    vcor_path = Path(f"{ld_output_prefix}.vcor")
    ld_members, pair_count = parse_vcor(
        vcor_path, set(lead_ids), bim_variants, args.ld_r2_threshold
    )

    log("stage 4/7: first define primary region, then process exploratory region")
    seed_regions_by_scope: dict[str, list[dict]] = {
        PRIMARY_SCOPE: [],
        EXPLORATORY_SCOPE: [],
    }
    seed_region_by_candidate_id = {}
    for row in exact_rows:
        seed_region = make_seed_region(
            row,
            ld_members,
            bim_variants,
            chrom_bounds,
            args.fallback_window_bp,
        )
        seed_regions_by_scope[row["analysis_scope"]].append(seed_region)
        seed_region_by_candidate_id[row["candidate_row_id"]] = seed_region

    primary_regions = merge_regions(seed_regions_by_scope[PRIMARY_SCOPE])
    assign_locus_ids(primary_regions, PRIMARY_SCOPE)
    exploratory_merged = merge_regions(seed_regions_by_scope[EXPLORATORY_SCOPE])
    retained_exploratory = []
    exploratory_forced_primary: dict[str, dict] = {}
    for region in exploratory_merged:
        primary_match = choose_primary_overlap(region, primary_regions)
        if primary_match is None:
            retained_exploratory.append(region)
        else:
            for component_row in region["component_rows"]:
                exploratory_forced_primary[component_row["candidate_row_id"]] = primary_match
    assign_locus_ids(retained_exploratory, EXPLORATORY_SCOPE)
    all_regions = primary_regions + retained_exploratory
    require(all_regions, "no generate any locus")

    exact_component_assignment: dict[str, dict] = {}
    for region in all_regions:
        component_traits = {row["trait"] for row in region["component_rows"]}
        component_scopes = {row["analysis_scope"] for row in region["component_rows"]}
        require(len(component_traits) == 1, f"locus unexpected inclusion of multiple trait: {region['locus_id']}")
        require(len(component_scopes) == 1, f"locus unexpected inclusion of multiple candidate branch: {region['locus_id']}")
        require(
            all(
                region["analysis_start"] <= row["_pos"] <= region["analysis_end"]
                for row in region["component_rows"]
            ),
            f"locus the analysis interval does not contain all component lead: {region['locus_id']}",
        )
        for component_row in region["component_rows"]:
            exact_component_assignment[component_row["candidate_row_id"]] = region
    exact_component_assignment.update(exploratory_forced_primary)
    require(len(exact_component_assignment) == len(exact_rows), "not all exact candidate matches obtained locus assignment")

    for index, left in enumerate(all_regions):
        for right in all_regions[index + 1 :]:
            if left["representative_lead"]["trait"] != right["representative_lead"]["trait"]:
                continue
            require(
                not core_intervals_overlap(left, right),
                f"final locus height LD core regions overlap: {left['locus_id']} {right['locus_id']}",
            )

    log("stage 5/7: create all candidate, missing candidate and manuscript ROI of assignment audit")
    assignments = []
    locus_assignment_rows: defaultdict[str, list[dict]] = defaultdict(list)
    for row in compatibility_rows:
        exact = row["ld_compatibility_status"] == EXACT_STATUS
        exact_region = exact_component_assignment.get(row["candidate_row_id"])

        matching_regions = (
            [exact_region]
            if exact_region is not None
            else find_containing_loci(row, all_regions)
        )
        matching_locus_ids = "|".join(
            region["locus_id"] for region in matching_regions
        )
        matching_locus_count = len(matching_regions)
        is_ambiguous_missing = not exact and matching_locus_count > 1
        if is_ambiguous_missing:
            assignment_status = "MISSING_CANDIDATE_AMBIGUOUS_OVERLAPPING_ANALYSIS_WINDOWS"
            assigned_locus_id = matching_locus_ids
            assigned_scope = "|".join(
                sorted({region["analysis_scope"] for region in matching_regions})
            )
            is_representative_lead = "no"
            assignment_ambiguity = "overlapping_analysis_windows"
        elif not matching_regions:
            assignment_status = "MISSING_CANDIDATE_UNRESOLVED_NO_PRIMARY_LD_LOCUS"
            assigned_locus_id = ""
            assigned_scope = ""
            is_representative_lead = "no"
            assignment_ambiguity = "none"
        else:
            region = matching_regions[0]
            assigned_locus_id = region["locus_id"]
            assigned_scope = region["analysis_scope"]
            representative = region["representative_lead"]
            is_representative_lead = (
                "yes" if row["candidate_row_id"] == representative["candidate_row_id"] else "no"
            )
            if exact and row["analysis_scope"] == PRIMARY_SCOPE:
                assignment_status = "PRIMARY_LEAD" if is_representative_lead == "yes" else "ABSORBED_BY_PRIMARY_LOCUS"
            elif exact and row["analysis_scope"] == EXPLORATORY_SCOPE and assigned_scope == PRIMARY_SCOPE:
                assignment_status = "EXPLORATORY_LD_REGION_ABSORBED_BY_PRIMARY_LOCUS"
            elif exact and row["analysis_scope"] == EXPLORATORY_SCOPE:
                assignment_status = "EXPLORATORY_LEAD" if is_representative_lead == "yes" else "ABSORBED_BY_EXPLORATORY_LOCUS"
            elif assigned_scope == PRIMARY_SCOPE:
                assignment_status = "MISSING_CANDIDATE_WITHIN_PRIMARY_LOCUS"
            else:
                assignment_status = "MISSING_CANDIDATE_WITHIN_EXPLORATORY_LOCUS"
            assignment_ambiguity = "none"
            locus_assignment_rows[assigned_locus_id].append(row)
        assignments.append(
            {
                **row,
                "assigned_locus_id": assigned_locus_id,
                "assigned_locus_scope": assigned_scope,
                "locus_assignment_status": assignment_status,
                "candidate_is_final_locus_representative_lead": is_representative_lead,
                "matching_locus_ids": matching_locus_ids,
                "matching_locus_count": matching_locus_count,
                "assignment_ambiguity": assignment_ambiguity,
            }
        )

    pnas_by_locus: defaultdict[str, list[str]] = defaultdict(list)
    for paper_row in pnas_rows:
        paper_chrom = normalize_chrom(paper_row["chrom"])
        paper_start, paper_end, _ = pnas_roi_bounds(paper_row)
        for region in all_regions:
            if region["representative_lead"]["trait"] != paper_row["trait"] or region["chrom"] != paper_chrom:
                continue
            if interval_overlap(region["analysis_start"], region["analysis_end"], paper_start, paper_end) > 0:
                pnas_by_locus[region["locus_id"]].append(paper_row["table1_locus_id"])

    analysis_window_overlap_ids: defaultdict[str, set[str]] = defaultdict(set)
    analysis_window_overlap_rows = []
    for index, left in enumerate(all_regions):
        for right in all_regions[index + 1 :]:
            if left["representative_lead"]["trait"] != right["representative_lead"]["trait"]:
                continue
            overlap = interval_overlap(
                left["analysis_start"],
                left["analysis_end"],
                right["analysis_start"],
                right["analysis_end"],
            )
            if left["chrom"] != right["chrom"] or overlap == 0:
                continue
            analysis_window_overlap_ids[left["locus_id"]].add(right["locus_id"])
            analysis_window_overlap_ids[right["locus_id"]].add(left["locus_id"])
            analysis_window_overlap_rows.append(
                {
                    "trait": left["representative_lead"]["trait"],
                    "chrom": left["chrom"],
                    "left_locus_id": left["locus_id"],
                    "left_scope": left["analysis_scope"],
                    "left_region_source": left["region_source"],
                    "right_locus_id": right["locus_id"],
                    "right_scope": right["analysis_scope"],
                    "right_region_source": right["region_source"],
                    "analysis_window_overlap_bp": overlap,
                    "core_interval_overlap_bp": interval_overlap(
                        left["core_start"],
                        left["core_end"],
                        right["core_start"],
                        right["core_end"],
                    ),
                    "interpretation": "SEPARATE_LD_CORES_WITH_OVERLAPPING_ANALYSIS_WINDOWS",
                }
            )

    locus_rows = []
    lead_rows = []
    for region in sorted(
        all_regions,
        key=lambda item: (
            0 if item["analysis_scope"] == PRIMARY_SCOPE else 1,
            item["representative_lead"]["trait"],
            chrom_sort_key(item["chrom"]),
            item["analysis_start"],
        ),
    ):
        lead = region["representative_lead"]
        assigned = locus_assignment_rows[region["locus_id"]]
        component_ids = sorted(
            {row["final_ld_variant_id"] for row in region["component_rows"]}
        )
        lead_member_count = len(ld_members[lead["final_ld_variant_id"]])
        locus_row = {
            "locus_id": region["locus_id"],
            "trait": lead["trait"],
            "trait_label": trait_by_id[lead["trait"]]["trait_label"],
            "analysis_type": lead["analysis_type"],
            "analysis_scope": region["analysis_scope"],
            "chrom": region["chrom"],
            "core_start": region["core_start"],
            "core_end": region["core_end"],
            "analysis_start": region["analysis_start"],
            "analysis_end": region["analysis_end"],
            "analysis_length_bp": region["analysis_end"] - region["analysis_start"] + 1,
            "region_source": region["region_source"],
            "component_region_sources": "|".join(sorted(region["component_region_sources"])),
            "merge_basis": region["merge_basis"],
            "merge_event_count": region["merge_event_count"],
            "ld_core_overlap_merge_count": region["ld_core_overlap_merge_count"],
            "no_proxy_core_overlap_merge_count": region["no_proxy_core_overlap_merge_count"],
            "analysis_window_overlaps_other_locus": (
                "yes" if analysis_window_overlap_ids[region["locus_id"]] else "no"
            ),
            "overlapping_analysis_window_locus_ids": "|".join(
                sorted(analysis_window_overlap_ids[region["locus_id"]])
            ),
            "ld_r2_threshold": format_number(args.ld_r2_threshold),
            "fallback_window_bp": args.fallback_window_bp,
            "lead_candidate_row_id": lead["candidate_row_id"],
            "lead_variant_id": lead["final_ld_variant_id"],
            "lead_pos": lead["_pos"],
            "lead_allele1": lead["allele1"],
            "lead_allele0": lead["allele0"],
            "lead_p_lrt": lead["p_lrt"],
            "lead_p_wald": lead["p_wald"],
            "lead_p_score": lead["p_score"],
            "lead_beta": lead["beta"],
            "lead_se": lead["se"],
            "lead_af": lead["af"],
            "lead_n_miss": lead["n_miss"],
            "lead_p_lrt_tier": lead["p_lrt_tier"],
            "lead_high_ld_member_count_including_self": lead_member_count,
            "lead_high_ld_proxy_count": lead_member_count - 1,
            "component_seed_count": len(component_ids),
            "component_seed_variant_ids": "|".join(component_ids),
            "fallback_component_count": region["fallback_component_count"],
            "assigned_candidate_rows": len(assigned),
            "assigned_missing_candidate_rows": sum(
                row["ld_compatibility_status"] != EXACT_STATUS for row in assigned
            ),
            "pnas_table1_locus_ids": "|".join(sorted(pnas_by_locus[region["locus_id"]])),
        }
        locus_rows.append(locus_row)
        lead_rows.append(
            {
                field: locus_row[field]
                for field in (
                    "locus_id",
                    "trait",
                    "trait_label",
                    "analysis_type",
                    "analysis_scope",
                    "chrom",
                    "lead_candidate_row_id",
                    "lead_variant_id",
                    "lead_pos",
                    "lead_allele1",
                    "lead_allele0",
                    "lead_p_lrt",
                    "lead_p_lrt_tier",
                    "component_seed_count",
                    "component_seed_variant_ids",
                    "region_source",
                    "merge_basis",
                    "analysis_window_overlaps_other_locus",
                    "analysis_start",
                    "analysis_end",
                )
            }
        )

    pnas_reconciliation = []
    for paper_row in pnas_rows:
        paper_chrom = normalize_chrom(paper_row["chrom"])
        paper_pos = parse_int(paper_row["lead_pos"], "lead_pos")
        paper_alleles = allele_pair(paper_row["table1_allele_a"], paper_row["table1_allele_b"])
        paper_candidates = [
            row
            for row in compatibility_rows
            if row["trait"] == paper_row["trait"]
            and row["_chrom"] == paper_chrom
            and row["_pos"] == paper_pos
            and row["_candidate_alleles"] == paper_alleles
        ]
        require(len(paper_candidates) <= 1, f"manuscript lead matching multiple candidate: {paper_row['table1_locus_id']}")
        paper_candidate = paper_candidates[0] if paper_candidates else None
        roi_start, roi_end, roi_defined = pnas_roi_bounds(paper_row)
        matching_regions = []
        roi_overlap_intervals = []
        for region in all_regions:
            if region["representative_lead"]["trait"] != paper_row["trait"] or region["chrom"] != paper_chrom:
                continue
            overlap = interval_overlap(region["analysis_start"], region["analysis_end"], roi_start, roi_end)
            if overlap > 0:
                matching_regions.append(region)
                roi_overlap_intervals.append(
                    (
                        max(region["analysis_start"], roi_start),
                        min(region["analysis_end"], roi_end),
                    )
                )
        roi_length = roi_end - roi_start + 1

        overlap_total = interval_union_length(roi_overlap_intervals)
        require(
            0 <= overlap_total <= roi_length,
            f"PNAS ROI covered length exceeds span: {paper_row['table1_locus_id']} "
            f"overlap={overlap_total} roi_length={roi_length}",
        )
        overlap_fraction = overlap_total / roi_length
        pnas_reconciliation.append(
            {
                "table1_locus_id": paper_row["table1_locus_id"],
                "trait": paper_row["trait"],
                "chrom": paper_chrom,
                "paper_lead_pos": paper_pos,
                "paper_lead_alleles": paper_alleles,
                "paper_p_displayed": paper_row["table1_p_displayed"],
                "raw_p_lrt": paper_row["raw_p_lrt"],
                "proximal_candidate_gene": paper_row["proximal_candidate_gene"],
                "official_roi_defined": "yes" if roi_defined else "no",
                "official_roi_start": paper_row["official_roi_start"],
                "official_roi_end": paper_row["official_roi_end"],
                "paper_lead_ld_compatibility_status": "" if paper_candidate is None else paper_candidate["ld_compatibility_status"],
                "paper_lead_removal_reason": "" if paper_candidate is None else paper_candidate["removal_reason"],
                "matched_current_locus_ids": "|".join(sorted(region["locus_id"] for region in matching_regions)),
                "matched_current_locus_count": len(matching_regions),
                "overlap_bp": overlap_total,
                "overlap_fraction_of_official_roi": format_number(overlap_fraction),
                "reconciliation_status": "MATCHED_CURRENT_LOCUS" if matching_regions else "NO_CURRENT_LOCUS_OVERLAP",
            }
        )

    missing_assignments = [
        row for row in assignments if row["ld_compatibility_status"] != EXACT_STATUS
    ]
    region_by_locus_id = {region["locus_id"]: region for region in all_regions}
    require(len(region_by_locus_id) == len(all_regions), "final locus_id is not unique")
    non_finemappable_rows = []
    for row in missing_assignments:
        candidate_key = (
            row["trait"],
            row["_chrom"],
            row["_pos"],
            row["_candidate_alleles"],
        )
        paper_row = pnas_by_candidate_key.get(candidate_key)
        paper_roi_start = ""
        paper_roi_end = ""
        paper_roi_defined = False
        if paper_row is not None:
            paper_roi_start_value, paper_roi_end_value, paper_roi_defined = pnas_roi_bounds(
                paper_row
            )
            if paper_roi_defined:
                paper_roi_start = paper_roi_start_value
                paper_roi_end = paper_roi_end_value
        chrom_min, chrom_max = chrom_bounds[row["_chrom"]]
        assigned_locus_id = row["assigned_locus_id"]
        matching_regions = [
            region_by_locus_id[locus_id]
            for locus_id in row["matching_locus_ids"].split("|")
            if locus_id
        ]
        if row["assignment_ambiguity"] == "overlapping_analysis_windows":
            analysis_start = min(region["analysis_start"] for region in matching_regions)
            analysis_end = max(region["analysis_end"] for region in matching_regions)
            region_source = "AMBIGUOUS_OVERLAPPING_ANALYSIS_WINDOWS"
            regional_analysis_status = (
                "EXISTING_LOCUS_MATCH_IS_AMBIGUOUS_DO_NOT_ASSIGN_UNIQUELY"
            )
        elif assigned_locus_id:
            region = matching_regions[0]
            analysis_start = region["analysis_start"]
            analysis_end = region["analysis_end"]
            region_source = "ASSIGNED_EXISTING_LD_LOCUS"
            regional_analysis_status = "EXISTING_LOCUS_FINE_MAPPING_WITHOUT_THIS_CANDIDATE"
        else:
            analysis_start = max(chrom_min, row["_pos"] - args.fallback_window_bp)
            analysis_end = min(chrom_max, row["_pos"] + args.fallback_window_bp)
            region_source = "SUMMARY_CANDIDATE_PLUS_MINUS_WINDOW_NO_LD_GENOTYPE"
            regional_analysis_status = "PHYSICAL_WINDOW_REQUIRES_VARIANT_LD_INTERSECTION"

        removal_reason = row["removal_reason"]
        if "HWE" in removal_reason:
            recommended_handling = "REPORT_AND_TEST_HWE_SENSITIVITY_PANEL"
        elif "MAF" in removal_reason:
            recommended_handling = "REPORT_AND_TEST_MAF_SENSITIVITY_PANEL"
        elif row["ld_compatibility_status"] == "NOT_IN_ORIGINAL_LD":
            recommended_handling = "REPORT_ONLY_NO_LD_GENOTYPE_AVAILABLE"
        else:
            recommended_handling = "REPORT_AND_REVIEW_UPSTREAM_LD_EXCLUSION"

        non_finemappable_rows.append(
            {
                "candidate_row_id": row["candidate_row_id"],
                "trait": row["trait"],
                "trait_label": trait_by_id[row["trait"]]["trait_label"],
                "analysis_type": row["analysis_type"],
                "analysis_scope": row["analysis_scope"],
                "p_lrt_tier": row["p_lrt_tier"],
                "chrom": row["_chrom"],
                "candidate_pos": row["_pos"],
                "candidate_alleles": row["_candidate_alleles"],
                "candidate_rs": row["rs"],
                "candidate_p_lrt": row["p_lrt"],
                "candidate_beta": row["beta"],
                "candidate_se": row["se"],
                "ld_compatibility_status": row["ld_compatibility_status"],
                "removal_reason": removal_reason,
                "candidate_in_ld_matrix": "no",
                "candidate_pip_eligible": "no",
                "assigned_locus_id": assigned_locus_id,
                "assigned_locus_scope": row["assigned_locus_scope"],
                "matching_locus_ids": row["matching_locus_ids"],
                "matching_locus_count": row["matching_locus_count"],
                "assignment_ambiguity": row["assignment_ambiguity"],
                "analysis_start": analysis_start,
                "analysis_end": analysis_end,
                "analysis_length_bp": analysis_end - analysis_start + 1,
                "region_source": region_source,
                "regional_analysis_status": regional_analysis_status,
                "pnas_table1_locus_id": "" if paper_row is None else paper_row["table1_locus_id"],
                "pnas_official_roi_defined": "yes" if paper_roi_defined else "no",
                "pnas_official_roi_start": paper_roi_start,
                "pnas_official_roi_end": paper_roi_end,
                "pnas_roi_role": "audit_only_not_region_override",
                "recommended_handling": recommended_handling,
            }
        )
    require(
        len(non_finemappable_rows) == len(missing_assignments),
        "the table of candidate regions not directly fine-mappable does not cover all missing candidates",
    )

    log("stage 6/7: write locus, candidate assignment, manuscript review and lineage file")
    locus_fields = [
        "locus_id", "trait", "trait_label", "analysis_type", "analysis_scope",
        "chrom", "core_start", "core_end", "analysis_start", "analysis_end",
        "analysis_length_bp", "region_source", "component_region_sources", "merge_basis",
        "merge_event_count", "ld_core_overlap_merge_count",
        "no_proxy_core_overlap_merge_count", "analysis_window_overlaps_other_locus",
        "overlapping_analysis_window_locus_ids",
        "ld_r2_threshold", "fallback_window_bp",
        "lead_candidate_row_id", "lead_variant_id", "lead_pos", "lead_allele1",
        "lead_allele0", "lead_p_lrt", "lead_p_wald", "lead_p_score", "lead_beta",
        "lead_se", "lead_af", "lead_n_miss", "lead_p_lrt_tier",
        "lead_high_ld_member_count_including_self", "lead_high_ld_proxy_count",
        "component_seed_count", "component_seed_variant_ids", "fallback_component_count",
        "assigned_candidate_rows", "assigned_missing_candidate_rows", "pnas_table1_locus_ids",
    ]
    lead_fields = [
        "locus_id", "trait", "trait_label", "analysis_type", "analysis_scope", "chrom",
        "lead_candidate_row_id", "lead_variant_id", "lead_pos", "lead_allele1",
        "lead_allele0", "lead_p_lrt", "lead_p_lrt_tier", "component_seed_count",
        "component_seed_variant_ids", "region_source", "merge_basis",
        "analysis_window_overlaps_other_locus", "analysis_start", "analysis_end",
    ]
    assignment_fields = compatibility_columns + [
        "assigned_locus_id",
        "assigned_locus_scope",
        "locus_assignment_status",
        "candidate_is_final_locus_representative_lead",
        "matching_locus_ids",
        "matching_locus_count",
        "assignment_ambiguity",
    ]
    pnas_fields_out = [
        "table1_locus_id", "trait", "chrom", "paper_lead_pos", "paper_lead_alleles",
        "paper_p_displayed", "raw_p_lrt", "proximal_candidate_gene", "official_roi_defined",
        "official_roi_start", "official_roi_end", "paper_lead_ld_compatibility_status",
        "paper_lead_removal_reason", "matched_current_locus_ids",
        "matched_current_locus_count", "overlap_bp", "overlap_fraction_of_official_roi",
        "reconciliation_status",
    ]
    write_tsv(args.output_dir / "all_locus_regions.tsv", locus_fields, locus_rows)
    write_tsv(
        args.output_dir / "primary_locus_regions.tsv",
        locus_fields,
        [row for row in locus_rows if row["analysis_scope"] == PRIMARY_SCOPE],
    )
    write_tsv(
        args.output_dir / "exploratory_locus_regions.tsv",
        locus_fields,
        [row for row in locus_rows if row["analysis_scope"] == EXPLORATORY_SCOPE],
    )
    write_tsv(
        args.output_dir / "primary_lead_variants.tsv",
        lead_fields,
        [row for row in lead_rows if row["analysis_scope"] == PRIMARY_SCOPE],
    )
    write_tsv(
        args.output_dir / "exploratory_lead_variants.tsv",
        lead_fields,
        [row for row in lead_rows if row["analysis_scope"] == EXPLORATORY_SCOPE],
    )
    merge_audit_rows = []
    for row in sorted(
        exact_rows,
        key=lambda item: (
            item["trait"],
            chrom_sort_key(item["_chrom"]),
            item["_pos"],
            item["candidate_row_id"],
        ),
    ):
        seed_region = seed_region_by_candidate_id[row["candidate_row_id"]]
        final_region = exact_component_assignment[row["candidate_row_id"]]
        merge_audit_rows.append(
            {
                "candidate_row_id": row["candidate_row_id"],
                "trait": row["trait"],
                "candidate_scope": row["analysis_scope"],
                "seed_variant_id": row["final_ld_variant_id"],
                "chrom": row["_chrom"],
                "seed_pos": row["_pos"],
                "seed_core_start": seed_region["core_start"],
                "seed_core_end": seed_region["core_end"],
                "seed_analysis_start": seed_region["analysis_start"],
                "seed_analysis_end": seed_region["analysis_end"],
                "seed_region_source": seed_region["region_source"],
                "final_locus_id": final_region["locus_id"],
                "final_locus_scope": final_region["analysis_scope"],
                "final_core_start": final_region["core_start"],
                "final_core_end": final_region["core_end"],
                "final_analysis_start": final_region["analysis_start"],
                "final_analysis_end": final_region["analysis_end"],
                "final_merge_basis": final_region["merge_basis"],
                "absorbed_exploratory_into_primary": (
                    "yes"
                    if row["analysis_scope"] == EXPLORATORY_SCOPE
                    and final_region["analysis_scope"] == PRIMARY_SCOPE
                    else "no"
                ),
            }
        )
    write_tsv(
        args.output_dir / "locus_merge_audit.tsv",
        [
            "candidate_row_id", "trait", "candidate_scope", "seed_variant_id",
            "chrom", "seed_pos", "seed_core_start", "seed_core_end",
            "seed_analysis_start", "seed_analysis_end", "seed_region_source",
            "final_locus_id", "final_locus_scope", "final_core_start",
            "final_core_end", "final_analysis_start", "final_analysis_end",
            "final_merge_basis", "absorbed_exploratory_into_primary",
        ],
        merge_audit_rows,
    )
    write_tsv(
        args.output_dir / "locus_region_source_audit.tsv",
        [
            "locus_id", "trait", "analysis_scope", "chrom", "core_start", "core_end",
            "analysis_start", "analysis_end", "region_source",
            "component_region_sources", "fallback_component_count", "merge_basis",
            "analysis_window_overlaps_other_locus",
            "overlapping_analysis_window_locus_ids",
        ],
        locus_rows,
    )
    write_tsv(
        args.output_dir / "overlapping_analysis_windows_audit.tsv",
        [
            "trait", "chrom", "left_locus_id", "left_scope",
            "left_region_source", "right_locus_id", "right_scope",
            "right_region_source", "analysis_window_overlap_bp",
            "core_interval_overlap_bp", "interpretation",
        ],
        analysis_window_overlap_rows,
    )
    write_tsv(args.output_dir / "candidate_locus_assignment.tsv", assignment_fields, assignments)
    write_tsv(
        args.output_dir / "missing_candidate_locus_resolution.tsv",
        assignment_fields,
        missing_assignments,
    )
    write_tsv(
        args.output_dir / "non_finemappable_candidate_regions.tsv",
        [
            "candidate_row_id", "trait", "trait_label", "analysis_type",
            "analysis_scope", "p_lrt_tier", "chrom", "candidate_pos",
            "candidate_alleles", "candidate_rs", "candidate_p_lrt", "candidate_beta",
            "candidate_se", "ld_compatibility_status", "removal_reason",
            "candidate_in_ld_matrix", "candidate_pip_eligible", "assigned_locus_id",
            "assigned_locus_scope", "matching_locus_ids", "matching_locus_count",
            "assignment_ambiguity", "analysis_start", "analysis_end",
            "analysis_length_bp", "region_source", "regional_analysis_status",
            "pnas_table1_locus_id", "pnas_official_roi_defined",
            "pnas_official_roi_start", "pnas_official_roi_end", "pnas_roi_role",
            "recommended_handling",
        ],
        non_finemappable_rows,
    )
    write_tsv(
        args.output_dir / "exploratory_candidates_absorbed_by_primary.tsv",
        assignment_fields,
        [
            row
            for row in assignments
            if row["analysis_scope"] == EXPLORATORY_SCOPE
            and row["assigned_locus_scope"] == PRIMARY_SCOPE
        ],
    )
    write_tsv(
        args.output_dir / "pnas_table1_locus_reconciliation.tsv",
        pnas_fields_out,
        pnas_reconciliation,
    )

    proxy_rows = []
    for lead_id in lead_ids:
        lead_variant = bim_variants[lead_id]
        for member_id, r2_value in sorted(
            ld_members[lead_id].items(),
            key=lambda item: (bim_variants[item[0]].pos, item[0]),
        ):
            member = bim_variants[member_id]
            proxy_rows.append(
                {
                    "lead_variant_id": lead_id,
                    "lead_chrom": lead_variant.chrom,
                    "lead_pos": lead_variant.pos,
                    "member_variant_id": member_id,
                    "member_chrom": member.chrom,
                    "member_pos": member.pos,
                    "unphased_r2": format_number(r2_value),
                    "is_self": "yes" if member_id == lead_id else "no",
                }
            )
    write_tsv_gz(
        args.output_dir / "lead_ld_proxy_variants.tsv.gz",
        [
            "lead_variant_id", "lead_chrom", "lead_pos", "member_variant_id",
            "member_chrom", "member_pos", "unphased_r2", "is_self",
        ],
        proxy_rows,
    )

    parameter_rows = [
        ["ld_statistic", "unphased_hardcall_dosage_r2"],
        ["ld_r2_threshold", format_number(args.ld_r2_threshold)],
        ["ld_window_kb", args.ld_window_kb],
        ["ld_window_variant_count", args.ld_window_variant_count],
        ["fallback_window_bp", args.fallback_window_bp],
        ["scope_priority", "primary_then_exploratory"],
        ["paper_roi_role", "audit_only_not_region_override"],
        ["missing_candidate_seed_policy", "never_impute_or_substitute_ld_seed"],
        ["missing_candidate_regional_policy", "existing_locus_else_physical_window"],
        ["bed_magic_hex", bed_magic_hex],
        ["bed_actual_size_bytes", bed_actual_size],
        ["bed_expected_size_bytes", bed_expected_size],
        ["plink2_version", version],
    ]
    write_tsv(
        args.output_dir / "pipeline_parameters.tsv",
        ["parameter", "value"],
        ({"parameter": key, "value": value} for key, value in parameter_rows),
    )
    lineage_rows = [
        ["compatibility_table", args.compatibility_table, args.audit_producer],
        ["trait_metadata", args.trait_metadata, args.candidate_producer],
        ["pnas_table1_loci", args.pnas_table1_loci, args.candidate_producer],
        ["step16_completion", args.step16_completion, args.audit_producer],
        ["ld_prefix", args.ld_prefix, args.ld_producer],
        ["step15_completion", args.step15_completion, args.ld_producer],
        ["plink2_binary", args.plink2_bin, "local_PLINK2_installation"],
    ]
    write_tsv(
        args.output_dir / "input_lineage.tsv",
        ["role", "path", "producer"],
        (
            {"role": role, "path": path, "producer": producer}
            for role, path, producer in lineage_rows
        ),
    )

    log("stage 7/7: check terminal-state invariants and write the completion marker")
    assignment_status_counts = Counter(row["locus_assignment_status"] for row in assignments)
    unresolved_missing = assignment_status_counts[
        "MISSING_CANDIDATE_UNRESOLVED_NO_PRIMARY_LD_LOCUS"
    ]
    ambiguous_missing = assignment_status_counts[
        "MISSING_CANDIDATE_AMBIGUOUS_OVERLAPPING_ANALYSIS_WINDOWS"
    ]
    fallback_loci = sum(row["fallback_component_count"] > 0 for row in locus_rows)
    no_proxy_core_merged_loci = sum(
        row["no_proxy_core_overlap_merge_count"] > 0 for row in locus_rows
    )
    non_finemappable_source_counts = Counter(
        row["region_source"] for row in non_finemappable_rows
    )
    unmatched_pnas = sum(
        row["reconciliation_status"] == "NO_CURRENT_LOCUS_OVERLAP"
        for row in pnas_reconciliation
    )
    overall_status = (
        "PASS_WITH_WARNINGS"
        if unresolved_missing > 0
        or ambiguous_missing > 0
        or fallback_loci > 0
        or unmatched_pnas > 0
        else "PASS"
    )
    summary_rows = [
        ["overall_status", overall_status, "PASS", "LD locus definition and audit completed."],
        ["candidate_rows", len(compatibility_rows), "PASS", "Step 1.6 compatibility records."],
        ["exact_candidate_rows", len(exact_rows), "PASS", "Eligible as final main-LD seeds."],
        ["unique_ld_seed_variant_ids", len(lead_ids), "PASS", "BIM IDs deduplicated across traits."],
        ["final_ld_samples", sample_count, "PASS", "All samples are coded as founders."],
        ["final_ld_variants", bim_row_count, "PASS", "Final main-LD BIM."],
        ["bed_magic_hex", bed_magic_hex, "PASS", "The PLINK SNP-major header should be 6c1b01."],
        ["bed_actual_size_bytes", bed_actual_size, "PASS", "Consistent with the size expected from FAM/BIM."],
        ["bed_expected_size_bytes", bed_expected_size, "PASS", "3+ceil(N/4)*M"],
        ["plink2_high_ld_pair_rows", pair_count, "INFO", "Pairs in VCOR with r2 reaching the threshold."],
        ["primary_locus_count", len(primary_regions), "PASS", "Primary regions defined first."],
        ["exploratory_locus_count", len(retained_exploratory), "PASS", "After excluding primary-region coverage."],
        ["fallback_locus_count", fallback_loci, "WARN" if fallback_loci else "PASS", "At least one component has no r2 proxy."],
        ["no_proxy_lead_position_merged_locus_count", no_proxy_core_merged_loci, "WARN" if no_proxy_core_merged_loci else "PASS", "Same lead position, but at least one component has no r2 proxy."],
        ["overlapping_analysis_window_pair_count", len(analysis_window_overlap_rows), "WARN" if analysis_window_overlap_rows else "PASS", "Distinct LD cores remain separate, although alternative analysis windows may intersect."],
        ["missing_candidate_rows", len(missing_assignments), "INFO", "Not forcibly restored to the main LD panel."],
        ["non_finemappable_candidate_region_rows", len(non_finemappable_rows), "PASS", "An individual disposition record is retained for every missing candidate."],
        ["ambiguous_missing_candidate_rows", ambiguous_missing, "WARN" if ambiguous_missing else "PASS", "Retain all matches to overlapping analysis windows without forcing a unique assignment."],
        ["non_finemappable_assigned_existing_locus_rows", non_finemappable_source_counts["ASSIGNED_EXISTING_LD_LOCUS"], "INFO", "The region is analyzable, but this candidate is not eligible for a PIP."],
        ["non_finemappable_pnas_lead_rows", sum(bool(row["pnas_table1_locus_id"]) for row in non_finemappable_rows), "INFO", "Published ROI are used only for auditing and do not override the current regions."],
        ["non_finemappable_physical_window_rows", non_finemappable_source_counts["SUMMARY_CANDIDATE_PLUS_MINUS_WINDOW_NO_LD_GENOTYPE"], "INFO", "Preset physical window when no LD seed is available."],
        ["unresolved_missing_candidate_rows", unresolved_missing, "WARN" if unresolved_missing else "PASS", "Does not fall within any current locus."],
        ["pnas_table1_loci", len(pnas_rows), "INFO", "For review only."],
        ["pnas_table1_unmatched_loci", unmatched_pnas, "WARN" if unmatched_pnas else "PASS", "No overlap with a current region."],
        ["nonfounder_samples", nonfounder_count, "PASS", "PLINK2 LD does not silently drop samples."],
        ["step15_status", "COMPLETED", "PASS", str(args.step15_completion)],
        ["step16_status", step16_completion["overall_status"], "PASS", str(args.step16_completion)],
    ]
    for status_name, count in sorted(assignment_status_counts.items()):
        summary_rows.append([f"assignment_count::{status_name}", count, "INFO", "Candidate row count."])
    for summary_name in ("locus_definition_summary.tsv", "ld_locus_definition_summary.tsv"):
        write_tsv(
            args.output_dir / summary_name,
            ["metric", "value", "status", "details"],
            (
                {"metric": metric, "value": value, "status": status, "details": details}
                for metric, value, status, details in summary_rows
            ),
        )
    completion_path = args.output_dir / "COMPLETED.ok"
    completion_path.write_text(
        "\n".join(
            [
                f"overall_status\t{overall_status}",
                f"candidate_rows\t{len(compatibility_rows)}",
                f"non_finemappable_candidate_rows\t{len(non_finemappable_rows)}",
                f"ambiguous_missing_candidate_rows\t{ambiguous_missing}",
                f"primary_locus_count\t{len(primary_regions)}",
                f"exploratory_locus_count\t{len(retained_exploratory)}",
                f"output_dir\t{args.output_dir}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    log(
        f"completed: status={overall_status}, primary locus={len(primary_regions)}, "
        f"exploratory locus={len(retained_exploratory)}, unresolved missing candidates={unresolved_missing}, "
        f"missing candidates with ambiguous assignment to multiple windows={ambiguous_missing}"
    )
    log("task completed / done")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"ERROR: {exc}")
        sys.exit(1)
