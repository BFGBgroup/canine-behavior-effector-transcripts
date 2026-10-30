#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


AUTOSOME_MIN = 1
AUTOSOME_MAX = 38
MISSING_VALUES = frozenset({"", ".", "NA", "NaN", "nan", "N/A"})

LOCUS_REQUIRED_FIELDS = (
    "locus_id",
    "trait",
    "trait_label",
    "analysis_type",
    "analysis_scope",
    "chrom",
    "core_start",
    "core_end",
    "analysis_start",
    "analysis_end",
    "region_source",
    "lead_variant_id",
    "lead_pos",
    "lead_allele1",
    "lead_allele0",
    "lead_p_lrt",
    "lead_p_lrt_tier",
)
TRAIT_CONTRACT_REQUIRED_FIELDS = (
    "trait",
    "expected_analysis_type",
    "expected_record_count",
    "input_file",
)
TRAIT_METADATA_REQUIRED_FIELDS = (
    "trait",
    "trait_label",
    "analysis_type",
    "n_main",
    "n_sensitivity",
    "locus_selection_p",
    "fine_mapping_z",
)
GWAS_REQUIRED_FIELDS = (
    "chr",
    "rs",
    "ps",
    "n_miss",
    "allele1",
    "allele0",
    "af",
    "beta",
    "se",
    "p_wald",
    "p_lrt",
    "p_score",
)
NON_FINEMAPPABLE_REQUIRED_FIELDS = (
    "candidate_row_id",
    "trait",
    "analysis_scope",
    "candidate_rs",
    "ld_compatibility_status",
    "candidate_pip_eligible",
    "assigned_locus_id",
    "regional_analysis_status",
    "recommended_handling",
)

AUDIT_FIELDS = (
    "locus_id",
    "trait",
    "trait_label",
    "analysis_type",
    "analysis_scope",
    "chrom",
    "pos",
    "gwas_variant_id",
    "gwas_allele1",
    "gwas_allele0",
    "gwas_af_for_allele1",
    "beta",
    "se",
    "z_for_gwas_allele1",
    "p_wald",
    "p_lrt",
    "p_score",
    "n_main",
    "n_miss",
    "n_obs",
    "missing_rate",
    "ld_variant_id",
    "ld_bim_order",
    "ld_a1",
    "ld_a2",
    "ld_a1_frequency",
    "gwas_frequency_aligned_to_ld_a1",
    "frequency_difference_label_orientation",
    "frequency_difference_swapped_orientation",
    "ambiguous_at_cg",
    "z_for_ld_a1",
    "z_orientation",
    "harmonization_status",
    "main_analysis_include",
    "strict_no_ambiguous_include",
    "status_detail",
)

HARMONIZED_FIELDS = (
    "variant_index",
    *AUDIT_FIELDS,
)

ORDER_FIELDS = (
    "variant_index",
    "locus_id",
    "analysis_scope",
    "trait",
    "chrom",
    "pos",
    "ld_variant_id",
    "ld_a1",
    "ld_a2",
    "z_for_ld_a1",
    "n_obs",
    "ambiguous_at_cg",
    "strict_no_ambiguous_include",
)

SUMMARY_FIELDS = (
    "locus_id",
    "trait",
    "trait_label",
    "analysis_type",
    "analysis_scope",
    "chrom",
    "analysis_start",
    "analysis_end",
    "region_source",
    "lead_variant_id",
    "lead_pos",
    "lead_p_lrt",
    "lead_p_lrt_tier",
    "gwas_rows_in_region",
    "harmonized_variant_count",
    "strict_non_ambiguous_variant_count",
    "dropped_variant_count",
    "no_ld_position_count",
    "allele_set_mismatch_count",
    "missing_rate_filtered_count",
    "ambiguous_frequency_resolved_count",
    "ambiguous_frequency_dropped_count",
    "nonambiguous_frequency_warning_count",
    "z_flipped_count",
    "n_locus_median",
    "n_obs_min",
    "n_obs_max",
    "lead_in_harmonized_variants",
    "locus_input_status",
)

MANIFEST_FIELDS = (
    "locus_id",
    "trait",
    "analysis_scope",
    "locus_input_status",
    "harmonized_variant_count",
    "n_locus_median",
    "locus_directory",
    "harmonized_variants",
    "variant_order",
    "allele_harmonization_audit",
    "dropped_variants",
    "locus_metadata",
)


@dataclass(frozen=True)
class TraitInfo:
    trait: str
    trait_label: str
    analysis_type: str
    n_main: int
    n_sensitivity: str
    input_file: Path
    expected_record_count: int


@dataclass(frozen=True)
class Locus:
    locus_id: str
    trait: str
    trait_label: str
    analysis_type: str
    analysis_scope: str
    chrom: str
    core_start: int
    core_end: int
    analysis_start: int
    analysis_end: int
    region_source: str
    lead_variant_id: str
    lead_pos: int
    lead_allele1: str
    lead_allele0: str
    lead_p_lrt: float
    lead_p_lrt_tier: str


@dataclass(frozen=True)
class BimRecord:
    chrom: str
    variant_id: str
    pos: int
    a1: str
    a2: str
    order: int


@dataclass(frozen=True)
class FrequencyRecord:
    chrom: str
    variant_id: str
    pos: int
    ref: str
    alt: str
    alt_frequency: float
    observation_count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Harmonize canine behavioral GWAS summary statistics with the final LD panel by locus, including coordinates, alleles, z-score orientation, and sample size."
    )
    parser.add_argument("--locus-regions", required=True, type=Path)
    parser.add_argument("--non-finemappable-candidates", required=True, type=Path)
    parser.add_argument("--step17-completion", required=True, type=Path)
    parser.add_argument("--trait-input-contract", required=True, type=Path)
    parser.add_argument("--trait-metadata", required=True, type=Path)
    parser.add_argument("--ld-prefix", required=True, type=Path)
    parser.add_argument("--ld-frequency", required=True, type=Path)
    parser.add_argument("--step15-completion", required=True, type=Path)
    parser.add_argument("--step11-producer", required=True, type=Path)
    parser.add_argument("--step12-producer", required=True, type=Path)
    parser.add_argument("--step15-producer", required=True, type=Path)
    parser.add_argument("--step17-producer", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-missing-rate", type=float, default=0.05)
    parser.add_argument("--ambiguous-frequency-tolerance", type=float, default=0.10)
    parser.add_argument("--ambiguous-frequency-separation", type=float, default=0.10)
    parser.add_argument("--min-variants-for-finemapping", type=int, default=2)
    parser.add_argument("--progress-every", type=int, default=100000)
    return parser.parse_args()


def log(message: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S ]")
    print(f"{stamp} {message}", flush=True)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def parse_int(value: object, label: str) -> int:
    text = clean(value)
    require(text not in MISSING_VALUES, f"{label} missing")
    try:
        number = float(text)
    except ValueError as exc:
        raise RuntimeError(f"{label} is not numeric: {text}") from exc
    require(math.isfinite(number) and number.is_integer(), f"{label} is not an integer: {text}")
    return int(number)


def parse_float(value: object, label: str) -> float:
    text = clean(value)
    require(text not in MISSING_VALUES, f"{label} missing")
    try:
        number = float(text)
    except ValueError as exc:
        raise RuntimeError(f"{label} is not numeric: {text}") from exc
    require(math.isfinite(number), f"{label} is not finite numeric value: {text}")
    return number


def maybe_float(value: object) -> float | None:
    text = clean(value)
    if text in MISSING_VALUES:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def maybe_int(value: object) -> int | None:
    text = clean(value)
    if text in MISSING_VALUES:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


def normalize_chrom(value: object) -> str:
    text = clean(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    require(text.isdigit(), f"chromosome is not a numeric autosomal identifier: {value}")
    chrom = int(text)
    require(AUTOSOME_MIN <= chrom <= AUTOSOME_MAX, f"chromosome exceeds 1-38: {value}")
    return str(chrom)


def normalize_allele(value: object) -> str:
    return clean(value).upper()


def allele_pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def is_palindromic_snp(a: str, b: str) -> bool:
    if len(a) != 1 or len(b) != 1:
        return False
    return frozenset((a, b)) in (
        frozenset(("A", "T")),
        frozenset(("C", "G")),
    )


def format_float(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def read_tsv(path: Path, required_fields: Sequence[str]) -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file() and path.stat().st_size > 0, f"input does not exist or is empty: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        missing = [field for field in required_fields if field not in header]
        require(not missing, f"{path} missing fields: {','.join(missing)}")
        return header, [dict(row) for row in reader]


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_gzip_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def read_completion_status(
    path: Path,
    allowed: set[str],
    *,
    fallback_completion_key: str | None = None,
) -> str:
    require(path.is_file() and path.stat().st_size > 0, f"completion marker does not exist or is empty: {path}")
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 1)
            if len(fields) == 2:
                values[fields[0]] = fields[1]
    status = values.get("overall_status", "")
    if not status and fallback_completion_key is not None:
        require(
            clean(values.get(fallback_completion_key)) != "",
            f"upstream completion marker missing {fallback_completion_key}: {path}",
        )
        return "COMPLETED"
    require(status in allowed, f"upstream completion status invalid: {path} -> {status}")
    return status


def load_traits(contract_path: Path, metadata_path: Path) -> dict[str, TraitInfo]:
    _, contract_rows = read_tsv(contract_path, TRAIT_CONTRACT_REQUIRED_FIELDS)
    _, metadata_rows = read_tsv(metadata_path, TRAIT_METADATA_REQUIRED_FIELDS)
    contracts: dict[str, dict[str, str]] = {}
    for row in contract_rows:
        trait = clean(row["trait"])
        require(trait and trait not in contracts, f"trait_input_contract in trait duplicate or is empty: {trait}")
        contracts[trait] = row
    metadata: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        trait = clean(row["trait"])
        require(trait and trait not in metadata, f"trait_metadata in trait duplicate or is empty: {trait}")
        metadata[trait] = row
    require(set(contracts) == set(metadata), "Trait sets differ between the contract and metadata")
    require(len(contracts) == 14, f"actual requirement 14 traits, current {len(contracts)}")

    result: dict[str, TraitInfo] = {}
    for trait in sorted(contracts):
        contract = contracts[trait]
        meta = metadata[trait]
        analysis_type = clean(meta["analysis_type"])
        require(analysis_type in {"continuous", "binary"}, f"unknown analysis_type: {trait} -> {analysis_type}")
        require(
            clean(contract["expected_analysis_type"]) == analysis_type,
            f"trait type mismatch: {trait}",
        )
        require(clean(meta["locus_selection_p"]) == "p_lrt", f"{trait} not used p_lrt definition locus")
        require(clean(meta["fine_mapping_z"]) == "beta/se", f"{trait} not used beta/se definition z")
        input_file = Path(clean(contract["input_file"]))
        require(input_file.is_absolute(), f"GWAS input is not an absolute path: {input_file}")
        require(input_file.is_file() and input_file.stat().st_size > 0, f"GWAS input does not exist or is empty: {input_file}")
        n_main = parse_int(meta["n_main"], f"{trait}.n_main")
        expected_record_count = parse_int(
            contract["expected_record_count"], f"{trait}.expected_record_count"
        )
        require(n_main > 0 and expected_record_count > 0, f"{trait} of N or expected row count invalid")
        result[trait] = TraitInfo(
            trait=trait,
            trait_label=clean(meta["trait_label"]),
            analysis_type=analysis_type,
            n_main=n_main,
            n_sensitivity=clean(meta["n_sensitivity"]),
            input_file=input_file,
            expected_record_count=expected_record_count,
        )
    return result


def load_loci(path: Path, traits: Mapping[str, TraitInfo]) -> list[Locus]:
    _, rows = read_tsv(path, LOCUS_REQUIRED_FIELDS)
    loci: list[Locus] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        locus_id = clean(row["locus_id"])
        require(locus_id and locus_id not in seen, f"{path}:{row_number} locus_id duplicate or is empty: {locus_id}")
        require(re.fullmatch(r"[A-Za-z0-9_.-]+", locus_id) is not None, f"locus_id contains unsafe characters: {locus_id}")
        seen.add(locus_id)
        trait = clean(row["trait"])
        require(trait in traits, f"locus use unknown trait: {trait}")
        info = traits[trait]
        scope = clean(row["analysis_scope"])
        require(scope in {"primary", "exploratory"}, f"unknown analysis_scope: {locus_id} -> {scope}")
        require(clean(row["trait_label"]) == info.trait_label, f"locus trait_label mismatch: {locus_id}")
        require(clean(row["analysis_type"]) == info.analysis_type, f"locus analysis_type mismatch: {locus_id}")
        chrom = normalize_chrom(row["chrom"])
        core_start = parse_int(row["core_start"], f"{locus_id}.core_start")
        core_end = parse_int(row["core_end"], f"{locus_id}.core_end")
        analysis_start = parse_int(row["analysis_start"], f"{locus_id}.analysis_start")
        analysis_end = parse_int(row["analysis_end"], f"{locus_id}.analysis_end")
        lead_pos = parse_int(row["lead_pos"], f"{locus_id}.lead_pos")
        require(1 <= analysis_start <= core_start <= core_end <= analysis_end, f"locus boundary order invalid: {locus_id}")
        require(analysis_start <= lead_pos <= analysis_end, f"locus the analysis interval does not contain lead: {locus_id}")
        loci.append(
            Locus(
                locus_id=locus_id,
                trait=trait,
                trait_label=info.trait_label,
                analysis_type=info.analysis_type,
                analysis_scope=scope,
                chrom=chrom,
                core_start=core_start,
                core_end=core_end,
                analysis_start=analysis_start,
                analysis_end=analysis_end,
                region_source=clean(row["region_source"]),
                lead_variant_id=clean(row["lead_variant_id"]),
                lead_pos=lead_pos,
                lead_allele1=normalize_allele(row["lead_allele1"]),
                lead_allele0=normalize_allele(row["lead_allele0"]),
                lead_p_lrt=parse_float(row["lead_p_lrt"], f"{locus_id}.lead_p_lrt"),
                lead_p_lrt_tier=clean(row["lead_p_lrt_tier"]),
            )
        )
    require(loci, "no read to locus")
    scope_counts = Counter(locus.analysis_scope for locus in loci)
    require(scope_counts == Counter({"exploratory": 87, "primary": 36}), f"1.7 actual locus count does not match: {dict(scope_counts)}")
    return sorted(loci, key=lambda item: (item.trait, int(item.chrom), item.analysis_start, item.locus_id))


def load_bim(prefix: Path) -> tuple[dict[tuple[str, int], BimRecord], dict[str, BimRecord], int]:
    bed = Path(f"{prefix}.bed")
    bim = Path(f"{prefix}.bim")
    fam = Path(f"{prefix}.fam")
    for path in (bed, bim, fam):
        require(path.is_file() and path.stat().st_size > 0, f"LD prefix file does not exist or is empty: {path}")
    with fam.open("r", encoding="utf-8") as handle:
        sample_count = sum(1 for line in handle if line.strip())
    require(sample_count > 0, "final LD FAM no sample")

    by_position: dict[tuple[str, int], BimRecord] = {}
    by_id: dict[str, BimRecord] = {}
    with bim.open("r", encoding="utf-8") as handle:
        for order, line in enumerate(handle, start=1):
            fields = line.split()
            require(len(fields) == 6, f"{bim}:{order} does not have six columns BIM")
            chrom = normalize_chrom(fields[0])
            variant_id = fields[1]
            pos = parse_int(fields[3], f"{bim}:{order}.POS")
            a1 = normalize_allele(fields[4])
            a2 = normalize_allele(fields[5])
            require(a1 and a2 and a1 != a2 and a1 != "0" and a2 != "0", f"BIM allele invalid: {variant_id}")
            record = BimRecord(chrom, variant_id, pos, a1, a2, order)
            key = (chrom, pos)
            require(key not in by_position, f"final LD contains duplicate CHR+POS: {key}")
            require(variant_id not in by_id, f"final LD contains duplicate variant ID: {variant_id}")
            by_position[key] = record
            by_id[variant_id] = record

    expected_bed_size = 3 + math.ceil(sample_count / 4) * len(by_id)
    with bed.open("rb") as handle:
        magic = handle.read(3).hex()
    require(magic == "6c1b01", f"BED is not SNP-major format: magic={magic}")
    require(bed.stat().st_size == expected_bed_size, "BED size differs from the value expected from FAM/BIM")
    return by_position, by_id, sample_count


def load_frequency(path: Path, bim_by_id: Mapping[str, BimRecord]) -> dict[str, FrequencyRecord]:
    require(path.is_file() and path.stat().st_size > 0, f"LD frequency file does not exist or is empty: {path}")
    result: dict[str, FrequencyRecord] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"#CHROM", "POS", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"}
        missing = required - set(reader.fieldnames or [])
        require(not missing, f"LD frequency file missing fields: {','.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            variant_id = clean(row["ID"])
            require(variant_id in bim_by_id, f"LD frequency contains BIM outside variant: {variant_id}")
            require(variant_id not in result, f"LD frequency variant duplicate: {variant_id}")
            bim = bim_by_id[variant_id]
            record = FrequencyRecord(
                chrom=normalize_chrom(row["#CHROM"]),
                variant_id=variant_id,
                pos=parse_int(row["POS"], f"{path}:{row_number}.POS"),
                ref=normalize_allele(row["REF"]),
                alt=normalize_allele(row["ALT"]),
                alt_frequency=parse_float(row["ALT_FREQS"], f"{path}:{row_number}.ALT_FREQS"),
                observation_count=parse_int(row["OBS_CT"], f"{path}:{row_number}.OBS_CT"),
            )
            require(0.0 <= record.alt_frequency <= 1.0, f"LD frequency out of bounds: {variant_id}")
            require(record.observation_count > 0, f"LD frequency OBS_CT invalid: {variant_id}")
            require(
                (record.chrom, record.pos, record.alt, record.ref)
                == (bim.chrom, bim.pos, bim.a1, bim.a2),
                f"LD frequency REF/ALT differs from BIM A2/A1: {variant_id}",
            )
            result[variant_id] = record
    require(set(result) == set(bim_by_id), "Variant sets differ between LD frequencies and the final BIM")
    return result


def blank_audit_row(locus: Locus, trait: TraitInfo, row: Mapping[str, str], chrom: str, pos: int) -> dict[str, object]:
    return {
        "locus_id": locus.locus_id,
        "trait": locus.trait,
        "trait_label": locus.trait_label,
        "analysis_type": locus.analysis_type,
        "analysis_scope": locus.analysis_scope,
        "chrom": chrom,
        "pos": pos,
        "gwas_variant_id": clean(row.get("rs")),
        "gwas_allele1": normalize_allele(row.get("allele1")),
        "gwas_allele0": normalize_allele(row.get("allele0")),
        "gwas_af_for_allele1": clean(row.get("af")),
        "beta": clean(row.get("beta")),
        "se": clean(row.get("se")),
        "z_for_gwas_allele1": "",
        "p_wald": clean(row.get("p_wald")),
        "p_lrt": clean(row.get("p_lrt")),
        "p_score": clean(row.get("p_score")),
        "n_main": trait.n_main,
        "n_miss": clean(row.get("n_miss")),
        "n_obs": "",
        "missing_rate": "",
        "ld_variant_id": "",
        "ld_bim_order": "",
        "ld_a1": "",
        "ld_a2": "",
        "ld_a1_frequency": "",
        "gwas_frequency_aligned_to_ld_a1": "",
        "frequency_difference_label_orientation": "",
        "frequency_difference_swapped_orientation": "",
        "ambiguous_at_cg": "no",
        "z_for_ld_a1": "",
        "z_orientation": "",
        "harmonization_status": "",
        "main_analysis_include": "no",
        "strict_no_ambiguous_include": "no",
        "status_detail": "",
    }


def set_drop(audit: dict[str, object], status: str, detail: str) -> dict[str, object]:
    audit["harmonization_status"] = status
    audit["status_detail"] = detail
    return audit


def harmonize_row(
    locus: Locus,
    trait: TraitInfo,
    row: Mapping[str, str],
    chrom: str,
    pos: int,
    bim_by_position: Mapping[tuple[str, int], BimRecord],
    frequencies: Mapping[str, FrequencyRecord],
    max_missing_rate: float,
    ambiguous_tolerance: float,
    ambiguous_separation: float,
) -> dict[str, object]:
    audit = blank_audit_row(locus, trait, row, chrom, pos)
    a1 = normalize_allele(row.get("allele1"))
    a0 = normalize_allele(row.get("allele0"))
    if not a1 or not a0 or a1 == a0:
        return set_drop(audit, "INVALID_GWAS_ALLELES", "allele1/allele0 missing or identical.")

    bim = bim_by_position.get((chrom, pos))
    if bim is None:
        return set_drop(audit, "NO_LD_POSITION", "No matching CHR+POS in the final LD panel.")
    audit.update(
        {
            "ld_variant_id": bim.variant_id,
            "ld_bim_order": bim.order,
            "ld_a1": bim.a1,
            "ld_a2": bim.a2,
        }
    )
    if allele_pair(a1, a0) != allele_pair(bim.a1, bim.a2):
        return set_drop(audit, "ALLELE_SET_MISMATCH", "Matching position but different unordered allele sets.")

    beta = maybe_float(row.get("beta"))
    se = maybe_float(row.get("se"))
    if beta is None or se is None or se <= 0:
        return set_drop(audit, "INVALID_EFFECT_STATISTICS", "beta or se missing/invalid, or se<=0.")
    p_values = [maybe_float(row.get(field)) for field in ("p_wald", "p_lrt", "p_score")]
    if any(value is None or value < 0 or value > 1 for value in p_values):
        return set_drop(audit, "INVALID_P_VALUE", "p_wald/p_lrt/p_score missing or out of bounds.")
    af = maybe_float(row.get("af"))
    if af is None or af < 0 or af > 1:
        return set_drop(audit, "INVALID_GWAS_AF", "GWAS allele1 frequency missing or out of bounds.")
    n_miss = maybe_int(row.get("n_miss"))
    if n_miss is None or n_miss < 0 or n_miss >= trait.n_main:
        return set_drop(audit, "INVALID_N_MISS", "n_miss missing, negative, or not smaller than the trait N.")

    n_obs = trait.n_main - n_miss
    missing_rate = n_miss / trait.n_main
    z_gwas = beta / se
    audit.update(
        {
            "z_for_gwas_allele1": format_float(z_gwas),
            "n_miss": n_miss,
            "n_obs": n_obs,
            "missing_rate": format_float(missing_rate),
        }
    )
    if missing_rate > max_missing_rate:
        return set_drop(
            audit,
            "VARIANT_MISSING_RATE_GT_THRESHOLD",
            f"n_miss/n_main={missing_rate:.6g}>{max_missing_rate:.6g}",
        )

    frequency = frequencies[bim.variant_id]
    ld_a1_frequency = frequency.alt_frequency
    if a1 == bim.a1:
        z_ld = z_gwas
        aligned_af = af
        orientation = "GWAS_ALLELE1_EQUALS_LD_A1"
    elif a1 == bim.a2:
        z_ld = -z_gwas
        aligned_af = 1.0 - af
        orientation = "GWAS_ALLELE1_EQUALS_LD_A2_Z_FLIPPED_TO_LD_A1"
    else:
        return set_drop(audit, "EFFECT_ALLELE_NOT_IN_LD_PAIR", "GWAS allele1 is not in the LD allele pair.")

    label_difference = abs(aligned_af - ld_a1_frequency)
    swapped_difference = abs((1.0 - aligned_af) - ld_a1_frequency)
    ambiguous = is_palindromic_snp(a1, a0)
    audit.update(
        {
            "ld_a1_frequency": format_float(ld_a1_frequency),
            "gwas_frequency_aligned_to_ld_a1": format_float(aligned_af),
            "frequency_difference_label_orientation": format_float(label_difference),
            "frequency_difference_swapped_orientation": format_float(swapped_difference),
            "ambiguous_at_cg": "yes" if ambiguous else "no",
            "z_for_ld_a1": format_float(z_ld),
            "z_orientation": orientation,
        }
    )

    if ambiguous:
        if label_difference > ambiguous_tolerance:
            return set_drop(
                audit,
                "AMBIGUOUS_FREQUENCY_LABEL_ORIENTATION_MISMATCH",
                f"Label-orientation frequency difference {label_difference:.6g}>{ambiguous_tolerance:.6g}",
            )
        if swapped_difference <= ambiguous_tolerance:
            return set_drop(
                audit,
                "AMBIGUOUS_FREQUENCY_BOTH_ORIENTATIONS_PLAUSIBLE",
                "Both label and swapped orientations are within the frequency tolerance.",
            )
        if swapped_difference - label_difference < ambiguous_separation:
            return set_drop(
                audit,
                "AMBIGUOUS_FREQUENCY_NOT_SEPARATED",
                f"Separation between the two orientation-specific frequency differences is less than {ambiguous_separation:.6g}",
            )
        status = "HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED"
        detail = "For this A/T or C/G site, GWAS allele1 and LD A1 frequencies uniquely support the label orientation."
    else:
        status = "HARMONIZED_NONAMBIGUOUS"
        detail = (
            f"Non-ambiguous site with frequency difference >{ambiguous_tolerance:.6g}; retained and flagged for LD mismatch diagnostics."
            if label_difference > ambiguous_tolerance
            else "Coordinates, allele set, effect direction, and sample missingness passed checks."
        )

    audit.update(
        {
            "harmonization_status": status,
            "main_analysis_include": "yes",
            "strict_no_ambiguous_include": "no" if ambiguous else "yes",
            "status_detail": detail,
        }
    )
    return audit


def locus_output_dir(output_dir: Path, locus: Locus) -> Path:
    return output_dir / locus.analysis_scope / locus.trait / locus.locus_id


def write_locus_outputs(
    output_dir: Path,
    locus: Locus,
    audits: Sequence[Mapping[str, object]],
    min_variants: int,
    frequency_warning_threshold: float,
) -> tuple[dict[str, object], dict[str, object]]:
    locus_dir = locus_output_dir(output_dir, locus)
    locus_dir.mkdir(parents=True, exist_ok=False)
    harmonized = [dict(row) for row in audits if row["main_analysis_include"] == "yes"]
    harmonized.sort(key=lambda row: (int(row["ld_bim_order"]), int(row["pos"])))
    for index, row in enumerate(harmonized, start=1):
        row["variant_index"] = index
    dropped = [row for row in audits if row["main_analysis_include"] != "yes"]

    harmonized_path = locus_dir / "harmonized_variants.tsv.gz"
    order_path = locus_dir / "variant_order.tsv"
    audit_path = locus_dir / "allele_harmonization_audit.tsv.gz"
    dropped_path = locus_dir / "dropped_variants.tsv.gz"
    metadata_path = locus_dir / "locus_metadata.tsv"
    write_gzip_tsv(harmonized_path, HARMONIZED_FIELDS, harmonized)
    write_tsv(order_path, ORDER_FIELDS, harmonized)
    write_gzip_tsv(audit_path, AUDIT_FIELDS, audits)
    write_gzip_tsv(dropped_path, AUDIT_FIELDS, dropped)

    status_counts = Counter(str(row["harmonization_status"]) for row in audits)
    n_obs = [int(row["n_obs"]) for row in harmonized]
    n_locus = statistics.median(n_obs) if n_obs else None
    lead_included = any(str(row["ld_variant_id"]) == locus.lead_variant_id for row in harmonized)
    if not harmonized:
        locus_status = "NO_HARMONIZED_VARIANTS"
    elif len(harmonized) < min_variants:
        locus_status = "SINGLE_MARKER_NOT_RESOLUTION_ASSESSABLE"
    elif not lead_included:
        locus_status = "READY_WITHOUT_TRIGGER_LEAD"
    else:
        locus_status = "READY_FOR_GENOTYPE_MATRIX"

    metadata_rows = [
        ("locus_id", locus.locus_id),
        ("trait", locus.trait),
        ("trait_label", locus.trait_label),
        ("analysis_type", locus.analysis_type),
        ("analysis_scope", locus.analysis_scope),
        ("chrom", locus.chrom),
        ("analysis_start", locus.analysis_start),
        ("analysis_end", locus.analysis_end),
        ("region_source", locus.region_source),
        ("lead_variant_id", locus.lead_variant_id),
        ("lead_pos", locus.lead_pos),
        ("lead_p_lrt", format_float(locus.lead_p_lrt)),
        ("lead_p_lrt_tier", locus.lead_p_lrt_tier),
        ("harmonized_variant_count", len(harmonized)),
        ("strict_non_ambiguous_variant_count", sum(row["strict_no_ambiguous_include"] == "yes" for row in harmonized)),
        ("n_locus_median", format_float(float(n_locus)) if n_locus is not None else ""),
        ("lead_in_harmonized_variants", "yes" if lead_included else "no"),
        ("locus_input_status", locus_status),
    ]
    write_tsv(metadata_path, ("key", "value"), ({"key": key, "value": value} for key, value in metadata_rows))

    summary = {
        "locus_id": locus.locus_id,
        "trait": locus.trait,
        "trait_label": locus.trait_label,
        "analysis_type": locus.analysis_type,
        "analysis_scope": locus.analysis_scope,
        "chrom": locus.chrom,
        "analysis_start": locus.analysis_start,
        "analysis_end": locus.analysis_end,
        "region_source": locus.region_source,
        "lead_variant_id": locus.lead_variant_id,
        "lead_pos": locus.lead_pos,
        "lead_p_lrt": format_float(locus.lead_p_lrt),
        "lead_p_lrt_tier": locus.lead_p_lrt_tier,
        "gwas_rows_in_region": len(audits),
        "harmonized_variant_count": len(harmonized),
        "strict_non_ambiguous_variant_count": sum(row["strict_no_ambiguous_include"] == "yes" for row in harmonized),
        "dropped_variant_count": len(dropped),
        "no_ld_position_count": status_counts["NO_LD_POSITION"],
        "allele_set_mismatch_count": status_counts["ALLELE_SET_MISMATCH"],
        "missing_rate_filtered_count": status_counts["VARIANT_MISSING_RATE_GT_THRESHOLD"],
        "ambiguous_frequency_resolved_count": status_counts["HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED"],
        "ambiguous_frequency_dropped_count": sum(
            count for status, count in status_counts.items() if status.startswith("AMBIGUOUS_FREQUENCY_")
        ),
        "nonambiguous_frequency_warning_count": sum(
            row["harmonization_status"] == "HARMONIZED_NONAMBIGUOUS"
            and maybe_float(row["frequency_difference_label_orientation"]) is not None
            and float(row["frequency_difference_label_orientation"])
            > frequency_warning_threshold
            for row in harmonized
        ),
        "z_flipped_count": sum("Z_FLIPPED" in str(row["z_orientation"]) for row in harmonized),
        "n_locus_median": format_float(float(n_locus)) if n_locus is not None else "",
        "n_obs_min": min(n_obs) if n_obs else "",
        "n_obs_max": max(n_obs) if n_obs else "",
        "lead_in_harmonized_variants": "yes" if lead_included else "no",
        "locus_input_status": locus_status,
    }
    manifest = {
        "locus_id": locus.locus_id,
        "trait": locus.trait,
        "analysis_scope": locus.analysis_scope,
        "locus_input_status": locus_status,
        "harmonized_variant_count": len(harmonized),
        "n_locus_median": summary["n_locus_median"],
        "locus_directory": locus_dir,
        "harmonized_variants": harmonized_path,
        "variant_order": order_path,
        "allele_harmonization_audit": audit_path,
        "dropped_variants": dropped_path,
        "locus_metadata": metadata_path,
    }
    return summary, manifest


def process_all_traits(
    args: argparse.Namespace,
    traits: Mapping[str, TraitInfo],
    loci: Sequence[Locus],
    bim_by_position: Mapping[tuple[str, int], BimRecord],
    frequencies: Mapping[str, FrequencyRecord],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    loci_by_trait_chrom: dict[tuple[str, str], list[Locus]] = defaultdict(list)
    audits_by_locus: dict[str, list[dict[str, object]]] = {locus.locus_id: [] for locus in loci}
    for locus in loci:
        loci_by_trait_chrom[(locus.trait, locus.chrom)].append(locus)
    for key in loci_by_trait_chrom:
        loci_by_trait_chrom[key].sort(key=lambda locus: (locus.analysis_start, locus.analysis_end, locus.locus_id))

    for trait_index, trait_name in enumerate(sorted({locus.trait for locus in loci}), start=1):
        trait = traits[trait_name]
        log(f"scan GWAS summary {trait_index}/{len(traits)}: {trait_name}")
        observed_rows = 0
        relevant_keys: set[tuple[str, int, tuple[str, str]]] = set()
        with gzip.open(trait.input_file, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [field for field in GWAS_REQUIRED_FIELDS if field not in (reader.fieldnames or [])]
            require(not missing, f"{trait.input_file} missing fields: {','.join(missing)}")
            for row_number, row in enumerate(reader, start=2):
                observed_rows += 1
                if args.progress_every > 0 and observed_rows % args.progress_every == 0:
                    log(f"{trait_name} scanned {observed_rows} row")
                try:
                    chrom = normalize_chrom(row.get("chr"))
                except RuntimeError as exc:
                    raise RuntimeError(f"{trait.input_file}:{row_number} {exc}") from exc
                pos = maybe_int(row.get("ps"))
                require(pos is not None and pos > 0, f"{trait.input_file}:{row_number} ps invalid")
                matching_loci = [
                    locus
                    for locus in loci_by_trait_chrom.get((trait_name, chrom), ())
                    if locus.analysis_start <= pos <= locus.analysis_end
                ]
                if not matching_loci:
                    continue
                a1 = normalize_allele(row.get("allele1"))
                a0 = normalize_allele(row.get("allele0"))
                key = (chrom, pos, allele_pair(a1, a0))
                require(key not in relevant_keys, f"{trait.input_file}:{row_number} locus within region GWAS duplicate key: {key}")
                relevant_keys.add(key)
                for locus in matching_loci:
                    audits_by_locus[locus.locus_id].append(
                        harmonize_row(
                            locus,
                            trait,
                            row,
                            chrom,
                            pos,
                            bim_by_position,
                            frequencies,
                            args.max_missing_rate,
                            args.ambiguous_frequency_tolerance,
                            args.ambiguous_frequency_separation,
                        )
                    )
        require(
            observed_rows == trait.expected_record_count,
            f"{trait_name}: actual row count {observed_rows} differs from the Step 1.1 expected count {trait.expected_record_count}",
        )
        log(f"{trait_name} scan completed: total row count={observed_rows}; within region unique GWAS key={len(relevant_keys)}")

    summaries: list[dict[str, object]] = []
    manifests: list[dict[str, object]] = []
    for locus_index, locus in enumerate(loci, start=1):
        audits = sorted(
            audits_by_locus[locus.locus_id],
            key=lambda row: (int(row["pos"]), str(row["gwas_variant_id"])),
        )
        summary, manifest = write_locus_outputs(
            args.output_dir,
            locus,
            audits,
            args.min_variants_for_finemapping,
            args.ambiguous_frequency_tolerance,
        )
        summaries.append(summary)
        manifests.append(manifest)
        if locus_index % 20 == 0 or locus_index == len(loci):
            log(f"per locus output progress: {locus_index}/{len(loci)}")
    return summaries, manifests


def carry_forward_non_finemappable(
    input_path: Path,
    output_path: Path,
    locus_ids: set[str],
) -> tuple[int, int, int]:
    header, rows = read_tsv(input_path, NON_FINEMAPPABLE_REQUIRED_FIELDS)
    require(len(rows) == 74, f"1.7 cannot direct fine-mapping candidate expected 74 row, current {len(rows)}")
    seen: set[str] = set()
    assigned = 0
    unresolved = 0
    output_rows: list[dict[str, object]] = []
    for row in rows:
        candidate_id = clean(row["candidate_row_id"])
        require(candidate_id and candidate_id not in seen, f"cannot fine-mapping candidate_row_id duplicate: {candidate_id}")
        seen.add(candidate_id)
        require(clean(row["candidate_pip_eligible"]) == "no", f"candidate is error label as PIP eligible: {candidate_id}")
        assigned_locus = clean(row["assigned_locus_id"])
        if assigned_locus:
            require(assigned_locus in locus_ids, f"cannot fine-mapping candidate points to unknown locus: {candidate_id}")
            assigned += 1
            disposition = "LOCUS_CONTEXT_ONLY_NO_VARIANT_PIP"
        else:
            unresolved += 1
            disposition = "OUTSIDE_MODELED_LOCI_AUDIT_ONLY"
        output = dict(row)
        output["step18_disposition"] = disposition
        output_rows.append(output)
    require((assigned, unresolved) == (31, 43), f"74 missing-candidate strata and 1.7 result does not match: assigned={assigned}; unresolved={unresolved}")
    write_tsv(output_path, (*header, "step18_disposition"), output_rows)
    return len(rows), assigned, unresolved


def write_root_outputs(
    args: argparse.Namespace,
    loci: Sequence[Locus],
    summaries: Sequence[Mapping[str, object]],
    manifests: Sequence[Mapping[str, object]],
    traits: Mapping[str, TraitInfo],
    sample_count: int,
    non_finemappable_counts: tuple[int, int, int],
    upstream_statuses: Mapping[str, str],
) -> None:
    write_tsv(args.output_dir / "harmonization_summary.tsv", SUMMARY_FIELDS, summaries)
    write_tsv(args.output_dir / "locus_input_manifest.tsv", MANIFEST_FIELDS, manifests)

    status_counts = Counter(str(row["locus_input_status"]) for row in summaries)
    scope_counts = Counter(str(row["analysis_scope"]) for row in summaries)
    total_harmonized = sum(int(row["harmonized_variant_count"]) for row in summaries)
    total_strict = sum(int(row["strict_non_ambiguous_variant_count"]) for row in summaries)
    total_dropped = sum(int(row["dropped_variant_count"]) for row in summaries)
    total_ambiguous_resolved = sum(int(row["ambiguous_frequency_resolved_count"]) for row in summaries)
    total_ambiguous_dropped = sum(int(row["ambiguous_frequency_dropped_count"]) for row in summaries)
    warning_loci = sum(status != "READY_FOR_GENOTYPE_MATRIX" for status in (row["locus_input_status"] for row in summaries))
    overall_status = "PASS_WITH_WARNINGS" if warning_loci or non_finemappable_counts[0] else "PASS"
    require(total_harmonized > 0, "all locus all no available harmonize variant")

    overall_rows: list[dict[str, object]] = [
        {"metric": "overall_status", "value": overall_status, "status": "PASS", "details": "Per-locus GWAS-LD allele harmonization completed."},
        {"metric": "trait_count", "value": len(traits), "status": "PASS", "details": "Traits with complete GWAS summary statistics."},
        {"metric": "locus_count", "value": len(loci), "status": "PASS", "details": "1.7 all_locus_regions"},
        {"metric": "primary_locus_count", "value": scope_counts["primary"], "status": "PASS", "details": "Primary interpretation branch."},
        {"metric": "exploratory_locus_count", "value": scope_counts["exploratory"], "status": "PASS", "details": "Exploratory interpretation branch."},
        {"metric": "ld_reference_samples", "value": sample_count, "status": "PASS", "details": "Independent samples in the final GRLS+PILOT panel."},
        {"metric": "harmonized_locus_variant_rows", "value": total_harmonized, "status": "PASS", "details": "The same variant may enter multiple loci through overlapping windows."},
        {"metric": "strict_non_ambiguous_locus_variant_rows", "value": total_strict, "status": "INFO", "details": "Subsequent sensitivity branch excluding all A/T and C/G variants."},
        {"metric": "dropped_locus_variant_rows", "value": total_dropped, "status": "INFO", "details": "See the per-locus dropped_variants files."},
        {"metric": "ambiguous_frequency_resolved_rows", "value": total_ambiguous_resolved, "status": "INFO", "details": "A/T or C/G variants for which frequencies uniquely support the label orientation."},
        {"metric": "ambiguous_frequency_dropped_rows", "value": total_ambiguous_dropped, "status": "WARN" if total_ambiguous_dropped else "PASS", "details": "A/T or C/G variants unresolved by frequency information."},
        {"metric": "non_finemappable_candidate_rows", "value": non_finemappable_counts[0], "status": "WARN", "details": "Will not receive a main-LD PIP."},
        {"metric": "non_finemappable_within_locus_rows", "value": non_finemappable_counts[1], "status": "INFO", "details": "Retained only as regional context."},
        {"metric": "non_finemappable_outside_locus_rows", "value": non_finemappable_counts[2], "status": "INFO", "details": "Audit records outside regions."},
    ]
    for status, count in sorted(status_counts.items()):
        overall_rows.append(
            {
                "metric": f"locus_input_status::{status}",
                "value": count,
                "status": "PASS" if status == "READY_FOR_GENOTYPE_MATRIX" else "WARN",
                "details": "Number of loci.",
            }
        )
    write_tsv(
        args.output_dir / "harmonization_overall_summary.tsv",
        ("metric", "value", "status", "details"),
        overall_rows,
    )

    lineage_rows = [
        {"role": "locus_regions", "input_path": args.locus_regions, "producer_script": args.step17_producer},
        {"role": "non_finemappable_candidates", "input_path": args.non_finemappable_candidates, "producer_script": args.step17_producer},
        {"role": "trait_input_contract", "input_path": args.trait_input_contract, "producer_script": args.step11_producer},
        {"role": "trait_metadata", "input_path": args.trait_metadata, "producer_script": args.step12_producer},
        {"role": "final_ld_bed", "input_path": Path(f"{args.ld_prefix}.bed"), "producer_script": args.step15_producer},
        {"role": "final_ld_bim", "input_path": Path(f"{args.ld_prefix}.bim"), "producer_script": args.step15_producer},
        {"role": "final_ld_fam", "input_path": Path(f"{args.ld_prefix}.fam"), "producer_script": args.step15_producer},
        {"role": "final_ld_frequency", "input_path": args.ld_frequency, "producer_script": args.step15_producer},
    ]
    for trait in sorted(traits):
        lineage_rows.append(
            {"role": f"gwas_summary::{trait}", "input_path": traits[trait].input_file, "producer_script": "external_PNAS_GWAS_summary_statistics"}
        )
    write_tsv(
        args.output_dir / "input_lineage.tsv",
        ("role", "input_path", "producer_script"),
        lineage_rows,
    )
    parameter_rows = [
        ("max_missing_rate", args.max_missing_rate),
        ("ambiguous_frequency_tolerance", args.ambiguous_frequency_tolerance),
        ("ambiguous_frequency_separation", args.ambiguous_frequency_separation),
        ("min_variants_for_finemapping", args.min_variants_for_finemapping),
        ("variant_match_rule", "same_CHR_POS_and_same_unordered_allele_pair_no_complement_guess"),
        ("z_definition", "beta_over_se_oriented_to_final_BIM_A1"),
        ("primary_locus_count", scope_counts["primary"]),
        ("exploratory_locus_count", scope_counts["exploratory"]),
        ("step15_status", upstream_statuses["step15"]),
        ("step17_status", upstream_statuses["step17"]),
    ]
    write_tsv(
        args.output_dir / "pipeline_parameters.tsv",
        ("parameter", "value"),
        ({"parameter": key, "value": value} for key, value in parameter_rows),
    )
    with (args.output_dir / "COMPLETED.ok").open("w", encoding="utf-8") as handle:
        handle.write(f"overall_status\t{overall_status}\n")
        handle.write(f"trait_count\t{len(traits)}\n")
        handle.write(f"locus_count\t{len(loci)}\n")
        handle.write(f"harmonized_locus_variant_rows\t{total_harmonized}\n")
        handle.write(f"output_dir\t{args.output_dir}\n")


def validate_args(args: argparse.Namespace) -> None:
    require(args.output_dir.is_absolute(), f"output directory must be an absolute path: {args.output_dir}")
    require(not args.output_dir.exists(), f"refuse to overwrite existing output directory: {args.output_dir}")
    require(0 <= args.max_missing_rate < 1, "max_missing_rate must be in[0,1)within")
    require(0 <= args.ambiguous_frequency_tolerance <= 0.5, "ambiguous_frequency_tolerance must be in[0,0.5]within")
    require(0 <= args.ambiguous_frequency_separation <= 1, "ambiguous_frequency_separation must be in[0,1]within")
    require(args.min_variants_for_finemapping >= 2, "min_variants_for_finemapping at least 2")
    for path in (
        args.locus_regions,
        args.non_finemappable_candidates,
        args.step17_completion,
        args.trait_input_contract,
        args.trait_metadata,
        args.ld_frequency,
        args.step15_completion,
        args.step11_producer,
        args.step12_producer,
        args.step15_producer,
        args.step17_producer,
    ):
        require(path.is_absolute(), f"input must be an absolute path: {path}")
        require(path.is_file() and path.stat().st_size > 0, f"input does not exist or is empty: {path}")


def main() -> int:
    args = parse_args()
    validate_args(args)
    log("stage 1/6: check upstream completion status, trait contract and locus region")
    upstream_statuses = {
        "step15": read_completion_status(
            args.step15_completion,
            {"COMPLETED", "PASS", "PASS_WITH_WARNINGS"},
            fallback_completion_key="completed_at",
        ),
        "step17": read_completion_status(args.step17_completion, {"PASS", "PASS_WITH_WARNINGS"}),
    }
    traits = load_traits(args.trait_input_contract, args.trait_metadata)
    loci = load_loci(args.locus_regions, traits)

    log("stage 2/6: check final LD BIM/FAM/BED and A1 frequency")
    bim_by_position, bim_by_id, sample_count = load_bim(args.ld_prefix)
    frequencies = load_frequency(args.ld_frequency, bim_by_id)
    log(f"final LD check completed: sample={sample_count}; variant={len(bim_by_id)}")

    log("stage 3/6: create this run independent output directory")
    args.output_dir.mkdir(parents=True, exist_ok=False)

    log("stage 4/6: scan 14 complete GWAS summary and per locus harmonize")
    summaries, manifests = process_all_traits(
        args, traits, loci, bim_by_position, frequencies
    )

    log("stage 5/6: continue 74 items cannot direct obtain PIP candidate of audit status")
    non_finemappable_counts = carry_forward_non_finemappable(
        args.non_finemappable_candidates,
        args.output_dir / "non_finemappable_candidate_carry_forward.tsv",
        {locus.locus_id for locus in loci},
    )

    log("stage 6/6: write overview, manifest, parameter, input provenance and completion marker")
    write_root_outputs(
        args,
        loci,
        summaries,
        manifests,
        traits,
        sample_count,
        non_finemappable_counts,
        upstream_statuses,
    )
    log(f"output directory: {args.output_dir}")
    log("task completed / done")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"error: {exc}")
        raise
