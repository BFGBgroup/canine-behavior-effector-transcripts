#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import shlex
import subprocess
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_GROUP_COUNTS = {"GRLS": 2846, "PILOT": 178, "Oldies": 200}
EXPECTED_TOTAL_SAMPLES = 3224
EXPECTED_TOTAL_VARIANTS = 913984
EXPECTED_AUTOSOME_VARIANTS = 884031
AUTOSOME_RE = re.compile(r"^chr([1-9]|[12][0-9]|3[0-8])$")
MISSING_VALUES = {".", "NA", "nan", "-nan"}


class PipelineError(RuntimeError):
    pass


def log(message: str) -> None:
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S} ] {message}", flush=True)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PipelineError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run 1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.py; see the workflow documentation.')
    parser.add_argument("--input-prefix", required=True, type=Path)
    parser.add_argument("--step11-dir", required=True, type=Path)
    parser.add_argument("--step12-dir", required=True, type=Path)
    parser.add_argument("--plink2-bin", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--mind", type=float, default=0.05)
    parser.add_argument("--geno", type=float, default=0.05)
    parser.add_argument("--maf", type=float, default=0.05)
    parser.add_argument("--hwe", type=float, default=1e-6)
    parser.add_argument("--king-relatedness", type=float, default=0.177)
    parser.add_argument("--prune-window", default="500kb")
    parser.add_argument("--prune-step", type=int, default=1)
    parser.add_argument("--prune-r2", type=float, default=0.2)
    parser.add_argument("--max-final-qc-iterations", type=int, default=10)
    return parser.parse_args()


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file(), f"missing TSV file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"TSV no column header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def read_whitespace_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file(), f"missing PLINK 2 result: {path}")
    with path.open("r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    require(lines, f"PLINK 2 result is empty: {path}")
    raw_header = lines[0].split()
    header = [field.lstrip("#") for field in raw_header]
    require(len(header) == len(set(header)), f"duplicate column headers after normalization: {path}: {header}")
    rows: list[dict[str, str]] = []
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split()
        require(
            len(fields) == len(header),
            f"{path}:{line_number} column count {len(fields)}, expected {len(header)}",
        )
        rows.append(dict(zip(header, fields)))
    return header, rows


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_id_pairs(path: Path, ids: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for fid, iid in sorted(set(ids)):
            handle.write(f"{fid}\t{iid}\n")


def write_variant_ids(path: Path, ids: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for variant_id in sorted(set(ids)):
            handle.write(f"{variant_id}\n")


def finite_float(value: str, context: str) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise PipelineError(f"{context} is not numeric: {value}") from exc
    require(math.isfinite(number), f"{context} is not finite numeric value: {value}")
    return number


def run_command(command: Sequence[object], label: str) -> None:
    command_text = " ".join(shlex.quote(str(part)) for part in command)
    log(f"RUN {label}: {command_text}")
    completed = subprocess.run([str(part) for part in command], check=False)
    require(completed.returncode == 0, f"{label} failed, exit code {completed.returncode}")


def read_plink2_version(plink2_bin: Path) -> str:
    completed = subprocess.run(
        [str(plink2_bin), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    require(completed.returncode == 0, "cannot read PLINK 2 version")
    version = (completed.stdout or completed.stderr).strip()
    require(version.startswith("PLINK v2."), f"is not PLINK 2: {version}")
    return version


def plink2_base(
    args: argparse.Namespace, input_prefix: Path, output_prefix: Path
) -> list[object]:
    return [
        args.plink2_bin,
        "--bfile",
        input_prefix,
        "--dog",
        "--allow-extra-chr",
        "--threads",
        args.threads,
        "--out",
        output_prefix,
    ]


def plink2_pfile_base(
    args: argparse.Namespace, input_prefix: Path, output_prefix: Path
) -> list[object]:
    return [
        args.plink2_bin,
        "--pfile",
        input_prefix,
        "--dog",
        "--allow-extra-chr",
        "--threads",
        args.threads,
        "--out",
        output_prefix,
    ]


def make_bed(
    args: argparse.Namespace,
    input_prefix: Path,
    output_prefix: Path,
    *,
    keep: Path | None = None,
    remove: Path | None = None,
    extract: Path | None = None,
    exclude: Path | None = None,
    update_chr: Path | None = None,
) -> None:
    command_output_prefix = (
        Path(f"{output_prefix}.sorted_pgen")
        if update_chr is not None
        else output_prefix
    )
    command = plink2_base(args, input_prefix, command_output_prefix)
    if keep is not None:
        command.extend(["--keep", keep])
    if remove is not None and remove.stat().st_size > 0:
        command.extend(["--remove", remove])
    if extract is not None:
        command.extend(["--extract", extract])
    if exclude is not None and exclude.stat().st_size > 0:
        command.extend(["--exclude", exclude])
    if update_chr is not None:

        sorted_pgen_prefix = command_output_prefix
        command.extend(
            ["--update-chr", update_chr, "1", "2", "--sort-vars", "--make-pgen"]
        )
        run_command(command, f"PLINK 2 sorted PGEN {output_prefix.name}")
        for suffix in (".pgen", ".pvar", ".psam"):
            path = Path(f"{sorted_pgen_prefix}{suffix}")
            require(
                path.is_file() and path.stat().st_size > 0,
                f"not generated valid of sort PGEN file: {path}",
            )
        conversion_command = plink2_pfile_base(
            args, sorted_pgen_prefix, output_prefix
        )
        conversion_command.append("--make-bed")
        run_command(
            conversion_command,
            f"PLINK 2 sorted PGEN to BED {output_prefix.name}",
        )
    else:
        command.append("--make-bed")
        run_command(command, f"PLINK 2 make-bed {output_prefix.name}")
    for suffix in (".bed", ".bim", ".fam"):
        path = Path(f"{output_prefix}{suffix}")
        require(path.is_file() and path.stat().st_size > 0, f"not generated valid file: {path}")


def run_missing(args: argparse.Namespace, input_prefix: Path, output_prefix: Path) -> None:
    command = plink2_base(args, input_prefix, output_prefix)
    command.extend(
        [
            "--missing",
            "scols=fid,nmiss,nobs,fmiss",
            "vcols=chrom,pos,ref,alt,provref,nmiss,nobs,fmiss",
        ]
    )
    run_command(command, f"PLINK 2 missingness {output_prefix.name}")


def run_freq(args: argparse.Namespace, input_prefix: Path, output_prefix: Path) -> None:
    command = plink2_base(args, input_prefix, output_prefix)
    command.extend(["--freq", "cols=chrom,pos,ref,alt,provref,altfreq,nobs"])
    run_command(command, f"PLINK 2 frequency {output_prefix.name}")


def run_hardy(args: argparse.Namespace, input_prefix: Path, output_prefix: Path) -> None:
    command = plink2_base(args, input_prefix, output_prefix)
    command.extend(["--hardy", "cols=chrom,pos,ref,alt,provref,ax,gcounts,p"])
    run_command(command, f"PLINK 2 Hardy-Weinberg {output_prefix.name}")


def prefix_counts(prefix: Path) -> tuple[int, int]:
    fam = Path(f"{prefix}.fam")
    bim = Path(f"{prefix}.bim")
    require(fam.is_file() and bim.is_file(), f"PLINK prefix incomplete: {prefix}")
    with fam.open("r", encoding="utf-8") as handle:
        sample_count = sum(1 for line in handle if line.strip())
    with bim.open("r", encoding="utf-8") as handle:
        variant_count = sum(1 for line in handle if line.strip())
    return sample_count, variant_count


def read_fam_ids(prefix: Path) -> set[tuple[str, str]]:
    fam = Path(f"{prefix}.fam")
    require(fam.is_file(), f"missing FAM: {fam}")
    samples: set[tuple[str, str]] = set()
    with fam.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            require(len(fields) == 6, f"{fam}:{line_number} does not have six columns FAM")
            sample = (fields[0], fields[1])
            require(sample not in samples, f"{fam}:{line_number} sample duplicate: {sample}")
            samples.add(sample)
    return samples


def read_bim_ids(prefix: Path) -> set[str]:
    bim = Path(f"{prefix}.bim")
    require(bim.is_file(), f"missing BIM: {bim}")
    variants: set[str] = set()
    with bim.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            require(len(fields) == 6, f"{bim}:{line_number} does not have six columns BIM")
            variant_id = fields[1]
            require(variant_id not in variants, f"{bim}:{line_number} variant ID duplicate: {variant_id}")
            variants.add(variant_id)
    return variants


def read_id_pairs(path: Path, *, allow_empty: bool = False) -> set[tuple[str, str]]:
    require(path.is_file(), f"missing sample ID file: {path}")
    samples: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            require(len(fields) == 2, f"{path}:{line_number} must be FID/IID two columns")
            sample = (fields[0], fields[1])
            require(sample not in samples, f"{path}:{line_number} sample duplicate: {sample}")
            samples.add(sample)
    require(allow_empty or samples, f"sample ID file is empty: {path}")
    return samples


def validate_prerequisites(args: argparse.Namespace) -> str:
    require(args.threads >= 1, "--threads must be greater than 0")
    require(0.0 <= args.mind < 1.0, "--mind must be located in [0, 1)")
    require(0.0 <= args.geno < 1.0, "--geno must be located in [0, 1)")
    require(0.0 < args.maf <= 0.5, "--maf must be located in (0, 0.5]")
    require(0.0 < args.hwe <= 1.0, "--hwe must be located in (0, 1]")
    require(0.0 < args.king_relatedness < 0.5, "--king-relatedness must be located in (0, 0.5)")
    require(args.prune_step >= 1, "--prune-step must be greater than 0")
    require(0.0 < args.prune_r2 < 1.0, "--prune-r2 must be located in (0, 1)")
    require(args.max_final_qc_iterations >= 1, "--max-final-qc-iterations must be greater than 0")
    if str(args.prune_window).strip().lower().endswith("kb"):
        require(args.prune_step == 1, "kb window by PLINK 2 requirement can only use step=1")
    require(args.plink2_bin.is_file(), f"PLINK 2 does not exist: {args.plink2_bin}")
    for suffix in (".bed", ".bim", ".fam"):
        path = Path(f"{args.input_prefix}{suffix}")
        require(path.is_file() and path.stat().st_size > 0, f"input missing or is empty: {path}")

    step11_summary = args.step11_dir / "validation_summary.tsv"
    step11_complete = args.step11_dir / "COMPLETED.ok"
    step12_summary = args.step12_dir / "manifest_summary.tsv"
    step12_complete = args.step12_dir / "COMPLETED.ok"
    for path in (step11_summary, step11_complete, step12_summary, step12_complete):
        require(path.is_file() and path.stat().st_size > 0, f"upstream not completed or file is empty: {path}")

    _, step11_rows = read_tsv(step11_summary)
    step11_values = {row.get("metric", ""): row.get("value", "") for row in step11_rows}
    require(
        step11_values.get("overall_status") in {"PASS", "PASS_WITH_WARNINGS"},
        f"Step 1.1 status unexpected: {step11_values.get('overall_status')}",
    )
    require(step11_values.get("ld_samples") == str(EXPECTED_TOTAL_SAMPLES), "Step 1.1 sample count change")
    require(step11_values.get("ld_variants") == str(EXPECTED_TOTAL_VARIANTS), "Step 1.1 variant count change")

    _, step12_rows = read_tsv(step12_summary)
    step12_values = {row.get("metric", ""): row.get("value", "") for row in step12_rows}
    require(step12_values.get("overall_status") == "PASS", "Step 1.2 failed")

    manifest_path = args.step11_dir / "input_file_manifest.tsv"
    header, manifest_rows = read_tsv(manifest_path)
    require({"path", "sha256"}.issubset(header), f"Step 1.1 manifest header mismatch: {header}")
    recorded_hash = {row["path"]: row["sha256"] for row in manifest_rows}
    for suffix in (".bed", ".bim", ".fam"):
        path = Path(f"{args.input_prefix}{suffix}")
        require(str(path) in recorded_hash, f"Step 1.1 manifest not recorded: {path}")
        require(sha256sum(path) == recorded_hash[str(path)], f"Input file has changed since Step 1.1: {path}")
    return read_plink2_version(args.plink2_bin)


def parse_release_inputs(
    args: argparse.Namespace,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]], dict[str, tuple[str, str, str, str]]]:
    fam_path = Path(f"{args.input_prefix}.fam")
    keep_ids: list[tuple[str, str]] = []
    group_remove_ids: list[tuple[str, str]] = []
    group_counts: Counter[str] = Counter()
    seen_samples: set[tuple[str, str]] = set()
    with fam_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            require(len(fields) == 6, f"{fam_path}:{line_number} does not have six columns FAM")
            fid, iid = fields[0], fields[1]
            sample = (fid, iid)
            require(sample not in seen_samples, f"duplicate sample: {fid}/{iid}")
            seen_samples.add(sample)
            group_counts[fid] += 1
            if fid in {"GRLS", "PILOT"}:
                keep_ids.append(sample)
            elif fid == "Oldies":
                group_remove_ids.append(sample)
            else:
                raise PipelineError(f"unknown value found FID group: {fid}")
    require(len(seen_samples) == EXPECTED_TOTAL_SAMPLES, "official release total samples change")
    require(dict(group_counts) == EXPECTED_GROUP_COUNTS, f"FID group counts have changed: {dict(group_counts)}")
    require(len(keep_ids) == 3024 and len(group_remove_ids) == 200, "retain/exclude unexpected sample count")

    bim_path = Path(f"{args.input_prefix}.bim")
    variant_info: dict[str, tuple[str, str, str, str]] = {}
    seen_variant_ids: set[str] = set()
    with bim_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.split()
            require(len(fields) == 6, f"{bim_path}:{line_number} does not have six columns BIM")
            chrom, variant_id, _, bp, allele1, allele2 = fields
            require(variant_id not in seen_variant_ids, f"duplicate variant ID: {variant_id}")
            seen_variant_ids.add(variant_id)
            match = AUTOSOME_RE.fullmatch(chrom)
            if match:
                require(allele1 != "0" and allele2 != "0", f"autosome variant contains missing allele: {variant_id}")
                variant_info[variant_id] = (match.group(1), bp, allele1, allele2)
    require(len(seen_variant_ids) == EXPECTED_TOTAL_VARIANTS, "official release total variant count has changed")
    require(len(variant_info) == EXPECTED_AUTOSOME_VARIANTS, "autosome variant count change")
    return keep_ids, group_remove_ids, variant_info


def parse_missingness(prefix: Path) -> tuple[dict[tuple[str, str], float], dict[str, float]]:
    smiss_header, smiss_rows = read_whitespace_table(Path(f"{prefix}.smiss"))
    vmiss_header, vmiss_rows = read_whitespace_table(Path(f"{prefix}.vmiss"))
    require({"FID", "IID", "F_MISS"}.issubset(smiss_header), f"SMISS header mismatch: {smiss_header}")
    require({"ID", "F_MISS"}.issubset(vmiss_header), f"VMISS header mismatch: {vmiss_header}")
    sample_missing: dict[tuple[str, str], float] = {}
    for row in smiss_rows:
        sample = (row["FID"], row["IID"])
        require(sample not in sample_missing, f"SMISS sample duplicate: {sample}")
        sample_missing[sample] = finite_float(row["F_MISS"], "sample F_MISS")
    variant_missing: dict[str, float] = {}
    for row in vmiss_rows:
        variant_id = row["ID"]
        require(variant_id not in variant_missing, f"VMISS variant duplicate: {variant_id}")
        variant_missing[variant_id] = finite_float(row["F_MISS"], "variant F_MISS")
    return sample_missing, variant_missing


def parse_maf(prefix: Path) -> dict[str, float | None]:
    header, rows = read_whitespace_table(Path(f"{prefix}.afreq"))
    require({"ID", "ALT_FREQS"}.issubset(header), f"AFREQ header mismatch: {header}")
    result: dict[str, float | None] = {}
    for row in rows:
        variant_id = row["ID"]
        require(variant_id not in result, f"AFREQ variant duplicate: {variant_id}")
        value = row["ALT_FREQS"]
        if value in MISSING_VALUES:
            result[variant_id] = None
            continue
        require("," not in value, f"input expected two allele, but ALT_FREQS contains multiple value: {variant_id}")
        alt_frequency = finite_float(value, "ALT frequency")
        require(0.0 <= alt_frequency <= 1.0, f"ALT frequency out of bounds: {variant_id}")
        result[variant_id] = min(alt_frequency, 1.0 - alt_frequency)
    return result


def parse_hwe(prefix: Path) -> dict[str, float | None]:
    header, rows = read_whitespace_table(Path(f"{prefix}.hardy"))
    require({"ID", "P"}.issubset(header), f"HARDY header mismatch: {header}")
    result: dict[str, float | None] = {}
    for row in rows:
        variant_id = row["ID"]
        require(variant_id not in result, f"HARDY variant duplicate: {variant_id}")
        value = row["P"]
        result[variant_id] = None if value in MISSING_VALUES else finite_float(value, "HWE P")
    return result


def count_prefix_groups(prefix: Path) -> Counter[str]:
    groups: Counter[str] = Counter()
    with Path(f"{prefix}.fam").open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                groups[line.split()[0]] += 1
    return groups


def validate_final_reference(
    args: argparse.Namespace,
    final_prefix: Path,
    final_report_prefix: Path,
    original_variant_info: dict[str, tuple[str, str, str, str]],
) -> dict[str, object]:
    final_sample_missing, final_variant_missing = parse_missingness(final_report_prefix)
    final_maf = parse_maf(final_report_prefix)
    final_hwe = parse_hwe(final_report_prefix)
    require(final_sample_missing and final_variant_missing, "final reference set no sample or variant")
    require(max(final_sample_missing.values()) <= args.mind, "final reference set still has sample missingness exceeds threshold")
    require(max(final_variant_missing.values()) <= args.geno, "final reference set still has variant missingness exceeds threshold")
    require(all(value is not None and value >= args.maf for value in final_maf.values()), "the final reference set still contains low MAF variant")
    require(all(value is not None and value >= args.hwe for value in final_hwe.values()), "final reference set still has HWE variants failing the criteria")

    groups = count_prefix_groups(final_prefix)
    require(set(groups).issubset({"GRLS", "PILOT"}), f"final reference set contains error FID: {dict(groups)}")
    require("Oldies" not in groups, "final reference set still contains Oldies")

    final_variant_count = 0
    seen_variants: set[str] = set()
    with Path(f"{final_prefix}.bim").open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            require(len(fields) == 6, f"final BIM row {line_number} does not have six columns")
            chrom, variant_id, _, bp, allele1, allele2 = fields
            require(chrom.isdigit() and 1 <= int(chrom) <= 38, f"final BIM non-autosomal: {chrom}")
            require(variant_id not in seen_variants, f"final BIM duplicate variant ID: {variant_id}")
            seen_variants.add(variant_id)
            require(variant_id in original_variant_info, f"final BIM encountered input outside variant: {variant_id}")
            expected = original_variant_info[variant_id]
            require((chrom, bp, allele1, allele2) == expected, f"final BIM coordinates or A1/A2 change: {variant_id}")
            final_variant_count += 1

    sample_count, counted_variants = prefix_counts(final_prefix)
    require(counted_variants == final_variant_count, "final BIM count mismatch")
    require(set(final_variant_missing) == seen_variants, "Variant sets differ between the final VMISS and BIM")
    require(set(final_maf) == seen_variants, "Variant sets differ between the final AFREQ and BIM")
    require(set(final_hwe) == seen_variants, "Variant sets differ between the final HARDY and BIM")

    bed_path = Path(f"{final_prefix}.bed")
    with bed_path.open("rb") as handle:
        bed_magic = handle.read(3)
    require(bed_magic == b"\x6c\x1b\x01", "final BED is not PLINK SNP-major format")
    expected_bed_size = 3 + final_variant_count * ((sample_count + 3) // 4)
    actual_bed_size = bed_path.stat().st_size
    require(actual_bed_size == expected_bed_size, f"final BED size does not match: actual={actual_bed_size}, expected={expected_bed_size}")
    return {
        "samples": sample_count,
        "variants": final_variant_count,
        "groups": dict(sorted(groups.items())),
        "max_sample_missingness": max(final_sample_missing.values()),
        "max_variant_missingness": max(final_variant_missing.values()),
        "min_maf": min(value for value in final_maf.values() if value is not None),
        "min_hwe_p": min(value for value in final_hwe.values() if value is not None),
        "bed_size_bytes": actual_bed_size,
    }


def build_stable_final_reference(
    args: argparse.Namespace,
    input_prefix: Path,
    final_prefix: Path,
    intermediate: Path,
    original_variant_info: dict[str, tuple[str, str, str, str]],
    joint_round: int,
) -> tuple[dict[str, object], set[tuple[str, str]], dict[str, list[str]], int, list[list[object]]]:
    current_prefix = input_prefix
    all_sample_removals: set[tuple[str, str]] = set()
    all_variant_reasons: dict[str, list[str]] = defaultdict(list)
    audit_rows: list[list[object]] = []
    stable_iteration: int | None = None
    round_label = f"joint_round_{joint_round:02d}"

    for iteration in range(1, args.max_final_qc_iterations + 1):
        report_prefix = intermediate / f"10_{round_label}_qc_iteration_{iteration:02d}_reports"
        run_missing(args, current_prefix, report_prefix)
        run_freq(args, current_prefix, report_prefix)
        run_hardy(args, current_prefix, report_prefix)
        sample_missing, variant_missing = parse_missingness(report_prefix)
        maf = parse_maf(report_prefix)
        hwe = parse_hwe(report_prefix)
        expected_variants = set(variant_missing)
        require(set(maf) == expected_variants, "Final QC: variant sets differ between AFREQ and VMISS")
        require(set(hwe) == expected_variants, "Final QC: variant sets differ between HARDY and VMISS")

        sample_removals = {sample for sample, value in sample_missing.items() if value > args.mind}
        variant_reasons: dict[str, list[str]] = defaultdict(list)
        for variant_id, missing_value in variant_missing.items():
            if missing_value > args.geno:
                variant_reasons[variant_id].append("GENO")
            maf_value = maf[variant_id]
            if maf_value is None or maf_value < args.maf:
                variant_reasons[variant_id].append("MAF")
            hwe_value = hwe[variant_id]
            if hwe_value is None or hwe_value < args.hwe:
                variant_reasons[variant_id].append("HWE")

        current_counts = prefix_counts(current_prefix)
        audit_rows.append(
            [joint_round, iteration, current_prefix, current_counts[0], current_counts[1], len(sample_removals), len(variant_reasons)]
        )
        if not sample_removals and not variant_reasons:
            stable_iteration = iteration
            make_bed(args, current_prefix, final_prefix)
            break

        all_sample_removals.update(sample_removals)
        for variant_id, reasons in variant_reasons.items():
            for reason in reasons:
                if reason not in all_variant_reasons[variant_id]:
                    all_variant_reasons[variant_id].append(reason)
        sample_remove_path = intermediate / f"10_{round_label}_qc_iteration_{iteration:02d}.sample_remove.txt"
        variant_remove_path = intermediate / f"10_{round_label}_qc_iteration_{iteration:02d}.variant_remove.txt"
        write_id_pairs(sample_remove_path, sample_removals)
        write_variant_ids(variant_remove_path, variant_reasons)
        next_prefix = intermediate / f"10_{round_label}_qc_iteration_{iteration:02d}_filtered"
        make_bed(
            args,
            current_prefix,
            next_prefix,
            remove=sample_remove_path,
            exclude=variant_remove_path,
        )
        next_counts = prefix_counts(next_prefix)
        require(next_counts[0] > 0 and next_counts[1] > 0, "iteration QC after no available sample or variant")
        require(next_counts[0] <= current_counts[0] and next_counts[1] <= current_counts[1], "iteration QC count unexpected increase")
        current_prefix = next_prefix

    require(stable_iteration is not None, f"final QC in {args.max_final_qc_iterations} rounds without stabilization")
    validation_prefix = intermediate / f"11_{round_label}_final_validation"
    run_missing(args, final_prefix, validation_prefix)
    run_freq(args, final_prefix, validation_prefix)
    run_hardy(args, final_prefix, validation_prefix)
    validation = validate_final_reference(args, final_prefix, validation_prefix, original_variant_info)
    return validation, all_sample_removals, all_variant_reasons, stable_iteration, audit_rows


def main() -> int:
    args = parse_args()
    plink2_version = validate_prerequisites(args)
    require(not args.output_dir.exists(), f"refuse to overwrite existing output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    intermediate = args.output_dir / "intermediate"
    intermediate.mkdir()

    log("parse and check official PLINK release")
    keep_ids, group_remove_ids, original_variant_info = parse_release_inputs(args)
    initial_keep = args.output_dir / "ld_sample_keep.initial_grls_pilot.txt"
    group_remove = args.output_dir / "ld_sample_remove.group_oldies.txt"
    autosome_extract = intermediate / "ld_autosome_extract.ids"
    chromosome_update = intermediate / "ld_chromosome_update.txt"
    write_id_pairs(initial_keep, keep_ids)
    write_id_pairs(group_remove, group_remove_ids)
    write_variant_ids(autosome_extract, original_variant_info)
    with chromosome_update.open("w", encoding="utf-8") as handle:
        for variant_id in sorted(original_variant_info):
            handle.write(f"{original_variant_info[variant_id][0]}\t{variant_id}\n")

    stage_counts: list[list[object]] = []
    input_counts = prefix_counts(args.input_prefix)
    require(input_counts == (EXPECTED_TOTAL_SAMPLES, EXPECTED_TOTAL_VARIANTS), f"input count unexpected: {input_counts}")
    stage_counts.append(["input_official_release", *input_counts])

    stage1 = intermediate / "01_grls_pilot_autosome"
    make_bed(
        args,
        args.input_prefix,
        stage1,
        keep=initial_keep,
        extract=autosome_extract,
        update_chr=chromosome_update,
    )
    stage1_counts = prefix_counts(stage1)
    require(stage1_counts == (3024, EXPECTED_AUTOSOME_VARIANTS), f"GRLS+PILOT autosome count unexpected: {stage1_counts}")
    stage_counts.append(["grls_pilot_autosome", *stage1_counts])

    preliminary_missing_prefix = intermediate / "02_preliminary_missingness"
    run_missing(args, stage1, preliminary_missing_prefix)
    preliminary_sample_missing, preliminary_variant_missing = parse_missingness(preliminary_missing_prefix)
    preliminary_sample_remove = {sample for sample, value in preliminary_sample_missing.items() if value > args.mind}
    preliminary_variant_remove = {variant_id for variant_id, value in preliminary_variant_missing.items() if value > args.geno}
    sample_missing_path = args.output_dir / "ld_sample_remove.preliminary_missingness.txt"
    variant_missing_path = args.output_dir / "ld_variant_remove.preliminary_missingness.txt"
    write_id_pairs(sample_missing_path, preliminary_sample_remove)
    write_variant_ids(variant_missing_path, preliminary_variant_remove)

    stage2 = intermediate / "03_preliminary_missingness_qc"
    make_bed(args, stage1, stage2, remove=sample_missing_path, exclude=variant_missing_path)
    stage2_counts = prefix_counts(stage2)
    require(stage2_counts[0] > 0 and stage2_counts[1] > 0, "missingness pre-cleaning after no data")
    stage_counts.append(["preliminary_missingness_qc", *stage2_counts])

    preliminary_freq_prefix = intermediate / "04_preliminary_frequency"
    run_freq(args, stage2, preliminary_freq_prefix)
    preliminary_maf = parse_maf(preliminary_freq_prefix)
    preliminary_maf_remove = {
        variant_id for variant_id, value in preliminary_maf.items() if value is None or value < args.maf
    }
    preliminary_maf_path = args.output_dir / "ld_variant_remove.preliminary_maf.txt"
    write_variant_ids(preliminary_maf_path, preliminary_maf_remove)

    preliminary_qc = intermediate / "05_preliminary_maf_qc"
    make_bed(args, stage2, preliminary_qc, exclude=preliminary_maf_path)
    preliminary_qc_counts = prefix_counts(preliminary_qc)
    require(preliminary_qc_counts[0] > 0 and preliminary_qc_counts[1] > 0, "pre-cleaning after no data")
    stage_counts.append(["preliminary_maf_qc", *preliminary_qc_counts])

    prune_prefix = intermediate / "06_ld_prune_for_king"
    prune_command = plink2_base(args, preliminary_qc, prune_prefix)
    prune_command.extend(
        [
            "--indep-pairwise",
            args.prune_window,
            args.prune_step,
            args.prune_r2,
            "--indep-order",
            "1",
        ]
    )
    run_command(prune_command, "PLINK 2 LD pruning for KING")
    prune_in = Path(f"{prune_prefix}.prune.in")
    require(prune_in.is_file() and prune_in.stat().st_size > 0, "LD pruning not generated available marker")
    with prune_in.open("r", encoding="utf-8") as handle:
        pruned_variants = {line.strip() for line in handle if line.strip()}
    require(pruned_variants.issubset(read_bim_ids(preliminary_qc)), "prune.in contains pre-cleaning data external variants")

    king_prefix = intermediate / "07_king_primary"
    king_command = plink2_base(args, preliminary_qc, king_prefix)
    king_command.extend(
        ["--extract", prune_in, "--king-cutoff", args.king_relatedness, "--no-id-header"]
    )
    run_command(king_command, "PLINK 2 KING-robust primary pruning")
    king_keep_source = Path(f"{king_prefix}.king.cutoff.in.id")
    king_remove_source = Path(f"{king_prefix}.king.cutoff.out.id")
    king_kept = read_id_pairs(king_keep_source)
    king_removals = read_id_pairs(king_remove_source, allow_empty=True)
    preliminary_samples = read_fam_ids(preliminary_qc)
    require(not (king_kept & king_removals), "KING retain with remove sample set overlap")
    require(king_kept | king_removals == preliminary_samples, "KING retained/removed sets do not match the FAM before filtering")
    king_initial_keep_path = args.output_dir / "ld_plink2_king_primary_sample_keep.initial_king.txt"
    king_initial_remove_path = args.output_dir / "ld_plink2_king_primary_sample_remove.initial_king.txt"
    write_id_pairs(king_initial_keep_path, king_kept)
    write_id_pairs(king_initial_remove_path, king_removals)

    post_related = intermediate / "08_king_primary_post_relatedness"
    make_bed(args, preliminary_qc, post_related, remove=king_initial_remove_path)
    post_related_counts = prefix_counts(post_related)
    require(post_related_counts[0] == len(king_kept), "FAM sample count after KING relatedness filtering differs from the retained-set count")
    stage_counts.append(["king_primary_post_relatedness", *post_related_counts])

    final_prefix = args.output_dir / "ld_reference_qc_plink2_king_primary"
    joint_input_prefix = post_related
    post_qc_king_removals: set[tuple[str, str]] = set()
    final_sample_removals: set[tuple[str, str]] = set()
    final_variant_reasons: dict[str, list[str]] = defaultdict(list)
    final_qc_audit_rows: list[list[object]] = []
    stable_qc_iterations: list[int] = []
    final_king_kept: set[tuple[str, str]] = set()
    joint_stable_round: int | None = None

    for joint_round in range(1, args.max_final_qc_iterations + 1):
        validation, round_sample_removals, round_variant_reasons, stable_iteration, round_audit = (
            build_stable_final_reference(
                args,
                joint_input_prefix,
                final_prefix,
                intermediate,
                original_variant_info,
                joint_round,
            )
        )
        final_sample_removals.update(round_sample_removals)
        for variant_id, reasons in round_variant_reasons.items():
            for reason in reasons:
                if reason not in final_variant_reasons[variant_id]:
                    final_variant_reasons[variant_id].append(reason)
        final_qc_audit_rows.extend(round_audit)
        stable_qc_iterations.append(stable_iteration)

        final_variants = read_bim_ids(final_prefix)
        final_king_variants = pruned_variants & final_variants
        require(final_king_variants, "final reference set no available in KING reviewed LD-pruned marker")
        final_king_extract = intermediate / f"12_joint_round_{joint_round:02d}_king_extract.ids"
        write_variant_ids(final_king_extract, final_king_variants)
        final_king_prefix = intermediate / f"12_joint_round_{joint_round:02d}_king_validation"
        final_king_command = plink2_base(args, final_prefix, final_king_prefix)
        final_king_command.extend(
            ["--extract", final_king_extract, "--king-cutoff", args.king_relatedness, "--no-id-header"]
        )
        run_command(final_king_command, f"PLINK 2 joint round {joint_round} KING-robust validation")
        final_king_kept = read_id_pairs(Path(f"{final_king_prefix}.king.cutoff.in.id"))
        round_king_removed = read_id_pairs(
            Path(f"{final_king_prefix}.king.cutoff.out.id"), allow_empty=True
        )
        require(not (final_king_kept & round_king_removed), "final KING retain with remove sample set overlap")
        require(
            final_king_kept | round_king_removed == read_fam_ids(final_prefix),
            "Final KING retained/removed sets do not match the FAM",
        )
        if not round_king_removed:
            joint_stable_round = joint_round
            break

        post_qc_king_removals.update(round_king_removed)
        round_remove_path = intermediate / f"13_joint_round_{joint_round:02d}_king_remove.txt"
        write_id_pairs(round_remove_path, round_king_removed)
        next_joint_input = intermediate / f"13_joint_round_{joint_round:02d}_post_king"
        make_bed(args, final_prefix, next_joint_input, remove=round_remove_path)
        next_counts = prefix_counts(next_joint_input)
        require(next_counts[0] > 0 and next_counts[1] > 0, "KING no samples or variants remain after joint iteration")
        require(next_counts[0] < validation["samples"], "KING sample removal was reported but the sample count did not decrease")
        joint_input_prefix = next_joint_input

    require(
        joint_stable_round is not None,
        f"KING and variants QC in {args.max_final_qc_iterations} joint iterations without stabilization",
    )
    require(final_king_kept == read_fam_ids(final_prefix), "Final KING retained set does not match the FAM")

    all_king_removals = king_removals | post_qc_king_removals
    write_id_pairs(args.output_dir / "ld_plink2_king_primary_sample_keep.txt", final_king_kept)
    write_id_pairs(args.output_dir / "ld_plink2_king_primary_sample_remove.relatedness.txt", all_king_removals)
    write_id_pairs(
        args.output_dir / "ld_plink2_king_primary_sample_remove.post_qc_relatedness.txt",
        post_qc_king_removals,
    )
    write_id_pairs(args.output_dir / "ld_plink2_king_primary_sample_remove.final_qc.txt", final_sample_removals)
    write_variant_ids(args.output_dir / "ld_plink2_king_primary_variant_remove.final_qc.txt", final_variant_reasons)
    write_tsv(
        args.output_dir / "ld_plink2_king_primary_variant_remove.final_qc_reasons.tsv",
        ["variant_id", "reasons"],
        [[variant_id, ";".join(reasons)] for variant_id, reasons in sorted(final_variant_reasons.items())],
    )
    write_tsv(
        args.output_dir / "ld_final_qc_iteration_counts.tsv",
        ["joint_round", "qc_iteration", "input_prefix", "sample_count", "variant_count", "sample_removals", "variant_removals"],
        final_qc_audit_rows,
    )
    final_counts = prefix_counts(final_prefix)
    stage_counts.append(["plink2_king_primary_final_qc", *final_counts])

    write_tsv(args.output_dir / "ld_qc_stage_counts.tsv", ["stage", "sample_count", "variant_count"], stage_counts)
    write_tsv(
        args.output_dir / "pipeline_parameters.tsv",
        ["parameter", "value"],
        [
            ["input_prefix", args.input_prefix],
            ["sample_groups_retained", "GRLS;PILOT"],
            ["sample_group_excluded", "Oldies"],
            ["chromosomes", "1-38"],
            ["mind", args.mind],
            ["geno", args.geno],
            ["maf", args.maf],
            ["hwe", args.hwe],
            ["prune_window", args.prune_window],
            ["prune_step", args.prune_step],
            ["prune_r2", args.prune_r2],
            ["indep_order", "1 (imitate PLINK 1.x)"],
            ["relatedness_method", "PLINK2_KING_robust"],
            ["king_relatedness_rule", f"kinship>{args.king_relatedness}"],
            ["king_selection_rule", "PLINK2_king_cutoff_greedy_iterated_with_variant_QC_until_joint_stability"],
            ["reference_allele_handling", "input_A2_is_provisional_REF; no_maj_ref; no_real_ref_assertion"],
            ["allele_order_validation", "final_BIM_exactly_matches_input_A1_A2"],
            ["plink2_bin", args.plink2_bin],
            ["plink2_version", plink2_version],
            ["threads", args.threads],
        ],
    )
    write_tsv(
        args.output_dir / "input_lineage.tsv",
        ["role", "path", "producer_or_evidence"],
        [
            ["official_plink_bed", f"{args.input_prefix}.bed", args.step11_dir / "input_file_manifest.tsv"],
            ["official_plink_bim", f"{args.input_prefix}.bim", args.step11_dir / "input_file_manifest.tsv"],
            ["official_plink_fam", f"{args.input_prefix}.fam", args.step11_dir / "input_file_manifest.tsv"],
            ["step_1_1_validation", args.step11_dir / "validation_summary.tsv", args.step11_dir / "COMPLETED.ok"],
            ["step_1_2_manifest", args.step12_dir / "manifest_summary.tsv", args.step12_dir / "COMPLETED.ok"],
        ],
    )

    summary_rows = [
        ["overall_status", "PASS", "PASS", "all PLINK 2/KING primary reference invariants satisfied"],
        ["input_samples", input_counts[0], "PASS", "official merged release"],
        ["input_variants", input_counts[1], "PASS", "official merged release"],
        ["initial_keep_samples", len(keep_ids), "PASS", "GRLS+PILOT"],
        ["excluded_oldies_samples", len(group_remove_ids), "PASS", "FID=Oldies"],
        ["autosome_variants", len(original_variant_info), "PASS", "chr1-chr38"],
        ["preliminary_missingness_sample_removals", len(preliminary_sample_remove), "INFO", f"F_MISS>{args.mind}"],
        ["preliminary_missingness_variant_removals", len(preliminary_variant_remove), "INFO", f"F_MISS>{args.geno}"],
        ["preliminary_maf_variant_removals", len(preliminary_maf_remove), "INFO", f"MAF<{args.maf}"],
        ["ld_pruned_variants_for_king", len(pruned_variants), "PASS", f"{args.prune_window}/{args.prune_step}/{args.prune_r2}; indep-order=1"],
        ["king_relatedness_threshold", args.king_relatedness, "INFO", "KING-robust kinship coefficient"],
        ["king_initial_sample_removals", len(king_removals), "INFO", "initial PLINK 2 --king-cutoff greedy selection"],
        ["king_post_qc_sample_removals", len(post_qc_king_removals), "INFO", "additional removals after the KING marker set changed during variant QC"],
        ["king_relatedness_sample_removals", len(all_king_removals), "INFO", "union across all KING/QC joint rounds"],
        ["king_qc_joint_rounds_to_stability", joint_stable_round, "PASS", "KING and variant QC jointly stable"],
        ["final_qc_iterations_to_stability", sum(stable_qc_iterations), "PASS", "total QC iterations across KING/QC joint rounds"],
        ["final_qc_sample_removals", len(final_sample_removals), "INFO", f"F_MISS>{args.mind}"],
        ["final_qc_variant_removals", len(final_variant_reasons), "INFO", "GENO/MAF/HWE union across iterations"],
        ["final_samples", validation["samples"], "PASS", str(validation["groups"])],
        ["final_variants", validation["variants"], "PASS", "autosomal PLINK 2/KING primary LD reference"],
        ["final_bed_size_bytes", validation["bed_size_bytes"], "PASS", "matches SNP-major theoretical size"],
        ["final_max_sample_missingness", f"{validation['max_sample_missingness']:.12g}", "PASS", f"<={args.mind}"],
        ["final_max_variant_missingness", f"{validation['max_variant_missingness']:.12g}", "PASS", f"<={args.geno}"],
        ["final_min_maf", f"{validation['min_maf']:.12g}", "PASS", f">={args.maf}"],
        ["final_min_hwe_p", f"{validation['min_hwe_p']:.12g}", "PASS", f">={args.hwe}"],
        ["final_allele_order_preserved", "yes", "PASS", "final BIM matches input CHROM/POS/A1/A2"],
        ["final_unresolved_king_relatives", 0, "PASS", f"no kinship>{args.king_relatedness} pair"],
    ]
    write_tsv(args.output_dir / "ld_reference_qc_summary.tsv", ["metric", "value", "status", "details"], summary_rows)
    with (args.output_dir / "COMPLETED.ok").open("w", encoding="utf-8") as handle:
        handle.write(f"completed_at\t{datetime.now():%Y-%m-%d %H:%M:%S}\n")
        handle.write(f"primary_prefix\t{final_prefix}\n")
        handle.write(f"primary_samples\t{validation['samples']}\n")
        handle.write(f"primary_variants\t{validation['variants']}\n")
        handle.write(f"relatedness_method\tPLINK2_KING_robust\n")
        handle.write(f"king_threshold\t{args.king_relatedness}\n")

    log(f"PASS primary_prefix={final_prefix}")
    log(f"PASS samples={validation['samples']} variants={validation['variants']}")
    log("task completed / done")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PipelineError as exc:
        log(f"ERROR {exc}")
        raise SystemExit(1) from exc
