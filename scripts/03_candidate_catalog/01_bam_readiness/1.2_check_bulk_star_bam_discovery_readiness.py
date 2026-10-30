#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


PROJECT_ROOT = Path("/path/to/project/legacy/five-dogs")
DEFAULT_SAMPLE_LIST = PROJECT_ROOT / "data/statistics/five_dogs_fastq_pair_stats.tsv"
DEFAULT_STAR_ROOT = PROJECT_ROOT / "5.cnv/1.star-v5"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "9.new-gene/1.2_check_bulk_star_bam_discovery_readiness-v3"
DEFAULT_SAMTOOLS = Path("/path/to/software/miniconda3/bin/samtools")
DEFAULT_CHROM_SIZES = Path(
    "/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"
)

READINESS_FIELDS = [
    "animal_id",
    "tissue",
    "sample_tissue",
    "raw_data_name",
    "star_sample_dir",
    "align_star_dir",
    "raw_star_bam",
    "raw_star_bam_exists",
    "raw_star_bam_size_bytes",
    "raw_star_bam_size_status",
    "raw_star_bam_index_exists",
    "raw_star_bam_index_path",
    "star_log_final",
    "star_log_final_exists",
    "star_input_reads",
    "star_uniquely_mapped_reads_number",
    "star_uniquely_mapped_reads_percent",
    "star_multimapping_percent",
    "star_unmapped_too_short_percent",
    "star_total_splices",
    "sj_out_tab",
    "sj_out_tab_exists",
    "sj_out_line_count",
    "samtools_quickcheck_status",
    "samtools_quickcheck_message",
    "bam_sort_order",
    "is_genome_coordinate_bam",
    "sq_contig_count",
    "chrom_sizes",
    "chrom_sizes_contig_count",
    "bam_chrom_size_match_status",
    "bam_sq_missing_in_chrom_sizes_count",
    "bam_sq_missing_in_chrom_sizes_examples",
    "chrom_sizes_missing_in_bam_sq_count",
    "chrom_sizes_missing_in_bam_sq_examples",
    "bam_sq_length_mismatch_count",
    "bam_sq_length_mismatch_examples",
    "sampled_alignment_count",
    "has_spliced_cigar_N_in_sample",
    "spliced_alignment_sample_count",
    "spliced_alignment_with_XS_count",
    "has_XS_tag_in_spliced_reads",
    "xs_tag_status",
    "splice_evidence_status",
    "raw_bam_content_usable_for_discovery",
    "raw_bam_ready_for_region_query",
    "raw_bam_needs_index",
    "can_use_for_novel_transcript_discovery",
    "recommended_for_stringtie_discovery",
    "recommended_for_stringtie3_discovery",
    "cannot_use_reason",
    "warning_reason",
]

SUMMARY_FIELDS = ["metric", "value", "note"]
MANIFEST_FIELDS = [
    "sample_list",
    "star_root",
    "output_root",
    "samtools",
    "chrom_sizes",
    "max_alignments",
    "workers",
    "min_bam_size_bytes",
    "command_timeout_sec",
    "readiness_table",
    "summary_table",
    "manifest_table",
]


@dataclass(frozen=True)
class SampleRecord:
    animal_id: str
    tissue: str
    sample_tissue: str
    raw_data_name: str


@dataclass(frozen=True)
class SampleReadiness:
    row: Dict[str, object]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check whether the five-dog bulk RNA-seq STAR coordinate-sorted BAMs are suitable for novel transcript discovery."
    )
    parser.add_argument("--sample-list", type=Path, default=DEFAULT_SAMPLE_LIST)
    parser.add_argument("--star-root", type=Path, default=DEFAULT_STAR_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--samtools", type=Path, default=DEFAULT_SAMTOOLS)
    parser.add_argument("--chrom-sizes", type=Path, default=DEFAULT_CHROM_SIZES)
    parser.add_argument("--workers", type=int, default=19)
    parser.add_argument("--max-alignments", type=int, default=100_000)
    parser.add_argument("--min-bam-size-bytes", type=int, default=10_000_000)
    parser.add_argument("--command-timeout-sec", type=int, default=300)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%F %T",
    )


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return text if text else "NA"


def check_fields(fieldnames: Sequence[str] | None, required: Sequence[str], source: Path) -> None:
    if fieldnames is None:
        raise ValueError(f"cannot read column header: {source}")
    missing = [field for field in required if field not in fieldnames]
    if missing:
        raise ValueError(f"{source} missing required fields: {','.join(missing)}")


def read_sample_list(path: Path) -> List[SampleRecord]:
    records: List[SampleRecord] = []
    seen: set[str] = set()
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check_fields(reader.fieldnames, ["animal_id", "term", "raw_data_name"], path)
        for row_number, row in enumerate(reader, start=2):
            animal_id = clean(row.get("animal_id"))
            tissue = clean(row.get("term"))
            raw_data_name = clean(row.get("raw_data_name"))
            if animal_id == "NA" or tissue == "NA":
                raise ValueError(f"{path} row {row_number} is missing animal_id or term")
            sample_tissue = f"{animal_id}_{tissue}"
            if sample_tissue in seen:
                raise ValueError(f"{path} Duplicate sample_tissue at row {row_number}: {sample_tissue}")
            seen.add(sample_tissue)
            records.append(
                SampleRecord(
                    animal_id=animal_id,
                    tissue=tissue,
                    sample_tissue=sample_tissue,
                    raw_data_name=raw_data_name,
                )
            )
    if not records:
        raise ValueError(f"sample list is empty: {path}")
    return records


def align_star_dir(star_root: Path, sample: SampleRecord) -> Path:
    return star_root / sample.animal_id / sample.tissue / "align_star"


def star_bam_path(star_root: Path, sample: SampleRecord) -> Path:
    return align_star_dir(star_root, sample) / f"{sample.sample_tissue}.Aligned.sortedByCoord.out.bam"


def star_log_final_path(star_root: Path, sample: SampleRecord) -> Path:
    return align_star_dir(star_root, sample) / f"{sample.sample_tissue}.Log.final.out"


def sj_out_path(star_root: Path, sample: SampleRecord) -> Path:
    return align_star_dir(star_root, sample) / f"{sample.sample_tissue}.SJ.out.tab"


def bam_index_candidates(bam_path: Path) -> List[Path]:
    return [
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
        Path(f"{bam_path}.csi"),
        bam_path.with_suffix(".csi"),
    ]


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def run_command(cmd: Sequence[str], timeout_sec: int) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(cmd),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_sec,
        check=False,
    )


def parse_star_log(path: Path) -> Dict[str, str]:
    metrics = {
        "star_input_reads": "NA",
        "star_uniquely_mapped_reads_number": "NA",
        "star_uniquely_mapped_reads_percent": "NA",
        "star_multimapping_percent": "NA",
        "star_unmapped_too_short_percent": "NA",
        "star_total_splices": "NA",
    }
    if not path.exists():
        return metrics
    key_map = {
        "Number of input reads": "star_input_reads",
        "Uniquely mapped reads number": "star_uniquely_mapped_reads_number",
        "Uniquely mapped reads %": "star_uniquely_mapped_reads_percent",
        "% of reads mapped to multiple loci": "star_multimapping_percent",
        "% of reads unmapped: too short": "star_unmapped_too_short_percent",
        "Number of splices: Total": "star_total_splices",
    }
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = [part.strip() for part in line.rstrip("\n").split("|", 1)]
            if key in key_map:
                metrics[key_map[key]] = clean(value).rstrip("%")
    return metrics


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    chrom_sizes: Dict[str, int] = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path} row {row_number} has an invalid .fai/chromosome-size format")
            chrom = clean(parts[0])
            if chrom in chrom_sizes:
                raise ValueError(f"{path} row {row_number} contains a duplicate chromosome name: {chrom}")
            try:
                chrom_sizes[chrom] = int(parts[1])
            except ValueError as exc:
                raise ValueError(f"{path} row {row_number} chromosome length is not an integer: {parts[1]}") from exc
    if not chrom_sizes:
        raise ValueError(f"chrom sizes file is empty: {path}")
    return chrom_sizes


def parse_header(header_text: str) -> Tuple[str, Dict[str, int]]:
    sort_order = "NA"
    sq_lengths: Dict[str, int] = {}
    for line in header_text.splitlines():
        if line.startswith("@HD"):
            for part in line.split("\t"):
                if part.startswith("SO:"):
                    sort_order = part.removeprefix("SO:")
        elif line.startswith("@SQ"):
            name = "NA"
            length: int | None = None
            for part in line.split("\t"):
                if part.startswith("SN:"):
                    name = part.removeprefix("SN:")
                elif part.startswith("LN:"):
                    try:
                        length = int(part.removeprefix("LN:"))
                    except ValueError:
                        length = None
            if name != "NA" and length is not None:
                sq_lengths[name] = length
    return sort_order, sq_lengths


def format_examples(values: Sequence[str], limit: int = 10) -> str:
    if not values:
        return "NA"
    return "|".join(values[:limit])


def compare_chrom_sizes(
    bam_sq_lengths: Mapping[str, int],
    chrom_sizes: Mapping[str, int],
) -> Dict[str, object]:
    if not bam_sq_lengths:
        return {
            "bam_chrom_size_match_status": "not_tested",
            "bam_sq_missing_in_chrom_sizes_count": "NA",
            "bam_sq_missing_in_chrom_sizes_examples": "NA",
            "chrom_sizes_missing_in_bam_sq_count": "NA",
            "chrom_sizes_missing_in_bam_sq_examples": "NA",
            "bam_sq_length_mismatch_count": "NA",
            "bam_sq_length_mismatch_examples": "NA",
        }

    bam_names = set(bam_sq_lengths)
    chrom_names = set(chrom_sizes)
    missing_in_chrom_sizes = sorted(bam_names - chrom_names)
    missing_in_bam = sorted(chrom_names - bam_names)
    length_mismatches = sorted(
        name for name in bam_names & chrom_names if bam_sq_lengths[name] != chrom_sizes[name]
    )
    status = "match"
    if missing_in_chrom_sizes or missing_in_bam or length_mismatches:
        status = "mismatch"

    return {
        "bam_chrom_size_match_status": status,
        "bam_sq_missing_in_chrom_sizes_count": len(missing_in_chrom_sizes),
        "bam_sq_missing_in_chrom_sizes_examples": format_examples(missing_in_chrom_sizes),
        "chrom_sizes_missing_in_bam_sq_count": len(missing_in_bam),
        "chrom_sizes_missing_in_bam_sq_examples": format_examples(missing_in_bam),
        "bam_sq_length_mismatch_count": len(length_mismatches),
        "bam_sq_length_mismatch_examples": format_examples(length_mismatches),
    }


def count_lines(path: Path) -> int | str:
    if not path.exists():
        return "NA"
    total = 0
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for _line in handle:
            total += 1
    return total


def int_or_none(value: object) -> int | None:
    text = clean(value).replace(",", "")
    if text == "NA":
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def inspect_alignment_sample(
    samtools: Path,
    bam_path: Path,
    max_alignments: int,
    timeout_sec: int,
) -> Tuple[int, int, int, int, str]:
    cmd = [str(samtools), "view", str(bam_path)]
    sampled = 0
    spliced_sampled = 0
    spliced_with_xs = 0
    timed_out = 0
    start_time = time.monotonic()
    proc = subprocess.Popen(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.stdout is not None
    assert proc.stderr is not None
    try:
        for line in proc.stdout:
            if time.monotonic() - start_time > timeout_sec:
                timed_out = 1
                break
            if sampled >= max_alignments:
                break
            if not line or line.startswith("@"):
                continue
            sampled += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 5 and "N" in parts[5]:
                spliced_sampled += 1
                if any(field.startswith("XS:A:") for field in parts[11:]):
                    spliced_with_xs += 1
                if spliced_sampled >= 20 and spliced_with_xs > 0:
                    break
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=5)
        stderr = clean(proc.stderr.read())
    except Exception:
        if proc.poll() is None:
            proc.kill()
        raise

    has_spliced = int(spliced_sampled > 0)
    if timed_out:
        return sampled, spliced_sampled, spliced_with_xs, has_spliced, f"alignment_sample_timeout_after_{timeout_sec}s"
    if sampled == 0 and proc.returncode not in (0, -15):
        return sampled, spliced_sampled, spliced_with_xs, has_spliced, stderr
    return sampled, spliced_sampled, spliced_with_xs, has_spliced, "NA"


def inspect_sample(task: Tuple[SampleRecord, Path, Path, Path, Mapping[str, int], int, int, int]) -> SampleReadiness:
    sample, star_root, samtools, chrom_sizes_path, chrom_sizes, max_alignments, min_bam_size_bytes, timeout_sec = task
    sample_dir = star_root / sample.animal_id / sample.tissue
    align_dir = align_star_dir(star_root, sample)
    bam_path = star_bam_path(star_root, sample)
    log_path = star_log_final_path(star_root, sample)
    sj_path = sj_out_path(star_root, sample)
    bam_exists = int(bam_path.exists())
    bam_size = bam_path.stat().st_size if bam_exists else None
    size_status = "not_found"
    if bam_size is not None:
        size_status = "ok" if bam_size >= min_bam_size_bytes else "too_small"

    index_path = first_existing(bam_index_candidates(bam_path)) if bam_exists else None
    star_metrics = parse_star_log(log_path)
    sj_line_count = count_lines(sj_path)
    quickcheck_status = "not_tested"
    quickcheck_message = "NA"
    sort_order = "NA"
    sq_count: int | None = None
    sq_lengths: Dict[str, int] = {}
    sampled = 0
    has_spliced = 0
    spliced_sampled = 0
    spliced_with_xs = 0
    has_xs_tag = 0
    xs_tag_status = "not_tested"
    reasons: List[str] = []
    warnings: List[str] = []

    if not bam_exists:
        reasons.append("raw_star_bam_missing")
    else:
        quickcheck = run_command([str(samtools), "quickcheck", "-v", str(bam_path)], timeout_sec)
        quickcheck_status = "pass" if quickcheck.returncode == 0 else "fail"
        quickcheck_message = clean(quickcheck.stderr)
        if quickcheck.returncode != 0:
            reasons.append("samtools_quickcheck_fail")

        header = run_command([str(samtools), "view", "-H", str(bam_path)], timeout_sec)
        if header.returncode == 0:
            sort_order, sq_lengths = parse_header(header.stdout)
            sq_count = len(sq_lengths)
        else:
            reasons.append("samtools_header_fail")

        sampled, spliced_sampled, spliced_with_xs, has_spliced, sample_message = inspect_alignment_sample(
            samtools, bam_path, max_alignments, timeout_sec
        )
        if sample_message != "NA":
            warnings.append(f"samtools_alignment_sample_warning:{sample_message}")
        has_xs_tag = int(spliced_with_xs > 0)
        if spliced_sampled > 0 and has_xs_tag:
            xs_tag_status = "detected"
        elif spliced_sampled > 0:
            xs_tag_status = "missing"
            warnings.append("xs_tag_missing_in_sampled_spliced_reads")
        elif sampled > 0:
            xs_tag_status = "not_tested_no_spliced_read_in_sample"
        else:
            xs_tag_status = "not_tested_no_alignment_sample"

    if bam_exists and size_status != "ok":
        reasons.append(f"bam_size_{size_status}")
    is_coordinate = int(sort_order == "coordinate" and (sq_count or 0) > 0)
    if bam_exists and not is_coordinate:
        reasons.append("not_coordinate_sorted_genome_bam")
    chrom_size_check = compare_chrom_sizes(sq_lengths, chrom_sizes)
    if chrom_size_check["bam_chrom_size_match_status"] == "mismatch":
        warnings.append("chrom_sizes_mismatch")
    elif chrom_size_check["bam_chrom_size_match_status"] == "not_tested":
        warnings.append("chrom_sizes_not_tested")
    if not log_path.exists():
        warnings.append("star_log_final_missing")
    if not sj_path.exists():
        warnings.append("sj_out_tab_missing")

    sj_has_rows = isinstance(sj_line_count, int) and sj_line_count > 0
    total_splices = int_or_none(star_metrics["star_total_splices"])
    splice_evidence = has_spliced or sj_has_rows or (total_splices is not None and total_splices > 0)
    if splice_evidence:
        splice_status = "detected"
    else:
        splice_status = "not_detected"
        warnings.append("splice_evidence_not_detected")
    if splice_status == "detected" and xs_tag_status.startswith("not_tested"):
        warnings.append(xs_tag_status)

    needs_index = int(bam_exists and index_path is None)
    content_usable = int(not reasons)
    ready_for_region_query = int(content_usable and not needs_index)
    can_use = content_usable
    recommended_for_stringtie = int(
        bam_exists
        and quickcheck_status == "pass"
        and is_coordinate
        and size_status == "ok"
        and splice_status == "detected"
    )
    recommended_for_stringtie3 = int(recommended_for_stringtie and has_xs_tag)

    row = {
        "animal_id": sample.animal_id,
        "tissue": sample.tissue,
        "sample_tissue": sample.sample_tissue,
        "raw_data_name": sample.raw_data_name,
        "star_sample_dir": sample_dir,
        "align_star_dir": align_dir,
        "raw_star_bam": bam_path,
        "raw_star_bam_exists": bam_exists,
        "raw_star_bam_size_bytes": bam_size if bam_size is not None else "NA",
        "raw_star_bam_size_status": size_status,
        "raw_star_bam_index_exists": int(index_path is not None),
        "raw_star_bam_index_path": index_path if index_path is not None else "NA",
        "star_log_final": log_path,
        "star_log_final_exists": int(log_path.exists()),
        "star_input_reads": star_metrics["star_input_reads"],
        "star_uniquely_mapped_reads_number": star_metrics["star_uniquely_mapped_reads_number"],
        "star_uniquely_mapped_reads_percent": star_metrics["star_uniquely_mapped_reads_percent"],
        "star_multimapping_percent": star_metrics["star_multimapping_percent"],
        "star_unmapped_too_short_percent": star_metrics["star_unmapped_too_short_percent"],
        "star_total_splices": star_metrics["star_total_splices"],
        "sj_out_tab": sj_path,
        "sj_out_tab_exists": int(sj_path.exists()),
        "sj_out_line_count": sj_line_count,
        "samtools_quickcheck_status": quickcheck_status,
        "samtools_quickcheck_message": quickcheck_message,
        "bam_sort_order": sort_order,
        "is_genome_coordinate_bam": is_coordinate,
        "sq_contig_count": sq_count if sq_count is not None else "NA",
        "chrom_sizes": chrom_sizes_path,
        "chrom_sizes_contig_count": len(chrom_sizes),
        **chrom_size_check,
        "sampled_alignment_count": sampled,
        "has_spliced_cigar_N_in_sample": has_spliced,
        "spliced_alignment_sample_count": spliced_sampled,
        "spliced_alignment_with_XS_count": spliced_with_xs,
        "has_XS_tag_in_spliced_reads": has_xs_tag,
        "xs_tag_status": xs_tag_status,
        "splice_evidence_status": splice_status,
        "raw_bam_content_usable_for_discovery": content_usable,
        "raw_bam_ready_for_region_query": ready_for_region_query,
        "raw_bam_needs_index": needs_index,
        "can_use_for_novel_transcript_discovery": can_use,
        "recommended_for_stringtie_discovery": recommended_for_stringtie,
        "recommended_for_stringtie3_discovery": recommended_for_stringtie3,
        "cannot_use_reason": "|".join(reasons) if reasons else "NA",
        "warning_reason": "|".join(warnings) if warnings else "NA",
    }
    return SampleReadiness(row=row)


def output_paths(outdir: Path) -> Dict[str, Path]:
    return {
        "readiness": outdir / "bam_discovery_readiness.tsv",
        "summary": outdir / "bam_discovery_readiness_summary.tsv",
        "manifest": outdir / "bam_discovery_readiness_manifest.tsv",
    }


def open_output(path: Path, overwrite: bool):
    if path.exists() and not overwrite:
        raise FileExistsError(f"output file already exists, if confirm overwrite please plus --overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", encoding="utf-8", newline="")


def check_inputs(args: argparse.Namespace) -> None:
    if not args.sample_list.exists():
        raise FileNotFoundError(f"cannot find sample list: {args.sample_list}")
    if not args.star_root.exists():
        raise FileNotFoundError(f"cannot find STAR root directory: {args.star_root}")
    if not args.samtools.exists():
        raise FileNotFoundError(f"cannot find samtools: {args.samtools}")
    if not args.chrom_sizes.exists():
        raise FileNotFoundError(f"cannot find chrom sizes/.fai file: {args.chrom_sizes}")
    if args.workers < 1:
        raise ValueError("--workers must be >= 1")
    if args.max_alignments < 1:
        raise ValueError("--max-alignments must be >= 1")
    if args.min_bam_size_bytes < 0:
        raise ValueError("--min-bam-size-bytes must be >= 0")
    if args.command_timeout_sec < 1:
        raise ValueError("--command-timeout-sec must be >= 1")


def check_outputs(paths: Mapping[str, Path], overwrite: bool) -> None:
    if overwrite:
        return
    existing = [str(path) for path in paths.values() if path.exists()]
    if existing:
        raise FileExistsError("output file already exists, if confirm overwrite please plus --overwrite:\n" + "\n".join(existing))


def write_readiness_table(path: Path, rows: Sequence[Mapping[str, object]], overwrite: bool) -> None:
    with open_output(path, overwrite) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=READINESS_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in READINESS_FIELDS})


def count_rows(rows: Sequence[Mapping[str, object]], field: str, value: object) -> int:
    return sum(1 for row in rows if row.get(field) == value)


def build_summary_rows(samples: Sequence[SampleRecord], rows: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    return [
        {
            "metric": "sample_rows",
            "value": len(samples),
            "note": "Number of animal_id × tissue rows in the sample manifest.",
        },
        {
            "metric": "raw_star_bam_exists_count",
            "value": count_rows(rows, "raw_star_bam_exists", 1),
            "note": "Number of samples with Aligned.sortedByCoord.out.bam located.",
        },
        {
            "metric": "raw_star_bam_index_exists_count",
            "value": count_rows(rows, "raw_star_bam_index_exists", 1),
            "note": "Number of samples with a .bai/.csi index available for region queries.",
        },
        {
            "metric": "samtools_quickcheck_pass_count",
            "value": count_rows(rows, "samtools_quickcheck_status", "pass"),
            "note": "Number of BAM files passing samtools quickcheck.",
        },
        {
            "metric": "coordinate_sorted_bam_count",
            "value": count_rows(rows, "is_genome_coordinate_bam", 1),
            "note": "Number of BAM files whose headers indicate coordinate sorting and include @SQ contigs.",
        },
        {
            "metric": "sj_out_tab_exists_count",
            "value": count_rows(rows, "sj_out_tab_exists", 1),
            "note": "Number of samples with a STAR SJ.out.tab file.",
        },
        {
            "metric": "bam_chrom_size_match_count",
            "value": count_rows(rows, "bam_chrom_size_match_status", "match"),
            "note": "Number of samples whose BAM @SQ contig names and lengths match chrom sizes/.fai.",
        },
        {
            "metric": "bam_chrom_size_mismatch_count",
            "value": count_rows(rows, "bam_chrom_size_match_status", "mismatch"),
            "note": "Number of samples with BAM @SQ contig names or lengths inconsistent with chrom sizes/.fai.",
        },
        {
            "metric": "splice_evidence_detected_count",
            "value": count_rows(rows, "splice_evidence_status", "detected"),
            "note": "Number of samples with splice evidence from at least one of CIGAR N, SJ.out.tab, or STAR total splices.",
        },
        {
            "metric": "xs_tag_detected_count",
            "value": count_rows(rows, "has_XS_tag_in_spliced_reads", 1),
            "note": "Number of samples with XS tags detected in sampled spliced reads.",
        },
        {
            "metric": "can_use_for_novel_transcript_discovery_count",
            "value": count_rows(rows, "can_use_for_novel_transcript_discovery", 1),
            "note": "Number of samples whose BAM content is suitable for downstream coverage/junction/StringTie/Scallop analysis.",
        },
        {
            "metric": "recommended_for_stringtie_discovery_count",
            "value": count_rows(rows, "recommended_for_stringtie_discovery", 1),
            "note": "Number of samples with an existing BAM, passing quickcheck, coordinate sorting, matching sizes, and detected splice evidence.",
        },
        {
            "metric": "recommended_for_stringtie3_discovery_count",
            "value": count_rows(rows, "recommended_for_stringtie3_discovery", 1),
            "note": "Number of samples meeting the recommended StringTie conditions and additionally showing XS tags in sampled spliced reads.",
        },
        {
            "metric": "raw_bam_ready_for_region_query_count",
            "value": count_rows(rows, "raw_bam_ready_for_region_query", 1),
            "note": "Number of samples with usable BAM content and an existing index.",
        },
        {
            "metric": "raw_bam_needs_index_count",
            "value": count_rows(rows, "raw_bam_needs_index", 1),
            "note": "Number of samples with an existing BAM but no index.",
        },
        {
            "metric": "unusable_sample_count",
            "value": len(rows) - count_rows(rows, "can_use_for_novel_transcript_discovery", 1),
            "note": "Number of unusable samples due to missing BAM files or quickcheck/header/coordinate/size issues.",
        },
    ]


def write_summary_table(path: Path, rows: Sequence[Mapping[str, object]], overwrite: bool) -> None:
    with open_output(path, overwrite) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=SUMMARY_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_manifest(path: Path, args: argparse.Namespace, paths: Mapping[str, Path], overwrite: bool) -> None:
    row = {
        "sample_list": args.sample_list,
        "star_root": args.star_root,
        "output_root": args.outdir,
        "samtools": args.samtools,
        "chrom_sizes": args.chrom_sizes,
        "max_alignments": args.max_alignments,
        "workers": args.workers,
        "min_bam_size_bytes": args.min_bam_size_bytes,
        "command_timeout_sec": args.command_timeout_sec,
        "readiness_table": paths["readiness"],
        "summary_table": paths["summary"],
        "manifest_table": paths["manifest"],
    }
    with open_output(path, overwrite) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=MANIFEST_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    configure_logging()
    args = parse_args()
    check_inputs(args)
    paths = output_paths(args.outdir)
    check_outputs(paths, args.overwrite)

    logging.info("START five-dogs bulk STAR BAM discovery readiness")
    logging.info("sample_list=%s", args.sample_list)
    logging.info("star_root=%s", args.star_root)
    logging.info("outdir=%s", args.outdir)
    logging.info("samtools=%s", args.samtools)
    logging.info("chrom_sizes=%s", args.chrom_sizes)
    logging.info("Stage 1/4: loading sample list and chrom sizes")
    samples = read_sample_list(args.sample_list)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    logging.info("Stage 1/4 done: samples=%s chrom_size_contigs=%s", len(samples), len(chrom_sizes))

    tasks = [
        (
            sample,
            args.star_root,
            args.samtools,
            args.chrom_sizes,
            chrom_sizes,
            args.max_alignments,
            args.min_bam_size_bytes,
            args.command_timeout_sec,
        )
        for sample in samples
    ]
    workers = min(args.workers, len(tasks))
    logging.info("Stage 2/4: inspecting STAR BAM files with workers=%s", workers)
    rows_by_sample: Dict[str, Dict[str, object]] = {}
    completed = 0
    if workers == 1:
        for task in tasks:
            readiness = inspect_sample(task)
            rows_by_sample[str(readiness.row["sample_tissue"])] = readiness.row
            completed += 1
            logging.info("checked %s/%s: %s", completed, len(tasks), readiness.row["sample_tissue"])
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_sample = {executor.submit(inspect_sample, task): task[0].sample_tissue for task in tasks}
            for future in as_completed(future_to_sample):
                sample_tissue = future_to_sample[future]
                try:
                    readiness = future.result()
                except Exception as exc:
                    raise RuntimeError(f"STAR BAM readiness check failed: {sample_tissue}") from exc
                rows_by_sample[str(readiness.row["sample_tissue"])] = readiness.row
                completed += 1
                logging.info("checked %s/%s: %s", completed, len(tasks), readiness.row["sample_tissue"])
    rows = [rows_by_sample[sample.sample_tissue] for sample in samples]
    logging.info("Stage 2/4 done: inspected=%s", len(rows))

    logging.info("Stage 3/4: building summary")
    summary_rows = build_summary_rows(samples, rows)
    logging.info("Stage 3/4 done")

    logging.info("Stage 4/4: writing outputs")
    write_readiness_table(paths["readiness"], rows, args.overwrite)
    write_summary_table(paths["summary"], summary_rows, args.overwrite)
    write_manifest(paths["manifest"], args, paths, args.overwrite)
    logging.info(
        "FINISH samples=%s usable=%s region_ready=%s stringtie_recommended=%s stringtie3_recommended=%s outdir=%s",
        len(rows),
        count_rows(rows, "can_use_for_novel_transcript_discovery", 1),
        count_rows(rows, "raw_bam_ready_for_region_query", 1),
        count_rows(rows, "recommended_for_stringtie_discovery", 1),
        count_rows(rows, "recommended_for_stringtie3_discovery", 1),
        args.outdir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
