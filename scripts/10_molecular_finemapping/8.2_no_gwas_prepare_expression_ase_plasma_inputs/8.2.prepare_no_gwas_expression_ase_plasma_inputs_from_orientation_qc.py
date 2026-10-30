#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3


from __future__ import annotations

import argparse
import csv
import bisect
import math
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


ORIENTATION_BLOCK_QC = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.1_no_gwas_plasma_as_only/20260821_165620/orientation_block_qc.tsv'
)
ORIENTATION_COMPLETION_MARKER = ORIENTATION_BLOCK_QC.parent / "execution_completion.marker"
TASK_FILE = Path(
    '/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv'
)
PHASE_BCF_ROOT = Path(
    '/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.10_reconcile_dog11_whatsHap_shapeit5/reconciled_bcf'
)
BCFTOOLS_BIN = Path(
    "/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"
)


EXPRESSION_MATRIX = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.2.4.normalized_expression/20260821_142351/bulk_rna.cnv_expression.manual_new_gene.log2_sizefactor_norm_counts.tsv'
)
GENE_ANNOTATION = Path(
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.2.4.normalized_expression/20260821_142351/bulk_rna.cnv_expression.manual_new_gene.gene_annotation.tsv'
)
EXPRESSION_COMPLETION_MARKER = EXPRESSION_MATRIX.parent / "execution_completion.marker"

CIS_WINDOW_BP = 100000
MIN_TOTAL_COUNT = 2
MIN_SELECTED_SAMPLE_COUNT = 2
MIN_COMMON_MARKERS = 2


EXPRESSION_SOURCE = "manual_new_gene_log2_sizefactor_norm_counts"
CHROMOSOMES = tuple("chr{0}".format(index) for index in range(1, 39))


@dataclass(frozen=True)
class SampleBlock:
    source_block_id: str
    sample: str
    tissue: str
    feature_id: str
    feature_primary_id: str
    feature_label: str
    feature_class: str
    chrom: str
    phase_set: str
    variant_count: int
    hap1_count: int
    hap2_count: int
    total_count: int
    min_pos: int
    max_pos: int
    orientation_status: str
    orientation_evidence_summary: str


@dataclass(frozen=True)
class Marker:
    chrom: str
    pos: int
    ref: str
    alt: str
    hap_a: Tuple[Optional[int], ...]
    hap_b: Tuple[Optional[int], ...]
    local_ps: Tuple[str, ...]


@dataclass(frozen=True)
class GeneAnnotation:
    gene_id: str
    gene_name: str
    seqname: str
    start: int
    end: int
    strand: str


def log(message: str) -> None:
    stamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S ]")
    print("{0} {1}".format(stamp, message), flush=True)


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(clean(value))
    except (TypeError, ValueError):
        return default


def safe_float(value: object) -> Optional[float]:
    try:
        parsed = float(clean(value))
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def format_float(value: float) -> str:
    return "{0:.8g}".format(value)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", clean(value))
    return cleaned.strip("_") or "NA"


def feature_class_for(feature_id: str) -> str:
    text = clean(feature_id)
    if not text or text in {"NA", "."} or text.startswith("(NA)"):
        return "unannotated_like"
    return "gene_like"


def primary_feature_id(feature_id: str) -> str:
    return clean(feature_id)


def sample_column(block: SampleBlock) -> str:
    return "{0}_{1}".format(block.sample, block.tissue)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare five-dog no-GWAS PLASMA expression+ASE joint inputs."
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Pre-created, empty timestamped output directory.",
    )
    return parser.parse_args()


def require_success_marker(path: Path, source_label: str) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(
            "{0} completion marker is missing or empty: {1}".format(source_label, path)
        )
    normalized_lines = {
        line.strip().replace("\t", "=", 1)
        for line in path.read_text().splitlines()
        if line.strip()
    }
    if "execution_status=success" not in normalized_lines:
        raise ValueError(
            "{0} completion marker does not report success: {1}".format(source_label, path)
        )


def parse_gt_to_haps(gt: str, alt: str) -> Optional[Tuple[int, int]]:
    text = clean(gt)
    if "|" not in text or "," in clean(alt):
        return None
    left, right = text.split("|", 1)
    if left not in {"0", "1"} or right not in {"0", "1"}:
        return None
    return int(left), int(right)


def normalize_local_ps(value: object) -> str:
    text = clean(value)
    return "NA" if text in {"", ".", "NA"} else text


def phase_bcf_path(chrom: str) -> Path:
    return PHASE_BCF_ROOT / "{0}.dog11.ref1519.whatsHap_reconciled.phased.bcf".format(chrom)


def group_key(block: SampleBlock) -> Tuple[str, str, str]:
    return (block.tissue, block.feature_id, block.chrom)


def no_gwas_block_id(tissue: str, feature_id: str, chrom: str) -> str:
    return "no_gwas__{0}__{1}__{2}".format(
        safe_name(tissue),
        safe_name(feature_id),
        safe_name(chrom),
    )


def load_task_manifest(path: Path) -> Tuple[List[str], List[str]]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError("five-dog task manifest is missing or empty: {0}".format(path))

    animals: List[str] = []
    sample_columns: List[str] = []
    seen_animals: Set[str] = set()
    seen_columns: Set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"animal_id", "term", "raw_data_name"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError("five-dog task manifest missing required columns: {0}".format(",".join(missing)))
        for line_number, row in enumerate(reader, 2):
            animal = clean(row.get("animal_id"))
            tissue = clean(row.get("term"))
            raw_data_name = clean(row.get("raw_data_name"))
            if not animal or not tissue or not raw_data_name:
                raise ValueError("five-dog task manifest has empty required value at line {0}".format(line_number))
            column = "{0}_{1}".format(animal, tissue)
            if column in seen_columns:
                raise ValueError("Five-dog task manifest contains a duplicate animal+tissue entry: {0}".format(column))
            seen_columns.add(column)
            sample_columns.append(column)
            if animal not in seen_animals:
                seen_animals.add(animal)
                animals.append(animal)
    if len(animals) != 5 or len(sample_columns) != 19:
        raise ValueError(
            "five-dog task manifest expected 5 animals and 19 animal+tissue rows; found {0} and {1}".format(
                len(animals), len(sample_columns)
            )
        )
    return animals, sample_columns


def load_gene_annotation(path: Path) -> Dict[str, GeneAnnotation]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError("gene annotation file is missing or empty: {0}".format(path))

    annotations: Dict[str, GeneAnnotation] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_id", "gene_name", "seqname", "start", "end", "strand"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError("gene annotation missing required columns: {0}".format(",".join(missing)))
        for line_number, row in enumerate(reader, 2):
            gene_id = clean(row.get("gene_id"))
            gene_name = clean(row.get("gene_name"))
            seqname = clean(row.get("seqname"))
            start = safe_int(row.get("start"))
            end = safe_int(row.get("end"))
            strand = clean(row.get("strand"))
            if not gene_id or not gene_name or not seqname or start <= 0 or end < start or strand not in {"+", "-"}:
                raise ValueError("invalid gene annotation at line {0}: {1}".format(line_number, gene_id or "NA"))
            if gene_id in annotations:
                raise ValueError("duplicate gene_id in gene annotation: {0}".format(gene_id))
            annotations[gene_id] = GeneAnnotation(
                gene_id=gene_id,
                gene_name=gene_name,
                seqname=seqname,
                start=start,
                end=end,
                strand=strand,
            )
    return annotations


def load_expression_matrix(path: Path) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError("expression matrix is missing or empty: {0}".format(path))

    expression: Dict[str, Dict[str, float]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("expression matrix must contain gene_id column: {0}".format(path))
        sample_columns = [field for field in reader.fieldnames if field != "gene_id"]
        if not sample_columns:
            raise ValueError("expression matrix has no sample columns: {0}".format(path))
        for row in reader:
            gene_id = clean(row.get("gene_id"))
            if not gene_id:
                continue
            values: Dict[str, float] = {}
            for column in sample_columns:
                value = safe_float(row.get(column))
                if value is not None:
                    values[column] = value
            expression[gene_id] = values
    return expression, sample_columns


def resolve_expression_gene_id(feature_id: str, expression: Dict[str, Dict[str, float]]) -> str:
    gene_id = clean(feature_id)
    return gene_id if gene_id in expression else ""


def read_orientation_blocks(
    path: Path,
    allowed_sample_columns: Set[str],
) -> Tuple[List[SampleBlock], List[Dict[str, object]]]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError("orientation block QC file is missing or empty: {0}".format(path))

    blocks: List[SampleBlock] = []
    inventory_rows: List[Dict[str, object]] = []
    required = {
        "block_id",
        "sample",
        "tissue",
        "feature_id",
        "feature_label",
        "gene_id",
        "chrom",
        "phase_set",
        "variant_count",
        "checked_variant_count",
        "hap1_count",
        "hap2_count",
        "total_count",
        "min_pos",
        "max_pos",
        "orientation_status",
        "orientation_evidence_summary",
        "aggregation_status",
    }

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError("orientation_block_qc.tsv missing required columns: {0}".format(",".join(missing)))

        for row_index, row in enumerate(reader, 1):
            feature_id = clean(row.get("feature_id"))
            feature_class = feature_class_for(feature_id)
            hap1_count = safe_int(row.get("hap1_count"))
            hap2_count = safe_int(row.get("hap2_count"))
            total_count = safe_int(row.get("total_count"))
            variant_count = safe_int(row.get("variant_count"))
            checked_variant_count = safe_int(row.get("checked_variant_count"))
            chrom = clean(row.get("chrom"))
            phase_set = clean(row.get("phase_set"))
            min_pos = safe_int(row.get("min_pos"))
            max_pos = safe_int(row.get("max_pos"))
            orientation_status = clean(row.get("orientation_status"))
            aggregation_status = clean(row.get("aggregation_status"))
            sample = clean(row.get("sample"))
            tissue = clean(row.get("tissue"))

            exclude_reasons: List[str] = []
            if feature_id != clean(row.get("gene_id")):
                exclude_reasons.append("feature_id_gene_id_mismatch")
            if orientation_status != "consistent":
                exclude_reasons.append("orientation_not_consistent")
            if aggregation_status != "matched_1.5":
                exclude_reasons.append("aggregation_not_matched_1.5")
            if checked_variant_count != variant_count:
                exclude_reasons.append("not_all_variants_checked")
            if hap1_count <= 0 or hap2_count <= 0:
                exclude_reasons.append("haplotype_zero_count")
            if total_count < MIN_TOTAL_COUNT:
                exclude_reasons.append("total_below_min")
            if not chrom or chrom == "mixed":
                exclude_reasons.append("invalid_chrom")
            if not phase_set:
                exclude_reasons.append("missing_phase_set")
            elif phase_set != chrom:
                exclude_reasons.append("phase_scope_not_chromosome")
            if min_pos <= 0 or max_pos <= 0 or min_pos > max_pos:
                exclude_reasons.append("invalid_position_span")
            if total_count != hap1_count + hap2_count:
                exclude_reasons.append("haplotype_total_count_mismatch")
            if feature_class == "unannotated_like":
                exclude_reasons.append("unannotated_feature_not_eligible_for_joint_input")
            if "{0}_{1}".format(sample, tissue) not in allowed_sample_columns:
                exclude_reasons.append("sample_tissue_not_in_five_dog_manifest")

            sample_allowed = 0 if exclude_reasons else 1
            feature_primary_id = primary_feature_id(feature_id)
            feature_label = clean(row.get("feature_label")) or feature_id
            inventory_rows.append(
                {
                    "sample": sample,
                    "tissue": tissue,
                    "feature_id": feature_id,
                    "feature_primary_id": feature_primary_id,
                    "feature_label": feature_label,
                    "feature_class": feature_class,
                    "chrom": chrom,
                    "phase_scope": phase_set,
                    "hap1_count": hap1_count,
                    "hap2_count": hap2_count,
                    "total_count": total_count,
                    "variant_count": variant_count,
                    "checked_variant_count": checked_variant_count,
                    "min_pos": min_pos if min_pos > 0 else "",
                    "max_pos": max_pos if max_pos > 0 else "",
                    "region_source": "ase_variant_span",
                    "count_source": "gene_id_chromosome_aggregated",
                    "orientation_status": orientation_status,
                    "aggregation_status": aggregation_status,
                    "plasma_sample_block_allowed": sample_allowed,
                    "exclude_reason": ";".join(exclude_reasons),
                    "source_orientation_block_id": clean(row.get("block_id")) or "row_{0}".format(row_index),
                }
            )
            if exclude_reasons:
                continue

            blocks.append(
                SampleBlock(
                    source_block_id=clean(row.get("block_id")) or "row_{0}".format(row_index),
                    sample=sample,
                    tissue=tissue,
                    feature_id=feature_id,
                    feature_primary_id=feature_primary_id,
                    feature_label=feature_label,
                    feature_class=feature_class,
                    chrom=chrom,
                    phase_set=phase_set,
                    variant_count=variant_count,
                    hap1_count=hap1_count,
                    hap2_count=hap2_count,
                    total_count=total_count,
                    min_pos=min_pos,
                    max_pos=max_pos,
                    orientation_status=orientation_status,
                    orientation_evidence_summary=clean(row.get("orientation_evidence_summary")),
                )
            )
    return blocks, inventory_rows


def choose_unique_sample_blocks(blocks: Sequence[SampleBlock]) -> Tuple[List[SampleBlock], str]:
    by_sample: Dict[str, SampleBlock] = {}
    for block in sorted(blocks, key=lambda item: (item.sample, item.source_block_id)):
        previous = by_sample.get(block.sample)
        if previous is None:
            by_sample[block.sample] = block
            continue
        raise ValueError(
            "duplicate 8.1 block for sample+tissue+GeneID+chrom: {0}; {1}".format(
                previous.source_block_id, block.source_block_id
            )
        )
    return [by_sample[sample] for sample in sorted(by_sample)], ""


def build_candidate_groups(blocks: Sequence[SampleBlock]) -> List[Tuple[str, Tuple[str, str, str], List[SampleBlock], str]]:
    grouped: Dict[Tuple[str, str, str], List[SampleBlock]] = defaultdict(list)
    for block in blocks:
        grouped[group_key(block)].append(block)

    candidates: List[Tuple[str, Tuple[str, str, str], List[SampleBlock], str]] = []
    for key, group_blocks in sorted(grouped.items()):
        unique_blocks, dedup_note = choose_unique_sample_blocks(group_blocks)
        if len(unique_blocks) < MIN_SELECTED_SAMPLE_COUNT:
            continue
        tissue, feature_id, chrom = key
        block_id = no_gwas_block_id(tissue, feature_id, chrom)
        candidates.append((block_id, key, unique_blocks, dedup_note))
    return candidates


def candidate_window(
    blocks: Sequence[SampleBlock],
    annotation: GeneAnnotation,
) -> Tuple[int, int, int, int, int, int]:
    ase_min_pos = min(block.min_pos for block in blocks)
    ase_max_pos = max(block.max_pos for block in blocks)
    cis_start = max(1, annotation.start - CIS_WINDOW_BP)
    cis_end = annotation.end + CIS_WINDOW_BP
    return annotation.start, annotation.end, ase_min_pos, ase_max_pos, cis_start, cis_end


def validate_phase_bcf_inputs(animals: Sequence[str]) -> Dict[str, Path]:
    if not BCFTOOLS_BIN.exists() or not BCFTOOLS_BIN.is_file():
        raise FileNotFoundError("bcftools is missing: {0}".format(BCFTOOLS_BIN))

    paths: Dict[str, Path] = {}
    expected_samples = set(animals)
    for chrom in CHROMOSOMES:
        path = phase_bcf_path(chrom)
        index_path = Path("{0}.csi".format(path))
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError("final phased BCF is missing or empty: {0}".format(path))
        if not index_path.exists() or index_path.stat().st_size == 0:
            raise FileNotFoundError("final phased BCF index is missing or empty: {0}".format(index_path))
        result = subprocess.run(
            [str(BCFTOOLS_BIN), "query", "-l", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        bcf_samples = {line.strip() for line in result.stdout.splitlines() if line.strip()}
        missing_samples = sorted(expected_samples - bcf_samples)
        if missing_samples:
            raise ValueError("{0} missing five-dog samples: {1}".format(path, ",".join(missing_samples)))
        paths[chrom] = path
    return paths


def scan_bcf_markers(path: Path, chrom: str, animals: Sequence[str]) -> Tuple[List[Marker], List[int]]:
    query_format = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%PS]\\n"
    command = [
        str(BCFTOOLS_BIN),
        "query",
        "-s",
        ",".join(animals),
        "-f",
        query_format,
        str(path),
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if process.stdout is None:
        raise RuntimeError("failed to capture bcftools stdout: {0}".format(path))

    markers: List[Marker] = []
    expected_fields = 4 + 2 * len(animals)
    previous_pos = 0
    seen_marker_keys: Set[Tuple[int, str, str]] = set()
    for line_number, line in enumerate(process.stdout, 1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != expected_fields:
            process.kill()
            raise ValueError(
                "unexpected bcftools query field count at {0}: line {1}; found {2}, expected {3}".format(
                    path, line_number, len(fields), expected_fields
                )
            )
        record_chrom, pos_text, ref, alt = fields[:4]
        if record_chrom != chrom:
            process.kill()
            raise ValueError("unexpected chromosome in {0}: {1}".format(path, record_chrom))
        if "," in alt or len(ref) != 1 or len(alt) != 1:
            continue
        pos = safe_int(pos_text)
        if pos <= 0:
            process.kill()
            raise ValueError("invalid BCF position at {0}: line {1}".format(path, line_number))
        marker_key = (pos, ref, alt)
        if pos < previous_pos or marker_key in seen_marker_keys:
            process.kill()
            raise ValueError("BCF markers are not position-sorted with unique alleles: {0}; {1}".format(path, marker_key))
        previous_pos = pos
        seen_marker_keys.add(marker_key)

        hap_a_values: List[Optional[int]] = []
        hap_b_values: List[Optional[int]] = []
        local_ps_values: List[str] = []
        for sample_index in range(len(animals)):
            gt = fields[4 + 2 * sample_index]
            local_ps = fields[5 + 2 * sample_index]
            hap_pair = parse_gt_to_haps(gt, alt)
            hap_a_values.append(hap_pair[0] if hap_pair is not None else None)
            hap_b_values.append(hap_pair[1] if hap_pair is not None else None)
            local_ps_values.append(normalize_local_ps(local_ps))
        markers.append(
            Marker(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                hap_a=tuple(hap_a_values),
                hap_b=tuple(hap_b_values),
                local_ps=tuple(local_ps_values),
            )
        )

    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError("bcftools query failed: {0}; rc={1}; {2}".format(path, return_code, stderr.strip()))
    positions = [marker.pos for marker in markers]
    return markers, positions


def common_informative_markers(
    blocks: Sequence[SampleBlock],
    markers: Sequence[Marker],
    marker_positions: Sequence[int],
    sample_indices: Dict[str, int],
    cis_start: int,
    cis_end: int,
) -> Tuple[List[Marker], str]:
    left = bisect.bisect_left(marker_positions, cis_start)
    right = bisect.bisect_right(marker_positions, cis_end)
    selected_indices = [sample_indices[block.sample] for block in blocks]
    kept: List[Marker] = []
    missing_genotype_count = 0
    monomorphic_count = 0
    for marker in markers[left:right]:
        haplotypes: List[int] = []
        complete = True
        for sample_index in selected_indices:
            hap_a = marker.hap_a[sample_index]
            hap_b = marker.hap_b[sample_index]
            if hap_a is None or hap_b is None:
                complete = False
                break
            haplotypes.extend([hap_a, hap_b])
        if not complete:
            missing_genotype_count += 1
            continue
        alt_count = sum(haplotypes)
        if alt_count == 0 or alt_count == len(haplotypes):
            monomorphic_count += 1
            continue
        kept.append(marker)

    if not kept:
        return [], "no_common_informative_markers"
    details: List[str] = []
    if missing_genotype_count:
        details.append("missing_or_unphased_gt_excluded={0}".format(missing_genotype_count))
    if monomorphic_count:
        details.append("monomorphic_in_selected_dogs_excluded={0}".format(monomorphic_count))
    status = "consistent" if not details else "consistent_with_" + ";".join(details)
    return kept, status


def write_plasma_input(
    block_id: str,
    blocks: Sequence[SampleBlock],
    markers: Sequence[Marker],
    sample_indices: Dict[str, int],
    expression_gene_id: str,
    total_exp_by_sample: Dict[str, float],
    out_root: Path,
) -> Path:
    out_dir = out_root / "plasma_inputs" / block_id
    if out_dir.exists():
        raise FileExistsError("PLASMA block output already exists; refusing to overwrite: {0}".format(out_dir))
    out_dir.mkdir(parents=True, exist_ok=True)

    hap_a_rows: List[str] = []
    hap_b_rows: List[str] = []
    counts_a_rows: List[str] = []
    counts_b_rows: List[str] = []
    total_exp_rows: List[str] = []
    sample_rows: List[Dict[str, object]] = []
    marker_rows: List[Dict[str, object]] = []

    for block in blocks:
        sample_index = sample_indices[block.sample]
        hap_a_rows.append(" ".join(str(marker.hap_a[sample_index]) for marker in markers))
        hap_b_rows.append(" ".join(str(marker.hap_b[sample_index]) for marker in markers))
        counts_a_rows.append(str(block.hap1_count))
        counts_b_rows.append(str(block.hap2_count))
        total_exp = total_exp_by_sample[block.sample]
        total_exp_rows.append(format_float(total_exp))
        sample_rows.append(
            {
                "sample": block.sample,
                "tissue": block.tissue,
                "feature_id": block.feature_id,
                "feature_primary_id": block.feature_primary_id,
                "feature_label": block.feature_label,
                "feature_class": block.feature_class,
                "expression_gene_id": expression_gene_id,
                "expression_source": EXPRESSION_SOURCE,
                "total_exp": format_float(total_exp),
                "chrom": block.chrom,
                "phase_scope": block.phase_set,
                "count_source": "gene_id_chromosome_aggregated",
                "orientation_status": block.orientation_status,
                "orientation_evidence_summary": block.orientation_evidence_summary,
                "hap1_count": block.hap1_count,
                "hap2_count": block.hap2_count,
                "total_count": block.total_count,
                "source_orientation_block_id": block.source_block_id,
            }
        )

    for marker in markers:
        local_ps_by_sample = ",".join(
            "{0}:{1}".format(block.sample, marker.local_ps[sample_indices[block.sample]]) for block in blocks
        )
        marker_rows.append(
            {
                "chrom": marker.chrom,
                "pos": marker.pos,
                "marker_id": "{0}:{1}_{2}>{3}".format(marker.chrom, marker.pos, marker.ref, marker.alt),
                "ref": marker.ref,
                "alt": marker.alt,
                "phase_scope": marker.chrom,
                "local_ps_by_sample": local_ps_by_sample,
            }
        )

    (out_dir / "hap_A.txt").write_text("\n".join(hap_a_rows) + "\n")
    (out_dir / "hap_B.txt").write_text("\n".join(hap_b_rows) + "\n")
    (out_dir / "counts_A.txt").write_text("\n".join(counts_a_rows) + "\n")
    (out_dir / "counts_B.txt").write_text("\n".join(counts_b_rows) + "\n")
    (out_dir / "total_exp.txt").write_text("\n".join(total_exp_rows) + "\n")
    write_tsv(
        out_dir / "samples.tsv",
        [
            "sample",
            "tissue",
            "feature_id",
            "feature_primary_id",
            "feature_label",
            "feature_class",
            "expression_gene_id",
            "expression_source",
            "total_exp",
            "chrom",
            "phase_scope",
            "count_source",
            "orientation_status",
            "orientation_evidence_summary",
            "hap1_count",
            "hap2_count",
            "total_count",
            "source_orientation_block_id",
        ],
        sample_rows,
    )
    write_tsv(
        out_dir / "markers.tsv",
        [
            "chrom",
            "pos",
            "marker_id",
            "ref",
            "alt",
            "phase_scope",
            "local_ps_by_sample",
        ],
        marker_rows,
    )
    (out_dir / "README.txt").write_text(
        "This is a non-GWAS PLASMA expression + ASE joint input.\n"
        "Gene/feature block is the molecular phenotype unit; phase_scope stores the chromosome-wide scope as chrN.\n"
        "Variant/SNP is the fine-mapping target.\n"
        "The result is for candidate cis-regulatory variant catalogue.\n"
        "It does not imply phenotype causal variants or causal genes.\n"
        "counts_A/counts_B are written from HAP1/HAP2 only when orientation_status=consistent.\n"
        "PLASMA does not require cross-sample absolute HAP1/HAP2 direction; each sample row only needs hap_A/hap_B aligned to counts_A/counts_B.\n"
        "Local BCF PS is provenance only and is not used to restrict chromosome-wide phased markers.\n"
        "total_exp.txt is written from the normalized expression matrix and matches samples.tsv row order.\n"
        "Run PLASMA without --as_only and with --total_exp_path total_exp.txt.\n"
    )
    return out_dir


def build_outputs(out_dir: Path) -> None:
    log("[start] prepare no-GWAS PLASMA expression + ASE joint inputs")
    log("[config] ORIENTATION_BLOCK_QC={0}".format(ORIENTATION_BLOCK_QC))
    log("[config] TASK_FILE={0}".format(TASK_FILE))
    log("[config] PHASE_BCF_ROOT={0}".format(PHASE_BCF_ROOT))
    log("[config] BCFTOOLS_BIN={0}".format(BCFTOOLS_BIN))
    log("[config] EXPRESSION_MATRIX={0}".format(EXPRESSION_MATRIX))
    log("[config] GENE_ANNOTATION={0}".format(GENE_ANNOTATION))
    log("[config] OUT_DIR={0}".format(out_dir))
    log("[config] CIS_WINDOW_BP={0}".format(CIS_WINDOW_BP))
    log("[config] MIN_TOTAL_COUNT={0}".format(MIN_TOTAL_COUNT))
    log("[config] MIN_SELECTED_SAMPLE_COUNT={0}".format(MIN_SELECTED_SAMPLE_COUNT))
    log("[config] MIN_COMMON_MARKERS={0}".format(MIN_COMMON_MARKERS))
    log("[config] unannotated ASE features are inventory-only and excluded from joint input")
    log("[config] EXPRESSION_SOURCE={0}".format(EXPRESSION_SOURCE))

    if not out_dir.is_absolute():
        raise ValueError("--out-dir must be an absolute path: {0}".format(out_dir))
    if not out_dir.exists() or not out_dir.is_dir():
        raise FileNotFoundError("--out-dir must be a pre-created directory: {0}".format(out_dir))
    if any(out_dir.iterdir()):
        raise FileExistsError("--out-dir is not empty; refusing to overwrite: {0}".format(out_dir))

    log("[stage:provenance] validating upstream completion markers")
    require_success_marker(ORIENTATION_COMPLETION_MARKER, "8.1 orientation QC")
    require_success_marker(EXPRESSION_COMPLETION_MARKER, "1.2.4 normalized expression")

    log("[stage:manifest] loading five-dog task manifest")
    animals, expected_expression_columns = load_task_manifest(TASK_FILE)
    allowed_sample_columns = set(expected_expression_columns)
    sample_indices = {sample: index for index, sample in enumerate(animals)}
    log("[stage:manifest] animals={0}; animal+tissue rows={1}".format(len(animals), len(expected_expression_columns)))

    log("[stage:expression] loading exact GeneID annotation and expression matrix")
    gene_annotations = load_gene_annotation(GENE_ANNOTATION)
    expression_matrix, expression_sample_columns = load_expression_matrix(EXPRESSION_MATRIX)
    if expression_sample_columns != expected_expression_columns:
        missing_columns = sorted(set(expected_expression_columns) - set(expression_sample_columns))
        extra_columns = sorted(set(expression_sample_columns) - set(expected_expression_columns))
        raise ValueError(
            "expression columns do not exactly match the five-dog task manifest; missing={0}; extra={1}".format(
                ",".join(missing_columns) or "none", ",".join(extra_columns) or "none"
            )
        )
    annotation_gene_ids = set(gene_annotations)
    expression_gene_ids = set(expression_matrix)
    if annotation_gene_ids != expression_gene_ids:
        raise ValueError(
            "gene annotation and expression GeneID sets differ; annotation_only={0}; expression_only={1}".format(
                len(annotation_gene_ids - expression_gene_ids), len(expression_gene_ids - annotation_gene_ids)
            )
        )
    log(
        "[stage:expression] expression genes={0}; annotation genes={1}; sample columns={2}".format(
            len(expression_matrix),
            len(gene_annotations),
            len(expression_sample_columns),
        )
    )

    log("[stage:BCF-check] validating 38 indexed final multi-sample BCF files")
    bcf_paths = validate_phase_bcf_inputs(animals)
    log("[stage:BCF-check] valid BCF files={0}; target animals={1}".format(len(bcf_paths), len(animals)))

    log("[stage:read] loading orientation block QC")
    sample_blocks, inventory_rows = read_orientation_blocks(ORIENTATION_BLOCK_QC, allowed_sample_columns)
    log("[stage:read] allowed sample blocks={0}; inventory rows={1}".format(len(sample_blocks), len(inventory_rows)))

    log("[stage:group] grouping sample blocks by tissue+GeneID+chromosome")
    candidates = build_candidate_groups(sample_blocks)
    unknown_chromosomes = sorted({key[2] for _, key, _, _ in candidates} - set(CHROMOSOMES))
    if unknown_chromosomes:
        raise ValueError("candidate groups contain unsupported chromosomes: {0}".format(",".join(unknown_chromosomes)))
    log("[stage:group] candidate groups with selected_sample_count >= {0}: {1}".format(MIN_SELECTED_SAMPLE_COUNT, len(candidates)))


    plasma_root = out_dir / "plasma_inputs"
    plasma_root.mkdir(parents=True, exist_ok=True)

    inventory_fields = [
        "sample",
        "tissue",
        "feature_id",
        "feature_primary_id",
        "feature_label",
        "feature_class",
        "chrom",
        "phase_scope",
        "hap1_count",
        "hap2_count",
        "total_count",
        "variant_count",
        "checked_variant_count",
        "min_pos",
        "max_pos",
        "region_source",
        "count_source",
        "orientation_status",
        "aggregation_status",
        "plasma_sample_block_allowed",
        "exclude_reason",
        "source_orientation_block_id",
    ]
    log("[stage:write] writing ase_feature_block_inventory.tsv")
    write_tsv(out_dir / "ase_feature_block_inventory.tsv", inventory_fields, inventory_rows)

    manifest_rows: List[Dict[str, object]] = []
    written = 0
    processed = 0
    candidates_by_chrom: Dict[str, List[Tuple[str, Tuple[str, str, str], List[SampleBlock], str]]] = defaultdict(list)
    for candidate in candidates:
        candidates_by_chrom[candidate[1][2]].append(candidate)

    for chrom_index, chrom in enumerate(CHROMOSOMES, 1):
        chrom_candidates = candidates_by_chrom.get(chrom, [])
        if not chrom_candidates:
            continue
        log(
            "[stage:BCF {0}/{1}] reading chromosome-wide phased markers; chrom={2}; candidate_groups={3}".format(
                chrom_index, len(CHROMOSOMES), chrom, len(chrom_candidates)
            )
        )
        chrom_markers, chrom_marker_positions = scan_bcf_markers(bcf_paths[chrom], chrom, animals)
        log("[stage:BCF {0}/{1}] chrom={2}; biallelic_records={3}".format(chrom_index, len(CHROMOSOMES), chrom, len(chrom_markers)))

        for block_id, key, blocks, dedup_note in chrom_candidates:
            processed += 1
            if processed == 1 or processed % 1000 == 0 or processed == len(candidates):
                log("[stage:matrix] candidate groups {0}/{1}; written={2}".format(processed, len(candidates), written))
            tissue, feature_id, _ = key
            annotation = gene_annotations.get(feature_id)
            annotation_status = "ok"
            feature_start: object = ""
            feature_end: object = ""
            feature_strand = ""
            ase_min_pos = min(block.min_pos for block in blocks)
            ase_max_pos = max(block.max_pos for block in blocks)
            cis_start: object = ""
            cis_end: object = ""
            marker_status = "not_scanned_missing_gene_annotation"
            selected_markers: List[Marker] = []
            ase_span_status = "not_checked"

            if annotation is None:
                annotation_status = "gene_annotation_missing"
            elif annotation.seqname != chrom:
                annotation_status = "gene_annotation_chrom_mismatch"
                marker_status = "not_scanned_gene_annotation_chrom_mismatch"
            else:
                feature_start, feature_end, ase_min_pos, ase_max_pos, cis_start, cis_end = candidate_window(
                    blocks, annotation
                )
                feature_strand = annotation.strand
                ase_span_status = (
                    "inside_gene_annotation"
                    if annotation.start <= ase_min_pos <= ase_max_pos <= annotation.end
                    else "outside_gene_annotation"
                )
                selected_markers, marker_status = common_informative_markers(
                    blocks,
                    chrom_markers,
                    chrom_marker_positions,
                    sample_indices,
                    int(cis_start),
                    int(cis_end),
                )

            expression_gene_id = resolve_expression_gene_id(feature_id, expression_matrix)
            total_exp_by_sample: Dict[str, float] = {}
            expression_status = "ok"
            if not expression_gene_id:
                expression_status = "expression_gene_missing"
            else:
                expression_values = expression_matrix[expression_gene_id]
                missing_columns: List[str] = []
                for block in blocks:
                    column = sample_column(block)
                    value = expression_values.get(column)
                    if value is None:
                        missing_columns.append(column)
                    else:
                        total_exp_by_sample[block.sample] = value
                if missing_columns:
                    expression_status = "expression_sample_missing={0}".format(",".join(missing_columns[:10]))

            exclude_reasons: List[str] = []
            if annotation_status != "ok":
                exclude_reasons.append(annotation_status)
            if ase_span_status == "outside_gene_annotation":
                exclude_reasons.append("ase_span_outside_gene_annotation")
            if expression_status != "ok":
                exclude_reasons.append(expression_status)
            if len(selected_markers) < MIN_COMMON_MARKERS:
                exclude_reasons.append("common_informative_marker_count_below_min")

            plasma_allowed = 0 if exclude_reasons else 1
            input_dir = ""
            if plasma_allowed:
                input_dir = str(
                    write_plasma_input(
                        block_id,
                        blocks,
                        selected_markers,
                        sample_indices,
                        expression_gene_id,
                        total_exp_by_sample,
                        out_dir,
                    )
                )
                written += 1

            selected_samples = ",".join(block.sample for block in blocks)
            selected_phase_scopes = ",".join("{0}:{1}".format(block.sample, block.phase_set) for block in blocks)
            total_exp_values = ",".join(
                format_float(total_exp_by_sample[block.sample]) for block in blocks if block.sample in total_exp_by_sample
            )
            manifest_rows.append(
                {
                    "block_id": block_id,
                    "tissue": tissue,
                    "feature_id": feature_id,
                    "feature_primary_id": primary_feature_id(feature_id),
                    "feature_label": blocks[0].feature_label,
                    "feature_class": feature_class_for(feature_id),
                    "expression_gene_id": expression_gene_id,
                    "expression_source": EXPRESSION_SOURCE,
                    "expression_status": expression_status,
                    "total_exp_values": total_exp_values,
                    "chrom": chrom,
                    "feature_start": feature_start,
                    "feature_end": feature_end,
                    "feature_strand": feature_strand,
                    "ase_min_pos": ase_min_pos,
                    "ase_max_pos": ase_max_pos,
                    "ase_span_status": ase_span_status,
                    "cis_window_start": cis_start,
                    "cis_window_end": cis_end,
                    "cis_window_bp": CIS_WINDOW_BP,
                    "region_source": "gene_annotation_gene_span",
                    "annotation_status": annotation_status,
                    "selected_sample_count": len(blocks),
                    "selected_samples": selected_samples,
                    "selected_phase_scopes": selected_phase_scopes,
                    "dedup_note": dedup_note,
                    "common_marker_count": len(selected_markers),
                    "common_marker_positions": ",".join(str(marker.pos) for marker in selected_markers),
                    "total_reads": sum(block.total_count for block in blocks),
                    "min_haplotype_count": min(min(block.hap1_count, block.hap2_count) for block in blocks),
                    "orientation_status": "consistent",
                    "orientation_evidence_count": sum(block.variant_count for block in blocks),
                    "marker_allele_status": marker_status,
                    "matrix_status": "ready" if plasma_allowed else "excluded",
                    "exclude_reason": ";".join(exclude_reasons),
                    "plasma_input_allowed": plasma_allowed,
                    "input_dir": input_dir,
                }
            )

    manifest_rows.sort(
        key=lambda row: (
            -safe_int(row.get("plasma_input_allowed")),
            -safe_int(row.get("selected_sample_count")),
            -safe_int(row.get("common_marker_count")),
            -safe_int(row.get("total_reads")),
            clean(row.get("block_id")),
        )
    )
    manifest_fields = [
        "block_id",
        "tissue",
        "feature_id",
        "feature_primary_id",
        "feature_label",
        "feature_class",
        "expression_gene_id",
        "expression_source",
        "expression_status",
        "total_exp_values",
        "chrom",
        "feature_start",
        "feature_end",
        "feature_strand",
        "ase_min_pos",
        "ase_max_pos",
        "ase_span_status",
        "cis_window_start",
        "cis_window_end",
        "cis_window_bp",
        "region_source",
        "annotation_status",
        "selected_sample_count",
        "selected_samples",
        "selected_phase_scopes",
        "dedup_note",
        "common_marker_count",
        "common_marker_positions",
        "total_reads",
        "min_haplotype_count",
        "orientation_status",
        "orientation_evidence_count",
        "marker_allele_status",
        "matrix_status",
        "exclude_reason",
        "plasma_input_allowed",
        "input_dir",
    ]
    log("[stage:write] writing no_gwas_plasma_candidate_blocks.tsv")
    write_tsv(out_dir / "no_gwas_plasma_candidate_blocks.tsv", manifest_fields, manifest_rows)

    log("[summary] inventory_rows={0}".format(len(inventory_rows)))
    log("[summary] allowed_sample_blocks={0}".format(len(sample_blocks)))
    log("[summary] candidate_groups={0}".format(len(candidates)))
    log("[summary] plasma_inputs_written={0}".format(written))
    log("[summary] matrix_status={0}".format(dict(Counter(row["matrix_status"] for row in manifest_rows))))
    log("[output] {0}".format(out_dir / "ase_feature_block_inventory.tsv"))
    log("[output] {0}".format(out_dir / "no_gwas_plasma_candidate_blocks.tsv"))
    log("[output] {0}".format(plasma_root))
    log("[done] prepare no-GWAS PLASMA expression + ASE joint inputs")


if __name__ == "__main__":
    arguments = parse_args()
    build_outputs(arguments.out_dir)
