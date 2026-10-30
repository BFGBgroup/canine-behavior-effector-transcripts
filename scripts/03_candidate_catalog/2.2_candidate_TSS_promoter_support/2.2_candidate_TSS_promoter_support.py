#!/usr/bin/python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

for thread_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[thread_key] = "1"


PROJECT_ROOT = Path("/path/to/project/analysis_20260709/5-dogs/9.new-gene")
STEP21_ROOT = PROJECT_ROOT / "2.1_transcript_first_preASE_candidate_assessment"
STEP16_ROOT = Path(
    '/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6'
)
REFERENCE_ROOT = Path("/path/to/project/legacy/five-dogs/9.new-gene/reference/TSS")
CANFAM_ROOT = Path("/path/to/project/legacy/reference")

DEFAULT_RETAINED_TSV = STEP21_ROOT / "2.1_retained_transcript_candidates.preASE.tsv"
DEFAULT_SUMMARY_TSV = STEP21_ROOT / "2.1_candidate_transcript_summary.preASE.tsv"
DEFAULT_TRANSCRIPT_EXON_BED = STEP16_ROOT / "transcript_exon.bed"
DEFAULT_CHROM_SIZES = CANFAM_ROOT / "GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"
DEFAULT_DOGA_CTSS_BED = (
    REFERENCE_ROOT
    / "3.DoGA-CanFam4-to-CanFam3-liftover/results/doga_CTSS_brain_CNS.CanFam3.bed"
)
DEFAULT_DOGA_PROMOTER_BED = (
    REFERENCE_ROOT
    / "3.DoGA-CanFam4-to-CanFam3-liftover/results/doga_promoters_all.CanFam3.bed"
)
DEFAULT_BARKBASE_ATAC_BED = (
    REFERENCE_ROOT / "4.BarkBase-ATAC-seq/results/barkbase_ATAC_all_selected.merged.bed"
)
DEFAULT_EPICDOG_PROMOTER_BED = (
    REFERENCE_ROOT
    / "5.EpiC-Dog-chromatin-states/results/epicdog_brain_CNS_promoter_like_states.bed"
)
DEFAULT_EPICDOG_ACTIVE_BED = (
    REFERENCE_ROOT
    / "5.EpiC-Dog-chromatin-states/results/epicdog_brain_CNS_active_H3K4me3_H3K27ac_states.bed"
)
DEFAULT_EPICDOG_SUPER_BED = (
    REFERENCE_ROOT / "5.EpiC-Dog-chromatin-states/results/epicdog_brain_CNS_super_enhancers.bed"
)
DEFAULT_OUTDIR = PROJECT_ROOT / "2.2_candidate_TSS_promoter_support"

PRIMARY_CHROMS = {*(f"chr{index}" for index in range(1, 39)), "chrX"}
RETAINED_REQUIRED_FIELDS = ["query_transcript_id", "preASE_retained_for_downstream"]
SUMMARY_REQUIRED_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "chr",
    "transcript_strand",
    "transcript_span_start0",
    "transcript_span_end0",
    "transcript_span_length",
    "transcript_exonic_length",
    "transcript_exon_count",
    "class_code",
    "novelty_category",
    "FEELnc_biotype",
    "FEELnc_codpot_category",
    "FEELnc_lncRNA_class",
    "candidate_gene_strand",
    "gene_locus_filter_reason",
    "preASE_retained_for_downstream",
]
OUTPUT_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "class_code",
    "novelty_category",
    "FEELnc_biotype",
    "FEELnc_codpot_category",
    "FEELnc_lncRNA_class",
    "candidate_gene_strand",
    "gene_locus_filter_reason",
    "preASE_retained_for_downstream",
    "transcript_chr",
    "transcript_strand",
    "transcript_exon_count",
    "transcript_span_start0",
    "transcript_span_end0",
    "transcript_span_length",
    "transcript_exonic_length",
    "tss_inferred",
    "reference_assessable",
    "not_assessed_reason",
    "tss0",
    "tss_1bp_start0",
    "tss_1bp_end0",
    "tss_pm500_start0",
    "tss_pm500_end0",
    "promoter_1kb_100bp_start0",
    "promoter_1kb_100bp_end0",
    "coordinate_overlap_group",
    "doga_ctss_tss1bp_overlap_count",
    "doga_promoter_tss_pm500_overlap_count",
    "doga_promoter_promoter_overlap_count",
    "barkbase_atac_promoter_overlap_count",
    "epicdog_promoter_like_promoter_overlap_count",
    "epicdog_active_regulatory_promoter_overlap_count",
    "epicdog_super_enhancer_promoter_overlap_count",
    "ctss_exact_supported",
    "ctss_within_50bp",
    "ctss_within_100bp",
    "ctss_within_500bp",
    "nearest_ctss_match_status",
    "nearest_same_strand_ctss_positions0",
    "nearest_same_strand_ctss_signed_distances_bp",
    "nearest_same_strand_ctss_abs_distance_bp",
    "nearest_ctss_tied_position_count",
    "tss_inside_same_strand_doga_promoter",
    "doga_promoter_tss1bp_overlap_count",
    "nearest_promoter_match_status",
    "nearest_same_strand_doga_promoter_distance_bp",
    "nearest_same_strand_doga_promoter_intervals_json",
    "nearest_promoter_tied_interval_count",
    "doga_promoter_within_100bp",
    "doga_promoter_within_500bp",
    "doga_promoter_within_1000bp",
    "doga_promoter_neighborhood_supported",
    "accessible_context_supported",
    "promoter_like_context_supported",
    "active_regulatory_context_supported",
    "epicdog_super_enhancer_context",
    "has_tss_position_support",
    "has_cage_proximity_support_500bp",
    "has_promoter_context_support",
    "support_sources",
    "support_level",
    "support_interpretation",
]
QC_FIELDS = ["metric", "value", "note"]
MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]

SUPPORT_FLAGS = [
    "ctss_exact_supported", "ctss_within_50bp", "ctss_within_100bp", "ctss_within_500bp",
    "tss_inside_same_strand_doga_promoter", "doga_promoter_within_100bp",
    "doga_promoter_within_500bp", "doga_promoter_within_1000bp",
    "doga_promoter_neighborhood_supported", "accessible_context_supported",
    "promoter_like_context_supported", "active_regulatory_context_supported",
    "epicdog_super_enhancer_context", "has_tss_position_support", "has_promoter_context_support",
    "has_cage_proximity_support_500bp",
]


@dataclass
class TranscriptModel:
    chroms: Set[str] = field(default_factory=set)
    strands: Set[str] = field(default_factory=set)
    min_start: int | None = None
    max_end: int | None = None
    exon_count: int = 0


@dataclass
class IntervalTrack:
    label: str
    starts: Dict[Tuple[str, str], List[int]]
    ends: Dict[Tuple[str, str], List[int]]
    input_rows: int
    kept_rows: int
    skipped_rows: int
    positions: Dict[Tuple[str, str], List[int]] = field(default_factory=dict)
    named_intervals: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = field(default_factory=dict)
    site_counts: Dict[Tuple[str, str], Counter] = field(default_factory=dict)


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s ] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        force=True,
    )


def fail(message: str) -> None:
    raise RuntimeError(message)


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip().strip("\r")
    return text if text else "NA"


def is_missing(value: object) -> bool:
    return clean(value) in {"NA", ".", "None", "none", "nan", "NaN"}


def to_int(value: object, label: str) -> int:
    try:
        return int(clean(value))
    except ValueError as exc:
        raise RuntimeError(f"{label} is not an integer: {value!r}") from exc


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"{label} does not exist or is empty: {path}")


def require_fields(path: Path, fieldnames: Sequence[str] | None, required: Sequence[str]) -> None:
    if fieldnames is None:
        fail(f"{path} missing TSV header")
    missing = [field for field in required if field not in fieldnames]
    if missing:
        fail(f"{path} missing required fields: {','.join(missing)}")


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    chrom_sizes: Dict[str, int] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                fail(f"{path}:{line_number} .fai must contain at least two columns")
            chrom_sizes[fields[0]] = to_int(fields[1], f"{path}:{line_number} chrom_size")
    if not chrom_sizes:
        fail(f"{path} could not read chromosome size")
    return chrom_sizes


def read_retained_ids(path: Path) -> Set[str]:
    retained: Set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, RETAINED_REQUIRED_FIELDS)
        for line_number, row in enumerate(reader, start=2):
            transcript_id = clean(row.get("query_transcript_id"))
            if is_missing(transcript_id):
                fail(f"{path}:{line_number} missing query_transcript_id")
            if clean(row.get("preASE_retained_for_downstream")) != "1":
                fail(f"{path}:{line_number} retained table contains non-retained transcript: {transcript_id}")
            if transcript_id in retained:
                fail(f"{path}:{line_number} query_transcript_id duplicate: {transcript_id}")
            retained.add(transcript_id)
    if not retained:
        fail(f"{path} contains no retained transcripts")
    return retained


def read_summary(path: Path, retained_ids: Set[str]) -> Dict[str, Dict[str, str]]:
    summary: Dict[str, Dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, SUMMARY_REQUIRED_FIELDS)
        for line_number, row in enumerate(reader, start=2):
            transcript_id = clean(row.get("query_transcript_id"))
            if transcript_id not in retained_ids:
                continue
            if transcript_id in summary:
                fail(f"{path}:{line_number} query_transcript_id duplicate: {transcript_id}")
            if clean(row.get("preASE_retained_for_downstream")) != "1":
                fail(f"{path}:{line_number} Summary contains a non-retained transcript: {transcript_id}")
            summary[transcript_id] = {field: clean(row.get(field)) for field in SUMMARY_REQUIRED_FIELDS}
    missing = retained_ids.difference(summary)
    if missing:
        fail(f"{path} is missing summary records for {len(missing)} retained transcripts, for example {sorted(missing)[0]}")
    return summary


def read_models(path: Path, retained_ids: Set[str]) -> Dict[str, TranscriptModel]:
    models: Dict[str, TranscriptModel] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                fail(f"{path}:{line_number} transcript_exon.bed less than 8 column")
            chrom, start_text, end_text, transcript_id, _gene_id, strand = fields[:6]
            if transcript_id not in retained_ids:
                continue
            start = to_int(start_text, f"{path}:{line_number} exon start")
            end = to_int(end_text, f"{path}:{line_number} exon end")
            if start < 0 or end <= start:
                fail(f"{path}:{line_number} exon coordinates invalid: {start}-{end}")
            model = models.setdefault(transcript_id, TranscriptModel())
            model.chroms.add(chrom)
            model.strands.add(strand)
            model.min_start = start if model.min_start is None else min(model.min_start, start)
            model.max_end = end if model.max_end is None else max(model.max_end, end)
            model.exon_count += 1
    return models


def read_track(path: Path, label: str, stranded: bool, positional: str = "") -> IntervalTrack:
    intervals: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    named_intervals = defaultdict(set)
    input_rows = kept_rows = skipped_rows = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            input_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                fail(f"{path}:{line_number} BED less than 3 column")
            chrom = fields[0]
            start = to_int(fields[1], f"{path}:{line_number} start")
            end = to_int(fields[2], f"{path}:{line_number} end")
            if start < 0 or end <= start:
                fail(f"{path}:{line_number} BED coordinates invalid: {start}-{end}")
            if chrom not in PRIMARY_CHROMS:
                skipped_rows += 1
                continue
            strand = "*"
            if stranded:
                if len(fields) < 6 or fields[5] not in {"+", "-"}:
                    skipped_rows += 1
                    continue
                strand = fields[5]
            intervals[(chrom, strand)].append((start, end))
            if positional == "promoter":
                named_intervals[(chrom, strand)].add((start, end, fields[3] if len(fields) > 3 else "."))
            kept_rows += 1

    starts: Dict[Tuple[str, str], List[int]] = {}
    ends: Dict[Tuple[str, str], List[int]] = {}
    for key, values in intervals.items():


        values.sort(key=lambda item: (item[0], item[1]))
        key_starts = [start for start, _end in values]
        key_ends = sorted(end for _start, end in values)
        starts[key] = key_starts
        ends[key] = key_ends
    if not kept_rows:
        fail(f"{path} no retain primary chromosome BED record")
    return IntervalTrack(label, starts, ends, input_rows, kept_rows, skipped_rows,
                         named_intervals={key: sorted(values) for key, values in named_intervals.items()})


def read_ctss(path: Path, label: str, chrom_sizes: Mapping[str, int]) -> IntervalTrack:

    counts = defaultdict(Counter)
    input_rows = kept_rows = skipped_rows = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            input_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] not in {"+", "-"}:
                fail(f"{path}:{line_number} CTSS must contain six columns and a valid +/- strand")
            chrom, strand = fields[0], fields[5]
            start = to_int(fields[1], f"{path}:{line_number} start")
            end = to_int(fields[2], f"{path}:{line_number} end")
            if start < 0 or end - start != 1:
                fail(f"{path}:{line_number} CTSS must be a 1 bp BED interval with a nonnegative start")
            if chrom not in PRIMARY_CHROMS:
                skipped_rows += 1
            else:
                if chrom not in chrom_sizes or end > chrom_sizes[chrom]:
                    fail(f"{path}:{line_number} CTSS exceeds CanFam3.1 chromosome span")
                counts[(chrom, strand)][start] += 1
                kept_rows += 1
            if input_rows % 5000000 == 0:
                logging.info("%s read %d CTSS records", label, input_rows)
    if not kept_rows:
        fail(f"{path} no retain primary chromosome CTSS")
    return IntervalTrack(label, {}, {}, input_rows, kept_rows, skipped_rows,
                         positions={key: sorted(values) for key, values in counts.items()},
                         site_counts=dict(counts))


def nearest_ctss(track: IntervalTrack, chrom: str, strand: str, tss0: int) -> Dict[str, object]:

    positions = track.positions.get((chrom, strand), [])
    if not positions:
        return {**{key: "NA" for key in (
                    "nearest_same_strand_ctss_positions0",
                    "nearest_same_strand_ctss_signed_distances_bp",
                    "nearest_same_strand_ctss_abs_distance_bp",
                    "nearest_ctss_tied_position_count")},
                "nearest_ctss_match_status": "no_same_strand_reference",
                **{flag: "no" for flag in ("ctss_exact_supported", "ctss_within_50bp", "ctss_within_100bp", "ctss_within_500bp")}}
    i = bisect_left(positions, tss0)
    candidates = positions[max(0, i-1):min(len(positions), i+1)]
    distance = min(abs(p-tss0) for p in candidates)
    tied = [p for p in candidates if abs(p-tss0) == distance]
    direction = 1 if strand == "+" else -1
    return {
        "nearest_ctss_match_status": "found",
        "nearest_same_strand_ctss_positions0": ";".join(map(str, tied)),
        "nearest_same_strand_ctss_signed_distances_bp": ";".join(str(direction*(p-tss0)) for p in tied),
        "nearest_same_strand_ctss_abs_distance_bp": distance,
        "nearest_ctss_tied_position_count": len(tied),
        "ctss_exact_supported": yes_no(distance == 0),
        **{f"ctss_within_{cutoff}bp": yes_no(distance <= cutoff) for cutoff in (50, 100, 500)},
    }


def nearest_promoter(track: IntervalTrack, chrom: str, strand: str, tss0: int) -> Dict[str, object]:

    intervals = track.named_intervals.get((chrom, strand), [])
    if not intervals:
        return {"nearest_promoter_match_status": "no_same_strand_reference",
                **{f"doga_promoter_within_{cutoff}bp": "no" for cutoff in (100, 500, 1000)}}
    distances = [max(start-tss0, tss0-(end-1), 0) for start, end, _ in intervals]
    distance = min(distances)
    tied = [record for record, d in zip(intervals, distances) if d == distance]
    return {
        "nearest_promoter_match_status": "found",
        "nearest_same_strand_doga_promoter_distance_bp": distance,

        "nearest_same_strand_doga_promoter_intervals_json": json.dumps(
            [{"chrom": chrom, "start0": s, "end0": e, "strand": strand, "source_name": name}
             for s, e, name in tied], ensure_ascii=False, separators=(",", ":")),
        "nearest_promoter_tied_interval_count": len(tied),
        **{f"doga_promoter_within_{cutoff}bp": yes_no(distance <= cutoff) for cutoff in (100, 500, 1000)},
    }


def count_overlaps(track: IntervalTrack, chrom: str, strand: str, start: int, end: int) -> int:
    key = (chrom, strand if strand in {"+", "-"} else "*")
    if track.site_counts:
        positions = track.positions.get(key, [])
        counts = track.site_counts.get(key, {})
        return sum(counts[p] for p in positions[bisect_left(positions, start):bisect_left(positions, end)])
    starts = track.starts.get(key, [])
    ends = track.ends.get(key, [])
    if not starts:
        return 0

    return max(0, bisect_left(starts, end) - bisect_right(ends, start))


def empty_support(reason: str, coordinate_group: str = "not_assessed") -> Dict[str, object]:
    return {
        "tss_inferred": "no",
        "reference_assessable": "no",
        "not_assessed_reason": reason,
        "tss0": "NA",
        "tss_1bp_start0": "NA",
        "tss_1bp_end0": "NA",
        "tss_pm500_start0": "NA",
        "tss_pm500_end0": "NA",
        "promoter_1kb_100bp_start0": "NA",
        "promoter_1kb_100bp_end0": "NA",
        "coordinate_overlap_group": coordinate_group,
        "doga_ctss_tss1bp_overlap_count": "NA",
        "doga_promoter_tss_pm500_overlap_count": "NA",
        "doga_promoter_promoter_overlap_count": "NA",
        "barkbase_atac_promoter_overlap_count": "NA",
        "epicdog_promoter_like_promoter_overlap_count": "NA",
        "epicdog_active_regulatory_promoter_overlap_count": "NA",
        "epicdog_super_enhancer_promoter_overlap_count": "NA",
        **{flag: "NA" for flag in SUPPORT_FLAGS},
        "support_sources": "NA",
        "support_level": "not_assessed",
        "support_interpretation": reason,
    }


def derive_support(
    transcript_id: str,
    summary: Mapping[str, str],
    model: TranscriptModel | None,
    chrom_sizes: Mapping[str, int],
    tracks: Mapping[str, IntervalTrack],
) -> Dict[str, object]:
    base: Dict[str, object] = {
        "query_transcript_id": transcript_id,
        "query_gene_id": summary["query_gene_id"],
        "class_code": summary["class_code"],
        "novelty_category": summary["novelty_category"],
        "FEELnc_biotype": summary["FEELnc_biotype"],
        "FEELnc_codpot_category": summary["FEELnc_codpot_category"],
        "FEELnc_lncRNA_class": summary["FEELnc_lncRNA_class"],
        "candidate_gene_strand": summary["candidate_gene_strand"],
        "gene_locus_filter_reason": summary["gene_locus_filter_reason"],
        "preASE_retained_for_downstream": summary["preASE_retained_for_downstream"],
        "transcript_chr": summary["chr"],
        "transcript_strand": summary["transcript_strand"],
        "transcript_exon_count": summary["transcript_exon_count"],
        "transcript_span_start0": summary["transcript_span_start0"],
        "transcript_span_end0": summary["transcript_span_end0"],
        "transcript_span_length": summary["transcript_span_length"],
        "transcript_exonic_length": summary["transcript_exonic_length"],
    }
    if model is None:
        base.update(empty_support("Retained transcript missing from transcript_exon_bed."))
        return base
    if len(model.chroms) != 1 or len(model.strands) != 1:
        base.update(empty_support("Exons of the same transcript span multiple chromosomes or strands."))
        return base
    chrom = next(iter(model.chroms))
    strand = next(iter(model.strands))
    base.update(
        {
            "transcript_chr": chrom,
            "transcript_strand": strand,
            "transcript_exon_count": model.exon_count,
            "transcript_span_start0": model.min_start,
            "transcript_span_end0": model.max_end,
            "transcript_span_length": model.max_end - model.min_start if model.min_start is not None and model.max_end is not None else "NA",
        }
    )
    if strand not in {"+", "-"} or model.min_start is None or model.max_end is None:
        base.update(empty_support("Transcript strand is not +/- or the exon model is empty."))
        return base
    tss0 = model.min_start if strand == "+" else model.max_end - 1
    base.update({"tss_inferred": "yes", "tss0": tss0,
                 "tss_1bp_start0": tss0, "tss_1bp_end0": tss0+1})
    if chrom not in chrom_sizes or chrom not in PRIMARY_CHROMS:
        reason = "Reference length unavailable for this chromosome." if chrom not in chrom_sizes else "Non-primary chromosome; overlap with canine TSS/promoter tracks was not assessed."
        group = "missing_chromosome_size" if chrom not in chrom_sizes else "alt_or_unplaced_contig"
        base.update(empty_support(reason, group))
        base.update({"tss_inferred": "yes", "tss0": tss0,
                     "tss_1bp_start0": tss0, "tss_1bp_end0": tss0+1})
        return base

    chrom_size = chrom_sizes[chrom]
    if not 0 <= tss0 < chrom_size:
        fail(f"{transcript_id} infer TSS exceeds chromosome span: {chrom}:{tss0}")
    promoter_start, promoter_end = (
        (max(0, tss0 - 1000), min(chrom_size, tss0 + 101))
        if strand == "+"
        else (max(0, tss0 - 100), min(chrom_size, tss0 + 1001))
    )
    tss_1bp = (tss0, tss0 + 1)
    tss_pm500 = (max(0, tss0 - 500), min(chrom_size, tss0 + 501))
    counts = {
        "doga_ctss_tss1bp_overlap_count": count_overlaps(tracks["doga_ctss"], chrom, strand, *tss_1bp),
        "doga_promoter_tss1bp_overlap_count": count_overlaps(tracks["doga_promoter"], chrom, strand, *tss_1bp),
        "doga_promoter_tss_pm500_overlap_count": count_overlaps(tracks["doga_promoter"], chrom, strand, *tss_pm500),
        "doga_promoter_promoter_overlap_count": count_overlaps(tracks["doga_promoter"], chrom, strand, promoter_start, promoter_end),
        "barkbase_atac_promoter_overlap_count": count_overlaps(tracks["barkbase_atac"], chrom, "*", promoter_start, promoter_end),
        "epicdog_promoter_like_promoter_overlap_count": count_overlaps(tracks["epicdog_promoter_like"], chrom, "*", promoter_start, promoter_end),
        "epicdog_active_regulatory_promoter_overlap_count": count_overlaps(tracks["epicdog_active"], chrom, "*", promoter_start, promoter_end),
        "epicdog_super_enhancer_promoter_overlap_count": count_overlaps(tracks["epicdog_super"], chrom, "*", promoter_start, promoter_end),
    }
    doga_ctss = counts["doga_ctss_tss1bp_overlap_count"] > 0
    doga_promoter = counts["doga_promoter_tss_pm500_overlap_count"] > 0 or counts["doga_promoter_promoter_overlap_count"] > 0
    barkbase = counts["barkbase_atac_promoter_overlap_count"] > 0
    epic_promoter = counts["epicdog_promoter_like_promoter_overlap_count"] > 0
    epic_active = counts["epicdog_active_regulatory_promoter_overlap_count"] > 0
    epic_super = counts["epicdog_super_enhancer_promoter_overlap_count"] > 0
    inside_promoter = counts["doga_promoter_tss1bp_overlap_count"] > 0

    positional = doga_ctss or inside_promoter


    promoter_match = nearest_promoter(tracks["doga_promoter"], chrom, strand, tss0)
    proximity = promoter_match["doga_promoter_within_500bp"]
    chromatin = barkbase or epic_promoter
    sources = []
    if doga_ctss:
        sources.append("DoGA_CTSS_TSS_1bp")
    if doga_promoter:
        sources.append("DoGA_promoter_neighborhood")
    if inside_promoter:
        sources.append("TSS_inside_DoGA_promoter")
    if proximity == "yes":
        sources.append("DoGA_robust_promoter_within_500bp")
    if barkbase:
        sources.append("BarkBase_ATAC_promoter")
    if epic_promoter:
        sources.append("EpiC_Dog_promoter_like")
    if epic_active:
        sources.append("EpiC_Dog_active_regulatory_states_1_7")
    if epic_super:
        sources.append("EpiC_Dog_super_enhancer_context")
    if positional and chromatin:
        level = "position_and_promoter_context"
        interpretation = "Positional evidence supports the candidate start, and the promoter window overlaps ATAC or promoter-like states; this does not validate all exact boundaries."
    elif positional:
        level = "position_support_only"
        interpretation = "The start exactly overlaps a CTSS or lies within a same-strand DoGA promoter, without selected promoter-context evidence."
    elif chromatin:
        level = "promoter_context_only"
        interpretation = "The promoter window overlaps BarkBase ATAC or EpiC Dog states 1-4, without selected positional evidence."
    else:
        level = "no_position_or_promoter_context"
        interpretation = "No selected positional/promoter-context evidence; inspect nearby CTSS, promoter, active, or super-enhancer evidence separately."
    base.update(
        {
            "tss_inferred": "yes",
            "reference_assessable": "yes",
            "not_assessed_reason": "NA",
            "tss0": tss0,
            "tss_1bp_start0": tss_1bp[0],
            "tss_1bp_end0": tss_1bp[1],
            "tss_pm500_start0": tss_pm500[0],
            "tss_pm500_end0": tss_pm500[1],
            "promoter_1kb_100bp_start0": promoter_start,
            "promoter_1kb_100bp_end0": promoter_end,
            "coordinate_overlap_group": "primary_chr",
            **counts,
            **nearest_ctss(tracks["doga_ctss"], chrom, strand, tss0),
            **promoter_match,
            "tss_inside_same_strand_doga_promoter": yes_no(inside_promoter),
            "doga_promoter_neighborhood_supported": yes_no(doga_promoter),
            "accessible_context_supported": yes_no(barkbase),
            "promoter_like_context_supported": yes_no(epic_promoter),
            "active_regulatory_context_supported": yes_no(epic_active),
            "epicdog_super_enhancer_context": yes_no(epic_super),
            "has_tss_position_support": yes_no(positional),
            "has_cage_proximity_support_500bp": proximity,
            "has_promoter_context_support": yes_no(chromatin),
            "support_sources": "|".join(sources) if sources else "NA",
            "support_level": level,
            "support_interpretation": interpretation,
        }
    )
    return base


def check_primary_support(row: Mapping[str, object]) -> None:

    if row["reference_assessable"] != "yes":
        return
    transcript = row["query_transcript_id"]
    flags = ["ctss_exact_supported", "ctss_within_50bp", "ctss_within_100bp", "ctss_within_500bp"]
    for names in (flags,
                  ["tss_inside_same_strand_doga_promoter", "doga_promoter_within_100bp",
                   "doga_promoter_within_500bp", "doga_promoter_within_1000bp"]):
        values = [row[name] == "yes" for name in names]
        if any(a and not b for a, b in zip(values, values[1:])):
            fail(f"{transcript} cumulative thresholds are not monotonic: {names}")
    restricted_distance = row.get("nearest_same_strand_ctss_abs_distance_bp", "NA")
    if (row["doga_ctss_tss1bp_overlap_count"] > 0) != (restricted_distance == 0):
        fail(f"{transcript}: CTSS overlap and nearest exact-match results disagree")
    if (row["doga_promoter_tss1bp_overlap_count"] > 0) != (row.get("nearest_same_strand_doga_promoter_distance_bp") == 0):
        fail(f"{transcript}: promoter overlap and nearest inside-promoter results disagree")


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean(row.get(field)) for field in fields})


def write_bed(path: Path, rows: Iterable[Sequence[object]]) -> None:
    with path.open("w") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def manifest_row(item: str, path: Path, note: str) -> Dict[str, object]:
    return {"item": item, "path": str(path), "exists": int(path.exists()), "size_bytes": path.stat().st_size if path.exists() else "NA", "note": note}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Transcript-first dog TSS/promoter support annotation.")
    parser.add_argument("--retained-transcript-tsv", type=Path, default=DEFAULT_RETAINED_TSV)
    parser.add_argument("--transcript-summary-tsv", type=Path, default=DEFAULT_SUMMARY_TSV)
    parser.add_argument("--transcript-exon-bed", type=Path, default=DEFAULT_TRANSCRIPT_EXON_BED)
    parser.add_argument("--chrom-sizes", type=Path, default=DEFAULT_CHROM_SIZES)
    parser.add_argument("--doga-ctss-bed", type=Path, default=DEFAULT_DOGA_CTSS_BED)
    parser.add_argument("--doga-promoter-bed", type=Path, default=DEFAULT_DOGA_PROMOTER_BED)
    parser.add_argument("--barkbase-atac-bed", type=Path, default=DEFAULT_BARKBASE_ATAC_BED)
    parser.add_argument("--epicdog-promoter-like-bed", type=Path, default=DEFAULT_EPICDOG_PROMOTER_BED)
    parser.add_argument("--epicdog-active-bed", type=Path, default=DEFAULT_EPICDOG_ACTIVE_BED)
    parser.add_argument("--epicdog-super-enhancer-bed", type=Path, default=DEFAULT_EPICDOG_SUPER_BED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S"))
    return parser.parse_args()


def main() -> None:
    setup_logging()
    args = parse_args()
    if not args.outdir.is_absolute():
        fail("--outdir must be an absolute path")
    if len(args.run_id) != 15 or args.run_id[8] != "_" or not (args.run_id[:8]+args.run_id[9:]).isdigit():
        fail("--run-id must use YYYYMMDD_HHMMSS format")
    args.outdir = args.outdir / args.run_id
    start_time = datetime.now().strftime("%F %T")
    output_paths = {
        "support": args.outdir / "all_candidate_transcript_TSS_promoter_support.tsv",
        "tss_1bp": args.outdir / "all_candidate_TSS_1bp.bed",
        "tss_pm500": args.outdir / "all_candidate_TSS_pm500bp.bed",
        "promoter": args.outdir / "all_candidate_promoter_1kb_100bp.bed",
        "qc": args.outdir / "TSS_promoter_support_qc_summary.tsv",
        "manifest": args.outdir / "TSS_promoter_support_manifest.tsv",
        "counts": args.outdir / "TSS_support_counts.tsv",
        "policy": args.outdir / "TSS_support_definitions.txt",
        "site_qc": args.outdir / "CTSS_reference_site_qc.tsv",
    }
    for path, label in [
        (args.retained_transcript_tsv, "Step 2.1 retained transcript table."),
        (args.transcript_summary_tsv, "Step 2.1 transcript summary"),
        (args.transcript_exon_bed, "Step 1.6 transcript exon BED"),
        (args.chrom_sizes, "CanFam3.1 .fai"),
        (args.doga_ctss_bed, "DoGA CTSS"),
        (args.doga_promoter_bed, "DoGA promoter"),
        (args.barkbase_atac_bed, "BarkBase ATAC"),
        (args.epicdog_promoter_like_bed, "EpiC Dog promoter-like"),
        (args.epicdog_active_bed, "EpiC Dog active H3K4me3/H3K27ac"),
        (args.epicdog_super_enhancer_bed, "EpiC Dog super-enhancer"),
    ]:
        require_file(path, label)
    if args.outdir.exists():
        fail(f"run directory already exists, refuse to overwrite: {args.outdir}")

    logging.info("python=%s version=%s concurrent_tasks=1 threads_per_task=1", sys.executable, sys.version.split()[0])
    logging.info("output=%s; current CTSS default reference by CNS enrichment promoter region filter, is not brain sample exclusive CTSS", args.outdir)
    logging.info("all_cage_sensitivity_analysis=not_performed; full analysis was not performed for CAGE reference sensitivity analysis")
    logging.info("[stage 1/5] read and validate Step 2.1 retained transcript input")
    retained_ids = read_retained_ids(args.retained_transcript_tsv)
    summary = read_summary(args.transcript_summary_tsv, retained_ids)
    logging.info("retained transcripts=%d", len(retained_ids))

    logging.info("[stage 2/5] read exon model and CanFam3.1 chromosome sizes")
    models = read_models(args.transcript_exon_bed, retained_ids)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)

    logging.info("[stage 3/5] read dog TSS/promoter with chromatin track")
    tracks = {
        "doga_ctss": read_ctss(args.doga_ctss_bed, "DoGA_CTSS", chrom_sizes),
        "doga_promoter": read_track(args.doga_promoter_bed, "DoGA_promoter", stranded=True, positional="promoter"),
        "barkbase_atac": read_track(args.barkbase_atac_bed, "BarkBase_ATAC", stranded=False),
        "epicdog_promoter_like": read_track(args.epicdog_promoter_like_bed, "EpiC_Dog_promoter_like", stranded=False),
        "epicdog_active": read_track(args.epicdog_active_bed, "EpiC_Dog_active_regulatory_states_1_7", stranded=False),
        "epicdog_super": read_track(args.epicdog_super_enhancer_bed, "EpiC_Dog_super_enhancer", stranded=False),
    }

    logging.info("[stage 4/5] infer transcript TSS: assess exact-position evidence, same-strand robust promoter proximity <=500 bp, and contextual evidence; by default, CTSS is a region-restricted subset and promoters comprise the full robust set")
    support_rows: List[Dict[str, object]] = []
    for index, transcript_id in enumerate(sorted(retained_ids), start=1):
        support_rows.append(derive_support(transcript_id, summary[transcript_id], models.get(transcript_id), chrom_sizes, tracks))
        if index % 1000 == 0:
            logging.info("completed %d/%d transcripts", index, len(retained_ids))
    for row in support_rows:
        check_primary_support(row)

    tss_1bp_rows = []
    tss_pm500_rows = []
    promoter_rows = []
    for row in support_rows:
        if row["tss_inferred"] != "yes":
            continue
        common = [row["transcript_chr"], row["query_transcript_id"], 0, row["transcript_strand"]]
        tss_1bp_rows.append([common[0], row["tss_1bp_start0"], row["tss_1bp_end0"], common[1], common[2], common[3]])
        if row["reference_assessable"] != "yes":
            continue
        tss_pm500_rows.append([common[0], row["tss_pm500_start0"], row["tss_pm500_end0"], common[1], common[2], common[3]])
        promoter_rows.append([common[0], row["promoter_1kb_100bp_start0"], row["promoter_1kb_100bp_end0"], common[1], common[2], common[3]])

    logging.info("[stage 5/5] write transcript support table, BED, QC with manifest")

    args.outdir.mkdir(parents=True, exist_ok=False)
    write_tsv(output_paths["support"], OUTPUT_FIELDS, support_rows)
    track = tracks["doga_ctss"]
    unique = sum(map(len, track.positions.values()))
    site_rows = [{"reference": "region_restricted_CAGE", "input_record_count": track.input_rows,
                  "kept_primary_record_count": track.kept_rows,
                  "nonprimary_record_count": track.skipped_rows,
                  "unique_chr_pos_strand_site_count": unique,
                  "duplicate_record_count": track.kept_rows - unique}]
    write_tsv(output_paths["site_qc"], ["reference", "input_record_count", "kept_primary_record_count",
              "nonprimary_record_count", "unique_chr_pos_strand_site_count", "duplicate_record_count"], site_rows)
    write_bed(output_paths["tss_1bp"], tss_1bp_rows)
    write_bed(output_paths["tss_pm500"], tss_pm500_rows)
    write_bed(output_paths["promoter"], promoter_rows)
    level_counts = Counter(clean(row["support_level"]) for row in support_rows)
    qc_rows: List[Dict[str, object]] = [
        {"metric": "retained_transcript_count", "value": len(retained_ids), "note": "Number of Step 2.1 retained query_transcript_id values."},
        {"metric": "summary_transcript_count", "value": len(summary), "note": "Number of retained transcripts in the Step 2.1 transcript summary."},
        {"metric": "exon_model_count", "value": len(models), "note": "Number of models in transcript_exon.bed matching retained transcripts."},
        {"metric": "tss_inferred_count", "value": sum(row["tss_inferred"] == "yes" for row in support_rows), "note": "TSS inferred from transcript models, including unevaluated non-primary chromosomes."},
        {"metric": "reference_assessable_count", "value": sum(row["reference_assessable"] == "yes" for row in support_rows), "note": "Primary denominator for all support proportions."},
        {"metric": "not_assessed_count", "value": sum(row["reference_assessable"] != "yes" for row in support_rows), "note": "Unevaluated records are not counted as unsupported."},
        {"metric": "output_transcript_rows", "value": len(support_rows), "note": "One row per retained query_transcript_id."},
        {"metric": "all_candidate_TSS_1bp_bed_rows", "value": len(tss_1bp_rows), "note": "Number of TSS 1bp BED rows."},
        {"metric": "all_candidate_promoter_bed_rows", "value": len(promoter_rows), "note": "Number of promoter BED rows."},
    ]
    for level, count in sorted(level_counts.items()):
        qc_rows.append({"metric": f"support_level_{level}_count", "value": count, "note": "transcript-level TSS/promoter support"})
    for track in tracks.values():
        qc_rows.extend(
            [
                {"metric": f"{track.label}_input_rows", "value": track.input_rows, "note": "Number of BED records scanned."},
                {"metric": f"{track.label}_kept_primary_rows", "value": track.kept_rows, "note": "Number of primary-chromosome records retained."},
                {"metric": f"{track.label}_skipped_rows", "value": track.skipped_rows,
                 "note": "Number of non-primary-chromosome records; invalid CTSS strands raise an error and are not counted as skipped."
                 if track is tracks["doga_ctss"] else "Number of records with non-primary chromosomes or invalid strands."},
            ]
        )

    for reason, n in Counter(row["not_assessed_reason"] for row in support_rows if row["reference_assessable"] != "yes").items():
        qc_rows.append({"metric": "not_assessed_reason", "value": n, "note": reason})
    write_tsv(output_paths["qc"], QC_FIELDS, qc_rows)
    assessed = [row for row in support_rows if row["reference_assessable"] == "yes"]
    summary_rows = []
    for flag in SUPPORT_FLAGS:
        n = sum(row[flag] == "yes" for row in assessed)
        summary_rows.append({"metric": flag, "supported_transcripts": n,
                             "assessable_transcripts": len(assessed), "total_transcripts": len(support_rows),
                             "not_assessed": len(support_rows)-len(assessed),
                             "percent": 100*n/len(assessed) if assessed else "NA"})
    write_tsv(output_paths["counts"], ["metric", "supported_transcripts", "assessable_transcripts",
                                      "total_transcripts", "not_assessed", "percent"], summary_rows)
    with output_paths["policy"].open("x", encoding="utf-8") as handle:
        handle.write("""schema_version=6
unit=retained transcript; candidate selection remains Step 2.1
reference_assembly=CanFam3.1; reference chromosomes=chr1-38,chrX
tss=minimum exon start on +; maximum exon end minus 1 on - (BED 0-based half-open)
tss_inferred is distinct from reference_assessable; unassessed evidence=NA
CTSS nearest search=whole chromosome, same strand, distinct positions; not restricted to 500 bp
CTSS signed distances=(reference-candidate) on +, (candidate-reference) on -; positive=downstream
CTSS tied positions and signed distances are parallel semicolon-separated lists
CTSS thresholds=0,50,100,500 bp (cumulative, not mutually exclusive)
CTSS overlap_count counts input records, including repeated sample positions; not independent sites or read depth
Promoter distance=min distance to a covered genomic base; inside=0; BED end is excluded
Promoter nearest intervals=JSON list of objects with chrom,start0,end0,strand,source_name; all ties retained
Read the result with a TSV/CSV parser first (standard CSV quoting), then JSON-decode the intervals field; source names may contain semicolons
Source promoter names may contain original CanFam4 coordinates; use output interval coordinates for CanFam3.1
Promoter near thresholds=100,500,1000 bp, cumulative
Promoter neighborhood=TSS +/-500 OR strand-oriented upstream1000/downstream100 window
The two inclusive-base windows span 1001 and 1101 bases before chromosome-boundary clipping
ATAC and EpiC overlap=unstranded, at least 1 overlapping base in the promoter window
has_tss_position_support=exact CTSS OR candidate TSS inside same-strand DoGA promoter
has_cage_proximity_support_500bp=doga_promoter_within_500bp; same chromosome and strand; minimum distance to a covered promoter base <=500 bp, including peak overlap
recommended_neighborhood_metric=has_cage_proximity_support_500bp; CAGE-derived transcription-initiation proximity support, not precise TSS validation
sensitivity_comparison=doga_promoter_within_100bp; retain strict has_tss_position_support separately
Proximity uses the robust promoter input only, not the region-restricted CTSS input; unassessed=NA
Strict position and proximity flags are overlapping evidence sets, not mutually exclusive groups; do not sum their counts
has_promoter_context_support=ATAC OR EpiC states1-4
support_level and support_interpretation summarize strict has_tss_position_support and has_promoter_context_support only; the 500bp proximity flag is reported separately
support_sources includes DoGA_robust_promoter_within_500bp when the proximity flag is yes; a source label does not imply strict position support
METHOD REFERENCES:
Hardwick et al. 2019, Frontiers in Genetics, doi:10.3389/fgene.2019.00309; novel-transcript TSS within 500 bp of a same-strand FANTOM CAGE robust peak; adapted here to dog DoGA robust peaks, with explicit nearest-covered-base distance
Hu et al. 2022, Nature Communications, doi:10.1038/s41467-022-34568-z; CAGE/chromatin evidence within 500 nt around TSS; contextual support for a neighborhood criterion, not a claim of exact pipeline replication
DoGA CAGE resource: Takan et al., The DoGA Consortium Atlas of Canine Enhancers and Promoters Across Tissues and Development, doi:10.64898/2026.03.24.713649 (preprint); local CAGE metadata match resource scale, not a file-version identity check
DoGA 2024 doi:10.1038/s41467-024-52798-1 primarily describes STRT2; do not label the local CAGE reference as that STRT2 dataset
active_regulatory_context=EpiC states1-7 (includes enhancer-like states5-7); not raw H3K4me3/H3K27ac peaks
DEFAULT REFERENCE PROVENANCE (custom paths require independent provenance review):
PRIMARY CTSS REFERENCE: region-restricted DoGA CAGE CTSS
CAGE-derived; all 114 local CTSS files initially pooled; unstranded genomic overlap with 13 selected CNS/brain/pituitary-related enriched promoter groups
Not selected by CTSS library tissue metadata; then lifted CanFam4 -> CanFam3.1; candidate matching requires same strand
ctss_reference_scope=region_restricted_DoGA_CAGE
all_cage_sensitivity_analysis=not_performed
promoter_reference_scope=all_DoGA_robust_promoters_successfully_lifted
default_reference_combination=region_restricted_CTSS_plus_all_robust_promoters
This run expands the promoter reference only; it is not an all-CTSS analysis.
The effect of promoter-region reference restriction on CTSS support rates was not quantified.
unique_chr_pos_strand_site_count covers kept primary records only; duplicate_record_count=kept_primary_record_count-unique_chr_pos_strand_site_count
Record count != independent CTSS sites != read depth; no new score/read-count cutoff
Primary CTSS/promoter cumulative-threshold and overlap/nearest consistency checks are retained as requested; failures occur before any success marker
DoGA CTSS and DoGA promoter evidence are related CAGE-derived annotations and should not be described as statistically independent validation sources.
Missing CTSS support means no support in the selected reference, not absence of transcription initiation.
DoGA promoters=all robust promoters from _Promoter_robust.bed, successfully lifted to CanFam3.1
DoGA CTSS=all-sample CTSS selected by unstranded genomic overlap with the original 13 brain/CNS/pituitary enriched promoter groups before liftover; not reselected using the all-promoter input
Thus brain_CNS CTSS is not a sample-metadata-selected brain CTSS set; CTSS and promoter evidence are not independent
DoGA CTSS and promoter coordinates were lifted from CanFam4 to CanFam3.1; only mapped records are available
BarkBase ATAC=Adult5 occipital cortex and Adult4/5/6 pituitary peak union; any selected peak, no replicate-count threshold
EpiC brain/CNS=CL cerebellum and CR cerebrum; promoter-like=states1-4; source coordinates CanFam3.1
No new score/read-count threshold is applied here; no claim of independent donor validation
Missing overlap means no selected-reference support, not absence of transcription
concurrent_tasks=1; threads_per_task=1
Execution completion is not a post-run quality audit
""")
    end_time = datetime.now().strftime("%F %T")
    manifest_rows = [
        manifest_row("retained_transcript_tsv_input", args.retained_transcript_tsv, "Step 2.1 retained transcript candidates"),
        manifest_row("transcript_summary_tsv_input", args.transcript_summary_tsv, "Step 2.1 one-row-per-retained-transcript summary"),
        manifest_row("transcript_exon_bed_input", args.transcript_exon_bed, "five-dogs Step 1.6 transcript exon BED"),
        manifest_row("chrom_sizes_input", args.chrom_sizes, "CanFam3.1 FASTA index"),
        manifest_row("doga_ctss_bed_input", args.doga_ctss_bed, "Default reference: CTSS from all samples filtered by CNS-enriched promoter regions and then lifted over; not a brain-sample-only set."),
        manifest_row("doga_promoter_bed_input", args.doga_promoter_bed, "Default reference: all DoGA robust promoters successfully converted to CanFam3.1, not only a brain/CNS-enriched subset; verify provenance separately for custom paths."),
        manifest_row("barkbase_atac_bed_input", args.barkbase_atac_bed, "BarkBase selected merged ATAC regions"),
        manifest_row("epicdog_promoter_like_bed_input", args.epicdog_promoter_like_bed, "EpiC Dog brain/CNS promoter-like states"),
        manifest_row("epicdog_active_bed_input", args.epicdog_active_bed, "EpiC Dog states 1-7 active regulatory context, including enhancer-like states."),
        manifest_row("epicdog_super_enhancer_bed_input", args.epicdog_super_enhancer_bed, "EpiC Dog super-enhancer context"),
        manifest_row("transcript_support_output", output_paths["support"], "Step 2.2 transcript-level TSS/promoter support"),
        manifest_row("tss_1bp_bed_output", output_paths["tss_1bp"], "candidate transcript TSS 1bp BED"),
        manifest_row("tss_pm500_bed_output", output_paths["tss_pm500"], "candidate transcript TSS plus/minus 500bp BED"),
        manifest_row("promoter_bed_output", output_paths["promoter"], "candidate transcript strand-aware promoter BED"),
        manifest_row("qc_output", output_paths["qc"], "Step 2.2 QC summary"),
        manifest_row("support_counts_output", output_paths["counts"], "Numerator, denominator, unevaluated count, and percentage for each metric."),
        manifest_row("definitions_output", output_paths["policy"], "Updated metric definitions, reference limitations, and distance rules."),
        {"item": "run_time", "path": f"{start_time} -> {end_time}", "exists": "NA", "size_bytes": "NA", "note": "script runtime"},
    ]
    manifest_rows.append(manifest_row("ctss_site_qc_output", output_paths["site_qc"], "Restricted CTSS record count and unique-site count."))
    write_tsv(output_paths["manifest"], MANIFEST_FIELDS, manifest_rows)
    logging.info("compute completed: %d transcript rows -> %s", len(support_rows), output_paths["support"])


if __name__ == "__main__":
    main()
