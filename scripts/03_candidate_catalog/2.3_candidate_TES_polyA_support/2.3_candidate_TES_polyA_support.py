#!/usr/bin/python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import os

for thread_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[thread_key] = "1"

import argparse
import csv
import logging
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, List, Mapping, Sequence, Set, Tuple


PROJECT_ROOT = Path("/path/to/project/analysis_20260709/5-dogs/9.new-gene")
STEP21_ROOT = PROJECT_ROOT / "2.1_transcript_first_preASE_candidate_assessment"
STEP16_EXON_BED = Path("/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/transcript_exon.bed")
REFERENCE_ROOT = Path("/path/to/project/legacy")
GENOME_FASTA = REFERENCE_ROOT / "reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
CHROM_SIZES = Path(str(GENOME_FASTA) + ".fai")
TES_ROOT = REFERENCE_ROOT / "five-dogs/9.new-gene/reference/TES"
ANIMAL_ROOT = TES_ROOT / "1.Animal-APAdb-dog-PAS/results-v2"
GSE_ROOT = TES_ROOT / "2.GEO-GSE30198-PolyA-Seq/results-v1"

DEFAULT_RETAINED_TSV = STEP21_ROOT / "2.1_retained_transcript_candidates.preASE.tsv"
DEFAULT_SUMMARY_TSV = STEP21_ROOT / "2.1_candidate_transcript_summary.preASE.tsv"
DEFAULT_OUTDIR = PROJECT_ROOT / "2.3_candidate_TES_polyA_support"
PRIMARY_CHROMS = {*(f"chr{index}" for index in range(1, 39)), "chrX"}
CANONICAL_PAS_MOTIFS = ("AATAAA", "ATTAAA")
VARIANT_PAS_MOTIFS = ("TATAAA", "AGTAAA", "AAGAAA", "AATATA", "AATACA", "CATAAA", "GATAAA")
ALL_PAS_MOTIFS = CANONICAL_PAS_MOTIFS + VARIANT_PAS_MOTIFS


EXTENDED_PAS_MOTIFS = ALL_PAS_MOTIFS + ("AATGAA", "AATAGA")
SPLICED_FIELDS = (
    "spliced_pas_scan_available", "spliced_pas_scan_failure_reason",
    "spliced_terminal_50nt_sequence", "spliced_terminal_50nt_length",
    "spliced_pas_selected9_found", "spliced_pas_extended11_found",
    "spliced_pas_canonical_found", "spliced_pas_variant_found",
    "spliced_pas_hit_list", "spliced_pas_best_motif",
    "spliced_pas_best_distance_to_TES_nt", "downstream_20nt_sequence",
    "downstream_20nt_A_percent", "downstream_20nt_assessable",
)
RADIUS100_FLAGS = tuple(f"{prefix}_within_100bp" for prefix in
    ("polyaseq", "animal_apadb", "external", "brain_CNS", "brain_or_pituitary", "nonbrain"))

RETAINED_REQUIRED_FIELDS = ("query_transcript_id", "preASE_retained_for_downstream")
SUMMARY_REQUIRED_FIELDS = (
    "query_transcript_id", "query_gene_id", "chr", "transcript_strand",
    "transcript_span_start0", "transcript_span_end0", "transcript_span_length",
    "transcript_exonic_length", "transcript_exon_count", "class_code", "novelty_category",
    "FEELnc_biotype", "FEELnc_codpot_category", "FEELnc_lncRNA_class",
    "candidate_gene_strand", "gene_locus_filter_reason", "preASE_retained_for_downstream",
)
OUTPUT_FIELDS = (
    "query_transcript_id", "query_gene_id", "class_code", "novelty_category", "FEELnc_biotype",
    "FEELnc_codpot_category", "FEELnc_lncRNA_class", "candidate_gene_strand", "gene_locus_filter_reason",
    "preASE_retained_for_downstream", "transcript_chr", "transcript_strand", "transcript_exon_count",
    "transcript_span_start0", "transcript_span_end0", "transcript_span_length", "transcript_exonic_length",
    "tes_inferred", "tes_failure_reason", "not_assessed_reason", "tes0", "tes_1bp_start0", "tes_1bp_end0",
    "tes_pm50_start0", "tes_pm50_end0", "tes_pm500_start0", "tes_pm500_end0",
    "threeprime_pm1kb_start0", "threeprime_pm1kb_end0", "coordinate_overlap_group",
    "reference_assessable", "animal_apadb_all_TES_1bp_overlap_count",
    "animal_apadb_all_TES_pm50_overlap_count", "animal_apadb_all_TES_pm500_overlap_count",
    "animal_apadb_brain_CNS_TES_1bp_overlap_count", "animal_apadb_brain_CNS_TES_pm50_overlap_count",
    "animal_apadb_brain_CNS_TES_pm500_overlap_count", "animal_apadb_brain_CNS_plus_pituitary_TES_1bp_overlap_count",
    "animal_apadb_brain_CNS_plus_pituitary_TES_pm50_overlap_count",
    "animal_apadb_brain_CNS_plus_pituitary_TES_pm500_overlap_count", "gse30198_brain_TES_1bp_overlap_count",
    "gse30198_brain_TES_pm50_overlap_count", "gse30198_brain_TES_pm500_overlap_count",
    "gse30198_all_TES_1bp_overlap_count", "gse30198_all_TES_pm50_overlap_count",
    "gse30198_all_TES_pm500_overlap_count", "gse30198_nonbrain_TES_1bp_overlap_count",
    "gse30198_nonbrain_TES_pm50_overlap_count", "gse30198_nonbrain_TES_pm500_overlap_count",
    "pas_motif_scan_available", "pas_motif_scan_failure_reason", "canonical_pas_motif_found",
    "variant_pas_motif_found", "best_pas_motif", "best_pas_motif_class",
    "best_pas_motif_distance_to_TES_bp", "pas_motif_hit_list", "support_sources", "support_category",
    "support_interpretation",
    "external_exact_supported", "external_within_50bp", "external_within_500bp",
    "brain_CNS_within_50bp", "brain_CNS_within_500bp",
    "brain_or_pituitary_within_50bp", "brain_or_pituitary_within_500bp",
    "nonbrain_within_50bp", "nonbrain_within_500bp", "any_pas_motif_found",
    "polyaseq_within_50bp", "polyaseq_within_500bp",
    "animal_apadb_within_50bp", "animal_apadb_within_500bp",
    "source_support_class_50bp", "source_support_class_500bp",
    "terminal_exon_start0", "terminal_exon_end0", "terminal_exon_status",
    "motif_window_crosses_terminal_exon_boundary",
)
MOTIF_FIELDS = (
    "query_transcript_id", "transcript_chr", "transcript_strand", "tes_inferred", "tes0",
    "pas_motif_scan_available", "pas_motif_scan_failure_reason", "motif_scan_window_start0",
    "motif_scan_window_end0", "motif_scan_window_orientation", "motif_scan_sequence_length",
    "canonical_pas_motif_found", "variant_pas_motif_found", "best_pas_motif", "best_pas_motif_class",
    "best_pas_motif_distance_to_TES_bp", "pas_motif_hit_list",
    "terminal_exon_start0", "terminal_exon_end0", "terminal_exon_status",
    "motif_window_crosses_terminal_exon_boundary",
)
QC_FIELDS = ("metric", "value", "note")
MANIFEST_FIELDS = ("item", "path", "exists", "size_bytes", "note")
TRACK_NAMES = (
    "animal_apadb_all", "animal_apadb_brain_CNS", "animal_apadb_brain_CNS_plus_pituitary",
    "gse30198_brain", "gse30198_all", "gse30198_nonbrain",
)


OUTPUT_FIELDS += ("tes_pm100_start0", "tes_pm100_end0", "source_support_class_100bp",
                  "has_relaxed_3prime_evidence") + RADIUS100_FLAGS + SPLICED_FIELDS + tuple(
    f"{name}_TES_pm100_overlap_count" for name in TRACK_NAMES)
MOTIF_FIELDS += SPLICED_FIELDS


@dataclass
class TranscriptModel:
    chroms: Set[str] = field(default_factory=set)
    strands: Set[str] = field(default_factory=set)
    min_start: int | None = None
    max_end: int | None = None
    exon_count: int = 0
    exons: Set[Tuple[int, int]] = field(default_factory=set)


@dataclass(frozen=True)
class IntervalTrack:
    starts: Dict[Tuple[str, str], List[int]]
    ends: Dict[Tuple[str, str], List[int]]
    input_rows: int
    kept_rows: int
    skipped_rows: int


@dataclass(frozen=True)
class FastaIndexRecord:
    length: int
    offset: int
    line_bases: int
    line_width: int


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s ] %(message)s", datefmt="%Y-%m-%d %H:%M:%S", stream=sys.stderr)


def fail(message: str) -> None:
    raise RuntimeError(message)


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return text if text else "NA"


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"{label} does not exist, is not regular file or is empty: {path}")


def to_int(value: str, label: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise RuntimeError(f"{label} is not an integer: {value!r}") from exc


def require_fields(path: Path, fieldnames: Sequence[str] | None, required: Sequence[str]) -> None:
    if fieldnames is None:
        fail(f"{path} missing TSV header")
    missing = [field for field in required if field not in fieldnames]
    if missing:
        fail(f"{path} missing required fields: {','.join(missing)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate TES/polyA support for transcript-first five-dog candidates.")
    parser.add_argument("--retained-transcript-tsv", type=Path, default=DEFAULT_RETAINED_TSV)
    parser.add_argument("--transcript-summary-tsv", type=Path, default=DEFAULT_SUMMARY_TSV)
    parser.add_argument("--transcript-exon-bed", type=Path, default=STEP16_EXON_BED)
    parser.add_argument("--chrom-sizes", type=Path, default=CHROM_SIZES)
    parser.add_argument("--genome-fasta", type=Path, default=GENOME_FASTA)
    parser.add_argument("--animal-apadb-all-bed", type=Path, default=ANIMAL_ROOT / "animal_apadb_dog_PAS_all.CanFam3.bed")
    parser.add_argument("--animal-apadb-brain-cns-bed", type=Path, default=ANIMAL_ROOT / "animal_apadb_dog_PAS_brain_CNS.CanFam3.bed")
    parser.add_argument("--animal-apadb-brain-cns-plus-pituitary-bed", type=Path, default=ANIMAL_ROOT / "animal_apadb_dog_PAS_brain_CNS_plus_pituitary.CanFam3.bed")
    parser.add_argument("--gse30198-brain-bed", type=Path, default=GSE_ROOT / "GSE30198_dog_brain_polyA_sites.CanFam3.bed")
    parser.add_argument("--gse30198-all-bed", type=Path, default=GSE_ROOT / "GSE30198_dog_polyA_sites.CanFam3.bed")
    parser.add_argument("--gse30198-nonbrain-bed", type=Path, default=GSE_ROOT / "GSE30198_dog_nonbrain_polyA_sites.CanFam3.bed")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR / datetime.now().strftime("%Y%m%d_%H%M%S"))

    return parser.parse_args()


def read_retained_ids(path: Path) -> Set[str]:
    retained: Set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, RETAINED_REQUIRED_FIELDS)
        for line_number, row in enumerate(reader, start=2):
            transcript_id = clean(row.get("query_transcript_id"))
            if transcript_id == "NA" or row.get("preASE_retained_for_downstream", "").strip() != "1":
                fail(f"{path}:{line_number} is not a valid retained transcript")
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
            if row.get("preASE_retained_for_downstream", "").strip() != "1":
                fail(f"{path}:{line_number} retained transcript marker is not 1: {transcript_id}")
            summary[transcript_id] = {field: clean(row.get(field)) for field in SUMMARY_REQUIRED_FIELDS}
    missing = retained_ids.difference(summary)
    if missing:
        fail(f"{path} is missing {len(missing)} retained transcripts, for example {sorted(missing)[0]}")
    return summary


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with path.open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not row:
                continue
            if len(row) < 2:
                fail(f"{path}:{line_number} less than two columns")
            sizes[row[0]] = to_int(row[1], f"{path}:{line_number} chromosome length")
    if not sizes:
        fail(f"{path} no chromosome size")
    return sizes


def read_models(path: Path, retained_ids: Set[str], chrom_sizes: Mapping[str, int]) -> Dict[str, TranscriptModel]:
    models: Dict[str, TranscriptModel] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                fail(f"{path}:{line_number} transcript_exon.bed less than 8 column")
            chrom, start_text, end_text, transcript_id, _gene_id, strand = fields[:6]
            if transcript_id not in retained_ids:
                continue
            start, end = to_int(start_text, f"{path}:{line_number} exon start"), to_int(end_text, f"{path}:{line_number} exon end")
            if start < 0 or end <= start or (chrom in chrom_sizes and end > chrom_sizes[chrom]):
                fail(f"{path}:{line_number} exon coordinates or contig invalid")
            model = models.setdefault(transcript_id, TranscriptModel())
            model.chroms.add(chrom)
            model.strands.add(strand)
            model.min_start = start if model.min_start is None else min(model.min_start, start)
            model.max_end = end if model.max_end is None else max(model.max_end, end)
            model.exon_count += 1
            model.exons.add((start, end))
    return models


def read_track(path: Path, label: str, chrom_sizes: Mapping[str, int]) -> IntervalTrack:
    intervals: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    input_rows = kept_rows = skipped_rows = 0
    with path.open(newline="") as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not fields or fields[0].startswith("#") or fields[0].startswith(("track", "browser")):
                continue
            if len(fields) < 6:
                fail(f"{path}:{line_number} {label} BED has fewer than 6 columns; cannot evaluate strand-specific overlap")
            input_rows += 1
            chrom, strand = fields[0], fields[5]
            start, end = to_int(fields[1], f"{path}:{line_number} start"), to_int(fields[2], f"{path}:{line_number} end")
            if start < 0 or end <= start:
                fail(f"{path}:{line_number} BED coordinates invalid: {start}-{end}")
            if strand not in {"+", "-"}:
                fail(f"{path}:{line_number} invalid reference strand: {strand}")
            if end - start != 1:
                fail(f"{path}:{line_number} this analysis requires single-base PAS site")
            if chrom not in chrom_sizes or chrom not in PRIMARY_CHROMS:
                skipped_rows += 1
                continue
            if end > chrom_sizes[chrom]:
                fail(f"{path}:{line_number} BED end exceeds chromosome size: {chrom}:{end}")
            intervals[(chrom, strand)].append((start, end))
            kept_rows += 1
    if not kept_rows:
        fail(f"{path} no available in primary chromosome of strand specific {label} record")
    starts, ends = {}, {}
    for key, values in intervals.items():
        values.sort(key=lambda item: (item[0], item[1]))
        starts[key] = [start for start, _end in values]
        ends[key] = sorted(end for _start, end in values)
    return IntervalTrack(starts, ends, input_rows, kept_rows, skipped_rows)


def count_overlaps(track: IntervalTrack, chrom: str, strand: str, start: int, end: int) -> int:
    starts, ends = track.starts.get((chrom, strand), []), track.ends.get((chrom, strand), [])
    return max(0, bisect_left(starts, end) - bisect_right(ends, start)) if starts else 0


def read_fasta_index(path: Path) -> Dict[str, FastaIndexRecord]:
    index: Dict[str, FastaIndexRecord] = {}
    with path.open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(row) < 5:
                fail(f"{path}:{line_number} FASTA index less than 5 column")
            record = FastaIndexRecord(*(to_int(row[column], f"{path}:{line_number} FAI column {column + 1}") for column in range(1, 5)))
            if record.length <= 0 or record.line_bases <= 0 or record.line_width < record.line_bases:
                fail(f"{path}:{line_number} FASTA index parameter invalid")
            index[row[0]] = record
    return index


def fetch_fasta_sequence(handle: BinaryIO, fai: Mapping[str, FastaIndexRecord], chrom: str, start0: int, end0: int) -> str:
    record = fai.get(chrom)
    if record is None or start0 < 0 or end0 < start0 or end0 > record.length:
        fail(f"FASTA extract coordinates invalid: {chrom}:{start0}-{end0}")
    chunks: List[bytes] = []
    position = start0
    while position < end0:
        line_index, offset_in_line = divmod(position, record.line_bases)
        take = min(record.line_bases - offset_in_line, end0 - position)
        handle.seek(record.offset + line_index * record.line_width + offset_in_line)
        chunks.append(handle.read(take))
        position += take
    return b"".join(chunks).decode("ascii").upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


EVIDENCE_FLAGS = (
    "polyaseq_within_50bp", "polyaseq_within_500bp",
    "animal_apadb_within_50bp", "animal_apadb_within_500bp",
    "external_exact_supported", "external_within_50bp", "external_within_500bp",
    "brain_CNS_within_50bp", "brain_CNS_within_500bp",
    "brain_or_pituitary_within_50bp", "brain_or_pituitary_within_500bp",
    "nonbrain_within_50bp", "nonbrain_within_500bp",
)

EVIDENCE_FLAGS += RADIUS100_FLAGS


def empty_support(base: Dict[str, object], reason: str, group: str = "not_assessed") -> Dict[str, object]:

    base.setdefault("tes_inferred", "0")
    base.setdefault("tes_failure_reason", reason if base["tes_inferred"] == "0" else "NA")
    for key in ("tes0", "tes_1bp_start0", "tes_1bp_end0"):
        base.setdefault(key, "NA")
    base.update({
        "not_assessed_reason": reason,
        "tes_pm50_start0": "NA", "tes_pm50_end0": "NA",
        "tes_pm500_start0": "NA", "tes_pm500_end0": "NA",
        "threeprime_pm1kb_start0": "NA", "threeprime_pm1kb_end0": "NA",
        "coordinate_overlap_group": group, "reference_assessable": "no",
        "pas_motif_scan_available": "0", "pas_motif_scan_failure_reason": reason,
        "canonical_pas_motif_found": "NA", "variant_pas_motif_found": "NA",
        "any_pas_motif_found": "NA", "best_pas_motif": "NA",
        "best_pas_motif_class": "NA", "best_pas_motif_distance_to_TES_bp": "NA",
        "pas_motif_hit_list": "NA", "support_sources": "NA",
        "support_category": "not_assessed", "support_interpretation": reason,
        "source_support_class_50bp": "not_assessed", "source_support_class_500bp": "not_assessed",
        "terminal_exon_start0": "NA", "terminal_exon_end0": "NA",
        "terminal_exon_status": "not_assessed", "motif_window_crosses_terminal_exon_boundary": "NA",
    })
    base.update({key: "NA" for key in SPLICED_FIELDS})
    base.update({"spliced_pas_scan_available": "0", "spliced_pas_scan_failure_reason": reason,
                 "downstream_20nt_assessable": "no", "tes_pm100_start0": "NA",
                 "tes_pm100_end0": "NA", "source_support_class_100bp": "not_assessed",
                 "has_relaxed_3prime_evidence": "NA"})
    for flag in EVIDENCE_FLAGS:
        base[flag] = "NA"
    for track_name in TRACK_NAMES:
        for window in ("TES_1bp", "TES_pm50", "TES_pm100", "TES_pm500"):
            base[f"{track_name}_{window}_overlap_count"] = "NA"
    return base


def derive_motif(handle: BinaryIO, fai: Mapping[str, FastaIndexRecord], chrom: str, strand: str, tes0: int) -> Dict[str, object]:
    chrom_size = fai[chrom].length
    if strand == "+":
        start0, end0, orientation = max(0, tes0 - 50), max(0, tes0 - 9), "genomic_plus"
        sequence = fetch_fasta_sequence(handle, fai, chrom, start0, end0)
    else:
        start0, end0, orientation = min(chrom_size, tes0 + 10), min(chrom_size, tes0 + 51), "reverse_complement"
        sequence = reverse_complement(fetch_fasta_sequence(handle, fai, chrom, start0, end0))


    scan_failure = ("scan_window_shorter_than_6bp" if len(sequence) < 6 else
                    "scan_window_contains_non_ACGT" if set(sequence) - set("ACGT") else "NA")
    if scan_failure != "NA":
        return {
            "pas_motif_scan_available": "0", "pas_motif_scan_failure_reason": scan_failure,
            "motif_scan_window_start0": start0, "motif_scan_window_end0": end0,
            "motif_scan_window_orientation": orientation, "motif_scan_sequence_length": len(sequence),
            "canonical_pas_motif_found": "NA", "variant_pas_motif_found": "NA",
            "best_pas_motif": "NA", "best_pas_motif_class": "NA",
            "best_pas_motif_distance_to_TES_bp": "NA", "pas_motif_hit_list": "NA",
        }
    hits: List[Tuple[int, str, str, int]] = []
    for motif in ALL_PAS_MOTIFS:
        category = "canonical" if motif in CANONICAL_PAS_MOTIFS else "variant"
        offset = sequence.find(motif)
        while offset >= 0:
            distance = tes0 - (start0 + offset + len(motif) - 1) if strand == "+" else end0 - offset - len(motif) - tes0
            hits.append((0 if category == "canonical" else 1, motif, category, distance))
            offset = sequence.find(motif, offset + 1)
    hits.sort(key=lambda item: (item[0], abs(item[3] - 20), item[1], item[3]))
    best = hits[0] if hits else None
    return {
        "pas_motif_scan_available": "1", "pas_motif_scan_failure_reason": "NA", "motif_scan_window_start0": start0,
        "motif_scan_window_end0": end0, "motif_scan_window_orientation": orientation, "motif_scan_sequence_length": len(sequence),
        "canonical_pas_motif_found": "1" if any(hit[2] == "canonical" for hit in hits) else "0",
        "variant_pas_motif_found": "1" if any(hit[2] == "variant" for hit in hits) else "0",
        "best_pas_motif": best[1] if best else "NA", "best_pas_motif_class": best[2] if best else "NA",
        "best_pas_motif_distance_to_TES_bp": best[3] if best else "NA",
        "pas_motif_hit_list": "|".join(f"{motif}:{category}:{distance}bp" for _rank, motif, category, distance in hits) or "NA",
    }


def derive_spliced_motif(handle: BinaryIO, fai: Mapping[str, FastaIndexRecord],
                         chrom: str, strand: str, tes0: int,
                         model: TranscriptModel) -> Dict[str, object]:

    result: Dict[str, object] = {key: "NA" for key in SPLICED_FIELDS}
    result.update({"spliced_pas_scan_available": "0", "downstream_20nt_assessable": "no"})

    start, end = ((tes0 + 1, min(fai[chrom].length, tes0 + 21)) if strand == "+"
                  else (max(0, tes0 - 20), tes0))
    downstream = fetch_fasta_sequence(handle, fai, chrom, start, end)
    if strand == "-":
        downstream = reverse_complement(downstream)
    result["downstream_20nt_sequence"] = downstream or "NA"
    if len(downstream) == 20 and not set(downstream) - set("ACGT"):
        result.update({"downstream_20nt_assessable": "yes",
                       "downstream_20nt_A_percent": 100 * downstream.count("A") / 20})
    exons = sorted(model.exons)

    if any(left[1] > right[0] for left, right in zip(exons, exons[1:])):
        result["spliced_pas_scan_failure_reason"] = "overlapping_exon_model"
        return result
    chunks, remaining = [], 50

    for start, end in (reversed(exons) if strand == "+" else exons):
        take = min(remaining, end - start)
        fragment = fetch_fasta_sequence(handle, fai, chrom,
                                        end - take if strand == "+" else start,
                                        end if strand == "+" else start + take)
        chunks.append(fragment if strand == "+" else reverse_complement(fragment))
        remaining -= take
        if remaining == 0:
            break
    sequence = "".join(reversed(chunks))
    result.update({"spliced_terminal_50nt_sequence": sequence or "NA",
                   "spliced_terminal_50nt_length": len(sequence)})
    reason = ("scan_window_shorter_than_6nt" if len(sequence) < 6 else
              "scan_window_contains_non_ACGT" if set(sequence) - set("ACGT") else "NA")
    result["spliced_pas_scan_failure_reason"] = reason
    if reason != "NA":
        return result
    hits = []
    for motif in EXTENDED_PAS_MOTIFS:
        offset = sequence.find(motif)
        while offset >= 0:
            distance = len(sequence) - offset - len(motif)
            hits.append((motif, distance))
            offset = sequence.find(motif, offset + 1)
    hits.sort(key=lambda hit: (hit[0] not in CANONICAL_PAS_MOTIFS,
                              abs(hit[1] - 20), hit[0], hit[1]))
    result.update({
        "spliced_pas_scan_available": "1",
        "spliced_pas_selected9_found": "yes" if any(m in ALL_PAS_MOTIFS for m, d in hits) else "no",
        "spliced_pas_extended11_found": "yes" if hits else "no",
        "spliced_pas_canonical_found": "yes" if any(m in CANONICAL_PAS_MOTIFS for m, d in hits) else "no",
        "spliced_pas_variant_found": "yes" if any(m not in CANONICAL_PAS_MOTIFS for m, d in hits) else "no",
        "spliced_pas_hit_list": "|".join(f"{m}:{d}nt" for m, d in hits) or "NA",
        "spliced_pas_best_motif": hits[0][0] if hits else "NA",
        "spliced_pas_best_distance_to_TES_nt": hits[0][1] if hits else "NA",
    })
    return result


def annotate_evidence(row: Dict[str, object]) -> None:
    def supported(names: Sequence[str], window: str) -> str:
        return "yes" if any(int(row[f"{name}_{window}_overlap_count"]) > 0 for name in names) else "no"


    all_names = ("animal_apadb_all", "gse30198_all")
    brain_names = ("animal_apadb_brain_CNS", "gse30198_brain")
    brain_pituitary_names = ("animal_apadb_brain_CNS_plus_pituitary", "gse30198_brain")
    row["external_exact_supported"] = supported(all_names, "TES_1bp")
    for radius in (50, 100, 500):
        window = f"TES_pm{radius}"

        polyaseq = supported(("gse30198_all",), window)
        animal = supported(("animal_apadb_all",), window)
        row[f"polyaseq_within_{radius}bp"] = polyaseq
        row[f"animal_apadb_within_{radius}bp"] = animal
        row[f"source_support_class_{radius}bp"] = (
            "both_sources" if polyaseq == "yes" and animal == "yes" else
            "polyaseq_only" if polyaseq == "yes" else
            "animal_apadb_only" if animal == "yes" else "neither")
        row[f"external_within_{radius}bp"] = supported(all_names, window)
        row[f"brain_CNS_within_{radius}bp"] = supported(brain_names, window)
        row[f"brain_or_pituitary_within_{radius}bp"] = supported(brain_pituitary_names, window)
        row[f"nonbrain_within_{radius}bp"] = supported(("gse30198_nonbrain",), window)
    motif_ok = row["pas_motif_scan_available"] == "1"
    motif = row["canonical_pas_motif_found"] == "1" or row["variant_pas_motif_found"] == "1"
    row["any_pas_motif_found"] = ("yes" if motif else "no") if motif_ok else "NA"
    sources = [name for name in TRACK_NAMES if int(row[f"{name}_TES_pm500_overlap_count"]) > 0]
    if row["canonical_pas_motif_found"] == "1":
        sources.append("canonical_PAS_motif")
    if row["variant_pas_motif_found"] == "1":
        sources.append("variant_PAS_motif")
    row["support_sources"] = "|".join(sources) or "NA"

    if row["external_exact_supported"] == "yes":
        category = "external_exact"
    elif row["external_within_50bp"] == "yes":
        category = "external_1_to_50bp"
    elif row["external_within_100bp"] == "yes":
        category = "external_51_to_100bp"
    elif row["external_within_500bp"] == "yes":
        category = "external_101_to_500bp"
    elif not motif_ok:
        category = "no_external_within_500bp_motif_not_assessed"
    elif motif:
        category = "motif_only_within_500bp_policy"
    else:
        category = "no_selected_support"


    relaxed_motif = row["spliced_pas_extended11_found"]
    row["has_relaxed_3prime_evidence"] = (
        "yes" if row["external_within_500bp"] == "yes" or relaxed_motif == "yes"
        else "no" if relaxed_motif == "no" else "NA")
    row["support_category"] = category
    row["support_interpretation"] = "Descriptive evidence from same-strand site distances and PAS motifs; this does not validate actual cleavage sites."


def derive_support(transcript_id: str, summary: Mapping[str, str], model: TranscriptModel | None, chrom_sizes: Mapping[str, int], tracks: Mapping[str, IntervalTrack], fasta_handle: BinaryIO, fai: Mapping[str, FastaIndexRecord]) -> Dict[str, object]:
    row: Dict[str, object] = {field: summary[field] for field in SUMMARY_REQUIRED_FIELDS if field != "chr"}
    row["query_transcript_id"] = transcript_id
    if model is None:
        row.update({"transcript_chr": "NA", "transcript_strand": "NA", "transcript_exon_count": 0, "transcript_span_start0": "NA", "transcript_span_end0": "NA", "transcript_span_length": "NA"})
        return empty_support(row, "transcript_exon_bed_missing")
    if len(model.chroms) != 1 or len(model.strands) != 1 or model.min_start is None or model.max_end is None:
        row.update({"transcript_chr": "|".join(sorted(model.chroms)), "transcript_strand": "|".join(sorted(model.strands)), "transcript_exon_count": model.exon_count, "transcript_span_start0": model.min_start if model.min_start is not None else "NA", "transcript_span_end0": model.max_end if model.max_end is not None else "NA", "transcript_span_length": "NA"})
        return empty_support(row, "inconsistent_transcript_exon_model")
    chrom, strand = next(iter(model.chroms)), next(iter(model.strands))
    row.update({"transcript_chr": chrom, "transcript_strand": strand, "transcript_exon_count": model.exon_count, "transcript_span_start0": model.min_start, "transcript_span_end0": model.max_end, "transcript_span_length": model.max_end - model.min_start})
    if strand not in {"+", "-"}:
        return empty_support(row, "invalid_transcript_strand")
    tes0 = model.max_end - 1 if strand == "+" else model.min_start
    row.update({"tes_inferred": "1", "tes_failure_reason": "NA", "tes0": tes0,
                "tes_1bp_start0": tes0, "tes_1bp_end0": tes0 + 1})
    if chrom not in chrom_sizes:
        return empty_support(row, "missing_chromosome_size", "missing_chromosome_size")
    if chrom not in fai:
        return empty_support(row, "missing_fasta_contig", "missing_fasta_contig")
    if chrom not in PRIMARY_CHROMS:
        return empty_support(row, "non_primary_chromosome_not_assessed", "alt_or_unplaced_contig")
    chrom_size = chrom_sizes[chrom]
    if not 0 <= tes0 < chrom_size:
        fail(f"{transcript_id} TES exceeds chromosome span: {chrom}:{tes0}")
    windows = {"TES_1bp": (tes0, tes0 + 1), "TES_pm50": (max(0, tes0 - 50), min(chrom_size, tes0 + 51)), "TES_pm500": (max(0, tes0 - 500), min(chrom_size, tes0 + 501))}
    row.update({"tes_inferred": "1", "tes_failure_reason": "NA", "tes0": tes0, "tes_1bp_start0": windows["TES_1bp"][0], "tes_1bp_end0": windows["TES_1bp"][1], "tes_pm50_start0": windows["TES_pm50"][0], "tes_pm50_end0": windows["TES_pm50"][1], "tes_pm500_start0": windows["TES_pm500"][0], "tes_pm500_end0": windows["TES_pm500"][1], "threeprime_pm1kb_start0": max(0, tes0 - 1000), "threeprime_pm1kb_end0": min(chrom_size, tes0 + 1001), "coordinate_overlap_group": "primary_chr", "reference_assessable": "yes", "not_assessed_reason": "NA"})
    windows["TES_pm100"] = (max(0, tes0 - 100), min(chrom_size, tes0 + 101))
    row.update({"tes_pm100_start0": windows["TES_pm100"][0],
                "tes_pm100_end0": windows["TES_pm100"][1]})
    for track_name, track in tracks.items():
        for window_name, (start0, end0) in windows.items():
            row[f"{track_name}_{window_name}_overlap_count"] = count_overlaps(track, chrom, strand, start0, end0)
    row.update(derive_motif(fasta_handle, fai, chrom, strand, tes0))
    row.update(derive_spliced_motif(fasta_handle, fai, chrom, strand, tes0, model))


    terminal_exons = [exon for exon in model.exons
                      if (exon[1] == model.max_end if strand == "+" else exon[0] == model.min_start)]
    row.update({"terminal_exon_start0": "NA", "terminal_exon_end0": "NA",
                "terminal_exon_status": "ambiguous_terminal_exon",
                "motif_window_crosses_terminal_exon_boundary": "NA"})
    if len(terminal_exons) == 1:
        terminal_start, terminal_end = terminal_exons[0]
        contained = (terminal_start <= row["motif_scan_window_start0"] and
                     row["motif_scan_window_end0"] <= terminal_end)
        row.update({"terminal_exon_start0": terminal_start, "terminal_exon_end0": terminal_end,
                    "terminal_exon_status": "unique_terminal_exon",
                    "motif_window_crosses_terminal_exon_boundary": "no" if contained else "yes"})
    annotate_evidence(row)
    return row


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean(row.get(field)) for field in fields})


def write_bed(path: Path, rows: Iterable[Sequence[object]]) -> int:
    count = 0
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def write_evidence_summaries(outputs: Mapping[str, Path], rows: Sequence[Mapping[str, object]]) -> None:
    reference_rows = [r for r in rows if r["reference_assessable"] == "yes"]
    motif_rows = [r for r in rows if r["pas_motif_scan_available"] == "1"]
    joint_rows = [r for r in reference_rows if r["pas_motif_scan_available"] == "1"]
    summary = []
    for field in EVIDENCE_FLAGS + ("canonical_pas_motif_found", "variant_pas_motif_found", "any_pas_motif_found"):
        motif_field = field in ("canonical_pas_motif_found", "variant_pas_motif_found", "any_pas_motif_found")
        eligible = motif_rows if motif_field else reference_rows
        positive = "1" if field in ("canonical_pas_motif_found", "variant_pas_motif_found") else "yes"
        n = sum(r[field] == positive for r in eligible)
        summary.append({"metric": field, "supported_transcripts": n,
                        "assessable_transcripts": len(eligible), "total_transcripts": len(rows),
                        "not_assessed": len(rows)-len(eligible),
                        "percent": 100*n/len(eligible) if eligible else "NA"})
    for field in ("spliced_pas_selected9_found", "spliced_pas_extended11_found",
                  "spliced_pas_canonical_found", "spliced_pas_variant_found",
                  "has_relaxed_3prime_evidence"):
        eligible = [r for r in rows if r[field] != "NA"]
        n = sum(r[field] == "yes" for r in eligible)
        summary.append({"metric": field, "supported_transcripts": n,
                        "assessable_transcripts": len(eligible), "total_transcripts": len(rows),
                        "not_assessed": len(rows) - len(eligible),
                        "percent": 100*n/len(eligible) if eligible else "NA"})
    write_tsv(outputs["evidence_summary"],
              ("metric", "supported_transcripts", "assessable_transcripts", "total_transcripts", "not_assessed", "percent"), summary)
    joint = []
    for threshold, field in (("exact", "external_exact_supported"), ("within_50bp", "external_within_50bp"),
                             ("within_100bp", "external_within_100bp"), ("within_500bp", "external_within_500bp")):
        for ext, motif, category in ((True, True, "both"), (True, False, "external_only"),
                                     (False, True, "motif_only"), (False, False, "neither")):
            n = sum((r[field] == "yes") == ext and (r["any_pas_motif_found"] == "yes") == motif for r in joint_rows)
            joint.append({"external_threshold": threshold, "category": category,
                          "transcripts": n, "joint_assessable_transcripts": len(joint_rows),
                          "not_jointly_assessed": len(rows)-len(joint_rows),
                          "percent": 100*n/len(joint_rows) if joint_rows else "NA"})
    write_tsv(outputs["joint_summary"],
              ("external_threshold", "category", "transcripts", "joint_assessable_transcripts", "not_jointly_assessed", "percent"), joint)
    spliced_joint = []
    eligible = [r for r in reference_rows if r["spliced_pas_scan_available"] == "1"]
    for radius in (50, 100, 500):
        for ext, motif, category in ((True, True, "both"), (True, False, "external_only"),
                                     (False, True, "motif_only"), (False, False, "neither")):
            n = sum((r[f"external_within_{radius}bp"] == "yes") == ext and
                    (r["spliced_pas_extended11_found"] == "yes") == motif for r in eligible)
            spliced_joint.append({"external_threshold_bp": radius,
                "motif_definition": "spliced_terminal_50nt_extended11", "category": category,
                "transcripts": n, "joint_assessable_transcripts": len(eligible),
                "not_jointly_assessed": len(rows)-len(eligible),
                "percent": 100*n/len(eligible) if eligible else "NA"})
    write_tsv(outputs["spliced_joint_summary"],
              ("external_threshold_bp", "motif_definition", "category", "transcripts",
               "joint_assessable_transcripts", "not_jointly_assessed", "percent"), spliced_joint)

    source_summary = []
    for radius in (50, 100, 500):
        for category in ("both_sources", "polyaseq_only", "animal_apadb_only", "neither"):
            n = sum(r[f"source_support_class_{radius}bp"] == category for r in reference_rows)
            source_summary.append({"distance_threshold_bp": radius, "category": category,
                                   "transcripts": n, "assessable_transcripts": len(reference_rows),
                                   "not_assessed": len(rows)-len(reference_rows),
                                   "percent": 100*n/len(reference_rows) if reference_rows else "NA"})
    write_tsv(outputs["source_summary"],
              ("distance_threshold_bp", "category", "transcripts", "assessable_transcripts", "not_assessed", "percent"), source_summary)
    outputs["policy"].write_text(
        "schema_version=3\n"
        "Unit=retained transcript; original query IDs are retained, not gene-level counts.\n"
        "End=inferred strand-specific outermost exon boundary, not a validated cleavage site.\n"
        "Reference scope=chr1-38 and chrX; inferred ends on other contigs are retained but not assessed.\n"
        "TES 1-bp BED includes all inferred endpoints; window BEDs include only reference-assessable transcripts.\n"
        "External union=Animal-APAdb all OR GSE30198 all; same strand; 1-bp reference sites.\n"
        "Distances=exact, <=50 bp, <=100 bp, <=500 bp; cumulative. Primary relaxed proximity threshold=500 bp; 50 bp close support and 100 bp sensitivity are retained; not precise TES validation.\n"
        "polyaseq_within_50bp/100bp/500bp=GSE30198 all; animal_apadb_within_50bp/100bp/500bp=Animal-APAdb all; derived from existing same-strand overlap counts.\n"
        "source_support_class_50bp/100bp/500bp=both_sources|polyaseq_only|animal_apadb_only|neither; unassessed=not_assessed.\n"
        "Both sources may hit different sites in a window; their intersection is not two independent experimental validations of one cleavage site.\n"
        "GSE30198 is a dedicated PolyA-seq experimental reference; Animal-APAdb representative sites are primarily RNA-seq/DaPars2 predictions, with DAPA/QAPA usage annotations. External union is not uniformly experimental evidence.\n"
        "Tissue flags are separate from distance. Brain/CNS excludes pituitary; brain_or_pituitary includes it.\n"
        "GSE30198 nonbrain=selected kidney/testis data, not all nonbrain tissues. Nested tracks are not independent assays.\n"
        "Animal-APAdb=CanFam3.1 representative PAS with selected DAPA/QAPA tissue usage annotations.\n"
        "GSE30198=PolyA-Seq, canFam2 lifted to CanFam3.1 upstream according to reference provenance files.\n"
        "Motif=genomic sequence at relative bases -50..-10, reverse complemented on minus strand; full hexamer must fit.\n"
        "Motif distance=distance of motif downstream-most base to inferred TES. Best=canonical first, then closest to 20 bp.\n"
        "Canonical=AATAAA|ATTAAA; variants=TATAAA|AGTAAA|AAGAAA|AATATA|AATACA|CATAAA|GATAAA.\n"
        "Motif scan is primary-contig-only, not exon-spliced; short terminal exons may include intronic sequence.\n"
        "motif_window_crosses_terminal_exon_boundary=yes when the actual genomic scan interval is not wholly contained in the unique terminal exon; no otherwise; NA if unassessed or terminal exon ambiguous.\n"
        "Terminal exon is the distinct exon ending at TES on +, or starting at TES on -; multiple distinct boundary choices are ambiguous. The boundary flag does not filter candidates or change motif detection.\n"
        "Genomic motif fields retain the selected9 panel and genomic -50..-10 window. Spliced fields independently scan the final up-to-50 transcript nt, including TES base at distance 0; negative strand is reverse-complemented; overlapping exon models are unassessed. Extended11=selected9+AATGAA+AATAGA (Beaudoing 2000); not an exhaustive or dog-validated panel.\n"
        "Unscanned motif=NA; genomic canonical/variant use 1/0, genomic any and spliced motif flags use yes/no. Fewer than six available bases is not assessable.\n"
        "Windows containing any non-ACGT base are excluded in full from motif assessment, including locally recognizable hits; external evidence is unaffected.\n"
        "Joint denominator=both reference and motif assessable; counts are provided separately at each external threshold.\n"
        "support_category=nearest distance band of external evidence, then motif-only within the 500-bp policy; not a validation grade.\n"
        "support_sources=tracks hit within 500 bp and all motif classes found; not an independent-assay count.\n"
        "Downstream20 A percent is reported in transcript orientation only for 20 ACGT bases; no internal-priming filter is applied. Motif presence or external proximity does not establish full-length structure.\n"
        "has_relaxed_3prime_evidence=external_within_500bp OR spliced_pas_extended11_found; positive if either positive, negative only if both assessable and negative, otherwise NA. Denominator may differ from complete-case joint summaries.\n"
        "support_category/support_sources and the original joint table use genomic selected9 motifs; spliced joint table is separate. No candidate is removed.\n"
        "METHOD REFERENCES:\n"
        "SQANTI2 official README https://github.com/Magdoll/SQANTI2: historical within_polya_site uses 100 bp; this is a sensitivity threshold, not a claim that all current SQANTI3 versions use identical semantics.\n"
        "Singh et al. 2018, Nature Communications, doi:10.1038/s41467-018-04112-z: 500 nt associates assembled structures with experimental IpA ends in a constrained workflow; this relaxed dog proximity annotation is an adaptation, not exact-end validation or a reproduction of all their filters.\n"
        "Derti et al. 2012, Genome Research, doi:10.1101/gr.132563.111 (GSE30198): dedicated PolyA-seq with upstream internal-priming filtering; same-strand peaks clustered within 30 bases. This does not validate each candidate endpoint.\n"
        "Jin et al. 2021, NAR, doi:10.1093/nar/gkaa778: Animal-APAdb uses RNA-seq/DaPars2 and cluster-median representative sites; not a dedicated 3prime sequencing experiment.\n"
        "Xia et al. 2014, Nature Communications, doi:10.1038/ncomms6274: DaPars simulation evaluation uses a 50-bp positional tolerance with an additional usage criterion; supports a close-distance convention, not validation of this workflow.\n"
        "Beaudoing et al. 2000, Genome Research, doi:10.1101/gr.10.7.1001: canonical and variant polyadenylation signals; the selected 9 motifs are a subset, not an exact reproduction of the published panel.\n"
        "Pardo-Palacios et al. 2024, Nature Methods, doi:10.1038/s41592-024-02229-2: SQANTI3 examines motifs in the final 50 bp of transcript sequence; this genomic -50..-10 scan is different.\n"
        "Summary and QC tables are descriptive; no postflight biological acceptance test is part of execution status.\n",
        encoding="utf-8")


def main() -> None:
    setup_logging()
    args = parse_args()
    fasta_fai = Path(str(args.genome_fasta) + ".fai")
    inputs = (("retained_transcript_tsv", args.retained_transcript_tsv), ("transcript_summary_tsv", args.transcript_summary_tsv), ("transcript_exon_bed", args.transcript_exon_bed), ("chrom_sizes", args.chrom_sizes), ("genome_fasta", args.genome_fasta), ("genome_fasta_index", fasta_fai), ("animal_apadb_all_bed", args.animal_apadb_all_bed), ("animal_apadb_brain_cns_bed", args.animal_apadb_brain_cns_bed), ("animal_apadb_brain_cns_plus_pituitary_bed", args.animal_apadb_brain_cns_plus_pituitary_bed), ("gse30198_brain_bed", args.gse30198_brain_bed), ("gse30198_all_bed", args.gse30198_all_bed), ("gse30198_nonbrain_bed", args.gse30198_nonbrain_bed))
    for label, path in inputs:
        require_file(path, label)
    outputs = {"support": args.outdir / "all_candidate_transcript_TES_polyA_support.tsv", "motif": args.outdir / "candidate_TES_PAS_motif_support.tsv", "tes_1bp": args.outdir / "all_candidate_TES_1bp.bed", "tes_pm50": args.outdir / "all_candidate_TES_pm50bp.bed", "tes_pm100": args.outdir / "all_candidate_TES_pm100bp.bed", "tes_pm500": args.outdir / "all_candidate_TES_pm500bp.bed", "threeprime": args.outdir / "all_candidate_3prime_pm1kb.bed", "qc": args.outdir / "TES_polyA_support_qc_summary.tsv", "manifest": args.outdir / "TES_polyA_support_manifest.tsv"}
    outputs.update({"evidence_summary": args.outdir / "TES_polyA_evidence_summary.tsv",
                    "source_summary": args.outdir / "TES_polyA_source_support_summary.tsv",
                    "spliced_joint_summary": args.outdir / "TES_polyA_spliced_joint_support.tsv",
                    "joint_summary": args.outdir / "TES_polyA_joint_support.tsv",
                    "policy": args.outdir / "TES_polyA_support_policy.txt"})
    if not args.outdir.is_absolute():
        fail("--outdir must be an absolute path")
    args.outdir.mkdir(parents=True, exist_ok=False)
    logging.info("concurrent_tasks=1 threads_per_task=1(serial settings retained)")
    logging.info("[stage 1/5] read and validate Step 2.1 retained transcript input")
    retained_ids = read_retained_ids(args.retained_transcript_tsv)
    summary = read_summary(args.transcript_summary_tsv, retained_ids)
    logging.info("retained transcripts=%d", len(retained_ids))
    logging.info("[stage 2/5] read exon model, CanFam3.1 sizes with FASTA index")
    chrom_sizes, fasta_index = read_chrom_sizes(args.chrom_sizes), read_fasta_index(fasta_fai)
    models = read_models(args.transcript_exon_bed, retained_ids, chrom_sizes)
    logging.info("[stage 3/5] read dog PAS/polyA strand specific track")
    tracks = {"animal_apadb_all": read_track(args.animal_apadb_all_bed, "Animal-APAdb all", chrom_sizes), "animal_apadb_brain_CNS": read_track(args.animal_apadb_brain_cns_bed, "Animal-APAdb brain/CNS", chrom_sizes), "animal_apadb_brain_CNS_plus_pituitary": read_track(args.animal_apadb_brain_cns_plus_pituitary_bed, "Animal-APAdb brain/CNS plus pituitary", chrom_sizes), "gse30198_brain": read_track(args.gse30198_brain_bed, "GSE30198 brain", chrom_sizes), "gse30198_all": read_track(args.gse30198_all_bed, "GSE30198 all", chrom_sizes), "gse30198_nonbrain": read_track(args.gse30198_nonbrain_bed, "GSE30198 nonbrain", chrom_sizes)}
    logging.info("[stage 4/5] infer strand specific TES, respectively annotation PolyA-seq/Animal-APAdb; respectively scan genome/splicing terminal PAS, and annotation downstream A content")
    rows: List[Dict[str, object]] = []
    with args.genome_fasta.open("rb") as fasta_handle:
        for index, transcript_id in enumerate(sorted(retained_ids), start=1):
            rows.append(derive_support(transcript_id, summary[transcript_id], models.get(transcript_id), chrom_sizes, tracks, fasta_handle, fasta_index))
            if index % 1000 == 0:
                logging.info("completed %d/%d transcripts", index, len(retained_ids))
    logging.info("[stage 5/5] write TES/polyA table, BED, QC with manifest")
    write_tsv(outputs["support"], OUTPUT_FIELDS, rows)
    write_tsv(outputs["motif"], MOTIF_FIELDS, rows)
    bed_rows = {"tes_1bp": [], "tes_pm50": [], "tes_pm100": [], "tes_pm500": [], "threeprime": []}
    for row in rows:
        if row["tes_inferred"] != "1":
            continue
        shared = (row["transcript_chr"], row["query_transcript_id"], 0, row["transcript_strand"])
        for output_key, start_field, end_field in (("tes_1bp", "tes_1bp_start0", "tes_1bp_end0"), ("tes_pm50", "tes_pm50_start0", "tes_pm50_end0"), ("tes_pm100", "tes_pm100_start0", "tes_pm100_end0"), ("tes_pm500", "tes_pm500_start0", "tes_pm500_end0"), ("threeprime", "threeprime_pm1kb_start0", "threeprime_pm1kb_end0")):
            if output_key != "tes_1bp" and row["reference_assessable"] != "yes":
                continue
            bed_rows[output_key].append((shared[0], row[start_field], row[end_field], shared[1], shared[2], shared[3]))
    bed_counts = {key: write_bed(outputs[key], records) for key, records in bed_rows.items()}
    level_counts = Counter(clean(row["support_category"]) for row in rows)
    qc_rows: List[Dict[str, object]] = [{"metric": "retained_transcript_count", "value": len(retained_ids), "note": "Number of Step 2.1 retained query_transcript_id values."}, {"metric": "summary_transcript_count", "value": len(summary), "note": "Number of retained transcripts in the Step 2.1 summary."}, {"metric": "exon_model_count", "value": len(models), "note": "Number of models in transcript_exon.bed matching retained transcripts."}, {"metric": "tes_inferred_count", "value": sum(row["tes_inferred"] == "1" for row in rows), "note": "Number of transcripts with successfully inferred strand-aware TES."}, {"metric": "tes_unavailable_count", "value": sum(row["tes_inferred"] != "1" for row in rows), "note": "Number of transcripts without an inferred TES."}, {"metric": "output_transcript_rows", "value": len(rows), "note": "Number of output rows, one per retained query_transcript_id."}]
    qc_rows.extend({"metric": f"support_category_{level}_count", "value": count, "note": "transcript-level TES/polyA support"} for level, count in sorted(level_counts.items()))
    qc_rows.extend({"metric": f"{key}_rows", "value": count, "note": "Number of TES BED rows."} for key, count in sorted(bed_counts.items()))
    qc_rows.extend({"metric": f"{name}_{metric}", "value": value, "note": "PAS/polyA BED reading statistics."} for name, track in tracks.items() for metric, value in (("input_rows", track.input_rows), ("kept_primary_rows", track.kept_rows), ("skipped_rows", track.skipped_rows)))
    qc_rows.extend([
        {"metric": "reference_assessable_count", "value": sum(r["reference_assessable"] == "yes" for r in rows), "note": "Evaluability against reference tracks, distinct from whether endpoints can be inferred."},
        {"metric": "motif_assessable_count", "value": sum(r["pas_motif_scan_available"] == "1" for r in rows), "note": "Denominator for the motif metric itself."},
        {"metric": "spliced_motif_assessable_count", "value": sum(r["spliced_pas_scan_available"] == "1" for r in rows), "note": "Separate denominator for motifs in the terminal 50 nt of spliced transcripts."},
        {"metric": "downstream_20nt_assessable_count", "value": sum(r["downstream_20nt_assessable"] == "yes" for r in rows), "note": "Complete downstream 20 nt containing only ACGT; no A-content filter is applied."},
    ])
    for boundary_state in ("yes", "no", "NA"):
        qc_rows.append({"metric": f"motif_window_crosses_terminal_exon_boundary_{boundary_state}",
                        "value": sum(r["motif_window_crosses_terminal_exon_boundary"] == boundary_state for r in rows),
                        "note": "Whether the actual genomic scanning window extends beyond the terminal exon; a descriptive count, not a candidate filter."})

    for field in ("not_assessed_reason", "pas_motif_scan_failure_reason", "spliced_pas_scan_failure_reason"):
        reasons = Counter(str(r[field]) for r in rows if r[field] != "NA")
        qc_rows.extend({"metric": f"{field}:{reason}", "value": count,
                        "note": "Reason for non-evaluability; this does not imply absence of biological support."}
                       for reason, count in sorted(reasons.items()))
    write_tsv(outputs["qc"], QC_FIELDS, qc_rows)
    write_evidence_summaries(outputs, rows)
    end_time = datetime.now().strftime("%F %T")
    manifest_rows = [{"item": item, "path": str(path), "exists": int(path.exists()), "size_bytes": path.stat().st_size if path.exists() else "NA", "note": "input"} for item, path in inputs]

    manifest_rows.extend({"item": item, "path": str(path), "exists": int(path.exists()), "size_bytes": path.stat().st_size if path.exists() else "NA", "note": "output"} for item, path in outputs.items() if item != "manifest")
    manifest_rows.extend(({"item": "canonical_pas_motifs", "path": "NA", "exists": "NA", "size_bytes": "NA", "note": "|".join(CANONICAL_PAS_MOTIFS)}, {"item": "variant_pas_motifs", "path": "NA", "exists": "NA", "size_bytes": "NA", "note": "|".join(VARIANT_PAS_MOTIFS)}, {"item": "run_end_time", "path": "NA", "exists": "NA", "size_bytes": "NA", "note": end_time}))
    manifest_rows.extend([
        {"item": "relaxed_proximity_bp", "path": "NA", "exists": "NA", "size_bytes": "NA", "note": "500"},
        {"item": "spliced_extended11_motifs", "path": "NA", "exists": "NA", "size_bytes": "NA", "note": "|".join(EXTENDED_PAS_MOTIFS)},
        {"item": "motif_extension_source", "path": "https://doi.org/10.1101/gr.10.7.1001", "exists": "NA", "size_bytes": "NA", "note": "Selected9 plus AATGAA and AATAGA; fixed selected11-v1 panel"},
    ])
    write_tsv(outputs["manifest"], MANIFEST_FIELDS, manifest_rows)
    logging.info("main program compute completed: %d transcript rows -> %s", len(rows), outputs["support"])
    logging.info("result completeness, statistical denominator and biological quality requires separately check; from batch based on exit status record execution completed.")


if __name__ == "__main__":
    main()
