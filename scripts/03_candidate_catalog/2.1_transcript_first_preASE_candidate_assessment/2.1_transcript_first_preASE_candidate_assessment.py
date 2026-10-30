#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple


PROJECT_ROOT = Path("/path/to/project/analysis_20260709/5-dogs/9.new-gene")
INPUT_ROOT = Path("/path/to/project/legacy/five-dogs/9.new-gene")
REFERENCE_ROOT = Path("/path/to/project/legacy/reference")
STEP16_ROOT = INPUT_ROOT / "1.6_gffcompare_stringtie3_merged_transcriptome-v6"
STEP165C_ROOT = INPUT_ROOT / "1.6.5_feelnc_transcript_annotation/1.6.5c_run_feelnc_on_all_query-v4"
DEFAULT_CLASSIFICATION_TSV = STEP16_ROOT / "1.6_gffcompare_transcript_classification.tsv"
DEFAULT_TRANSCRIPT_EXON_BED = STEP16_ROOT / "transcript_exon.bed"
DEFAULT_FEELNC_STATUS_TSV = STEP165C_ROOT / "qc/1.6.5c_all_query_feelnc_module_status.tsv"
DEFAULT_REFERENCE_GTF = REFERENCE_ROOT / "GCF_000002285.3_CanFam3.1_genomic_chr.gtf"
DEFAULT_OUTDIR = PROJECT_ROOT / "2.1_transcript_first_preASE_candidate_assessment"
REFERENCE_BIN_SIZE = 100000
DEFAULT_CONTAINMENT_OVERLAP_FRACTION = 0.95

CLASSIFICATION_REQUIRED_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "reference_gene_id",
    "reference_transcript_id",
    "class_code",
    "novelty_category",
    "num_exons",
    "FPKM",
    "TPM",
    "cov",
    "length",
    "major_iso_id",
    "ref_match_len",
    "tmap_row_number",
]

FEELNC_REQUIRED_FIELDS = [
    "query_transcript_id",
    "query_gene_id",
    "class_code",
    "transcript_length",
    "exon_count",
    "single_exon_flag",
    "FEELnc_filter_status",
    "FEELnc_CPS",
    "FEELnc_RF_label",
    "FEELnc_codpot_category",
    "FEELnc_noORF_flag",
    "FEELnc_biotype",
    "FEELnc_classifier_best_isBest",
    "FEELnc_lncRNA_class",
    "FEELnc_lncRNA_subtype",
    "FEELnc_lncRNA_location",
    "FEELnc_partner_gene",
    "FEELnc_partner_transcript",
    "FEELnc_partner_direction",
    "FEELnc_partner_distance",
    "FEELnc_consistency_flags",
]

TRANSCRIPT_GENE_FILTER_FIELDS = [
    "preASE_retained_by_class_code",
    "candidate_gene_chr",
    "candidate_gene_start0",
    "candidate_gene_end0",
    "candidate_gene_span_length",
    "candidate_gene_strand",
    "candidate_gene_transcript_count",
    "candidate_gene_same_strand_contains_known_gene_count",
    "candidate_gene_same_strand_within_known_gene_count",
    "candidate_gene_same_strand_equal_known_gene_count",
    "candidate_gene_same_strand_partial_overlap_known_gene_count",
    "candidate_gene_same_strand_overlap_known_gene_count",
    "candidate_gene_same_strand_contains_known_gene_id_list",
    "candidate_gene_same_strand_within_known_gene_id_list",
    "candidate_gene_same_strand_equal_known_gene_id_list",
    "candidate_gene_same_strand_partial_overlap_known_gene_id_list",
    "candidate_gene_same_strand_overlap_known_gene_id_list",
    "candidate_gene_opposite_strand_complete_overlap_known_gene_count",
    "candidate_gene_opposite_strand_partial_overlap_known_gene_count",
    "candidate_gene_opposite_strand_overlap_known_gene_count",
    "candidate_gene_opposite_strand_complete_overlap_known_gene_id_list",
    "candidate_gene_opposite_strand_partial_overlap_known_gene_id_list",
    "candidate_gene_opposite_strand_overlap_known_gene_id_list",
    "gene_locus_filter_reason",
    "preASE_retained_for_downstream",
]

ASSESSMENT_FIELDS = [
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
    "reference_gene_id",
    "reference_transcript_id",
    "major_iso_id",
    "ref_match_len",
    "tmap_row_number",
    "classification_num_exons",
    "classification_length",
    "FPKM",
    "TPM",
    "cov",
    "FEELnc_query_gene_id",
    "FEELnc_class_code",
    "FEELnc_transcript_length",
    "FEELnc_exon_count",
    "FEELnc_single_exon_flag",
    "FEELnc_filter_status",
    "FEELnc_CPS",
    "FEELnc_RF_label",
    "FEELnc_codpot_category",
    "FEELnc_noORF_flag",
    "FEELnc_biotype",
    "FEELnc_classifier_best_isBest",
    "FEELnc_lncRNA_class",
    "FEELnc_lncRNA_subtype",
    "FEELnc_lncRNA_location",
    "FEELnc_partner_gene",
    "FEELnc_partner_transcript",
    "FEELnc_partner_direction",
    "FEELnc_partner_distance",
    "FEELnc_consistency_flags",
    *TRANSCRIPT_GENE_FILTER_FIELDS,
]

TRANSCRIPT_SUMMARY_FIELDS = ASSESSMENT_FIELDS

GENE_SUMMARY_FIELDS = [
    "query_gene_id",
    "candidate_gene_chr",
    "candidate_gene_start0",
    "candidate_gene_end0",
    "candidate_gene_span_length",
    "candidate_gene_strand",
    "candidate_gene_transcript_strand_list",
    "preliminary_transcript_count",
    "retained_transcript_count",
    "gene_locus_retained_for_downstream",
    "gene_locus_filter_reason",
    "gene_same_strand_contains_known_gene_count",
    "gene_same_strand_within_known_gene_count",
    "gene_same_strand_equal_known_gene_count",
    "gene_same_strand_partial_overlap_known_gene_count",
    "gene_same_strand_overlap_known_gene_count",
    "gene_same_strand_contains_known_gene_id_list",
    "gene_same_strand_within_known_gene_id_list",
    "gene_same_strand_equal_known_gene_id_list",
    "gene_same_strand_partial_overlap_known_gene_id_list",
    "gene_same_strand_overlap_known_gene_id_list",
    "gene_opposite_strand_complete_overlap_known_gene_count",
    "gene_opposite_strand_partial_overlap_known_gene_count",
    "gene_opposite_strand_overlap_known_gene_count",
    "gene_opposite_strand_complete_overlap_known_gene_id_list",
    "gene_opposite_strand_partial_overlap_known_gene_id_list",
    "gene_opposite_strand_overlap_known_gene_id_list",
    "chr_list",
    "transcript_strand_list",
    "class_code_list",
    "novelty_category_list",
    "FEELnc_class_code_list",
    "FEELnc_codpot_category_list",
    "FEELnc_biotype_list",
    "query_transcript_id_list",
]

QC_FIELDS = ["metric", "category", "value", "note"]
MANIFEST_FIELDS = ["item", "path", "exists", "size_bytes", "note"]


@dataclass
class TranscriptModel:
    query_transcript_id: str
    query_gene_id: str
    chrom: str
    strand: str
    span_start0: int
    span_end0: int
    exon_count: int
    exonic_length: int

    @property
    def span_length(self) -> int:
        return max(0, self.span_end0 - self.span_start0)


@dataclass(frozen=True)
class KnownGene:
    gene_id: str
    gene_name: str
    gene_type: str
    chrom: str
    strand: str
    start0: int
    end0: int


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S ",
        stream=sys.stderr,
    )


def clean(value: object) -> str:
    text = "" if value is None else str(value)
    text = text.strip()
    return text if text else "NA"


def merge_values(values: Iterable[object]) -> str:
    merged = sorted({clean(value) for value in values if clean(value) != "NA"})
    return "|".join(merged) if merged else "NA"


def parse_int(value: object, default: int = 0) -> int:
    text = clean(value)
    if text == "NA":
        return default
    try:
        return int(float(text))
    except ValueError:
        return default


def list_field(values: Iterable[str]) -> str:
    merged = sorted({clean(value) for value in values if clean(value) != "NA"})
    return "|".join(merged) if merged else "NA"


def require_readable(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} does not exist or is empty: {path}")


def require_fields(path: Path, fieldnames: Sequence[str] | None, required: Sequence[str]) -> None:
    available = list(fieldnames or [])
    missing = [field for field in required if field not in available]
    if missing:
        raise RuntimeError(f"{path} missing required fields: {', '.join(missing)}")


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(fieldnames),
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean(row.get(field)) for field in fieldnames})


def parse_gtf_attributes(attribute_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attribute_text.rstrip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            attrs[item] = ""
            continue
        key, value = item.split(" ", 1)
        value = value.strip()
        if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
            value = value[1:-1]
        attrs[key] = value
    return attrs


def read_classification(path: Path) -> List[Dict[str, str]]:
    logging.info("read five-dogs Step 1.6 gffcompare transcript classification: %s", path)
    rows: List[Dict[str, str]] = []
    seen: Set[str] = set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, CLASSIFICATION_REQUIRED_FIELDS)
        for row_number, row in enumerate(reader, start=1):
            transcript_id = clean(row.get("query_transcript_id"))
            if transcript_id == "NA":
                raise RuntimeError(f"{path} row {row_number} is missing query_transcript_id")
            if transcript_id in seen:
                raise RuntimeError(f"{path} query_transcript_id duplicate: {transcript_id}")
            seen.add(transcript_id)
            rows.append({field: clean(value) for field, value in row.items()})
            if row_number % 100000 == 0:
                logging.info("read transcript classification rows: %s", row_number)
    logging.info("Step 1.6 transcript classification read completed: %s transcripts", len(rows))
    return rows


def read_feelnc_status(path: Path) -> Dict[str, Dict[str, str]]:
    logging.info("read five-dogs Step 1.6.5c FEELnc module status: %s", path)
    rows: Dict[str, Dict[str, str]] = {}
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(path, reader.fieldnames, FEELNC_REQUIRED_FIELDS)
        for row_number, row in enumerate(reader, start=1):
            transcript_id = clean(row.get("query_transcript_id"))
            if transcript_id == "NA":
                raise RuntimeError(f"{path} row {row_number} is missing query_transcript_id")
            if transcript_id in rows:
                raise RuntimeError(f"{path} query_transcript_id duplicate: {transcript_id}")
            rows[transcript_id] = {field: clean(value) for field, value in row.items()}
            if row_number % 100000 == 0:
                logging.info("read FEELnc status rows: %s", row_number)
    logging.info("Step 1.6.5c FEELnc status read completed: %s transcripts", len(rows))
    return rows


def read_transcript_models(path: Path) -> Dict[str, TranscriptModel]:
    logging.info("read five-dogs Step 1.6 transcript exon BED: %s", path)
    models: MutableMapping[str, TranscriptModel] = {}
    with path.open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise RuntimeError(f"{path} row {row_number} has fewer than 8 BED columns")
            chrom, start_text, end_text, transcript_id, gene_id, strand = fields[:6]
            start0 = parse_int(start_text, default=-1)
            end0 = parse_int(end_text, default=-1)
            if start0 < 0 or end0 <= start0:
                raise RuntimeError(f"{path} row {row_number} has invalid coordinates: {line}")
            exon_length = end0 - start0
            transcript_id = clean(transcript_id)
            gene_id = clean(gene_id)
            chrom = clean(chrom)
            strand = clean(strand)
            if "NA" in {transcript_id, gene_id, chrom, strand}:
                raise RuntimeError(f"{path} row {row_number} is missing transcript/gene/chrom/strand: {line}")
            if transcript_id not in models:
                models[transcript_id] = TranscriptModel(
                    query_transcript_id=transcript_id,
                    query_gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    span_start0=start0,
                    span_end0=end0,
                    exon_count=1,
                    exonic_length=exon_length,
                )
            else:
                model = models[transcript_id]
                if (model.query_gene_id, model.chrom, model.strand) != (gene_id, chrom, strand):
                    raise RuntimeError(
                        f"{path} transcript of exon metadata mismatch: transcript={transcript_id} "
                        f"expected=({model.query_gene_id},{model.chrom},{model.strand}) "
                        f"observed=({gene_id},{chrom},{strand})"
                    )
                model.span_start0 = min(model.span_start0, start0)
                model.span_end0 = max(model.span_end0, end0)
                model.exon_count += 1
                model.exonic_length += exon_length
            if row_number % 1000000 == 0:
                logging.info("read transcript exon BED rows: %s", row_number)
    logging.info("Step 1.6 transcript exon BED read completed: %s transcripts", len(models))
    return dict(models)


def read_reference_genes_gtf(path: Path, allowed_chroms: Set[str]) -> List[KnownGene]:
    logging.info("Read the CanFam3.1 reference annotation and aggregate gene spans by gene_id: %s", path)


    grouped: Dict[Tuple[str, str], Dict[str, object]] = {}
    transcript_to_gene_id: Dict[str, str] = {}
    reference_rows_seen = 0
    exon_gene_id_fallback_rows = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line_number % 1000000 == 0:
                logging.info("scanned CanFam3.1 GTF rows: %s, aggregated_genes=%s", line_number, len(grouped))
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise RuntimeError(f"{path} row {line_number} does not have 9 GTF columns")
            chrom, _source, feature, start_text, end_text, _score, strand, _frame, attr_text = fields
            if feature not in {"transcript", "exon"}:
                continue
            reference_rows_seen += 1
            if allowed_chroms and chrom not in allowed_chroms:
                continue
            start1 = parse_int(start_text, default=-1)
            end1 = parse_int(end_text, default=-1)
            if start1 <= 0 or end1 < start1:
                raise RuntimeError(f"{path} row {line_number} has invalid reference coordinates: {line.rstrip()}")
            attrs = parse_gtf_attributes(attr_text)
            gene_id = clean(attrs.get("gene_id"))
            transcript_id = clean(attrs.get("transcript_id"))
            if feature == "transcript" and transcript_id != "NA" and gene_id != "NA":
                transcript_to_gene_id[transcript_id] = gene_id
            if gene_id == "NA" and feature == "exon" and transcript_id != "NA":
                gene_id = transcript_to_gene_id.get(transcript_id, "NA")
                if gene_id != "NA":
                    exon_gene_id_fallback_rows += 1
            if gene_id == "NA":
                raise RuntimeError(
                    f"{path} row {line_number} reference feature is missing gene_id, "
                    f"and transcript_id cannot be mapped: transcript_id={transcript_id}"
                )
            key = (clean(chrom), gene_id)
            current = grouped.get(key)
            if current is None:
                grouped[key] = {
                    "gene_id": gene_id,
                    "gene_name": clean(attrs.get("gene_name")),
                    "gene_type": clean(attrs.get("gene_type", attrs.get("gene_biotype"))),
                    "chrom": clean(chrom),
                    "strands": {clean(strand)},
                    "start0": start1 - 1,
                    "end0": end1,
                }
            else:
                current["start0"] = min(int(current["start0"]), start1 - 1)
                current["end0"] = max(int(current["end0"]), end1)
                cast_strands = current["strands"]
                if isinstance(cast_strands, set):
                    cast_strands.add(clean(strand))
                if clean(current.get("gene_name")) == "NA" and clean(attrs.get("gene_name")) != "NA":
                    current["gene_name"] = clean(attrs.get("gene_name"))

    genes: List[KnownGene] = []
    for current in grouped.values():
        strands = current["strands"] if isinstance(current["strands"], set) else set()
        genes.append(
            KnownGene(
                gene_id=clean(current["gene_id"]),
                gene_name=clean(current["gene_name"]),
                gene_type=clean(current["gene_type"]),
                chrom=clean(current["chrom"]),
                strand=merge_values(strands),
                start0=int(current["start0"]),
                end0=int(current["end0"]),
            )
        )
    genes.sort(key=lambda gene: (gene.chrom, gene.start0, gene.end0, gene.gene_id))
    logging.info(
        "CanFam3.1 gene span aggregation completed: source_rows=%s exon_gene_id_fallback_rows=%s genes=%s",
        reference_rows_seen,
        exon_gene_id_fallback_rows,
        len(genes),
    )
    return genes


def validate_linked_inputs(
    classification_rows: Sequence[Mapping[str, str]],
    feelnc_by_transcript: Mapping[str, Mapping[str, str]],
    models_by_transcript: Mapping[str, TranscriptModel],
) -> None:

    classification_by_transcript = {
        clean(row.get("query_transcript_id")): row for row in classification_rows
    }
    expected_ids = set(classification_by_transcript)
    for label, observed_ids in [
        ("FEELnc status", set(feelnc_by_transcript)),
        ("transcript exon BED", set(models_by_transcript)),
    ]:
        missing = sorted(expected_ids - observed_ids)
        extra = sorted(observed_ids - expected_ids)
        if missing or extra:
            raise RuntimeError(
                f"query_transcript_id sets differ between {label} and the classification table: "
                f"missing={len(missing)} extra={len(extra)} "
                f"missing_examples={missing[:5]} extra_examples={extra[:5]}"
            )

    mismatches: List[str] = []
    for transcript_id, source in classification_by_transcript.items():
        source_gene_id = clean(source.get("query_gene_id"))
        source_class_code = clean(source.get("class_code"))
        feelnc = feelnc_by_transcript[transcript_id]
        model = models_by_transcript[transcript_id]
        checks = [
            ("classification query_gene_id", source_gene_id != "NA"),
            ("classification/FEELnc query_gene_id", source_gene_id == clean(feelnc.get("query_gene_id"))),
            ("classification/BED query_gene_id", source_gene_id == model.query_gene_id),
            ("classification/FEELnc class_code", source_class_code == clean(feelnc.get("class_code"))),
            ("classification/BED exon_count", parse_int(source.get("num_exons"), -1) == model.exon_count),
            ("classification/BED transcript_length", parse_int(source.get("length"), -1) == model.exonic_length),
            ("FEELnc/BED exon_count", parse_int(feelnc.get("exon_count"), -1) == model.exon_count),
            ("FEELnc/BED transcript_length", parse_int(feelnc.get("transcript_length"), -1) == model.exonic_length),
        ]
        failed_labels = [label for label, passed in checks if not passed]
        if failed_labels:
            mismatches.append(f"{transcript_id}:{','.join(failed_labels)}")
            if len(mismatches) == 5:
                break
    if mismatches:
        raise RuntimeError("transcript input join field mismatch, example: " + "; ".join(mismatches))
    logging.info("three sets transcript input main key, gene, class, exon count and length checks passed")


def build_reference_bins(genes: Sequence[KnownGene]) -> Dict[str, Dict[int, List[int]]]:
    bins: Dict[str, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))
    for index, gene in enumerate(genes):
        start_bin = gene.start0 // REFERENCE_BIN_SIZE
        end_bin = max(gene.start0, gene.end0 - 1) // REFERENCE_BIN_SIZE
        for bin_id in range(start_bin, end_bin + 1):
            bins[gene.chrom][bin_id].append(index)
    return bins


def overlapping_known_gene_indices(
    chrom: str,
    start0: int,
    end0: int,
    reference_genes: Sequence[KnownGene],
    reference_bins: Mapping[str, Mapping[int, Sequence[int]]],
) -> List[int]:
    seq_bins = reference_bins.get(chrom)
    if not seq_bins:
        return []
    start_bin = start0 // REFERENCE_BIN_SIZE
    end_bin = max(start0, end0 - 1) // REFERENCE_BIN_SIZE
    candidates: Set[int] = set()
    for bin_id in range(start_bin, end_bin + 1):
        candidates.update(seq_bins.get(bin_id, []))
    overlaps = [
        index
        for index in candidates
        if reference_genes[index].start0 < end0 and reference_genes[index].end0 > start0
    ]
    overlaps.sort(key=lambda index: (reference_genes[index].start0, reference_genes[index].end0, reference_genes[index].gene_id))
    return overlaps


def model_fields(model: TranscriptModel | None) -> Dict[str, str]:
    if model is None:
        return {
            "chr": "NA",
            "transcript_strand": "NA",
            "transcript_span_start0": "NA",
            "transcript_span_end0": "NA",
            "transcript_span_length": "NA",
            "transcript_exonic_length": "NA",
            "transcript_exon_count": "NA",
        }
    return {
        "chr": model.chrom,
        "transcript_strand": model.strand,
        "transcript_span_start0": str(model.span_start0),
        "transcript_span_end0": str(model.span_end0),
        "transcript_span_length": str(model.span_length),
        "transcript_exonic_length": str(model.exonic_length),
        "transcript_exon_count": str(model.exon_count),
    }


def feelnc_fields(feelnc: Mapping[str, str]) -> Dict[str, str]:
    return {
        "FEELnc_query_gene_id": clean(feelnc.get("query_gene_id")),
        "FEELnc_class_code": clean(feelnc.get("class_code")),
        "FEELnc_transcript_length": clean(feelnc.get("transcript_length")),
        "FEELnc_exon_count": clean(feelnc.get("exon_count")),
        "FEELnc_single_exon_flag": clean(feelnc.get("single_exon_flag")),
        "FEELnc_filter_status": clean(feelnc.get("FEELnc_filter_status")),
        "FEELnc_CPS": clean(feelnc.get("FEELnc_CPS")),
        "FEELnc_RF_label": clean(feelnc.get("FEELnc_RF_label")),
        "FEELnc_codpot_category": clean(feelnc.get("FEELnc_codpot_category")),
        "FEELnc_noORF_flag": clean(feelnc.get("FEELnc_noORF_flag")),
        "FEELnc_biotype": clean(feelnc.get("FEELnc_biotype")),
        "FEELnc_classifier_best_isBest": clean(feelnc.get("FEELnc_classifier_best_isBest")),
        "FEELnc_lncRNA_class": clean(feelnc.get("FEELnc_lncRNA_class")),
        "FEELnc_lncRNA_subtype": clean(feelnc.get("FEELnc_lncRNA_subtype")),
        "FEELnc_lncRNA_location": clean(feelnc.get("FEELnc_lncRNA_location")),
        "FEELnc_partner_gene": clean(feelnc.get("FEELnc_partner_gene")),
        "FEELnc_partner_transcript": clean(feelnc.get("FEELnc_partner_transcript")),
        "FEELnc_partner_direction": clean(feelnc.get("FEELnc_partner_direction")),
        "FEELnc_partner_distance": clean(feelnc.get("FEELnc_partner_distance")),
        "FEELnc_consistency_flags": clean(feelnc.get("FEELnc_consistency_flags")),
    }


def empty_gene_filter_fields(reason: str = "exact_known_transcript_not_used_for_gene_locus") -> Dict[str, str]:
    return {
        "candidate_gene_chr": "NA",
        "candidate_gene_start0": "NA",
        "candidate_gene_end0": "NA",
        "candidate_gene_span_length": "NA",
        "candidate_gene_strand": "NA",
        "candidate_gene_transcript_count": "0",
        "candidate_gene_same_strand_contains_known_gene_count": "0",
        "candidate_gene_same_strand_within_known_gene_count": "0",
        "candidate_gene_same_strand_equal_known_gene_count": "0",
        "candidate_gene_same_strand_partial_overlap_known_gene_count": "0",
        "candidate_gene_same_strand_overlap_known_gene_count": "0",
        "candidate_gene_same_strand_contains_known_gene_id_list": "NA",
        "candidate_gene_same_strand_within_known_gene_id_list": "NA",
        "candidate_gene_same_strand_equal_known_gene_id_list": "NA",
        "candidate_gene_same_strand_partial_overlap_known_gene_id_list": "NA",
        "candidate_gene_same_strand_overlap_known_gene_id_list": "NA",
        "candidate_gene_opposite_strand_complete_overlap_known_gene_count": "0",
        "candidate_gene_opposite_strand_partial_overlap_known_gene_count": "0",
        "candidate_gene_opposite_strand_overlap_known_gene_count": "0",
        "candidate_gene_opposite_strand_complete_overlap_known_gene_id_list": "NA",
        "candidate_gene_opposite_strand_partial_overlap_known_gene_id_list": "NA",
        "candidate_gene_opposite_strand_overlap_known_gene_id_list": "NA",
        "gene_locus_filter_reason": reason,
    }


def gene_filter_fields(gene_row: Mapping[str, str]) -> Dict[str, str]:
    return {
        "candidate_gene_chr": clean(gene_row.get("candidate_gene_chr")),
        "candidate_gene_start0": clean(gene_row.get("candidate_gene_start0")),
        "candidate_gene_end0": clean(gene_row.get("candidate_gene_end0")),
        "candidate_gene_span_length": clean(gene_row.get("candidate_gene_span_length")),
        "candidate_gene_strand": clean(gene_row.get("candidate_gene_strand")),
        "candidate_gene_transcript_count": clean(gene_row.get("preliminary_transcript_count")),
        "candidate_gene_same_strand_contains_known_gene_count": clean(gene_row.get("gene_same_strand_contains_known_gene_count")),
        "candidate_gene_same_strand_within_known_gene_count": clean(gene_row.get("gene_same_strand_within_known_gene_count")),
        "candidate_gene_same_strand_equal_known_gene_count": clean(gene_row.get("gene_same_strand_equal_known_gene_count")),
        "candidate_gene_same_strand_partial_overlap_known_gene_count": clean(gene_row.get("gene_same_strand_partial_overlap_known_gene_count")),
        "candidate_gene_same_strand_overlap_known_gene_count": clean(gene_row.get("gene_same_strand_overlap_known_gene_count")),
        "candidate_gene_same_strand_contains_known_gene_id_list": clean(gene_row.get("gene_same_strand_contains_known_gene_id_list")),
        "candidate_gene_same_strand_within_known_gene_id_list": clean(gene_row.get("gene_same_strand_within_known_gene_id_list")),
        "candidate_gene_same_strand_equal_known_gene_id_list": clean(gene_row.get("gene_same_strand_equal_known_gene_id_list")),
        "candidate_gene_same_strand_partial_overlap_known_gene_id_list": clean(gene_row.get("gene_same_strand_partial_overlap_known_gene_id_list")),
        "candidate_gene_same_strand_overlap_known_gene_id_list": clean(gene_row.get("gene_same_strand_overlap_known_gene_id_list")),
        "candidate_gene_opposite_strand_complete_overlap_known_gene_count": clean(gene_row.get("gene_opposite_strand_complete_overlap_known_gene_count")),
        "candidate_gene_opposite_strand_partial_overlap_known_gene_count": clean(gene_row.get("gene_opposite_strand_partial_overlap_known_gene_count")),
        "candidate_gene_opposite_strand_overlap_known_gene_count": clean(gene_row.get("gene_opposite_strand_overlap_known_gene_count")),
        "candidate_gene_opposite_strand_complete_overlap_known_gene_id_list": clean(gene_row.get("gene_opposite_strand_complete_overlap_known_gene_id_list")),
        "candidate_gene_opposite_strand_partial_overlap_known_gene_id_list": clean(gene_row.get("gene_opposite_strand_partial_overlap_known_gene_id_list")),
        "candidate_gene_opposite_strand_overlap_known_gene_id_list": clean(gene_row.get("gene_opposite_strand_overlap_known_gene_id_list")),
        "gene_locus_filter_reason": clean(gene_row.get("gene_locus_filter_reason")),
    }


def build_transcript_assessment_rows(
    classification_rows: Sequence[Mapping[str, str]],
    feelnc_by_transcript: Mapping[str, Mapping[str, str]],
    models_by_transcript: Mapping[str, TranscriptModel],
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], Counter]:
    logging.info("Stage 1/5: build transcript-level assessment, and flag class_code != '=' initial screening result")
    metrics: Counter = Counter()
    all_rows: List[Dict[str, str]] = []
    preliminary_rows: List[Dict[str, str]] = []
    for index, source in enumerate(classification_rows, start=1):
        transcript_id = clean(source.get("query_transcript_id"))
        model = models_by_transcript.get(transcript_id)
        feelnc = feelnc_by_transcript.get(transcript_id, {})
        retained_by_class = "0" if clean(source.get("class_code")) == "=" else "1"
        row: Dict[str, str] = {
            "query_transcript_id": transcript_id,
            "query_gene_id": clean(source.get("query_gene_id")),
            "class_code": clean(source.get("class_code")),
            "novelty_category": clean(source.get("novelty_category")),
            "reference_gene_id": clean(source.get("reference_gene_id")),
            "reference_transcript_id": clean(source.get("reference_transcript_id")),
            "major_iso_id": clean(source.get("major_iso_id")),
            "ref_match_len": clean(source.get("ref_match_len")),
            "tmap_row_number": clean(source.get("tmap_row_number")),
            "classification_num_exons": clean(source.get("num_exons")),
            "classification_length": clean(source.get("length")),
            "FPKM": clean(source.get("FPKM")),
            "TPM": clean(source.get("TPM")),
            "cov": clean(source.get("cov")),
            "preASE_retained_by_class_code": retained_by_class,
            "preASE_retained_for_downstream": "0",
        }
        row.update(model_fields(model))
        row.update(feelnc_fields(feelnc))
        row.update(empty_gene_filter_fields())
        all_rows.append(row)
        metrics["all_transcript_rows"] += 1
        metrics[f"all_class_code_{row['class_code']}"] += 1
        if not feelnc:
            metrics["missing_feelnc_status_rows"] += 1
        elif clean(feelnc.get("class_code")) != clean(source.get("class_code")):
            metrics["feelnc_class_code_mismatch_rows"] += 1
        if model is None:
            metrics["missing_exon_model_rows"] += 1
        if retained_by_class == "1":
            preliminary_rows.append(row)
            metrics["preliminary_retained_transcript_rows"] += 1
            metrics[f"preliminary_retained_class_code_{row['class_code']}"] += 1
        else:
            metrics["exact_known_transcript_excluded_rows"] += 1
        if index % 100000 == 0:
            logging.info("completed transcript assessment rows: %s", index)
    metrics["preliminary_unique_transcript_count"] = len({row["query_transcript_id"] for row in preliminary_rows})
    metrics["preliminary_unique_gene_locus_count"] = len({row["query_gene_id"] for row in preliminary_rows})
    return all_rows, preliminary_rows, metrics


def classify_gene_filter_reason(
    retained: bool,
    contains_count: int,
    within_count: int,
    equal_count: int,
    same_strand_partial_count: int,
    opposite_strand_complete_count: int,
    opposite_strand_partial_count: int,
) -> str:
    if retained:
        if opposite_strand_complete_count > 0:
            return "retained_opposite_strand_complete_overlap"
        if same_strand_partial_count > 0:
            return "retained_same_strand_partial_overlap"
        if opposite_strand_partial_count > 0:
            return "retained_opposite_strand_partial_overlap"
        return "retained_no_known_gene_overlap"
    if equal_count > 0:
        return "excluded_same_strand_equal_known_gene_interval"
    if contains_count > 0 and within_count > 0:
        return "excluded_same_strand_new_contains_and_within_known_gene"
    if contains_count > 0:
        return "excluded_same_strand_new_contains_known_gene"
    if within_count > 0:
        return "excluded_same_strand_new_within_known_gene"
    return "excluded_invalid_gene_locus"


def build_gene_locus_assessment(
    preliminary_rows: Sequence[Mapping[str, str]],
    reference_genes: Sequence[KnownGene],
    reference_bins: Mapping[str, Mapping[int, Sequence[int]]],
    containment_overlap_fraction: float,
) -> Tuple[List[Dict[str, str]], List[Dict[str, str]], Counter]:
    logging.info(
        "Stage 2/5: aggregate gene loci by query_gene_id; exclude near-complete containment only for same-strand overlap fraction >= %.3f",
        containment_overlap_fraction,
    )
    metrics: Counter = Counter()
    grouped: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in preliminary_rows:
        grouped[clean(row.get("query_gene_id"))].append(row)

    all_gene_rows: List[Dict[str, str]] = []
    retained_gene_rows: List[Dict[str, str]] = []
    for gene_id in sorted(grouped):
        rows = grouped[gene_id]
        transcript_ids = {clean(row.get("query_transcript_id")) for row in rows}
        chroms = {clean(row.get("chr")) for row in rows if clean(row.get("chr")) != "NA"}

        strands = {clean(row.get("transcript_strand")) for row in rows if clean(row.get("transcript_strand")) != "NA"}
        candidate_strand = next(iter(strands)) if len(strands) == 1 else "NA"
        if candidate_strand not in {"+", "-"}:
            candidate_strand = "NA"
        starts = [parse_int(row.get("transcript_span_start0"), default=-1) for row in rows]
        ends = [parse_int(row.get("transcript_span_end0"), default=-1) for row in rows]
        valid_coords = [coord for coord in zip(starts, ends) if coord[0] >= 0 and coord[1] > coord[0]]

        metrics["preliminary_gene_locus_count"] += 1
        if len(chroms) == 1 and valid_coords and candidate_strand != "NA":
            chrom = next(iter(chroms))
            start0 = min(start for start, _end in valid_coords)
            end0 = max(end for _start, end in valid_coords)
            overlap_indices = overlapping_known_gene_indices(chrom, start0, end0, reference_genes, reference_bins)
            candidate_gene_length = end0 - start0
            same_strand_contains_ids: Set[str] = set()
            same_strand_within_ids: Set[str] = set()
            same_strand_equal_ids: Set[str] = set()
            same_strand_partial_ids: Set[str] = set()
            same_strand_overlap_ids: Set[str] = set()
            opposite_strand_complete_ids: Set[str] = set()
            opposite_strand_partial_ids: Set[str] = set()
            opposite_strand_overlap_ids: Set[str] = set()
            for known_index in overlap_indices:
                known = reference_genes[known_index]
                overlap_bp = max(0, min(end0, known.end0) - max(start0, known.start0))
                known_gene_length = known.end0 - known.start0
                known_coverage = overlap_bp / known_gene_length if known_gene_length > 0 else 0.0
                candidate_coverage = overlap_bp / candidate_gene_length if candidate_gene_length > 0 else 0.0
                new_contains_known = known_coverage >= containment_overlap_fraction
                new_within_known = candidate_coverage >= containment_overlap_fraction
                equal_interval = start0 == known.start0 and end0 == known.end0

                if known.strand == candidate_strand:
                    same_strand_overlap_ids.add(known.gene_id)
                    if new_contains_known:
                        same_strand_contains_ids.add(known.gene_id)
                    if new_within_known:
                        same_strand_within_ids.add(known.gene_id)
                    if equal_interval:
                        same_strand_equal_ids.add(known.gene_id)
                    if not new_contains_known and not new_within_known:
                        same_strand_partial_ids.add(known.gene_id)
                elif known.strand in {"+", "-"}:
                    opposite_strand_overlap_ids.add(known.gene_id)
                    if new_contains_known or new_within_known:
                        opposite_strand_complete_ids.add(known.gene_id)
                    else:
                        opposite_strand_partial_ids.add(known.gene_id)

            contains_count = len(same_strand_contains_ids)
            within_count = len(same_strand_within_ids)
            equal_count = len(same_strand_equal_ids)
            same_strand_partial_count = len(same_strand_partial_ids)
            same_strand_overlap_count = len(same_strand_overlap_ids)
            opposite_strand_complete_count = len(opposite_strand_complete_ids)
            opposite_strand_partial_count = len(opposite_strand_partial_ids)
            opposite_strand_overlap_count = len(opposite_strand_overlap_ids)
            retained = contains_count == 0 and within_count == 0
            reason = classify_gene_filter_reason(
                retained,
                contains_count,
                within_count,
                equal_count,
                same_strand_partial_count,
                opposite_strand_complete_count,
                opposite_strand_partial_count,
            )
        else:
            chrom = merge_values(chroms)
            start0 = -1
            end0 = -1
            same_strand_contains_ids = set()
            same_strand_within_ids = set()
            same_strand_equal_ids = set()
            same_strand_partial_ids = set()
            same_strand_overlap_ids = set()
            opposite_strand_complete_ids = set()
            opposite_strand_partial_ids = set()
            opposite_strand_overlap_ids = set()
            contains_count = 0
            within_count = 0
            equal_count = 0
            same_strand_partial_count = 0
            same_strand_overlap_count = 0
            opposite_strand_complete_count = 0
            opposite_strand_partial_count = 0
            opposite_strand_overlap_count = 0
            retained = False
            if candidate_strand == "NA":
                reason = "excluded_mixed_or_noncanonical_strand"
                metrics["gene_locus_mixed_or_noncanonical_strand_count"] += 1
            else:
                reason = "excluded_mixed_or_missing_gene_coordinates"
                metrics["gene_locus_mixed_or_missing_coordinate_count"] += 1

        row = {
            "query_gene_id": gene_id,
            "candidate_gene_chr": chrom,
            "candidate_gene_start0": str(start0) if start0 >= 0 else "NA",
            "candidate_gene_end0": str(end0) if end0 > start0 else "NA",
            "candidate_gene_span_length": str(max(0, end0 - start0)) if start0 >= 0 and end0 > start0 else "NA",
            "candidate_gene_strand": candidate_strand,
            "candidate_gene_transcript_strand_list": merge_values(row.get("transcript_strand") for row in rows),
            "preliminary_transcript_count": str(len(transcript_ids)),
            "retained_transcript_count": str(len(transcript_ids) if retained else 0),
            "gene_locus_retained_for_downstream": "1" if retained else "0",
            "gene_locus_filter_reason": reason,
            "gene_same_strand_contains_known_gene_count": str(contains_count),
            "gene_same_strand_within_known_gene_count": str(within_count),
            "gene_same_strand_equal_known_gene_count": str(equal_count),
            "gene_same_strand_partial_overlap_known_gene_count": str(same_strand_partial_count),
            "gene_same_strand_overlap_known_gene_count": str(same_strand_overlap_count),
            "gene_same_strand_contains_known_gene_id_list": list_field(same_strand_contains_ids),
            "gene_same_strand_within_known_gene_id_list": list_field(same_strand_within_ids),
            "gene_same_strand_equal_known_gene_id_list": list_field(same_strand_equal_ids),
            "gene_same_strand_partial_overlap_known_gene_id_list": list_field(same_strand_partial_ids),
            "gene_same_strand_overlap_known_gene_id_list": list_field(same_strand_overlap_ids),
            "gene_opposite_strand_complete_overlap_known_gene_count": str(opposite_strand_complete_count),
            "gene_opposite_strand_partial_overlap_known_gene_count": str(opposite_strand_partial_count),
            "gene_opposite_strand_overlap_known_gene_count": str(opposite_strand_overlap_count),
            "gene_opposite_strand_complete_overlap_known_gene_id_list": list_field(opposite_strand_complete_ids),
            "gene_opposite_strand_partial_overlap_known_gene_id_list": list_field(opposite_strand_partial_ids),
            "gene_opposite_strand_overlap_known_gene_id_list": list_field(opposite_strand_overlap_ids),
            "chr_list": merge_values(row.get("chr") for row in rows),
            "transcript_strand_list": merge_values(row.get("transcript_strand") for row in rows),
            "class_code_list": merge_values(row.get("class_code") for row in rows),
            "novelty_category_list": merge_values(row.get("novelty_category") for row in rows),
            "FEELnc_class_code_list": merge_values(row.get("FEELnc_class_code") for row in rows),
            "FEELnc_codpot_category_list": merge_values(row.get("FEELnc_codpot_category") for row in rows),
            "FEELnc_biotype_list": merge_values(row.get("FEELnc_biotype") for row in rows),
            "query_transcript_id_list": merge_values(transcript_ids),
        }
        all_gene_rows.append(row)
        metrics[f"gene_locus_filter_reason_{reason}"] += 1
        if contains_count > 0:
            metrics["gene_locus_same_strand_new_contains_known_gene_count"] += 1
        if within_count > 0:
            metrics["gene_locus_same_strand_new_within_known_gene_count"] += 1
        if equal_count > 0:
            metrics["gene_locus_same_strand_equal_known_gene_count"] += 1
        if same_strand_partial_count > 0 and contains_count == 0 and within_count == 0:
            metrics["gene_locus_same_strand_partial_overlap_only_count"] += 1
        if opposite_strand_complete_count > 0:
            metrics["gene_locus_opposite_strand_complete_overlap_count"] += 1
        if opposite_strand_partial_count > 0:
            metrics["gene_locus_opposite_strand_partial_overlap_count"] += 1
        if retained:
            retained_gene_rows.append(row)
            metrics["retained_gene_locus_count"] += 1

    logging.info(
        "gene-level containment filter completed: preliminary_gene_loci=%s retained_gene_loci=%s",
        len(all_gene_rows),
        len(retained_gene_rows),
    )
    return all_gene_rows, retained_gene_rows, metrics


def apply_gene_filter_to_transcripts(
    all_rows: Sequence[MutableMapping[str, str]],
    gene_rows_by_id: Mapping[str, Mapping[str, str]],
    retained_gene_ids: Set[str],
) -> Tuple[List[Dict[str, str]], Counter]:
    logging.info("Stage 3/5: assign transcript-level downstream flags using the retained gene loci")
    metrics: Counter = Counter()
    retained_rows: List[Dict[str, str]] = []
    for row in all_rows:
        gene_id = clean(row.get("query_gene_id"))
        retained_by_class = clean(row.get("preASE_retained_by_class_code")) == "1"
        gene_row = gene_rows_by_id.get(gene_id)
        if gene_row:
            row.update(gene_filter_fields(gene_row))
        final_retained = retained_by_class and gene_id in retained_gene_ids
        row["preASE_retained_for_downstream"] = "1" if final_retained else "0"
        if final_retained:
            retained_rows.append({field: clean(row.get(field)) for field in ASSESSMENT_FIELDS})
            metrics["retained_transcript_rows"] += 1
            metrics[f"retained_class_code_{clean(row.get('class_code'))}"] += 1
        elif retained_by_class:
            metrics["transcripts_excluded_by_gene_locus_filter"] += 1
    metrics["retained_unique_transcript_count"] = len({row["query_transcript_id"] for row in retained_rows})
    metrics["retained_unique_gene_locus_count"] = len({row["query_gene_id"] for row in retained_rows})
    return retained_rows, metrics


def build_transcript_summary(retained_rows: Sequence[Mapping[str, str]]) -> List[Dict[str, str]]:
    logging.info("Stage 4/5: aggregate the retained-transcript summary")
    return [{field: clean(row.get(field)) for field in TRANSCRIPT_SUMMARY_FIELDS} for row in retained_rows]


def qc_rows(
    args: argparse.Namespace,
    classification_rows: Sequence[Mapping[str, str]],
    feelnc_by_transcript: Mapping[str, Mapping[str, str]],
    models_by_transcript: Mapping[str, TranscriptModel],
    reference_genes: Sequence[KnownGene],
    all_rows: Sequence[Mapping[str, str]],
    preliminary_rows: Sequence[Mapping[str, str]],
    all_gene_rows: Sequence[Mapping[str, str]],
    retained_gene_rows: Sequence[Mapping[str, str]],
    retained_rows: Sequence[Mapping[str, str]],
    transcript_rows: Sequence[Mapping[str, str]],
    metrics: Mapping[str, int],
) -> List[Dict[str, str]]:
    expected_ok = args.expected_transcript_count <= 0 or len(classification_rows) == args.expected_transcript_count
    pass_status = (
        "PASS"
        if expected_ok
        and len(reference_genes) > 0
        and len(all_rows) == len(classification_rows)
        and len(all_rows) == metrics.get("all_transcript_rows", -1)
        and len(preliminary_rows) == metrics.get("preliminary_retained_transcript_rows", -1)
        and len(all_gene_rows) == metrics.get("preliminary_gene_locus_count", -1)
        and len(retained_gene_rows) == metrics.get("retained_gene_locus_count", -1)
        and len(retained_rows) == metrics.get("retained_transcript_rows", -1)
        and len(transcript_rows) == metrics.get("retained_unique_transcript_count", -1)
        and len(retained_gene_rows) == metrics.get("retained_unique_gene_locus_count", -1)
        and metrics.get("missing_feelnc_status_rows", 0) == 0
        and metrics.get("missing_exon_model_rows", 0) == 0
        else "FAIL"
    )
    rows = [
        {
            "metric": "expected_transcript_count",
            "category": "all",
            "value": str(args.expected_transcript_count),
            "note": "Derived by the batch launcher from the Step 1.6 classification row count; 0 disables the check.",
        },
        {
            "metric": "classification_transcript_count",
            "category": "all",
            "value": str(len(classification_rows)),
            "note": str(args.classification_tsv),
        },
        {
            "metric": "feelnc_status_transcript_count",
            "category": "all",
            "value": str(len(feelnc_by_transcript)),
            "note": str(args.feelnc_status_tsv),
        },
        {
            "metric": "exon_model_transcript_count",
            "category": "all",
            "value": str(len(models_by_transcript)),
            "note": str(args.transcript_exon_bed),
        },
        {
            "metric": "transcript_input_link_validation",
            "category": "all",
            "value": "PASS",
            "note": "Transcript IDs, gene IDs, classes, exon counts, and lengths agree across classification, FEELnc, and BED inputs.",
        },
        {
            "metric": "reference_gtf_aggregated_gene_count",
            "category": "all",
            "value": str(len(reference_genes)),
            "note": f"{args.reference_gtf}; transcript/exon records aggregated by gene_id to define known-gene spans.",
        },
        {
            "metric": "containment_overlap_fraction",
            "category": "parameter",
            "value": f"{args.containment_overlap_fraction:.6g}",
            "note": "For same-strand overlaps only, remove the gene locus if known_coverage=overlap/known_gene_length or candidate_coverage=overlap/candidate_gene_length reaches this threshold.",
        },
        {
            "metric": "containment_strand_rule",
            "category": "parameter",
            "value": "same_strand_only",
            "note": "Retain and separately label complete/partial antisense overlaps; exclude candidate loci with mixed strands or strands other than + / -.",
        },
        {
            "metric": "all_transcript_rows",
            "category": "all",
            "value": str(len(all_rows)),
            "note": "One row per Step 1.6 query_transcript_id; retain GffCompare/FEELnc annotation fields.",
        },
        {
            "metric": "preliminary_retained_transcript_rows",
            "category": "all",
            "value": str(len(preliminary_rows)),
            "note": "Initial transcript filtering excludes only gffcompare class_code == '='.",
        },
        {
            "metric": "all_gene_locus_rows",
            "category": "all",
            "value": str(len(all_gene_rows)),
            "note": "One row per query_gene_id/StringTie locus aggregated from transcripts with class_code != '='.",
        },
        {
            "metric": "retained_gene_locus_rows",
            "category": "all",
            "value": str(len(retained_gene_rows)),
            "note": "Remove only same-strand new-contains-known, new-within-known, and equal-interval genes whose overlap coverage reaches the threshold; retain and annotate antisense overlaps.",
        },
        {
            "metric": "retained_transcript_rows",
            "category": "all",
            "value": str(len(retained_rows)),
            "note": "Transcripts belonging to retained genes with class_code != '='.",
        },
        {
            "metric": "candidate_transcript_summary_rows",
            "category": "all",
            "value": str(len(transcript_rows)),
            "note": "One row per final retained query_transcript_id.",
        },
        {
            "metric": "candidate_gene_locus_summary_rows",
            "category": "all",
            "value": str(len(retained_gene_rows)),
            "note": "One row per final retained query_gene_id/StringTie locus.",
        },
    ]
    for key in sorted(metrics):
        if (
            key.startswith("all_class_code_")
            or key.startswith("preliminary_retained_class_code_")
            or key.startswith("retained_class_code_")
            or key.startswith("gene_locus_filter_reason_")
            or key
            in {
                "exact_known_transcript_excluded_rows",
                "missing_feelnc_status_rows",
                "missing_exon_model_rows",
                "feelnc_class_code_mismatch_rows",
                "gene_locus_same_strand_new_contains_known_gene_count",
                "gene_locus_same_strand_new_within_known_gene_count",
                "gene_locus_same_strand_equal_known_gene_count",
                "gene_locus_same_strand_partial_overlap_only_count",
                "gene_locus_opposite_strand_complete_overlap_count",
                "gene_locus_opposite_strand_partial_overlap_count",
                "gene_locus_mixed_or_noncanonical_strand_count",
                "gene_locus_mixed_or_missing_coordinate_count",
                "transcripts_excluded_by_gene_locus_filter",
                "retained_unique_transcript_count",
                "retained_unique_gene_locus_count",
            }
        ):
            rows.append(
                {
                    "metric": key,
                    "category": "distribution",
                    "value": str(metrics[key]),
                    "note": "Step 2.1 transcript-first statistics for same-strand near-complete containment and antisense overlaps.",
                }
            )
    rows.append(
        {
            "metric": "pass_status",
            "category": "all",
            "value": pass_status,
            "note": "PASS indicates only consistency of input joins, row counts, and aggregation; it does not establish biological correctness.",
        }
    )
    return rows


def manifest_row(item: str, path: Path, note: str) -> Dict[str, str]:
    return {
        "item": item,
        "path": str(path),
        "exists": "1" if path.exists() else "0",
        "size_bytes": str(path.stat().st_size) if path.exists() and path.is_file() else "NA",
        "note": note,
    }


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Five-dogs Step 2.1 transcript-first gene-level near-complete containment filter.")
    parser.add_argument("--classification-tsv", type=Path, default=DEFAULT_CLASSIFICATION_TSV)
    parser.add_argument("--transcript-exon-bed", type=Path, default=DEFAULT_TRANSCRIPT_EXON_BED)
    parser.add_argument("--feelnc-status-tsv", type=Path, default=DEFAULT_FEELNC_STATUS_TSV)
    parser.add_argument("--reference-gtf", type=Path, default=DEFAULT_REFERENCE_GTF)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--containment-overlap-fraction", type=float, default=DEFAULT_CONTAINMENT_OVERLAP_FRACTION)
    parser.add_argument("--expected-transcript-count", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def ensure_output_paths(args: argparse.Namespace) -> Dict[str, Path]:
    paths = {
        "all_assessment": args.outdir / "2.1_all_transcript_preASE_assessment.tsv",
        "retained_candidates": args.outdir / "2.1_retained_transcript_candidates.preASE.tsv",
        "transcript_summary": args.outdir / "2.1_candidate_transcript_summary.preASE.tsv",
        "all_gene_locus_assessment": args.outdir / "2.1_all_gene_locus_preASE_assessment.tsv",
        "retained_gene_locus_candidates": args.outdir / "2.1_retained_gene_locus_candidates.preASE.tsv",
        "gene_locus_summary": args.outdir / "2.1_candidate_gene_locus_summary.preASE.tsv",
        "qc_summary": args.outdir / "2.1_transcript_first_preASE_candidate_assessment.qc_summary.tsv",
        "manifest": args.outdir / "2.1_transcript_first_preASE_candidate_assessment.manifest.tsv",
    }
    if not args.overwrite:
        existing = [str(path) for path in paths.values() if path.exists()]
        if existing:
            raise RuntimeError("output file already exists; such as requires overwrite please explicitly provided --overwrite: " + ",".join(existing))
    args.outdir.mkdir(parents=True, exist_ok=True)
    return paths


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    setup_logging()
    logging.info("Step 2.1 transcript-first gene-level near-complete containment filter start")
    if args.containment_overlap_fraction <= 0 or args.containment_overlap_fraction > 1:
        raise RuntimeError(
            f"--containment-overlap-fraction must be within (0, 1]: {args.containment_overlap_fraction}"
        )
    for path, label in [
        (args.classification_tsv, "five-dogs Step 1.6 gffcompare transcript classification"),
        (args.transcript_exon_bed, "five-dogs Step 1.6 transcript exon BED"),
        (args.feelnc_status_tsv, "five-dogs Step 1.6.5c FEELnc module status"),
        (args.reference_gtf, "CanFam3.1 reference GTF"),
    ]:
        require_readable(path, label)
    output_paths = ensure_output_paths(args)

    classification_rows = read_classification(args.classification_tsv)
    if args.expected_transcript_count > 0 and len(classification_rows) != args.expected_transcript_count:
        raise RuntimeError(
            f"Step 1.6 classification transcript count does not equal expected-transcript-count: "
            f"classification={len(classification_rows)} expected={args.expected_transcript_count}"
        )
    feelnc_by_transcript = read_feelnc_status(args.feelnc_status_tsv)
    models_by_transcript = read_transcript_models(args.transcript_exon_bed)
    validate_linked_inputs(classification_rows, feelnc_by_transcript, models_by_transcript)

    all_rows, preliminary_rows, transcript_metrics = build_transcript_assessment_rows(
        classification_rows,
        feelnc_by_transcript,
        models_by_transcript,
    )
    allowed_chroms = {clean(row.get("chr")) for row in preliminary_rows if clean(row.get("chr")) != "NA"}
    reference_genes = read_reference_genes_gtf(args.reference_gtf, allowed_chroms)
    reference_bins = build_reference_bins(reference_genes)
    all_gene_rows, retained_gene_rows, gene_metrics = build_gene_locus_assessment(
        preliminary_rows,
        reference_genes,
        reference_bins,
        args.containment_overlap_fraction,
    )
    gene_rows_by_id = {clean(row.get("query_gene_id")): row for row in all_gene_rows}
    retained_gene_ids = {clean(row.get("query_gene_id")) for row in retained_gene_rows}
    retained_rows, final_metrics = apply_gene_filter_to_transcripts(
        all_rows,
        gene_rows_by_id,
        retained_gene_ids,
    )
    metrics = transcript_metrics + gene_metrics + final_metrics
    transcript_rows = build_transcript_summary(retained_rows)
    qc = qc_rows(
        args,
        classification_rows,
        feelnc_by_transcript,
        models_by_transcript,
        reference_genes,
        all_rows,
        preliminary_rows,
        all_gene_rows,
        retained_gene_rows,
        retained_rows,
        transcript_rows,
        metrics,
    )

    logging.info("Stage 5/5: write Step 2.1 output table")
    logging.info("write all transcript table: %s", output_paths["all_assessment"])
    write_tsv(output_paths["all_assessment"], ASSESSMENT_FIELDS, all_rows)
    logging.info("write retained transcript table: %s", output_paths["retained_candidates"])
    write_tsv(output_paths["retained_candidates"], ASSESSMENT_FIELDS, retained_rows)
    logging.info("write retained transcript summary: %s", output_paths["transcript_summary"])
    write_tsv(output_paths["transcript_summary"], TRANSCRIPT_SUMMARY_FIELDS, transcript_rows)
    logging.info("write all gene/locus assessment: %s", output_paths["all_gene_locus_assessment"])
    write_tsv(output_paths["all_gene_locus_assessment"], GENE_SUMMARY_FIELDS, all_gene_rows)
    logging.info("write retained gene/locus candidates: %s", output_paths["retained_gene_locus_candidates"])
    write_tsv(output_paths["retained_gene_locus_candidates"], GENE_SUMMARY_FIELDS, retained_gene_rows)
    logging.info("write downstream-compatible retained gene/locus summary: %s", output_paths["gene_locus_summary"])
    write_tsv(output_paths["gene_locus_summary"], GENE_SUMMARY_FIELDS, retained_gene_rows)
    logging.info("write QC summary: %s", output_paths["qc_summary"])
    write_tsv(output_paths["qc_summary"], QC_FIELDS, qc)

    manifest = [
        manifest_row("classification_tsv_input", args.classification_tsv, "five-dogs Step 1.6 gffcompare transcript classification"),
        manifest_row("transcript_exon_bed_input", args.transcript_exon_bed, "five-dogs Step 1.6 transcript exon BED"),
        manifest_row("feelnc_status_tsv_input", args.feelnc_status_tsv, "five-dogs Step 1.6.5c all-query FEELnc module status"),
        manifest_row("reference_gtf_input", args.reference_gtf, "CanFam3.1 GTF; transcript/exon records aggregated by gene_id to define known-gene spans."),
        manifest_row("all_transcript_preASE_assessment_tsv", output_paths["all_assessment"], "Step 2.1 transcript table with class-code, same-strand filter and opposite-strand annotation fields"),
        manifest_row("retained_transcript_candidates_tsv", output_paths["retained_candidates"], "Step 2.1 final retained transcripts: retained gene, class_code != '=' and no same-strand complete containment"),
        manifest_row("candidate_transcript_summary_tsv", output_paths["transcript_summary"], "Step 2.1 one-row-per-retained-transcript table"),
        manifest_row("all_gene_locus_assessment_tsv", output_paths["all_gene_locus_assessment"], "Step 2.1 all preliminary gene loci with same-strand containment and opposite-strand overlap annotation"),
        manifest_row("retained_gene_locus_candidates_tsv", output_paths["retained_gene_locus_candidates"], "Step 2.1 retained gene loci after same-strand near-complete containment filter"),
        manifest_row("candidate_gene_locus_summary_tsv", output_paths["gene_locus_summary"], "Step 2.1 retained gene loci; downstream-compatible filename"),
        manifest_row("qc_summary_tsv", output_paths["qc_summary"], "Step 2.1 QC summary"),
        {
            "item": "manifest_tsv",
            "path": str(output_paths["manifest"]),
            "exists": "1",
            "size_bytes": "NA",
            "note": "Step 2.1 manifest; size is finalized after writing",
        },
    ]
    logging.info("write manifest: %s", output_paths["manifest"])
    write_tsv(output_paths["manifest"], MANIFEST_FIELDS, manifest)

    pass_row = next((row for row in qc if row["metric"] == "pass_status"), None)
    if pass_row and pass_row["value"] != "PASS":
        raise RuntimeError(f"Step 2.1 QC failed: {pass_row}")
    logging.info("Step 2.1 transcript-first gene-level near-complete containment filter task completed / done")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except Exception as exc:
        logging.exception("Step 2.1 transcript-first gene-level near-complete containment filter failed: %s", exc)
        raise SystemExit(1)
