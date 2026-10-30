#!/usr/bin/python3.10
# -*- coding: utf-8 -*-


from __future__ import annotations

import csv
import fcntl
import logging
import shutil
import sys
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Mapping, MutableMapping, Sequence, Tuple


SCRIPT_DIR = Path(
    '/path/to/repository/scripts/03_candidate_catalog/candidate_gtf'
)
LOG_DIR = SCRIPT_DIR / "log"
MASTER_LOG = LOG_DIR / "1.1.build_five_only_new_gene_gtf.log"
SUCCESS_TSV = LOG_DIR / "1.1.build_five_only_new_gene_gtf.success.tsv"
FAILED_TSV = LOG_DIR / "1.1.build_five_only_new_gene_gtf.failed.tsv"
LOCK_FILE = LOG_DIR / "1.1.build_five_only_new_gene_gtf.lock"

BASE_GTF = Path(
    '/path/to/project/legacy/five-dogs/5.cnv/0.gtf/GCF_000002285.3_CanFam3.1_genomic_chr.gene_id_patched.gtf'
)
REFERENCE_FASTA_FAI = Path(
    '/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai'
)
GENE_TABLE = Path(
    '/path/to/project/analysis_20260709/5-dogs/9.new-gene/2.1_transcript_first_preASE_candidate_assessment/2.1_candidate_gene_locus_summary.preASE.tsv'
)
TRANSCRIPT_TABLE = Path(
    '/path/to/project/analysis_20260709/5-dogs/9.new-gene/2.1_transcript_first_preASE_candidate_assessment/2.1_candidate_transcript_summary.preASE.tsv'
)
SOURCE_GTF = Path(
    '/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/gffcmp.stringtie3.annotated.gtf'
)
OUTPUT_ROOT = Path(
    '/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf'
)

EXPECTED_GENE_COUNT = 3828
EXPECTED_TRANSCRIPT_COUNT = 4870
FINAL_GENE_PREFIX = "DOG5_NEWGENE_"
SOURCE_DATASET = "five-dogs"
SELECTION_STAGE = "2.1_transcript_first_preASE_candidate_assessment"

STATUS_FIELDS = [
    "task_id",
    "run_id",
    "status",
    "exit_code",
    "started_at",
    "finished_at",
    "output_path",
    "task_log",
    "reason",
]

GENE_REQUIRED_FIELDS = [
    "query_gene_id",
    "candidate_gene_chr",
    "candidate_gene_start0",
    "candidate_gene_end0",
    "candidate_gene_span_length",
    "candidate_gene_strand",
    "candidate_gene_transcript_strand_list",
    "retained_transcript_count",
    "gene_locus_retained_for_downstream",
    "gene_locus_filter_reason",
    "query_transcript_id_list",
]

TRANSCRIPT_REQUIRED_FIELDS = [
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
    "FEELnc_codpot_category",
    "FEELnc_biotype",
    "candidate_gene_chr",
    "candidate_gene_start0",
    "candidate_gene_end0",
    "candidate_gene_strand",
    "candidate_gene_transcript_count",
    "gene_locus_filter_reason",
    "preASE_retained_for_downstream",
]

TRANSCRIPT_MAP_FIELDS = [
    "final_gene_id",
    "final_transcript_id",
    "source_dataset",
    "selection_stage",
    "original_gene_id",
    "original_transcript_id",
    "chr",
    "strand",
    "gene_start_0based",
    "gene_end_0based",
    "transcript_start_0based",
    "transcript_end_0based",
    "class_code",
    "novelty_category",
    "reference_gene_id",
    "reference_transcript_id",
    "FEELnc_biotype",
    "FEELnc_codpot_category",
    "gene_locus_filter_reason",
    "transcript_row_count",
    "exon_row_count",
    "gtf_row_count",
]


@dataclass(frozen=True)
class GeneRecord:
    original_gene_id: str
    chrom: str
    start0: int
    end0: int
    strand: str
    retained_transcript_count: int
    transcript_ids: Tuple[str, ...]
    filter_reason: str
    final_gene_id: str


@dataclass
class TranscriptRecord:
    original_transcript_id: str
    original_gene_id: str
    chrom: str
    start0: int
    end0: int
    strand: str
    expected_exon_count: int
    class_code: str
    novelty_category: str
    reference_gene_id: str
    reference_transcript_id: str
    feelnc_biotype: str
    feelnc_codpot_category: str
    gene_filter_reason: str
    final_gene_id: str = ""
    final_transcript_id: str = ""
    gtf_rows: List[Tuple[List[str], OrderedDict[str, str]]] = field(default_factory=list)
    transcript_row_count: int = 0
    exon_row_count: int = 0


def now_text() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clean_status_value(value: object) -> str:
    return str(value).replace("\t", " ").replace("\n", " ").strip() or "NA"


def write_status_table(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean_status_value(row.get(field, "NA")) for field in STATUS_FIELDS})


def setup_logging() -> None:
    formatter = logging.Formatter(
        fmt="[%(asctime)s ] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(MASTER_LOG, mode="w", encoding="utf-8")
    stream_handler = logging.StreamHandler(sys.stdout)
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    logging.basicConfig(
        level=logging.INFO,
        handlers=[file_handler, stream_handler],
        force=True,
    )


def require_input_file(path: Path, label: str) -> None:
    if not path.is_absolute():
        raise RuntimeError(f"{label} is not an absolute path: {path}")
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} does not exist or is empty: {path}")


def read_tsv(path: Path, label: str) -> Tuple[List[str], List[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"{label} missing header: {path}")
        return list(reader.fieldnames), list(reader)


def require_columns(header: Sequence[str], required: Sequence[str], label: str) -> None:
    missing = [field for field in required if field not in header]
    if missing:
        raise RuntimeError(f"{label} missing required columns: {','.join(missing)}")


def parse_int(value: str, label: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{label} is not an integer: {value}") from exc


def split_pipe(value: str) -> Tuple[str, ...]:
    if value in {"", "NA"}:
        return tuple()
    return tuple(part for part in value.split("|") if part and part != "NA")


def chromosome_sort_key(chrom: str) -> Tuple[int, int, str]:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return (0, int(suffix), chrom)
        if suffix == "X":
            return (0, 39, chrom)
    return (1, 0, chrom)


def load_genes() -> Tuple[List[GeneRecord], dict[str, GeneRecord]]:
    header, rows = read_tsv(GENE_TABLE, "five dogs Step 2.1 gene table")
    require_columns(header, GENE_REQUIRED_FIELDS, "five dogs Step 2.1 gene table")
    if len(rows) != EXPECTED_GENE_COUNT:
        raise RuntimeError(
            f"Five-dog Step 2.1 gene row count differs from the checked input: observed={len(rows)} "
            f"expected={EXPECTED_GENE_COUNT}"
        )

    parsed: List[Tuple[dict[str, str], int, int, int, Tuple[str, ...]]] = []
    seen_gene_ids: set[str] = set()
    for row in rows:
        gene_id = row["query_gene_id"]
        if not gene_id or gene_id == "NA":
            raise RuntimeError("five dogs Step 2.1 gene table contains empty query_gene_id")
        if gene_id in seen_gene_ids:
            raise RuntimeError(f"five dogs Step 2.1 gene table contains duplicate query_gene_id: {gene_id}")
        seen_gene_ids.add(gene_id)

        if row["gene_locus_retained_for_downstream"] != "1":
            raise RuntimeError(f"candidate gene is not Step 2.1 retain record: {gene_id}")
        if not row["gene_locus_filter_reason"].startswith("retained_"):
            raise RuntimeError(f"The candidate gene filter_reason is not retained: {gene_id}")
        if row["candidate_gene_strand"] not in {"+", "-"}:
            raise RuntimeError(f"candidate gene strand is not+or-: {gene_id}")
        if row["candidate_gene_transcript_strand_list"] != row["candidate_gene_strand"]:
            raise RuntimeError(f"Candidate gene strand differs from its transcript strand list: {gene_id}")

        start0 = parse_int(row["candidate_gene_start0"], f"{gene_id}.candidate_gene_start0")
        end0 = parse_int(row["candidate_gene_end0"], f"{gene_id}.candidate_gene_end0")
        span = parse_int(row["candidate_gene_span_length"], f"{gene_id}.candidate_gene_span_length")
        retained_count = parse_int(row["retained_transcript_count"], f"{gene_id}.retained_transcript_count")
        transcript_ids = split_pipe(row["query_transcript_id_list"])
        if start0 < 0 or end0 <= start0 or span != end0 - start0:
            raise RuntimeError(f"candidate gene coordinates or span invalid: {gene_id}")
        if not transcript_ids or retained_count != len(transcript_ids):
            raise RuntimeError(f"The candidate gene retained-transcript count does not match its ID list: {gene_id}")
        if len(set(transcript_ids)) != len(transcript_ids):
            raise RuntimeError(f"The candidate gene transcript ID list contains duplicates: {gene_id}")
        parsed.append((row, start0, end0, retained_count, transcript_ids))

    parsed.sort(
        key=lambda item: (
            chromosome_sort_key(item[0]["candidate_gene_chr"]),
            item[1],
            item[2],
            item[0]["query_gene_id"],
        )
    )
    genes: List[GeneRecord] = []
    for serial, (row, start0, end0, retained_count, transcript_ids) in enumerate(parsed, start=1):
        genes.append(
            GeneRecord(
                original_gene_id=row["query_gene_id"],
                chrom=row["candidate_gene_chr"],
                start0=start0,
                end0=end0,
                strand=row["candidate_gene_strand"],
                retained_transcript_count=retained_count,
                transcript_ids=transcript_ids,
                filter_reason=row["gene_locus_filter_reason"],
                final_gene_id=f"{FINAL_GENE_PREFIX}{serial:05d}",
            )
        )
    return genes, {gene.original_gene_id: gene for gene in genes}


def load_transcripts(
    genes: Sequence[GeneRecord],
    genes_by_id: Mapping[str, GeneRecord],
) -> Tuple[List[TranscriptRecord], dict[str, TranscriptRecord]]:
    header, rows = read_tsv(TRANSCRIPT_TABLE, "five dogs Step 2.1 transcript table")
    require_columns(header, TRANSCRIPT_REQUIRED_FIELDS, "five dogs Step 2.1 transcript table")
    if len(rows) != EXPECTED_TRANSCRIPT_COUNT:
        raise RuntimeError(
            f"Five-dog Step 2.1 transcript row count differs from the checked input: observed={len(rows)} "
            f"expected={EXPECTED_TRANSCRIPT_COUNT}"
        )

    transcripts_by_gene: MutableMapping[str, List[TranscriptRecord]] = defaultdict(list)
    seen_transcript_ids: set[str] = set()
    for row in rows:
        transcript_id = row["query_transcript_id"]
        gene_id = row["query_gene_id"]
        if not transcript_id or transcript_id == "NA":
            raise RuntimeError("five dogs Step 2.1 transcript table contains empty query_transcript_id")
        if transcript_id in seen_transcript_ids:
            raise RuntimeError(f"five dogs Step 2.1 transcript table contains duplicate ID: {transcript_id}")
        seen_transcript_ids.add(transcript_id)
        gene = genes_by_id.get(gene_id)
        if gene is None:
            raise RuntimeError(f"transcript references gene in the table does not exist of query_gene_id: {transcript_id} {gene_id}")
        if transcript_id not in gene.transcript_ids:
            raise RuntimeError(f"transcript absent from gene of query_transcript_id_list: {transcript_id}")
        if row["preASE_retained_for_downstream"] != "1":
            raise RuntimeError(f"transcript is not Step 2.1 retain record: {transcript_id}")
        if row["class_code"] == "=" or not row["class_code"]:
            raise RuntimeError(f"candidate transcript contains empty class_code or exact reference matching=: {transcript_id}")
        if row["transcript_strand"] not in {"+", "-"}:
            raise RuntimeError(f"transcript strand is not+or-: {transcript_id}")

        start0 = parse_int(row["transcript_span_start0"], f"{transcript_id}.transcript_span_start0")
        end0 = parse_int(row["transcript_span_end0"], f"{transcript_id}.transcript_span_end0")
        span = parse_int(row["transcript_span_length"], f"{transcript_id}.transcript_span_length")
        exon_count = parse_int(row["transcript_exon_count"], f"{transcript_id}.transcript_exon_count")
        candidate_gene_start0 = parse_int(row["candidate_gene_start0"], f"{transcript_id}.candidate_gene_start0")
        candidate_gene_end0 = parse_int(row["candidate_gene_end0"], f"{transcript_id}.candidate_gene_end0")
        candidate_gene_count = parse_int(
            row["candidate_gene_transcript_count"],
            f"{transcript_id}.candidate_gene_transcript_count",
        )
        if start0 < 0 or end0 <= start0 or span != end0 - start0 or exon_count < 1:
            raise RuntimeError(f"candidate transcript coordinates, span or exon count invalid: {transcript_id}")
        if row["chr"] != gene.chrom or row["candidate_gene_chr"] != gene.chrom:
            raise RuntimeError(f"Transcript and gene chromosomes differ: {transcript_id}")
        if row["transcript_strand"] != gene.strand or row["candidate_gene_strand"] != gene.strand:
            raise RuntimeError(f"Transcript and gene strands differ: {transcript_id}")
        if candidate_gene_start0 != gene.start0 or candidate_gene_end0 != gene.end0:
            raise RuntimeError(f"Candidate gene interval in the transcript table differs from the gene table: {transcript_id}")
        if start0 < gene.start0 or end0 > gene.end0:
            raise RuntimeError(f"transcript interval exceeds belonging to candidate gene: {transcript_id}")
        if candidate_gene_count != gene.retained_transcript_count:
            raise RuntimeError(f"Candidate transcript count in the transcript table differs from the gene table: {transcript_id}")
        if row["gene_locus_filter_reason"] != gene.filter_reason:
            raise RuntimeError(f"Transcript and gene filter_reason values differ: {transcript_id}")

        transcripts_by_gene[gene_id].append(
            TranscriptRecord(
                original_transcript_id=transcript_id,
                original_gene_id=gene_id,
                chrom=row["chr"],
                start0=start0,
                end0=end0,
                strand=row["transcript_strand"],
                expected_exon_count=exon_count,
                class_code=row["class_code"],
                novelty_category=row["novelty_category"] or "NA",
                reference_gene_id=row["reference_gene_id"] or "NA",
                reference_transcript_id=row["reference_transcript_id"] or "NA",
                feelnc_biotype=row["FEELnc_biotype"] or "NA",
                feelnc_codpot_category=row["FEELnc_codpot_category"] or "NA",
                gene_filter_reason=row["gene_locus_filter_reason"],
            )
        )

    transcripts: List[TranscriptRecord] = []
    for gene in genes:
        gene_transcripts = sorted(
            transcripts_by_gene.get(gene.original_gene_id, []),
            key=lambda record: record.original_transcript_id,
        )
        observed_ids = {record.original_transcript_id for record in gene_transcripts}
        if observed_ids != set(gene.transcript_ids):
            raise RuntimeError(f"ID sets differ between the gene and transcript tables: {gene.original_gene_id}")
        for serial, record in enumerate(gene_transcripts, start=1):
            record.final_gene_id = gene.final_gene_id
            record.final_transcript_id = f"{gene.final_gene_id}.T{serial:03d}"
            transcripts.append(record)

    if len(transcripts) != EXPECTED_TRANSCRIPT_COUNT:
        raise RuntimeError("gene-transcript after joining transcript count changed")
    return transcripts, {record.original_transcript_id: record for record in transcripts}


def parse_gtf_attributes(text: str) -> OrderedDict[str, str]:
    attributes: OrderedDict[str, str] = OrderedDict()
    for raw_item in text.strip().rstrip(";").split(";"):
        item = raw_item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        value = value.strip()
        if value.startswith('"') and value.endswith('"'):
            value = value[1:-1]
        attributes[key] = value
    return attributes


def clean_attribute_value(value: object) -> str:
    text = str(value).replace('"', "'").replace(";", ",").replace("\t", " ").replace("\n", " ").strip()
    return text or "NA"


def format_gtf_attributes(attributes: Mapping[str, object]) -> str:
    return " ".join(
        f'{key} "{clean_attribute_value(value)}";'
        for key, value in attributes.items()
    )


def collect_source_gtf_rows(transcripts_by_id: Mapping[str, TranscriptRecord]) -> None:
    with SOURCE_GTF.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            if columns[1] != "StringTie" or columns[2] not in {"transcript", "exon"}:
                continue
            attributes = parse_gtf_attributes(columns[8])
            transcript_id = attributes.get("transcript_id", "")
            record = transcripts_by_id.get(transcript_id)
            if record is None:
                continue
            gene_id = attributes.get("gene_id", "")
            if gene_id != record.original_gene_id:
                raise RuntimeError(
                    f"Source GTF gene_id differs from Step 2.1: line={line_number} transcript={transcript_id}"
                )
            if columns[0] != record.chrom or columns[6] != record.strand:
                raise RuntimeError(
                    f"Source GTF chromosome or strand differs from Step 2.1: line={line_number} transcript={transcript_id}"
                )
            row_start0 = parse_int(columns[3], f"source_gtf.line_{line_number}.start") - 1
            row_end0 = parse_int(columns[4], f"source_gtf.line_{line_number}.end")
            if columns[2] == "transcript":
                if row_start0 != record.start0 or row_end0 != record.end0:
                    raise RuntimeError(
                        f"Source GTF transcript interval differs from Step 2.1: line={line_number} transcript={transcript_id}"
                    )
                record.transcript_row_count += 1
            else:
                if row_start0 < record.start0 or row_end0 > record.end0 or row_end0 <= row_start0:
                    raise RuntimeError(
                        f"source GTF exon interval exceeds transcript: line={line_number} transcript={transcript_id}"
                    )
                record.exon_row_count += 1
            record.gtf_rows.append((columns, attributes))
            if line_number % 250000 == 0:
                logging.info("scanned source GTF row count=%d", line_number)

    missing = [
        record.original_transcript_id
        for record in transcripts_by_id.values()
        if record.transcript_row_count != 1
        or record.exon_row_count != record.expected_exon_count
    ]
    if missing:
        raise RuntimeError(
            "subset candidate transcript no parse to one records transcript and expected exon count, example="
            + ",".join(sorted(missing)[:10])
        )


def inspect_base_gtf() -> None:
    with BASE_GTF.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            if FINAL_GENE_PREFIX in columns[8]:
                raise RuntimeError(f"baseline GTF already contains {FINAL_GENE_PREFIX}: line={line_number}")


def inspect_reference_sequences(candidate_chromosomes: Iterable[str]) -> None:
    reference_sequences: set[str] = set()
    with REFERENCE_FASTA_FAI.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2 or not columns[0]:
                raise RuntimeError(f"reference FASTA FAI format unexpected: line={line_number}")
            reference_sequences.add(columns[0])
    missing = sorted(set(candidate_chromosomes) - reference_sequences)
    if missing:
        raise RuntimeError("candidate gene chromosome/contig does not exist in CanFam3.1 FASTA FAI: " + ",".join(missing[:10]))


def make_new_gtf_line(
    columns: Sequence[str],
    old_attributes: Mapping[str, str],
    record: TranscriptRecord,
) -> str:
    new_attributes: OrderedDict[str, object] = OrderedDict()
    new_attributes["gene_id"] = record.final_gene_id
    new_attributes["transcript_id"] = record.final_transcript_id
    new_attributes["gene_name"] = record.final_gene_id
    new_attributes["source_dataset"] = SOURCE_DATASET
    new_attributes["selection_stage"] = SELECTION_STAGE
    new_attributes["original_gene_id"] = record.original_gene_id
    new_attributes["original_transcript_id"] = record.original_transcript_id
    new_attributes["original_gtf_source"] = columns[1]
    new_attributes["class_code"] = record.class_code
    new_attributes["novelty_category"] = record.novelty_category
    new_attributes["reference_gene_id"] = record.reference_gene_id
    new_attributes["reference_transcript_id"] = record.reference_transcript_id
    new_attributes["FEELnc_biotype"] = record.feelnc_biotype
    new_attributes["FEELnc_codpot_category"] = record.feelnc_codpot_category
    new_attributes["gene_locus_filter_reason"] = record.gene_filter_reason
    if "xloc" in old_attributes:
        new_attributes["xloc"] = old_attributes["xloc"]
    if columns[2] == "exon" and "exon_number" in old_attributes:
        new_attributes["exon_number"] = old_attributes["exon_number"]

    new_columns = list(columns)
    new_columns[1] = "ManualNewGene"
    new_columns[8] = format_gtf_attributes(new_attributes)
    return "\t".join(new_columns)


def write_new_only_gtf(path: Path, transcripts: Sequence[TranscriptRecord]) -> int:
    row_count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# five-dogs-only Step 2.1 retained candidate gene-loci\n")
        handle.write("# source rows: StringTie transcript/exon records from five-dogs gffcompare GTF\n")
        for record in transcripts:
            for columns, old_attributes in record.gtf_rows:
                handle.write(make_new_gtf_line(columns, old_attributes, record))
                handle.write("\n")
                row_count += 1
    return row_count


def write_combined_gtf(new_only_gtf: Path, combined_gtf: Path) -> None:
    with combined_gtf.open("w", encoding="utf-8", newline="") as output_handle:
        with BASE_GTF.open("r", encoding="utf-8") as base_handle:
            shutil.copyfileobj(base_handle, output_handle)
        output_handle.write("\n# appended five-dogs-only Step 2.1 retained candidate gene-loci\n")
        with new_only_gtf.open("r", encoding="utf-8") as new_handle:
            shutil.copyfileobj(new_handle, output_handle)


def write_transcript_map(path: Path, transcripts: Sequence[TranscriptRecord], genes_by_id: Mapping[str, GeneRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRANSCRIPT_MAP_FIELDS, delimiter="\t")
        writer.writeheader()
        for record in transcripts:
            gene = genes_by_id[record.original_gene_id]
            writer.writerow(
                {
                    "final_gene_id": record.final_gene_id,
                    "final_transcript_id": record.final_transcript_id,
                    "source_dataset": SOURCE_DATASET,
                    "selection_stage": SELECTION_STAGE,
                    "original_gene_id": record.original_gene_id,
                    "original_transcript_id": record.original_transcript_id,
                    "chr": record.chrom,
                    "strand": record.strand,
                    "gene_start_0based": gene.start0,
                    "gene_end_0based": gene.end0,
                    "transcript_start_0based": record.start0,
                    "transcript_end_0based": record.end0,
                    "class_code": record.class_code,
                    "novelty_category": record.novelty_category,
                    "reference_gene_id": record.reference_gene_id,
                    "reference_transcript_id": record.reference_transcript_id,
                    "FEELnc_biotype": record.feelnc_biotype,
                    "FEELnc_codpot_category": record.feelnc_codpot_category,
                    "gene_locus_filter_reason": record.gene_filter_reason,
                    "transcript_row_count": record.transcript_row_count,
                    "exon_row_count": record.exon_row_count,
                    "gtf_row_count": len(record.gtf_rows),
                }
            )


def write_summary(path: Path, metrics: Mapping[str, object]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, value])


def run_workflow(run_id: str, run_dir: Path) -> Path:
    logging.info("Stage 1/5: check five-dog-specific input file")
    for path, label in [
        (BASE_GTF, "Reference GTF"),
        (REFERENCE_FASTA_FAI, "CanFam3.1 FASTA FAI"),
        (GENE_TABLE, "Five-dog Step 2.1 gene table"),
        (TRANSCRIPT_TABLE, "Five-dog Step 2.1 transcript table"),
        (SOURCE_GTF, "Five-dog gffcompare annotated GTF"),
    ]:
        require_input_file(path, label)
        logging.info("input=%s", path)

    logging.info("Stage 2/5: load and join five dogs Step 2.1 gene/transcript candidate")
    genes, genes_by_id = load_genes()
    transcripts, transcripts_by_id = load_transcripts(genes, genes_by_id)
    logging.info("candidate_gene_count=%d", len(genes))
    logging.info("candidate_transcript_count=%d", len(transcripts))

    logging.info("Stage 3/5: check reference sequences/baseline GTF and extract five dogs StringTie transcript/exon structure")
    inspect_reference_sequences(gene.chrom for gene in genes)
    inspect_base_gtf()
    collect_source_gtf_rows(transcripts_by_id)

    logging.info("Stage 4/5: create this run timestamp output directory and write GTF with provenance table")
    run_dir.mkdir(parents=True, exist_ok=False)
    new_only_gtf = run_dir / "five_only_new_gene_loci.new_only.gtf"
    combined_gtf = run_dir / "GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf"
    transcript_map = run_dir / "five_only_new_gene_loci.transcript_map.tsv"
    summary_tsv = run_dir / "five_only_new_gene_loci.summary.tsv"
    completion_marker = run_dir / "execution_completion.marker"

    new_only_row_count = write_new_only_gtf(new_only_gtf, transcripts)
    write_combined_gtf(new_only_gtf, combined_gtf)
    write_transcript_map(transcript_map, transcripts, genes_by_id)
    metrics: OrderedDict[str, object] = OrderedDict(
        [
            ("run_id", run_id),
            ("source_dataset", SOURCE_DATASET),
            ("selection_stage", SELECTION_STAGE),
            ("base_gtf", BASE_GTF),
            ("reference_fasta_fai", REFERENCE_FASTA_FAI),
            ("candidate_gene_table", GENE_TABLE),
            ("candidate_transcript_table", TRANSCRIPT_TABLE),
            ("source_gtf", SOURCE_GTF),
            ("candidate_gene_count", len(genes)),
            ("candidate_transcript_count", len(transcripts)),
            ("appended_gtf_row_count", new_only_row_count),
            ("final_gene_id_prefix", FINAL_GENE_PREFIX),
            ("new_only_gtf", new_only_gtf),
            ("combined_gtf", combined_gtf),
            ("transcript_map", transcript_map),
        ]
    )
    write_summary(summary_tsv, metrics)

    logging.info("Stage 5/5: check this run main output path")
    logging.info("new_only_gtf=%s", new_only_gtf)
    logging.info("combined_gtf=%s", combined_gtf)
    logging.info("transcript_map=%s", transcript_map)
    logging.info("summary=%s", summary_tsv)
    return completion_marker


def main() -> int:
    if len(sys.argv) != 1:
        print("this script does not accept command-line arguments, please direct run fixed five dogs input version.", file=sys.stderr)
        return 2

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with LOCK_FILE.open("a+", encoding="utf-8") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(f"existing the same step task in run, no new task was started: {LOCK_FILE}", file=sys.stderr)
            return 3

        setup_logging()
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = OUTPUT_ROOT / run_id
        started_at = now_text()
        task_id = "build_five_only_new_gene_gtf"
        write_status_table(SUCCESS_TSV, [])
        write_status_table(FAILED_TSV, [])

        logging.info("task started")
        logging.info("run_id=%s", run_id)
        logging.info("output_dir=%s", run_dir)
        logging.info("master_log=%s", MASTER_LOG)
        try:
            completion_marker = run_workflow(run_id, run_dir)
            finished_at = now_text()
            write_status_table(
                SUCCESS_TSV,
                [
                    {
                        "task_id": task_id,
                        "run_id": run_id,
                        "status": "success",
                        "exit_code": 0,
                        "started_at": started_at,
                        "finished_at": finished_at,
                        "output_path": run_dir,
                        "task_log": MASTER_LOG,
                        "reason": "main_workflow_completed",
                    }
                ],
            )
            completion_marker.write_text(
                f"run_id={run_id}\nexecution_status=success\n",
                encoding="utf-8",
            )
        except BaseException as exc:
            finished_at = now_text()
            logging.exception("task failed: %s", exc)
            write_status_table(SUCCESS_TSV, [])
            write_status_table(
                FAILED_TSV,
                [
                    {
                        "task_id": task_id,
                        "run_id": run_id,
                        "status": "failed",
                        "exit_code": 1,
                        "started_at": started_at,
                        "finished_at": finished_at,
                        "output_path": run_dir,
                        "task_log": MASTER_LOG,
                        "reason": exc,
                    }
                ],
            )
            logging.info("execution_status=failed")
            return 1

        logging.info("execution_completion_marker=%s", completion_marker)
        logging.info("task completed / done")
        logging.info("execution_status=success")
        logging.info("reminder: execution completion does not establish output completeness, format, counts and biological quality already passed, please check separately.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
