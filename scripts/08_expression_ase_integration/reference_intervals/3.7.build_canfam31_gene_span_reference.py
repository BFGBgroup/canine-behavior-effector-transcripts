#!/usr/bin/env python3


import csv
import gzip
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from urllib.parse import unquote


GTF_FILE = "/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf"
GFF_FILE = "/path/to/project/legacy/five-dogs/7.Allelome.LINK/reference/GCF_000002285.3_CanFam3.1_genomic.gff.gz"
OUT_ROOT = "/path/to/project/analysis_20260709/5-dogs/3.ASE-process/3.7_build_canfam31_gene_span_reference"
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 39))
AUTOSOME_ORDER = {chrom: index for index, chrom in enumerate(AUTOSOMES, start=1)}

GENE_SPAN_BED_OUT = os.path.join(OUT_ROOT, "CanFam3.1.autosomes.gene_span.bed")
GENE_SPAN_ANNOTATION_OUT = os.path.join(
    OUT_ROOT, "CanFam3.1.autosomes.gene_span.annotation.tsv"
)
EXON_UNION_BED_OUT = os.path.join(OUT_ROOT, "CanFam3.1.autosomes.exon_union.bed")
EXON_UNION_ANNOTATION_OUT = os.path.join(
    OUT_ROOT, "CanFam3.1.autosomes.exon_union.annotation.tsv"
)
GENE_SPAN_OVERLAP_OUT = os.path.join(
    OUT_ROOT, "CanFam3.1.autosomes.gene_span_overlaps.tsv"
)
QC_OUT = os.path.join(OUT_ROOT, "CanFam3.1.autosomes.gene_span.qc.tsv")
OUTPUTS = (
    GENE_SPAN_BED_OUT,
    GENE_SPAN_ANNOTATION_OUT,
    EXON_UNION_BED_OUT,
    EXON_UNION_ANNOTATION_OUT,
    GENE_SPAN_OVERLAP_OUT,
    QC_OUT,
)

GTF_ATTRIBUTE_PATTERN = re.compile(r'(\S+)\s+"([^"]*)"')

def log(message):
    print(f"[{datetime.now():%F %T} ] {message}", file=sys.stderr, flush=True)


def strip_gene_version(gene_id):
    return re.sub(r"\.\d+$", "", gene_id)


def parse_gff_attributes(attribute_text):
    values = {}
    for item in attribute_text.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        values[key] = unquote(value)
    return values


def parse_gtf_attributes(attribute_text):
    return dict(GTF_ATTRIBUTE_PATTERN.findall(attribute_text))


def ncbi_gene_id(values):
    for item in values.get("Dbxref", "").split(","):
        if item.startswith("GeneID:"):
            return item.split(":", 1)[1]
    return "NA"


def require_paths():
    if not os.path.isfile(GTF_FILE):
        raise FileNotFoundError(f"GTF not found: {GTF_FILE}")
    if not os.path.isfile(GFF_FILE):
        raise FileNotFoundError(f"GFF3 not found: {GFF_FILE}")

    existing_outputs = [path for path in OUTPUTS if os.path.exists(path)]
    if existing_outputs:
        raise FileExistsError(
            "Refusing to overwrite existing gene-region output(s): "
            + ", ".join(existing_outputs)
        )


def parse_interval(fields, line_number, invalid_examples):
    try:
        start1 = int(fields[3])
        end1 = int(fields[4])
    except ValueError:
        invalid_examples.append(f"line {line_number}: non-integer coordinate")
        return None

    if start1 < 1 or end1 < start1:
        invalid_examples.append(
            f"line {line_number}: invalid interval {start1}-{end1}"
        )
        return None
    return start1, end1


def open_gff():
    if GFF_FILE.endswith(".gz"):
        return gzip.open(GFF_FILE, "rt", encoding="utf-8")
    return open(GFF_FILE, "rt", encoding="utf-8")


def load_autosome_seqids():
    seqid_to_chrom = {}
    with open_gff() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "region":
                continue
            values = parse_gff_attributes(fields[8])
            chromosome = values.get("chromosome", "")
            chrom = f"chr{chromosome}"
            if chrom in AUTOSOME_ORDER:
                seqid_to_chrom[fields[0]] = chrom

    missing = sorted(set(AUTOSOMES) - set(seqid_to_chrom.values()))
    if missing:
        raise ValueError(
            "GFF3 region features do not define all expected autosomes: "
            + ", ".join(missing)
        )
    return seqid_to_chrom


def resolve_top_genes(feature_id, parents_by_id, top_gene_ids, cache, visiting):
    if feature_id in cache:
        return cache[feature_id]
    if feature_id in top_gene_ids:
        cache[feature_id] = {feature_id}
        return cache[feature_id]
    if feature_id in visiting:
        raise ValueError(f"Cycle detected in GFF3 Parent hierarchy: {feature_id}")

    visiting.add(feature_id)
    resolved = set()
    for parent_id in parents_by_id.get(feature_id, set()):
        resolved.update(
            resolve_top_genes(
                parent_id, parents_by_id, top_gene_ids, cache, visiting
            )
        )
    visiting.remove(feature_id)
    cache[feature_id] = resolved
    return resolved


def load_gff_hierarchy(seqid_to_chrom, invalid_examples, qc):
    top_gene_records = {}
    parents_by_id = defaultdict(set)
    feature_types = {}

    with open_gff() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                qc["nonstandard_field_count_rows"] += 1
                continue

            chrom = seqid_to_chrom.get(fields[0])
            if chrom is None:
                qc["gff_non_autosome_feature_rows"] += 1
                continue
            feature = fields[2]
            values = parse_gff_attributes(fields[8])
            feature_id = values.get("ID", "").strip()
            parent_ids = {
                parent.strip()
                for parent in values.get("Parent", "").split(",")
                if parent.strip()
            }
            strand = fields[6]
            interval = parse_interval(fields, line_number, invalid_examples)
            if interval is None:
                continue
            start1, end1 = interval

            if feature_id:
                previous_type = feature_types.setdefault(feature_id, feature)
                if previous_type != feature:
                    invalid_examples.append(
                        f"line {line_number}: feature ID {feature_id} has types "
                        f"{previous_type!r} and {feature!r}"
                    )
                parents_by_id[feature_id].update(parent_ids)

            if feature in {"gene", "pseudogene"} and not parent_ids:
                qc[f"autosome_top_level_{feature}_rows"] += 1
                if not feature_id:
                    invalid_examples.append(
                        f"line {line_number}: top-level {feature} missing ID"
                    )
                    continue
                if strand not in {"+", "-"}:
                    invalid_examples.append(
                        f"line {line_number}: top-level {feature_id} has invalid strand={strand!r}"
                    )
                    continue
                gene_symbol = (
                    values.get("gene", "").strip()
                    or values.get("Name", "").strip()
                    or feature_id
                )
                gene_symbol_source = (
                    "gene" if values.get("gene", "").strip() else
                    "Name" if values.get("Name", "").strip() else
                    "gene_id_fallback"
                )
                record = {
                    "GeneID": feature_id,
                    "GeneID_base": strip_gene_version(feature_id),
                    "NCBI_GeneID": ncbi_gene_id(values),
                    "GeneSymbol": gene_symbol,
                    "GeneSymbolSource": gene_symbol_source,
                    "GeneFeatureType": feature,
                    "GeneBiotype": values.get("gene_biotype", "NA") or "NA",
                    "AnnotationSource": fields[1],
                    "GeneStrand": strand,
                    "CHROM": chrom,
                    "GFF_START1": start1,
                    "GFF_END1": end1,
                }
                if (
                    feature_id in top_gene_records
                    and top_gene_records[feature_id] != record
                ):
                    invalid_examples.append(
                        f"line {line_number}: conflicting top-level gene ID={feature_id}"
                    )
                else:
                    top_gene_records[feature_id] = record

    if not top_gene_records:
        raise ValueError(f"No autosomal top-level genes found in GFF3: {GFF_FILE}")
    return top_gene_records, parents_by_id


def resolve_one_top_gene(
    feature_id, line_number, parents_by_id, top_gene_ids, cache, invalid_examples
):
    top_genes = resolve_top_genes(
        feature_id, parents_by_id, top_gene_ids, cache, set()
    )
    if len(top_genes) != 1:
        invalid_examples.append(
            f"GTF line {line_number}: transcript_id={feature_id!r} maps to "
            f"{len(top_genes)} top-level genes"
        )
        return None
    return next(iter(top_genes))


def load_gene_regions():
    seqid_to_chrom = load_autosome_seqids()
    qc = defaultdict(int)
    invalid_examples = []
    top_gene_records, parents_by_id = load_gff_hierarchy(
        seqid_to_chrom, invalid_examples, qc
    )

    top_gene_ids = set(top_gene_records)
    cache = {}
    gene_records = {}
    transcript_ids_by_gene = defaultdict(set)
    source_gtf_gene_ids_by_gene = defaultdict(set)
    remapped_source_gene_ids = set()
    raw_exon_rows = []
    exon_rows = defaultdict(list)

    with open(GTF_FILE, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                qc["nonstandard_field_count_rows"] += 1
                continue
            feature = fields[2]
            if feature not in {"transcript", "exon"}:
                continue
            chrom = fields[0]
            if chrom not in AUTOSOME_ORDER:
                qc[f"gtf_non_autosome_{feature}_rows"] += 1
                continue
            qc[f"autosome_{feature}_rows"] += 1

            interval = parse_interval(fields, line_number, invalid_examples)
            if interval is None:
                continue
            start1, end1 = interval
            strand = fields[6]
            values = parse_gtf_attributes(fields[8])
            transcript_id = values.get("transcript_id", "").strip()
            source_gene_id = values.get("gene_id", "").strip()
            source_gene_symbol = values.get("gene_name", "").strip()
            if not transcript_id:
                invalid_examples.append(
                    f"GTF line {line_number}: {feature} missing transcript_id"
                )
                continue
            if strand not in {"+", "-"}:
                invalid_examples.append(
                    f"GTF line {line_number}: invalid strand={strand!r}"
                )
                continue

            if feature == "exon":
                raw_exon_rows.append(
                    {
                        "TranscriptID": transcript_id,
                        "START1": start1,
                        "END1": end1,
                        "CHROM": chrom,
                        "GeneStrand": strand,
                        "LineNumber": line_number,
                    }
                )
                continue

            gene_id = resolve_one_top_gene(
                transcript_id,
                line_number,
                parents_by_id,
                top_gene_ids,
                cache,
                invalid_examples,
            )
            if gene_id is None:
                continue
            metadata = top_gene_records[gene_id]
            if metadata["CHROM"] != chrom or metadata["GeneStrand"] != strand:
                invalid_examples.append(
                    f"GTF line {line_number}: transcript_id={transcript_id!r} "
                    "has CHROM/strand inconsistent with its top-level gene"
                )
                continue
            if source_gene_symbol and source_gene_symbol != metadata["GeneSymbol"]:
                qc["autosome_transcript_rows_symbol_differs_from_top_gene"] += 1

            record = gene_records.get(gene_id)
            if record is None:
                record = {
                    **metadata,
                    "START1": start1,
                    "END1": end1,
                    "TranscriptN": 0,
                    "SourceGTFGeneIDN": 0,
                    "SourceGTFGeneIDs": "",
                }
                gene_records[gene_id] = record
            else:
                record["START1"] = min(record["START1"], start1)
                record["END1"] = max(record["END1"], end1)

            transcript_ids_by_gene[gene_id].add(transcript_id)
            if source_gene_id:
                source_gtf_gene_ids_by_gene[gene_id].add(source_gene_id)
                if source_gene_id != gene_id:
                    remapped_source_gene_ids.add(source_gene_id)

    for row in raw_exon_rows:
        gene_id = resolve_one_top_gene(
            row["TranscriptID"],
            row["LineNumber"],
            parents_by_id,
            top_gene_ids,
            cache,
            invalid_examples,
        )
        if gene_id is None:
            continue
        if gene_id not in gene_records:
            invalid_examples.append(
                f"GTF line {row['LineNumber']}: exon maps to {gene_id}, "
                "but no matching transcript row was found"
            )
            continue
        exon_rows[gene_id].append(row)

    for gene_id, record in gene_records.items():
        record["TranscriptN"] = len(transcript_ids_by_gene.get(gene_id, set()))
        source_ids = sorted(source_gtf_gene_ids_by_gene.get(gene_id, set()))
        record["SourceGTFGeneIDN"] = len(source_ids)
        record["SourceGTFGeneIDs"] = "|".join(source_ids)

    exon_union = build_exon_union(gene_records, exon_rows, invalid_examples)
    if invalid_examples:
        preview = "; ".join(invalid_examples[:10])
        raise ValueError(
            "Cannot build an unambiguous top-level gene reference; "
            f"{len(invalid_examples)} invalid/conflicting rows. {preview}"
        )

    qc["autosome_gene_n"] = len(gene_records)
    qc["autosome_gff_top_level_gene_n"] = len(top_gene_records)
    qc["autosome_gff_top_level_gene_without_gtf_transcript_n"] = (
        len(top_gene_records) - len(gene_records)
    )
    qc["autosome_gtf_gene_ids_remapped_to_top_gene_n"] = len(
        remapped_source_gene_ids
    )
    qc["autosome_genes_with_multiple_source_gtf_gene_ids_n"] = sum(
        record["SourceGTFGeneIDN"] > 1 for record in gene_records.values()
    )
    qc["autosome_gene_feature_n"] = sum(
        record["GeneFeatureType"] == "gene" for record in gene_records.values()
    )
    qc["autosome_pseudogene_feature_n"] = sum(
        record["GeneFeatureType"] == "pseudogene"
        for record in gene_records.values()
    )
    qc["autosome_plus_strand_gene_n"] = sum(
        record["GeneStrand"] == "+" for record in gene_records.values()
    )
    qc["autosome_minus_strand_gene_n"] = sum(
        record["GeneStrand"] == "-" for record in gene_records.values()
    )
    qc["autosome_gene_symbol_fallback_n"] = sum(
        record["GeneSymbolSource"] == "gene_id_fallback"
        for record in gene_records.values()
    )
    qc["autosome_gene_with_exon_n"] = sum(
        bool(blocks) for blocks in exon_union.values()
    )
    qc["autosome_gene_without_exon_n"] = sum(
        not blocks for blocks in exon_union.values()
    )
    qc["autosome_exon_union_block_n"] = sum(
        len(blocks) for blocks in exon_union.values()
    )
    qc["autosome_exon_union_bp"] = sum(
        block["END1"] - block["START1"] + 1
        for blocks in exon_union.values()
        for block in blocks
    )
    return gene_records, exon_union, qc


def build_exon_union(gene_records, exon_rows, invalid_examples):
    exon_union = {}

    for gene_id, record in gene_records.items():
        rows = exon_rows.get(gene_id, [])
        if not rows:

            record["ExonRowN"] = 0
            record["ExonUnionBlockN"] = 0
            record["ExonUnionBp"] = 0
            exon_union[gene_id] = []
            continue

        intervals = []
        for row in rows:
            conflicts = []
            if row["CHROM"] != record["CHROM"]:
                conflicts.append("CHROM")
            if row["GeneStrand"] != record["GeneStrand"]:
                conflicts.append("GeneStrand")
            if conflicts:
                invalid_examples.append(
                    f"line {row['LineNumber']}: exon gene_id={gene_id} has conflicting "
                    + ",".join(conflicts)
                )
                continue
            intervals.append((row["START1"], row["END1"]))

        intervals.sort()
        if intervals:
            if (
                intervals[0][0] < record["GFF_START1"]
                or max(end1 for _, end1 in intervals) > record["GFF_END1"]
            ):
                invalid_examples.append(
                    f"descendant exon lies outside official top-level span: {gene_id}"
                )

            record["START1"] = min(record["START1"], intervals[0][0])
            record["END1"] = max(
                record["END1"], max(end1 for _, end1 in intervals)
            )

        merged = []
        for start1, end1 in intervals:
            if not merged or start1 > merged[-1]["END1"] + 1:
                merged.append(
                    {"START1": start1, "END1": end1, "SourceExonRowN": 1}
                )
            else:
                merged[-1]["END1"] = max(merged[-1]["END1"], end1)
                merged[-1]["SourceExonRowN"] += 1

        record["ExonRowN"] = len(rows)
        record["ExonUnionBlockN"] = len(merged)
        record["ExonUnionBp"] = sum(
            block["END1"] - block["START1"] + 1 for block in merged
        )
        exon_union[gene_id] = merged

    return exon_union


def ordered_gene_records(gene_records):
    return sorted(
        gene_records.values(),
        key=lambda record: (
            AUTOSOME_ORDER[record["CHROM"]],
            record["START1"],
            record["END1"],
            record["GeneID"],
        ),
    )


def classify_span_overlap(gene1, gene2):
    if gene1["START1"] == gene2["START1"] and gene1["END1"] == gene2["END1"]:
        return "identical_span"
    if gene1["START1"] <= gene2["START1"] and gene1["END1"] >= gene2["END1"]:
        return "gene1_contains_gene2"
    if gene2["START1"] <= gene1["START1"] and gene2["END1"] >= gene1["END1"]:
        return "gene2_contains_gene1"
    return "partial_overlap"


def build_gene_span_overlaps(ordered_records, qc):
    overlap_rows = []
    active_by_chrom = defaultdict(list)
    overlapping_gene_ids = set()

    for gene2 in ordered_records:
        chrom = gene2["CHROM"]
        active = [
            gene1
            for gene1 in active_by_chrom[chrom]
            if gene1["END1"] >= gene2["START1"]
        ]

        for gene1 in active:
            overlap_start1 = max(gene1["START1"], gene2["START1"])
            overlap_end1 = min(gene1["END1"], gene2["END1"])
            if overlap_start1 > overlap_end1:
                continue

            strand_relationship = (
                "same_strand"
                if gene1["GeneStrand"] == gene2["GeneStrand"]
                else "opposite_strand"
            )
            overlap_type = classify_span_overlap(gene1, gene2)
            overlap_rows.append(
                {
                    "CHROM": chrom,
                    "Gene1ID": gene1["GeneID"],
                    "Gene1ID_base": gene1["GeneID_base"],
                    "Gene1Symbol": gene1["GeneSymbol"],
                    "Gene1Strand": gene1["GeneStrand"],
                    "Gene1Start1": gene1["START1"],
                    "Gene1End1": gene1["END1"],
                    "Gene2ID": gene2["GeneID"],
                    "Gene2ID_base": gene2["GeneID_base"],
                    "Gene2Symbol": gene2["GeneSymbol"],
                    "Gene2Strand": gene2["GeneStrand"],
                    "Gene2Start1": gene2["START1"],
                    "Gene2End1": gene2["END1"],
                    "OverlapStart1": overlap_start1,
                    "OverlapEnd1": overlap_end1,
                    "OverlapBp": overlap_end1 - overlap_start1 + 1,
                    "StrandRelationship": strand_relationship,
                    "SpanOverlapType": overlap_type,
                }
            )
            overlapping_gene_ids.update((gene1["GeneID"], gene2["GeneID"]))
            qc[f"gene_span_overlap_{strand_relationship}_pair_n"] += 1
            qc[f"gene_span_overlap_{overlap_type}_pair_n"] += 1

        active.append(gene2)
        active_by_chrom[chrom] = active

    qc["gene_span_overlap_pair_n"] = len(overlap_rows)
    qc["gene_span_overlapping_gene_n"] = len(overlapping_gene_ids)
    return overlap_rows


def write_gene_span_outputs(ordered_records):
    log(f"[stage] writing gene-span BED: {GENE_SPAN_BED_OUT}")
    with open(GENE_SPAN_BED_OUT, "xt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for record in ordered_records:
            writer.writerow(
                (
                    record["CHROM"],
                    record["START1"] - 1,
                    record["END1"],
                    record["GeneID"],
                    0,
                    record["GeneStrand"],
                )
            )

    fields = (
        "GeneID",
        "GeneID_base",
        "NCBI_GeneID",
        "GeneSymbol",
        "GeneSymbolSource",
        "GeneFeatureType",
        "GeneBiotype",
        "AnnotationSource",
        "SourceGTFGeneIDN",
        "SourceGTFGeneIDs",
        "GeneStrand",
        "CHROM",
        "GFF_START1",
        "GFF_END1",
        "START1",
        "END1",
        "START0",
        "END0",
        "SpanBp",
        "TranscriptN",
        "ExonRowN",
        "ExonUnionBlockN",
        "ExonUnionBp",
    )
    log(f"[stage] writing gene-span annotation: {GENE_SPAN_ANNOTATION_OUT}")
    with open(
        GENE_SPAN_ANNOTATION_OUT, "xt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in ordered_records:
            writer.writerow(
                {
                    **record,
                    "START0": record["START1"] - 1,
                    "END0": record["END1"],
                    "SpanBp": record["END1"] - record["START1"] + 1,
                }
            )


def write_exon_union_outputs(ordered_records, exon_union):
    fields = (
        "ExonBlockID",
        "GeneID",
        "GeneID_base",
        "GeneSymbol",
        "GeneStrand",
        "CHROM",
        "BlockIndex",
        "START1",
        "END1",
        "START0",
        "END0",
        "BlockBp",
        "SourceExonRowN",
    )

    log(f"[stage] writing exon-union BED: {EXON_UNION_BED_OUT}")
    log(f"[stage] writing exon-union annotation: {EXON_UNION_ANNOTATION_OUT}")
    with open(EXON_UNION_BED_OUT, "xt", encoding="utf-8", newline="") as bed_handle, open(
        EXON_UNION_ANNOTATION_OUT, "xt", encoding="utf-8", newline=""
    ) as annotation_handle:
        bed_writer = csv.writer(bed_handle, delimiter="\t", lineterminator="\n")
        annotation_writer = csv.DictWriter(
            annotation_handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        annotation_writer.writeheader()

        for record in ordered_records:
            blocks = exon_union[record["GeneID"]]
            for block_index, block in enumerate(blocks, start=1):
                block_id = (
                    f"{record['GeneID_base']}|exon_union|{block_index:04d}"
                )
                bed_writer.writerow(
                    (
                        record["CHROM"],
                        block["START1"] - 1,
                        block["END1"],
                        block_id,
                        0,
                        record["GeneStrand"],
                    )
                )
                annotation_writer.writerow(
                    {
                        "ExonBlockID": block_id,
                        "GeneID": record["GeneID"],
                        "GeneID_base": record["GeneID_base"],
                        "GeneSymbol": record["GeneSymbol"],
                        "GeneStrand": record["GeneStrand"],
                        "CHROM": record["CHROM"],
                        "BlockIndex": block_index,
                        "START1": block["START1"],
                        "END1": block["END1"],
                        "START0": block["START1"] - 1,
                        "END0": block["END1"],
                        "BlockBp": block["END1"] - block["START1"] + 1,
                        "SourceExonRowN": block["SourceExonRowN"],
                    }
                )


def write_overlap_output(overlap_rows):
    fields = (
        "CHROM",
        "Gene1ID",
        "Gene1ID_base",
        "Gene1Symbol",
        "Gene1Strand",
        "Gene1Start1",
        "Gene1End1",
        "Gene2ID",
        "Gene2ID_base",
        "Gene2Symbol",
        "Gene2Strand",
        "Gene2Start1",
        "Gene2End1",
        "OverlapStart1",
        "OverlapEnd1",
        "OverlapBp",
        "StrandRelationship",
        "SpanOverlapType",
    )
    log(f"[stage] writing gene-span overlaps: {GENE_SPAN_OVERLAP_OUT}")
    with open(GENE_SPAN_OVERLAP_OUT, "xt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(overlap_rows)


def write_qc(qc):
    qc_rows = (
        ("input_gtf", GTF_FILE),
        ("input_gff3", GFF_FILE),
        ("autosomes", ",".join(AUTOSOMES)),
        (
            "gene_span_definition",
            "min-max of GTF transcript/exon records after GFF3 top-level mapping",
        ),
        ("gene_identity", "GFF3 top-level ID; NCBI GeneID retained separately"),
        (
            "coordinate_contract",
            "GTF/GFF3/TSV=1-based-inclusive; BED=0-based-half-open",
        ),
        ("autosome_top_level_gene_rows", qc["autosome_top_level_gene_rows"]),
        (
            "autosome_top_level_pseudogene_rows",
            qc["autosome_top_level_pseudogene_rows"],
        ),
        ("autosome_transcript_rows", qc["autosome_transcript_rows"]),
        ("autosome_exon_rows", qc["autosome_exon_rows"]),
        (
            "autosome_transcript_rows_symbol_differs_from_top_gene",
            qc["autosome_transcript_rows_symbol_differs_from_top_gene"],
        ),
        (
            "autosome_gff_top_level_gene_n",
            qc["autosome_gff_top_level_gene_n"],
        ),
        (
            "autosome_gff_top_level_gene_without_gtf_transcript_n",
            qc["autosome_gff_top_level_gene_without_gtf_transcript_n"],
        ),
        (
            "autosome_gtf_gene_ids_remapped_to_top_gene_n",
            qc["autosome_gtf_gene_ids_remapped_to_top_gene_n"],
        ),
        (
            "autosome_genes_with_multiple_source_gtf_gene_ids_n",
            qc["autosome_genes_with_multiple_source_gtf_gene_ids_n"],
        ),
        ("autosome_gene_n", qc["autosome_gene_n"]),
        ("autosome_gene_feature_n", qc["autosome_gene_feature_n"]),
        (
            "autosome_pseudogene_feature_n",
            qc["autosome_pseudogene_feature_n"],
        ),
        ("autosome_plus_strand_gene_n", qc["autosome_plus_strand_gene_n"]),
        ("autosome_minus_strand_gene_n", qc["autosome_minus_strand_gene_n"]),
        ("autosome_gene_with_exon_n", qc["autosome_gene_with_exon_n"]),
        ("autosome_gene_without_exon_n", qc["autosome_gene_without_exon_n"]),
        ("autosome_exon_union_block_n", qc["autosome_exon_union_block_n"]),
        ("autosome_exon_union_bp", qc["autosome_exon_union_bp"]),
        (
            "gtf_non_autosome_transcript_rows",
            qc["gtf_non_autosome_transcript_rows"],
        ),
        ("gtf_non_autosome_exon_rows", qc["gtf_non_autosome_exon_rows"]),
        ("gff_non_autosome_feature_rows", qc["gff_non_autosome_feature_rows"]),
        ("autosome_gene_symbol_fallback_n", qc["autosome_gene_symbol_fallback_n"]),
        ("nonstandard_field_count_rows", qc["nonstandard_field_count_rows"]),
        ("gene_span_overlap_pair_n", qc["gene_span_overlap_pair_n"]),
        ("gene_span_overlapping_gene_n", qc["gene_span_overlapping_gene_n"]),
        (
            "gene_span_overlap_same_strand_pair_n",
            qc["gene_span_overlap_same_strand_pair_n"],
        ),
        (
            "gene_span_overlap_opposite_strand_pair_n",
            qc["gene_span_overlap_opposite_strand_pair_n"],
        ),
        (
            "gene_span_overlap_identical_span_pair_n",
            qc["gene_span_overlap_identical_span_pair_n"],
        ),
        (
            "gene_span_overlap_gene1_contains_gene2_pair_n",
            qc["gene_span_overlap_gene1_contains_gene2_pair_n"],
        ),
        (
            "gene_span_overlap_gene2_contains_gene1_pair_n",
            qc["gene_span_overlap_gene2_contains_gene1_pair_n"],
        ),
        (
            "gene_span_overlap_partial_overlap_pair_n",
            qc["gene_span_overlap_partial_overlap_pair_n"],
        ),
        ("gene_span_bed_output", GENE_SPAN_BED_OUT),
        ("gene_span_annotation_output", GENE_SPAN_ANNOTATION_OUT),
        ("exon_union_bed_output", EXON_UNION_BED_OUT),
        ("exon_union_annotation_output", EXON_UNION_ANNOTATION_OUT),
        ("gene_span_overlap_output", GENE_SPAN_OVERLAP_OUT),
    )

    log(f"[stage] writing QC table: {QC_OUT}")
    with open(QC_OUT, "xt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerows(qc_rows)


def main():
    log("[start] build CanFam3.1 autosomal gene-region reference")
    require_paths()
    log(f"[stage] reading GFF3 top-level gene hierarchy from: {GFF_FILE}")
    log(f"[stage] reading GTF transcript/exon coordinates from: {GTF_FILE}")
    gene_records, exon_union, qc = load_gene_regions()
    ordered_records = ordered_gene_records(gene_records)
    log(
        f"[stage] validated {qc['autosome_gene_n']} top-level genes, "
        f"including {qc['autosome_pseudogene_feature_n']} pseudogenes; "
        f"transcript_rows={qc['autosome_transcript_rows']}; "
        f"exon_rows={qc['autosome_exon_rows']}; "
        f"remapped_gtf_gene_ids="
        f"{qc['autosome_gtf_gene_ids_remapped_to_top_gene_n']}"
    )
    log("[stage] classifying pairwise gene-span overlaps")
    overlap_rows = build_gene_span_overlaps(ordered_records, qc)

    os.makedirs(OUT_ROOT, exist_ok=True)
    write_gene_span_outputs(ordered_records)
    write_exon_union_outputs(ordered_records, exon_union)
    write_overlap_output(overlap_rows)
    write_qc(qc)

    log(
        f"[summary] genes={qc['autosome_gene_n']}; "
        f"exon_union_blocks={qc['autosome_exon_union_block_n']}; "
        f"gene_span_overlap_pairs={qc['gene_span_overlap_pair_n']}"
    )
    log(f"[summary] output_root={OUT_ROOT}")
    log("[done] CanFam3.1 autosomal gene-region reference completed")


if __name__ == "__main__":
    main()
