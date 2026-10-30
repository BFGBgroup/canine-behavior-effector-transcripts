#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3


import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime


ASE_ROOT = "/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.6bulk-ase-count"
PHASE_ROOT = "/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.10_reconcile_dog11_whatsHap_shapeit5/reconciled_bcf"

REFERENCE_ROOT = "/path/to/project/analysis_20260709/5-dogs/3.ASE-process/3.7_build_canfam31_gene_span_reference"
GENE_SPAN_ANNOTATION = os.path.join(
    REFERENCE_ROOT, "CanFam3.1.autosomes.gene_span.annotation.tsv"
)
EXON_UNION_ANNOTATION = os.path.join(
    REFERENCE_ROOT, "CanFam3.1.autosomes.exon_union.annotation.tsv"
)
MANUAL_NEW_GTF = "/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf"
MANUAL_NEW_GTF_COMPLETION_MARKER = "/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/execution_completion.marker"
BCFTOOLS_BIN = "/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"
CHROMOSOMES = tuple(f"chr{i}" for i in range(1, 39))
AUTOSOME_SET = {str(i) for i in range(1, 39)}
BIN_SIZE = 100000
EXPECTED_MANUAL_GENE_N = 3828
EXPECTED_MANUAL_AUTOSOMAL_GENE_N = 3528
EXPECTED_MANUAL_NON_AUTOSOMAL_GENE_N = 300
EXPECTED_MANUAL_SELECTION_STAGE = "2.1_transcript_first_preASE_candidate_assessment"
GTF_ATTRIBUTE_PATTERN = re.compile(r'(\S+)\s+"([^"]*)"')


def log(message):
    print(f"[{datetime.now():%F %T} ] {message}", file=sys.stderr, flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Integrate ASEReadCounter output, global phasing, and CanFam3.1 gene intervals for one animal and tissue."
    )
    parser.add_argument("animal_id", help="Animal ID, for example Ace_396226.")
    parser.add_argument("term", help="Tissue name, for example cerebellum.")
    parser.add_argument(
        "--output-root",
        required=True,
        help="Absolute output directory for this run; the batch launcher supplies the timestamped directory.",
    )
    args = parser.parse_args()
    if not os.path.isabs(args.output_root):
        parser.error("--output-root must be an absolute path")
    return args


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def normalize_chrom(chrom):
    text = str(chrom or "").strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return text


def canonical_autosome(chrom):
    normalized = normalize_chrom(chrom)
    if normalized not in AUTOSOME_SET:
        return ""
    return f"chr{int(normalized)}"


def variant_key(chrom, pos, ref, alt):
    return (
        normalize_chrom(chrom),
        str(pos).strip(),
        str(ref).strip().upper(),
        str(alt).strip().upper(),
    )


def variant_label(chrom, pos, ref, alt):
    return f"{chrom}:{pos}_{ref}>{alt}"


def require_columns(header, required, path):
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"{path} missing required columns: {', '.join(missing)}")


def parse_nonnegative_int(value, column, path, line_number):
    try:
        parsed = int(str(value).strip())
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number} of {column} is not an integer: {value!r}"
        ) from error
    if parsed < 0:
        raise ValueError(f"{path}:{line_number} of {column} is negative: {parsed}")
    return parsed


def parse_ase_row(row, path, line_number):
    chrom = canonical_autosome(row["contig"])
    try:
        pos = int(str(row["position"]).strip())
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number} of position is not an integer: {row['position']!r}"
        ) from error
    if pos < 1:
        raise ValueError(f"{path}:{line_number} of position invalid: {pos}")

    ref = str(row["refAllele"] or "").strip().upper()
    alt = str(row["altAllele"] or "").strip().upper()
    if not ref or not alt or ref == alt:
        raise ValueError(
            f"{path}:{line_number} of REF/ALT invalid: REF={ref!r}, ALT={alt!r}"
        )

    ref_count = parse_nonnegative_int(row["refCount"], "refCount", path, line_number)
    alt_count = parse_nonnegative_int(row["altCount"], "altCount", path, line_number)
    total_count = parse_nonnegative_int(
        row["totalCount"], "totalCount", path, line_number
    )
    if ref_count + alt_count != total_count:
        raise ValueError(
            f"{path}:{line_number} count mismatch: refCount={ref_count}, "
            f"altCount={alt_count}, totalCount={total_count}"
        )

    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "total_count": total_count,
    }


def build_paths(animal_id, term, output_root):
    ase_table = os.path.join(
        ASE_ROOT,
        animal_id,
        term,
        f"{animal_id}_{term}.ASEReadCounter.table",
    )
    phased_bcfs = [
        os.path.join(
            PHASE_ROOT,
            f"{chrom}.dog11.ref1519.whatsHap_reconciled.phased.bcf",
        )
        for chrom in CHROMOSOMES
    ]
    out_dir = os.path.join(output_root, animal_id, term)
    out_all = os.path.join(out_dir, f"{animal_id}_{term}.ase_haplotype.tsv")
    out_qc = os.path.join(out_dir, f"{animal_id}_{term}.ase_haplotype.qc.tsv")
    return ase_table, phased_bcfs, out_all, out_qc


def load_gene_reference(gene_span_path, exon_union_path):
    span_required = {
        "GeneID",
        "GeneID_base",
        "NCBI_GeneID",
        "GeneSymbol",
        "GeneFeatureType",
        "GeneBiotype",
        "AnnotationSource",
        "SourceGTFGeneIDN",
        "SourceGTFGeneIDs",
        "GeneStrand",
        "CHROM",
        "START1",
        "END1",
        "TranscriptN",
    }
    exon_required = {
        "GeneID",
        "GeneID_base",
        "GeneSymbol",
        "GeneStrand",
        "CHROM",
        "START1",
        "END1",
    }
    gene_metadata = {}
    span_bins = defaultdict(list)
    exon_bins = defaultdict(list)

    with open(gene_span_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"gene-span annotation no column header: {gene_span_path}")
        require_columns(reader.fieldnames, span_required, gene_span_path)

        for line_number, row in enumerate(reader, start=2):
            gene_id = str(row["GeneID"] or "").strip()
            chrom = canonical_autosome(row["CHROM"])
            strand = str(row["GeneStrand"] or "").strip()
            try:
                start = int(row["START1"])
                end = int(row["END1"])
                transcript_n = int(row["TranscriptN"])
                source_gtf_gene_id_n = int(row["SourceGTFGeneIDN"])
            except ValueError as error:
                raise ValueError(
                    f"{gene_span_path}:{line_number} coordinates, TranscriptN "
                    "or SourceGTFGeneIDN invalid"
                ) from error

            if not gene_id or gene_id in gene_metadata:
                raise ValueError(
                    f"{gene_span_path}:{line_number} of GeneID missing or duplicate: {gene_id!r}"
                )
            if not chrom or strand not in {"+", "-"}:
                raise ValueError(
                    f"{gene_span_path}:{line_number} of CHROM/GeneStrand invalid"
                )
            if (
                start < 1
                or end < start
                or transcript_n < 1
                or source_gtf_gene_id_n < 1
            ):
                raise ValueError(
                    f"{gene_span_path}:{line_number} of gene interval, TranscriptN "
                    "or SourceGTFGeneIDN invalid"
                )

            source_gtf_gene_ids = str(row["SourceGTFGeneIDs"] or "").strip()
            source_ids = source_gtf_gene_ids.split("|") if source_gtf_gene_ids else []
            metadata = {
                "GeneID": gene_id,
                "GeneID_base": str(row["GeneID_base"] or "").strip(),
                "NCBI_GeneID": str(row["NCBI_GeneID"] or "").strip(),
                "GeneSymbol": str(row["GeneSymbol"] or "").strip(),
                "GeneFeatureType": str(row["GeneFeatureType"] or "").strip(),
                "GeneBiotype": str(row["GeneBiotype"] or "").strip(),
                "AnnotationSource": str(row["AnnotationSource"] or "").strip(),
                "SourceGTFGeneIDN": source_gtf_gene_id_n,
                "SourceGTFGeneIDs": source_gtf_gene_ids,
                "GeneStrand": strand,
                "CHROM": chrom,
                "START1": start,
                "END1": end,
            }
            if (
                not metadata["GeneID_base"]
                or not metadata["NCBI_GeneID"]
                or not metadata["GeneSymbol"]
                or metadata["GeneFeatureType"] not in {"gene", "pseudogene"}
                or not metadata["GeneBiotype"]
                or not metadata["AnnotationSource"]
                or len(source_ids) != source_gtf_gene_id_n
                or len(source_ids) != len(set(source_ids))
            ):
                raise ValueError(
                    f"{gene_span_path}:{line_number} top-level gene metadata invalid"
                )
            gene_metadata[gene_id] = metadata
            for bin_id in range(start // BIN_SIZE, end // BIN_SIZE + 1):
                span_bins[(chrom, bin_id)].append(gene_id)

    exon_block_n = 0
    with open(exon_union_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"exon-union annotation no column header: {exon_union_path}")
        require_columns(reader.fieldnames, exon_required, exon_union_path)

        for line_number, row in enumerate(reader, start=2):
            gene_id = str(row["GeneID"] or "").strip()
            metadata = gene_metadata.get(gene_id)
            if metadata is None:
                raise ValueError(
                    f"{exon_union_path}:{line_number} of GeneID not in gene span in: {gene_id!r}"
                )
            chrom = canonical_autosome(row["CHROM"])
            strand = str(row["GeneStrand"] or "").strip()
            try:
                start = int(row["START1"])
                end = int(row["END1"])
            except ValueError as error:
                raise ValueError(
                    f"{exon_union_path}:{line_number} coordinates invalid"
                ) from error

            if (
                chrom != metadata["CHROM"]
                or strand != metadata["GeneStrand"]
                or str(row["GeneID_base"] or "").strip() != metadata["GeneID_base"]
                or str(row["GeneSymbol"] or "").strip() != metadata["GeneSymbol"]
            ):
                raise ValueError(
                    f"{exon_union_path}:{line_number}: exon union differs from the gene-span metadata: {gene_id}"
                )
            if start < metadata["START1"] or end > metadata["END1"] or end < start:
                raise ValueError(
                    f"{exon_union_path}:{line_number} of exon block exceeds gene span: {gene_id}"
                )

            for bin_id in range(start // BIN_SIZE, end // BIN_SIZE + 1):
                exon_bins[(chrom, bin_id)].append((start, end, gene_id))
            exon_block_n += 1

    for key in span_bins:
        span_bins[key].sort(
            key=lambda gene_id: (
                gene_metadata[gene_id]["START1"],
                gene_metadata[gene_id]["END1"],
                gene_id,
            )
        )
    for key in exon_bins:
        exon_bins[key].sort(key=lambda record: (record[0], record[1], record[2]))

    log(
        f"[info] load gene region reference: genes={len(gene_metadata)}, "
        f"exon_union_blocks={exon_block_n}"
    )
    return gene_metadata, span_bins, exon_bins


def parse_gtf_attributes(attribute_text):
    return dict(GTF_ATTRIBUTE_PATTERN.findall(attribute_text))


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def add_manual_new_gene_reference(
    gtf_path,
    gene_metadata,
    span_bins,
    exon_bins,
):

    manual_genes = {}
    gtf_rows = 0

    with open(gtf_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gtf_path}:{line_number} is not 9 column GTF record")

            seqname, _source, feature, start_text, end_text, _score, strand, _frame, attribute_text = fields
            if feature not in {"transcript", "exon"}:
                raise ValueError(
                    f"{gtf_path}:{line_number} contains non-transcript/exon feature: {feature!r}"
                )
            attributes = parse_gtf_attributes(attribute_text)
            required_attributes = {
                "gene_id",
                "transcript_id",
                "gene_name",
                "source_dataset",
                "selection_stage",
                "original_gene_id",
            }
            missing = sorted(required_attributes - set(attributes))
            if missing:
                raise ValueError(
                    f"{gtf_path}:{line_number} missing GTF attribute: {','.join(missing)}"
                )

            gene_id = attributes["gene_id"].strip()
            gene_name = attributes["gene_name"].strip()
            transcript_id = attributes["transcript_id"].strip()
            source_dataset = attributes["source_dataset"].strip()
            selection_stage = attributes["selection_stage"].strip()
            original_gene_id = attributes["original_gene_id"].strip()
            if (
                not gene_id.startswith("DOG5_NEWGENE_")
                or gene_name != gene_id
                or not transcript_id
                or source_dataset != "five-dogs"
                or selection_stage != EXPECTED_MANUAL_SELECTION_STAGE
                or not original_gene_id
                or strand not in {"+", "-"}
            ):
                raise ValueError(
                    f"{gtf_path}:{line_number} of manual new gene identifier, source or strand invalid"
                )
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as error:
                raise ValueError(
                    f"{gtf_path}:{line_number} coordinates is not an integer"
                ) from error
            if start < 1 or end < start:
                raise ValueError(f"{gtf_path}:{line_number} coordinates invalid: {start}-{end}")

            chrom = canonical_autosome(seqname)
            record = manual_genes.setdefault(
                gene_id,
                {
                    "GeneSymbol": gene_name,
                    "Seqname": seqname,
                    "CHROM": chrom,
                    "GeneStrand": strand,
                    "START1": start,
                    "END1": end,
                    "TranscriptIDs": set(),
                    "SourceGTFGeneIDs": set(),
                    "Exons": [],
                },
            )
            if (
                record["GeneSymbol"] != gene_name
                or record["Seqname"] != seqname
                or record["CHROM"] != chrom
                or record["GeneStrand"] != strand
            ):
                raise ValueError(
                    f"{gtf_path}:{line_number} of {gene_id} across seqname, strand or gene_name mismatch"
                )
            record["START1"] = min(record["START1"], start)
            record["END1"] = max(record["END1"], end)
            record["TranscriptIDs"].add(transcript_id)
            record["SourceGTFGeneIDs"].add(
                f"{source_dataset}:{original_gene_id}"
            )
            if feature == "exon":
                record["Exons"].append((start, end))
            gtf_rows += 1

    if len(manual_genes) != EXPECTED_MANUAL_GENE_N:
        raise ValueError(
            f"manual new gene expected {EXPECTED_MANUAL_GENE_N} items, observed {len(manual_genes)}"
        )

    autosomal_gene_n = 0
    non_autosomal_gene_n = 0
    exon_union_block_n = 0
    for gene_id, record in sorted(manual_genes.items()):
        if not record["TranscriptIDs"] or not record["Exons"]:
            raise ValueError(f"manual new gene missing transcript or exon: {gene_id}")
        if not record["CHROM"]:
            non_autosomal_gene_n += 1
            continue
        if gene_id in gene_metadata:
            raise ValueError(f"manual new gene ID with known gene conflict: {gene_id}")

        source_ids = sorted(record["SourceGTFGeneIDs"])
        metadata = {
            "GeneID": gene_id,
            "GeneID_base": gene_id,
            "NCBI_GeneID": "NA",
            "GeneSymbol": record["GeneSymbol"],
            "GeneFeatureType": "manual_new_gene",
            "GeneBiotype": "manual_new_gene",
            "AnnotationSource": "ManualNewGene",
            "SourceGTFGeneIDN": len(source_ids),
            "SourceGTFGeneIDs": "|".join(source_ids),
            "GeneStrand": record["GeneStrand"],
            "CHROM": record["CHROM"],
            "START1": record["START1"],
            "END1": record["END1"],
        }
        gene_metadata[gene_id] = metadata
        for bin_id in range(metadata["START1"] // BIN_SIZE, metadata["END1"] // BIN_SIZE + 1):
            span_bins[(metadata["CHROM"], bin_id)].append(gene_id)

        for start, end in merge_intervals(record["Exons"]):
            for bin_id in range(start // BIN_SIZE, end // BIN_SIZE + 1):
                exon_bins[(metadata["CHROM"], bin_id)].append(
                    (start, end, gene_id)
                )
            exon_union_block_n += 1
        autosomal_gene_n += 1

    if autosomal_gene_n != EXPECTED_MANUAL_AUTOSOMAL_GENE_N:
        raise ValueError(
            "manual autosome new gene expected"
            f"{EXPECTED_MANUAL_AUTOSOMAL_GENE_N} items, observed {autosomal_gene_n}"
        )
    if non_autosomal_gene_n != EXPECTED_MANUAL_NON_AUTOSOMAL_GENE_N:
        raise ValueError(
            "manual non-autosomal new gene expected"
            f"{EXPECTED_MANUAL_NON_AUTOSOMAL_GENE_N} items, observed {non_autosomal_gene_n}"
        )

    for key in span_bins:
        span_bins[key].sort(
            key=lambda gene_id: (
                gene_metadata[gene_id]["START1"],
                gene_metadata[gene_id]["END1"],
                gene_id,
            )
        )
    for key in exon_bins:
        exon_bins[key].sort(key=lambda record: (record[0], record[1], record[2]))

    stats = {
        "manual_gtf_rows": gtf_rows,
        "manual_gene_total": len(manual_genes),
        "manual_gene_autosomal": autosomal_gene_n,
        "manual_gene_non_autosomal_skipped": non_autosomal_gene_n,
        "manual_exon_union_blocks": exon_union_block_n,
    }
    log(
        "[info] load manual new gene reference: "
        f"total={stats['manual_gene_total']}, autosomal={autosomal_gene_n}, "
        f"non_autosomal_skipped={non_autosomal_gene_n}, "
        f"exon_union_blocks={exon_union_block_n}"
    )
    return stats


def annotate_variant(chrom, pos, label, gene_metadata, span_bins, exon_bins):
    bin_key = (chrom, pos // BIN_SIZE)
    span_gene_ids = {
        gene_id
        for gene_id in span_bins.get(bin_key, [])
        if gene_metadata[gene_id]["START1"] <= pos <= gene_metadata[gene_id]["END1"]
    }
    exon_gene_ids = {
        gene_id
        for start, end, gene_id in exon_bins.get(bin_key, [])
        if start <= pos <= end
    }
    if not exon_gene_ids.issubset(span_gene_ids):
        raise ValueError(f"{chrom}:{pos} of exon match does not belong to gene-span match")

    ordered_gene_ids = sorted(
        span_gene_ids,
        key=lambda gene_id: (
            gene_metadata[gene_id]["START1"],
            gene_metadata[gene_id]["END1"],
            gene_id,
        ),
    )
    gene_hit_n = len(ordered_gene_ids)
    if gene_hit_n == 0:
        gene_overlap_type = "INTERGENIC"
    elif gene_hit_n == 1:
        gene_overlap_type = "UNIQUE_GENE"
    else:
        strands = {gene_metadata[gene_id]["GeneStrand"] for gene_id in ordered_gene_ids}
        gene_overlap_type = (
            "MULTI_GENE_SAME_STRAND"
            if len(strands) == 1
            else "MULTI_GENE_OPPOSITE_STRAND"
        )

    hits = []
    for gene_id in ordered_gene_ids:
        metadata = gene_metadata[gene_id]
        feature = "exon_union" if gene_id in exon_gene_ids else "intron"
        hits.append(
            {
                **metadata,
                "GeneFeature": feature,
            }
        )

    if not hits:

        unannotated_id = f"(NA){label}"
        hits = [
            {
                "GeneID": unannotated_id,
                "GeneID_base": unannotated_id,
                "NCBI_GeneID": "NA",
                "GeneSymbol": unannotated_id,
                "GeneFeatureType": "unannotated",
                "GeneBiotype": "unannotated",
                "AnnotationSource": "unannotated",
                "SourceGTFGeneIDN": 0,
                "SourceGTFGeneIDs": "NA",
                "GeneStrand": "NA",
                "GeneFeature": "intergenic",
            }
        ]

    all_ids = "|".join(ordered_gene_ids) or hits[0]["GeneID"]

    summary = {
        "AllGeneIDs": all_ids,
        "GeneHitN": gene_hit_n,
        "GeneOverlapType": gene_overlap_type,
    }
    return hits, summary


def collect_ase_variant_keys(ase_table):
    required = [
        "contig",
        "position",
        "refAllele",
        "altAllele",
        "refCount",
        "altCount",
        "totalCount",
    ]
    keys = set()
    input_rows = 0
    eligible_rows = 0
    skipped_non_autosome = 0
    skipped_zero_total = 0

    with open_text(ase_table) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"ASE table no column header: {ase_table}")
        require_columns(reader.fieldnames, required, ase_table)

        for line_number, row in enumerate(reader, start=2):
            input_rows += 1
            parsed = parse_ase_row(row, ase_table, line_number)
            if not parsed["chrom"]:
                skipped_non_autosome += 1
                continue
            if parsed["total_count"] == 0:
                skipped_zero_total += 1
                continue
            key = variant_key(
                parsed["chrom"], parsed["pos"], parsed["ref"], parsed["alt"]
            )
            if key in keys:
                raise ValueError(f"{ase_table}:{line_number} encountered duplicate ASE variant: {key}")
            keys.add(key)
            eligible_rows += 1

    log(
        f"[info] scan ASE table: input_rows={input_rows}, eligible_rows={eligible_rows}, "
        f"skipped_non_autosome={skipped_non_autosome}, "
        f"skipped_zero_total={skipped_zero_total}"
    )
    return keys


def validate_bcf_sample(bcf_path, animal_id):
    result = subprocess.run(
        [BCFTOOLS_BIN, "query", "-l", bcf_path],
        check=True,
        text=True,
        capture_output=True,
    )
    samples = result.stdout.splitlines()
    if animal_id not in samples:
        raise ValueError(f"Sample {animal_id} is absent from the BCF: {bcf_path}")


def load_target_bcf_records(bcf_paths, animal_id, target_keys):
    vcf_map = {}
    matched_positions = set()
    total_records = 0
    matched_records = 0
    target_chromosomes = {key[0] for key in target_keys}
    target_positions = {(key[0], key[1]) for key in target_keys}

    for index, bcf_path in enumerate(bcf_paths, start=1):
        chromosome = CHROMOSOMES[index - 1]
        if normalize_chrom(chromosome) not in target_chromosomes:
            continue

        command = [
            BCFTOOLS_BIN,
            "query",
            "-s",
            animal_id,
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%PS]\n",
            bcf_path,
        ]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
        if process.stdout is None:
            raise RuntimeError(f"cannot read bcftools output: {bcf_path}")

        chromosome_records = 0
        chromosome_matches = 0
        for line in process.stdout:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                process.kill()
                raise ValueError(f"bcftools query row format unexpected: {line.rstrip()}")

            chrom, pos, ref, alts, gt, ps = fields
            chromosome_records += 1
            normalized_position = (normalize_chrom(chrom), str(pos).strip())
            if normalized_position in target_positions:
                matched_positions.add(normalized_position)
            for alt_index, alt in enumerate(alts.split(","), start=1):
                key = variant_key(chrom, pos, ref, alt)
                if key not in target_keys:
                    continue
                if key in vcf_map:
                    process.kill()
                    raise ValueError(f"phase BCF in encountered duplicate exact variant: {key}")
                vcf_map[key] = {
                    "gt": str(gt).strip(),
                    "ps": str(ps).strip(),
                    "alt_index": str(alt_index),
                }
                chromosome_matches += 1

        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"bcftools query failed: {bcf_path}, rc={return_code}")

        total_records += chromosome_records
        matched_records += chromosome_matches
        log(
            f"[phase {index}/{len(bcf_paths)}] {chromosome}: "
            f"records={chromosome_records}, matched_ASE={chromosome_matches}"
        )

    log(
        f"[info] phase BCF scan completed: records={total_records}, "
        f"matched_ASE={matched_records}, matched_positions={len(matched_positions)}"
    )
    return vcf_map, matched_positions


def normalize_ps(ps):
    text = str(ps or "").strip()
    return "NA" if text in {"", ".", "NA"} else text


def count_for_allele(allele_code, alt_index, ref_count, alt_count):
    if allele_code == "0":
        return ref_count
    if allele_code == alt_index:
        return alt_count
    return None


def base_for_allele(allele_code, alt_index, ref, alt):
    if allele_code == "0":
        return ref
    if allele_code == alt_index:
        return alt
    return ""


def assign_haplotype_counts(
    vcf_record,
    coordinate_in_bcf,
    chrom,
    ref,
    alt,
    ref_count,
    alt_count,
):
    if vcf_record is None:
        if coordinate_in_bcf:
            return "", "", "", "", "", "NA", "NA", "REF_ALT_MISMATCH"
        return "", "", "", "", "", "NA", "NA", "NOT_IN_PHASE_BCF"

    gt = vcf_record["gt"]
    ps = normalize_ps(vcf_record["ps"])
    alt_index = vcf_record["alt_index"]
    if not gt or gt in {".", "./.", ".|."}:
        return "", "", "", "", gt, ps, "NA", "MISSING_GT"
    if "|" not in gt:
        return "", "", "", "", gt, ps, "NA", "UNPHASED_GT"

    alleles = gt.split("|")
    if len(alleles) != 2:
        return "", "", "", "", gt, ps, "NA", "INVALID_GT"
    allele1, allele2 = alleles
    if allele1 == allele2:
        return "", "", "", "", gt, ps, "NA", "HOMOZYGOUS_GT"

    hap1_count = count_for_allele(allele1, alt_index, ref_count, alt_count)
    hap2_count = count_for_allele(allele2, alt_index, ref_count, alt_count)
    if hap1_count is None or hap2_count is None:
        return "", "", "", "", gt, ps, "NA", "ALLELE_INDEX_MISMATCH"

    hap1_allele = base_for_allele(allele1, alt_index, ref, alt)
    hap2_allele = base_for_allele(allele2, alt_index, ref, alt)
    return (
        hap1_allele,
        hap1_count,
        hap2_allele,
        hap2_count,
        gt,
        ps,
        chrom,
        "PHASED_HET_EXACT_MATCH",
    )


def output_header(extra_ase_columns):
    base = [
        "Variant",
        "GeneSymbol",
        "GeneID",
        "GeneID_base",
        "NCBI_GeneID",
        "GeneFeatureType",
        "GeneBiotype",
        "AnnotationSource",
        "SourceGTFGeneIDN",
        "SourceGTFGeneIDs",
        "GeneFeature",
        "GeneStrand",
        "GeneOverlapType",
        "AllGeneIDs",
        "GeneHitN",
        "SAMPLE",
        "TISSUE",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "VCF_GT",
        "PS",
        "PHASE_SET",
        "PHASE_STATUS",
        "HAP1_ALLELE",
        "HAP1_COUNT",
        "HAP2_ALLELE",
        "HAP2_COUNT",
        "REF_COUNT",
        "ALT_COUNT",
        "TOTAL_COUNT",
    ]
    return base + extra_ase_columns


def write_output(
    animal_id,
    term,
    ase_table,
    out_all,
    gene_metadata,
    span_bins,
    exon_bins,
    vcf_map,
    vcf_positions,
):
    required = [
        "contig",
        "position",
        "refAllele",
        "altAllele",
        "refCount",
        "altCount",
        "totalCount",
    ]
    stats = defaultdict(int)
    temp_output = f"{out_all}.tmp.{os.getpid()}"

    try:
        with open_text(ase_table) as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"ASE table no column header: {ase_table}")
            require_columns(reader.fieldnames, required, ase_table)
            extra_ase_columns = [x for x in reader.fieldnames if x not in required]
            header = output_header(extra_ase_columns)

            with open(temp_output, "w", encoding="utf-8", newline="") as output_handle:
                writer = csv.DictWriter(
                    output_handle,
                    fieldnames=header,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()

                for line_number, row in enumerate(reader, start=2):
                    stats["input_rows"] += 1
                    parsed = parse_ase_row(row, ase_table, line_number)
                    if not parsed["chrom"]:
                        stats["skipped_non_autosome_rows"] += 1
                        continue
                    if parsed["total_count"] == 0:
                        stats["skipped_zero_total_rows"] += 1
                        continue
                    stats["eligible_autosomal_positive_count_rows"] += 1

                    chrom = parsed["chrom"]
                    pos = parsed["pos"]
                    ref = parsed["ref"]
                    alt = parsed["alt"]
                    ref_count = parsed["ref_count"]
                    alt_count = parsed["alt_count"]
                    total_count = parsed["total_count"]
                    key = variant_key(chrom, pos, ref, alt)

                    (
                        hap1_allele,
                        hap1_count,
                        hap2_allele,
                        hap2_count,
                        gt,
                        ps,
                        phase_set,
                        status,
                    ) = assign_haplotype_counts(
                        vcf_map.get(key),
                        (normalize_chrom(chrom), str(pos)) in vcf_positions,
                        chrom,
                        ref,
                        alt,
                        ref_count,
                        alt_count,
                    )
                    stats[f"phase_status_{status}"] += 1

                    label = variant_label(chrom, pos, ref, alt)
                    hits, annotation = annotate_variant(
                        chrom, pos, label, gene_metadata, span_bins, exon_bins
                    )
                    stats[f"gene_overlap_{annotation['GeneOverlapType']}_variants"] += 1

                    for hit in hits:
                        out_row = {
                            "Variant": label,
                            "GeneSymbol": hit["GeneSymbol"],
                            "GeneID": hit["GeneID"],
                            "GeneID_base": hit["GeneID_base"],
                            "NCBI_GeneID": hit["NCBI_GeneID"],
                            "GeneFeatureType": hit["GeneFeatureType"],
                            "GeneBiotype": hit["GeneBiotype"],
                            "AnnotationSource": hit["AnnotationSource"],
                            "SourceGTFGeneIDN": hit["SourceGTFGeneIDN"],
                            "SourceGTFGeneIDs": hit["SourceGTFGeneIDs"],
                            "GeneFeature": hit["GeneFeature"],
                            "GeneStrand": hit["GeneStrand"],
                            "GeneOverlapType": annotation["GeneOverlapType"],
                            "AllGeneIDs": annotation["AllGeneIDs"],
                            "GeneHitN": annotation["GeneHitN"],
                            "SAMPLE": animal_id,
                            "TISSUE": term,
                            "CHROM": chrom,
                            "POS": pos,
                            "REF": ref,
                            "ALT": alt,
                            "VCF_GT": gt,
                            "PS": ps,
                            "PHASE_SET": phase_set,
                            "PHASE_STATUS": status,
                            "HAP1_ALLELE": hap1_allele,
                            "HAP1_COUNT": hap1_count,
                            "HAP2_ALLELE": hap2_allele,
                            "HAP2_COUNT": hap2_count,
                            "REF_COUNT": ref_count,
                            "ALT_COUNT": alt_count,
                            "TOTAL_COUNT": total_count,
                        }
                        for column in extra_ase_columns:
                            out_row[column] = row.get(column, "")
                        writer.writerow(out_row)
                        stats["output_variant_gene_rows"] += 1
                        stats[f"gene_feature_{hit['GeneFeature']}_rows"] += 1
                        stats[
                            f"annotation_source_{hit['AnnotationSource']}_rows"
                        ] += 1

        os.replace(temp_output, out_all)
    except Exception:
        if os.path.exists(temp_output):
            os.unlink(temp_output)
        raise

    return stats


def write_qc(
    out_qc,
    animal_id,
    term,
    ase_table,
    out_all,
    manual_reference_stats,
    stats,
):
    temp_qc = f"{out_qc}.tmp.{os.getpid()}"
    rows = [
        ("animal_id", animal_id),
        ("term", term),
        ("ase_input", ase_table),
        ("phase_root", PHASE_ROOT),
        ("gene_span_annotation", GENE_SPAN_ANNOTATION),
        ("exon_union_annotation", EXON_UNION_ANNOTATION),
        ("manual_new_gene_gtf", MANUAL_NEW_GTF),
        ("manual_new_gene_gtf_completion_marker", MANUAL_NEW_GTF_COMPLETION_MARKER),
        ("output_table", out_all),
        ("analysis_scope", "chr1-chr38; totalCount>0"),
        ("count_unit", "GATK_ASEReadCounter_site_allele_count"),
        ("variant_match_key", "CHROM+POS+REF+ALT exact match"),
        (
            "unannotated_id_contract",
            "GeneID/GeneID_base/GeneSymbol=(NA)<Variant> when GeneHitN=0",
        ),
        (
            "unannotated_metadata_contract",
            "NCBI_GeneID/SourceGTFGeneIDs=NA; GeneFeatureType/GeneBiotype=unannotated; SourceGTFGeneIDN=0",
        ),
        (
            "gene_annotation_contract",
            "3.7 known genes retain NCBI metadata; manual genes use GeneFeatureType/GeneBiotype=manual_new_gene and NCBI_GeneID=NA",
        ),
    ]
    rows.extend((key, manual_reference_stats[key]) for key in sorted(manual_reference_stats))
    rows.extend((key, stats[key]) for key in sorted(stats))

    try:
        with open(temp_qc, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("metric", "value"))
            writer.writerows(rows)
        os.replace(temp_qc, out_qc)
    except Exception:
        if os.path.exists(temp_qc):
            os.unlink(temp_qc)
        raise


def main():
    args = parse_args()
    animal_id = args.animal_id
    term = args.term
    ase_table, phased_bcfs, out_all, out_qc = build_paths(
        animal_id, term, args.output_root
    )

    required_files = [
        ase_table,
        GENE_SPAN_ANNOTATION,
        EXON_UNION_ANNOTATION,
        MANUAL_NEW_GTF,
        MANUAL_NEW_GTF_COMPLETION_MARKER,
    ]
    required_files.extend(phased_bcfs)
    required_files.extend(f"{path}.csi" for path in phased_bcfs)
    for path in required_files:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    if not os.path.isfile(BCFTOOLS_BIN) or not os.access(BCFTOOLS_BIN, os.X_OK):
        raise FileNotFoundError(f"bcftools not executable: {BCFTOOLS_BIN}")
    existing_outputs = [path for path in (out_all, out_qc) if os.path.exists(path)]
    if existing_outputs:
        raise FileExistsError(
            "refuse to overwrite existing 1.3 result: " + ", ".join(existing_outputs)
        )

    with open(MANUAL_NEW_GTF_COMPLETION_MARKER, "rt", encoding="utf-8") as handle:
        marker_lines = {line.strip() for line in handle if line.strip()}
    if "execution_status=success" not in marker_lines:
        raise ValueError(
            "five-dog-specific manual new gene GTF of completion marker missing execution_status=success: "
            f"{MANUAL_NEW_GTF_COMPLETION_MARKER}"
        )

    os.makedirs(os.path.dirname(out_all), exist_ok=True)
    log(f"[start] animal_id={animal_id} term={term}")
    log(f"[input] ASE_TABLE={ase_table}")
    log(f"[input] PHASE_ROOT={PHASE_ROOT}")
    log(f"[input] GENE_SPAN_ANNOTATION={GENE_SPAN_ANNOTATION}")
    log(f"[input] EXON_UNION_ANNOTATION={EXON_UNION_ANNOTATION}")
    log(f"[input] MANUAL_NEW_GTF={MANUAL_NEW_GTF}")
    log(
        "[input] MANUAL_NEW_GTF_COMPLETION_MARKER="
        f"{MANUAL_NEW_GTF_COMPLETION_MARKER}"
    )
    log(f"[output] OUT_ALL={out_all}")
    log(f"[output] OUT_QC={out_qc}")

    log("[stage 1/6] verify animal sample exists in phase BCF")
    validate_bcf_sample(phased_bcfs[0], animal_id)
    log("[stage 2/6] scan ASE site and restrict autosome valid count")
    target_keys = collect_ase_variant_keys(ase_table)
    log("[stage 3/6] load target site of chromosome level phase")
    vcf_map, vcf_positions = load_target_bcf_records(
        phased_bcfs, animal_id, target_keys
    )
    log("[stage 4/6] load CanFam3.1 known gene gene span with exon union")
    gene_metadata, span_bins, exon_bins = load_gene_reference(
        GENE_SPAN_ANNOTATION, EXON_UNION_ANNOTATION
    )
    log("[stage 5/6] incremental load five-dog-specific manual new gene")
    manual_reference_stats = add_manual_new_gene_reference(
        MANUAL_NEW_GTF,
        gene_metadata,
        span_bins,
        exon_bins,
    )
    log("[stage 6/6] write ASE-haplotype details and QC")
    stats = write_output(
        animal_id,
        term,
        ase_table,
        out_all,
        gene_metadata,
        span_bins,
        exon_bins,
        vcf_map,
        vcf_positions,
    )
    write_qc(
        out_qc,
        animal_id,
        term,
        ase_table,
        out_all,
        manual_reference_stats,
        stats,
    )

    for key in sorted(stats):
        log(f"[summary] {key}={stats[key]}")
    log(f"[done] task completed / done; OUT_ALL={out_all}; OUT_QC={out_qc}")


if __name__ == "__main__":
    main()
