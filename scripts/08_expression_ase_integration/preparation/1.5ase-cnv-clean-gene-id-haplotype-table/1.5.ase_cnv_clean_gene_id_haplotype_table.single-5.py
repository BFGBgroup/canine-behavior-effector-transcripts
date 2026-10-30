#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3


import argparse
import csv
import os
import sys
import tempfile
from collections import OrderedDict
from datetime import datetime


EXPECTED_PHASE_STATUS = "PHASED_HET_EXACT_MATCH"
ALLOWED_CNV_STATUSES = {
    "OK",
    "NO_STEP7_SEGMENT_RESULT",
    "NO_STEP7_SAMPLE_RESULT",
    "INCOMPLETE_STEP7_SCALE_SET",
}
EXPECTED_PHASE_SET_SCOPE = "chromosome"
COUNT_UNIT = "ASEReadCounter_site_allele_support_sum"
AGGREGATION_METHOD = "sum_distinct_variant_allele_support_within_gene_id"
ALLOWED_CNV_FILTER_ACTIONS = {
    "keep",
    "flag_no_cnv_consensus",
    "flag_soft_cnv_overlap",
    "flag_strong_cnv_overlap",
    "flag_cnv_conflict",
    "flag_not_evaluable",
}
CNV_ACTION_COUNT_COLUMNS = {
    "keep": "keep_variant_count",
    "flag_no_cnv_consensus": "flag_no_cnv_consensus_variant_count",
    "flag_soft_cnv_overlap": "flag_soft_cnv_overlap_variant_count",
    "flag_strong_cnv_overlap": "flag_strong_cnv_overlap_variant_count",
    "flag_cnv_conflict": "flag_cnv_conflict_variant_count",
    "flag_not_evaluable": "flag_not_evaluable_variant_count",
}
GROUP_METADATA_FIELDS = {
    "gene_id_base": "GeneID_base",
    "ncbi_gene_id": "NCBI_GeneID",
    "gene_feature_type": "GeneFeatureType",
    "gene_biotype": "GeneBiotype",
    "annotation_source": "AnnotationSource",
    "source_gtf_gene_ids": "SourceGTFGeneIDs",
    "gene_strand": "GeneStrand",
}
TRUE_VALUES = {"TRUE", "True", "true", "1", "YES", "Yes", "yes", "Y", "y"}
FALSE_VALUES = {"FALSE", "False", "false", "0", "NO", "No", "no", "N", "n"}

REQUIRED_COLS = [
    "Variant",
    "GeneSymbol",
    "GeneID",
    "GeneID_base",
    "NCBI_GeneID",
    "GeneFeatureType",
    "GeneBiotype",
    "AnnotationSource",
    "SourceGTFGeneIDs",
    "GeneFeature",
    "GeneStrand",
    "GeneHitN",
    "SAMPLE",
    "TISSUE",
    "CHROM",
    "VCF_GT",
    "PHASE_SET",
    "PHASE_STATUS",
    "HAP1_ALLELE",
    "HAP1_COUNT",
    "HAP2_ALLELE",
    "HAP2_COUNT",
    "TOTAL_COUNT",
    "cnv_result_status",
    "cnv_evaluable",
    "cnv_filter_action",
    "cnv_clean_pass",
]

OUTPUT_HEADER = [
    "group_id",
    "record_type",
    "gene_id",
    "gene_symbol",
    "gene_id_base",
    "ncbi_gene_id",
    "gene_feature_type",
    "gene_biotype",
    "annotation_source",
    "source_gtf_gene_ids",
    "gene_strand",
    "chromosome",
    "sample",
    "tissue",
    "variant_count",
    "unambiguous_variant_count",
    "multi_gene_variant_count",
    "keep_variant_count",
    "flag_no_cnv_consensus_variant_count",
    "flag_soft_cnv_overlap_variant_count",
    "flag_strong_cnv_overlap_variant_count",
    "flag_cnv_conflict_variant_count",
    "flag_not_evaluable_variant_count",
    "hap1_count_sum",
    "hap2_count_sum",
    "total_count_sum",
    "unambiguous_hap1_count_sum",
    "unambiguous_hap2_count_sum",
    "multi_gene_hap1_count_sum",
    "multi_gene_hap2_count_sum",
    "count_unit",
    "aggregation_method",
]


def log(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp} ] {message}", file=sys.stderr, flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate a step 1.4 CNV-clean ASE long-format table by GeneID+CHROM+SAMPLE+TISSUE."
    )
    parser.add_argument("animal_id", help="Animal ID, for example Ace_396226.")
    parser.add_argument("term", help="Tissue name, for example cerebellum.")
    parser.add_argument(
        "--input-tsv",
        required=True,
        help="Absolute path to the per-sample ase_haplotype.cnv_clean.tsv produced by step 1.4.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Absolute output directory for this task.",
    )
    return parser.parse_args()


def build_paths(animal_id, term, input_tsv, out_dir):
    sample_id = f"{animal_id}_{term}"
    input_file = os.path.abspath(input_tsv)
    output_dir = os.path.abspath(out_dir)
    output_file = os.path.join(output_dir, f"{sample_id}.cnv_clean.gene_id_haplotype.tsv")
    qc_file = os.path.join(output_dir, f"{sample_id}.cnv_clean.gene_id_haplotype.qc.tsv")
    return input_file, output_file, qc_file


def clean_value(value, default="NA"):
    text = "" if value is None else str(value).strip()
    return text if text else default


def parse_bool(value, column, row_number):
    text = clean_value(value, default="")
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"row {row_number} {column} is not a valid boolean value: {text!r}")


def parse_nonnegative_int(value, column, row_number):
    text = clean_value(value, default="")
    try:
        number = int(text)
    except ValueError as exc:
        raise ValueError(f"row {row_number} {column} is not an integer: {text!r}") from exc
    if number < 0:
        raise ValueError(f"row {row_number} {column} cannot be negative: {number}")
    return number


def parse_gene_hit_n(value, row_number):
    return parse_nonnegative_int(value, "GeneHitN", row_number)


def require_columns(fieldnames, input_file):
    missing = [column for column in REQUIRED_COLS if column not in fieldnames]
    if missing:
        raise ValueError(f"input file missing required columns: {input_file}; missing={missing}")


def add_token(store, value):
    token = clean_value(value)
    if token not in {"NA", "nan", "None"}:
        store[token] = None


def is_pseudo_gene_id(gene_id):
    return gene_id.startswith("(NA)")


def validate_annotation_contract(row, row_number):
    variant = clean_value(row["Variant"])
    gene_id = clean_value(row["GeneID"])
    gene_id_base = clean_value(row["GeneID_base"])
    gene_symbol = clean_value(row["GeneSymbol"])
    ncbi_gene_id = clean_value(row["NCBI_GeneID"])
    feature_type = clean_value(row["GeneFeatureType"])
    gene_biotype = clean_value(row["GeneBiotype"])
    annotation_source = clean_value(row["AnnotationSource"])
    source_gtf_gene_ids = clean_value(row["SourceGTFGeneIDs"])
    gene_strand = clean_value(row["GeneStrand"])
    gene_hit_n = parse_gene_hit_n(row["GeneHitN"], row_number)

    if gene_hit_n == 0:
        expected = f"(NA){variant}"
        observed = {
            "GeneID": gene_id,
            "GeneID_base": gene_id_base,
            "GeneSymbol": gene_symbol,
        }
        invalid = {key: value for key, value in observed.items() if value != expected}
        if invalid:
            raise ValueError(
                f"row {row_number} contains an unannotated synthetic ID that violates the expected format: expected={expected}; observed={invalid}"
            )
        if clean_value(row["GeneFeature"]) != "intergenic":
            raise ValueError(f"row {row_number} GeneHitN=0 but GeneFeature is not intergenic")
        if (
            ncbi_gene_id != "NA"
            or feature_type != "unannotated"
            or gene_biotype != "unannotated"
            or annotation_source != "unannotated"
            or source_gtf_gene_ids != "NA"
            or gene_strand != "NA"
        ):
            raise ValueError(f"row {row_number} contains an unannotated record inconsistent with the required gene metadata")
        return "unannotated_variant", gene_hit_n

    if is_pseudo_gene_id(gene_id):
        raise ValueError(f"row {row_number} GeneHitN>0 but uses unannotated synthetic GeneID: {gene_id}")

    if feature_type == "manual_new_gene":
        if (
            not gene_id.startswith("DOG5_NEWGENE_")
            or gene_id_base != gene_id
            or gene_symbol != gene_id
            or ncbi_gene_id != "NA"
            or gene_biotype != "manual_new_gene"
            or annotation_source != "ManualNewGene"
            or source_gtf_gene_ids == "NA"
            or gene_strand not in {"+", "-"}
        ):
            raise ValueError(f"row {row_number} contains manual new-gene metadata inconsistent with the expected format: {gene_id}")
        return "manual_new_gene", gene_hit_n

    if (
        feature_type not in {"gene", "pseudogene"}
        or ncbi_gene_id == "NA"
        or gene_biotype in {"NA", "unannotated", "manual_new_gene"}
        or annotation_source in {"NA", "unannotated", "ManualNewGene"}
        or source_gtf_gene_ids == "NA"
        or gene_strand not in {"+", "-"}
    ):
        raise ValueError(f"row {row_number} contains known-gene metadata inconsistent with the expected format: {gene_id}")
    return "known_gene", gene_hit_n


def new_group(row, record_type):
    group = {
        "record_type": record_type,
        "gene_id": clean_value(row["GeneID"]),
        "gene_symbols": OrderedDict(),
        "chromosome": clean_value(row["CHROM"]),
        "sample": clean_value(row["SAMPLE"]),
        "tissue": clean_value(row["TISSUE"]),
        "variants": OrderedDict(),
        "input_row_count": 0,
        "deduplicated_row_count": 0,
    }
    for output_field, input_field in GROUP_METADATA_FIELDS.items():
        group[output_field] = clean_value(row[input_field])
    return group


def add_row_to_group(group, row, row_number, gene_hit_n):
    variant = clean_value(row["Variant"])
    hap1_count = parse_nonnegative_int(row["HAP1_COUNT"], "HAP1_COUNT", row_number)
    hap2_count = parse_nonnegative_int(row["HAP2_COUNT"], "HAP2_COUNT", row_number)
    total_count = parse_nonnegative_int(row["TOTAL_COUNT"], "TOTAL_COUNT", row_number)
    if hap1_count + hap2_count != total_count:
        raise ValueError(
            f"row {row_number} HAP1_COUNT+HAP2_COUNT does not equal TOTAL_COUNT: "
            f"{hap1_count}+{hap2_count}!={total_count}"
        )
    cnv_filter_action = clean_value(row["cnv_filter_action"])
    evidence = (hap1_count, hap2_count, gene_hit_n, cnv_filter_action)

    for output_field, input_field in GROUP_METADATA_FIELDS.items():
        current_value = clean_value(row[input_field])
        if group[output_field] != current_value:
            raise ValueError(
                f"row {row_number} has conflicting {input_field} values for the same GeneID: "
                f"GeneID={group['gene_id']}; previous={group[output_field]!r}; "
                f"current={current_value!r}"
            )

    group["input_row_count"] += 1
    add_token(group["gene_symbols"], row["GeneSymbol"])

    if variant in group["variants"]:
        previous_evidence = group["variants"][variant]
        if previous_evidence != evidence:
            raise ValueError(
                f"row {row_number} has conflicting counts or GeneHitN values for the same GeneID+Variant: "
                f"GeneID={group['gene_id']}; Variant={variant}; "
                f"previous={previous_evidence}; current={evidence}"
            )
        group["deduplicated_row_count"] += 1
        return

    group["variants"][variant] = evidence


def group_to_output(group):
    gene_id = group["gene_id"]
    chromosome = group["chromosome"]
    sample = group["sample"]
    tissue = group["tissue"]
    gene_symbols = list(group["gene_symbols"])
    if len(gene_symbols) != 1:
        raise ValueError(
            f"the same GeneID corresponding multiple GeneSymbol: GeneID={gene_id}; symbols={gene_symbols}"
        )

    unambiguous_evidence = [evidence for evidence in group["variants"].values() if evidence[2] <= 1]
    multi_evidence = [evidence for evidence in group["variants"].values() if evidence[2] > 1]
    unambiguous_hap1_sum = sum(evidence[0] for evidence in unambiguous_evidence)
    unambiguous_hap2_sum = sum(evidence[1] for evidence in unambiguous_evidence)
    multi_hap1_sum = sum(evidence[0] for evidence in multi_evidence)
    multi_hap2_sum = sum(evidence[1] for evidence in multi_evidence)
    hap1_sum = unambiguous_hap1_sum + multi_hap1_sum
    hap2_sum = unambiguous_hap2_sum + multi_hap2_sum

    action_counts = {
        output_column: sum(
            evidence[3] == action for evidence in group["variants"].values()
        )
        for action, output_column in CNV_ACTION_COUNT_COLUMNS.items()
    }

    return {
        "group_id": f"{gene_id}|{chromosome}|{sample}|{tissue}",
        "record_type": group["record_type"],
        "gene_id": gene_id,
        "gene_symbol": gene_symbols[0],
        **{field: group[field] for field in GROUP_METADATA_FIELDS},
        "chromosome": chromosome,
        "sample": sample,
        "tissue": tissue,
        "variant_count": len(group["variants"]),
        "unambiguous_variant_count": len(unambiguous_evidence),
        "multi_gene_variant_count": len(multi_evidence),
        **action_counts,
        "hap1_count_sum": hap1_sum,
        "hap2_count_sum": hap2_sum,
        "total_count_sum": hap1_sum + hap2_sum,
        "unambiguous_hap1_count_sum": unambiguous_hap1_sum,
        "unambiguous_hap2_count_sum": unambiguous_hap2_sum,
        "multi_gene_hap1_count_sum": multi_hap1_sum,
        "multi_gene_hap2_count_sum": multi_hap2_sum,
        "count_unit": COUNT_UNIT,
        "aggregation_method": AGGREGATION_METHOD,
    }


def write_tsv_atomic(path, fieldnames, rows):
    output_dir = os.path.dirname(path)
    os.makedirs(output_dir, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            dir=output_dir,
            prefix=f".{os.path.basename(path)}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = handle.name
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, path)
    except Exception:
        if temporary_path and os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise


def write_qc_atomic(path, metrics):
    rows = [{"metric": key, "value": value} for key, value in metrics.items()]
    write_tsv_atomic(path, ["metric", "value"], rows)


def aggregate_one_file(input_file, output_file, qc_file, animal_id, term):
    if not os.path.isfile(input_file) or os.path.getsize(input_file) == 0:
        raise FileNotFoundError(f"1.4 input missing or is empty: {input_file}")
    for formal_output in (output_file, qc_file):
        if os.path.exists(formal_output):
            raise FileExistsError(f"Refusing to overwrite an existing result; output already exists: {formal_output}")

    groups = OrderedDict()
    counters = OrderedDict(
        input_rows=0,
        cnv_clean_rows=0,
        cnv_status_ok_rows=0,
        cnv_status_not_evaluable_rows=0,
        cnv_evaluable_rows=0,
        cnv_not_evaluable_rows=0,
        phase_exact_rows=0,
        phased_het_gt_rows=0,
        skipped_phase_not_exact=0,
        skipped_invalid_phased_gt=0,
    )

    log(f"[stage] read and verify 1.4 input: {input_file}")
    with open(input_file, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"input file no column header: {input_file}")
        reader.fieldnames = [str(name).strip() for name in reader.fieldnames]
        require_columns(reader.fieldnames, input_file)

        for row_number, raw_row in enumerate(reader, start=2):
            counters["input_rows"] += 1
            row = {
                str(key).strip(): ("" if value is None else str(value).strip())
                for key, value in raw_row.items()
            }

            if clean_value(row["SAMPLE"]) != animal_id or clean_value(row["TISSUE"]) != term:
                raise ValueError(
                    f"row {row_number} sample information does not match the task: "
                    f"SAMPLE={row['SAMPLE']}; TISSUE={row['TISSUE']}; "
                    f"expected={animal_id}/{term}"
                )

            cnv_filter_action = clean_value(row["cnv_filter_action"])
            if cnv_filter_action not in ALLOWED_CNV_FILTER_ACTIONS:
                raise ValueError(
                    f"row {row_number} cnv_filter_action is not an allowed value for the Step 1.4 known_cnv_removed clean table: "
                    f"{cnv_filter_action!r}"
                )
            if not parse_bool(row["cnv_clean_pass"], "cnv_clean_pass", row_number):
                raise ValueError(
                    f"row {row_number} clean table contains cnv_clean_pass=FALSE: "
                    f"cnv_filter_action={cnv_filter_action}"
                )
            counters["cnv_clean_rows"] += 1
            action_counter = f"retained_cnv_action_{cnv_filter_action}_rows"
            counters[action_counter] = counters.get(action_counter, 0) + 1

            cnv_result_status = clean_value(row["cnv_result_status"])
            if cnv_result_status not in ALLOWED_CNV_STATUSES:
                raise ValueError(
                    f"row {row_number} cnv_result_status is not an allowed status: "
                    f"{cnv_result_status!r}"
                )
            cnv_evaluable = parse_bool(
                row["cnv_evaluable"], "cnv_evaluable", row_number
            )
            if cnv_result_status == "OK":
                if not cnv_evaluable or cnv_filter_action == "flag_not_evaluable":
                    raise ValueError(
                        f"row {row_number} OK status is inconsistent with cnv_evaluable/action"
                    )
                counters["cnv_status_ok_rows"] += 1
                counters["cnv_evaluable_rows"] += 1
            else:
                if cnv_evaluable or cnv_filter_action != "flag_not_evaluable":
                    raise ValueError(
                        f"row {row_number} non-evaluable status is inconsistent with cnv_evaluable/action"
                    )
                counters["cnv_status_not_evaluable_rows"] += 1
                counters["cnv_not_evaluable_rows"] += 1

            if clean_value(row["PHASE_STATUS"]) != EXPECTED_PHASE_STATUS:
                counters["skipped_phase_not_exact"] += 1
                continue
            counters["phase_exact_rows"] += 1

            if clean_value(row["VCF_GT"]) not in {"0|1", "1|0"}:
                counters["skipped_invalid_phased_gt"] += 1
                continue
            counters["phased_het_gt_rows"] += 1

            chromosome = clean_value(row["CHROM"])
            phase_set = clean_value(row["PHASE_SET"])
            if chromosome == "NA" or phase_set != chromosome:
                raise ValueError(
                    f"row {row_number} has inconsistent chromosome-wide phase fields: CHROM={chromosome}; PHASE_SET={phase_set}"
                )

            if clean_value(row["HAP1_ALLELE"]) == "NA" or clean_value(row["HAP2_ALLELE"]) == "NA":
                raise ValueError(f"row {row_number} strict phase record is missing HAP1_ALLELE or HAP2_ALLELE")

            record_type, gene_hit_n = validate_annotation_contract(row, row_number)
            gene_id = clean_value(row["GeneID"])
            group_key = (gene_id, chromosome, animal_id, term)
            if group_key not in groups:
                groups[group_key] = new_group(row, record_type)
            elif groups[group_key]["record_type"] != record_type:
                raise ValueError(f"Conflicting record_type values within the group for GeneID: {gene_id}")

            add_row_to_group(groups[group_key], row, row_number, gene_hit_n)

            if counters["input_rows"] % 50000 == 0:
                log(
                    f"[progress] read {counters['input_rows']} row, "
                    f"strict phase and CNV-clean={counters['phased_het_gt_rows']} row"
                )

    log(f"[stage] completed aggregation and construct output record: groups={len(groups)}")
    output_rows = [
        group_to_output(groups[key])
        for key in sorted(groups, key=lambda item: (item[1], item[0], item[2], item[3]))
    ]
    known_gene_groups = sum(row["record_type"] == "known_gene" for row in output_rows)
    manual_new_gene_groups = sum(
        row["record_type"] == "manual_new_gene" for row in output_rows
    )
    unannotated_groups = sum(row["record_type"] == "unannotated_variant" for row in output_rows)
    deduplicated_rows = sum(group["deduplicated_row_count"] for group in groups.values())

    metrics = OrderedDict(
        animal_id=animal_id,
        tissue=term,
        input_tsv=input_file,
        output_tsv=output_file,
        qc_tsv=qc_file,
        aggregation_key="GeneID+CHROM+SAMPLE+TISSUE",
        phase_status_required=EXPECTED_PHASE_STATUS,
        vcf_gt_required="0|1_or_1|0",
        cnv_result_status_allowed=",".join(sorted(ALLOWED_CNV_STATUSES)),
        cnv_clean_pass_required="TRUE",
        cnv_evaluable_contract="OK=TRUE;non_OK=FALSE_and_flag_not_evaluable",
        cnv_input_policy="known_cnv_removed:exclude_only_single_state_consensus_cnv_at_7_of_9_scales",
        retained_cnv_filter_actions="keep,flag_no_cnv_consensus,flag_soft_cnv_overlap,flag_strong_cnv_overlap,flag_cnv_conflict,flag_not_evaluable",
        phase_set_scope=EXPECTED_PHASE_SET_SCOPE,
        local_ps_location="retained_in_1.4_input_not_repeated_in_gene_id_table",
        unannotated_id_contract="GeneID/GeneID_base/GeneSymbol=(NA)<Variant> when GeneHitN=0",
        count_unit=COUNT_UNIT,
        aggregation_method=AGGREGATION_METHOD,
        **counters,
        output_groups=len(output_rows),
        known_gene_groups=known_gene_groups,
        manual_new_gene_groups=manual_new_gene_groups,
        unannotated_variant_groups=unannotated_groups,
        deduplicated_same_gene_variant_rows=deduplicated_rows,
        output_variant_references=sum(row["variant_count"] for row in output_rows),
        output_unambiguous_variant_references=sum(
            row["unambiguous_variant_count"] for row in output_rows
        ),
        output_multi_gene_variant_references=sum(
            row["multi_gene_variant_count"] for row in output_rows
        ),
        output_hap1_count_sum=sum(row["hap1_count_sum"] for row in output_rows),
        output_hap2_count_sum=sum(row["hap2_count_sum"] for row in output_rows),
        output_unambiguous_hap1_count_sum=sum(
            row["unambiguous_hap1_count_sum"] for row in output_rows
        ),
        output_unambiguous_hap2_count_sum=sum(
            row["unambiguous_hap2_count_sum"] for row in output_rows
        ),
        output_multi_gene_hap1_count_sum=sum(
            row["multi_gene_hap1_count_sum"] for row in output_rows
        ),
        output_multi_gene_hap2_count_sum=sum(
            row["multi_gene_hap2_count_sum"] for row in output_rows
        ),
    )

    log(f"[stage] atomic write aggregation result: {output_file}")
    write_tsv_atomic(output_file, OUTPUT_HEADER, output_rows)
    try:
        log(f"[stage] atomic write QC summary: {qc_file}")
        write_qc_atomic(qc_file, metrics)
    except Exception:
        if os.path.exists(output_file):
            os.unlink(output_file)
        raise
    return metrics


def main():
    args = parse_args()
    input_file, output_file, qc_file = build_paths(
        args.animal_id,
        args.term,
        args.input_tsv,
        args.out_dir,
    )

    log(f"[start] animal_id={args.animal_id} term={args.term}")
    metrics = aggregate_one_file(
        input_file,
        output_file,
        qc_file,
        args.animal_id,
        args.term,
    )
    log(
        f"[summary] input_rows={metrics['input_rows']} "
        f"phased_het_gt_rows={metrics['phased_het_gt_rows']} "
        f"output_groups={metrics['output_groups']}"
    )
    log("[done] 1.5 GeneID haplotype-count aggregation completed; task completed / done")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"[ERROR] run failed: {exc}")
        raise
