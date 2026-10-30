#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3


import argparse
import csv
import os
import re
import subprocess
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime


TASK_FILE = (
    '/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv'
)
LONG_ROOT = (
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.4ase-cnv-qc-filter/20260821_154914'
)
GENE_ROOT = (
    '/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.5ase-cnv-clean-gene-id-haplotype-table/20260821_161805'
)
PHASE_ROOT = (
    '/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.10_reconcile_dog11_whatsHap_shapeit5/reconciled_bcf'
)
BCFTOOLS_BIN = "/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"
CHROMOSOMES = tuple("chr{0}".format(index) for index in range(1, 39))
EXPECTED_PHASE_STATUS = "PHASED_HET_EXACT_MATCH"
EXPECTED_COUNT_UNIT = "ASEReadCounter_site_allele_support_sum"
EXPECTED_AGGREGATION_METHOD = "sum_distinct_variant_allele_support_within_gene_id"
LONG_COMPLETION_MARKER = os.path.join(LONG_ROOT, "execution_completion.marker")
GENE_COMPLETION_MARKER = os.path.join(GENE_ROOT, "execution_completion.marker")

LONG_REQUIRED = {
    "Variant", "GeneSymbol", "GeneID", "SAMPLE", "TISSUE", "CHROM", "POS",
    "REF", "ALT", "VCF_GT", "PS", "PHASE_SET", "PHASE_STATUS",
    "HAP1_ALLELE", "HAP1_COUNT", "HAP2_ALLELE", "HAP2_COUNT",
    "REF_COUNT", "ALT_COUNT", "TOTAL_COUNT", "cnv_filter_action", "cnv_clean_pass",
}
GENE_REQUIRED = {
    "group_id", "record_type", "gene_id", "gene_symbol", "chromosome", "sample",
    "tissue", "variant_count", "hap1_count_sum", "hap2_count_sum", "total_count_sum",
    "count_unit", "aggregation_method",
}


def log(message):
    stamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S ]")
    print("{0} {1}".format(stamp, message), flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check ASE haplotype-count orientation across all chromosomes in the five dogs and generate upstream QC tables for PLASMA."
    )
    parser.add_argument("--out-dir", required=True, help="Absolute output directory for this run.")
    args = parser.parse_args()
    if not os.path.isabs(args.out_dir):
        parser.error("--out-dir must be an absolute path")
    return args


def clean(value):
    return "" if value is None else str(value).strip()


def normalize_ps(value):
    text = clean(value)
    return "NA" if text in {"", ".", "NA"} else text


def parse_nonnegative_int(value, field, context):
    text = clean(value)
    try:
        parsed = int(text)
    except ValueError as exc:
        raise ValueError("{0}: {1} is not an integer: {2!r}".format(context, field, text)) from exc
    if parsed < 0:
        raise ValueError("{0}: {1} cannot be negative: {2}".format(context, field, parsed))
    return parsed


def require_columns(path, fieldnames, required):
    missing = sorted(required - set(fieldnames or []))
    if missing:
        raise ValueError("{0} missing fields: {1}".format(path, ",".join(missing)))


def safe_block_id(gene_id, chrom, sample, tissue):
    raw = "{0}__{1}__{2}__{3}".format(gene_id, chrom, sample, tissue)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw).strip("_")


def long_path(animal_id, tissue):
    return os.path.join(
        LONG_ROOT,
        animal_id,
        tissue,
        "{0}_{1}.ase_haplotype.cnv_clean.tsv".format(animal_id, tissue),
    )


def gene_path(animal_id, tissue):
    return os.path.join(
        GENE_ROOT,
        animal_id,
        tissue,
        "{0}_{1}.cnv_clean.gene_id_haplotype.tsv".format(animal_id, tissue),
    )


def bcf_path(chrom):
    return os.path.join(
        PHASE_ROOT,
        "{0}.dog11.ref1519.whatsHap_reconciled.phased.bcf".format(chrom),
    )


def read_tasks():
    if not os.path.isfile(TASK_FILE):
        raise FileNotFoundError("sample task list does not exist: {0}".format(TASK_FILE))
    tasks = []
    seen = set()
    with open(TASK_FILE, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(TASK_FILE, reader.fieldnames, {"animal_id", "term"})
        for line_number, row in enumerate(reader, start=2):
            animal_id = clean(row.get("animal_id"))
            tissue = clean(row.get("term"))
            if not animal_id or not tissue:
                raise ValueError("{0}: row {1} has an empty animal or tissue field".format(TASK_FILE, line_number))
            key = (animal_id, tissue)
            if key in seen:
                raise ValueError("sample task duplicate: {0}/{1}".format(animal_id, tissue))
            seen.add(key)
            tasks.append(key)
    if not tasks:
        raise ValueError("sample task list no data row: {0}".format(TASK_FILE))
    return tasks


def validate_inputs(tasks):
    if not os.path.isfile(BCFTOOLS_BIN) or not os.access(BCFTOOLS_BIN, os.X_OK):
        raise FileNotFoundError("bcftools not executable: {0}".format(BCFTOOLS_BIN))
    missing = []
    for animal_id, tissue in tasks:
        for path in (long_path(animal_id, tissue), gene_path(animal_id, tissue)):
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                missing.append(path)
    for chrom in CHROMOSOMES:
        path = bcf_path(chrom)
        for required in (path, path + ".csi"):
            if not os.path.isfile(required) or os.path.getsize(required) == 0:
                missing.append(required)
    if missing:
        raise FileNotFoundError("missing or is empty of input:\n" + "\n".join(missing))
    for marker in (LONG_COMPLETION_MARKER, GENE_COMPLETION_MARKER):
        if not os.path.isfile(marker):
            raise FileNotFoundError("upstream completion marker does not exist: {0}".format(marker))
        with open(marker) as handle:
            marker_lines = {line.strip() for line in handle if line.strip()}
        if "execution_status=success" not in marker_lines:
            raise ValueError("the upstream batch is not marked as successful: {0}".format(marker))


def validate_bcf_samples(tasks):
    first_bcf = bcf_path(CHROMOSOMES[0])
    result = subprocess.run(
        [BCFTOOLS_BIN, "query", "-l", first_bcf],
        check=True,
        text=True,
        capture_output=True,
    )
    bcf_samples = set(result.stdout.splitlines())
    animals = sorted({animal_id for animal_id, _ in tasks})
    missing = [animal_id for animal_id in animals if animal_id not in bcf_samples]
    if missing:
        raise ValueError("final BCF missing five dogs sample: {0}".format(",".join(missing)))
    return animals, len(bcf_samples)


def new_variant_row(row, source_path, line_number):
    context = "{0}: row {1}".format(source_path, line_number)
    sample = clean(row["SAMPLE"])
    tissue = clean(row["TISSUE"])
    chrom = clean(row["CHROM"])
    gene_id = clean(row["GeneID"])
    gene_symbol = clean(row["GeneSymbol"])
    pos = parse_nonnegative_int(row["POS"], "POS", context)
    return {
        "block_id": safe_block_id(gene_id, chrom, sample, tissue),
        "sample": sample,
        "tissue": tissue,
        "feature_id": gene_id,
        "feature_label": gene_symbol,
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "chrom": chrom,
        "pos": pos,
        "phase_set": clean(row["PHASE_SET"]),
        "source_ps": normalize_ps(row["PS"]),
        "variant_id": clean(row["Variant"]),
        "ase_ref": clean(row["REF"]),
        "ase_alt": clean(row["ALT"]),
        "ase_gt": clean(row["VCF_GT"]),
        "ase_hap1_allele": clean(row["HAP1_ALLELE"]),
        "ase_hap2_allele": clean(row["HAP2_ALLELE"]),
        "ase_hap1_count": clean(row["HAP1_COUNT"]),
        "ase_hap2_count": clean(row["HAP2_COUNT"]),
        "ase_total_count": clean(row["TOTAL_COUNT"]),
        "ref_count": clean(row["REF_COUNT"]),
        "alt_count": clean(row["ALT_COUNT"]),
        "phase_status": clean(row["PHASE_STATUS"]),
        "cnv_filter_action": clean(row["cnv_filter_action"]),
        "vcf_gt": "",
        "vcf_ps": "",
        "vcf_ref": "",
        "vcf_alt": "",
        "vcf_hapA_allele": "",
        "vcf_hapB_allele": "",
        "orientation_status": "unknown",
        "evidence_source": "phase_not_exact",
        "source_long_path": source_path,
        "source_line_number": line_number,
    }


def read_long_tables(tasks):
    rows = []
    target_keys = defaultdict(set)
    phase_statuses = Counter()
    for task_index, (animal_id, tissue) in enumerate(tasks, start=1):
        path = long_path(animal_id, tissue)
        log("[stage:1.4 {0}/{1}] {2}".format(task_index, len(tasks), path))
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(path, reader.fieldnames, LONG_REQUIRED)
            for line_number, row in enumerate(reader, start=2):
                if clean(row["SAMPLE"]) != animal_id or clean(row["TISSUE"]) != tissue:
                    raise ValueError("{0}: row {1} SAMPLE/TISSUE does not match the task path".format(path, line_number))
                if clean(row["cnv_clean_pass"]).upper() != "TRUE":
                    raise ValueError("{0}: row {1} clean table contains cnv_clean_pass!=TRUE".format(path, line_number))
                parsed = new_variant_row(row, path, line_number)
                rows.append(parsed)
                phase_statuses[parsed["phase_status"]] += 1
                if parsed["phase_status"] == EXPECTED_PHASE_STATUS:
                    target_keys[animal_id].add(
                        (parsed["chrom"], parsed["pos"], parsed["ase_ref"], parsed["ase_alt"])
                    )
    return rows, target_keys, phase_statuses


def load_bcf_records(animals, target_keys):
    records = {}
    duplicate_keys = set()
    positions_by_animal_chrom = defaultdict(lambda: defaultdict(set))
    target_chromosomes = set()
    for animal_id in animals:
        for chrom, pos, _, _ in target_keys[animal_id]:
            positions_by_animal_chrom[animal_id][chrom].add(pos)
            target_chromosomes.add(chrom)

    for chrom_index, chrom in enumerate(CHROMOSOMES, start=1):
        if chrom not in target_chromosomes:
            continue
        path = bcf_path(chrom)
        union_positions = set().union(
            *(positions_by_animal_chrom[animal_id].get(chrom, set()) for animal_id in animals)
        )
        log(
            "[stage:BCF chrom {0}/38] {1}; animals={2}; target_positions={3}".format(
                chrom_index, chrom, len(animals), len(union_positions)
            )
        )
        command = [
            BCFTOOLS_BIN, "query", "-s", ",".join(animals),
            "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%GT\t%PS]\n", path,
        ]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
        if process.stdout is None:
            raise RuntimeError("cannot read bcftools output: {0}".format(path))
        expected_field_count = 4 + 3 * len(animals)
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != expected_field_count:
                process.kill()
                raise ValueError("bcftools query row format unexpected: {0}".format(line.rstrip()))
            bcf_chrom, pos_text, ref, alts_text = fields[:4]
            pos = int(pos_text)
            if pos not in union_positions:
                continue
            for offset in range(4, expected_field_count, 3):
                animal_id, gt, ps = fields[offset:offset + 3]
                if pos not in positions_by_animal_chrom[animal_id].get(chrom, set()):
                    continue
                for alt_index, alt in enumerate(alts_text.split(","), start=1):
                    key = (animal_id, bcf_chrom, pos, ref, alt)
                    if key in records:
                        duplicate_keys.add(key)
                        continue
                    records[key] = {
                        "gt": clean(gt),
                        "ps": normalize_ps(ps),
                        "ref": ref,
                        "alt": alt,
                        "alts": alts_text.split(","),
                        "alt_index": alt_index,
                    }
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError("bcftools query failed: {0}; rc={1}".format(path, return_code))
    return records, duplicate_keys


def allele_from_code(code, ref, alts):
    if code == "0":
        return ref
    try:
        index = int(code)
    except ValueError:
        return ""
    if 1 <= index <= len(alts):
        return alts[index - 1]
    return ""


def classify_exact_variant(row, bcf_records, duplicate_keys):
    context = "{0}: row {1}".format(row["source_long_path"], row["source_line_number"])
    if row["phase_set"] != row["chrom"]:
        row["evidence_source"] = "phase_set_not_chromosome"
        row["orientation_status"] = "conflict"
        return
    key = (row["sample"], row["chrom"], row["pos"], row["ase_ref"], row["ase_alt"])
    if key in duplicate_keys:
        row["evidence_source"] = "duplicate_exact_bcf_variant"
        row["orientation_status"] = "conflict"
        return
    record = bcf_records.get(key)
    if record is None:
        row["evidence_source"] = "exact_variant_missing_from_final_bcf"
        row["orientation_status"] = "unknown"
        return
    row.update(
        {
            "vcf_gt": record["gt"],
            "vcf_ps": record["ps"],
            "vcf_ref": record["ref"],
            "vcf_alt": record["alt"],
        }
    )
    gt_parts = record["gt"].split("|")
    if len(gt_parts) != 2 or gt_parts[0] == gt_parts[1]:
        row["evidence_source"] = "final_bcf_not_phased_heterozygote"
        row["orientation_status"] = "conflict"
        return
    hap_a = allele_from_code(gt_parts[0], record["ref"], record["alts"])
    hap_b = allele_from_code(gt_parts[1], record["ref"], record["alts"])
    row["vcf_hapA_allele"] = hap_a
    row["vcf_hapB_allele"] = hap_b
    ref_count = parse_nonnegative_int(row["ref_count"], "REF_COUNT", context)
    alt_count = parse_nonnegative_int(row["alt_count"], "ALT_COUNT", context)
    hap1_count = parse_nonnegative_int(row["ase_hap1_count"], "HAP1_COUNT", context)
    hap2_count = parse_nonnegative_int(row["ase_hap2_count"], "HAP2_COUNT", context)
    total_count = parse_nonnegative_int(row["ase_total_count"], "TOTAL_COUNT", context)
    count_by_allele = {record["ref"]: ref_count, record["alt"]: alt_count}
    checks = [
        row["ase_gt"] == record["gt"],
        row["source_ps"] == record["ps"],
        row["ase_hap1_allele"] == hap_a,
        row["ase_hap2_allele"] == hap_b,
        hap_a in count_by_allele,
        hap_b in count_by_allele,
        hap1_count == count_by_allele.get(hap_a),
        hap2_count == count_by_allele.get(hap_b),
        hap1_count + hap2_count == total_count,
    ]
    if all(checks):
        row["orientation_status"] = "consistent"
        row["evidence_source"] = "final_reconciled_bcf_exact_match"
    else:
        row["orientation_status"] = "conflict"
        row["evidence_source"] = "final_bcf_or_count_contract_mismatch"


def classify_variants(rows, bcf_records, duplicate_keys):
    for index, row in enumerate(rows, start=1):
        if row["phase_status"] == EXPECTED_PHASE_STATUS:
            classify_exact_variant(row, bcf_records, duplicate_keys)
        if index % 50000 == 0:
            log("[stage:variant_qc] already verify {0} row".format(index))


def read_gene_tables(tasks):
    gene_rows = []
    seen_keys = set()
    for task_index, (animal_id, tissue) in enumerate(tasks, start=1):
        path = gene_path(animal_id, tissue)
        log("[stage:1.5 {0}/{1}] {2}".format(task_index, len(tasks), path))
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(path, reader.fieldnames, GENE_REQUIRED)
            for line_number, row in enumerate(reader, start=2):
                if clean(row["sample"]) != animal_id or clean(row["tissue"]) != tissue:
                    raise ValueError("{0}: row {1} sample/tissue does not match the task path".format(path, line_number))
                if clean(row["count_unit"]) != EXPECTED_COUNT_UNIT:
                    raise ValueError("{0}: row {1} count_unit does not satisfy 1.5 contract".format(path, line_number))
                if clean(row["aggregation_method"]) != EXPECTED_AGGREGATION_METHOD:
                    raise ValueError("{0}: row {1} aggregation_method does not satisfy 1.5 contract".format(path, line_number))
                key = (clean(row["gene_id"]), clean(row["chromosome"]), animal_id, tissue)
                if key in seen_keys:
                    raise ValueError("1.5 encountered duplicate aggregation key: {0}".format(key))
                seen_keys.add(key)
                copied = dict(row)
                copied["source_gene_path"] = path
                copied["source_line_number"] = line_number
                copied["group_key"] = key
                gene_rows.append(copied)
    return gene_rows


def group_exact_variants(rows):
    grouped = defaultdict(OrderedDict)
    for row in rows:
        if row["phase_status"] != EXPECTED_PHASE_STATUS:
            continue
        key = (row["gene_id"], row["chrom"], row["sample"], row["tissue"])
        variant_id = row["variant_id"]
        evidence = (
            row["ase_hap1_count"], row["ase_hap2_count"], row["ase_total_count"],
            row["orientation_status"], row["evidence_source"], row["pos"], row["source_ps"],
        )
        previous = grouped[key].get(variant_id)
        if previous is not None and previous != evidence:
            raise ValueError("the same GeneID+Variant encountered conflict evidence: {0}; {1}".format(key, variant_id))
        grouped[key][variant_id] = evidence
    return grouped


def build_block_rows(gene_rows, grouped):
    blocks = []
    gene_keys = {row["group_key"] for row in gene_rows}
    extra_group_keys = sorted(set(grouped) - gene_keys)
    if extra_group_keys:
        raise ValueError(
            "1.4 the strict-phase aggregation key is absent from 1.5 in, before 5 items: {0}".format(extra_group_keys[:5])
        )
    for row in gene_rows:
        key = row["group_key"]
        variants = grouped.get(key, OrderedDict())
        statuses = Counter(value[3] for value in variants.values())
        sources = Counter(value[4] for value in variants.values())
        hap1_sum = sum(int(value[0]) for value in variants.values())
        hap2_sum = sum(int(value[1]) for value in variants.values())
        positions = [value[5] for value in variants.values()]
        source_ps_values = sorted({value[6] for value in variants.values() if value[6]})
        expected_variant_count = int(row["variant_count"])
        expected_hap1 = int(row["hap1_count_sum"])
        expected_hap2 = int(row["hap2_count_sum"])
        expected_total = int(row["total_count_sum"])
        aggregate_ok = (
            len(variants) == expected_variant_count
            and hap1_sum == expected_hap1
            and hap2_sum == expected_hap2
            and hap1_sum + hap2_sum == expected_total
        )
        checked = statuses.get("consistent", 0) + statuses.get("conflict", 0)
        if aggregate_ok and checked == expected_variant_count and statuses.get("conflict", 0) == 0:
            orientation = "consistent"
        elif not variants or checked == 0:
            orientation = "unknown"
        else:
            orientation = "conflict"
        source_summary = ";".join(
            "{0}:{1}".format(name, sources[name]) for name in sorted(sources)
        )
        blocks.append(
            {
                "block_id": safe_block_id(row["gene_id"], row["chromosome"], row["sample"], row["tissue"]),
                "sample": row["sample"],
                "tissue": row["tissue"],
                "feature_id": row["gene_id"],
                "feature_label": row["gene_symbol"],
                "gene_id": row["gene_id"],
                "gene_symbol": row["gene_symbol"],
                "record_type": row["record_type"],
                "chrom": row["chromosome"],
                "phase_set": row["chromosome"],
                "source_ps_values": "|".join(source_ps_values),
                "variant_count": expected_variant_count,
                "checked_variant_count": checked,
                "consistent_count": statuses.get("consistent", 0),
                "swapped_count": 0,
                "conflict_count": statuses.get("conflict", 0),
                "unknown_count": max(0, expected_variant_count - checked),
                "orientation_status": orientation,
                "orientation_evidence_summary": source_summary,
                "aggregation_status": "matched_1.5" if aggregate_ok else "mismatch_1.5",
                "hap1_count": expected_hap1,
                "hap2_count": expected_hap2,
                "total_count": expected_total,
                "min_pos": min(positions) if positions else "",
                "max_pos": max(positions) if positions else "",
                "source_gene_path": row["source_gene_path"],
            }
        )
    return blocks


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    log("[start] five dogs all chromosome HAP1/HAP2 orientation QC")
    log("[config] TASK_FILE={0}".format(TASK_FILE))
    log("[config] LONG_ROOT={0}".format(LONG_ROOT))
    log("[config] GENE_ROOT={0}".format(GENE_ROOT))
    log("[config] PHASE_ROOT={0}".format(PHASE_ROOT))
    log("[config] BCFTOOLS_BIN={0}".format(BCFTOOLS_BIN))
    log("[config] OUT_DIR={0}".format(out_dir))

    log("[stage:input] read sample list and verify all absolute path")
    tasks = read_tasks()
    validate_inputs(tasks)
    animals, bcf_sample_count = validate_bcf_samples(tasks)
    log("[summary:input] tasks={0}; animals={1}; BCF_samples={2}".format(len(tasks), len(animals), bcf_sample_count))

    log("[stage:1.4] read CNV-clean site-level ASE evidence")
    variant_rows, target_keys, phase_statuses = read_long_tables(tasks)
    log("[summary:1.4] rows={0}; phase_status={1}".format(len(variant_rows), dict(phase_statuses)))

    log("[stage:BCF] cross-check final 38 chromosome multiple sample phase BCF")
    bcf_records, duplicate_keys = load_bcf_records(animals, target_keys)
    log("[summary:BCF] matched_exact_records={0}; duplicate_keys={1}".format(len(bcf_records), len(duplicate_keys)))

    log("[stage:variant_qc] verify GT, local PS, allele orientation and count")
    classify_variants(variant_rows, bcf_records, duplicate_keys)

    log("[stage:1.5] read GeneID+chromosome aggregation table and check site summary")
    gene_rows = read_gene_tables(tasks)
    grouped = group_exact_variants(variant_rows)
    block_rows = build_block_rows(gene_rows, grouped)

    variant_out = os.path.join(out_dir, "orientation_variant_qc.tsv")
    block_out = os.path.join(out_dir, "orientation_block_qc.tsv")
    existing_outputs = [path for path in (variant_out, block_out) if os.path.exists(path)]
    if existing_outputs:
        raise FileExistsError("output already exists, refuse to overwrite: {0}".format(",".join(existing_outputs)))
    os.makedirs(out_dir, exist_ok=True)

    variant_fields = [
        "block_id", "sample", "tissue", "feature_id", "feature_label", "gene_id", "gene_symbol",
        "chrom", "pos", "phase_set", "source_ps", "variant_id", "phase_status",
        "ase_ref", "ase_alt", "ase_gt", "ase_hap1_allele", "ase_hap2_allele",
        "ase_hap1_count", "ase_hap2_count", "ase_total_count", "ref_count", "alt_count",
        "vcf_gt", "vcf_ps", "vcf_ref", "vcf_alt", "vcf_hapA_allele", "vcf_hapB_allele",
        "cnv_filter_action", "orientation_status", "evidence_source", "source_long_path",
    ]
    block_fields = [
        "block_id", "sample", "tissue", "feature_id", "feature_label", "gene_id", "gene_symbol",
        "record_type", "chrom", "phase_set", "source_ps_values", "variant_count",
        "checked_variant_count", "consistent_count", "swapped_count", "conflict_count", "unknown_count",
        "orientation_status", "orientation_evidence_summary", "aggregation_status",
        "hap1_count", "hap2_count", "total_count", "min_pos", "max_pos", "source_gene_path",
    ]
    log("[stage:write] write per site and per GeneID-chromosome QC")
    write_tsv(variant_out, variant_rows, variant_fields)
    write_tsv(block_out, block_rows, block_fields)

    variant_status = Counter(row["orientation_status"] for row in variant_rows)
    block_status = Counter(row["orientation_status"] for row in block_rows)
    aggregation_status = Counter(row["aggregation_status"] for row in block_rows)
    log("[summary] variant_orientation={0}".format(dict(variant_status)))
    log("[summary] block_orientation={0}".format(dict(block_status)))
    log("[summary] aggregation_status={0}".format(dict(aggregation_status)))
    log("[output] {0}".format(variant_out))
    log("[output] {0}".format(block_out))
    log("[done] five dogs all chromosome HAP1/HAP2 orientation QC task completed / done")


if __name__ == "__main__":
    main()
