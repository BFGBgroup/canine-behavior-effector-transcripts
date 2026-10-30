#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3


import argparse
import bisect
import csv
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


STEP7_ROOT = Path("/path/to/project/analysis_20260609/5-dogs/4.cnv-rna/8.casper_dog_keep_locus_genes-v2")

ABNORMAL_STATES = {"amp", "del", "cnloh"}
STEP7_SCALE_PAIRS = tuple(f"cs{cnv}_ls{loh}" for cnv in range(1, 4) for loh in range(1, 4))
STEP7_SCALE_SET = set(STEP7_SCALE_PAIRS)
SCALE_ORDER = {scale: idx for idx, scale in enumerate(STEP7_SCALE_PAIRS)}
VARIANT_RE = re.compile(
    r"^(?:chr)?([0-9]+|X|Y|M|MT):([0-9]+)_([^>]+)>(.+)$",
    re.IGNORECASE,
)
UNANNOTATED_PREFIX = "(NA)"
HARD_FILTER_ACTIONS = {
    "exclude_cnv_overlap",
}

CNV_COLUMNS = [
    "cnv_sample_id",
    "cnv_result_status",
    "cnv_evaluable",
    "cnv_variant_coordinate_status",
    "cnv_variant_chr_normalized",
    "cnv_variant_pos",
    "cnv_step7_overlap_any",
    "cnv_step7_overlap_states",
    "cnv_step7_overlap_event_scales",
    "cnv_step7_overlap_segments",
    "cnv_step7_overlap_n_segments",
    "cnv_step7_available_scales",
    "cnv_step7_consensus_required_support",
    "cnv_step7_consensus_state",
    "cnv_step7_consensus_support",
    "cnv_step7_consensus_fraction",
    "cnv_step7_consensus_scales",
    "cnv_step7_state_support_counts",
    "cnv_step7_best_state",
    "cnv_step7_best_support",
    "cnv_step7_conflict_states",
    "cnv_step7_baf_decisions",
    "cnv_evidence_source",
    "cnv_filter_action",
    "cnv_filter_reason",
    "cnv_clean_pass",
]


def log(message):
    print(f"[{datetime.now():%F %T} ] {message}", flush=True)


def normalize_chr(value):
    text = str(value or "").strip()
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)
    match = re.match(r"^([0-9]+|X|Y|M|MT)", text, flags=re.IGNORECASE)
    if not match:
        return ""
    chrom = match.group(1).upper()
    return "M" if chrom == "MT" else chrom


def parse_int(value):
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def validate_annotation_contract(row, input_row_number, ase_tsv):
    variant = str(row.get("Variant", "") or "").strip()
    gene_hit_n = parse_int(row.get("GeneHitN"))
    if gene_hit_n is None or gene_hit_n < 0:
        raise ValueError(
            f"1.3 GeneHitN is invalid at input row {input_row_number}: "
            f"GeneHitN={row.get('GeneHitN', '')}; Variant={variant}; input={ase_tsv}"
        )

    id_fields = ("GeneSymbol", "GeneID", "GeneID_base", "AllGeneIDs")
    if gene_hit_n == 0:
        expected_id = f"{UNANNOTATED_PREFIX}{variant}"
        mismatches = [
            f"{field}={str(row.get(field, '') or '').strip()!r}"
            for field in id_fields
            if str(row.get(field, "") or "").strip() != expected_id
        ]
        if mismatches:
            raise ValueError(
                f"1.3 unannotated pseudo-ID contract mismatch at input row "
                f"{input_row_number}: expected={expected_id!r}; "
                f"observed={','.join(mismatches)}; input={ase_tsv}"
            )
        if (
            str(row.get("GeneFeature", "") or "").strip() != "intergenic"
            or str(row.get("GeneStrand", "") or "").strip() != "NA"
            or str(row.get("GeneOverlapType", "") or "").strip() != "INTERGENIC"
            or str(row.get("NCBI_GeneID", "") or "").strip() != "NA"
            or str(row.get("GeneFeatureType", "") or "").strip() != "unannotated"
            or str(row.get("GeneBiotype", "") or "").strip() != "unannotated"
            or parse_int(row.get("SourceGTFGeneIDN")) != 0
            or str(row.get("SourceGTFGeneIDs", "") or "").strip() != "NA"
        ):
            raise ValueError(
                f"1.3 unannotated annotation contract mismatch at input row "
                f"{input_row_number}: GeneFeature={row.get('GeneFeature', '')}; "
                f"GeneStrand={row.get('GeneStrand', '')}; "
                f"GeneOverlapType={row.get('GeneOverlapType', '')}; "
                f"NCBI_GeneID={row.get('NCBI_GeneID', '')}; "
                f"GeneFeatureType={row.get('GeneFeatureType', '')}; "
                f"GeneBiotype={row.get('GeneBiotype', '')}; "
                f"SourceGTFGeneIDN={row.get('SourceGTFGeneIDN', '')}; "
                f"SourceGTFGeneIDs={row.get('SourceGTFGeneIDs', '')}; "
                f"Variant={variant}; input={ase_tsv}"
            )
        return "unannotated"

    pseudo_fields = [
        field
        for field in id_fields
        if str(row.get(field, "") or "").strip().startswith(UNANNOTATED_PREFIX)
    ]
    if pseudo_fields:
        raise ValueError(
            f"1.3 annotated row contains unannotated pseudo-ID at input row "
            f"{input_row_number}: fields={','.join(pseudo_fields)}; "
            f"Variant={variant}; input={ase_tsv}"
        )

    ncbi_gene_id = str(row.get("NCBI_GeneID", "") or "").strip()
    feature_type = str(row.get("GeneFeatureType", "") or "").strip()
    gene_biotype = str(row.get("GeneBiotype", "") or "").strip()
    annotation_source = str(row.get("AnnotationSource", "") or "").strip()
    source_gene_id_n = parse_int(row.get("SourceGTFGeneIDN"))
    source_gene_ids = str(row.get("SourceGTFGeneIDs", "") or "").strip()

    if feature_type == "manual_new_gene":
        gene_id = str(row.get("GeneID", "") or "").strip()
        gene_id_base = str(row.get("GeneID_base", "") or "").strip()
        gene_symbol = str(row.get("GeneSymbol", "") or "").strip()
        if (
            not gene_id.startswith("DOG5_NEWGENE_")
            or gene_id_base != gene_id
            or gene_symbol != gene_id
            or ncbi_gene_id != "NA"
            or gene_biotype != "manual_new_gene"
            or annotation_source != "ManualNewGene"
            or source_gene_id_n is None
            or source_gene_id_n < 1
            or not source_gene_ids
            or source_gene_ids == "NA"
        ):
            raise ValueError(
                f"1.3 manual new-gene metadata mismatch at input row "
                f"{input_row_number}: GeneID={gene_id!r}; "
                f"GeneID_base={gene_id_base!r}; GeneSymbol={gene_symbol!r}; "
                f"NCBI_GeneID={ncbi_gene_id!r}; GeneBiotype={gene_biotype!r}; "
                f"AnnotationSource={annotation_source!r}; "
                f"SourceGTFGeneIDN={row.get('SourceGTFGeneIDN', '')!r}; "
                f"SourceGTFGeneIDs={source_gene_ids!r}; Variant={variant}; input={ase_tsv}"
            )
        return "manual_new_gene"

    if (
        not ncbi_gene_id
        or ncbi_gene_id == "NA"
        or feature_type not in {"gene", "pseudogene"}
        or not gene_biotype
        or gene_biotype == "unannotated"
        or source_gene_id_n is None
        or source_gene_id_n < 1
        or not source_gene_ids
        or source_gene_ids == "NA"
    ):
        raise ValueError(
            f"1.3 annotated known-gene metadata mismatch at input row "
            f"{input_row_number}: NCBI_GeneID={ncbi_gene_id!r}; "
            f"GeneFeatureType={feature_type!r}; GeneBiotype={gene_biotype!r}; "
            f"SourceGTFGeneIDN={row.get('SourceGTFGeneIDN', '')!r}; "
            f"SourceGTFGeneIDs={source_gene_ids!r}; Variant={variant}; input={ase_tsv}"
        )
    return "known_gene"


def sort_scales(scales):
    return sorted(
        (str(scale) for scale in scales if str(scale or "").strip()),
        key=lambda item: (SCALE_ORDER.get(item, 999), item),
    )


def join_unique(values):
    out = []
    seen = set()
    for value in values:
        text = str(value or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return "|".join(out)


def validate_variant_identity(row, input_row_number, ase_tsv):
    variant = str(row.get("Variant", "") or "").strip()
    chrom = normalize_chr(row.get("CHROM", ""))
    pos = parse_int(row.get("POS"))
    ref = str(row.get("REF", "") or "").strip().upper()
    alt = str(row.get("ALT", "") or "").strip().upper()
    match = VARIANT_RE.fullmatch(variant)
    if not match or not chrom or pos is None or not ref or not alt:
        raise ValueError(
            f"1.3 variant identity is incomplete at input row {input_row_number}: "
            f"Variant={variant!r}; CHROM={row.get('CHROM', '')!r}; "
            f"POS={row.get('POS', '')!r}; REF={ref!r}; ALT={alt!r}; input={ase_tsv}"
        )

    variant_chrom = normalize_chr(match.group(1))
    variant_pos = parse_int(match.group(2))
    variant_ref = match.group(3).strip().upper()
    variant_alt = match.group(4).strip().upper()
    if (variant_chrom, variant_pos, variant_ref, variant_alt) != (chrom, pos, ref, alt):
        raise ValueError(
            f"1.3 Variant disagrees with CHROM/POS/REF/ALT at input row "
            f"{input_row_number}: Variant={variant!r}; "
            f"CHROM/POS/REF/ALT={chrom}:{pos}_{ref}>{alt}; input={ase_tsv}"
        )
    return chrom, pos, "OK"


def load_step7_segments(step7_root, tissue, sample_id):
    segment_tsv = step7_root / tissue / "segments.tsv"
    stats = Counter()
    sample_scales = set()
    segments_by_chr = defaultdict(list)

    if not segment_tsv.is_file() or segment_tsv.stat().st_size == 0:
        stats["missing_step7_segments_tsv"] += 1
        return {}, sample_scales, stats

    with segment_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "chr", "start", "end", "finalCNVState", "scale_pair"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Step7 segments missing columns {','.join(missing)}: {segment_tsv}")

        for row in reader:
            stats["segment_rows_total"] += 1
            if str(row.get("ID", "")).strip() != sample_id:
                continue

            stats["segment_rows_for_sample"] += 1
            scale = str(row.get("scale_pair", "")).strip()
            if scale:
                sample_scales.add(scale)

            state = str(row.get("finalCNVState", "")).strip().lower()
            if state not in ABNORMAL_STATES:
                continue

            chrom = normalize_chr(row.get("chr", ""))
            start = parse_int(row.get("start"))
            end = parse_int(row.get("end"))
            if not chrom or start is None or end is None:
                stats["segment_rows_bad_coordinate"] += 1
                continue
            if end < start:
                start, end = end, start

            segment = {
                "start": start,
                "end": end,
                "state": state,
                "scale": scale,
                "event_scale": str(row.get("event_scale", "")).strip(),
                "baf_decision": str(row.get("baf_decision", "")).strip(),
                "label": (
                    f"{scale}:{row.get('chr', '')}:{start}-{end}:{state}:"
                    f"{row.get('baf_decision', '')}:{row.get('event_scale', '')}"
                ),
            }
            segments_by_chr[chrom].append(segment)
            stats["abnormal_segments_loaded"] += 1

    indexed = {}
    for chrom, segments in segments_by_chr.items():
        segments.sort(key=lambda item: (item["start"], item["end"]))
        indexed[chrom] = {
            "segments": segments,
            "starts": [item["start"] for item in segments],
        }

    return indexed, sample_scales, stats


def find_step7_overlaps(chrom, pos, segments_by_chr):
    if not chrom or pos is None:
        return []
    bucket = segments_by_chr.get(chrom)
    if not bucket:
        return []

    starts = bucket["starts"]
    segments = bucket["segments"]
    idx = bisect.bisect_right(starts, pos)
    overlaps = []
    seen = set()
    for segment in segments[:idx]:
        if segment["end"] < pos:
            continue
        key = segment["label"]
        if key in seen:
            continue
        seen.add(key)
        overlaps.append(segment)
    return overlaps


def format_state_support(state_to_scales):
    parts = []
    for state in sorted(state_to_scales):
        scales = sort_scales(state_to_scales[state])
        parts.append(f"{state}:{len(scales)}[{','.join(scales)}]")
    return ";".join(parts)


def summarize_step7(overlaps, sample_scales, consensus_fraction):
    required_support = int(math.ceil(consensus_fraction * len(STEP7_SCALE_PAIRS)))
    state_to_scales = defaultdict(set)
    for segment in overlaps:
        scale = segment.get("scale", "")
        state = segment.get("state", "")
        if scale in STEP7_SCALE_SET and state in ABNORMAL_STATES:
            state_to_scales[state].add(scale)

    best_state = ""
    best_support = 0
    if state_to_scales:
        best_state, best_scales = sorted(
            state_to_scales.items(),
            key=lambda item: (-len(item[1]), item[0]),
        )[0]
        best_support = len(best_scales)

    consensus_states = sorted(
        state for state, scales in state_to_scales.items()
        if len(scales) >= required_support
    )
    is_conflict = len(consensus_states) > 1
    consensus_state = "" if is_conflict or not consensus_states else consensus_states[0]
    consensus_scales = state_to_scales.get(consensus_state, set()) if consensus_state else set()
    consensus_support = len(consensus_scales)
    consensus_fraction_observed = ""
    if consensus_support > 0:
        consensus_fraction_observed = f"{consensus_support / len(STEP7_SCALE_PAIRS):.6g}"

    return {
        "available_scales": "|".join(sort_scales(sample_scales)),
        "required_support": required_support,
        "consensus_state": consensus_state,
        "consensus_support": consensus_support,
        "consensus_fraction": consensus_fraction_observed,
        "consensus_scales": "|".join(sort_scales(consensus_scales)),
        "state_support_counts": format_state_support(state_to_scales),
        "best_state": best_state,
        "best_support": best_support,
        "conflict_states": "|".join(consensus_states) if is_conflict else "",
        "is_conflict": is_conflict,
    }


def decide_action(evaluable, step7_overlaps, step7_summary):
    if not evaluable:
        return "flag_not_evaluable", "Step7 segment result is not available for this sample"

    if step7_overlaps and step7_summary["consensus_state"]:
        return (
            "exclude_cnv_overlap",
            f"Step7 finalCNVState={step7_summary['consensus_state']} reaches "
            f"{step7_summary['consensus_support']}/{len(STEP7_SCALE_PAIRS)} "
            f"multi-scale consensus",
        )

    if step7_overlaps:
        if step7_summary["is_conflict"]:
            return "flag_cnv_conflict", "Step7 has conflicting consensus CNV states"
        if step7_summary["best_support"] >= 5:
            return "flag_strong_cnv_overlap", "Step7 finalCNVState support is 5-6/9"
        if step7_summary["best_support"] >= 3:
            return "flag_soft_cnv_overlap", "Step7 finalCNVState support is 3-4/9"
        return "flag_no_cnv_consensus", "Step7 abnormal segment overlap has weak multi-scale support"

    return "keep", "No Step7 consensus CNV overlap"


def evidence_source(step7_overlaps, step7_summary):
    sources = []
    if step7_summary["consensus_state"]:
        sources.append("Step7_segment_consensus")
    elif step7_overlaps:
        sources.append("Step7_segment_overlap")
    return "|".join(sources)


def clean_keep(row, clean_mode):
    if clean_mode == "strict":
        return row.get("cnv_filter_action") == "keep"
    return row.get("cnv_filter_action") not in HARD_FILTER_ACTIONS


def annotate_ase(args):
    sample_id = f"{args.animal_id}_{args.tissue}"
    ase_tsv = args.ase_tsv
    out_dir = args.out_dir
    full_tsv = out_dir / f"{sample_id}.ase_haplotype.cnv_annotated.tsv"
    clean_tsv = out_dir / f"{sample_id}.ase_haplotype.cnv_clean.tsv"
    summary_tsv = out_dir / f"{sample_id}.ase_haplotype.cnv_qc_summary.tsv"
    full_tmp = full_tsv.with_name(f"{full_tsv.name}.tmp")
    clean_tmp = clean_tsv.with_name(f"{clean_tsv.name}.tmp")
    summary_tmp = summary_tsv.with_name(f"{summary_tsv.name}.tmp")

    log(f"[START] 1.4 ASE variant-gene CNV QC: sample_id={sample_id}")
    log(f"[INFO] ase_tsv={ase_tsv}")
    log(f"[INFO] step7_root={args.step7_root}")
    log(f"[INFO] out_dir={out_dir}")
    log(f"[INFO] consensus_fraction={args.consensus_fraction} clean_mode={args.clean_mode}")
    log("[INFO] phase_contract=PHASE_SET matches CHROM; PS is local WhatsHap provenance")

    if not ase_tsv.is_file():
        raise FileNotFoundError(f"ASE input not found: {ase_tsv}")
    existing_outputs = [
        path for path in (full_tsv, clean_tsv, summary_tsv) if path.exists()
    ]
    if existing_outputs:
        raise FileExistsError(
            "refuse to overwrite existing 1.4 main result: "
            + ", ".join(str(path) for path in existing_outputs)
        )

    out_dir.mkdir(parents=True, exist_ok=True)

    log("[STAGE] loading Step7 segment evidence")
    segments_by_chr, sample_scales, step7_stats = load_step7_segments(args.step7_root, args.tissue, sample_id)
    missing_scales = STEP7_SCALE_SET - sample_scales
    unexpected_scales = sample_scales - STEP7_SCALE_SET
    if step7_stats.get("missing_step7_segments_tsv"):
        result_status = "NO_STEP7_SEGMENT_RESULT"
    elif not step7_stats.get("segment_rows_for_sample", 0):
        result_status = "NO_STEP7_SAMPLE_RESULT"
    elif missing_scales or unexpected_scales:
        result_status = "INCOMPLETE_STEP7_SCALE_SET"
    else:
        result_status = "OK"
    evaluable = result_status == "OK"
    log(f"[INFO] step7_segment_rows_for_sample={step7_stats.get('segment_rows_for_sample', 0)}")
    log(f"[INFO] step7_abnormal_segments_loaded={step7_stats.get('abnormal_segments_loaded', 0)}")
    log(f"[INFO] step7_scales={','.join(sort_scales(sample_scales))}")
    log(f"[INFO] step7_missing_scales={','.join(sort_scales(missing_scales)) or 'NONE'}")
    log(f"[INFO] cnv_result_status={result_status}")

    action_row_counts = Counter()
    action_variant_counts = Counter()
    input_rows = 0
    clean_rows = 0
    seen_variants = set()
    clean_unique_variant_n = 0
    phase_contract_stats = Counter()
    annotation_contract_stats = Counter()
    previous_variant_key = None
    previous_cnv_fields = None
    previous_clean_pass = False
    previous_action = ""

    log("[STAGE] streaming 1.3 full ASE variant-gene rows")
    try:
        with ase_tsv.open("r", encoding="utf-8", newline="") as in_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"ASE input is empty: {ase_tsv}")

            required_columns = {
                "Variant", "GeneSymbol", "GeneID", "GeneID_base", "GeneFeature",
                "GeneStrand", "GeneHitN", "GeneOverlapType", "AllGeneIDs",
                "NCBI_GeneID", "GeneFeatureType", "GeneBiotype",
                "AnnotationSource", "SourceGTFGeneIDN", "SourceGTFGeneIDs",
                "SAMPLE", "TISSUE", "CHROM", "POS", "REF", "ALT",
                "PS", "PHASE_SET", "PHASE_STATUS",
            }
            missing_columns = sorted(required_columns - set(reader.fieldnames))
            if missing_columns:
                raise ValueError(
                    f"1.3 ASE input missing required columns "
                    f"{','.join(missing_columns)}: {ase_tsv}"
                )

            out_fields = list(reader.fieldnames)
            for column in CNV_COLUMNS:
                if column not in out_fields:
                    out_fields.append(column)

            with full_tmp.open("w", encoding="utf-8", newline="") as full_handle, clean_tmp.open("w", encoding="utf-8", newline="") as clean_handle:
                full_writer = csv.DictWriter(full_handle, delimiter="\t", fieldnames=out_fields, lineterminator="\n")
                clean_writer = csv.DictWriter(clean_handle, delimiter="\t", fieldnames=out_fields, lineterminator="\n")
                full_writer.writeheader()
                clean_writer.writeheader()

                for row in reader:
                    input_rows += 1
                    annotation_class = validate_annotation_contract(
                        row, input_rows, ase_tsv
                    )
                    annotation_contract_stats[f"{annotation_class}_rows"] += 1
                    row_animal = str(row.get("SAMPLE", "") or "").strip()
                    row_tissue = str(row.get("TISSUE", "") or "").strip()
                    if row_animal != args.animal_id or row_tissue != args.tissue:
                        raise ValueError(
                            f"1.3 SAMPLE/TISSUE mismatch at input row {input_rows}: "
                            f"observed={row_animal!r}/{row_tissue!r}; "
                            f"expected={args.animal_id!r}/{args.tissue!r}; input={ase_tsv}"
                        )
                    chrom, pos, variant_status = validate_variant_identity(
                        row, input_rows, ase_tsv
                    )
                    phase_status = str(row.get("PHASE_STATUS", "")).strip()
                    phase_set = str(row.get("PHASE_SET", "")).strip()
                    phase_set_chrom = normalize_chr(phase_set)
                    ps = str(row.get("PS", "")).strip()
                    phase_contract_stats[f"phase_status_{phase_status or 'MISSING'}_rows"] += 1

                    if phase_status == "PHASED_HET_EXACT_MATCH":
                        phase_contract_stats["phase_eligible_rows"] += 1
                        if not chrom or phase_set_chrom != chrom:
                            phase_contract_stats["phase_set_chrom_mismatch_rows"] += 1
                            raise ValueError(
                                f"Chromosome-wide PHASE_SET mismatch at input row {input_rows}: "
                                f"CHROM={row.get('CHROM', '')}; PHASE_SET={phase_set}; "
                                f"Variant={row.get('Variant', '')}; input={ase_tsv}"
                            )
                        if ps in {"", ".", "NA"}:
                            phase_contract_stats["phase_eligible_without_local_ps_rows"] += 1
                        else:
                            phase_contract_stats["phase_eligible_with_local_ps_rows"] += 1
                    elif phase_set not in {"", ".", "NA"}:
                        phase_contract_stats["phase_ineligible_with_phase_set_rows"] += 1
                        raise ValueError(
                            f"Phase-ineligible row has non-empty PHASE_SET at input row {input_rows}: "
                            f"PHASE_STATUS={phase_status}; PHASE_SET={phase_set}; "
                            f"Variant={row.get('Variant', '')}; input={ase_tsv}"
                        )

                    variant_key = (
                        row_animal,
                        row_tissue,
                        chrom,
                        "" if pos is None else str(pos),
                        str(row.get("REF", "")).strip(),
                        str(row.get("ALT", "")).strip(),
                    )

                    if variant_key != previous_variant_key:
                        if variant_key in seen_variants:
                            raise ValueError(
                                f"1.3 input has non-contiguous rows for variant {variant_key}: {ase_tsv}"
                            )
                        seen_variants.add(variant_key)
                        step7_overlaps = find_step7_overlaps(
                            chrom, pos, segments_by_chr
                        )
                        step7_summary = summarize_step7(
                            step7_overlaps, sample_scales, args.consensus_fraction
                        )
                        action, reason = decide_action(
                            evaluable,
                            step7_overlaps,
                            step7_summary,
                        )
                        clean_pass = clean_keep(
                            {"cnv_filter_action": action}, args.clean_mode
                        )
                        previous_cnv_fields = {
                            "cnv_sample_id": sample_id,
                            "cnv_result_status": result_status,
                            "cnv_evaluable": "TRUE" if evaluable else "FALSE",
                            "cnv_variant_coordinate_status": variant_status,
                            "cnv_variant_chr_normalized": chrom,
                            "cnv_variant_pos": "" if pos is None else str(pos),
                            "cnv_step7_overlap_any": "TRUE" if step7_overlaps else "FALSE",
                            "cnv_step7_overlap_states": join_unique(segment["state"] for segment in step7_overlaps),
                            "cnv_step7_overlap_event_scales": join_unique(segment["event_scale"] for segment in step7_overlaps),
                            "cnv_step7_overlap_segments": join_unique(segment["label"] for segment in step7_overlaps),
                            "cnv_step7_overlap_n_segments": str(len(step7_overlaps)),
                            "cnv_step7_available_scales": step7_summary["available_scales"],
                            "cnv_step7_consensus_required_support": str(step7_summary["required_support"]),
                            "cnv_step7_consensus_state": step7_summary["consensus_state"],
                            "cnv_step7_consensus_support": str(step7_summary["consensus_support"]),
                            "cnv_step7_consensus_fraction": step7_summary["consensus_fraction"],
                            "cnv_step7_consensus_scales": step7_summary["consensus_scales"],
                            "cnv_step7_state_support_counts": step7_summary["state_support_counts"],
                            "cnv_step7_best_state": step7_summary["best_state"],
                            "cnv_step7_best_support": str(step7_summary["best_support"]),
                            "cnv_step7_conflict_states": step7_summary["conflict_states"],
                            "cnv_step7_baf_decisions": join_unique(segment["baf_decision"] for segment in step7_overlaps),
                            "cnv_evidence_source": evidence_source(step7_overlaps, step7_summary),
                            "cnv_filter_action": action,
                            "cnv_filter_reason": reason,
                            "cnv_clean_pass": "TRUE" if clean_pass else "FALSE",
                        }
                        previous_variant_key = variant_key
                        previous_clean_pass = clean_pass
                        previous_action = action
                        action_variant_counts[action] += 1
                        if clean_pass:
                            clean_unique_variant_n += 1

                    annotated = dict(row)
                    annotated.update(previous_cnv_fields)
                    action_row_counts[previous_action] += 1
                    full_writer.writerow(annotated)
                    if previous_clean_pass:
                        clean_writer.writerow(annotated)
                        clean_rows += 1

                    if input_rows % args.progress_every == 0:
                        log(f"[PROGRESS] processed_rows={input_rows} clean_rows={clean_rows}")

        log("[STAGE] writing summary")
        metrics = [
            ("animal_id", args.animal_id),
            ("tissue", args.tissue),
            ("sample_id", sample_id),
            ("ase_tsv", str(ase_tsv)),
            ("step7_root", str(args.step7_root)),
            ("full_tsv", str(full_tsv)),
            ("clean_tsv", str(clean_tsv)),
            ("summary_tsv", str(summary_tsv)),
            ("clean_mode", args.clean_mode),
            (
                "clean_policy",
                "strict_keep_only"
                if args.clean_mode == "strict"
                else (
                    "exclude_only_single_state_consensus_cnv_at_"
                    f"{math.ceil(args.consensus_fraction * len(STEP7_SCALE_PAIRS))}"
                    f"_of_{len(STEP7_SCALE_PAIRS)}_scales"
                ),
            ),
            (
                "flag_not_evaluable_policy",
                "excluded_in_strict;retained_with_flag_in_known_cnv_removed",
            ),
            ("consensus_fraction", args.consensus_fraction),
            ("phase_set_scope", "chromosome"),
            ("ps_role", "WhatsHap_local_provenance"),
            (
                "unannotated_id_contract",
                "GeneID/GeneID_base/GeneSymbol/AllGeneIDs=(NA)<Variant> when GeneHitN=0",
            ),
            ("annotation_contract_known_gene_rows", annotation_contract_stats.get("known_gene_rows", 0)),
            ("annotation_contract_manual_new_gene_rows", annotation_contract_stats.get("manual_new_gene_rows", 0)),
            ("annotation_contract_unannotated_rows", annotation_contract_stats.get("unannotated_rows", 0)),
            ("phase_eligible_rows", phase_contract_stats.get("phase_eligible_rows", 0)),
            ("phase_eligible_with_local_ps_rows", phase_contract_stats.get("phase_eligible_with_local_ps_rows", 0)),
            ("phase_eligible_without_local_ps_rows", phase_contract_stats.get("phase_eligible_without_local_ps_rows", 0)),
            ("phase_set_chrom_mismatch_rows", phase_contract_stats.get("phase_set_chrom_mismatch_rows", 0)),
            ("phase_ineligible_with_phase_set_rows", phase_contract_stats.get("phase_ineligible_with_phase_set_rows", 0)),
            ("step7_expected_scales", ",".join(STEP7_SCALE_PAIRS)),
            ("step7_scales_used", ",".join(sort_scales(sample_scales))),
            ("step7_missing_scales", ",".join(sort_scales(missing_scales))),
            ("step7_unexpected_scales", ",".join(sort_scales(unexpected_scales))),
            ("cnv_result_status", result_status),
            ("input_variant_gene_rows", input_rows),
            ("full_variant_gene_rows", input_rows),
            ("clean_variant_gene_rows", clean_rows),
            ("unique_variants", len(seen_variants)),
            ("clean_unique_variants", clean_unique_variant_n),
            ("step7_segment_rows_total", step7_stats.get("segment_rows_total", 0)),
            ("step7_segment_rows_for_sample", step7_stats.get("segment_rows_for_sample", 0)),
            ("step7_abnormal_segments_loaded", step7_stats.get("abnormal_segments_loaded", 0)),
            ("step7_segment_rows_bad_coordinate", step7_stats.get("segment_rows_bad_coordinate", 0)),
        ]
        for status, count in sorted(phase_contract_stats.items()):
            if status.startswith("phase_status_"):
                metrics.append((status, count))
        for action, count in sorted(action_row_counts.items()):
            metrics.append((f"cnv_filter_action_rows.{action}", count))
        for action, count in sorted(action_variant_counts.items()):
            metrics.append((f"cnv_filter_action_variants.{action}", count))

        with summary_tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["metric", "value"])
            writer.writerows(metrics)

        full_tmp.replace(full_tsv)
        clean_tmp.replace(clean_tsv)
        summary_tmp.replace(summary_tsv)
    except Exception:
        for temp_path in (full_tmp, clean_tmp, summary_tmp):
            if temp_path.exists():
                temp_path.unlink()
        raise

    log("[DONE] 1.4 ASE-CNV QC complete / done")
    log(f"[DONE] input_rows={input_rows} clean_rows={clean_rows}")
    log(f"[DONE] full_tsv={full_tsv}")
    log(f"[DONE] clean_tsv={clean_tsv}")
    log(f"[DONE] summary_tsv={summary_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate five-dogs 1.3 full ASE variant-gene rows with Step7 segment CNV evidence."
    )
    parser.add_argument("animal_id")
    parser.add_argument("tissue")
    parser.add_argument("--ase-tsv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--step7-root", type=Path, default=STEP7_ROOT)
    parser.add_argument("--consensus-fraction", type=float, default=0.75)
    parser.add_argument("--clean-mode", choices=["known_cnv_removed", "strict"], default="known_cnv_removed")
    parser.add_argument("--progress-every", type=int, default=50000)
    args = parser.parse_args()

    if not 0 < args.consensus_fraction <= 1:
        raise ValueError("--consensus-fraction must be > 0 and <= 1")
    annotate_ase(args)


if __name__ == "__main__":
    main()
