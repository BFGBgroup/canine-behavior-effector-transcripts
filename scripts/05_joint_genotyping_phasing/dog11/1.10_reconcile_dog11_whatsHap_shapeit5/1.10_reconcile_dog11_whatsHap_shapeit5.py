#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass
class BlockStats:
    sample_id: str
    chromosome: str
    ps: int
    start: int
    end: int
    whatshap_phased_hets: int = 0
    matched_hets: int = 0
    same_orientation: int = 0
    flipped_orientation: int = 0
    genotype_mismatches: int = 0
    orientation: str = "pending"

    @property
    def concordant_after_orientation(self) -> int:
        if self.orientation in {"keep", "flip"}:
            return max(self.same_orientation, self.flipped_orientation)
        return 0

    @property
    def discordant_after_orientation(self) -> int:
        if self.orientation in {"keep", "flip"}:
            return min(self.same_orientation, self.flipped_orientation)
        return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Orient per-sample WhatsHap PS blocks to the multisample SHAPEIT5 scaffold and output BCFs corrected using local reads."
    )
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--shapeit-bcf", required=True, type=Path)
    parser.add_argument("--whatshap-root", required=True, type=Path)
    parser.add_argument("--output-bcf", required=True, type=Path)
    parser.add_argument("--block-qc", required=True, type=Path)
    parser.add_argument("--corrections", required=True, type=Path)
    parser.add_argument("--sample-summary", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    return parser.parse_args()


def variant_key(record: pysam.VariantRecord) -> tuple[int, str, tuple[str, ...]]:
    return record.pos, record.ref, tuple(record.alts or ())


def format_gt(gt: tuple[int | None, ...] | None) -> str:
    if gt is None:
        return "."
    return "|".join("." if allele is None else str(allele) for allele in gt)


def is_phased_heterozygote(sample: pysam.VariantRecordSample) -> bool:
    gt = sample.get("GT")
    return bool(
        sample.phased
        and gt is not None
        and len(gt) == 2
        and gt[0] is not None
        and gt[1] is not None
        and gt[0] != gt[1]
    )


def choose_orientation(block: BlockStats) -> str:
    anchors = block.same_orientation + block.flipped_orientation
    if anchors == 0:
        return "no_anchor"
    if block.same_orientation == block.flipped_orientation:
        return "ambiguous_tie"
    if block.same_orientation > block.flipped_orientation:
        return "keep"
    return "flip"


def load_whatshap(
    root: Path,
    sample_ids: list[str],
    chromosome: str,
) -> tuple[
    dict[str, dict[tuple[int, str, tuple[str, ...]], tuple[tuple[int, int], int]]],
    dict[str, dict[int, BlockStats]],
]:
    variants_by_sample: dict[
        str, dict[tuple[int, str, tuple[str, ...]], tuple[tuple[int, int], int]]
    ] = {}
    blocks_by_sample: dict[str, dict[int, BlockStats]] = {}

    for sample_id in sample_ids:
        path = root / sample_id / f"{sample_id}.{chromosome}.whatshap.phased.vcf.gz"
        if not path.is_file() or not Path(f"{path}.tbi").is_file():
            raise FileNotFoundError(f"WhatsHap VCF or index does not exist: {path}")

        sample_variants: dict[
            tuple[int, str, tuple[str, ...]], tuple[tuple[int, int], int]
        ] = {}
        sample_blocks: dict[int, BlockStats] = {}

        with pysam.VariantFile(path) as vcf:
            input_samples = list(vcf.header.samples)
            if input_samples != [sample_id]:
                raise ValueError(
                    f"WhatsHap sample does not match: {path}; expected={sample_id}; actual={input_samples}"
                )
            for record in vcf.fetch(chromosome):
                if len(record.alts or ()) != 1:
                    continue
                call = record.samples[sample_id]
                ps_value = call.get("PS")
                if not is_phased_heterozygote(call) or ps_value is None:
                    continue
                gt = call.get("GT")
                assert gt is not None and gt[0] is not None and gt[1] is not None
                ps = int(ps_value)
                key = variant_key(record)
                if key in sample_variants:
                    raise ValueError(f"WhatsHap encountered duplicate variant: {sample_id} {chromosome} {key}")
                sample_variants[key] = ((int(gt[0]), int(gt[1])), ps)

                block = sample_blocks.get(ps)
                if block is None:
                    block = BlockStats(sample_id, chromosome, ps, record.pos, record.pos)
                    sample_blocks[ps] = block
                block.start = min(block.start, record.pos)
                block.end = max(block.end, record.pos)
                block.whatshap_phased_hets += 1

        variants_by_sample[sample_id] = sample_variants
        blocks_by_sample[sample_id] = sample_blocks

    return variants_by_sample, blocks_by_sample


def compare_with_shapeit(
    shapeit_path: Path,
    chromosome: str,
    sample_ids: list[str],
    variants_by_sample: dict[
        str, dict[tuple[int, str, tuple[str, ...]], tuple[tuple[int, int], int]]
    ],
    blocks_by_sample: dict[str, dict[int, BlockStats]],
) -> int:
    record_count = 0
    with pysam.VariantFile(shapeit_path) as shapeit:
        if list(shapeit.header.samples) != sample_ids:
            raise ValueError("SHAPEIT5 sample order in repeated reading when changed")
        for record in shapeit.fetch(chromosome):
            record_count += 1
            key = variant_key(record)
            for sample_id in sample_ids:
                whatshap_value = variants_by_sample[sample_id].get(key)
                if whatshap_value is None:
                    continue
                whatshap_gt, ps = whatshap_value
                block = blocks_by_sample[sample_id][ps]
                call = record.samples[sample_id]
                shapeit_gt = call.get("GT")
                if not is_phased_heterozygote(call) or shapeit_gt is None:
                    block.genotype_mismatches += 1
                    continue
                shapeit_pair = (int(shapeit_gt[0]), int(shapeit_gt[1]))
                if shapeit_pair == whatshap_gt:
                    block.same_orientation += 1
                    block.matched_hets += 1
                elif shapeit_pair == tuple(reversed(whatshap_gt)):
                    block.flipped_orientation += 1
                    block.matched_hets += 1
                else:
                    block.genotype_mismatches += 1

    for sample_blocks in blocks_by_sample.values():
        for block in sample_blocks.values():
            block.orientation = choose_orientation(block)
    return record_count


def write_reconciled_bcf(
    shapeit_path: Path,
    chromosome: str,
    output_path: Path,
    corrections_path: Path,
    sample_ids: list[str],
    variants_by_sample: dict[
        str, dict[tuple[int, str, tuple[str, ...]], tuple[tuple[int, int], int]]
    ],
    blocks_by_sample: dict[str, dict[int, BlockStats]],
) -> tuple[int, dict[str, dict[str, int]]]:
    per_sample: dict[str, dict[str, int]] = {
        sample_id: defaultdict(int) for sample_id in sample_ids
    }

    with pysam.VariantFile(shapeit_path) as source:
        header = source.header.copy()
        if "PS" not in header.formats:
            header.formats.add(
                "PS",
                number=1,
                type="Integer",
                description="WhatsHap phase set oriented to the SHAPEIT5 chromosome-wide scaffold",
            )

        with pysam.VariantFile(output_path, "wb", header=header) as output, gzip.open(
            corrections_path, "wt", encoding="utf-8", newline=""
        ) as corrections_handle:
            correction_writer = csv.writer(corrections_handle, delimiter="\t", lineterminator="\n")
            correction_writer.writerow(
                [
                    "sample_id",
                    "chromosome",
                    "position",
                    "ref",
                    "alt",
                    "ps",
                    "shapeit_gt_before",
                    "whatshap_oriented_gt",
                    "block_orientation",
                    "block_same_sites",
                    "block_flipped_sites",
                ]
            )

            output_records = 0
            for record in source.fetch(chromosome):
                output_records += 1
                record.translate(header)
                key = variant_key(record)
                for sample_id in sample_ids:
                    whatshap_value = variants_by_sample[sample_id].get(key)
                    if whatshap_value is None:
                        continue
                    whatshap_gt, ps = whatshap_value
                    block = blocks_by_sample[sample_id][ps]
                    if block.orientation not in {"keep", "flip"}:
                        continue

                    oriented_gt = (
                        whatshap_gt
                        if block.orientation == "keep"
                        else tuple(reversed(whatshap_gt))
                    )
                    call = record.samples[sample_id]
                    shapeit_gt = call.get("GT")
                    if not is_phased_heterozygote(call) or shapeit_gt is None:
                        continue
                    shapeit_pair = (int(shapeit_gt[0]), int(shapeit_gt[1]))
                    if shapeit_pair not in {oriented_gt, tuple(reversed(oriented_gt))}:
                        continue

                    call["PS"] = ps
                    per_sample[sample_id]["ps_annotated_sites"] += 1
                    if shapeit_pair != oriented_gt:
                        call["GT"] = oriented_gt
                        call.phased = True
                        if "PP" in header.formats:
                            call["PP"] = None
                        per_sample[sample_id]["corrected_sites"] += 1
                        correction_writer.writerow(
                            [
                                sample_id,
                                chromosome,
                                record.pos,
                                record.ref,
                                record.alts[0],
                                ps,
                                format_gt(shapeit_pair),
                                format_gt(oriented_gt),
                                block.orientation,
                                block.same_orientation,
                                block.flipped_orientation,
                            ]
                        )
                output.write(record)

    return output_records, per_sample


def write_block_qc(
    path: Path,
    sample_ids: list[str],
    blocks_by_sample: dict[str, dict[int, BlockStats]],
) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "chromosome",
                "ps",
                "start",
                "end",
                "whatshap_phased_hets",
                "matched_hets",
                "same_orientation",
                "flipped_orientation",
                "genotype_mismatches",
                "orientation",
                "concordant_after_orientation",
                "discordant_after_orientation",
                "concordance_after_orientation",
            ]
        )
        for sample_id in sample_ids:
            for ps in sorted(blocks_by_sample[sample_id]):
                block = blocks_by_sample[sample_id][ps]
                anchors = block.same_orientation + block.flipped_orientation
                concordance = (
                    f"{block.concordant_after_orientation / anchors:.6f}"
                    if anchors and block.orientation in {"keep", "flip"}
                    else "NA"
                )
                writer.writerow(
                    [
                        block.sample_id,
                        block.chromosome,
                        block.ps,
                        block.start,
                        block.end,
                        block.whatshap_phased_hets,
                        block.matched_hets,
                        block.same_orientation,
                        block.flipped_orientation,
                        block.genotype_mismatches,
                        block.orientation,
                        block.concordant_after_orientation,
                        block.discordant_after_orientation,
                        concordance,
                    ]
                )


def summarize_sample(
    sample_id: str,
    chromosome: str,
    blocks: dict[int, BlockStats],
    write_counts: dict[str, int],
) -> dict[str, int | str]:
    values: dict[str, int | str] = {
        "sample_id": sample_id,
        "chromosome": chromosome,
        "whatshap_phased_hets": sum(block.whatshap_phased_hets for block in blocks.values()),
        "ps_blocks": len(blocks),
        "matched_hets": sum(block.matched_hets for block in blocks.values()),
        "genotype_mismatches": sum(block.genotype_mismatches for block in blocks.values()),
        "oriented_blocks": sum(block.orientation in {"keep", "flip"} for block in blocks.values()),
        "keep_blocks": sum(block.orientation == "keep" for block in blocks.values()),
        "flip_blocks": sum(block.orientation == "flip" for block in blocks.values()),
        "ambiguous_tie_blocks": sum(block.orientation == "ambiguous_tie" for block in blocks.values()),
        "no_anchor_blocks": sum(block.orientation == "no_anchor" for block in blocks.values()),
        "discordant_after_orientation": sum(
            block.discordant_after_orientation for block in blocks.values()
        ),
        "ps_annotated_sites": write_counts.get("ps_annotated_sites", 0),
        "corrected_sites": write_counts.get("corrected_sites", 0),
    }
    return values


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    chromosome = args.chromosome
    if not chromosome.startswith("chr"):
        raise ValueError(f"Chromosome names must start with chr: {chromosome}")
    if not args.shapeit_bcf.is_file() or not Path(f"{args.shapeit_bcf}.csi").is_file():
        raise FileNotFoundError(f"SHAPEIT5 BCF or index does not exist: {args.shapeit_bcf}")

    with pysam.VariantFile(args.shapeit_bcf) as shapeit:
        sample_ids = list(shapeit.header.samples)
    if len(sample_ids) != 11:
        raise ValueError(f"SHAPEIT5 sample count is not 11: {len(sample_ids)}")

    variants_by_sample, blocks_by_sample = load_whatshap(
        args.whatshap_root, sample_ids, chromosome
    )
    shapeit_records = compare_with_shapeit(
        args.shapeit_bcf,
        chromosome,
        sample_ids,
        variants_by_sample,
        blocks_by_sample,
    )
    output_records, per_sample_write_counts = write_reconciled_bcf(
        args.shapeit_bcf,
        chromosome,
        args.output_bcf,
        args.corrections,
        sample_ids,
        variants_by_sample,
        blocks_by_sample,
    )
    if output_records != shapeit_records:
        raise ValueError(
            f"Output record count differs from the SHAPEIT5 input: input={shapeit_records}; output={output_records}"
        )

    write_block_qc(args.block_qc, sample_ids, blocks_by_sample)
    sample_rows = [
        summarize_sample(
            sample_id,
            chromosome,
            blocks_by_sample[sample_id],
            per_sample_write_counts[sample_id],
        )
        for sample_id in sample_ids
    ]
    sample_fields = list(sample_rows[0])
    write_tsv(args.sample_summary, sample_fields, sample_rows)

    metric_row: dict[str, object] = {
        "chromosome": chromosome,
        "shapeit_records": shapeit_records,
        "output_records": output_records,
        "output_samples": len(sample_ids),
    }
    for field in sample_fields[2:]:
        metric_row[field] = sum(int(row[field]) for row in sample_rows)
    write_tsv(args.metrics, list(metric_row), [metric_row])
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        raise
