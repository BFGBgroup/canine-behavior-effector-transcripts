#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import annotations

import argparse
import csv
import logging
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


PROJECT_ROOT = Path("/path/to/project/legacy/five-dogs")
STEP12_ROOT = PROJECT_ROOT / "9.new-gene/1.2_check_bulk_star_bam_discovery_readiness-v3"

DEFAULT_READINESS_TSV = STEP12_ROOT / "bam_discovery_readiness.tsv"
DEFAULT_GUIDE_GTF = Path("/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf")
DEFAULT_STRINGTIE = Path("/path/to/software/legacy_envs/stringtie3/bin/stringtie")
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "9.new-gene/1.5_stringtie3_bulk_transcript_discovery-v6"

READINESS_REQUIRED_FIELDS = [
    "animal_id",
    "tissue",
    "sample_tissue",
    "raw_star_bam",
    "raw_star_bam_exists",
    "raw_star_bam_index_exists",
    "recommended_for_stringtie3_discovery",
    "has_XS_tag_in_spliced_reads",
]

MANIFEST_FIELDS = [
    "sample",
    "tissue",
    "sample_tissue",
    "bam_path",
    "strandness",
    "stringtie_strand_option",
    "stringtie_nascent_option",
    "stringtie_label",
    "output_gtf",
    "output_gene_abund",
    "log_file",
]

STATUS_FIELDS = [
    "sample",
    "tissue",
    "sample_tissue",
    "bam_path",
    "gtf_path",
    "gene_abund_path",
    "log_file",
    "status",
    "return_code",
    "runtime_sec",
    "transcript_count",
    "exon_count",
    "error_message",
]


@dataclass(frozen=True)
class Sample:
    sample: str
    tissue: str
    sample_tissue: str
    bam_path: Path
    label: str
    output_gtf: Path
    output_gene_abund: Path
    log_file: Path


@dataclass
class SampleStatus:
    sample: str
    tissue: str
    sample_tissue: str
    bam_path: str
    gtf_path: str
    gene_abund_path: str
    log_file: str
    status: str
    return_code: str
    runtime_sec: str
    transcript_count: str
    exon_count: str
    error_message: str

    def as_row(self) -> Dict[str, str]:
        return {
            "sample": self.sample,
            "tissue": self.tissue,
            "sample_tissue": self.sample_tissue,
            "bam_path": self.bam_path,
            "gtf_path": self.gtf_path,
            "gene_abund_path": self.gene_abund_path,
            "log_file": self.log_file,
            "status": self.status,
            "return_code": self.return_code,
            "runtime_sec": self.runtime_sec,
            "transcript_count": self.transcript_count,
            "exon_count": self.exon_count,
            "error_message": self.error_message,
        }


def clean(value: object) -> str:
    text = "" if value is None else str(value).strip()
    return text if text else "NA"


def fail(message: str) -> None:
    raise RuntimeError(message)


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


def require_fields(path: Path, fieldnames: Sequence[str], required: Sequence[str]) -> None:
    missing = [field for field in required if field not in fieldnames]
    if missing:
        fail(f"{path} missing required fields: {','.join(missing)}")


def sanitize_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_]+", "_", value.strip())
    label = re.sub(r"_+", "_", label).strip("_")
    if not label:
        fail(f"Cannot generate a StringTie label from sample_tissue: {value!r}")
    if not re.match(r"^[A-Za-z]", label):
        label = f"st_{label}"
    return label


def resolve_executable(path: Path) -> Path:
    text = str(path)
    if "/" in text:
        if not path.exists():
            fail(f"StringTie executable does not exist: {path}")
        if not path.is_file():
            fail(f"StringTie path is not a regular file: {path}")
        if not path.stat().st_mode & 0o111:
            fail(f"StringTie file not executable: {path}")
        return path
    resolved = shutil.which(text)
    if not resolved:
        fail(f"Cannot find StringTie in PATH: {text}")
    return Path(resolved)


def make_paths(outdir: Path) -> Dict[str, Path]:
    per_sample_dir = outdir / "stringtie3_per_sample"
    return {
        "per_sample_dir": per_sample_dir,
        "input_manifest": outdir / "stringtie3_input_manifest.tsv",
        "status": outdir / "stringtie3_per_sample_status.tsv",
        "mergelist": outdir / "stringtie3_mergelist.txt",
        "merged_gtf": outdir / "merged.stringtie3.gtf",
        "merge_log": outdir / "merged.stringtie3.log",
        "qc": outdir / "stringtie3_discovery_qc_summary.tsv",
        "manifest": outdir / "stringtie3_discovery_manifest.tsv",
    }


def prepare_output(outdir: Path, paths: Mapping[str, Path], overwrite: bool, resume: bool) -> None:
    if outdir.exists() and not outdir.is_dir():
        fail(f"output path already exists but is not a directory: {outdir}")
    if outdir.exists() and not overwrite and not resume:
        protected = [
            paths["input_manifest"],
            paths["status"],
            paths["mergelist"],
            paths["merged_gtf"],
            paths["qc"],
            paths["manifest"],
        ]
        existing = [str(path) for path in protected if path.exists()]
        if existing:
            fail("output directory existing main result file; such as requires continue please use --resume, such as requires rerun please use --overwrite: " + "|".join(existing))
    outdir.mkdir(parents=True, exist_ok=True)
    paths["per_sample_dir"].mkdir(parents=True, exist_ok=True)


def read_samples(readiness_tsv: Path, per_sample_dir: Path) -> Tuple[List[Sample], int]:
    if not readiness_tsv.exists():
        fail(f"readiness TSV does not exist: {readiness_tsv}")
    samples: List[Sample] = []
    total_rows = 0
    with readiness_tsv.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(readiness_tsv, reader.fieldnames or [], READINESS_REQUIRED_FIELDS)
        for row in reader:
            total_rows += 1
            if clean(row.get("recommended_for_stringtie3_discovery")) != "1":
                continue
            if clean(row.get("raw_star_bam_exists")) != "1":
                continue
            if clean(row.get("raw_star_bam_index_exists")) != "1":
                continue
            if clean(row.get("has_XS_tag_in_spliced_reads")) != "1":
                continue
            sample = clean(row.get("animal_id"))
            tissue = clean(row.get("tissue"))
            sample_tissue = clean(row.get("sample_tissue"))
            bam_path = Path(clean(row.get("raw_star_bam")))
            if not bam_path.exists():
                fail(f"readiness in recommended sample of BAM does not exist: sample_tissue={sample_tissue} bam={bam_path}")
            label = sanitize_label(sample_tissue)
            samples.append(
                Sample(
                    sample=sample,
                    tissue=tissue,
                    sample_tissue=sample_tissue,
                    bam_path=bam_path,
                    label=label,
                    output_gtf=per_sample_dir / f"{sample_tissue}.stringtie3.gtf",
                    output_gene_abund=per_sample_dir / f"{sample_tissue}.stringtie3.gene_abund.tsv",
                    log_file=per_sample_dir / f"{sample_tissue}.stringtie3.log",
                )
            )
    if not samples:
        fail(f"no recommended_for_stringtie3_discovery == 1 and XS tag available of sample: {readiness_tsv}")
    duplicated = [name for name, count in count_values(sample.sample_tissue for sample in samples).items() if count > 1]
    if duplicated:
        fail("sample_tissue is not unique: " + "|".join(sorted(duplicated)))
    return samples, total_rows


def count_values(values: Iterable[str]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts


def count_gtf_features(path: Path) -> Tuple[int, int]:
    transcript_count = 0
    exon_count = 0
    if not path.exists():
        return transcript_count, exon_count
    with path.open("r") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            if fields[2] == "transcript":
                transcript_count += 1
            elif fields[2] == "exon":
                exon_count += 1
    return transcript_count, exon_count


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: clean(row.get(field, "NA")) for field in fieldnames})


def write_input_manifest(path: Path, samples: Sequence[Sample], strand_option: str, nascent_option: str) -> None:
    strandness = "unstranded_or_unknown" if not strand_option else strand_option.lstrip("-")
    nascent_mode = "off" if not nascent_option else nascent_option
    rows = [
        {
            "sample": sample.sample,
            "tissue": sample.tissue,
            "sample_tissue": sample.sample_tissue,
            "bam_path": str(sample.bam_path),
            "strandness": strandness,
            "stringtie_strand_option": strand_option or "NA",
            "stringtie_nascent_option": nascent_mode,
            "stringtie_label": sample.label,
            "output_gtf": str(sample.output_gtf),
            "output_gene_abund": str(sample.output_gene_abund),
            "log_file": str(sample.log_file),
        }
        for sample in samples
    ]
    write_tsv(path, MANIFEST_FIELDS, rows)


def run_one_sample(
    sample: Sample,
    stringtie: Path,
    guide_gtf: Path,
    threads: int,
    strand_option: str,
    nascent_option: str,
    overwrite: bool,
    resume: bool,
) -> SampleStatus:
    start = time.time()
    if resume and sample.output_gtf.exists() and sample.output_gene_abund.exists():
        transcript_count, exon_count = count_gtf_features(sample.output_gtf)
        if transcript_count > 0 and exon_count > 0:
            return SampleStatus(
                sample=sample.sample,
                tissue=sample.tissue,
                sample_tissue=sample.sample_tissue,
                bam_path=str(sample.bam_path),
                gtf_path=str(sample.output_gtf),
                gene_abund_path=str(sample.output_gene_abund),
                log_file=str(sample.log_file),
                status="success_reused",
                return_code="0",
                runtime_sec=f"{time.time() - start:.3f}",
                transcript_count=str(transcript_count),
                exon_count=str(exon_count),
                error_message="NA",
            )
    if not overwrite and not resume and (sample.output_gtf.exists() or sample.output_gene_abund.exists()):
        return SampleStatus(
            sample=sample.sample,
            tissue=sample.tissue,
            sample_tissue=sample.sample_tissue,
            bam_path=str(sample.bam_path),
            gtf_path=str(sample.output_gtf),
            gene_abund_path=str(sample.output_gene_abund),
            log_file=str(sample.log_file),
            status="failed",
            return_code="NA",
            runtime_sec=f"{time.time() - start:.3f}",
            transcript_count="NA",
            exon_count="NA",
            error_message="per-sample output exists; use --resume or --overwrite",
        )
    cmd = [
        str(stringtie),
        str(sample.bam_path),
        "-p",
        str(threads),
        "-G",
        str(guide_gtf),
        "-o",
        str(sample.output_gtf),
        "-A",
        str(sample.output_gene_abund),
        "-l",
        sample.label,
    ]
    if strand_option:
        cmd.append(strand_option)
    if nascent_option:
        cmd.append(nascent_option)
    try:
        with sample.log_file.open("w") as log_handle:
            log_handle.write("[cmd] " + " ".join(cmd) + "\n")
            log_handle.flush()
            result = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=False)
    except Exception as exc:
        return SampleStatus(
            sample=sample.sample,
            tissue=sample.tissue,
            sample_tissue=sample.sample_tissue,
            bam_path=str(sample.bam_path),
            gtf_path=str(sample.output_gtf),
            gene_abund_path=str(sample.output_gene_abund),
            log_file=str(sample.log_file),
            status="failed",
            return_code="NA",
            runtime_sec=f"{time.time() - start:.3f}",
            transcript_count="NA",
            exon_count="NA",
            error_message=str(exc).replace("\t", " "),
        )
    transcript_count, exon_count = count_gtf_features(sample.output_gtf)
    status = "success" if result.returncode == 0 and transcript_count > 0 and exon_count > 0 else "failed"
    error_message = "NA" if status == "success" else "StringTie failed or produced empty GTF"
    return SampleStatus(
        sample=sample.sample,
        tissue=sample.tissue,
        sample_tissue=sample.sample_tissue,
        bam_path=str(sample.bam_path),
        gtf_path=str(sample.output_gtf),
        gene_abund_path=str(sample.output_gene_abund),
        log_file=str(sample.log_file),
        status=status,
        return_code=str(result.returncode),
        runtime_sec=f"{time.time() - start:.3f}",
        transcript_count=str(transcript_count),
        exon_count=str(exon_count),
        error_message=error_message,
    )


def write_status(path: Path, statuses: Sequence[SampleStatus]) -> None:
    write_tsv(path, STATUS_FIELDS, [status.as_row() for status in sorted(statuses, key=lambda item: item.sample_tissue)])


def run_samples(
    samples: Sequence[Sample],
    stringtie: Path,
    guide_gtf: Path,
    threads_per_sample: int,
    workers: int,
    strand_option: str,
    nascent_option: str,
    overwrite: bool,
    resume: bool,
    status_path: Path,
) -> List[SampleStatus]:
    logging.info("start per sample StringTie3 assembly: samples=%s workers=%s threads_per_sample=%s", len(samples), workers, threads_per_sample)
    statuses: List[SampleStatus] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                run_one_sample,
                sample,
                stringtie,
                guide_gtf,
                threads_per_sample,
                strand_option,
                nascent_option,
                overwrite,
                resume,
            ): sample
            for sample in samples
        }
        for index, future in enumerate(as_completed(futures), start=1):
            sample = futures[future]
            try:
                status = future.result()
            except Exception as exc:
                status = SampleStatus(
                    sample=sample.sample,
                    tissue=sample.tissue,
                    sample_tissue=sample.sample_tissue,
                    bam_path=str(sample.bam_path),
                    gtf_path=str(sample.output_gtf),
                    gene_abund_path=str(sample.output_gene_abund),
                    log_file=str(sample.log_file),
                    status="failed",
                    return_code="NA",
                    runtime_sec="NA",
                    transcript_count="NA",
                    exon_count="NA",
                    error_message=str(exc).replace("\t", " "),
                )
            statuses.append(status)
            write_status(status_path, statuses)
            logging.info("sample completed collect: %s status=%s progress=%s/%s", status.sample_tissue, status.status, index, len(samples))
    return statuses


def write_mergelist(path: Path, statuses: Sequence[SampleStatus]) -> List[Path]:
    success_status = {"success", "success_reused"}
    gtf_paths = [Path(status.gtf_path) for status in sorted(statuses, key=lambda item: item.sample_tissue) if status.status in success_status]
    with path.open("w") as handle:
        for gtf_path in gtf_paths:
            handle.write(str(gtf_path) + "\n")
    return gtf_paths


def run_merge(stringtie: Path, guide_gtf: Path, merge_threads: int, mergelist: Path, merged_gtf: Path, merge_log: Path) -> int:
    cmd = [
        str(stringtie),
        "--merge",
        "-p",
        str(merge_threads),
        "-G",
        str(guide_gtf),
        "-o",
        str(merged_gtf),
        str(mergelist),
    ]
    with merge_log.open("w") as log_handle:
        log_handle.write("[cmd] " + " ".join(cmd) + "\n")
        log_handle.flush()
        result = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=False)
    return int(result.returncode)


def write_qc(
    path: Path,
    input_row_count: int,
    samples: Sequence[Sample],
    statuses: Sequence[SampleStatus],
    gtf_paths: Sequence[Path],
    merged_gtf: Path,
    merge_return_code: str,
    expected_sample_count: int,
) -> None:
    success_status = {"success", "success_reused"}
    successful = [status for status in statuses if status.status in success_status]
    failed = [status for status in statuses if status.status not in success_status]
    merged_transcripts, merged_exons = count_gtf_features(merged_gtf)
    rows = [
        {"metric": "input_readiness_rows", "value": input_row_count, "note": "bam_discovery_readiness.tsv data rows"},
        {"metric": "selected_sample_count", "value": len(samples), "note": "recommended_for_stringtie3_discovery == 1 and XS tag detected"},
        {"metric": "expected_sample_count", "value": expected_sample_count, "note": "formal PASS expected sample count"},
        {"metric": "successful_sample_count", "value": len(successful), "note": "per-sample StringTie status success/success_reused"},
        {"metric": "failed_sample_count", "value": len(failed), "note": "per-sample StringTie status failed"},
        {"metric": "per_sample_gtf_count", "value": len(gtf_paths), "note": "GTF files written to mergelist"},
        {"metric": "merge_return_code", "value": merge_return_code, "note": "StringTie --merge return code"},
        {"metric": "merged_gtf_exists", "value": "1" if merged_gtf.exists() else "0", "note": str(merged_gtf)},
        {"metric": "merged_gtf_transcript_count", "value": merged_transcripts, "note": "feature == transcript"},
        {"metric": "merged_gtf_exon_count", "value": merged_exons, "note": "feature == exon"},
        {
            "metric": "pass_status",
            "value": "PASS" if len(successful) == expected_sample_count and not failed and merge_return_code == "0" and merged_transcripts > 0 else "FAIL",
            "note": "PASS requires all expected samples and non-empty merged GTF",
        },
    ]
    write_tsv(path, ["metric", "value", "note"], rows)


def write_run_manifest(path: Path, args: argparse.Namespace, stringtie: Path, paths: Mapping[str, Path], samples: Sequence[Sample]) -> None:
    rows = [
        {"key": "readiness_tsv", "value": str(args.readiness_tsv)},
        {"key": "guide_gtf", "value": str(args.guide_gtf)},
        {"key": "stringtie", "value": str(stringtie)},
        {"key": "outdir", "value": str(args.outdir)},
        {"key": "input_manifest", "value": str(paths["input_manifest"])},
        {"key": "per_sample_status", "value": str(paths["status"])},
        {"key": "mergelist", "value": str(paths["mergelist"])},
        {"key": "merged_gtf", "value": str(paths["merged_gtf"])},
        {"key": "qc_summary", "value": str(paths["qc"])},
        {"key": "workers", "value": str(args.workers)},
        {"key": "threads_per_sample", "value": str(args.threads_per_sample)},
        {"key": "merge_threads", "value": str(args.merge_threads)},
        {"key": "strand_option", "value": args.strand_option or "NA"},
        {"key": "nascent_option", "value": args.nascent_option or "NA"},
        {"key": "expected_sample_count", "value": str(args.expected_sample_count)},
        {"key": "overwrite", "value": str(args.overwrite)},
        {"key": "resume", "value": str(args.resume)},
        {"key": "wasp_bam_used", "value": "no"},
        {"key": "stringtie_e_mode_used", "value": "no"},
        {"key": "assembly_scope", "value": "whole_genome"},
        {"key": "sample_tissue_list", "value": "|".join(sample.sample_tissue for sample in samples)},
    ]
    write_tsv(path, ["key", "value"], rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Step 1.5: run whole-genome StringTie3 transcript discovery for five-dogs bulk RNA raw STAR BAMs."
    )
    parser.add_argument("--readiness-tsv", type=Path, default=DEFAULT_READINESS_TSV)
    parser.add_argument("--guide-gtf", type=Path, default=DEFAULT_GUIDE_GTF)
    parser.add_argument("--stringtie", type=Path, default=DEFAULT_STRINGTIE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--threads-per-sample", type=int, default=8)
    parser.add_argument("--merge-threads", type=int, default=8)
    parser.add_argument("--strand-option", choices=["", "--rf", "--fr"], default="")
    parser.add_argument("--nascent-option", choices=["", "-N", "--nasc"], default="")
    parser.add_argument("--expected-sample-count", type=int, default=19)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging()
    logging.info("Step 1.5 StringTie3 transcript discovery start")
    if args.workers < 1:
        fail("--workers must be >= 1")
    if args.threads_per_sample < 1:
        fail("--threads-per-sample must be >= 1")
    if args.merge_threads < 1:
        fail("--merge-threads must be >= 1")
    if args.resume and args.overwrite:
        fail("--resume and --overwrite cannot be used together")
    if not args.guide_gtf.exists():
        fail(f"guide GTF does not exist: {args.guide_gtf}")
    stringtie = resolve_executable(args.stringtie)
    paths = make_paths(args.outdir)
    prepare_output(args.outdir, paths, overwrite=args.overwrite, resume=args.resume)

    logging.info("read Step 1.2 v3 readiness: %s", args.readiness_tsv)
    samples, input_row_count = read_samples(args.readiness_tsv, paths["per_sample_dir"])
    logging.info("filter to StringTie3 input sample: %s / input_rows=%s", len(samples), input_row_count)
    if len(samples) != args.expected_sample_count:
        fail(f"recommended sample count does not equal expected_sample_count: selected={len(samples)} expected={args.expected_sample_count}")

    logging.info("write StringTie3 input manifest: %s", paths["input_manifest"])
    write_input_manifest(paths["input_manifest"], samples, args.strand_option, args.nascent_option)
    write_run_manifest(paths["manifest"], args, stringtie, paths, samples)

    statuses = run_samples(
        samples=samples,
        stringtie=stringtie,
        guide_gtf=args.guide_gtf,
        threads_per_sample=args.threads_per_sample,
        workers=args.workers,
        strand_option=args.strand_option,
        nascent_option=args.nascent_option,
        overwrite=args.overwrite,
        resume=args.resume,
        status_path=paths["status"],
    )
    write_status(paths["status"], statuses)
    failed = [status for status in statuses if status.status not in {"success", "success_reused"}]
    gtf_paths = write_mergelist(paths["mergelist"], statuses)
    merge_return_code = "NA"
    if failed:
        logging.error("Some samples failed; StringTie --merge will not run: failed=%s", len(failed))
        write_qc(paths["qc"], input_row_count, samples, statuses, gtf_paths, paths["merged_gtf"], merge_return_code, args.expected_sample_count)
        return 1
    if len(gtf_paths) != args.expected_sample_count:
        logging.error("mergelist GTF count does not equal expected_sample_count: gtf=%s expected=%s", len(gtf_paths), args.expected_sample_count)
        write_qc(paths["qc"], input_row_count, samples, statuses, gtf_paths, paths["merged_gtf"], merge_return_code, args.expected_sample_count)
        return 1

    logging.info("run StringTie --merge: mergelist=%s", paths["mergelist"])
    merge_return_code = str(run_merge(stringtie, args.guide_gtf, args.merge_threads, paths["mergelist"], paths["merged_gtf"], paths["merge_log"]))
    write_qc(paths["qc"], input_row_count, samples, statuses, gtf_paths, paths["merged_gtf"], merge_return_code, args.expected_sample_count)
    merged_transcripts, merged_exons = count_gtf_features(paths["merged_gtf"])
    if merge_return_code != "0" or merged_transcripts == 0 or merged_exons == 0:
        logging.error("StringTie --merge failed: rc=%s transcripts=%s exons=%s", merge_return_code, merged_transcripts, merged_exons)
        return 1
    logging.info("Step 1.5 StringTie3 transcript discovery completed: PASS merged_transcripts=%s merged_exons=%s", merged_transcripts, merged_exons)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        logging.error("%s", exc)
        raise SystemExit(1)
