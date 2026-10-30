#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
from pathlib import Path

CSV_PATH = Path("/path/to/project/legacy/five-dogs/data/metadata/SraRunTable.csv")
ROOT = Path("/path/to/project/legacy/five-dogs/data/raw-rna-data")
OUT_DIR = Path("/path/to/project/legacy/five-dogs/scripts/1.preparation/1")
OUT_DIR.mkdir(parents=True, exist_ok=True)

VALID_DIRS = {
    "Cerebellum",
    "FrontalCortex",
    "PituitaryGland",
    "occipitalcortex",
}

def load_mapping(csv_path: Path):
    run_to_sample = {}
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)

        need = {"Run", "Sample Name"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV missing required columns: {missing}; actual columns is: {reader.fieldnames}")

        for row in reader:
            run = (row.get("Run") or "").strip()
            sample_name = (row.get("Sample Name") or "").strip()
            if run:
                run_to_sample[run] = sample_name

    return run_to_sample

def parse_run_and_suffix(filename: str):

    name = filename

    for suffix in [".fastq.gz", ".fq.gz"]:
        if name.endswith(suffix):
            stem = name[:-len(suffix)]
            if "_" in stem:
                run, tail = stem.split("_", 1)
                return run, "_" + tail + suffix
            return stem, suffix

    return None, None

def tsv_writer(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "brain_region",
            "run",
            "sample_name",
            "old_path",
            "new_path",
            "status",
            "note",
        ])
        w.writerows(rows)

def main():
    parser = argparse.ArgumentParser(
        description="Preview or apply FASTQ renaming under five-dogs raw-rna-data using the Run -> Sample Name mapping in SraRunTable.csv. Preview is the default; --apply is required to change files."
    )
    parser.add_argument("--apply", action="store_true", help="Apply the renaming operations.")
    parser.add_argument("--only-run", nargs="*", default=None, help="Process only the specified runs, for example --only-run SRR5889308 SRR5889305.")
    parser.add_argument("--only-dir", nargs="*", default=None, help="Process only the specified brain-region directories, for example --only-dir Cerebellum.")
    args = parser.parse_args()

    run_to_sample = load_mapping(CSV_PATH)

    preview_rows = []
    apply_rows = []

    dirs = [p for p in sorted(ROOT.iterdir()) if p.is_dir() and p.name in VALID_DIRS]

    if args.only_dir:
        wanted_dirs = set(args.only_dir)
        dirs = [d for d in dirs if d.name in wanted_dirs]

    wanted_runs = set(args.only_run) if args.only_run else None

    for d in dirs:
        brain_region = d.name

        for fp in sorted(d.iterdir()):
            if not fp.is_file():
                continue

            run, suffix = parse_run_and_suffix(fp.name)
            if not run:
                preview_rows.append([
                    brain_region, "", "", str(fp), "", "SKIP", "Filename is not in a supported FASTQ format"
                ])
                continue

            if wanted_runs and run not in wanted_runs:
                continue

            sample_name = run_to_sample.get(run, "")
            if not sample_name:
                preview_rows.append([
                    brain_region, run, "", str(fp), "", "SKIP", "No Sample Name found for the Run in the CSV"
                ])
                continue

            new_name = f"{sample_name}{suffix}"
            new_fp = d / new_name

            status = "PREVIEW"
            note = ""
            if new_fp.exists() and new_fp != fp:
                status = "CONFLICT"
                note = "Target filename already exists"

            preview_rows.append([
                brain_region, run, sample_name, str(fp), str(new_fp), status, note
            ])

            if args.apply:
                if fp == new_fp:
                    apply_rows.append([
                        brain_region, run, sample_name, str(fp), str(new_fp), "SKIP", "Filename does not require renaming"
                    ])
                    continue

                if new_fp.exists():
                    apply_rows.append([
                        brain_region, run, sample_name, str(fp), str(new_fp), "SKIP", "Target file already exists; operation skipped"
                    ])
                    continue

                try:
                    fp.rename(new_fp)
                    apply_rows.append([
                        brain_region, run, sample_name, str(fp), str(new_fp), "RENAMED", ""
                    ])
                except Exception as e:
                    apply_rows.append([
                        brain_region, run, sample_name, str(fp), str(new_fp), "ERROR", str(e)
                    ])

    preview_tsv = OUT_DIR / "five_dogs_fastq_rename.preview.tsv"
    tsv_writer(preview_tsv, preview_rows)
    print(f"[OK] preview file generated: {preview_tsv}")

    if args.apply:
        apply_tsv = OUT_DIR / "five_dogs_fastq_rename.apply.tsv"
        tsv_writer(apply_tsv, apply_rows)
        print(f"[OK] execution result file generated: {apply_tsv}")
        renamed_n = sum(1 for r in apply_rows if r[5] == "RENAMED")
        print(f"[OK] actual completed of rename record count: {renamed_n}")
    else:
        print("[INFO] preview mode is active, no files were modified.")
        print("[INFO] after verifying correctness, use --apply apply the operations.")

if __name__ == "__main__":
    main()
