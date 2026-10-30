#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
from pathlib import Path

CSV_PATH = Path("/path/to/project/legacy/five-dogs/data/metadata/SraRunTable.csv")
ROOT = Path("/path/to/project/legacy/five-dogs/data/WGS")
OUT_DIR = Path("/path/to/project/legacy/five-dogs/scripts/1.preparation/4")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def load_mapping(csv_path: Path):
    run_to_sample = {}
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        need = {"Run", "Sample Name"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SraRunTable.csv missing required columns: {missing}; actual columns is: {reader.fieldnames}")

        for row in reader:
            run = (row.get("Run") or "").strip()
            sample_name = (row.get("Sample Name") or "").strip()
            if run:
                run_to_sample[run] = sample_name

    return run_to_sample

def sample_name_to_animal_id(sample_name: str) -> str:

    s = (sample_name or "").strip()
    parts = s.split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return s

def parse_run_and_suffix(filename: str):

    for suffix in [".fastq.gz", ".fq.gz"]:
        if filename.endswith(suffix):
            stem = filename[:-len(suffix)]
            if "_" not in stem:
                return None, None
            run, rest = stem.split("_", 1)
            return run, "_" + rest + suffix
    return None, None

def write_tsv(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "run",
            "sample_name",
            "animal_id",
            "old_path",
            "new_path",
            "status",
            "note",
        ])
        w.writerows(rows)

def main():
    parser = argparse.ArgumentParser(
        description="Preview or apply FASTQ renaming in the WGS directory using the Run -> Sample Name mapping in SraRunTable.csv. Preview is the default; --apply is required to change files."
    )
    parser.add_argument("--apply", action="store_true", help="Apply the renaming operations.")
    parser.add_argument("--only-run", nargs="*", default=None,
                        help="Process only the specified runs, for example --only-run SRR9001402 SRR9001404.")
    args = parser.parse_args()

    run_to_sample = load_mapping(CSV_PATH)
    wanted_runs = set(args.only_run) if args.only_run else None

    preview_rows = []
    apply_rows = []

    files = sorted([p for p in ROOT.iterdir() if p.is_file()])

    for fp in files:
        name = fp.name


        run, suffix = parse_run_and_suffix(name)
        if not run:
            continue

        if wanted_runs and run not in wanted_runs:
            continue

        sample_name = run_to_sample.get(run, "")
        if not sample_name:
            preview_rows.append([
                run, "", "", str(fp), "", "SKIP", "No Sample Name found for this Run in the CSV"
            ])
            continue

        animal_id = sample_name_to_animal_id(sample_name)
        new_name = f"{animal_id}{suffix}"
        new_fp = ROOT / new_name

        status = "PREVIEW"
        note = ""
        if new_fp.exists() and new_fp != fp:
            status = "CONFLICT"
            note = "Target file already exists"

        preview_rows.append([
            run, sample_name, animal_id, str(fp), str(new_fp), status, note
        ])

        if args.apply:
            if fp == new_fp:
                apply_rows.append([
                    run, sample_name, animal_id, str(fp), str(new_fp), "SKIP", "Filename does not require renaming"
                ])
                continue

            if new_fp.exists():
                apply_rows.append([
                    run, sample_name, animal_id, str(fp), str(new_fp), "SKIP", "Target file already exists; operation skipped"
                ])
                continue

            try:
                fp.rename(new_fp)
                apply_rows.append([
                    run, sample_name, animal_id, str(fp), str(new_fp), "RENAMED", ""
                ])
            except Exception as e:
                apply_rows.append([
                    run, sample_name, animal_id, str(fp), str(new_fp), "ERROR", str(e)
                ])

    preview_tsv = OUT_DIR / "wgs_fastq_rename.preview.tsv"
    write_tsv(preview_tsv, preview_rows)
    print(f"[OK] preview file generated: {preview_tsv}")

    if args.apply:
        apply_tsv = OUT_DIR / "wgs_fastq_rename.apply.tsv"
        write_tsv(apply_tsv, apply_rows)
        print(f"[OK] execution result file generated: {apply_tsv}")
        renamed_n = sum(1 for r in apply_rows if r[5] == "RENAMED")
        print(f"[OK] actual completed of rename record count: {renamed_n}")
    else:
        print("[INFO] preview mode is active, no files were modified.")
        print("[INFO] after verifying correctness, add --apply apply the operations.")

if __name__ == "__main__":
    main()
