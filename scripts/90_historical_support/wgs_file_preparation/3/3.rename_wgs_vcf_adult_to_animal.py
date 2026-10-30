#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import argparse
import re
from pathlib import Path

TABLES1_CSV = Path("/path/to/project/legacy/five-dogs/data/metadata/TableS1.csv")
SRARUN_CSV  = Path("/path/to/project/legacy/five-dogs/data/metadata/SraRunTable.csv")
ROOT        = Path("/path/to/project/legacy/five-dogs/data/WGS-VCF")
OUT_DIR     = Path("/path/to/project/legacy/five-dogs/scripts/1.preparation/3")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def load_srarun_mapping(csv_path: Path):
    run_to_sample = {}
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        need = {"Run", "Sample Name"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SraRunTable.csv missing required columns: {missing}; actual columns: {reader.fieldnames}")

        for row in reader:
            run = (row.get("Run") or "").strip()
            sample_name = (row.get("Sample Name") or "").strip()
            if run:
                run_to_sample[run] = sample_name
    return run_to_sample

def parse_tables1_wgs_section(csv_path: Path):

    rows = []
    in_section = False
    header = None

    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        for raw in reader:
            row = [x.strip() for x in raw]

            if not any(row):
                if in_section and header is not None:
                    break
                continue

            first = row[0] if len(row) > 0 else ""

            if first == "A. Whole genome sequencing":
                in_section = True
                header = None
                continue

            if not in_section:
                continue

            if header is None:
                if "Individual" in row and "SRA Accession" in row:
                    header = row
                continue

            record = {}
            for i, col in enumerate(header):
                record[col] = row[i].strip() if i < len(row) else ""

            individual = record.get("Individual", "")
            sra = record.get("SRA Accession", "")

            if individual and sra:
                rows.append({
                    "Individual": individual,
                    "SRA Accession": sra
                })

    return rows

def sample_name_to_animal_id(sample_name: str) -> str:

    s = (sample_name or "").strip()
    parts = s.split("_")
    if len(parts) >= 2:
        return parts[0] + "_" + parts[1]
    return s

def individual_to_old_basename(individual: str) -> str:

    m = re.fullmatch(r"Adult(\d+)", individual)
    if not m:
        raise ValueError(f"cannot from Individual parse filename: {individual}")
    n = m.group(1)
    return f"Adult_{n}.WGS.vcf.gz"

def write_tsv(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "individual",
            "sra_accession",
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
        description="Use TableS1.csv + SraRunTable.csv to preview or rename Adult_1~5.WGS.vcf.gz as animal_id.WGS.vcf.gz. Preview is the default; --apply is required to change files."
    )
    parser.add_argument("--apply", action="store_true", help="Apply the renaming operations.")
    parser.add_argument("--only-individual", nargs="*", default=None,
                        help="Process only the specified Individual values, for example --only-individual Adult1 Adult2.")
    args = parser.parse_args()

    run_to_sample = load_srarun_mapping(SRARUN_CSV)
    table_rows = parse_tables1_wgs_section(TABLES1_CSV)

    if args.only_individual:
        wanted = set(args.only_individual)
        table_rows = [r for r in table_rows if r["Individual"] in wanted]

    preview_rows = []
    apply_rows = []

    for rec in table_rows:
        individual = rec["Individual"]
        sra = rec["SRA Accession"]

        sample_name = run_to_sample.get(sra, "")
        if not sample_name:
            preview_rows.append([
                individual, sra, "", "", "", "", "SKIP",
                "No Sample Name found for this SRA Accession in SraRunTable.csv"
            ])
            continue

        animal_id = sample_name_to_animal_id(sample_name)

        old_name = individual_to_old_basename(individual)
        new_name = f"{animal_id}.WGS.vcf.gz"

        old_fp = ROOT / old_name
        new_fp = ROOT / new_name

        status = "PREVIEW"
        note = ""
        if not old_fp.exists():
            status = "SKIP"
            note = "Source file does not exist"
        elif new_fp.exists() and new_fp != old_fp:
            status = "CONFLICT"
            note = "Target file already exists"

        preview_rows.append([
            individual, sra, sample_name, animal_id,
            str(old_fp), str(new_fp), status, note
        ])

        if args.apply:
            if not old_fp.exists():
                apply_rows.append([
                    individual, sra, sample_name, animal_id,
                    str(old_fp), str(new_fp), "SKIP", "Source file does not exist"
                ])
                continue

            if old_fp == new_fp:
                apply_rows.append([
                    individual, sra, sample_name, animal_id,
                    str(old_fp), str(new_fp), "SKIP", "Filename does not require renaming"
                ])
                continue

            if new_fp.exists():
                apply_rows.append([
                    individual, sra, sample_name, animal_id,
                    str(old_fp), str(new_fp), "SKIP", "Target file already exists; operation skipped"
                ])
                continue

            try:
                old_fp.rename(new_fp)
                apply_rows.append([
                    individual, sra, sample_name, animal_id,
                    str(old_fp), str(new_fp), "RENAMED", ""
                ])
            except Exception as e:
                apply_rows.append([
                    individual, sra, sample_name, animal_id,
                    str(old_fp), str(new_fp), "ERROR", str(e)
                ])


            for ext in [".tbi", ".csi"]:
                old_idx = Path(str(old_fp) + ext)
                new_idx = Path(str(new_fp) + ext)
                if old_idx.exists():
                    if new_idx.exists():
                        apply_rows.append([
                            individual, sra, sample_name, animal_id,
                            str(old_idx), str(new_idx), "SKIP", f"Target index file already exists; operation skipped"
                        ])
                    else:
                        try:
                            old_idx.rename(new_idx)
                            apply_rows.append([
                                individual, sra, sample_name, animal_id,
                                str(old_idx), str(new_idx), "RENAMED", "index"
                            ])
                        except Exception as e:
                            apply_rows.append([
                                individual, sra, sample_name, animal_id,
                                str(old_idx), str(new_idx), "ERROR", f"index: {e}"
                            ])

    preview_tsv = OUT_DIR / "wgs_vcf_adult_rename.preview.tsv"
    write_tsv(preview_tsv, preview_rows)
    print(f"[OK] preview file generated: {preview_tsv}")

    if args.apply:
        apply_tsv = OUT_DIR / "wgs_vcf_adult_rename.apply.tsv"
        write_tsv(apply_tsv, apply_rows)
        print(f"[OK] execution result file generated: {apply_tsv}")
        renamed_n = sum(1 for r in apply_rows if r[6] == "RENAMED")
        print(f"[OK] actual completed of rename record count: {renamed_n}")
    else:
        print("[INFO] preview mode is active, no files were modified.")
        print("[INFO] after verifying correctness, add --apply apply the operations.")

if __name__ == "__main__":
    main()
