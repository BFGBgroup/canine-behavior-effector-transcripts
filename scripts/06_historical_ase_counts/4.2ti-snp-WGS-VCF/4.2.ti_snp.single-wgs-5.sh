#!/bin/bash
set -euo pipefail


REF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"

VCF_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_gatk_norm"
OUT_ROOT="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.2WGS_ti_snp"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <animal_id>" >&2
  exit 1
fi

ANIMAL_ID="$1"

IN_VCF="${VCF_ROOT}/${ANIMAL_ID}/${ANIMAL_ID}.filtered.PASS.snp.norm.vcf.gz"
IN_TBI="${IN_VCF}.tbi"

SAMPLE_OUT_DIR="${OUT_ROOT}/${ANIMAL_ID}"
mkdir -p "$SAMPLE_OUT_DIR"

OUT_VCF="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.minimal.het.snp.vcf.gz"
OUT_TBI="${OUT_VCF}.tbi"
OUT_VCF_TXT="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.minimal.het.snp.vcf"
COUNT_TXT="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.minimal.het.snp.counts.txt"

if ! command -v bcftools >/dev/null 2>&1; then
  echo "[ERROR] bcftools not found in PATH" >&2
  exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
  echo "[ERROR] tabix not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$REF" ]]; then
  echo "[ERROR] REF not found: $REF" >&2
  exit 1
fi

if [[ ! -f "${REF}.fai" ]]; then
  echo "[ERROR] FASTA index not found: ${REF}.fai" >&2
  exit 1
fi

if [[ ! -f "$IN_VCF" ]]; then
  echo "[ERROR] input VCF not found: $IN_VCF" >&2
  exit 1
fi

if [[ ! -f "$IN_TBI" ]]; then
  echo "[ERROR] input VCF index not found: $IN_TBI" >&2
  exit 1
fi

if [[ -s "$OUT_VCF" && -s "$OUT_TBI" && -s "$OUT_VCF_TXT" && -s "$COUNT_TXT" ]]; then
  echo "[skip] ${ANIMAL_ID} already finished."
  exit 0
fi

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[info] REF=$REF"
echo "[info] IN_VCF=$IN_VCF"
echo "[info] OUTDIR=$SAMPLE_OUT_DIR"

echo "[step] extract minimal het biallelic SNP VCF"
bcftools view \
  -m2 -M2 \
  -v snps \
  -i 'GT="het" && N_ALT=1' \
  -Oz \
  -o "$OUT_VCF" \
  "$IN_VCF"

tabix -f -p vcf "$OUT_VCF"

echo "[step] export plain vcf"
gunzip -c "$OUT_VCF" > "$OUT_VCF_TXT"

echo "[step] summarize counts"
{
  echo -n "input_total="
  bcftools view -H "$IN_VCF" | wc -l
  echo -n "minimal_het_total="
  bcftools view -H "$OUT_VCF" | wc -l
} > "$COUNT_TXT"

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[done] OUT_VCF=$OUT_VCF"
echo "[done] OUT_VCF_TXT=$OUT_VCF_TXT"
echo "[done] COUNT_TXT=$COUNT_TXT"
