#!/bin/bash
set -euo pipefail


REF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"

VCF_ROOT="/path/to/project/doga_six_dogs/2.wgs-fastq/2.3WGS_gatk"
OUT_ROOT="/path/to/project/doga_six_dogs/2.wgs-fastq/2.4WGS_gatk_norm"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <animal_id>" >&2
  exit 1
fi

ANIMAL_ID="$1"

IN_VCF="${VCF_ROOT}/${ANIMAL_ID}/${ANIMAL_ID}.filtered.snp.vcf.gz"
IN_TBI="${IN_VCF}.tbi"

SAMPLE_OUT_DIR="${OUT_ROOT}/${ANIMAL_ID}"
mkdir -p "$SAMPLE_OUT_DIR"

NORM_VCF="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.filtered.snp.norm.vcf.gz"
NORM_TBI="${NORM_VCF}.tbi"

PASS_NORM_VCF="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.filtered.PASS.snp.norm.vcf.gz"
PASS_NORM_TBI="${PASS_NORM_VCF}.tbi"

COUNT_TXT="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.norm.counts.txt"

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

if [[ -s "$NORM_VCF" && -s "$NORM_TBI" \
   && -s "$PASS_NORM_VCF" && -s "$PASS_NORM_TBI" \
   && -s "$COUNT_TXT" ]]; then
  echo "[skip] ${ANIMAL_ID} already finished."
  exit 0
fi

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[info] REF=$REF"
echo "[info] IN_VCF=$IN_VCF"
echo "[info] OUTDIR=$SAMPLE_OUT_DIR"

echo "[step] normalize filtered SNP VCF"
bcftools norm \
  -f "$REF" \
  -m -both \
  -Oz \
  -o "$NORM_VCF" \
  "$IN_VCF"

tabix -f -p vcf "$NORM_VCF"

echo "[step] extract PASS SNP from normalized VCF"
bcftools view \
  -f PASS \
  -v snps \
  -Oz \
  -o "$PASS_NORM_VCF" \
  "$NORM_VCF"

tabix -f -p vcf "$PASS_NORM_VCF"

echo "[step] summarize counts"
{
  echo -n "input_total="
  bcftools view -H "$IN_VCF" | wc -l
  echo -n "norm_total="
  bcftools view -H "$NORM_VCF" | wc -l
  echo -n "norm_pass_total="
  bcftools view -H "$PASS_NORM_VCF" | wc -l
} > "$COUNT_TXT"

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[done] NORM_VCF=$NORM_VCF"
echo "[done] PASS_NORM_VCF=$PASS_NORM_VCF"
echo "[done] COUNT_TXT=$COUNT_TXT"
