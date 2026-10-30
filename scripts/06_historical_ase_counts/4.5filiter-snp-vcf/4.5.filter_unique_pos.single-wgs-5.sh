#!/bin/bash
set -euo pipefail


IN_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.2WGS_ti_snp"
OUT_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.5filiter-snp-vcf"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <animal_id>" >&2
  exit 1
fi

ANIMAL_ID="$1"

IN_DIR="${IN_BASE}/${ANIMAL_ID}"
OUT_DIR="${OUT_BASE}/${ANIMAL_ID}"
mkdir -p "$OUT_DIR"

IN_VCF="${IN_DIR}/${ANIMAL_ID}.minimal.het.snp.vcf.gz"
IN_TBI="${IN_VCF}.tbi"

OUT_VCF="${OUT_DIR}/${ANIMAL_ID}.minimal.het.snp.unique_pos_only.vcf.gz"
OUT_TBI="${OUT_VCF}.tbi"
DUP_POS_TSV="${OUT_DIR}/${ANIMAL_ID}.minimal.het.snp.dup_positions.tsv"
COUNT_TXT="${OUT_DIR}/${ANIMAL_ID}.minimal.het.snp.unique_pos_only.counts.txt"

if ! command -v bcftools >/dev/null 2>&1; then
  echo "[ERROR] bcftools not found in PATH" >&2
  exit 1
fi

if ! command -v bgzip >/dev/null 2>&1; then
  echo "[ERROR] bgzip not found in PATH" >&2
  exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
  echo "[ERROR] tabix not found in PATH" >&2
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

if [[ -s "$OUT_VCF" && -s "$OUT_TBI" && -s "$DUP_POS_TSV" && -s "$COUNT_TXT" ]]; then
  echo "[skip] ${ANIMAL_ID} already finished."
  exit 0
fi

TMP_DIR="$(mktemp -d "${OUT_DIR}/tmp.unique_pos_filter.${ANIMAL_ID}.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

HEADER_VCF="${TMP_DIR}/header.vcf"
BODY_VCF="${TMP_DIR}/body.vcf"
BODY_FILTERED_VCF="${TMP_DIR}/body.filtered.vcf"
POS_ALL="${TMP_DIR}/all.pos.tsv"
DUP_KEYS="${TMP_DIR}/dup.keys.tsv"

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[info] IN_VCF=$IN_VCF"
echo "[info] OUT_VCF=$OUT_VCF"

echo "[step] collect all CHROM+POS"
bcftools query -f '%CHROM\t%POS\n' "$IN_VCF" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  > "$POS_ALL"

echo "[step] find duplicated CHROM+POS"
awk '
BEGIN{
  FS=OFS="\t";
  print "CHROM","POS","N_RECORDS";
}
NR==1{
  prev_chr=$1;
  prev_pos=$2;
  cnt=1;
  next;
}
{
  if($1==prev_chr && $2==prev_pos){
    cnt++;
  } else {
    if(cnt>1){
      print prev_chr, prev_pos, cnt;
    }
    prev_chr=$1;
    prev_pos=$2;
    cnt=1;
  }
}
END{
  if(NR>0 && cnt>1){
    print prev_chr, prev_pos, cnt;
  }
}
' "$POS_ALL" > "$DUP_POS_TSV"

tail -n +2 "$DUP_POS_TSV" | cut -f1,2 > "$DUP_KEYS" || true

echo "[step] extract header and body"
bcftools view -h "$IN_VCF" > "$HEADER_VCF"
bcftools view -H "$IN_VCF" > "$BODY_VCF"

echo "[step] remove all duplicated-position records"
if [[ -s "$DUP_KEYS" ]]; then
  awk '
  BEGIN{FS=OFS="\t"}
  NR==FNR{
    bad[$1 FS $2]=1;
    next;
  }
  !(($1 FS $2) in bad){
    print
  }
  ' "$DUP_KEYS" "$BODY_VCF" > "$BODY_FILTERED_VCF"
else
  cp "$BODY_VCF" "$BODY_FILTERED_VCF"
fi

echo "[step] write output VCF"
cat "$HEADER_VCF" "$BODY_FILTERED_VCF" | bgzip -c > "$OUT_VCF"
tabix -f -p vcf "$OUT_VCF"

echo "[step] summarize counts"
INPUT_TOTAL=$(bcftools view -H "$IN_VCF" | wc -l)
OUTPUT_TOTAL=$(bcftools view -H "$OUT_VCF" | wc -l)
DUP_POS_N=$(tail -n +2 "$DUP_POS_TSV" | wc -l)
REMOVED_TOTAL=$((INPUT_TOTAL - OUTPUT_TOTAL))

REMAIN_DUP_POS_N=$(
  bcftools query -f '%CHROM\t%POS\n' "$OUT_VCF" \
  | awk '
    {k=$1 FS $2; n[k]++}
    END{
      c=0;
      for(k in n){
        if(n[k]>1) c++;
      }
      print c;
    }
  '
)

{
  echo "input_vcf=${IN_VCF}"
  echo "output_vcf=${OUT_VCF}"
  echo "input_total=${INPUT_TOTAL}"
  echo "output_total=${OUTPUT_TOTAL}"
  echo "duplicated_position_count=${DUP_POS_N}"
  echo "removed_record_total=${REMOVED_TOTAL}"
  echo "remaining_duplicated_position_count=${REMAIN_DUP_POS_N}"
} > "$COUNT_TXT"

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[done] OUT_VCF=$OUT_VCF"
echo "[done] DUP_POS_TSV=$DUP_POS_TSV"
echo "[done] COUNT_TXT=$COUNT_TXT"
