#!/bin/bash
set -euo pipefail


THREADS="${THREADS:-12}"
LIMIT_BAM_SORT_RAM="${LIMIT_BAM_SORT_RAM:-120000000000}"
FORCE_RERUN="${FORCE_RERUN:-0}"

RAW_BASE="/path/to/project/legacy/five-dogs/3.bulk-rna-qc/3.1qc"
STAR_INDEX_BASE="/path/to/project/legacy/five-dogs/5.cnv/1.1star-index"
OUT_BASE="/path/to/project/legacy/five-dogs/5.cnv/1.star-v5"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <animal_id> <term>" >&2
  exit 1
fi

ANIMAL_ID="$1"
TERM="$2"

R1="${RAW_BASE}/${ANIMAL_ID}/${TERM}/01.clean_fastp/${ANIMAL_ID}_${TERM}_1.clean.fastq.gz"
R2="${RAW_BASE}/${ANIMAL_ID}/${TERM}/01.clean_fastp/${ANIMAL_ID}_${TERM}_2.clean.fastq.gz"

STAR_INDEX_DIR="${STAR_INDEX_BASE}/${ANIMAL_ID}/${TERM}/STARindex"

SAMPLE_OUT_DIR="${OUT_BASE}/${ANIMAL_ID}/${TERM}/align_star"
mkdir -p "${SAMPLE_OUT_DIR}"

OUT_PREFIX="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_${TERM}."
OUT_BAM="${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
OUT_BAI="${OUT_BAM}.bai"
OUT_COUNTS="${OUT_PREFIX}ReadsPerGene.out.tab"
OUT_SJ="${OUT_PREFIX}SJ.out.tab"
STAR_STD_LOG="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_${TERM}.STAR.stdout.log"
STAR_MAIN_LOG="${OUT_PREFIX}Log.out"
STAR_FINAL_LOG="${OUT_PREFIX}Log.final.out"
STAR_PROGRESS_LOG="${OUT_PREFIX}Log.progress.out"

RG_ID="${ANIMAL_ID}_${TERM}"
RG_SM="${ANIMAL_ID}_${TERM}"
RG_LB="${ANIMAL_ID}_${TERM}"
RG_PL="ILLUMINA"
RG_PU="${ANIMAL_ID}_${TERM}"

if ! command -v STAR >/dev/null 2>&1; then
  echo "[ERROR] STAR not found in PATH" >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools not found in PATH" >&2
  exit 1
fi

if [[ ! -f "${R1}" ]]; then
  echo "[ERROR] R1 not found: ${R1}" >&2
  exit 1
fi

if [[ ! -f "${R2}" ]]; then
  echo "[ERROR] R2 not found: ${R2}" >&2
  exit 1
fi

if [[ ! -d "${STAR_INDEX_DIR}" || ! -s "${STAR_INDEX_DIR}/SA" ]]; then
  echo "[ERROR] STAR index not found or incomplete: ${STAR_INDEX_DIR}" >&2
  echo "[ERROR] Please build index first with scripts/1.1star-index." >&2
  exit 1
fi

if [[ "${FORCE_RERUN}" != "1" && -s "${OUT_BAM}" && -s "${OUT_BAI}" && -s "${OUT_COUNTS}" ]]; then
  echo "[skip] ${ANIMAL_ID} ${TERM} already finished. Set FORCE_RERUN=1 to rebuild with current STAR parameters."
  exit 0
fi

if [[ "${FORCE_RERUN}" == "1" ]]; then
  echo "[force] FORCE_RERUN=1, removing old STAR outputs under the same output prefix before rebuilding."
  rm -f \
    "${OUT_BAM}" \
    "${OUT_BAI}" \
    "${OUT_COUNTS}" \
    "${OUT_SJ}" \
    "${STAR_STD_LOG}" \
    "${STAR_MAIN_LOG}" \
    "${STAR_FINAL_LOG}" \
    "${STAR_PROGRESS_LOG}"
fi

echo "[start] $(date '+%F %T') animal_id=${ANIMAL_ID} term=${TERM}"
echo "[info] THREADS=${THREADS}"
echo "[info] LIMIT_BAM_SORT_RAM=${LIMIT_BAM_SORT_RAM}"
echo "[info] FORCE_RERUN=${FORCE_RERUN}"
echo "[info] R1=${R1}"
echo "[info] R2=${R2}"
echo "[info] STAR_INDEX_DIR=${STAR_INDEX_DIR}"
echo "[info] OUT_PREFIX=${OUT_PREFIX}"
echo "[info] STAR strand field=--outSAMstrandField intronMotif"
echo "[info] RG_ID=${RG_ID}"
echo "[info] RG_SM=${RG_SM}"
echo "[info] RG_LB=${RG_LB}"
echo "[info] RG_PL=${RG_PL}"
echo "[info] RG_PU=${RG_PU}"

ulimit -n 65535

STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX_DIR}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_PREFIX}" \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM "${LIMIT_BAM_SORT_RAM}" \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM \
  --outSAMattrRGline ID:"${RG_ID}" SM:"${RG_SM}" LB:"${RG_LB}" PL:"${RG_PL}" PU:"${RG_PU}" \
  --quantMode GeneCounts \
  > "${STAR_STD_LOG}" 2>&1

samtools index -@ "${THREADS}" "${OUT_BAM}"

echo "[done] $(date '+%F %T') animal_id=${ANIMAL_ID} term=${TERM}"
echo "[done] BAM=${OUT_BAM}"
echo "[done] BAM index=${OUT_BAI}"
echo "[done] Counts=${OUT_COUNTS}"
echo "[done] STAR stdout/stderr log=${STAR_STD_LOG}"
echo "[done] STAR main log=${OUT_PREFIX}Log.out"
echo "[done] STAR final log=${OUT_PREFIX}Log.final.out"
