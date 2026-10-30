#!/usr/bin/env bash
set -euo pipefail


if [[ $# -lt 2 ]]; then
    echo "usage: bash $0 <animal_id> <term>"
    exit 1
fi

ANIMAL_ID="$1"
TERM="$2"

THREADS=10
FORCE_REBUILD=0

RAW_BASE="/path/to/project/legacy/five-dogs/3.bulk-rna-qc/3.1qc"
OUT_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.1star-index"

GENOME_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
GTF_FILE="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf"

RAW_DIR="${RAW_BASE}/${ANIMAL_ID}/${TERM}/01.clean_fastp"
FQ1_PATH="${RAW_DIR}/${ANIMAL_ID}_${TERM}_1.clean.fastq.gz"
FQ2_PATH="${RAW_DIR}/${ANIMAL_ID}_${TERM}_2.clean.fastq.gz"

[[ -d "${RAW_BASE}" ]]     || { echo "[ERROR] cannot find RAW_BASE: ${RAW_BASE}"; exit 1; }
[[ -f "${GENOME_FASTA}" ]] || { echo "[ERROR] cannot find GENOME_FASTA: ${GENOME_FASTA}"; exit 1; }
[[ -f "${GTF_FILE}" ]]     || { echo "[ERROR] cannot find GTF_FILE: ${GTF_FILE}"; exit 1; }
[[ -d "${RAW_DIR}" ]]      || { echo "[ERROR] cannot find original directory: ${RAW_DIR}"; exit 1; }
[[ -f "${FQ1_PATH}" ]]     || { echo "[ERROR] cannot find R1: ${FQ1_PATH}"; exit 1; }
[[ -f "${FQ2_PATH}" ]]     || { echo "[ERROR] cannot find R2: ${FQ2_PATH}"; exit 1; }

SAMPLE_DIR="${OUT_BASE}/${ANIMAL_ID}/${TERM}"
STAR_INDEX_DIR="${SAMPLE_DIR}/STARindex"
STAT_FILE="${SAMPLE_DIR}/${ANIMAL_ID}_${TERM}.readlen.tsv"
LOG_FILE="${SAMPLE_DIR}/${ANIMAL_ID}_${TERM}.STAR_genomeGenerate.log"

mkdir -p "${SAMPLE_DIR}"

check_one_fastq() {
    local fq="$1"
    zcat "$fq" | awk '
    NR%4==2{
        l=length($0)
        n++
        cnt[l]++
        if(min==0 || l<min) min=l
        if(l>max) max=l
    }
    END{
        uniq=0
        for(i in cnt) uniq++
        printf "%d\t%d\t%d\t%d\n", n, min, max, uniq
    }'
}

R1_STAT="$(check_one_fastq "${FQ1_PATH}")"
R2_STAT="$(check_one_fastq "${FQ2_PATH}")"

READS1="$(echo "${R1_STAT}" | cut -f1)"
MIN1="$(echo "${R1_STAT}"   | cut -f2)"
MAX1="$(echo "${R1_STAT}"   | cut -f3)"
UNIQ1="$(echo "${R1_STAT}"  | cut -f4)"

READS2="$(echo "${R2_STAT}" | cut -f1)"
MIN2="$(echo "${R2_STAT}"   | cut -f2)"
MAX2="$(echo "${R2_STAT}"   | cut -f3)"
UNIQ2="$(echo "${R2_STAT}"  | cut -f4)"

MAX_LEN="${MAX1}"
if (( MAX2 > MAX_LEN )); then
    MAX_LEN="${MAX2}"
fi

SJDB_OVERHANG=$((MAX_LEN - 1))

echo "=============================="
echo "[INFO] build per-sample STAR index"
echo "[INFO] ANIMAL_ID       : ${ANIMAL_ID}"
echo "[INFO] TERM            : ${TERM}"
echo "[INFO] RAW_DIR         : ${RAW_DIR}"
echo "[INFO] FQ1             : ${FQ1_PATH}"
echo "[INFO] FQ2             : ${FQ2_PATH}"
echo "[INFO] READS1          : ${READS1}"
echo "[INFO] READS2          : ${READS2}"
echo "[INFO] MIN1/MAX1/UNIQ1 : ${MIN1}/${MAX1}/${UNIQ1}"
echo "[INFO] MIN2/MAX2/UNIQ2 : ${MIN2}/${MAX2}/${UNIQ2}"
echo "[INFO] MAX_READ_LENGTH : ${MAX_LEN}"
echo "[INFO] sjdbOverhang    : ${SJDB_OVERHANG}"
echo "[INFO] SAMPLE_DIR      : ${SAMPLE_DIR}"
echo "[INFO] INDEX_DIR       : ${STAR_INDEX_DIR}"
echo "[INFO] LOG_FILE        : ${LOG_FILE}"
echo "=============================="

cat > "${STAT_FILE}" <<EOF
animal_id	term	raw_dir	fq1_path	fq2_path	reads1	reads2	min_len1	max_len1	uniq_len_count1	min_len2	max_len2	uniq_len_count2	max_read_length	sjdbOverhang
${ANIMAL_ID}	${TERM}	${RAW_DIR}	${FQ1_PATH}	${FQ2_PATH}	${READS1}	${READS2}	${MIN1}	${MAX1}	${UNIQ1}	${MIN2}	${MAX2}	${UNIQ2}	${MAX_LEN}	${SJDB_OVERHANG}
EOF

if [[ -s "${STAR_INDEX_DIR}/SA" && "${FORCE_REBUILD}" -eq 0 ]]; then
    echo "[SKIP] Index already exists; skipping: ${STAR_INDEX_DIR}"
    exit 0
fi

rm -rf "${STAR_INDEX_DIR}"
mkdir -p "${STAR_INDEX_DIR}"

STAR \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF_FILE}" \
  --sjdbOverhang "${SJDB_OVERHANG}" \
  > "${LOG_FILE}" 2>&1

echo "[OK] ${ANIMAL_ID} ${TERM} STAR genomeGenerate completed"
echo "[OK] statistics file: ${STAT_FILE}"
echo "[OK] log file: ${LOG_FILE}"
