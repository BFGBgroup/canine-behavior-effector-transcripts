#!/usr/bin/bash
set -euo pipefail


PYTHON_BIN="/usr/bin/python3"
SCRIPT="/path/to/repository/scripts/03_candidate_catalog/2.3_candidate_TES_polyA_support/2.3_candidate_TES_polyA_support.py"
PROJECT_ROOT="/path/to/project/analysis_20260709/5-dogs/9.new-gene"
STEP21_ROOT="${PROJECT_ROOT}/2.1_transcript_first_preASE_candidate_assessment"
OUTPUT_ROOT="${PROJECT_ROOT}/2.3_candidate_TES_polyA_support"
LOG_DIR="/path/to/repository/scripts/03_candidate_catalog/2.3_candidate_TES_polyA_support/log"

RETAINED_TSV="${STEP21_ROOT}/2.1_retained_transcript_candidates.preASE.tsv"
SUMMARY_TSV="${STEP21_ROOT}/2.1_candidate_transcript_summary.preASE.tsv"
EXON_BED="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/transcript_exon.bed"
GENOME_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
CHROM_SIZES="${GENOME_FASTA}.fai"
TES_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/reference/TES"
ANIMAL_ROOT="${TES_ROOT}/1.Animal-APAdb-dog-PAS/results-v2"
GSE_ROOT="${TES_ROOT}/2.GEO-GSE30198-PolyA-Seq/results-v1"
RUN_ID="$(/usr/bin/date '+%Y%m%d_%H%M%S')"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/${RUN_ID}.master.log"
SUCCESS_FILE="${LOG_DIR}/${RUN_ID}.success.tsv"
FAILED_FILE="${LOG_DIR}/${RUN_ID}.failed.tsv"
/usr/bin/mkdir -p "$LOG_DIR"

(set -o noclobber; : > "$MASTER_LOG") || exit 1
exec >> "$MASTER_LOG" 2>&1
printf 'task\tstatus\toutput_dir\n' > "$SUCCESS_FILE"
printf 'task\tstatus\texit_code\toutput_dir\n' > "$FAILED_FILE"
finish() {
  status=$?
  if (( status != 0 )); then
    printf 'TES_polyA\tfailed\t%s\t%s\n' "$status" "$RUN_DIR" >> "$FAILED_FILE"
    printf '[%s ] total=1 success=0 failed=1 skipped=0 execution_status=failed exit_code=%s\n' "$(/usr/bin/date '+%F %T')" "$status"
  fi
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
printf '[%s ] START concurrent_tasks=1 threads_per_task=1(serial settings retained) output=%s\n' "$(/usr/bin/date '+%F %T')" "$RUN_DIR"
printf 'close_support_bp=50 sensitivity_bp=100 primary_relaxed_proximity_bp=500\nsource_evidence=PolyA-seq_experimental_reference_and_Animal-APAdb_inferred_reference_reported_separately\nmotif_window=genomic_-50_to_-10_selected9_and_spliced_terminal50nt_selected9_extended11\n'

for required_file in "$PYTHON_BIN" "$SCRIPT" "$RETAINED_TSV" "$SUMMARY_TSV" "$EXON_BED" "$GENOME_FASTA" "$CHROM_SIZES" \
  "${ANIMAL_ROOT}/animal_apadb_dog_PAS_all.CanFam3.bed" \
  "${ANIMAL_ROOT}/animal_apadb_dog_PAS_brain_CNS.CanFam3.bed" \
  "${ANIMAL_ROOT}/animal_apadb_dog_PAS_brain_CNS_plus_pituitary.CanFam3.bed" \
  "${GSE_ROOT}/GSE30198_dog_brain_polyA_sites.CanFam3.bed" \
  "${GSE_ROOT}/GSE30198_dog_polyA_sites.CanFam3.bed" \
  "${GSE_ROOT}/GSE30198_dog_nonbrain_polyA_sites.CanFam3.bed"; do
  if [[ ! -s "$required_file" ]]; then
    echo "[ERROR] required input missing or empty: $required_file" >&2
    exit 1
  fi
done

"$PYTHON_BIN" "$SCRIPT" \
  --retained-transcript-tsv "$RETAINED_TSV" \
  --transcript-summary-tsv "$SUMMARY_TSV" \
  --transcript-exon-bed "$EXON_BED" \
  --chrom-sizes "$CHROM_SIZES" \
  --genome-fasta "$GENOME_FASTA" \
  --animal-apadb-all-bed "${ANIMAL_ROOT}/animal_apadb_dog_PAS_all.CanFam3.bed" \
  --animal-apadb-brain-cns-bed "${ANIMAL_ROOT}/animal_apadb_dog_PAS_brain_CNS.CanFam3.bed" \
  --animal-apadb-brain-cns-plus-pituitary-bed "${ANIMAL_ROOT}/animal_apadb_dog_PAS_brain_CNS_plus_pituitary.CanFam3.bed" \
  --gse30198-brain-bed "${GSE_ROOT}/GSE30198_dog_brain_polyA_sites.CanFam3.bed" \
  --gse30198-all-bed "${GSE_ROOT}/GSE30198_dog_polyA_sites.CanFam3.bed" \
  --gse30198-nonbrain-bed "${GSE_ROOT}/GSE30198_dog_nonbrain_polyA_sites.CanFam3.bed" \
  --outdir "$RUN_DIR"
printf 'execution_status=success\n' > "${RUN_DIR}/execution_completion.marker"
printf 'TES_polyA\tsuccess\t%s\n' "$RUN_DIR" >> "$SUCCESS_FILE"
printf '[%s ] total=1 success=1 failed=0 skipped=0\n' "$(/usr/bin/date '+%F %T')"
printf '[%s ] task completed / done execution_status=success\n' "$(/usr/bin/date '+%F %T')"
printf 'result completeness, format, statistical denominators and biological quality require separate review; successful execution does not establish output quality.\n'
