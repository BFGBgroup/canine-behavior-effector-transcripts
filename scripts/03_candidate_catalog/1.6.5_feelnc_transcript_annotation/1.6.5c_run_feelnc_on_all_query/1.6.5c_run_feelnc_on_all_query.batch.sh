#!/bin/bash
set -euo pipefail


SCRIPT="/path/to/repository/scripts/03_candidate_catalog/1.6.5_feelnc_transcript_annotation/1.6.5c_run_feelnc_on_all_query/1.6.5c_collect_feelnc_module_outputs.py"
PYTHON_BIN="/usr/bin/python3"

QUERY_GTF="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5b_prepare_query_transcripts/query_transcripts/stringtie3_gffcompare_all_query_transcripts.gtf"
QUERY_BASIC_INFO="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5b_prepare_query_transcripts/query_transcripts/stringtie3_gffcompare_all_query_transcript_basic_info.tsv"
QUERY_QC="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5b_prepare_query_transcripts/query_transcripts/1.6.5b_prepare_query_transcripts.qc_summary.tsv"
CODPOT_TRAINING_GTF="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5a/training_sets/CanFam3.1_refseq_protein_coding_transcripts.from_patched_CDS.gtf"
CODPOT_TRAINING_QC="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5a/training_sets/CanFam3.1_refseq_protein_coding_training_from_patched_CDS.qc_summary.tsv"
REFERENCE_PROTEIN_CODING_GTF="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5a/training_sets/CanFam3.1_refseq_protein_coding_transcripts.from_patched_CDS.gtf"
GENOME_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"

OUTPUT_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5c_run_feelnc_on_all_query-v4"
FILTER_DIR="${OUTPUT_ROOT}/filter"
CODPOT_DIR="${OUTPUT_ROOT}/codpot"
CLASSIFIER_DIR="${OUTPUT_ROOT}/classifier"
QC_DIR="${OUTPUT_ROOT}/qc"
LOG_DIR="${OUTPUT_ROOT}/logs-v1"

CODPOT_TRAINING_CLEAN_GTF="${CODPOT_DIR}/CanFam3.1_refseq_protein_coding_transcripts.from_patched_CDS.no_chrM_for_codpot.gtf"
FILTER_GTF="${FILTER_DIR}/all_query.feelnc_filter_pass.gtf"
FILTER_LOG="${FILTER_DIR}/all_query.feelnc_filter.log"
CODPOT_OUTNAME="all_query.feelnc_codpot"
CODPOT_RF="${CODPOT_DIR}/${CODPOT_OUTNAME}_RF.txt"
CODPOT_RF_SUMMARY="${CODPOT_DIR}/${CODPOT_OUTNAME}_RF_summary.txt"
CODPOT_LNCRNA_GTF="${CODPOT_DIR}/${CODPOT_OUTNAME}.lncRNA.gtf"
CODPOT_MRNA_GTF="${CODPOT_DIR}/${CODPOT_OUTNAME}.mRNA.gtf"
CODPOT_TUCP_GTF="${CODPOT_DIR}/${CODPOT_OUTNAME}.TUCp.gtf"
CODPOT_NOORF_GTF="${CODPOT_DIR}/${CODPOT_OUTNAME}.noORF.gtf"
CLASSIFIER_INPUT_GTF="${CLASSIFIER_DIR}/all_query.feelnc_classifier_input_lncRNA_noORF.gtf"
CLASSIFIER_CLASSES="${CLASSIFIER_DIR}/all_query.feelnc_classifier.classes.txt"
CLASSIFIER_LOG="${CLASSIFIER_DIR}/all_query.feelnc_classifier.log"
STATUS_TSV="${QC_DIR}/1.6.5c_all_query_feelnc_module_status.tsv"
SUMMARY_TSV="${QC_DIR}/1.6.5c_feelnc_module_summary.tsv"
MANIFEST="${QC_DIR}/1.6.5c_feelnc_module_manifest.tsv"

CONDA_ENV_PATH="/path/to/software/environments/stringtie3"
CONDA_SH="/path/to/software/legacy_miniconda3/etc/profile.d/conda.sh"
ENV_BIN="${CONDA_ENV_PATH}/bin"
FEELNC_HOME="${CONDA_ENV_PATH}/share/FEELnc"

THREADS=50
FILTER_MONOEX="-1"
FILTER_MIN_SIZE=200
CODPOT_MODE="shuffle"
CODPOT_SPETHRES="0.93,0.93"
CODPOT_SEED=1234
CLASSIFIER_WINDOW=10000
CLASSIFIER_MAXWINDOW=100000
OVERWRITE=0

if [[ "$#" -ne 0 ]]; then
  echo "[ERROR] This batch does not accept positional arguments." >&2
  echo "Usage: bash $0" >&2
  exit 2
fi

if [[ "$OVERWRITE" != "0" && "$OVERWRITE" != "1" ]]; then
  echo "[ERROR] OVERWRITE must be 0 or 1, got: $OVERWRITE" >&2
  exit 2
fi

if [[ -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
  conda activate "$CONDA_ENV_PATH"
else
  echo "[WARN] conda profile script not found: $CONDA_SH; continue with explicit ENV_BIN/FEELNC_HOME PATH" >&2
fi

export FEELNCPATH="$FEELNC_HOME"
export PERL5LIB="${FEELNC_HOME}/lib:${PERL5LIB:-}"
export PATH="${ENV_BIN}:${FEELNC_HOME}/scripts:${FEELNC_HOME}/utils:${FEELNC_HOME}/bin/LINUX:${PATH}"

for path in "$SCRIPT" "$QUERY_GTF" "$QUERY_BASIC_INFO" "$QUERY_QC" "$CODPOT_TRAINING_GTF" "$CODPOT_TRAINING_QC" "$REFERENCE_PROTEIN_CODING_GTF" "$GENOME_FASTA"; do
  if [[ ! -s "$path" ]]; then
    echo "[ERROR] required input missing or empty: $path" >&2
    exit 1
  fi
done

for tool in perl Rscript FEELnc_filter.pl FEELnc_codpot.pl FEELnc_classifier.pl fasta_ushuffle KmerInShort; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] required executable not found in PATH: $tool" >&2
    exit 1
  fi
done

QUERY_PASS_STATUS=$(awk -F '\t' '$1=="pass_status"{print $2}' "$QUERY_QC")
TRAINING_PASS_STATUS=$(awk -F '\t' '$1=="pass_status"{print $2}' "$CODPOT_TRAINING_QC")
if [[ "$QUERY_PASS_STATUS" != "PASS" ]]; then
  echo "[ERROR] 1.6.5b query QC is not PASS: $QUERY_QC" >&2
  exit 1
fi
if [[ "$TRAINING_PASS_STATUS" != "PASS" ]]; then
  echo "[ERROR] 1.6.5a training QC is not PASS: $CODPOT_TRAINING_QC" >&2
  exit 1
fi

EXPECTED_OUTPUTS=(
  "$FILTER_GTF"
  "$FILTER_LOG"
  "$CODPOT_TRAINING_CLEAN_GTF"
  "$CODPOT_RF"
  "$CODPOT_RF_SUMMARY"
  "$CODPOT_LNCRNA_GTF"
  "$CODPOT_MRNA_GTF"
  "$CODPOT_TUCP_GTF"
  "$CODPOT_NOORF_GTF"
  "$CLASSIFIER_INPUT_GTF"
  "$CLASSIFIER_CLASSES"
  "$CLASSIFIER_LOG"
  "$STATUS_TSV"
  "$SUMMARY_TSV"
  "$MANIFEST"
)
if [[ "$OVERWRITE" == "0" ]]; then
  for path in "${EXPECTED_OUTPUTS[@]}"; do
    if [[ -e "$path" ]]; then
      echo "[ERROR] output already exists; set OVERWRITE=1 only after confirmation: $path" >&2
      exit 1
    fi
  done
fi

mkdir -p "$FILTER_DIR" "$CODPOT_DIR" "$CLASSIFIER_DIR" "$QC_DIR" "$LOG_DIR"

MASTER_LOG="${LOG_DIR}/1.6.5c_run_feelnc_on_all_query.batch.master.log"
FILTER_STDERR="${LOG_DIR}/1.6.5c_feelnc_filter.stderr.log"
CODPOT_STDERR="${LOG_DIR}/1.6.5c_feelnc_codpot.stderr.log"
CLASSIFIER_STDERR="${LOG_DIR}/1.6.5c_feelnc_classifier.stderr.log"
COLLECT_STDERR="${LOG_DIR}/1.6.5c_collect_feelnc_module_outputs.stderr.log"
SUCCESS_LIST="${LOG_DIR}/1.6.5c_run_feelnc_on_all_query.success.tsv"
FAILED_LIST="${LOG_DIR}/1.6.5c_run_feelnc_on_all_query.failed.tsv"

printf "task\toutdir\trc\tlog\n" > "$SUCCESS_LIST"
printf "task\toutdir\trc\tlog\n" > "$FAILED_LIST"

{
  echo "[batch] $(date '+%F %T') START"
  echo "[batch] QUERY_GTF=$QUERY_GTF"
  echo "[batch] QUERY_BASIC_INFO=$QUERY_BASIC_INFO"
  echo "[batch] CODPOT_TRAINING_GTF=$CODPOT_TRAINING_GTF"
  echo "[batch] CODPOT_TRAINING_CLEAN_GTF=$CODPOT_TRAINING_CLEAN_GTF"
  echo "[batch] REFERENCE_PROTEIN_CODING_GTF=$REFERENCE_PROTEIN_CODING_GTF"
  echo "[batch] GENOME_FASTA=$GENOME_FASTA"
  echo "[batch] OUTPUT_ROOT=$OUTPUT_ROOT"
  echo "[batch] FEELNC_HOME=$FEELNC_HOME"
  echo "[batch] PATH=$(command -v FEELnc_filter.pl);$(command -v FEELnc_codpot.pl);$(command -v FEELnc_classifier.pl);$(command -v fasta_ushuffle);$(command -v KmerInShort)"
  echo "[batch] THREADS=$THREADS"
  echo "[batch] FILTER_MONOEX=$FILTER_MONOEX"
  echo "[batch] FILTER_MIN_SIZE=$FILTER_MIN_SIZE"
  echo "[batch] CODPOT_MODE=$CODPOT_MODE"
  echo "[batch] CODPOT_SPETHRES=$CODPOT_SPETHRES"
  echo "[batch] CODPOT_SEED=$CODPOT_SEED"
  echo "[batch] CLASSIFIER_WINDOW=$CLASSIFIER_WINDOW"
  echo "[batch] CLASSIFIER_MAXWINDOW=$CLASSIFIER_MAXWINDOW"
  echo "[batch] OVERWRITE=$OVERWRITE"
} > "$MASTER_LOG"

echo "[batch] $(date '+%F %T') prepare codpot training GTF without chrM records" | tee -a "$MASTER_LOG"
awk 'BEGIN{FS=OFS="\t"} /^#/ {next} $1!="chrM" {print}' "$CODPOT_TRAINING_GTF" > "$CODPOT_TRAINING_CLEAN_GTF"
CODPOT_TRAINING_CLEAN_TX=$(awk -F '\t' '$3=="transcript"{n++} END{print n+0}' "$CODPOT_TRAINING_CLEAN_GTF")
CODPOT_TRAINING_EXCLUDED_CHRM_TX=$(awk -F '\t' '$1=="chrM" && $3=="transcript"{n++} END{print n+0}' "$CODPOT_TRAINING_GTF")
if [[ ! -s "$CODPOT_TRAINING_CLEAN_GTF" || "$CODPOT_TRAINING_CLEAN_TX" -lt 100 ]]; then
  echo "[ERROR] codpot clean training GTF is unexpectedly small: $CODPOT_TRAINING_CLEAN_GTF transcripts=$CODPOT_TRAINING_CLEAN_TX" >&2
  exit 1
fi
echo "[batch] $(date '+%F %T') codpot clean training transcripts=$CODPOT_TRAINING_CLEAN_TX excluded_chrM_transcripts=$CODPOT_TRAINING_EXCLUDED_CHRM_TX" | tee -a "$MASTER_LOG"

set +e
FEELnc_filter.pl \
  --infile="$QUERY_GTF" \
  --mRNAfile="$REFERENCE_PROTEIN_CODING_GTF" \
  --monoex="$FILTER_MONOEX" \
  --size="$FILTER_MIN_SIZE" \
  --proc="$THREADS" \
  --outlog="$FILTER_LOG" \
  > "$FILTER_GTF" 2> "$FILTER_STDERR"
RC_FILTER=$?
set -e
if [[ "$RC_FILTER" -ne 0 ]]; then
  printf "FEELnc_filter\t%s\t%s\t%s\n" "$FILTER_DIR" "$RC_FILTER" "$FILTER_STDERR" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FEELnc_filter failed rc=$RC_FILTER" | tee -a "$MASTER_LOG"
  exit "$RC_FILTER"
fi
printf "FEELnc_filter\t%s\t%s\t%s\n" "$FILTER_DIR" "$RC_FILTER" "$FILTER_STDERR" >> "$SUCCESS_LIST"
echo "[batch] $(date '+%F %T') FEELnc_filter success" | tee -a "$MASTER_LOG"

CODPOT_CMD=(
  FEELnc_codpot.pl
  --infile="$QUERY_GTF"
  --mRNAfile="$CODPOT_TRAINING_CLEAN_GTF"
  --genome="$GENOME_FASTA"
  --mode="$CODPOT_MODE"
  --outname="$CODPOT_OUTNAME"
  --outdir="$CODPOT_DIR"
  --seed="$CODPOT_SEED"
  --processor="$THREADS"
)
if [[ -n "$CODPOT_SPETHRES" ]]; then
  CODPOT_CMD+=(--spethres="$CODPOT_SPETHRES")
fi

set +e
"${CODPOT_CMD[@]}" > "${LOG_DIR}/1.6.5c_feelnc_codpot.stdout.log" 2> "$CODPOT_STDERR"
RC_CODPOT=$?
set -e
if [[ "$RC_CODPOT" -ne 0 ]]; then
  printf "FEELnc_codpot\t%s\t%s\t%s\n" "$CODPOT_DIR" "$RC_CODPOT" "$CODPOT_STDERR" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FEELnc_codpot failed rc=$RC_CODPOT" | tee -a "$MASTER_LOG"
  exit "$RC_CODPOT"
fi
printf "FEELnc_codpot\t%s\t%s\t%s\n" "$CODPOT_DIR" "$RC_CODPOT" "$CODPOT_STDERR" >> "$SUCCESS_LIST"
echo "[batch] $(date '+%F %T') FEELnc_codpot success" | tee -a "$MASTER_LOG"

: > "$CLASSIFIER_INPUT_GTF"
if [[ -s "$CODPOT_LNCRNA_GTF" ]]; then
  cat "$CODPOT_LNCRNA_GTF" >> "$CLASSIFIER_INPUT_GTF"
fi
if [[ -s "$CODPOT_NOORF_GTF" ]]; then
  cat "$CODPOT_NOORF_GTF" >> "$CLASSIFIER_INPUT_GTF"
fi

if [[ ! -s "$CLASSIFIER_INPUT_GTF" ]]; then
  echo "[ERROR] no lncRNA/noORF GTF available for FEELnc_classifier: $CLASSIFIER_INPUT_GTF" >&2
  exit 1
fi

set +e
FEELnc_classifier.pl \
  --lncrna="$CLASSIFIER_INPUT_GTF" \
  --mrna="$REFERENCE_PROTEIN_CODING_GTF" \
  --window="$CLASSIFIER_WINDOW" \
  --maxwindow="$CLASSIFIER_MAXWINDOW" \
  --log="$CLASSIFIER_LOG" \
  > "$CLASSIFIER_CLASSES" 2> "$CLASSIFIER_STDERR"
RC_CLASSIFIER=$?
set -e
if [[ "$RC_CLASSIFIER" -ne 0 ]]; then
  printf "FEELnc_classifier\t%s\t%s\t%s\n" "$CLASSIFIER_DIR" "$RC_CLASSIFIER" "$CLASSIFIER_STDERR" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FEELnc_classifier failed rc=$RC_CLASSIFIER" | tee -a "$MASTER_LOG"
  exit "$RC_CLASSIFIER"
fi
printf "FEELnc_classifier\t%s\t%s\t%s\n" "$CLASSIFIER_DIR" "$RC_CLASSIFIER" "$CLASSIFIER_STDERR" >> "$SUCCESS_LIST"
echo "[batch] $(date '+%F %T') FEELnc_classifier success" | tee -a "$MASTER_LOG"

COLLECT_CMD=(
  "$PYTHON_BIN"
  "$SCRIPT"
  --query-basic-info "$QUERY_BASIC_INFO"
  --query-gtf "$QUERY_GTF"
  --filter-gtf "$FILTER_GTF"
  --filter-log "$FILTER_LOG"
  --codpot-rf "$CODPOT_RF"
  --codpot-rf-summary "$CODPOT_RF_SUMMARY"
  --codpot-lncrna-gtf "$CODPOT_LNCRNA_GTF"
  --codpot-mrna-gtf "$CODPOT_MRNA_GTF"
  --codpot-tucp-gtf "$CODPOT_TUCP_GTF"
  --codpot-noorf-gtf "$CODPOT_NOORF_GTF"
  --classifier-input-gtf "$CLASSIFIER_INPUT_GTF"
  --classifier-classes "$CLASSIFIER_CLASSES"
  --classifier-log "$CLASSIFIER_LOG"
  --status-tsv "$STATUS_TSV"
  --summary-tsv "$SUMMARY_TSV"
  --manifest "$MANIFEST"
)

if [[ "$OVERWRITE" == "1" ]]; then
  COLLECT_CMD+=(--overwrite)
fi

set +e
"${COLLECT_CMD[@]}" > "${LOG_DIR}/1.6.5c_collect_feelnc_module_outputs.stdout.log" 2> "$COLLECT_STDERR"
RC_COLLECT=$?
set -e
if [[ "$RC_COLLECT" -ne 0 ]]; then
  printf "collect_feelnc_outputs\t%s\t%s\t%s\n" "$QC_DIR" "$RC_COLLECT" "$COLLECT_STDERR" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') collect failed rc=$RC_COLLECT" | tee -a "$MASTER_LOG"
  exit "$RC_COLLECT"
fi
printf "collect_feelnc_outputs\t%s\t%s\t%s\n" "$QC_DIR" "$RC_COLLECT" "$COLLECT_STDERR" >> "$SUCCESS_LIST"

{
  echo "[batch] $(date '+%F %T') FINISH success"
  echo "[batch] FILTER_GTF=$FILTER_GTF"
  echo "[batch] CODPOT_RF=$CODPOT_RF"
  echo "[batch] CLASSIFIER_CLASSES=$CLASSIFIER_CLASSES"
  echo "[batch] STATUS_TSV=$STATUS_TSV"
  echo "[batch] SUMMARY_TSV=$SUMMARY_TSV"
} | tee -a "$MASTER_LOG"
