#!/usr/bin/env bash
# run.sh - GSE66494 pipeline: CKD vs Control gene signature.
# Run from anywhere: the script moves to its own folder (the repo root).
set -Eeuo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

ACCESSION="GSE66494"
BASE_FTP="https://ftp.ncbi.nlm.nih.gov/geo"
URL_RAW="$BASE_FTP/series/GSE66nnn/$ACCESSION/suppl/${ACCESSION}_RAW.tar"
URL_MATRIX="$BASE_FTP/series/GSE66nnn/$ACCESSION/matrix/${ACCESSION}_series_matrix.txt.gz"
URL_ANNOT="$BASE_FTP/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz"

RAW_DIR="data/raw"
SAMPLE_DIR="$RAW_DIR/$ACCESSION"
PROC_DIR="data/processed"
RAW_TAR="$RAW_DIR/${ACCESSION}_RAW.tar"
MATRIX_GZ="$RAW_DIR/${ACCESSION}_series_matrix.txt.gz"
ANNOT_GZ="$RAW_DIR/GPL6480.annot.gz"
CHECKSUMS="checksums.sha256"
EXPECTED_SAMPLES=61

OUT_DIR="results"
SEED=42
FORCE=0
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="logs/run_${RUN_ID}.log"
TMP_DIR="tmp/run_$$"

usage() {
  cat << USAGE
Usage: ./run.sh [-o OUTPUT_DIR] [-s SEED] [-f] [-h]

Runs the GSE66494 pipeline (CKD vs Control gene signature):
download + checksums -> labels -> gene matrix -> signature analysis -> archive.

Options:
  -o DIR   output directory for the results (default: results)
  -s SEED  random seed of the train/test split (default: 42)
  -f       force: rebuild the processed matrix even if it already exists
  -h       show this help and exit

Outputs (in the output directory):
  selected_genes.csv, sample_scores.csv, metrics.json,
  figures/score_boxplot.png, figures/roc_curve.png,
  archives/results_<date>.tar.gz
A timestamped log is written to logs/.
USAGE
}

# ---------- 1. options (before anything is created) ----------
while getopts ":ho:s:f" opt; do
  case "$opt" in
    h) usage; exit 0 ;;
    o) OUT_DIR="$OPTARG" ;;
    s) SEED="$OPTARG" ;;
    f) FORCE=1 ;;
    :)  echo "ERREUR: option -$OPTARG needs a value (see ./run.sh -h)" >&2; exit 2 ;;
    \?) echo "ERREUR: unknown option -$OPTARG (see ./run.sh -h)" >&2; exit 2 ;;
  esac
done
if ! [[ "$SEED" =~ ^[0-9]+$ ]]; then
  echo "ERREUR: the seed (-s) must be a non-negative integer, got '$SEED'" >&2
  exit 2
fi

# ---------- 2. log, temp folder, trap ----------
mkdir -p logs tmp "$TMP_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { echo "ERREUR: $*" >&2; exit 1; }

cleanup() { rm -rf "$TMP_DIR"; }
trap cleanup EXIT
trap 'log "Interrupted, cleaning temporary files"; exit 130' INT TERM
trap 'echo "ERREUR: a command failed at line $LINENO of run.sh" >&2' ERR

# ---------- 3. functions ----------
check_file() {
  [[ -f "$1" ]] || die "required file missing: $1"
}

check_tools() {
  local tool
  for tool in curl tar sha256sum python3 zcat; do
    command -v "$tool" > /dev/null 2>&1 || die "tool not installed: $tool"
  done
  python3 -c "import numpy, pandas, scipy, sklearn, matplotlib" 2> /dev/null \
    || die "Python packages missing: activate the venv (source .venv/bin/activate) and run: pip install -r requirements.txt"
}

download_file() {  # usage: download_file URL DEST
  local url="$1" dest="$2" part
  part="$TMP_DIR/$(basename "$dest")"
  if [[ -f "$dest" ]]; then
    log "Already present, skipping download: $dest"
    return 0
  fi
  log "Downloading $url"
  curl -fsSL --retry 3 -o "$part" "$url" || die "download failed: $url"
  mkdir -p "$(dirname "$dest")"
  mv -f "$part" "$dest"
  chmod 444 "$dest"
}

verify_checksums() {
  check_file "$CHECKSUMS"
  log "Verifying SHA256 checksums"
  sha256sum -c "$CHECKSUMS" || die "checksum verification failed: a raw file is corrupted or changed"
}

count_samples() {
  if [[ -d "$SAMPLE_DIR" ]]; then
    find "$SAMPLE_DIR" -maxdepth 1 -name 'GSM*.txt.gz' | wc -l
  else
    echo 0
  fi
}

extract_raw() {
  local n
  n="$(count_samples)"
  if [[ "$n" -eq "$EXPECTED_SAMPLES" ]]; then
    log "Samples already extracted ($n files)"
    return 0
  fi
  log "Extracting $RAW_TAR"
  mkdir -p "$SAMPLE_DIR"
  tar -xf "$RAW_TAR" -C "$SAMPLE_DIR"
  n="$(count_samples)"
  [[ "$n" -eq "$EXPECTED_SAMPLES" ]] || die "expected $EXPECTED_SAMPLES sample files after extraction, found $n"
  chmod 444 "$SAMPLE_DIR"/*.txt.gz
}

series_field() {  # pipe chain: one value per line from a series-matrix header line
  zcat "$MATRIX_GZ" | grep "^!Sample_$1" | tr '\t' '\n' | sed 's/"//g' | tail -n +2
}

build_labels() {
  log "Building labels from the series matrix"
  check_file "$MATRIX_GZ"
  mkdir -p "$PROC_DIR"
  series_field geo_accession > "$TMP_DIR/gsms.txt"
  series_field title > "$TMP_DIR/titles.txt"
  paste "$TMP_DIR/gsms.txt" "$TMP_DIR/titles.txt" \
    | awk -F'\t' 'BEGIN{OFS="\t"; print "gsm","title","study_set","condition"}
                  {split($2,a,"_"); print $1,$2,tolower(a[1]),($2 ~ /CKD/ ? "CKD" : "Control")}' \
    > "$TMP_DIR/labels.tsv"
  local n
  n="$(( $(wc -l < "$TMP_DIR/labels.tsv") - 1 ))"
  [[ "$n" -eq "$EXPECTED_SAMPLES" ]] || die "expected $EXPECTED_SAMPLES labelled samples, found $n"
  mv -f "$TMP_DIR/labels.tsv" "$PROC_DIR/labels.tsv"
  log "Label counts (study set, condition):"
  tail -n +2 "$PROC_DIR/labels.tsv" | cut -f3,4 | sort | uniq -c
}

build_matrix() {
  local matrix="$PROC_DIR/expression_matrix.tsv"
  if [[ "$FORCE" -eq 0 && -s "$matrix" ]]; then
    log "Gene matrix already exists, skipping (use -f to rebuild)"
    return 0
  fi
  log "Building the gene expression matrix (takes about a minute)"
  if ! python3 scripts/prepare_matrix.py 2> "$TMP_DIR/prepare.err"; then
    cat "$TMP_DIR/prepare.err" >&2
    die "matrix preparation failed"
  fi
}

run_analysis() {
  log "Running the signature analysis (seed $SEED)"
  python3 scripts/analyze_signature.py --out-dir "$OUT_DIR" --seed "$SEED" \
    || die "signature analysis failed"
}

archive_results() {
  local f
  for f in selected_genes.csv sample_scores.csv metrics.json \
           figures/score_boxplot.png figures/roc_curve.png; do
    check_file "$OUT_DIR/$f"
  done
  mkdir -p "$OUT_DIR/archives"
  local archive="$OUT_DIR/archives/results_${RUN_ID}.tar.gz"
  tar -czf "$archive" -C "$OUT_DIR" selected_genes.csv sample_scores.csv metrics.json \
    figures/score_boxplot.png figures/roc_curve.png
  log "Results archived: $archive"
}

# ---------- 4. main ----------
main() {
  log "Pipeline start ($ACCESSION), log: $LOG_FILE"
  check_tools
  check_file scripts/prepare_matrix.py
  check_file scripts/analyze_signature.py
  check_file "$CHECKSUMS"
  download_file "$URL_RAW" "$RAW_TAR"
  download_file "$URL_MATRIX" "$MATRIX_GZ"
  download_file "$URL_ANNOT" "$ANNOT_GZ"
  verify_checksums
  extract_raw
  build_labels
  build_matrix
  run_analysis
  archive_results
  log "Pipeline finished successfully"
}

main
