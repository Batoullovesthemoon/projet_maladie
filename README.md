# Projet Maladie - Group [GROUP NUMBER]

**Enquête sur une maladie: a gene-expression signature of chronic kidney disease (CKD).**

Members: [MEMBERS]

Dataset: NCBI GEO **GSE66494** (kidney biopsies, Agilent GPL6480 arrays, 61 samples, CKD vs control).
One command, `./run.sh`, downloads the data, checks its integrity, builds a gene
expression matrix, selects a 40-gene signature on a training set, evaluates it on a
test set, and archives the results.

## Requirements

Tested on Ubuntu 26.04 (WSL2) with Python 3.14. Needs `bash`, `curl`, `tar`, `gzip`
(`zcat`), `sha256sum`, `git`, `python3` with `venv`, about 1 GB of free disk space,
and an internet connection for the first run.

## Installation (clean Ubuntu)

    sudo apt update && sudo apt install -y git curl python3 python3-venv python3-pip
    git clone <REPOSITORY_URL> projet_maladie
    cd projet_maladie
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt

## Usage

    ./run.sh -h              # help
    ./run.sh                 # full pipeline with default settings
    ./run.sh -o my_results   # write results to another folder
    ./run.sh -s 7            # use another random seed for the train/test split
    ./run.sh -f              # force: rebuild the processed gene matrix

| Option | Meaning | Default |
|---|---|---|
| `-h` | show help and exit | |
| `-o DIR` | output directory for the results | `results` |
| `-s SEED` | random seed of the train/test split | `42` |
| `-f` | rebuild the processed matrix even if it exists | off |

Run it from the repository root, with the virtual environment active.
Exit codes: `0` success, `1` error (missing file, failed download, failed checksum, ...),
`2` wrong option, `130` interrupted with Ctrl+C.

**Duration:** the first run downloads about 160 MB and takes several minutes
(about 10 minutes on a slow connection). Later runs reuse the downloaded data
and take under a minute.

## What the pipeline does

1. Downloads `GSE66494_RAW.tar`, the series matrix and the GPL6480 annotation from the
   NCBI FTP server into `data/raw/` (skipped if the files are already there).
2. Verifies the three files with `sha256sum -c checksums.sha256` and makes them read-only.
3. Extracts the 61 sample files.
4. Builds `data/processed/labels.tsv` (CKD / Control for each sample) from the series
   matrix with a chain of piped shell commands.
5. `scripts/prepare_matrix.py`: builds the genes x samples matrix.
6. `scripts/analyze_signature.py`: selects the signature, scores samples, evaluates, writes
   results and figures.
7. Archives the results in `results/archives/` and writes a timestamped log in `logs/`.

## Method in short

- **Quality control:** control spots are removed. Arrays with fewer than 40% of probes
  detected above background are excluded (3 of 61 samples: GSM1623315, GSM1623352,
  GSM1623356, listed in `data/processed/excluded_samples.tsv`). Probes detected in
  fewer than 20% of the remaining samples are removed. Values are log2-transformed and
  quantile-normalized.
- **Probe to gene:** probes without a gene symbol, or with several genes, are dropped.
  When several probes measure the same gene, the most variable one is kept
  (`data/processed/gene_probe_map.tsv`).
- **Split:** stratified 70/30 train/test split of the 58 samples (same CKD/Control
  proportions in both).
- **Everything learned from data uses the training samples only:** Welch t-test per gene,
  Benjamini-Hochberg adjusted p-values, log2 fold change, selection of the top 20 up and
  top 20 down genes (ranked by adjusted p-value, with |log2 fold change| >= 0.5),
  the mean and standard deviation of each gene (for z-scores), and the threshold.
- **Score:** `score = mean(z of up genes) - mean(z of down genes)`.
  A sample is predicted CKD if its score is at or above the threshold.
- **Threshold:** chosen on the training set with Youden's index.
- **Evaluation on the test set:** AUC, balanced accuracy, sensitivity and specificity
  (CKD is the positive class).

## What you should get (default seed 42)

Split: train 34 CKD + 6 Control, test 16 CKD + 2 Control. 20 up and 20 down genes,
threshold about -0.213. On the test set: AUC 1.00, balanced accuracy 0.97,
sensitivity 0.94, specificity 1.00. The test set contains only 2 controls, so these
numbers are very uncertain (see Limitations).

## Repository tree

    projet_maladie/
    |-- run.sh                  entry point of the pipeline
    |-- requirements.txt        Python dependencies
    |-- checksums.sha256        SHA256 fingerprints of the raw files
    |-- SOURCES.md              accession, URLs, download date, checksums
    |-- crontab.example         example of an automatic weekly run
    |-- README.md               this file
    |-- scripts/                prepare_matrix.py, analyze_signature.py
    |-- tests/                  pytest tests
    |-- data/
    |   |-- raw/                downloaded data (read-only, not in Git)
    |   `-- processed/          labels and gene matrix (not in Git)
    |-- results/                results of a run (not in Git)
    |   |-- figures/
    |   `-- archives/
    |-- logs/                   timestamped logs (not in Git)
    `-- tmp/                    temporary files, cleaned by the trap in run.sh

## Files produced

| File | Description |
|---|---|
| `data/processed/labels.tsv` | sample, title, study set, condition (CKD / Control) |
| `data/processed/excluded_samples.tsv` | arrays excluded by quality control |
| `data/processed/expression_matrix.tsv` | genes x samples matrix (log2, quantile-normalized) |
| `data/processed/gene_probe_map.tsv` | probe kept for each gene |
| `results/selected_genes.csv` | the 20 up and 20 down genes with log2 fold change, p-values and training mean/std |
| `results/sample_scores.csv` | score, train/test split and predicted label of each sample |
| `results/metrics.json` | parameters, sample counts, threshold, train and test metrics |
| `results/figures/score_boxplot.png` | signature score by group, with the threshold |
| `results/figures/roc_curve.png` | ROC curves (train and test) |
| `results/archives/results_<date>_<time>.tar.gz` | archive of the five result files above |
| `logs/run_<date>_<time>.log` | timestamped log of the run |

With `-o DIR`, the `results/` files go to `DIR` instead.

## Tests

    python -m pytest tests -v

9 tests on a small synthetic dataset: statistics, gene selection, score formula,
threshold, metrics, the whole script, the error on a missing file, and a check that
the test samples do not influence the selected genes or the threshold (no data leakage).

## Automatic run (cron)

`crontab.example` runs the pipeline every Monday at 03:00 and appends the output to
`logs/cron.log`. To install it: `crontab -e`, paste the last line of the file, save,
and check with `crontab -l`. Change the project path in that line if your
repository is not in `~/projet_maladie`.

## Reproducibility

- The split uses a fixed seed (`-s`), so the same seed gives the same results.
- The downloaded files are checked against `checksums.sha256` at every run.
- Raw data is read-only, and `run.sh` can be rerun without any manual cleanup.
- A `trap` removes the temporary files if the run is interrupted.

## Work split

| Member | Work |
|---|---|
| [MEMBER 1] | [e.g. data download, checksums, run.sh] |
| [MEMBER 2] | [e.g. matrix preparation, quality control] |
| [MEMBER 3] | [e.g. signature analysis, tests, report] |

## Troubleshooting

- `Permission denied` when running `./run.sh`: run `chmod +x run.sh scripts/*.py`.
- `ERREUR: Python packages missing`: activate the environment
  (`source .venv/bin/activate`) and run `pip install -r requirements.txt`.
- `ERREUR: tool not installed: ...`: install it, for example `sudo apt install curl`.
- `ERREUR: download failed`: check your internet connection and rerun. Files already
  downloaded are not downloaded again.
- `checksum verification failed`: a raw file is corrupted. Delete that file in
  `data/raw/` (it is read-only, so use `rm -f`) and rerun, so it is downloaded again.

## Limitations

- Only 8 of the 58 samples are controls, so the test set has 2 controls: the test
  metrics have a very wide uncertainty.
- CKD arrays have a lower overall signal than control arrays. Normalization removes the
  global difference, but a technical effect cannot be completely separated from the
  disease with this dataset.
- Quantile normalization and the choice of one probe per gene use all samples (without
  labels), including the test samples.
- Probes without a gene symbol (6056) or with several genes (29) are dropped.

## Sources

See `SOURCES.md` (accession GSE66494, platform GPL6480, URLs, dates, checksums).
