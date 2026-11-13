Projet Maladie – La maladie rénale chronique

The disease survey: can a small gene signature tell chronic kidney disease (CKD) apart from healthy kidney tissue?

Members: [CHOULI Amel El Batoul 232331673901,ARABI Ghalia Bahria 212132027986, MAOUDJ Fella 23233179110, BENHADID Nafissa Racha 222231649203, BENCHOUCH Ines 222231440315]

We work with the public dataset GSE66494 from NCBI GEO: kidney biopsies measured on Agilent GPL6480 arrays, 61 samples in total (CKD patients and controls).

The whole project runs with one command, ./run.sh. It downloads the data, checks that the files are intact, builds a gene expression matrix, picks a 40-gene signature on a training set, tests it on samples it has never seen, and archives the results.
What you need

We developed and tested everything on Ubuntu 26.04 (WSL2) with Python 3.14. You need bash, curl, tar, gzip (for zcat), sha256sum, git and python3 with venv. Plan for about 1 GB of free disk space and an internet connection for the first run.
Installation (clean Ubuntu)

sudo apt update && sudo apt install -y git curl python3 python3-venv python3-pip
git clone https://github.com/Batoullovesthemoon/projet_maladie.git projet_maladie
cd projet_maladie
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

How to run it

Always run from the repository root, with the virtual environment active:

./run.sh -h              # show the help
./run.sh                 # full pipeline, default settings
./run.sh -o my_results   # write the results in another folder
./run.sh -s 7            # use another random seed for the train/test split
./run.sh -f              # force a rebuild of the processed gene matrix

Option	What it does	Default
-h	shows the help and exits	
-o DIR	output folder for the results	results
-s SEED	random seed for the train/test split	42
-f	rebuilds the processed matrix even if it already exists	off

Exit codes: 0 success, 1 error (missing file, failed download, failed checksum...), 2 wrong option, 130 stopped with Ctrl+C.

How long does it take? The first run downloads about 160 MB, so count several minutes (around 10 on a slow connection). After that the data is reused and a run takes less than a minute.
What the pipeline does, step by step

    Downloads GSE66494_RAW.tar, the series matrix and the GPL6480 annotation from the NCBI FTP server into data/raw/. Files that are already there are not downloaded again.
    Checks the three files with sha256sum -c checksums.sha256, then makes them read-only.
    Extracts the 61 sample files.
    Builds data/processed/labels.tsv (CKD or Control for each sample) from the series matrix, using a chain of piped shell commands.
    scripts/prepare_matrix.py builds the genes x samples matrix.
    scripts/analyze_signature.py selects the signature, scores the samples, evaluates it and writes the results and figures.
    Archives the results in results/archives/ and writes a timestamped log in logs/.

The method in a nutshell

Cleaning the data

    Control spots are removed. Arrays where fewer than 40% of the probes are detected above background are excluded. This removes 3 of the 61 samples (GSM1623315, GSM1623352, GSM1623356); they are listed in data/processed/excluded_samples.tsv.
    Probes detected in fewer than 20% of the remaining samples are dropped. The values are then log2-transformed and quantile-normalized.
    Probes with no gene symbol, or with several genes, are dropped. When several probes measure the same gene, we keep the most variable one (data/processed/gene_probe_map.tsv).

Splitting the samples

    The 58 remaining samples are split 70/30 into train and test, keeping the same CKD/Control proportions in both.

Building the signature (training samples only)

    For each gene: Welch t-test, Benjamini-Hochberg adjusted p-value and log2 fold change.
    We keep the 20 most significant up-regulated and the 20 most significant down-regulated genes, with |log2 fold change| >= 0.5.
    The mean and standard deviation of each gene (used for z-scores) and the decision threshold also come from the training samples only.

Scoring and evaluating

    score = mean(z of up genes) - mean(z of down genes). A sample is predicted CKD if its score is at or above the threshold.
    The threshold is chosen on the training set with Youden's index.
    On the test set we report AUC, balanced accuracy, sensitivity and specificity, with CKD as the positive class.

What you should get (seed 42)

Train: 34 CKD + 6 Control. Test: 16 CKD + 2 Control. The signature has 20 up and 20 down genes and the threshold is about -0.213. On the test set we get AUC 1.00, balanced accuracy 0.97, sensitivity 0.94 and specificity 1.00.

Please read these numbers with caution: the test set has only 2 controls, so they are very uncertain (see Limitations).
Repository tree

projet_maladie/
|-- run.sh                  entry point of the pipeline
|-- requirements.txt        Python dependencies
|-- checksums.sha256        SHA256 fingerprints of the raw files
|-- SOURCES.md              accession, URLs, download dates, checksums
|-- crontab.example         example of an automatic weekly run
|-- README.md               this file
|-- scripts/                prepare_matrix.py, analyze_signature.py
|-- tests/                  test_analyze_signature.py
|-- data/
|   |-- raw/                downloaded data (read-only, not in Git)
|   `-- processed/          labels and gene matrix (not in Git)
|-- results/                results of a run (not in Git)
|   |-- figures/
|   `-- archives/
|-- logs/                   timestamped logs (not in Git)
`-- tmp/                    temporary files, cleaned up by the trap in run.sh

Files produced
File	Content
data/processed/labels.tsv	sample, title, study set, condition (CKD / Control)
data/processed/excluded_samples.tsv	arrays removed by quality control
data/processed/expression_matrix.tsv	genes x samples matrix (log2, quantile-normalized)
data/processed/gene_probe_map.tsv	probe kept for each gene
results/selected_genes.csv	the 20 up and 20 down genes, with log2 fold change, p-values, training mean and std
results/sample_scores.csv	score, train/test split and predicted label of every sample
results/metrics.json	parameters, sample counts, threshold, train and test metrics
results/figures/score_boxplot.png	signature score by group, with the threshold
results/figures/roc_curve.png	ROC curves (train and test)
results/archives/results_<date>_<time>.tar.gz	archive of the five result files above
logs/run_<date>_<time>.log	timestamped log of the run

With -o DIR, the files normally in results/ are written to DIR instead.
Tests

python -m pytest tests -v

There are 9 tests, run on a small synthetic dataset. They cover the statistics, the gene selection, the score formula, the threshold, the metrics, the whole script end to end, the error on a missing file, and a data-leakage check: scrambling the test samples must not change the selected genes or the threshold.
Automatic run (cron)

crontab.example runs the pipeline every Monday at 03:00 and appends the output to logs/cron.log. To use it, run crontab -e, paste the last line of the file, save, and check with crontab -l. If your copy of the repository is not in ~/projet_maladie, change the path in that line.
Reproducibility

    The split uses a fixed seed (-s), so the same seed always gives the same results.
    The downloaded files are checked against checksums.sha256 on every run.
    Raw data is read-only, and run.sh can be rerun at any time without cleaning anything by hand.
    If a run is interrupted, a trap removes the temporary files.
If something goes wrong

    Permission denied when running ./run.sh: run chmod +x run.sh scripts/*.py.
    ERREUR: Python packages missing: activate the environment (source .venv/bin/activate) and run pip install -r requirements.txt.
    ERREUR: tool not installed: ...: install the tool, for example sudo apt install curl.
    ERREUR: download failed: check your internet connection and run again. Files that were already downloaded are kept.
    checksum verification failed: a raw file is corrupted. Delete it in data/raw/ (it is read-only, so use rm -f) and run again; it will be downloaded again.

Limitations

    Only 8 of the 58 samples are controls, so the test set contains just 2 of them. The test metrics are therefore very uncertain.
    CKD arrays have a lower overall signal than control arrays. Normalization removes the global difference, but with this dataset we cannot fully separate a technical effect from the disease itself.
    Quantile normalization and the choice of one probe per gene use all the samples (without their labels), including the test samples.
    We drop probes without a gene symbol (6056) and probes matching several genes (29).

Sources

See SOURCES.md for the accession (GSE66494), the platform (GPL6480), the URLs, the download dates and the checksums.