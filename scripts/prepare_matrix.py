#!/usr/bin/env python3
"""Build a probes x samples log2 expression matrix from Agilent FE files.

Steps: read files -> drop failed arrays -> filter probes -> quantile normalize.
"""
import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

COLS = ["ProbeName", "ControlType", "gProcessedSignal", "gIsWellAboveBG"]


def die(msg):
    print(f"ERREUR: {msg}", file=sys.stderr)
    sys.exit(1)


def find_header_line(path):
    """Return the 0-based line number of the FEATURES header line."""
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith("FEATURES"):
                return i
    die(f"no FEATURES header line in {path}")


def read_sample(path):
    """Read one sample: log2 signal, 'detected' flag per probe, array quality."""
    skip = find_header_line(path)
    df = pd.read_csv(path, sep="\t", skiprows=skip, usecols=COLS)
    df = df[df["ControlType"] == 0]  # keep real probes, drop control spots
    frac_detected = df["gIsWellAboveBG"].mean()  # quality of the whole array
    grp = df.groupby("ProbeName")  # some probes are printed several times
    sig = grp["gProcessedSignal"].mean()
    sig = sig.where(sig > 0)  # log2 needs positive values
    detected = grp["gIsWellAboveBG"].mean() >= 0.5
    return np.log2(sig), detected, frac_detected


def quantile_normalize(df):
    """Give every column (sample) the same distribution of values."""
    ranks = df.rank(method="first").astype(int) - 1  # rank inside each column
    mean_sorted = np.sort(df.values, axis=0).mean(axis=1)  # average distribution
    return pd.DataFrame(mean_sorted[ranks.values], index=df.index, columns=df.columns)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", default="data/raw/GSE66494")
    p.add_argument("--labels", default="data/processed/labels.tsv")
    p.add_argument("--out", default="data/processed/expression_matrix.tsv")
    p.add_argument("--min-detected", type=float, default=0.2,
                   help="keep probes detected above background in at least "
                        "this fraction of samples (default 0.2)")
    p.add_argument("--min-sample-detected", type=float, default=0.4,
                   help="exclude arrays where fewer than this fraction of "
                        "probes are above background (default 0.4)")
    args = p.parse_args()

    raw_dir, labels_path, out = Path(args.raw_dir), Path(args.labels), Path(args.out)
    if not raw_dir.is_dir():
        die(f"raw directory not found: {raw_dir}")
    if not labels_path.is_file():
        die(f"labels file not found: {labels_path}")

    labels = pd.read_csv(labels_path, sep="\t")
    expr, det, fracs = {}, {}, {}
    for gsm in labels["gsm"]:
        files = sorted(raw_dir.glob(f"{gsm}_*.txt.gz"))
        if len(files) != 1:
            die(f"expected 1 file for {gsm} in {raw_dir}, found {len(files)}")
        expr[gsm], det[gsm], fracs[gsm] = read_sample(files[0])
        print(f"read {gsm}", file=sys.stderr)

    expr = pd.DataFrame(expr)  # rows = probes, columns = samples
    det = pd.DataFrame(det)
    n_start = expr.shape[0]

    # 1. drop failed arrays and record why
    frac = pd.Series(fracs, name="frac_detected")
    bad = frac[frac < args.min_sample_detected].sort_index()
    out.parent.mkdir(parents=True, exist_ok=True)
    bad.to_csv(out.parent / "excluded_samples.tsv", sep="\t", index_label="gsm")
    expr = expr.drop(columns=bad.index)
    det = det.drop(columns=bad.index)

    # 2. filter probes (on the samples that remain)
    keep = (det.mean(axis=1) >= args.min_detected) & (expr.isna().mean(axis=1) <= 0.1)
    expr = expr[keep]
    n_filled = int(expr.isna().sum().sum())
    expr = expr.apply(lambda row: row.fillna(row.median()), axis=1)

    # 3. quantile normalization
    expr = quantile_normalize(expr)
    expr.to_csv(out, sep="\t", index_label="ProbeName", float_format="%.4f")

    meds = expr.median()
    print(f"samples excluded        : {len(bad)} {list(bad.index)}")
    print(f"probes before filtering : {n_start}")
    print(f"probes kept             : {expr.shape[0]}")
    print(f"missing values filled   : {n_filled}")
    print(f"samples kept            : {expr.shape[1]}")
    print(f"sample medians (log2)   : min {meds.min():.2f}, max {meds.max():.2f}")
    print(f"written: {out}")


if __name__ == "__main__":
    main()
