#!/usr/bin/env python3
"""Build a genes x samples log2 expression matrix from Agilent FE files.

Steps: read files -> drop failed arrays -> filter probes -> quantile normalize
-> map probes to gene symbols (GPL6480) -> keep the most variable probe per gene.
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


def read_annotation(path):
    """Read the GPL6480 table: probe ID -> gene symbol."""
    skip = None
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith("!platform_table_begin"):
                skip = i + 1  # the column names are on the next line
                break
    if skip is None:
        die(f"no !platform_table_begin line in {path}")
    ann = pd.read_csv(path, sep="\t", skiprows=skip, usecols=["ID", "Gene symbol"],
                      dtype=str, quoting=3)
    ann = ann[~ann["ID"].str.startswith("!", na=False)]  # drop '!platform_table_end'
    ann = ann.rename(columns={"Gene symbol": "gene"})
    ann["gene"] = ann["gene"].str.strip()
    return ann


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", default="data/raw/GSE66494")
    p.add_argument("--labels", default="data/processed/labels.tsv")
    p.add_argument("--annot", default="data/raw/GPL6480.annot.gz")
    p.add_argument("--out", default="data/processed/expression_matrix.tsv")
    p.add_argument("--min-detected", type=float, default=0.2,
                   help="keep probes detected above background in at least "
                        "this fraction of samples (default 0.2)")
    p.add_argument("--min-sample-detected", type=float, default=0.4,
                   help="exclude arrays where fewer than this fraction of "
                        "probes are above background (default 0.4)")
    args = p.parse_args()

    raw_dir, labels_path = Path(args.raw_dir), Path(args.labels)
    annot_path, out = Path(args.annot), Path(args.out)
    if not raw_dir.is_dir():
        die(f"raw directory not found: {raw_dir}")
    if not labels_path.is_file():
        die(f"labels file not found: {labels_path}")
    if not annot_path.is_file():
        die(f"annotation file not found: {annot_path}")

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
    n_probes = expr.shape[0]

    # 3. quantile normalization
    expr = quantile_normalize(expr)

    # 4. probes -> genes: drop probes without a clear gene symbol, then keep
    #    the most variable probe of each gene (uses no labels)
    ann = read_annotation(annot_path)
    gene_of = ann.drop_duplicates("ID").set_index("ID")["gene"]
    gene = expr.index.to_series().map(gene_of)
    ambiguous = gene.str.contains("///", na=False)  # several genes in one probe
    ok = gene.notna() & (gene != "") & ~ambiguous
    n_no_gene = int((~(gene.notna() & (gene != ""))).sum())
    expr, gene = expr[ok], gene[ok]
    order = expr.var(axis=1).sort_values(ascending=False).index
    best = gene.loc[order].drop_duplicates()  # first = most variable probe per gene
    expr_gene = expr.loc[best.index]
    expr_gene.index = best.values
    expr_gene = expr_gene.sort_index()

    expr_gene.to_csv(out, sep="\t", index_label="Gene", float_format="%.4f")
    map_df = best.rename("Gene").rename_axis("ProbeName").reset_index()
    map_df.sort_values("Gene").to_csv(out.parent / "gene_probe_map.tsv",
                                      sep="\t", index=False)

    meds = expr.median()
    print(f"samples excluded        : {len(bad)} {list(bad.index)}")
    print(f"probes before filtering : {n_start}")
    print(f"probes kept             : {n_probes}")
    print(f"missing values filled   : {n_filled}")
    print(f"probes without gene     : {n_no_gene}")
    print(f"probes with 2+ genes    : {int(ambiguous.sum())}")
    print(f"genes kept (1 probe each): {expr_gene.shape[0]}")
    print(f"samples kept            : {expr_gene.shape[1]}")
    print(f"sample medians (log2)   : min {meds.min():.2f}, max {meds.max():.2f}")
    print(f"written: {out}")


if __name__ == "__main__":
    main()
