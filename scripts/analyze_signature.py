#!/usr/bin/env python3
"""Find a CKD vs Control gene signature on a train set, evaluate it on a test set.

Everything learned from the data (selected genes, z-score means/stds, threshold)
is computed on the TRAIN samples only. The TEST samples are only scored.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import (balanced_accuracy_score, confusion_matrix,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import train_test_split

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def die(msg):
    print(f"ERREUR: {msg}", file=sys.stderr)
    sys.exit(1)


def welch_bh(expr_train, y_train):
    """Welch t-test per gene (CKD vs Control), BH-adjusted p-values, log2 fold change."""
    ckd = expr_train[y_train.index[y_train == 1]]
    ctl = expr_train[y_train.index[y_train == 0]]
    _, p = stats.ttest_ind(ckd.values, ctl.values, axis=1, equal_var=False)
    p = np.nan_to_num(p, nan=1.0)
    padj = stats.false_discovery_control(p, method="bh")
    return pd.DataFrame({"log2fc": ckd.mean(axis=1) - ctl.mean(axis=1),
                         "pvalue": p, "padj": padj}, index=expr_train.index)


def select_genes(res, top_n, min_lfc):
    """Top N up and top N down genes, ranked by adjusted p-value, with |log2FC| >= min_lfc."""
    up = res[res["log2fc"] >= min_lfc].sort_values(["padj", "pvalue"]).head(top_n)
    down = res[res["log2fc"] <= -min_lfc].sort_values(["padj", "pvalue"]).head(top_n)
    return up.index.tolist(), down.index.tolist()


def signature_score(expr, up, down, mu, sd):
    """score = mean(z of up genes) - mean(z of down genes), z uses TRAIN mean/std."""
    genes = up + down
    z = expr.loc[genes].sub(mu[genes], axis=0).div(sd[genes], axis=0)
    return z.loc[up].mean() - z.loc[down].mean()


def youden_threshold(y, score):
    """Threshold maximizing sensitivity + specificity - 1 (Youden's index)."""
    fpr, tpr, thr = roc_curve(y, score)
    return float(thr[np.argmax(tpr - fpr)])


def evaluate(y, score, thr):
    """Metrics with CKD as the positive class."""
    pred = (score >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {"auc": float(roc_auc_score(y, score)),
            "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
            "sensitivity": float(tp / (tp + fn)),
            "specificity": float(tn / (tn + fp)),
            "tp": int(tp), "fn": int(fn), "tn": int(tn), "fp": int(fp)}


def make_figures(table, y_train, s_train, y_test, s_test, thr, fig_dir):
    """Boxplot of scores by group and ROC curves."""
    fig, ax = plt.subplots(figsize=(5, 5))
    rng = np.random.default_rng(0)
    groups = ["Control", "CKD"]
    ax.boxplot([table.loc[table.condition == g, "score"] for g in groups],
               positions=[0, 1], widths=0.5, showfliers=False)
    for i, g in enumerate(groups):
        for split, marker, color in (("train", "o", "tab:gray"), ("test", "^", "tab:red")):
            d = table[(table.condition == g) & (table.split == split)]
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(d)), d.score,
                       marker=marker, color=color, alpha=0.8, label=split if i == 0 else None)
    ax.axhline(thr, ls="--", color="k", lw=1, label="threshold (train, Youden)")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(groups)
    ax.set_ylabel("signature score")
    ax.set_title("Signature score by group")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "score_boxplot.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    for y, s, name, style in ((y_train, s_train, "train", "--"), (y_test, s_test, "test", "-")):
        fpr, tpr, _ = roc_curve(y, s)
        ax.plot(fpr, tpr, style, label=f"{name} (AUC = {roc_auc_score(y, s):.2f})")
    ax.plot([0, 1], [0, 1], ":", color="gray")
    ax.set_xlabel("1 - specificity")
    ax.set_ylabel("sensitivity")
    ax.set_title("ROC curve (positive class = CKD)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "roc_curve.png", dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", default="data/processed/expression_matrix.tsv")
    p.add_argument("--labels", default="data/processed/labels.tsv")
    p.add_argument("--out-dir", default="results")
    p.add_argument("--test-size", type=float, default=0.3,
                   help="fraction of samples kept for the test set (default 0.3)")
    p.add_argument("--top-n", type=int, default=20,
                   help="number of up and of down genes in the signature (default 20)")
    p.add_argument("--min-lfc", type=float, default=0.5,
                   help="minimum |log2 fold change| to select a gene (default 0.5)")
    p.add_argument("--seed", type=int, default=42, help="random seed of the split")
    args = p.parse_args()

    matrix_path, labels_path, out_dir = Path(args.matrix), Path(args.labels), Path(args.out_dir)
    if not matrix_path.is_file():
        die(f"expression matrix not found: {matrix_path}")
    if not labels_path.is_file():
        die(f"labels file not found: {labels_path}")
    if not 0 < args.test_size < 1:
        die("--test-size must be between 0 and 1")

    expr = pd.read_csv(matrix_path, sep="\t", index_col=0)
    labels = pd.read_csv(labels_path, sep="\t", index_col="gsm")
    missing = [c for c in expr.columns if c not in labels.index]
    if missing:
        die(f"samples without a label: {missing[:5]}")
    labels = labels.loc[expr.columns]
    y = (labels["condition"] == "CKD").astype(int)

    # 1. stratified split: same CKD/Control proportions in train and test
    try:
        train_ids, test_ids = train_test_split(expr.columns.to_list(), test_size=args.test_size,
                                               stratify=y.values, random_state=args.seed)
    except ValueError as e:
        die(f"cannot split the samples: {e}")
    y_train, y_test = y[train_ids], y[test_ids]
    if y_test.nunique() < 2 or y_train.nunique() < 2:
        die("train or test set has only one class; change --seed or --test-size")

    # 2. everything below is learned on TRAIN only
    expr_train = expr[train_ids]
    res = welch_bh(expr_train, y_train)
    up, down = select_genes(res, args.top_n, args.min_lfc)
    if not up or not down:
        die("no up or no down gene passes the filter; lower --min-lfc")
    if len(up) < args.top_n or len(down) < args.top_n:
        print(f"ATTENTION: only {len(up)} up and {len(down)} down genes pass the filter",
              file=sys.stderr)
    mu, sd = expr_train.mean(axis=1), expr_train.std(axis=1, ddof=1)

    scores = signature_score(expr, up, down, mu, sd)
    thr = youden_threshold(y_train, scores[train_ids])

    # 3. evaluation
    train_metrics = evaluate(y_train, scores[train_ids], thr)
    test_metrics = evaluate(y_test, scores[test_ids], thr)

    # 4. outputs
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    sel = pd.concat([res.loc[up].assign(direction="up"), res.loc[down].assign(direction="down")])
    sel["train_mean"], sel["train_std"] = mu[sel.index], sd[sel.index]
    sel.index.name = "gene"
    sel[["direction", "log2fc", "pvalue", "padj", "train_mean", "train_std"]].to_csv(
        out_dir / "selected_genes.csv", float_format="%.6g")

    train_set = set(train_ids)
    table = pd.DataFrame({"gsm": expr.columns, "title": labels["title"].values,
                          "condition": labels["condition"].values,
                          "split": ["train" if c in train_set else "test" for c in expr.columns],
                          "score": scores.values})
    table["predicted"] = np.where(table["score"] >= thr, "CKD", "Control")
    table.to_csv(out_dir / "sample_scores.csv", index=False, float_format="%.6g")

    metrics = {"parameters": {"test_size": args.test_size, "top_n": args.top_n,
                              "min_lfc": args.min_lfc, "seed": args.seed},
               "n_train": {"CKD": int((y_train == 1).sum()), "Control": int((y_train == 0).sum())},
               "n_test": {"CKD": int((y_test == 1).sum()), "Control": int((y_test == 0).sum())},
               "n_up_genes": len(up), "n_down_genes": len(down),
               "threshold": thr, "train": train_metrics, "test": test_metrics}
    with open(out_dir / "metrics.json", "w") as fh:
        json.dump(metrics, fh, indent=2)

    make_figures(table, y_train, scores[train_ids], y_test, scores[test_ids], thr, fig_dir)

    print(f"train: {metrics['n_train']} | test: {metrics['n_test']}")
    print(f"genes selected: {len(up)} up, {len(down)} down | threshold: {thr:.3f}")
    for name, m in (("train", train_metrics), ("test", test_metrics)):
        print(f"{name:5s} AUC {m['auc']:.2f} | balanced acc {m['balanced_accuracy']:.2f} | "
              f"sensitivity {m['sensitivity']:.2f} | specificity {m['specificity']:.2f}")
    print(f"written: {out_dir}/selected_genes.csv, sample_scores.csv, metrics.json, figures/")


if __name__ == "__main__":
    main()
