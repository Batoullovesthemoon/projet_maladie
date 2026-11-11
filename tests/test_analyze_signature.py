"""Automated tests of scripts/analyze_signature.py on a tiny synthetic dataset."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_signature as sig  # noqa: E402


@pytest.fixture
def toy():
    """60 genes x 30 samples (20 CKD, 10 Control). 5 genes are higher in CKD,
    5 are lower, 50 are pure noise."""
    rng = np.random.default_rng(0)
    samples = [f"S{i:02d}" for i in range(30)]
    y = pd.Series([1] * 20 + [0] * 10, index=samples)
    genes = ([f"UP{i}" for i in range(5)] + [f"DOWN{i}" for i in range(5)]
             + [f"N{i}" for i in range(50)])
    expr = pd.DataFrame(rng.normal(8, 0.3, (len(genes), len(samples))),
                        index=genes, columns=samples)
    expr.loc[[g for g in genes if g.startswith("UP")], y.index[y == 1]] += 2.0
    expr.loc[[g for g in genes if g.startswith("DOWN")], y.index[y == 1]] -= 2.0
    return expr, y


def write_toy(folder, expr, y, matrix_name="matrix.tsv"):
    expr.to_csv(folder / matrix_name, sep="\t", index_label="Gene")
    labels = pd.DataFrame({"gsm": y.index, "title": y.index, "study_set": "discovery",
                           "condition": np.where(y == 1, "CKD", "Control")})
    labels.to_csv(folder / "labels.tsv", sep="\t", index=False)


def run_script(matrix, labels, out_dir):
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "analyze_signature.py"),
         "--matrix", str(matrix), "--labels", str(labels),
         "--out-dir", str(out_dir), "--top-n", "3"],
        capture_output=True, text=True)


def test_welch_bh_ranks_planted_genes_first(toy):
    expr, y = toy
    res = sig.welch_bh(expr, y)
    ups = [g for g in expr.index if g.startswith("UP")]
    downs = [g for g in expr.index if g.startswith("DOWN")]
    noise = [g for g in expr.index if g.startswith("N")]
    assert (res.loc[ups, "log2fc"] > 1).all()
    assert (res.loc[downs, "log2fc"] < -1).all()
    assert res.loc[ups + downs, "padj"].max() < res.loc[noise, "padj"].min()


def test_select_genes_takes_top_up_and_down(toy):
    expr, y = toy
    up, down = sig.select_genes(sig.welch_bh(expr, y), top_n=3, min_lfc=0.5)
    assert len(up) == 3 and all(g.startswith("UP") for g in up)
    assert len(down) == 3 and all(g.startswith("DOWN") for g in down)


def test_select_genes_ignores_small_fold_change():
    res = pd.DataFrame({"log2fc": [0.1, -0.1, 2.0, -2.0],
                        "pvalue": [1e-9, 1e-9, 1e-6, 1e-6],
                        "padj": [1e-9, 1e-9, 1e-6, 1e-6]},
                       index=["g1", "g2", "g3", "g4"])
    up, down = sig.select_genes(res, top_n=5, min_lfc=0.5)
    assert up == ["g3"] and down == ["g4"]


def test_signature_score_uses_given_train_stats():
    expr = pd.DataFrame({"s1": [3.0, 2.0], "s2": [5.0, 0.0]}, index=["U", "D"])
    mu = pd.Series({"U": 4.0, "D": 1.0})
    sd = pd.Series({"U": 1.0, "D": 2.0})
    score = sig.signature_score(expr, ["U"], ["D"], mu, sd)
    assert np.allclose(score.values, [-1.5, 1.5])


def test_youden_threshold_separates_two_groups():
    y = np.array([0, 0, 0, 1, 1, 1])
    score = np.array([-2.0, -1.0, -1.5, 1.0, 2.0, 3.0])
    thr = sig.youden_threshold(y, score)
    assert (score[y == 1] >= thr).all()
    assert (score[y == 0] < thr).all()


def test_evaluate_known_values():
    y = np.array([1, 1, 1, 1, 0, 0])
    score = np.array([0.9, 0.8, 0.7, 0.1, 0.2, -0.5])
    m = sig.evaluate(y, score, 0.5)
    assert (m["tp"], m["fn"], m["tn"], m["fp"]) == (3, 1, 2, 0)
    assert m["sensitivity"] == pytest.approx(0.75)
    assert m["specificity"] == pytest.approx(1.0)
    assert m["balanced_accuracy"] == pytest.approx(0.875)
    assert m["auc"] == pytest.approx(0.875)


def test_script_end_to_end(tmp_path, toy):
    expr, y = toy
    write_toy(tmp_path, expr, y)
    r = run_script(tmp_path / "matrix.tsv", tmp_path / "labels.tsv", tmp_path / "out")
    assert r.returncode == 0, r.stderr
    for name in ("selected_genes.csv", "sample_scores.csv", "metrics.json",
                 "figures/score_boxplot.png", "figures/roc_curve.png"):
        assert (tmp_path / "out" / name).is_file(), name
    metrics = json.loads((tmp_path / "out" / "metrics.json").read_text())
    assert metrics["test"]["auc"] > 0.9


def test_script_fails_clearly_on_missing_file(tmp_path):
    r = run_script(tmp_path / "nope.tsv", tmp_path / "nope_labels.tsv", tmp_path / "out")
    assert r.returncode == 1
    assert "ERREUR" in r.stderr


def test_test_samples_do_not_change_selection(tmp_path, toy):
    """No leakage: scrambling the TEST samples must not change genes or threshold."""
    expr, y = toy
    write_toy(tmp_path, expr, y)
    out1 = tmp_path / "out1"
    assert run_script(tmp_path / "matrix.tsv", tmp_path / "labels.tsv", out1).returncode == 0
    scores = pd.read_csv(out1 / "sample_scores.csv")
    test_ids = scores.loc[scores["split"] == "test", "gsm"].tolist()
    assert test_ids

    expr2 = expr.copy()
    expr2[test_ids] = np.random.default_rng(1).normal(8, 3, (expr.shape[0], len(test_ids)))
    write_toy(tmp_path, expr2, y, matrix_name="matrix2.tsv")
    out2 = tmp_path / "out2"
    assert run_script(tmp_path / "matrix2.tsv", tmp_path / "labels.tsv", out2).returncode == 0

    pd.testing.assert_frame_equal(pd.read_csv(out1 / "selected_genes.csv"),
                                  pd.read_csv(out2 / "selected_genes.csv"))
    t1 = json.loads((out1 / "metrics.json").read_text())["threshold"]
    t2 = json.loads((out2 / "metrics.json").read_text())["threshold"]
    assert t1 == pytest.approx(t2)
