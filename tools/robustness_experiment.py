"""
Degradation-robustness experiment for the SARS-CoV-2 strain classifier.

Trains a Random Forest once on clean genomes (data/train), then measures test
accuracy as the test genomes are progressively degraded by four mechanisms that
mimic real sequencing failure modes: truncation, N-injection, base
substitutions, and fragment dropout. Produces an accuracy-vs-degradation curve
for each mechanism.

Usage (from project root, with the venv python):
  .venv/Scripts/python tools/robustness_experiment.py
  .venv/Scripts/python tools/robustness_experiment.py --data-dir data --kmer-size 6

Outputs:
  robustness_results.json   structured results (one entry per degradation level)
  robustness_curves.png     2x2 accuracy-vs-severity plot
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

# Allow running as a script from the project root.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.degradation import DEGRADATIONS  # noqa: E402
from src.training import build_kmer_matrix, load_genome_dataset  # noqa: E402


def _accuracy(y_true: list[str], y_pred: list[str]) -> float:
    correct = sum(1 for t, p in zip(y_true, y_pred) if t == p)
    return correct / len(y_true) if y_true else 0.0


def run_experiment(data_dir: str = "data", k: int = 6, seed: int = 42) -> dict:
    from sklearn.ensemble import RandomForestClassifier  # type: ignore

    print(f"[*] Loading genome dataset from {data_dir} ...")
    train_genomes, test_genomes, all_labels = load_genome_dataset(data_dir, k=k)
    if not train_genomes or not test_genomes:
        raise SystemExit("[!] Missing train/ or test/ data. Run main.py to build the split first.")
    print(f"  Train: {len(train_genomes)} genomes, Test: {len(test_genomes)} genomes")
    print(f"  Strains: {sorted(all_labels)}")

    # Fit the classifier once on clean training genomes.
    print(f"\n[*] Training Random Forest on clean genomes (k={k}) ...")
    X_train, y_train, kmer_vocab = build_kmer_matrix(train_genomes, k=k)
    clf = RandomForestClassifier(n_estimators=100, max_depth=20, random_state=seed)
    clf.fit(X_train, y_train)
    print(f"  Feature dimension: {len(kmer_vocab)} k-mers")

    # Clean-test baseline using the frozen vocabulary.
    X_test_clean, y_test, _ = build_kmer_matrix(test_genomes, k=k, kmer_vocab=kmer_vocab)
    baseline_acc = _accuracy(y_test, list(clf.predict(X_test_clean)))
    print(f"  Clean test accuracy (baseline): {baseline_acc:.4f}")

    results: dict = {
        "data_dir": str(Path(data_dir).resolve()),
        "kmer_size": k,
        "seed": seed,
        "test_samples": len(test_genomes),
        "baseline_accuracy": baseline_acc,
        "degradations": {},
    }

    # Sweep each degradation across its severity grid. A fresh RandomState per
    # (degradation, level) keeps every point independent yet reproducible.
    for name, spec in DEGRADATIONS.items():
        fn = spec["fn"]
        param = spec["param"]
        levels = spec["levels"]
        severity_of = spec["severity"]
        print(f"\n[*] Degradation: {name} (varying {param})")
        points = []
        for level in levels:
            rng = np.random.RandomState(seed)
            degraded = []
            for g in test_genomes:
                seq = fn(g["sequence"], level, rng)
                degraded.append({"accession": g["accession"], "sequence": seq, "label": g["label"]})
            X_deg, y_deg, _ = build_kmer_matrix(degraded, k=k, kmer_vocab=kmer_vocab)
            # build_kmer_matrix drops genomes shorter than k or empty; align labels.
            acc = _accuracy(y_deg, list(clf.predict(X_deg))) if X_deg else 0.0
            usable = len(X_deg)
            points.append({
                "level": level,
                "severity": severity_of(level),
                "accuracy": acc,
                "usable_genomes": usable,
            })
            print(f"  {param}={level:<6} severity={severity_of(level):.3f}  acc={acc:.4f}  ({usable}/{len(test_genomes)} usable)")
        results["degradations"][name] = {
            "param": param,
            "xlabel": spec["xlabel"],
            "points": points,
        }

    return results


def plot_results(results: dict, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    degradations = results["degradations"]
    baseline = results["baseline_accuracy"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (name, spec) in zip(axes.flat, degradations.items()):
        pts = spec["points"]
        sev = [p["severity"] for p in pts]
        acc = [p["accuracy"] for p in pts]
        ax.plot(sev, acc, marker="o", linewidth=2, color="#2b6cb0")
        ax.axhline(baseline, linestyle="--", color="gray", alpha=0.7, label=f"clean baseline ({baseline:.3f})")
        ax.set_title(f"{name}")
        ax.set_xlabel(spec["xlabel"] + "  (severity ->)")
        ax.set_ylabel("test accuracy")
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="lower left", fontsize=8)
    # Hide any unused axes if there are fewer than 4 degradations.
    for ax in axes.flat[len(degradations):]:
        ax.axis("off")
    fig.suptitle("Strain classifier robustness to genome degradation", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=120)
    print(f"\n[OK] Saved plot to {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Genome degradation robustness experiment")
    ap.add_argument("--data-dir", "-d", default="data", help="Path to data dir with train/ and test/")
    ap.add_argument("--kmer-size", "-k", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--results-out", default="robustness_results.json")
    ap.add_argument("--plot-out", default="robustness_curves.png")
    args = ap.parse_args()

    results = run_experiment(data_dir=args.data_dir, k=args.kmer_size, seed=args.seed)

    results_path = Path(args.results_out)
    results_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\n[OK] Saved results to {results_path}")

    plot_results(results, Path(args.plot_out))


if __name__ == "__main__":
    main()
