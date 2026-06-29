"""
SPIKE (throwaway de-risk experiment) — gate before building the full robustness study.

Tests two load-bearing assumptions cheaply, on ONE corruption type (substitution
noise) at THREE severities, for BOTH feature representations:

  1. STRUCTURE: does macro-F1 decline as a structured slope (a finding) rather
     than a cliff or a flat line (not a finding)?
  2. DIVERGENCE: do the two representations -- k-mer frequencies vs. mutation
     (substitution-site) features -- fail DIFFERENTLY? The mechanistic prediction
     is that substitution noise shreds k-mers (one SNP corrupts k overlapping
     k-mers) but barely touches mutation features (spurious mutations scatter
     across ~30kb and miss the small informative site set).

If we see slope + divergence -> green light for the full harness.
If we see a cliff or glued-together lines -> stop and rethink before investing.

This is deliberately minimal and not production code. Run from project root:
  .venv/Scripts/python tools/spike_representation_robustness.py
"""

from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.degradation import add_substitutions  # noqa: E402
from src.training import build_kmer_matrix, load_genome_dataset  # noqa: E402

# Restored Nextclade binary + dataset live in the Windows temp dir.
_TEMP = Path(tempfile.gettempdir())
NEXTCLADE_BIN = _TEMP / "nextclade.exe"
NEXTCLADE_DATASET = _TEMP / "nextclade" / "sars-cov-2-dataset"

SUB_LEVELS = [0.0, 0.01, 0.05]   # substitution rates: clean, mild, severe
K = 6
SEED = 42
MIN_SITE_COUNT = 2  # keep substitution sites seen in >= this many training genomes


def macro_f1(y_true: list[str], y_pred: list[str]) -> float:
    from sklearn.metrics import f1_score  # type: ignore
    return float(f1_score(y_true, y_pred, average="macro", zero_division=0))


# ----------------------------- mutation features -----------------------------

def _write_fasta(genomes: list[dict], path: Path) -> None:
    with path.open("w", encoding="ascii") as fh:
        for g in genomes:
            fh.write(f">{g['accession']}\n{g['sequence']}\n")


def run_nextclade(genomes: list[dict], tag: str, workdir: Path) -> dict[str, set[str]]:
    """
    Run Nextclade on the given genomes and return {accession -> set of
    nucleotide substitutions like 'C241T'}. Missing/failed sequences map to an
    empty set.
    """
    fasta = workdir / f"{tag}.fasta"
    out_tsv = workdir / f"{tag}.tsv"
    _write_fasta(genomes, fasta)
    cmd = [
        str(NEXTCLADE_BIN), "run",
        "--input-dataset", str(NEXTCLADE_DATASET),
        "--output-tsv", str(out_tsv),
        "--silent",
        str(fasta),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    subs: dict[str, set[str]] = {g["accession"]: set() for g in genomes}
    with out_tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = row.get("seqName", "").strip()
            raw = (row.get("substitutions") or "").strip()
            if name in subs and raw:
                subs[name] = {tok for tok in raw.split(",") if tok}
    return subs


def build_mutation_matrix(
    subs_by_acc: dict[str, set[str]],
    order: list[str],
    vocab: list[str] | None = None,
    min_count: int = MIN_SITE_COUNT,
) -> tuple[list[list[int]], list[str]]:
    """Binary substitution-site features, vocab fixed from the training set."""
    if vocab is None:
        counts: dict[str, int] = {}
        for acc in order:
            for s in subs_by_acc.get(acc, set()):
                counts[s] = counts.get(s, 0) + 1
        vocab = sorted(s for s, c in counts.items() if c >= min_count)
    index = {s: i for i, s in enumerate(vocab)}
    X = []
    for acc in order:
        row = [0] * len(vocab)
        for s in subs_by_acc.get(acc, set()):
            j = index.get(s)
            if j is not None:
                row[j] = 1
        X.append(row)
    return X, vocab


def main() -> None:
    if not NEXTCLADE_BIN.exists() or not NEXTCLADE_DATASET.exists():
        raise SystemExit(f"[!] Nextclade not found at {NEXTCLADE_BIN} / {NEXTCLADE_DATASET}")

    from sklearn.ensemble import RandomForestClassifier  # type: ignore

    print("[*] Loading genomes...")
    train_genomes, test_genomes, labels = load_genome_dataset("data", k=K)
    print(f"  Train: {len(train_genomes)}, Test: {len(test_genomes)}, strains: {sorted(labels)}")
    y_train = [g["label"] for g in train_genomes]
    y_test = [g["label"] for g in test_genomes]

    workdir = Path(tempfile.mkdtemp(prefix="spike_"))
    print(f"  Scratch dir: {workdir}")

    # ---- Representation 1: k-mers ----
    print("\n[*] k-mer representation: training on clean genomes...")
    Xtr_k, ytr_k, kvocab = build_kmer_matrix(train_genomes, k=K)
    rf_k = RandomForestClassifier(n_estimators=100, max_depth=20, random_state=SEED)
    rf_k.fit(Xtr_k, ytr_k)

    # ---- Representation 2: mutation features (Nextclade on CLEAN train) ----
    print("[*] mutation representation: running Nextclade on clean train genomes...")
    train_subs = run_nextclade(train_genomes, "train_clean", workdir)
    train_order = [g["accession"] for g in train_genomes]
    Xtr_m, mvocab = build_mutation_matrix(train_subs, train_order)
    print(f"  mutation feature dimension: {len(mvocab)} substitution sites")
    rf_m = RandomForestClassifier(n_estimators=100, max_depth=20, random_state=SEED)
    rf_m.fit(Xtr_m, y_train)

    test_order = [g["accession"] for g in test_genomes]

    print("\n[*] Sweeping substitution noise on the TEST set:\n")
    print(f"  {'sub_rate':>9} | {'kmer_F1':>8} | {'mut_F1':>8} | {'gap':>6}")
    print("  " + "-" * 40)
    results = []
    for level in SUB_LEVELS:
        rng = np.random.RandomState(SEED)
        degraded = [
            {"accession": g["accession"], "label": g["label"],
             "sequence": add_substitutions(g["sequence"], level, rng)}
            for g in test_genomes
        ]
        # k-mer eval (frozen training vocab)
        Xk, _, _ = build_kmer_matrix(degraded, k=K, kmer_vocab=kvocab)
        f1_k = macro_f1(y_test, list(rf_k.predict(Xk)))
        # mutation eval (fresh Nextclade on the degraded sequences, frozen vocab)
        deg_subs = run_nextclade(degraded, f"test_sub_{level}", workdir)
        Xm, _ = build_mutation_matrix(deg_subs, test_order, vocab=mvocab)
        f1_m = macro_f1(y_test, list(rf_m.predict(Xm)))
        gap = f1_k - f1_m
        results.append((level, f1_k, f1_m, gap))
        print(f"  {level:>9.3f} | {f1_k:>8.4f} | {f1_m:>8.4f} | {gap:>+6.3f}")

    # ---- Verdict ----
    f1k = [r[1] for r in results]
    f1m = [r[2] for r in results]
    kmer_drop = f1k[0] - f1k[-1]
    mut_drop = f1m[0] - f1m[-1]
    max_gap = max(abs(r[3]) for r in results)

    print("\n[*] SPIKE VERDICT")
    print(f"  k-mer macro-F1 drop (clean->severe): {kmer_drop:+.3f}")
    print(f"  mutation macro-F1 drop (clean->severe): {mut_drop:+.3f}")
    print(f"  max divergence between representations: {max_gap:.3f}")
    structured = kmer_drop > 0.05 and f1k[1] < f1k[0]  # mid point between clean and severe
    diverges = max_gap > 0.05
    print(f"  STRUCTURE (k-mers decline with a mid-point, not a cliff): {'YES' if structured else 'NO'}")
    print(f"  DIVERGENCE (representations differ by >0.05 F1): {'YES' if diverges else 'NO'}")
    if structured and diverges:
        print("  => GREEN LIGHT: build the full harness.")
    else:
        print("  => CAUTION: assumption(s) not met; rethink before the full build.")


if __name__ == "__main__":
    main()
