"""
Training module for SARS-CoV-2 genome strain classification.

Supports:
- k-mer + Random Forest
- k-mer + MLP (sklearn)
- wandb integration for monitoring
"""

import zipfile
from collections import Counter
from pathlib import Path

from Bio import SeqIO  # type: ignore[import-untyped]
from Bio.Seq import Seq  # type: ignore[import-untyped]

from src.visualization import get_strain_from_accession

# Optional wandb - no-op if unavailable
def _wandb_log(data, step=None):
    try:
        import wandb  # type: ignore[import-untyped]
        if wandb.run is not None:
            wandb.log(data, step=step)
    except Exception:
        pass


def load_sequence_from_zip(zip_path: Path) -> str | None:
    """
    Extract genomic sequence from an NCBI virus genome zip file.
    Looks for ncbi_dataset/data/genomic.fna or any .fna file.
    """
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            # Try standard path first
            candidates = [
                "ncbi_dataset/data/genomic.fna",
                "genomic.fna",
                "ncbi_dataset/data/genomic.fa",
            ]
            for name in candidates:
                if name in zf.namelist():
                    with zf.open(name) as f:
                        for record in SeqIO.parse(f, "fasta"):
                            return str(record.seq).upper()
                    return None
            # Fallback: find any .fna or .fa file (prefer genomic)
            fa_names = [n for n in zf.namelist() if n.endswith(".fna") or n.endswith(".fa")]
            genomic = [n for n in fa_names if "genomic" in n.lower()]
            for name in (genomic if genomic else fa_names):
                with zf.open(name) as f:
                    for record in SeqIO.parse(f, "fasta"):
                        return str(record.seq).upper()
    except Exception:
        pass
    return None


def extract_kmers(sequence: str, k: int = 6) -> Counter[str]:
    """Extract k-mer counts from a DNA sequence."""
    seq = Seq(sequence)
    # Get canonical k-mers (lexicographically smaller of k-mer and reverse complement)
    kmers = []
    for i in range(len(seq) - k + 1):
        kmer = str(seq[i : i + k])
        if "N" in kmer or len(kmer) != k:
            continue
        rc = str(Seq(kmer).reverse_complement())
        kmers.append(min(kmer, rc))
    return Counter(kmers)


def get_all_kmers(vocab: set[str]) -> list[str]:
    """Return sorted k-mer vocabulary for consistent feature ordering."""
    return sorted(vocab)


def build_kmer_matrix(genomes: list[dict], k: int = 6, kmer_vocab: list[str] | None = None) -> tuple[list[list[float]], list[str], list[str]]:
    """
    Build k-mer count matrix from genome sequences.
    Returns (X, y, kmer_vocab) where X is normalized counts.
    """
    kmer_counts: list[Counter[str]] = []
    labels: list[str] = []

    for item in genomes:
        seq = item["sequence"]
        if not seq or len(seq) < k:
            continue
        kmer_counts.append(extract_kmers(seq, k))
        labels.append(item["label"])

    # Build vocabulary from all k-mers seen
    all_kmers: set[str] = set()
    for c in kmer_counts:
        all_kmers.update(c.keys())

    if kmer_vocab is None:
        kmer_vocab = get_all_kmers(all_kmers)

    # Build matrix (normalize by total k-mer count per genome)
    X = []
    for c in kmer_counts:
        total = sum(c.values())
        row = [c.get(kmer, 0) / total if total > 0 else 0 for kmer in kmer_vocab]
        X.append(row)

    return X, labels, kmer_vocab


def load_genome_dataset(data_dir: str | Path, k: int = 6) -> tuple[list[dict], list[dict], list[str]]:
    """
    Load genomes from train/ and test/ subdirectories.
    Returns (train_genomes, test_genomes, strain_labels) where each genome dict has
    'accession', 'sequence', 'label'.
    """
    data_path = Path(data_dir)
    train_dir = data_path / "train"
    test_dir = data_path / "test"

    def load_from_dir(d: Path) -> list[dict]:
        genomes = []
        for zip_path in sorted(d.glob("*.zip")):
            accession = zip_path.stem
            label = get_strain_from_accession(accession)
            seq = load_sequence_from_zip(zip_path)
            if seq:
                genomes.append({"accession": accession, "sequence": seq, "label": label})
        return genomes

    train_genomes = load_from_dir(train_dir) if train_dir.exists() else []
    test_genomes = load_from_dir(test_dir) if test_dir.exists() else []
    return train_genomes, test_genomes, list(set(g["label"] for g in train_genomes + test_genomes))


def train_random_forest(
    X_train: list[list[float]],
    y_train: list[str],
    X_test: list[list[float]],
    y_test: list[str],
    n_estimators: int = 100,
    max_depth: int | None = 20,
    random_state: int = 42,
    use_wandb: bool = True,
):
    """Train Random Forest classifier and return fitted model + test accuracy."""
    from sklearn.ensemble import RandomForestClassifier  # type: ignore
    from sklearn.metrics import accuracy_score, classification_report  # type: ignore

    clf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        random_state=random_state,
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    acc = accuracy_score(y_test, y_pred)

    if use_wandb:
        _wandb_log({
            "model": "random_forest",
            "test_accuracy": acc,
            "train_samples": len(y_train),
            "test_samples": len(y_test),
        })

    print("\n[✓] Random Forest training complete")
    print(f"  Test accuracy: {acc:.4f}")
    print(classification_report(y_test, y_pred))
    return clf, acc


def train_mlp(
    X_train: list[list[float]],
    y_train: list[str],
    X_test: list[list[float]],
    y_test: list[str],
    label_encoder: "dict[str, int]",
    hidden_layer_sizes: tuple[int, ...] = (128, 64),
    max_iter: int = 200,
    random_state: int = 42,
    use_wandb: bool = True,
):
    """Train MLP classifier (sklearn) and return fitted model + test accuracy."""
    from sklearn.metrics import accuracy_score, classification_report  # type: ignore
    from sklearn.neural_network import MLPClassifier  # type: ignore

    y_train_int = [label_encoder[l] for l in y_train]
    y_test_int = [label_encoder[l] for l in y_test]

    clf = MLPClassifier(
        hidden_layer_sizes=hidden_layer_sizes,
        max_iter=max_iter,
        random_state=random_state,
        early_stopping=True,
        validation_fraction=0.1,
    )
    inv_encoder = {v: k for k, v in label_encoder.items()}
    clf.fit(X_train, y_train_int)

    # Log MLP training and validation curves to wandb
    if use_wandb and hasattr(clf, "loss_curve_"):
        for step, loss in enumerate(clf.loss_curve_):
            log_data = {"train/loss": loss}
            if hasattr(clf, "validation_scores_") and clf.validation_scores_ is not None and step < len(clf.validation_scores_):
                log_data["val/accuracy"] = clf.validation_scores_[step]
            _wandb_log(log_data, step=step)

    y_pred_int = clf.predict(X_test)
    y_pred = [inv_encoder[v] for v in y_pred_int]
    acc = accuracy_score(y_test_int, y_pred_int)

    if use_wandb:
        _wandb_log({
            "model": "mlp",
            "test_accuracy": acc,
            "train_samples": len(y_train),
            "test_samples": len(y_test),
        })

    print("\n[✓] MLP training complete")
    print(f"  Test accuracy: {acc:.4f}")
    print(classification_report(y_test, y_pred))
    return clf, acc


def run_training(
    data_dir: str | Path,
    model: str = "random_forest",
    k: int = 6,
    use_wandb: bool = True,
    **kwargs,
):
    """
    Main training entry point. Loads data, trains model, evaluates.
    """
    print(f"\n[*] Loading genome dataset from {data_dir}...")
    train_genomes, test_genomes, all_labels = load_genome_dataset(data_dir, k=k)

    if not train_genomes:
        print("[!] No training data found. Run main.py first to download and split data.")
        return None
    if not test_genomes:
        print("[!] No test data found. Run main.py first to download and split data.")
        return None

    print(f"  Train: {len(train_genomes)} genomes, Test: {len(test_genomes)} genomes")
    print(f"  Strains: {sorted(all_labels)}")

    # Build k-mer features
    print(f"\n[*] Extracting {k}-mers and building feature matrix...")
    X_train, y_train, kmer_vocab = build_kmer_matrix(train_genomes, k=k)
    X_test, y_test, _ = build_kmer_matrix(test_genomes, k=k, kmer_vocab=kmer_vocab)

    print(f"  Feature dimension: {len(kmer_vocab)} k-mers")
    print(f"  Train X: {len(X_train)} samples, Test X: {len(X_test)} samples")

    # Label encoding for MLP
    label_encoder = {l: i for i, l in enumerate(sorted(all_labels))}

    if use_wandb:
        try:
            import wandb  # type: ignore[import-untyped]
            if wandb.run is None:
                wandb.init(project="covid-genome-classification", config={
                    "model": model,
                    "k": k,
                    "train_samples": len(train_genomes),
                    "test_samples": len(test_genomes),
                    "kmer_vocab_size": len(kmer_vocab),
                    **kwargs,
                })
        except Exception:
            use_wandb = False

    if model == "random_forest":
        clf, acc = train_random_forest(X_train, y_train, X_test, y_test, use_wandb=use_wandb, **kwargs)
    elif model == "mlp":
        clf, acc = train_mlp(
            X_train, y_train, X_test, y_test, label_encoder,
            use_wandb=use_wandb, **kwargs
        )
    else:
        raise ValueError(f"Unknown model: {model}. Use 'random_forest' or 'mlp'.")

    if use_wandb:
        try:
            import wandb  # type: ignore[import-untyped]
            wandb.finish()
        except Exception:
            pass

    return clf


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Train strain classifier on existing data")
    ap.add_argument("--data-dir", "-d", default="data", help="Path to data dir with train/ and test/")
    ap.add_argument("--model", "-m", default="random_forest", choices=["random_forest", "mlp"])
    ap.add_argument("--kmer-size", "-k", type=int, default=6)
    ap.add_argument("--no-wandb", action="store_true", help="Disable wandb logging")
    args = ap.parse_args()
    run_training(
        data_dir=args.data_dir,
        model=args.model,
        k=args.kmer_size,
        use_wandb=not args.no_wandb,
    )
