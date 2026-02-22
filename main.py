import argparse
import os
from pathlib import Path

from src.visualization import (
    split_dataset,
    visualize_dataset,
    visualize_dataset_composition,
    visualize_train_test_split,
)
from tools.download_data import download_dataset_balanced, query_ncbi_genome_count

# Optional wandb integration - no-ops if wandb unavailable or init fails
_wandb_run = None


def _wandb_init(config):
    global _wandb_run
    try:
        import wandb  # type: ignore[import-untyped]
        run = wandb.init(project="covid-genome-classification", config=config)
        _wandb_run = run if run is not None else None
        return _wandb_run is not None
    except Exception:
        _wandb_run = None
        return False


def _wandb_log(data):
    if _wandb_run is None:
        return
    try:
        import wandb  # type: ignore[import-untyped]
        wandb.log(data)
    except Exception:
        pass


def _wandb_log_image(key, path):
    if _wandb_run is None or not path or not os.path.isfile(path):
        return
    try:
        import wandb  # type: ignore[import-untyped]
        wandb.log({key: wandb.Image(path)})
    except Exception:
        pass


def _wandb_finish():
    global _wandb_run
    if _wandb_run is None:
        return
    try:
        import wandb  # type: ignore[import-untyped]
        wandb.finish()
    except Exception:
        pass
    finally:
        _wandb_run = None


def main():
    parser = argparse.ArgumentParser(description="Download and visualize SARS-CoV-2 genome dataset")
    parser.add_argument(
        "--full-dataset",
        action="store_true",
        help="Download full dataset (all available accessions). Default: partial dataset (0.05 GB)",
    )
    parser.add_argument(
        "--size-gb",
        type=float,
        default=0.05,
        help="Target dataset size in GB for partial dataset (default: 0.05)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data",
        help="Output directory for downloaded data (default: data)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for deterministic sampling (default: 42)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel download workers (default: 4, max: 16)",
    )
    parser.add_argument(
        "--query-size",
        action="store_true",
        help="Query NCBI for total count and estimated size of all available genomes (human host, complete-only)",
    )

    args = parser.parse_args()

    # Initialize wandb (optional - no-op if unavailable)
    _wandb_init(vars(args))

    # Validate workers parameter (cap at reasonable maximum)
    MAX_WORKERS = 16
    if args.workers < 1:
        print(f"[!] Invalid workers value ({args.workers}). Using minimum: 1")
        args.workers = 1
    elif args.workers > MAX_WORKERS:
        print(f"[!] Workers value ({args.workers}) exceeds maximum. Using: {MAX_WORKERS}")
        args.workers = MAX_WORKERS

    # If query-size flag is set, just query and exit
    if args.query_size:
        result = query_ncbi_genome_count(virus_name="sars-cov-2", host="human", complete_only=True)
        if result:
            print("\n" + "=" * 60)
            print("NCBI SARS-CoV-2 Dataset Information")
            print("=" * 60)
            print(f"Total available genomes (human host, complete-only): {result['total_count']:,}")
            print("\nEstimated dataset size:")
            print(f"  - {result['estimated_size_mb']:,.0f} MB")
            print(f"  - {result['estimated_size_gb']:,.2f} GB")
            print(f"  - {result['estimated_size_tb']:,.2f} TB")
            print("\nNote: This is an estimate based on ~6MB per genome")
            print("      (includes FASTA, GFF, CDS, protein files, and metadata)")
            print("=" * 60)
        _wandb_finish()
        return

    try:
        # Determine dataset size based on toggle
        if args.full_dataset:
            size_gb = None  # None means full dataset
            print("[*] Full dataset mode: downloading all available accessions")
        else:
            size_gb = args.size_gb
            print(f"[*] Partial dataset mode: target size {size_gb} GB")

        # Download dataset (will use Drive if available)
        dataset_path = download_dataset_balanced(
            virus_name="sars-cov-2", output_dir=args.output_dir, size_gb=size_gb, seed=args.seed, workers=args.workers
        )

        # Log download metrics to wandb
        try:
            n_files = len(list(Path(dataset_path).glob("*.zip")))
            _wandb_log({"dataset/total_genomes": n_files, "dataset/path": dataset_path})
        except Exception:
            pass

        # Visualize downloaded data
        visualize_dataset(dataset_path)

        # Log visualization outputs
        _wandb_log_image("download_summary", "download_summary.png")
        for col in ["collection_date", "host", "geo_loc_name"]:
            _wandb_log_image(f"dataset/{col}", f"{col.lower().replace(' ', '_')}_plot.png")

        # Split into train/test datasets (use the same path as download)
        print("\n[*] Splitting dataset into train/test...")
        train_files, test_files = split_dataset(dataset_path, test_ratio=0.2, random_seed=args.seed)

        _wandb_log({"split/train_count": len(train_files), "split/test_count": len(test_files)})

        # Visualize train/test split
        print("\n[*] Creating train/test split visualization...")
        visualize_train_test_split(dataset_path)
        _wandb_log_image("split/train_test_split", "train_test_split.png")

        # Visualize dataset composition
        print("\n[*] Creating detailed dataset composition visualization...")
        visualize_dataset_composition(dataset_path)
        _wandb_log_image("split/dataset_composition", "dataset_composition.png")
    finally:
        _wandb_finish()


if __name__ == "__main__":
    main()
