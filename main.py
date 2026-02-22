import argparse

from src.training import run_training
from src.visualization import (
    split_dataset,
    visualize_dataset,
    visualize_dataset_composition,
    visualize_train_test_split,
)
from tools.download_data import download_dataset_balanced, query_ncbi_genome_count


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
    parser.add_argument(
        "--train",
        action="store_true",
        help="Run strain classification training after data pipeline (k-mer + Random Forest or MLP)",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="random_forest",
        choices=["random_forest", "mlp"],
        help="Model for training: random_forest or mlp (default: random_forest)",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=6,
        help="k-mer size for sequence features (default: 6)",
    )

    args = parser.parse_args()

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
        return

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

    # Visualize downloaded data
    visualize_dataset(dataset_path)

    # Split into train/test datasets (use the same path as download)
    print("\n[*] Splitting dataset into train/test...")
    train_files, test_files = split_dataset(dataset_path, test_ratio=0.2, random_seed=args.seed)

    # Visualize train/test split
    print("\n[*] Creating train/test split visualization...")
    visualize_train_test_split(dataset_path)

    # Visualize dataset composition
    print("\n[*] Creating detailed dataset composition visualization...")
    visualize_dataset_composition(dataset_path)

    # Optional training
    if args.train:
        print("\n" + "=" * 60)
        print("Training strain classification model")
        print("=" * 60)
        run_training(
            data_dir=dataset_path,
            model=args.model,
            k=args.kmer_size,
            use_wandb=True,
        )


if __name__ == "__main__":
    main()
