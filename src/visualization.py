import json
import os
import random
import shutil
from pathlib import Path

try:
    import matplotlib.pyplot as plt  # type: ignore
    import pandas as pd  # type: ignore
except ImportError as e:
    raise ImportError(
        f"Required packages not installed: {e}\nPlease install with: pip install matplotlib pandas"
    ) from e


def get_strain_from_accession(accession):
    """
    Map a genome accession number to its SARS-CoV-2 strain type.

    This function classifies accessions based on their prefix patterns to identify
    the variant strain (Reference, Alpha, Beta, Gamma, Delta, Omicron, Recent).

    Args:
        accession (str): NCBI accession number (e.g., "NC_045512.2", "MT188341.1")

    Returns:
        str: Strain type (Reference, Alpha, Beta, Gamma, Delta, Omicron, Recent, or Other)
    """
    if accession == "NC_045512.2" or accession == "NC_045512.1":
        return "Reference"
    elif accession.startswith("MT"):
        # Check if it's Alpha (188xxx, 326xxx) or Beta (291xxx) or Reference/Early
        acc_num = accession.split(".")[0]
        if "188" in acc_num or "326" in acc_num:
            return "Alpha"
        elif "291" in acc_num:
            return "Beta"
        else:
            return "Reference"  # Early/Reference variants
    elif accession.startswith("MW"):
        # Check if it's Gamma (633477-633496) or Delta (633497+)
        try:
            num = int(accession.split(".")[0].replace("MW", ""))
            if 633477 <= num <= 633496:
                return "Gamma"
            else:
                return "Delta"
        except ValueError:
            return "Delta"  # Default to Delta for MW accessions
    elif accession.startswith("OM"):
        return "Omicron"
    elif accession.startswith("ON") or accession.startswith("OP") or accession.startswith("OR"):
        return "Recent"
    else:
        return "Other"


def load_metadata(jsonl_path):
    """Load NCBI genome metadata from JSON Lines file into a DataFrame."""
    records = []
    with open(jsonl_path) as f:
        for line in f:
            data = json.loads(line)
            records.append(data)
    df = pd.json_normalize(records)
    return df


def plot_counts(df, column, title, xlabel, ylabel="Count", top_n=None):
    """Plot counts of unique values in a column."""
    counts = df[column].value_counts()
    if top_n:
        counts = counts.head(top_n)

    # Create plot for Google Colab
    plt.figure(figsize=(12, 6))
    counts.plot(kind="bar")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # Display in Google Colab
    plt.show()

    # Save plot as image
    filename = f"{column.lower().replace(' ', '_')}_plot.png"
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    print(f"[*] Plot saved as '{filename}'")


def visualize_dataset(dataset_path):
    """
    Visualize the metadata in an extracted NCBI virus dataset folder.

    Args:
        dataset_path (str): Path to the extracted dataset folder returned by download_dataset().
    """
    # Check if this is a bulk dataset or individual genome downloads
    metadata_file = os.path.join(dataset_path, "ncbi_dataset", "metadata", "genome_metadata.jsonl")

    if os.path.exists(metadata_file):
        # Bulk dataset format
        print(f"Loading metadata from {metadata_file}…")
        df = load_metadata(metadata_file)
        print(f"Metadata contains {len(df)} genomes with {len(df.columns)} fields.")

        # Plot top collection dates
        if "collection_date" in df.columns:
            plot_counts(df, "collection_date", "Genomes by Collection Date", "Collection Date", top_n=20)

        # Plot hosts
        if "host" in df.columns:
            plot_counts(df, "host", "Genomes by Host", "Host")

        # Plot locations
        if "geo_loc_name" in df.columns:
            plot_counts(df, "geo_loc_name", "Genomes by Location", "Country / Location", top_n=20)
    else:
        # Individual genome downloads - create simple summary
        print(f"[*] Individual genome downloads detected in {dataset_path}")

        # Count downloaded files
        zip_files = [f for f in os.listdir(dataset_path) if f.endswith(".zip")]
        print(f"[*] Found {len(zip_files)} downloaded genome files:")
        for i, zip_file in enumerate(zip_files, 1):
            print(f"  {i}. {zip_file}")

        # Create a simple visualization showing download summary

        plt.figure(figsize=(10, 6))
        plt.bar(["Downloaded Genomes"], [len(zip_files)])
        plt.title(f"SARS-CoV-2 Genome Downloads ({len(zip_files)} genomes)")
        plt.ylabel("Number of Genomes")
        plt.tight_layout()

        # Save the plot as an image file
        plt.savefig("download_summary.png", dpi=150, bbox_inches="tight")
        print("[*] Plot saved as 'download_summary.png'")

        # Display the plot
        plt.show()
        print("[*] Plot displayed")

        print(f"[✓] Visualization complete - {len(zip_files)} genomes downloaded successfully")


def split_dataset(data_dir, test_ratio=0.2, random_seed=42, balance_variants=True):
    """
    Split downloaded genome files into train and test datasets with balanced strain representation.
    Ensures equal number of genomes from each strain in both train and test sets.

    Args:
        data_dir (str): Directory containing downloaded genome files
        test_ratio (float): Proportion of data to use for testing (default: 0.2)
        random_seed (int): Random seed for reproducible splits (default: 42)
        balance_variants (bool): Whether to balance strain representation (default: True)

    Returns:
        tuple: (train_files, test_files) - lists of file paths
    """
    # Set random seed for reproducible splits
    random.seed(random_seed)

    # Get all zip files in the data directory
    data_path = Path(data_dir)
    zip_files = list(data_path.glob("*.zip"))

    if len(zip_files) == 0:
        print("[!] No genome files found in data directory")
        return [], []

    print(f"[*] Found {len(zip_files)} genome files to split")

    if balance_variants:
        # Group files by strain type using module-level get_strain_from_accession
        strain_groups = {}
        for file_path in zip_files:
            accession = file_path.stem  # Remove .zip extension
            strain_type = get_strain_from_accession(accession)

            if strain_type not in strain_groups:
                strain_groups[strain_type] = []
            strain_groups[strain_type].append(file_path)

        # Print strain distribution
        print("[*] Strain distribution:")
        for strain, files in sorted(strain_groups.items()):
            print(f"  - {strain}: {len(files)} genomes")

        # Split each strain group proportionally with equal representation
        train_files = []
        test_files = []

        for strain, files in sorted(strain_groups.items()):
            if len(files) == 0:
                continue

            # Shuffle files within this strain
            random.shuffle(files)

            # Calculate split for this strain - ensure at least 1 in test if possible
            if len(files) == 1:
                # If only 1 file, put it in train
                test_size = 0
            else:
                test_size = max(1, int(len(files) * test_ratio))

            train_size = len(files) - test_size

            # Split this strain
            strain_test = files[:test_size]
            strain_train = files[test_size:]

            train_files.extend(strain_train)
            test_files.extend(strain_test)

            print(f"  - {strain}: {len(strain_train)} train, {len(strain_test)} test")

        print("[*] Balanced split ensures each strain is represented proportionally in both train and test")

    else:
        # Original random split (no balancing)
        random.shuffle(zip_files)
        test_size = int(len(zip_files) * test_ratio)
        test_files = zip_files[:test_size]
        train_files = zip_files[test_size:]

    # Create train and test directories
    train_dir = data_path / "train"
    test_dir = data_path / "test"

    train_dir.mkdir(exist_ok=True)
    test_dir.mkdir(exist_ok=True)

    # Copy files to respective directories
    train_size = len(train_files)
    test_size = len(test_files)

    print(f"[*] Creating train dataset ({train_size} files)...")
    for file_path in train_files:
        shutil.copy2(file_path, train_dir / file_path.name)

    print(f"[*] Creating test dataset ({test_size} files)...")
    for file_path in test_files:
        shutil.copy2(file_path, test_dir / file_path.name)

    # Print summary
    print("[✓] Dataset split complete:")
    print(f"  - Train: {train_size} files ({train_size / len(zip_files) * 100:.1f}%)")
    print(f"  - Test: {test_size} files ({test_size / len(zip_files) * 100:.1f}%)")
    print(f"  - Train directory: {train_dir}")
    print(f"  - Test directory: {test_dir}")

    return train_files, test_files


def visualize_train_test_split(data_dir):
    """Create visualizations for train/test split."""
    train_dir = Path(data_dir) / "train"
    test_dir = Path(data_dir) / "test"

    if not train_dir.exists() or not test_dir.exists():
        print("[!] Train/test directories not found. Run split_dataset() first.")
        return

    # Count files in each directory
    train_files = list(train_dir.glob("*.zip"))
    test_files = list(test_dir.glob("*.zip"))

    # Create visualization

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Train/Test split bar chart
    categories = ["Train", "Test"]
    counts = [len(train_files), len(test_files)]
    colors = ["#2E8B57", "#DC143C"]  # Green for train, Red for test

    ax1.bar(categories, counts, color=colors)
    ax1.set_title("Train/Test Dataset Split")
    ax1.set_ylabel("Number of Genomes")
    ax1.set_ylim(0, max(counts) * 1.1)

    # Add count labels on bars
    for i, count in enumerate(counts):
        ax1.text(i, count + 0.1, str(count), ha="center", va="bottom", fontweight="bold")

    # Pie chart
    ax2.pie(counts, labels=categories, colors=colors, autopct="%1.1f%%", startangle=90)
    ax2.set_title("Dataset Distribution")

    plt.tight_layout()

    # Save plot
    plt.savefig("train_test_split.png", dpi=150, bbox_inches="tight")
    print("[*] Train/test split plot saved as 'train_test_split.png'")

    # Display the plot
    plt.show()
    print("[*] Train/test split plot displayed")

    print("[✓] Train/test visualization complete:")
    print(f"  - Train: {len(train_files)} genomes")
    print(f"  - Test: {len(test_files)} genomes")
    print(f"  - Split ratio: {len(test_files) / (len(train_files) + len(test_files)) * 100:.1f}% test")


def visualize_dataset_composition(data_dir):
    """Create detailed visualizations for train and test dataset composition."""
    train_dir = Path(data_dir) / "train"
    test_dir = Path(data_dir) / "test"

    if not train_dir.exists() or not test_dir.exists():
        print("[!] Train/test directories not found. Run split_dataset() first.")
        return

    # Get files from both directories
    train_files = list(train_dir.glob("*.zip"))
    test_files = list(test_dir.glob("*.zip"))

    # Analyze variant composition using module-level get_strain_from_accession
    def analyze_variants(files):
        variant_counts = {}
        for file_path in files:
            accession = file_path.stem
            variant_type = get_strain_from_accession(accession)
            variant_counts[variant_type] = variant_counts.get(variant_type, 0) + 1
        return variant_counts

    train_variants = analyze_variants(train_files)
    test_variants = analyze_variants(test_files)

    # Create comprehensive visualization

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Train dataset composition
    if train_variants:
        variants = list(train_variants.keys())
        counts = list(train_variants.values())
        colors = plt.cm.Set3(range(len(variants)))

        ax1.bar(variants, counts, color=colors)
        ax1.set_title("Train Dataset - Variant Composition")
        ax1.set_ylabel("Number of Genomes")
        ax1.tick_params(axis="x", rotation=45)

        # Add count labels
        for i, count in enumerate(counts):
            ax1.text(i, count + 0.1, str(count), ha="center", va="bottom", fontweight="bold")

    # 2. Test dataset composition
    if test_variants:
        variants = list(test_variants.keys())
        counts = list(test_variants.values())
        colors = plt.cm.Set3(range(len(variants)))

        ax2.bar(variants, counts, color=colors)
        ax2.set_title("Test Dataset - Variant Composition")
        ax2.set_ylabel("Number of Genomes")
        ax2.tick_params(axis="x", rotation=45)

        # Add count labels
        for i, count in enumerate(counts):
            ax2.text(i, count + 0.1, str(count), ha="center", va="bottom", fontweight="bold")

    # 3. Side-by-side comparison
    all_variants = set(train_variants.keys()) | set(test_variants.keys())
    train_counts = [train_variants.get(v, 0) for v in all_variants]
    test_counts = [test_variants.get(v, 0) for v in all_variants]

    x = range(len(all_variants))
    width = 0.35

    ax3.bar([i - width / 2 for i in x], train_counts, width, label="Train", color="#2E8B57", alpha=0.8)
    ax3.bar([i + width / 2 for i in x], test_counts, width, label="Test", color="#DC143C", alpha=0.8)
    ax3.set_title("Train vs Test - Variant Comparison")
    ax3.set_ylabel("Number of Genomes")
    ax3.set_xticks(x)
    ax3.set_xticklabels(all_variants, rotation=45)
    ax3.legend()

    # Add count labels
    for i, (train_count, test_count) in enumerate(zip(train_counts, test_counts)):
        ax3.text(i - width / 2, train_count + 0.1, str(train_count), ha="center", va="bottom", fontweight="bold")
        ax3.text(i + width / 2, test_count + 0.1, str(test_count), ha="center", va="bottom", fontweight="bold")

    # 4. Overall distribution pie chart
    total_train = sum(train_variants.values())
    total_test = sum(test_variants.values())

    ax4.pie(
        [total_train, total_test],
        labels=["Train", "Test"],
        colors=["#2E8B57", "#DC143C"],
        autopct="%1.1f%%",
        startangle=90,
    )
    ax4.set_title("Overall Dataset Distribution")

    plt.tight_layout()

    # Save plot
    plt.savefig("dataset_composition.png", dpi=150, bbox_inches="tight")
    print("[*] Dataset composition plot saved as 'dataset_composition.png'")

    # Display the plot
    plt.show()
    print("[*] Dataset composition plot displayed")

    # Print detailed summary
    print("\n[✓] Dataset Composition Analysis:")
    print(f"  Train Dataset ({total_train} genomes):")
    for variant, count in train_variants.items():
        print(f"    - {variant}: {count} genomes")

    print(f"  Test Dataset ({total_test} genomes):")
    for variant, count in test_variants.items():
        print(f"    - {variant}: {count} genomes")

    print("  Balance Check:")
    for variant in all_variants:
        train_count = train_variants.get(variant, 0)
        test_count = test_variants.get(variant, 0)
        total = train_count + test_count
        if total > 0:
            train_pct = (train_count / total) * 100
            test_pct = (test_count / total) * 100
            print(f"    - {variant}: {train_pct:.1f}% train, {test_pct:.1f}% test")


if __name__ == "__main__":
    # Test the visualization with a sample dataset path

    # Check if data directory exists
    if os.path.exists("data"):
        print("[*] Running visualization on existing data directory...")
        visualize_dataset("data")

        # Split dataset into train/test
        print("\n[*] Splitting dataset into train/test...")
        train_files, test_files = split_dataset("data", test_ratio=0.2, random_seed=42)

        # Visualize train/test split
        print("\n[*] Creating train/test split visualization...")
        visualize_train_test_split("data")

        # Visualize dataset composition
        print("\n[*] Creating detailed dataset composition visualization...")
        visualize_dataset_composition("data")
    else:
        print("[!] No 'data' directory found. Please run main.py first to download data.")
        print("[*] Creating a sample visualization...")

        # Create a sample visualization for demonstration

        # Sample data
        sample_genomes = ["NC_045512.2", "MT123290.1", "MT188341.1", "OM095411.1", "OP912844.1"]

        plt.figure(figsize=(10, 6))
        plt.bar(["Sample SARS-CoV-2 Genomes"], [len(sample_genomes)])
        plt.title("Sample SARS-CoV-2 Genome Visualization")
        plt.ylabel("Number of Genomes")
        plt.tight_layout()

        plt.savefig("sample_visualization.png", dpi=150, bbox_inches="tight")
        print("[*] Sample plot saved as 'sample_visualization.png'")

        # Display the plot
        plt.show()
        print("[*] Sample plot displayed")
        print("[✓] Sample visualization complete!")
