#!/usr/bin/env python3
"""
src/advanced_visualizations.py
Advanced genomic visualizations for SARS-CoV-2 analysis.

Includes:
- GC Content Distribution
- Nucleotide Composition (Stacked Bar)
- Sequence Length Distribution
- Phylogenetic Tree
- Mutation Heatmap
- Temporal Evolution
"""

import random
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap
except ImportError as e:
    raise ImportError(
        f"Required packages not installed: {e}\nPlease install with: pip install matplotlib pandas numpy"
    ) from e

# Try to import optional dependencies
try:
    import seaborn as sns  # noqa: F401

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    print("[!] Seaborn not installed. Some visualizations will use matplotlib fallback.")

try:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("[!] Biopython not installed. Phylogenetic tree will use simplified visualization.")


# =============================================================================
# VARIANT DEFINITIONS
# =============================================================================

VARIANTS = {
    "Reference": {"color": "#1f77b4", "emergence": "2019-12"},
    "Alpha": {"color": "#ff7f0e", "emergence": "2020-09"},
    "Beta": {"color": "#2ca02c", "emergence": "2020-05"},
    "Gamma": {"color": "#d62728", "emergence": "2020-11"},
    "Delta": {"color": "#9467bd", "emergence": "2021-04"},
    "Omicron": {"color": "#8c564b", "emergence": "2021-11"},
    "Recent": {"color": "#e377c2", "emergence": "2022-06"},
}

# Spike protein mutation positions (simplified - key mutations per variant)
SPIKE_MUTATIONS = {
    "Reference": [],
    "Alpha": [69, 70, 144, 501, 570, 681, 716, 982, 1118],
    "Beta": [80, 215, 241, 242, 243, 417, 484, 501, 614, 701],
    "Gamma": [18, 20, 26, 138, 190, 417, 484, 501, 614, 655, 1027, 1176],
    "Delta": [19, 142, 156, 157, 158, 452, 478, 614, 681, 950],
    "Omicron": [
        67,
        69,
        70,
        95,
        142,
        143,
        144,
        211,
        212,
        339,
        371,
        373,
        375,
        417,
        440,
        446,
        477,
        478,
        484,
        493,
        496,
        498,
        501,
        505,
        547,
        614,
        655,
        679,
        681,
        764,
        796,
        856,
        954,
        969,
        981,
    ],
    "Recent": [
        67,
        69,
        70,
        95,
        142,
        143,
        144,
        211,
        339,
        371,
        373,
        375,
        417,
        440,
        446,
        452,
        478,
        484,
        486,
        490,
        493,
        498,
        501,
        505,
        614,
        655,
        679,
        681,
        764,
        796,
        954,
        969,
        981,
        1020,
    ],
}


# =============================================================================
# SAMPLE DATA GENERATION (for demonstration)
# =============================================================================


def generate_sample_sequence(length=29903, gc_bias=0.0):
    """Generate a sample genome sequence with optional GC bias."""
    # Base probabilities (roughly matching SARS-CoV-2: ~38% GC)
    base_gc = 0.38 + gc_bias
    base_at = 1 - base_gc

    weights = {
        "A": base_at * 0.5,
        "T": base_at * 0.5,
        "G": base_gc * 0.5,
        "C": base_gc * 0.5,
    }

    bases = list(weights.keys())
    probs = list(weights.values())

    return "".join(random.choices(bases, weights=probs, k=length))


def generate_sample_data(n_per_variant=5):
    """Generate sample genomic data for all variants."""
    data = []

    # GC content biases per variant (slight variations)
    gc_biases = {
        "Reference": 0.0,
        "Alpha": -0.01,
        "Beta": 0.005,
        "Gamma": -0.005,
        "Delta": 0.01,
        "Omicron": -0.02,
        "Recent": -0.015,
    }

    # Length variations (SARS-CoV-2 is ~29,903 bp)
    base_length = 29903

    for variant, info in VARIANTS.items():
        emergence = datetime.strptime(info["emergence"], "%Y-%m")

        for i in range(n_per_variant):
            # Add some variation to collection date
            days_offset = random.randint(0, 180)
            collection_date = emergence + timedelta(days=days_offset)

            # Generate sequence with variant-specific GC bias
            gc_bias = gc_biases.get(variant, 0.0) + random.uniform(-0.005, 0.005)
            length = base_length + random.randint(-50, 50)
            sequence = generate_sample_sequence(length, gc_bias)

            data.append(
                {
                    "accession": f"{variant[:2].upper()}{100000 + i}",
                    "variant": variant,
                    "sequence": sequence,
                    "length": length,
                    "collection_date": collection_date,
                    "gc_content": calculate_gc_content(sequence),
                }
            )

    return pd.DataFrame(data)


def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence) * 100


def calculate_nucleotide_composition(sequence):
    """Calculate nucleotide composition of a sequence."""
    if not sequence:
        return {"A": 0, "T": 0, "G": 0, "C": 0}

    total = len(sequence)
    return {
        "A": sequence.count("A") / total * 100,
        "T": sequence.count("T") / total * 100,
        "G": sequence.count("G") / total * 100,
        "C": sequence.count("C") / total * 100,
    }


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================


def plot_gc_content_distribution(df, output_path=None, show=True):
    """
    Plot GC content distribution across variants.

    Args:
        df: DataFrame with 'variant' and 'gc_content' columns
        output_path: Path to save the figure
        show: Whether to display the plot

    Returns:
        Figure object
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Box plot by variant
    ax1 = axes[0]
    variants = list(VARIANTS.keys())
    colors = [VARIANTS[v]["color"] for v in variants]

    data_by_variant = [df[df["variant"] == v]["gc_content"].values for v in variants]

    bp = ax1.boxplot(data_by_variant, labels=variants, patch_artist=True)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax1.set_xlabel("Variant", fontsize=12)
    ax1.set_ylabel("GC Content (%)", fontsize=12)
    ax1.set_title("GC Content Distribution by Variant", fontsize=14, fontweight="bold")
    ax1.tick_params(axis="x", rotation=45)
    ax1.grid(axis="y", alpha=0.3)

    # Right: Histogram overlay
    ax2 = axes[1]
    for variant in variants:
        variant_data = df[df["variant"] == variant]["gc_content"]
        ax2.hist(variant_data, bins=15, alpha=0.5, label=variant, color=VARIANTS[variant]["color"], density=True)

    ax2.set_xlabel("GC Content (%)", fontsize=12)
    ax2.set_ylabel("Density", fontsize=12)
    ax2.set_title("GC Content Distribution (Overlaid)", fontsize=14, fontweight="bold")
    ax2.legend(loc="upper right")
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"[*] GC content plot saved to: {output_path}")

    if show:
        plt.show()

    return fig


def plot_nucleotide_composition(df, output_path=None, show=True):
    """
    Plot nucleotide composition as stacked bar chart.

    Args:
        df: DataFrame with 'variant' and 'sequence' columns
        output_path: Path to save the figure
        show: Whether to display the plot

    Returns:
        Figure object
    """
    # Calculate average composition per variant
    composition_data = []

    for variant in VARIANTS.keys():
        variant_seqs = df[df["variant"] == variant]["sequence"]

        avg_comp = {"A": 0, "T": 0, "G": 0, "C": 0}
        for seq in variant_seqs:
            comp = calculate_nucleotide_composition(seq)
            for base in avg_comp:
                avg_comp[base] += comp[base]

        n = len(variant_seqs)
        if n > 0:
            for base in avg_comp:
                avg_comp[base] /= n

        avg_comp["variant"] = variant
        composition_data.append(avg_comp)

    comp_df = pd.DataFrame(composition_data)

    fig, ax = plt.subplots(figsize=(12, 7))

    x = np.arange(len(comp_df))
    width = 0.6

    # Stack the bars
    base_colors = {"A": "#e41a1c", "T": "#377eb8", "G": "#4daf4a", "C": "#ff7f00"}

    bottom = np.zeros(len(comp_df))
    for base in ["A", "T", "G", "C"]:
        values = comp_df[base].values
        ax.bar(x, values, width, bottom=bottom, label=base, color=base_colors[base], alpha=0.85)

        # Add percentage labels in the middle of each segment
        for i, (val, bot) in enumerate(zip(values, bottom)):
            if val > 3:  # Only label if segment is large enough
                ax.text(
                    i,
                    bot + val / 2,
                    f"{val:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                    color="white",
                )

        bottom += values

    ax.set_xlabel("Variant", fontsize=12)
    ax.set_ylabel("Nucleotide Composition (%)", fontsize=12)
    ax.set_title("Nucleotide Composition by Variant", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(comp_df["variant"], rotation=45, ha="right")
    ax.legend(title="Nucleotide", loc="upper right")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"[*] Nucleotide composition plot saved to: {output_path}")

    if show:
        plt.show()

    return fig


def plot_sequence_length_distribution(df, output_path=None, show=True):
    """
    Plot sequence length distribution.

    Args:
        df: DataFrame with 'variant' and 'length' columns
        output_path: Path to save the figure
        show: Whether to display the plot

    Returns:
        Figure object
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Violin plot
    ax1 = axes[0]
    variants = list(VARIANTS.keys())

    data_by_variant = [df[df["variant"] == v]["length"].values for v in variants]
    colors = [VARIANTS[v]["color"] for v in variants]

    parts = ax1.violinplot(data_by_variant, positions=range(len(variants)), showmeans=True, showmedians=True)

    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)

    ax1.set_xticks(range(len(variants)))
    ax1.set_xticklabels(variants, rotation=45, ha="right")
    ax1.set_xlabel("Variant", fontsize=12)
    ax1.set_ylabel("Sequence Length (bp)", fontsize=12)
    ax1.set_title("Sequence Length Distribution (Violin Plot)", fontsize=14, fontweight="bold")
    ax1.axhline(y=29903, color="red", linestyle="--", alpha=0.5, label="Reference (29,903 bp)")
    ax1.legend()
    ax1.grid(axis="y", alpha=0.3)

    # Right: Histogram
    ax2 = axes[1]
    ax2.hist(df["length"], bins=30, color="steelblue", alpha=0.7, edgecolor="black")
    ax2.axvline(x=29903, color="red", linestyle="--", linewidth=2, label="Reference (29,903 bp)")
    ax2.axvline(
        x=df["length"].mean(), color="green", linestyle="-", linewidth=2, label=f"Mean ({df['length'].mean():,.0f} bp)"
    )
    ax2.set_xlabel("Sequence Length (bp)", fontsize=12)
    ax2.set_ylabel("Count", fontsize=12)
    ax2.set_title("Overall Sequence Length Distribution", fontsize=14, fontweight="bold")
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"[*] Sequence length plot saved to: {output_path}")

    if show:
        plt.show()

    return fig


def plot_phylogenetic_tree(df, output_path=None, show=True):
    """
    Plot a phylogenetic tree showing variant relationships.
    Uses Biopython if available, otherwise creates a simplified visualization.

    Args:
        df: DataFrame with variant information
        output_path: Path to save the figure
        show: Whether to display the plot

    Returns:
        Figure object
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    if HAS_BIOPYTHON:
        # Create a distance matrix based on mutation differences
        variants = list(VARIANTS.keys())

        # Calculate pairwise distances based on spike mutations
        matrix = []
        for i, v1 in enumerate(variants):
            row = []
            for j, v2 in enumerate(variants[: i + 1]):
                if i == j:
                    row.append(0)
                else:
                    # Distance = number of unique mutations between variants
                    m1 = set(SPIKE_MUTATIONS[v1])
                    m2 = set(SPIKE_MUTATIONS[v2])
                    distance = len(m1.symmetric_difference(m2)) / 100
                    row.append(distance)
            matrix.append(row)

        # Create distance matrix
        dm = DistanceMatrix(variants, matrix)

        # Build tree using UPGMA
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)

        # Draw tree
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title(
            "Phylogenetic Tree of SARS-CoV-2 Variants\n(Based on Spike Protein Mutations)",
            fontsize=14,
            fontweight="bold",
        )
    else:
        # Simplified tree visualization without Biopython
        # Create a manual tree structure based on known evolutionary relationships

        # Define tree structure (simplified)
        tree_structure = {
            "Reference": {"x": 0, "y": 3, "parent": None},
            "Alpha": {"x": 1, "y": 5, "parent": "Reference"},
            "Beta": {"x": 1, "y": 4, "parent": "Reference"},
            "Gamma": {"x": 1, "y": 3, "parent": "Reference"},
            "Delta": {"x": 1, "y": 2, "parent": "Reference"},
            "Omicron": {"x": 2, "y": 1, "parent": "Reference"},
            "Recent": {"x": 3, "y": 0.5, "parent": "Omicron"},
        }

        # Draw connections
        for _variant, info in tree_structure.items():
            if info["parent"]:
                parent_info = tree_structure[info["parent"]]
                # Draw horizontal then vertical line
                ax.plot(
                    [parent_info["x"], info["x"] - 0.3, info["x"]],
                    [parent_info["y"], parent_info["y"], info["y"]],
                    color="gray",
                    linewidth=2,
                    zorder=1,
                )

        # Draw nodes
        for variant, info in tree_structure.items():
            color = VARIANTS[variant]["color"]
            ax.scatter(info["x"], info["y"], s=500, c=color, zorder=2, edgecolors="black", linewidth=2)
            ax.annotate(
                variant, (info["x"], info["y"]), fontsize=11, fontweight="bold", ha="center", va="center", color="white"
            )

            # Add emergence date
            emergence = VARIANTS[variant]["emergence"]
            ax.annotate(emergence, (info["x"] + 0.15, info["y"] - 0.2), fontsize=9, color="gray")

        ax.set_xlim(-0.5, 4)
        ax.set_ylim(-0.5, 6)
        ax.set_title(
            "Phylogenetic Tree of SARS-CoV-2 Variants\n(Simplified Evolutionary Relationships)",
            fontsize=14,
            fontweight="bold",
        )
        ax.axis("off")

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"[*] Phylogenetic tree saved to: {output_path}")

    if show:
        plt.show()

    return fig


def plot_mutation_heatmap(output_path=None, show=True):
    """
    Plot mutation heatmap showing spike protein mutations across variants.

    Args:
        output_path: Path to save the figure
        show: Whether to display the plot

    Returns:
        Figure object
    """
    # Create mutation matrix
    variants = list(VARIANTS.keys())

    # Get all unique mutation positions
    all_positions = set()
    for mutations in SPIKE_MUTATIONS.values():
        all_positions.update(mutations)
    all_positions = sorted(all_positions)

    # Create binary matrix
    matrix = np.zeros((len(variants), len(all_positions)))
    for i, variant in enumerate(variants):
        for j, pos in enumerate(all_positions):
            if pos in SPIKE_MUTATIONS[variant]:
                matrix[i, j] = 1

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 8))

    # Create custom colormap
    colors = ["#f7f7f7", "#d73027"]
    cmap = LinearSegmentedColormap.from_list("mutation", colors)

    # Plot heatmap
    im = ax.imshow(matrix, cmap=cmap, aspect="auto", interpolation="nearest")

    # Set labels
    ax.set_yticks(range(len(variants)))
    ax.set_yticklabels(variants, fontsize=11)

    # Show fewer x-tick labels to avoid crowding
    step = max(1, len(all_positions) // 20)
    ax.set_xticks(range(0, len(all_positions), step))
    ax.set_xticklabels([all_positions[i] for i in range(0, len(all_positions), step)], rotation=90, fontsize=9)

    ax.set_xlabel("Spike Protein Position", fontsize=12)
    ax.set_ylabel("Variant", fontsize=12)
    ax.set_title(
        "Spike Protein Mutation Heatmap\n(Key mutations defining each variant)", fontsize=14, fontweight="bold"
    )

    # Add mutation count annotations on the right
    for i, variant in enumerate(variants):
        count = len(SPIKE_MUTATIONS[variant])
        ax.annotate(
            f"{count} mutations",
            xy=(len(all_positions) + 0.5, i),
            fontsize=10,
            va="center",
            color=VARIANTS[variant]["color"],
            fontweight="bold",
        )

    # Adjust layout to make room for annotations
    plt.subplots_adjust(right=0.85)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.5, label="Mutation Present")
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["No", "Yes"])

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"[*] Mutation heatmap saved to: {output_path}")

    if show:
        plt.show()

    return fig


def plot_temporal_evolution(df, output_path=None, show=True):
    """
    Plot temporal evolution of variants.

    Args:
        df: DataFrame with 'variant' and 'collection_date' columns
        output_path: Path to save the figure
        show: Whether to display the plot

    Returns:
        Figure object
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Top: Timeline showing variant emergence
    ax1 = axes[0]

    variants = list(VARIANTS.keys())
    y_positions = list(range(len(variants)))

    for i, variant in enumerate(variants):
        emergence = datetime.strptime(VARIANTS[variant]["emergence"], "%Y-%m")
        color = VARIANTS[variant]["color"]

        # Draw emergence point
        ax1.scatter(emergence, i, s=300, c=color, zorder=3, edgecolors="black", linewidth=2)

        # Draw timeline extending to "now"
        end_date = datetime(2024, 1, 1)
        ax1.plot([emergence, end_date], [i, i], color=color, linewidth=4, alpha=0.5, zorder=2)

        # Add label
        ax1.annotate(variant, (emergence, i + 0.3), fontsize=11, fontweight="bold", color=color)

    ax1.set_yticks(y_positions)
    ax1.set_yticklabels([""] * len(variants))
    ax1.set_xlabel("Date", fontsize=12)
    ax1.set_title("SARS-CoV-2 Variant Emergence Timeline", fontsize=14, fontweight="bold")
    ax1.grid(axis="x", alpha=0.3)
    ax1.set_xlim(datetime(2019, 10, 1), datetime(2024, 3, 1))

    # Bottom: Stacked area chart of sample counts over time
    ax2 = axes[1]

    # Group by month and variant
    df["month"] = df["collection_date"].dt.to_period("M")
    monthly_counts = df.groupby(["month", "variant"]).size().unstack(fill_value=0)

    # Ensure all variants are present
    for variant in variants:
        if variant not in monthly_counts.columns:
            monthly_counts[variant] = 0

    monthly_counts = monthly_counts[variants]  # Reorder columns

    # Convert period to datetime for plotting
    x = monthly_counts.index.to_timestamp()

    # Create stacked area chart
    colors = [VARIANTS[v]["color"] for v in variants]
    ax2.stackplot(x, monthly_counts.T, labels=variants, colors=colors, alpha=0.8)

    ax2.set_xlabel("Date", fontsize=12)
    ax2.set_ylabel("Sample Count", fontsize=12)
    ax2.set_title("Temporal Distribution of Samples by Variant", fontsize=14, fontweight="bold")
    ax2.legend(loc="upper left", ncol=4)
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"[*] Temporal evolution plot saved to: {output_path}")

    if show:
        plt.show()

    return fig


# =============================================================================
# MAIN FUNCTION
# =============================================================================


def generate_all_visualizations(output_dir, n_samples_per_variant=10, show=False):
    """
    Generate all visualizations with sample data.

    Args:
        output_dir: Directory to save figures
        n_samples_per_variant: Number of sample sequences per variant
        show: Whether to display plots

    Returns:
        List of generated file paths
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"[*] Generating sample data ({n_samples_per_variant} samples per variant)...")
    df = generate_sample_data(n_samples_per_variant)
    print(f"[*] Generated {len(df)} total samples")

    generated_files = []

    # 1. GC Content Distribution
    print("\n[1/6] Generating GC Content Distribution...")
    gc_path = output_path / "gc_content_distribution.png"
    plot_gc_content_distribution(df, output_path=gc_path, show=show)
    generated_files.append(gc_path)

    # 2. Nucleotide Composition
    print("\n[2/6] Generating Nucleotide Composition...")
    nuc_path = output_path / "nucleotide_composition.png"
    plot_nucleotide_composition(df, output_path=nuc_path, show=show)
    generated_files.append(nuc_path)

    # 3. Sequence Length Distribution
    print("\n[3/6] Generating Sequence Length Distribution...")
    len_path = output_path / "sequence_length_distribution.png"
    plot_sequence_length_distribution(df, output_path=len_path, show=show)
    generated_files.append(len_path)

    # 4. Phylogenetic Tree
    print("\n[4/6] Generating Phylogenetic Tree...")
    tree_path = output_path / "phylogenetic_tree.png"
    plot_phylogenetic_tree(df, output_path=tree_path, show=show)
    generated_files.append(tree_path)

    # 5. Mutation Heatmap
    print("\n[5/6] Generating Mutation Heatmap...")
    heatmap_path = output_path / "mutation_heatmap.png"
    plot_mutation_heatmap(output_path=heatmap_path, show=show)
    generated_files.append(heatmap_path)

    # 6. Temporal Evolution
    print("\n[6/6] Generating Temporal Evolution...")
    temporal_path = output_path / "temporal_evolution.png"
    plot_temporal_evolution(df, output_path=temporal_path, show=show)
    generated_files.append(temporal_path)

    print(f"\n[*] All visualizations saved to: {output_path}")
    return generated_files


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate advanced genomic visualizations")
    parser.add_argument("--output", "-o", default="./visualizations", help="Output directory for figures")
    parser.add_argument("--samples", "-n", type=int, default=10, help="Number of samples per variant")
    parser.add_argument("--show", action="store_true", help="Display plots interactively")

    args = parser.parse_args()

    generate_all_visualizations(args.output, args.samples, args.show)
