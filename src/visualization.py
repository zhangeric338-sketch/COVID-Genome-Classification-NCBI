import json
import os
import matplotlib.pyplot as plt
import pandas as pd

def load_metadata(jsonl_path):
    """Load NCBI genome metadata from JSON Lines file into a DataFrame."""
    records = []
    with open(jsonl_path, "r") as f:
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
    counts.plot(kind="bar", figsize=(12,6))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

def visualize_dataset(dataset_path):
    """
    Visualize the metadata in an extracted NCBI virus dataset folder.
    
    Args:
        dataset_path (str): Path to the extracted dataset folder returned by download_dataset().
    """
    # Path to metadata file
    metadata_file = os.path.join(dataset_path, "ncbi_dataset", "metadata", "genome_metadata.jsonl")
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

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
