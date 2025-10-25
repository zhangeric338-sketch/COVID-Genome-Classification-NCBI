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
    plt.savefig(filename, dpi=150, bbox_inches='tight')
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
        zip_files = [f for f in os.listdir(dataset_path) if f.endswith('.zip')]
        print(f"[*] Found {len(zip_files)} downloaded genome files:")
        for i, zip_file in enumerate(zip_files, 1):
            print(f"  {i}. {zip_file}")
        
        # Create a simple visualization showing download summary
        import matplotlib.pyplot as plt
        
        # Ensure plots display in Google Colab
        plt.figure(figsize=(10, 6))
        plt.bar(['Downloaded Genomes'], [len(zip_files)])
        plt.title(f'SARS-CoV-2 Genome Downloads ({len(zip_files)} genomes)')
        plt.ylabel('Number of Genomes')
        plt.tight_layout()
        
        # Display the plot in Google Colab
        plt.show()
        
        # Also save the plot as an image file
        plt.savefig('download_summary.png', dpi=150, bbox_inches='tight')
        print(f"[*] Plot saved as 'download_summary.png'")
        
        print(f"[✓] Visualization complete - {len(zip_files)} genomes downloaded successfully")
