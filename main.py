from src.download import download_dataset
from src.visualization import visualize_dataset

def main():
    """
    Download viral genome datasets from NCBI and visualize the metadata.
    """
    # Download dataset (you can pass virus_name/output_dir if your download_dataset supports it)
    dataset_path = download_dataset(virus_name="sars-cov-2", output_dir="data")
    
    # Visualize the downloaded dataset
    visualize_dataset(dataset_path)

if __name__ == "__main__":
    main()
