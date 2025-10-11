from src.download import download_dataset_balanced
from src.visualization import visualize_dataset

def main():
    dataset_path = download_dataset_balanced(
        virus_name="sars-cov-2",
        output_dir="data",
        size_gb=8,
        seed=42,
        workers=4
    )
    visualize_dataset(dataset_path)

if __name__ == "__main__":
    main()
