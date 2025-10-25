from tools.download_data import download_small_subset
from src.visualization import visualize_dataset

def main():
    dataset_path = download_small_subset("data")
    visualize_dataset(dataset_path)

if __name__ == "__main__":
    main()
