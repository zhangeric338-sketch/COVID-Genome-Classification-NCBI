from tools.download_data import download_dataset_balanced
from src.visualization import visualize_dataset, split_dataset, visualize_train_test_split, visualize_dataset_composition

def main():
    # Download dataset
    dataset_path = download_dataset_balanced(
        virus_name="sars-cov-2",
        output_dir="data",
        size_gb=0.05,
        seed=42,
        workers=4
    )
    
    # Visualize downloaded data
    visualize_dataset(dataset_path)
    
    # Split into train/test datasets
    print("\n[*] Splitting dataset into train/test...")
    train_files, test_files = split_dataset("data", test_ratio=0.2, random_seed=42)
    
    # Visualize train/test split
    print("\n[*] Creating train/test split visualization...")
    visualize_train_test_split("data")
    
    # Visualize dataset composition
    print("\n[*] Creating detailed dataset composition visualization...")
    visualize_dataset_composition("data")

if __name__ == "__main__":
    main()
