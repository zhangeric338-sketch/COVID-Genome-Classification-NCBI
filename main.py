from src.download import download_dataset


def main():
    """
    Download viral genome datasets from NCBI.
    """

    dataset_path = download_dataset()


if __name__ == "__main__":
    main()
