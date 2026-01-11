import os
import subprocess
import zipfile


def download_dataset(virus_name: str = "sars-cov-2", output_dir: str = "data") -> str:
    """
    Download a virus genome dataset from NCBI using the datasets CLI.
    Only downloads if the dataset hasn't been downloaded before.

    Args:
        virus_name (str): Taxon name or ID (e.g., "sars-cov-2", "monkeypox").
        output_dir (str): Directory where files will be stored.

    Returns:
        str: Path to the extracted dataset folder.
    """
    os.makedirs(output_dir, exist_ok=True)
    zip_path = os.path.join(output_dir, f"{virus_name}.zip")
    extract_path = os.path.join(output_dir, virus_name)

    # If extracted folder already exists, skip download
    if os.path.exists(extract_path):
        print(f"✅ Dataset already exists at {extract_path}, skipping download.")
        return extract_path

    # Command: datasets download virus genome taxon <virus> --filename <zip>
    cmd = ["datasets", "download", "virus", "genome", "taxon", virus_name, "--filename", zip_path]

    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Extract
    print(f"Extracting {zip_path} to {extract_path}")
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(extract_path)

    print(f"✅ Dataset downloaded and extracted to {extract_path}")
    return extract_path
