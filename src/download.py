import subprocess
import os
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
    cmd = [
        "datasets", "download", "virus", "genome",
        "taxon", virus_name,
        "--filename", zip_path
    ]

    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Extract
    print(f"Extracting {zip_path} to {extract_path}")
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(extract_path)

    print(f"✅ Dataset downloaded and extracted to {extract_path}")
    return extract_path

def download_dataset_balanced(
    virus_name="sars-cov-2", output_dir="data", size_gb=8, seed=42, workers=4
):
    """
    Top-level function: fetch metadata, sample balanced by variant, download subset.
    Returns path to downloaded genomes folder.
    """
    os.makedirs(output_dir, exist_ok=True)
    genomes_dir = Path(output_dir) / virus_name

    # If already exists, skip
    if genomes_dir.exists():
        print(f"✅ Genomes folder already exists at {genomes_dir}, skipping download.")
        return str(genomes_dir)

    with tempfile.TemporaryDirectory() as tmp:
        print("[*] Fetching metadata...")
        mpath = fetch_metadata(tmp)
        print("[*] Parsing metadata...")
        df = parse_metadata(mpath)
        print(f"[*] Found {len(df)} entries, across {df['lineage'].nunique()} variants.")
        print("[*] Sampling balanced subset...")
        chosen = sample_balanced_by_variant(df, size_gb, seed)
        print(f"[*] Selected {len(chosen)} genomes in subset (~{size_gb} GB target).")
        print("[*] Downloading selected genomes...")
        download_all(chosen, genomes_dir, workers=workers)
        print("[✓] Download complete.")
    return str(genomes_dir)
