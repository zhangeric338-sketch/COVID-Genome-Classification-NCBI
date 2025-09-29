#!/usr/bin/env python3
"""
tools/download_data.py
Subset and download SARS-CoV-2 assemblies (complete, human host) from NCBI Datasets.

- Starts with metadata query
- Deterministically samples assemblies to limit download size (~5-10 GB)
- Supports parallel download
- Stores data in Google Drive or a local folder

Requirements:
    pip install pandas tqdm joblib
    Install NCBI Datasets CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

Usage:
    python tools/download_data.py --workers 4 --seed 42 --size-gb 8 --output /content/drive/MyDrive/sarscov2_data
"""
import argparse
import os
import random
import subprocess
import tempfile
from pathlib import Path
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm


def run_cmd(cmd):
    """Run a shell command and raise if it fails."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout


def fetch_metadata(tmpdir):
    """Fetch SARS-CoV-2 assembly metadata from NCBI Datasets."""
    metadata_path = Path(tmpdir) / "sarscov2_metadata.jsonl"
    cmd = (
        "datasets summary virus genome taxon sars-cov-2 "
        "--assembly-level complete "
        "--host 'Homo sapiens' "
        f"--as-json-lines > {metadata_path}"
    )
    run_cmd(cmd)
    return metadata_path


def parse_metadata(metadata_path):
    """Convert JSONL metadata to DataFrame."""
    df = pd.read_json(metadata_path, lines=True)
    # Flatten for accession IDs and estimated size
    records = []
    for _, row in df.iterrows():
        accession = row["accession"]
        size_bytes = row.get("seq_length", 0) * 2  # fasta + annotation rough estimate
        records.append({"accession": accession, "size_bytes": size_bytes})
    return pd.DataFrame(records)


def sample_accessions(df, target_size_gb, seed):
    """Sample accessions deterministically until ~target_size_gb is reached."""
    rng = random.Random(seed)
    accs = df.sample(frac=1, random_state=seed).to_dict("records")
    total, chosen = 0, []
    target_bytes = target_size_gb * (1024**3)
    for rec in accs:
        if total + rec["size_bytes"] > target_bytes:
            break
        chosen.append(rec["accession"])
        total += rec["size_bytes"]
    return chosen


def download_one(acc, outdir):
    """Download one accession with datasets CLI."""
    outdir = Path(outdir)
    cmd = (
        f"datasets download virus genome accession {acc} "
        f"--include genome,protein,cds,annotation_report "
        f"--filename {outdir}/{acc}.zip --no-progressbar"
    )
    run_cmd(cmd)
    return acc


def download_all(accessions, outdir, workers):
    """Download all accessions in parallel."""
    os.makedirs(outdir, exist_ok=True)
    Parallel(n_jobs=workers)(
        delayed(download_one)(acc, outdir) for acc in tqdm(accessions, desc="Downloading")
    )


def main(args):
    with tempfile.TemporaryDirectory() as tmpdir:
        print("[*] Fetching metadata...")
        metadata_path = fetch_metadata(tmpdir)
        df = parse_metadata(metadata_path)
        print(f"[*] {len(df)} complete assemblies found")

        chosen = sample_accessions(df, args.size_gb, args.seed)
        print(f"[*] Selected {len(chosen)} accessions (~{args.size_gb} GB target)")

        print(f"[*] Downloading data to {args.output}")
        download_all(chosen, args.output, args.workers)

        print("[✓] Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SARS-CoV-2 assemblies from NCBI")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling")
    parser.add_argument("--size-gb", type=int, default=8, help="Approximate dataset size to download")
    parser.add_argument("--output", type=str, required=True, help="Output directory (Google Drive mount)")
    args = parser.parse_args()
    main(args)

