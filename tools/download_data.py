#!/usr/bin/env python3
"""
tools/download_data.py
Download a balanced subset of SARS-CoV-2 (complete, human) assemblies from NCBI,
such that total download ≈ size_gb, and the data is balanced across variants.

Usage:
    python tools/download_data.py --workers 4 --seed 42 --size-gb 8 --output /path
"""

import argparse
import os
import subprocess
import tempfile
import random
import json
from pathlib import Path
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm

def run_cmd(cmd):
    """Run a shell command and raise on failure."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout

def fetch_metadata(tmpdir):
    """Fetch SARS-CoV-2 metadata and save as JSONL."""
    metadata_path = Path(tmpdir) / "sarscov2_metadata.jsonl"
    cmd = f"datasets summary virus genome taxon sars-cov-2 --as-json-lines > {metadata_path}"
    run_cmd(cmd)
    return metadata_path

def parse_metadata(metadata_path):
    """Parse JSONL metadata and return DataFrame with accession, lineage, seq_length."""
    records = []
    with open(metadata_path, "r") as f:
        for line in f:
            obj = json.loads(line)

            # Keep only complete genomes
            if obj.get("completeness") != "COMPLETE":
                continue

            # Keep only human host
            host_name = obj.get("host", {}).get("organism_name")
            if host_name != "Homo sapiens":
                continue

            accession = obj.get("accession")
            seq_len = obj.get("length", 0)
            lineage = obj.get("virus", {}).get("pangolin_classification", "Unknown")

            records.append({
                "accession": accession,
                "lineage": lineage,
                "seq_length": seq_len
            })

    df = pd.DataFrame(records)
    df = df.dropna(subset=["accession"])
    df = df[df["seq_length"] > 0]
    return df

def sample_balanced_by_variant(df, total_gb, seed):
    """Sample accessions so total ~total_gb, distributing bytes across variants."""
    random.seed(seed)
    target_bytes = total_gb * (1024 ** 3)
    variants = df["lineage"].unique().tolist()
    n_variants = len(variants)
    if n_variants == 0:
        raise ValueError("No variants found in metadata")

    bytes_per_variant = target_bytes / n_variants
    chosen = []

    for v in variants:
        df_v = df[df["lineage"] == v].sample(frac=1, random_state=seed)
        cumulative = 0
        for _, row in df_v.iterrows():
            size_est = row["seq_length"] * 2  # genome + annotation rough estimate
            if cumulative + size_est > bytes_per_variant:
                break
            chosen.append(row["accession"])
            cumulative += size_est

    random.shuffle(chosen)
    return chosen

def download_one(accession, outdir):
    """Download one accession via datasets CLI."""
    outdir = Path(outdir)
    os.makedirs(outdir, exist_ok=True)
    cmd = (
        f"datasets download virus genome accession {accession} "
        f"--include genome,protein,cds,annotation_report "
        f"--filename {outdir}/{accession}.zip --no-progressbar"
    )
    run_cmd(cmd)
    return accession

def download_all(accessions, outdir, workers):
    """Download all accessions in parallel."""
    os.makedirs(outdir, exist_ok=True)
    Parallel(n_jobs=workers)(
        delayed(download_one)(acc, outdir) for acc in tqdm(accessions, desc="Downloading")
    )

def download_dataset_balanced(virus_name="sars-cov-2", output_dir="data", size_gb=8, seed=42, workers=4):
    """Top-level function: fetch metadata, sample balanced, download subset."""
    os.makedirs(output_dir, exist_ok=True)
    genomes_dir = Path(output_dir) / virus_name

    if genomes_dir.exists():
        print(f"✅ Genomes folder already exists at {genomes_dir}, skipping download.")
        return str(genomes_dir)

    with tempfile.TemporaryDirectory() as tmp:
        print("[*] Fetching metadata...")
        mpath = fetch_metadata(tmp)
        print("[*] Parsing metadata...")
        df = parse_metadata(mpath)
        print(f"[*] Found {len(df)} complete human SARS-CoV-2 genomes across {df['lineage'].nunique()} variants.")
        print("[*] Sampling balanced subset...")
        chosen = sample_balanced_by_variant(df, size_gb, seed)
        print(f"[*] Selected {len(chosen)} genomes (~{size_gb} GB target).")
        print("[*] Downloading selected genomes...")
        download_all(chosen, genomes_dir, workers=workers)
        print("[✓] Download complete.")

    return str(genomes_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SARS-CoV-2 assemblies from NCBI")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling")
    parser.add_argument("--size-gb", type=int, default=8, help="Approximate dataset size to download")
    parser.add_argument("--output", type=str, required=True, help="Output directory (Google Drive mount)")
    args = parser.parse_args()
    download_dataset_balanced(
        virus_name="sars-cov-2",
        output_dir=args.output,
        size_gb=args.size_gb,
        seed=args.seed,
        workers=args.workers
    )
