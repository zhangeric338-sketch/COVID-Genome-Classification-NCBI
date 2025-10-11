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
    # We use shell=True sometimes for redirections, etc.
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout

def fetch_metadata(tmpdir):
    """
    Fetch metadata for SARS-CoV-2 (complete, human) using NCBI Datasets.
    Saves to JSONL in tmpdir and returns path.
    """
    metadata_path = Path(tmpdir) / "sarscov2_metadata.jsonl"
    # Use datasets summary or download command to get metadata lines, with filtering flags
    # Note: using `datasets summary ... --as-json-lines` to just pull metadata
    cmd = (
        "datasets summary virus genome taxon sars-cov-2 "
        "--assembly-level complete "
        "--host 'Homo sapiens' "
        "--as-json-lines > " + str(metadata_path)
    )
    run_cmd(cmd)
    return metadata_path

def parse_metadata(metadata_path):
    """
    Parse metadata JSONL and return a DataFrame with columns:
    accession, lineage (or variant), seq_length
    """
    records = []
    with open(metadata_path, "r") as f:
        for line in f:
            obj = json.loads(line)
            # Some metadata fields; adapt based on what the JSON has
            accession = obj.get("accession")
            seq_len = obj.get("seq_length", 0)
            # Try different fields for lineage / variant
            lineage = None
            # Some metadata have `virus_pangolin` or `lineage` or similar
            if "virus_pangolin" in obj:
                lineage = obj["virus_pangolin"]
            elif "lineage" in obj:
                lineage = obj["lineage"]
            else:
                # Maybe in nested metadata
                lineage = obj.get("isolate_lineage") or obj.get("virus", {}).get("lineage")
            if lineage is None:
                lineage = "Unknown"
            records.append({
                "accession": accession,
                "lineage": lineage,
                "seq_length": seq_len
            })
    df = pd.DataFrame(records)
    # Drop entries missing accession or seq_length = 0
    df = df.dropna(subset=["accession"])
    df = df[df["seq_length"] > 0]
    return df

def sample_balanced_by_variant(df, total_gb, seed):
    """
    From df, sample accessions so that total bytes ~ total_gb,
    distributing budget equally (or nearly so) across variants.
    Returns list of accessions.
    """
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
            size_est = row["seq_length"] * 2  # rough: genome + annotation
            if cumulative + size_est > bytes_per_variant:
                break
            chosen.append(row["accession"])
            cumulative += size_est
    # Optionally shuffle the final list
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
    """Download the list of accessions in parallel."""
    os.makedirs(outdir, exist_ok=True)
    Parallel(n_jobs=workers)(
        delayed(download_one)(acc, outdir) for acc in tqdm(accessions, desc="Downloading")
    )

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
