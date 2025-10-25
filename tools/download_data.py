#!/usr/bin/env python3
"""
tools/download_data.py
Download a small subset of SARS-CoV-2 genomes from NCBI (~50MB) to Google Drive.

Usage:
    python tools/download_data.py [--output /path/to/google/drive/ncbi_data]
"""

import argparse
import os
import subprocess
from pathlib import Path

def run_cmd(cmd):
    """Run a shell command and raise on failure."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout

def download_small_subset(output_dir):
    """Download a small subset of SARS-CoV-2 genomes (~50MB)."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Pre-selected small set of SARS-CoV-2 accessions for ~50MB download
    # These are well-known complete genomes that should total around 50MB
    accessions = [
        "NC_045512.2",  # Reference genome (Wuhan-Hu-1)
        "MT123290.1",   # Early variant
        "MT188341.1",   # Alpha variant example
        "MW633477.1",   # Delta variant example
        "OM095411.1",   # Omicron variant example
        "ON563414.1",   # Recent variant
        "OP912844.1",   # Recent variant
        "OR064389.1",   # Recent variant
    ]
    
    print(f"[*] Downloading {len(accessions)} SARS-CoV-2 genomes to {output_path}")
    print(f"[*] Expected total size: ~50MB")
    
    for i, accession in enumerate(accessions, 1):
        print(f"[{i}/{len(accessions)}] Downloading {accession}...")
        cmd = (
            f"datasets download virus genome accession {accession} "
            f"--include genome,annotation_report "
            f"--host human --assembly-level complete "
            f"--filename {output_path}/{accession}.zip --no-progressbar"
        )
        try:
            run_cmd(cmd)
            print(f"    ✓ Downloaded {accession}")
        except RuntimeError as e:
            print(f"    ✗ Failed to download {accession}: {e}")
    
    print(f"[✓] Download complete. Files saved to: {output_path}")
    return str(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download small SARS-CoV-2 dataset from NCBI")
    parser.add_argument("--output", type=str, 
                       default="G:/My Drive/ncbi_data",  # Default to Google Drive
                       help="Output directory (default: G:/My Drive/ncbi_data)")
    args = parser.parse_args()
    
    # Check if output directory exists, if not try common Google Drive paths
    if not os.path.exists(args.output):
        common_paths = [
            "G:/My Drive/ncbi_data",
            "C:/Users/zhouf/Google Drive/ncbi_data", 
            "D:/Google Drive/ncbi_data",
            "./ncbi_data"  # Fallback to local directory
        ]
        
        for path in common_paths:
            if os.path.exists(os.path.dirname(path)) or path == "./ncbi_data":
                args.output = path
                print(f"[*] Using output directory: {args.output}")
                break
        else:
            print(f"[!] Warning: Could not find Google Drive. Using local directory: ./ncbi_data")
            args.output = "./ncbi_data"
    
    download_small_subset(args.output)
