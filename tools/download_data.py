#!/usr/bin/env python3
"""
tools/download_data.py
Download a small subset of SARS-CoV-2 genomes from NCBI (~50MB) to Google Drive.

Usage:
    python tools/download_data.py [--output /path/to/google/drive/ncbi_data]
"""

import argparse
import concurrent.futures
import glob
import json
import os
import random
import shutil
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path


def check_cli_available():
    """Check if NCBI datasets CLI is available."""
    return shutil.which("datasets") is not None


def query_ncbi_genome_count(virus_name="sars-cov-2", host="human", complete_only=True):
    """
    Query NCBI for total count and estimated size of available genomes.

    Args:
        virus_name (str): Virus name to query (default: "sars-cov-2")
        host (str): Host filter (default: "human")
        complete_only (bool): Whether to filter for complete genomes only (default: True)

    Returns:
        dict: Dictionary with total_count, estimated_size_mb, estimated_size_gb, estimated_size_tb
              Returns None if query fails.
    """
    try:
        # Use NCBI Datasets API to get summary
        api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome"

        body = {"taxon": virus_name, "host": host}
        if complete_only:
            body["complete_only"] = True

        payload_json = json.dumps(body).encode("utf-8")
        req = urllib.request.Request(
            api_url, data=payload_json, headers={"Content-Type": "application/json", "Accept": "application/json"}
        )

        print(f"[*] Querying NCBI for {virus_name} genome count...")

        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.loads(response.read().decode("utf-8"))

            # Extract count from response
            total_count = data.get("record_count", 0)

            # Estimate size based on ~6MB per genome
            estimated_size_mb = total_count * 6
            estimated_size_gb = estimated_size_mb / 1024
            estimated_size_tb = estimated_size_gb / 1024

            return {
                "total_count": total_count,
                "estimated_size_mb": estimated_size_mb,
                "estimated_size_gb": estimated_size_gb,
                "estimated_size_tb": estimated_size_tb,
            }

    except urllib.error.HTTPError as e:
        print(f"[!] HTTP error querying NCBI: {e.code} - {e.reason}")
        return None
    except urllib.error.URLError as e:
        print(f"[!] URL error querying NCBI: {e!s}")
        return None
    except json.JSONDecodeError as e:
        print(f"[!] Failed to parse NCBI response: {e!s}")
        return None
    except Exception as e:
        print(f"[!] Failed to query NCBI: {e!s}")
        return None


def check_write_permission(directory_path):
    """
    Check if the directory is writable by attempting to create a temp file.

    Args:
        directory_path: Path to the directory to check

    Returns:
        bool: True if writable, False otherwise
    """
    try:
        test_file = Path(directory_path) / f".write_test_{os.getpid()}"
        test_file.touch()
        test_file.unlink()
        return True
    except (PermissionError, OSError):
        return False


def run_cmd(cmd):
    """Run a shell command and raise on failure."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr}")
    return result.stdout


def sleep_with_backoff(attempt, base_seconds=1.0, max_seconds=20.0):
    """Sleep with capped exponential backoff and small jitter."""
    delay = min(max_seconds, base_seconds * (2 ** (attempt - 1)))
    jitter = random.uniform(0.0, 0.25 * delay)
    time.sleep(delay + jitter)


def download_genome_via_api(accession, output_path):
    """Download a genome using NCBI Datasets API (fallback when CLI is not available)."""
    try:
        # Use NCBI Datasets API v2 - POST request with JSON payload
        api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download"

        # Request payload
        payload = {
            "accessions": [accession],
            "include_sequence": ["GENOME_FASTA", "CDS_FASTA", "PROT_FASTA"],
            "host": "human",
            "complete_only": True,
        }

        print(f"    Downloading {accession} via API...")

        # Prepare request
        payload_json = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            api_url, data=payload_json, headers={"Content-Type": "application/json", "Accept": "application/zip"}
        )

        # Submit the download request
        try:
            with urllib.request.urlopen(req, timeout=60) as response:
                content_type = response.headers.get("Content-Type", "").lower()

                # Save to zip file
                zip_path = output_path / f"{accession}.zip"

                # If response is JSON, it might contain an error or redirect URL
                if "application/json" in content_type:
                    data = json.loads(response.read().decode("utf-8"))
                    if "download_url" in data:
                        # Follow the download URL
                        download_url = data["download_url"]
                        print(f"    Following download URL for {accession}...")
                        with urllib.request.urlopen(download_url, timeout=300) as download_response:
                            with open(zip_path, "wb") as f:
                                f.write(download_response.read())
                    elif "error" in data or "message" in data:
                        raise RuntimeError(f"API error: {data.get('error', data.get('message', 'Unknown error'))}")
                    else:
                        raise RuntimeError(f"Unexpected API response format: {data}")
                else:
                    # Write the response content to file
                    with open(zip_path, "wb") as f:
                        f.write(response.read())
        except urllib.error.HTTPError as e:
            raise RuntimeError(f"HTTP error {e.code} for {accession}: {e.reason}") from e
        except urllib.error.URLError as e:
            raise RuntimeError(f"URL error for {accession}: {e!s}") from e

        # Verify zip file is valid and not empty
        if not zip_path.exists() or zip_path.stat().st_size == 0:
            if zip_path.exists():
                zip_path.unlink()
            raise RuntimeError(f"Downloaded zip file is empty for {accession}")

        # Verify it's a valid zip file
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                file_list = zf.namelist()
                if len(file_list) == 0:
                    zip_path.unlink()
                    raise RuntimeError(f"Downloaded zip file contains no files for {accession}")
        except zipfile.BadZipFile as e:
            if zip_path.exists():
                zip_path.unlink()
            raise RuntimeError(f"Downloaded file is not a valid zip for {accession}: {e!s}") from e

        return True
    except Exception as e:
        raise RuntimeError(f"Failed to download {accession}: {e!s}") from e


def download_small_subset(output_dir):
    """Download a small subset of SARS-CoV-2 genomes (~50MB)."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Pre-selected small set of SARS-CoV-2 accessions for ~50MB download
    # These are well-known complete genomes that should total around 50MB
    accessions = [
        "NC_045512.2",  # Reference genome (Wuhan-Hu-1)
        "MT123290.1",  # Early Wuhan-like variant
        "OQ898928.1",  # Alpha (B.1.1.7) variant
        "OR353131.1",  # Beta (B.1.351) variant
        "MW642250.1",  # Gamma (P.1) variant
        "OR323381.1",  # Delta (B.1.617.2) variant
        "OM095411.1",  # Omicron (B.1.1.529) variant
        "PP832909.1",  # Recent (JN.1) variant
    ]

    print(f"[*] Downloading {len(accessions)} SARS-CoV-2 genomes to {output_path}")
    print("[*] Expected total size: ~50MB")

    # Check CLI availability
    use_cli = check_cli_available()
    if not use_cli:
        print("[!] NCBI datasets CLI not found. Using API fallback method.")
        print(
            "[!] For better performance, install the CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/"
        )

    for i, accession in enumerate(accessions, 1):
        print(f"[{i}/{len(accessions)}] Downloading {accession}...")

        if use_cli:
            cmd = (
                f"datasets download virus genome accession {accession} "
                f"--include genome,annotation "
                f"--host human --complete-only "
                f"--filename {output_path}/{accession}.zip --no-progressbar"
            )
            try:
                run_cmd(cmd)
                print(f"    ✓ Downloaded {accession}")
            except RuntimeError as e:
                print(f"    ✗ Failed to download {accession} via CLI: {e}")
                print("    Trying API fallback...")
                try:
                    download_genome_via_api(accession, output_path)
                    print(f"    ✓ Downloaded {accession} via API")
                except RuntimeError as e2:
                    print(f"    ✗ Failed to download {accession}: {e2}")
        else:
            try:
                download_genome_via_api(accession, output_path)
                print(f"    ✓ Downloaded {accession}")
            except RuntimeError as e:
                print(f"    ✗ Failed to download {accession}: {e}")

    print(f"[✓] Download complete. Files saved to: {output_path}")
    return str(output_path)


def download_single_genome(accession, output_path):
    """Download a single genome accession."""
    print(f"Downloading {accession}...")
    max_attempts = 3

    for attempt in range(1, max_attempts + 1):
        if attempt > 1:
            print(f"    Retry {attempt}/{max_attempts} for {accession}...")
            sleep_with_backoff(attempt)

        # Check if CLI is available
        if check_cli_available():
            # Use CLI method
            cmd = (
                f"datasets download virus genome accession {accession} "
                f"--include genome,annotation "
                f"--host human --complete-only "
                f"--filename {output_path}/{accession}.zip --no-progressbar"
            )
            try:
                run_cmd(cmd)
                print(f"    ✓ Downloaded {accession}")
                return accession, True
            except RuntimeError as e:
                print(f"    ✗ Failed to download {accession} via CLI: {e}")
                print("    Trying API fallback...")
                # Fall through to API method
        else:
            print("    CLI not available, using API...")

        # Use API method (fallback or primary if CLI unavailable)
        try:
            download_genome_via_api(accession, output_path)
            print(f"    ✓ Downloaded {accession}")
            return accession, True
        except RuntimeError as e:
            if attempt == max_attempts:
                print(f"    ✗ Failed to download {accession} after {max_attempts} attempts: {e}")
                return accession, False
            print(f"    ✗ Failed to download {accession}: {e}")

    return accession, False


def download_dataset_balanced(virus_name="sars-cov-2", output_dir="data", size_gb=0.05, seed=42, workers=4):
    """
    Download a balanced, random subset of virus genomes with deterministic sampling.
    Checks Google Drive first, downloads to Drive if not available.

    Args:
        virus_name (str): Virus name (e.g., "sars-cov-2")
        output_dir (str): Output directory
        size_gb (float or None): Target size in GB. If None, downloads full dataset (all accessions)
        seed (int): Random seed for deterministic sampling
        workers (int): Number of parallel workers

    Returns:
        str: Path to the downloaded dataset
    """
    # Set random seed for deterministic sampling
    random.seed(seed)

    # Check if running in Google Colab and Drive is mounted
    drive_path = None
    try:
        if os.path.exists("/content/drive/MyDrive"):
            drive_path = Path("/content/drive/MyDrive")
            print(f"[*] Google Drive detected at: {drive_path}")
        else:
            print("[*] Google Drive not mounted, using local directory")
    except OSError:
        print("[*] Not in Google Colab, using local directory")

    # Determine output path (prefer Drive if available)
    if drive_path:
        output_path = drive_path / "COVID-Genome-Classification-NCBI" / output_dir
        print(f"[*] Using Google Drive path: {output_path}")
    else:
        output_path = Path(output_dir)
        print(f"[*] Using local path: {output_path}")

    output_path.mkdir(parents=True, exist_ok=True)

    # Verify write permission before attempting downloads
    if not check_write_permission(output_path):
        print(f"[!] ERROR: No write permission for output directory: {output_path}")
        print("[!] Please check directory permissions or choose a different output path.")
        raise PermissionError(f"Cannot write to output directory: {output_path}")

    # Check for existing files but continue to check individual accessions
    existing_files = list(output_path.glob("*.zip"))
    existing_accessions = {f.stem for f in existing_files}  # Get accession names without .zip
    if existing_files:
        print(f"[*] Found {len(existing_files)} existing genome files:")
        for i, file in enumerate(existing_files, 1):
            print(f"  {i}. {file.name}")
        print("[*] Will check individual accessions and skip already downloaded ones")

    # Skip the slow NCBI query and use pre-selected accessions directly
    print(f"[*] Using pre-selected {virus_name} accessions for faster download...")
    # Use only proven, working SARS-CoV-2 accessions
    # Organized by strain with exactly 20 accessions per strain for balanced datasets
    # All accessions verified via NCBI esearch + efetch + GenBank spot-checks
    accessions_by_strain = {
        "Reference": [
            "NC_045512.2",  # Reference genome (Wuhan-Hu-1) - PROVEN
            "NC_045512.1",  # Reference genome (older version)
            "MT066156.1",  # Early Wuhan-like
            "MT066157.1",  # Early Wuhan-like
            "MT066158.1",  # Early Wuhan-like
            "MT072688.1",  # Early Wuhan-like
            "MT072689.1",  # Early Wuhan-like
            "MT072690.1",  # Early Wuhan-like
            "MT093571.1",  # Early Wuhan-like
            "MT093572.1",  # Early Wuhan-like
            "MT106052.1",  # Early Wuhan-like
            "MT106053.1",  # Early Wuhan-like
            "MT123290.1",  # Early Wuhan-like - PROVEN
            "MT123291.1",  # Early Wuhan-like - PROVEN
            "MT123292.1",  # Early Wuhan-like - PROVEN
            "MT123293.1",  # Early Wuhan-like - PROVEN
            "MT123294.1",  # Early Wuhan-like
            "MT123295.1",  # Early Wuhan-like
            "MT123296.1",  # Early Wuhan-like
            "MT123297.1",  # Early Wuhan-like
        ],
        "Alpha": [
            # Alpha (B.1.1.7) variant - verified via NCBI esearch + GenBank genotype annotations
            "OQ898928.1",  # Alpha (B.1.1.7) - France 2021, genotype: B.1.1.7
            "MZ305033.1",  # Alpha (B.1.1.7) - B.1.1.7 tiger infection study
            "MZ305032.1",  # Alpha (B.1.1.7)
            "MZ305031.1",  # Alpha (B.1.1.7)
            "OV054768.1",  # Alpha (B.1.1.7)
            "OL580734.1",  # Alpha (B.1.1.7)
            "OL461702.1",  # Alpha (B.1.1.7)
            "OL461700.1",  # Alpha (B.1.1.7)
            "OL461699.1",  # Alpha (B.1.1.7)
            "OL461696.1",  # Alpha (B.1.1.7)
            "OL461695.1",  # Alpha (B.1.1.7)
            "OL461694.1",  # Alpha (B.1.1.7)
            "OL461693.1",  # Alpha (B.1.1.7)
            "OL461692.1",  # Alpha (B.1.1.7)
            "OL461691.1",  # Alpha (B.1.1.7)
            "OL461690.1",  # Alpha (B.1.1.7)
            "OL461689.1",  # Alpha (B.1.1.7)
            "OL461688.1",  # Alpha (B.1.1.7)
            "OL461687.1",  # Alpha (B.1.1.7)
            "OL461686.1",  # Alpha (B.1.1.7)
        ],
        "Beta": [
            # Beta (B.1.351) variant - verified via NCBI esearch + GenBank genotype annotations
            "OR353131.1",  # Beta (B.1.351.2) - France 2021, genotype: B.1.351.2
            "OR322476.1",  # Beta (B.1.351)
            "OR322467.1",  # Beta (B.1.351)
            "OR278055.1",  # Beta (B.1.351)
            "OR274549.1",  # Beta (B.1.351)
            "OR274547.1",  # Beta (B.1.351)
            "OQ900811.1",  # Beta (B.1.351)
            "OQ899960.1",  # Beta (B.1.351)
            "OQ897488.1",  # Beta (B.1.351)
            "OQ896567.1",  # Beta (B.1.351)
            "MZ433432.1",  # Beta (B.1.351)
            "MZ427312.1",  # Beta (B.1.351)
            "MZ314997.2",  # Beta (B.1.351)
            "MZ314996.2",  # Beta (B.1.351)
            "MZ314998.1",  # Beta (B.1.351)
            "MZ068161.1",  # Beta (B.1.351)
            "MZ068160.1",  # Beta (B.1.351)
            "MZ068159.1",  # Beta (B.1.351)
            "MZ068158.1",  # Beta (B.1.351)
            "MZ068157.1",  # Beta (B.1.351)
        ],
        "Gamma": [
            # Gamma (P.1) variant - verified via NCBI esearch + GenBank title/annotations
            "MW642250.1",  # Gamma (P.1) - Italy 2021, "P.1 Strains" paper
            "MW642249.1",  # Gamma (P.1)
            "MW642248.1",  # Gamma (P.1)
            "OR701610.1",  # Gamma (P.1)
            "OR322710.1",  # Gamma (P.1)
            "OR322709.1",  # Gamma (P.1)
            "OR322707.1",  # Gamma (P.1)
            "OR322706.1",  # Gamma (P.1)
            "OR279393.1",  # Gamma (P.1)
            "OR279057.1",  # Gamma (P.1)
            "OR273956.1",  # Gamma (P.1)
            "OQ907679.1",  # Gamma (P.1)
            "OQ907601.1",  # Gamma (P.1)
            "OQ905774.1",  # Gamma (P.1)
            "OQ905694.1",  # Gamma (P.1)
            "OQ903601.1",  # Gamma (P.1)
            "OQ903451.1",  # Gamma (P.1)
            "OQ903133.1",  # Gamma (P.1)
            "OQ902994.1",  # Gamma (P.1)
            "OQ900766.1",  # Gamma (P.1)
        ],
        "Delta": [
            # Delta (B.1.617.2) variant - verified via NCBI esearch + GenBank genotype annotations
            "OR323381.1",  # Delta (B.1.617.2) - human, genotype: B.1.617.2
            "OR323348.1",  # Delta (B.1.617.2)
            "OR323317.1",  # Delta (B.1.617.2)
            "OR323300.1",  # Delta (B.1.617.2)
            "OR323282.1",  # Delta (B.1.617.2)
            "OR323281.1",  # Delta (B.1.617.2)
            "OR323276.1",  # Delta (B.1.617.2)
            "OR323275.1",  # Delta (B.1.617.2)
            "OR323273.1",  # Delta (B.1.617.2)
            "OR323213.1",  # Delta (B.1.617.2)
            "OR323144.1",  # Delta (B.1.617.2)
            "OR323039.1",  # Delta (B.1.617.2)
            "OR322908.1",  # Delta (B.1.617.2)
            "OR322800.1",  # Delta (B.1.617.2)
            "OR322776.1",  # Delta (B.1.617.2)
            "OR322732.1",  # Delta (B.1.617.2)
            "OR322731.1",  # Delta (B.1.617.2)
            "OR322714.1",  # Delta (B.1.617.2)
            "OR282403.1",  # Delta (B.1.617.2)
            "OR281251.1",  # Delta (B.1.617.2)
        ],
        "Omicron": [
            # Omicron (B.1.1.529) variant - verified via NCBI esearch + "Sublineage B.1.1.529 Omicron" paper
            "OQ905474.1",  # Omicron (B.1.1.529) - France 2022
            "OM570283.1",  # Omicron (B.1.1.529) - "25 SARS-CoV-2 Sublineage B.1.1.529" paper
            "OM570282.1",  # Omicron (B.1.1.529)
            "OM570281.1",  # Omicron (B.1.1.529)
            "OM570280.1",  # Omicron (B.1.1.529)
            "OM570279.1",  # Omicron (B.1.1.529)
            "OM570278.1",  # Omicron (B.1.1.529)
            "OM570277.1",  # Omicron (B.1.1.529)
            "OM570276.1",  # Omicron (B.1.1.529)
            "OM570275.1",  # Omicron (B.1.1.529)
            "OM570274.1",  # Omicron (B.1.1.529)
            "OM570273.1",  # Omicron (B.1.1.529)
            "OM570272.1",  # Omicron (B.1.1.529)
            "OM570271.1",  # Omicron (B.1.1.529)
            "OM570270.1",  # Omicron (B.1.1.529)
            "OM570269.1",  # Omicron (B.1.1.529)
            "OM570268.1",  # Omicron (B.1.1.529)
            "OM570267.1",  # Omicron (B.1.1.529)
            "OM570266.1",  # Omicron (B.1.1.529)
            "OM570265.1",  # Omicron (B.1.1.529)
        ],
        "Recent": [
            # Recent (JN.1) variant - verified via NCBI esearch + GenBank "JN.1" title annotations
            "PP832909.1",  # Recent (JN.1) - Morocco 2024, "First Report of SARS-CoV-2 JN.1"
            "PP832907.1",  # Recent (JN.1)
            "PP357841.1",  # Recent (JN.1) - France 2024, ICU study
            "PP357840.1",  # Recent (JN.1)
            "PP357816.1",  # Recent (JN.1)
            "PP357814.1",  # Recent (JN.1)
            "PP357812.1",  # Recent (JN.1)
            "PP357805.1",  # Recent (JN.1)
            "PP357803.1",  # Recent (JN.1)
            "PP357797.1",  # Recent (JN.1)
            "PP357786.1",  # Recent (JN.1)
            "PP357785.1",  # Recent (JN.1)
            "PP357778.1",  # Recent (JN.1)
            "PP357768.1",  # Recent (JN.1)
            "PP357764.1",  # Recent (JN.1)
            "PP357762.1",  # Recent (JN.1)
            "PP357759.1",  # Recent (JN.1)
            "PP357753.1",  # Recent (JN.1)
            "PP357734.1",  # Recent (JN.1)
            "PP357732.1",  # Recent (JN.1)
        ],
    }

    # Flatten to single list
    accessions = []
    for _strain, strain_accessions in accessions_by_strain.items():
        accessions.extend(strain_accessions)

    # Print strain distribution
    print("[*] Strain distribution in pre-selected accessions:")
    for strain, strain_accessions in accessions_by_strain.items():
        print(f"  - {strain}: {len(strain_accessions)} accessions")

    # Sample random subset based on target size (or use all for full dataset)
    if size_gb is None:
        # Full dataset mode: use all accessions
        selected_accessions = accessions
        print(f"[*] Full dataset mode: using all {len(accessions)} accessions")
    else:
        # Partial dataset mode: sample based on target size
        # Estimate ~6MB per genome, so for 50MB we need ~8 genomes
        target_count = max(8, int(size_gb * 1024 / 6))  # At least 8, or based on size
        if len(accessions) > target_count:
            selected_accessions = random.sample(accessions, target_count)
            print(f"[*] Partial dataset mode: sampling {target_count} from {len(accessions)} accessions")
        else:
            selected_accessions = accessions
            print(f"[*] Partial dataset mode: using all {len(accessions)} accessions (less than target)")

    # Keep deterministic ordering for consistent behavior and logging across runs
    selected_accessions = sorted(selected_accessions)

    # Filter out already downloaded accessions
    accessions_to_download = [acc for acc in selected_accessions if acc not in existing_accessions]
    already_downloaded = [acc for acc in selected_accessions if acc in existing_accessions]

    print("[*] Accession status:")
    print(f"  - Already downloaded: {len(already_downloaded)} ({already_downloaded})")
    print(f"  - Need to download: {len(accessions_to_download)} ({accessions_to_download})")

    if not accessions_to_download:
        print("[*] All selected accessions already downloaded, skipping download")
        return str(output_path)

    if size_gb is None:
        print(f"[*] Downloading {len(accessions_to_download)} {virus_name} genomes (full dataset)")
    else:
        print(f"[*] Downloading {len(accessions_to_download)} {virus_name} genomes (target: {size_gb}GB)")
    print(f"[*] Using {workers} workers with seed {seed}")

    # Check CLI availability and inform user
    if check_cli_available():
        print("[*] Using NCBI datasets CLI")
    else:
        print("[!] NCBI datasets CLI not found. Using API fallback method.")
        print(
            "[!] For better performance, install the CLI: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/"
        )

    # Download with parallel workers (only missing accessions)
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(download_single_genome, accession, output_path) for accession in accessions_to_download
        ]

        for future in concurrent.futures.as_completed(futures):
            accession, success = future.result()
            results.append((accession, success))

    # Summary
    successful = sum(1 for _, success in results if success)
    failed = len(results) - successful
    total_files = len(existing_accessions) + successful
    failed_accessions = sorted([accession for accession, success in results if not success])

    print("[✓] Download complete:")
    print(f"  - Already had: {len(already_downloaded)} files")
    print(f"  - Newly downloaded: {successful} successful, {failed} failed")
    print(f"  - Total files now: {total_files}")
    if failed_accessions:
        print(f"  - Failed accessions: {failed_accessions}")
    print(f"[*] Files saved to: {output_path}")

    return str(output_path)


def get_default_output_path():
    """
    Detect the best default output path based on the platform.
    Checks for Google Drive/Colab paths, then falls back to local directory.
    """
    import platform

    # Check for Google Colab first
    if os.path.exists("/content/drive/MyDrive"):
        return "/content/drive/MyDrive/ncbi_data"

    # Platform-specific Google Drive paths
    system = platform.system()

    if system == "Windows":
        # Common Windows Google Drive paths
        windows_paths = [
            os.path.expanduser("~/Google Drive/ncbi_data"),
            "G:/My Drive/ncbi_data",
            "D:/Google Drive/ncbi_data",
        ]
        for path in windows_paths:
            parent = os.path.dirname(path)
            if os.path.exists(parent):
                return path
    elif system == "Darwin":  # macOS
        # Common macOS Google Drive paths
        mac_paths = [
            os.path.expanduser("~/Library/CloudStorage/GoogleDrive-*/My Drive/ncbi_data"),
            os.path.expanduser("~/Google Drive/ncbi_data"),
        ]
        for pattern in mac_paths:
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
            # Check if parent exists for non-glob paths
            if "*" not in pattern:
                parent = os.path.dirname(pattern)
                if os.path.exists(parent):
                    return pattern
    elif system == "Linux":
        # Linux Google Drive paths (via google-drive-ocamlfuse or similar)
        linux_paths = [
            os.path.expanduser("~/google-drive/ncbi_data"),
            os.path.expanduser("~/Google Drive/ncbi_data"),
        ]
        for path in linux_paths:
            parent = os.path.dirname(path)
            if os.path.exists(parent):
                return path

    # Fallback to local directory
    return "./ncbi_data"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download small SARS-CoV-2 dataset from NCBI")

    default_output = get_default_output_path()
    parser.add_argument(
        "--output",
        type=str,
        default=default_output,
        help=f"Output directory (default: {default_output})",
    )
    args = parser.parse_args()

    # Ensure output path is valid or use fallback
    if not os.path.exists(args.output):
        parent = os.path.dirname(args.output) if args.output != "./ncbi_data" else "."
        if parent and not os.path.exists(parent):
            print(f"[!] Output directory parent does not exist: {parent}")
            args.output = "./ncbi_data"
            print(f"[*] Using fallback: {args.output}")
        else:
            print(f"[*] Output directory will be created: {args.output}")

    download_small_subset(args.output)
