#!/usr/bin/env python3
"""
tools/download_genomes.py
Bulk download SARS-CoV-2 genomes from NCBI using Datasets CLI (preferred)
or API fallback. Reads accessions from accessions.json produced by fetch_accessions.py.

Outputs per-accession ZIP files compatible with the existing training pipeline.

Usage:
    python tools/download_genomes.py [--accessions-file tools/accessions.json] [--output-dir data]
"""

import argparse
import io
import json
import os
import random
import shutil
import subprocess
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def check_cli_available():
    """Check if NCBI datasets CLI is available."""
    return shutil.which("datasets") is not None


def install_cli_colab():
    """Attempt to install NCBI datasets CLI in Google Colab."""
    if not os.path.exists("/content"):
        return False
    print("[*] Detected Google Colab. Installing NCBI datasets CLI...")
    try:
        subprocess.run(
            "curl -o /usr/local/bin/datasets "
            "https://ftp.ncbi.nlm.nih.gov/datasets/cli/v2/linux-amd64/datasets "
            "&& chmod +x /usr/local/bin/datasets",
            shell=True,
            check=True,
            capture_output=True,
        )
        print("[✓] NCBI datasets CLI installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[!] Failed to install CLI: {e}")
        return False


def load_accessions(accessions_file):
    """Load accessions from JSON file produced by fetch_accessions.py."""
    with open(accessions_file) as f:
        data = json.load(f)

    all_accessions = []
    strain_map = {}
    for strain_name, strain_data in data["strains"].items():
        for acc in strain_data["accessions"]:
            all_accessions.append(acc)
            strain_map[acc] = strain_name

    return sorted(all_accessions), strain_map, data


def download_bulk_cli(accessions, output_path, batch_size=200):
    """
    Download genomes in bulk using NCBI datasets CLI, then repackage
    the multi-FASTA output into per-accession ZIPs.

    Returns:
        tuple: (successful_accessions, failed_accessions)
    """
    successful = []
    failed = []

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(accessions) + batch_size - 1) // batch_size
        print(f"\n[*] Batch {batch_num}/{total_batches}: downloading {len(batch)} genomes via CLI...")

        # Write batch accessions to temp file for --inputfile
        acc_file = output_path / f"_batch_{batch_num}.txt"
        with open(acc_file, "w") as f:
            f.write("\n".join(batch) + "\n")

        bulk_zip = output_path / f"_batch_{batch_num}.zip"

        cmd = (
            f"datasets download virus genome accession "
            f"--inputfile {acc_file} "
            f"--include genome "
            f"--filename {bulk_zip} --no-progressbar"
        )

        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
            if result.returncode != 0:
                print(f"  [!] CLI download failed: {result.stderr}")
                failed.extend(batch)
                continue
        except subprocess.TimeoutExpired:
            print(f"  [!] CLI download timed out for batch {batch_num}")
            failed.extend(batch)
            continue

        # Repackage bulk ZIP into per-accession ZIPs
        batch_successful, batch_failed = repackage_bulk_zip(bulk_zip, batch, output_path)
        successful.extend(batch_successful)
        failed.extend(batch_failed)

        # Cleanup temp files
        if acc_file.exists():
            acc_file.unlink()
        if bulk_zip.exists():
            bulk_zip.unlink()

    return successful, failed


def repackage_bulk_zip(bulk_zip_path, expected_accessions, output_path):
    """
    Extract a bulk CLI ZIP (multi-record genomic.fna) into per-accession ZIPs.

    Returns:
        tuple: (successful_accessions, failed_accessions)
    """
    successful = []
    found_accessions = set()

    try:
        with zipfile.ZipFile(bulk_zip_path, "r") as zf:
            # Find the genomic.fna file
            fna_path = None
            for name in zf.namelist():
                if name.endswith("genomic.fna"):
                    fna_path = name
                    break

            if not fna_path:
                print("  [!] No genomic.fna found in bulk ZIP")
                return [], expected_accessions

            with zf.open(fna_path) as fna_file:
                fna_text = io.TextIOWrapper(fna_file)
                if SeqIO is None:
                    # Fallback: parse FASTA manually
                    records = _parse_fasta_manual(fna_text.read())
                else:
                    records = list(SeqIO.parse(fna_text, "fasta"))

            for record in records:
                # Extract accession from FASTA header
                # Header format: >accession.version description
                if SeqIO is not None:
                    acc = record.id
                    seq_str = str(record.seq)
                    description = record.description
                else:
                    acc = record["id"]
                    seq_str = record["seq"]
                    description = record["description"]

                found_accessions.add(acc)

                # Skip if already exists
                per_acc_zip = output_path / f"{acc}.zip"
                if per_acc_zip.exists():
                    successful.append(acc)
                    continue

                # Create per-accession ZIP
                fasta_content = f">{description}\n{seq_str}\n"
                buf = io.BytesIO()
                with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as out_zip:
                    out_zip.writestr("ncbi_dataset/data/genomic.fna", fasta_content)
                buf.seek(0)
                with open(per_acc_zip, "wb") as f:
                    f.write(buf.read())

                successful.append(acc)

    except (zipfile.BadZipFile, Exception) as e:
        print(f"  [!] Error processing bulk ZIP: {e}")

    failed = [acc for acc in expected_accessions if acc not in found_accessions]
    if failed:
        print(f"  [!] {len(failed)} accessions not found in bulk download")

    print(f"  [✓] Repackaged {len(successful)} genomes into per-accession ZIPs")
    return successful, failed


def _parse_fasta_manual(text):
    """Simple FASTA parser fallback when BioPython is not available."""
    records = []
    current_header = None
    current_seq = []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                acc = current_header.split()[0]
                records.append({
                    "id": acc,
                    "description": current_header,
                    "seq": "".join(current_seq),
                })
            current_header = line[1:]
            current_seq = []
        elif line:
            current_seq.append(line)
    if current_header:
        acc = current_header.split()[0]
        records.append({
            "id": acc,
            "description": current_header,
            "seq": "".join(current_seq),
        })
    return records


def download_single_api(accession, output_path):
    """Download a single genome via NCBI Datasets API v2."""
    zip_path = output_path / f"{accession}.zip"
    if zip_path.exists():
        return accession, True

    api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download"
    payload = {
        "accessions": [accession],
        "include_sequence": ["GENOME_FASTA"],
    }

    max_attempts = 3
    for attempt in range(1, max_attempts + 1):
        try:
            payload_json = json.dumps(payload).encode("utf-8")
            req = urllib.request.Request(
                api_url,
                data=payload_json,
                headers={"Content-Type": "application/json", "Accept": "application/zip"},
            )
            with urllib.request.urlopen(req, timeout=60) as response:
                content = response.read()

            # Save and validate
            with open(zip_path, "wb") as f:
                f.write(content)

            with zipfile.ZipFile(zip_path, "r") as zf:
                if not zf.namelist():
                    zip_path.unlink()
                    raise RuntimeError("Empty ZIP")

            return accession, True

        except Exception as e:
            if attempt < max_attempts:
                delay = min(20, 1.0 * (2 ** (attempt - 1)))
                jitter = random.uniform(0, 0.25 * delay)
                time.sleep(delay + jitter)
            else:
                if zip_path.exists():
                    zip_path.unlink()
                print(f"  [!] Failed {accession} after {max_attempts} attempts: {e}")
                return accession, False

    return accession, False


def download_api_fallback(accessions, output_path, workers=4):
    """Download genomes one at a time via API with parallel workers."""
    import concurrent.futures

    print(f"\n[*] Downloading {len(accessions)} genomes via API ({workers} workers)...")

    successful = []
    failed = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_single_api, acc, output_path): acc
            for acc in accessions
        }
        done_count = 0
        for future in concurrent.futures.as_completed(futures):
            acc, success = future.result()
            done_count += 1
            if success:
                successful.append(acc)
            else:
                failed.append(acc)
            if done_count % 50 == 0:
                print(f"  Progress: {done_count}/{len(accessions)}")

    return successful, failed


def download_genomes(accessions_file, output_dir, workers=4, batch_size=200):
    """
    Main download pipeline: reads accessions.json, downloads genomes,
    writes strain_manifest.json to output directory.
    """
    print(f"[*] Loading accessions from: {accessions_file}")
    all_accessions, strain_map, metadata = load_accessions(accessions_file)
    print(f"[*] Total accessions to download: {len(all_accessions)}")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Check which accessions already exist
    existing = {f.stem for f in output_path.glob("*.zip") if not f.name.startswith("_")}
    to_download = [acc for acc in all_accessions if acc not in existing]
    already_have = len(all_accessions) - len(to_download)

    if already_have > 0:
        print(f"[*] Already downloaded: {already_have}")
    if not to_download:
        print("[✓] All accessions already downloaded")
    else:
        print(f"[*] Need to download: {len(to_download)}")

        # Try CLI first
        use_cli = check_cli_available()
        if not use_cli and os.path.exists("/content"):
            use_cli = install_cli_colab()

        if use_cli:
            print("[*] Using NCBI datasets CLI for bulk download")
            successful, failed = download_bulk_cli(to_download, output_path, batch_size)
        else:
            print("[!] NCBI datasets CLI not found. Using API fallback (slower).")
            successful, failed = [], to_download

        # Retry failed accessions via API
        if failed:
            print(f"\n[*] Retrying {len(failed)} failed accessions via API...")
            api_success, api_failed = download_api_fallback(failed, output_path, workers)
            successful.extend(api_success)
            failed = api_failed

        # Final summary
        print(f"\n[✓] Download complete:")
        print(f"  Already had: {already_have}")
        print(f"  Newly downloaded: {len(successful)}")
        print(f"  Failed: {len(failed)}")
        if failed:
            print(f"  Failed accessions: {failed[:20]}{'...' if len(failed) > 20 else ''}")

    # Write strain manifest to data directory
    manifest_path = output_path / "strain_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(strain_map, f, indent=2)
    print(f"[✓] Strain manifest saved to: {manifest_path}")

    # Validate final state
    final_zips = {f.stem for f in output_path.glob("*.zip") if not f.name.startswith("_")}
    matched = final_zips & set(all_accessions)
    print(f"[✓] Final dataset: {len(matched)}/{len(all_accessions)} accessions available")

    return str(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Bulk download SARS-CoV-2 genomes from NCBI"
    )
    parser.add_argument(
        "--accessions-file",
        type=str,
        default="tools/accessions.json",
        help="Path to accessions.json from fetch_accessions.py (default: tools/accessions.json)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data",
        help="Output directory for genome ZIPs (default: data)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Parallel workers for API fallback (default: 4)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="Accessions per CLI batch (default: 200)",
    )

    args = parser.parse_args()

    if not os.path.exists(args.accessions_file):
        print(f"[!] Accessions file not found: {args.accessions_file}")
        print("[!] Run fetch_accessions.py first to generate it.")
        return

    download_genomes(
        accessions_file=args.accessions_file,
        output_dir=args.output_dir,
        workers=args.workers,
        batch_size=args.batch_size,
    )


if __name__ == "__main__":
    main()
