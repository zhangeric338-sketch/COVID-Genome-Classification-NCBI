#!/usr/bin/env python3
"""
tools/download_genomes.py
Bulk download SARS-CoV-2 genomes from NCBI using accessions.json
produced by fetch_accessions.py.

Uses the download engine (CLI primary, Entrez fallback) for reliable downloads.

Usage:
    python tools/download_genomes.py [--accessions-file tools/accessions.json] [--output-dir data]
"""

import argparse
import json
import os
from pathlib import Path

from tools.download_engine import run_download


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


def download_genomes(accessions_file, output_dir, workers=4, batch_size=200):
    """
    Main download pipeline: reads accessions.json, downloads genomes,
    writes strain_manifest.json to output directory.
    """
    print(f"[*] Loading accessions from: {accessions_file}")
    all_accessions, strain_map, _metadata = load_accessions(accessions_file)
    print(f"[*] Total accessions to download: {len(all_accessions)}")

    # Invert strain_map to accessions_by_strain
    accessions_by_strain = {}
    for acc, strain in strain_map.items():
        accessions_by_strain.setdefault(strain, []).append(acc)

    output_path = Path(output_dir)
    dataset_path = run_download(
        accessions_by_strain=accessions_by_strain,
        output_dir=str(output_path),
        batch_size=batch_size,
        max_entrez_workers=min(workers, 3),
    )

    # Write strain manifest (accession -> strain) for visualization/training
    manifest_path = output_path / "strain_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(strain_map, f, indent=2)
    print(f"[OK] Strain manifest saved to: {manifest_path}")

    return dataset_path


def main():
    parser = argparse.ArgumentParser(description="Bulk download SARS-CoV-2 genomes from NCBI")
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
        help="Parallel workers for Entrez fallback (default: 4)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Accessions per CLI batch (default: 50)",
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
