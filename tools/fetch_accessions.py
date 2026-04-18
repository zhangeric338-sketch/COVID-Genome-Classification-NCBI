#!/usr/bin/env python3
"""
tools/fetch_accessions.py
Query NCBI via Entrez E-utilities to find ~1000 balanced SARS-CoV-2 accessions
(complete genomes, human host) across 7 strain categories.

Usage:
    python tools/fetch_accessions.py --email you@example.com [--per-strain 143] [--seed 42]
"""

import argparse
import json
import os
import random
import time
from pathlib import Path

from Bio import Entrez


# Strain definitions: name -> list of Entrez search terms for Pango lineages
# Each strain may have multiple lineage patterns to search for
STRAIN_QUERIES = {
    "Reference": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND ("Wuhan-Hu-1"[Title] OR "WIV04"[Title])',
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND 2020/01/01:2020/03/31[Collection Date] AND "A"[Lineage]',
    ],
    "Alpha": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND "B.1.1.7"[All Fields]',
    ],
    "Beta": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND "B.1.351"[All Fields]',
    ],
    "Gamma": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND "P.1"[All Fields]',
    ],
    "Delta": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND "B.1.617.2"[All Fields]',
    ],
    "Omicron": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND "B.1.1.529"[All Fields]',
    ],
    "Recent": [
        'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
        ' AND "JN.1"[All Fields]',
    ],
}


def search_strain(strain_name, queries, max_results=1000):
    """
    Search NCBI Nucleotide for accessions matching a strain's queries.

    Args:
        strain_name: Name of the strain (for logging)
        queries: List of Entrez search query strings
        max_results: Maximum IDs to retrieve per query

    Returns:
        set of NCBI nucleotide IDs
    """
    all_ids = set()
    for query in queries:
        print(f"  Searching: {query[:80]}...")
        time.sleep(0.34)  # Rate limit: ~3 requests/sec without API key
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        count = int(record["Count"])
        ids = record["IdList"]
        print(f"  Found {count} total, retrieved {len(ids)} IDs")
        all_ids.update(ids)
    print(f"  [{strain_name}] Total unique IDs: {len(all_ids)}")
    return all_ids


def fetch_accessions_for_ids(id_list):
    """
    Convert NCBI nucleotide IDs to accession.version strings via efetch.

    Args:
        id_list: List of NCBI IDs (strings)

    Returns:
        list of accession strings (e.g., ["NC_045512.2", "MT123290.1"])
    """
    accessions = []
    batch_size = 200
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i : i + batch_size]
        time.sleep(0.34)
        handle = Entrez.efetch(db="nucleotide", id=batch, rettype="acc", retmode="text")
        text = handle.read()
        handle.close()
        batch_accessions = [line.strip() for line in text.strip().split("\n") if line.strip()]
        accessions.extend(batch_accessions)
    return accessions


def fetch_all_accessions(per_strain=143, seed=42):
    """
    Search NCBI for balanced SARS-CoV-2 accessions across all strains.

    Args:
        per_strain: Target number of accessions per strain
        seed: Random seed for reproducible sampling

    Returns:
        dict with strain -> list of accession strings
    """
    random.seed(seed)
    results = {}

    for strain_name, queries in STRAIN_QUERIES.items():
        print(f"\n[*] Searching for {strain_name} strain...")
        ids = search_strain(strain_name, queries)

        if not ids:
            print(f"  [!] WARNING: No results found for {strain_name}")
            results[strain_name] = {
                "query": queries,
                "total_found": 0,
                "sampled": 0,
                "accessions": [],
            }
            continue

        # Sample IDs if we have more than needed
        id_list = list(ids)
        if len(id_list) > per_strain:
            sampled_ids = random.sample(id_list, per_strain)
        else:
            sampled_ids = id_list
            if len(id_list) < per_strain:
                print(
                    f"  [!] WARNING: Only {len(id_list)} results for {strain_name}, "
                    f"target was {per_strain}"
                )

        # Convert IDs to accession strings
        print(f"  Fetching accession names for {len(sampled_ids)} IDs...")
        accessions = fetch_accessions_for_ids(sampled_ids)

        results[strain_name] = {
            "query": queries,
            "total_found": len(ids),
            "sampled": len(accessions),
            "accessions": sorted(accessions),
        }
        print(f"  [{strain_name}] Sampled {len(accessions)} accessions")

    return results


def save_results(results, output_dir, seed, email):
    """Save accession results to JSON and flat text file."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Build summary
    total = sum(r["sampled"] for r in results.values())
    summary = {
        "generated": time.strftime("%Y-%m-%d"),
        "seed": seed,
        "email": email,
        "total": total,
        "per_strain_target": max(r["sampled"] for r in results.values()) if results else 0,
        "strains": results,
    }

    # Save JSON
    json_path = output_path / "accessions.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[OK] Saved accession metadata to: {json_path}")

    # Save flat accession list
    txt_path = output_path / "accession_list.txt"
    all_accessions = []
    for strain_data in results.values():
        all_accessions.extend(strain_data["accessions"])
    with open(txt_path, "w") as f:
        f.write("\n".join(sorted(all_accessions)) + "\n")
    print(f"[OK] Saved flat accession list to: {txt_path} ({len(all_accessions)} accessions)")

    # Save strain manifest (accession -> strain mapping)
    manifest = {}
    for strain_name, strain_data in results.items():
        for acc in strain_data["accessions"]:
            manifest[acc] = strain_name
    manifest_path = output_path / "strain_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[OK] Saved strain manifest to: {manifest_path}")

    return json_path, txt_path, manifest_path


def print_summary(results):
    """Print a summary table of the search results."""
    print("\n" + "=" * 60)
    print("Accession Search Summary")
    print("=" * 60)
    total = 0
    for strain_name, data in results.items():
        found = data["total_found"]
        sampled = data["sampled"]
        total += sampled
        print(f"  {strain_name:12s}: {sampled:4d} sampled  (from {found:,} found)")
    print("-" * 60)
    print(f"  {'TOTAL':12s}: {total:4d} accessions")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch balanced SARS-CoV-2 accessions from NCBI via Entrez"
    )
    parser.add_argument(
        "--email",
        type=str,
        default=os.environ.get("ENTREZ_EMAIL", ""),
        help="Email for NCBI Entrez (required by NCBI). Can also set ENTREZ_EMAIL env var.",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key for higher rate limits (optional). Can also set NCBI_API_KEY env var.",
    )
    parser.add_argument(
        "--per-strain",
        type=int,
        default=143,
        help="Target accessions per strain (default: 143, for ~1001 total across 7 strains)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible sampling (default: 42)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="tools",
        help="Directory to save accessions.json and accession_list.txt (default: tools)",
    )

    args = parser.parse_args()

    if not args.email:
        parser.error(
            "Email is required by NCBI Entrez. Use --email or set ENTREZ_EMAIL env var."
        )

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
        print(f"[*] Using API key (10 requests/sec limit)")
    else:
        print("[*] No API key set (3 requests/sec limit). Set NCBI_API_KEY for faster queries.")

    print(f"[*] Searching NCBI for {args.per_strain} accessions per strain (7 strains)")
    print(f"[*] Target total: ~{args.per_strain * 7} accessions")
    print(f"[*] Random seed: {args.seed}")

    results = fetch_all_accessions(per_strain=args.per_strain, seed=args.seed)
    print_summary(results)
    save_results(results, args.output_dir, args.seed, args.email)


if __name__ == "__main__":
    main()
