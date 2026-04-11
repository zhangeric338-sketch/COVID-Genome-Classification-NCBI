#!/usr/bin/env python3
"""
tools/download_data.py
Download a balanced subset of SARS-CoV-2 genomes from NCBI.

Uses the download engine (CLI primary, Entrez fallback) for reliable downloads
with persistent state, rate limiting, and automatic retry.

Usage:
    python tools/download_data.py [--output ./data]
"""

import argparse
import json
import random
import urllib.error
import urllib.request
from pathlib import Path

from tools.download_engine import run_download


def query_ncbi_genome_count(virus_name="sars-cov-2", host="human", complete_only=True):
    """
    Query NCBI for total count and estimated size of available genomes.

    Returns:
        dict with total_count, estimated_size_mb/gb/tb, or None on failure.
    """
    try:
        api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome"
        body = {"taxon": virus_name, "host": host}
        if complete_only:
            body["complete_only"] = True

        payload_json = json.dumps(body).encode("utf-8")
        req = urllib.request.Request(
            api_url,
            data=payload_json,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )

        print(f"[*] Querying NCBI for {virus_name} genome count...")
        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.loads(response.read().decode("utf-8"))
            total_count = data.get("record_count", 0)
            estimated_size_mb = total_count * 6
            return {
                "total_count": total_count,
                "estimated_size_mb": estimated_size_mb,
                "estimated_size_gb": estimated_size_mb / 1024,
                "estimated_size_tb": estimated_size_mb / 1024 / 1024,
            }
    except (urllib.error.HTTPError, urllib.error.URLError, json.JSONDecodeError, Exception) as e:
        print(f"[!] Failed to query NCBI: {e}")
        return None


# ---------------------------------------------------------------------------
# Curated accessions — 7 strains x 20 = 140 verified genomes
# ---------------------------------------------------------------------------

ACCESSIONS_BY_STRAIN = {
    "Reference": [
        "NC_045512.2",  # Reference genome (Wuhan-Hu-1)
        "NC_045512.1",  # Reference genome (older version)
        "MT066156.1",
        "MT066157.1",
        "MT066158.1",
        "MT072688.1",
        "MT072689.1",
        "MT072690.1",
        "MT093571.1",
        "MT093572.1",
        "MT106052.1",
        "MT106053.1",
        "MT123290.1",
        "MT123291.1",
        "MT123292.1",
        "MT123293.1",
        "MT123294.1",
        "MT123295.1",
        "MT123296.1",
        "MT123297.1",
    ],
    "Alpha": [
        "OQ898928.1",  # B.1.1.7 - France 2021
        "MZ305033.1",
        "MZ305032.1",
        "MZ305031.1",
        "OV054768.1",
        "OL580734.1",
        "OL461702.1",
        "OL461700.1",
        "OL461699.1",
        "OL461696.1",
        "OL461695.1",
        "OL461694.1",
        "OL461693.1",
        "OL461692.1",
        "OL461691.1",
        "OL461690.1",
        "OL461689.1",
        "OL461688.1",
        "OL461687.1",
        "OL461686.1",
    ],
    "Beta": [
        "OR353131.1",  # B.1.351.2 - France 2021
        "OR322476.1",
        "OR322467.1",
        "OR278055.1",
        "OR274549.1",
        "OR274547.1",
        "OQ900811.1",
        "OQ899960.1",
        "OQ897488.1",
        "OQ896567.1",
        "MZ433432.1",
        "MZ427312.1",
        "MZ314997.2",
        "MZ314996.2",
        "MZ314998.1",
        "MZ068161.1",
        "MZ068160.1",
        "MZ068159.1",
        "MZ068158.1",
        "MZ068157.1",
    ],
    "Gamma": [
        "MW642250.1",  # P.1 - Italy 2021
        "MW642249.1",
        "MW642248.1",
        "OR701610.1",
        "OR322710.1",
        "OR322709.1",
        "OR322707.1",
        "OR322706.1",
        "OR279393.1",
        "OR279057.1",
        "OR273956.1",
        "OQ907679.1",
        "OQ907601.1",
        "OQ905774.1",
        "OQ905694.1",
        "OQ903601.1",
        "OQ903451.1",
        "OQ903133.1",
        "OQ902994.1",
        "OQ900766.1",
    ],
    "Delta": [
        "OR323381.1",  # B.1.617.2
        "OR323348.1",
        "OR323317.1",
        "OR323300.1",
        "OR323282.1",
        "OR323281.1",
        "OR323276.1",
        "OR323275.1",
        "OR323273.1",
        "OR323213.1",
        "OR323144.1",
        "OR323039.1",
        "OR322908.1",
        "OR322800.1",
        "OR322776.1",
        "OR322732.1",
        "OR322731.1",
        "OR322714.1",
        "OR282403.1",
        "OR281251.1",
    ],
    "Omicron": [
        "OQ905474.1",  # B.1.1.529 - France 2022
        "OM570283.1",
        "OM570282.1",
        "OM570281.1",
        "OM570280.1",
        "OM570279.1",
        "OM570278.1",
        "OM570277.1",
        "OM570276.1",
        "OM570275.1",
        "OM570274.1",
        "OM570273.1",
        "OM570272.1",
        "OM570271.1",
        "OM570270.1",
        "OM570269.1",
        "OM570268.1",
        "OM570267.1",
        "OM570266.1",
        "OM570265.1",
    ],
    "Recent": [
        "PP832909.1",  # JN.1 - Morocco 2024
        "PP832907.1",
        "PP357841.1",
        "PP357840.1",
        "PP357816.1",
        "PP357814.1",
        "PP357812.1",
        "PP357805.1",
        "PP357803.1",
        "PP357797.1",
        "PP357786.1",
        "PP357785.1",
        "PP357778.1",
        "PP357768.1",
        "PP357764.1",
        "PP357762.1",
        "PP357759.1",
        "PP357753.1",
        "PP357734.1",
        "PP357732.1",
    ],
}


def download_dataset_balanced(virus_name="sars-cov-2", output_dir="data", size_gb=0.05, seed=42, workers=4):
    """
    Download a balanced, random subset of virus genomes with deterministic sampling.

    Args:
        virus_name: Virus name (for logging only)
        output_dir: Output directory
        size_gb: Target size in GB (None = full dataset)
        seed: Random seed for deterministic sampling
        workers: Max parallel workers for Entrez fallback

    Returns:
        str: Path to the downloaded dataset
    """
    random.seed(seed)
    output_path = Path(output_dir)

    # Flatten all accessions
    all_accessions = []
    for strain_accs in ACCESSIONS_BY_STRAIN.values():
        all_accessions.extend(strain_accs)

    print(f"[*] Strain distribution ({len(all_accessions)} total):")
    for strain, accs in ACCESSIONS_BY_STRAIN.items():
        print(f"  - {strain}: {len(accs)}")

    # Sample based on target size (or use all)
    if size_gb is None:
        selected = set(all_accessions)
        print(f"[*] Full dataset mode: using all {len(selected)} accessions")
    else:
        target_count = max(8, int(size_gb * 1024 / 6))
        if len(all_accessions) > target_count:
            selected = set(random.sample(all_accessions, target_count))
            print(f"[*] Sampling {target_count} accessions (target: {size_gb} GB)")
        else:
            selected = set(all_accessions)

    # Build filtered accessions_by_strain for the selected subset
    selected_by_strain = {}
    for strain, accs in ACCESSIONS_BY_STRAIN.items():
        filtered = [a for a in accs if a in selected]
        if filtered:
            selected_by_strain[strain] = filtered

    return run_download(
        accessions_by_strain=selected_by_strain,
        output_dir=str(output_path),
        max_entrez_workers=min(workers, 3),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SARS-CoV-2 dataset from NCBI")
    parser.add_argument("--output", type=str, default="data", help="Output directory (default: data)")
    parser.add_argument("--size-gb", type=float, default=None, help="Target size in GB (default: full dataset)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    download_dataset_balanced(output_dir=args.output, size_gb=args.size_gb, seed=args.seed)
