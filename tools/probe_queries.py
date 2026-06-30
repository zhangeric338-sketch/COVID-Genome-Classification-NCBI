#!/usr/bin/env python3
"""Probe Beta/Gamma queries - they're rare in GenBank under human+complete filters."""
import os
import sys
import time
from Bio import Entrez

Entrez.email = sys.argv[1] if len(sys.argv) > 1 else "probe@example.com"
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")

BASE_STRICT = 'txid2697049[Organism] AND "complete genome"[Title] AND "Homo sapiens"[Host]'
BASE_NOHOST = 'txid2697049[Organism] AND "complete genome"[Title]'
BASE_LOOSE = 'txid2697049[Organism]'

TESTS = {
    "Beta strict B.1.351": f'{BASE_STRICT} AND "B.1.351"[All Fields]',
    "Beta no-host": f'{BASE_NOHOST} AND "B.1.351"[All Fields]',
    "Beta loose": f'{BASE_LOOSE} AND "B.1.351"[All Fields]',
    "Beta sublineages": f'{BASE_STRICT} AND (B.1.351[All Fields] OR B.1.351.1[All Fields] OR B.1.351.2[All Fields] OR B.1.351.3[All Fields])',
    "Beta date window": f'{BASE_STRICT} AND 2020/10/01:2021/09/30[PDAT] AND ("B.1.351"[All Fields] OR beta[Title])',

    "Gamma strict P.1": f'{BASE_STRICT} AND "P.1"[All Fields]',
    "Gamma no-host": f'{BASE_NOHOST} AND "P.1"[All Fields]',
    "Gamma loose": f'{BASE_LOOSE} AND "P.1"[All Fields]',
    "Gamma sublineages": f'{BASE_STRICT} AND (P.1[All Fields] OR P.1.1[All Fields] OR P.1.2[All Fields])',
    "Gamma date window": f'{BASE_STRICT} AND 2020/12/01:2021/09/30[PDAT] AND ("P.1"[All Fields] OR gamma[Title])',
}

for label, q in TESTS.items():
    time.sleep(0.12)
    try:
        h = Entrez.esearch(db="nucleotide", term=q, retmax=0)
        r = Entrez.read(h); h.close()
        count = int(r["Count"])
        flag = " OK" if count >= 143 else ""
        print(f"  {count:8d}  {label}{flag}")
    except Exception as e:
        print(f"  ERROR  {label}: {e}")
