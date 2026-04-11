#!/usr/bin/env python3
"""
tools/download_engine.py
Robust genome download engine with persistent state, rate limiting,
and automatic fallback from NCBI CLI to Entrez efetch.

Used by download_data.py and download_genomes.py — not called directly.
"""

import io
import json
import os
import random
import shutil
import subprocess
import threading
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


# ---------------------------------------------------------------------------
# Download manifest — persistent JSON state for resume / crash recovery
# ---------------------------------------------------------------------------


class DownloadManifest:
    """Track per-accession download status in a JSON file."""

    def __init__(self, path):
        self._path = Path(path)
        self._lock = threading.Lock()
        self._data = {"version": 1, "accessions": {}}

    # -- persistence --

    @classmethod
    def load(cls, output_dir):
        manifest_path = Path(output_dir) / ".download_manifest.json"
        m = cls(manifest_path)
        if manifest_path.exists():
            with open(manifest_path) as f:
                m._data = json.load(f)
        return m

    def save(self):
        with self._lock:
            self._data["updated"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            tmp = self._path.with_suffix(".tmp")
            with open(tmp, "w") as f:
                json.dump(self._data, f, indent=2)
            os.replace(tmp, self._path)

    # -- bulk init --

    def set_accessions(self, accessions_by_strain):
        """Register accessions with their strains, preserving existing state."""
        for strain, accs in accessions_by_strain.items():
            for acc in accs:
                if acc not in self._data["accessions"]:
                    self._data["accessions"][acc] = {
                        "status": "pending",
                        "strain": strain,
                        "method": None,
                        "error": None,
                        "error_type": None,
                        "attempts": 0,
                    }
        self.save()

    # -- status updates --

    def mark_downloaded(self, acc, method, size_bytes=None):
        with self._lock:
            self._data["accessions"][acc].update(
                {
                    "status": "downloaded",
                    "method": method,
                    "size_bytes": size_bytes,
                    "error": None,
                    "error_type": None,
                }
            )
            self._data["accessions"][acc]["attempts"] += 1
        self.save()

    def mark_failed(self, acc, error_msg, error_type="transient"):
        with self._lock:
            self._data["accessions"][acc].update(
                {
                    "status": "failed",
                    "error": str(error_msg)[:200],
                    "error_type": error_type,
                }
            )
            self._data["accessions"][acc]["attempts"] += 1
        self.save()

    def mark_invalid(self, acc, reason):
        with self._lock:
            self._data["accessions"][acc].update(
                {
                    "status": "invalid",
                    "error": str(reason)[:200],
                    "error_type": "permanent",
                }
            )
            self._data["accessions"][acc]["attempts"] += 1
        self.save()

    def reset_to_pending(self, acc):
        with self._lock:
            self._data["accessions"][acc]["status"] = "pending"
            self._data["accessions"][acc]["error"] = None
            self._data["accessions"][acc]["error_type"] = None
        self.save()

    # -- queries --

    def get_pending(self):
        """Accessions that need downloading: pending or transient failures."""
        return [
            acc
            for acc, info in self._data["accessions"].items()
            if info["status"] == "pending" or (info["status"] == "failed" and info.get("error_type") != "permanent")
        ]

    def get_downloaded(self):
        return [acc for acc, info in self._data["accessions"].items() if info["status"] == "downloaded"]

    def get_strain(self, acc):
        return self._data["accessions"].get(acc, {}).get("strain")

    def summary(self):
        counts = {"pending": 0, "downloaded": 0, "failed": 0, "invalid": 0}
        for info in self._data["accessions"].values():
            counts[info["status"]] = counts.get(info["status"], 0) + 1
        return counts


# ---------------------------------------------------------------------------
# Rate limiter — thread-safe token bucket for Entrez API
# ---------------------------------------------------------------------------


class RateLimiter:
    """Simple thread-safe rate limiter."""

    def __init__(self, requests_per_second=3.0):
        self._lock = threading.Lock()
        self._min_interval = 1.0 / requests_per_second
        self._last_request = 0.0

    def acquire(self):
        with self._lock:
            now = time.monotonic()
            wait = self._min_interval - (now - self._last_request)
            if wait > 0:
                time.sleep(wait)
            self._last_request = time.monotonic()


# ---------------------------------------------------------------------------
# ZIP validation
# ---------------------------------------------------------------------------


def validate_genome_zip(zip_path, min_sequence_length=25000):
    """
    Validate a downloaded genome ZIP.

    Returns (is_valid, reason_if_invalid).
    SARS-CoV-2 is ~29,903 bp; 25,000 catches truncated downloads.
    """
    zip_path = Path(zip_path)

    if not zip_path.exists():
        return False, "file does not exist"
    if zip_path.stat().st_size == 0:
        return False, "file is empty"

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            # Find a FASTA file
            fasta_name = None
            for name in zf.namelist():
                if name.endswith((".fna", ".fa")):
                    fasta_name = name
                    break
            if not fasta_name:
                return False, "no FASTA file in ZIP"

            with zf.open(fasta_name) as f:
                text = f.read().decode("utf-8", errors="replace")

            # Parse first sequence
            seq_lines = []
            found_header = False
            for line in text.split("\n"):
                line = line.strip()
                if line.startswith(">"):
                    if found_header:
                        break  # only need first record
                    found_header = True
                    continue
                if found_header and line:
                    seq_lines.append(line)

            if not found_header:
                return False, "no FASTA header found"

            seq_len = sum(len(s) for s in seq_lines)
            if seq_len < min_sequence_length:
                return False, f"sequence too short ({seq_len} bp, need {min_sequence_length})"

    except zipfile.BadZipFile:
        return False, "corrupt ZIP file"
    except Exception as e:
        return False, f"validation error: {e}"

    return True, None


# ---------------------------------------------------------------------------
# FASTA parsing (no BioPython needed)
# ---------------------------------------------------------------------------


def _parse_fasta(text):
    """Parse FASTA text into list of {id, description, seq} dicts."""
    records = []
    header = None
    seq_parts = []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append(
                    {
                        "id": header.split()[0],
                        "description": header,
                        "seq": "".join(seq_parts),
                    }
                )
            header = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line)
    if header is not None:
        records.append(
            {
                "id": header.split()[0],
                "description": header,
                "seq": "".join(seq_parts),
            }
        )
    return records


# ---------------------------------------------------------------------------
# CLI transport — bulk download via NCBI datasets CLI
# ---------------------------------------------------------------------------


def _check_cli_available():
    return shutil.which("datasets") is not None


def _repackage_bulk_zip(bulk_zip_path, expected_accessions, output_path):
    """
    Extract multi-record genomic.fna from a CLI bulk ZIP into per-accession ZIPs.

    Returns (successful, failed) lists.
    """
    successful = []
    found = set()

    try:
        with zipfile.ZipFile(bulk_zip_path, "r") as zf:
            fna_path = None
            for name in zf.namelist():
                if name.endswith("genomic.fna"):
                    fna_path = name
                    break
            if not fna_path:
                return [], list(expected_accessions)

            with zf.open(fna_path) as fna_file:
                text = io.TextIOWrapper(fna_file).read()

            if SeqIO is not None:
                # BioPython path
                from io import StringIO

                records = [
                    {"id": r.id, "description": r.description, "seq": str(r.seq)}
                    for r in SeqIO.parse(StringIO(text), "fasta")
                ]
            else:
                records = _parse_fasta(text)

            for rec in records:
                acc = rec["id"]
                found.add(acc)
                per_zip = output_path / f"{acc}.zip"
                if per_zip.exists():
                    successful.append(acc)
                    continue

                fasta_content = f">{rec['description']}\n{rec['seq']}\n"
                buf = io.BytesIO()
                with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as out_zip:
                    out_zip.writestr("ncbi_dataset/data/genomic.fna", fasta_content)
                buf.seek(0)
                with open(per_zip, "wb") as f:
                    f.write(buf.read())
                successful.append(acc)

    except Exception as e:
        print(f"  [!] Error processing bulk ZIP: {e}")

    failed = [acc for acc in expected_accessions if acc not in found]
    return successful, failed


def download_via_cli(accessions, output_path, manifest, batch_size=50):
    """
    Download genomes via NCBI datasets CLI with adaptive batch sizing.
    On batch failure, bisects and retries each half once.
    """
    if not accessions:
        return

    def _run_batch(batch, label=""):
        acc_file = output_path / f"_batch{label}.txt"
        bulk_zip = output_path / f"_batch{label}.zip"

        try:
            with open(acc_file, "w") as f:
                f.write("\n".join(batch) + "\n")

            cmd = [
                "datasets",
                "download",
                "virus",
                "genome",
                "accession",
                "--inputfile",
                str(acc_file),
                "--include",
                "genome",
                "--filename",
                str(bulk_zip),
                "--no-progressbar",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode != 0:
                return None, result.stderr.strip()
            return bulk_zip, None

        except subprocess.TimeoutExpired:
            return None, "timeout"
        except Exception as e:
            return None, str(e)
        finally:
            if acc_file.exists():
                acc_file.unlink()

    total = len(accessions)
    for i in range(0, total, batch_size):
        batch = accessions[i : i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size
        print(f"  [CLI] Batch {batch_num}/{total_batches}: {len(batch)} genomes...")

        zip_path, err = _run_batch(batch, label=f"_{batch_num}")

        if err:
            # Bisect: try each half once before falling through to Entrez
            if len(batch) > 1:
                mid = len(batch) // 2
                for half_label, half in [("a", batch[:mid]), ("b", batch[mid:])]:
                    hz, herr = _run_batch(half, label=f"_{batch_num}{half_label}")
                    if herr:
                        print(f"  [!] Half-batch failed: {herr}")
                        for acc in half:
                            manifest.mark_failed(acc, herr, "transient")
                    else:
                        _process_bulk_zip(hz, half, output_path, manifest)
                        if hz.exists():
                            hz.unlink()
            else:
                manifest.mark_failed(batch[0], err, "transient")
            continue

        _process_bulk_zip(zip_path, batch, output_path, manifest)
        if zip_path.exists():
            zip_path.unlink()


def _process_bulk_zip(zip_path, batch, output_path, manifest):
    """Repackage a bulk ZIP and validate each resulting per-accession ZIP."""
    successful, failed = _repackage_bulk_zip(zip_path, batch, output_path)

    for acc in successful:
        per_zip = output_path / f"{acc}.zip"
        valid, reason = validate_genome_zip(per_zip)
        if valid:
            manifest.mark_downloaded(acc, "cli", per_zip.stat().st_size)
        else:
            print(f"  [!] {acc}: invalid after repackage ({reason})")
            per_zip.unlink(missing_ok=True)
            manifest.mark_failed(acc, reason, "transient")

    for acc in failed:
        manifest.mark_failed(acc, "not found in bulk ZIP", "transient")


# ---------------------------------------------------------------------------
# Entrez efetch transport — stable fallback
# ---------------------------------------------------------------------------


def _classify_http_error(e):
    """Classify an HTTP error as transient or permanent."""
    if isinstance(e, urllib.error.HTTPError):
        if e.code in (429, 500, 502, 503, 504):
            return "transient"
        if e.code == 404:
            return "permanent"
    if isinstance(e, (urllib.error.URLError, TimeoutError, ConnectionError)):
        return "transient"
    return "unknown"


def _backoff_sleep(attempt, base=1.0, cap=30.0):
    delay = min(cap, base * (2 ** (attempt - 1)))
    jitter = random.uniform(0, 0.25 * delay)
    time.sleep(delay + jitter)


def download_single_entrez(acc, output_path, rate_limiter, manifest):
    """Download one genome via Entrez efetch and wrap in a per-accession ZIP."""
    zip_path = output_path / f"{acc}.zip"
    max_attempts = 3

    api_key = os.environ.get("NCBI_API_KEY", "")
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = f"db=nucleotide&id={acc}&rettype=fasta&retmode=text"
    if api_key:
        params += f"&api_key={api_key}"
    url = f"{base_url}?{params}"

    for attempt in range(1, max_attempts + 1):
        if attempt > 1:
            _backoff_sleep(attempt)

        rate_limiter.acquire()

        try:
            req = urllib.request.Request(url, headers={"Accept": "text/plain"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                body = resp.read().decode("utf-8", errors="replace")

            # Check for NCBI error responses
            if not body.strip() or not body.strip().startswith(">"):
                if "Nothing has been found" in body or "not found" in body.lower():
                    manifest.mark_invalid(acc, "accession not found on NCBI")
                    return
                if attempt == max_attempts:
                    manifest.mark_failed(acc, "empty or non-FASTA response", "transient")
                    return
                continue

            # Parse and validate sequence length
            records = _parse_fasta(body)
            if not records or len(records[0]["seq"]) < 25000:
                if not records:
                    reason = "no FASTA records in response"
                else:
                    reason = f"sequence too short ({len(records[0]['seq'])} bp)"
                if attempt == max_attempts:
                    manifest.mark_failed(acc, reason, "transient")
                    return
                continue

            # Wrap as ZIP matching training pipeline's expected format
            fasta_content = f">{records[0]['description']}\n{records[0]['seq']}\n"
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.writestr("ncbi_dataset/data/genomic.fna", fasta_content)
            buf.seek(0)
            with open(zip_path, "wb") as f:
                f.write(buf.read())

            manifest.mark_downloaded(acc, "entrez", zip_path.stat().st_size)
            return

        except urllib.error.HTTPError as e:
            err_type = _classify_http_error(e)
            if err_type == "permanent":
                manifest.mark_invalid(acc, f"HTTP {e.code}: {e.reason}")
                return
            if attempt == max_attempts:
                manifest.mark_failed(acc, f"HTTP {e.code}: {e.reason}", err_type)
                return

        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if attempt == max_attempts:
                manifest.mark_failed(acc, str(e), "transient")
                return

        except Exception as e:
            if attempt == max_attempts:
                manifest.mark_failed(acc, str(e), "unknown")
                return


def download_via_entrez(accessions, output_path, manifest, max_workers=3):
    """Download remaining accessions via Entrez efetch with rate limiting."""
    if not accessions:
        return

    import concurrent.futures

    rate = 10.0 if os.environ.get("NCBI_API_KEY") else 3.0
    limiter = RateLimiter(rate)
    print(f"  [Entrez] Downloading {len(accessions)} genomes ({max_workers} workers, {rate:.0f} req/s)...")

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(download_single_entrez, acc, output_path, limiter, manifest): acc for acc in accessions}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            future.result()  # propagate any uncaught exception
            done += 1
            if done % 25 == 0 or done == len(accessions):
                print(f"  [Entrez] Progress: {done}/{len(accessions)}")


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_download(accessions_by_strain, output_dir, batch_size=50, max_entrez_workers=3):
    """
    Download genomes: CLI primary, Entrez fallback, persistent state.

    Args:
        accessions_by_strain: {"Reference": ["NC_045512.2", ...], ...}
        output_dir: Path to write per-accession ZIP files
        batch_size: CLI batch size (default 50)
        max_entrez_workers: Parallel workers for Entrez fallback (default 3)

    Returns:
        str: Path to output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load or create persistent manifest
    manifest = DownloadManifest.load(output_dir)
    manifest.set_accessions(accessions_by_strain)

    total = sum(len(v) for v in accessions_by_strain.values())
    print(f"[*] Download manifest: {total} accessions across {len(accessions_by_strain)} strains")

    # Phase 0: validate existing ZIPs — delete corrupt ones
    already_downloaded = manifest.get_downloaded()
    revalidated = 0
    for acc in already_downloaded:
        zip_path = output_path / f"{acc}.zip"
        if not zip_path.exists():
            manifest.reset_to_pending(acc)
            revalidated += 1
            continue
        valid, reason = validate_genome_zip(zip_path)
        if not valid:
            print(f"  [!] Removing corrupt ZIP: {acc} ({reason})")
            zip_path.unlink(missing_ok=True)
            manifest.reset_to_pending(acc)
            revalidated += 1
    if revalidated:
        print(f"[*] Reset {revalidated} corrupt/missing files for re-download")

    # Check what still needs downloading
    pending = manifest.get_pending()
    counts = manifest.summary()
    print(f"[*] Status: {counts['downloaded']} downloaded, {counts['invalid']} invalid, {len(pending)} pending")

    if not pending:
        print("[OK] All accessions already downloaded")
        _print_summary(manifest)
        return str(output_path)

    # Phase 1: CLI bulk download
    if _check_cli_available():
        print(f"\n[*] Phase 1: CLI bulk download ({len(pending)} accessions, batch_size={batch_size})")
        download_via_cli(pending, output_path, manifest, batch_size)
    else:
        print("[!] NCBI datasets CLI not found — skipping to Entrez fallback")

    # Phase 2: Entrez fallback for anything still pending
    still_pending = manifest.get_pending()
    if still_pending:
        print(f"\n[*] Phase 2: Entrez efetch fallback ({len(still_pending)} accessions)")
        download_via_entrez(still_pending, output_path, manifest, max_entrez_workers)

    _print_summary(manifest)
    return str(output_path)


def _print_summary(manifest):
    counts = manifest.summary()
    total = sum(counts.values())
    print(f"\n{'=' * 50}")
    print("Download Summary")
    print(f"{'=' * 50}")
    print(f"  Downloaded : {counts['downloaded']}/{total}")
    print(f"  Failed     : {counts['failed']}")
    print(f"  Invalid    : {counts['invalid']} (permanently unavailable)")
    print(f"  Pending    : {counts['pending']}")
    print(f"{'=' * 50}")
