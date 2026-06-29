"""
Sequence degradation transforms for robustness experiments.

The strain classifier is trained on clean, complete ~29,900 bp genomes and
scores ~99.5% on equally clean test genomes. Real clinical samples are rarely
that clean: sequencing produces partial genomes, ambiguous bases (N), base-call
errors, and sparse read coverage. These transforms simulate each failure mode
so we can measure how far a genome can degrade before classification breaks.

Every transform is deterministic given a numpy RandomState, so a sweep is
reproducible. None of these mutate the input string.
"""

from __future__ import annotations

import numpy as np

_BASES = ("A", "C", "G", "T")


def truncate_sequence(sequence: str, keep_fraction: float, rng: np.random.RandomState) -> str:
    """
    Keep a contiguous window covering `keep_fraction` of the genome, starting
    at a random offset. Simulates incomplete sequencing where only part of the
    genome is recovered.
    """
    if keep_fraction >= 1.0:
        return sequence
    keep_fraction = max(keep_fraction, 0.0)
    keep_len = int(len(sequence) * keep_fraction)
    if keep_len <= 0:
        return ""
    max_start = len(sequence) - keep_len
    start = int(rng.randint(0, max_start + 1)) if max_start > 0 else 0
    return sequence[start : start + keep_len]


def inject_ns(sequence: str, n_fraction: float, rng: np.random.RandomState) -> str:
    """
    Replace a random `n_fraction` of positions with 'N'. Simulates low-coverage
    or ambiguous base calls. Because k-mers containing N are dropped downstream,
    this thins the usable feature signal rather than corrupting it.
    """
    if n_fraction <= 0.0:
        return sequence
    n_fraction = min(n_fraction, 1.0)
    seq = np.frombuffer(sequence.encode("ascii"), dtype="S1").copy()
    mask = rng.random_sample(len(seq)) < n_fraction
    seq[mask] = b"N"
    return seq.tobytes().decode("ascii")


def add_substitutions(sequence: str, sub_rate: float, rng: np.random.RandomState) -> str:
    """
    Mutate a random `sub_rate` of positions to a different base (A/C/G/T).
    Simulates sequencing base-call errors. Unlike N-injection, these positions
    stay valid bases, so they actively corrupt k-mers rather than removing them.
    Existing non-ACGT characters (e.g. N) are left untouched.
    """
    if sub_rate <= 0.0:
        return sequence
    sub_rate = min(sub_rate, 1.0)
    seq = list(sequence)
    n = len(seq)
    n_sub = int(n * sub_rate)
    if n_sub <= 0:
        return sequence
    positions = rng.choice(n, size=n_sub, replace=False)
    for pos in positions:
        original = seq[pos]
        if original not in _BASES:
            continue
        alt = _BASES[rng.randint(0, 4)]
        while alt == original:
            alt = _BASES[rng.randint(0, 4)]
        seq[pos] = alt
    return "".join(seq)


def fragment_sequence(
    sequence: str,
    keep_fraction: float,
    rng: np.random.RandomState,
    fragment_length: int = 300,
) -> str:
    """
    Split the genome into non-overlapping `fragment_length` fragments (like
    short sequencing reads), keep a random `keep_fraction` of them, and
    concatenate the survivors in original order. Simulates sparse read coverage
    where only some regions of the genome are observed.
    """
    if keep_fraction >= 1.0:
        return sequence
    keep_fraction = max(keep_fraction, 0.0)
    fragments = [sequence[i : i + fragment_length] for i in range(0, len(sequence), fragment_length)]
    if not fragments:
        return ""
    keep_mask = rng.random_sample(len(fragments)) < keep_fraction
    kept = [frag for frag, keep in zip(fragments, keep_mask) if keep]
    return "".join(kept)


# Registry consumed by the experiment runner. Each entry maps a degradation
# name to (transform, severity_levels, x_axis_label). Severity grids run from
# "clean" to "severe" so the resulting curves are directly comparable.
DEGRADATIONS = {
    "truncation": {
        "fn": truncate_sequence,
        "param": "keep_fraction",
        "levels": [1.0, 0.75, 0.5, 0.25, 0.1, 0.05],
        "xlabel": "fraction of genome kept",
        # smaller keep_fraction = more degraded; plot as fraction removed
        "severity": lambda v: 1.0 - v,
    },
    "n_injection": {
        "fn": inject_ns,
        "param": "n_fraction",
        "levels": [0.0, 0.05, 0.1, 0.2, 0.4, 0.6],
        "xlabel": "fraction of bases set to N",
        "severity": lambda v: v,
    },
    "substitution": {
        "fn": add_substitutions,
        "param": "sub_rate",
        "levels": [0.0, 0.005, 0.01, 0.02, 0.05, 0.1],
        "xlabel": "base substitution rate",
        "severity": lambda v: v,
    },
    "fragmentation": {
        "fn": fragment_sequence,
        "param": "keep_fraction",
        "levels": [1.0, 0.5, 0.25, 0.1, 0.05],
        "xlabel": "fraction of 300bp fragments kept",
        "severity": lambda v: 1.0 - v,
    },
}
