# Design Doc — Degradation Robustness of SARS-CoV-2 Strain Classifiers

**Status:** Draft v1 · **Owner:** Eric Zhang · **Last updated:** 2026-06-28
**Branch:** `feature/degradation-robustness`
**Target venue:** Journal of Emerging Investigators (JEI), or arXiv preprint

---

## 0. TL;DR

A SARS-CoV-2 strain classifier hits ~99.5% clean accuracy — a trivial result, because the
strains are genetically far apart. This study asks the question that the clean number hides:
**what happens when the input genome is corrupted, and does the answer depend on how we
represent the genome?**

We subject a Random Forest classifier to four controlled, parametric corruptions
(truncation, ambiguous-base injection, substitution noise, fragmentation) and measure how
macro-F1 degrades, under **two feature representations**:

1. **k-mer frequency features** (swept across multiple k), and
2. **mutation / substitution-site features** (derived from Nextclade).

The core hypothesis is that these two representations **fail differently** — k-mers are
fragile to point corruption (one substitution damages *k* overlapping k-mers) while mutation
features watch a small set of informative sites and should be comparatively robust to
substitution but fragile to truncation. We confirm the breaking points, isolate the
mechanisms, and add a single realistic-sequencing-error condition to show the controlled
findings are not an artifact of crude uniform noise.

---

## 1. Motivation & Problem Statement

### 1.1 Why the clean result is not interesting
The 7 SARS-CoV-2 strains in this dataset (Alpha, Beta, Delta, Gamma, Omicron, Recent,
Reference) are separated by many lineage-defining mutations. Distinguishing them from clean,
complete genomes is close to a solved problem — ~99.5% accuracy reflects the *task* being
easy, not the *model* being good. A paper that reports only this number has no finding.

### 1.2 The real question
Real genomic data is rarely clean: assemblies are truncated, basecalls produce ambiguous
`N`s, sequencing introduces substitution errors, and low-coverage assemblies fragment. A
useful, publishable question is therefore:

> How gracefully does strain classification degrade as the genome is corrupted, where are the
> breaking points, and does the **feature representation** change the answer?

This reframes a trivial classification task into a study of **robustness and representation**,
which is both more honest and more interesting.

---

## 2. Goals & Success Criteria

### 2.1 Goals
- **G1.** Quantify how macro-F1 degrades for each of 4 corruption types across a severity sweep.
- **G2.** Compare two feature representations (k-mer vs. mutation) and show *whether and how*
  they fail differently — including identifying each one's breaking point.
- **G3.** Explain the *mechanism* behind the difference (k-mer fragility vs. site sparsity),
  not just report curves.
- **G4.** Show robustness findings survive under one **realistic** sequencing-error simulator,
  not only under hand-rolled uniform noise.
- **G5.** A small **robustness-transfer** probe: does training on some corruptions help on a
  *held-out* corruption type?

### 2.2 Definition of success (falsifiable)
A single, self-contained paper that:
- reports degradation curves for 4 corruptions × 2 representations × a k-sweep,
- states each representation's breaking point with a clear mechanistic explanation, and
- is clean enough to submit to JEI without an obvious "you didn't control for X" rejection.

Success is **not** a leaderboard number; it is a defensible finding the author can fully
explain and own.

### 2.3 Primary metric
**Macro-averaged F1.** The dataset is imbalanced (Beta is the minority class, ~23 samples),
so raw accuracy would be dominated by majority classes and would hide collapse on rare
strains. Macro-F1 weights every strain equally. Accuracy is reported as a secondary metric.

---

## 3. Background & Prior Art

### 3.1 The closest prior work
**"Benchmarking machine learning robustness in Covid-19 genome sequence classification"**
(Scientific Reports, 2023; arXiv:2207.08898; PMC10010240) already benchmarks SARS-CoV-2
lineage-classifier robustness across three representations (k-mer, PSSM, minimizer) under
realistic sequencing-error simulators.

**Honest consequence:** the *headline* idea — "test classifier robustness to genome
corruption" — is **not novel.** Any design doc that ignores this gets desk-rejected. We cite
it directly and position against it.

### 3.2 Our surviving delta (what is actually new here)
1. **Mutation / substitution-site features** as a representation. The 2023 benchmark does not
   include a biology-informed mutation-site representation. This is the strongest novelty:
   it ties the robustness question to *what virologists actually look at* (lineage-defining
   mutations), and the mechanistic prediction (robust to substitution, fragile to truncation)
   is specific and testable.
2. **Controlled, per-axis parametric corruption.** Rather than only realistic mixed-error
   simulators, we isolate one corruption mechanism at a time and sweep its severity, so we can
   attribute a breaking point to a *specific* mechanism.
3. **Robustness-transfer probe.** A small leave-one-corruption-out experiment asking whether
   robustness *generalizes* across corruption types — framed as a proof-of-concept, not a full
   benchmark.

### 3.3 Forced scope consequences from the prior art
- **Varying k is mandatory, not optional.** The literature reports k-mers as relatively robust
  at coarse taxonomy / small k (k≈3); our spike shows k=6 collapsing under 1–5% substitution.
  If we do not sweep k, "k-mers are fragile" reads as a k=6 artifact and a reviewer kills it.
- **A realistic-simulator condition is effectively required** (one condition) to preempt the
  "why only crude uniform noise?" objection.

---

## 4. Data

- **Source:** NCBI GenBank SARS-CoV-2 complete genomes (human host).
- **Size:** 986 genomes, split **791 train / 195 test**, balanced-by-strain split, seed=42.
- **Classes (7):** Alpha, Beta, Delta, Gamma, Omicron, Recent, Reference.
- **Imbalance:** Beta is the minority class (~23 samples) — the reason for macro-F1.
- **Provenance:** training records a SHA-256 hash over the sorted train/test accession lists,
  plus an identical-sequence leakage check across train/test.

> **NOTE (precondition):** The working dataset and the BioPython Windows loader fix currently
> live on `fix/download-issues`, **not** on `main`. See §9 (Preconditions).

---

## 5. Method

### 5.1 Classifier
Random Forest (100 trees, `max_depth=20`, `random_state=42`). Deliberately **not** deep
learning (see Non-Goals §8): RF on interpretable features is explainable, reproducible on a
laptop, and fully ownable by a solo author. The classifier is held fixed across all
conditions so that differences are attributable to representation and corruption, not to model
tuning.

### 5.2 Representations
| | k-mer frequencies | Mutation features |
|---|---|---|
| **Feature** | normalized counts of canonical k-mers | binary indicators over substitution sites (e.g. `C241T`) |
| **Source** | sequence directly | Nextclade `substitutions` column |
| **Dimension** | grows with k | ~hundreds of informative sites (min site count ≥ 2) |
| **k handling** | swept (see §5.4) | n/a |
| **Predicted weakness** | point corruption (substitution, N) | truncation / fragmentation (loses sites) |

Both vocabularies are **frozen on the clean training set** and reused at evaluation, so we
never leak test-time information into the feature space.

### 5.3 Corruptions (4 axes, seeded & reproducible)
Implemented in `src/degradation.py`, each a deterministic function of a `numpy.RandomState`:

1. **Truncation** — keep a leading fraction of the genome.
2. **N-injection** — replace a fraction of bases with `N` (ambiguous basecalls).
3. **Substitution noise** — flip a fraction of bases to a *different* base (sequencing error).
4. **Fragmentation** — keep a fraction of the genome as short fragments (low-coverage assembly).

Each corruption is swept across a severity grid. A fresh `RandomState(seed)` per
(corruption, level) keeps every point independent yet reproducible.

### 5.4 The experiment matrix
```
representations (2) × corruptions (4) × severities (~6) × k-values (~4: k∈{3,4,5,6,7})
                    × { stress-test, robustness-transfer(small) }
```
Mostly cheap compute. The cost is **narrative coherence**, not CPU — see §7.

### 5.5 Experiments
- **E1 — Stress test (the spine).** Train on clean genomes; evaluate on degraded test genomes
  across all corruptions × severities × k. Produces the core degradation curves and breaking
  points. Answers G1–G3.
- **E2 — Robustness-transfer probe (small).** One leave-one-corruption-out run: augment
  training with some corruption types, hold one out, test on the held-out type. Answers G5 as
  a proof-of-concept (not a full 4-fold matrix).
- **E3 — Realistic-simulator validation (one condition).** Re-run the key comparison under a
  realistic read simulator (e.g. InSilicoSeq / PBSIM-style error model) to show the controlled
  findings hold. Answers G4.

### 5.6 Mechanistic prediction (the thing we are testing)
A single substitution corrupts *k* overlapping k-mers, so k-mer features should degrade fast
under substitution and worse at larger k. Mutation features watch a small set of informative
sites scattered across ~30 kb, so random substitutions mostly miss them → comparative
robustness to substitution, but fragility to truncation/fragmentation (which delete sites
wholesale). E1–E3 are designed to confirm or refute this.

---

## 6. De-risking: the spike (already run)

Before committing to the full harness, a throwaway **spike**
(`tools/spike_representation_robustness.py`) tested the two load-bearing assumptions on one
corruption (substitution) at three severities, for both representations:

- **Assumption 1 (structure):** macro-F1 declines as a structured slope, not a cliff or a flat
  line. → **confirmed.**
- **Assumption 2 (divergence):** the two representations fail *differently*. → **confirmed**;
  k-mer (k=6) macro-F1 collapsed under 1–5% substitution while mutation features stayed near
  baseline.

> **TODO:** insert the exact spike F1 figures (clean / mild / severe for each representation)
> from a clean re-run once the foundation in §9 is resolved. Do **not** paraphrase numbers from
> memory in the paper.

Verdict at spike time: **green light**, with the forced addition of a k-sweep (§3.3).

---

## 7. Risks & Mitigations

| Risk | Severity | Mitigation |
|---|---|---|
| **Scope sprawl** — the matrix is large and the *narrative* can balloon into an unfinishable paper | High | Stress test is the spine; transfer is one small experiment; simulator is one condition. Build incrementally (§10). |
| **Novelty challenge** — 2023 benchmark already exists | High | Cite it; lead with mutation-feature representation + per-axis mechanism isolation as the delta. |
| **"k=6 is cherry-picked"** | High | k-sweep is mandatory, baked into E1. |
| **"Uniform noise is unrealistic"** | Medium | E3 realistic-simulator condition. |
| **Nextclade re-runs on every degraded set are slow** | Medium | Cache; run mutation eval only at the severities that matter; parallelize. |
| **Class imbalance distorts results** | Medium | Macro-F1 primary; report per-class F1 for Beta. |
| **Reproducibility** | Low | Fixed seeds, frozen vocab, dataset hash, leakage check already in `training.py`. |

---

## 8. Non-Goals (explicitly out of scope)

1. **No full Pango lineage classification.** The 215-class sublineage problem is unrealistic at
   this data scale. Task is strain-level (7 classes).
2. **No deep learning / CNNs / transformers.** RF on interpretable features is a deliberate
   scoping decision, not a limitation to apologize for.
3. **No new data collection.** We use the existing 791/195 split; no scraping or rebalancing by
   gathering more Beta samples.
4. **No chasing clean accuracy.** The ~99.5% clean number is the starting point, not the
   contribution.
5. **No production tool / deployment.** This is a controlled study of mechanism, not a hardened
   classifier for labs.

---

## 9. Preconditions (MUST resolve before the first code PR)

These are blocking for the **build**, not for this design:

- **P1 — Branch base.** `main` lacks (a) the BioPython `TextIOWrapper` fix (so genomes do not
  load on Windows) and (b) the `load_genome_dataset(..., k=...)` signature the spike calls.
  The working pipeline + dataset live on `fix/download-issues`. The robustness work must be
  rebased onto a base that actually loads data. Decision required:
  (i) reconcile `fix/download-issues` → `main` first, then rebase; or
  (ii) rebase this branch onto `fix/download-issues`.
- **P2 — Spike compatibility.** `tools/spike_representation_robustness.py` calls
  `load_genome_dataset("data", k=K)`, which only exists on `fix`. Must be reconciled with
  whatever base is chosen in P1.
- **P3 — Dataset in git.** `fix` commits ~1000 genome `.zip` files into version control while
  `.gitignore` ignores `data/`. Decide on a storage strategy (manifest + download script, or a
  Zenodo/release asset) before this repo becomes the paper's public artifact.

---

## 10. Implementation Plan (incremental — do NOT build the whole matrix at once)

1. **Resolve §9 preconditions** (branch base + spike compatibility).
2. **PR #1 — stress-test spine, single k.** Clean-train → degraded-eval for all 4 corruptions
   at one k, both representations. Produces first real curves. *(This is the first authorized
   build step.)*
3. **PR #2 — k-sweep.** Add k ∈ {3,4,5,6,7} to the stress test.
4. **PR #3 — robustness-transfer probe.** One leave-one-corruption-out run.
5. **PR #4 — realistic-simulator condition.**
6. **PR #5 — figures, tables, and writeup.**

Each PR is independently reviewable and leaves the repo in a working state.

---

## 11. Open Questions
- Which realistic simulator (InSilicoSeq vs. PBSIM3) best matches the assembly failure modes we
  care about, with the least integration cost?
- For mutation features under truncation: do we re-run Nextclade on the truncated sequence
  (realistic) or mask sites (faster)? Realistic is correct but slower — confirm in PR #3.
- Final k grid: is k=7 worth the dimensionality blow-up, or stop at k=6?
