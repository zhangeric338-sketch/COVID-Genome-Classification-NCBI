# CLAUDE.md

Guidance for future coding agents working in this repository.

## Project Overview
- Purpose: classify COVID-related viral genomes from NCBI data.
- Main entry point: `main.py`.
- Core package code: `src/`.
- Utility scripts: `tools/`.

## Directory Layout
- `src/`: model, preprocessing, and pipeline modules.
- `tools/`: data and maintenance utilities.
- `data/`: downloaded and generated datasets/artifacts.

## Run Data Locations
- Base run data directory is controlled by `--output-dir` (default: `data`).
- Downloaded genomes are stored as `*.zip` files directly under `<output-dir>/`.
- Strain manifest is typically at `<output-dir>/strain_manifest.json` (fallbacks may use `tools/accessions.json`).
- Train/test split output is written to:
  - `<output-dir>/train/*.zip`
  - `<output-dir>/test/*.zip`
- Accessions source file defaults to `tools/accessions.json` when present or when fetched.
- Visualization images from `main.py` are saved in the repository root (current working directory), including:
  - `download_summary.png`
  - `train_test_split.png`
  - `dataset_composition.png`
  - `<column>_plot.png` (for metadata plots such as collection date/host/location)
- Training consumes `<output-dir>/train` and `<output-dir>/test`; model metrics/provenance are printed and optionally logged to Weights & Biases.

## Required Workflow
1. Read `README.md` and `main.py` before making substantial changes.
2. Prefer editing existing files over creating new files.
3. Keep changes scoped and minimal.
4. Run targeted validation after edits.

## Data Download Safety
- Default to partial dataset only (`--size-gb 0.05`).
- Never use full dataset unless explicitly requested with `--full-dataset`.
- Use deterministic seed `42` unless told otherwise.
- Check for existing files before downloading.
- Warn before large downloads (>1 GB).

## Coding Conventions
- Follow existing Python style and naming.
- Keep logic readable and avoid overengineering.
- Add concise comments only where logic is non-obvious.
- Avoid broad refactors unless requested.

## Validation Checklist
- For `main.py` flow changes, run: `python main.py --help` and relevant mode checks.
- For download-path edits, verify with a small partial run.
- Confirm no empty or obviously invalid output artifacts are produced.

## Agent Notes
- Do not remove user changes unrelated to the task.
- If repository state is dirty, isolate edits to requested scope.
- Keep dependency changes explicit and minimal.
