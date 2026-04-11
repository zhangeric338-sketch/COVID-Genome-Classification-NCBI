.PHONY: setup lint format check data train clean help

PYTHON  ?= python3.12
VENV    := .venv
PIP     := $(VENV)/bin/pip
PY      := $(VENV)/bin/python
RUFF    := $(VENV)/bin/ruff

SEED    ?= 42
WORKERS ?= 4
MODEL   ?= random_forest
KMER    ?= 6

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?##' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

setup: $(VENV)/bin/activate ## Create venv and install dependencies

$(VENV)/bin/activate: requirements.txt
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip setuptools
	$(PIP) install -r requirements.txt
	touch $(VENV)/bin/activate

lint: ## Run ruff linter
	$(RUFF) check .

format: ## Run ruff formatter
	$(RUFF) format .

check: ## Run lint + format check (CI-friendly, no edits)
	$(RUFF) check .
	$(RUFF) format --check .

data: ## Download balanced dataset via Entrez pipeline
	$(PY) main.py --fetch-accessions --seed $(SEED) --workers $(WORKERS)

train: ## Run training (MODEL=random_forest|mlp, KMER=6)
	$(PY) main.py --fetch-accessions --train --model $(MODEL) --kmer-size $(KMER) --seed $(SEED) --workers $(WORKERS)

clean: ## Remove generated data, caches, and plots
	rm -rf data/ __pycache__ src/__pycache__ tools/__pycache__
	rm -f tools/accessions.json
	rm -f *.png
