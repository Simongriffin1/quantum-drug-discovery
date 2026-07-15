# PeptideForge monorepo Makefile
# Reference CURSOR_PROJECT_CONTEXT.md for non-negotiables.

.PHONY: help install install-core install-backend install-frontend install-benchmarks schemas test test-core test-backend test-benchmarks lint lint-core typecheck typecheck-core docker-up docker-down docker-ps oracle-validity benchmark-report ci

PYTHON ?= python3.11
POETRY ?= poetry
COMPOSE := docker compose -f infra/docker-compose.yml

help:
	@echo "PeptideForge targets:"
	@echo "  make install            - install core + benchmarks + backend + frontend"
	@echo "  make test               - run Python + frontend unit/contract tests"
	@echo "  make oracle-validity    - run OpenMM MM-GBSA Spearman + MLflow log"
	@echo "  make benchmark-report   - regenerate P12 credibility markdown report"
	@echo "  make ci                 - lint + typecheck + test (local CI parity)"

install: install-benchmarks install-core install-backend install-frontend

install-benchmarks:
	cd benchmarks && $(POETRY) env use $(PYTHON) && $(POETRY) install

install-core:
	cd core && $(POETRY) env use $(PYTHON) && $(POETRY) install
	cd core && $(POETRY) run pip install --only-binary=:all: 'openmm>=8.1' 'mlflow>=2.13' 'pyscf>=2.4' 'pennylane>=0.36' pdbfixer || true

install-backend:
	cd backend && $(POETRY) env use $(PYTHON) && $(POETRY) install

install-frontend:
	cd frontend && npm install

schemas:
	cd core && $(POETRY) run python -m peptideforge.contracts.export_schemas

test: test-benchmarks test-core test-backend
	cd frontend && npm test

test-benchmarks:
	cd benchmarks && $(POETRY) run pytest

test-core:
	cd core && $(POETRY) run pytest -m "not slow"

test-backend:
	cd backend && $(POETRY) run pytest

lint: lint-core
	cd frontend && npm run lint

lint-core:
	cd core && $(POETRY) run ruff check src tests && $(POETRY) run black --check src tests
	cd benchmarks && $(POETRY) run ruff check peptideforge_benchmarks tests && $(POETRY) run black --check peptideforge_benchmarks tests

typecheck: typecheck-core
	cd frontend && npm run typecheck

typecheck-core:
	cd core && $(POETRY) run mypy src
	cd benchmarks && $(POETRY) run mypy peptideforge_benchmarks

docker-up:
	$(COMPOSE) up -d

docker-down:
	$(COMPOSE) down

docker-ps:
	$(COMPOSE) ps

oracle-validity:
	cd core && $(POETRY) run python -m peptideforge.oracles.run_oracle_validity --mlflow

benchmark-report:
	cd core && $(POETRY) run python -m peptideforge.reports.generate_benchmark_report --also-json \
		--mlflow-uri sqlite:///$$(cd .. && pwd)/benchmarks/fixtures/structures/oracle_validity_mlflow.db

ci: lint typecheck test
