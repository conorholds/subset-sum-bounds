PYTHON ?= python3

.PHONY: check test verify search

check: test verify

test:
	$(PYTHON) -B research/filter-budgets/test_package.py
	$(PYTHON) -B -m unittest discover -s tests -v

verify:
	$(PYTHON) -B scripts/check_project.py --replay

search:
	@mkdir -p .context
	$(PYTHON) -B research/filter-budgets/search.py --start combined --parameters combined --degree 14 --steps 4000 --output .context/candidate.json
