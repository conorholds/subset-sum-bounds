"""Run the explicit tunable-multiplier model and its failure controls."""
import importlib.util
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / 'research/chernoff-multipliers/test_package.py'
SPEC = importlib.util.spec_from_file_location('chernoff_multiplier_controls', PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def load_tests(loader, tests, pattern):
    return loader.loadTestsFromModule(MODULE)
