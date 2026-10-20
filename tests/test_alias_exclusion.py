"""Include the finite alias-exclusion audit in the project test router."""
import importlib.util
from pathlib import Path

PATH = Path(__file__).resolve().parents[1]/'research/alias-exclusion/test_package.py'
SPEC = importlib.util.spec_from_file_location('alias_exclusion_controls', PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def load_tests(loader, tests, pattern):
    return loader.loadTestsFromModule(MODULE)
