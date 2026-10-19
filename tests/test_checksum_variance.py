"""Route the isolated checksum audit's exact tests through make check."""
import importlib.util
from pathlib import Path
import sys

PACKAGE = Path(__file__).resolve().parents[1]/'research/checksum-variance'


def load_tests(loader, tests, pattern):
    sys.path.insert(0, str(PACKAGE))
    spec = importlib.util.spec_from_file_location('checksum_variance_tests', PACKAGE/'test_package.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return loader.loadTestsFromModule(module)
