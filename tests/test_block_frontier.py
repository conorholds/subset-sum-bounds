"""Include the portable block-frontier package controls in make check."""
import importlib.util
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]/'research/block-count-frontier/test_package.py'
SPEC = importlib.util.spec_from_file_location('block_frontier_controls', PACKAGE)
controls = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(controls)
BlockBoundTests = controls.BlockBoundTests
FrontierTests = controls.FrontierTests
