"""Run checksum controls without sharing another package's audit module."""
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT/'research/checksum-variance/test_package.py'


class ChecksumVarianceTests(unittest.TestCase):
    def test_package_in_an_isolated_interpreter(self):
        subprocess.run([sys.executable, '-B', str(PATH)], cwd=ROOT,
                       check=True, timeout=120)
