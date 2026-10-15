"""Regression checks for stale evidence, selected-result claims, and package discovery."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from check_project import CORE, check_project, check_receipt, check_selection, read_json


class ProjectTests(unittest.TestCase):
    def setUp(self):
        self.selected = read_json(ROOT/'certificates/selected-result.json')
        self.witness_path = ROOT/self.selected['witnesses']
        self.receipt = read_json(ROOT/self.selected['certificate'])

    def test_project_evidence_and_links_are_consistent(self):
        result = check_project(ROOT)
        self.assertEqual(result['status'], 'PASS')
        self.assertGreaterEqual(result['packages'], 1)
        self.assertGreaterEqual(result['witnesses'], 4)
        self.assertFalse(result['replayed'])

    def test_selected_exponent_cannot_drift_from_its_witness(self):
        selected = copy.deepcopy(self.selected)
        selected['alpha'] = '0.489'
        with self.assertRaisesRegex(ValueError, 'selected exponent'):
            check_selection(ROOT, selected)

    def test_stale_witness_and_checker_hashes_fail(self):
        for field in ['witnesses_sha256', 'package_sha256']:
            receipt = copy.deepcopy(self.receipt)
            if field == 'witnesses_sha256':
                receipt[field] = '0'*64
            else:
                receipt[field]['certify.py'] = '0'*64
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'stale'):
                check_receipt(ROOT, self.witness_path, receipt)

    def test_incomplete_receipt_coverage_is_rejected(self):
        for mutate in [lambda r: r['results'].pop(),
                       lambda r: r['package_sha256'].pop('SOURCE.json')]:
            receipt = copy.deepcopy(self.receipt)
            mutate(receipt)
            with self.assertRaises(ValueError):
                check_receipt(ROOT, self.witness_path, receipt)

    def test_replay_configuration_is_recovered_and_validated(self):
        self.assertEqual(check_receipt(ROOT, self.witness_path, self.receipt), (128, 50))
        receipt = copy.deepcopy(self.receipt)
        receipt['results'][0]['moments']['grid'] = 63
        with self.assertRaisesRegex(ValueError, 'replay precision'):
            check_receipt(ROOT, self.witness_path, receipt)

    def test_new_packages_are_discovered_and_bound(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            shutil.copytree(ROOT/'research', root/'research')
            shutil.copytree(ROOT/'certificates', root/'certificates')
            shutil.copytree(ROOT/'docs', root/'docs')
            for path in ROOT.glob('*.md'):
                shutil.copyfile(path, root/path.name)
            for name in ['LICENSE', 'NOTICE', 'CITATION.cff']:
                shutil.copyfile(ROOT/name, root/name)
            before = check_project(root)
            addition = Path(tempfile.mkdtemp(prefix='new-package-', dir=root/'research'))
            shutil.copytree(root/CORE, addition, dirs_exist_ok=True)
            result = check_project(root)
            added_witnesses = len(read_json(addition/'witnesses.json')['witnesses'])
            self.assertEqual(result['packages'], before['packages']+1)
            self.assertEqual(result['witnesses'], before['witnesses']+added_witnesses)
            path = addition/'witnesses.json'
            data = read_json(path)
            data['witnesses'][0]['parameters']['alpha'] = '0.4'
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, 'stale witness'):
                check_project(root)


if __name__ == '__main__':
    unittest.main()
