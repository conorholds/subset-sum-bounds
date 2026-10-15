"""Source controls, independent scalar quadrature, and rejection tests."""
import copy
from decimal import Decimal
import json
from pathlib import Path
import tempfile
import unittest

from mpmath import iv

from certify import budget_margins, certify_moments, enclosure, q
from filters import direct_nodes, moment_logs
from search import search
from verify import load_witnesses


class BudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.witnesses = {w['id']: w for w in load_witnesses()}
        iv.dps = 50

    def test_all_retained_budget_margins_positive(self):
        for name, witness in self.witnesses.items():
            with self.subTest(name=name):
                self.assertTrue(all(x.a > 0 for x in budget_margins(witness['parameters']).values()))

    def test_pair_cost_blocks_unearned_exponent(self):
        parameters = copy.deepcopy(self.witnesses['combined']['parameters'])
        parameters['tau'] = '0.489'
        with self.assertRaisesRegex(ValueError, 'filtered_lists|visited_pairs'):
            budget_margins(parameters)

    def test_representations_require_coverage(self):
        parameters = copy.deepcopy(self.witnesses['combined']['parameters'])
        parameters['local_deficiency'] = '0.001'
        with self.assertRaisesRegex(ValueError, 'coverage'):
            budget_margins(parameters)

    def test_disjoint_block_boundary_is_rejected(self):
        parameters = copy.deepcopy(self.witnesses['combined']['parameters'])
        parameters['rho'] = '0.0625'
        with self.assertRaisesRegex(ValueError, 'block_space'):
            budget_margins(parameters)

    def test_final_alpha_must_absorb_checksum_overhead(self):
        parameters = copy.deepcopy(self.witnesses['combined']['parameters'])
        parameters['alpha'] = parameters['tau']
        with self.assertRaisesRegex(ValueError, 'assembly_checksum'):
            budget_margins(parameters)

    def test_unknown_budget_fields_fail(self):
        parameters = dict(self.witnesses['combined']['parameters'], unpriced='0')
        with self.assertRaisesRegex(ValueError, 'unexpected'):
            budget_margins(parameters)

    def test_exact_input_required(self):
        for value in [0.1, True, 'NaN', 'Infinity', 'not-a-number']:
            with self.subTest(value=value), self.assertRaises(ValueError):
                q(value)

    def test_output_enclosures_round_outward(self):
        interval = q('0.12345678901234567890123456789')
        left, right = enclosure(interval)
        self.assertLessEqual(Decimal(left), Decimal('0.12345678901234567890123456789'))
        self.assertGreaterEqual(Decimal(right), Decimal('0.12345678901234567890123456789'))

    def test_empty_or_duplicate_witnesses_are_rejected(self):
        source = self.witnesses['source']
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'witnesses.json'
            for witnesses in [[], [source, source]]:
                path.write_text(json.dumps({'witnesses': witnesses}))
                with self.subTest(witnesses=len(witnesses)), self.assertRaises(ValueError):
                    load_witnesses(path)


class QuadratureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.witnesses = {w['id']: w for w in load_witnesses()}

    def test_published_scalar_nodes(self):
        p, y = direct_nodes(self.witnesses['source']['design'])
        # Two separated rows of the manuscript's independently printed table.
        self.assertLess(abs(p[0]-1.7255739), 4e-6)
        self.assertLess(abs(y[0]-2.3814527), 2e-5)
        self.assertLess(abs(p[32]-1.0678994), 4e-6)
        self.assertLess(abs(y[32]-.1742636), 2e-5)

    def test_source_moments_and_grid_refinement(self):
        witness = self.witnesses['source']
        coarse = moment_logs(witness['design'], witness['parameters'], 64)
        fine = moment_logs(witness['design'], witness['parameters'], 256)
        self.assertLess(abs(fine[0]-.5149740483510736), 1e-12)
        self.assertLess(abs(fine[1]-.5148514512729696), 1e-12)
        self.assertLess(max(abs(a-b) for a, b in zip(coarse, fine)), 1e-10)

    def test_filter_ablation_improves_both_moments(self):
        original, new = self.witnesses['source'], self.witnesses['filter_only']
        a = moment_logs(original['design'], original['parameters'], 256)
        b = moment_logs(new['design'], new['parameters'], 256)
        self.assertTrue(all(y < x for x, y in zip(a, b)))

    def test_new_filter_is_not_sufficient_without_rebalancing(self):
        # The final assembly still uses the old alpha when only filters change.
        self.assertEqual(self.witnesses['source']['parameters'], self.witnesses['filter_only']['parameters'])
        self.assertLess(Decimal(self.witnesses['combined']['parameters']['alpha']),
                        Decimal(self.witnesses['parameters_only']['parameters']['alpha']))

    def test_uniform_filter_fails_the_rigorous_moment_budget(self):
        witness = self.witnesses['source']
        uniform = {'r': ['0']*6, 's': ['0']*6, 'kappa': '1.0315'}
        with self.assertRaisesRegex(ValueError, 'continuous moment'):
            certify_moments(uniform, witness['parameters'], points=64)

    def test_search_is_bounded_and_can_extend_an_existing_filter(self):
        parameters = self.witnesses['combined']['parameters']
        a = search(parameters, steps=40, degree=12, start='combined')
        b = search(parameters, steps=40, degree=12, start='combined')
        self.assertEqual(a['evaluations'], 41)
        self.assertEqual(a['design'], b['design'])
        self.assertLessEqual(a['candidate_max_log2_moment'], a['starting_max_log2_moment'])

    def test_search_rejects_truncation_and_excessive_budget(self):
        parameters = self.witnesses['combined']['parameters']
        for options in [{'degree': 6, 'start': 'combined'}, {'steps': 100001}]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                search(parameters, **options)


if __name__ == '__main__':
    unittest.main(verbosity=2)
