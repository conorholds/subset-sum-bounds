"""Independent derivative controls, bounded-search checks, and retained-budget tests."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import optimize
import retune
import compare
from certify import budget_margins
from filters import moment_logs


class DerivativeTests(unittest.TestCase):
    def test_all_derivatives_against_the_independent_shared_scorer(self):
        for degree in [10, 14]:
            design = json.loads((HERE/f'discovery-degree{degree}.json').read_text())['design']
            retained = np.array([*map(float, design['r']), *map(float, design['s']), float(design['kappa'])])
            for x in [optimize.initial_vector(degree), retained]:
                for points in [128, 256]:
                    logs, jacobian = optimize.moments_and_jacobian(x, points=points)
                    np.testing.assert_allclose(logs, moment_logs(optimize.unpack(x), optimize.PARAMETERS, points), atol=2e-12, rtol=0)
                    for index in range(len(x)):
                        step = np.zeros(len(x))
                        step[index] = .00001
                        high = np.array(moment_logs(optimize.unpack(x+step), optimize.PARAMETERS, points))
                        low = np.array(moment_logs(optimize.unpack(x-step), optimize.PARAMETERS, points))
                        with self.subTest(degree=degree, points=points, coefficient=index):
                            np.testing.assert_allclose(jacobian[:, index], (high-low)/.00002, atol=2e-8, rtol=0)

    def test_uniform_density_has_known_moments_and_zero_first_derivatives(self):
        logs, jacobian = optimize.moments_and_jacobian([0., 0., 1.])
        expected = [1+8.5*(float(optimize.PARAMETERS['eta'])-float(optimize.PARAMETERS['tau'])), 1.]
        np.testing.assert_allclose(logs, expected, atol=1e-13, rtol=0)
        np.testing.assert_allclose(jacobian, 0., atol=1e-13, rtol=0)

    def test_invalid_inputs_fail_explicitly(self):
        for vector in [[], [0., 1.], [0., 0., 0.], [float('nan'), 0., 1.]]:
            with self.subTest(vector=vector), self.assertRaises(ValueError):
                optimize.moments_and_jacobian(vector)
        for points in [63, 65, True]:
            with self.subTest(points=points), self.assertRaises(ValueError):
                optimize.moments_and_jacobian([0., 0., 1.], points=points)
        with self.assertRaises(ValueError):
            optimize.moments_and_jacobian([0., 0., 1.], {'eta': '.5', 'tau': '.5'})

    def test_smoothing_bounds_and_derivative_with_a_tied_and_extreme_pair(self):
        jacobian = np.array([[2., -1.], [-1., 3.]])
        value, gradient = optimize.smooth_max(np.array([.5, .5]), jacobian, .001)
        self.assertAlmostEqual(value, .5+.001*np.log(2))
        np.testing.assert_allclose(gradient, [.5, 1.])
        value, gradient = optimize.smooth_max(np.array([.5, -1000.]), jacobian, .001)
        self.assertEqual(value, .5)
        np.testing.assert_allclose(gradient, jacobian[0])
        for temperature in [0., -1., float('nan'), float('inf')]:
            with self.subTest(temperature=temperature), self.assertRaises(ValueError):
                optimize.smooth_max(np.array([.5, .5]), jacobian, temperature)


class SearchTests(unittest.TestCase):
    def test_evaluation_caps_include_initial_and_serialized_candidates(self):
        for limit in [1, 2, 3, 40]:
            with self.subTest(limit=limit), patch.object(optimize, 'moments_and_jacobian', wraps=optimize.moments_and_jacobian) as scorer:
                result = optimize.optimize(max_evaluations=limit, iterations=8)
                self.assertEqual(result['evaluations'], scorer.call_count)
                self.assertLessEqual(result['evaluations'], limit)
                self.assertLessEqual(max(result['candidate_log2_moments']), max(result['starting_log2_moments']))

    def test_descent_is_reproducible_and_retains_an_improved_unsmoothed_objective(self):
        first = optimize.optimize(max_evaluations=80, iterations=8)
        second = optimize.optimize(max_evaluations=80, iterations=8)
        self.assertEqual(first['design'], second['design'])
        self.assertEqual(first['evaluations'], second['evaluations'])
        self.assertLess(max(first['candidate_log2_moments']), max(first['starting_log2_moments'])-1e-6)
        np.testing.assert_allclose(first['candidate_log2_moments'], moment_logs(first['design'], optimize.PARAMETERS), atol=2e-12, rtol=0)

    def test_campaign_limits_cannot_be_expanded(self):
        for arguments in [{'degree': 9}, {'degree': 15}, {'degree': True}, {'max_evaluations': 0},
                          {'max_evaluations': 12001}, {'iterations': 0}, {'iterations': 201}]:
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                optimize.optimize(**arguments)
        for steps in [-1, 4001, True]:
            with self.subTest(steps=steps), self.assertRaises(ValueError):
                compare.compare(steps)

    def test_retained_evidence_is_bound_to_the_actual_code_and_inputs(self):
        for degree in [10, 14]:
            result = json.loads((HERE/f'discovery-degree{degree}.json').read_text())
            self.assertEqual(result['optimizer_sha256'], hashlib.sha256((HERE/'optimize.py').read_bytes()).hexdigest())
            self.assertEqual(result['baseline_sha256'], hashlib.sha256((HERE/'baseline.json').read_bytes()).hexdigest())
            for name, digest in result['inputs_sha256'].items():
                self.assertEqual(digest, hashlib.sha256((optimize.CORE/name).read_bytes()).hexdigest())
            np.testing.assert_allclose(result['candidate_log2_moments'], moment_logs(result['design'], optimize.PARAMETERS), atol=2e-12, rtol=0)

    def test_mutation_controls_use_the_same_inputs_and_replay_their_objectives(self):
        control = json.loads((HERE/'mutation-control.json').read_text())
        self.assertEqual(control['search_sha256'], hashlib.sha256((optimize.CORE/'search.py').read_bytes()).hexdigest())
        self.assertEqual(control['baseline_sha256'], hashlib.sha256((HERE/'baseline.json').read_bytes()).hexdigest())
        for result in control['results']:
            self.assertEqual(result['parameters'], optimize.PARAMETERS)
            self.assertEqual(result['start'], 'combined')
            self.assertEqual(result['seed'], 20261015)
            self.assertEqual(result['evaluations'], 4001)
            self.assertAlmostEqual(result['candidate_max_log2_moment'], max(moment_logs(result['design'], optimize.PARAMETERS)), places=12)
        initial = max(moment_logs(optimize.BASELINE['witness']['design'], optimize.PARAMETERS))
        for result in compare.compare(steps=0)['results']:
            self.assertEqual(result['evaluations'], 1)
            self.assertAlmostEqual(result['candidate_max_log2_moment'], initial, places=12)


class RetuningTests(unittest.TestCase):
    def test_retuning_receipt_and_ablations_match_the_retained_search_inputs(self):
        witnesses = {w['id']: w for w in json.loads((HERE/'witnesses.json').read_text())['witnesses']}
        report = json.loads((HERE/'retuning.json').read_text())
        self.assertEqual(report['retuner_sha256'], hashlib.sha256((HERE/'retune.py').read_bytes()).hexdigest())
        for result in report['results']:
            degree = result['degree']
            path = HERE/f'discovery-degree{degree}.json'
            design = json.loads(path.read_text())['design']
            self.assertEqual(result['input_sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
            for suffix in ['filter_only', 'retuned']:
                self.assertEqual(witnesses[f'degree{degree}_{suffix}']['design'], design)
            self.assertEqual(witnesses[f'degree{degree}_filter_only']['parameters'], optimize.PARAMETERS)
            self.assertEqual(witnesses[f'degree{degree}_retuned']['parameters'], result['parameters'])

    def test_all_retained_budgets_and_key_reserves_survive(self):
        witnesses = json.loads((HERE/'witnesses.json').read_text())['witnesses']
        baseline = witnesses[0]
        self.assertEqual(dict(baseline, id='block17_retuned'), optimize.BASELINE['witness'])
        for witness in witnesses:
            with self.subTest(witness=witness['id']):
                p = witness['parameters']
                margins = budget_margins(p)
                self.assertEqual(p['blocks'], 17)
                for key in ['epsilon', 'guard']:
                    self.assertEqual(Decimal(p[key]), Decimal(baseline['parameters'][key]))
                self.assertEqual(Decimal(p['eta'])-Decimal(p['moment_ceiling']), Decimal('.0001'))
                self.assertEqual(Decimal(p['deficiency'])-(Decimal('.5')-Decimal(p['tau'])), Decimal('.00002'))
                self.assertGreater(17*Decimal(p['rho'])*Decimal(p['local_deficiency'])-Decimal(p['deficiency']), Decimal('.000026'))
                self.assertEqual(Decimal('.5')-Decimal(p['local_deficiency'])-Decimal(p['gamma']), Decimal('.00003'))
                self.assertEqual(Decimal(p['preliminary'])-(Decimal('.5')-Decimal(p['deficiency'])), Decimal('.00001'))
                self.assertEqual(Decimal(p['alpha'])-Decimal(p['tau']), Decimal('.00001'))
                self.assertGreater(float(margins['filtered_lists'].a), .000016)
                self.assertGreater(float(margins['visited_pairs'].a), .000010)
                if witness['id'].endswith('_retuned'):
                    self.assertLess(Decimal(p['alpha']), Decimal(baseline['parameters']['alpha']))
                    replayed = retune.retune(witness['design'])
                    self.assertEqual(replayed['parameters'], p)
                    self.assertEqual(replayed['evaluations'], 43)
        degree14 = next(w for w in witnesses if w['id'] == 'degree14_filter_only')
        self.assertTrue(all(a < b for a, b in zip(moment_logs(degree14['design'], optimize.PARAMETERS),
                                                 moment_logs(baseline['design'], optimize.PARAMETERS))))

    def test_retuning_rejects_a_uniform_filter_and_unclipped_slice_violation(self):
        with self.assertRaisesRegex(ValueError, 'not bracketed'):
            retune.retune({'r': ['0']*10, 's': ['0']*10, 'kappa': '1'})
        with self.assertRaisesRegex(ValueError, 'outside'):
            retune.parameters_at(.55)

    def test_a_more_aggressive_cap_fails_the_shared_budget_checker(self):
        witness = json.loads((HERE/'witnesses.json').read_text())['witnesses'][-1]
        parameters = dict(witness['parameters'], tau='0.4896')
        with self.assertRaisesRegex(ValueError, 'filtered_lists'):
            budget_margins(parameters)


if __name__ == '__main__':
    unittest.main(verbosity=2)
