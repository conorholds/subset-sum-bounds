"""Failure controls, independent quadrature, and exact domination examples."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import multiplier_search as search
import interval_check as check
import frontier_bound
import joint_search
import odd_probe
from certify import certify as legacy_certify
from filters import scores


def candidates():
    return json.loads((HERE / 'candidates.json').read_text())['candidates']


class DominationTests(unittest.TestCase):
    def test_exact_domination_in_both_branches_including_rounding_extremes(self):
        delta, gamma, epsilon = Fraction(7, 250), Fraction(23, 50), Fraction(1, 10**8)
        ratio = delta / gamma
        lam1, lam2 = map(Fraction, candidates()[-1]['multipliers'])
        for gap in [epsilon / 100, epsilon, Fraction(1, 1000), Fraction(1, 5)]:
            rounded_u = delta - gap
            true_u = rounded_u + 3 * epsilon
            self.assertLessEqual(-rounded_u, -delta + lam1 * (delta - true_u) + 30 * epsilon)
            for alpha in [Fraction(1), Fraction(999, 1000), Fraction(1, 3)]:
                rounded_w = gamma / alpha
                rounded_u = ratio * rounded_w - gap
                true_u, true_v = rounded_u + 3 * epsilon, rounded_w - 4 * epsilon
                self.assertLessEqual(-alpha * rounded_u,
                                     -delta + lam2 * (ratio * true_v - true_u) + 30 * epsilon)

    def test_slopes_below_one_have_exact_counterexamples_without_rounding(self):
        delta, gamma, lam = Fraction(7, 250), Fraction(23, 50), Fraction(99, 100)
        gap = Fraction(1, 1000)
        u = delta - gap
        self.assertGreater(-u, -delta + lam * (delta - u))
        alpha = (1 + lam) / 2
        w = gamma / alpha
        u = delta / gamma * w - gap
        self.assertGreater(-alpha * u, -delta + lam * (delta / gamma * w - u))

    def test_list_pair_combination_is_a_lower_bound_for_arbitrary_densities(self):
        check.iv.dps = 50
        h = -check.q('.25') * check.log2(check.q('.25')) - check.q('.75') * check.log2(check.q('.75'))
        for eta in ['.51', '.51798721', '.53']:
            floor, derivative = frontier_bound.cost_floor(check.q(eta))
            self.assertLess(float(derivative.b), 0)
            for rho in ['.001', '.04', '.0577721', '.06', '.12']:
                left = (.5 + (h - 1 + check.q('.00003')) * check.q(rho)
                        + check.q('.500046') / 17 + check.q('.000016')) / (1 + check.q(1) / 17)
                pair = (1 + (2 * h - check.q('1.5') + check.q('.00003')) * check.q(rho)
                        + check.q('.500046') / 17 - check.q(eta) + check.q('.00002')) / (1 + check.q(1) / 17)
                self.assertGreater(float(max(left.a, pair.a)), float(floor.b))


class NumericalTests(unittest.TestCase):
    def test_receipt_requires_exact_coverage_status_multipliers_and_precision(self):
        original = json.loads((HERE / 'CERTIFICATE.json').read_text())
        self.assertEqual(check.validate_receipt(original), (128, 50))
        mutations = []
        bad = deepcopy(original); bad['results'].pop(); mutations.append(bad)
        bad = deepcopy(original); bad['results'][1] = deepcopy(bad['results'][0]); mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['status'] = 'FAIL'; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['multipliers'] = ['1', '1']; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['moments']['grid'] = True; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['moments']['interval_digits'] = 20; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['moments']['grid'] = 256; mutations.append(bad)
        bad = deepcopy(original); bad['inputs_sha256'] = {}; mutations.append(bad)
        for index, receipt in enumerate(mutations):
            with self.subTest(index=index), self.assertRaises(ValueError):
                check.validate_receipt(receipt)

    def test_interval_replay_and_independent_fft_quadrature_agree(self):
        for candidate in candidates():
            result = check.certify_candidate(candidate)
            self.assertEqual(result['status'], 'PASS')
            p = candidate['parameters']
            eta, tau = float(p['eta']), float(p['tau'])
            for points in [256, 512]:
                u, v, m2 = scores(candidate['design'], points)
                directions = [eta - tau - u, (eta - tau) / (2 * tau - eta) * v - u]
                for direction, lam, interval in zip(directions, candidate['multipliers'],
                                                     result['moments']['log2_moment_intervals']):
                    sampled = search.log_moment(direction, m2, float(lam))
                    self.assertLess(float(interval[0]), sampled)
                    self.assertLess(sampled, float(interval[1]))

    def test_fixed_multiplier_checker_rejects_the_improved_parameters(self):
        candidate = candidates()[-1]
        with self.assertRaisesRegex(ValueError, 'continuous moment'):
            legacy_certify({key: candidate[key] for key in ['id', 'design', 'parameters']})

    def test_invalid_multipliers_and_insufficient_rounding_budgets_fail(self):
        for multipliers, message in [(['.999', '8.3'], 'at least one'),
                                     (['8.5', '.999'], 'at least one'),
                                     (['10', '8.3'], 'domination_one'),
                                     (['8.5', '10'], 'domination_two'),
                                     ([8.5, '8.3'], 'exact decimal'),
                                     (['NaN', '8.3'], 'finite exact')]:
            candidate = dict(candidates()[-1], multipliers=multipliers)
            with self.subTest(multipliers=multipliers), self.assertRaisesRegex(ValueError, message):
                check.certify_candidate(candidate)

    def test_more_aggressive_cost_budget_and_uniform_filters_fail(self):
        candidate = deepcopy(candidates()[-1])
        candidate['parameters']['tau'] = '.4896'
        with self.assertRaisesRegex(ValueError, 'filtered_lists'):
            check.certify_candidate(candidate)
        candidate = deepcopy(candidates()[-1])
        candidate['design'] = {'r': ['0'], 's': ['0'], 'kappa': '1'}
        with self.assertRaisesRegex(ValueError, 'continuous moment'):
            check.certify_candidate(candidate)

    def test_analytic_multiplier_derivative_matches_independent_finite_differences(self):
        candidate = candidates()[-1]
        result = check.enclose_moments(candidate['design'], candidate['parameters'],
                                      candidate['multipliers'], points=256, digits=70)
        eta, tau = map(float, [candidate['parameters'][k] for k in ['eta', 'tau']])
        u, _, m2 = scores(candidate['design'], 512)
        direction, lam = eta - tau - u, float(candidate['multipliers'][0])
        derivative = (search.log_moment(direction, m2, lam + .0001)
                      - search.log_moment(direction, m2, lam - .0001)) / .0002
        bounds = list(map(float, result['first_multiplier_derivative']))
        self.assertLess(abs(derivative - sum(bounds) / 2), 1e-9)

    def test_frontier_lower_bound_is_positive_and_near_the_retained_witness(self):
        bound = frontier_bound.prove_bound()
        self.assertGreater(Decimal(bound['first_moment_residual_lower_bound'][0]), 0)
        floor = Decimal(bound['alpha_lower_bound'][0])
        fixed_candidate = next(row for row in candidates() if row['id'] == 'multipliers_retuned')
        retained = Decimal(fixed_candidate['parameters']['alpha'])
        self.assertGreater(retained, floor)
        self.assertLess(retained - floor, Decimal('.00000013'))


class SearchTests(unittest.TestCase):
    def test_convex_search_handles_both_endpoint_minima(self):
        m2 = np.ones(64)
        low = search.minimize_multiplier(np.ones(64), m2, 9., iterations=12)
        high = search.minimize_multiplier(-np.ones(64), m2, 9., iterations=12)
        self.assertEqual(low['multiplier'], 1.)
        self.assertEqual(high['multiplier'], 9.)
        self.assertEqual(low['evaluations'], 2)

    def test_search_limits_are_bounded(self):
        for kwargs in [{'points': 1024}, {'points': True}, {'bisections': 46},
                       {'bisections': 0}, {'multiplier_iterations': 61}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                search.search(**kwargs)
        with self.assertRaises(ValueError):
            search.minimize_multiplier(np.ones(64), np.ones(64), 1.)
        for eta in [float('inf'), .55]:
            with self.assertRaises(ValueError):
                search.slice_parameters(eta)

    def test_search_and_certificate_evidence_bind_exact_inputs(self):
        discovery = json.loads((HERE / 'discovery.json').read_text())
        self.assertEqual(discovery['search_sha256'], hashlib.sha256((HERE / 'multiplier_search.py').read_bytes()).hexdigest())
        self.assertEqual(discovery['predecessor_sha256'], hashlib.sha256((HERE / 'predecessor.json').read_bytes()).hexdigest())
        receipt = json.loads((HERE / 'CERTIFICATE.json').read_text())
        self.assertEqual(receipt['inputs_sha256'], check.input_hashes())
        bound = json.loads((HERE / 'BOUND.json').read_text())
        self.assertEqual(bound['inputs_sha256'], check.input_hashes())
        self.assertEqual(bound['bounder_sha256'], hashlib.sha256((HERE / 'frontier_bound.py').read_bytes()).hexdigest())
        predecessor = search.PREDECESSOR['witness']
        for candidate in candidates():
            if not candidate['id'].startswith('joint_'):
                self.assertEqual(candidate['design'], predecessor['design'])
            for key in ['guard', 'epsilon', 'blocks']:
                self.assertEqual(candidate['parameters'][key], predecessor['parameters'][key])
            p = candidate['parameters']
            self.assertEqual(Decimal(p['eta']) - Decimal(p['moment_ceiling']), Decimal('.0001'))
            self.assertEqual(Decimal(p['alpha']) - Decimal(p['tau']), Decimal('.00001'))
        self.assertLess(Decimal(candidates()[-1]['parameters']['alpha']), Decimal(predecessor['parameters']['alpha']))


class JointFrontierTests(unittest.TestCase):
    def test_all_implicit_frontier_derivatives_against_independent_perturbations(self):
        vectors = [joint_search.initial_vector()]
        for degree in [14, 20]:
            report = json.loads((HERE / f'joint-degree{degree}.json').read_text())
            candidate = report['candidate']
            design = candidate['design']
            vectors.append(np.array([*map(float, design['r']), *map(float, design['s']),
                                     float(design['kappa']), *map(float, candidate['multipliers'])]))
        for vector in vectors:
            for points in [128, 256]:
                _, jacobian = joint_search.frontier_and_jacobian(vector, points)
                for index in range(len(vector)):
                    step = np.zeros(len(vector))
                    step[index] = .00001
                    high = joint_search.frontier_and_jacobian(vector + step, points)[0]
                    low = joint_search.frontier_and_jacobian(vector - step, points)[0]
                    with self.subTest(variables=len(vector), points=points, coefficient=index):
                        np.testing.assert_allclose(jacobian[:, index], (high - low) / .00002, atol=2e-9, rtol=0)

    def test_joint_evaluation_caps_include_rounded_retention(self):
        from unittest.mock import patch
        for limit in [1, 2, 3, 8]:
            with self.subTest(limit=limit), patch.object(joint_search, 'frontier_and_jacobian',
                                                       wraps=joint_search.frontier_and_jacobian) as scorer:
                result = joint_search.optimize(iterations=2, max_evaluations=limit)
                self.assertEqual(scorer.call_count, result['evaluations'])
                self.assertLessEqual(result['evaluations'], limit)
                self.assertGreaterEqual(min(result['candidate_boundaries']), min(result['initial_boundaries']) - 1e-10)

    def test_joint_inputs_and_all_exact_retained_candidates_are_bound(self):
        retained = {row['id']: row for row in candidates()}
        for degree in [14, 20]:
            report = json.loads((HERE / f'joint-degree{degree}.json').read_text())
            for name, digest in report['inputs_sha256'].items():
                self.assertEqual(digest, hashlib.sha256((HERE / name).read_bytes()).hexdigest())
            for name, digest in report['shared_sha256'].items():
                self.assertEqual(digest, hashlib.sha256((search.CORE / name).read_bytes()).hexdigest())
            self.assertEqual(report['candidate'], retained[f'joint_degree{degree}_retuned'])
            self.assertLessEqual(report['evaluations'], report['limits']['evaluations'])
            self.assertGreater(min(report['candidate_boundaries']), min(report['initial_boundaries']))
            if report['start']:
                path = HERE / report['start']['path']
                self.assertEqual(report['start']['sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
        floor = Decimal(json.loads((HERE / 'BOUND.json').read_text())['alpha_lower_bound'][0])
        self.assertLess(Decimal(retained['joint_degree14_retuned']['parameters']['alpha']), floor)

    def test_invalid_joint_inputs_fail_and_embedding_preserves_the_function(self):
        for vector in [[], [0.] * 6, [0., 0., 1., float('nan'), 8.], [0., 0., 1., .9, 8.]]:
            with self.subTest(vector=vector), self.assertRaises(ValueError):
                joint_search.frontier_and_jacobian(vector)
        for arguments in [{'degree': 15}, {'degree': True}, {'iterations': 0},
                          {'max_evaluations': 4001}, {'max_evaluations': 0}]:
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                joint_search.optimize(**arguments)
        low = joint_search.initial_vector(degree=14, start=HERE / 'joint-degree14.json')
        high = joint_search.initial_vector(degree=20, start=HERE / 'joint-degree14.json')
        np.testing.assert_allclose(joint_search.frontier_and_jacobian(low)[0],
                                   joint_search.frontier_and_jacobian(high)[0], atol=1e-14, rtol=0)

    def test_original_source_start_reaches_the_same_bounded_frontier(self):
        report = json.loads((HERE / 'joint-source-start.json').read_text())
        for name, digest in report['inputs_sha256'].items():
            self.assertEqual(digest, hashlib.sha256((HERE / name).read_bytes()).hexdigest())
        start = HERE / report['start']['path']
        self.assertEqual(report['start']['sha256'], hashlib.sha256(start.read_bytes()).hexdigest())
        original = json.loads(start.read_text())
        core_source = search.CORE / 'witnesses.json'
        self.assertEqual(original['source_sha256'], hashlib.sha256(core_source.read_bytes()).hexdigest())
        self.assertEqual(original['candidate'], next(w for w in json.loads(core_source.read_text())['witnesses'] if w['id'] == 'source'))
        default = json.loads((HERE / 'joint-degree20.json').read_text())
        self.assertLess(abs(min(report['candidate_boundaries']) - min(default['candidate_boundaries'])), 1e-9)

    def test_odd_probe_replays_the_retained_nonimproving_steps_and_reflection(self):
        control = json.loads((HERE / 'odd-control.json').read_text())
        for name, digest in control['inputs_sha256'].items():
            self.assertEqual(digest, hashlib.sha256((HERE / name).read_bytes()).hexdigest())
        direction = np.array(control['direction'])
        baseline = odd_probe.boundaries(np.zeros(len(direction)))
        np.testing.assert_allclose(baseline, control['origin'], atol=1e-14, rtol=0)
        for trial in control['trials']:
            plus = odd_probe.boundaries(trial['amplitude'] * direction)
            minus = odd_probe.boundaries(-trial['amplitude'] * direction)
            np.testing.assert_allclose(plus, minus, atol=1e-14, rtol=0)
            np.testing.assert_allclose(plus, trial['boundaries'], atol=1e-14, rtol=0)
            self.assertLess(min(plus), min(baseline))


if __name__ == '__main__':
    unittest.main(verbosity=2)
