"""Independent finite oracles and rejection controls for fractional coverage."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
from itertools import combinations
import json
from math import comb, prod
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


audit = load('fractional_coverage_exact_audit', 'audit.py')
check = load('fractional_coverage_interval_tests', 'interval_check.py')
search = load('fractional_coverage_search_tests', 'search.py')
bound = load('fractional_coverage_bound_tests', 'frontier_bound.py')


def candidates():
    return json.loads((HERE/'candidates.json').read_text())['candidates']


def brute_sums(weights):
    return {sum(weights[j] for j in range(len(weights)) if mask & (1 << j))
            for mask in range(1 << len(weights))}


class ProjectionTests(unittest.TestCase):
    def test_canonical_projections_are_injective_for_signed_and_repeated_weights(self):
        for weights in [[1, 1, 2], [-2, -1, 0, 1, 2], [0, 0, 0], [1, 2, 4, 8]]:
            representatives, _ = audit.canonical_masks(weights)
            self.assertEqual(set(representatives), brute_sums(weights))
            for mask in range(1 << len(weights)):
                projection = {value & mask for value in representatives.values()}
                projected_sums = {sum(weights[i] for i in range(len(weights)) if value & (1 << i))
                                  for value in projection}
                self.assertEqual(len(projection), len(projected_sums))

    def test_arbitrary_representatives_do_not_supply_the_projection_argument(self):
        weights = [1, 1, 2]
        representatives = [0, 1, 3, 6, 7]  # One valid representative of each sum 0,...,4.
        self.assertEqual({sum(weights[i] for i in range(3) if mask & (1 << i))
                          for mask in representatives}, set(range(5)))
        projection = {mask & 3 for mask in representatives}
        self.assertEqual(len(projection), 4)
        self.assertEqual(len(brute_sums(weights[:2])), 3)

    def test_the_uniform_representative_law_is_not_assumed_independent(self):
        representatives, _ = audit.canonical_masks([1, 1])
        masks = list(representatives.values())
        first_zero = Fraction(sum(not mask & 1 for mask in masks), len(masks))
        second_one = Fraction(sum(bool(mask & 2) for mask in masks), len(masks))
        joint = Fraction(sum(not mask & 1 and bool(mask & 2) for mask in masks), len(masks))
        self.assertGreater(first_zero*second_one, joint)

    def test_exact_uniform_cover_products_on_independent_enumerations(self):
        for weights in [[1, 1, 1, 1], [-2, 0, 1, 2], [0, 1, 2, 4]]:
            n = len(weights)
            whole = len(brute_sums(weights))
            for q in range(1, n+1):
                product = prod(len(brute_sums([weights[i] for i in indexes]))
                               for indexes in combinations(range(n), q))
                self.assertLessEqual(whole**comb(n-1, q-1), product)

    def test_full_exact_audit_matches_the_committed_receipt(self):
        receipt = json.loads((HERE/'AUDIT.json').read_text())
        replay = audit.run(receipt['limits']['maximum_n'])
        self.assertEqual(replay, receipt)
        n = receipt['limits']['maximum_n']
        self.assertEqual(receipt['cardinality_controls']['weight_multisets'], comb(n+5, 5))
        expected = sum(comb(k+4, 4)*(1 << k) for k in range(n+1))
        self.assertEqual(receipt['cardinality_controls']['canonical_projections'], expected)

    def test_superadditivity_alone_does_not_prove_the_new_inequality(self):
        result = audit.superadditivity_counterexample()
        self.assertGreater(Fraction(result['small_expectation']), Fraction(result['claimed_fractional_upper']))
        self.assertEqual(3//3, 1)
        self.assertEqual(5//3, 1)

    def test_zero_and_power_family_attains_the_fractional_coefficient(self):
        result = audit.sharpness_controls()
        self.assertEqual(result['hypergeometric_cases'], 818)
        for weights in [[0, 0, 1, 2], [0, 1, 2, 4], [0, 0, 0, 0]]:
            n, zeros = len(weights), weights.count(0)
            for q in range(n+1):
                mean = sum((Fraction(q-(len(brute_sums([weights[i] for i in indexes])).bit_length()-1))
                            for indexes in combinations(range(n), q)), Fraction(0))/comb(n, q)
                self.assertEqual(mean, Fraction(q*zeros, n))

    def test_positive_transformed_example_by_independent_full_half_enumeration(self):
        toy = audit.positive_toy()
        weights, n = toy['transformed_weights'], toy['n']
        half = n//2
        actual_product = prod(len(brute_sums([weights[i] for i in indexes]))
                              for indexes in combinations(range(n), half))
        d = Fraction(toy['d'])
        self.assertGreaterEqual(actual_product**d.denominator,
                                1 << ((half*d.denominator-d.numerator*n)*comb(n, half)))
        delta, m = Fraction(toy['delta']), toy['m']
        for indexes in combinations(range(n), m//2):
            size = len(brute_sums([weights[i] for i in indexes]))
            self.assertGreaterEqual(size**delta.denominator, 1 << (m//2*delta.denominator-delta.numerator*m))
        self.assertLess(Fraction(toy['old_coverage_margin']), 0)
        self.assertEqual(toy['uniform_useful_restriction_probability'], '1')

    def test_hypergeometric_mgf_is_bounded_by_the_binomial_mgf(self):
        for n in range(2, 15):
            for small in range(n+1):
                p = Fraction(small, n)
                for q in range(1, n+1):
                    for t in [Fraction(1, 4), Fraction(1, 2), Fraction(2)]:
                        mgf = sum((Fraction(comb(small, k)*comb(n-small, q-k), comb(n, q))*t**k
                                   for k in range(max(0, q-(n-small)), min(q, small)+1)), Fraction(0))
                        self.assertLessEqual(mgf, (1-p+p*t)**q)

    def test_below_threshold_obstruction_has_compatible_guarded_sizes(self):
        p = candidates()[1]['parameters']
        result = audit.below_threshold_obstruction(p)
        n = result['first_compatible_n_after_cutoffs']
        rho, d, delta = [Fraction(result[key]) for key in ['rho', 'd', 'smaller_delta']]
        self.assertEqual((2*d*n).denominator, 1)
        self.assertEqual((rho*n/2).denominator, 1)
        self.assertGreater(Fraction(result['base_two_probability_exponent_per_n']), 0)
        self.assertGreater(Fraction(result['exponent_at_first_compatible_n']), 300)
        self.assertGreaterEqual(n, result['bit_guard_n_cutoff'])
        q0 = result['q_logarithm_cutoff']
        self.assertGreaterEqual((d-delta)*q0, q0.bit_length())

    def test_probability_bounds_cover_zero_blocks_zero_deficiency_and_odd_n(self):
        self.assertEqual(audit.useful_probability({'rho': '.1', 'deficiency': '0', 'local_deficiency': '.1'}, 2), 0)
        self.assertEqual(audit.useful_probability({'rho': '.4', 'deficiency': '0', 'local_deficiency': '.1'}, 5), 1)

    def test_component_cutoff_preserves_both_parities_and_larger_sizes(self):
        for row in candidates():
            p = row['parameters']
            cutoff = audit.half_probability_cutoff(p)
            for n in [cutoff, cutoff+1, cutoff+2, 2*cutoff, 10*cutoff+1]:
                self.assertGreaterEqual(audit.useful_probability(p, n), Fraction(1, 2))

    def test_exponential_tail_keeps_floors_and_gives_nontrivial_high_confidence(self):
        p = candidates()[1]['parameters']
        self.assertGreater(audit.exponential_failure_exponent(p, 10**12), 170)
        self.assertGreater(audit.exponential_failure_exponent(p, 10**12+1), 170)
        self.assertEqual(audit.exponential_failure_exponent(
            {'rho': '.1', 'deficiency': '0', 'local_deficiency': '.1'}, 2), 0)

    def test_invalid_or_unbounded_finite_campaigns_are_rejected(self):
        for n in [-1, 9, True]:
            with self.assertRaises(ValueError):
                audit.cardinality_controls(n)
        for n in [4, 33, True]:
            with self.assertRaises(ValueError):
                audit.probability_controls(n)
        with self.assertRaises(ValueError):
            audit.supports([1]*9)
        with self.assertRaises(ValueError):
            audit.positive_family(65, 10)
        with self.assertRaises(ValueError):
            audit.useful_probability({'rho': '.1', 'deficiency': '.2', 'local_deficiency': '.2'}, 100)


class ModelTests(unittest.TestCase):
    def test_source_and_exact_predecessor_bindings_match_the_retained_files(self):
        snapshot = json.loads((HERE/'predecessor.json').read_text())
        source = HERE.parent/'chernoff-multipliers/candidates.json'
        self.assertEqual(snapshot['sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertIn(snapshot['candidate'], json.loads(source.read_text())['candidates'])
        manifest = HERE.parent/'filter-budgets/SOURCE.json'
        self.assertEqual(snapshot['source_manifest_sha256'], hashlib.sha256(manifest.read_bytes()).hexdigest())
        self.assertEqual(json.loads(manifest.read_text())['commit'], '1581d4cbdd69a0faf34b2089614fe2a814316a6a')

    def test_common_budget_margins_are_unchanged(self):
        predecessor = json.loads((HERE/'predecessor.json').read_text())['candidate']
        check.CHERNOFF.iv.dps = 50
        old = check.CHERNOFF.budget_margins(predecessor['parameters'])
        new = check.budget.budget_margins(candidates()[0]['parameters'])
        for key in set(old)-{'coverage', 'block_space'}:
            self.assertEqual(check.CHERNOFF.enclosure(old[key]), check.CHERNOFF.enclosure(new[key]))

    def test_fixed_control_changes_neither_filters_nor_checksum_multipliers(self):
        before = json.loads((HERE/'predecessor.json').read_text())['candidate']
        fixed = candidates()[1]
        self.assertEqual(fixed['design'], before['design'])
        self.assertEqual(fixed['multipliers'], before['multipliers'])
        self.assertEqual(Decimal(before['parameters']['alpha'])-Decimal(fixed['parameters']['alpha']), Decimal('.00000599'))
        for row in candidates():
            self.assertEqual(row['design'], before['design'])

    def test_old_integer_packing_checker_rejects_every_possible_J(self):
        candidate = candidates()[1]
        rho = Fraction(candidate['parameters']['rho'])
        maximum = audit.ceil(1/rho)-1
        self.assertEqual(maximum, 17)
        for j in range(1, maximum+2):
            old = deepcopy(candidate)
            old['parameters']['blocks'] = j
            with self.subTest(j=j), self.assertRaisesRegex(ValueError, 'coverage|block_space'):
                check.CHERNOFF.certify_candidate(old)

    def test_new_model_rejects_legacy_schema_missing_fields_and_nonpositive_coverage(self):
        original = candidates()[1]
        for key, value in [('blocks', 17), ('local_deficiency', original['parameters']['deficiency']),
                           ('guard', '0'), ('rho', '.5'), ('eta', float(original['parameters']['eta']))]:
            candidate = deepcopy(original)
            candidate['parameters'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                check.certify_candidate(candidate)
        bad = deepcopy(original)
        bad['parameters'].pop('gamma')
        with self.assertRaises(ValueError):
            check.certify_candidate(bad)

    def test_invalid_multipliers_cost_budgets_and_uniform_filters_fail(self):
        original = candidates()[1]
        for multipliers in [['.99', '8.3'], ['10', '8.3'], ['8.5', '10'], ['NaN', '8.3'], [8.5, '8.3']]:
            with self.subTest(multipliers=multipliers), self.assertRaises(ValueError):
                check.certify_candidate(dict(original, multipliers=multipliers))
        bad = deepcopy(original)
        bad['parameters']['tau'] = '.4896'
        with self.assertRaisesRegex(ValueError, 'filtered_lists'):
            check.certify_candidate(bad)
        bad = deepcopy(original)
        bad['design'] = {'r': ['0'], 's': ['0'], 'kappa': '1'}
        with self.assertRaisesRegex(ValueError, 'continuous moment'):
            check.certify_candidate(bad)

    def test_interval_replays_match_independent_fft_quadrature(self):
        receipt = json.loads((HERE/'CERTIFICATE.json').read_text())
        points, digits = check.validate_receipt(receipt)
        fresh = check.replay(points, digits)
        self.assertEqual(len(fresh['results']), 3)
        for candidate, row in zip(candidates(), fresh['results']):
            p = candidate['parameters']
            eta, tau = float(p['eta']), float(p['tau'])
            u, v, m2 = search.previous.scores(candidate['design'], 512)
            directions = [eta-tau-u, (eta-tau)/(2*tau-eta)*v-u]
            for direction, lam, interval in zip(directions, candidate['multipliers'], row['moments']['log2_moment_intervals']):
                sampled = search.previous.log_moment(direction, m2, float(lam))
                self.assertLess(float(interval[0]), sampled)
                self.assertLess(sampled, float(interval[1]))

    def test_all_retained_reserves_are_at_least_the_stated_minima(self):
        for row in candidates()[1:]:
            p = {key: Fraction(value) for key, value in row['parameters'].items()}
            self.assertGreaterEqual(p['deficiency']-(Fraction(1, 2)-p['tau']), Fraction(1, 50000))
            self.assertGreaterEqual(p['rho']*(p['local_deficiency']-p['deficiency']), Fraction(1, 625000))
            self.assertGreaterEqual(Fraction(1, 2)-p['local_deficiency']-p['gamma'], Fraction(3, 100000))
            self.assertEqual(p['eta']-p['moment_ceiling'], Fraction(1, 10000))
            self.assertEqual(p['alpha']-p['tau'], Fraction(1, 100000))
            self.assertGreater(Fraction(1, 625000), Fraction(26, 1000000*17))

    def test_receipt_rejects_stale_inputs_wrong_model_missing_rows_and_precision_changes(self):
        receipt = json.loads((HERE/'CERTIFICATE.json').read_text())
        mutations = []
        bad = deepcopy(receipt); bad['model'] = 'chernoff-multipliers-v1'; mutations.append(bad)
        bad = deepcopy(receipt); bad['inputs_sha256'] = {}; mutations.append(bad)
        bad = deepcopy(receipt); bad['results'].pop(); mutations.append(bad)
        bad = deepcopy(receipt); bad['results'][1] = deepcopy(bad['results'][0]); mutations.append(bad)
        bad = deepcopy(receipt); bad['results'][0]['status'] = 'FAIL'; mutations.append(bad)
        bad = deepcopy(receipt); bad['results'][0]['multipliers'] = ['1', '1']; mutations.append(bad)
        bad = deepcopy(receipt); bad['results'][0]['moments']['grid'] = True; mutations.append(bad)
        bad = deepcopy(receipt); bad['results'][0]['moments']['grid'] = 256; mutations.append(bad)
        bad = deepcopy(receipt); bad['results'][0]['moments']['interval_digits'] = 20; mutations.append(bad)
        for bad in mutations:
            with self.assertRaises(ValueError):
                check.validate_receipt(bad)

    def test_bounded_fixed_and_tuned_searches_reproduce_the_retained_candidates(self):
        start = json.loads((HERE/'predecessor.json').read_text())['candidate']
        for tune, name in [(False, 'fixed-search.json'), (True, 'retuned-search.json')]:
            replay = search.search(start, tune)
            receipt = json.loads((HERE/name).read_text())
            self.assertEqual(replay['candidate'], receipt['candidate'])
            self.assertEqual(replay['moment_evaluations'], receipt['moment_evaluations'])
            self.assertEqual(replay['inputs_sha256'], receipt['inputs_sha256'])

    def test_search_limits_and_bad_parameter_slices_fail(self):
        start = candidates()[0]
        for kwargs in [{'points': 64}, {'points': True}, {'bisections': 46},
                       {'bisections': True}, {'tune_multipliers': 1}, {'author': ''}]:
            with self.assertRaises(ValueError):
                search.search(start, **kwargs)
        for eta in [float('nan'), .4, .8]:
            with self.assertRaises(ValueError):
                search.slice_parameters(eta)


class FullFrontierTests(unittest.TestCase):
    def test_full_fixed_filter_bound_and_independent_tilted_derivatives(self):
        result = bound.prove_bound()
        self.assertGreater(Decimal(result['second_moment_residual_lower_bound'][0]), 0)
        self.assertGreater(Decimal(result['larger_eta_monotonicity_margin'][0]), 0)
        alpha_gap = Decimal(candidates()[1]['parameters']['alpha'])-Decimal(result['alpha_lower_bound'][0])
        self.assertGreater(alpha_gap, 0)
        self.assertLess(alpha_gap, Decimal('.00000012'))
        eta, tau = float(result['eta_cutoff']), float(result['tau_probe'])
        candidate = json.loads((HERE/'predecessor.json').read_text())['candidate']
        u, v, m2 = search.previous.scores(candidate['design'], 1024)
        direction = (eta-tau)/(2*tau-eta)*v-u
        lam = float(result['support_multiplier'])
        _, derivative = search.previous.log_moment_and_derivative(direction, m2, lam)
        self.assertLess(float(result['lambda_derivative_enclosure'][0]), derivative)
        self.assertLess(derivative, float(result['lambda_derivative_enclosure'][1]))
        weights = m2*np.exp2(lam*direction)
        self.assertGreater(float(weights@v/weights.sum()), float(result['tilted_v_lower_bound'][0]))

    def test_rational_cost_floor_matches_the_balanced_slice_and_bounds_other_densities(self):
        bound.iv.dps = 70
        for eta in ['.406', '.5', '.51799028', '.55', '.7']:
            floor, derivative = bound.cost_floor(bound.q(eta))
            self.assertLess(float(derivative.b), 0)
            h = -bound.q('.25')*bound.budget.log2(bound.q('.25'))-bound.q('.75')*bound.budget.log2(bound.q('.75'))
            for rho in ['.001', '.03', '.05778', '.1', '.49']:
                r = bound.q(rho)
                left = (bound.q('.5000176')+(h-bound.q('.5')+bound.q('.00005'))*r)/(1+r)
                pair = (bound.q('1.0000216')+(2*h-1+bound.q('.00005'))*r-bound.q(eta))/(1+r)
                self.assertGreaterEqual(max(float(left.b), float(pair.b)), float(floor.a))
        for eta in [.5178, .51799, .5181]:
            row = search.slice_parameters(eta)
            floor, _ = bound.cost_floor(bound.q(str(eta)))
            self.assertLess(abs(row['tau']-float(floor.a)), 1e-15)

    def test_bound_rejects_bad_precision_and_an_unproved_earlier_cutoff(self):
        for kwargs in [{'points': 64}, {'points': 129}, {'digits': 40}]:
            with self.assertRaises(ValueError):
                bound.prove_bound(**kwargs)
        with self.assertRaisesRegex(ValueError, 'infeasibility not certified'):
            bound.prove_bound(eta='.51798', points=128, digits=50)


if __name__ == '__main__':
    unittest.main()
