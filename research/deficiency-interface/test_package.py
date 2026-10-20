"""Exact controls for the audited probability interface and retained evidence."""
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit

INPUTS = json.loads((HERE/'parameters.json').read_text())
WITNESSES = [w for source in INPUTS['sources'] for w in source['parameters']]


class ExactSupportTests(unittest.TestCase):
    def test_signed_and_colliding_supports_against_direct_examples(self):
        self.assertEqual(audit.subset_sums([0, 1, 1, -2]), {-2, -1, 0, 1, 2})
        self.assertEqual(audit.subset_sums([]), {0})
        self.assertEqual(audit.subset_sums([1, 2, 4]), set(range(8)))

    def test_exhaustive_control_counts_match_independent_combinatorial_totals(self):
        counts = audit.support_controls(3)
        weights = lambda n: comb(n+4, 4)
        self.assertEqual(counts['weight_multisets'], sum(weights(n) for n in range(4)))
        self.assertEqual(counts['supports'], sum(weights(n)*2**n for n in range(4)))
        self.assertEqual(counts['sign_translations'], counts['supports'])
        self.assertEqual(counts['disjoint_unions'], sum(weights(n)*3**n for n in range(4)))
        self.assertEqual(counts['increments'], sum(weights(n)*n*2**(n-1) for n in range(1, 4)))
        self.assertEqual(counts['exchanges'], sum(weights(n)*n*(n-1)*2**(n-2) for n in range(2, 4)))

    def test_flipping_and_balancing_give_uniform_blocks_for_every_toy_alias(self):
        self.assertEqual(audit.sampling_controls(4), sum((n//2)*2**n for n in range(2, 5)))

    def test_oracle_limits_and_inexact_inputs_fail(self):
        for operation in [lambda: audit.support_controls(7), lambda: audit.sampling_controls(7),
                          lambda: audit.variance_controls(11), lambda: audit.subset_sums([0]*17),
                          lambda: audit.subset_sums([True]), lambda: audit.exact(.1)]:
            with self.assertRaises(ValueError):
                operation()


class ProbabilityTests(unittest.TestCase):
    def test_displayed_rational_bounds_round_in_the_claimed_direction(self):
        for value in [Fraction(1, 3), Fraction(-2, 7), Fraction(0), Fraction(1)]:
            lower = audit.outward(value)
            upper = audit.outward(value, upper=True)
            self.assertEqual(audit.exact(lower['exact']), value)
            self.assertLessEqual(audit.exact(lower['outward_decimal']), value)
            self.assertGreaterEqual(audit.exact(upper['outward_decimal']), value)

    def test_variance_on_actual_deficiencies_of_a_zero_and_power_family(self):
        weights = [0]*4+[1, 2, 4, 8]
        values = []
        for selected in combinations(weights, 4):
            size = len(audit.subset_sums(selected))
            self.assertEqual(size & (size-1), 0)
            values.append(4-(size.bit_length()-1))
        mean = Fraction(sum(values), len(values))
        variance = sum(((value-mean)**2 for value in values), Fraction(0))/len(values)
        self.assertEqual(variance, Fraction(4, 7))
        self.assertLessEqual(variance, audit.coupling_variance_bound(8, 4))
        self.assertLessEqual(audit.coupling_variance_bound(8, 4), audit.variance_bound(8, 4))

    def test_complementary_variance_bounds_and_deterministic_endpoints(self):
        self.assertEqual(audit.coupling_variance_bound(4, 2), Fraction(13, 36))
        for n in range(2, 12):
            self.assertEqual(audit.variance_bound(n, 0), 0)
            self.assertEqual(audit.variance_bound(n, n), 0)
            for q in range(n+1):
                self.assertEqual(audit.variance_bound(n, q), audit.variance_bound(n, n-q))
                self.assertLessEqual(audit.coupling_variance_bound(n, q), audit.variance_bound(n, q))
        for n, q in [(0, 0), (2, -1), (2, 3), (True, 0)]:
            with self.subTest(n=n, q=q), self.assertRaises(ValueError):
                audit.variance_bound(n, q)

    def test_finite_useful_bound_against_an_exact_low_deficiency_instance(self):
        parameters = {'preliminary': '3/10', 'deficiency': '13/50', 'rho': '2/5',
                      'blocks': 2, 'local_deficiency': '2/5', 'gamma': '1/20'}
        # A half has expected deficiency four, below (13/50)*16.
        weights = [0]*8+[1 << index for index in range(8)]
        m = 2*((Fraction(2, 5)*16)//2)
        useful = total = 0
        for selected in combinations(weights, m//2):
            support = len(audit.subset_sums(selected))
            deficiency = m//2-(support.bit_length()-1)
            useful += deficiency <= Fraction(2, 5)*m
            total += 1
        actual = Fraction(useful, total)
        bound = audit.useful_probability(parameters, 16)
        self.assertEqual(actual, Fraction(9, 10))
        self.assertEqual(bound, Fraction(2, 15))
        self.assertLessEqual(bound, actual)

    def test_finite_preliminary_bound_against_a_positive_high_deficiency_family(self):
        parameters = {'preliminary': '7/20', 'deficiency': '1/5', 'rho': '2/5',
                      'blocks': 2, 'local_deficiency': '3/10', 'gamma': '1/10'}
        n, half = 128, 64
        # The positive-weight family in the proof has mean deficiency >26>128/5.
        self.assertLess(33, 2**6)
        self.assertGreater(26, Fraction(n, 5))
        cap = 1 << audit.ceil(Fraction(7, 20)*n)
        overflow = sum(comb(half, count)**2 for count in range(half+1)
                       if max((1 << count)*(half-count+1), (1 << (half-count))*(count+1)) > cap)
        actual = Fraction(overflow, comb(n, half))
        bound = audit.preliminary_continue_upper(parameters, n)
        self.assertGreater(actual, 0)
        self.assertLessEqual(actual, bound)
        self.assertLess(bound, 1)

    def test_continuation_can_coexist_with_high_expected_deficiency(self):
        example = audit.continuation_counterexample()
        self.assertEqual(example['uniform_half_cap_continue_probability'], '1/143')
        self.assertEqual(example['cap'], 128)
        self.assertTrue(all(weight > 0 for weight in example['original_weights']))
        self.assertLess(6**5, 2**13)
        for witness in WITNESSES:
            family = audit.scalable_counterexample(witness['parameters'])
            self.assertGreater(family['expected_deficiency_strict_lower'], audit.exact(witness['parameters']['deficiency'])*family['n'])
            self.assertLess(family['cap_exponent'], family['one_half_support_exponent'])
            self.assertGreater(audit.exact(family['uniform_half_cap_continue_probability_lower']), 0)

    def test_all_pinned_parameters_have_positive_gaps_and_effective_cutoffs(self):
        for witness in WITNESSES:
            with self.subTest(witness=witness['id']):
                parameters = witness['parameters']
                report = audit.parameter_report(parameters)
                self.assertTrue(all(audit.exact(value) > 0 for value in report['gaps'].values()))
                cutoff = report['n_useful_probability_at_least_half']
                self.assertGreaterEqual(audit.useful_probability(parameters, cutoff), Fraction(1, 2))
                self.assertGreater(audit.useful_probability(parameters, 10**12), Fraction(99, 100))
                self.assertLess(audit.preliminary_continue_upper(parameters, 10**12), Fraction(1, 500))
                rho = audit.exact(parameters['rho'])
                n = report['n_absorbing_cardinality_loss']
                m = 2*((rho*n)//2)
                self.assertGreaterEqual(m, report['dyadic_m_absorbing_cardinality_loss'])
                self.assertGreaterEqual(audit.gaps(parameters)['representations']*m, m.bit_length())

    def test_missing_strict_gaps_and_invalid_sizes_are_rejected(self):
        baseline = WITNESSES[-1]['parameters']
        for changed in [dict(baseline, preliminary=str(Fraction(1, 2)-audit.exact(baseline['deficiency']))),
                        dict(baseline, rho='1'), dict(baseline, local_deficiency='0'),
                        dict(baseline, rho=str(Fraction(1, baseline['blocks']))),
                        dict(baseline, local_deficiency=str(audit.exact(baseline['deficiency'])/(baseline['blocks']*audit.exact(baseline['rho'])))),
                        dict(baseline, gamma=str(Fraction(1, 2)-audit.exact(baseline['local_deficiency']))),
                        dict(baseline, blocks=True)]:
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                audit.gaps(changed)
        for n in [0, -1, True]:
            with self.subTest(n=n), self.assertRaises(ValueError):
                audit.useful_probability(baseline, n)
        self.assertEqual(audit.useful_probability(baseline, 1), 0)
        self.assertEqual(audit.preliminary_continue_upper(baseline, 1), 1)
        tiny = dict(baseline, gamma=str(Fraction(1, 2)-audit.exact(baseline['local_deficiency'])-Fraction(1, 10**30)))
        with self.assertRaisesRegex(ValueError, 'bounded dyadic'):
            audit.parameter_report(tiny)


class EvidenceTests(unittest.TestCase):
    def test_committed_exact_receipt_replays_without_numerical_dependencies(self):
        receipt = json.loads((HERE/'AUDIT.json').read_text())
        for name, digest in receipt['sha256'].items():
            self.assertEqual(digest, hashlib.sha256((HERE/name).read_bytes()).hexdigest())
        self.assertEqual(receipt['source_manifest_sha256'], hashlib.sha256((audit.CORE/'SOURCE.json').read_bytes()).hexdigest())
        self.assertEqual(receipt, audit.run())
        self.assertEqual(len(receipt['parameters']), 8)
        self.assertEqual(receipt['sampling_cases'], 300)
        self.assertEqual(receipt['variance_cases'], 501)


if __name__ == '__main__':
    unittest.main(verbosity=2)
