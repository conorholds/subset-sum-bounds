"""Exact thresholds, rounding, likelihood controls and failure paths."""
from copy import deepcopy
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb, isqrt
from pathlib import Path
import unittest

import coverage_audit
import coverage_bounds as bounds
import implementation_bounds as implementation

HERE = Path(__file__).resolve().parent


def candidate():
    return json.loads((HERE.parent/'checksum-variance/FINITE-INPUT.json').read_text())['candidate']


class ProfileTests(unittest.TestCase):
    def test_outward_rounding_preserves_order_and_half_epsilon_error(self):
        epsilon = Fraction(1, 8)
        for d, da, de in product([Fraction(0), Fraction(1, 64), Fraction(1, 3), Fraction(2)], repeat=3):
            actual = (d, d-da, d-de)
            for errors in product([-epsilon/8, Fraction(0), epsilon/8], repeat=3):
                approx = tuple(x+e for x, e in zip(actual, errors))
                rounded = bounds.round_profile(approx, epsilon)
                self.assertGreaterEqual(rounded[0], actual[0])
                self.assertLessEqual(rounded[1], actual[1])
                self.assertLessEqual(rounded[2], actual[2])
                self.assertTrue(all(abs(x-y) <= epsilon/2 for x, y in zip(rounded, actual)))
                self.assertGreaterEqual(rounded[0], 0)
                self.assertLessEqual(max(rounded[1:]), rounded[0])
                self.assertEqual(rounded, bounds.round_profile(approx, epsilon))

    def test_grid_counts_include_both_extreme_approximation_errors(self):
        c = bounds.constants(candidate())
        eps, step = c['epsilon'], c['grid_step']
        absolute = [c['likelihood_width'], c['high_width'], c['low_width']]
        indices = [set(), set(), set()]
        for sign, error in product([-1, 1], [-eps/8, eps/8]):
            for index, width in enumerate(absolute):
                actual = width*sign if index else max(Fraction(0), width*sign)
                approximation = actual+error
                rounded = (step*bounds.ceil((approximation+eps/8)/step) if not index
                           else step*((approximation-eps/8)//step))
                indices[index].add(int(rounded/step))
        for count, attained in zip(c['grid_counts'], indices):
            self.assertLessEqual(max(attained)-min(attained)+1, count)

    def test_bad_rounding_inputs_and_broken_order_promise_fail(self):
        for values, epsilon in [((0, 0, 0), Fraction(1, 8)),
                                 ((Fraction(0),)*2, Fraction(1, 8)),
                                 ([Fraction(0)]*3, Fraction(1, 8)),
                                 ((Fraction(0),)*3, Fraction(0)),
                                 ((Fraction(-1), Fraction(2), Fraction(0)), Fraction(1, 8))]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                bounds.round_profile(values, epsilon)


class PartitionTests(unittest.TestCase):
    def test_balanced_partition_reconstructs_every_size_and_limits_spread(self):
        for n in range(1, 201):
            for q in range(1, n+1):
                p = bounds.balanced(n, q)
                small, large, count = p['minimum'], p['maximum'], p['larger_blocks']
                self.assertEqual((q-count)*small+count*large, n)
                self.assertIn(large-small, [0, 1])

    def test_square_root_rule_at_squares_and_adjacent_values(self):
        for root in range(1, 101):
            for n in [root*root, root*root+1, (root+1)**2-1]:
                p = bounds.square_root_partition(n)
                self.assertEqual(p['q'], isqrt(n)+(isqrt(n)**2 != n))
                self.assertGreaterEqual(p['minimum'], isqrt(n)-1)

    def test_minimum_block_rule_eventually_recovers_the_source_partition(self):
        c = bounds.constants(candidate())
        for power in [122, 123, 160, 256, 512, 1024]:
            for offset in [-1, 0, 1]:
                n = (1 << power)+offset
                p = bounds.minimum_block_partition(c, n)
                self.assertEqual(p, bounds.square_root_partition(n))
                self.assertGreaterEqual(p['minimum'], bounds.required_block_size(c, n))

    def test_minimum_block_rule_at_its_boundary_and_between_square_roots(self):
        c = bounds.constants(candidate())
        d0 = c['minimum_block_size']
        with self.assertRaises(ValueError):
            bounds.minimum_block_partition(c, d0-1)
        for n in [d0, d0+1, 2*d0-1, 2*d0, 3*d0-1, 3*d0, 1 << 61, 1 << 80, 1 << 120]:
            p = bounds.minimum_block_partition(c, n)
            target = bounds.required_block_size(c, n)
            self.assertGreaterEqual(p['minimum'], target)
            self.assertLessEqual(p['q'], bounds.square_root_partition(n)['q'])
            self.assertLessEqual(p['maximum'], max(2*target, bounds.square_root_partition(n)['q']))

    def test_logarithmic_and_pattern_density_minimums_can_each_be_binding(self):
        # Algebraic controls with smaller constants exercise branches that
        # the retained filter reaches only beyond the calculator's size cap.
        c = deepcopy(bounds.constants(candidate()))
        c.update(epsilon=Fraction(1, 8), minimum_block_size=3, pattern_bits_per_block=3)
        n = 256
        self.assertEqual(bounds.required_block_size(c, n), 60)
        self.assertGreaterEqual(bounds.minimum_block_partition(c, n)['minimum'], 60)
        c['pattern_bits_per_block'] = 300
        n = 2048
        self.assertEqual(bounds.required_block_size(c, n), 1204)
        self.assertEqual(bounds.minimum_block_partition(c, n)['q'], 1)
        self.assertEqual(bounds.minimum_block_partition(c, 1 << 40), bounds.square_root_partition(1 << 40))

    def test_invalid_sizes_counts_and_excessive_integers_fail(self):
        for n, q in [(0, 1), (True, 1), (1, False), (2, 3), (-1, 1), (1.0, 1), (1 << 4096, 1)]:
            with self.subTest(n=str(n)[:20]), self.assertRaises(ValueError):
                bounds.balanced(n, q)
        for value in [0, -1, True, Fraction(1), 1 << 4096]:
            with self.assertRaises(ValueError):
                bounds.ceil_log2(value)


class CutoffTests(unittest.TestCase):
    def test_exact_constants_and_independent_coarser_envelopes(self):
        c = bounds.constants(candidate())
        self.assertEqual(c['minimum_block_size'], 2006670879202202015)
        self.assertEqual(c['grid_counts'], [6772500498, 7026149213, 6991195439])
        self.assertEqual(c['pattern_bits_per_block'], 99)
        # Independent coarser envelopes need only small rational arithmetic.
        eps = Fraction(1, 10**8)
        kappa, r, s = Fraction(1035, 1000), Fraction(4244, 1000), Fraction(4223, 1000)
        coarser = [bounds.ceil(Fraction(7, 10)*(2*(r+s))**2/eps**2),
                   bounds.ceil(Fraction(21, 20)*(2*kappa*r)**2/eps**2),
                   bounds.ceil(Fraction(21, 20)*(2*kappa*s)**2/eps**2)]
        self.assertTrue(all(x <= y for x, y in zip(c['concentration_blocks'], coarser)))
        self.assertLess(max(coarser), 1 << 61)
        self.assertGreater(c['minimum_block_size'], 1 << 60)
        # exp(7/10)>2 from four terms, with no floating-point logarithm.
        x = Fraction(7, 10)
        self.assertGreater(1+x+x*x/2+x*x*x/6, 2)

    def test_both_all_n_cutoffs_and_previous_dyadic_failures(self):
        c = bounds.constants(candidate())
        for kind, expected in [('square-root', 122), ('minimum-block', 61)]:
            result = bounds.uniform_cutoff(c, kind)
            self.assertEqual(result['power'], expected)
            self.assertTrue(all(result['uniform_conditions'].values()))
            self.assertFalse(all(result['previous_power_conditions'].values()))
            self.assertEqual(result['at_cutoff']['status'], 'PASS')

    def test_point_bounds_after_dyadic_cutoff_and_partition_transitions(self):
        c = bounds.constants(candidate())
        d0 = c['minimum_block_size']
        for n in [(1 << 61)+1, (1 << 61)+1234567, 3*d0-1, 3*d0, 3*d0+1,
                  1 << 80, (1 << 100)+1, 1 << 121, 1 << 122, 1 << 512]:
            p = bounds.minimum_block_partition(c, n)
            self.assertEqual(bounds.point_certificate(c, n, p)['status'], 'PASS')

    def test_each_probability_failure_is_detected(self):
        c = bounds.constants(candidate())
        n = 1 << 61
        failed = bounds.point_certificate(c, n, bounds.square_root_partition(n))
        self.assertFalse(failed['checks']['concentration'])
        modified = deepcopy(c)
        modified['pattern_bits_per_block'] = n
        self.assertFalse(bounds.point_certificate(modified, n, bounds.balanced(n, 1))['checks']['caps'])
        # Deliberately undersized blocks cannot clear the coverage exponent.
        self.assertFalse(bounds.point_certificate(c, n, bounds.balanced(n, n))['checks']['coverage'])

    def test_exact_ideal_checksum_composition(self):
        result = bounds.ideal_checksum_conditions(candidate(), 61)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['sampling_exponent'], Fraction('0.48943092'))
        self.assertEqual(result['sampling_decay'], Fraction('0.00002488'))
        self.assertEqual(result['combined_failure_at_most'], '2^-127')
        self.assertEqual(bounds.ideal_checksum_conditions(candidate(), 61, 1024)['status'], 'PASS')
        self.assertEqual(bounds.ideal_checksum_conditions(candidate(), 8)['status'], 'FAIL')

    def test_finite_moments_alone_do_not_validate_checksum_domination(self):
        for multipliers in [['0.5', '8'], ['11', '8'], ['8', '16']]:
            invalid = deepcopy(candidate())
            invalid['multipliers'] = multipliers
            # The transfer calculator may bound moments for these slopes,
            # but the checksum composition must reject their reserve loss.
            with self.assertRaises(ValueError):
                bounds.ideal_checksum_conditions(invalid, 61)
            with self.assertRaises(ValueError):
                implementation.finite_conditions(invalid, 61)

    def test_bad_candidate_precision_probability_and_partition_schemas_fail(self):
        c = bounds.constants(candidate())
        invalid = deepcopy(candidate())
        invalid['parameters']['epsilon'] = '.0000000000000000000000000000001'
        with self.assertRaises(ValueError):
            bounds.constants(invalid)
        invalid = deepcopy(candidate())
        invalid['design']['r'][0] = 0.1
        with self.assertRaises(ValueError):
            bounds.constants(invalid)
        for kind, bits in [('other', 128), ('square-root', 0), ('minimum-block', True)]:
            with self.assertRaises(ValueError):
                bounds.uniform_cutoff(c, kind, bits)
        for power in [0, 513, True]:
            with self.assertRaises(ValueError):
                bounds.uniform_conditions(c, power, 'square-root')
        for p in [{'q': 1}, dict(bounds.balanced(2, 1), minimum=1)]:
            with self.assertRaises(ValueError):
                bounds.point_certificate(c, 2, p)
        for degree in [0, True, 1025]:
            with self.assertRaises(ValueError):
                bounds.ideal_checksum_conditions(candidate(), 61, degree)

    def test_huge_input_certificate_uses_exponents_not_huge_powers(self):
        c = bounds.constants(candidate())
        n = 1 << 2048
        p = bounds.minimum_block_partition(c, n)
        self.assertEqual(bounds.point_certificate(c, n, p)['status'], 'PASS')


class ProbabilityControls(unittest.TestCase):
    def test_digit_prefactors_cover_nonintegral_budgets_and_are_not_dropped(self):
        for denominator in [1, 2, 5, 10]:
            for numerator in range(denominator, 20*denominator+1):
                budget = Fraction(numerator, denominator)
                for digit_scale in [Fraction(1), budget/2, budget]:
                    if digit_scale < 1:
                        continue
                    base = 1 << bounds.ceil_log2(bounds.ceil(digit_scale))
                    for modulus in range(2, int(4*budget)+1):
                        high = (modulus+base-1)//base
                        self.assertLessEqual(high, 5*budget/digit_scale)
                        self.assertLessEqual(base, 2*digit_scale)
                        self.assertLessEqual(high*base, 6*budget)
        # Dropping the high-side rounding prefactor can underestimate a list.
        self.assertGreater(Fraction(5), 4*Fraction(83, 10)/8)
        # Three length-16 rows represent 48 pairs for a budget of 8.3.
        self.assertGreater(48, 5*Fraction(83, 10))

    def test_exact_binomial_tails_obey_the_dyadic_hoeffding_consequence(self):
        checks = 0
        for n in range(1, 17):
            for numerator, denominator in [(1, 4), (1, 2), (3, 4)]:
                p = Fraction(numerator, denominator)
                for threshold in range(1, n+1):
                    deviation = Fraction(threshold, n)-p
                    if deviation <= 0:
                        continue
                    tail = sum((comb(n, j)*p**j*(1-p)**(n-j) for j in range(threshold, n+1)), Fraction(0))
                    exponent = 2*n*deviation**2
                    # exp(-x)<=2^-x; clear the rational exponent exactly.
                    self.assertLessEqual(tail**exponent.denominator, Fraction(1, 1 << exponent.numerator))
                    checks += 1
        self.assertEqual(checks, 220)

    def test_halving_the_true_width_produces_a_false_tail_bound(self):
        actual = Fraction(sum(comb(8, j) for j in range(6, 9)), 1 << 8)
        self.assertEqual(actual, Fraction(37, 256))
        self.assertGreater(actual, Fraction(1, 16))

    def test_exact_likelihood_and_typical_event_receipt(self):
        result = coverage_audit.run('conorholds')
        self.assertEqual(result, json.loads((HERE/'AUDIT.json').read_text()))
        self.assertEqual(result['cases'], 9)
        self.assertGreaterEqual(Fraction(result['minimum_typical_mass_lower_4096']), Fraction(1, 2))
        self.assertEqual(result['small_block_typical_mass'], '0')

    def test_exhaustive_small_labeled_law_matches_the_count_recurrence(self):
        for d in range(1, 7):
            alpha = Fraction(1, 2)
            original, tilted, denominator = coverage_audit.distributions(alpha, (d,))[d]
            # Six equally likely original labels: three stars, two low
            # f-values, and one high f-value. The tilted six-label law
            # has one low value and two high values instead.
            for calculated, labels in [(original, [0, 0, 0, -1, -1, 1]),
                                       (tilted, [0, 0, 0, -1, 1, 1])]:
                literal = [0]*(2*d+1)
                for coordinates in product(labels, repeat=d):
                    literal[sum(coordinates)+d] += 1
                self.assertEqual(tuple(literal), calculated)
                self.assertEqual(sum(literal), denominator)

    def test_small_block_and_invalid_laws_fail_meaningfully(self):
        _, tilted, denominator = coverage_audit.distributions(Fraction(1, 2), (1,))[1]
        self.assertEqual(coverage_audit.typical_mass(Fraction(1, 2), 1, Fraction(1, 32), tilted, denominator), 0)
        for alpha, sizes in [(Fraction(0), (1,)), (0.5, (1,)), (Fraction(1, 65), (1,)),
                              (Fraction(1, 2), (513,)), (Fraction(1), (True,)), (Fraction(1), ())]:
            with self.assertRaises(ValueError):
                coverage_audit.distributions(alpha, sizes)
        for operation in [bounds.report, coverage_audit.run]:
            with self.assertRaises(ValueError):
                operation('not-the-approved-handle')


class BindingTests(unittest.TestCase):
    def test_bound_receipt_replays_exactly(self):
        self.assertEqual(bounds.report('conorholds'), json.loads((HERE/'BOUNDS.json').read_text()))

    def test_predecessor_and_source_bindings_remain_exact(self):
        manifest = json.loads((HERE/'INPUTS.json').read_text())
        for name, expected in manifest['predecessor']['sha256'].items():
            self.assertEqual(hashlib.sha256((HERE.parent/'checksum-variance'/name).read_bytes()).hexdigest(), expected)
        preceding = json.loads((HERE.parent/'checksum-variance/INPUTS.json').read_text())
        self.assertEqual(manifest['source_commit'], '1581d4cbdd69a0faf34b2089614fe2a814316a6a')
        # Compare the new source hashes to the complete, hash-bound
        # seven-section predecessor manifest, rather than to this code's output.
        files = {Path(item['path']).name: item['sha256']
                 for item in preceding['mathematical_source']['files']}
        for name, expected in manifest['source_sections'].items():
            self.assertEqual(expected, files[name])
        self.assertEqual(manifest['source_commit'], preceding['mathematical_source']['commit'])


class ImplementationTests(unittest.TestCase):
    def test_independent_coarser_envelope_recovers_the_overhead_power(self):
        n = 1 << 93
        root = isqrt(n)+(isqrt(n)**2 != n)
        # a<17, D0<2^61 and g=99 suffice independently of the
        # implementation's rational expression for the envelope.
        upper = 215*root+1700000001*93+87+34*(1 << 61)+183600000000
        self.assertLess(upper, Fraction(n, 10**8))
        # Already the constant term exceeds the preceding power's budget.
        lower_constant = 2*Fraction(169, 10)*2000000000000000000
        self.assertGreater(lower_constant, Fraction(1 << 92, 10**8))

    def test_uniform_table_pattern_overhead_has_an_explicit_cutoff(self):
        c = bounds.constants(candidate())
        result = implementation.overhead_cutoff(c)
        self.assertEqual(result['power'], 93)
        self.assertLessEqual(result['envelope_at_cutoff_at_most'], c['epsilon']*result['n'])
        for power in [93, 94, 122, 256, 512]:
            for offset in [-1, 0, 1]:
                n = (1 << power)+offset
                if n < result['n']:
                    continue
                p = bounds.minimum_block_partition(c, n)
                actual_upper = (c['pattern_bits_per_block']*p['q']+bounds.ceil_log2(p['q'])+1
                                +(c['d_bound']+4*c['epsilon'])*p['maximum'])
                self.assertLessEqual(actual_upper, c['epsilon']*n)

    def test_smaller_probability_cutoff_does_not_pay_the_table_scan_overhead(self):
        c = bounds.constants(candidate())
        n = 1 << 61
        p = bounds.minimum_block_partition(c, n)
        self.assertGreater((c['d_bound']+4*c['epsilon'])*p['maximum'], c['epsilon']*n)

    def test_finite_error_budget_at_both_component_cutoffs(self):
        for power in [61, 93]:
            result = implementation.finite_conditions(candidate(), power)
            self.assertEqual(result['status'], 'PASS')
            self.assertEqual(result['combined_failure_at_most'], '2^-128')
            self.assertEqual(result['target_degree'], 2)
            self.assertEqual(result['per_target_total_error_at_most'], '1/8')
            self.assertGreaterEqual(result['coupling_failure_bits_at_least_at_cutoff'], 130)
            self.assertGreaterEqual(result['evaluation_error_bits_at_least_at_cutoff'], 5)
            self.assertGreaterEqual(result['rounding_error_bits_at_least_at_cutoff'], 5)

    def test_every_finite_error_prefactor_matches_the_written_accounting(self):
        n = 1 << 61
        result = implementation.finite_conditions(candidate(), 61)
        self.assertEqual(result['partition']['q'], 1)
        c = -(-(Fraction('0.48945580')+20*Fraction(1, 10**8))*n//1)
        self.assertEqual(result['evaluation_error_bits_at_least_at_cutoff'],
                         2*n*n-2*n-3*c-99-62)
        self.assertEqual(result['rounding_error_bits_at_least_at_cutoff'], n*n-c-99-1)
        self.assertEqual(result['uniform_rejections_per_call'], n**3)
        self.assertEqual(result['coupling_failure_bits_at_least_at_cutoff'],
                         n**3-result['ordinary_uniform_call_count_log2_at_most'])

    def test_linear_precision_and_rejection_schedules_keep_every_allowance(self):
        for power in [61, 93, 122, 256]:
            result = implementation.finite_conditions(candidate(), power)
            n = result['n']
            self.assertEqual(result['proposed_evaluation_error_bits_at_least'], 5)
            self.assertEqual(result['proposed_rounding_error_bits_at_least'], 5)
            self.assertEqual(result['proposed_coupling_failure_bits_at_least'], 130+power)
            self.assertLess(result['proposed_unit_phase_accuracy_bits'], 2*n*n)
            self.assertLess(result['proposed_contribution_fractional_bits'], n*n)
            self.assertLess(result['proposed_uniform_rejection_cap'], n**3)
            self.assertLessEqual(result['proposed_unit_phase_accuracy_bits'], 105*n+5)
            self.assertLessEqual(result['proposed_contribution_fractional_bits'], 100*n+6)
            self.assertLessEqual(result['proposed_uniform_rejection_cap'], 119*n+133)

    def test_exact_weighted_product_supports_the_linear_precision_budget(self):
        # A finite pessimistic all-real product: Q=2, C=4, weight=C^2.
        n, c, pattern_bits = 3, 2, 1
        p = 2*n+3*c+pattern_bits+bounds.ceil_log2(n+1)+5
        r = c+pattern_bits+6
        delta = Fraction(1, 1 << p)
        one_product_error = (1+delta)*(2+delta)**n-2**n
        all_terms = (1 << (pattern_bits+c))*Fraction(1 << (2*c))*one_product_error
        all_rounding = Fraction(1 << (pattern_bits+c+1), 1 << r)
        self.assertLessEqual(all_terms, Fraction(1, 32))
        self.assertLessEqual(all_rounding, Fraction(1, 32))
        self.assertLessEqual(all_terms+all_rounding, Fraction(1, 16))
        coarse = Fraction(1, 1 << 10)
        self.assertGreater((1 << (pattern_bits+3*c))*((1+coarse)*(2+coarse)**n-2**n), Fraction(1, 16))

    def test_failure_allocations_and_large_target_lists_are_checked(self):
        self.assertEqual(implementation.finite_conditions(candidate(), 61, 1024, 512)['status'], 'PASS')
        for power, degree, bits in [(0, 1, 130), (513, 1, 130), (61, True, 130),
                                    (61, 1025, 130), (61, 1, 1), (61, 1, True)]:
            with self.assertRaises(ValueError):
                implementation.finite_conditions(candidate(), power, degree, bits)
        with self.assertRaises(ValueError):
            implementation.report('not-the-approved-handle')

    def test_rejection_sampler_has_exact_accepted_distribution(self):
        # Exhaust every proposal stream for small non-power moduli.
        for modulus, attempts in [(3, 4), (5, 3), (6, 3), (7, 3)]:
            base = 1 << bounds.ceil_log2(modulus)
            outcomes = [0]*modulus
            failures = 0
            for proposals in product(range(base), repeat=attempts):
                accepted = next((value for value in proposals if value < modulus), None)
                if accepted is None:
                    failures += 1
                else:
                    outcomes[accepted] += 1
            self.assertEqual(len(set(outcomes)), 1)
            self.assertEqual(failures, (base-modulus)**attempts)
            self.assertLessEqual(Fraction(failures, base**attempts), Fraction(1, 1 << attempts))
            self.assertEqual(sum(outcomes)+failures, base**attempts)

    def test_success_conditioning_biases_an_adaptive_draw_and_three_terms_need_budget(self):
        # A fair bit chooses between no second call and one capped draw
        # below three. One rejected proposal aborts. Conditioning on the
        # entire execution succeeding changes the initial fair bit's law.
        success = Fraction(1, 2)+Fraction(1, 2)*Fraction(3, 4)
        bit_one_given_success = Fraction(1, 2)*Fraction(3, 4)/success
        self.assertEqual(bit_one_given_success, Fraction(3, 7))
        self.assertNotEqual(bit_one_given_success, Fraction(1, 2))
        total = 3*Fraction(1, 1 << 130)
        self.assertLess(total, Fraction(1, 1 << 128))
        self.assertGreater(total, Fraction(1, 1 << 129))

    def test_phase_product_and_complex_rounding_bounds_are_not_omitted(self):
        # The leading phase also needs approximation: n+1, not n factors.
        n, delta = 3, Fraction(1, 100)
        exact_product = Fraction(2)**n
        approximate = (1+delta)*(2+delta)**n
        self.assertLessEqual(abs(approximate-exact_product), (n+1)*delta*(2+delta)**n)
        # Rounding both coordinates down can have norm exceeding one grid step.
        self.assertGreater(2*Fraction(9, 10)**2, 1)
        self.assertLessEqual(2*Fraction(9, 10)**2, 4)

    def test_explicit_phase_evaluator_accuracy_budget_is_rational(self):
        for p in range(1, 65):
            m = 2*(p+8)+256
            # The factorial lower bound uses the final half of its factors.
            taylor_upper = Fraction(128, m)**(m//2)
            self.assertLessEqual(taylor_upper, Fraction(1, 1 << (p+136)))
            angle = Fraction(1, 1 << (2*p+13))
            rounding = Fraction(1, 1 << (p+7))
            self.assertLess(angle+taylor_upper+rounding, Fraction(1, 1 << p))
        self.assertEqual((Fraction(1, 2)+Fraction(1, 3))/(1-Fraction(1, 6)), 1)

    def test_implementation_receipt_replays_exactly(self):
        self.assertEqual(implementation.report('conorholds'), json.loads((HERE/'IMPLEMENTATION.json').read_text()))


if __name__ == '__main__':
    unittest.main()
