"""Exact happy paths, failure controls and source-bound checksum identities."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb, cos, log2, pi, sin
from pathlib import Path
import unittest

import audit
import finite_circle
from oracle import MatchIndex, Scalar, Tables, character, mean, variance
import reference

HERE = Path(__file__).resolve().parent


def all_star(modulus=4, residues=(1,), phases=(0,), base=2):
    return Tables(modulus, residues, phases, base, [(i,) for i in range(len(residues))],
                  [[(None,)] for _ in residues], [1]*modulus, [1]*modulus,
                  [1]*len(residues), [1]*len(residues))


def repeated_table():
    return Tables(4, [1], [0], 2, [(0,)], [[(None,), (None,), (0,)]],
                  [4, 2, 1, 1], [1, 3, 2, 2], [1], [1])


class CharacterTests(unittest.TestCase):
    def test_exact_characters_and_conjugation_on_all_four_circles(self):
        for modulus in [2, 3, 4, 6]:
            zero, one = Scalar(modulus), Scalar(modulus, 1)
            for a, b in product(range(-2*modulus, 2*modulus), repeat=2):
                self.assertEqual(character(modulus, a)*character(modulus, b), character(modulus, a+b))
                self.assertEqual(character(modulus, a).conjugate(), character(modulus, -a))
                self.assertEqual(character(modulus, a).norm2(), 1)
            for frequency in range(modulus):
                self.assertEqual(sum((character(modulus, frequency*k) for k in range(modulus)), zero),
                                 one*modulus if frequency == 0 else zero)

    def test_nonreal_variance_uses_conjugate_squares(self):
        values = [Scalar(4, 1), Scalar(4, 0, 1)]
        self.assertEqual(variance(values, 4), Fraction(1, 2))
        center = mean(values, 4)
        incorrect = sum(((value-center)*(value-center) for value in values), Scalar(4))/2
        self.assertNotEqual(incorrect, Scalar(4, variance(values, 4)))

    def test_fourier_identity_matches_subset_expansion_for_composite_moduli(self):
        for modulus in [2, 3, 4, 6]:
            for residues, phases in [([0, 0, 0], [0, 1, 0]), ([1, 1, 0], [1, 0, 1]),
                                     ([modulus-1, 1, 1], [0, 1, modulus-1])]:
                table = all_star(modulus, residues, phases)
                expected = reference.subset_checksums(modulus, residues, phases)
                for target in range(modulus):
                    calculated = sum((table.fourier_term(k, target) for k in range(modulus)), Scalar(modulus))/modulus
                    self.assertEqual(calculated, expected[target])
                    self.assertEqual(table.direct_checksum(target), expected[target])

    def test_invalid_exact_arithmetic_is_rejected(self):
        for args in [(5,), (True,), (4, 0.1), (2, 0, 1)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                Scalar(*args)
        for call in [lambda: Scalar(4, 1)+Scalar(3, 1), lambda: Scalar(4, 1)/0,
                     lambda: character(4, 1.0), lambda: character(4, 10**6+1)]:
            with self.assertRaises(ValueError):
                call()


class MatchTests(unittest.TestCase):
    def test_group_prefix_rank_bijection_including_unmatched_keys(self):
        high = [(0, (0,)), (1, (0,)), (2, (1,)), (3, (3,))]
        low = [(0, (0,)), (1, (0,)), (2, (0,)), (0, (1,)), (1, (2,))]
        index = MatchIndex(high, low)
        expected = sorted((h, low_digit, key) for h, key in high for low_digit, other_key in low if key == other_key)
        self.assertEqual(index.total, 7)
        self.assertEqual(sorted(index.match(rank) for rank in range(index.total)), expected)
        self.assertEqual(index.match(5), (1, 2, (0,)))
        self.assertEqual(index.match(6), (2, 0, (1,)))
        self.assertEqual(index.enumerate(), ('PASS', expected))

    def test_large_cartesian_products_are_indexed_without_enumeration(self):
        count = 1024
        index = MatchIndex(((i, (0,)) for i in range(count)), ((i, (0,)) for i in range(count)), cap=count)
        self.assertEqual(index.total, count**2)
        self.assertEqual(len(index.groups), 1)
        self.assertEqual(index.match(count**2-1), (count-1, count-1, (0,)))
        self.assertEqual(index.enumerate(cap=2), ('match_cap', [(0, 0, (0,)), (0, 1, (0,))]))

    def test_empty_matches_and_list_enumeration_caps_fail_before_an_extra_record(self):
        empty = MatchIndex([(0, (0,))], [(0, (1,))])
        self.assertEqual(empty.total, 0)
        self.assertEqual(empty.enumerate(cap=0), ('PASS', []))
        with self.assertRaises(ValueError):
            empty.match(0)
        for high, low, cap, status in [([(0, (0,))], [], 0, 'high_cap'),
                                      ([], [(0, (0,))], 0, 'low_cap'),
                                      ([(0, (0,)), (1, (0,))], [], 1, 'high_cap')]:
            index = MatchIndex(high, low, cap)
            self.assertEqual(index.status, status)
            self.assertEqual(index.total, 0)
            with self.assertRaises(ValueError):
                index.enumerate()

    def test_invalid_entries_duplicate_indices_ranks_and_caps_fail(self):
        for high, low, cap in [([(0, (0,)), (0, (0,))], [], 2),
                               ([(True, (0,))], [], 1), ([(0, [0])], [], 1),
                               ([(0, (-1,))], [], 1), ([], [], True), ([], [], 65537)]:
            with self.subTest(high=high, cap=cap), self.assertRaises(ValueError):
                MatchIndex(high, low, cap)
        index = MatchIndex([(0, (0,))], [(0, (0,))])
        for rank in [-1, 1, True, Fraction(0)]:
            with self.assertRaises(ValueError):
                index.match(rank)

    def test_literal_reference_matches_grouped_emission_for_multicoordinate_blocks(self):
        table = Tables(6, [0, 1, 5], [1, 0, 5], 4, [(0, 2), (1,)],
                       [[(None, None), (0, 1), (0, 1)], [(None,), (5,)]],
                       [1, 2, 3, 4, 5, 6], [6, 2, 5, 1, 4, 3],
                       [Fraction(2, 3), 1], [1, Fraction(4, 5)])
        index = table.index()
        literal = reference.matches(table)
        self.assertEqual(sorted(index.match(rank) for rank in range(index.total)), literal)
        counts = Counter((h, low) for h, low, _ in literal)
        for high, low in product(table.high_range, table.low_range):
            self.assertEqual(table.multiplicity(high, low), counts[(high, low)])

    def test_zero_entries_and_invalid_acceptance_inputs_are_rejected(self):
        table = Tables(4, [1], [0], 2, [(0,)], [[(0,)]], [1]*4, [1]*4, [2], [2])
        self.assertEqual(table.index().total, 0)
        self.assertEqual(table.distribution(0, [True]*4), [])
        arguments = dict(modulus=4, residues=[1], phases=[0], base=2, blocks=[(0,)],
                         labels=[[(None,)]], f=[1]*4, g=[1]*4, high_thresholds=[1], low_thresholds=[1])
        for mutation in [{'modulus': 5}, {'phases': [4]}, {'base': 3}, {'blocks': [(0,), (0,)]},
                         {'labels': [[(4,)]]}, {'labels': [[]]}, {'f': [1, 1, 0, 1]},
                         {'g': [1.0]*4}, {'high_thresholds': [0]}, {'residues': [True]}]:
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                Tables(**(arguments | mutation))
        with self.assertRaises(ValueError):
            all_star().distribution(0, [1]*4)


class SamplingTests(unittest.TestCase):
    def test_original_phase_and_shift_law_has_unit_joint_score_moment(self):
        for modulus in [2, 3, 4, 6]:
            f = [Fraction(modulus*x, sum(range(1, modulus+1))) for x in range(1, modulus+1)]
            g = list(reversed(f))
            for offset_high, offset_low in product(range(modulus), repeat=2):
                joint = sum((f[(phase+offset_high-shift) % modulus]*g[(offset_low-shift) % modulus]
                             for phase, shift in product(range(modulus), repeat=2)), Fraction(0))/modulus**2
                self.assertEqual(joint, 1)
                for tested_probability in [Fraction(1, 4), Fraction(1, 2), Fraction(1)]:
                    self.assertEqual(1-tested_probability+tested_probability*joint, 1)

    def test_two_coordinate_original_moment_factorizes_without_input_randomness(self):
        for modulus in [2, 3, 4, 6]:
            f = [Fraction(modulus*x, sum(range(1, modulus+1))) for x in range(1, modulus+1)]
            g = list(reversed(f))
            # Draw a star with probability 1/3; every other label is uniform.
            laws = []
            for offset_high, offset_low in [(0, 0), (modulus-1, 1)]:
                law = []
                for phase in range(modulus):
                    law.append((Fraction(1, 3*modulus), Fraction(1)))
                    for shift in range(modulus):
                        value = f[(phase+offset_high-shift) % modulus]*g[(offset_low-shift) % modulus]
                        law.append((Fraction(2, 3*modulus**2), value))
                self.assertEqual(sum((p for p, _ in law), Fraction(0)), 1)
                laws.append(law)
            second = sum((p*q*x*y for p, x in laws[0] for q, y in laws[1]), Fraction(0))
            self.assertEqual(second, 1)

    def test_exact_inverse_multiplicity_second_moment_and_direct_weight(self):
        table = repeated_table()
        index = table.index()
        values = table.distribution(0, [True]*4)
        self.assertEqual([table.multiplicity(k//2, k%2) for k in range(4)], [2, 3, 2, 2])
        self.assertEqual(mean(values, 4), table.direct_checksum(0))
        self.assertEqual(variance(values, 4), Fraction(17, 16))
        direct = sum((value/index.total for value in values), Scalar(4))
        self.assertEqual(direct, Scalar(4, 1))
        second = sum((value.norm2() for value in values), Fraction(0))/len(values)
        formula = Fraction(index.total, 16)*sum((table.fourier_term(k, 0).norm2()/table.multiplicity(k//2, k%2)
                                               for k in range(4)), Fraction(0))
        self.assertEqual(second, formula)

    def test_independent_repeated_draws_divide_complex_variance_by_sample_count(self):
        values = repeated_table().distribution(0, [True]*4)
        center = mean(values, 4)
        for samples in [1, 2, 3]:
            squared = sum(((sum(draw, Scalar(4))/samples-center).norm2()
                           for draw in product(values, repeat=samples)), Fraction(0))/len(values)**samples
            self.assertEqual(squared, variance(values, 4)/samples)

    def test_conditional_second_moment_upper_bound_can_be_sharp(self):
        table = all_star()
        values = table.distribution(2, [True]*4)
        index = table.index()
        bound = Fraction(index.total, 16)*sum((table.fourier_term(k, 2).norm2() for k in range(4)), Fraction(0))
        self.assertEqual(table.direct_checksum(2), Scalar(4))
        self.assertTrue(all(table.multiplicity(k//2, k%2) == 1 for k in range(4)))
        self.assertEqual(variance(values, 4), bound)
        self.assertEqual(bound, 2)

    def test_invalid_frequency_and_wrong_owner_matches_have_zero_value(self):
        table = all_star(3, [1], [0], 2)
        values = table.distribution(0, [True, False, False])
        self.assertEqual(len(values), 4)
        self.assertEqual(sum(value != Scalar(3) for value in values), 1)
        self.assertEqual(mean(values, 3), table.fourier_term(0, 0)/3)
        all_values = table.distribution(0, [True]*3)
        self.assertEqual(mean(all_values, 3), table.direct_checksum(0))
        # Resampling until a valid match changes the law and the estimator.
        # The out-of-range match must consume a draw and contribute zero.
        valid_only = all_values[:3]
        self.assertEqual(mean(valid_only, 3), Scalar(3, Fraction(4, 3)))
        self.assertNotEqual(mean(valid_only, 3), table.direct_checksum(0))
        # Unattained patterns may have raw matches. Waiting for an owned
        # frequency would never finish; a fixed number of zero draws does.
        unattained = table.distribution(0, [False]*3)
        self.assertEqual(len(unattained), 4)
        self.assertTrue(all(value == Scalar(3) for value in unattained))
        # Counting the out-of-range digit pair as frequency zero adds J(0)/P.
        wrong = list(all_values)
        wrong[-1] = table.fourier_term(0, 0)*Fraction(4, 3)
        self.assertNotEqual(mean(wrong, 3), table.direct_checksum(0))

    def test_three_literal_sampling_biases_are_detected(self):
        control = audit.multiplicity_control()
        self.assertEqual(control['correct_mean'], ['1', '0'])
        self.assertEqual(control['omitted_multiplicity_mean'], ['9/4', '1/4'])
        self.assertEqual(control['inconsistent_label_deduplication_mean'], ['15/8', '-1/8'])
        self.assertEqual(control['uniform_group_instead_of_match_mean'], ['15/16', '3/16'])

    def test_fresh_pattern_draws_add_variances_and_reused_rank_can_fail(self):
        control = audit.independent_pattern_control()
        self.assertEqual(control['independent_mean_squared_error'], '3/2')
        self.assertEqual(control['reused_rank_mean_squared_error'], '2')
        table = repeated_table()
        distributions = [table.distribution(0, [k % 3 == j for k in range(4)]) for j in range(3)]
        actual = reference.product_mean_square(distributions, 4)
        self.assertEqual(actual, sum((variance(x, 4) for x in distributions), Fraction(0)))

    def test_different_pattern_table_environments_and_sample_space_sizes(self):
        first = repeated_table().distribution(0, [True, False, False, False])
        second = all_star().distribution(0, [False, True, True, True])
        self.assertEqual((len(first), len(second)), (9, 4))
        self.assertEqual(mean(first, 4)+mean(second, 4), Scalar(4, 1))
        squared = reference.product_mean_square([first, second], 4)
        self.assertEqual(squared, variance(first, 4)+variance(second, 4))

    def test_conditioning_on_coverage_changes_phase_moments_and_failure_can_be_biased(self):
        control = audit.phase_conditioning_control()
        self.assertEqual(control['coverage_probability'], '1/4')
        self.assertEqual(control['unconditional_amplitude_second_moment'], '2')
        self.assertEqual(control['coverage_conditioned_amplitude_second_moment'], '4')
        self.assertEqual(control['unconditional_joint_score_moment'], '1')
        self.assertEqual(control['coverage_conditioned_joint_score_moment'], '9/4')
        self.assertEqual(control['uncovered_nonzero_bias'], [{'phase': 0, 'shift': 1, 'bias': ['-2', '0']}])


class FiniteCircleTests(unittest.TestCase):
    def candidate(self):
        return json.loads((HERE/'FINITE-INPUT.json').read_text())['candidate']

    def test_exact_effective_cutoff_and_its_receipt(self):
        report = finite_circle.report('conorholds')
        self.assertEqual(report, json.loads((HERE/'FINITE.json').read_text()))
        bounds = report['cutoff']
        self.assertEqual(bounds['exponent_amplitude_integer'], 162)
        self.assertEqual(bounds['minimum_modulus_power'], 351)
        self.assertEqual(bounds['sufficient_input_size'], 681)
        self.assertEqual(Fraction(bounds['finite_ceiling']), Fraction('0.51841070'))
        self.assertLessEqual(Fraction(bounds['total_log_error_constant'])/(1 << 351), Fraction('0.00004976'))

    def test_coarser_independent_rational_envelopes_recover_the_same_modulus_power(self):
        candidate = self.candidate()
        design = candidate['design']
        # Recompute deliberately looser envelopes, without the cutoff helpers.
        amp = sum((abs(Fraction(x)) for key in ['r', 'r_sin', 's', 's_sin'] for x in design[key]), Fraction(0))
        rd = sum((j*abs(Fraction(x)) for key in ['r', 'r_sin'] for j, x in enumerate(design[key], 1)), Fraction(0))
        total_derivative = rd+sum((j*abs(Fraction(x)) for key in ['s', 's_sin']
                                   for j, x in enumerate(design[key], 1)), Fraction(0))
        self.assertLess(amp, Fraction('8.47'))
        self.assertLess(rd, Fraction('13.92'))
        self.assertLess(total_derivative, Fraction('28.09'))
        self.assertTrue(1 <= Fraction(design['kappa']) <= Fraction('1.035'))
        p = {key: Fraction(value) for key, value in candidate['parameters'].items()}
        a = (p['eta']-p['tau'])/(2*p['tau']-p['eta'])
        self.assertLess(a, Fraction('0.0631'))
        self.assertLess(p['eta']-p['tau'], Fraction('0.0291'))
        self.assertLess(Fraction(candidate['multipliers'][0]), Fraction('8.389'))
        self.assertLess(Fraction(candidate['multipliers'][1]), Fraction('8.199'))
        b, kappa, lam1, lam2, a = map(Fraction, ['8.47', '1.035', '8.389', '8.199', '0.0631'])
        upper_u, upper_v = 2*kappa*b, 4*b
        self.assertLess(lam1*(Fraction('0.0291')+upper_u), 162)
        self.assertLess(lam2*(a*upper_v+upper_u), 162)
        td = Fraction('13.92')*(1+2*b)
        ud, vd = Fraction('0.035')*td+Fraction('13.92'), td+2*Fraction('13.92')
        self.assertLess(max(lam1*ud, lam2*(a*vd+ud)), 330)
        ue = Fraction('0.035')*(1+2*b)+1+kappa
        ve = (1+2*b)+2+kappa
        self.assertLess(8*Fraction('28.09')*max(lam1*ue, lam2*(a*ve+ue)), 7350)
        independent_constant = 7350+2644*(1 << 325)
        self.assertLess(Fraction(independent_constant, 1 << 351), Fraction('0.00004976'))
        self.assertGreater(Fraction(independent_constant, 1 << 350), Fraction('0.00004976'))
        self.assertGreaterEqual(p['eta']*681-2, 351)
        self.assertLess(p['eta']*680-2, 351)

    def test_source_snapshot_and_finite_normalizer_coefficients_are_bound(self):
        finite_input = json.loads((HERE/'FINITE-INPUT.json').read_text())
        inputs = json.loads((HERE/'INPUTS.json').read_text())
        prior = next(x for x in inputs['parameter_sources'] if x['path'] == finite_input['path'])
        self.assertEqual(finite_input['sha256'], prior['sha256'])
        self.assertEqual(finite_input['commit'], prior['commit'])
        snapshot = inputs['parameter_witnesses'][-1]
        candidate = finite_input['candidate']
        self.assertEqual(snapshot['multipliers'], candidate['multipliers'])
        self.assertEqual(snapshot['parameters'], {key: candidate['parameters'][key] for key in snapshot['parameters']})
        # Separate floating control of the raw-normalizer cancellation, not
        # part of the exact cutoff certificate or a quadrature proof.
        design, kappa = candidate['design'], float(candidate['design']['kappa'])
        for modulus in [8, 16]:
            values = {}
            for key in ['r', 's']:
                values[key] = [sum(float(c)*cos(j*2*pi*k/modulus)+float(d)*sin(j*2*pi*k/modulus)
                                   for j, (c, d) in enumerate(zip(design[key], design[key+'_sin']), 1))
                               for k in range(modulus)]
            rawf, rawg = [[2**x for x in values[key]] for key in ['r', 's']]
            f = [x/(sum(rawf)/modulus) for x in rawf]
            g = [x/(sum(rawg)/modulus) for x in rawg]
            fs = [x**kappa/(sum(y**kappa for y in f)/modulus) for x in f]
            gs = [x**kappa/(sum(y**kappa for y in g)/modulus) for x in g]
            z = log2(sum(2**(kappa*x) for x in values['r'])/modulus)+log2(sum(2**(kappa*x) for x in values['s'])/modulus)
            for difference in range(modulus):
                weights = [f[(difference+w) % modulus]*g[w] for w in range(modulus)]
                kernel = [x/(sum(weights)/modulus) for x in weights]
                d = sum(x*log2(x) for x in kernel)/modulus
                a = sum(kernel[w]*log2(fs[(difference+w) % modulus]) for w in range(modulus))/modulus
                e = sum(kernel[w]*log2(gs[w]) for w in range(modulus))/modulus
                t = sum(kernel[w]*(values['r'][(difference+w) % modulus]+values['s'][w]) for w in range(modulus))/modulus
                l = log2(sum(rawf[(difference+w) % modulus]*rawg[w] for w in range(modulus))/modulus)
                self.assertAlmostEqual(a+e-d, (kappa-1)*t+l-z, places=12)
                self.assertAlmostEqual(2*d-a-e, (2-kappa)*t-2*l+z, places=12)

    def test_reflection_and_smaller_remaining_reserve(self):
        candidate = self.candidate()
        original = finite_circle.cutoff(candidate)
        reflected = deepcopy(candidate)
        # Preserve exact finite decimal strings when reflecting coefficients.
        for key in ['r_sin', 's_sin']:
            reflected['design'][key] = [x[1:] if x.startswith('-') else '-'+x for x in candidate['design'][key]]
        self.assertEqual(finite_circle.cutoff(reflected), original)
        tighter = deepcopy(candidate)
        tighter['parameters']['moment_ceiling'] = '0.51841070'
        smaller = finite_circle.cutoff(tighter)
        self.assertGreater(smaller['minimum_modulus_power'], original['minimum_modulus_power'])

    def test_invalid_or_unbounded_cutoff_inputs_are_rejected(self):
        candidate = self.candidate()
        for field, value in [('kappa', '0'), ('kappa', '5'), ('kappa', '1e1000000000'),
                             ('r_sin', []), ('r', [1.0]*20), ('r', ['100']*20)]:
            bad = deepcopy(candidate)
            bad['design'][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                finite_circle.cutoff(bad)
        for changed in [{'eta': '0.4'}, {'moment_ceiling': '0.51846094'}, {'tau': '0.34'}]:
            bad = deepcopy(candidate)
            bad['parameters'].update(changed)
            with self.assertRaises(ValueError):
                finite_circle.cutoff(bad)
        bad = deepcopy(candidate)
        bad['multipliers'] = ['1000000', '8']
        with self.assertRaises(ValueError):
            finite_circle.cutoff(bad)
        bad = deepcopy(candidate)
        bad['parameters'].update(eta='1000000', tau='800000', moment_ceiling='999999')
        with self.assertRaises(ValueError):
            finite_circle.cutoff(bad)
        with self.assertRaises(ValueError):
            finite_circle.report('unapproved')


class BudgetAndEvidenceTests(unittest.TestCase):
    def test_uniform_multiplicity_binomial_tail_against_exact_probabilities(self):
        probabilities = [Fraction(1, 16), Fraction(1, 8), Fraction(1, 4), Fraction(1, 2),
                         Fraction(3, 4), Fraction(7, 8), Fraction(15, 16)]
        for rows, probability in product(range(1, 65), probabilities):
            expectation = rows*probability
            tail = sum((Fraction(comb(rows, accepted))*probability**accepted*(1-probability)**(rows-accepted)
                        for accepted in range(rows+1) if accepted < expectation/2), Fraction(0))
            # Raise to an integer power so the exponential bound is checked
            # with exact integers, without a floating logarithm or root.
            power = 8*probability.denominator
            self.assertLessEqual(tail.numerator**power*(1 << (rows*probability.numerator)),
                                 tail.denominator**power)

    def test_all_sixteen_snapshots_have_positive_exact_variance_reserve(self):
        rows = audit.budget_profiles()
        self.assertEqual(len(rows), 16)
        self.assertTrue(all(Fraction(row['variance_decay']) > 0 for row in rows))
        best = rows[-1]
        self.assertEqual(Fraction(best['variance_decay']), Fraction('0.00004976'))
        self.assertEqual(Fraction(best['finite_ceiling']), Fraction('0.51841070'))
        self.assertTrue(all(row['pattern_factor_in_mse'] == 1 and row['prefactor_upper_bound'] == 8 for row in rows))

    def test_reduced_sample_exponent_preserves_a_positive_exact_error_reserve(self):
        inputs = json.loads((HERE/'INPUTS.json').read_text())
        for row in inputs['parameter_witnesses']:
            p = {key: Fraction(value) for key, value in row['parameters'].items()}
            reserve = p['eta']-p['moment_ceiling']-48*p['epsilon']
            finite_ceiling = p['moment_ceiling']+reserve/2
            sigma = p['tau']-reserve/4
            self.assertGreater(sigma, 0)
            self.assertLess(sigma, p['tau'])
            self.assertEqual(sigma+(p['eta']-p['tau'])-finite_ceiling-48*p['epsilon'], reserve/4)
            if row['id'] == 'asymmetric-filters/asymmetric_degree20_retuned':
                self.assertEqual(sigma, Fraction('0.48943092'))
                self.assertEqual(reserve/4, Fraction('0.00002488'))

    def test_invalid_checksum_domain_rounding_and_zero_reserve_fail(self):
        base = json.loads((HERE/'INPUTS.json').read_text())['parameter_witnesses'][-1]
        for changed in [{'eta': '0.4'}, {'eta': '0.99'}, {'tau': '0'}, {'epsilon': '0'},
                        {'moment_ceiling': base['parameters']['eta']}, {'epsilon': '0.001'},
                        {'tau': 0.4894558}]:
            bad = deepcopy(base)
            bad['parameters'].update(changed)
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                audit.budget_profile(bad)
        for slopes in [['0.99', '8'], ['8', '0.99'], ['10.0001', '8'], ['8', '10'], ['8']]:
            bad = deepcopy(base)
            bad['multipliers'] = slopes
            with self.subTest(slopes=slopes), self.assertRaises(ValueError):
                audit.budget_profile(bad)

    def test_frozen_mathematical_source_matches_the_existing_manifest(self):
        inputs = json.loads((HERE/'INPUTS.json').read_text())
        manifest = HERE.parent/'filter-budgets/SOURCE.json'
        self.assertEqual(inputs['mathematical_source'], json.loads(manifest.read_text()))
        self.assertEqual(inputs['source_manifest_sha256'], hashlib.sha256(manifest.read_bytes()).hexdigest())
        sources = {row['path'].removeprefix('research/'): row for row in inputs['parameter_sources']}
        self.assertEqual(len(sources), 4)
        self.assertTrue(all(row['source'] in sources for row in inputs['parameter_witnesses']))

    def test_retained_receipts_bind_current_code_and_exact_controls(self):
        receipt = json.loads((HERE/'AUDIT.json').read_text())
        self.assertEqual(receipt['status'], 'PASS')
        self.assertEqual(receipt['author'], 'conorholds')
        self.assertEqual(receipt['inputs_sha256'], audit.input_hashes())
        for name in ['multiplicity_control', 'independent_pattern_control', 'phase_conditioning_control']:
            self.assertEqual(receipt[name], getattr(audit, name)())
        self.assertEqual(receipt['budget_profiles'], audit.budget_profiles())
        self.assertEqual(receipt['finite_campaign']['configurations'], 163120)
        self.assertEqual(receipt['finite_campaign']['fourier_target_controls'], 655120)

    def test_standard_campaign_replays_the_retained_exhaustive_controls(self):
        receipt = json.loads((HERE/'STANDARD.json').read_text())
        self.assertEqual(receipt['inputs_sha256'], audit.input_hashes())
        self.assertEqual(audit.finite_campaign(), receipt['finite_campaign'])
        self.assertEqual(receipt['finite_campaign']['configurations'], 3120)
        with self.assertRaises(ValueError):
            audit.finite_campaign(1)
        with self.assertRaises(ValueError):
            audit.run('unapproved')


if __name__ == '__main__':
    unittest.main(verbosity=2)
