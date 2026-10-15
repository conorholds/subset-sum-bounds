"""Controls for the block-count proof, discrete search, and retained reserves."""
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block_bound
import frontier
from certify import budget_margins


class BlockBoundTests(unittest.TestCase):
    def test_exact_entropy_bounds_and_elementary_integer_controls(self):
        lower, upper = block_bound.entropy_bounds()
        self.assertLess(3**75, 2**119)
        self.assertGreater(3**12, 2**19)
        self.assertGreater(lower, Fraction('0.81'))
        self.assertLess(upper, Fraction('0.8125'))
        self.assertLess(upper-lower, Fraction(1, 10**23))
        coarse = block_bound.logarithm_bounds(Fraction(2), 4)
        fine = block_bound.logarithm_bounds(Fraction(2), 24)
        self.assertLess(coarse[0], fine[0])
        self.assertLess(fine[1], coarse[1])

    def test_eighteen_blocks_are_excluded_but_seventeen_are_not(self):
        self.assertGreater(block_bound.alpha_floor(18)[0], Fraction('0.49'))
        self.assertLess(block_bound.alpha_floor(17)[1], Fraction('0.4897358'))
        self.assertLess(block_bound.alpha_floor(18)[1], block_bound.alpha_floor(19)[0])
        self.assertEqual(block_bound.certificate(), json.loads((HERE/'BLOCK-BOUND.json').read_text()))

    def test_invalid_block_counts_and_logarithm_budgets_fail(self):
        for blocks in [0, -1, True, 17.5]:
            with self.subTest(blocks=blocks), self.assertRaises(ValueError):
                block_bound.alpha_floor(blocks)
        for argument, terms in [(Fraction(1), 24), (Fraction(3), 24),
                                (Fraction(2), 0), (Fraction(2), 65), (2.0, 24)]:
            with self.subTest(argument=argument, terms=terms), self.assertRaises(ValueError):
                block_bound.logarithm_bounds(argument, terms)


class FrontierTests(unittest.TestCase):
    def test_balanced_density_and_clipped_boundary(self):
        rho = frontier.balanced_rho(.51798, 17)
        self.assertAlmostEqual(.5+(frontier.ENTROPY-.5)*rho+.00001+.00001-.000016, .51798)
        self.assertAlmostEqual(frontier.balanced_rho(.55, 17), .999/17)
        self.assertIsNone(frontier.screen(.500004, 17, lambda eta: .49))
        self.assertIsNone(frontier.screen(.51798, 17, lambda eta: float('inf')))

    def test_opposite_density_perturbations_fail_independent_budget_checks(self):
        design = next(w for w in frontier.load_witnesses() if w['id'] == 'combined')['design']
        balanced = frontier.screen(.5179812, 17, lambda eta: .48)
        good = frontier.witness(balanced, design)
        budget_margins(good['parameters'])
        for change, failure in [(-.001, 'filtered_lists'), (.001, 'visited_pairs')]:
            row = dict(balanced, rho=balanced['rho']+change)
            parameters = frontier.witness(row, design)['parameters']
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, failure):
                budget_margins(parameters)

    def test_retained_witnesses_replay_rounding_and_preserve_key_reserves(self):
        retained = json.loads((HERE/'witnesses.json').read_text())
        discovery = json.loads((HERE/'discovery.json').read_text())
        self.assertEqual(frontier.retained_witnesses(discovery), retained)
        baseline = next(w for w in frontier.load_witnesses() if w['id'] == 'combined')
        baseline_margins = budget_margins(baseline['parameters'])
        for witness in retained['witnesses']:
            self.assertEqual(witness['design'], baseline['design'])
            p = witness['parameters']
            margin = budget_margins(p)
            self.assertLess(Decimal(p['alpha']), Decimal(baseline['parameters']['alpha']))
            self.assertEqual(Decimal(p['epsilon']), Decimal('.00000001'))
            self.assertEqual(Decimal(p['guard']), Decimal('.00001'))
            self.assertEqual(Decimal(p['eta'])-Decimal(p['moment_ceiling']), Decimal('.0001'))
            self.assertEqual(Decimal('.5')-Decimal(p['local_deficiency'])-Decimal(p['gamma']), Decimal('.00003'))
            self.assertEqual(Decimal(p['preliminary'])-(Decimal('.5')-Decimal(p['deficiency'])), Decimal('.00001'))
            for name in ['coverage', 'filtered_lists', 'visited_pairs']:
                self.assertGreater(margin[name].a, baseline_margins[name].b)

    def test_each_block_count_is_refined_reproducibly_with_a_fixed_budget(self):
        forward = frontier.search(blocks=(16, 17), coarse_points=41, refinement_points=21)
        reverse = frontier.search(blocks=(17, 16), coarse_points=41, refinement_points=21)
        self.assertEqual(sorted(forward['frontier'], key=lambda row: row['blocks']),
                         sorted(reverse['frontier'], key=lambda row: row['blocks']))
        self.assertEqual(forward['parameter_screens'], 2*(41+2*21))
        self.assertLessEqual(forward['moment_evaluations'], (41+2*2*21)*41)
        self.assertLess(forward['frontier'][1]['tau'], forward['frontier'][0]['tau'])
        self.assertTrue(all(.510 <= row['eta'] <= .550 for row in forward['frontier']))

    def test_invalid_search_ranges_cannot_expand_the_campaign(self):
        for arguments in [{'blocks': (18,)}, {'blocks': (17, 17)}, {'blocks': (True,)},
                          {'coarse_points': 402}, {'refinement_points': 202}]:
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                frontier.search(**arguments)
        with self.assertRaises(ValueError):
            frontier.balanced_rho(float('nan'), 17)


if __name__ == '__main__':
    unittest.main(verbosity=2)
