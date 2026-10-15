"""Controls for the pilot discovery replay; certificates remain a separate check."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('initial_retuning_replay', ROOT/'research/initial-filter-retuning/replay.py')
replay = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(replay)


class InitialRetuningTests(unittest.TestCase):
    def test_source_screen_and_infeasible_controls(self):
        source = next(w for w in replay.load_witnesses() if w['id'] == 'source')['design']
        sc = replay.scores(source)
        result = replay.parameter_screen(sc, .51786, .05735, 17)
        self.assertIsNotNone(result)
        self.assertLess(result['tau'], .49)
        self.assertLess(max(result['moment_logs']), result['eta']-.0001)
        self.assertIsNone(replay.parameter_screen(sc, .51786, 1/16, 16))
        uniform = replay.scores({'r': ['0']*6, 's': ['0']*6, 'kappa': '1'})
        self.assertIsNone(replay.parameter_screen(uniform, .51786, .05735, 17))

    def test_coordinate_search_is_reproducible_and_bounded(self):
        source = next(w for w in replay.load_witnesses() if w['id'] == 'source')['design']
        r, s = list(map(float, source['r'])), list(map(float, source['s']))
        before = max(replay.moments(replay.filter_scores(r, s, float(source['kappa'])), .51772, .48996))
        a = replay.coordinate_search(r, s, float(source['kappa']), 2)
        b = replay.coordinate_search(r, s, float(source['kappa']), 2)
        self.assertEqual(a, b)
        self.assertLessEqual(a['objective'], before)
        self.assertLessEqual(a['evaluations'], 1+2*13*2)
        with self.assertRaises(ValueError):
            replay.coordinate_search(r, s, float(source['kappa']), 66)

    def test_simplex_solves_a_known_quadratic_with_bounded_calls(self):
        calls = 0

        def objective(x):
            nonlocal calls
            calls += 1
            return float(np.sum((x-np.array([1., -2.]))**2))

        point, value, counted = replay.nelder_mead(objective, [2., -3.], 200)
        self.assertLess(np.linalg.norm(point-np.array([1., -2.])), 1e-8)
        self.assertLess(value, 1e-16)
        self.assertEqual(calls, counted)
        self.assertLessEqual(counted, 3+200*4)
        for start, iterations in [([], 1), ([float('nan')], 1), ([1.], -1), ([1.], 3501)]:
            with self.subTest(start=start, iterations=iterations), self.assertRaises(ValueError):
                replay.nelder_mead(objective, start, iterations)


if __name__ == '__main__':
    unittest.main()
