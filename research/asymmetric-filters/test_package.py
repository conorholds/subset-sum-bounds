"""Orientation, continuous-bound, derivative, and changed-model controls."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import ceil
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


spectral = load('asymmetric_spectral_controls', HERE/'spectral.py')
check = load('asymmetric_interval_controls', HERE/'interval_check.py')
search = load('asymmetric_search_controls', HERE/'search.py')
density = load('asymmetric_density_controls', HERE/'density_probe.py')
old = load('asymmetric_old_interval_control', HERE.parent/'chernoff-multipliers/interval_check.py')


def candidates():
    return json.loads((HERE/'candidates.json').read_text())['candidates']


def rotate(design, phase, both=True):
    result = deepcopy(design)
    for key in (['r', 's'] if both else ['s']):
        c, d = np.array(list(map(float, design[key]))), np.array(list(map(float, design[key+'_sin'])))
        x = np.arange(1, len(c)+1)*phase
        result[key] = [f'{v:.12f}' for v in c*np.cos(x)+d*np.sin(x)]
        result[key+'_sin'] = [f'{v:.12f}' for v in -c*np.sin(x)+d*np.cos(x)]
    return result


def direct_scores(design, points=256):
    x = np.arange(points)*2*np.pi/points
    angles = np.outer(np.arange(1, len(design['r'])+1), x)
    values = []
    for key in ['r', 's']:
        values.append(np.array(list(map(float, design[key])))@np.cos(angles)
                      +np.array(list(map(float, design[key+'_sin'])))@np.sin(angles))
    r, s = values
    shifts = (np.arange(points)[:, None]+np.arange(points)[None, :]) % points
    raw = np.exp2(r[shifts]+s[None, :])
    conditional = raw/raw.sum(axis=1)[:, None]
    kappa = float(design['kappa'])
    test_f, test_g = np.exp2(kappa*r), np.exp2(kappa*s)
    test_f /= test_f.mean(); test_g /= test_g.mean()
    d = (conditional*np.log2(points*conditional)).sum(axis=1)
    a = (conditional*np.log2(test_f[shifts])).sum(axis=1)
    e = (conditional*np.log2(test_g[None, :])).sum(axis=1)
    return a+e-d, 2*d-a-e, d, a, e


class FourierTests(unittest.TestCase):
    def test_zero_sine_embedding_recovers_the_previous_continuous_model(self):
        base = candidates()[0]
        prior = json.loads((HERE/'predecessor.json').read_text())['candidate']
        new = check.certify_candidate(base)
        previous = old.enclose_moments(prior['design'], prior['parameters'], prior['multipliers'])
        for a, b in zip(new['moments']['log2_moment_intervals'], previous['log2_moment_intervals']):
            self.assertLessEqual(Decimal(a[0]), Decimal(b[1]))
            self.assertLessEqual(Decimal(b[0]), Decimal(a[1]))

    def test_fft_orientation_and_entropy_scores_match_an_independent_kernel(self):
        for candidate in [candidates()[1], candidates()[-1]]:
            u, v, _ = spectral.scores(candidate['design'], 256)
            direct_u, direct_v, d, a, e = direct_scores(candidate['design'])
            np.testing.assert_allclose(u, direct_u, atol=2e-13, rtol=0)
            np.testing.assert_allclose(v, direct_v, atol=2e-13, rtol=0)
            self.assertGreaterEqual(min(d.min(), (d-a).min(), (d-e).min()), -1e-13)
            self.assertGreater(np.max(np.abs(u-u[(-np.arange(256)) % 256])), .001)

    def test_common_rotation_preserves_scores_but_relative_rotation_changes_them(self):
        design = candidates()[-1]['design']
        u, v, _ = spectral.scores(design, 512)
        ru, rv, _ = spectral.scores(rotate(design, .713), 512)
        np.testing.assert_allclose(u, ru, atol=2e-11, rtol=0)
        np.testing.assert_allclose(v, rv, atol=2e-11, rtol=0)
        wrong, _, _ = spectral.scores(rotate(design, .713, both=False), 512)
        self.assertGreater(np.max(abs(u-wrong)), .001)

    def test_common_reflection_preserves_both_integrated_moments(self):
        candidate = candidates()[-1]
        reflected = deepcopy(candidate['design'])
        for key in ['r_sin', 's_sin']:
            reflected[key] = [str(-Decimal(v)) for v in reflected[key]]
        u, v, m2 = spectral.scores(candidate['design'], 512)
        ru, rv, _ = spectral.scores(reflected, 512)
        reverse = (-np.arange(512)) % 512
        np.testing.assert_allclose(u[reverse], ru, atol=1e-13, rtol=0)
        np.testing.assert_allclose(v[reverse], rv, atol=1e-13, rtol=0)

    def test_exact_rotation_invariant_proves_the_filter_is_not_shifted_even(self):
        design = candidates()[-1]['design']
        c1, d1, c2, d2 = [Fraction(v) for v in [design['r'][0], design['r_sin'][0], design['r'][1], design['r_sin'][1]]]
        invariant = 2*c2*c1*d1-d2*(c1*c1-d1*d1)
        self.assertNotEqual(invariant, 0)
        self.assertGreater(abs(invariant), Fraction(1, 5))

    def test_analytic_moment_and_eta_derivatives_match_finite_differences(self):
        x = np.concatenate((spectral.vector_from_design(candidates()[-1]['design']),
                            np.array(list(map(float, candidates()[-1]['multipliers'])))))
        for points in [128, 256]:
            state = spectral.score_state(x[:-2], points)
            eta = float(candidates()[-1]['parameters']['eta'])
            _, jacobian, slopes = search.at_eta(state, x[-2:], eta)
            for j in range(len(x)):
                lo, hi = x.copy(), x.copy()
                lo[j] -= 1e-5; hi[j] += 1e-5
                a = search.at_eta(spectral.score_state(lo[:-2], points), lo[-2:], eta)[0]
                b = search.at_eta(spectral.score_state(hi[:-2], points), hi[-2:], eta)[0]
                np.testing.assert_allclose((b-a)/2e-5, jacobian[:, j], atol=1e-7, rtol=0)
            a = search.at_eta(state, x[-2:], eta-1e-6)[0]
            b = search.at_eta(state, x[-2:], eta+1e-6)[0]
            np.testing.assert_allclose((b-a)/2e-6, slopes, atol=1e-7, rtol=0)

    def test_implicit_frontier_derivatives_and_common_rotation_null_direction(self):
        candidate = candidates()[-1]
        x = np.concatenate((spectral.vector_from_design(candidate['design']),
                            np.array(list(map(float, candidate['multipliers'])))))
        _, jacobian = search.frontier_and_jacobian(x)
        for j in range(len(x)):
            lo, hi = x.copy(), x.copy(); lo[j] -= 1e-5; hi[j] += 1e-5
            a = search.frontier_and_jacobian(lo)[0]
            b = search.frontier_and_jacobian(hi)[0]
            np.testing.assert_allclose((b-a)/2e-5, jacobian[:, j], atol=1e-7, rtol=0)
        n = len(candidate['design']['r'])
        frequency = np.arange(1, n+1)
        gauge = np.concatenate((frequency*x[n:2*n], -frequency*x[:n],
                                frequency*x[3*n:4*n], -frequency*x[2*n:3*n], [0., 0., 0.]))
        np.testing.assert_allclose(jacobian@gauge, 0, atol=1e-12, rtol=0)

    def test_shifted_squared_norm_formula_and_zero_free_disk_without_evenness(self):
        candidate = candidates()[-1]
        row = check.certify_candidate(candidate)
        design = candidate['design']
        points, beta = 512, .12
        x = np.arange(points)*2*np.pi/points
        j = np.arange(1, len(design['r'])+1)[:, None]
        values = []
        for key in ['r', 's']:
            c, d = np.array(list(map(float, design[key])))[:, None], np.array(list(map(float, design[key+'_sin'])))[:, None]
            real = (c*np.cos(j*x)+d*np.sin(j*x)).sum(axis=0)
            shifted = (c*np.cos(j*(x+1j*beta))+d*np.sin(j*(x+1j*beta))).sum(axis=0)
            doubled = (2*np.cosh(j*beta)*(c*np.cos(j*x)+d*np.sin(j*x))).sum(axis=0)
            np.testing.assert_allclose(abs(np.exp2(shifted))**2, np.exp2(doubled), atol=1e-12, rtol=1e-12)
            values.append((c, d, np.exp2(real).mean()))
        observed = 0.
        for difference in x[::8]:
            factors = []
            for (c, d, z), angle in zip(values, [x+difference+1j*beta, x-1j*beta]):
                poly = (c*np.cos(j*angle)+d*np.sin(j*angle)).sum(axis=0)
                factors.append(np.exp2(poly)/z)
            observed = max(observed, abs(np.mean(factors[0]*factors[1])-1))
        self.assertLess(observed, float(row['moments']['zero_free_theta_upper']))


class CertificateTests(unittest.TestCase):
    def test_all_four_candidates_replay_and_match_finer_independent_grids(self):
        receipt = json.loads((HERE/'CERTIFICATE.json').read_text())
        points, digits = check.validate_receipt(receipt)
        fresh = check.replay(points, digits)
        self.assertEqual(len(fresh['results']), 4)
        for candidate, row in zip(candidates(), fresh['results']):
            eta, tau = [float(candidate['parameters'][key]) for key in ['eta', 'tau']]
            for grid in [256, 512, 1024]:
                u, v, m2 = spectral.scores(candidate['design'], grid)
                for direction, lam, interval in zip([eta-tau-u, (eta-tau)/(2*tau-eta)*v-u],
                                                    candidate['multipliers'], row['moments']['log2_moment_intervals']):
                    value = search.tune.log_moment(direction, m2, float(lam))
                    self.assertLess(float(interval[0]), value)
                    self.assertLess(value, float(interval[1]))

    def test_fractional_coverage_and_all_prescribed_reserves_are_present(self):
        candidate = candidates()[-1]
        p = {key: Fraction(value) for key, value in candidate['parameters'].items()}
        self.assertGreater(17*p['rho'], 1)
        self.assertGreaterEqual(p['rho']*(p['local_deficiency']-p['deficiency']), Fraction(1, 625000))
        self.assertEqual(p['eta']-p['moment_ceiling'], Fraction(1, 10000))
        self.assertEqual(p['alpha']-p['tau'], Fraction(1, 100000))
        self.assertGreaterEqual(p['deficiency']-(Fraction(1, 2)-p['tau']), Fraction(1, 50000))
        maximum = ceil(1/p['rho'])-1
        self.assertEqual(maximum, 16)
        for count in range(1, maximum+2):
            parameters = dict(candidate['parameters'], blocks=count)
            with self.assertRaisesRegex(ValueError, 'coverage|block_space'):
                old.budget_margins(parameters)

    def test_legacy_schema_and_invalid_exact_designs_fail_closed(self):
        original = candidates()[-1]
        bad_designs = []
        bad = deepcopy(original['design']); bad.pop('r_sin'); bad_designs.append(bad)
        bad = deepcopy(original['design']); bad['s_sin'].pop(); bad_designs.append(bad)
        bad = deepcopy(original['design']); bad['r_sin'][0] = .1; bad_designs.append(bad)
        bad = deepcopy(original['design']); bad['s_sin'][0] = 'NaN'; bad_designs.append(bad)
        bad = deepcopy(original['design']); bad['kappa'] = '0'; bad_designs.append(bad)
        for design in bad_designs:
            with self.assertRaises(ValueError):
                check.certify_candidate(dict(original, design=design))
        with self.assertRaises(ValueError):
            old.enclose_moments(original['design'], original['parameters'], original['multipliers'])

    def test_invalid_multipliers_inadequate_quadrature_and_uniform_filters_fail(self):
        original = candidates()[-1]
        for multipliers in [['.99', '8.2'], ['10', '8.2'], ['8.4', '10'], ['NaN', '8.2']]:
            with self.assertRaises(ValueError):
                check.certify_candidate(dict(original, multipliers=multipliers))
        with self.assertRaisesRegex(ValueError, 'continuous moment'):
            check.certify_candidate(original, points=64)
        design = {'r': ['0'], 's': ['0'], 'r_sin': ['0'], 's_sin': ['0'], 'kappa': '1'}
        with self.assertRaisesRegex(ValueError, 'continuous moment'):
            check.certify_candidate(dict(original, design=design))

    def test_receipt_rejects_missing_rows_wrong_model_status_and_stale_dependencies(self):
        original = json.loads((HERE/'CERTIFICATE.json').read_text())
        mutations = []
        bad = deepcopy(original); bad['results'].pop(); mutations.append(bad)
        bad = deepcopy(original); bad['results'][1] = deepcopy(bad['results'][0]); mutations.append(bad)
        bad = deepcopy(original); bad['model'] = 'fractional-coverage-v1'; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['status'] = 'FAIL'; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['multipliers'] = ['1', '1']; mutations.append(bad)
        bad = deepcopy(original); bad['results'][0]['moments']['grid'] = True; mutations.append(bad)
        bad = deepcopy(original); bad['inputs_sha256']['fractional']['budget.py'] = '0'*64; mutations.append(bad)
        for receipt in mutations:
            with self.assertRaises(ValueError):
                check.validate_receipt(receipt)


class SearchTests(unittest.TestCase):
    def test_density_projection_and_predecessor_bindings_match_exact_inputs(self):
        discovery = json.loads((HERE/'density-discovery.json').read_text())
        start = json.loads((HERE/'start.json').read_text())
        self.assertEqual(density.project_to_fourier(discovery['selected']), start)
        self.assertEqual(discovery['search_sha256'], hashlib.sha256((HERE/'density_probe.py').read_bytes()).hexdigest())
        snapshot = json.loads((HERE/'predecessor.json').read_text())
        self.assertEqual(snapshot['sha256'], hashlib.sha256((HERE.parent/'fractional-coverage/candidates.json').read_bytes()).hexdigest())
        self.assertEqual(snapshot['candidate'], json.loads((HERE.parent/'fractional-coverage/candidates.json').read_text())['candidates'][1])

    def test_uniform_density_fixed_point_is_a_failure_control_not_an_optimum_claim(self):
        row = density.run_one('uniform_control', np.ones(density.N), np.ones(density.N), .4245, iterations=2)
        self.assertEqual(row['accepted'], 0)
        self.assertGreater(max(row['log_moments']), 1)

    def test_density_update_decreases_the_scalar_objective(self):
        h, g = density.initializations()['split_and_single_peaks']
        state = density.optimize_kernel(h, g, .4245)
        before = density.objective(state, .4245)
        row = density.run_one('bounded_control', h, g, .4245, iterations=2)
        self.assertLess(row['objective'], before)
        self.assertLessEqual(row['largest_kernel_update_increase'], 1e-12)

    def test_search_limits_and_invalid_fourier_inputs_fail(self):
        for kwargs in [{'iterations': 0}, {'max_evaluations': 1}, {'degree': 21}, {'author': ''}]:
            with self.assertRaises(ValueError):
                search.optimize(**kwargs)
        for vector in [[], [1, 2, 3], [1, 2, 3, 4, 0], [float('nan')]*5]:
            with self.assertRaises(ValueError):
                spectral.score_state(vector)
        with self.assertRaises(ValueError):
            density.density(np.array([0., 1.]))
        with self.assertRaises(ValueError):
            density.project_to_fourier({}, degree=64)
        report = search.optimize(iterations=1, max_evaluations=2)
        self.assertLessEqual(report['evaluations'], 2)

    def test_retained_search_receipt_has_fresh_code_and_start_bindings(self):
        report = json.loads((HERE/'joint-search.json').read_text())
        self.assertEqual(report['candidate'], candidates()[-1])
        self.assertEqual(report['inputs_sha256']['search.py'], hashlib.sha256((HERE/'search.py').read_bytes()).hexdigest())
        self.assertEqual(report['inputs_sha256']['spectral.py'], hashlib.sha256((HERE/'spectral.py').read_bytes()).hexdigest())
        self.assertEqual(report['start']['sha256'], hashlib.sha256((HERE/'start.json').read_bytes()).hexdigest())
        self.assertLessEqual(report['evaluations'], report['limits']['evaluations'])


if __name__ == '__main__':
    unittest.main()
