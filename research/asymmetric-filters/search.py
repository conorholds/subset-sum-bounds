"""Exploratory cosine/sine filter search on the fractional-coverage frontier.

The analytic score derivatives adapt minimax-filters/optimize.py from PR 13.
This experiment changes the objective to the two implicit eta boundaries.
"""
import json
import hashlib
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'chernoff-multipliers'))
import multiplier_search as tune

from spectral import score_state, vector_from_design
import importlib.util
SPEC = importlib.util.spec_from_file_location('asymmetric_fractional_search', HERE.parent/'fractional-coverage/search.py')
FRACTIONAL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FRACTIONAL)


def at_eta(state, multipliers, eta):
    row = FRACTIONAL.slice_parameters(eta)
    tau = row['tau']
    tau_prime = (tune.H-1+.0000324)/((tune.H-.5)*(1+row['rho'])**2)
    delta, gamma = eta - tau, 2 * tau - eta
    ratio = delta / gamma
    ratio_prime = ((1 - tau_prime) * gamma - delta * (2 * tau_prime - 1)) / gamma**2
    directions = [delta - state['u'], ratio * state['v'] - state['u']]
    logs, jacobian, eta_slopes = [], [], []
    for branch in range(2):
        lam = multipliers[branch]
        exponent = lam * directions[branch]
        shift = np.max(exponent)
        weights = state['m2'] * np.exp2(exponent - shift)
        total = weights.sum()
        weights /= total
        logs.append(shift + np.log2(total / len(weights)))
        coefficient_gradient = (-lam * state['du'] if branch == 0
                                else lam * (ratio * state['dv'] - state['du'])) @ weights
        multiplier_gradient = np.zeros(2)
        multiplier_gradient[branch] = directions[branch] @ weights
        jacobian.append(np.concatenate((coefficient_gradient, multiplier_gradient)))
        eta_slopes.append(lam * (1 - tau_prime) if branch == 0
                          else lam * ratio_prime * (state['v'] @ weights))
    return np.array(logs), np.array(jacobian), np.array(eta_slopes)


def frontier_and_jacobian(vector, points=128, bisections=40):
    vector = np.asarray(vector, dtype=float)
    if (vector.ndim != 1 or len(vector) < 7 or (len(vector)-3) % 4
            or not np.all(np.isfinite(vector)) or np.min(vector[-2:]) < 1):
        raise ValueError('finite paired coefficients, positive kappa, and two multipliers at least one required')
    if type(bisections) is not int or not 1 <= bisections <= 40:
        raise ValueError('frontier bisections must be from 1 through 40')
    state = score_state(vector[:-2], points)
    multipliers = vector[-2:]
    boundaries, jacobian = [], []
    for branch in range(2):
        lo, hi = .5175, .5195
        low = at_eta(state, multipliers, lo)[0][branch] - lo + .0001
        high = at_eta(state, multipliers, hi)[0][branch] - hi + .0001
        if not low < 0 < high:
            raise ValueError('joint frontier outside the bounded eta bracket')
        for _ in range(bisections):
            mid = (lo + hi) / 2
            residual = at_eta(state, multipliers, mid)[0][branch] - mid + .0001
            if residual < 0:
                lo = mid
            else:
                hi = mid
        eta = (lo + hi) / 2
        _, gradient, eta_slopes = at_eta(state, multipliers, eta)
        if eta_slopes[branch] <= 1:
            raise ValueError('nonmonotone frontier branch')
        boundaries.append(eta)
        jacobian.append(-gradient[branch] / (eta_slopes[branch] - 1))
    return np.array(boundaries), np.array(jacobian)


def initial_vector(degree=20, start=None):
    if type(degree) is not int or degree not in [14, 20, 32]:
        raise ValueError('degree must be 14, 20, or 32')
    baseline = json.loads((Path(start) if start else HERE/'start.json').read_text())
    if 'candidate' in baseline:
        baseline = baseline['candidate']
    design = baseline['design']
    old_degree = len(design['r'])
    if not 0 < old_degree <= degree or len(design['s']) != old_degree:
        raise ValueError('start degree must fit the target degree')
    arrays = [design['r'], design.get('r_sin', ['0']*old_degree),
              design['s'], design.get('s_sin', ['0']*old_degree)]
    if any(len(v) != old_degree for v in arrays):
        raise ValueError('all start coefficient arrays must have the same length')
    padded = [list(map(float, array))+[0.]*(degree-old_degree) for array in arrays]
    return np.array([v for array in padded for v in array]+[float(design['kappa']),
                     *map(float, baseline['multipliers'])])


def smooth_objective(boundaries, jacobian, temperature):
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError('finite positive smoothing temperature required')
    losses = -1000 * boundaries
    weights = np.exp((losses - max(losses)) / temperature)
    value = max(losses) + temperature * np.log(weights.sum())
    gradient = -1000 * (weights @ jacobian) / weights.sum()
    return value, gradient


def optimize(iterations=100, max_evaluations=4000, degree=20, start=None, author='conorholds'):
    if type(iterations) is not int or not 1 <= iterations <= 200:
        raise ValueError('iteration limit must be from 1 through 200')
    if type(max_evaluations) is not int or not 2 <= max_evaluations <= 4000:
        raise ValueError('evaluation limit must be from 2 through 4000')
    if not isinstance(author, str) or not author:
        raise ValueError('explicit nonempty author required')
    started = time.perf_counter()
    start_path = Path(start) if start else HERE/'start.json'
    x = initial_vector(degree, start_path)
    low, high = x - .35, x + .35
    low[-3:], high[-3:] = [.95, 7., 7.], [1.12, 9.99, 9.15]
    evaluations = 1
    boundaries, jacobian = frontier_and_jacobian(x)
    initial = boundaries.copy()
    best, best_boundaries = x.copy(), boundaries.copy()
    stages, rejected = [], 0
    for temperature in [.01, .001, .0001, .00001]:
        value, gradient = smooth_objective(boundaries, jacobian, temperature)
        inverse, accepted = np.eye(len(x)), 0
        for _ in range(iterations):
            if evaluations >= max_evaluations - 1:
                break
            direction = -inverse @ gradient
            if gradient @ direction >= 0:
                direction, inverse = -gradient, np.eye(len(x))
            direction /= max(1., np.max(np.abs(direction)) / .05)
            found = False
            step = 1.
            for _ in range(24):
                if evaluations >= max_evaluations - 1:
                    break
                trial = np.clip(x + step * direction, low, high)
                displacement = trial - x
                if np.max(np.abs(displacement)) < 1e-12 or gradient @ displacement >= 0:
                    break
                evaluations += 1
                try:
                    trial_boundaries, trial_jacobian = frontier_and_jacobian(trial)
                except ValueError as error:
                    if str(error) != 'joint frontier outside the bounded eta bracket':
                        raise
                    rejected += 1
                    step /= 2
                    continue
                trial_value, trial_gradient = smooth_objective(trial_boundaries, trial_jacobian, temperature)
                if min(trial_boundaries) > min(best_boundaries):
                    best, best_boundaries = trial.copy(), trial_boundaries.copy()
                if trial_value <= value + 1e-4 * (gradient @ displacement):
                    found = True
                    break
                step /= 2
            if not found:
                break
            change = trial_gradient - gradient
            curvature = displacement @ change
            if curvature > 1e-14:
                transform = np.eye(len(x)) - np.outer(displacement, change) / curvature
                inverse = transform @ inverse @ transform.T + np.outer(displacement, displacement) / curvature
            else:
                inverse = np.eye(len(x))
            x, boundaries, jacobian = trial, trial_boundaries, trial_jacobian
            value, gradient = trial_value, trial_gradient
            accepted += 1
        stages.append({'temperature': temperature, 'accepted': accepted, 'evaluations': evaluations,
                       'eta_boundary': float(min(boundaries))})
    degree = (len(best)-3)//4
    design = {key: [f'{v:.12f}' for v in best[i*degree:(i+1)*degree]]
              for i,key in enumerate(['r', 'r_sin', 's', 's_sin'])}
    design['kappa'] = f'{best[-3]:.12f}'
    multipliers = [f'{v:.12f}' for v in best[-2:]]
    serialized = np.concatenate((vector_from_design(design), np.array(list(map(float,multipliers)))))
    if not np.array_equal(serialized, best):
        best_boundaries, _ = frontier_and_jacobian(serialized)
        evaluations += 1
    candidate = {'id': f'asymmetric_degree{degree}_retuned', 'design': design, 'multipliers': multipliers,
                 'parameters': FRACTIONAL.exact_parameters(FRACTIONAL.slice_parameters(min(best_boundaries)))}
    report = {'status': 'EXPLORATORY', 'author': author, 'initial_boundaries': initial.tolist(),
              'candidate_boundaries': best_boundaries.tolist(), 'design': design,
              'multipliers': multipliers, 'candidate': candidate,
              'evaluations': evaluations, 'rejected_outside_bracket': rejected, 'stages': stages,
              'limits': {'iterations_per_stage': iterations, 'evaluations': max_evaluations,
                         'coefficient_displacement': '.35', 'line_search_trials': 24,
                         'eta': ['.5175', '.5195'], 'frontier_bisections': 40, 'grid': 128,
                         'kappa': ['.95', '1.12'], 'multipliers': [['7', '9.99'], ['7', '9.15']]},
              'inputs_sha256': {'search.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                'spectral.py': hashlib.sha256((HERE/'spectral.py').read_bytes()).hexdigest(),
                                'fractional_search.py': hashlib.sha256((HERE.parent/'fractional-coverage/search.py').read_bytes()).hexdigest(),
                                'predecessor.json': hashlib.sha256((HERE/'predecessor.json').read_bytes()).hexdigest()},
              'shared_sha256': {name: hashlib.sha256((tune.CORE / name).read_bytes()).hexdigest()
                                for name in ['filters.py', 'SOURCE.json', 'requirements.txt']},
              'start': {'path': start_path.name, 'sha256': hashlib.sha256(start_path.read_bytes()).hexdigest()},
              'degree': degree, 'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
              'seconds': time.perf_counter() - started}
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--degree', type=int, choices=[14, 20, 32], default=20)
    parser.add_argument('--author', required=True)
    parser.add_argument('--start', type=Path)
    args = parser.parse_args()
    result = optimize(degree=args.degree, start=args.start, author=args.author)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['status', 'initial_boundaries', 'candidate_boundaries',
                                                 'evaluations', 'seconds']}))
