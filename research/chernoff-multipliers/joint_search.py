"""Exploratory joint filter/multiplier search for the actual balanced frontier.

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
import multiplier_search as tune

def score_state(vector, points=128):
    """Differentiate normalizers, correlations and both log2 moments analytically."""
    vector = np.asarray(vector, dtype=float)
    if (vector.ndim != 1 or len(vector) < 3 or len(vector) % 2 != 1
            or not np.all(np.isfinite(vector)) or vector[-1] <= 0):
        raise ValueError('finite paired coefficient vectors and positive kappa required')
    if type(points) is not int or points < 64 or points % 2:
        raise ValueError('even quadrature grid of at least 64 required')
    degree = (len(vector)-1)//2
    r, s, k = vector[:degree], vector[degree:2*degree], vector[-1]
    x = np.arange(points)*(2*np.pi/points)
    cosines = np.cos(np.outer(np.arange(1, degree+1), x))
    rv, sv = r@cosines, s@cosines
    ln2 = np.log(2.)

    def correlate(left, right):
        return np.fft.ifft(np.fft.fft(left)*np.conj(np.fft.fft(right))).real/points

    with np.errstate(over='raise', invalid='raise', divide='raise'):
        f, g = np.exp2(rv), np.exp2(sv)
        p = correlate(f, g)
        y = correlate(f*rv, g)+correlate(f, g*sv)
        if np.min(p) <= 0:
            raise ValueError('quadrature lost positive convolution')
        z = np.log2(f.mean())+np.log2(g.mean())
        fk, gk = np.exp2(k*rv), np.exp2(k*sv)
        h = np.log2(fk.mean())+np.log2(gk.mean())-k*z
        ell, b = np.log2(p)-z, y/p-z
        df, dg = ln2*f*cosines, ln2*g*cosines
        dp = np.concatenate((correlate(df, g), correlate(f, dg), np.zeros((1, points))))
        dy = np.concatenate((correlate(f*cosines*(1+ln2*rv), g)+correlate(df, g*sv),
                             correlate(f*rv, dg)+correlate(f, g*cosines*(1+ln2*sv)),
                             np.zeros((1, points))))
        dz = np.concatenate(((f*cosines).mean(axis=1)/f.mean(),
                             (g*cosines).mean(axis=1)/g.mean(), [0.]))
        dh = k*np.concatenate(((fk*cosines).mean(axis=1)/fk.mean(),
                                (gk*cosines).mean(axis=1)/gk.mean(), [0.]))-k*dz
        dh[-1] = (fk*rv).mean()/fk.mean()+(gk*sv).mean()/gk.mean()-z
        dell = dp/(ln2*p)-dz[:, None]
        db = dy/p-y*dp/p**2-dz[:, None]
        du = dell+(k-1)*db-dh[:, None]
        dv = (2-k)*db-2*dell+dh[:, None]
        du[-1] += b
        dv[-1] -= b
        u, v = ell+(k-1)*b-h, (2-k)*b-2*ell+h
        return {'u': u, 'v': v, 'm2': 2+2*np.cos(x), 'du': du, 'dv': dv}


def at_eta(state, multipliers, eta):
    row = tune.slice_parameters(eta)
    tau = row['tau']
    tau_prime = (tune.H - 1 + .00003) / ((tune.H - .5) * (1 + 1 / 17))
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
    if (vector.ndim != 1 or len(vector) < 5 or len(vector) % 2 != 1
            or not np.all(np.isfinite(vector)) or np.min(vector[-2:]) < 1):
        raise ValueError('finite paired coefficients, positive kappa, and two multipliers at least one required')
    if type(bisections) is not int or not 1 <= bisections <= 40:
        raise ValueError('frontier bisections must be from 1 through 40')
    state = score_state(vector[:-2], points)
    multipliers = vector[-2:]
    boundaries, jacobian = [], []
    for branch in range(2):
        lo, hi = .5175, .51828
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


def initial_vector(degree=14, start=None):
    if type(degree) is not int or degree not in [14, 20]:
        raise ValueError('degree must be 14 or 20')
    baseline = json.loads(Path(start).read_text())['candidate'] if start else tune.PREDECESSOR['witness']
    design = baseline['design']
    if len(design['r']) != len(design['s']) or not 0 < len(design['r']) <= degree:
        raise ValueError('start degree must fit the target degree')
    sc = tune.scores(design)
    multipliers = baseline.get('multipliers')
    if multipliers is None:
        multipliers = [f'{row["multiplier"]:.9f}' for row in tune.tune_at(sc, baseline['parameters'])]
    padding = [0.] * (degree - len(design['r']))
    return np.array([*map(float, design['r']), *padding, *map(float, design['s']), *padding,
                     float(design['kappa']), *map(float, multipliers)])


def smooth_objective(boundaries, jacobian, temperature):
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError('finite positive smoothing temperature required')
    losses = -1000 * boundaries
    weights = np.exp((losses - max(losses)) / temperature)
    value = max(losses) + temperature * np.log(weights.sum())
    gradient = -1000 * (weights @ jacobian) / weights.sum()
    return value, gradient


def optimize(iterations=100, max_evaluations=4000, degree=14, start=None):
    if type(iterations) is not int or not 1 <= iterations <= 200:
        raise ValueError('iteration limit must be from 1 through 200')
    if type(max_evaluations) is not int or not 1 <= max_evaluations <= 4000:
        raise ValueError('evaluation limit must be from 1 through 4000')
    started = time.perf_counter()
    x = initial_vector(degree, start)
    low, high = x - .35, x + .35
    low[-3:], high[-3:] = [.95, 7.5, 7.5], [1.12, 9.5, 9.]
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
    degree = (len(best) - 3) // 2
    design = {'r': [f'{v:.9f}' for v in best[:degree]],
              's': [f'{v:.9f}' for v in best[degree:2*degree]], 'kappa': f'{best[-3]:.9f}'}
    multipliers = [f'{v:.9f}' for v in best[-2:]]
    serialized = np.array([*map(float, design['r']), *map(float, design['s']),
                           float(design['kappa']), *map(float, multipliers)])
    if not np.array_equal(serialized, best):
        best_boundaries, _ = frontier_and_jacobian(serialized)
        evaluations += 1
    candidate = {'id': f'joint_degree{degree}_retuned', 'design': design, 'multipliers': multipliers,
                 'parameters': tune.exact_parameters(tune.slice_parameters(min(best_boundaries)))}
    report = {'status': 'EXPLORATORY', 'initial_boundaries': initial.tolist(),
              'candidate_boundaries': best_boundaries.tolist(), 'design': design,
              'multipliers': multipliers, 'candidate': candidate,
              'evaluations': evaluations, 'rejected_outside_bracket': rejected, 'stages': stages,
              'limits': {'iterations_per_stage': iterations, 'evaluations': max_evaluations,
                         'coefficient_displacement': '.35', 'line_search_trials': 24,
                         'eta': ['.5175', '.51828'], 'frontier_bisections': 40, 'grid': 128},
              'inputs_sha256': {'joint_search.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                'multiplier_search.py': hashlib.sha256((HERE / 'multiplier_search.py').read_bytes()).hexdigest(),
                                'predecessor.json': hashlib.sha256((HERE / 'predecessor.json').read_bytes()).hexdigest()},
              'shared_sha256': {name: hashlib.sha256((tune.CORE / name).read_bytes()).hexdigest()
                                for name in ['filters.py', 'SOURCE.json', 'requirements.txt']},
              'start': {'path': Path(start).name, 'sha256': hashlib.sha256(Path(start).read_bytes()).hexdigest()} if start else None,
              'degree': degree, 'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
              'seconds': time.perf_counter() - started}
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--degree', type=int, choices=[14, 20], default=14)
    parser.add_argument('--start', type=Path)
    args = parser.parse_args()
    result = optimize(degree=args.degree, start=args.start)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ['status', 'initial_boundaries', 'candidate_boundaries',
                                                 'evaluations', 'seconds']}))
