"""Bounded analytic-gradient minimax continuation for the existing cosine family."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent/'filter-budgets'
sys.path.insert(0, str(CORE))
from verify import load_witnesses

BASELINE = json.loads((HERE/'baseline.json').read_text())
PARAMETERS = BASELINE['witness']['parameters']
TEMPERATURES = (.001, .0001, .00001, .000001)


def unpack(vector):
    degree = (len(vector)-1)//2
    return {'r': [f'{v:.9f}' for v in vector[:degree]],
            's': [f'{v:.9f}' for v in vector[degree:2*degree]],
            'kappa': f'{vector[-1]:.9f}'}


def initial_vector(degree):
    if type(degree) is not int or not 10 <= degree <= 14:
        raise ValueError('degree must be an integer from 10 through 14')
    design = next(w for w in load_witnesses() if w['id'] == 'combined')['design']
    return np.array([*map(float, design['r']), *([0.]*(degree-10)),
                     *map(float, design['s']), *([0.]*(degree-10)), float(design['kappa'])])


def moments_and_jacobian(vector, parameters=PARAMETERS, points=128):
    """Differentiate normalizers, correlations and both log2 moments analytically."""
    vector = np.asarray(vector, dtype=float)
    if (vector.ndim != 1 or len(vector) < 3 or len(vector) % 2 != 1
            or not np.all(np.isfinite(vector)) or vector[-1] <= 0):
        raise ValueError('finite paired coefficient vectors and positive kappa required')
    if type(points) is not int or points < 64 or points % 2:
        raise ValueError('even quadrature grid of at least 64 required')
    eta, tau = float(parameters['eta']), float(parameters['tau'])
    if not 0 < tau < eta < 2*tau:
        raise ValueError('require 0 < tau < eta < 2 tau')
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
        a = (eta-tau)/(2*tau-eta)
        exponents = np.array([8.5*(eta-tau-u), 8.3*(a*v-u)])
        shifted = np.exp2(exponents-exponents.max(axis=1)[:, None])*(2+2*np.cos(x))
        logs = exponents.max(axis=1)+np.log2(shifted.mean(axis=1))
        weights = shifted/shifted.sum(axis=1)[:, None]
        jacobian = np.array([-8.5*(du@weights[0]), 8.3*((a*dv-du)@weights[1])])
    if not np.all(np.isfinite([*logs, *jacobian.ravel()])):
        raise ValueError('nonfinite moments or derivatives')
    return logs, jacobian


def smooth_max(logs, jacobian, temperature):
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError('finite positive smoothing temperature required')
    weights = np.exp((logs-max(logs))/temperature)
    return max(logs)+temperature*np.log(weights.sum()), (weights@jacobian)/weights.sum()


def optimize(degree=10, max_evaluations=12000, iterations=200):
    if type(max_evaluations) is not int or not 1 <= max_evaluations <= 12000:
        raise ValueError('evaluation budget must be from 1 through 12000')
    if type(iterations) is not int or not 1 <= iterations <= 200:
        raise ValueError('stage iteration budget must be from 1 through 200')
    x = initial_vector(degree)
    lo, hi = x-.35, x+.35
    lo[-1], hi[-1] = .95, 1.12
    started = time.perf_counter()
    evaluations = 0
    initial_logs, jacobian = moments_and_jacobian(x)
    evaluations += 1
    logs = initial_logs.copy()
    best, best_logs = x.copy(), initial_logs.copy()
    stages = []
    for temperature in TEMPERATURES:
        if evaluations >= max_evaluations-1:
            break
        value, gradient = smooth_max(logs, jacobian, temperature)
        inverse = np.eye(len(x))
        accepted = 0
        for iteration in range(iterations):
            if evaluations >= max_evaluations-1:
                break
            direction = -inverse@gradient
            if gradient@direction >= 0:
                direction = -gradient
                inverse = np.eye(len(x))
            direction /= max(1., np.max(np.abs(direction))/.1)
            step = 1.
            found = False
            for _ in range(24):
                if evaluations >= max_evaluations-1:
                    break
                trial = np.clip(x+step*direction, lo, hi)
                displacement = trial-x
                if np.max(np.abs(displacement)) < 1e-12 or gradient@displacement >= 0:
                    break
                trial_logs, trial_jacobian = moments_and_jacobian(trial)
                evaluations += 1
                trial_value, trial_gradient = smooth_max(trial_logs, trial_jacobian, temperature)
                if max(trial_logs) < max(best_logs):
                    best, best_logs = trial.copy(), trial_logs.copy()
                if trial_value <= value+1e-4*(gradient@displacement):
                    found = True
                    break
                step /= 2
            if not found:
                break
            change = trial_gradient-gradient
            curvature = displacement@change
            if curvature > 1e-14:
                transform = np.eye(len(x))-np.outer(displacement, change)/curvature
                inverse = transform@inverse@transform.T+np.outer(displacement, displacement)/curvature
            else:
                inverse = np.eye(len(x))
            x, logs, jacobian, value, gradient = trial, trial_logs, trial_jacobian, trial_value, trial_gradient
            accepted += 1
        stages.append({'temperature': temperature, 'accepted_steps': accepted,
                       'max_log2_moment': float(max(logs)), 'evaluations_so_far': evaluations})
        if evaluations >= max_evaluations-1:
            break
    design = unpack(best)
    serialized = np.array([*map(float, design['r']), *map(float, design['s']), float(design['kappa'])])
    if not np.array_equal(serialized, best):
        best_logs, _ = moments_and_jacobian(serialized)
        evaluations += 1
    return {'status': 'EXPLORATORY', 'scope': 'Bounded floating-point search; separate interval replay required.',
            'degree': degree, 'parameters': PARAMETERS, 'temperatures': list(TEMPERATURES),
            'limits': {'evaluations': max_evaluations, 'iterations_per_stage': iterations,
                       'line_search_trials': 24, 'coefficient_displacement': '.35', 'kappa': ['.95', '1.12']},
            'evaluations': evaluations, 'stages': stages, 'starting_log2_moments': initial_logs.tolist(),
            'candidate_log2_moments': best_logs.tolist(), 'design': design,
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'optimizer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'inputs_sha256': {name: hashlib.sha256((CORE/name).read_bytes()).hexdigest()
                              for name in ['filters.py', 'verify.py', 'certify.py', 'SOURCE.json', 'witnesses.json']},
            'baseline_sha256': hashlib.sha256((HERE/'baseline.json').read_bytes()).hexdigest(),
            'seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=10)
    parser.add_argument('--evaluations', type=int, default=12000)
    parser.add_argument('--iterations', type=int, default=200)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = optimize(args.degree, args.evaluations, args.iterations)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['degree', 'evaluations', 'starting_log2_moments',
                                                  'candidate_log2_moments', 'seconds']}))


if __name__ == '__main__':
    main()
