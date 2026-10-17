"""Bounded convex multiplier tuning on the retained 17-block parameter slice.

The filter coefficients are unchanged from PR 13. Floating-point search is
exploratory; interval_check.py supplies the separate numerical acceptance path.
"""
import argparse
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / 'filter-budgets'
sys.path.insert(0, str(CORE))
from filters import scores

PREDECESSOR = json.loads((HERE / 'predecessor.json').read_text())
H = -.25 * np.log2(.25) - .75 * np.log2(.75)


def slice_parameters(eta):
    """Same strictly positive reserves and balanced costs as PRs 12 and 13."""
    if not np.isfinite(eta):
        raise ValueError('finite eta required')
    rho = (eta - .500004) / (H - .5)
    if not 0 < 17 * rho < .999:
        raise ValueError('outside the unclipped 17-block slice')
    tau = (.5 + (H - 1 + .00003) * rho + .500046 / 17 + .000016) / (1 + 1 / 17)
    return {'eta': eta, 'tau': tau, 'rho': rho}


def log_moment(direction, m2, multiplier):
    return log_moment_and_derivative(direction, m2, multiplier)[0]


def log_moment_and_derivative(direction, m2, multiplier):
    exponent = multiplier * direction
    shift = np.max(exponent)
    weights = m2 * np.exp2(exponent - shift)
    return (float(shift + np.log2(np.mean(weights))),
            float(np.sum(weights * direction) / np.sum(weights)))


def minimize_multiplier(direction, m2, upper, iterations=60):
    """Bisect the monotone derivative of a convex log-sum-exp."""
    if type(iterations) is not int or not 1 <= iterations <= 60:
        raise ValueError('multiplier iterations must be from 1 through 60')
    if not np.isfinite(upper) or upper <= 1:
        raise ValueError('upper multiplier must be finite and exceed one')
    lo, hi = 1., upper
    low_value, low_gradient = log_moment_and_derivative(direction, m2, lo)
    high_value, high_gradient = log_moment_and_derivative(direction, m2, hi)
    if low_gradient >= 0:
        return {'multiplier': lo, 'log2_moment': low_value, 'evaluations': 2}
    if high_gradient <= 0:
        return {'multiplier': hi, 'log2_moment': high_value, 'evaluations': 2}
    evaluations = 2
    for _ in range(iterations):
        mid = (lo + hi) / 2
        _, derivative = log_moment_and_derivative(direction, m2, mid)
        evaluations += 1
        if derivative > 0:
            hi = mid
        else:
            lo = mid
    multiplier = (lo + hi) / 2
    value = log_moment(direction, m2, multiplier)
    evaluations += 1
    return {'multiplier': multiplier, 'log2_moment': value, 'evaluations': evaluations}


def tune_at(sc, parameters, iterations=60):
    u, v, m2 = sc
    eta, tau = float(parameters['eta']), float(parameters['tau'])
    if not 0 < tau < eta < 2 * tau:
        raise ValueError('require 0 < tau < eta < 2 tau')
    delta, ratio = eta - tau, (eta - tau) / (2 * tau - eta)
    # Keep at least .01 epsilon of the original 30 epsilon domination reserve.
    upper = [29.99 / 3, 29.99 / (4 * ratio + 3)]
    return [minimize_multiplier(direction, m2, limit, iterations)
            for direction, limit in zip([delta - u, ratio * v - u], upper)]


def exact_parameters(row):
    eta = Decimal(str(row['eta'])).quantize(Decimal('.00000001'), rounding=ROUND_FLOOR)
    recalculated = slice_parameters(float(eta))
    rho = Decimal(str(recalculated['rho'])).quantize(Decimal('.00000001'), rounding=ROUND_HALF_EVEN)
    tau = Decimal(str(recalculated['tau'])).quantize(Decimal('.00000001'), rounding=ROUND_CEILING) + Decimal('.0000001')
    deficiency = Decimal('.5') - tau + Decimal('.00002')
    local = ((deficiency + Decimal('.000026')) / (17 * rho)).quantize(Decimal('.0000000001'), rounding=ROUND_CEILING)
    return {'eta': str(eta), 'tau': str(tau), 'rho': str(rho), 'blocks': 17,
            'deficiency': str(deficiency), 'local_deficiency': str(local),
            'gamma': str(Decimal('.5') - local - Decimal('.00003')),
            'guard': '.00001', 'preliminary': str(tau - Decimal('.00001')),
            'alpha': str(tau + Decimal('.00001')), 'epsilon': '.00000001',
            'moment_ceiling': str(eta - Decimal('.0001'))}


def search(points=256, bisections=45, multiplier_iterations=60):
    if type(points) is not int or points not in [128, 256, 512]:
        raise ValueError('grid must be 128, 256, or 512')
    if type(bisections) is not int or not 1 <= bisections <= 45:
        raise ValueError('frontier bisections must be from 1 through 45')
    started = time.perf_counter()
    predecessor = PREDECESSOR['witness']
    sc = scores(predecessor['design'], points)
    evaluations = 0

    def evaluate(eta):
        nonlocal evaluations
        row = slice_parameters(eta)
        result = tune_at(sc, row, multiplier_iterations)
        evaluations += sum(x['evaluations'] for x in result)
        return max(x['log2_moment'] for x in result) - eta + .0001

    lo, hi = .5179, .5181
    if not evaluate(lo) < 0 < evaluate(hi):
        raise ValueError('retained frontier bracket has no sign change')
    for _ in range(bisections):
        mid = (lo + hi) / 2
        if evaluate(mid) < 0:
            lo = mid
        else:
            hi = mid
    parameters = exact_parameters(slice_parameters(lo))
    tuned = tune_at(sc, parameters, multiplier_iterations)
    evaluations += sum(x['evaluations'] for x in tuned)
    multipliers = [f'{x["multiplier"]:.9f}' for x in tuned]
    baseline = dict(predecessor, id='predecessor_degree14', multipliers=['8.5', '8.3'])
    fixed_tuning = tune_at(sc, predecessor['parameters'], multiplier_iterations)
    evaluations += sum(x['evaluations'] for x in fixed_tuning)
    fixed = dict(predecessor, id='multipliers_only', multipliers=[f'{x["multiplier"]:.9f}' for x in fixed_tuning])
    candidate = {'id': 'multipliers_retuned', 'design': predecessor['design'],
                 'parameters': parameters, 'multipliers': multipliers}
    candidates = {'model': 'chernoff-multipliers-v1', 'candidates': [baseline, fixed, candidate]}
    report = {'status': 'EXPLORATORY', 'author': 'conorholds',
              'predecessor_sha256': hashlib.sha256((HERE / 'predecessor.json').read_bytes()).hexdigest(),
              'search_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'inputs_sha256': {name: hashlib.sha256((CORE / name).read_bytes()).hexdigest()
                                for name in ['filters.py', 'certify.py', 'SOURCE.json']},
              'limits': {'points': points, 'frontier_bisections': bisections,
                         'multiplier_iterations': multiplier_iterations, 'eta': ['.5179', '.5181'],
                         'domination_rounding_budget': '29.99'},
              'frontier_bracket': [lo, hi], 'moment_evaluations': evaluations,
              'parameters': parameters, 'multipliers': multipliers,
              'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
              'seconds': time.perf_counter() - started}
    return candidates, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--candidates-output', type=Path, required=True)
    args = parser.parse_args()
    candidates, report = search()
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    args.candidates_output.write_text(json.dumps(candidates, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ['status', 'moment_evaluations', 'parameters', 'multipliers', 'seconds']}))


if __name__ == '__main__':
    main()
