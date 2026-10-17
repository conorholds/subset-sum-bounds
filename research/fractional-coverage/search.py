"""Bounded fixed-filter retuning after removing the integer packing loss."""
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
sys.path.insert(0, str(HERE.parent / 'chernoff-multipliers'))
import multiplier_search as previous

H = previous.H
COVERAGE = .0000016


def slice_parameters(eta):
    if not np.isfinite(eta):
        raise ValueError('finite eta required')
    rho = (eta-.500004)/(H-.5)
    if not 0 < rho < .5:
        raise ValueError('outside the fractional-coverage density slice')
    tau = (.5+(H-.5+.00005)*rho+COVERAGE+.000016)/(1+rho)
    return {'eta': eta, 'tau': tau, 'rho': rho}


def exact_parameters(row):
    eta = Decimal(str(row['eta'])).quantize(Decimal('.00000001'), rounding=ROUND_FLOOR)
    fresh = slice_parameters(float(eta))
    rho = Decimal(str(fresh['rho'])).quantize(Decimal('.00000001'), rounding=ROUND_HALF_EVEN)
    tau = Decimal(str(fresh['tau'])).quantize(Decimal('.00000001'), rounding=ROUND_CEILING)+Decimal('.0000001')
    d = Decimal('.5')-tau+Decimal('.00002')
    delta = (d+Decimal('.0000016')/rho).quantize(Decimal('.0000000001'), rounding=ROUND_CEILING)
    return {'eta': str(eta), 'tau': str(tau), 'rho': str(rho),
            'deficiency': str(d), 'local_deficiency': str(delta),
            'gamma': str(Decimal('.5')-delta-Decimal('.00003')),
            'guard': '.00001', 'preliminary': str(tau-Decimal('.00001')),
            'alpha': str(tau+Decimal('.00001')), 'epsilon': '.00000001',
            'moment_ceiling': str(eta-Decimal('.0001'))}


def search(candidate, tune_multipliers=True, points=256, bisections=45, author='conorholds'):
    if type(points) is not int or points not in [128, 256, 512]:
        raise ValueError('grid must be 128, 256, or 512')
    if type(bisections) is not int or not 1 <= bisections <= 45:
        raise ValueError('bisections must be from 1 through 45')
    if type(tune_multipliers) is not bool:
        raise ValueError('multiplier tuning flag must be boolean')
    if not isinstance(author, str) or not author:
        raise ValueError('explicit nonempty author required')
    started = time.perf_counter()
    sc = previous.scores(candidate['design'], points)
    evaluations = 0

    def evaluate(eta):
        nonlocal evaluations
        row = slice_parameters(eta)
        if tune_multipliers:
            result = previous.tune_at(sc, row)
            evaluations += sum(r['evaluations'] for r in result)
        else:
            delta = row['eta']-row['tau']
            a = delta/(2*row['tau']-row['eta'])
            result = [{'multiplier': float(lam),
                       'log2_moment': previous.log_moment(direction, sc[2], float(lam))}
                      for direction, lam in zip([delta-sc[0], a*sc[1]-sc[0]], candidate['multipliers'])]
            evaluations += 2
        return max(r['log2_moment'] for r in result)-eta+.0001, result

    lo, hi = .5175, .51828
    if not evaluate(lo)[0] < 0 < evaluate(hi)[0]:
        raise ValueError('fractional frontier outside the bounded eta bracket')
    for _ in range(bisections):
        mid = (lo+hi)/2
        if evaluate(mid)[0] < 0:
            lo = mid
        else:
            hi = mid
    parameters = exact_parameters(slice_parameters(lo))
    if tune_multipliers:
        tuned = previous.tune_at(sc, parameters)
        evaluations += sum(r['evaluations'] for r in tuned)
        multipliers = [f'{r["multiplier"]:.12f}' for r in tuned]
    else:
        multipliers = candidate['multipliers']
    retained = {'id': 'fractional_retuned' if tune_multipliers else 'fractional_fixed_multipliers',
                'design': candidate['design'], 'multipliers': multipliers, 'parameters': parameters}
    return {'status': 'EXPLORATORY', 'author': author, 'candidate': retained,
            'frontier_bracket': [lo, hi], 'moment_evaluations': evaluations,
            'limits': {'grid': points, 'eta': ['.5175', '.51828'], 'bisections': bisections,
                       'multiplier_bisections': 60, 'coverage_reserve_per_n': '.0000016',
                       'tune_multipliers': tune_multipliers},
            'inputs_sha256': {'search.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                              'budget.py': hashlib.sha256((HERE/'budget.py').read_bytes()).hexdigest(),
                              'predecessor.json': hashlib.sha256((HERE/'predecessor.json').read_bytes()).hexdigest(),
                              'multiplier_search.py': hashlib.sha256((HERE.parent/'chernoff-multipliers/multiplier_search.py').read_bytes()).hexdigest(),
                              'filters.py': hashlib.sha256((HERE.parent/'filter-budgets/filters.py').read_bytes()).hexdigest(),
                              'SOURCE.json': hashlib.sha256((HERE.parent/'filter-budgets/SOURCE.json').read_bytes()).hexdigest()},
            'start_sha256': hashlib.sha256(json.dumps(candidate, sort_keys=True).encode()).hexdigest(),
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'seconds': time.perf_counter()-started}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predecessor', type=Path, default=HERE/'predecessor.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fixed-multipliers', action='store_true')
    parser.add_argument('--author', required=True)
    args = parser.parse_args()
    report = search(json.loads(args.predecessor.read_text())['candidate'], not args.fixed_multipliers,
                    author=args.author)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ['status', 'frontier_bracket', 'moment_evaluations', 'seconds']}))
