"""Exploratory discrete-density relaxation; no continuous or algorithm claim."""
import argparse
import hashlib
import platform
import json
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent/'fractional-coverage'
START = json.loads((PACKAGE/'candidates.json').read_text())['candidates'][1]
N = 128
X = np.arange(N)*2*np.pi/N
SHIFTS = (np.arange(N)[:, None]+np.arange(N)[None, :]) % N
M2 = 2+2*np.cos(X)
ETA, TAU = [float(START['parameters'][key]) for key in ['eta', 'tau']]
DELTA, RATIO = ETA-TAU, (ETA-TAU)/(2*TAU-ETA)
LAM = np.array(list(map(float, START['multipliers'])))
LN2 = np.log(2.)


def density(values):
    values = np.asarray(values)
    if not np.all(np.isfinite(values)) or values.min() <= 0:
        raise ValueError('strictly positive finite density required')
    return values/values.mean()


def initializations():
    design = START['design']
    cs = np.cos(np.outer(np.arange(1, len(design['r'])+1), X))
    k = float(design['kappa'])
    baseline = tuple(density(np.exp2(k*np.array(list(map(float, design[key])))@cs)) for key in ['r', 's'])
    wrap = lambda x: (x+np.pi)%(2*np.pi)-np.pi
    return {
        'retained_power_filters': baseline,
        'uniform_fixed_point_control': (np.ones(N), np.ones(N)),
        'opposed_single_peaks': (density(np.exp2(2*np.cos(X))), density(np.exp2(2*np.cos(X-np.pi)))),
        'split_and_single_peaks': (density(np.exp2(1.5*np.cos(2*X)+.2*np.cos(X))), density(np.exp2(2*np.cos(X+.3)))),
        'poisson_shifted': (density((1-.75**2)/(1-2*.75*np.cos(X)+.75**2)),
                            density((1-.9**2)/(1-2*.9*np.cos(X-.7)+.9**2))),
        'triangular_supports': (density(1e-4+np.maximum(0, 1-np.abs(wrap(X))/(np.pi*.7))),
                                density(1e-4+np.maximum(0, 1-np.abs(wrap(X-1.2))/(np.pi*.45)))),
    }


def state(h, g, theta):
    values = np.log2(h)[SHIFTS]+np.log2(g)[None, :]
    raw = np.exp2(theta[:, None]*values)
    p = raw.mean(axis=1)
    conditional = raw/raw.sum(axis=1)[:, None]
    b = (conditional*values).sum(axis=1)
    d = theta*b-np.log2(p)
    s = b
    exponents = np.array([LAM[0]*(DELTA+d-s), LAM[1]*((1+2*RATIO)*d-(1+RATIO)*s)])
    return {'K': conditional, 'D': d, 'S': s,
            'integrands': M2*np.exp2(exponents), 'theta': theta}


def optimize_kernel(h, g, beta):
    lo, hi = np.full(N, (1+RATIO)/(1+2*RATIO)), np.ones(N)
    for _ in range(36):
        theta = (lo+hi)/2
        st = state(h, g, theta)
        # M2 cancels; use the two exponents without M2 at the zero node.
        d, s = st['D'], st['S']
        first = np.exp2(LAM[0]*(DELTA+d-s))
        second = np.exp2(LAM[1]*((1+2*RATIO)*d-(1+RATIO)*s))
        derivative = (beta*LAM[0]*first*(theta-1)
                      +(1-beta)*LAM[1]*second*((1+2*RATIO)*theta-(1+RATIO)))
        hi = np.where(derivative > 0, theta, hi)
        lo = np.where(derivative > 0, lo, theta)
    return state(h, g, (lo+hi)/2)


def objective(st, beta):
    return beta*st['integrands'][0].mean()+(1-beta)*st['integrands'][1].mean()


def fixed_kernel_objective(h, g, st, beta):
    scores = (st['K']*(np.log2(h)[SHIFTS]+np.log2(g)[None, :])).sum(axis=1)
    values = [LAM[0]*(DELTA+st['D']-scores),
              LAM[1]*((1+2*RATIO)*st['D']-(1+RATIO)*scores)]
    return beta*np.mean(M2*np.exp2(values[0]))+(1-beta)*np.mean(M2*np.exp2(values[1]))


def run_one(name, h, g, beta, iterations=1200):
    started = time.perf_counter()
    initial = optimize_kernel(h, g, beta)
    st, old = initial, objective(initial, beta)
    accepted, evaluated, largest_increase = 0, 1, 0.
    for iteration in range(iterations):
        weights = (beta*LAM[0]*st['integrands'][0]
                   +(1-beta)*LAM[1]*(1+RATIO)*st['integrands'][1])/N
        joint = weights[:, None]*st['K']
        mass = weights.sum()
        target_h = N*np.bincount(SHIFTS.ravel(), weights=joint.ravel(), minlength=N)/mass
        target_g = N*joint.sum(axis=0)/mass
        found, step = False, 1.
        for _ in range(24):
            nh, ng = density((1-step)*h+step*target_h), density((1-step)*g+step*target_g)
            trial = fixed_kernel_objective(nh, ng, st, beta)
            evaluated += 1
            if trial <= old-1e-14:
                found = True
                break
            step /= 2
        if not found:
            break
        h, g = nh, ng
        st = optimize_kernel(h, g, beta)
        current = objective(st, beta)
        largest_increase = max(largest_increase, current-old)
        if current > old+1e-12:
            raise ValueError('kernel update increased the scalar objective')
        accepted += 1
        old = current
        if iteration % 200 == 199:
            print(json.dumps({'name': name, 'iteration': iteration+1, 'log_moments': np.log2(st['integrands'].mean(axis=1)).tolist()}), flush=True)
    return {'name': name, 'status': 'EXPLORATORY', 'beta': beta, 'accepted': accepted,
            'iterations_limit': iterations, 'value_evaluations': evaluated,
            'initial_log_moments': np.log2(initial['integrands'].mean(axis=1)).tolist(),
            'log_moments': np.log2(st['integrands'].mean(axis=1)).tolist(),
            'objective': old, 'largest_kernel_update_increase': largest_increase,
            'test_h': [f'{v:.14g}' for v in h], 'test_g': [f'{v:.14g}' for v in g],
            'kernel_temperature': [f'{v:.14g}' for v in st['theta']],
            'seconds': time.perf_counter()-started}


def project_to_fourier(row, degree=20, theta='.96735'):
    if type(degree) is not int or not 1 <= degree < N//2:
        raise ValueError('projection degree must be between one and half the grid minus one')
    theta_value = float(theta)
    if not np.isfinite(theta_value) or theta_value <= 0:
        raise ValueError('positive finite projection temperature required')
    design = {}
    for key, field in [('r', 'test_h'), ('s', 'test_g')]:
        values = np.array(list(map(float, row[field])))
        coefficients = np.fft.rfft(np.log2(density(values)))/N
        design[key] = [f'{2*theta_value*v.real:.12f}' for v in coefficients[1:degree+1]]
        design[key+'_sin'] = [f'{-2*theta_value*v.imag:.12f}' for v in coefficients[1:degree+1]]
    design['kappa'] = f'{1/theta_value:.12f}'
    return {'id': 'projected_density_start', 'design': design,
            'parameters': START['parameters'], 'multipliers': START['multipliers']}


def run(author='conorholds'):
    if not isinstance(author, str) or not author:
        raise ValueError('explicit nonempty author required')
    beta = .4245
    reports = []
    for name, (h, g) in initializations().items():
        report = run_one(name, h, g, beta)
        reports.append(report)
        print(json.dumps({key: report[key] for key in ['name', 'accepted', 'log_moments', 'seconds']}), flush=True)
    selected = min(reports, key=lambda row: max(row['log_moments']))
    summary = [{key: value for key,value in row.items()
                if key not in ['test_h', 'test_g', 'kernel_temperature']} for row in reports]
    return {'status': 'EXPLORATORY', 'author': author,
            'scope': 'Discrete-grid scalarized relaxation; no continuous certificate or full-budget claim.',
            'grid': N, 'limits': {'initializations': 6, 'iterations_per_initialization': 1200,
                                  'kernel_bisections': 36, 'backtracking_trials': 24},
            'candidate_sha256': hashlib.sha256((PACKAGE/'candidates.json').read_bytes()).hexdigest(),
            'search_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'reports': summary, 'selected': selected,
            'projected_start': project_to_fourier(selected)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.author)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'selected': result['selected']['name'],
                      'log_moments': result['selected']['log_moments']}))
