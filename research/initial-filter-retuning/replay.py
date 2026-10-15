"""Replay the bounded pilot discovery campaign. Output is exploratory, not a certificate."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE.parent/'filter-budgets'
sys.path.insert(0, str(CORE))
from filters import scores, moment_logs_from_scores
from verify import load_witnesses

ENTROPY = -.25*math.log2(.25)-.75*math.log2(.75)
SEED = 20261015


def filter_scores(r, s, kappa):
    return scores({'r': r, 's': s, 'kappa': kappa})


def moments(sc, eta, tau):
    return np.asarray(moment_logs_from_scores(sc, {'eta': eta, 'tau': tau}))


def parameter_screen(sc, eta, rho, blocks):
    """Historical floating-point screen; final witnesses require the interval checker."""
    if type(blocks) is not int or blocks <= 0 or not 0 < rho < .5 or blocks*rho >= 1:
        return None
    # Fixed exploratory reserves: d=.5-tau+.00002, then two .00002 gaps.
    tau = max((.5+.50002/blocks+(ENTROPY-1+.00004)*rho)/(1+1/blocks),
              (1+.50002/blocks+(2*ENTROPY-1.5+.00004)*rho-eta+.00001)/(1+1/blocks),
              rho+.25, eta/2)+.000001
    if tau >= .5 or eta <= tau:
        return None
    if max(moments(sc, eta, .49999)) >= eta-.0001:
        return None
    lo, hi = tau, .49999
    if max(moments(sc, eta, lo)) > eta-.0001:
        for _ in range(32):
            mid = (lo+hi)/2
            if max(moments(sc, eta, mid)) < eta-.0001:
                hi = mid
            else:
                lo = mid
        tau = hi+.000001
    d = .5-tau+.00002
    delta = d/(blocks*rho)+.00002
    gamma = .5-delta-.00002
    if gamma <= 0 or -(rho*math.log2(rho)+(1-rho)*math.log2(1-rho))+rho+.00001 >= eta:
        return None
    return dict(eta=eta, tau=tau, rho=rho, blocks=blocks, deficiency=d,
                local_deficiency=delta, gamma=gamma, guard=.00001,
                preliminary=.5-d+.000005, moment_logs=moments(sc, eta, tau).tolist())


def scan(sc):
    best = None
    for blocks in range(11, 23):
        for fill in [.93, .97, .995]:
            rho = fill/blocks
            for eta in np.linspace(.510, .550, 81):
                row = parameter_screen(sc, float(eta), rho, blocks)
                if row and (best is None or row['tau'] < best['tau']):
                    best = row
    if best is None:
        return None
    for scale in [.0005, .00005]:
        anchor = best.copy()
        for eta in np.linspace(anchor['eta']-scale, anchor['eta']+scale, 21):
            for rho in np.linspace(anchor['rho']-scale, min(anchor['rho']+scale, .9999/anchor['blocks']), 11):
                row = parameter_screen(sc, float(eta), float(rho), anchor['blocks'])
                if row and row['tau'] < best['tau']:
                    best = row
    return best


def coordinate_search(r, s, kappa, sweeps):
    if len(r) != len(s) or len(r) not in {6, 10} or type(sweeps) is not int or not 0 <= sweeps <= 65:
        raise ValueError('pilot replay requires degree 6 or 10 and at most 65 sweeps')
    rng = np.random.default_rng(SEED)
    x = np.array([*r, *s, kappa], dtype=float)
    width = len(r)

    def objective(z):
        if not .95 < z[-1] < 1.15:
            return float('inf')
        sc = filter_scores(z[:width], z[width:2*width], z[-1])
        return float(max(moments(sc, .51772, .48996)))

    value = objective(x)
    evaluations = 1
    for sweep in range(sweeps):
        step = .035*.94**sweep
        changed = False
        for j in rng.permutation(len(x)):
            for sign in [-1, 1]:
                trial = x.copy()
                trial[j] += sign*(step if j != len(x)-1 else step/8)
                candidate = objective(trial)
                evaluations += 1
                if candidate < value:
                    x, value, changed = trial, candidate, True
        if not changed and sweep > 25:
            break
    return dict(r=x[:width].tolist(), s=x[width:2*width].tolist(), kappa=float(x[-1]),
                objective=value, evaluations=evaluations)


def nelder_mead(objective, start, iterations):
    start = np.asarray(start, dtype=float)
    if start.ndim != 1 or not len(start) or not np.all(np.isfinite(start)):
        raise ValueError('a nonempty finite starting vector is required')
    if type(iterations) is not int or not 0 <= iterations <= 3500:
        raise ValueError('pilot replay allows at most 3500 simplex iterations')
    simplex = [start.copy()]
    for j in range(len(start)):
        x = start.copy()
        x[j] += .025 if j < len(start)-1 else .0015
        simplex.append(x)
    values = [objective(x) for x in simplex]
    calls = len(values)
    for _ in range(iterations):
        order = np.argsort(values)
        simplex, values = [simplex[i] for i in order], [values[i] for i in order]
        center = np.mean(simplex[:-1], axis=0)
        reflected = center+(center-simplex[-1])
        rv = objective(reflected)
        calls += 1
        if values[0] <= rv < values[-2]:
            simplex[-1], values[-1] = reflected, rv
            continue
        if rv < values[0]:
            expanded = center+2*(reflected-center)
            ev = objective(expanded)
            calls += 1
            simplex[-1], values[-1] = (expanded, ev) if ev < rv else (reflected, rv)
            continue
        contracted = center+.5*((reflected if rv < values[-1] else simplex[-1])-center)
        cv = objective(contracted)
        calls += 1
        if cv < min(rv, values[-1]):
            simplex[-1], values[-1] = contracted, cv
            continue
        for j in range(1, len(simplex)):
            simplex[j] = simplex[0]+.5*(simplex[j]-simplex[0])
            values[j] = objective(simplex[j])
            calls += 1
    best = int(np.argmin(values))
    return simplex[best], values[best], calls


def campaign(quick=False):
    started = time.perf_counter()
    witnesses = {w['id']: w for w in load_witnesses()}
    source = witnesses['source']['design']
    r, s = list(map(float, source['r'])), list(map(float, source['s']))
    kappa = float(source['kappa'])
    sweeps, iterations, rounds = (2, 80, 1) if quick else (65, 3500, 3)
    initial = {'baseline_parameters': scan(filter_scores(r, s, kappa)), 'families': []}
    refinement = []
    for degree in [6, 10]:
        row = coordinate_search(r+[0.]*(degree-6), s+[0.]*(degree-6), kappa, sweeps)
        row['degree'] = degree
        row['parameters'] = scan(filter_scores(row['r'], row['s'], row['kappa']))
        if row['parameters'] is None:
            raise ValueError('coordinate stage found no feasible screened parameters')
        initial['families'].append(row)
    for old in initial['families']:
        degree = old['degree']
        x = np.array([*old['r'], *old['s'], old['kappa']])
        parameters = old['parameters']
        for round_number in range(rounds):
            eta, tau = parameters['eta'], parameters['tau']

            def objective(z):
                if not .8 < z[-1] < 1.3 or max(abs(z[:-1])) > 4:
                    return 100.
                logs = moments(filter_scores(z[:degree], z[degree:2*degree], z[-1]), eta, tau)
                return float(np.logaddexp(*(logs/.00002))*.00002)

            x, value, calls = nelder_mead(objective, x, iterations)
            parameters = scan(filter_scores(x[:degree], x[degree:2*degree], x[-1]))
            if parameters is None:
                raise ValueError('simplex stage found no feasible screened parameters')
            refinement.append(dict(degree=degree, round=round_number,
                                   r=x[:degree].tolist(), s=x[degree:2*degree].tolist(),
                                   kappa=float(x[-1]), at_parameters=dict(eta=eta, tau=tau),
                                   objective=value, evaluations=calls, parameters=parameters))
    last = next(row for row in reversed(refinement) if row['degree'] == 10)
    rounded = {'r': [f'{v:.6f}' for v in last['r']], 's': [f'{v:.6f}' for v in last['s']],
               'kappa': f'{last["kappa"]:.6f}'}
    return {'status': 'EXPLORATORY',
            'scope': 'Discovery replay only. Exact retained witnesses require the separate interval certificate.',
            'seed': SEED, 'quick': quick,
            'bounds': {'coordinate_sweeps': sweeps, 'simplex_iterations': iterations,
                       'refinement_rounds_per_degree': rounds, 'grid': 128},
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'inputs_sha256': {name: hashlib.sha256((CORE/name).read_bytes()).hexdigest()
                              for name in ['witnesses.json', 'filters.py', 'SOURCE.json']},
            'replay_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'initial': initial, 'refinement': {'results': refinement},
            'rounded_degree10_design': rounded,
            'matches_retained_design': rounded == witnesses['combined']['design'],
            'seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='short smoke campaign; does not reproduce the retained design')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = campaign(args.quick)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in {'initial', 'refinement'}}, sort_keys=True))


if __name__ == '__main__':
    main()
