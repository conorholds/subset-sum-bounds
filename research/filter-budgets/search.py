"""Bounded deterministic-seed filter search; output requires separate certification."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from filters import moment_logs
from verify import load_witnesses


def search(parameters, steps=4000, degree=10, seed=20261015, start='source'):
    if not 6 <= degree <= 24 or not 0 <= steps <= 100000:
        raise ValueError('require degree in [6,24] and steps in [0,100000]')
    source = next(w for w in load_witnesses() if w['id'] == start)['design']
    width = len(source['r'])
    if degree < width:
        raise ValueError('degree cannot truncate the starting filter')
    initial = np.array([*map(float, source['r']), *([0.]*(degree-width)),
                        *map(float, source['s']), *([0.]*(degree-width)), float(source['kappa'])])

    def design(x):
        return {'r': [f'{v:.9f}' for v in x[:degree]],
                's': [f'{v:.9f}' for v in x[degree:2*degree]], 'kappa': f'{x[-1]:.9f}'}

    def objective(x):
        return max(moment_logs(design(x), parameters))

    rng = np.random.default_rng(seed)
    best = initial.copy()
    before = value = objective(best)
    started = time.perf_counter()
    for iteration in range(steps):
        trial = best.copy()
        count = 1 if iteration % 3 else 3
        indices = rng.choice(len(trial), size=count, replace=False)
        scale = .035*(.2**(iteration/max(1, steps)))
        trial[indices] += rng.normal(size=count)*scale
        trial[-1] = np.clip(trial[-1], .8, 1.3)
        candidate = objective(trial)
        if candidate < value:
            best, value = trial, candidate
    return {'status': 'EXPLORATORY', 'scope': 'Floating-point quadrature; no proof or global exponent claim.',
            'seed': seed, 'degree': degree, 'evaluations': steps+1, 'start': start,
            'parameters': parameters, 'starting_max_log2_moment': before,
            'candidate_max_log2_moment': value, 'design': design(best),
            'seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', type=int, default=4000)
    parser.add_argument('--degree', type=int, default=10)
    parser.add_argument('--seed', type=int, default=20261015)
    parser.add_argument('--start', choices=['source', 'combined'], default='source')
    parser.add_argument('--parameters', choices=['source', 'parameters_only', 'combined'], default='source')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    parameters = next(w for w in load_witnesses() if w['id'] == args.parameters)['parameters']
    result = search(parameters, args.steps, args.degree, args.seed, args.start)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
