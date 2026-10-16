"""Bounded exploratory test of odd Fourier directions; no interval certificate."""
import json
import argparse
import hashlib
from pathlib import Path
import sys
import time
import numpy as np

PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE))
import multiplier_search as tune

candidate = json.loads((PACKAGE / 'joint-degree20.json').read_text())['candidate']
design = candidate['design']
degree, points = len(design['r']), 256
x = np.arange(points) * (2 * np.pi / points)
indices = np.arange(1, degree + 1)
cosines, sines = np.cos(np.outer(indices, x)), np.sin(np.outer(indices, x))
rv0, sv0 = np.array(design['r'], float) @ cosines, np.array(design['s'], float) @ cosines
kappa = float(design['kappa'])
lam1, lam2 = map(float, candidate['multipliers'])
m2 = 2 + 2 * np.cos(x)
tau_prime = (tune.H - 1 + .00003) / ((tune.H - .5) * (1 + 1 / 17))
tau_constant = tune.slice_parameters(.51799)['tau'] - tau_prime * .51799
calls = 0


def boundaries(odd):
    global calls
    calls += 1
    rv, sv = rv0 + odd[:degree] @ sines, sv0 + odd[degree:] @ sines
    f, g = np.exp2(rv), np.exp2(sv)
    ff, gg = np.fft.fft(f), np.fft.fft(g)
    p = np.fft.ifft(ff * np.conj(gg)).real / points
    y = np.fft.ifft(np.fft.fft(f * rv) * np.conj(gg) + ff * np.conj(np.fft.fft(g * sv))).real / points
    z = np.log2(f.mean()) + np.log2(g.mean())
    ell, b = np.log2(p) - z, y / p - z
    h = np.log2(np.exp2(kappa * rv).mean()) + np.log2(np.exp2(kappa * sv).mean()) - kappa * z
    u, v = ell + (kappa - 1) * b - h, (2 - kappa) * b - 2 * ell + h
    constant = tune.log_moment(-u, m2, lam1)
    first = (lam1 * tau_constant - constant - .0001) / (lam1 * (1 - tau_prime) - 1)
    lo, hi = .5175, .51828
    for _ in range(40):
        eta = (lo + hi) / 2
        tau = tau_constant + tau_prime * eta
        ratio = (eta - tau) / (2 * tau - eta)
        residual = tune.log_moment(ratio * v - u, m2, lam2) - eta + .0001
        if residual < 0:
            lo = eta
        else:
            hi = eta
    return np.array([first, (lo + hi) / 2])


def main(output):
    global calls
    calls = 0
    started = time.perf_counter()
    dimension, step = 2 * degree, .001
    zero = np.zeros(dimension)
    origin = boundaries(zero)
    hessians = np.zeros((2, dimension, dimension))
    for i in range(dimension):
        e = zero.copy(); e[i] = step
        hessians[:, i, i] = 2 * (boundaries(e) - origin) / step**2
        for j in range(i):
            plus, minus = e.copy(), e.copy()
            plus[j], minus[j] = step, -step
            value = (boundaries(plus) - boundaries(minus)) / (2 * step**2)
            hessians[:, i, j] = hessians[:, j, i] = value
    def top(t):
        values, vectors = np.linalg.eigh(t * hessians[0] + (1 - t) * hessians[1])
        return values[-1], vectors[:, -1]
    lo, hi = 0., 1.
    for _ in range(45):
        a, b = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if top(a)[0] < top(b)[0]: hi = b
        else: lo = a
    weight = (lo + hi) / 2
    eigenvalue, direction = top(weight)
    trials = []
    for amplitude in [.001, .003, .01, .03, .1]:
        result = boundaries(amplitude * direction)
        trials.append({'amplitude': amplitude, 'boundaries': result.tolist(),
                       'worst_eta_change': float(min(result) - min(origin))})
    report = {'status': 'EXPLORATORY', 'points': points, 'degree': degree, 'hessian_step': step,
              'origin': origin.tolist(), 'weight': weight, 'largest_weighted_eigenvalue': float(eigenvalue),
              'direction': direction.tolist(),
              'direction_curvatures': [float(direction @ h @ direction) for h in hessians],
              'trials': trials, 'calls': calls,
              'inputs_sha256': {name: hashlib.sha256((PACKAGE / name).read_bytes()).hexdigest()
                                for name in ['odd_probe.py', 'multiplier_search.py', 'joint-degree20.json']}, 'seconds': time.perf_counter() - started}
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)
