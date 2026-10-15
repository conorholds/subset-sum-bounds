"""Exploratory circle quadrature for OpenAI family 138. Not a proof checker."""
import math

import numpy as np


def scores(design, points=128):
    """Return the u, v and M^2 arrays of filtering.tex on a uniform grid."""
    r = np.asarray(design['r'], dtype=float)
    s = np.asarray(design['s'], dtype=float)
    if r.ndim != 1 or len(r) != len(s) or len(r) == 0:
        raise ValueError('coefficient arrays must have equal positive length')
    kappa = float(design['kappa'])
    if not np.all(np.isfinite([*r, *s, kappa])) or kappa <= 0:
        raise ValueError('finite coefficients and positive kappa required')
    x = np.arange(points) * (2 * np.pi / points)
    cosines = np.cos(np.outer(np.arange(1, len(r) + 1), x))
    rv, sv = r @ cosines, s @ cosines
    f, g = np.exp2(rv), np.exp2(sv)
    ff, gg = np.fft.fft(f), np.fft.fft(g)
    p = np.fft.ifft(ff * np.conj(gg)).real / points
    y = np.fft.ifft(np.fft.fft(f * rv) * np.conj(gg)
                    + ff * np.conj(np.fft.fft(g * sv))).real / points
    if np.min(p) <= 0 or not np.all(np.isfinite([*p, *y])):
        raise ValueError('quadrature lost positive finite convolution')
    z = np.log2(f.mean()) + np.log2(g.mean())
    ell, b = np.log2(p) - z, y / p - z
    h = np.log2(np.exp2(kappa * rv).mean()) + np.log2(np.exp2(kappa * sv).mean()) - kappa * z
    return ell + (kappa - 1) * b - h, (2 - kappa) * b - 2 * ell + h, 2 + 2 * np.cos(x)


def moment_logs(design, parameters, points=128):
    return moment_logs_from_scores(scores(design, points), parameters)


def moment_logs_from_scores(sc, parameters):
    u, v, m2 = sc
    eta, tau = float(parameters['eta']), float(parameters['tau'])
    if not tau < eta < 2 * tau:
        raise ValueError('require tau < eta < 2 tau')
    a = (eta-tau) / (2*tau-eta)
    return np.log2([(m2*np.exp2(8.5*(eta-tau-u))).mean(),
                   (m2*np.exp2(8.3*(a*v-u))).mean()]).tolist()


def direct_nodes(design, points=64):
    """Separate scalar O(N^2) control for the FFT correlation convention."""
    r, s = list(map(float, design['r'])), list(map(float, design['s']))
    rv = [sum(c*math.cos(2*math.pi*j*k/points) for k, c in enumerate(r, 1)) for j in range(points)]
    sv = [sum(c*math.cos(2*math.pi*j*k/points) for k, c in enumerate(s, 1)) for j in range(points)]
    f, g = [2**v for v in rv], [2**v for v in sv]
    p = [sum(f[(i+j) % points]*g[i] for i in range(points))/points for j in range(points)]
    y = [sum(f[(i+j) % points]*g[i]*(rv[(i+j) % points]+sv[i]) for i in range(points))/points
         for j in range(points)]
    return p, y
