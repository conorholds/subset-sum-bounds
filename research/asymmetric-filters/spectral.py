"""Exploratory Fourier scores and analytic derivatives with cosine and sine modes."""
import numpy as np


def vector_from_design(design):
    if set(design) != {'r', 's', 'r_sin', 's_sin', 'kappa'}:
        raise ValueError('explicit cosine/sine arrays and kappa required')
    arrays = [design[key] for key in ['r', 'r_sin', 's', 's_sin']]
    if not arrays[0] or len({len(array) for array in arrays}) != 1:
        raise ValueError('four nonempty equal-length arrays required')
    return np.array([float(v) for array in arrays for v in array]+[float(design['kappa'])])


def score_state(vector, points=128):
    vector = np.asarray(vector, dtype=float)
    if (vector.ndim != 1 or len(vector) < 5 or (len(vector)-1) % 4
            or not np.all(np.isfinite(vector)) or vector[-1] <= 0):
        raise ValueError('finite cosine/sine arrays and positive kappa required')
    if type(points) is not int or points < 64 or points % 2:
        raise ValueError('even grid of at least 64 required')
    degree = (len(vector)-1)//4
    r, s, kappa = vector[:2*degree], vector[2*degree:4*degree], vector[-1]
    x = np.arange(points)*(2*np.pi/points)
    angles = np.outer(np.arange(1, degree+1), x)
    basis = np.concatenate((np.cos(angles), np.sin(angles)))
    rv, sv = r@basis, s@basis
    ln2 = np.log(2.)

    def corr(left, right):
        return np.fft.ifft(np.fft.fft(left)*np.conj(np.fft.fft(right))).real/points

    with np.errstate(over='raise', invalid='raise', divide='raise'):
        f, g = np.exp2(rv), np.exp2(sv)
        p = corr(f, g)
        if p.min() <= 0:
            raise ValueError('quadrature lost positive convolution')
        y = corr(f*rv, g)+corr(f, g*sv)
        z = np.log2(f.mean())+np.log2(g.mean())
        fk, gk = np.exp2(kappa*rv), np.exp2(kappa*sv)
        h = np.log2(fk.mean())+np.log2(gk.mean())-kappa*z
        ell, b = np.log2(p)-z, y/p-z
        df, dg = ln2*f*basis, ln2*g*basis
        dp = np.concatenate((corr(df, g), corr(f, dg), np.zeros((1, points))))
        dy = np.concatenate((corr(f*basis*(1+ln2*rv), g)+corr(df, g*sv),
                             corr(f*rv, dg)+corr(f, g*basis*(1+ln2*sv)), np.zeros((1, points))))
        dz = np.concatenate(((f*basis).mean(axis=1)/f.mean(),
                             (g*basis).mean(axis=1)/g.mean(), [0.]))
        dh = kappa*np.concatenate(((fk*basis).mean(axis=1)/fk.mean(),
                                   (gk*basis).mean(axis=1)/gk.mean(), [0.]))-kappa*dz
        dh[-1] = (fk*rv).mean()/fk.mean()+(gk*sv).mean()/gk.mean()-z
        dell, db = dp/(ln2*p)-dz[:, None], dy/p-y*dp/p**2-dz[:, None]
        du, dv = dell+(kappa-1)*db-dh[:, None], (2-kappa)*db-2*dell+dh[:, None]
        du[-1] += b
        dv[-1] -= b
        return {'u': ell+(kappa-1)*b-h, 'v': (2-kappa)*b-2*ell+h,
                'du': du, 'dv': dv, 'm2': 2+2*np.cos(x)}


def scores(design, points=256):
    state = score_state(vector_from_design(design), points)
    return state['u'], state['v'], state['m2']
