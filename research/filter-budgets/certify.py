"""Outward-rounded moment bounds and a conditional asymptotic budget audit.

The analytic-trapezoid and split-contour arguments are rederived in README.md.
This checks those numerical sufficient conditions, not the full Subset Sum
algorithm or every proof in the source manuscript.
"""
from decimal import Decimal, InvalidOperation, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import time

from mpmath import iv


def require(condition, message):
    if not condition:
        raise ValueError(message)


def q(value):
    require(isinstance(value, (str, int)) and not isinstance(value, bool),
            'exact decimal strings or integers required')
    try:
        finite = Decimal(value).is_finite()
    except InvalidOperation as error:
        raise ValueError('finite exact decimal required') from error
    require(finite, 'finite exact parameters required')
    return iv.mpf(value)


def upper(x):
    return x.b


def lower(x):
    return x.a


def log2(x):
    return iv.ln(x)/iv.ln(2)


def exp2(x):
    return iv.exp(x*iv.ln(2))


def cosh(x):
    return (iv.exp(x)+iv.exp(-x))/2


def endpoint_decimal(x, side, precision=18):
    sign, mantissa, exponent, _ = x._mpi_[side]
    value = Fraction((-1 if sign else 1)*mantissa)
    value = value * 2**exponent if exponent >= 0 else value / 2**(-exponent)
    with localcontext() as ctx:
        ctx.prec = precision
        ctx.rounding = ROUND_FLOOR if side == 0 else ROUND_CEILING
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def enclosure(x):
    return [endpoint_decimal(x, 0), endpoint_decimal(x, 1)]


def budget_margins(parameters):
    expected = {'eta', 'tau', 'rho', 'blocks', 'deficiency', 'local_deficiency',
                'gamma', 'guard', 'preliminary', 'alpha', 'epsilon', 'moment_ceiling'}
    require(set(parameters) == expected, 'unexpected or missing budget parameters')
    p = {k: q(v) for k, v in parameters.items() if k != 'blocks'}
    blocks = parameters['blocks']
    require(type(blocks) is int and blocks > 0, 'positive integer block count required')
    eta, tau, rho = p['eta'], p['tau'], p['rho']
    d, delta, gamma = p['deficiency'], p['local_deficiency'], p['gamma']
    require(0 < lower(tau) and upper(tau) < lower(eta) and upper(eta) < lower(2*tau),
            'require 0 < tau < eta < 2 tau')
    require(0 < lower(rho) and upper(rho) < q('0.5'), 'require 0 < rho < 1/2')
    require(0 < lower(d) and upper(d) < q('0.5') and 0 < lower(delta)
            and 0 < lower(gamma), 'positive deficiency and representation parameters required')
    require(0 < lower(p['guard']) and 0 < lower(p['epsilon']), 'positive guard and epsilon required')
    hquarter = -q('.25')*log2(q('.25'))-q('.75')*log2(q('.75'))
    hrho = -rho*log2(rho)-(1-rho)*log2(1-rho)
    a = (eta-tau)/(2*tau-eta)
    # The original 30 epsilon score-domination allowance remains adequate.
    rounding = q('30')-q('8.3')*(4*a+3)
    margins = {
        'preliminary': p['preliminary']-(q('.5')-d),
        'block_space': 1-blocks*rho,
        'coverage': blocks*rho*delta-d,
        'representations': q('.5')-delta-gamma,
        'sparse_distinctions': eta-(hrho+rho+p['guard']),
        'generation': tau-(rho+q('.25')),
        'filtered_lists': tau-(q('.5')+(hquarter-q('.5')-gamma)*rho),
        'visited_pairs': tau-(1+(2*hquarter-1-gamma)*rho-eta+p['guard']),
        'rounding_allowance': rounding,
        'moment_reserve': eta-p['moment_ceiling']-48*p['epsilon'],
        'assembly_preliminary': p['alpha']-p['preliminary'],
        'assembly_checksum': p['alpha']-(tau+20*p['epsilon']),
        'assembly_primes': p['alpha']-eta/2,
        'below_meet_in_middle': q('.5')-p['alpha'],
    }
    for name, margin in margins.items():
        require(lower(margin) > 0, 'nonpositive margin: '+name)
    return margins


def certify_moments(design, parameters, points=128, digits=50, beta='0.12', sigma='0.3'):
    require(type(points) is int and points >= 64 and points % 2 == 0, 'even grid of at least 64 required')
    require(type(digits) is int and digits >= 30, 'at least 30 decimal interval digits required')
    iv.dps = digits
    require(set(design) == {'r', 's', 'kappa'} and isinstance(design['r'], list)
            and isinstance(design['s'], list), 'unexpected filter design schema')
    r, s = list(map(q, design['r'])), list(map(q, design['s']))
    require(len(r) == len(s) and len(r) > 0, 'equal nonempty coefficient arrays required')
    kappa, beta, sigma = q(design['kappa']), q(beta), q(sigma)
    require(lower(kappa) > 0 and lower(beta) > 0 and lower(sigma) > 0, 'positive analytic parameters required')
    eta, tau = q(parameters['eta']), q(parameters['tau'])
    require(upper(tau) < lower(eta) and upper(eta) < lower(2*tau), 'invalid moment budgets')
    angles = [2*iv.pi*j/points for j in range(points)]
    cosines = [[iv.cos(k*x) for x in angles] for k in range(1, len(r)+1)]
    rv = [sum((c*cosines[k][j] for k, c in enumerate(r)), q(0)) for j in range(points)]
    sv = [sum((c*cosines[k][j] for k, c in enumerate(s)), q(0)) for j in range(points)]
    f, g = list(map(exp2, rv)), list(map(exp2, sv))

    def strip_sum(cs, width):
        return upper(sum((abs(c)*cosh((j+1)*width) for j, c in enumerate(cs)), q(0)))

    rs, ss = strip_sum(r, sigma), strip_sum(s, sigma)
    denominator = iv.exp(points*sigma)-1

    def with_inner_error(value, bound):
        error = upper(2*bound/denominator)
        return value + iv.mpf([-error.b, error.b])

    zf = with_inner_error(sum(f)/points, exp2(rs))
    zg = with_inner_error(sum(g)/points, exp2(ss))
    zfk = with_inner_error(sum(exp2(kappa*x) for x in rv)/points, exp2(kappa*rs))
    zgk = with_inner_error(sum(exp2(kappa*x) for x in sv)/points, exp2(kappa*ss))
    require(all(lower(z) > 0 for z in [zf, zg, zfk, zgk]), 'inner error lost positive normalizers')
    zlog = log2(zf)+log2(zg)
    h = log2(zfk)+log2(zgk)-kappa*zlog
    delta, a = eta-tau, (eta-tau)/(2*tau-eta)
    first, second = [], []
    node_zero = None
    for j in range(points):
        products = [f[(i+j) % points]*g[i] for i in range(points)]
        p = with_inner_error(sum(products)/points, exp2(rs+ss))
        y = with_inner_error(sum(products[i]*(rv[(i+j) % points]+sv[i]) for i in range(points))/points,
                             (rs+ss)*exp2(rs+ss))
        require(lower(p) > 0, 'inner error lost positive convolution')
        ell, b = log2(p)-zlog, y/p-zlog
        u, v = ell+(kappa-1)*b-h, (2-kappa)*b-2*ell+h
        m2 = 2+2*iv.cos(angles[j])
        first.append(m2*exp2(q('8.5')*(delta-u)))
        second.append(m2*exp2(q('8.3')*(a*v-u)))
        if j == 0:
            node_zero = {'P': enclosure(p), 'Y': enclosure(y)}

    # Shifted L2 norms give a zero-free strip for Q via split contours.
    shifted = []
    for cs, z in [(r, zf), (s, zg)]:
        doubled = [2*c*cosh((j+1)*beta) for j, c in enumerate(cs)]
        values = [exp2(sum((c*cosines[k][j] for k, c in enumerate(doubled)), q(0)))
                  for j in range(points)]
        integral = with_inner_error(sum(values)/points, exp2(strip_sum(doubled, sigma)))
        shifted.append(upper(integral/(z*z)))
    uf, ug = shifted
    theta = upper(iv.sqrt((uf-1)*(ug-1)))
    require(upper(theta) < 1, 'split-contour bound does not certify a zero-free strip')
    qmin = 1-theta
    bmax = upper(iv.sqrt(uf*ug)*(strip_sum(r, beta)+strip_sum(s, beta)+upper(abs(zlog)))/qmin)
    lminus = upper(-log2(qmin))
    e1 = upper(q('8.5')*(delta+lminus+abs(kappa-1)*bmax+upper(h)))
    e2 = upper(q('8.3')*((1+2*a)*lminus+abs(a*(2-kappa)-(kappa-1))*bmax+(1+a)*upper(h)))
    outer_prefactor = 2+2*cosh(2*beta)
    errors = [upper(2*outer_prefactor*exp2(e)/(iv.exp(points*2*beta)-1)) for e in [e1, e2]]
    moments = [sum(values)/points+iv.mpf([-err.b, err.b])
               for values, err in zip([first, second], errors)]
    ceiling = q(parameters['moment_ceiling'])
    for moment in moments:
        require(lower(moment) > 0 and upper(moment) < lower(exp2(ceiling)),
                'continuous moment does not meet the supplied ceiling')
    return {'grid': points, 'interval_digits': digits, 'beta': enclosure(beta), 'sigma': enclosure(sigma),
            'moment_intervals': [enclosure(v) for v in moments],
            'log2_moment_intervals': [enclosure(log2(v)) for v in moments],
            'outer_error_upper': [endpoint_decimal(v, 1) for v in errors],
            'zero_free_theta_upper': endpoint_decimal(theta, 1), 'node_zero': node_zero}


def certify(witness, points=128, digits=50):
    started = time.perf_counter()
    require(type(digits) is int and digits >= 30, 'at least 30 decimal interval digits required')
    iv.dps = digits
    margins = budget_margins(witness['parameters'])
    moments = certify_moments(witness['design'], witness['parameters'], points, digits)
    return {'id': witness['id'], 'status': 'PASS',
            'scope': 'Certified filter moments and conditional sufficient-budget inequalities; not a full algorithm proof.',
            'margins': {k: enclosure(v) for k, v in margins.items()}, 'moments': moments,
            'seconds': time.perf_counter()-started}
