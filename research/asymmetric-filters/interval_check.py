"""Continuous interval bounds for real cosine/sine filters and fractional coverage.

Adapted from filter-budgets/certify.py: the two Chernoff multipliers are
explicit exact inputs in both integrands and analytic remainder bounds.
The original frozen checker and all its receipts remain unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import mpmath
from mpmath import iv

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / 'filter-budgets'
sys.path.insert(0, str(CORE))
from certify import (cosh, enclosure, endpoint_decimal, exp2,
                     log2, lower, q, require, upper)

import importlib.util
SPEC = importlib.util.spec_from_file_location('asymmetric_fractional_budget', HERE.parent / 'fractional-coverage/budget.py')
FRACTIONAL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FRACTIONAL)
budget_margins = FRACTIONAL.budget_margins


def enclose_moments(design, parameters, multipliers, points=128, digits=50, beta='0.12', sigma='0.3'):
    require(type(points) is int and points >= 64 and points % 2 == 0, 'even grid of at least 64 required')
    require(type(digits) is int and digits >= 30, 'at least 30 decimal interval digits required')
    iv.dps = digits
    lam1, lam2 = map(q, multipliers)
    require(lower(lam1) >= 1 and lower(lam2) >= 1, 'domination multipliers must be at least one')
    require(set(design) == {'r', 's', 'r_sin', 's_sin', 'kappa'} and isinstance(design['r'], list)
            and isinstance(design['s'], list), 'unexpected filter design schema')
    require(isinstance(design['r_sin'], list) and isinstance(design['s_sin'], list), 'sine arrays required')
    r, s, rsin, ssin = [list(map(q, design[key])) for key in ['r', 's', 'r_sin', 's_sin']]
    require(len(r) == len(s) == len(rsin) == len(ssin) and len(r) > 0, 'equal nonempty coefficient arrays required')
    rh, sh = list(zip(r, rsin)), list(zip(s, ssin))
    kappa, beta, sigma = q(design['kappa']), q(beta), q(sigma)
    require(lower(kappa) > 0 and lower(beta) > 0 and lower(sigma) > 0, 'positive analytic parameters required')
    eta, tau = q(parameters['eta']), q(parameters['tau'])
    require(upper(tau) < lower(eta) and upper(eta) < lower(2*tau), 'invalid moment budgets')
    angles = [2*iv.pi*j/points for j in range(points)]
    cosines = [[iv.cos(k*x) for x in angles] for k in range(1, len(r)+1)]
    sines = [[iv.sin(k*x) for x in angles] for k in range(1, len(r)+1)]
    rv = [sum((c*cosines[k][j]+d*sines[k][j] for k, (c,d) in enumerate(rh)), q(0)) for j in range(points)]
    sv = [sum((c*cosines[k][j]+d*sines[k][j] for k, (c,d) in enumerate(sh)), q(0)) for j in range(points)]
    f, g = list(map(exp2, rv)), list(map(exp2, sv))

    def strip_sum(cs, width):
        return upper(sum((iv.sqrt(c*c+d*d)*cosh((j+1)*width) for j, (c,d) in enumerate(cs)), q(0)))

    rs, ss = strip_sum(rh, sigma), strip_sum(sh, sigma)
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
    first, second, first_derivative = [], [], []
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
        first.append(m2*exp2(lam1*(delta-u)))
        first_derivative.append(first[-1]*(delta-u))
        second.append(m2*exp2(lam2*(a*v-u)))
        if j == 0:
            node_zero = {'P': enclosure(p), 'Y': enclosure(y)}

    # Shifted L2 norms give a zero-free strip for Q via split contours.
    shifted = []
    for cs, z in [(rh, zf), (sh, zg)]:
        doubled = [(2*c*cosh((j+1)*beta), 2*d*cosh((j+1)*beta)) for j, (c,d) in enumerate(cs)]
        values = [exp2(sum((c*cosines[k][j]+d*sines[k][j] for k, (c,d) in enumerate(doubled)), q(0)))
                  for j in range(points)]
        integral = with_inner_error(sum(values)/points, exp2(strip_sum(doubled, sigma)))
        shifted.append(upper(integral/(z*z)))
    uf, ug = shifted
    theta = upper(iv.sqrt((uf-1)*(ug-1)))
    require(upper(theta) < 1, 'split-contour bound does not certify a zero-free strip')
    qmin = 1-theta
    bmax = upper(iv.sqrt(uf*ug)*(strip_sum(rh, beta)+strip_sum(sh, beta)+upper(abs(zlog)))/qmin)
    lminus = upper(-log2(qmin))
    e1 = upper(lam1*(delta+lminus+abs(kappa-1)*bmax+upper(h)))
    e2 = upper(lam2*((1+2*a)*lminus+abs(a*(2-kappa)-(kappa-1))*bmax+(1+a)*upper(h)))
    outer_prefactor = 2+2*cosh(2*beta)
    errors = [upper(2*outer_prefactor*exp2(e)/(iv.exp(points*2*beta)-1)) for e in [e1, e2]]
    moments = [sum(values)/points+iv.mpf([-err.b, err.b])
               for values, err in zip([first, second], errors)]
    # d/d lambda log2 I1 = integral(c1 M^2 2^(lambda c1)) / I1.
    # |log2 Q| <= -log2(1-theta) follows from its disk power series.
    c1_bound = upper(delta+lminus+abs(kappa-1)*bmax+upper(abs(h)))
    derivative_error = upper(errors[0]*c1_bound)
    derivative_numerator = sum(first_derivative)/points + iv.mpf([-derivative_error.b, derivative_error.b])
    ceiling = q(parameters['moment_ceiling'])
    for moment in moments:
        require(lower(moment) > 0 and upper(moment) < lower(exp2(ceiling)),
                'continuous moment does not meet the supplied ceiling')
    return {'grid': points, 'interval_digits': digits, 'beta': enclosure(beta), 'sigma': enclosure(sigma),
            'moment_intervals': [enclosure(v) for v in moments],
            'log2_moment_intervals': [enclosure(log2(v)) for v in moments],
            'outer_error_upper': [endpoint_decimal(v, 1) for v in errors],
            'first_multiplier_derivative': enclosure(derivative_numerator/moments[0]),
            'zero_free_theta_upper': endpoint_decimal(theta, 1), 'node_zero': node_zero}


def certify_candidate(candidate, points=128, digits=50):
    started = time.perf_counter()
    require(type(digits) is int and digits >= 30, 'at least 30 interval digits required')
    iv.dps = digits
    require(isinstance(candidate, dict) and set(candidate) == {'id', 'design', 'parameters', 'multipliers'},
            'unexpected or missing candidate fields')
    require(isinstance(candidate['id'], str) and candidate['id'], 'nonempty candidate id required')
    multipliers = candidate['multipliers']
    require(isinstance(multipliers, list) and len(multipliers) == 2, 'two exact multipliers required')
    lam1, lam2 = map(q, multipliers)
    require(lower(lam1) >= 1 and lower(lam2) >= 1, 'domination multipliers must be at least one')
    parameters = candidate['parameters']
    margins = budget_margins(parameters)
    eta, tau = q(parameters['eta']), q(parameters['tau'])
    ratio = (eta - tau) / (2 * tau - eta)
    margins['domination_one'] = 30 - 3 * lam1
    margins['domination_two'] = 30 - (4 * ratio + 3) * lam2
    for name in ['domination_one', 'domination_two']:
        require(lower(margins[name]) > 0, 'nonpositive margin: ' + name)
    moments = enclose_moments(candidate['design'], parameters, multipliers, points, digits)
    return {'id': candidate['id'], 'status': 'PASS', 'multipliers': multipliers,
            'margins': {key: enclosure(value) for key, value in margins.items()},
            'moments': moments, 'seconds': time.perf_counter() - started}


def input_hashes():
    sources = {'package': (HERE, ['interval_check.py', 'spectral.py', 'search.py', 'candidates.json', 'predecessor.json']),
               'fractional': (HERE.parent/'fractional-coverage', ['budget.py']),
               'shared': (CORE, ['certify.py', 'SOURCE.json', 'requirements.txt'])}
    return {group: {name: hashlib.sha256((path/name).read_bytes()).hexdigest() for name in names}
            for group, (path,names) in sources.items()}


def replay(points=128, digits=50, author='conorholds'):
    require(isinstance(author, str) and author, 'explicit nonempty author required')
    data = json.loads((HERE / 'candidates.json').read_text())
    require(set(data) == {'model', 'candidates'} and data['model'] == 'asymmetric-filters-v1',
            'unrecognized model')
    require(isinstance(data['candidates'], list) and data['candidates'], 'nonempty candidates required')
    ids = [candidate['id'] for candidate in data['candidates']]
    require(len(set(ids)) == len(ids), 'duplicate candidate ids')
    return {'status': 'PASS', 'model': data['model'], 'author': author,
            'scope': 'Cosine/sine moments and fractional-coverage conditional budgets; independent algorithm review pending.',
            'inputs_sha256': input_hashes(),
            'runtime': {'python': platform.python_version(), 'mpmath': mpmath.__version__},
            'results': [certify_candidate(candidate, points, digits) for candidate in data['candidates']]}


def validate_receipt(receipt):
    require(isinstance(receipt, dict) and receipt.get('status') == 'PASS'
            and receipt.get('model') == 'asymmetric-filters-v1', 'unrecognized or failing receipt')
    require(receipt.get('inputs_sha256') == input_hashes(), 'receipt inputs have changed')
    candidates = json.loads((HERE / 'candidates.json').read_text())['candidates']
    expected = {candidate['id']: candidate['multipliers'] for candidate in candidates}
    rows = receipt.get('results')
    require(isinstance(rows, list) and len(rows) == len(expected), 'incomplete receipt coverage')
    require(all(isinstance(row, dict) and isinstance(row.get('id'), str) for row in rows),
            'invalid receipt rows')
    require(len({row['id'] for row in rows}) == len(rows)
            and {row['id'] for row in rows} == set(expected), 'receipt candidate coverage mismatch')
    configurations = set()
    for row in rows:
        require(row.get('status') == 'PASS' and row.get('multipliers') == expected[row['id']],
                'receipt status or multipliers mismatch')
        moments = row.get('moments')
        require(isinstance(moments, dict), 'missing receipt moments')
        points, digits = moments.get('grid'), moments.get('interval_digits')
        require(type(points) is int and points >= 64 and points % 2 == 0
                and type(digits) is int and digits >= 30, 'invalid receipt precision')
        configurations.add((points, digits))
    require(len(configurations) == 1, 'one receipt configuration required')
    return configurations.pop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--author', required=True)
    parser.add_argument('--points', type=int)
    parser.add_argument('--digits', type=int)
    parser.add_argument('--check-receipt', type=Path)
    args = parser.parse_args()
    points, digits = 128, 50
    if args.check_receipt:
        receipt = json.loads(args.check_receipt.read_text())
        points, digits = validate_receipt(receipt)
    if args.points is not None:
        points = args.points
    if args.digits is not None:
        digits = args.digits
    report = replay(points, digits, args.author)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'model': report['model'],
                      'results': [{'id': row['id'], 'log2_moments': row['moments']['log2_moment_intervals'],
                                   'seconds': round(row['seconds'], 3)} for row in report['results']]}))


if __name__ == '__main__':
    main()
