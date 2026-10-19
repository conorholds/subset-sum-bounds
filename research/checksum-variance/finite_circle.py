"""Exact conservative finite-circle cutoff from elementary Lipschitz bounds.

This binds an existing continuous moment ceiling; it does not re-prove that
ceiling. Every retained numerical input is a rational decimal string.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def ceil(value):
    return -(-value.numerator//value.denominator)


def exact(value):
    require(isinstance(value, str) and len(value) <= 128
            and re.fullmatch(r'-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)', value) is not None,
            'bounded exact decimal strings required')
    return Fraction(value)


def polynomial_bounds(cosines, sines):
    require(isinstance(cosines, list) and isinstance(sines, list)
            and 1 <= len(cosines) <= 32 and len(cosines) == len(sines),
            'one to thirty-two matched real harmonics required')
    amplitudes = [abs(exact(c))+abs(exact(s)) for c, s in zip(cosines, sines)]
    return sum(amplitudes, Fraction(0)), sum((j*x for j, x in enumerate(amplitudes, 1)), Fraction(0))


def cutoff(candidate):
    require(isinstance(candidate, dict) and set(candidate) == {'id', 'design', 'parameters', 'multipliers'},
            'invalid exact candidate schema')
    design = candidate['design']
    require(isinstance(design, dict) and set(design) == {'r', 'r_sin', 's', 's_sin', 'kappa'},
            'full real Fourier schema required')
    r, rd = polynomial_bounds(design['r'], design['r_sin'])
    s, sd = polynomial_bounds(design['s'], design['s_sin'])
    b, derivative = r+s, rd+sd
    kappa = exact(design['kappa'])
    require(0 < kappa <= 4 and b <= 32 and derivative <= 512,
            'cutoff calculator domain: kappa at most four, amplitude at most 32, derivative sum at most 512')
    p = candidate['parameters']
    eta, tau, epsilon, chi = [exact(p[key]) for key in ['eta', 'tau', 'epsilon', 'moment_ceiling']]
    require(0 < tau < eta < 2*tau and eta <= 1 and epsilon > 0 and 0 < chi < eta,
            'invalid checksum domain')
    reserve = eta-chi-48*epsilon
    require(reserve > 0, 'positive remaining variance reserve required')
    require(isinstance(candidate['multipliers'], list) and len(candidate['multipliers']) == 2,
            'two exact multiplier strings required')
    lam1, lam2 = map(exact, candidate['multipliers'])
    require(0 < lam1 <= 16 and 0 < lam2 <= 16, 'cutoff calculator slopes must be positive and at most sixteen')
    delta, ratio = eta-tau, (eta-tau)/(2*tau-eta)
    require(ratio <= 1, 'cutoff calculator requires a at most one')
    # Uniform absolute bounds for continuous u and v.
    u_bound = b*(abs(kappa-1)+1+kappa)
    v_bound = b*(abs(2-kappa)+2+kappa)
    exponent_bound = ceil(max(lam1*(delta+u_bound), lam2*(ratio*v_bound+u_bound)))
    # Uniform x-derivative bound. Only r(x+w) depends on x.
    t_derivative = rd*(1+2*b)
    u_derivative = abs(kappa-1)*t_derivative+rd
    v_derivative = abs(2-kappa)*t_derivative+2*rd
    exponent_derivative = max(lam1*u_derivative, lam2*(ratio*v_derivative+u_derivative))
    # Uniform finite-inner score discrepancy at the common outer grid nodes.
    u_error = abs(kappa-1)*(1+2*b)+1+kappa
    v_error = abs(2-kappa)*(1+2*b)+2+kappa
    inner_log_error = 8*derivative*max(lam1*u_error, lam2*(ratio*v_error+u_error))
    # C bounds the relative error of the continuous outer integrand's grid.
    outer_relative_error = (1 << (2*exponent_bound+1))*(2+4*exponent_derivative)
    total_log_error = inner_log_error+2*outer_relative_error
    required_modulus = max(2, ceil(16*derivative), ceil(2*total_log_error/reserve))
    power = (required_modulus-1).bit_length()
    modulus = 1 << power
    n_cutoff = ceil(Fraction(power+2)/eta)
    require(Fraction(total_log_error, modulus) <= reserve/2, 'dyadic cutoff misses moment reserve')
    require(Fraction(8*derivative, modulus) <= Fraction(1, 2), 'inner cell distortion exceeds its domain')
    require(eta*n_cutoff-2 >= power, 'input-size conversion misses modulus cutoff')
    fields = {'r_amplitude': r, 's_amplitude': s, 'r_derivative': rd, 's_derivative': sd,
              'exponent_amplitude_integer': exponent_bound, 'exponent_derivative': exponent_derivative,
              'inner_log_error_constant': inner_log_error, 'outer_relative_error_constant': outer_relative_error,
              'total_log_error_constant': total_log_error, 'remaining_variance_reserve': reserve,
              'finite_ceiling': chi+reserve/2, 'minimum_modulus_power': power,
              'sufficient_input_size': n_cutoff}
    return {key: str(value) if isinstance(value, Fraction) else value for key, value in fields.items()}


def report(author):
    require(author == 'conorholds', 'approved author handle required')
    data = json.loads((HERE/'FINITE-INPUT.json').read_text())
    return {'author': author, 'status': 'PASS',
            'scope': 'Effective finite-circle transfer of the pinned continuous ceiling; excludes coverage and full algorithm cutoffs.',
            'inputs_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                              for name in ['finite_circle.py', 'FINITE-INPUT.json']},
            'candidate_id': data['candidate']['id'], 'cutoff': cutoff(data['candidate'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = report(args.author)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'candidate': result['candidate_id'],
                      'minimum_modulus_power': result['cutoff']['minimum_modulus_power'],
                      'sufficient_input_size': result['cutoff']['sufficient_input_size']}))


if __name__ == '__main__':
    main()
