"""A fixed-filter lower bound over all parameters in the stated reserve class.

Eliminate rho between the two costs, then support the second log moment by
a convex plane in (lambda*a, lambda). Derivatives use rigorous convex secants,
not unvalidated finite differences.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import budget
import importlib.util
SPEC = importlib.util.spec_from_file_location('fractional_frontier_checker', HERE/'interval_check.py')
check = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(check)
sys.path.insert(0, str(budget.CHERNOFF))
import multiplier_search as tune

q, require, lower, upper = budget.q, budget.require, budget.lower, budget.upper
iv = check.CHERNOFF.iv
enclosure, endpoint_decimal = check.CHERNOFF.enclosure, check.CHERNOFF.endpoint_decimal


def cost_floor(eta):
    h = -q('.25')*budget.log2(q('.25'))-q('.75')*budget.log2(q('.75'))
    a, b = h-q('.5')+q('.00005'), 2*h-1+q('.00005')
    c1, c2 = q('.5000176'), q('1.0000216')
    denominator = eta+c1+b-c2-a
    require(lower(denominator) > 0, 'positive cost-floor denominator required')
    numerator = a*eta+c1*b-a*c2
    derivative_numerator = a*(c1+b-c2-a)-c1*b+a*c2
    require(upper(derivative_numerator) < 0, 'cost floor must strictly decrease')
    return numerator/denominator, derivative_numerator/denominator**2


def prove_bound(eta='0.517990280', points=256, digits=70, author='conorholds'):
    require(type(points) is int and points >= 128 and points % 2 == 0,
            'even bound grid of at least 128 required')
    require(type(digits) is int and digits >= 50, 'at least 50 bound digits required')
    require(isinstance(author, str) and author, 'explicit nonempty author required')
    iv.dps = digits
    eta0 = q(eta)
    tau, slope = cost_floor(eta0)
    tau_probe = endpoint_decimal(tau, 1, precision=45)
    parameters = {'eta': eta, 'tau': tau_probe, 'moment_ceiling': '2'}
    predecessor = json.loads((HERE/'predecessor.json').read_text())['candidate']
    design = predecessor['design']
    sc = tune.scores(design, 512)
    multiplier = f'{tune.tune_at(sc, parameters)[1]["multiplier"]:.12f}'
    lam, h = q(multiplier), q('.0000001')

    def evaluate(p, l):
        result = check.CHERNOFF.enclose_moments(design, p, ['8.5', l], points, digits)
        return iv.mpf(result['log2_moment_intervals'][1]), result

    value, moments = evaluate(parameters, multiplier)
    minus_lambda = endpoint_decimal(lam-h, 0, precision=45)
    plus_lambda = endpoint_decimal(lam+h, 1, precision=45)
    minus, minus_moments = evaluate(parameters, minus_lambda)
    plus, plus_moments = evaluate(parameters, plus_lambda)
    left_derivative = (value-minus)/(lam-q(minus_lambda))
    right_derivative = (plus-value)/(q(plus_lambda)-lam)
    derivative_abs = max(upper(abs(left_derivative)), upper(abs(right_derivative)))
    distance = max(upper(lam-1), upper(10-lam))

    # Increasing tau lowers a. A backward secant of the convex function in a
    # is a lower bound on lambda*E_tilt[v] at the supporting point.
    tau_minus = endpoint_decimal(q(tau_probe)+q('.00001'), 1, precision=45)
    smaller_a_parameters = dict(parameters, tau=tau_minus)
    smaller_value, smaller_moments = evaluate(smaller_a_parameters, multiplier)
    a0 = (eta0-q(tau_probe))/(2*q(tau_probe)-eta0)
    aminus = (eta0-q(tau_minus))/(2*q(tau_minus)-eta0)
    require(lower(a0-aminus) > 0, 'positive a-secant displacement required')
    v_lower = lower((lower(value)-upper(smaller_value))/(lam*(a0-aminus)))
    monotone = v_lower*q(tau_probe)/(2*q(tau_probe)-eta0)**2-1
    residual = lower(value)-distance*derivative_abs-eta0+q('.0001')
    require(lower(v_lower) > 0 and lower(monotone) > 0,
            'second-moment plane does not control all larger eta')
    require(lower(residual) > 0, 'second-moment infeasibility not certified at this cutoff')
    return {'status': 'PASS', 'author': author,
            'scope': 'Fixed predecessor filters, all eta/tau/rho and multipliers, with the listed minimum reserves.',
            'minimum_reserves': {'deficiency': '.00002', 'fractional_coverage_per_n': '.0000016',
                                 'representation': '.00003', 'filtered_lists': '.000016',
                                 'visited_pairs': '.000010', 'guard': '.00001',
                                 'moment': '.0001', 'assembly': '.00001'},
            'eta_cutoff': eta, 'cost_floor_tau': enclosure(tau), 'tau_probe': tau_probe,
            'cost_floor_derivative': enclosure(slope), 'alpha_lower_bound': enclosure(tau+q('.00001')),
            'support_multiplier': multiplier, 'allowed_second_multiplier': ['1', '10'],
            'lambda_secant_step': '.0000001',
            'lambda_derivative_enclosure': [endpoint_decimal(left_derivative, 0), endpoint_decimal(right_derivative, 1)],
            'tilted_v_lower_bound': enclosure(v_lower), 'larger_eta_monotonicity_margin': enclosure(monotone),
            'second_moment_residual_lower_bound': enclosure(residual),
            'support_moments': moments,
            'secant_moments': {'minus_lambda': minus_moments, 'plus_lambda': plus_moments,
                               'smaller_a': smaller_moments},
            'inputs_sha256': check.input_hashes(),
            'bounder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = prove_bound(author=args.author)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in ['status', 'alpha_lower_bound',
                                                 'second_moment_residual_lower_bound',
                                                 'larger_eta_monotonicity_margin']}))
