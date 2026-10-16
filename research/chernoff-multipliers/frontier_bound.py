"""Certified lower bound for the fixed-filter, fixed-reserve parameter frontier.

Convexity in the first multiplier and a positive combination of the list/pair
inequalities give a bound over all eta, tau and rho, not only searched grids.
"""
import argparse
import hashlib
import json
from pathlib import Path

import interval_check as check
import multiplier_search as search

HERE = Path(__file__).resolve().parent


def cost_floor(eta):
    """Necessary tau bound from eliminating rho between two affine costs."""
    h = -check.q('.25') * check.log2(check.q('.25')) - check.q('.75') * check.log2(check.q('.75'))
    slope_left, slope_pair = h - 1 + check.q('.00003'), 2 * h - check.q('1.5') + check.q('.00003')
    check.require(check.upper(slope_left) < 0 < check.lower(slope_pair), 'opposite cost slopes required')
    pair_weight = -slope_left / (slope_pair - slope_left)
    left_weight = 1 - pair_weight
    constant = check.q('.500046') / 17
    left_cost = check.q('.5') + constant + check.q('.000016')
    pair_cost = 1 + constant - eta + check.q('.00001') + check.q('.00001')
    floor = (left_weight * left_cost + pair_weight * pair_cost) / (1 + check.q(1) / 17)
    derivative = -pair_weight / (1 + check.q(1) / 17)
    return floor, derivative


def prove_bound(eta='0.517987190', points=256, digits=70):
    check.require(type(points) is int and points >= 128 and points % 2 == 0,
                  'even bound grid of at least 128 required')
    check.require(type(digits) is int and digits >= 50, 'at least 50 bound digits required')
    check.iv.dps = digits
    eta_exact = check.q(eta)
    tau, slope = cost_floor(eta_exact)
    # Increasing tau reduces the first moment, so this outward upper endpoint
    # gives a conservative lower-bound test at the exact affine cost floor.
    tau_upper = check.endpoint_decimal(tau, 1, precision=45)
    parameters = {'eta': eta, 'tau': tau_upper, 'moment_ceiling': '2'}
    design = search.PREDECESSOR['witness']['design']
    sc = search.scores(design, 512)
    tuned = search.tune_at(sc, parameters)[0]
    multiplier = f'{tuned["multiplier"]:.9f}'
    moments = check.enclose_moments(design, parameters, [multiplier, '8.3'], points, digits)
    value = check.iv.mpf(moments['log2_moment_intervals'][0])
    gradient = check.iv.mpf(moments['first_multiplier_derivative'])
    lam = check.q(multiplier)
    distance = max(check.upper(lam - 1), check.upper(10 - lam))
    # Log moment convexity: f(l) >= f(l0) + (l-l0) f'(l0).
    tangent_floor = check.lower(value) - distance * check.upper(abs(gradient))
    residual = tangent_floor - eta_exact + check.q('.0001')
    check.require(check.lower(residual) > 0, 'first-moment infeasibility not certified at this cutoff')
    check.require(check.upper(slope) < 0, 'cost floor must decrease with eta')
    alpha_floor = tau + check.q('.00001')
    h = -check.q('.25') * check.log2(check.q('.25')) - check.q('.75') * check.log2(check.q('.75'))
    large_block_floor = check.q('.5') - (1 - h) / 19
    check.require(check.lower(large_block_floor) > check.upper(alpha_floor),
                  'large block counts are not excluded by their elementary lower bound')
    return {'status': 'PASS', 'author': 'conorholds',
            'scope': 'Fixed predecessor filters and the explicitly retained reserves; all positive integer block counts.',
            'eta_cutoff': eta, 'cost_floor_tau': check.enclosure(tau),
            'cost_floor_derivative': check.enclosure(slope),
            'alpha_lower_bound': check.enclosure(alpha_floor),
            'tau_lower_bound_for_18_or_more_blocks': check.enclosure(large_block_floor),
            'tangent_multiplier': multiplier, 'allowed_first_multiplier': ['1', '10'],
            'tangent_log2_moment': moments['log2_moment_intervals'][0],
            'tangent_derivative': moments['first_multiplier_derivative'],
            'first_moment_residual_lower_bound': check.enclosure(residual),
            'moments': moments, 'inputs_sha256': check.input_hashes(),
            'bounder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = prove_bound()
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in ['status', 'eta_cutoff', 'alpha_lower_bound',
                                                 'first_moment_residual_lower_bound']}))


if __name__ == '__main__':
    main()
