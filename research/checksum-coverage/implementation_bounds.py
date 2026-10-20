"""Explicit checksum-only draw coupling, rounding and exponential overhead.

The bounds assume the source's deterministic score evaluator and phase
approximations. They do not implement a full estimator or bound its complete
polynomial operation count, alias extraction, prime sampling or final assembly.
"""
import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path

import coverage_bounds as coverage

HERE = Path(__file__).resolve().parent


def overhead_envelope(c):
    """H(n)<=a sqrt(n)+b log2(n)+constant for the minimum-block rule."""
    eps, g = c['epsilon'], c['pattern_bits_per_block']
    table_slope = c['d_bound']+4*eps
    return {'sqrt_coefficient': 2*g+table_slope,
            'log_coefficient': 1+table_slope/eps,
            'constant': 2+5*table_slope+2*table_slope*c['minimum_block_size']
                        +table_slope*(g+9)/eps}


def overhead_cutoff(c):
    envelope = overhead_envelope(c)
    for power in range(2, 513):
        n = 1 << power
        root = isqrt(n)
        root += root*root != n
        upper = (envelope['sqrt_coefficient']*root
                 +envelope['log_coefficient']*power+envelope['constant'])
        if upper <= c['epsilon']*n:
            return {'power': power, 'n': n, 'envelope': envelope,
                    'envelope_at_cutoff_at_most': upper,
                    'allowed_overhead': c['epsilon']*n}
    raise ValueError('overhead cutoff outside the bounded calculator domain')


def finite_conditions(candidate, power, target_degree=2, failure_bits=130):
    coverage.require(type(power) is int and 1 <= power <= 512, 'dyadic power from one to 512 required')
    coverage.require(type(target_degree) is int and 1 <= target_degree <= 1024,
                     'integer target degree from one to 1024 required')
    coverage.require(type(failure_bits) is int and 2 <= failure_bits <= 512,
                     'failure bits must be an integer from two to 512')
    c = coverage.constants(candidate)
    n, g = 1 << power, c['pattern_bits_per_block']
    partition = coverage.minimum_block_partition(c, n)
    eps = c['epsilon']
    tau = coverage.exact(candidate['parameters']['tau'])
    reserve = c['eta']-coverage.exact(candidate['parameters']['moment_ceiling'])-48*eps
    sigma = tau-reserve/4
    table_slope = c['d_bound']+4*eps
    # Phases n; labels <=2*n*Q*Rmax; matches <=T*Q*Ndraw.
    # Rmax<=2^(table_slope*dmax+1), Ndraw<=2^(sigma*n+1).
    # For T<=n^d, their sum is <=7*n^d*2^(g*q+max(...)).
    calls_bits = (g*partition['q']+coverage.ceil(max(table_slope*partition['maximum'], sigma*n))
                  +3+target_degree*coverage.ceil_log2(n))
    scalar = coverage.ceil(max(table_slope, sigma))
    # A global all-history bound, using q,dmax<=n and ceil(log2 n)<=n.
    linear = g+scalar+target_degree
    coupling_checks = {'cubic_dominates_linear': n*n >= 2*linear,
                       'coupling_tail': n**3 >= 2*(failure_bits+3)}
    cap_exponent = coverage.ceil((tau+20*eps)*n)
    # Both bounds hold for each target on every nonaborted history.
    evaluation_bits = (2*n*n-2*n-3*cap_exponent-g*partition['q']
                       -coverage.ceil_log2(n+1))
    rounding_bits = n*n-cap_exponent-g*partition['q']-1
    # Proposed smaller budgets. They depend only on public n/parameters,
    # retain the same error reserves, and do not condition the phase law.
    linear_phase_bits = 2*n+3*cap_exponent+g*partition['q']+coverage.ceil_log2(n+1)+5
    linear_contribution_bits = cap_exponent+g*partition['q']+6
    linear_rejection_cap = calls_bits+failure_bits+coverage.ceil_log2(n)
    arithmetic_checks = {'cap_exponent_at_most_n': 0 < tau+20*eps <= 1,
                         'uniform_quadratic_bound': n >= max(4, 2*(g+6)),
                         'evaluation_error': evaluation_bits >= 5,
                         'rounding_error': rounding_bits >= 5}
    sampling = coverage.ideal_checksum_conditions(candidate, power, target_degree, failure_bits)
    event_checks = coverage.uniform_conditions(c, power, 'minimum-block', failure_bits)
    checks = {'coverage_caps': all(event_checks.values()), 'sampling': sampling['status'] == 'PASS',
              'coupling': all(coupling_checks.values()), 'arithmetic': all(arithmetic_checks.values())}
    return {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
            'power': power, 'n': n, 'target_degree': target_degree, 'partition': partition,
            'ordinary_uniform_call_count_log2_at_most': calls_bits,
            'uniform_rejections_per_call': n**3,
            'coupling_failure_bits_at_least_at_cutoff': n**3-calls_bits,
            'proposed_uniform_rejection_cap': linear_rejection_cap,
            'proposed_coupling_failure_bits_at_least': linear_rejection_cap-calls_bits,
            'coupling_uniform_conditions': coupling_checks,
            'arithmetic_uniform_conditions': arithmetic_checks,
            'evaluation_error_bits_at_least_at_cutoff': evaluation_bits,
            'rounding_error_bits_at_least_at_cutoff': rounding_bits,
            'proposed_unit_phase_accuracy_bits': linear_phase_bits,
            'proposed_contribution_fractional_bits': linear_contribution_bits,
            'proposed_evaluation_error_bits_at_least': (linear_phase_bits-2*n-3*cap_exponent
                                                      -g*partition['q']-coverage.ceil_log2(n+1)),
            'proposed_rounding_error_bits_at_least': linear_contribution_bits-cap_exponent-g*partition['q']-1,
            'per_target_arithmetic_error_at_most': '1/16',
            'per_target_total_error_at_most': '1/8',
            'each_random_failure_at_most': f'2^-{failure_bits}',
            'three_failure_terms': f'3*2^-{failure_bits}',
            'combined_failure_at_most': f'2^-{failure_bits-2}'}


def report(author):
    coverage.require(author == 'conorholds', 'approved author handle required')
    candidate = json.loads((coverage.PREDECESSOR/'FINITE-INPUT.json').read_text())['candidate']
    c = coverage.constants(candidate)
    probability = coverage.uniform_cutoff(c, 'minimum-block', 130)
    work = overhead_cutoff(c)
    power = max(probability['power'], work['power'])
    finite = finite_conditions(candidate, power)
    coverage.require(finite['status'] == 'PASS', 'finite checksum component conditions failed')
    return coverage.serializable({'author': author, 'status': 'PASS', 'candidate_id': candidate['id'],
                                  'scope': 'Checksum-only error and exponential overhead bounds, conditional on the retained filtering and deterministic evaluator proofs; no full operation-count or Subset Sum cutoff.',
                                  'probability_cutoff': finite_conditions(candidate, probability['power']),
                                  'table_and_pattern_overhead_cutoff': work,
                                  'combined_cutoff': finite,
                                  'inputs_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                                    for name in ['implementation_bounds.py', 'coverage_bounds.py', 'INPUTS.json',
                                                                 '../checksum-variance/FINITE-INPUT.json',
                                                                 '../checksum-variance/finite_circle.py']}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = report(args.author)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'],
                      'probability_power': result['probability_cutoff']['power'],
                      'probability_and_overhead_power': result['combined_cutoff']['power']}))


if __name__ == '__main__':
    main()
