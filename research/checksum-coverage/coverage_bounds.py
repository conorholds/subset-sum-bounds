"""Exact sufficient thresholds for checksum coverage and simultaneous caps.

This calculates bounds from a written argument, not actual random tables or
the continuous filtering certificate. Huge powers in probability estimates
are compared through their exponents and are never materialized.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREDECESSOR = HERE.parent/'checksum-variance'
_spec = importlib.util.spec_from_file_location('coverage_finite_input', PREDECESSOR/'finite_circle.py')
_finite = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_finite)
ceil, exact, require = _finite.ceil, _finite.exact, _finite.require


def positive_integer(value, label):
    require(type(value) is int and 1 <= value and value.bit_length() <= 4096,
            label+' must be a positive integer of at most 4096 bits')


def ceil_log2(value):
    positive_integer(value, 'logarithm argument')
    return (value-1).bit_length()


def constants(candidate):
    # Reuse the predecessor's exact input/domain validation and coefficient
    # bounds. Its moment cutoff is a separate dependency, not re-proved here.
    _finite.cutoff(candidate)
    design = candidate['design']
    r, _ = _finite.polynomial_bounds(design['r'], design['r_sin'])
    s, _ = _finite.polynomial_bounds(design['s'], design['s_sin'])
    kappa = exact(design['kappa'])
    epsilon = exact(candidate['parameters']['epsilon'])
    require(Fraction(1, 1 << 64) <= epsilon <= Fraction(1, 8),
            'coverage calculator requires 2^-64 <= epsilon <= 1/8')
    widths = [2*(r+s), 2*kappa*r, 2*kappa*s]
    # Allocate failure <1/4 to the upper likelihood tail and <1/8 to
    # each lower score tail. exp(7/10)>2 makes these rational constants
    # sufficient in the one-sided bounded-variable inequality.
    block_bounds = [ceil(Fraction(7, 10)*widths[0]**2/epsilon**2)]
    block_bounds += [ceil(Fraction(21, 20)*w**2/epsilon**2) for w in widths[1:]]
    d0 = max(1, *block_bounds)
    step = epsilon/4
    # Outward endpoint rounding: error <=epsilon/2. The upper counts
    # include boundary indices, even when a boundary is itself a grid point.
    grid_counts = [ceil((2*(r+s)+epsilon/2)/step)+1]
    grid_counts += [ceil(w/step)+ceil((w+epsilon/2)/step)+1 for w in widths[1:]]
    grid_bits = sum(ceil_log2(count) for count in grid_counts)
    return {'epsilon': epsilon, 'eta': exact(candidate['parameters']['eta']),
            'likelihood_width': widths[0], 'high_width': widths[1], 'low_width': widths[2],
            'd_bound': 2*(r+s)+epsilon/2, 'concentration_blocks': block_bounds,
            'minimum_block_size': d0, 'grid_step': step, 'grid_counts': grid_counts,
            'pattern_bits_per_block': grid_bits}


def round_profile(approximations, epsilon):
    """Round values known to be within epsilon/8 of the true averages.

    Correctness of that numerical approximation is a caller obligation.
    This routine does not evaluate the filter or certify those errors.
    """
    require(isinstance(approximations, tuple) and len(approximations) == 3
            and all(isinstance(x, Fraction) for x in approximations),
            'three exact rational approximations required')
    require(isinstance(epsilon, Fraction) and epsilon > 0, 'positive rational epsilon required')
    step = epsilon/4
    d, a, e = approximations
    result = (step*ceil((d+epsilon/8)/step),
              step*((a-epsilon/8)//step), step*((e-epsilon/8)//step))
    require(result[0] >= 0 and result[1] <= result[0] and result[2] <= result[0],
            'approximations violate the promised entropy-order interface')
    return result


def balanced(n, q):
    positive_integer(n, 'input size')
    positive_integer(q, 'block count')
    require(q <= n, 'block count exceeds input size')
    minimum, larger = divmod(n, q)
    return {'q': q, 'minimum': minimum, 'maximum': minimum+bool(larger),
            'larger_blocks': larger}


def square_root_partition(n):
    positive_integer(n, 'input size')
    root = isqrt(n)
    return balanced(n, root+(root*root != n))


def required_block_size(c, n):
    positive_integer(n, 'input size')
    epsilon = c['epsilon']
    return max(c['minimum_block_size'], ceil(Fraction(ceil_log2(n)+7)/(2*epsilon)),
               ceil(Fraction(c['pattern_bits_per_block']+1)/(2*epsilon)))


def minimum_block_partition(c, n):
    """Use the source's block count whenever its blocks are large enough."""
    target = required_block_size(c, n)
    require(target <= n, 'input too small for the sufficient minimum block size')
    q = min(square_root_partition(n)['q'], n//target)
    return balanced(n, q)


def point_certificate(c, n, partition, failure_bits=128):
    positive_integer(n, 'input size')
    require(type(failure_bits) is int and 1 <= failure_bits <= 512,
            'failure bits must be an integer from one to 512')
    require(isinstance(partition, dict) and set(partition) == set(balanced(n, 1)),
            'balanced partition schema required')
    require(partition == balanced(n, partition['q']), 'incorrect balanced partition')
    q, dmin, dmax = [partition[key] for key in ['q', 'minimum', 'maximum']]
    epsilon, g = c['epsilon'], c['pattern_bits_per_block']
    # qP exp(-2^r) <= 2^(ceil(log2 q)+eta*n+2-2^r).
    # Demand <=2^(-failure_bits-1), without constructing 2^r.
    r = 2*epsilon*dmin-4
    coverage_required = ceil_log2(ceil(c['eta']*n+2+ceil_log2(q)+failure_bits+1))
    # 13 < 16 gives the following lower bound on the cap decay in bits.
    cap_bits = 6*epsilon*n-(g+1)*q-4
    checks = {'concentration': dmin >= c['minimum_block_size'],
              'coverage': r >= coverage_required, 'caps': cap_bits >= failure_bits+1}
    return {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
            'n': n, 'partition': partition, 'coverage_outer_exponent': r,
            'coverage_required_exponent': coverage_required, 'cap_failure_bits_at_least': cap_bits,
            'pattern_overhead_bits': g*q,
            'per_digit_table_scan_bits_at_most': ceil_log2(q)+1+(c['d_bound']+4*epsilon)*dmax}


def uniform_conditions(c, power, kind, failure_bits=128):
    require(type(power) is int and 1 <= power <= 512, 'dyadic power from one to 512 required')
    require(type(failure_bits) is int and 1 <= failure_bits <= 512,
            'failure bits must be an integer from one to 512')
    require(kind in ['square-root', 'minimum-block'], 'unknown partition rule')
    n = 1 << power
    epsilon, g = c['epsilon'], c['pattern_bits_per_block']
    if kind == 'square-root':
        root = isqrt(n)
        return {'concentration': root-2 >= c['minimum_block_size'],
                'cap_density': root >= 2*(g+1)/epsilon,
                'cap_tail': 5*epsilon*n >= failure_bits+5,
                'coverage_growth': 2*epsilon*root >= 4,
                'coverage_initial': 2*epsilon*(root-2) >= power+7,
                'coverage_tail': 2*n-2 >= failure_bits+1}
    return {'concentration': n >= c['minimum_block_size'],
            'partition_logarithm': 2*epsilon*n >= power+8,
            'partition_growth': 2*epsilon*n >= 2,
            'cap_density': 2*epsilon*n >= g+1,
            'cap_tail': 4*epsilon*n >= failure_bits+5,
            'coverage_tail': n >= failure_bits+1}


def uniform_cutoff(c, kind, failure_bits=128):
    for power in range(1, 513):
        checks = uniform_conditions(c, power, kind, failure_bits)
        if all(checks.values()):
            n = 1 << power
            partition = (square_root_partition(n) if kind == 'square-root'
                         else minimum_block_partition(c, n))
            point = point_certificate(c, n, partition, failure_bits)
            require(point['status'] == 'PASS', 'uniform cutoff fails its point check')
            return {'power': power, 'n': n, 'uniform_conditions': checks, 'at_cutoff': point,
                    'previous_power_conditions': uniform_conditions(c, power-1, kind, failure_bits)
                    if power > 1 else None}
    raise ValueError('no sufficient cutoff in the bounded calculator domain')


def serializable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [serializable(item) for item in value]
    return value


def ideal_checksum_conditions(candidate, power, target_degree=1, failure_bits=128):
    """Sufficient all-n conditions for composing the preceding ideal MSE bound.

    Targets <=n^target_degree and tolerance 1/16 give a sampling-failure
    bound 2048*n^target_degree*2^(-s*n/4). This excludes arithmetic,
    bounded-draw coupling, and independent review of the continuous ceiling.
    """
    require(type(power) is int and 1 <= power <= 512, 'dyadic power from one to 512 required')
    require(type(target_degree) is int and 1 <= target_degree <= 1024,
            'integer target degree from one to 1024 required')
    require(type(failure_bits) is int and 1 <= failure_bits <= 512,
            'failure bits must be an integer from one to 512')
    finite = _finite.cutoff(candidate)
    eta = exact(candidate['parameters']['eta'])
    tau = exact(candidate['parameters']['tau'])
    ratio = (eta-tau)/(2*tau-eta)
    lam1, lam2 = map(exact, candidate['multipliers'])
    require(lam1 >= 1 and lam2 >= 1 and 3*lam1 <= 30 and (4*ratio+3)*lam2 <= 30,
            'checksum domination slopes exceed their retained rounding reserves')
    reserve = Fraction(finite['remaining_variance_reserve'])
    decay = reserve/4
    n = 1 << power
    checks = {'finite_circle': n >= finite['sufficient_input_size'],
              'sampling_initial': decay*n >= target_degree*(power+1)+failure_bits+11,
              'sampling_growth': decay*n >= 2*target_degree,
              'positive_draw_exponent': exact(candidate['parameters']['tau'])-reserve/4 > 0}
    return {'checks': checks, 'status': 'PASS' if all(checks.values()) else 'FAIL',
            'target_degree': target_degree, 'tolerance': '1/16',
            'sampling_exponent': exact(candidate['parameters']['tau'])-reserve/4,
            'sampling_decay': decay, 'sampling_failure_at_most': f'2^-{failure_bits}',
            'combined_failure_at_most': f'2^-{failure_bits-1}'}


def report(author):
    require(author == 'conorholds', 'approved author handle required')
    snapshot = json.loads((PREDECESSOR/'FINITE-INPUT.json').read_text())
    c = constants(snapshot['candidate'])
    partitions = {kind: uniform_cutoff(c, kind) for kind in ['square-root', 'minimum-block']}
    ideal = ideal_checksum_conditions(snapshot['candidate'], partitions['minimum-block']['power'])
    require(ideal['status'] == 'PASS', 'ideal checksum composition misses its component cutoff')
    result = {'author': author, 'status': 'PASS', 'candidate_id': snapshot['candidate']['id'],
              'scope': 'Coverage, multiplicity and simultaneous caps, composed with the preceding ideal checksum MSE; excludes finite arithmetic and a practical runtime cutoff.',
              'failure_probability_at_most': '2^-128', 'constants': c,
              'partitions': partitions, 'ideal_checksum': ideal,
              'inputs_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                for name in ['coverage_bounds.py', 'INPUTS.json',
                                             '../checksum-variance/finite_circle.py',
                                             '../checksum-variance/FINITE-INPUT.json']}}
    return serializable(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = report(args.author)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'candidate': result['candidate_id'],
                      'cutoff_powers': {key: value['power'] for key, value in result['partitions'].items()}}))


if __name__ == '__main__':
    main()
