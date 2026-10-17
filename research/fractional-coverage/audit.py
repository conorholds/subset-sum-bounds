"""Exact finite controls for fractional subset-sum coverage and its cutoffs."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
from math import comb, lcm, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(value):
    require(type(value) in (str, int), 'exact strings or integers required')
    return Fraction(value)


def ceil(value):
    return -(-value.numerator//value.denominator)


def supports(weights):
    require(len(weights) <= 8 and all(type(value) is int for value in weights),
            'finite oracle accepts at most eight integer weights')
    data = [{0}]
    for mask in range(1, 1 << len(weights)):
        bit = mask & -mask
        before = data[mask ^ bit]
        data.append(before | {value+weights[bit.bit_length()-1] for value in before})
    return data


def canonical_masks(weights):
    """One least binary mask for each sum, under a restriction-consistent order."""
    supports(weights)  # Validate the deliberately bounded finite domain.
    sums, first = [0], {0: 0}
    for mask in range(1, 1 << len(weights)):
        bit = mask & -mask
        value = sums[mask ^ bit]+weights[bit.bit_length()-1]
        sums.append(value)
        first.setdefault(value, mask)
    return first, sums


def cardinality_controls(max_n=7):
    require(type(max_n) is int and 0 <= max_n <= 8, 'cardinality campaign is capped at n=8')
    counts = Counter()
    for n in range(max_n+1):
        for weights in combinations_with_replacement(range(-2, 3), n):
            counts['weight_multisets'] += 1
            sets = supports(weights)
            canonical, sums = canonical_masks(weights)
            products = [1]*(n+1)
            for mask, support in enumerate(sets):
                size = mask.bit_count()
                products[size] *= len(support)
                projection = {representative & mask for representative in canonical.values()}
                require(len({sums[value] for value in projection}) == len(projection),
                        'canonical projection is not sum-injective')
                require(len(projection) <= len(support), 'projection exceeds its sumset')
                counts['canonical_projections'] += 1
            for q in range(1, n):
                # Exact form of E log|Sigma(V_q)|/q >= E log|Sigma(V_{q+1})|/(q+1).
                require(products[q]**(n-q) >= products[q+1]**q,
                        'normalized restriction inequality failed')
                counts['normalized_adjacent_inequalities'] += 1
            for q in range(1, n+1):
                exponent = comb(n-1, q-1)
                require(len(sets[-1])**exponent <= products[q], 'uniform-cover product failed')
                counts['uniform_cover_inequalities'] += 1
    return dict(counts)


def superadditivity_counterexample():
    # f(U)=floor(|U|/3) has all the source's elementary deficiency properties,
    # yet E f(V_3)=1 > (3/5) E f(V_5)=3/5.
    for a in range(9):
        require(0 <= (a+1)//3-a//3 <= 1, 'counterexample increment failed')
        for b in range(9-a):
            require((a+b)//3 >= a//3+b//3, 'counterexample superadditivity failed')
    require(Fraction(1) > Fraction(3, 5), 'counterexample was not strict')
    return {'function': 'floor(|U|/3)', 'h': 5, 'q': 3,
            'small_expectation': '1', 'claimed_fractional_upper': '3/5',
            'meaning': 'The new sumset lemma does not follow from superadditivity and exchange bounds alone.'}


def sharpness_controls():
    cases = 0
    for n in range(1, 13):
        for zeros in range(n+1):
            for q in range(n+1):
                mean = sum((Fraction(k*comb(zeros, k)*comb(n-zeros, q-k), comb(n, q))
                            for k in range(max(0, q-(n-zeros)), min(q, zeros)+1)), Fraction(0))
                require(mean == Fraction(q*zeros, n), 'sharpness mean failed')
                cases += 1
    # Check the exact support formula independently for all submasks of small examples.
    support_cases = 0
    for n in range(1, 9):
        for zeros in range(n+1):
            weights = [0]*zeros+[1 << j for j in range(n-zeros)]
            for mask, values in enumerate(supports(weights)):
                nonzero_count = (mask >> zeros).bit_count()
                require(len(values) == 1 << nonzero_count, 'zero/powers support formula failed')
                support_cases += 1
    return {'hypergeometric_cases': cases, 'exact_support_cases': support_cases,
            'family': 'zero weights and distinct powers of two',
            'expected_deficiency': 'q*zero_count/n', 'normalized_inequality': 'equality'}


def positive_family(n, small):
    require(type(n) is int and 1 <= n <= 64 and type(small) is int and 0 <= small <= n,
            'positive finite family requires 1 <= n <= 64 and a valid small-weight count')
    n0 = 1
    while n0 < 12*n:
        n0 *= 2
    require(n0 < 24*n, 'isolation range failed')
    original = [1]*small+[(n+1)*(1 << j) for j in range(n-small)]
    c = 1+n*n0
    return n0, original, [c*value+1 for value in original]


def positive_support_controls():
    count = 0
    for n in range(1, 9):
        for small in range(n+1):
            _, _, weights = positive_family(n, small)
            for mask, values in enumerate(supports(weights)):
                large = (mask >> small).bit_count()
                little = (mask & ((1 << small)-1)).bit_count()
                require(len(values) == (1 << large)*(little+1), 'positive support formula failed')
                count += 1
    return count


def positive_toy():
    n, half, small, large, m, q = 14, 7, 7, 7, 4, 2
    d, delta, rho = Fraction(3, 25), Fraction(13, 100), Fraction(2, 7)
    n0, original, transformed = positive_family(n, small)
    half_product = prod(((1 << k)*(half-k+1))**(comb(large, k)*comb(small, half-k))
                        for k in range(half+1))
    # Exact comparison E D(V_half)<=d*n, with no approximate logarithms.
    require(half_product**d.denominator >=
            1 << ((half*d.denominator-d.numerator*n)*comb(n, half)),
            'positive toy does not meet the low-deficiency hypothesis')
    for k in range(q+1):
        size = (1 << k)*(q-k+1)
        require(size**delta.denominator >= 1 << (q*delta.denominator-delta.numerator*m),
                'positive toy restriction is not useful')
    j = ceil(1/rho)-1
    require(j*rho*delta-d < 0 and delta > d, 'positive toy does not distinguish the models')
    return {'n': n, 'm': m, 'rho': str(rho), 'd': str(d), 'delta': str(delta),
            'N0': n0, 'tie_weights': 'all one', 'original_weights': original,
            'transformed_weights': transformed, 'uniform_useful_restriction_probability': '1',
            'largest_strictly_fitting_J': j, 'old_coverage_margin': str(j*rho*delta-d),
            'scope': 'A finite coverage-interface control, not a complete algorithm parameter tuple.'}


def below_threshold_obstruction(parameters):
    rho, d = [exact(parameters[key]) for key in ['rho', 'deficiency']]
    require(0 < rho < Fraction(1, 2) and 0 < d < Fraction(1, 2), 'valid density and deficiency required')
    delta = d/2
    gap = d-delta
    # For q>=2^k, log2(q+1)<=gap*q; the same inequality persists by doubling.
    k = 0
    while gap*(1 << k) < k+1:
        k += 1
        require(k <= 128, 'obstruction cutoff exceeds bounded dyadic search')
    q0 = 1 << k
    guard = exact(parameters['guard'])
    require(guard > 0, 'positive bit guard required')
    guard_k = 1
    while guard*(1 << guard_k) < guard_k+2:
        guard_k += 1
        require(guard_k <= 128, 'bit-guard cutoff exceeds bounded dyadic search')
    # z=2*d*n and q=rho*n/2 must be integers; m=rho*n is then exactly even.
    multiple = lcm((2*d).denominator, 2*rho.denominator, 2)
    cutoff = max(ceil(2*q0/rho), 1 << guard_k)
    first = multiple*ceil(Fraction(cutoff, multiple))
    exponent = rho*gap**2
    require(exponent > 0 and 2*d*first == int(2*d*first)
            and rho*first/2 == int(rho*first/2), 'obstruction compatibility failed')
    return {'rho': str(rho), 'd': str(d), 'smaller_delta': str(delta),
            'compatible_n_multiple': multiple, 'first_compatible_n_after_cutoffs': first,
            'q_logarithm_cutoff': q0, 'bit_guard_n_cutoff': 1 << guard_k,
            'base_two_probability_exponent_per_n': str(exponent),
            'exponent_at_first_compatible_n': str(exponent*first),
            'bound': 'Pr[D(S)<=delta*m | balance] <= 2^(-exponent_per_n*n)',
            'scope': 'Positive legal transformed family; obstruction to a uniform useful-block guarantee below d, not to the full algorithm.'}


def restriction_profile(parameters, n):
    require(type(n) is int and n >= 2, 'integer n >= 2 required')
    rho, d, delta = [exact(parameters[key]) for key in ['rho', 'deficiency', 'local_deficiency']]
    require(0 < rho < Fraction(1, 2) and 0 <= d < delta, 'positive fractional coverage gap required')
    half = n//2
    m = 2*((rho*n)//2)
    require(m//2 <= half, 'restriction does not fit inside a half')
    mean_upper = Fraction(m, 2*half)*d*n
    separation = delta*m-mean_upper
    return m, mean_upper, separation


def useful_probability(parameters, n):
    m, mean_upper, separation = restriction_profile(parameters, n)
    if not m:
        return Fraction(0)
    delta = exact(parameters['local_deficiency'])
    markov = max(Fraction(0), 1-mean_upper/(delta*m))
    # Conditional martingale values lie in an interval of length at most one;
    # Popoviciu's bound gives Var D(V_q) <= q/4, with q=m/2.
    cantelli = (separation**2/(Fraction(m, 8)+separation**2)
                if separation > 0 else Fraction(0))
    return max(markov, cantelli)


def exponential_failure_exponent(parameters, n):
    """A lower exponent e with useful-block failure <= 2^(-e)."""
    m, _, separation = restriction_profile(parameters, n)
    return 4*separation**2/m if m and separation > 0 else Fraction(0)


def half_probability_cutoff(parameters):
    rho, d, delta = [exact(parameters[key]) for key in ['rho', 'deficiency', 'local_deficiency']]
    require(0 < rho < Fraction(1, 2) and 0 <= d < delta, 'positive fractional coverage gap required')
    gap = delta-d
    return max(2, ceil(4/rho+4*d/gap), ceil(1/(2*rho*gap**2)))


def probability_controls(max_n=32):
    """Compare the bounds with exact hypergeometric probabilities, including floors."""
    require(type(max_n) is int and 5 <= max_n <= 32, 'probability campaign requires 5 <= n <= 32')
    cases = 0
    for n in range(5, max_n+1):
        for m in range(2, (n-1)//2+1, 2):
            q = m//2
            for zeros in range(n+1):
                # D(S) equals its zero count for zeros and distinct powers.
                # This d gives equality in the finite low-deficiency hypothesis.
                d = Fraction((n//2)*zeros, n*n)
                for extra in [Fraction(1, 4*n), Fraction(1, n), Fraction(1, 20), Fraction(1, 4)]:
                    delta = d+extra
                    actual = sum((Fraction(comb(zeros, k)*comb(n-zeros, q-k), comb(n, q))
                                  for k in range(max(0, q-(n-zeros)), min(q, zeros)+1)
                                  if k <= delta*m), Fraction(0))
                    parameters = {'rho': str(Fraction(m, n)), 'deficiency': str(d),
                                  'local_deficiency': str(delta)}
                    require(useful_probability(parameters, n) <= actual,
                            'component lower bound exceeds exact useful probability')
                    exponent = exponential_failure_exponent(parameters, n)
                    bits = exponent.numerator//exponent.denominator
                    require(1-Fraction(1, 1 << bits) <= actual,
                            'exponential concentration exceeds exact useful probability')
                    cases += 1
    return cases


def parameter_report(candidate):
    p = candidate['parameters']
    rho, d, delta = [exact(p[key]) for key in ['rho', 'deficiency', 'local_deficiency']]
    jmax = ceil(1/rho)-1
    cutoff = half_probability_cutoff(p)
    require(useful_probability(p, cutoff) >= Fraction(1, 2), 'half-probability cutoff failed')
    require(useful_probability(p, cutoff+1) >= Fraction(1, 2), 'odd/even cutoff control failed')
    return {'id': candidate['id'], 'fractional_coverage_gap_per_n': str(rho*(delta-d)),
            'largest_strictly_fitting_integer_J': jmax,
            'integer_packing_coverage_margin_at_Jmax': str(jmax*rho*delta-d),
            'n_useful_probability_at_least_half': cutoff,
            'useful_probability_at_sizes': {str(n): str(useful_probability(p, n))
                                           for n in [10**8, 10**10, 10**12]},
            'useful_failure_dyadic_upper_at_sizes': {
                str(n): '2^(-'+str((exponential_failure_exponent(p, n)//1))+')'
                for n in [10**8, 10**10, 10**12]}}


def input_hashes():
    names = ['audit.py', 'budget.py', 'candidates.json', 'predecessor.json']
    return {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names}


def run(max_n=7, author='conorholds'):
    require(isinstance(author, str) and author, 'explicit nonempty author required')
    return {'status': 'PASS', 'author': author,
            'scope': 'Exact finite controls and component cutoffs; the written all-size proof needs review.',
            'inputs_sha256': input_hashes(), 'limits': {'maximum_n': max_n, 'probability_maximum_n': 32,
                                                      'integer_weights': [-2, -1, 0, 1, 2]},
            'cardinality_controls': cardinality_controls(max_n),
            'superadditivity_counterexample': superadditivity_counterexample(),
            'sharpness': sharpness_controls(),
            'positive_support_controls': positive_support_controls(),
            'positive_toy': positive_toy(),
            'exact_probability_controls': probability_controls(),
            'below_threshold_obstruction': below_threshold_obstruction(
                json.loads((HERE/'candidates.json').read_text())['candidates'][1]['parameters']),
            'parameters': [parameter_report(row) for row in json.loads((HERE/'candidates.json').read_text())['candidates']]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-n', type=int, default=7)
    args = parser.parse_args()
    report = run(args.max_n, args.author)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': report['status'], 'cardinality_controls': report['cardinality_controls'],
                      'sharpness': report['sharpness']}))
