"""Exact finite controls for the preliminary-search/representation proof interface."""
import argparse
from collections import Counter
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORE = HERE.parent/'filter-budgets'
SOURCE_COMMIT = '1581d4cbdd69a0faf34b2089614fe2a814316a6a'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(value):
    require(type(value) in (str, int), 'exact strings or integers required')
    return Fraction(value)


def ceil(value):
    return -(-value.numerator//value.denominator)


def outward(value, upper=False):
    with localcontext() as context:
        context.prec = 12
        context.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        decimal = str(Decimal(value.numerator)/Decimal(value.denominator))
    return {'exact': str(value), 'outward_decimal': decimal}


def coupling_variance_bound(n, q):
    """The smaller of two complementary martingale sums, for small controls."""
    require(type(n) is int and 1 <= n <= 64 and type(q) is int and 0 <= q <= n,
            'coupling controls require 1 <= n <= 64 and 0 <= q <= n')
    if q in (0, n):
        return Fraction(0)
    return min(Fraction((n-size)**2, 4)*sum((Fraction(1, k*k) for k in range(n-size, n)), Fraction(0))
               for size in (q, n-q))


def variance_bound(n, q):
    """Fast closed rational upper bound for any exchange-1-Lipschitz function."""
    require(type(n) is int and n >= 1 and type(q) is int and 0 <= q <= n,
            'positive integer n and a size from 0 through n required')
    if q in (0, n):
        return Fraction(0)
    return min(Fraction(1, 4)+Fraction((n-size)*(size-1), 4*(n-1))
               for size in (q, n-q))


def subset_sums(weights):
    require(len(weights) <= 16 and all(type(w) is int for w in weights),
            'at most 16 integer weights required by the finite oracle')
    support = {0}
    for weight in weights:
        support |= {value+weight for value in support}
    return support


def support_controls(max_n):
    require(type(max_n) is int and 0 <= max_n <= 6, 'support enumeration is capped at n=6')
    counts = Counter()
    for n in range(max_n+1):
        for weights in combinations_with_replacement(range(-2, 3), n):
            counts['weight_multisets'] += 1
            support = [{0}]
            for mask in range(1, 1 << n):
                bit = mask & -mask
                previous = support[mask ^ bit]
                support.append(previous | {value+weights[bit.bit_length()-1] for value in previous})
            whole = (1 << n)-1
            for mask, values in enumerate(support):
                counts['supports'] += 1
                require(1 <= len(values) <= 1 << mask.bit_count(), 'invalid support cardinality')
                for index in range(n):
                    bit = 1 << index
                    if not mask & bit:
                        require(len(values) <= len(support[mask | bit]) <= 2*len(values), 'deficiency increment failed')
                        counts['increments'] += 1
                    else:
                        for other in range(n):
                            exchange = 1 << other
                            if not mask & exchange:
                                size = len(support[(mask ^ bit) | exchange])
                                require(size <= 2*len(values) and len(values) <= 2*size, 'exchange Lipschitz bound failed')
                                counts['exchanges'] += 1
                other = whole ^ mask
                while True:
                    require(len(support[mask | other]) <= len(values)*len(support[other]), 'superadditivity failed')
                    counts['disjoint_unions'] += 1
                    if other == 0:
                        break
                    other = (other-1) & (whole ^ mask)
                flipped = tuple(-weight if mask & (1 << index) else weight for index, weight in enumerate(weights))
                translate = sum(weight for index, weight in enumerate(weights) if mask & (1 << index))
                require(subset_sums(flipped) == {value-translate for value in support[whole]}, 'signed support translation failed')
                counts['sign_translations'] += 1
    return dict(counts)


def sampling_controls(max_n):
    require(type(max_n) is int and 2 <= max_n <= 6, 'sampling enumeration requires 2 <= n <= 6')
    cases = 0
    for n in range(2, max_n+1):
        for m in range(2, n+1, 2):
            for alias in range(1 << n):
                counts = Counter()
                for positions in combinations(range(n), m):
                    mask = sum(1 << index for index in positions)
                    for flip in range(1 << n):
                        selected = (alias ^ flip) & mask
                        if selected.bit_count() == m//2:
                            counts[selected] += 1
                require(len(counts) == comb(n, m//2) and len(set(counts.values())) == 1,
                        'balanced selected block is not uniform')
                require(sum(counts.values()) == comb(n, m)*comb(m, m//2)*(1 << (n-m)),
                        'balance-event probability failed')
                cases += 1
    return cases


def variance_controls(max_n=10):
    require(type(max_n) is int and 2 <= max_n <= 10, 'variance controls require 2 <= n <= 10')
    cases = 0
    for n in range(2, max_n+1):
        for q in range(n+1):
            require(coupling_variance_bound(n, q) <= variance_bound(n, q), 'closed variance majorant failed')
            require(variance_bound(n, q) <= Fraction(n+3, 16), 'uniform variance bound failed')
            for distinguished in range(n+1):
                values = [sum(index < distinguished for index in selected) for selected in combinations(range(n), q)]
                mean = Fraction(sum(values), len(values))
                variance = sum(((value-mean)**2 for value in values), Fraction(0))/len(values)
                require(variance <= coupling_variance_bound(n, q), 'exchange-function variance bound failed')
                cases += 1
    return cases


def continuation_counterexample():
    n, half, large = 16, 8, 6
    cap = 1 << ceil(Fraction(2, 5)*n)
    n0 = 256
    require(12*n <= n0 < 24*n, 'illegal isolation range')
    scale = 1+n*n0
    original = [1]*10+[(n+1)*(1 << index) for index in range(large)]
    transformed = [scale*weight+1 for weight in original]
    distribution = []
    failures = 0
    for count in range(large+1):
        multiplicity = comb(large, count)*comb(n-large, half-count)
        size = (1 << count)*(half-count+1)
        other_size = (1 << (large-count))*(half-(large-count)+1)
        example = transformed[:half-count]+transformed[10:10+count]
        require(len(example) == half and len(subset_sums(example)) == size, 'counterexample support formula failed')
        if max(size, other_size) > cap:
            failures += multiplicity
        distribution.append({'large_positions': count, 'halves': multiplicity, 'support_size': size})
    require(6**5 < 2**13, 'elementary logarithm bound failed')
    probability = Fraction(failures, comb(n, half))
    require(probability == Fraction(1, 143), 'continuation probability changed')
    return {'n': n, 'original_weights': original, 'tie_weights': [1]*n, 'isolation_N0': n0,
            'transformed_weights': transformed, 'preliminary_density': '2/5', 'deficiency_density': '3/20',
            'cap': cap, 'half_support_distribution': distribution,
            'expected_deficiency_argument': 'E d(V) >= 5-log2(6) > 12/5 = (3/20)*16, since 6**5 < 2**13.',
            'uniform_half_cap_continue_probability': str(probability)}


def scalable_counterexample(parameters):
    """A positive-weight family that also covers all retained p,d thresholds."""
    n = 128
    half = n//2
    log_upper = (n//4).bit_length()  # log2(n/4+1) <= this integer.
    mean_strict_lower = n//4-log_upper
    cap_exponent = ceil(exact(parameters['preliminary'])*n)
    require(mean_strict_lower > exact(parameters['deficiency'])*n and cap_exponent < half,
            'retained thresholds are outside the scalable counterexample')
    return {'n': n, 'expected_deficiency_strict_lower': mean_strict_lower,
            'cap_exponent': cap_exponent, 'one_half_support_exponent': half,
            'uniform_half_cap_continue_probability_lower': str(Fraction(2, comb(n, half)))}


def gaps(parameters):
    p = {key: exact(parameters[key]) for key in ['preliminary', 'deficiency', 'rho', 'local_deficiency', 'gamma']}
    j = parameters['blocks']
    require(type(j) is int and j > 0, 'positive integer block count required')
    require(0 < p['preliminary'] < Fraction(1, 2) and 0 < p['deficiency'] < Fraction(1, 2)
            and 0 < p['rho'] < Fraction(1, 2) and p['local_deficiency'] > 0 and p['gamma'] > 0,
            'invalid preliminary or representation parameter range')
    result = {'preliminary': p['preliminary']-Fraction(1, 2)+p['deficiency'],
              'block_space': 1-j*p['rho'],
              'coverage': j*p['rho']*p['local_deficiency']-p['deficiency'],
              'representations': Fraction(1, 2)-p['local_deficiency']-p['gamma']}
    require(all(value > 0 for value in result.values()), 'a required interface gap is not positive')
    return result


def useful_probability(parameters, n):
    gaps(parameters)
    require(type(n) is int and n >= 1, 'positive integer n required')
    rho, delta, d = (exact(parameters[key]) for key in ['rho', 'local_deficiency', 'deficiency'])
    j = parameters['blocks']
    m = 2*((rho*n)//2)
    require(j*(m//2) <= n//2, 'blocks do not fit in a half')
    if m == 0:
        return Fraction(0)
    mean_upper = d*n/j
    markov = max(Fraction(0), 1-mean_upper/(delta*m))
    separation = delta*m-mean_upper
    cantelli = (separation**2/(variance_bound(n, m//2)+separation**2)
                if separation > 0 else Fraction(0))
    return max(markov, cantelli)


def preliminary_continue_upper(parameters, n):
    gaps(parameters)
    require(type(n) is int and n >= 1, 'positive integer n required')
    cap_exponent = ceil(exact(parameters['preliminary'])*n)
    bound = Fraction(1, n)  # n*2**(-n*n) <= 1/n; avoids constructing an enormous denominator.
    for size in (n//2, n-n//2):
        separation = exact(parameters['deficiency'])*n-size+cap_exponent
        if separation <= 0:
            return Fraction(1)
        variance = variance_bound(n, size)
        bound += variance/(variance+separation**2)
    return min(Fraction(1), bound)


def parameter_report(parameters):
    margins = gaps(parameters)
    rho, delta = (exact(parameters[key]) for key in ['rho', 'local_deficiency'])
    j = parameters['blocks']
    floor_cutoff = ceil(4*j*delta/margins['coverage'])
    half_cutoff = max(floor_cutoff, ceil(rho*j*j/(2*margins['coverage']**2)))
    require(useful_probability(parameters, half_cutoff) >= Fraction(1, 2), 'half-probability cutoff failed')
    block_size = 2
    for _ in range(64):
        if margins['representations']*block_size >= block_size.bit_length():
            break
        block_size *= 2
    else:
        raise ValueError('cardinality-loss cutoff exceeds the bounded dyadic search')
    return {'gaps': {key: str(value) for key, value in margins.items()},
            'limiting_markov_lower': outward(margins['coverage']/(j*rho*delta)),
            'n_useful_probability_at_least_half': half_cutoff,
            'dyadic_m_absorbing_cardinality_loss': block_size,
            'n_absorbing_cardinality_loss': ceil(Fraction(block_size)/rho),
            'scalable_counterexample': scalable_counterexample(parameters),
            'finite_bounds': [{'n': n, 'useful_block_lower': outward(useful_probability(parameters, n)),
                               'high_deficiency_continue_upper': outward(preliminary_continue_upper(parameters, n), upper=True)}
                              for n in [10**6, 10**9, 10**12]]}


def run(max_n=6):
    require(type(max_n) is int and 2 <= max_n <= 6, 'campaign size must be from 2 through 6')
    source = json.loads((CORE/'SOURCE.json').read_text())
    require(source['commit'] == SOURCE_COMMIT, 'unexpected mathematical source pin')
    inputs = json.loads((HERE/'parameters.json').read_text())
    rows = [{'id': witness['id'], 'commit': entry['commit'], **parameter_report(witness['parameters'])}
            for entry in inputs['sources'] for witness in entry['parameters']]
    return {'status': 'PASS',
            'scope': 'Exact finite controls and rational consequences of a written proof audit; not a full algorithm proof or a numerical filter certificate.',
            'source_commit': SOURCE_COMMIT,
            'sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['audit.py', 'parameters.json']},
            'source_manifest_sha256': hashlib.sha256((CORE/'SOURCE.json').read_bytes()).hexdigest(),
            'limits': {'support_max_n': max_n, 'weights': [-2, -1, 0, 1, 2], 'sampling_max_n': max_n,
                       'variance_max_n': 10, 'cardinality_doublings': 64},
            'support_controls': support_controls(max_n), 'sampling_cases': sampling_controls(max_n),
            'variance_cases': variance_controls(), 'continuation_counterexample': continuation_counterexample(),
            'parameters': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-n', type=int, default=6)
    args = parser.parse_args()
    result = run(args.max_n)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'support_controls': result['support_controls'],
                      'sampling_cases': result['sampling_cases'], 'variance_cases': result['variance_cases'],
                      'parameters': len(result['parameters']),
                      'counterexample_cap_continue_probability': result['continuation_counterexample']['uniform_half_cap_continue_probability']}))


if __name__ == '__main__':
    main()
