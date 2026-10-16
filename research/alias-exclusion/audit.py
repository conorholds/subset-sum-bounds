"""Bounded exact alias controls and the delayed-exclusion cost obstruction."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CORE = HERE.parent/'filter-budgets'
sys.path.insert(0, str(CORE))
from certify import budget_margins, enclosure, iv, log2, q
from oracle import balanced_left_size, exhaustive_layer, layer, require


def masks(weights, target, large):
    return {mask for mask in range(1 << len(weights))
            if (value := sum(w for i, w in enumerate(weights) if mask & (1 << i))) != target
            and value % large == target % large}


def finite_campaign():
    cases = [([1, 1, 1, 1], 2, [0, 1], range(16), [3, 5, 7]),
             ([1, 2, 4, 8], 5, [0, 2], range(16), [3, 5, 7]),
             ([0, 0, -1, 2], 0, [1, 3], range(16), [3, 5, 7]),
             ([1, 1, 3, 3, 7, 7], 6, [1, 4], [0, 1, 7, 15, 21, 42, 62, 63], [3, 5, 7]),
             ([1, 2, 3, 5, 8, 13, 21, 34], 17, [0, 2, 4, 6], [0, 1, 7, 15, 85, 170, 254, 255], [7])]
    layers = completeness = visits = exclusions = overlapping = aliases = 0
    digest = hashlib.sha256()
    for weights, target, block, flips, primes in cases:
        n, m = len(weights), len(block)
        block_mask = sum(1 << index for index in block)
        for large in primes:
            actual = masks(weights, target, large)
            for flip in flips:
                expected = {mask for mask in actual if ((mask ^ flip) & block_mask).bit_count() == m//2}
                for small in [2, 3]:
                    recovered = set()
                    for k in range(m//2+1):
                        for carry in range(n+1):
                            for residue in range(small):
                                args = (weights, target, large, block, flip, k, carry, small, residue)
                                grouped = layer(*args)
                                direct = exhaustive_layer(*args)
                                require(grouped['status'] == 'PASS', 'unexpected toy cap')
                                for key in ['left_length', 'right_length', 'excluded_exact', 'aliases']:
                                    require(grouped[key] == direct[key], 'grouped/exhaustive mismatch: '+key)
                                require(sorted(grouped['visited_pairs']) == sorted(direct['visited_pairs']), 'visited pairs differ')
                                require(set(grouped['aliases']) <= actual, 'saved a false alias')
                                recovered.update(grouped['aliases'])
                                visits += len(grouped['visited_pairs'])
                                exclusions += grouped['excluded_exact']
                                overlapping += direct['excluded_exact_overlapping']
                                aliases += len(grouped['aliases'])
                                layers += 1
                                digest.update(json.dumps(grouped, sort_keys=True).encode())
                    require(recovered == expected, 'all-residue recovery failed for a balanced alias')
                    completeness += 1
    return {'layers': layers, 'balanced_recovery_controls': completeness, 'visited_pairs': visits,
            'exact_pairs_skipped': exclusions, 'skipped_exact_overlapping_pairs': overlapping,
            'saved_aliases_before_cross_layer_deduplication': aliases, 'grouped_trace_sha256': digest.hexdigest()}


def flip_equality_probability(n, overlap):
    require(type(n) is int and 2 <= n <= 512 and n % 2 == 0 and type(overlap) is int
            and 0 <= overlap <= n, 'require even 2 <= n <= 512 and a valid overlap count')
    favorable = 0
    for count in range(overlap+1):
        rest = (n+2*overlap)//2-3*count
        if 0 <= rest <= n-overlap:
            favorable += comb(overlap, count)*comb(n-overlap, rest)
    return Fraction(favorable, 1 << n)


def homogeneous_counts(n, m):
    require(type(n) is int and type(m) is int and 32 <= m <= n//2 and m % 32 == 0
            and n % 2 == 0 and n <= 512, 'type controls require m divisible by 32, 2m <= even n <= 512')
    overlap = m//16
    sizes = (overlap, 3*m//16, 3*m//16, 9*m//16)
    denominator = 1
    for size in sizes:
        denominator *= factorial(size)
    types = factorial(m)//denominator
    raw = (1 << (n-m))*types
    probability = flip_equality_probability(n, overlap)
    require(types*(m+1)**4*3**(3*m//2) >= 1 << (4*m), 'type-class entropy lower bound failed')
    require(probability >= Fraction(1, (overlap+1)*(n-overlap+1)), 'balanced-flip lower bound failed')
    return {'n': n, 'm': m, 'overlap': overlap, 'raw_overlapping_pairs': raw,
            'exact_integer_equality_probability_over_flips': str(probability),
            'expected_exact_overlapping_pairs_times_small_modulus': str(raw*probability)}


def homogeneous_toy():
    n, n0 = 8, 128
    require(12*n <= n0 < 24*n and n0 & (n0-1) == 0, 'invalid isolation range')
    weight = 2+n*n0
    skipped = overlapping = visits = 0
    for flip in range(256):
        for residue in range(3):
            args = ([weight]*8, 4*weight, 101, [0, 1, 2, 3], flip, 1, 0, 3, residue)
            grouped, direct = layer(*args), exhaustive_layer(*args)
            require(grouped['status'] == 'PASS' and not grouped['aliases'], 'toy equal weights produced aliases')
            visits += len(grouped['visited_pairs'])
            skipped += grouped['excluded_exact']
            overlapping += direct['excluded_exact_overlapping']
    expected = 64*flip_equality_probability(8, 1)/3
    require(Fraction(overlapping, 256*3) == expected == Fraction(7, 2), 'homogeneous toy expectation failed')
    require(visits == 0, 'exclusion should skip every equal-weight reduced match')
    return {'n': 8, 'm': 4, 'original_weights': [1]*8, 'original_target': 4,
            'tie_weights': [1]*8, 'isolation_N0': n0, 'transformed_weight': weight,
            'large_prime': 101, 'small_prime': 3, 'all_flips_and_residues': 768,
            'exact_pairs_skipped': skipped, 'exact_overlapping_pairs_skipped': overlapping,
            'visited_after_exclusion': visits,
            'expected_overlapping_visits_if_exclusion_delayed': str(expected)}


def cost_profiles():
    iv.dps = 50
    h = -q('.25')*log2(q('.25'))-q('.75')*log2(q('.75'))
    data = json.loads((HERE/'parameters.json').read_text())
    profiles = []
    for source in data['sources']:
        for witness in source['parameters']:
            parameters = witness['parameters']
            budget_margins(parameters)
            rho, gamma, eta, guard, tau = (q(parameters[key]) for key in ['rho', 'gamma', 'eta', 'guard', 'tau'])
            generation = rho+q('.25')
            lists = q('.5')+(h-q('.5')-gamma)*rho
            pairs = 1+(2*h-1-gamma)*rho-eta+guard
            delayed = 1+(2*h-1-gamma)*rho
            require(generation.b < tau.a and lists.b < tau.a and pairs.b < tau.a,
                    'a retained alias cost does not fit its main cap')
            require(delayed.a > 1, 'delayed-exclusion obstruction should exceed exponent one')
            profiles.append({'id': witness['id'], 'commit': source['commit'], 'tau': parameters['tau'],
                             'generation_exponent': enclosure(generation), 'expected_list_exponent': enclosure(lists),
                             'expected_nonexact_pair_exponent_with_guard': enclosure(pairs),
                             'delayed_exclusion_expected_overlap_exponent': enclosure(delayed)})
    return profiles


def run(author):
    require(author == 'conorholds', 'author must be conorholds')
    source = json.loads((CORE/'SOURCE.json').read_text())
    require(source['commit'] == '1581d4cbdd69a0faf34b2089614fe2a814316a6a', 'unexpected manuscript pin')
    balanced = 0
    for n in range(4, 11):
        for m in range(2, n//2+1, 2):
            for k in range(m//2+1):
                size = balanced_left_size(n, m, k)
                left, right = (1 << size)*comb(m, k), (1 << (n-m-size))*comb(m, m//2-k)
                require(max(left, right) <= 2*min(left, right), 'raw records not balanced')
                require(left*right == (1 << (n-m))*comb(m, k)*comb(m, m//2-k), 'raw product depends on balancing')
                balanced += 1
    controls = [homogeneous_counts(2*m, m) for m in range(32, 257, 32)]
    return {'status': 'PASS', 'author': author,
            'scope': 'Finite alias oracles and a proposed asymptotic cost proof; no full solver or exponent promotion.',
            'source_commit': source['commit'],
            'sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['oracle.py', 'audit.py', 'parameters.json']},
            'shared_inputs_sha256': {name: hashlib.sha256((CORE/name).read_bytes()).hexdigest()
                                     for name in ['certify.py', 'SOURCE.json']},
            'limits': {'oracle_max_n': 10, 'max_records_per_list': 16384, 'max_visited_pairs': 16384,
                       'campaign_weight_families': 5, 'campaign_max_n': 8, 'type_control_max_n': 512},
            'balanced_partition_controls': balanced, 'finite_campaign': finite_campaign(),
            'homogeneous_toy': homogeneous_toy(), 'homogeneous_type_controls': controls,
            'cost_profiles': cost_profiles()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True, choices=['conorholds'])
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run(args.author)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'finite_campaign': result['finite_campaign'],
                      'homogeneous_toy': result['homogeneous_toy'], 'cost_profiles': len(result['cost_profiles'])}))


if __name__ == '__main__':
    main()
