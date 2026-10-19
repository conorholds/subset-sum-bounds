"""Exact finite controls for the coverage change of measure.

On Z/3Z use f=(1/2,1/2,2), g=1 and kappa=1. A tested
coordinate has high score and log likelihood -1 or +1; low score
is zero. A star has both scores and log likelihood zero. Dynamic
programming exhausts the score distribution, including indexed labels.
This synthetic control is not the retained degree-twenty filter.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def distributions(alpha, sizes):
    if not isinstance(alpha, Fraction) or not 0 < alpha <= 1:
        raise ValueError('rational alpha in (0,1] required')
    if (not isinstance(sizes, tuple) or not sizes
            or any(type(d) is not int or not 1 <= d <= 512 for d in sizes)):
        raise ValueError('nonempty tuple of block sizes from one to 512 required')
    numerator, denominator = alpha.numerator, alpha.denominator
    if denominator > 64:
        raise ValueError('bounded probability denominator required')
    base = 3*denominator
    weights = [(2*numerator, 3*(denominator-numerator), numerator),
               (numerator, 3*(denominator-numerator), 2*numerator)]
    laws = [[1], [1]]
    wanted, result = set(sizes), {}
    for d in range(1, max(sizes)+1):
        for law_index, (negative, zero, positive) in enumerate(weights):
            updated = [0]*(len(laws[law_index])+2)
            for index, count in enumerate(laws[law_index]):
                updated[index] += negative*count
                updated[index+1] += zero*count
                updated[index+2] += positive*count
            laws[law_index] = updated
        if d in wanted:
            result[d] = (tuple(laws[0]), tuple(laws[1]), base**d)
    return result


def typical_mass(alpha, d, epsilon, tilted, denominator):
    mean = alpha*d/3
    return Fraction(sum(count for index, count in enumerate(tilted)
                        if abs(Fraction(index-d)-mean) <= epsilon*d), denominator)


def run(author):
    if author != 'conorholds':
        raise ValueError('approved author handle required')
    trace = hashlib.sha256()
    cases, identities, minimum = 0, 0, Fraction(1)
    sizes = (288, 320, 384)
    epsilon = Fraction(1, 8)
    # Direct outward rounding, independently of coverage_bounds.round_profile:
    # grid step 1/32, approximation error allowance 1/64.
    d_profile, a_profile, e_profile = Fraction(3, 8), Fraction(5, 16), Fraction(-1, 32)
    for alpha in [Fraction(1, 4), Fraction(1, 2), Fraction(1)]:
        for d, (original, tilted, denominator) in distributions(alpha, sizes).items():
            if sum(original) != denominator or sum(tilted) != denominator:
                raise AssertionError('distribution lost probability mass')
            for index, (p_count, q_count) in enumerate(zip(original, tilted)):
                score = index-d
                correct = (q_count == p_count*(1 << score) if score >= 0
                           else p_count == q_count*(1 << -score))
                if not correct:
                    raise AssertionError('incorrect coordinate-product likelihood law')
                identities += 1
            probability = typical_mass(alpha, d, epsilon, tilted, denominator)
            if probability < Fraction(1, 2):
                raise AssertionError('typical event has insufficient tilted mass')
            minimum = min(minimum, probability)
            accepted = Fraction(sum(count for index, count in enumerate(original)
                                    if Fraction(index-d)-epsilon*d >= (alpha*a_profile-3*epsilon)*d),
                                denominator)
            upper = (alpha*d_profile+2*epsilon)*d
            if upper.denominator != 1 or 0 < (alpha*e_profile-3*epsilon):
                raise AssertionError('invalid exact fixture threshold')
            if accepted < probability/(1 << int(upper)):
                raise AssertionError('change-of-measure acceptance lower bound failed')
            # Check the implication on every possible score, including a
            # worst-case -epsilon score-approximation allowance.
            for score in range(-d, d+1):
                if abs(Fraction(score)-alpha*d/3) <= epsilon*d:
                    if (Fraction(score)-epsilon*d < (alpha*a_profile-3*epsilon)*d
                            or score > upper):
                        raise AssertionError('typical score failed acceptance/likelihood implication')
            trace.update(json.dumps([str(alpha), d, original, tilted], separators=(',', ':')).encode())
            cases += 1
    _, tilted, denominator = distributions(Fraction(1, 2), (1,))[1]
    small = typical_mass(Fraction(1, 2), 1, Fraction(1, 32), tilted, denominator)
    if small != 0:
        raise AssertionError('small-block failure control no longer fails')
    return {'author': author, 'status': 'PASS', 'scope': 'Exact synthetic coverage-law controls, not a full checksum implementation.',
            'fixture': {'circle': 3, 'f': ['1/2', '1/2', '2'], 'g': ['1', '1', '1'],
                        'kappa': '1', 'alpha_values': ['1/4', '1/2', '1'], 'block_sizes': list(sizes),
                        'epsilon': str(epsilon)},
            'cases': cases, 'likelihood_identities': identities,
            'minimum_typical_mass_lower_4096': str(Fraction(minimum.numerator*4096//minimum.denominator, 4096)),
            'small_block_typical_mass': str(small), 'trace_sha256': trace.hexdigest(),
            'inputs_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                              for name in ['coverage_audit.py', 'coverage_bounds.py', 'INPUTS.json']}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run(args.author)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: result[key] for key in ['status', 'cases', 'likelihood_identities', 'trace_sha256']}))


if __name__ == '__main__':
    main()
