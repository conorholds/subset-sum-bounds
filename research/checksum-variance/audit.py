"""Bounded exact checksum experiments, failure controls, and symbolic budgets."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path

from oracle import MatchIndex, Scalar, Tables, character, mean, require, variance
import reference

HERE = Path(__file__).resolve().parent


def encoded(value):
    if isinstance(value, Scalar):
        return value.encoded()
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: encoded(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encoded(item) for item in value]
    return value


def multiplicity_control():
    table = Tables(4, [1], [0], 2, [(0,)], [[(None,), (None,), (0,)]],
                   [4, 2, 1, 1], [1, 3, 2, 2], [1], [1])
    index = table.index()
    values = table.distribution(0, [True]*4, index)
    actual = table.direct_checksum(0)
    require(index.total == 9 and mean(values, 4) == actual == Scalar(4, 1),
            'inverse-multiplicity control failed')
    omitted, inconsistent, uniform_groups = [], [], Scalar(4)
    for rank in range(index.total):
        h, low, _ = index.match(rank)
        k = 2*h+low
        omitted.append(table.fourier_term(k, 0)*Fraction(index.total, 4))
        # Invalid shortcut: retain repeated table entries but count unique labels in mu.
        unique = sum(table.passes(0, h, 0, r) and table.passes(1, low, 0, r)
                     for r in [0, 2])
        inconsistent.append(table.fourier_term(k, 0)*Fraction(index.total, 4*unique))
    for _, hs, ls in index.groups:
        for h, low in product(hs, ls):
            k = 2*h+low
            weight = Fraction(1, len(index.groups)*len(hs)*len(ls))
            uniform_groups += (table.fourier_term(k, 0)
                               *Fraction(index.total, 4*table.multiplicity(h, low))*weight)
    for wrong in [mean(omitted, 4), mean(inconsistent, 4), uniform_groups]:
        require(wrong != actual, 'intended biased control became unbiased')
    return encoded({'matches': index.total, 'multiplicities': [table.multiplicity(k//2, k%2) for k in range(4)],
                    'correct_mean': actual, 'one_draw_variance': variance(values, 4),
                    'omitted_multiplicity_mean': mean(omitted, 4),
                    'inconsistent_label_deduplication_mean': mean(inconsistent, 4),
                    'uniform_group_instead_of_match_mean': uniform_groups})


def independent_pattern_control():
    table = Tables(4, [1], [0], 2, [(0,)], [[(None,)]], [1]*4, [1]*4, [1], [1])
    first = table.distribution(2, [True, False, False, False])
    second = table.distribution(2, [False, True, True, True])
    independent = reference.product_mean_square([first, second], 4)
    shared = variance([a+b for a, b in zip(first, second)], 4)
    require(independent == variance(first, 4)+variance(second, 4) == Fraction(3, 2)
            and shared == 2, 'independent/shared draw control failed')
    return encoded({'own_frequency_sets': [[0], [1, 2, 3]], 'target': 2,
                    'independent_mean_squared_error': independent,
                    'reused_rank_mean_squared_error': shared,
                    'total_exact_checksum': table.direct_checksum(2)})


def phase_conditioning_control():
    probability = Fraction(0)
    joint_good = Fraction(0)
    amplitude_good = Fraction(0)
    missing_bias = []
    for phase, shift in product(range(2), repeat=2):
        table = Tables(2, [0], [phase], 1, [(0,)], [[(shift,)]], [3, 1], [3, 1], [1], [1])
        good = all(table.multiplicity(k, 0) > 0 for k in range(2))
        if good:
            probability += Fraction(1, 4)
            amplitude_good += table.fourier_term(0, 0).norm2()/4
            joint_good += table.f[(phase-shift) % 2]*table.g[-shift % 2]/4
        else:
            values = table.distribution(0, [True]*2)
            bias = mean(values, 2)-table.direct_checksum(0)
            if bias.norm2():
                missing_bias.append({'phase': phase, 'shift': shift, 'bias': bias})
    require(probability == Fraction(1, 4) and amplitude_good/probability == 4
            and joint_good/probability == Fraction(9, 4) and missing_bias,
            'phase-conditioning or uncovered-bias control failed')
    return encoded({'coverage_probability': probability,
                    'unconditional_amplitude_second_moment': Fraction(2),
                    'coverage_conditioned_amplitude_second_moment': amplitude_good/probability,
                    'unconditional_joint_score_moment': Fraction(1),
                    'coverage_conditioned_joint_score_moment': joint_good/probability,
                    'uncovered_nonzero_bias': missing_bias})


def check_configuration(table, counts, digest, joint=False):
    index = table.index()
    literal = reference.matches(table)
    require(index.status == 'PASS', 'unexpected small entry cap')
    compressed = sorted(index.match(rank) for rank in range(index.total))
    require(compressed == literal, 'compressed and whole-key references disagree')
    frequencies = Counter(h*table.base+low for h, low, _ in literal)
    for high, low in product(table.high_range, table.low_range):
        require(table.multiplicity(high, low) == frequencies[high*table.base+low],
                'retested product multiplicity differs from literal frequency count')
    expected = reference.subset_checksums(table.modulus, table.residues, table.phases)
    for target in range(table.modulus):
        fourier = sum((table.fourier_term(k, target) for k in range(table.modulus)), Scalar(table.modulus))/table.modulus
        require(fourier == expected[target] == table.direct_checksum(target), 'character orthogonality failed')
        counts['fourier_target_controls'] += 1
    # This deterministic partition mimics the required unique pattern ownership;
    # it is not the source's entropy grid or an asymptotic coverage experiment.
    owner = [sum((phase+k*residue) % table.modulus for phase, residue in zip(table.phases, table.residues)) % 2
             for k in range(table.modulus)]
    distributions = []
    total_mean = Scalar(table.modulus)
    for pattern in [0, 1]:
        own = [x == pattern for x in owner]
        values = table.distribution(0, own, index)
        center = mean(values, table.modulus)
        covered = sum((table.fourier_term(k, 0) for k in range(table.modulus)
                       if own[k] and frequencies[k]), Scalar(table.modulus))/table.modulus
        require(center == covered, 'sample mean differs from exactly covered contribution')
        second = sum((value.norm2() for value in values), Fraction(0))/len(values) if values else Fraction(0)
        formula = Fraction(index.total, table.modulus**2)*sum(
            (table.fourier_term(k, 0).norm2()/frequencies[k] for k in range(table.modulus)
             if own[k] and frequencies[k]), Fraction(0))
        require(second == formula, 'exact inverse-multiplicity second moment failed')
        counts['pattern_second_moments'] += 1
        total_mean += center
        distributions.append(values)
    good = all(frequencies[k] for k in range(table.modulus))
    counts['covered_configurations'] += good
    if good:
        require(total_mean == expected[0], 'covered complete estimate is biased')
    if joint:
        actual = reference.product_mean_square(distributions, table.modulus)
        require(actual == sum((variance(values, table.modulus) for values in distributions), Fraction(0)),
                'product-space errors do not add in variance')
        counts['independent_joint_draw_controls'] += 1
    counts['configurations'] += 1
    counts['literal_matches'] += len(literal)
    counts['out_of_range_matches'] += sum(k >= table.modulus for k in frequencies.elements())
    digest.update(json.dumps(encoded({'matches': literal, 'mean': total_mean,
                                     'second_moments': [variance(x, table.modulus) for x in distributions]}),
                             sort_keys=True, separators=(',', ':')).encode())


def finite_campaign(full=False):
    require(type(full) is bool, 'campaign selector must be boolean')
    counts = Counter()
    digest = hashlib.sha256()
    # Different exact circles, a nondividing digit base, repeated and all-star
    # tables, zero amplitudes, and positive rational asymmetric tests.
    for modulus in [2, 3, 4, 6]:
        residue_pairs = [(0, 0), (0, 1), (1, 1), (1, modulus-1)]
        for residues, phases, base in product(residue_pairs, product(range(modulus), repeat=2), [1, 2, 4, 8]):
            for shifts in [(None, None), (None, 0), (0, 1)]:
                labels = [[(shifts[0],), (shifts[1],)], [(None,), ((modulus-1),)]]
                table = Tables(modulus, residues, phases, base, [(0,), (1,)], labels,
                               list(range(1, modulus+1)), list(range(modulus, 0, -1)),
                               [1, Fraction(3, 4)], [Fraction(3, 4), 1])
                check_configuration(table, counts, digest, joint=counts['configurations'] % 17 == 0)
    if full:
        # Exhaust every residue vector, phase vector, and two one-coordinate
        # table rows per block on Z/4Z. Duplicate labels retain their indices.
        for residues, phases, flat in product(product(range(4), repeat=2),
                                              product(range(4), repeat=2),
                                              product([None, 0, 1, 2, 3], repeat=4)):
            labels = [[(flat[0],), (flat[1],)], [(flat[2],), (flat[3],)]]
            table = Tables(4, residues, phases, 2, [(0,), (1,)], labels,
                           [4, 2, 1, 1], [1, 3, 2, 2], [1, 1], [1, 1])
            check_configuration(table, counts, digest)
    return dict(counts, trace_sha256=digest.hexdigest())


def budget_profile(item):
    require(isinstance(item, dict) and set(item) == {'id', 'source', 'parameters', 'multipliers'}
            and isinstance(item['id'], str) and item['id'], 'invalid parameter snapshot')
    require(isinstance(item['parameters'], dict)
            and set(item['parameters']) == {'eta', 'tau', 'epsilon', 'moment_ceiling', 'alpha'}
            and all(isinstance(x, str) for x in item['parameters'].values()),
            'exact decimal-string checksum parameters required')
    require(isinstance(item['multipliers'], list) and len(item['multipliers']) == 2
            and all(isinstance(x, str) for x in item['multipliers']), 'two exact multiplier strings required')
    p = {key: Fraction(value) for key, value in item['parameters'].items()}
    eta, tau, epsilon, chi = [p[key] for key in ['eta', 'tau', 'epsilon', 'moment_ceiling']]
    require(0 < tau < eta < 2*tau and epsilon > 0 and 0 < chi < eta,
            'invalid checksum parameter domain')
    ratio = (eta-tau)/(2*tau-eta)
    lam1, lam2 = map(Fraction, item['multipliers'])
    require(lam1 >= 1 and lam2 >= 1 and 30 >= 3*lam1 and 30 >= (4*ratio+3)*lam2,
            'checksum domination or rounding budget fails')
    reserve = eta-chi-48*epsilon
    require(reserve > 0, 'no finite-circle variance reserve')
    # Half the remaining reserve is allowed for finite-circle convergence.
    finite_ceiling = chi+reserve/2
    require(eta-finite_ceiling-48*epsilon == reserve/2, 'variance reserve split failed')
    return encoded({'id': item['id'], 'continuous_ceiling': chi,
                    'finite_ceiling': finite_ceiling, 'variance_decay': reserve/2,
                    'first_rounding_margin': 30-3*lam1,
                    'second_rounding_margin': 30-(4*ratio+3)*lam2,
                    'pattern_factor_in_mse': 1, 'prefactor_upper_bound': 8})


def budget_profiles():
    inputs = json.loads((HERE/'INPUTS.json').read_text())
    rows = inputs['parameter_witnesses']
    require(len(rows) == 16 and len({row['id'] for row in rows}) == len(rows),
            'sixteen distinct source-bound snapshots required')
    return [budget_profile(item) for item in rows]


def input_hashes():
    return {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ['oracle.py', 'reference.py', 'audit.py', 'INPUTS.json']}


def run(author, full=False):
    require(author == 'conorholds', 'approved author handle required')
    require(type(full) is bool, 'campaign selector must be boolean')
    return {'author': author, 'status': 'PASS',
            'scope': 'Exact finite checksum controls and symbolic conditional variance budgets; no algorithmic record.',
            'campaign': 'full' if full else 'standard', 'inputs_sha256': input_hashes(),
            'finite_campaign': finite_campaign(full), 'multiplicity_control': multiplicity_control(),
            'independent_pattern_control': independent_pattern_control(),
            'phase_conditioning_control': phase_conditioning_control(), 'budget_profiles': budget_profiles()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.author, args.full)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'campaign': result['campaign'],
                      'finite_campaign': result['finite_campaign']}))


if __name__ == '__main__':
    main()
