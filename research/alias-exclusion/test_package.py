"""Independent finite identities, exact exclusion controls, and cap failures."""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit
import oracle


class LayerTests(unittest.TestCase):
    def test_extracts_real_aliases_with_carries_and_flips(self):
        weights = [1, 2, 4, 8]
        carried = oracle.layer(weights, 5, 3, [0, 2], 0, 0, 1, 2, 0)
        self.assertEqual(carried['aliases'], [11, 14])
        flipped = oracle.layer(weights, 5, 3, [0, 2], 1, 0, 0, 2, 0)
        self.assertIn(2, flipped['aliases'])
        for result in [carried, flipped]:
            self.assertEqual(result['status'], 'PASS')
            self.assertTrue(set(result['aliases']) <= audit.masks(weights, 5, 3))

    def test_exact_overlapping_matches_are_skipped_before_the_zero_pair_cap(self):
        args = ([1]*8, 4, 11, [0, 1, 2, 3], 0, 1, 0, 3, 0)
        direct = oracle.exhaustive_layer(*args)
        result = oracle.layer(*args, pair_cap=0)
        self.assertGreater(direct['excluded_exact_overlapping'], 0)
        self.assertEqual(result['excluded_exact'], direct['excluded_exact'])
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['visited_pairs'], [])
        self.assertEqual(result['aliases'], [])

    def test_caps_fail_before_an_extra_filtered_record_or_pair_is_saved(self):
        args = ([1, 2, 4, 8], 5, 3, [0, 2], 0, 0, 1, 2, 0)
        for cap in [0, 1]:
            with self.subTest(cap=cap):
                result = oracle.layer(*args, pair_cap=cap)
                self.assertEqual(result['status'], 'pair_cap')
                self.assertEqual(len(result['visited_pairs']), cap)
                self.assertTrue(set(result['aliases']) <= {11, 14})
                limited = oracle.layer(*args, list_cap=cap)
                self.assertEqual(limited['status'], 'left_cap')
                self.assertEqual(limited['left_length'], cap)
                self.assertEqual(limited['visited_pairs'], [])
        right = oracle.layer([1, 2, 4, 8], 5, 3, [0, 2], 0, 0, 0, 2, 1, list_cap=1)
        self.assertEqual(right['status'], 'right_cap')
        self.assertEqual(right['right_length'], 1)
        fits = oracle.layer(*args, list_cap=4, pair_cap=2)
        self.assertEqual(fits['status'], 'PASS')
        self.assertEqual(len(fits['visited_pairs']), 2)

    def test_one_shared_residue_costs_one_small_prime_factor(self):
        weights, target, large, block, flip, k, carry, small = [1, 2, 4, 8], 5, 3, [0, 2], 0, 0, 1, 3
        data = oracle.prepare(weights, target, large, block, flip, k, carry, small, 0)
        left = list(oracle.raw_records(data['left'], data['block'], k, data))
        right = list(oracle.raw_records(data['right'], data['block'], 1-k, data))
        eligible = sum(a[1]+b[1] == data['reduced_target'] and a[2]+b[2] != data['integer_target']
                       for a in left for b in right)
        results = [oracle.layer(weights, target, large, block, flip, k, carry, small, residue)
                   for residue in range(small)]
        self.assertGreater(eligible, 0)
        self.assertEqual(sum(row['left_length'] for row in results), len(left))
        self.assertEqual(sum(row['right_length'] for row in results), len(right))
        self.assertEqual(sum(len(row['visited_pairs']) for row in results), eligible)

    def test_fixed_integer_difference_does_not_depend_on_prime_or_carry(self):
        weights, target, block, flip, k = [10, 17, 23, 30], 19, [0, 2], 5, 0
        differences = None
        matches = 0
        for large in [3, 5, 7, 11]:
            for carry in range(5):
                data = oracle.prepare(weights, target, large, block, flip, k, carry, 2, 0)
                left = list(oracle.raw_records(data['left'], data['block'], k, data))
                right = list(oracle.raw_records(data['right'], data['block'], 1-k, data))
                current = {(a[0], b[0]): a[2]+b[2]-data['integer_target'] for a in left for b in right}
                if differences is None:
                    differences = current
                self.assertEqual(current, differences)
                for a, b in product(left, right):
                    if a[1]+b[1] == data['reduced_target']:
                        self.assertEqual(current[a[0], b[0]] % large, 0)
                        matches += bool(current[a[0], b[0]])
        self.assertGreater(matches, 0)

    def test_invalid_inputs_and_campaign_expansion_fail(self):
        args = dict(weights=[1, 2, 4, 8], target=5, large=3, block=[0, 2], flip=0, k=0,
                    carry=0, small=2, residue=0)
        for changed in [{'large': 9}, {'small': 1}, {'block': [0, 0]}, {'flip': True},
                        {'carry': -1}, {'residue': 2}, {'target': 1.5}, {'weights': [True, 2, 4, 8]},
                        {'list_cap': -1}, {'pair_cap': 16385}]:
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                oracle.layer(**(args | changed))
        for n, m, k in [(11, 2, 0), (8, 3, 0), (4, 4, 0), (4, 2, 2), (4, 2, True)]:
            with self.subTest(n=n, m=m, k=k), self.assertRaises(ValueError):
                oracle.balanced_left_size(n, m, k)


class ObstructionTests(unittest.TestCase):
    def test_flip_probability_matches_every_flip_in_small_multiplicity_models(self):
        for n in range(2, 11, 2):
            for overlap in range(n+1):
                count = 0
                for flip in range(1 << n):
                    a = (flip & ((1 << overlap)-1)).bit_count()
                    b = (flip >> overlap).bit_count()
                    count += 2*overlap-3*a+b == n//2
                with self.subTest(n=n, overlap=overlap):
                    self.assertEqual(Fraction(count, 1 << n), audit.flip_equality_probability(n, overlap))

    def test_equal_weight_toy_has_no_aliases_but_exact_overlapping_work_without_exclusion(self):
        result = audit.homogeneous_toy()
        self.assertEqual(result['transformed_weight'], 1026)
        self.assertEqual(result['large_prime'], 101)
        self.assertEqual(audit.masks([1026]*8, 4104, 101), set())
        self.assertEqual(result['visited_after_exclusion'], 0)
        self.assertEqual(result['exact_overlapping_pairs_skipped'], 2688)
        self.assertEqual(result['expected_overlapping_visits_if_exclusion_delayed'], '7/2')

    def test_type_class_and_flip_lower_bounds_are_exact(self):
        for m in range(32, 257, 32):
            result = audit.homogeneous_counts(2*m, m)
            self.assertGreater(result['raw_overlapping_pairs'], 1 << (2*m))
            probability = Fraction(result['exact_integer_equality_probability_over_flips'])
            self.assertGreaterEqual(probability, Fraction(1, (2*m+1)**2))
            self.assertGreater(Fraction(result['expected_exact_overlapping_pairs_times_small_modulus']), 0)
        for n, m in [(64, 16), (62, 32), (65, 32), (1024, 32)]:
            with self.subTest(n=n, m=m), self.assertRaises(ValueError):
                audit.homogeneous_counts(n, m)
        for n, overlap in [(3, 1), (514, 2), (8, -1), (8, True)]:
            with self.subTest(n=n, overlap=overlap), self.assertRaises(ValueError):
                audit.flip_equality_probability(n, overlap)


class EvidenceTests(unittest.TestCase):
    def test_exact_receipt_replays_all_independent_controls(self):
        receipt = json.loads((HERE/'AUDIT.json').read_text())
        self.assertEqual(receipt, audit.run('conorholds'))
        self.assertEqual(receipt['finite_campaign']['layers'], 9960)
        self.assertEqual(receipt['finite_campaign']['balanced_recovery_controls'], 352)
        self.assertEqual(receipt['balanced_partition_controls'], 23)
        with self.assertRaises(ValueError):
            audit.run('unapproved')

    def test_all_retained_costs_fit_and_the_obstruction_removes_the_prime_gain(self):
        profiles = audit.cost_profiles()
        self.assertEqual(len(profiles), 8)
        inputs = json.loads((HERE/'parameters.json').read_text())
        parameters = {w['id']: w['parameters'] for source in inputs['sources'] for w in source['parameters']}
        for row in profiles:
            p = parameters[row['id']]
            tau = Fraction(p['tau'])
            for key in ['generation_exponent', 'expected_list_exponent', 'expected_nonexact_pair_exponent_with_guard']:
                self.assertLess(Fraction(row[key][1]), tau)
            gap = Fraction(row['delayed_exclusion_expected_overlap_exponent'][0])-Fraction(row['expected_nonexact_pair_exponent_with_guard'][1])
            self.assertGreater(gap, Fraction(p['eta'])-Fraction(p['guard'])-Fraction(1, 10**15))
            self.assertGreater(Fraction(row['delayed_exclusion_expected_overlap_exponent'][0]), 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
