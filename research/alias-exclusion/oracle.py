"""Finite alias-layer oracle; raw generation is exhaustive and capped at n=10."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from itertools import combinations
from math import comb, isqrt


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prime(value):
    return (type(value) is int and 2 <= value <= 1009
            and all(value % divisor for divisor in range(2, isqrt(value)+1)))


def balanced_left_size(n, m, k):
    require(type(n) is int and 4 <= n <= 10 and type(m) is int and m >= 2
            and m % 2 == 0 and 2*m <= n and type(k) is int and 0 <= k <= m//2,
            'require 4 <= n <= 10, even 2 <= m <= n/2, and 0 <= k <= m/2')
    a, b = comb(m, k), comb(m, m//2-k)
    for size in range(n-m+1):
        left, right = (1 << size)*a, (1 << (n-m-size))*b
        if left <= 2*right and right <= 2*left:
            return size
    raise ValueError('no balanced ordinary partition')


def prepare(weights, target, large, block, flip, k, carry, small, residue):
    n, m = len(weights), len(block)
    size = balanced_left_size(n, m, k)
    require(all(type(w) is int for w in weights) and type(target) is int,
            'integer weights and target required')
    require(prime(large) and prime(small), 'primes at most 1009 required by the finite oracle')
    require(all(type(i) is int and 0 <= i < n for i in block) and len(set(block)) == m,
            'distinct block indices required')
    require(type(flip) is int and 0 <= flip < 1 << n and type(carry) is int and 0 <= carry <= n
            and type(residue) is int and 0 <= residue < small, 'invalid flip, carry, or residue')
    block = tuple(sorted(block))
    ordinary = [i for i in range(n) if i not in block]
    reduced = [w % large for w in weights]
    signed = [w*(-1 if flip & (1 << i) else 1) for i, w in enumerate(weights)]
    signed_reduced = [w*(-1 if flip & (1 << i) else 1) for i, w in enumerate(reduced)]
    integer_target = target-sum(w for i, w in enumerate(weights) if flip & (1 << i))
    reduced_target = target % large+carry*large-sum(w for i, w in enumerate(reduced) if flip & (1 << i))
    return {'n': n, 'm': m, 'block': block, 'left': ordinary[:size], 'right': ordinary[size:],
            'signed': signed, 'signed_reduced': signed_reduced,
            'integer_target': integer_target, 'reduced_target': reduced_target}


def raw_records(ordinary, block, cardinality, data):
    for mask in range(1 << len(ordinary)):
        outside = sum(1 << index for bit, index in enumerate(ordinary) if mask & (1 << bit))
        for selected in combinations(block, cardinality):
            record_mask = outside | sum(1 << index for index in selected)
            integer = sum(w for i, w in enumerate(data['signed']) if record_mask & (1 << i))
            reduced = sum(w for i, w in enumerate(data['signed_reduced']) if record_mask & (1 << i))
            yield (record_mask, reduced, integer)


def layer(weights, target, large, block, flip, k, carry, small, residue,
          list_cap=16384, pair_cap=16384):
    """Visit nonexact reduced matches via interval skipping, then test disjointness."""
    require(type(list_cap) is int and 0 <= list_cap <= 16384
            and type(pair_cap) is int and 0 <= pair_cap <= 16384, 'caps must be integers from 0 through 16384')
    data = prepare(weights, target, large, block, flip, k, carry, small, residue)
    result = {'status': 'PASS', 'left_length': 0, 'right_length': 0,
              'excluded_exact': 0, 'group_queries': 0, 'visited_pairs': [], 'aliases': []}
    lists = []
    for side, cardinality, desired in [('left', k, residue),
                                     ('right', data['m']//2-k, (data['reduced_target']-residue) % small)]:
        records = []
        for record in raw_records(data[side], data['block'], cardinality, data):
            if record[1] % small != desired:
                continue
            if len(records) == list_cap:
                result['status'] = side+'_cap'
                return result
            records.append(record)
            result[side+'_length'] += 1
        lists.append(records)
    groups = defaultdict(list)
    for record in lists[1]:
        groups[record[1]].append(record)
    keys = {}
    for reduced, records in groups.items():
        records.sort(key=lambda record: (record[2], record[0]))
        keys[reduced] = [record[2] for record in records]
    aliases = set()
    for left in lists[0]:
        needed = data['reduced_target']-left[1]
        records = groups.get(needed, [])
        integers = keys.get(needed, [])
        result['group_queries'] += 1
        forbidden = data['integer_target']-left[2]
        lo, hi = bisect_left(integers, forbidden), bisect_right(integers, forbidden)
        result['excluded_exact'] += hi-lo
        for start, stop in [(0, lo), (hi, len(records))]:
            for index in range(start, stop):
                if len(result['visited_pairs']) == pair_cap:
                    result['status'] = 'pair_cap'
                    result['aliases'] = sorted(aliases)
                    return result
                right = records[index]
                result['visited_pairs'].append((left[0], right[0]))
                if left[0] & right[0]:
                    continue
                mask = (left[0] | right[0]) ^ flip
                value = sum(w for i, w in enumerate(weights) if mask & (1 << i))
                require(value != target and value % large == target % large, 'invalid extracted alias')
                aliases.add(mask)
    result['aliases'] = sorted(aliases)
    return result


def exhaustive_layer(weights, target, large, block, flip, k, carry, small, residue):
    """Separate full-mask enumeration of every raw pair, with exact diagnostics."""
    data = prepare(weights, target, large, block, flip, k, carry, small, residue)
    block_mask = sum(1 << index for index in block)
    records = []
    for side, count, desired in [('left', k, residue),
                                 ('right', data['m']//2-k, (data['reduced_target']-residue) % small)]:
        allowed = block_mask | sum(1 << index for index in data[side])
        selected = []
        for mask in range(1 << data['n']):
            if mask & ~allowed or (mask & block_mask).bit_count() != count:
                continue
            integer = sum(w for i, w in enumerate(data['signed']) if mask & (1 << i))
            reduced = sum(w for i, w in enumerate(data['signed_reduced']) if mask & (1 << i))
            if reduced % small == desired:
                selected.append((mask, reduced, integer))
        records.append(selected)
    result = {'left_length': len(records[0]), 'right_length': len(records[1]),
              'excluded_exact': 0, 'excluded_exact_overlapping': 0, 'visited_pairs': [], 'aliases': []}
    aliases = set()
    for left in records[0]:
        for right in records[1]:
            if left[1]+right[1] != data['reduced_target']:
                continue
            if left[2]+right[2] == data['integer_target']:
                result['excluded_exact'] += 1
                result['excluded_exact_overlapping'] += bool(left[0] & right[0])
                continue
            result['visited_pairs'].append((left[0], right[0]))
            if not left[0] & right[0]:
                aliases.add((left[0] | right[0]) ^ flip)
    result['aliases'] = sorted(aliases)
    return result
