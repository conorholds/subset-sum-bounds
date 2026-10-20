"""Exact bounded checksum tables and compressed match indexing.

Finite controls use moduli 2, 3, 4, and 6, whose characters have quadratic
coordinates. These are small execution oracles, not a Subset Sum solver.
"""
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import prod


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, low, high):
    return type(value) is int and low <= value <= high


@dataclass(frozen=True)
class Scalar:
    """a+b*zeta, where zeta^2=t*zeta-1, in one fixed circle."""
    modulus: int
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self):
        require(self.modulus in (2, 3, 4, 6) and type(self.modulus) is int,
                'exact circle must have modulus 2, 3, 4, or 6')
        require(type(self.a) in (int, Fraction) and type(self.b) in (int, Fraction),
                'exact rational coordinates required')
        require(self.modulus != 2 or self.b == 0, 'the binary circle is real')
        object.__setattr__(self, 'a', Fraction(self.a))
        object.__setattr__(self, 'b', Fraction(self.b))

    @property
    def trace(self):
        return {2: 0, 3: -1, 4: 0, 6: 1}[self.modulus]

    def coerce(self, other):
        if type(other) in (int, Fraction):
            return Scalar(self.modulus, other)
        require(isinstance(other, Scalar) and other.modulus == self.modulus,
                'incompatible exact circles')
        return other

    def __add__(self, other):
        other = self.coerce(other)
        return Scalar(self.modulus, self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Scalar(self.modulus, -self.a, -self.b)

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __mul__(self, other):
        other = self.coerce(other)
        return Scalar(self.modulus, self.a*other.a-self.b*other.b,
                      self.a*other.b+self.b*other.a+self.trace*self.b*other.b)

    __rmul__ = __mul__

    def __truediv__(self, other):
        require(type(other) in (int, Fraction) and other != 0,
                'nonzero rational divisor required')
        return Scalar(self.modulus, self.a/other, self.b/other)

    def conjugate(self):
        return Scalar(self.modulus, self.a+self.trace*self.b, -self.b)

    def norm2(self):
        return self.a*self.a+self.trace*self.a*self.b+self.b*self.b

    def encoded(self):
        return [str(self.a), str(self.b)]


def character(modulus, exponent):
    require(integer(exponent, -10**6, 10**6), 'bounded integer character exponent required')
    one = Scalar(modulus, 1)
    root = Scalar(modulus, -1) if modulus == 2 else Scalar(modulus, 0, 1)
    value = one
    for _ in range(exponent % modulus):
        value = value*root
    return value


class MatchIndex:
    """Group equal table keys and index their products without listing pairs."""
    def __init__(self, high, low, cap=65536):
        require(integer(cap, 0, 65536), 'entry cap must be an integer at most 65536')
        self.status = 'PASS'
        sides = []
        for name, entries in [('high', high), ('low', low)]:
            saved = []
            for entry in entries:
                if len(saved) == cap:
                    self.status = name+'_cap'
                    self.high = self.low = self.groups = self.prefix = ()
                    self.total = 0
                    return
                digit, key = entry
                require(integer(digit, 0, 65535) and isinstance(key, tuple)
                        and len(key) <= 8 and all(integer(x, 0, 65535) for x in key),
                        'invalid digit or table-index key')
                saved.append((digit, key))
            require(len(set(saved)) == len(saved), 'duplicate digit/table entry')
            sides.append(tuple(saved))
        self.high, self.low = sides
        grouped = []
        for entries in sides:
            by_key = defaultdict(list)
            for digit, key in entries:
                by_key[key].append(digit)
            grouped.append(by_key)
        self.groups, self.prefix = [], []
        total = 0
        for key in sorted(set(grouped[0]) & set(grouped[1])):
            hs, ls = tuple(sorted(grouped[0][key])), tuple(sorted(grouped[1][key]))
            total += len(hs)*len(ls)
            self.groups.append((key, hs, ls))
            self.prefix.append(total)
        self.total = total

    def match(self, rank):
        require(self.status == 'PASS' and integer(rank, 0, self.total-1),
                'rank outside the nonaborted match index')
        group = bisect_right(self.prefix, rank)
        offset = rank-(self.prefix[group-1] if group else 0)
        key, hs, ls = self.groups[group]
        h, ell = divmod(offset, len(ls))
        return hs[h], ls[ell], key

    def enumerate(self, cap=65536):
        require(self.status == 'PASS' and integer(cap, 0, 65536), 'invalid enumeration cap')
        result = []
        for rank in range(min(cap, self.total)):
            result.append(self.match(rank))
        return ('PASS' if self.total <= cap else 'match_cap'), result


class Tables:
    """Literal small block tables with exact likelihood-product thresholds."""
    def __init__(self, modulus, residues, phases, base, blocks, labels, f, g,
                 high_thresholds, low_thresholds):
        require(type(modulus) is int and modulus in (2, 3, 4, 6), 'unsupported toy modulus')
        n = len(residues)
        require(1 <= n <= 6 and len(phases) == n
                and all(integer(x, 0, modulus-1) for x in (*residues, *phases)),
                'one to six valid residues and phases required')
        require(integer(base, 1, 8) and base & (base-1) == 0, 'base must be a power of two at most eight')
        require(1 <= len(blocks) <= n and sorted(i for block in blocks for i in block) == list(range(n))
                and all(block and all(type(i) is int for i in block) for block in blocks),
                'blocks must partition the coordinates')
        require(len(labels) == len(blocks) and all(1 <= len(rows) <= 4 for rows in labels),
                'one to four table rows per block required')
        for block, rows in zip(blocks, labels):
            require(all(len(row) == len(block) and all(x is None or integer(x, 0, modulus-1)
                        for x in row) for row in rows), 'invalid table labels')
        require(len(f) == len(g) == modulus and all(type(x) in (int, Fraction) and x > 0 for x in (*f, *g)),
                'positive exact test masses required')
        require(len(high_thresholds) == len(low_thresholds) == len(blocks)
                and all(type(x) in (int, Fraction) and x > 0 for x in (*high_thresholds, *low_thresholds)),
                'positive exact block thresholds required')
        self.modulus, self.n, self.base = modulus, n, base
        self.residues, self.phases = tuple(residues), tuple(phases)
        self.blocks = tuple(tuple(block) for block in blocks)
        self.labels = tuple(tuple(tuple(row) for row in rows) for rows in labels)
        self.f = tuple(Fraction(modulus*x, sum(f)) for x in f)
        self.g = tuple(Fraction(modulus*x, sum(g)) for x in g)
        self.thresholds = (tuple(map(Fraction, high_thresholds)), tuple(map(Fraction, low_thresholds)))
        self.high_range = range((modulus+base-1)//base)
        self.low_range = range(base)

    def passes(self, side, digit, block, row):
        require(side in (0, 1), 'invalid digit side')
        density = self.f if side == 0 else self.g
        terms = []
        for coordinate, label in zip(self.blocks[block], self.labels[block][row]):
            if label is None:
                continue
            vector = (self.phases[coordinate]+digit*self.base*self.residues[coordinate]
                      if side == 0 else -digit*self.residues[coordinate])
            terms.append(density[(vector-label) % self.modulus])
        return prod(terms, start=Fraction(1)) >= self.thresholds[side][block]

    def entries(self, side):
        for digit in (self.high_range if side == 0 else self.low_range):
            accepted = [tuple(r for r in range(len(rows)) if self.passes(side, digit, j, r))
                        for j, rows in enumerate(self.labels)]
            for key in product(*accepted):
                yield digit, key

    def index(self, cap=65536):
        return MatchIndex(self.entries(0), self.entries(1), cap)

    def multiplicity(self, high, low):
        return prod(sum(self.passes(0, high, j, r) and self.passes(1, low, j, r)
                        for r in range(len(rows))) for j, rows in enumerate(self.labels))

    def fourier_term(self, frequency, target):
        require(integer(frequency, 0, self.modulus-1) and integer(target, 0, self.modulus-1),
                'frequency and target must be valid residues')
        value = character(self.modulus, -frequency*target)
        for phase, residue in zip(self.phases, self.residues):
            value = value*(character(self.modulus, 0)+character(self.modulus, phase+frequency*residue))
        return value

    def direct_checksum(self, target):
        value = Scalar(self.modulus)
        for mask in range(1 << self.n):
            if sum(self.residues[i] for i in range(self.n) if mask & (1 << i)) % self.modulus == target:
                value += character(self.modulus, sum(self.phases[i] for i in range(self.n) if mask & (1 << i)))
        return value

    def distribution(self, target, own, index=None):
        """Exact one-draw values, including zero-valued invalid/wrong-pattern matches."""
        require(len(own) == self.modulus and all(type(x) is bool for x in own), 'boolean ownership map required')
        index = self.index() if index is None else index
        require(index.status == 'PASS' and index.total <= 65536, 'aborted or oversized toy distribution')
        values = []
        for rank in range(index.total):
            high, low, _ = index.match(rank)
            frequency = high*self.base+low
            value = Scalar(self.modulus)
            if frequency < self.modulus and own[frequency]:
                mu = self.multiplicity(high, low)
                require(mu > 0, 'sampled frequency has zero multiplicity')
                value = self.fourier_term(frequency, target)*Fraction(index.total, self.modulus*mu)
            values.append(value)
        return values


def mean(values, modulus):
    return sum(values, Scalar(modulus))/len(values) if values else Scalar(modulus)


def variance(values, modulus):
    center = mean(values, modulus)
    return sum(((value-center).norm2() for value in values), Fraction(0))/len(values) if values else Fraction(0)
