"""Literal full-key and full-subset references, separate from grouped emission."""
from itertools import product
from fractions import Fraction

from oracle import Scalar, character


def matches(tables):
    """Visit the whole table universe; this is only the small reference."""
    result = []
    for key in product(*(range(len(rows)) for rows in tables.labels)):
        for high in tables.high_range:
            for low in tables.low_range:
                accepted = True
                for j, row in enumerate(key):
                    score_high = score_low = Fraction(1)
                    for coordinate, shift in zip(tables.blocks[j], tables.labels[j][row]):
                        if shift is not None:
                            yh = tables.phases[coordinate]+high*tables.base*tables.residues[coordinate]
                            zl = -low*tables.residues[coordinate]
                            score_high *= tables.f[(yh-shift) % tables.modulus]
                            score_low *= tables.g[(zl-shift) % tables.modulus]
                    if score_high < tables.thresholds[0][j] or score_low < tables.thresholds[1][j]:
                        accepted = False
                        break
                if accepted:
                    result.append((high, low, key))
    return sorted(result)


def subset_checksums(modulus, residues, phases):
    """Dynamic subset expansion, without using the Fourier product."""
    values = [Scalar(modulus) for _ in range(modulus)]
    values[0] = Scalar(modulus, 1)
    for residue, phase in zip(residues, phases):
        previous = values
        values = list(previous)
        multiplier = character(modulus, phase)
        for old in range(modulus):
            values[(old+residue) % modulus] += previous[old]*multiplier
    return values


def product_mean_square(distributions, modulus):
    """Exact product probability space, without a covariance shortcut."""
    centers = [sum(values, Scalar(modulus))/len(values) if values else Scalar(modulus)
               for values in distributions]
    choices = [values if values else [Scalar(modulus)] for values in distributions]
    total = Fraction(0)
    count = 0
    for draw in product(*choices):
        error = sum((value-center for value, center in zip(draw, centers)), Scalar(modulus))
        total += error.norm2()
        count += 1
    return total/count
