"""Exact rational bounds for a necessary block-count inequality in the budget model."""
import argparse
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORE = HERE.parent/'filter-budgets'


def logarithm_bounds(x, terms=24):
    """Enclose ln(x) by the positive atanh series and a geometric tail bound."""
    if not isinstance(x, Fraction) or not 1 < x <= 2:
        raise ValueError('require a rational argument in (1,2]')
    if type(terms) is not int or not 1 <= terms <= 64:
        raise ValueError('require 1 through 64 series terms')
    t = (x-1)/(x+1)
    lower = 2*sum((t**(2*k+1)/(2*k+1) for k in range(terms)), Fraction(0))
    tail = 2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
    return lower, lower+tail


def entropy_bounds():
    l15, u15 = logarithm_bounds(Fraction(3, 2))
    l2, u2 = logarithm_bounds(Fraction(2))
    # log2(3) = 1 + ln(3/2)/ln(2); h(1/4) = 2 - (3/4) log2(3).
    return 2-Fraction(3, 4)*(1+u15/l2), 2-Fraction(3, 4)*(1+l15/u2)


def alpha_floor(blocks):
    if type(blocks) is not int or blocks <= 0:
        raise ValueError('positive integer block count required')
    lo, hi = entropy_bounds()
    return Fraction(1, 2)-(1-lo)/(blocks+1), Fraction(1, 2)-(1-hi)/(blocks+1)


def decimal_interval(bounds):
    output = []
    for value, rounding in zip(bounds, [ROUND_FLOOR, ROUND_CEILING], strict=True):
        with localcontext() as context:
            context.prec = 20
            context.rounding = rounding
            output.append(str(Decimal(value.numerator)/Decimal(value.denominator)))
    return output


def certificate():
    lower, upper = alpha_floor(18)
    if not lower > Fraction(49, 100) or not 3**75 < 2**119 or not 3**3 > 2**4:
        raise ValueError('block-count exclusion did not pass')
    return {
        'status': 'PASS',
        'scope': 'Necessary condition inside the conditional budget model, for every filter.',
        'formula': 'alpha > 1/2 - (1-h(1/4))/(J+1)',
        'series_terms': 24,
        'quarter_entropy_interval': decimal_interval(entropy_bounds()),
        'first_excluded_blocks': 18,
        'alpha_floor_at_18': decimal_interval((lower, upper)),
        'gap_above_0.49': decimal_interval((lower-Fraction(49, 100), upper-Fraction(49, 100))),
        'elementary_entropy_check': '3**75 < 2**119 implies h(1/4) > 0.81',
        'entropy_below_one_check': '3**3 > 2**4 implies h(1/4) < 1',
        'proof_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'model_checker_sha256': hashlib.sha256((CORE/'certify.py').read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
