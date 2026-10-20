"""Conditional budgets using fractional sumset coverage rather than J blocks.

The other budget inequalities are adapted unchanged from the frozen
filter-budgets/certify.py. This explicitly different model has its own receipt.
"""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CORE = HERE.parent / 'filter-budgets'
CHERNOFF = HERE.parent / 'chernoff-multipliers'
sys.path.insert(0, str(CORE))
from certify import q, require, lower, upper, log2


def budget_margins(parameters):
    expected = {'eta', 'tau', 'rho', 'deficiency', 'local_deficiency', 'gamma',
                'guard', 'preliminary', 'alpha', 'epsilon', 'moment_ceiling'}
    require(isinstance(parameters, dict) and set(parameters) == expected,
            'unexpected or missing fractional-coverage parameters')
    p = {key: q(value) for key, value in parameters.items()}
    eta, tau, rho = [p[key] for key in ['eta', 'tau', 'rho']]
    d, delta, gamma = [p[key] for key in ['deficiency', 'local_deficiency', 'gamma']]
    require(0 < lower(tau) and upper(tau) < lower(eta) and upper(eta) < lower(2*tau),
            'require 0 < tau < eta < 2 tau')
    require(0 < lower(rho) and upper(rho) < q('.5'), 'require 0 < rho < 1/2')
    require(0 < lower(d) and upper(d) < q('.5') and 0 < lower(delta) and 0 < lower(gamma),
            'positive deficiency and representation parameters required')
    require(lower(p['guard']) > 0 and lower(p['epsilon']) > 0, 'positive guard and epsilon required')
    h = -q('.25')*log2(q('.25'))-q('.75')*log2(q('.75'))
    hrho = -rho*log2(rho)-(1-rho)*log2(1-rho)
    a = (eta-tau)/(2*tau-eta)
    margins = {
        'preliminary': p['preliminary']-(q('.5')-d),
        'fractional_coverage': rho*(delta-d),
        'representations': q('.5')-delta-gamma,
        'sparse_distinctions': eta-(hrho+rho+p['guard']),
        'generation': tau-(rho+q('.25')),
        'filtered_lists': tau-(q('.5')+(h-q('.5')-gamma)*rho),
        'visited_pairs': tau-(1+(2*h-1-gamma)*rho-eta+p['guard']),
        'rounding_allowance': q('30')-q('8.3')*(4*a+3),
        'moment_reserve': eta-p['moment_ceiling']-48*p['epsilon'],
        'assembly_preliminary': p['alpha']-p['preliminary'],
        'assembly_checksum': p['alpha']-(tau+20*p['epsilon']),
        'assembly_primes': p['alpha']-eta/2,
        'below_meet_in_middle': q('.5')-p['alpha'],
    }
    for name, margin in margins.items():
        require(lower(margin) > 0, 'nonpositive margin: '+name)
    return margins
