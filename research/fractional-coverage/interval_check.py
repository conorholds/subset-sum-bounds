"""Replay the explicit fractional-coverage model with unchanged moment bounds."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import mpmath

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import budget

# A distinct module name avoids aliasing another research package's checker.
import importlib.util
SPEC = importlib.util.spec_from_file_location('fractional_coverage_chernoff_bounds',
                                              budget.CHERNOFF / 'interval_check.py')
CHERNOFF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHERNOFF)

MODEL = 'fractional-coverage-v1'


def certify_candidate(candidate, points=128, digits=50):
    started = time.perf_counter()
    q, require, lower = budget.q, budget.require, budget.lower
    require(type(digits) is int and digits >= 30, 'at least 30 interval digits required')
    CHERNOFF.iv.dps = digits
    require(isinstance(candidate, dict) and set(candidate) == {'id', 'design', 'parameters', 'multipliers'},
            'unexpected or missing candidate fields')
    require(isinstance(candidate['id'], str) and candidate['id'], 'nonempty candidate id required')
    multipliers = candidate['multipliers']
    require(isinstance(multipliers, list) and len(multipliers) == 2, 'two exact multipliers required')
    lam1, lam2 = map(q, multipliers)
    require(lower(lam1) >= 1 and lower(lam2) >= 1, 'domination multipliers must be at least one')
    margins = budget.budget_margins(candidate['parameters'])
    eta, tau = [q(candidate['parameters'][key]) for key in ['eta', 'tau']]
    ratio = (eta-tau)/(2*tau-eta)
    margins['domination_one'] = 30-3*lam1
    margins['domination_two'] = 30-(4*ratio+3)*lam2
    for name in ['domination_one', 'domination_two']:
        require(lower(margins[name]) > 0, 'nonpositive margin: '+name)
    moments = CHERNOFF.enclose_moments(candidate['design'], candidate['parameters'], multipliers, points, digits)
    return {'id': candidate['id'], 'status': 'PASS', 'multipliers': multipliers,
            'margins': {key: CHERNOFF.enclosure(value) for key, value in margins.items()},
            'moments': moments, 'seconds': time.perf_counter()-started}


def input_hashes():
    files = {
        'package': (HERE, ['budget.py', 'search.py', 'interval_check.py', 'candidates.json', 'predecessor.json']),
        'chernoff': (budget.CHERNOFF, ['interval_check.py', 'multiplier_search.py']),
        'shared': (budget.CORE, ['certify.py', 'filters.py', 'SOURCE.json', 'requirements.txt']),
    }
    return {group: {name: hashlib.sha256((path/name).read_bytes()).hexdigest() for name in names}
            for group, (path, names) in files.items()}


def replay(points=128, digits=50, author='conorholds'):
    data = json.loads((HERE/'candidates.json').read_text())
    budget.require(set(data) == {'model', 'candidates'} and data['model'] == MODEL, 'unrecognized model')
    budget.require(isinstance(data['candidates'], list) and data['candidates'], 'nonempty candidates required')
    ids = [row['id'] for row in data['candidates']]
    budget.require(len(ids) == len(set(ids)), 'duplicate candidate ids')
    budget.require(isinstance(author, str) and author, 'explicit nonempty author required')
    return {'status': 'PASS', 'model': MODEL, 'author': author,
            'scope': 'Fractional sumset coverage and continuous conditional budgets; independent algorithm review pending.',
            'inputs_sha256': input_hashes(),
            'runtime': {'python': platform.python_version(), 'mpmath': mpmath.__version__},
            'results': [certify_candidate(row, points, digits) for row in data['candidates']]}


def validate_receipt(receipt):
    require = budget.require
    require(isinstance(receipt, dict) and receipt.get('status') == 'PASS' and receipt.get('model') == MODEL,
            'unrecognized or failing receipt')
    require(receipt.get('inputs_sha256') == input_hashes(), 'receipt inputs have changed')
    candidates = json.loads((HERE/'candidates.json').read_text())['candidates']
    expected = {row['id']: row['multipliers'] for row in candidates}
    rows = receipt.get('results')
    require(isinstance(rows, list) and len(rows) == len(expected), 'incomplete receipt coverage')
    require(all(isinstance(row, dict) and isinstance(row.get('id'), str) for row in rows), 'invalid receipt rows')
    require(len({row['id'] for row in rows}) == len(rows) and {row['id'] for row in rows} == set(expected),
            'receipt candidate coverage mismatch')
    configurations = set()
    for row in rows:
        require(row.get('status') == 'PASS' and row.get('multipliers') == expected[row['id']],
                'receipt status or multipliers mismatch')
        moments = row.get('moments')
        require(isinstance(moments, dict), 'missing receipt moments')
        points, digits = moments.get('grid'), moments.get('interval_digits')
        require(type(points) is int and points >= 64 and points % 2 == 0
                and type(digits) is int and digits >= 30, 'invalid receipt precision')
        configurations.add((points, digits))
    require(len(configurations) == 1, 'one receipt configuration required')
    return configurations.pop()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--author', required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check-receipt', type=Path)
    parser.add_argument('--points', type=int)
    parser.add_argument('--digits', type=int)
    args = parser.parse_args()
    points, digits = 128, 50
    if args.check_receipt:
        points, digits = validate_receipt(json.loads(args.check_receipt.read_text()))
    if args.points is not None:
        points = args.points
    if args.digits is not None:
        digits = args.digits
    report = replay(points, digits, args.author)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': report['status'], 'model': MODEL,
                      'results': [{'id': row['id'], 'log2_moments': row['moments']['log2_moment_intervals'],
                                   'seconds': round(row['seconds'], 3)} for row in report['results']]}))
