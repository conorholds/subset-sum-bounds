"""Replay the existing mutation search at the same fixed frontier parameters."""
import argparse
import hashlib
import json
from pathlib import Path

import optimize
from search import search

HERE = Path(__file__).resolve().parent


def compare(steps=4000):
    if type(steps) is not int or not 0 <= steps <= 4000:
        raise ValueError('mutation-control steps must be from 0 through 4000')
    return {
        'status': 'EXPLORATORY',
        'scope': 'Same starts, degrees and fixed parameters; one fixed mutation seed is not a universal method comparison.',
        'search_sha256': hashlib.sha256((optimize.CORE/'search.py').read_bytes()).hexdigest(),
        'baseline_sha256': hashlib.sha256((HERE/'baseline.json').read_bytes()).hexdigest(),
        'results': [search(optimize.PARAMETERS, steps=steps, degree=degree, seed=20261015, start='combined')
                    for degree in [10, 14]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps([{'degree': row['degree'], 'evaluations': row['evaluations'],
                       'candidate_max_log2_moment': row['candidate_max_log2_moment']}
                      for row in result['results']]))


if __name__ == '__main__':
    main()
