"""Replay every retained sufficient-inequality witness using interval arithmetic."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import mpmath

from certify import certify

HERE = Path(__file__).resolve().parent


def load_witnesses(path=HERE/'witnesses.json'):
    data = json.loads(path.read_text())
    witnesses = data['witnesses']
    if not isinstance(witnesses, list) or not witnesses:
        raise ValueError('nonempty witness list required')
    for witness in witnesses:
        if set(witness) != {'id', 'design', 'parameters'} or not isinstance(witness['id'], str) or not witness['id']:
            raise ValueError('invalid witness schema')
    if len({w['id'] for w in witnesses}) != len(witnesses):
        raise ValueError('duplicate witness identifiers')
    return witnesses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--witnesses', type=Path, default=HERE/'witnesses.json')
    parser.add_argument('--points', type=int, default=128)
    parser.add_argument('--digits', type=int, default=50)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = [certify(w, args.points, args.digits) for w in load_witnesses(args.witnesses)]
    report = {'status': 'PASS', 'scope': 'Filter moments and rederived conditional budget model only.',
              'witnesses_sha256': hashlib.sha256(args.witnesses.read_bytes()).hexdigest(),
              'package_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                 for name in ['certify.py', 'verify.py', 'SOURCE.json', 'requirements.txt']},
              'runtime': {'python': platform.python_version(), 'mpmath': mpmath.__version__},
              'results': results}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': 'PASS', 'scope': report['scope'],
                      'results': [{'id': row['id'], 'status': row['status'],
                                   'log2_moments': row['moments']['log2_moment_intervals'],
                                   'seconds': round(row['seconds'], 3)} for row in results]}, sort_keys=True))


if __name__ == '__main__':
    main()
