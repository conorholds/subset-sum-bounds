"""Check evidence binding and selected-result consistency; optionally replay all witnesses."""
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CORE = Path('research/filter-budgets')
BOUND_FILES = {'certify.py', 'verify.py', 'SOURCE.json', 'requirements.txt'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project_file(root, relative):
    path = Path(relative)
    require(not path.is_absolute(), 'project paths must be relative')
    resolved = (root/path).resolve()
    require(resolved.is_relative_to(root.resolve()) and resolved.is_file(),
            'missing or escaping project file: '+str(relative))
    return resolved


def check_receipt(root, witness_path, report):
    data = read_json(witness_path)['witnesses']
    ids = [w['id'] for w in data]
    require(ids and len(ids) == len(set(ids)), 'nonempty, unique witness identifiers required')
    require(report['status'] == 'PASS', 'receipt did not pass')
    require(report['witnesses_sha256'] == digest(witness_path), 'stale witness hash')
    require(set(report['package_sha256']) == BOUND_FILES, 'incomplete checker binding')
    for name, expected in report['package_sha256'].items():
        require(digest(root/CORE/name) == expected, 'stale checker/source binding: '+name)
    rows = report['results']
    require(len(rows) == len(ids) and {r['id'] for r in rows} == set(ids),
            'receipt does not cover exactly the submitted witnesses')
    require(all(row['status'] == 'PASS' for row in rows), 'a witness did not pass')
    configurations = set()
    for row in rows:
        points, digits = row['moments']['grid'], row['moments']['interval_digits']
        require(type(points) is int and points >= 64 and points % 2 == 0
                and type(digits) is int and digits >= 30, 'invalid replay precision')
        configurations.add((points, digits))
    require(len(configurations) == 1, 'one replay configuration required per witness package')
    return configurations.pop()


def check_selection(root, selected):
    require(selected['schema_version'] == 1, 'unknown selected-result schema')
    require(selected['validation_level'] == 'interval-checked-conditional',
            'selected result must preserve its conditional numerical scope')
    require(selected['independent_review'] in {'pending', 'completed'}, 'unknown review status')
    for name in ['model', 'review', 'source_manifest', 'certificate']:
        project_file(root, selected[name])
    witness_path = project_file(root, selected['witnesses'])
    source_path = project_file(root, selected['source_manifest'])
    require(source_path == (root/CORE/'SOURCE.json').resolve(), 'source manifest must match the checker')
    require(read_json(source_path)['commit'] == selected['source_commit'], 'source commit mismatch')
    matching = [w for w in read_json(witness_path)['witnesses'] if w['id'] == selected['witness_id']]
    require(len(matching) == 1, 'selected witness is missing or ambiguous')
    alpha = selected['alpha']
    require(isinstance(alpha, str) and Decimal(alpha).is_finite(), 'exact selected alpha required')
    require(alpha == matching[0]['parameters']['alpha'] and 0 < Decimal(alpha) < Decimal('.5'),
            'selected exponent does not match the witness')
    report = read_json(project_file(root, selected['certificate']))
    check_receipt(root, witness_path, report)


def check_links(root):
    paths = [*root.glob('*.md'), *(root/'docs').rglob('*.md'),
             *(root/'research').rglob('*.md'), *(root/'.github').rglob('*.md')]
    for path in paths:
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            if target.startswith(('https://', 'http://', 'mailto:', '#')):
                continue
            relative = target.strip('<>').split('#', 1)[0]
            destination = (path.parent/relative).resolve()
            require(destination.is_relative_to(root.resolve()) and destination.is_file(),
                    f'broken local link in {path.relative_to(root)}: {target}')
    return len(paths)


def check_project(root, replay=False):
    selected = read_json(root/'certificates/selected-result.json')
    check_selection(root, selected)
    packages = sorted((root/'research').glob('*/witnesses.json'))
    require(packages, 'no numerical witness packages found')
    require(project_file(root, selected['witnesses']) in [p.resolve() for p in packages],
            'selected witness package is outside automatic replay discovery')
    count = 0
    for witness_path in packages:
        receipt_path = witness_path.with_name('CERTIFICATE.json')
        points, digits = check_receipt(root, witness_path, read_json(receipt_path))
        count += len(read_json(witness_path)['witnesses'])
        if replay:
            output = root/'.context/verification'/f'{witness_path.parent.name}.json'
            output.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run([sys.executable, '-B', str(root/CORE/'verify.py'),
                            '--witnesses', str(witness_path), '--output', str(output),
                            '--points', str(points), '--digits', str(digits)],
                           cwd=root, check=True)
            check_receipt(root, witness_path, read_json(output))
    documents = check_links(root)
    return {'status': 'PASS', 'packages': len(packages), 'witnesses': count,
            'documents': documents, 'replayed': replay,
            'alpha': selected['alpha'], 'independent_review': selected['independent_review'],
            'scope': 'Conditional numerical model; full algorithmic proof remains outside this check.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    try:
        result = check_project(ROOT, args.replay)
    except (ValueError, KeyError, TypeError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'FAIL: {error}\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
