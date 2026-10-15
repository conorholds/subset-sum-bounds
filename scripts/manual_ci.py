"""Bind owner-triggered verification to the exact commit being checked."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess

REPOSITORY = 'conorholds/subset-sum-bounds'
CHECK_NAME = 'Verify certificates'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def commit_sha(value):
    require(isinstance(value, str) and re.fullmatch(r'[0-9a-f]{40}', value),
            'an exact commit SHA is required')
    return value


def owner_context(env):
    require(env.get('GITHUB_REPOSITORY') == REPOSITORY
            and env.get('GITHUB_EVENT_NAME') == 'workflow_dispatch'
            and env.get('GITHUB_REF') == 'refs/heads/main'
            and env.get('GITHUB_ACTOR_ID') == '86416572'
            and env.get('GITHUB_TRIGGERING_ACTOR') == 'conorholds',
            'only conorholds may run or rerun verification from main')
    run_id, attempt = env.get('GITHUB_RUN_ID', ''), env.get('GITHUB_RUN_ATTEMPT', '')
    require(re.fullmatch(r'[1-9][0-9]*', run_id)
            and re.fullmatch(r'[1-9][0-9]*', attempt), 'invalid workflow run identity')
    return (f'manual-verification:{run_id}',
            f'https://github.com/{REPOSITORY}/actions/runs/{run_id}/attempts/{attempt}')


def github_api(method, path, payload=None):
    command = ['gh', 'api', '--method', method, f'repos/{REPOSITORY}/{path}']
    if payload is not None:
        command += ['--input', '-']
    result = subprocess.run(command, input=None if payload is None else json.dumps(payload),
                            text=True, capture_output=True, timeout=30, check=True)
    return json.loads(result.stdout)


def start(pr, env, api=github_api):
    identity, url = owner_context(env)
    if pr:
        require(re.fullmatch(r'[1-9][0-9]*', pr), 'PR must be a positive number')
        pull = api('GET', f'pulls/{pr}')
        require(pull['state'] == 'open' and pull['base']['ref'] == 'main'
                and pull['base']['repo']['full_name'] == REPOSITORY,
                'PR must be open and target this repository main')
        sha, target = commit_sha(pull['head']['sha']), f'PR #{pr}'
    else:
        sha, target = commit_sha(env.get('GITHUB_SHA')), 'main'
    check = api('POST', 'check-runs', {
        'name': CHECK_NAME, 'head_sha': sha, 'status': 'in_progress',
        'external_id': identity, 'details_url': url,
        'output': {'title': f'Verifying {target}',
                   'summary': f'Manually requested by conorholds for commit `{sha}`.'}})
    require(type(check['id']) is int and check['id'] > 0, 'invalid check identity')
    return {'sha': sha, 'check_id': str(check['id'])}


def finish(check_id, sha, result, env, api=github_api):
    identity, url = owner_context(env)
    require(re.fullmatch(r'[1-9][0-9]*', check_id), 'invalid check identity')
    sha = commit_sha(sha)
    require(result in {'success', 'failure', 'cancelled', 'skipped'}, 'invalid job result')
    path = f'check-runs/{check_id}'
    check = api('GET', path)
    require(check['head_sha'] == sha and check['name'] == CHECK_NAME
            and check['external_id'] == identity, 'check does not belong to this commit and run')
    conclusion = 'failure' if result == 'skipped' else result
    summary = (f'Commit: `{sha}`\n\nVerification: **{conclusion}**. '
               f'[Logs and interval receipts]({url}).')
    api('PATCH', path, {
        'status': 'completed', 'conclusion': conclusion, 'details_url': url,
        'completed_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'output': {'title': f'Manual verification: {conclusion}', 'summary': summary}})
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    begin = commands.add_parser('start')
    begin.add_argument('--pr', default='')
    complete = commands.add_parser('finish')
    complete.add_argument('--check-id', required=True)
    complete.add_argument('--sha', required=True)
    complete.add_argument('--result', required=True)
    args = parser.parse_args()
    try:
        if args.command == 'start':
            output = start(args.pr, os.environ)
            with Path(os.environ['GITHUB_OUTPUT']).open('a') as stream:
                stream.writelines(f'{key}={value}\n' for key, value in output.items())
            print(f'Verifying commit {output["sha"]}')
        else:
            summary = finish(args.check_id, args.sha, args.result, os.environ)
            with Path(os.environ['GITHUB_STEP_SUMMARY']).open('a') as stream:
                stream.write(summary+'\n')
    except subprocess.CalledProcessError as error:
        parser.exit(1, f'FAIL: GitHub API request failed: {error.stderr or error}\n')
    except (ValueError, KeyError, TypeError, OSError, subprocess.SubprocessError) as error:
        parser.exit(1, f'FAIL: {error}\n')


if __name__ == '__main__':
    main()
