"""Authorization, commit binding, and failure controls for manual verification."""
import copy
from contextlib import redirect_stdout, redirect_stderr
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import manual_ci as ci

MAIN_SHA = '1'*40
PR_SHA = '2'*40


class ManualCITests(unittest.TestCase):
    def setUp(self):
        self.env = {
            'GITHUB_REPOSITORY': ci.REPOSITORY, 'GITHUB_EVENT_NAME': 'workflow_dispatch',
            'GITHUB_REF': 'refs/heads/main', 'GITHUB_ACTOR_ID': '86416572',
            'GITHUB_TRIGGERING_ACTOR': 'conorholds', 'GITHUB_RUN_ID': '123',
            'GITHUB_RUN_ATTEMPT': '1', 'GITHUB_SHA': MAIN_SHA}
        self.pull = {'state': 'open', 'base': {
            'ref': 'main', 'repo': {'full_name': ci.REPOSITORY}}, 'head': {'sha': PR_SHA}}
        self.check = {'id': 7, 'name': ci.CHECK_NAME, 'head_sha': PR_SHA,
                      'external_id': 'manual-verification:123'}

    def test_main_is_pinned_to_the_dispatched_commit(self):
        api = Mock(return_value={'id': 7})
        output = ci.start('', self.env, api)
        self.assertEqual(output, {'sha': MAIN_SHA, 'check_id': '7'})
        api.assert_called_once()
        method, path, payload = api.call_args.args
        self.assertEqual((method, path), ('POST', 'check-runs'))
        self.assertEqual(payload['head_sha'], MAIN_SHA)
        self.assertEqual(payload['status'], 'in_progress')
        self.assertEqual(payload['name'], 'Verify certificates')

    def test_pr_check_uses_the_resolved_head_not_main(self):
        api = Mock(side_effect=[self.pull, {'id': 7}])
        output = ci.start('8', self.env, api)
        self.assertEqual(output['sha'], PR_SHA)
        self.assertEqual(api.call_args_list[0].args, ('GET', 'pulls/8'))
        payload = api.call_args_list[1].args[2]
        self.assertEqual(payload['head_sha'], PR_SHA)
        self.assertEqual(payload['external_id'], self.check['external_id'])

    def test_other_actors_events_and_branches_cannot_start_or_report(self):
        changes = {'GITHUB_REPOSITORY': 'someone/fork', 'GITHUB_EVENT_NAME': 'push',
                   'GITHUB_REF': 'refs/heads/contribution', 'GITHUB_ACTOR_ID': '999',
                   'GITHUB_TRIGGERING_ACTOR': 'another-maintainer',
                   'GITHUB_RUN_ID': '123\nsha=wrong', 'GITHUB_RUN_ATTEMPT': '0'}
        for field, value in changes.items():
            with self.subTest(field=field):
                env = dict(self.env, **{field: value})
                api = Mock()
                with self.assertRaises(ValueError):
                    ci.start('8', env, api)
                with self.assertRaises(ValueError):
                    ci.finish('7', PR_SHA, 'success', env, api)
                api.assert_not_called()

    def test_invalid_pr_inputs_never_reach_the_api(self):
        for value in ['0', '-1', '8; echo bad', '$(pwd)', '../8', '8\nsha=wrong', ' 8']:
            with self.subTest(value=value):
                api = Mock()
                with self.assertRaises(ValueError):
                    ci.start(value, self.env, api)
                api.assert_not_called()

    def test_closed_foreign_and_malformed_prs_cannot_create_checks(self):
        variants = []
        for field, value in [('state', 'closed'), ('ref', 'other'),
                             ('full_name', 'someone/other'), ('sha', 'not-a-sha')]:
            pull = copy.deepcopy(self.pull)
            if field == 'state':
                pull[field] = value
            elif field == 'ref':
                pull['base'][field] = value
            elif field == 'full_name':
                pull['base']['repo'][field] = value
            else:
                pull['head'][field] = value
            variants.append(pull)
        for pull in variants:
            with self.subTest(pull=pull):
                api = Mock(return_value=pull)
                with self.assertRaises(ValueError):
                    ci.start('8', self.env, api)
                api.assert_called_once_with('GET', 'pulls/8')
        api = Mock(side_effect=subprocess.CalledProcessError(1, ['gh', 'api']))
        with self.assertRaises(subprocess.CalledProcessError):
            ci.start('8', self.env, api)
        api.assert_called_once()

    def test_only_a_successful_verification_reports_success(self):
        for result, expected in [('success', 'success'), ('failure', 'failure'),
                                 ('cancelled', 'cancelled'), ('skipped', 'failure')]:
            with self.subTest(result=result):
                api = Mock(side_effect=[self.check, {}])
                summary = ci.finish('7', PR_SHA, result, self.env, api)
                method, path, payload = api.call_args.args
                self.assertEqual((method, path), ('PATCH', 'check-runs/7'))
                self.assertEqual(payload['conclusion'], expected)
                self.assertEqual(payload['status'], 'completed')
                self.assertIn(PR_SHA, summary)
                self.assertRegex(payload['completed_at'], r'^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$')
        api = Mock()
        with self.assertRaises(ValueError):
            ci.finish('7', PR_SHA, 'neutral', self.env, api)
        api.assert_not_called()

    def test_wrong_commit_check_or_run_cannot_be_completed(self):
        for field, value in [('head_sha', MAIN_SHA), ('name', 'unrelated'),
                             ('external_id', 'manual-verification:456')]:
            with self.subTest(field=field):
                api = Mock(return_value=dict(self.check, **{field: value}))
                with self.assertRaisesRegex(ValueError, 'does not belong'):
                    ci.finish('7', PR_SHA, 'success', self.env, api)
                api.assert_called_once_with('GET', 'check-runs/7')
        for check_id, sha in [('7/../8', PR_SHA), ('7', 'refs/heads/main')]:
            api = Mock()
            with self.assertRaises(ValueError):
                ci.finish(check_id, sha, 'success', self.env, api)
            api.assert_not_called()

    def test_owner_can_rerun_failed_jobs_without_re_resolving_the_commit(self):
        env = dict(self.env, GITHUB_RUN_ATTEMPT='2')
        api = Mock(side_effect=[self.check, {}])
        ci.finish('7', PR_SHA, 'success', env, api)
        self.assertTrue(api.call_args.args[2]['details_url'].endswith('/attempts/2'))

    def test_api_uses_structured_json_and_propagates_errors(self):
        payload = {'output': {'summary': 'literal $() and `text`\nsecond line'}}
        with patch.object(ci.subprocess, 'run', return_value=Mock(stdout='{"id":7}')) as run:
            self.assertEqual(ci.github_api('POST', 'check-runs', payload), {'id': 7})
            self.assertEqual(json.loads(run.call_args.kwargs['input']), payload)
            self.assertTrue(run.call_args.kwargs['check'])
            self.assertEqual(run.call_args.kwargs['timeout'], 30)
            self.assertNotIn('shell', run.call_args.kwargs)
        with patch.object(ci.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, ['gh'])):
            with self.assertRaises(subprocess.CalledProcessError):
                ci.github_api('GET', 'pulls/8')

    def test_cli_records_the_pinned_commit_and_reports_api_failures(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)/'outputs'
            summary = Path(folder)/'summary'
            env = dict(self.env, GITHUB_OUTPUT=str(output), GITHUB_STEP_SUMMARY=str(summary))
            responses = [Mock(stdout=json.dumps(self.pull)), Mock(stdout='{"id":7}'),
                         Mock(stdout=json.dumps(self.check)), Mock(stdout='{}')]
            with patch.dict(ci.os.environ, env, clear=True), \
                    patch.object(ci.subprocess, 'run', side_effect=responses), \
                    redirect_stdout(io.StringIO()):
                with patch.object(sys, 'argv', ['manual_ci.py', 'start', '--pr', '8']):
                    ci.main()
                with patch.object(sys, 'argv', ['manual_ci.py', 'finish', '--check-id', '7',
                                               '--sha', PR_SHA, '--result', 'success']):
                    ci.main()
            self.assertEqual(output.read_text(), f'sha={PR_SHA}\ncheck_id=7\n')
            self.assertIn(f'Commit: `{PR_SHA}`', summary.read_text())
            self.assertIn('**success**', summary.read_text())
            before = output.read_text()
            error = subprocess.CalledProcessError(1, ['gh'], stderr='HTTP 403: access denied')
            stderr = io.StringIO()
            with patch.dict(ci.os.environ, env, clear=True), \
                    patch.object(ci.subprocess, 'run', side_effect=error), \
                    patch.object(sys, 'argv', ['manual_ci.py', 'start', '--pr', '8']), \
                    redirect_stderr(stderr), self.assertRaises(SystemExit) as raised:
                ci.main()
            self.assertEqual(raised.exception.code, 1)
            self.assertIn('HTTP 403: access denied', stderr.getvalue())
            self.assertEqual(output.read_text(), before)


if __name__ == '__main__':
    unittest.main()
