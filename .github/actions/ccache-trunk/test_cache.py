#!/usr/bin/env python3
"""The same observable cache contract, exercised through both real backends."""

import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import cache
import oras
import prune
import transport


class TransportContract:
    """Only external commands are mocked; policy and backend code both execute."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root / 'cache with spaces'
        self.cache.mkdir()
        (self.cache / 'old').write_text('original')
        self.runner = self.root / 'runner'
        self.runner.mkdir()
        self.output = self.root / 'output'
        self.source = self.root / 'snapshot'
        self.source.mkdir()
        (self.source / 'object').write_text('complete object')
        os.utime(self.source / 'object', (1_600_000_000, 1_600_000_000))
        self.archive = self.root / 'snapshot.tar'
        oras.pack(self.source, self.archive)
        self.manifest = {'artifactType': oras.ARTIFACT_TYPE, 'layers': [{
            'mediaType': oras.LAYER_TYPE, 'digest': 'sha256:' + oras.sha256(self.archive),
            'size': self.archive.stat().st_size}]}
        self.shared_config = self.root / 'shared-docker-config'
        self.shared_config.mkdir()
        (self.shared_config / 'config.json').write_text('existing credentials')
        self.env = dict(
            CACHE_MODE='restore', CACHE_TRANSPORT=self.backend,
            CACHE_KEY='Ubuntu-22.04 X64 GCC MPICH DOUBLE Particles ON',
            CACHE_PATH=str(self.cache), GHCR_TOKEN='synthetic-test-token',
            RUNNER_OS='Linux', RUNNER_ARCH='X64', RUNNER_TEMP=str(self.runner),
            GITHUB_OUTPUT=str(self.output), GITHUB_REPOSITORY_OWNER='ExampleOwner',
            GITHUB_ACTOR='tester', GITHUB_SHA='0' * 40,
            CACHE_SOURCE='https://example.invalid/repo', CACHE_OWNER_TYPE='User',
            DOCKER_CONFIG=str(self.shared_config))
        self.calls = []
        self.failure = None

    def command(self, args, **kwargs):
        operation = args[1]
        if operation.endswith('prune.py'):
            phase = 'prune'
        else:
            phase = {'pull': 'fetch', 'manifest': 'fetch', 'cp': 'extract',
                     'blob': 'extract', 'rm': 'close', 'build': 'package'}.get(operation, operation)
        self.calls.append(phase)
        if phase == 'login':
            self.assertEqual(kwargs['input'], 'synthetic-test-token')
            self.assertNotIn('synthetic-test-token', args)
            if self.backend == 'docker':
                location = Path(kwargs['env']['DOCKER_CONFIG']) / 'config.json'
                self.assertNotEqual(location.parent, self.shared_config)
            else:
                location = Path(args[args.index('--registry-config') + 1])
            location.write_text('temporary credentials')
        if self.failure:
            self.failure(phase, args, kwargs)
        stdout = ''
        if operation == 'manifest':
            stdout = json.dumps(self.manifest)
        elif operation == 'blob':
            shutil.copyfile(self.archive, args[args.index('--output') + 1])
        elif operation == 'cp':
            shutil.copytree(self.source, args[-1], dirs_exist_ok=True)
        elif phase == 'prune':
            self.assertEqual(kwargs['env']['GH_TOKEN'], 'synthetic-test-token')
        return subprocess.CompletedProcess(args, 0, stdout, '')

    def action(self, mode='restore', **env):
        self.calls.clear()
        self.output.unlink(missing_ok=True)
        with patch.dict(os.environ, {**self.env, 'CACHE_MODE': mode, **env}), \
                patch.object(oras, 'install_oras', return_value=Path('oras')), \
                patch.object(transport.subprocess, 'run', side_effect=self.command), \
                patch.object(transport.time, 'sleep'), \
                contextlib.redirect_stdout(io.StringIO()) as out, \
                contextlib.redirect_stderr(io.StringIO()) as err:
            try:
                status = cache.main()
            finally:
                self.assertEqual(list(self.runner.iterdir()), [], 'temporary credentials leaked')
                self.assertEqual(list(self.root.glob('*.ghcr-*')), [], 'staging/backup leaked')
                self.assertEqual((self.shared_config / 'config.json').read_text(), 'existing credentials')
        return status, out.getvalue(), err.getvalue()

    def preserved(self):
        self.assertEqual(list(self.cache.iterdir()), [self.cache / 'old'])
        self.assertEqual((self.cache / 'old').read_text(), 'original')
        self.assertEqual(self.output.read_text(), 'hit=false\n')

    def test_complete_restore_replaces_cache(self):
        status, out, err = self.action()
        self.assertEqual((status, err), (0, ''))
        self.assertEqual((self.cache / 'object').read_text(), 'complete object')
        self.assertFalse((self.cache / 'old').exists())
        self.assertEqual((self.cache / 'object').stat().st_mtime, 1_600_000_000)
        self.assertEqual(self.output.read_text(), 'hit=false\nhit=true\n')
        self.assertIn('ubuntu-22.04-x64-gcc-mpich-double-particles-on', out)
        if self.backend == 'docker':
            self.assertEqual(self.calls[-1], 'close')

    def test_failed_restore_preserves_cache(self):
        phases = ['login', 'fetch', 'extract'] + (['create'] if self.backend == 'docker' else [])
        for failed in phases:
            with self.subTest(phase=failed):
                def fail(phase, args, kwargs):
                    if phase == failed:
                        raise subprocess.CalledProcessError(1, args, stderr='503 registry unavailable')
                self.failure = fail
                status, _, err = self.action()
                self.assertEqual(status, 0)
                self.assertIn('503 registry unavailable', err)
                self.preserved()
                if failed in ('login', 'fetch'):
                    self.assertEqual(self.calls.count(failed), 2)

    def test_partial_extraction_preserves_cache(self):
        if self.backend == 'docker':
            def partial(phase, args, kwargs):
                if phase == 'extract':
                    (Path(args[-1]) / 'partial-object').write_text('incomplete')
                    raise subprocess.CalledProcessError(1, args, stderr='partial copy')
            self.failure = partial
        else:
            with self.archive.open('ab') as stream:
                stream.write(b'corrupt download')
        status, _, err = self.action()
        self.assertEqual(status, 0)
        self.assertIn('::warning::', err)
        self.preserved()

    def test_misses_are_quiet_but_auth_and_server_failures_warn(self):
        for detail, quiet in [('manifest unknown', True), ('no such manifest', True),
                              ('manifest not found', True), ('403 Forbidden', False),
                              ('503 Service Unavailable: manifest unknown', False)]:
            with self.subTest(detail=detail):
                def fail(phase, args, kwargs):
                    if phase == 'fetch':
                        raise subprocess.CalledProcessError(1, args, stderr=detail)
                self.failure = fail
                status, out, err = self.action()
                self.assertEqual(status, 0)
                self.preserved()
                if quiet:
                    self.assertEqual(err, '')
                    self.assertIn('no cache at', out)
                else:
                    self.assertIn(detail, err)
                self.assertNotIn('extract', self.calls)

    def test_transient_timeout_recovers(self):
        def fail_once(phase, args, kwargs):
            if phase == 'fetch' and self.calls.count('fetch') == 1:
                raise subprocess.TimeoutExpired(args, kwargs['timeout'])
        self.failure = fail_once
        self.assertEqual(self.action()[0], 0)
        self.assertEqual(self.calls.count('fetch'), 2)
        self.assertTrue(self.output.read_text().endswith('hit=true\n'))

    def test_persistent_timeout_preserves_cache(self):
        def timeout(phase, args, kwargs):
            if phase == 'fetch':
                raise subprocess.TimeoutExpired(args, kwargs['timeout'])
        self.failure = timeout
        status, _, err = self.action(CACHE_TIMEOUT='7')
        self.assertEqual(status, 0)
        self.assertIn('timed out after 7s', err)
        self.assertEqual(self.calls.count('fetch'), 2)
        self.preserved()

    def test_anonymous_restore_with_no_existing_directory(self):
        shutil.rmtree(self.cache)
        status, _, err = self.action(GHCR_TOKEN='')
        self.assertEqual((status, err), (0, ''))
        self.assertNotIn('login', self.calls)
        self.assertEqual((self.cache / 'object').read_text(), 'complete object')

    def test_save_then_prune_and_failure_is_nonfatal(self):
        def fail_prune(phase, args, kwargs):
            if phase == 'prune':
                raise subprocess.CalledProcessError(1, args, stderr='403 pruning denied synthetic-test-token')
        self.failure = fail_prune
        status, out, err = self.action('save')
        self.assertEqual(status, 0)
        self.assertEqual(self.calls[-2:], ['push', 'prune'])
        self.assertIn('saved ghcr.io/', out)
        self.assertIn('403 pruning denied', err)
        self.assertNotIn('synthetic-test-token', err)
        self.assertIn('a later save will retry', err)

    def test_transient_push_retries_before_pruning(self):
        def fail_once(phase, args, kwargs):
            if phase == 'push' and self.calls.count('push') == 1:
                raise subprocess.CalledProcessError(1, args, stderr='503 unavailable')
        self.failure = fail_once
        status, out, err = self.action('save')
        self.assertEqual((status, err), (0, ''))
        self.assertEqual(self.calls[-3:], ['push', 'push', 'prune'])
        self.assertIn('saved ghcr.io/', out)

    def test_failed_save_never_prunes(self):
        for failed in ['login', 'push'] + (['package'] if self.backend == 'docker' else []):
            with self.subTest(phase=failed):
                def fail(phase, args, kwargs):
                    if phase == failed:
                        raise subprocess.CalledProcessError(1, args, stderr='503 unavailable')
                self.failure = fail
                status, _, err = self.action('save')
                self.assertEqual(status, 0)
                self.assertIn('::warning::', err)
                self.assertNotIn('prune', self.calls)

    def test_empty_save_skips_transport(self):
        (self.cache / 'old').unlink()
        self.assertEqual(self.action('save')[0], 0)
        self.assertEqual(self.calls, [])

    def test_invalid_inputs_fail_before_transport(self):
        for env in ({'CACHE_MODE': 'invalid'}, {'CACHE_KEY': ''}, {'CACHE_PATH': ''},
                    {'CACHE_PATH': str(Path.home())}, {'CACHE_PATH': self.root.anchor},
                    {'CACHE_TIMEOUT': '601'}, {'CACHE_TIMEOUT': '01'},
                    {'CACHE_MODE': 'save', 'GHCR_TOKEN': ''}):
            with self.subTest(env=env):
                status, _, err = self.action(**env)
                self.assertEqual(status, 2)
                self.assertIn('::error::', err)
                self.assertEqual(self.calls, [])

    def test_failed_promotion_rolls_back(self):
        rename = Path.rename
        def fail(path, target):
            if '.ghcr-' in path.name and not path.name.endswith('.previous'):
                raise OSError('simulated rename failure')
            return rename(path, target)
        with patch.object(Path, 'rename', fail):
            status, _, err = self.action()
        self.assertEqual(status, 0)
        self.assertIn('simulated rename failure', err)
        self.preserved()

    def test_interrupted_promotion_rolls_back_and_propagates(self):
        rename = Path.rename
        for interruption in (SystemExit(143), KeyboardInterrupt()):
            with self.subTest(interruption=type(interruption).__name__):
                def interrupt(path, target):
                    if '.ghcr-' in path.name and not path.name.endswith('.previous'):
                        raise interruption
                    return rename(path, target)
                with patch.object(Path, 'rename', interrupt), self.assertRaises(type(interruption)):
                    self.action()
                self.preserved()

    def test_interruption_just_after_moving_old_cache_rolls_back(self):
        rename = Path.rename
        def interrupt_after_rename(path, target):
            result = rename(path, target)
            if target.name.endswith('.previous'):
                raise SystemExit(143)
            return result
        with patch.object(Path, 'rename', interrupt_after_rename), self.assertRaises(SystemExit):
            self.action()
        self.preserved()


class DockerContractTests(TransportContract, unittest.TestCase):
    backend = 'docker'


class OrasContractTests(TransportContract, unittest.TestCase):
    backend = 'oras'


class RetentionTests(unittest.TestCase):
    def test_named_tags_and_sha_retention(self):
        versions = ('1\t2026-10-27T01:00:00Z\tactive-trunk\n'
                    '2\t2026-10-27T00:00:00Z\t\n'
                    '3\t2026-10-27T02:00:00Z\tother-trunk\n'
                    '4\t2026-10-26T00:00:00Z\tsha-old\n'
                    '5\t2026-10-27T03:00:00Z\tsha-recent\n'
                    '6\t2026-10-27T04:00:00Z\tsha-newest\n'
                    '7\t2026-10-25T00:00:00Z\tsha-older,release\n')
        for keep, deleted in [(0, {'2', '4', '5', '6'}), (1, {'2', '4', '5'}), (2, {'2', '4'})]:
            with self.subTest(keep=keep):
                self.assertEqual(set(prune.select_versions(versions, keep)), deleted)

    def test_concurrent_delete_404_is_harmless(self):
        response = subprocess.CompletedProcess(['gh'], 1, '', 'gh: Not Found (HTTP 404)')
        with patch.object(prune.subprocess, 'run', return_value=response), \
                patch.object(prune.time, 'monotonic', return_value=0):
            self.assertIsNone(prune.github_api(['--method', 'DELETE', 'example'], deadline=60))


if __name__ == '__main__':
    unittest.main()
