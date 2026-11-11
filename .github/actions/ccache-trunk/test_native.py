#!/usr/bin/env python3
"""Native transport failures, archive safety, and portable pruning; no network."""

import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import cache
import oras
import transport
import prune


class NativeTests(unittest.TestCase):
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
        self.source = self.root / 'source'
        (self.source / 'a').mkdir(parents=True)
        (self.source / 'a' / 'object').write_bytes(b'new object')
        os.utime(self.source / 'a' / 'object', (1_600_000_000, 1_600_000_000))
        self.archive = self.root / 'remote.tar'
        oras.pack(self.source, self.archive)
        self.manifest = {'artifactType': oras.ARTIFACT_TYPE, 'layers': [{
            'mediaType': oras.LAYER_TYPE, 'digest': 'sha256:' + oras.sha256(self.archive),
            'size': self.archive.stat().st_size}]}
        self.env = dict(CACHE_MODE='restore', CACHE_KEY='macOS-X64-Clang',
                        CACHE_PATH=str(self.cache), GHCR_TOKEN='synthetic-token',
                        RUNNER_OS='macOS', RUNNER_ARCH='X64', RUNNER_TEMP=str(self.runner),
                        GITHUB_OUTPUT=str(self.output), GITHUB_ACTOR='tester',
                        GITHUB_REPOSITORY_OWNER='Example', GITHUB_SHA='0' * 40,
                        CACHE_SOURCE='https://example.invalid/repo', CACHE_OWNER_TYPE='User')
        self.calls = []

    def command(self, args, **kwargs):
        self.calls.append((args, kwargs))
        if args[1] == 'login':
            self.assertEqual(kwargs['input'], 'synthetic-token')
            self.assertNotIn('synthetic-token', args)
            Path(args[args.index('--registry-config') + 1]).write_text('synthetic credentials')
        elif args[1:3] == ['manifest', 'fetch']:
            return json.dumps(self.manifest)
        elif args[1:3] == ['blob', 'fetch']:
            self.assertTrue(args[-1].endswith('@' + self.manifest['layers'][0]['digest']))
            shutil.copyfile(self.archive, args[args.index('--output') + 1])
        return ''

    def action(self, mode='restore', command=None, missing=(), **env):
        with patch.dict(os.environ, {**self.env, "CACHE_MODE": mode, **env}), \
                patch.object(oras, 'install_oras', return_value=Path('oras')), \
                patch.object(transport, 'run', side_effect=command or self.command), \
                contextlib.redirect_stderr(io.StringIO()) as err:
            for name in missing:
                os.environ.pop(name, None)
            status = cache.main()
        self.assertEqual(list(self.runner.iterdir()), [], 'temporary credentials leaked')
        self.assertFalse(list(self.root.glob('*.ghcr-*')), 'staging directory leaked')
        return status, err.getvalue()

    def preserved(self):
        self.assertEqual(list(self.cache.iterdir()), [self.cache / 'old'])
        self.assertEqual((self.cache / 'old').read_text(), 'original')
        self.assertEqual(self.output.read_text(), 'hit=false\n')

    def test_round_trip_and_anonymous_restore(self):
        status, errors = self.action(GHCR_TOKEN='')
        self.assertEqual((status, errors), (0, ''))
        self.assertEqual((self.cache / 'a/object').read_bytes(), b'new object')
        self.assertFalse((self.cache / 'old').exists())
        self.assertEqual((self.cache / 'a/object').stat().st_mtime, 1_600_000_000)
        self.assertEqual(self.output.read_text(), 'hit=false\nhit=true\n')
        self.assertFalse(any(args[1] == 'login' for args, _ in self.calls))

    def test_failed_transfer_preserves_cache(self):
        for failed in ('login', 'manifest', 'blob'):
            with self.subTest(failed=failed):
                self.output.unlink(missing_ok=True)
                def command(args, **kwargs):
                    if args[1] == failed:
                        raise RuntimeError('simulated timeout/unavailable registry')
                    return self.command(args, **kwargs)
                status, errors = self.action(command=command)
                self.assertEqual(status, 0)
                self.assertIn('::warning::', errors)
                self.preserved()

    def test_anonymous_restore_needs_no_publishing_metadata(self):
        status, errors = self.action(GHCR_TOKEN='', missing=(
            'GITHUB_ACTOR', 'CACHE_SOURCE', 'GITHUB_SHA', 'CACHE_OWNER_TYPE'))
        self.assertEqual((status, errors), (0, ''))
        self.assertEqual(self.output.read_text(), 'hit=false\nhit=true\n')
        self.assertEqual((self.cache / 'a/object').read_bytes(), b'new object')

    def test_missing_required_environment_fails_before_registry_access(self):
        cases = [('restore', name) for name in
                 ('GITHUB_REPOSITORY_OWNER', 'GITHUB_OUTPUT', 'RUNNER_TEMP', 'GITHUB_ACTOR')]
        cases += [('save', name) for name in ('CACHE_SOURCE', 'GITHUB_SHA', 'CACHE_OWNER_TYPE')]
        for mode, name in cases:
            with self.subTest(mode=mode, name=name):
                status, errors = self.action(mode, missing=(name,))
                self.assertEqual(status, 2)
                self.assertIn('::error::', errors)
                self.assertIn(name, errors)
                self.assertEqual(self.calls, [])

    def test_configuration_stays_consistent_throughout_transfer(self):
        alternate_output = self.root / 'other-output'

        def command(args, **kwargs):
            result = self.command(args, **kwargs)
            if args[1] == 'login':
                os.environ.update(GITHUB_OUTPUT=str(alternate_output),
                                  CACHE_SOURCE='https://example.invalid/changed',
                                  GITHUB_SHA='1' * 40, GITHUB_REPOSITORY_OWNER='Changed',
                                  CACHE_OWNER_TYPE='Organization', GHCR_TOKEN='changed-token')
            return result

        for mode in ('restore', 'save'):
            with self.subTest(mode=mode):
                self.calls.clear()
                self.output.unlink(missing_ok=True)
                status, errors = self.action(mode, command=command)
                self.assertEqual((status, errors), (0, ''))
                self.assertFalse(alternate_output.exists())
                if mode == 'restore':
                    self.assertEqual(self.output.read_text(), 'hit=false\nhit=true\n')
                else:
                    push = next(args for args, _ in self.calls if args[1] == 'push')
                    self.assertIn('org.opencontainers.image.source=' + self.env['CACHE_SOURCE'], push)
                    self.assertIn('org.opencontainers.image.revision=' + self.env['GITHUB_SHA'], push)
                    env = self.calls[-1][1]['env']
                    self.assertEqual(env['OWNER'], 'Example')
                    self.assertEqual(env['OWNER_TYPE'], 'User')
                    self.assertEqual(env['GH_TOKEN'], 'synthetic-token')

    def test_config_repr_does_not_expose_credentials(self):
        config = cache.read_config({**self.env, 'CACHE_MODE': 'save'})
        self.assertNotIn('synthetic-token', repr(config))

    def test_corrupt_blob_preserves_cache(self):
        self.archive.write_bytes(b'truncated download')
        self.assertEqual(self.action()[0], 0)
        self.preserved()

    def test_manifest_misses_are_quiet_but_other_failures_warn(self):
        ref = 'ghcr.io/example/erf-ccache:macos-x64-clang'
        cases = [
            ('manifest unknown', True),
            ('MANIFEST_UNKNOWN: manifest does not exist', True),
            (f'Error: {ref}: not found', True),
            ('registry: 403 Forbidden: permission denied', False),
            ('registry: 401 Unauthorized', False),
            ('registry: 503 Service Unavailable', False),
            ('Get https://ghcr.io/token: 404 Not Found', False),
            ('manifest unknown: 403 Forbidden', False),
            ('503 Service Unavailable: manifest unknown', False),
            ('Get ghcr.io/example/cache:manifest_unknown: 503 Service Unavailable', False),
        ]
        for detail, missing in cases:
            with self.subTest(detail=detail):
                self.output.unlink(missing_ok=True)
                self.calls.clear()

                def response(args, **kwargs):
                    if args[1:3] == ['manifest', 'fetch']:
                        raise subprocess.CalledProcessError(1, args, stderr=detail)
                    return subprocess.CompletedProcess(args, 0, self.command(args, **kwargs), '')

                with patch.object(transport.subprocess, 'run', side_effect=response), \
                        patch.object(transport.time, 'sleep'), \
                        contextlib.redirect_stdout(io.StringIO()) as output:
                    status, errors = self.action(command=transport.run)
                self.assertEqual(status, 0)
                self.preserved()
                self.assertFalse(any(args[1] == 'blob' for args, _ in self.calls))
                if missing:
                    self.assertEqual(errors, '')
                    self.assertIn('no cache at ' + ref, output.getvalue())
                else:
                    self.assertIn('::warning::', errors)
                    self.assertIn(detail, errors)

    def test_rename_failure_rolls_back(self):
        rename = Path.rename
        def fail_promotion(path, target):
            if '.ghcr-' in path.name and not path.name.endswith('.previous'):
                raise OSError('simulated rename failure')
            return rename(path, target)
        with patch.object(Path, 'rename', fail_promotion):
            self.assertEqual(self.action()[0], 0)
        self.preserved()

    def test_missing_blob_is_a_warning_and_preserves_cache(self):
        def response(args, **kwargs):
            if args[1:3] == ['blob', 'fetch']:
                raise subprocess.CalledProcessError(1, args, stderr='404 Not Found: blob unknown')
            return subprocess.CompletedProcess(args, 0, self.command(args, **kwargs), '')

        with patch.object(transport.subprocess, 'run', side_effect=response), \
                patch.object(transport.time, 'sleep'):
            status, errors = self.action(command=transport.run)
        self.assertEqual(status, 0)
        self.assertIn('fetch cache archive', errors)
        self.assertIn('404 Not Found: blob unknown', errors)
        self.preserved()

    def test_save_failure_does_not_prune(self):
        def fail_push(args, **kwargs):
            if args[1] == 'push':
                raise RuntimeError('simulated push failure')
            return self.command(args, **kwargs)
        self.assertEqual(self.action('save', command=fail_push)[0], 0)
        self.assertFalse(any('prune.py' in args[1] for args, _ in self.calls))

    def test_save_then_prune_and_prune_failure_is_nonfatal(self):
        def fail_prune(args, **kwargs):
            self.command(args, **kwargs)
            if args[1].endswith('prune.py'):
                raise RuntimeError('simulated prune timeout')
            return ''
        status, errors = self.action('save', command=fail_prune)
        self.assertEqual(status, 0)
        self.assertIn('a later save will retry', errors)
        self.assertIn('simulated prune timeout', errors)
        self.assertEqual(self.calls[-2][0][1], 'push')
        self.assertTrue(self.calls[-1][0][1].endswith('prune.py'))

    def test_prune_command_error_reaches_warning_without_credentials(self):
        def response(args, **kwargs):
            output = self.command(args, **kwargs)
            if args[1].endswith('prune.py'):
                raise subprocess.CalledProcessError(1, args, stderr=(
                    'erf-ccache: could not delete version 42: HTTP 403; synthetic-token'))
            return subprocess.CompletedProcess(args, 0, output, '')

        with patch.object(transport.subprocess, 'run', side_effect=response), \
                patch.object(transport.time, 'sleep'):
            status, errors = self.action('save', command=transport.run)
        self.assertEqual(status, 0)
        self.assertIn('exit status 1', errors)
        self.assertIn('could not delete version 42: HTTP 403', errors)
        self.assertIn('a later save will retry', errors)
        self.assertNotIn('synthetic-token', errors)

    def test_invalid_inputs_fail_and_empty_save_skips_transport(self):
        for env in ({'CACHE_MODE': 'oops'}, {'CACHE_KEY': ''}, {'CACHE_PATH': ''},
                    {'CACHE_PATH': str(Path.home())}, {'CACHE_PATH': self.root.anchor},
                    {'CACHE_TIMEOUT': '601'}, {'RUNNER_ARCH': 'unknown'},
                    {'CACHE_MODE': 'save', 'GHCR_TOKEN': ''}):
            self.assertEqual(self.action(**env)[0], 2, env)
        self.assertEqual(self.calls, [])
        (self.cache / 'old').unlink()
        self.assertEqual(self.action('save')[0], 0)
        self.assertEqual(self.calls, [])

    def test_bad_paths_links_and_special_files_rejected(self):
        cases = [(name, tarfile.REGTYPE) for name in
                 ('../escape', '/absolute', 'C:/drive', 'a\\..\\escape', 'a:stream', 'nul', 'a/../b')]
        cases += [('link', kind) for kind in (tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE)]
        for name, kind in cases:
            with self.subTest(name=name, kind=kind):
                with tarfile.open(self.archive, 'w') as bundle:
                    member = tarfile.TarInfo(name)
                    member.type = kind
                    member.linkname = '../escape' if kind in (tarfile.SYMTYPE, tarfile.LNKTYPE) else ''
                    bundle.addfile(member)
                self.manifest['layers'][0].update(digest='sha256:' + oras.sha256(self.archive),
                                                  size=self.archive.stat().st_size)
                self.output.unlink(missing_ok=True)
                self.assertEqual(self.action()[0], 0)
                self.preserved()
        self.assertFalse((self.root / 'escape').exists())

    def test_wrong_artifact_format_rejected(self):
        for manifest in ([], {}, {'artifactType': 'wrong', 'layers': []},
                         {'artifactType': oras.ARTIFACT_TYPE, 'layers': [None]}):
            with self.assertRaises(ValueError):
                oras.layer_descriptor(manifest)

    def test_oras_checksum_verified_before_extracting(self):
        def wrong_download(args, **kwargs):
            Path(args[args.index('--output') + 1]).write_bytes(b'wrong executable')
        with patch.object(transport, 'run', side_effect=wrong_download):
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                oras.install_oras(self.runner, ('macOS', 'ARM64'))
        self.assertFalse((self.runner / 'oras').exists())

    def test_timeout_and_transient_command_retry(self):
        timeout = subprocess.TimeoutExpired('oras', 1)
        good = subprocess.CompletedProcess(['oras'], 0, 'complete', '')
        with patch.object(subprocess, 'run', side_effect=[timeout, good]) as command, \
                patch.object(transport.time, 'sleep'):
            self.assertEqual(transport.run(['oras', 'pull'], timeout=1), 'complete')
            self.assertEqual(command.call_count, 2)
        with patch.object(subprocess, 'run', side_effect=timeout) as command, \
                patch.object(transport.time, 'sleep'):
            with self.assertRaises(RuntimeError):
                transport.run(['oras', 'pull'], timeout=1)
            self.assertEqual(command.call_count, 2)

    def test_unreadable_archive_preserves_cache(self):
        with patch.object(tarfile.TarFile, 'extractfile', return_value=None):
            status, errors = self.action()
        self.assertEqual(status, 0)
        self.assertIn('cache archive member is not a file', errors)
        self.preserved()

    def test_oras_executable_must_be_a_file(self):
        filename = oras.ORAS_RELEASES[('macOS', 'ARM64')][0]
        archive = self.runner / f'oras_{oras.ORAS_VERSION}_{filename}'
        with tarfile.open(archive, 'w') as bundle:
            member = tarfile.TarInfo('oras')
            member.type = tarfile.DIRTYPE
            bundle.addfile(member)
        releases = {('macOS', 'ARM64'): (filename, oras.sha256(archive))}
        with patch.object(oras, 'ORAS_RELEASES', releases), patch.object(transport, 'run'):
            with self.assertRaisesRegex(ValueError, 'ORAS executable is not a file'):
                oras.install_oras(self.runner, ('macOS', 'ARM64'))
        self.assertFalse((self.runner / 'oras').exists())


class PruneTests(unittest.TestCase):
    VERSIONS = ('1\t2026-10-28T03:00:00Z\tactive-trunk\n'
                '2\t2026-10-28T02:00:00Z\t\n'
                '3\t2026-10-28T01:00:00Z\t\n')

    def test_failed_deletion_does_not_block_later_versions(self):
        deleted = []

        def github_response(command, **kwargs):
            if 'DELETE' not in command:
                return subprocess.CompletedProcess(command, 0, self.VERSIONS, '')
            version = command[-1].rsplit('/', 1)[1]
            deleted.append(version)
            if version == '2':
                return subprocess.CompletedProcess(command, 1, '', 'gh: Forbidden (HTTP 403)')
            return subprocess.CompletedProcess(command, 0, '', '')

        with patch.dict(os.environ, OWNER='Example', OWNER_TYPE='User'), \
                patch.object(prune.subprocess, 'run', side_effect=github_response), \
                patch.object(prune.time, 'monotonic', return_value=0), \
                patch.object(prune.time, 'sleep'), \
                contextlib.redirect_stderr(io.StringIO()) as errors, \
                contextlib.redirect_stdout(io.StringIO()) as output:
            status = prune.main(['erf-ccache:0'])

        self.assertEqual(deleted, ['2', '2', '2', '3'])
        self.assertEqual(status, 1, 'partial cleanup must still report failure')
        self.assertIn('could not delete version 2', errors.getvalue())
        self.assertIn('version 3 deleted', output.getvalue())

    def test_failed_deletion_stops_at_total_deadline(self):
        now = [0]
        deleted = []

        def github_response(command, **kwargs):
            if 'DELETE' not in command:
                return subprocess.CompletedProcess(command, 0, self.VERSIONS, '')
            deleted.append(command[-1].rsplit('/', 1)[1])
            now[0] = 61
            raise subprocess.TimeoutExpired(command, kwargs['timeout'])

        with patch.dict(os.environ, OWNER='Example', OWNER_TYPE='User'), \
                patch.object(prune.subprocess, 'run', side_effect=github_response), \
                patch.object(prune.time, 'monotonic', side_effect=lambda: now[0]), \
                patch.object(prune.time, 'sleep'), \
                contextlib.redirect_stderr(io.StringIO()) as errors:
            status = prune.main(['erf-ccache:0'])

        self.assertEqual(deleted, ['2'], 'do not retry or delete more versions after the deadline')
        self.assertEqual(status, 1)
        self.assertIn('pruning timed out', errors.getvalue())

    def test_total_deadline_stops_requests(self):
        with patch.object(prune.time, 'monotonic', return_value=61), \
                patch.object(prune.subprocess, 'run') as command:
            with self.assertRaisesRegex(RuntimeError, 'timed out'):
                prune.github_api(['--paginate', 'example'], deadline=60)
            command.assert_not_called()


if __name__ == '__main__':
    unittest.main()
