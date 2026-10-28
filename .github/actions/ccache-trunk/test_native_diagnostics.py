#!/usr/bin/env python3
"""Command diagnostics use synthetic failures and credentials; no network."""

import contextlib
import io
import os
import subprocess
import unittest
from unittest.mock import patch

import transport


class CommandTests(unittest.TestCase):
    def test_invalid_retry_count_fails_before_running_command(self):
        for attempts in (0, -1):
            with self.subTest(attempts=attempts), patch.object(transport.subprocess, 'run') as command:
                with self.assertRaisesRegex(ValueError, 'attempts must be positive'):
                    transport.run(['oras', 'version'], attempts=attempts)
                command.assert_not_called()

    def failure(self, exception, **kwargs):
        with patch.object(transport.subprocess, 'run', side_effect=exception) as command, \
                patch.object(transport.time, 'sleep'):
            with self.assertRaises(RuntimeError) as caught:
                transport.run(['oras', 'manifest', 'fetch', 'private-command-argument'], **kwargs)
        self.assertEqual(command.call_count, kwargs.get('attempts', 2))
        message = str(caught.exception)
        self.assertNotIn('private-command-argument', message)
        return message

    def test_exit_status_and_stderr_are_reported(self):
        message = self.failure(subprocess.CalledProcessError(
            7, ['oras'], stderr='registry: 403 Forbidden: permission denied'))
        self.assertIn('exit status 7', message)
        self.assertIn('403 Forbidden: permission denied', message)

    def test_stdout_is_used_when_stderr_is_empty(self):
        message = self.failure(subprocess.CalledProcessError(
            1, ['oras'], output='registry: 503 Service Unavailable', stderr=''))
        self.assertIn('503 Service Unavailable', message)

    def test_timeout_includes_partial_byte_output(self):
        message = self.failure(subprocess.TimeoutExpired(
            ['oras'], 3, stderr=b'TLS handshake stalled\xff'), timeout=3)
        self.assertIn('timed out after 3s', message)
        self.assertIn('TLS handshake stalled', message)

    def test_missing_executable_is_identified(self):
        message = self.failure(FileNotFoundError(2, 'No such file or directory'))
        self.assertIn('FileNotFoundError', message)
        self.assertIn('No such file or directory', message)

    def test_credentials_and_multiline_output_are_sanitized(self):
        credentials = ['synthetic-ghcr-secret', 'synthetic-gh-secret', 'synthetic-github-secret',
                       'synthetic-stdin-secret', 'opaque-bearer', 'opaque-basic',
                       'opaque-password', 'opaque-auth', 'url-password', 'signed-query']
        detail = (' '.join(credentials[:4]) + '\nAuthorization: Bearer opaque-bearer\n'
                  'Authorization: Basic opaque-basic\npassword=opaque-password\n'
                  '"auth": "opaque-auth"\n'
                  'https://user:url-password@example.invalid/v2?signature=signed-query\n'
                  '::error::untrusted output')
        with patch.dict(os.environ, GHCR_TOKEN=credentials[0], GH_TOKEN=credentials[1],
                        GITHUB_TOKEN=credentials[2]):
            message = self.failure(subprocess.CalledProcessError(1, ['oras'], stderr=detail),
                                   input=credentials[3])
        for credential in credentials:
            self.assertNotIn(credential, message)
        self.assertIn('[redacted]', message)
        self.assertNotIn('\n', message)
        self.assertNotIn('\r', message)

    def test_redaction_precedes_output_limit(self):
        token = 'synthetic-secret-that-crosses-the-output-limit'
        with patch.dict(os.environ, GHCR_TOKEN=token):
            message = self.failure(subprocess.CalledProcessError(
                1, ['oras'], stderr='x' * 1000 + token + 'y' * 2000))
        self.assertNotIn(token[:12], message)
        self.assertIn('[redacted]', message)
        self.assertLess(len(message), 1200)
        self.assertIn('truncated', message)

    def test_child_environment_credentials_are_redacted(self):
        message = self.failure(subprocess.CalledProcessError(
            1, ['oras'], stderr='permission denied for synthetic-child-token'),
            env={'GH_TOKEN': 'synthetic-child-token'})
        self.assertIn('permission denied', message)
        self.assertNotIn('synthetic-child-token', message)
        self.assertIn('[redacted]', message)

    def test_warning_is_a_single_escaped_annotation(self):
        with patch.dict(os.environ, GHCR_TOKEN='synthetic-warning-token'), \
                contextlib.redirect_stderr(io.StringIO()) as output:
            transport.warning('failed: synthetic-warning-token\n::error::external%0Amessage')
        self.assertEqual(output.getvalue(),
                         '::warning::failed: [redacted] ::error::external%250Amessage\n')


if __name__ == '__main__':
    unittest.main()
