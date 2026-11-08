"""Real compiler-cache and GHCR checks for the manual cache-smoke workflow."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class SmokeTest:
    def __init__(self):
        temporary = Path(os.environ['RUNNER_TEMP']).resolve()
        self.cache = temporary / 'ccache-smoke'
        # Empty only this workflow's dedicated directory, never a user's cache.
        require(Path(os.environ['CCACHE_DIR']).resolve() == self.cache,
                'CCACHE_DIR must be RUNNER_TEMP/ccache-smoke')
        self.work = temporary / 'compiler-smoke'
        self.work.mkdir(exist_ok=True)
        self.windows = os.environ['RUNNER_OS'] == 'Windows'
        self.compiler = {'Windows': 'cl', 'macOS': 'clang++', 'Linux': 'g++'}[os.environ['RUNNER_OS']]
        self.source = self.work / 'smoke.cpp'
        self.obj = self.work / ('smoke.obj' if self.windows else 'smoke.o')
        self.exe = self.work / ('smoke.exe' if self.windows else 'smoke')
        # A failed publication must not pass by restoring a previous run's program.
        run = os.environ.get('GITHUB_RUN_ID', 'local')
        attempt = os.environ.get('GITHUB_RUN_ATTEMPT', '1')
        self.original = f'{run}-{attempt}-original'
        self.changed = f'{run}-{attempt}-changed'

    def run(self, *args):
        try:
            return subprocess.check_output(args, cwd=self.work, text=True,
                                           stderr=subprocess.STDOUT, timeout=120)
        except subprocess.CalledProcessError as exc:
            print(exc.output, file=sys.stderr)
            raise

    def write_source(self, value):
        self.source.write_text('#include <stdio.h>\n'
                               f'int main() {{ puts({json.dumps(value)}); return 0; }}\n')

    def compile_and_run(self, value):
        if self.windows:
            self.run('ccache', self.compiler, '/nologo', '/c', str(self.source), '/Fo' + str(self.obj))
            self.run(self.compiler, '/nologo', str(self.obj), '/Fe' + str(self.exe))
        else:
            self.run('ccache', self.compiler, '-c', str(self.source), '-o', str(self.obj))
            self.run(self.compiler, str(self.obj), '-o', str(self.exe))
        actual = self.run(str(self.exe)).strip()
        require(actual == value, f'expected program output {value}, got {actual}')

    def stats(self):
        return {key: int(value) for key, value in
                (line.split() for line in self.run('ccache', '--print-stats').splitlines())}

    def require_hit(self):
        stats = self.stats()
        hits = stats['direct_cache_hit'] + stats['preprocessed_cache_hit']
        require(hits == 1 and stats['cache_miss'] == 0, f'expected one cache hit: {stats}')
        print(self.run('ccache', '-s'))

    def populate(self):
        self.write_source(self.original)
        self.run('ccache', '-z')
        self.compile_and_run(self.original)
        require(self.stats()['cache_miss'] == 1, 'initial compilation must miss the cache')
        print(self.run('ccache', '-s'))

    def empty(self):
        shutil.rmtree(self.cache)
        self.obj.unlink()
        self.exe.unlink()

    def verify(self):
        self.run('ccache', '-z')
        self.compile_and_run(self.original)
        self.require_hit()
        self.write_source(self.changed)
        self.compile_and_run(self.changed)
        require(self.stats()['cache_miss'] == 1, 'changed source must miss the cache')
        print('PASS: restored cache hit; changed source rebuilt with the correct output')

    def verify_replacement(self):
        self.run('ccache', '-z')
        self.compile_and_run(self.changed)
        self.require_hit()
        print('PASS: replacement snapshot includes the changed-source compilation')

    def versions(self):
        route = {'User': 'users', 'Organization': 'orgs'}[os.environ['CACHE_OWNER_TYPE']]
        owner = os.environ['GITHUB_REPOSITORY_OWNER']
        endpoint = f'{route}/{owner}/packages/container/erf-ccache/versions?per_page=100'
        pages = json.loads(self.run('gh', 'api', '--paginate', '--slurp', endpoint))
        return [version for page in pages for version in page]

    def tagged_version(self, versions):
        tag = os.environ['CCACHE_TRUNK_KEY']
        matches = [version['id'] for version in versions
                   if tag in version['metadata']['container']['tags']]
        require(len(matches) == 1, f'expected one published version for {tag}, got {matches}')
        return matches[0]

    def record_version(self):
        version = self.tagged_version(self.versions())
        (self.work / 'published-version').write_text(str(version))
        print(f'Published GHCR version {version}')

    def check_pruning(self):
        previous = int((self.work / 'published-version').read_text())
        # Allow brief API listing lag after the successful DELETE.
        for attempt in range(6):
            versions = self.versions()
            current = self.tagged_version(versions)
            if current != previous and all(version['id'] != previous for version in versions):
                print(f'PASS: replaced GHCR version {previous} with {current}; old version pruned')
                return
            if attempt < 5:
                time.sleep(2)
        raise RuntimeError(f'GHCR version {previous} was not replaced and pruned')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('populate', 'empty', 'verify', 'verify-replacement',
                                           'record-version', 'check-pruning'))
    args = parser.parse_args()
    getattr(SmokeTest(), args.command.replace('-', '_'))()
