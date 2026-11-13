// This runs only from trusted development pushes. PR metadata is data, never code.
module.exports = async function cleanup({ github, context, core,
  sleep = ms => new Promise(resolve => setTimeout(resolve, ms)),
}) {
  const { owner, repo } = context.repo;
  const repository = `${owner}/${repo}`;
  const allowedPaths = new Set([
    '.github/workflows/ci.yml',
    '.github/workflows/cuda-ci.yml',
    '.github/workflows/docs.yml',
    '.github/workflows/gcc-rrtmgp.yml',
    '.github/workflows/gcc.yml',
    '.github/workflows/hip.yml',
    '.github/workflows/lint.yml',
    '.github/workflows/macos.yml',
    '.github/workflows/parent-project.yml',
    '.github/workflows/style.yml',
    '.github/workflows/sycl.yml',
    '.github/workflows/windows.yml',
  ]);
  const { data: sweep } = await github.rest.actions.getWorkflowRun({
    owner, repo, run_id: context.runId,
  });
  const cutoff = Date.parse(sweep.created_at);
  if (!Number.isFinite(cutoff)) throw new Error('Missing cleanup run creation time');

  const runs = new Map();
  for (const status of ['queued', 'pending', 'in_progress', 'waiting', 'requested']) {
    const listed = await github.paginate(github.rest.actions.listWorkflowRunsForRepo, {
      owner, repo, event: 'pull_request', status, per_page: 100,
    });
    for (const run of listed) {
      if (run.event === 'pull_request' && run.status !== 'completed'
          && allowedPaths.has(run.path) && Date.parse(run.created_at) <= cutoff
          && run.head_repository?.id && run.head_repository.owner?.login) {
        runs.set(run.id, run);
      }
    }
  }

  async function closedPRsOnly(run) {
    let associations = run.pull_requests || [];
    if (associations.length === 0) {
      // Fork runs often have an empty pull_requests array. Resolve by commit,
      // not branch name alone: different forks can use the same branch name.
      associations = await github.paginate(github.rest.repos.listPullRequestsAssociatedWithCommit, {
        owner, repo, commit_sha: run.head_sha, per_page: 100,
      });
    }
    let matched = false;
    for (const number of new Set(associations.map(pr => pr.number))) {
      // Re-read state on each sweep/verification; cached closed state can go stale.
      const { data: pr } = await github.rest.pulls.get({ owner, repo, pull_number: number });
      if (pr.base.repo.full_name !== repository || pr.head.repo?.id !== run.head_repository.id
          || pr.head.ref !== run.head_branch) continue;
      if (pr.state !== 'closed') return false;
      matched = true;
    }
    if (!matched) return false;
    // A closed PR's source branch may have been reused for a different open PR.
    const open = await github.paginate(github.rest.pulls.list, {
      owner, repo, state: 'open',
      head: `${run.head_repository.owner.login}:${run.head_branch}`, per_page: 100,
    });
    return !open.some(pr => pr.head.repo?.id === run.head_repository.id
      && pr.head.ref === run.head_branch);
  }

  async function unchangedAttempt(run) {
    const { data: current } = await github.rest.actions.getWorkflowRun({ owner, repo, run_id: run.id });
    return current.status !== 'completed'
      && (current.run_attempt || 1) === (run.run_attempt || 1);
  }

  let failures = 0;
  const cancelled = [];
  for (const run of runs.values()) {
    try {
      if (!await closedPRsOnly(run) || !await unchangedAttempt(run)) continue;
      await github.rest.actions.cancelWorkflowRun({ owner, repo, run_id: run.id });
      cancelled.push(run);
      core.info(`Requested cancellation: ${run.id} (${run.path})`);
    } catch (error) {
      if (error.status === 404 || error.status === 409) {
        core.info(`Skipped unavailable or finished run/PR: ${run.id}`);
      } else {
        core.warning(`Could not clean up run ${run.id}: ${error.message}`);
        failures++;
      }
    }
  }

  if (cancelled.length > 0) {
    // A normal cancellation can leave deeply queued runs stuck. Allow graceful
    // cleanup first, then force-cancel only if the PR is still closed and the
    // run has not been restarted as a newer attempt.
    await sleep(60000);
    for (const run of cancelled) {
      try {
        if (!await unchangedAttempt(run) || !await closedPRsOnly(run)) continue;
        await github.request('POST /repos/{owner}/{repo}/actions/runs/{run_id}/force-cancel', {
          owner, repo, run_id: run.id,
        });
        core.info(`Force-cancelled: ${run.id} (${run.path})`);
      } catch (error) {
        if (error.status === 404 || error.status === 409) {
          core.info(`Run/PR ${run.id} became unavailable or finished`);
        } else {
          core.warning(`Could not force-cancel run ${run.id}: ${error.message}`);
          failures++;
        }
      }
    }
  }
  core.info(`Examined ${runs.size} eligible runs; cancellation requests: ${cancelled.length}; errors: ${failures}`);
  if (failures > 0) core.setFailed(`${failures} cleanup operation(s) failed`);
};
