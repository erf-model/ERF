---
name: fortran-to-cpp-microphysics
description: ERF microphysics porting and validation playbook for WRF-style Fortran-to-C++ conversions, including bridge wiring, native AMReX GPU translation, WSM6/WDM6/Morrison/MYNN25 application notes, and rule-gated validation campaigns.
---

Use this skill when working on ERF microphysics ports, bridge-vs-native parity, rulebook extraction, or campaign-style validation decisions.

## Workflow
1. Read [`references/core-rules.md`](references/core-rules.md).
2. Read the target application note in `references/applications/`, then its
   scheme subdirectory `references/applications/<scheme>/` for the group map,
   campaign definition, lane ledgers, and any handoff note.
3. Read [`references/validation-operators.md`](references/validation-operators.md) for parity, retreat, restart, milestone, or ledger work.
4. Read [`references/process-operators.md`](references/process-operators.md) for evidence capture, branch-history use, or debug discipline.
5. Read [`references/coverage-matrix.md`](references/coverage-matrix.md) before changing the skill package or when justifying a rule/operator.
6. Read [`references/source-materials.md`](references/source-materials.md) or [`references/source-ledger.md`](references/source-ledger.md) when tracing provenance.

## Routing
- `core-rules.md`: reusable cross-scheme porting rules.
- `validation-operators.md`: formal validation lanes, ledgers, and decision procedures.
- `process-operators.md`: campaign habits and operator discipline.
- `applications/wsm6.md`: mature WSM6 source corpus and extraction guidance.
- `applications/wdm6.md`: current WDM6 implementation note and divergence list.
- `applications/<scheme>/`: scheme-local campaign definitions, group maps, lane
  ledgers, and campaign runners. Milestone ladders and fcompare protocol are
  scheme-specific by rule and live here, not in the operator files.

## Non-Negotiable Discipline
- Treat the tracked Fortran bridge path as oracle until native parity is proven for the active scope.
- Do not silently weaken rule language during extraction. Use the coverage matrix to justify any wording change.
- Treat TSV ledgers as source-of-truth validation infrastructure, not as optional notes.
- Use tracked assets only for formal validation outcomes.
- Prefer repo facts, branch history, and GitHub evidence over memory.
- **Do not reinvent. Grep the package and the tree first.** This is the most
  frequently violated rule in this skill and it has cost real time every time:
  an existing `orchestrator_mandatory` rule renamed as a new idea, a tracked
  campaign runner reimplemented as a scratch script, a `build_ci` feature the
  corpus explicitly forbids, and a Tier 2 print format invented while a
  conforming emit helper already sat in the same source file. Naming a rule is
  not applying it, and an empty-looking inventory in a note is not proof that
  the instrument is absent.
- **State the governing rule before using a method.** Before running prints,
  traces, or a retreat step, name the rule in `validation-operators.md`,
  `process-operators.md`, `core-rules.md`, or the source corpus that authorizes
  it. If no rule covers it, say so explicitly rather than proceeding on handoff
  precedent or a prior `decisions.tsv` row; neither is the rule package.
- **When a gap forces an invention, promote the missing rule up into THIS
  file**, not only into the reference where it belongs. Inventions happen
  because this compressed document is what gets read; a rule that exists only
  in a reference file did not prevent anything. Cite the evidence in
  `decisions.tsv` and, if the wording changes, in the coverage matrix.
- **A missing instrument is not a blocker.** If a scheme lacks a diagnostic
  tier, spot-add it narrowly at the site that needs it. Do not stand up a
  scheme-wide framework in advance, and do not treat its absence as permission
  to patch on weaker evidence.
- **The source is the oracle for structure and placement, not only values.**
  Before deciding where a computation lives, what gates it, or how a loop is
  split, quote the source construct: *where* it sits in the per-step sequence,
  *what* it actually gates on rather than a proxy, and *how* it is nested. One
  grep. Four WDM6 defects traced to skipping it, and each was confidently wrong
  rather than uncertain — so deliberation does not substitute for the check.
- **When one edit covers N sites, grep for the other N-1 afterwards.** An
  incomplete `ratio_s`/`ratio_g` fix survived a build, a Tier 2 verification and
  a milestone gate because the sibling was genuinely fixed and the comment
  claimed both. A comment asserting a fix is not evidence of one.

## Tool Expectations
- Use local repo inspection and git history by default.
- Use GitHub context for PR lineage, branch diffs, and review context.
- Use builds/tests when the task needs factual validation.
- Do not reach for external docs unless the task truly becomes a library/API question.

## Classification Rule
Every extracted source item must land in exactly one home:
- `core-rules`
- `validation-operators`
- `process-operators`
- `applications/<scheme>`
- `source-materials`
- `retired`

Use near-verbatim carry-forward for hard-won rules and operators unless there is a strong reason to normalize wording.
