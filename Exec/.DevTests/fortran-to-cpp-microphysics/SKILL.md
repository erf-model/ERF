---
name: fortran-to-cpp-microphysics
description: ERF microphysics porting and validation playbook for WRF-style Fortran-to-C++ conversions, including bridge wiring, native AMReX GPU translation, WSM6/WDM6/Morrison/MYNN25 application notes, and rule-gated validation campaigns.
---

Use this skill when working on ERF microphysics ports, bridge-vs-native parity, rulebook extraction, or campaign-style validation decisions.

## Workflow
1. Read [`references/core-rules.md`](references/core-rules.md).
2. Read the target application note in `references/applications/`.
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

## Non-Negotiable Discipline
- Treat the tracked Fortran bridge path as oracle until native parity is proven for the active scope.
- Do not silently weaken rule language during extraction. Use the coverage matrix to justify any wording change.
- Treat TSV ledgers as source-of-truth validation infrastructure, not as optional notes.
- Use tracked assets only for formal validation outcomes.
- Prefer repo facts, branch history, and GitHub evidence over memory.

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
