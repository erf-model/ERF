# WDM6 Diagnostic Tier Inventory

## Provenance

Verbatim recovery of an analysis produced on 2026-08-13 that was never written to
disk at the time. Recovered 2026-08-14 from the Claude Code session transcript
`~/.claude/projects/-home-jmsexton-codes-ERF/8c9b33e6-c1af-40ec-840a-deafdd07a0d6.jsonl`,
assistant message timestamped `2026-08-13T21:25:37.670Z`, and committed here so
the package's provenance chain does not depend on a file outside the repository.

Reproduce with: parse the JSONL line by line, `json.loads` each, and read
`message.content` blocks of `type == "text"`.

## Status of its claims, as of 2026-08-14

Carried near-verbatim below. Three notes on what has moved since:

- The WDM6 counts predate the phase-2 all-k conversion, which converted 14 groups
  (G8, G10a-e, G11, G12, G13a-f). The `142 indexed at kts` figure is lower now.
  The structural finding is unchanged: P0, G1a-G7, G13g and G14-G17 are still
  single-cell, and G9 is deliberately excluded as surface-scoped.
- **Its central prediction was borne out.** It argued the Tier 1 conversion alone
  "would have caught the step-2 ice defect, since `qi` routes to
  G8/G13b/G13e/G13f and the failure is at k>=90". The step-2 frontier measured at
  clean SHA `bc03ab8d1` reports G8 POST diverging on `qi` at exactly k=90.
- The two generalizable rules extracted from this file live in
  `validation-operators.md` under the Diagnostic Contract: measure coverage per
  tag rather than tag count, and treat a tag reading zero where its inputs cannot
  activate the guarded block as a null result. See `coverage-matrix.md` rows
  `EVID_TIER_INVENTORY`, `RULE_TAG_COVERAGE` and `RULE_NULL_TAG`.

---

Data gathered. Here's the comparison.

## Tier structure across schemes

| | **Tier 1** canonical tag | **Tier 1.5** block signature | **Tier 2** forensic |
|---|---|---|---|
| **WSM6**<br>`wsm6-mpi-cpu` | 17 sites. `WSM6-FORT_<NAME>`, gate `mpdbg_level>=1 .and. loop==1`, `do k=kts,kte`, one column `its`, positional `(A,I3,6E24.16)` | 48 sites. `WSM6-FORT PRE/POST-<G>`, gate `>=2`, loops **k and i**, labelled `' qv=',q(i,k)`, 7 fields | 79 F90 + 104 C++ calls to `wsm6_emit_diag_t2_blocksig_line`. Gate `>=2`, **selected k only** (1–4, 29, 30), one var per line, full `DIAG-T2` schema |
| **Morrison** | none | none | none live. Historical only: `d131c5fc`/`1720371b` used `write(10,…)` to a file, 24-field payloads, hardcoded cell lists, `#if 0`'d. `357ab3b6` added `MPDBG`, `0e3749e6` removed it |
| **WDM6**<br>current | **none** | 153 F90 tags / 157 C++ — but **142 indexed at `kts`**, so 1 cell each. 0 `loop==1` guards | **none** |

The surprise: WDM6 isn't under-instrumented by *count* — 153 Fortran tags versus WSM6's 48. It's under-instrumented in *coverage per tag*. Every one samples a single cell, and there's no tier structure at all: one flat surface doing the job of three.

## What WDM6 tags look like at each scale

| | shape | per-call volume | what it buys |
|---|---|---|---|
| **fewer** | current, `kts` only | 153 lines | nothing at k≠`kts`. Rule 37 returned `0.000e+00` for all 33 groups while the plotfile diverged at k=90 |
| **current + `loop==1`** | `kts`, loop-guarded | 153 ÷ `wdm6_loops` | cuts the 320,291-marker flood, still blind in z |
| **→ Tier 1 (proposed)** | `do k=kts,kte`, one column, `loop==1` | ~15,300 | sees k=90. This is the WSM6 template |
| **more — Tier 1.5** | add `do i=its,ite`, gate `>=2` | ~3,060,000 | full i–k plane; only for a hard divergence |
| **more — Tier 2** | selected k, one var per line, `DIAG-T2` schema | ~900 | surgical line-level attribution, not breadth |

Note Tier 2 is **smaller** than Tier 1, not bigger. It's narrow-and-deep: WSM6 restricts to six levels (`k-kts+1` ∈ 1–4, 29, 30) and emits one variable per line with full provenance. The ladder isn't "more output" — it's *breadth* (T1) → *plane* (T1.5) → *depth* (T2).

## Concrete shapes for WDM6

```fortran
! Tier 1  — target: all-k, one column, loop-guarded
if (debug_local >= 1 .and. loop == 1) then
  do k = kts, kte
    write(*,'(A,1X,I3,14(1X,ES24.16E3))') 'WDM6-FORT_PRE_G7', k, &
      qrs(i_dbg_local,k,1), ... t(i_dbg_local,k), ...
  enddo
endif

! Tier 1.5 — only for a hard divergence
if (debug_local >= 2) then
  do k = kts, kte ; do i = its, ite
    write(*,'(A,3(A,I4),7(A,E22.15))') 'WDM6-FORT PRE-G7', &
      ' i=',i,' j=',j_dbg_local,' k=',k,' qv=',q(i,k), ...

! Tier 2 — selected k, one var per line, full schema
if (debug_local >= 2 .and. j_dbg_local >= 0) then
  do k = kts, kte
    if (.not. (k_sel(k))) cycle
    call wdm6_emit_diag_t2_blocksig_line("G7","BLOCK_SIGNATURE","INCORE_FORTRAN", &
         "G7","input",loop,i_dbg_local,j_dbg_local,k-kts+1,k_raw,debug_local,"t",t(i_dbg_local,k))
```

WDM6 keeps one advantage over WSM6 throughout: it parameterizes the column as `i_dbg_local` where WSM6 hardcodes `its`.

## Reading

Action B as scoped is only the **Tier 1** row — 142 sites, three mechanical edits each. That alone would have caught the step-2 ice defect, since `qi` routes to G8/G13b/G13e/G13f and the failure is at k≥90.

Tier 1.5 mostly exists already in WDM6's current tags; they just need the `>= 2` gate and an `i` loop to become the real thing.

Tier 2 needs a `wdm6_emit_diag_t2_blocksig_line` helper that doesn't exist, plus the `micro_diag_mode` plumbing WDM6 lacks entirely. That's the largest piece and isn't needed for the ice defect.
