# WDM6 literal routing audit — Phase 0a. CLOSED 2026-08-17 at `05a3e1a52c1c`.

Purpose, from `core-rules.md` "Numerical Fidelity Rules": the audit **is the
precondition for the double-literal (Option B) comparison build to mean
anything.** An unrouted literal is already a true double on the C++ side, so it
agrees with the Fortran under `-fdefault-real-8` by accident and disagrees only
under the float configuration. Every unrouted literal silently removes the
discrimination Option B exists to provide.

So this audit closes only when every float32-lossy literal is either **routed**
through `wdm6_literal()` / `wdm6_default_real_pow()` or **proved a non-masker**.

## Method — reproducible, not recalled

A literal is **LOSSY** when `float32(x) != double(x)`; exact integers and dyadic
fractions (`2.`, `4.`, `24.`, `0.5`, `1.e3`) are excluded per `core-rules.md:67`.
A lossy literal is **ROUTED** when it is an argument of `wdm6_literal(...)` or
`wdm6_default_real_pow(...)`, or carries an `f` suffix.

The scan excludes comments, string literals (printf formats carry `%12.4E`), and
`#if 0` regions. **The `#if 0` exclusion is the whole reconciliation**: counting
dead code gives 123 lossy / 40 unrouted, and the three disabled blocks at
`ERF_AdvanceWDM6.cpp:4721-4800`, `:4811-4822` and `:4838-4917` hold exactly 24
lossy literals. They compile to nothing, so they cannot mask anything.

| File | lossy | routed | unrouted |
| --- | --- | --- | --- |
| `ERF_AdvanceWDM6.cpp` | 99 | 83 | 16 |
| `ERF_WDM6.H` | 23 | 23 | 0 |

Independently reproduced 2026-08-17 by a scanner written from the rule text
without reference to the prior count, matching the earlier Codex inventory on
every cell.

## The 16 unrouted literals are all non-maskers

Three disjoint classes, each with a structural argument. None is a
"probably fine" judgement about magnitude.

### Class 1 — unused helper (7 literals)

`wdm6_mean_droplet_diameter` (`:498`) and `wdm6_ccn_activation` (`:518`) have
**zero call sites** anywhere in the tree (`grep` over `Source/`, `Exec/`, all
`*.cpp *.H *.F90 *.h`, excluding `Submodules/`; the only hits are the two
definitions). Both are simplified formulations, not ports of a Fortran routine —
`wdm6_ccn_activation` is labelled "simplified Abdul-Razzak & Ghan 2000" and has
no Fortran counterpart to be faithful to.

| line | literal |
| --- | --- |
| 501 | `1.e-9` |
| 508 | `3.14159265359` |
| 531 | `1.e-8` |
| 535 | `0.0048` |
| 539 | `0.6` |
| 539 | `0.01` |
| 543 | `0.1` |

Non-masker because dead code contributes nothing to either leg's output. They
are compiled but never called, so no value they compute is observable.

### Class 2 — bridge-passed constant (5 literals)

| line | literal | C++ symbol |
| --- | --- | --- |
| 576 | `1.28` | `den0` → `mp_wdm6_init_c` (`:583`) |
| 625 | `273.15` | `t0c` → `mp_wdm6_run_c` (`:785`) |
| 628 | `1.0e-12` | `qmin` → `mp_wdm6_run_c` (`:785`) |
| 632 | `1.28` | `den0` → `mp_wdm6_run_c` (`:786`) |
| 637 | `610.78` | `psat` → `mp_wdm6_run_c` (`:786`) |

In the Fortran these four names are **`intent(in)` dummy arguments only** —
`ERF_module_mp_wdm6.F90:449`, `:451`, `:453`, `:3253`, `:4406`, `:4443-4444`.
Verified there is no local `PARAMETER` declaration and no raw literal
counterpart for `den0`, `qmin` or `psat` anywhere in the module.

So the value is defined once, in C++, and **the bridge hands that same double to
the Fortran**. There is no independent Fortran literal on the other side to
disagree with: whatever C++ computes, both legs use. Routing them would shift
both legs identically and change no comparison, in either configuration. They
carry zero discriminating power, which is exactly why they are not maskers.

**The `t0c` shadowing hazard is real but already closed.** The four slope
routines — `slope_wdm6`, `slope_rain`, `slope_snow`, `slope_graup` — each
declare their own `real(kind=kind_phys), PARAMETER :: t0c = 273.15`
(`F90:3393`, `:3475`, `:3520`, `:3568`), shadowing the passed-in dummy. That
local literal is unsuffixed, so the Fortran uses `float32(273.15)` inside those
routines and the caller's true double everywhere else — the Fortran is
internally inconsistent about `t0c`, and a faithful port must be too. The port
already handles this with a separately routed constant,
`wdm6_slope_t0c = wdm6_literal(273.15)` at `ERF_AdvanceWDM6.cpp:196`, used at
the six slope call sites (`:1461`, `:1754`, `:1834`, `:2113`, `:2607`, `:4458`)
and nowhere else. `t0c` at `:625` is therefore the caller's value, correctly
left as a true double. This is the `core-rules.md:62` "check for shadowed
constants on both sides" clause, discharged.

### Class 3 — structurally unreachable defensive floor (4 literals)

| line | literal | expression |
| --- | --- | --- |
| 2039 | `1.0e-30` | `xni_safe = max(xni_arr(i,j,k), 1.0e-30)` |
| 2062 | `1.0e-30` | `xni_s    = max(xni_arr(i,j,k), 1.0e-30)` |
| 3210 | `1.0e-30` | `xni_safe = max(xni_arr(i,j,k), 1.0e-30)` |
| 3654 | `1.0e-30` | `xni_safe = max(xni, 1.0e-30)` |

All four floor `xni` before the division `den*qi/xni_safe`. **`1.e-30` does not
appear anywhere in `ERF_module_mp_wdm6.F90`** — these are C++-only guards with
no Fortran counterpart.

The floor is unreachable. `xni_arr` has exactly two writers and both apply the
identical clamp into `[1e3, 1e6]`:

- `:1396` — `xni_arr = wdm6_xni_exact(...)`, whose body is
  `min(max(5.38e7*temp, 1.e3), 1.e6)` (`:134`). **Unconditional** within
  `ParallelFor(box, ...)` at `:1385`, so every cell of `box` is initialised.
- `:3203` — `min(max(5.38e7*temp, 1.e3), 1.e6)`, same clamp, conditional, and
  can only overwrite with another value from the same interval.

Both kernels iterate the same `box` as all four readers, and `:1396` precedes
every read in program order. So `xni >= 1.e3` at each site — **33 orders of
magnitude above the floor** — and `amrex::max` returns `xni` in every case. The
`1.0e-30` operand is never selected, so its float32-vs-double value is not
observable.

This survives the two clamp/gate rules deliberately, not by omission:
- `core-rules.md:72` (a literal in a clamp is a null result until something
  reaches it) — the argument here is not "no cell reached it in this run" but
  that the clamp on `xni` makes reaching it impossible for any input. Read from
  the source, not from a metric.
- `core-rules.md:82` (floor-and-gate pairs) — checked and does not apply: there
  is no companion `if (x > 1.0e-30)` gate, so no threshold is shared and no
  branch can diverge. The literal appears only as the floor operand.

## What was routed in this phase — `05a3e1a52c1c`

The commit routes the maskable remainder and closes two sites that had the
literal contract hardwired.

Two sites computed pi as `Real(4.0f * std::atan(1.0f))` — float32 **regardless of
the flag**, so they would have silently fought a `-fdefault-real-8` build. Both
now use `pi_wdm6_loc = m_pi_wdm6`, i.e. `wdm6_literal(4.0*atan(1.0))`
(`ERF_InitWDM6.cpp:203`), which follows the contract like the other 83.
Measured: under Option A both forms give `3.14159274101257324219E+00`
bit-for-bit, because scaling by 4 is exact in binary, so **the substitution is
value-preserving at the current gate configuration**; under Option B it becomes
the true double `3.14159265358979311600E+00`, which is the point.

`wdm6_default_real_pow(24.0, 0.3333333)` replaces
`Real(std::pow(24.0f, 0.3333333f))` at `:2638` and `:4420`. The Fortran is
`((24.)**(.3333333))` at `F90:1645` and `:3002` — both operands unsuffixed, so
the expression is evaluated in single precision **and its result is single**
before widening, per the whole-expression rule at `core-rules.md:55`. The helper
reproduces that: `Real(std::pow(float, float))` resolves to `powf`, whose result
is float. Value-preserving under Option A (`2.88449907302856445312E+00` both
ways), and correctly double under Option B.

Also routed: the `sr` denominator `1.e-12` at `:2195`. Note this is a genuine
unsuffixed Fortran literal (`F90:1391`, `:1393`) and is **not** the same object
as the `qmin = 1.0e-12` bridge constant at `:628`, despite the identical
numeric value — same number, different provenance, opposite treatment.

## Gate result

Re-gated the full serial A/B/C ladder at `05a3e1a52c1c`, campaign
`plotfile_parity_wdm6_bubble_v2s`, provenance repaired so the embedded
`26.06-384-g05a3e1a52c1c` matches the clean worktree SHA. Compared with
`Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex` (preferred when provenance
matters, being reproducible from tracked submodule source).

Prediction before running, from the bitwise equivalence above: **every figure
unchanged.** That is what a value-preserving routing change must produce, and
confirming it is what distinguishes "routed correctly" from "routed and
perturbed something".

| Milestone | recorded at `831690c21634` | v2s at `05a3e1a52c1c` | verdict |
| --- | --- | --- | --- |
| A (step 1) | PLOTFILE AGREE, 32/32 bitwise | PLOTFILE AGREE, 32/32 bitwise | unchanged |
| B (step 10) | `qt` 2.845847066e-16, 16/32 bitwise | `qt` 2.845847066e-16, 16/32 bitwise | unchanged, all ten digits |
| C (step 100) | `nn` 1.285451261e-11 at 128,3,81, 5/32 | see `campaign_executions.tsv` | — |

## Consequence for Option B

`ERF_AdvanceWDM6.cpp` and `ERF_WDM6.H` now carry **no lossy literal that can
mask an Option B discriminator**: 106 of 122 routed, and the 16 remaining are
dead code (7), bridge-passed and therefore shared by both legs (5), or an
unreachable floor operand (4). None can produce accidental agreement under
`-fdefault-real-8`.

Option B can now discriminate, so `plotfile_parity_wdm6_bubble_optionb_v1a` is
unblocked on this prerequisite.

**Still open, and not closed by this audit** (`ERF_WDM6.H:104-106`): `r0`,
`peaut` and `alpha_wdm6` have no matching Fortran parameter of the same name and
value. They are routed for uniformity, so they cannot mask; but their
*correspondence* is unverified, which is a separate question from routing and
remains open.
