
## Rule 5: WSM6-Specific Structural Differences [PARTIAL]

These are cases where the skill template needs extension beyond 
the Morrison baseline.

### 5a: 2D array indexing (its:ite, kts:kte)
Morrison was called column-by-column (KTS:KTE only).
WSM6 dimension(its:,:) processes multiple columns simultaneously.
AMReX kernel loop already handles (i,j,k) — the its/ite dimension 
maps naturally to the i-loop in ParallelFor. No special handling needed 
but the MicVar enum sizing must account for the full box, not a column.
[COMPLETE — resolved by AMReX ParallelFor pattern]

### 5b: Optional arguments
WSM6 has: snow, snowncv, graupel, graupelncv, rainprod2d, evapprod2d 
as Fortran optional.
[PENDING — need to check MYNN 2.5 for precedent on how ERF handles 
optional Fortran args. Candidates: nullable Array4, compile-time ifdef, 
boolean SolverChoice flag]

### 5c: errmsg/errflg
WSM6 returns error state via character(len=*) errmsg + integer errflg.
[PENDING — determine ERF convention. Likely exception or return code.]

### 5d: vrec/vsqrt from module_libmassv
WSM6 uses vectorized math from module_libmassv.
In C++ these become std::pow / std::sqrt inside ParallelFor 
(AMReX GPU kernels do not use libmassv).
[COMPLETE — straightforward substitution]

---

## MicVar_WSM6 Enum Draft [PARTIAL]

```cpp
namespace MicVar_WSM6 {
    enum {
        // independent variables
        rho = 0,   // density
        theta,     // liquid/ice water potential temperature
        tabs,      // temperature
        pres,      // pressure
        // moisture state
        qv,        // water vapor mixing ratio  (= q)
        qcl,       // cloud water               (= qc)
        qci,       // cloud ice                 (= qi)
        qpr,       // rain                      (= qr)
        qps,       // snow                      (= qs)
        qpg,       // graupel                   (= qg)
        // number concentrations — [PENDING: WSM6 is single-moment, 
        //   confirm no number conc prognostic vars]
        // surface accumulations
        rain_accum,
        snow_accum,
        graup_accum,
        //
        NumVars
    };
}
```

## WSM6Ind Enum Draft [PARTIAL]

```cpp
namespace WSM6Ind {
    enum {
        // state (mirrors Fortran inout args)
        t = 0,         // temperature
        q,             // water vapor
        qc,            // cloud water
        qi,            // cloud ice
        qr,            // rain
        qs,            // snow
        qg,            // graupel
        den,           // density (intent in, but needed in kernel)
        p,             // pressure
        delz,          // layer thickness
        // surface (1D in Fortran — handle separately)
        // rain, rainncv, sr, snow, snowncv, graupel, graupelncv
        // in-kernel derived (temperature dependent)
        xlv,           // latent heat of vaporization
        xls,           // latent heat of sublimation
        xlf,           // latent heat of fusion
        // [PENDING — full internal variable list from mp_wsm6_run body]
        NumInds
    };
}
```

---

## Pending Before WSM6 Implementation is Complete

1. **cice** — specific heat of ice, no Morrison analog, need value 
   and placement decision
2. **Optional arg pattern** — check MYNN 2.5 for ERF convention
3. **errmsg/errflg** — ERF error handling convention
4. **WSM6Ind complete** — need interior of mp_wsm6_run to capture 
   all internal working variables
5. **MYNN 2.5 cross-reference** — validate two-layer enum and 
   constants patterns match before trusting WSM6 application

---

## Files to Create for WSM6

- `Source/Microphysics/WSM6/ERF_WSM6.H`
- `Source/Microphysics/WSM6/ERF_InitWSM6.cpp`
- `Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp`

Pattern: copy Morrison file structure, swap enum names and variable 
mappings per this document.

---

## Rule 6: 3D Species-Indexed Arrays → Expanded Enum Entries [NEW - WSM6]

WSM6 has arrays dimensioned `(its:ite, kts:kte, 3)` where the third 
dimension indexes precip species: 1=rain, 2=snow, 3=graupel.

Variables: rh, qsat, rslope, rslope2, rslope3, rslopeb, qrs_tmp, 
           falk, fall, work1

In AMReX these cannot be a runtime-indexed third dimension in the 
FArrayBox component space cleanly. Resolution: expand into 3 separate 
enum entries per variable with suffix _r, _s, _g.

Example:
    rslope_r, rslope_s, rslope_g,
    rslope2_r, rslope2_s, rslope2_g,
    rslope3_r, rslope3_s, rslope3_g,
    rslopeb_r, rslopeb_s, rslopeb_g,
    // etc.

This adds ~30 entries to WSM6Ind but keeps the kernel access pattern 
uniform: wsm6_arr(i,j,k,WSM6Ind::rslope_r)

---

## Rule 7: 1D (its:ite) Arrays → Kernel-Local Scalars [NEW - WSM6]

Fortran arrays dimensioned only on (its:ite) are column-reduction or 
sedimentation sub-stepping variables. In AMReX each (i,j) column is 
an independent kernel thread, so these collapse to kernel-local scalars 
inside the ParallelFor lambda.

Variables affected:
    delqrs1, delqrs2, delqrs3, delqi   // sedimentation increments
    tstepsnow, tstepgraup               // sub-step timesteps
    dvec1, tvec1                        // vrec/vsqrt optimization 
                                        // (eliminated entirely — 
                                        //  replaced by std::sqrt)

Do NOT add these to WSM6Ind. Declare as Real inside the kernel.

---

## Rule 8: Integer and Logical Arrays → Kernel-Local int/bool [NEW - WSM6]

Fortran integer(its:ite) and logical(its:ite) arrays follow the same 
rule as Rule 7 — they become kernel-local typed variables.

    mstep, numdt  → int mstep, numdt inside ParallelFor
    flgcld        → bool flgcld inside ParallelFor

Do NOT store in the FArrayBox (Real-typed). Keep as kernel-local 
primitive types.

---

## WSM6Ind Enum — Updated [PARTIAL → MORE COMPLETE]

```cpp
namespace WSM6Ind {
    enum {
        // --- inout state (from argument list) ---
        t = 0,      // temperature
        q,          // water vapor mixing ratio
        qc,         // cloud water
        qi,         // cloud ice
        qr,         // rain
        qs,         // snow
        qg,         // graupel
        den,        // density (intent in, needed in kernel)
        p,          // pressure
        delz,       // layer thickness

        // --- sedimentation working arrays (its:ite,kts:kte) ---
        fallc, falkc, work1c, work2c, workr, worka,
        den_tmp, delz_tmp,

        // --- process rates (its:ite,kts:kte) ---
        pigen,      // ice nucleation
        pidep,      // ice deposition
        pcond,      // condensation
        prevp,      // rain evaporation
        psevp,      // snow evaporation
        pgevp,      // graupel evaporation
        psdep,      // snow deposition
        pgdep,      // graupel deposition
        praut,      // rain autoconversion
        psaut,      // snow autoconversion
        pgaut,      // graupel autoconversion
        piacr,      // rain-ice accretion
        pracw,      // rain-cloud accretion
        praci,      // rain-ice collection
        pracs,      // rain-snow collection
        psacw,      // snow-cloud accretion
        psaci,      // snow-ice accretion
        psacr,      // snow-rain accretion
        pgacw,      // graupel-cloud accretion
        pgaci,      // graupel-ice accretion
        pgacr,      // graupel-rain accretion
        pgacs,      // graupel-snow accretion
        paacw,      // accretion of cloud by all
        psmlt,      // snow melting
        pgmlt,      // graupel melting
        pseml,      // snow evap/melt latent
        pgeml,      // graupel evap/melt latent

        // --- derived thermodynamic (its:ite,kts:kte) ---
        qsum,
        xl,         // latent heat (mixed phase)
        cpm,        // effective cp (moist)
        work2,
        denfac,     // density factor
        xni,        // ice number concentration
        denqrs1,    // density * qr
        denqrs2,    // density * qs
        denqrs3,    // density * qg
        denqci,     // density * qi
        n0sfac,     // snow intercept factor

        // --- 3D species arrays expanded (its:ite,kts:kte,3→_r/_s/_g) ---
        rh_r,    rh_s,    rh_g,
        qsat_r,  qsat_s,  qsat_g,
        rslope_r,  rslope_s,  rslope_g,
        rslope2_r, rslope2_s, rslope2_g,
        rslope3_r, rslope3_s, rslope3_g,
        rslopeb_r, rslopeb_s, rslopeb_g,
        qrs_tmp_r, qrs_tmp_s, qrs_tmp_g,
        falk_r,  falk_s,  falk_g,
        fall_r,  fall_s,  fall_g,
        work1_r, work1_s, work1_g,

        NumInds
    };
}

## Rule 9: Module save Variables → Class Member Init [NEW - WSM6]

WSM6 has a separate `wsm6init` subroutine that pre-computes derived 
coefficients from module parameter constants and stores them in Fortran 
`save` variables. These persist across all calls to mp_wsm6_run.

In C++ these become:
- Private class member variables in ERF_WSM6.H
- Computed once in ERF_InitWSM6.cpp (called from Initialize(), 
  not from Advance())
- Passed into the kernel via capture in the ParallelFor lambda

Module `parameter` constants (lines 15-43) → `constexpr` in ERF_WSM6.H
Module `save` variables (lines 45-63) → class member Reals, 
  initialized in ERF_InitWSM6.cpp

[PENDING — wsm6init signature needed to confirm whether init 
requires runtime arguments or computes purely from parameters]

## Rule 9: Module save Variables → Class Member Init [COMPLETE]

mp_wsm6_init (= ERF's Initialize()) takes a strict subset of the 
mp_wsm6_run arguments: `den0, denr, dens, cl, cpv, hail_opt`

All map to already-resolved constants:
- den0  → m_rhosu (class member)
- denr  → m_rhow  (class member)  
- dens  → constexpr Real dens_snow = Real(100.0) (ERF_WSM6.H)
- cl    → m_cpw   (class member, = cliq)
- cpv   → Cp_v    (ERF_Constants.H)

This means ERF_InitWSM6.cpp can compute all save variables purely 
from ERF_Constants.H constexpr values + class members already set 
in Define(). No additional arguments needed at Initialize() time.

Port the exact Fortran `rgmma` loop in ERF_InitWSM6.cpp; do not use
`std::tgamma` during parity validation.

Full save variable list becomes private class members in ERF_WSM6.H,
computed once in Initialize(), captured by value in ParallelFor lambda.

---

## Rule 10: Runtime Physics Variant Scope [WSM6 current + future]

`hail_opt` is an integer flag that switches WSM6 between two physics 
regimes — hail vs graupel. It controls 5 save variables:

| Variable   | hail_opt=1 (hail) | hail_opt=0 (graupel) |
|------------|-------------------|----------------------|
| n0g        | 4.e4              | 4.e6                 |
| deng       | 700.              | 500.                 |
| avtg       | 285.0             | 330.0                |
| bvtg       | 0.8               | 0.8  (same)          |
| lamdagmax  | 2.e4              | 6.e4                 |

Current WSM6 parity campaign policy:
- Expose only `MoistureType::WSM6` (default graupel path).
- Set `m_hail_opt = false` and use `hail_opt = 0` in init bridge flow.
- Do not add `WSM6_Hail` yet.
- Do not add `WSM6_NoIce`.
- Do not add `SolverChoice::use_hail` or `SolverChoice::wsm6_hail_opt`.

Future extension policy (separate lane):
- If hail-mode support is needed later, add `MoistureType::WSM6_Hail`.
- Map `WSM6` and `WSM6_Hail` to the same WSM6 model class.
- Derive `m_hail_opt` from `sc.moisture_type` in `Define()`.
- Treat hail-mode validation as a separate campaign with separate
  acceptance criteria.
- Until that lane exists, all WSM6 acceptance criteria in this campaign
  are for default graupel mode only.

Scope precedent note:
- ERF Morrison exposes `Morrison` and `Morrison_NoIce` as user-visible
  `MoistureType` variants.
- Morrison's rimed-ice/hail switch (`morr_rimed_ice` / `IHAIL`) exists
  internally but is not exposed as `Morrison_Hail`; ERF currently runs
  graupel mode.
- WSM6 follows the same scope discipline for this campaign.

WSM6_Hail deferred:
- ERF Morrison currently hard-codes rimed/hail mode to graupel while
  exposing only `Morrison_NoIce` as a public variant.
- To match that scope discipline, the current WSM6 port validates
  default graupel mode only.
- Hail mode is a future extension, not part of the current acceptance
  criterion.

---

## ERF_WSM6.H Class Member List [PARTIAL → NEAR COMPLETE]

```cpp
// constexpr module parameters (from Fortran parameter declarations)
static constexpr amrex::Real dtcldcr    = Real(120.0);
static constexpr amrex::Real n0r        = Real(8.e6);
static constexpr amrex::Real avtr       = Real(841.9);
static constexpr amrex::Real bvtr       = Real(0.8);
static constexpr amrex::Real r0         = Real(0.8e-5);
static constexpr amrex::Real peaut      = Real(0.55);
static constexpr amrex::Real xncr       = Real(3.e8);
static constexpr amrex::Real xmyu       = Real(1.718e-5);
static constexpr amrex::Real avts       = Real(11.72);
static constexpr amrex::Real bvts       = Real(0.41);
static constexpr amrex::Real lamdarmax  = Real(8.e4);
static constexpr amrex::Real lamdasmax  = Real(1.e5);
static constexpr amrex::Real dicon      = Real(11.9);
static constexpr amrex::Real dimax      = Real(500.e-6);
static constexpr amrex::Real pfrz1      = Real(100.0);
static constexpr amrex::Real pfrz2      = Real(0.66);
static constexpr amrex::Real qcrmin     = Real(1.e-9);
static constexpr amrex::Real eacrc      = Real(1.0);
static constexpr amrex::Real dens_snow  = Real(100.0);
static constexpr amrex::Real qs0        = Real(6.e-4);
static constexpr amrex::Real n0smax     = Real(1.e11);
static constexpr amrex::Real n0s        = Real(2.e6);
static constexpr amrex::Real alpha      = Real(0.12);

// runtime physics variant (current campaign fixed to graupel mode)
bool m_hail_opt = false;

// hail_opt-dependent class members
amrex::Real m_n0g, m_deng, m_avtg, m_bvtg, m_lamdagmax;

// save variables computed in Initialize()
amrex::Real m_qc0, m_qck1, m_pidnc;
amrex::Real m_bvtr1, m_bvtr2, m_bvtr3, m_bvtr4, m_bvtr6;
amrex::Real m_g1pbr, m_g3pbr, m_g4pbr, m_g5pbro2, m_g6pbr;
amrex::Real m_pvtr, m_eacrr, m_pacrr;
amrex::Real m_precr1, m_precr2, m_roqimax;
amrex::Real m_bvts1, m_bvts2, m_bvts3, m_bvts4;
amrex::Real m_g1pbs, m_g3pbs, m_g4pbs, m_g5pbso2;
amrex::Real m_pvts, m_pacrs, m_precs1, m_precs2;
amrex::Real m_pidn0r, m_pidn0s, m_pacrc;
amrex::Real m_bvtg1, m_bvtg2, m_bvtg3, m_bvtg4;
amrex::Real m_g1pbg, m_g3pbg, m_g4pbg, m_g5pbgo2;
amrex::Real m_pvtg, m_pacrg, m_precg1, m_precg2, m_pidn0g;
amrex::Real m_rslopermax, m_rslopesmax, m_rslopegmax;
amrex::Real m_rsloperbmax, m_rslopesbmax, m_rslopegbmax;
amrex::Real m_rsloper2max, m_rslopes2max, m_rslopeg2max;
amrex::Real m_rsloper3max, m_rslopes3max, m_rslopeg3max;
amrex::Real m_xlv1, m_pi_wsm6;

// physical constants (class members per Rule 4)
amrex::Real m_rhosu, m_rhow, m_cpw;
amrex::Real m_ep_2, m_qsmall;
```

## Appendix A: Plotfile Parity Campaign — WSM6 SquallLine_2D

This appendix is WSM6-specific implementation guidance for the generic
Rules 35-38 defined in `fortran_to_cpp_microphysics_skill.md`.

### Scope gate

This appendix covers serial parity only:
- single-rank
- CPU path
- non-MPI executable mode (`ERF3d.gnu.DEBUG.ex`)

MPI and GPU campaign lanes are deferred until serial C++ parity is
verified and accepted.

### Fcompare binary

Repository-relative location:

```bash
Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex
```

### Executable and run policy

Run from: Exec/CanonicalTests/SquallLine_2D
Input file: inputs_moisture_WSM6
Use unique <campaign_id> root per campaign.

  Tag-frontier runs (Rules 30-34): ../../ERF3d.gnu.DEBUG.ex
  Plotfile short pair (Milestones A-B): ../../ERF3d.gnu.DEBUG.ex
  Plotfile long pair (Milestones C-I): ../../ERF3d.gnu.TEST.ex

Never run long plotfile campaigns with the DEBUG binary —
at ~38s/step in DEBUG, stop_time=9000 at dt=1.0 is impractical.
RELEASE and DEBUG may differ due to optimization and floating
point contraction — always validate Milestones A and B at
RELEASE before launching the long pair.

Build commands:
  DEBUG:   make -j8 USE_MPI=FALSE USE_NETCDF=FALSE \
                    USE_WSM6_FORT=TRUE DEBUG=TRUE
  RELEASE: make -j8 USE_MPI=FALSE USE_NETCDF=FALSE \
                    USE_WSM6_FORT=TRUE DEBUG=FALSE

Build type convention for run_id and notes:
  Encode build type in run_id token: _debug or _release
    e.g. 20260429T_plotparity_v1c_A_fortran_release
  Record in notes as build_type=debug or build_type=release
  Do not add build_type as a TSV column — notes token is enough.
  Rule: tag_frontier rows are always debug.
        Milestone A-B rows may be debug or release — record both
        if both are run, as separate rows.
        Milestone C-I rows must be release.

Set GFORTRAN_UNBUFFERED_ALL=y for Fortran-path runs only.
No mpirun for any run in this campaign.

### dt convention

Use `fixed_dt=1.0` for this campaign family so milestone times are
reproducible and aligned with existing WSM6 comparison practice.

### Milestone table (WSM6 SquallLine_2D)

| Milestone | Stop | Plot step | Active species / regime |
| --- | --- | --- | --- |
| A | step 1 | 1 | none |
| B | step 10 | 10 | none |
| C | t=300 s | 300 | qc, qr |
| D | t=500 s | 500 | qc, qr, qg |
| E | t=600 s | 600 | qc, qr, qg, qs, qi |
| F | t=3000 s | 3000 | all (paper comparison 1) |
| G | t=6000 s | 6000 | all (paper comparison 2) |
| H | t=9000 s | 9000 | all (paper comparison 3) |
| I | t=12000 s | 12000 | all (extended drift) |

With `fixed_dt=1.0` and `erf.plot_int_1=100`, milestones at
`t=300/500/600/3000/6000/9000` are exact plot steps.

### Run structure

- Short pair: milestones A-B (`max_step=10`, `erf.plot_int_1=1`,
  `erf.check_int=1`).
- Long pair: milestones C-I (`stop_time=12000`, `erf.plot_int_1=100`,
  `erf.check_int=100`).

### Command templates

Short run, Fortran path:

```bash
GFORTRAN_UNBUFFERED_ALL=y \
../../ERF3d.gnu.TEST.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=0 erf.microphysics_debug=0 \
  fixed_dt=1.0 max_step=10 \
  erf.plot_file_1=<campaign_id>/short_fortran/plt \
  erf.plot_int_1=1 \
  erf.check_file=<campaign_id>/short_fortran/chk \
  erf.check_int=1 \
  > <campaign_id>/log.short_fortran 2>&1
```

Short run, C++ path:

```bash
../../ERF3d.gnu.TEST.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=1 erf.microphysics_debug=0 \
  fixed_dt=1.0 max_step=10 \
  erf.plot_file_1=<campaign_id>/short_cpp/plt \
  erf.plot_int_1=1 \
  erf.check_file=<campaign_id>/short_cpp/chk \
  erf.check_int=1 \
  > <campaign_id>/log.short_cpp 2>&1
```

Long run, Fortran path:

```bash
GFORTRAN_UNBUFFERED_ALL=y \
../../ERF3d.gnu.TEST.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=0 erf.microphysics_debug=0 \
  fixed_dt=1.0 stop_time=12000 \
  erf.plot_file_1=<campaign_id>/long_fortran/plt \
  erf.plot_int_1=100 \
  erf.check_file=<campaign_id>/long_fortran/chk \
  erf.check_int=100 \
  > <campaign_id>/log.long_fortran 2>&1
```

Long run, C++ path:

```bash
../../ERF3d.gnu.TEST.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=1 erf.microphysics_debug=0 \
  fixed_dt=1.0 stop_time=12000 \
  erf.plot_file_1=<campaign_id>/long_cpp/plt \
  erf.plot_int_1=100 \
  erf.check_file=<campaign_id>/long_cpp/chk \
  erf.check_int=100 \
  > <campaign_id>/log.long_cpp 2>&1
```

### Per-milestone fcompare workflow

For matching plot step `N`:

```bash
Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex \
  <campaign_id>/<pair>_fortran/plt<N> <campaign_id>/<pair>_cpp/plt<N>
```

On first failing milestone (`N`, variable `<var>`):

```bash
Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex -z <var> \
  <fortran_plt_N> <cpp_plt_N>

Submodules/AMReX/Tools/Plotfile/fcompare.gnu.ex -d <var> \
  <fortran_plt_N> <cpp_plt_N>

# Archive default diff output immediately
mv diffs <campaign_id>/diffs_plt<N>_<var>
```

`-d` selects the diff variable; it does not take an output path argument.

### Restart reproducibility (Rule 36)

Restart from checkpoint at step `N-1` and advance exactly one step.
If `N-1` is not checkpointed, use nearest earlier `K` and run a bounded
refinement window to narrow first-fail substep `N*` before Rule 36 compare.
Use the same restart source path for both legs (default: Fortran checkpoint
path) to isolate implementation differences.

A typical ERF pattern is:

```bash
erf.restart=<campaign_id>/long_fortran/chk<N-1_padded> max_step=<N>
```

Use the same restart path on the C++ leg as well
(`erf.restart=<campaign_id>/long_fortran/chk<N-1_padded>`) and compare
regenerated `plt<N>` (or `plt<N*>` after refinement).
Record `restart_repro_status=MATCH|MISMATCH`.

### Bounded per-step refinement (coarse milestone failure)

When a long run uses coarse plot cadence (for example `erf.plot_int_1=100`)
and first fail appears at milestone step `N`, run a bounded restart refinement
from nearest checkpoint `K <= N-1`:

```bash
GFORTRAN_UNBUFFERED_ALL=y ../../ERF3d.gnu.TEST.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=0 erf.microphysics_debug=0 \
  erf.fixed_dt=1.0 amr.restart=<campaign_id>/long_fortran/chk<K_padded> \
  max_step=<N> stop_time=<N> \
  erf.plot_file_1=<campaign_id>/restart_stepwise/fortran/plt \
  erf.plot_int_1=1 \
  erf.check_file=<campaign_id>/restart_stepwise/fortran/chk \
  erf.check_int=1

../../ERF3d.gnu.TEST.ex inputs_moisture_WSM6 \
  erf.use_wsm6_cpp_answer=1 erf.microphysics_debug=0 \
  erf.fixed_dt=1.0 amr.restart=<campaign_id>/long_fortran/chk<K_padded> \
  max_step=<N> stop_time=<N> \
  erf.plot_file_1=<campaign_id>/restart_stepwise/cpp/plt \
  erf.plot_int_1=1 \
  erf.check_file=<campaign_id>/restart_stepwise/cpp/chk \
  erf.check_int=1
```

Then run explicit per-step paired `fcompare` from `K+1..N`, stop at earliest
failing substep `N*`, and use `N*` for tag retreat handoff.

### Soft reference baselines (early dry regime)

From terminal G16 dry-regime parity (`96b936d2`, 1-step):
- `qv` max abs about `5.20e-18`
- `t` max abs about `5.68e-14`
- `cpm` max abs about `4.55e-13`
- `xl` max abs about `4.66e-10`

Use as references for early dry milestones only; do not suppress active-regime
signal with tolerance flags.

### Variable-to-tag mapping (WSM6 manifest groups)

- `qr` -> NISLFV_R, PRAUT, PRACW, PREVP
- `qs` -> NISLFV_SG, PSACI, PRACS, PSAUT, PSEVP
- `qg` -> NISLFV_SG, PSACI, PRACS, PSAUT, PSEVP
- `qi` -> VICE, PRACI, PIDEP, PSAUT
- `qv`/`q` -> PREVP, PCOND, QSAT, QSAT2
- `t` -> MELT, PHASE, UPDATE, PCOND
- `fall_r/s/g` -> FALL
- `rho`/`den` -> SLOPE1, SLOPE2, SLOPE3

Canonical retreat order:
DENFAC -> QSAT -> SLOPE1 -> NISLFV_R -> NISLFV_SG -> FALL -> SLOPE2 ->
MELT -> VICE -> PRECIP -> PHASE -> SLOPE3 -> PRAUT -> PRACW -> PREVP ->
PRACI -> PSACI -> PRACS -> PSEML -> PIDEP -> PSAUT -> PSEVP -> UPDATE ->
QSAT2 -> PCOND

Use milestone active-species context to narrow the candidate set before
canonical-order retreat.
