
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

`rgmma` throughout init → std::tgamma() in C++

Full save variable list becomes private class members in ERF_WSM6.H,
computed once in Initialize(), captured by value in ParallelFor lambda.

---

## Rule 10: Runtime Physics Variant Flag → SolverChoice [NEW - WSM6]

`hail_opt` is an integer flag that switches WSM6 between two physics 
regimes — hail vs graupel. It controls 5 save variables:

| Variable   | hail_opt=1 (hail) | hail_opt=0 (graupel) |
|------------|-------------------|----------------------|
| n0g        | 4.e4              | 4.e6                 |
| deng       | 700.              | 500.                 |
| avtg       | 285.0             | 330.0                |
| bvtg       | 0.8               | 0.8  (same)          |
| lamdagmax  | 2.e4              | 6.e4                 |

In C++:
- Add `bool use_hail = false` to SolverChoice (parsed from input file)
- Pull in Define(): `m_hail_opt = sc.use_hail`
- Branch in Initialize() to set the 5 affected class members
- This is the only save variable that depends on a runtime flag — 
  all others compute purely from constants

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

// runtime physics variant (from hail_opt)
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
