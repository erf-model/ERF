# Fortran-to-C++ Microphysics Skill
# ERF/AMReX Framework Pattern
# Status: PARTIAL — validated on Morrison, pending MYNN 2.5 + WSM6

---

## Overview

This skill documents the repeatable pattern for porting WRF-style
Fortran microphysics schemes into ERF so that BOTH of the following
work simultaneously:

  Path A — Fortran Bridge: wire the existing Fortran source into ERF
  with minimal changes, behind a compile flag. This is always done
  first and serves as the validation ground truth.

  Path B — Native AMReX C++: replace the Fortran call with a
  ParallelFor GPU kernel using ERF-idiomatic patterns. Validated
  against Path A via amrex_fcompare.

Validated instances:
  Morrison: Path B complete (native C++ advance, Fortran init bridge)
  MYNN 2.5: Path B complete (GPU device struct strategy)
  WSM6: Path A complete (bridge smoke tested), Path B in progress

Rules below are tagged [Path A], [Path B], or [Both].

---


## Path A Rules: Fortran Bridge

### Path A — Rule 1: MoistureType Registration [COMPLETE]
Every new scheme requires exactly four wiring points:

1. Add entry to AMREX_ENUM(MoistureType,...) in ERF_DataStruct.H
   Example: WSM6, WSM6_Hail (variant flags become separate entries,
   not SolverChoice booleans — see Morrison_NoIce precedent)

2. Add Eulerian branch in ERF_EulerianMicrophysics.H constructor:
   } else if (a_model_type == MoistureType::WSM6) {
       SetModel<WSM6>();
   }

3. Add to modelType() in ERF_Microphysics.H Eulerian list:
   || (a_moisture_type == MoistureType::WSM6)

4. Create Source/Microphysics/<Scheme>/Make.package with:
   CEXE_sources += ERF_Init<Scheme>.cpp
   CEXE_sources += ERF_Advance<Scheme>.cpp
   CEXE_sources += ERF_Update<Scheme>.cpp
   CEXE_headers += ERF_<Scheme>.H

Validated: WSM6 (dc2968f8), Morrison, SAM, Kessler all follow
this exact pattern.

### Path A — Rule 2: Isohelper Fortran Bridge Pattern [COMPLETE]
The Fortran source is never modified. Instead:

1. Copy Fortran source into Source/Microphysics/<Scheme>/
   Rename with ERF_ prefix: ERF_module_mp_<scheme>.F90

2. Create ERF_module_mp_<scheme>_isohelper.F90 with bind(C) wrappers:
   subroutine mp_<scheme>_init_c(...) bind(C, name="mp_<scheme>_init_c")
     use iso_c_binding
     call mp_<scheme>_init(...)  ! calls original subroutine unchanged
   end subroutine

3. Create ERF_<Scheme>_Fortran_Interface.H with extern "C" declarations:
   extern "C" {
     void mp_<scheme>_init_c(...);
     void mp_<scheme>_run_c(...);
   }

4. Add conditional Fortran sources to Make.package:
   ifeq ($(USE_<SCHEME>_FORT),TRUE)
       F90EXE_sources += ERF_module_mp_<scheme>.F90
       F90EXE_sources += ERF_module_mp_<scheme>_isohelper.F90
   endif

The isohelper converts C-compatible types (c_double, c_int with VALUE)
to Fortran intrinsic types. errmsg/errflg from CCPP-style interfaces
are handled in the isohelper — on errflg /= 0, write and stop.
Never propagate errmsg/errflg into the C++ layer.

### Path A — Rule 3: Bridge Caller in Advance() [COMPLETE]
The C++ Advance() function gates on the compile flag:

   #ifdef ERF_USE_<SCHEME>_FORT
       // pack arrays, call Fortran, unpack
       mp_<scheme>_run_c(
           t_arr.dataPtr(), qv_arr.dataPtr(), ...,
           ims, ime, jms, jme, kms, kme,
           its, ite, jts, jte, kts, kte);
   #else
       amrex::Abort("<Scheme> native C++ not yet implemented");
   #endif

Array dimensions passed to Fortran:
  ims/ime/jms/jme/kms/kme → fab_box bounds (memory/halo extent)
  its/ite/jts/jte/kts/kte → box bounds (compute extent)
This matches WRF's memory vs. tile dimension convention exactly.

GNUmake development build:
   cd Exec/  (or Exec/<TestCase>/)
   make -j8 USE_<SCHEME>_FORT=TRUE USE_NETCDF=FALSE 2>&1 | tee make.ou

Smoke test: run 5 steps, confirm no abort, confirm plotfile written.

---

## Path B Rules: Native AMReX C++

(leave all existing Rule 1 through Rule 14 content exactly as-is
below this line — do not modify any of it)

---

## Rule 1: Two-Layer Enum Pattern [COMPLETE]

Every scheme produces two enums, never one.

### Layer 1: `MicVar_<Scheme>` — Persistent ERF State
- Lives in `ERF_<Scheme>.H`
- Uses ERF-idiomatic physics names (e.g. `qcl`, `qci`, `qpr`, `qps`, `qpg`)
- Backed by `amrex::Array<FabPtr, NumVars> mic_fab_vars` (one MultiFab per var)
- Initialized in `ERF_Init<Scheme>.cpp` via `make_shared<MultiFab>(...)`
- These are the variables ERF owns across timesteps (checkpoint/restart)

### Layer 2: `<Scheme>Ind` — Kernel Working Array
- Lives in `ERF_Advance<Scheme>.cpp`
- Uses Fortran-mirrored names 1:1 (e.g. `qc3d`, `qi3d`, `qr3d`, `t3d`)
- Backed by a single `FArrayBox <scheme>_fab(grown_box, NumInds)`
- Allocated fresh each Advance() call, not persistent
- Includes tendencies, slope parameters, derived quantities, 
  and temporaries that the Fortran kernel computed internally

### Why two layers:
ERF state uses clean physics names for readability and restartability.
The kernel working array mirrors Fortran exactly to make the conversion 
auditable and the Fortran source directly comparable line-by-line.

---

## Rule 2: Pack/Unpack Bridge Pattern [COMPLETE]

Between the two layers, there is always an explicit pack step before 
the kernel and an unpack step after.

### Pack (ERF state → kernel array):
```cpp
// From ERF_AdvanceMorrison.cpp:1071
morr_arr(i,j,k,MORRInd::qc3d)  = qcl_arr(i,j,k);
morr_arr(i,j,k,MORRInd::qi3d)  = qci_arr(i,j,k);
morr_arr(i,j,k,MORRInd::qni3d) = qps_arr(i,j,k);
morr_arr(i,j,k,MORRInd::qr3d)  = qpr_arr(i,j,k);
morr_arr(i,j,k,MORRInd::t3d)   = theta_arr(i,j,k) * pii_arr(i,j,k);
morr_arr(i,j,k,MORRInd::qv3d)  = qv_arr(i,j,k);
morr_arr(i,j,k,MORRInd::pres)  = pres_arr(i,j,k);
// ... etc
```

### Unpack (kernel array → ERF state):
Reverse mapping after kernel completes. Tendencies are applied to 
mic_fab_vars before returning from Advance().

### Skill action:
Given Fortran argument list, generate both directions of this mapping 
automatically using the name correspondence table.

---

## Rule 3: Entry Point Pattern [COMPLETE]

Every scheme implements the `NullMoist` interface:

```cpp
class <Scheme> : public NullMoist {
public:
    void Define(SolverChoice& sc) override;  // pull control flags only
    void Advance(const amrex::Real& dt_advance,
                 const SolverChoice& sc) override;
};
```

`Define()` pulls only control flags from SolverChoice 
(`ave_plane`, `rdOcp`, `use_shoc`). It does not pull physics constants.

---

## Rule 4: Physical Constants Resolution Hierarchy [COMPLETE]

WSM6 passes all physical constants as explicit Fortran arguments. 
The C++ equivalent resolves them via this priority order:

### Tier 1: ERF_Constants.H constexpr (prefer these)
```cpp
// Use when the value is universal and matches ERF conventions
constexpr amrex::Real R_d        = amrex::Real(287.0);
constexpr amrex::Real R_v        = amrex::Real(461.505);
constexpr amrex::Real Cp_d       = amrex::Real(1004.5);
constexpr amrex::Real Cp_v       = amrex::Real(1859.0);
constexpr amrex::Real CONST_GRAV = amrex::Real(9.81);
constexpr amrex::Real tmelt      = amrex::Real(273.15);  // = t0c
constexpr amrex::Real p_0        = amrex::Real(1.0e5);
```

### Tier 2: Class member variables (for derived or scheme-tuned values)
```cpp
// Set in Advance() preamble, not in Define()
m_ep_2   = R_d / R_v;          // = ep2 (note: Morrison uses Rd/Rv, 
                                //   WSM6 ep2 may differ — verify)
m_qsmall = Real(1.0E-14);      // = qmin
m_rhosu  = Real(85000.0)/...;  // = den0 reference density
m_rhow   = Real(997.0);        // = denr liquid water density
m_cpw    = Real(4187.0);       // = cliq
```

### Tier 3: In-kernel computation (for temperature-dependent quantities)
```cpp
// Computed inside ParallelFor as part of SchemeInd working array
// xls/xlv0/xlf0 analogs:
morr_arr(i,j,k,MORRInd::xxlv) = Real(3.1484E6) - Real(2370.0) * T;
morr_arr(i,j,k,MORRInd::xxls) = Real(3.15E6) - Real(2370.0)*T + Real(0.3337E6);
// xlf = xxls - xxlv (computed at point of use)

// psat analog: replaced by helper function, not a constant
evs = calc_saturation_vapor_pressure(T, 0);  // liquid
eis = calc_saturation_vapor_pressure(T, 1);  // ice
```

### SolverChoice: control flags only
```cpp
// Only these come from SolverChoice in Define():
m_rdOcp   = sc.rdOcp;      // R/cp ratio
m_axis    = sc.ave_plane;
m_do_cond = !sc.use_shoc;
```

### WSM6 constant mapping table:
| Fortran arg | C++ resolution         | Source                  |
|-------------|------------------------|-------------------------|
| g           | CONST_GRAV             | ERF_Constants.H         |
| cpd         | Cp_d                   | ERF_Constants.H         |
| cpv         | Cp_v                   | ERF_Constants.H         |
| rd          | R_d                    | ERF_Constants.H         |
| rv          | R_v                    | ERF_Constants.H         |
| t0c         | tmelt (273.15)         | ERF_Constants.H         |
| ep2         | m_ep_2 = R_d/R_v       | class member (verify)   |
| ep1         | R_v/R_d - 1.0          | derived (no Morrison analog) |
| qmin        | m_qsmall               | class member            |
| xls         | xxls (in-kernel)       | computed per-cell       |
| xlv0        | xxlv (in-kernel)       | computed per-cell       |
| xlf0        | xxls - xxlv            | computed per-cell       |
| cliq        | m_cpw                  | class member            |
| den0        | m_rhosu                | class member            |
| denr        | m_rhow                 | class member            |
| psat        | calc_saturation_vapor_pressure() | helper fn  |
| cice        | [PENDING — no Morrison analog] |                |

---
## Path B Rules (continued) — promoted from WSM6 implementation notes

- Rule 5: WSM6-Specific Structural Differences [PARTIAL]
- Rule 6: 3D Species-Indexed Arrays → Expanded Enum Entries [NEW - WSM6]
- Rule 7: 1D (its:ite) Arrays → Kernel-Local Scalars [NEW - WSM6]
- Rule 8: Integer and Logical Arrays → Kernel-Local int/bool [NEW - WSM6]

## Rule 9: Separate Coefficient Init [COMPLETE — confirmed Morrison, WSM6, MYNN]
Every scheme has a separate initialization subroutine distinct from 
Advance(). It pre-computes derived coefficients from physical constants 
and stores them as class members. Two sub-patterns exist:
- WSM6 style: derives coefficients from constexpr parameters using 
gamma functions and geometric series — purely mathematical, no 
runtime input needed beyond variant flag
- MYNN style: receives all physical constants as explicit arguments 
to init, stores them as module save variables for run-time access
Coefficient computation is done in a private initialize_coeffs() 
helper, called at the end of Init() after MultiFab allocation. 
No new virtual method is needed — Init() is already the correct 
setup hook in the NullMoist interface.

## Rule 10: Runtime Variant Flags [COMPLETE — confirmed WSM6, MYNN]
Schemes have integer or enum flags that switch between physics regimes.
- WSM6: hail_opt (0=graupel, 1=hail) controls 5 class members
- MYNN: bl_mynn_closure, bl_mynn_mixlength, bl_mynn_edmf etc.
In C++: add bool or enum to SolverChoice, pull in Define(), branch in 
initialize_coeffs() to set affected class members.

## Rule 11: Two Valid AMReX GPU Translation Strategies [COMPLETE]
Strategy 1 — Enum + FArrayBox working array (Morrison, WSM6):
  Use when: scheme is column-local, logic is sequential, Fortran 
  argument list maps cleanly to per-cell arrays.
  Pattern: SchemeInd enum + single FArrayBox per tile + pack/unpack bridge.

Strategy 2 — AMREX_GPU_DEVICE struct with inline methods (MYNN):
  Use when: scheme has complex reusable sub-calculations, 
  stability functions, or closure coefficients that benefit from 
  encapsulation. Coefficients live as struct members with defaults,
  init_coeffs() derives secondary values.

WSM6 follows Strategy 1.

Rules 9, 10, 11 confirmed generic via MYNN 2.5 cross-reference. 
Rules 1 and 2 confirmed on Morrison, pending WSM6 as second 
validation point.

## Rule 12: Fortran Bridge as Ground Truth During Native Port [COMPLETE]

When porting Fortran physics to native AMReX C++, maintain the 
Fortran bridge simultaneously behind a compile flag. The native 
C++ kernel goes in the #else branch.

```cpp
#ifdef ERF_USE_<SCHEME>_FORT
    // Fortran bridge path — validated, ground truth
    mp_<scheme>_run_c(...);
#else
    // Native AMReX GPU C++ kernel — under development
    // Must reproduce Fortran bridge output to machine precision
    amrex::Abort("<Scheme> native C++ kernel not yet implemented");
#endif
```

Build flag convention: `ERF_ENABLE_<SCHEME>_FORT=ON`
Make.package pattern:
```makefile
ifeq ($(USE_<SCHEME>_FORT),TRUE)
    F90EXE_sources += ERF_module_<scheme>.F90
    F90EXE_sources += ERF_module_<scheme>_isohelper.F90
endif
```

The isohelper F90 file provides `bind(C)` wrappers that convert 
C-compatible argument types to Fortran types and call the original 
physics subroutine unchanged.

Acceptance criterion: `amrex_fcompare` between native and Fortran 
bridge plotfiles should show differences at or below machine epsilon 
for dynamics variables early in the squall line test case.
Precipitation accumulations may be exact zero in dry test cases.

Morrison precedent: partial bridge — Fortran init only via 
`morr_two_moment_init_c`, native C++ advance.
WSM6 precedent: full bridge — both init and run gated behind 
`ERF_USE_WSM6_FORT`.

---

## Rule 13: 2D Surface Accumulation — box2d Slab Pattern [COMPLETE]

Surface accumulation variables (rain_accum, snow_accum, graupel_accum) 
are 2D fields stored in 3D MultiFabs at k=klo. Inside the MFIter tile 
loop they require a separate 2D slab FArrayBox.

```cpp
// Inside MFIter loop, after getting box and fab_box:
Box box2d(fab_box);
box2d.makeSlab(2, 0);  // collapse k dimension to single slab

FArrayBox rainncv_fab(box2d, 1);
FArrayBox snowncv_fab(box2d, 1);
FArrayBox graupelncv_fab(box2d, 1);

// Initialize from persistent state at klo:
ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
    rainacc_arr(i,j,0)    = rain_arr(i,j,klo);
    snowacc_arr(i,j,0)    = snow_arr(i,j,klo);
    graupacc_arr(i,j,0)   = graup_arr(i,j,klo);
    rainncv_arr(i,j,0)    = Real(0.0);  // per-step increment
    snowncv_arr(i,j,0)    = Real(0.0);
    graupelncv_arr(i,j,0) = Real(0.0);
});

// ... run kernel ...

// Write back to persistent state:
ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
    rain_arr(i,j,klo)  = rainacc_arr(i,j,0);
    snow_arr(i,j,klo)  = snowacc_arr(i,j,0);
    graup_arr(i,j,klo) = graupacc_arr(i,j,0);
});
```

The `_ncv` (non-convective) variants are per-step increments passed 
to the physics kernel but not stored persistently in mic_fab_vars.

Fortran optional args (snow, snowncv, graupel, graupelncv) that were 
conditionally written via present() become unconditional in C++ — 
just always allocate and pass the slab.

---

## Rule 14: Init() Signature — Required Members [COMPLETE]

Every scheme's Init() requires these members beyond the basics:

```cpp
void Init(const amrex::MultiFab& cons_in,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance,
          std::unique_ptr<amrex::MultiFab>& z_phys_nd,  // terrain nodes
          std::unique_ptr<amrex::MultiFab>& detJ_cc)     // terrain Jacobian
```

Required class members set in Init():
```cpp
amrex::MultiFab* m_z_phys_nd{nullptr};  // terrain node heights
amrex::MultiFab* m_detJ_cc{nullptr};    // cell-centered Jacobian
amrex::Geometry  m_geom;
amrex::Real      dt{0.0};
int              nlev{0}, zlo{0}, zhi{0};
```

`m_z_phys_nd` is used to compute per-cell `delz` from terrain:
```cpp
ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
    delz_arr(i,j,k) = (z_arr) ? 
        Real(0.25) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                     + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                     + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                     + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) )
        : dz_val;  // fallback to uniform spacing
});
```

This replaces Fortran's `delz` intent(in) argument.

---

## Development Workflow Notes [Both]

### GNUmake vs CMake
Use GNUmake for iterative development — no reconfigure step on
file changes, faster incremental builds.
Use CMake for CI, installation, and cross-platform builds.

GNUmake flags:
  USE_<SCHEME>_FORT=TRUE   — enable Fortran bridge objects
  USE_NETCDF=FALSE         — disable NetCDF (speeds up dev builds)
  
Build from Exec/ or Exec/<TestCase>/ directory.
Executable appears as ERF3d.gnu.TEST.MPI.ex (or similar).

### Acceptance Criterion for Path B
Run amrex_fcompare between Path A (Fortran bridge) and Path B
(native C++) plotfiles at the same timestep:
  amrex_fcompare plt_<scheme>_fort00120 plt_<scheme>_cpp00120

Expected: dynamics variables differ by <= machine epsilon.
Precipitation accumulations may be exact zero in dry test cases.
Any larger difference indicates a physics transcription error.

WSM6 validated values (step 120, squall line 2D):
  density:   9.1e-15 absolute, 8.1e-15 relative
  rhotheta:  1.3e-12 absolute, 3.9e-15 relative
  rain/snow/graupel accum: 0 (dry test, no precipitation)

---
---

## Path B — Rule 1 Addendum: Kernel Working Array Enum Placement

MORRInd lives in ERF_AdvanceMorrison.cpp (translation-unit scope).
WSM6Ind lives in ERF_WSM6.H (header scope). Both compile correctly.

Prefer .H placement when the enum needs to be visible to other
translation units or for forward reference. .cpp placement is
acceptable when the enum is truly local to the advance function.
Pick based on visibility needs, not convention.

---

## Path B — Rule 7 Extension: Fortran Scalar Locals [COMPLETE - WSM6]

Rule 7 previously covered only dimension(its:ite) 1D arrays.
WSM6 also has ~40 Fortran subroutine-scope scalars that are
undimensioned — these follow the same rule.

All Fortran subroutine-scope scalars become kernel-local declarations
inside the ParallelFor lambda. Do not add any to the SchemeInd enum.

Categories:
  Physics intermediates: cpmcal, xlcal, diffus, viscos, xka,
    venfac, conden, diffac, coeres, supsat, xmi, eacrs, satdt,
    qimax, diameter, xni0, roqi0, fallsum, fallsum_qsi, fallsum_qg,
    vt2i, vt2r, vt2s, vt2g, acrfac, egs, egi, xlwork2, factor,
    source, value, xlf, pfrzdtc, pfrzdtr, supice, alpha2,
    delta2, delta3, vt2ave, holdc, holdci, pvt, xlcal, supcol,
    supcolt, qdt, holdrr, holdrs, holdrg, temp
  Loop control: dtcld, loops, loop, ifsat, n, idim, kdim, mstepmax
  Generic temporaries: x, y, z, a, b, c, d, e
  Sub-step control: moved to Rule 15 below
  fpvs inlining temporaries: dldti, xb, xai, tr, xbi, xa,
    hvap, cvap, hsub, dldt, ttp
    (these replace the Fortran fpvs() saturation vapor pressure
    function — inline the polynomial directly rather than calling
    a helper function)

Declare all of these as Real or int or bool at the top of the
ParallelFor lambda before any physics computation.

---

## Path B — Rule 15: Sub-Stepping Loop — Outer Loop Pattern [COMPLETE - WSM6]

WRF schemes that use dtcld sub-stepping compute a global sub-step
count (loops) from the physics timestep before any column processing.
Because loops is uniform across all columns, the sub-step becomes
an outer C++ for loop wrapping multiple ParallelFor calls — NOT an
inner loop inside a single ParallelFor.

    int loops = std::max(static_cast<int>(std::round(dt / dtcldcr)), 1);
    Real dtcld = Real(dt) / loops;

    for (int loop = 1; loop <= loops; ++loop) {
        // Each Fortran do-k/do-i block inside do loop = a separate ParallelFor
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // initialization within loop iteration
        });
        amrex::Gpu::synchronize();  // required between dependent kernels
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // warm rain processes
        });
        amrex::Gpu::synchronize();
        // etc.
    }

Gpu::synchronize() is required between ParallelFor calls within one
loop iteration when later kernels read values written by earlier ones.

Contrast with Rule 15 (original draft): a single ParallelFor with an
inner for loop is only correct when the entire sub-step fits in one
kernel with no data dependencies requiring synchronization between
physics blocks. WSM6 has such dependencies (slope parameters computed
before process rates) so the outer-loop pattern is required.

mstep(its:ite) and flgcld(its:ite) are initialized once per loop
iteration and are never modified from their initial values (mstep=1,
flgcld=.true.) anywhere in mp_wsm6_run. They are vestigial
sub-step scaffolding. In C++ both become kernel-local scalars:
    int mstep = 1;
    bool flgcld = true;
Do not allocate a box2d FArrayBox for these variables.

## Path B — Rule 16: Fortran Statement Functions → GPU Free Functions [COMPLETE - WSM6]

Fortran statement functions (single-line inline expressions defined
in the declaration section before the first executable statement)
cannot become C++ lambdas inside AMREX_GPU_DEVICE kernels — nested
GPU lambda capture is not portable across CUDA/HIP/SYCL backends.

The established ERF pattern confirmed from Morrison (ERF_AdvanceMorrison.cpp
lines 254, 432, 444) is AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE.
Use this for all statement function equivalents called from ParallelFor.
Do not use constexpr for runtime-argument functions — constexpr is
reserved for compile-time mathematical constants (e.g. Morrison line 159:
constexpr Real xxx = Real(0.9189385332046727417803297)).

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real wsm6_cpmcal(Real x, Real qmin_arg, Real cpd_arg, Real cpv_arg) {
        return cpd_arg * (Real(1.0) - amrex::max(x, qmin_arg))
             + amrex::max(x, qmin_arg) * cpv_arg;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real wsm6_xlcal(Real x, Real xlv0_arg, Real xlv1_arg, Real t0c_arg) {
        return xlv0_arg - xlv1_arg * (x - t0c_arg);
    }

Annotation summary:
  Functions called from ParallelFor: AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Compile-time mathematical constants: constexpr Real
  Host-only helpers: AMREX_FORCE_INLINE (no GPU annotation needed)

WSM6 statement functions requiring this treatment:
  cpmcal, xlcal, diffus, viscos, xka, diffac, venfac, conden

Name prefix (wsm6_) avoids collision with any ERF global symbols.
All constants that the Fortran statement function captured from
module scope must become explicit function arguments.

## Path B — Rule 17: Loop-Invariant Scalars — CPU Precomputation [COMPLETE - WSM6]

Some Fortran schemes compute scalars outside the do-k/do-i loops
but inside the outer sub-step loop. These are loop-iteration-invariant
but not globally constant — they depend on constants available at
compile time or on the timestep.

In C++: compute these as plain Real scalars on the CPU before the
corresponding ParallelFor call. They are automatically captured by
value in the lambda.

WSM6 fpvs precomputation (appears twice per loop iteration —
once before the first qsat block, once before the second):

    // CPU side — before ParallelFor
    Real ttp   = static_cast<Real>(t0c) + Real(0.01);
    Real dldt  = cpv - cliq;
    Real xa    = -dldt / rv;
    Real xb    = xa + xlv0 / (rv * ttp);
    Real dldti = cpv - cice;
    Real xai   = -dldti / rv;
    Real xbi   = xai + xls / (rv * ttp);

    // GPU side — captured by value, used inside ParallelFor
    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        Real tr = ttp / wsm6_arr(i,j,k,WSM6Ind::t);
        Real qsat_liq = psat * std::exp(std::log(tr)*xa) * std::exp(xb*(Real(1.0)-tr));
        // ...
    });

Do not put these scalars in WSM6Ind. Do not recompute them inside
the ParallelFor. The compiler will place them in GPU constant/register
space automatically via lambda capture.

---

## Path B — Rule 5d Extension: libmassv Elimination [COMPLETE - WSM6]

Rule 5d previously noted that vrec/vsqrt become std::pow/std::sqrt.
More precisely:

dvec1 and tvec1 in WSM6 are vectorization-hint arrays used only
as arguments to vrec() and vsqrt() from module_libmassv. In the
GPU kernel both the arrays and the function calls disappear entirely.

Replace at point of use:
  vrec(y, x, n)  →  y = Real(1.0) / x
  vsqrt(y, x, n) →  y = std::sqrt(x)

Do not add dvec1 or tvec1 to WSM6Ind.
Do not include or call any module_libmassv function from the
#else native C++ kernel path — the isohelper and Fortran objects
are only compiled when USE_WSM6_FORT=TRUE.
---

## Path B — Rule 18: Cell-Local Internal Subroutines → Device Free Functions [COMPLETE - WSM6, Morrison]

Fortran internal subroutines with no vertical stencil (no k+1 or k-1
reads or writes) map directly to AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
free functions at file scope, callable from inside ParallelFor(box, ...).

Cross-scheme validation:
  Morrison: wrf_gamma, gamma_function, calc_saturation_vapor_pressure
    (ERF_AdvanceMorrison.cpp lines 254, 432, 444)
  WSM6: wsm6_slope_wsm6, wsm6_slope_rain, wsm6_slope_snow,
    wsm6_slope_graup, wsm6_lamdar, wsm6_lamdas, wsm6_lamdag
    (slope_rain/snow/graup confirmed cell-local, no stencil)
    (ERF_AdvanceWSM6.cpp — WSM6 Path B implementation)

Pattern:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real wsm6_lamdar(Real x, Real y, Real pidn0r_arg) {
        return std::sqrt(std::sqrt(pidn0r_arg / (x * y)));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void wsm6_slope_wsm6(Real qr, Real qs, Real qg, Real den,
                          Real denfac, Real t, Real n0sfac,
                          Real& rslope_r, Real& rslopeb_r, ...) {
        // all computation is cell-local, no array indexing
    }

Annotation decision tree:
  Called from ParallelFor(box,...) GPU kernel:
    → AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE (free function)
  Called from ParallelFor(box2d,...) column-sweep GPU kernel:
    → AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE (free function)
  Called only from host-side MFIter setup code or struct methods:
    → AMREX_FORCE_INLINE only, no GPU annotation needed
    (MYNN pattern: ERF_MYNNStruct.H lines 18, 95, 104, 112)
  Compile-time mathematical constant, not a function:
    → constexpr Real (Morrison line 159 pattern)

Internal statement functions of the Fortran subroutine (lamdar,
lamdas, lamdag inside slope_wsm6) get the same annotation as the
subroutine that contains them — they become separate free functions,
not nested lambdas.

Name all WSM6 device functions with wsm6_ prefix to avoid collision
with any ERF global symbols.

---

## Path B — Rule 19: Column-Serial Algorithms → box2d Column-Sweep Kernel [WSM6-first]

Fortran subroutines with serial vertical dependencies (cumulative
sums, order-dependent sweeps, convergence iterations, sequential
vertical searches) cannot be ParallelFor(box,...) kernels.

Diagnostic questions:
  1. Does any loop read k+1 or k-1 of an array it also writes?
     → serial vertical dependency
  2. Does any loop run top-to-bottom or bottom-to-top with a
     data-dependent update (wi(k) uses wi(k+1) just written)?
     → order-dependent, cannot be parallelized over k
  3. Is there a convergence iteration (go to / while) over the
     full column?
     → serial within the column

If any diagnostic is YES: use box2d column-sweep pattern.

AMReX translation — pre-allocated FArrayBox working arrays:
The ERF convention is to pre-allocate FArrayBoxes outside the
ParallelFor (in the MFIter loop body) and pass Array4 accessors
into the lambda. This is already established in ERF_AdvanceWSM6.cpp
for delz_fab, rainncv_fab, rainacc_fab etc.

For nislfv working arrays, follow the same pattern:

    // Before the box2d ParallelFor — allocate working arrays
    // sized to fab_box (km = fab_box.length(2))
    FArrayBox nislfv_work_fab(fab_box, NislfvWorkInd::NumInds);
    auto const& nw = nislfv_work_fab.array();

    ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
        // Each thread accesses its (i,j) column slice
        // No stack allocation needed — km is runtime, handled by Array4
        for (int k = klo; k <= khi; ++k) {
            nw(i,j,k,NislfvWorkInd::dz)  = delz_arr(i,j,k);
            nw(i,j,k,NislfvWorkInd::ww)  = wsm6_arr(i,j,k,WSM6Ind::work1_1);
            // ... fill remaining column arrays ...
        }
        // serial vertical sweep — all within this (i,j) thread
        // cumulative sums, stencil passes, interpolation
        // write results back to wsm6_arr
    });

GpuArray<Real, N> is used in ERF only for truly compile-time-fixed
small sizes (N = AMREX_SPACEDIM = 3, or NBCVAR_max for boundary
conditions). It is not appropriate for runtime column-length arrays.
Do not use GpuArray for nislfv working arrays.

Working arrays for nislfv_rain_plm (per column, km or km+1 sized):
  dz, ww, qq, wd, wa, was, den, denfac, tk, qn, qr, tmp, tmp1, tmp2, tmp3
  wi[km+1], zi[km+1], za[km+1], dza[km+1], qa[km+1], qmi[km+1], qpi[km+1]

These can be added as a local NislfvWorkInd enum (separate from WSM6Ind)
or folded into a dedicated FArrayBox with enough components.
The km+1 interface arrays require one extra level — allocate
the FArrayBox with domain box grown by 1 in z, or use a
separate interface FArrayBox.

// ERF_WSM6.H line 248 defines: constexpr int WSM6_MAX_LEVELS = 256;
// This was added in the WSM6 stub and is available for use.
// FArrayBox pre-allocation (Rule 19 primary approach) is preferred
// over stack arrays, but WSM6_MAX_LEVELS exists if needed.

WSM6 nislfv go-to convergence loop with iter=1:
The Fortran go to 100 with iter=1 runs exactly one convergence pass
(n starts at 1, 1<=1 is true, increments to 2, 2<=1 is false, exits).
In C++ this becomes a plain sequential block, not a loop:

    // Single convergence pass (iter=1 at call site)
    wsm6_slope_rain(qr, den, denfac, tk, tmp, tmp1, tmp2, tmp3, wa, km);
    for (int k = 0; k < km; ++k) {
        ww[k] = Real(0.5) * (wd[k] + wa[k]);
    }
    // No goto needed — iter=1 means exactly one pass

Alternative design — Morrison tendency-based sedimentation:
Morrison avoids column-serial sedimentation entirely by computing
fall speed parameters per cell (stored in MORRInd as ain, arn, asn,
acn, agn) and applying them as tendencies. This is simpler for GPU
and eliminates the need for Rule 19 entirely. When designing a new
scheme, prefer tendency-based sedimentation over semi-Lagrangian
PLM if GPU performance is a priority.

WSM6 nislfv subroutines requiring this treatment:
  nislfv_rain_plm  (~244 lines): rain sedimentation
  nislfv_rain_plm6 (~276 lines): combined snow+graupel sedimentation

Both have column-independent i_loop — columns are safe to parallelize
over. Serial dependency is only within each column's vertical sweep.

Validation status: WSM6-first. No prior ERF scheme uses this pattern.
Second validation pending (Thompson or NSSL if ported).

## Path B — Rule 20: Runtime Constants → [=]-Captured Locals in ParallelFor [Morrison + MYNN + WSM6 confirmed]

ERF-wide standard (confirmed Morrison, MYNN, WSM6):
GPU ParallelFor lambdas use [=] capture of local variables only.
Class member variables (true m_ members accessed via 'this') must
never appear directly inside AMREX_GPU_DEVICE lambdas.

MYNN (ERF_ComputeDiffusivityMYNN25.cpp line 76):
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    { ... }
Constants are local variables or Array4 handles declared in the
surrounding MFIter scope — no class member access inside kernel.

Two variants for how locals are populated before capture:

VARIANT A — cheap constants (Morrison pattern):
Declare and initialize directly as locals inside Advance().
m_ prefix is naming convention only — these are not class members.
    Real m_pi  = Real(3.14159...);
    Real m_rhow = Real(997.0);
    // ... ~300 lines of initialization ...
    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        // m_pi, m_rhow captured by value — GPU safe
    });
Use when: initialization is cheap arithmetic, no option branching.

VARIANT B — expensive constants (WSM6 pattern):
Store as true class members (m_ prefix), initialized once in
initialize_coeffs() called from Init(). Copy to l_ locals inside
Advance() before each ParallelFor.
    // WSM6::Advance(), before ParallelFor:
    auto l_pvtr        = m_pvtr;
    auto l_pvts        = m_pvts;
    auto l_pvtg        = m_pvtg;
    auto l_pidn0r      = m_pidn0r;
    auto l_pidn0s      = m_pidn0s;
    auto l_pidn0g      = m_pidn0g;
    auto l_rslopermax  = m_rslopermax;
    auto l_rslopesmax  = m_rslopesmax;
    auto l_rslopegmax  = m_rslopegmax;
    auto l_rsloperbmax = m_rsloperbmax;
    auto l_rslopesbmax = m_rslopesbmax;
    auto l_rslopegbmax = m_rslopegbmax;
    auto l_rsloper2max = m_rsloper2max;
    auto l_rslopes2max = m_rslopes2max;
    auto l_rslopeg2max = m_rslopeg2max;
    auto l_rsloper3max = m_rsloper3max;
    auto l_rslopes3max = m_rslopes3max;
    auto l_rslopeg3max = m_rslopeg3max;
    auto l_bvtg        = m_bvtg;
    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        // l_ locals captured by value — GPU safe
    });
Use when: initialization requires gamma functions, option-dependent
branching (hail_opt), or is otherwise expensive to recompute each
Advance() call.

Static constexpr class members (WSM6::bvtr, WSM6::bvts,
WSM6::qcrmin etc. — ERF_WSM6.H lines 51-73) are accessible
directly inside GPU lambdas without local copies — they have no
'this' dependency and are inlined by the compiler.

Device free functions (Rule 16, Rule 18) receive all runtime
constants as plain Real arguments:
    wsm6_slope_rain_cell(qr, den, denfac,
                          l_pidn0r, l_rslopermax, l_rsloperbmax,
                          l_rsloper2max, l_rsloper3max,
                          WSM6::bvtr, WSM6::avtr, WSM6::qcrmin,
                          rslope, rslopeb, rslope2, rslope3, vt);

Future scheme ports: prefer Variant A unless initialization cost
justifies once-per-Init computation.

---

## Rule 21: Exact Algorithm Porting for Init-Time Special Functions [WSM6 + Morrison - validation critical]

When a Fortran scheme uses a custom implementation of a mathematical
special function (gamma, erf, bessel, etc.) rather than an intrinsic,
the C++ port must replicate the exact algorithm — not substitute
std::tgamma, std::erf, or other library functions.

Reason: std:: library implementations and custom infinite-product or
series implementations agree to several significant digits but not
exactly. Any difference propagates into every coefficient derived from
that function (pvtr, pidn0r, pacrr, precr1/2, and all equivalents).
During one-to-one Fortran/C++ validation runs, these differences
appear as non-roundoff discrepancies that are indistinguishable from
logic bugs without careful diagnosis.

WSM6 case: rgmma(x) in ERF_module_mp_wsm6.F90 uses the Weierstrass
infinite product form to compute Γ(x). The C++ initialize_coeffs()
must port this exact loop, not call std::tgamma(x).
Verified: rgmma(4.8) = 17.837825 (Fortran comment), Γ(4.8) = 17.838.

Morrison case: gamma_function() and wrf_gamma() in
ERF_AdvanceMorrison.cpp — check whether these are exact ports of
Fortran intrinsics or custom algorithms. If custom, same rule applies.

Pattern — port the Fortran algorithm directly:

    // WSM6 rgmma: Weierstrass infinite product for Γ(x)
    // Mirrors ERF_module_mp_wsm6.F90 lines 1472-1495 exactly.
    // Do NOT replace with std::tgamma — validation runs require
    // bit-compatible output with the Fortran reference.
    auto rgmma = [](Real x) -> Real {
        // ... exact Weierstrass port, see ERF_InitWSM6.cpp ...
    };

Scope: applies to any special function computed by a custom loop
or series in the Fortran source. Does not apply to simple arithmetic
that happens to match a library function exactly.

WSM6 vs Morrison gamma implementations are NOT interchangeable:
  WSM6 rgmma: Weierstrass infinite product, 10000 iterations,
    CPU-only (used only in initialize_coeffs(), not GPU kernels)
  Morrison wrf_gamma: rational minimax approximation, ~20 operations,
    intended for GPU but currently missing AMREX_GPU_HOST_DEVICE
    annotation on wrf_gamma() itself (line 255) — GPU compilation
    bug, consistent with Morrison C++ being marked testing-only.

WSM6 does not need a GPU-callable gamma function: mp_wsm6_run uses
only precomputed coefficients from initialize_coeffs(), never calls
rgmma at runtime. Do not add GPU annotation to the rgmma lambda.

Cross-scheme validation note: WSM6 Fortran rgmma(4.8) and Morrison
wrf_gamma(4.8) will agree to ~6 significant digits but not exactly.
Do not use Morrison's wrf_gamma to initialize WSM6 coefficients —
this would break one-to-one Fortran/C++ validation for WSM6.

---

## Rule 22: Name New Enum Entries by Convention Survey Before Writing [process rule]

Before naming any new WSM6Ind or MicVar_WSM6 entry, grep all
existing C++ microphysics schemes in Source/Microphysics/ for
prior art on the same physical quantity. Do not invent names that
diverge from established codebase conventions.

Survey command:
    grep -rn "TARGET_NAME" /home/jmsexton/codes/ERF/Source/Microphysics/ \
      | grep -v ".F90\|.f90\|_fort\|FORT"

Naming priority order:
1. Match the WRF Fortran source name exactly where unambiguous
   (e.g. qr, qs, qg, qc, qi, t, den, denfac, xni, n0sfac)
2. For Fortran arrays with a third dimension index (e.g. qsat(i,k,1),
   qsat(i,k,2)), use the shortest unambiguous suffix that matches
   WRF convention — survey other schemes before assigning
3. Match Morrison MORRInd names for shared quantities
   (rain_accum, snow_accum, graup_accum already established)
4. If no prior art exists in the ERF codebase, WSM6 sets the
   convention — document the decision and rationale here

Work array naming: Fortran work1(i,k,3) that is dual-purpose
across a routine must become two separate named FABs in C++,
one per distinct use. Name each FAB after its physical role,
not after the Fortran scratch name.

Integer and logical per-column arrays (mstep, numdt, flgcld in
WSM6): these are different types from Real FABs and cannot go
in the main WSM6Ind enum. Use separate IArrayBox (integer) or
treat as local scalars inside ParallelFor if column-independence
allows. Document the decision per array.

### Rule 22 Addendum: Confirmed Naming Decisions for WSM6 Two-Phase Arrays

Survey of all C++ microphysics in Source/Microphysics/ confirmed
Morrison establishes the ERF-wide convention:
  w suffix = over liquid water  (qsatw, qvs in Morrison)
  i suffix = over ice           (qsati, qvi in Morrison)

WSM6 Fortran third-dimension index -> C++ name:
  qsat(i,k,1) -> qsatw    (saturation mixing ratio over liquid)
  qsat(i,k,2) -> qsati    (saturation mixing ratio over ice)
  rh(i,k,1)   -> rhw      (relative humidity over liquid)
  rh(i,k,2)   -> rhi      (relative humidity over ice)

Fortran work1(i,k,3) is dual-purpose across mp_wsm6_run:
  First use  (after slope_wsm6 call 1, lines 617-618):
    fall speed for rain/snow/graupel — kept as work1 with
    component index matching qrs_tmp(1/2/3)
  Second use (after slope_wsm6 call 2, lines 843-856):
    thermodynamic diffusion terms, overwriting work1 —
    in C++ these become two separate FABs:
      workdiffw   (liquid diffusion, Fortran work1(i,k,1))
      workdiffi   (ice diffusion,    Fortran work1(i,k,2))

rslope/rslopeb/rslope2/rslope3: dimensioned (i,k,3) with
rain(1)/snow(2)/graupel(3) — WSM6Ind entries use _r/s/g suffix,
consistent with Morrison lamr/lams/lamg convention.

mstep, numdt: integer per-column arrays — not in WSM6Ind Real enum.
Use local IArrayBox inside Advance() MFIter scope.
flgcld: logical per-column — treat as local IArrayBox (0/1) or
eliminate if column-independence allows scalar per ParallelFor thread.

### Rule 22 Addendum 2: WSM6Ind Final Confirmed State

After survey and correction, WSM6Ind naming is:

Two-component arrays (Fortran third dim = 2):
  rhw, rhi       — relative humidity over liquid/ice
                   Fortran: rh(i,k,1), rh(i,k,2)
  qsatw, qsati   — saturation mixing ratio over liquid/ice
                   Fortran: qsat(i,k,1), qsat(i,k,2)
  Matches Morrison convention: qsatw/qsati, qvs/qvi throughout
  ERF_Morrison.H and ERF_Morrison_Cloud/Precip.cpp

Three-component arrays (Fortran third dim = 3, 1=rain 2=snow 3=graupel):
  rslope_r/s/g, rslope2_r/s/g, rslope3_r/s/g, rslopeb_r/s/g
  qrs_tmp_r/s/g, falk_r/s/g, fall_r/s/g
  work1_r/s/g    — fall speeds from slope_wsm6 (first call)
  Suffix convention consistent with Morrison lamr/lams/lamg

Dual-use work1 split:
  work1_r/s/g   — fall speeds (Fortran work1 after first slope_wsm6)
  workdiffw     — liquid diffusion denominator (Fortran work1(i,k,1) after second slope_wsm6)
  workdiffi     — ice diffusion denominator    (Fortran work1(i,k,2) after second slope_wsm6)
  Fortran work1 is reused for two distinct purposes across the
  routine — C++ requires separate named FABs for clarity and
  to avoid overwriting live data.

Not in WSM6Ind (different types, handled separately):
  mstep  — integer per-column, use local IArrayBox in Advance()
  numdt  — integer per-column, use local IArrayBox in Advance()
  flgcld — logical per-column, use local IArrayBox (0/1) in Advance()
---

## Rule 23: Dead and Simplified Fortran Constructs — Omit or Inline in C++ Port

When porting Fortran to C++, identify variables that are initialized
but never read (dead code) or variables whose value is invariant
(always the same constant). Do not create FABs or IArrayBoxes for
these — omit or inline them.

WSM6 confirmed cases:

flgcld(i): logical array, set to .true. at mp_wsm6_run line 495,
never read anywhere in the 1469-line routine. Omit entirely in C++.
Do not create a FAB or IArrayBox for it.

mstep(i): integer array, initialized to 1 at line 494, used as
divisor at lines 681/696 (psmlt/pgmlt bounds), never modified.
The code that would update it (lines 498-502) is commented out.
Replace with literal 1 in C++ — no IArrayBox needed.

work1 third use (line 1429): Fortran reuses work1(i,k,1) as a
scratch variable for the condensation calculation immediately
before using it to compute pcond(i,k). In C++ this is a local
Real inside the ParallelFor lambda:
    Real workcond = wsm6_conden(t,q,qsatw,xl,cpm,...);
    Real work2loc = qc + workcond;
    pcond = ...;
Do not add a workcond FAB entry to WSM6Ind.

qsat computed twice inside do loop: lines 530-549 (first fpvs
expansion) and lines 1390-1420 (second fpvs expansion after state
update). Both write to the same qsatw/qsati FABs. In C++ these
are two separate ParallelFor passes to the same Array4 targets.
This is correct — the second pass overwrites with updated values.

General rule: before creating any FAB or IArrayBox for a Fortran
local variable, verify it is (a) read after being written and
(b) needed across multiple loop nests. Pure scratch within a
single (i,k) iteration becomes a local Real inside ParallelFor.

---

## Rule 24: nislfv Semi-Lagrangian Column Kernels

nislfv_rain_plm and nislfv_rain_plm6 are column-sequential
semi-Lagrangian sedimentation routines. Both have an outer
i_loop over columns with no cross-column dependencies.
All working arrays (dz, ww, qq, wi, zi, za, qn, etc.) are
column-local stack arrays dimensioned km.

C++ porting strategy:
- Port as ParallelFor over column index i only (not i,k)
- Local arrays become stack arrays inside the lambda,
  sized WSM6_MAX_LEVELS (already defined = 256 in ERF_WSM6.H)
- Fortran goto 100 with iter=1 becomes a single-pass loop
  (iter is always passed as 1 in all call sites in mp_wsm6_run)
- nislfv_rain_plm6 calls slope_snow and slope_graup internally
  on local 1D arrays — these must be AMREX_GPU_HOST_DEVICE
  device functions callable from inside the column kernel

nislfv_rain_plm6 processes snow (rql) and graupel (rql2) in
a single column pass using ist_loop (ist=1 snow, ist=2 graupel)
with shared terminal velocity wa averaged from both species.
precip1 receives snow flux, precip2 receives graupel flux.

Call sites in mp_wsm6_run:
  nislfv_rain_plm(idim,kdim,...,workr,denqrs1,delqrs1,dtcld,1,1)
    → rain sedimentation, iter=1
  nislfv_rain_plm6(idim,kdim,...,worka,denqrs2,denqrs3,
                   delqrs2,delqrs3,dtcld,1,1)
    → snow+graupel sedimentation, iter=1
  nislfv_rain_plm(idim,kdim,...,work1c,denqci,delqi,dtcld,1,0)
    → ice sedimentation, iter=0 (pure forward, no mean-wind iteration)

iter=0 path skips the mean-wind correction entirely.
iter=1 path runs one mean-wind correction iteration.

---

## Rule 25: Phase Structure for Large Routine Ports

When porting a Fortran routine with more than ~300 lines of physics
(WSM6 mp_wsm6_run is 1255 lines), the naive two-phase approach
(Fortran bridge → full C++ replacement) fails because:
- The full C++ kernel cannot be written, reviewed, and validated
  in a single context pass
- A half-implemented kernel produces no useful output for comparison
- Bugs in process A contaminate validation of process B

Required phase structure for large routine ports:

Phase 1 — Fortran bridge
  Compile and run guard: #ifdef ERF_USE_WSM6_FORT
  Establish working reference output.
  All subsequent phases must bit-reproduce or scientifically
  match this reference before proceeding.

Phase 2 — Deep analysis (read-only)
  Read full Fortran source. Build WSM6Ind. Build device functions.
  Write skill doc findings. Do not write any kernel code in this phase.
  Output: complete process inventory, naming decisions, type decisions,
  loop structure decisions, all documented in skill doc.

Phase 3 — Per-process diagnostic instrumentation
  Add write/print statements to the Fortran path at each physics
  process boundary. Target: one column, one timestep, stdout.
  Process boundaries in WSM6 (in order):
    denfac, qsat/rh computation
    slope_wsm6 call 1 → rslope/rslopeb/rslope2/rslope3/work1_r/s/g
    nislfv rain → qr, fall_r, delqrs1
    nislfv snow+graupel → qs, qg, fall_s, fall_g, delqrs2, delqrs3
    slope_wsm6 call 2 (after pimlt/pgfrz) → rslope update
    workdiffw, workdiffi, work2
    warm rain: praut, pracw, prevp
    cold rain T<T0: praci, piacr, psaci, pgaci, psacw, pgacw,
               paacw, pracs, psacr, pgacr, pgacs
    cold rain T<T0: pseml, pgeml, pidep, psdep, pgdep, pigen,
               psaut, pgaut
    T>=T0: psevp, pgevp
    mass conservation check and state update: q,qc,qi,qr,qs,qg,t
    second qsat computation
    pcond, final qc/t update

  Each diagnostic point prints: process name, i, k, value.
  This gives a per-process ground truth for Phase 4.

Phase 4 — Incremental C++ implementation with per-process validation
  Implement one process block at a time in the C++ #else path.
  For each process:
    a. Implement the C++ equivalent of that block
    b. Print the same diagnostic as Phase 3 from the C++ path
    c. Run both paths, diff the output for that process
    d. Green-light before proceeding to next process
  Order follows data dependencies — implement in Fortran
  execution order. Do not skip ahead.

Phase 5 — Full native C++ kernel
  Only reached when all Phase 4 process blocks are green-lit.
  Assemble into complete #else block.
  Remove diagnostic prints.
  Run full fcompare against Fortran reference.

## Rule 26: Context Budget Awareness for Large Ports

The WSM6 full #else kernel requires approximately:
  - nislfv_rain_plm column kernel: ~120 lines
  - nislfv_rain_plm6 column kernel: ~150 lines
  - slope_wsm6 ParallelFor: ~40 lines
  - denfac/qsat/rh init ParallelFor: ~30 lines
  - process rate ParallelFor (25 processes): ~300 lines
  - mass conservation + state update ParallelFor: ~120 lines
  - pcond + final update ParallelFor: ~30 lines
  Total: ~800 lines

800 lines cannot be written, reviewed, and debugged in a single
context pass alongside the Fortran source (1469 lines) plus
ERF_AdvanceWSM6.cpp (301 lines) plus skill doc (1083 lines).

Rule: Never attempt to write more than one Phase 4 process group
per context pass. One process group = one ParallelFor or one
column kernel. Write it, verify it compiles, then stop.
---

## Rule 27: Complete Ordered Process Inventory for mp_wsm6_run

The following is the canonical execution order of physics blocks
in mp_wsm6_run (Fortran lines shown for reference). Phase 4 must
implement in this order — no reordering, no skipping.

OUTSIDE the loops/dtcld minor timestep loop (execute once):
  A. Clamp qc/qr/qi/qs/qg >= 0          (lines 441-449)
  B. Compute cpm, xl                     (lines 455-460)
  C. Copy den, delz to temporaries       (lines 461-466)
  D. Initialize rainncv, sr, snowncv,
     graupelncv, tstepsnow, tstepgraup   (lines 471-479)
  E. Compute loops, dtcld                (lines 484-486)

INSIDE the loops/dtcld minor timestep loop (repeat loops times):

  Group 1 — Column setup (lines 493-549):
    G1a. mstep=1, flgcld=true (omit both — dead/invariant, Rule 23)
    G1b. Compute denfac via vrec/vsqrt   (lines 503-515)
    G1c. Compute qsatw, qsati, rhw, rhi  (lines 517-549)
         (first fpvs inline expansion)

  Group 2 — Initialize all process rate arrays to 0.0 (lines 555-594)
    prevp,psdep,pgdep,praut,psaut,pgaut,pracw,praci,piacr,
    psaci,psacw,pracs,psacr,pgacw,paacw,pgaci,pgacr,pgacs,
    pigen,pidep,pcond,psmlt,pgmlt,pseml,pgeml,psevp,pgevp,
    falk_r/s/g, fall_r/s/g, fallc, falkc, xni

  Group 3 — Ice number concentration xni (lines 598-604)

  Group 4 — First slope_wsm6 call (lines 610-618)
    Pack qrs_tmp_r/s/g from qr/qs/qg
    Call slope_wsm6 → rslope/rslopeb/rslope2/rslope3/work1_r/s/g

  Group 5 — Sedimentation (lines 620-653)
    G5a. workr, worka, denqrs1/2/3       (lines 620-635)
    G5b. nislfv_rain_plm  (rain)         (line 636-637)
    G5c. nislfv_rain_plm6 (snow+graupel) (line 638-639)
    G5d. Update qr/qs/qg from denqrs1/2/3, compute fall_r/s/g (lines 640-654)
    G5e. Slab fall at k=kts from delqrs1/2/3 (lines 650-654)

  Group 6 — Second slope_wsm6 call (lines 655-663)
    Pack qrs_tmp again, call slope_wsm6 (slope params updated
    after sedimentation has moved mass)

  Group 7 — Melting (T>T0 only): psmlt, pgmlt (lines 665-704)
    work2(i,k) = venfac(p,t,den) computed here
    psmlt and pgmlt applied immediately to qs/qg/qr/t
    (inline application, not deferred to mass conservation block)

  Group 8 — Ice crystal fallout via nislfv_rain_plm (lines 706-736)
    Vice computation → work1c
    denqci pack, nislfv_rain_plm call (iter=0)
    Update qi from denqci, fallc from delqi

  Group 9 — Surface precipitation accumulation (lines 738-772)
    fallsum, fallsum_qsi, fallsum_qg → rainncv, snowncv,
    graupelncv, tstepsnow, tstepgraup, sr

  Group 10 — Instantaneous phase changes (lines 774-830)
    G10a. pimlt: T>T0, I->C           (lines 783-787)
    G10b. pihmf: T<-40C, C->I         (lines 792-796)
    G10c. pihtf: T0>T>-40C, C->I      (lines 801-810)
    G10d. pgfrz: T<T0, R->G           (lines 815-829)
    (All four applied immediately, inline temperature updates)

  Group 11 — Third slope_wsm6 call (lines 836-844)
    Pack qrs_tmp again after phase changes

  Group 12 — n0sfac, work2 (venfac) (lines 845-870, inside k loop)

  Group 13 — Warm-rain processes (T>=T0 and T<T0 shared) (lines ~870-1062)
    praut: autoconversion C->R
    pracw: accretion C by R
    paacw: accretion C by S+G
    prevp: evaporation R->V
    Cold (supcol>0):
      praci, piacr: ice-rain interactions
      psaci, psacw: snow-cloud, snow-water accretion
      pgacw, pgaci: graupel-cloud, graupel-ice accretion
      pracs, psacr: snow-rain, rain-snow interactions
      pgacr, pgacs: graupel-rain, graupel-snow accretion
    T>=T0 melting:
      pseml: enhanced snow melt
      pgeml: enhanced graupel melt
    Cold (supcol>0) deposition/nucleation:
      pidep, psdep, pgdep, pigen, psaut, pgaut
    T>=T0 evaporation:
      psevp, pgevp

  Group 14 — Mass conservation check and state update (lines 1200-1388)
    Two branches: T<=T0 (cold) and T>T0 (warm)
    Each branch: scale rates to satisfy mass bounds,
    then apply all rates to q/qc/qi/qr/qs/qg/t in one update

  Group 15 — Second qsatw computation (lines 1390-1420)
    Second fpvs inline expansion, overwrites qsatw only
    (qsati not recomputed here)

  Group 16 — pcond: condensation/evaporation (lines 1427-1437)
    workcond = wsm6_conden(t,q,qsatw,xl,cpm)  (local Real)
    Updates q, qc, t inline

  Group 17 — Clamp qc/qi to 0 if <= qmin (lines 1444-1449)

END loops/dtcld loop

  Group 18 — rainprod2d, evapprod2d (optional, lines 1452-1459)

---

## Rule 28: FAB Allocation Strategy for WSM6Ind Working Arrays

All WSM6Ind entries are allocated as FArrayBox locals inside the
MFIter loop in Advance(), scoped to each tile. They are NOT
persistent MultiFabs — they are per-MFIter temporaries.

3D FABs (fab_box, 1 component each) — one per WSM6Ind entry:
  State (already in mic_fab_vars, accessed via Array4):
    t, q, qc, qi, qr, qs, qg, den, p — use mic_fab_vars directly
    delz — already allocated as FArrayBox delz_fab(fab_box,1)

  Sedimentation work (3D, fab_box):
    fallc, falkc, work1c, work2c, workr, worka, den_tmp, delz_tmp

  Process rates (3D, fab_box, one each):
    pigen, pidep, pcond, prevp, psevp, pgevp, psdep, pgdep,
    praut, psaut, pgaut, piacr, pracw, praci, pracs,
    psacw, psaci, psacr, pgacw, pgaci, pgacr, pgacs, paacw,
    psmlt, pgmlt, pseml, pgeml

  Derived thermodynamic (3D, fab_box):
    qsum, xl, cpm, work2, denfac, xni,
    denqrs1, denqrs2, denqrs3, denqci, n0sfac

  Multi-component (3D, fab_box):
    rhw, rhi                           (1 FAB each)
    qsatw, qsati                       (1 FAB each)
    rslope_r, rslope_s, rslope_g       (1 FAB each)
    rslope2_r, rslope2_s, rslope2_g    (1 FAB each)
    rslope3_r, rslope3_s, rslope3_g    (1 FAB each)
    rslopeb_r, rslopeb_s, rslopeb_g    (1 FAB each)
    qrs_tmp_r, qrs_tmp_s, qrs_tmp_g    (1 FAB each)
    falk_r, falk_s, falk_g             (1 FAB each)
    fall_r, fall_s, fall_g             (1 FAB each)
    work1_r, work1_s, work1_g          (1 FAB each)
    workdiffw, workdiffi               (1 FAB each)

box2d FABs (box2d = fab_box slabbed at k=0, 1 component each):
  rainncv, sr, snowncv, graupelncv     (already in current code)
  delqrs1, delqrs2, delqrs3, delqi     (1D per-column, new)
  tstepsnow, tstepgraup                (1D per-column, new)

NOT allocated as FABs (Rule 23):
  mstep  — replace with literal 1
  numdt  — not used in surviving code paths
  flgcld — dead code, omit

Column-stack arrays (local inside nislfv ParallelFor lambda):
  All nislfv local arrays (dz, ww, qq, wi, zi, za, qn, etc.)
  sized WSM6_MAX_LEVELS — allocated on stack inside lambda,
  not as FABs

Allocation pattern in code:
  FArrayBox NAME_fab(fab_box, 1);
  auto const& NAME_arr = NAME_fab.array();
  // then zero-initialize where Fortran initializes to 0:
  NAME_fab.setVal(Real(0.0));

## Rule 29: The loops/dtcld Minor Timestep Outer Loop

Fortran mp_wsm6_run lines 484-486, 488, 1450:
  loops = max(nint(delt/dtcldcr), 1)
  dtcld = delt/loops
  if (delt <= dtcldcr) dtcld = delt
  do loop = 1, loops
    ... all physics groups G1-G17 ...
  enddo

In C++ Advance(), this becomes a sequential for loop wrapping
all physics ParallelFor passes:

  const int loops = std::max(static_cast<int>(
                      std::round(dt / dtcldcr)), 1);
  const Real dtcld = dt / static_cast<Real>(loops);

  for (int loop = 0; loop < loops; ++loop) {
      // ParallelFor: G1 denfac
      // ParallelFor: G1 qsat/rh
      // ParallelFor: G2 zero process rates
      // ParallelFor: G3 xni
      // ParallelFor: G4 pack qrs_tmp, slope_wsm6
      // ParallelFor: G5 sedimentation (column kernels)
      // ... all groups in order ...
      // ParallelFor: G17 clamp qc/qi
  }

The for loop is NOT parallelized — dtcldcr=120s subdivision
is inherently sequential (each sub-step depends on prior state).
ParallelFor is only over spatial indices within each sub-step.

dtcldcr = 120.0 is already a static constexpr in WSM6 class.
---

## Rule 30: Phase 3 Diagnostic Instrumentation Protocol

Purpose: establish per-process ground truth from the Fortran path
for comparison against C++ process-by-process in Phase 4.

Target: single column i=ilo, j=jlo at timestep 1, all k levels.
Print to stdout in format parseable by diff or python script.

Fortran instrumentation pattern (add to mp_wsm6_run after each
process group, guarded by a compile flag):

  #ifdef WSM6_DIAG
    if (i.eq.its .and. loop.eq.1) then
      write(*,'(A,I3,6E16.8)') 'PRAUT ', k, &
        praut(i,k), qr(i,k), qc(i,k), t(i,k), q(i,k), den(i,k)
    endif
  #endif

C++ instrumentation pattern (add to #else block after each
process ParallelFor, for i==ilo, j==jlo only):

  if (wsm6_diag && loop==0) {
      for (int k=klo; k<=khi; ++k) {
          amrex::Print() << "PRAUT " << k << " "
              << praut_arr(ilo,jlo,k) << " "
              << qr_arr(ilo,jlo,k) << "\n";
      }
  }

Diagnostic tags (one per process group, Rule 27 order):
  DENFAC   after G1b
  QSAT     after G1c
  SLOPE1   after G4
  NISLFV_R after G5b
  NISLFV_SG after G5c
  FALL     after G5d-e
  SLOPE2   after G6
  MELT     after G7 (psmlt, pgmlt)
  VICE     after G8
  PRECIP   after G9
  PHASE    after G10 (pimlt,pihmf,pihtf,pgfrz)
  SLOPE3   after G11
  PRAUT    after praut computed in G13
  PRACW    after pracw
  PREVP    after prevp
  PRACI    after praci/piacr
  PSACI    after psaci/psacw/pgacw/pgaci
  PRACS    after pracs/psacr/pgacr/pgacs
  PSEML    after pseml/pgeml
  PIDEP    after pidep/psdep/pgdep/pigen
  PSAUT    after psaut/pgaut
  PSEVP    after psevp/pgevp
  UPDATE   after G14 state update (q,qc,qi,qr,qs,qg,t)
  QSAT2    after G15
  PCOND    after G16

Each diagnostic point is independent — add one, validate one,
then add the next. Never add all diagnostics at once.

Compile flag: -DWSM6_DIAG added to CMake for diagnostic builds only.
Remove all diagnostic prints before Phase 5 assembly.

## Rule 31: Runtime Fortran/C++ toggle — Morrison pattern

Use a runtime ParmParse query, not a build-time flag, to switch
between Fortran and C++ paths. Build-time flags produce the same
target name and require make clean to switch reliably.

Pattern (matches Morrison ERF_AdvanceMorrison.cpp lines 520-521):

  bool use_wsm6_cpp_answer = false;
  { amrex::ParmParse pp("erf");
    pp.query("use_wsm6_cpp_answer", use_wsm6_cpp_answer); }
  bool run_wsm6_fort = !use_wsm6_cpp_answer;

In inputs file:
  erf.use_wsm6_cpp_answer = 0   # Fortran path
  erf.use_wsm6_cpp_answer = 1   # native C++ path

Reference: ERF_AdvanceWSM6.cpp lines 672, 781

## Rule 32: C++ validation pipeline — debug-first, backtrace-before-prints

Order of operations for validating a new C++ physics path:

1. Build with DEBUG=TRUE. This gets -ftrapv, AMReX assertions,
   and Fortran -fbacktrace at no source cost. Do not add
   diagnostic prints before doing this.

2. Run the Fortran path first, save the plotfile before the
   next run overwrites it:
     cp -r plt_<case>00002 plt_<case>_fort_00002

3. Run the C++ path with identical inputs and step count.

4. On crash, read Backtrace.0 immediately. It gives file + line
   number without any added instrumentation. The format is:
     function_name at ERF_AdvanceWSM6.cpp:LINE
   Cross-reference that line against the Fortran source
   simultaneously (two subagents, one per file).

5. fcompare only after both paths run clean — it is a
   correctness gate, not a debug tool. A crash must be fixed
   before fcompare has meaning.

6. Diagnostic prints (microphysics_debug levels, PRE/POST
   reductions) are for understanding WHY something is wrong
   after Backtrace.0 tells you WHERE. See:
   Docs/wsm6_deep_dive/wsm6_diagnostic_sidecar_notes.md

Key finding from WSM6 port: Backtrace.0 is faster than any
diagnostic print strategy for locating a crash. Prints are
the scalpel for bisecting correctness divergence, not crashes.

Common crash sources in nislfv-style advection routines:
- dza(k) = za(k+1) - za(k) can be zero when arrival points
  collapse; always guard divisions by dza with a minimum.
- km passed as khi-klo instead of khi-klo+1 leaves garbage
  in the last array slot, causing downstream zero denominators.
- Column temporaries sized WSM6_MAX_LEVELS must be filled
  exactly km = khi-klo+1 elements, no more, no less.

## Rule 33: Pre-compile reflexive alignment pass — group-by-group

After writing a group in C++ but BEFORE compiling, show the
Fortran source for that group alongside the C++ (two subagents,
one per file, same line range) and scan for:

1. Statements OUTSIDE the main loop body: initializations
   before loops, sentinel assignments after loops, scalar
   resets between loops. These are the most commonly dropped
   lines in translation because they look like boilerplate.
   Motivating example: Fortran line 1866
     dza(km+1) = zi(km+1) - za(km+1)
   was dropped in the C++ translation of update_wind_and_state,
   leaving dza[km] uninitialized and producing NaN downstream.

2. Loop bounds: Fortran do k=1,km → C++ for(k=0;k<km;++k)
                Fortran do k=1,km+1 → C++ for(k=0;k<=km;++k)
   Arrays sized km+1 in Fortran need index [km] in C++.

3. Every array with size km+1 — verify the km+1 index is
   always written before it is read.

Translation bugs cluster at loop boundaries, not inside loop
bodies. The body gets transcribed carefully because it is
visually prominent. Before/after statements get dropped
because they look like boilerplate.

## Rule 34: Fortran control flow restructuring hazards

Before translating any Fortran group, scan for these patterns:

1. goto + labeled continue (backward jump = loop):
   Fortran:
     n=1
     100 continue
       ! body
       if (n.le.iter) then
         n=n+1
         go to 100
       endif
   Must become an explicit while/for loop in C++, NOT a
   one-shot if block. A one-shot if handles iter=0 and iter=1
   correctly but silently fails for iter>1.
   Always check ALL call sites for non-zero iter arguments
   even when the immediate test case passes.
   Motivating example: nislfv_rain_plm goto 100 loop,
   currently called with iter=0/1 only — silent for now.

2. cycle → continue, exit → break, but verify which enclosing
   construct they target. Fortran named loops (i_loop:,
   find_kb:) make targets explicit; C++ does not.

3. go to that jumps forward past an else block is an early
   exit — translate as a labeled break or restructure to
   if/else. Do not use C++ goto.
