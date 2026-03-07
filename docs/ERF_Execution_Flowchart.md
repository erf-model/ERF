# ERF Execution Flowchart — Isentropic Vortex (No Acoustic Substepping, Hybrid 5/6)

This document traces the complete execution path of ERF for the isentropic vortex test case. All formulas are extracted directly from the implementation source code.

**Simplifications active for this case:**
- Dry (no moisture: `MoistureType::None`)
- No terrain (`MeshType::ConstantDz`, `detJ = 1`, `ax = ay = az = 1`)
- No diffusion (`molec_diff_type = "None"`, `les_type = "None"`)
- No gravity (`use_gravity = false`, so buoyancy = 0)
- No Coriolis, no sponge, no Rayleigh damping
- Compressible (not anelastic)
- No acoustic substepping (`no_substepping = true` when `fixed_mri_dt_ratio` is small or forced)
- Blended 5th/6th order advection (`Blended_5th6th` with `upw_frac`)
- Single level (`max_level = 0`)
- Periodic in x,y; SlipWall in z

---

## 0. Prognostic Variables and Grid Layout

### Arakawa C-Grid Staggering

```
              y
              ^
              |
         +---------+-- v(i,j+1,k) --+---------+
         |         |                 |         |
         |         |                 |         |
    u(i,j,k)    [rho, Theta]    u(i+1,j,k)
         |       (i,j,k)            |         |
         |     CELL CENTER           |         |
         +---------+-- v(i,j,k) ----+---------+  --> x
                   |
                   w(i,j,k) on bottom z-face
                   w(i,j,k+1) on top z-face
```

| Variable | Symbol | Location | Index |
|----------|--------|----------|-------|
| Density | `rho` | Cell center | `Rho_comp = 0` |
| Conserved pot. temp. | `Theta = rho * theta` | Cell center | `RhoTheta_comp = 1` |
| x-momentum | `U = rho * u` | x-face (i+1/2, j, k) | `IntVars::xmom` |
| y-momentum | `V = rho * v` | y-face (i, j+1/2, k) | `IntVars::ymom` |
| z-momentum | `W = rho * w` | z-face (i, j, k+1/2) | `IntVars::zmom` |

### Physical Constants
*Source: `ERF_Constants.H`*

| Constant | Value | Symbol |
|----------|-------|--------|
| `R_d` | 287.0 J/(kg·K) | Dry air gas constant |
| `Cp_d` | 1004.5 J/(kg·K) | Specific heat at constant pressure |
| `p_0` | 1.0e5 Pa | Reference pressure |
| `Gamma` | 1.4 | c_p / c_v |
| `iGamma` | 1/1.4 | Inverse of Gamma |
| `rdOcp` | R_d / Cp_d ≈ 0.2857 | Kappa |

### Equation of State
*Source: `ERF_EOS.H:getPgivenRTh()`*

```
p = p_0 * ( R_d * (rho*theta) / p_0 )^Gamma
```

In code:
```cpp
// ERF_EOS.H line 83
return p_0 * std::pow(R_d * rhotheta * ip_0, Gamma);
```

### Exner Function
*Source: `ERF_EOS.H:getExnergivenRTh()`*

```
pi = ( R_d * (rho*theta) / p_0 )^(Gamma * rdOcp)
```

In code:
```cpp
// ERF_EOS.H line 161
return std::pow(R_d * rhotheta * ip_0, Gamma * rdOcp);
```

---

## 1. Top-Level Time Step: `ERF::Advance()`

*Source: `ERF_Advance.cpp`*

```
┌─────────────────────────────────────────────────────────────────┐
│                    ERF::Advance(lev, time, dt)                  │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. std::swap(vars_old, vars_new)                               │
│     // Previous step's "new" becomes this step's "old"          │
│                                                                 │
│  2. VelocityToMomentum(U_old, V_old, W_old, S_old,             │
│                         rU_old, rV_old, rW_old)                 │
│     // Convert fillpatched velocities to momenta                │
│     // INCLUDING ghost cells                                    │
│                                                                 │
│  3. Allocate source arrays:                                     │
│     - cc_source  (cell-centered, nvars components)              │
│     - xmom_source (x-face)                                     │
│     - ymom_source (y-face)                                     │
│     - zmom_source (z-face)                                     │
│     - buoyancy    (z-face)                                     │
│     All initialized to zero.                                    │
│                                                                 │
│  4. Build state_old = {S_old, rU_old, rV_old, rW_old}          │
│     Build state_new = {S_new, rU_new, rV_new, rW_new}          │
│                                                                 │
│  5. advance_dycore(lev, state_old, state_new, ...)              │
│     // >>> MAIN DYCORE CALL — see Section 2 <<<                 │
│                                                                 │
│  6. Post-dycore checks (NaN, negative theta)                    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Velocity-to-Momentum Conversion
*Source: `ERF_VelocityToMomentum.cpp`*

Density is averaged to faces by simple arithmetic mean:

```
U(i+1/2, j, k)  =  u(i+1/2, j, k)  *  0.5 * ( rho(i,j,k) + rho(i-1,j,k) )
V(i, j+1/2, k)  =  v(i, j+1/2, k)  *  0.5 * ( rho(i,j,k) + rho(i,j-1,k) )
W(i, j, k+1/2)  =  w(i, j, k+1/2)  *  0.5 * ( rho(i,j,k) + rho(i,j,k-1) )
```

---

## 2. Dycore Advance: `ERF::advance_dycore()`

*Source: `ERF_AdvanceDycore.cpp`*

```
┌─────────────────────────────────────────────────────────────────┐
│                   advance_dycore(level, ...)                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. Extract base state: r_hse, p_hse, pi_hse                   │
│                                                                 │
│  2. Allocate work arrays:                                       │
│     - S_prim   (primitive vars = cons / rho)                    │
│     - pi_stage (Exner function at cell centers)                 │
│     - fast_coeffs (for acoustic substepping — unused here)      │
│                                                                 │
│  3. [SKIPPED for isentropic vortex]:                            │
│     - Strain computation (l_use_diff = false)                   │
│     - Eddy viscosity (l_use_kturb = false)                      │
│                                                                 │
│  4. VelocityToMomentum on old state                             │
│     Copy old velocities into new velocities                     │
│                                                                 │
│  5. apply_bcs(state_old) — fill ghost cells                     │
│     cons_to_prim(state_old) — compute primitives                │
│                                                                 │
│  6. #include "ERF_TI_no_substep_fun.H"  — define no_substep    │
│     #include "ERF_TI_slow_rhs_pre.H"    — define slow_rhs_pre  │
│     #include "ERF_TI_slow_rhs_post.H"   — define slow_rhs_post │
│                                                                 │
│  7. Configure MRI integrator:                                   │
│     mri_integrator.set_slow_rhs_pre(slow_rhs_fun_pre)          │
│     mri_integrator.set_slow_rhs_post(slow_rhs_fun_post)        │
│     mri_integrator.set_no_substep(no_substep_fun)               │
│                                                                 │
│  8. mri_integrator.advance(state_old, state_new, time, dt)      │
│     // >>> RK3 LOOP — see Section 3 <<<                         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Primitive Variable Computation (`cons_to_prim`)
*Source: `ERF_TI_utils.H`*

```
theta(i,j,k) = (rho*theta)(i,j,k) / rho(i,j,k)

pi(i,j,k) = getExnergivenRTh( (rho*theta)(i,j,k) )
           = ( R_d * (rho*theta)(i,j,k) / p_0 )^(Gamma * rdOcp)
```

---

## 3. RK3 Integration Loop: `MRISplitIntegrator::advance()`

*Source: `ERF_MRI.H`*

The RK3 scheme advances the state vector **S** = {rho, rho*theta, U, V, W} from time `t^n` to `t^{n+1}`:

```
Stage 0:  S* = S^n   + (dt/3) * f(S^n)        [advances t^n  -> t^n + dt/3]
Stage 1:  S** = S^n  + (dt/2) * f(S*)         [advances t^n  -> t^n + dt/2]
Stage 2:  S^{n+1} = S^n + dt  * f(S**)        [advances t^n  -> t^n + dt  ]
```

where `f(S)` represents the full right-hand side evaluation.

### RK3 Loop Structure (No Substepping)

```
┌─────────────────────────────────────────────────────────────────────┐
│              MRISplitIntegrator::advance()                          │
│                                                                     │
│  Copy S_old -> S_new (initialization)                               │
│                                                                     │
│  FOR nrk = 0, 1, 2:                                                │
│  ┌─────────────────────────────────────────────────────────────┐    │
│  │                                                             │    │
│  │  Compute slow_dt for this stage:                            │    │
│  │    nrk=0: slow_dt = dt/3     (nsubsteps=1, dtau=dt/3)      │    │
│  │    nrk=1: slow_dt = dt/2     (nsubsteps=1, dtau=dt/2)      │    │
│  │    nrk=2: slow_dt = dt       (nsubsteps=1, dtau=dt)         │    │
│  │                                                             │    │
│  │  ┌───────────────────────────────────────────────┐          │    │
│  │  │ slow_rhs_pre(F_slow, S_old, S_new, nrk)      │          │    │
│  │  │  -> Evaluate RHS for {rho, rho*theta, U,V,W}  │          │    │
│  │  │  -> Stores result in F_slow                    │          │    │
│  │  │  -> See Section 4                              │          │    │
│  │  └───────────────────────────────────────────────┘          │    │
│  │                                                             │    │
│  │  ┌───────────────────────────────────────────────┐          │    │
│  │  │ no_substep(S_sum, S_old, F_slow, slow_dt)     │          │    │
│  │  │  -> Update fast variables:                     │          │    │
│  │  │     S_sum = S_old + slow_dt * F_slow           │          │    │
│  │  │  -> For rho, rho*theta, U, V, W               │          │    │
│  │  │  -> apply_bcs(S_sum) for fast vars             │          │    │
│  │  │  -> See Section 5                              │          │    │
│  │  └───────────────────────────────────────────────┘          │    │
│  │                                                             │    │
│  │  ┌───────────────────────────────────────────────┐          │    │
│  │  │ slow_rhs_post(F_slow, S_old, S_new, S_sum)    │          │    │
│  │  │  -> Evaluate RHS for slow scalars              │          │    │
│  │  │  -> Update: S_new = S_old + slow_dt * RHS      │          │    │
│  │  │  -> apply_bcs(S_new) for all vars              │          │    │
│  │  │  -> See Section 6                              │          │    │
│  │  └───────────────────────────────────────────────┘          │    │
│  │                                                             │    │
│  └─────────────────────────────────────────────────────────────┘    │
│  END FOR                                                            │
│                                                                     │
│  After nrk=2: S_new holds S^{n+1}                                  │
└─────────────────────────────────────────────────────────────────────┘
```

**Key: what lives where after each call:**

| Array | After `slow_rhs_pre` | After `no_substep` | After `slow_rhs_post` |
|-------|---------------------|--------------------|-----------------------|
| `F_slow` | RHS for rho, rho*theta, U, V, W | unchanged | RHS for slow scalars |
| `S_new` | current stage data | current stage data | updated stage data (all vars) |
| `S_sum` | — | updated fast vars | — |

---

## 4. Slow RHS Pre: `slow_rhs_fun_pre` → `erf_slow_rhs_pre()`

*Source: `ERF_TI_slow_rhs_pre.H`, `ERF_SlowRhsPre.cpp`*

This is the main RHS evaluation. It computes tendencies for the **fast variables** (rho, rho*theta) and **all three momenta** (U, V, W).

```
┌─────────────────────────────────────────────────────────────────────┐
│                       slow_rhs_fun_pre                              │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  IF nrk > 0:                                                        │
│     cons_to_prim(S_data)   // recompute theta, pi from S_data       │
│                                                                     │
│  A. make_sources(cc_src)                                            │
│     -> For isentropic vortex: cc_src = 0 (no external forcing)      │
│                                                                     │
│  B. make_gradp_pert(gradp)       >>> Section 4.1 <<<                │
│     -> Compute perturbation pressure gradient                       │
│                                                                     │
│  C. make_buoyancy(buoyancy)      >>> Section 4.2 <<<                │
│     -> For isentropic vortex: buoyancy = 0 (no gravity)             │
│                                                                     │
│  D. make_mom_sources(xmom_src, ymom_src, zmom_src)                  │
│     -> For isentropic vortex: all = 0 (no Coriolis, etc.)           │
│                                                                     │
│  E. erf_slow_rhs_pre(...)        >>> Section 4.3 <<<                │
│     -> Advection of {rho, rho*theta} on cell centers                │
│     -> Advection of {U, V, W} on faces                              │
│     -> Add pressure gradient + buoyancy + source terms              │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

### 4.1 Pressure Gradient: `make_gradp_pert()`

*Source: `ERF_MakeGradP.cpp`*

First, compute perturbation pressure at cell centers:

```
p'(i,j,k) = getPgivenRTh( (rho*theta)(i,j,k) ) - p0(i,j,k)
           = p_0 * ( R_d * (rho*theta)(i,j,k) / p_0 )^Gamma  -  p0(i,j,k)
```

Then take the gradient on faces (no terrain, so simple differences):

```
gpx(i,j,k) = (1/dx) * ( p'(i,j,k) - p'(i-1,j,k) )     [at x-face]
gpy(i,j,k) = (1/dy) * ( p'(i,j,k) - p'(i,j-1,k) )     [at y-face]
gpz(i,j,k) = (1/dz) * ( p'(i,j,k) - p'(i,j,k-1) )     [at z-face]
```

*Source: `ERF_MakeGradP.cpp` lines 130–203 (no-terrain branch)*

### 4.2 Buoyancy

*Source: `ERF_MakeBuoyancy.cpp`*

For the isentropic vortex with `use_gravity = false`:

```
buoyancy(i,j,k) = 0    for all (i,j,k)
```

If gravity were enabled with `buoyancy_type = 1` (dry compressible, `buoyancy_rhopert`):
```
rho'_hi = rho(i,j,k)   - rho0(i,j,k)
rho'_lo = rho(i,j,k-1) - rho0(i,j,k-1)
buoyancy(i,j,k) = g * 0.5 * (rho'_hi + rho'_lo)
```
*Source: `ERF_BuoyancyUtils.H` lines 117–126*

### 4.3 Core RHS Computation: `erf_slow_rhs_pre()`

*Source: `ERF_SlowRhsPre.cpp`*

#### 4.3.1 Contravariant Flux (Omega)

For no terrain, the vertical "Omega" flux is simply the z-momentum:

```
Omega(i,j,k) = W(i,j,k) = (rho*w)(i,j,k)
```

*Source: `ERF_SlowRhsPre.cpp` lines 414–418*

#### 4.3.2 Advection of Density (Continuity Equation)

*Source: `ERF_AdvectionSrcForState.cpp` via `AdvectionSrcForRho()`*

The continuity equation RHS is the negative divergence of mass flux. Since the mass fluxes are just the momenta on faces:

```
RHS_rho(i,j,k) = - [ (U(i+1,j,k) - U(i,j,k)) / dx
                    + (V(i,j+1,k) - V(i,j,k)) / dy
                    + (W(i,j,k+1) - W(i,j,k)) / dz ]
```

The `avg_xmom`, `avg_ymom`, `avg_zmom` arrays are set equal to the current momenta (these are later used for scalar advection).

#### 4.3.3 Advection of rho*theta (Potential Temperature Equation)

*Source: `ERF_AdvectionSrcForState.cpp` via `AdvectionSrcForScalars()`*

Using the **Blended 5th/6th order** (Wicker-Skamarock) scheme:

The advective flux at x-face `(i+1/2, j, k)` for the scalar `psi = theta` (primitive variable) advected by momentum `U`:

**6th-order centered component:**
```
F6(i+1/2) = U(i+1/2) / 60 * [ 37*(psi(i+1) + psi(i))
                              -  8*(psi(i+2) + psi(i-1))
                              +    (psi(i+3) + psi(i-2)) ]
```

**5th-order upwind dissipation:**
```
D5(i+1/2) = |U(i+1/2)| / 60 * [ 10*(psi(i+1) - psi(i))
                                -  5*(psi(i+2) - psi(i-1))
                                +    (psi(i+3) - psi(i-2)) ]
```

**Blended flux:**
```
F(i+1/2) = F6(i+1/2) - upw_frac * D5(i+1/2)
```

where `upw_frac` is the user-specified upwinding fraction:
- `upw_frac = 0`: pure 6th-order centered (non-dissipative)
- `upw_frac = 1`: full 5th-order upwind
- `upw_frac = 0.5`: half dissipation

The RHS for `rho*theta` is the negative divergence:
```
RHS_{rho*theta}(i,j,k) = - [ (Fx(i+1,j,k) - Fx(i,j,k)) / dx
                            + (Fy(i,j+1,k) - Fy(i,j,k)) / dy
                            + (Fz(i,j,k+1) - Fz(i,j,k)) / dz ]
```

*Source: `ERF_Discretization_Explained.md` Section 4; implementation in `ERF_AdvectionSrcForState.cpp`*

#### 4.3.4 Advection of Momentum

*Source: `ERF_AdvectionSrcForMom.cpp` via `AdvectionSrcForMom()`*

The momentum advection follows the same Blended 5th/6th order scheme but on staggered grids.

**x-momentum equation** (U lives at x-faces):

The x-flux of U is evaluated at **cell centers** (i,j,k):
- Advecting velocity: `u_c(i,j,k) = 0.5 * (u(i-1/2) + u(i+1/2))`
- Advected quantity: `U = rho*u` interpolated with WS5/6 stencil

The y-flux of U is evaluated at **(i+1/2, j+1/2, k)** cell edges:
- Advecting velocity: `v` averaged from y-faces to the edge
- Advected quantity: `U` interpolated in y

The z-flux of U is evaluated at **(i+1/2, j, k+1/2)** cell edges:
- Advecting velocity: `Omega` (= `W` for no terrain) averaged to the edge
- Advected quantity: `U` interpolated in z

Analogous for V (y-momentum) and W (z-momentum).

#### 4.3.5 Combining the Momentum RHS

*Source: `ERF_SlowRhsPre.cpp` lines 702–816*

After the advection source is computed, the pressure gradient, buoyancy, and external sources are added:

**x-momentum:**
```
RHS_U(i,j,k) = [advection_src_U]  +  (-gpx(i,j,k))  +  xmom_src(i,j,k)
```

**y-momentum:**
```
RHS_V(i,j,k) = [advection_src_V]  +  (-gpy(i,j,k))  +  ymom_src(i,j,k)
```

**z-momentum:**
```
RHS_W(i,j,k) = [advection_src_W]  +  (-gpz(i,j,k) + buoyancy(i,j,k))  +  zmom_src(i,j,k)
```

For the isentropic vortex: `buoyancy = 0`, `xmom_src = ymom_src = zmom_src = 0`.

#### 4.3.6 Adding Cell-Centered Sources to rho and rho*theta

*Source: `ERF_SlowRhsPre.cpp` lines 632–636*

```
RHS_rho(i,j,k)         += cc_src(i,j,k, Rho_comp)
RHS_{rho*theta}(i,j,k) += cc_src(i,j,k, RhoTheta_comp)
```

For the isentropic vortex: `cc_src = 0`.

---

## 5. No-Substep Update: `no_substep_fun`

*Source: `ERF_TI_no_substep_fun.H`*

When there is no acoustic substepping, the fast variables are updated with a simple forward Euler step using the slow RHS:

```
S_sum[rho](i,j,k)         = S_old[rho](i,j,k)         + slow_dt * F_slow[rho](i,j,k)
S_sum[rho*theta](i,j,k)   = S_old[rho*theta](i,j,k)   + slow_dt * F_slow[rho*theta](i,j,k)

S_sum[xmom](i,j,k) = S_old[xmom](i,j,k) + slow_dt * F_slow[xmom](i,j,k)
S_sum[ymom](i,j,k) = S_old[ymom](i,j,k) + slow_dt * F_slow[ymom](i,j,k)
S_sum[zmom](i,j,k) = S_old[zmom](i,j,k) + slow_dt * F_slow[zmom](i,j,k)
```

*Source: `ERF_TI_no_substep_fun.H` lines 113–176 (no-terrain, no-EB branch)*

After updating, boundary conditions are applied for the fast variables only:

```
apply_bcs(S_sum, fast_only=true)
```

Note: `S_sum` is the key output. With no substepping, `slow_dt` corresponds to the full RK stage time step.

---

## 6. Slow RHS Post: `slow_rhs_fun_post` → `erf_slow_rhs_post()`

*Source: `ERF_TI_slow_rhs_post.H`, `ERF_SlowRhsPost.cpp`*

This function handles the **slow variables** — scalars other than rho and rho*theta (e.g., RhoKE, RhoScalar, moisture variables). For the dry isentropic vortex with no TKE:

```
┌─────────────────────────────────────────────────────────────────┐
│                     slow_rhs_fun_post                           │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. Copy slow variables from S_new into S_data (cur_cons)       │
│     // So cur_cons has "new" slow vars + "current" fast vars    │
│                                                                 │
│  2. For each valid slow variable (RhoScalar, etc.):             │
│     a. Compute advection source using avg_xmom/ymom/zmom       │
│        (the time-averaged momenta from substepping, or just     │
│         the current momenta when no substepping)                │
│     b. [SKIPPED: Diffusion (l_use_diff = false)]                │
│                                                                 │
│  3. Update slow variables:                                      │
│     cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * RHS(i,j,k,n) │
│                                                                 │
│  4. Copy all of cur_cons -> new_cons                            │
│     Copy all of cur_mom  -> new_mom                             │
│                                                                 │
│  5. apply_bcs(S_new, fast_only=false)                           │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

For the isentropic vortex with `transport_scalar = true` and `RhoScalar_comp`:

```
RHS_scalar = -div(u * scalar) + cc_src[scalar]

S_new[RhoScalar](i,j,k) = S_old[RhoScalar](i,j,k) + slow_dt * RHS_scalar(i,j,k)
```

*Source: `ERF_SlowRhsPost.cpp` lines 551–563*

---

## 7. Complete RK3 Time Step — Putting It All Together

Let `f()` denote the combined RHS operator (slow_rhs_pre + no_substep + slow_rhs_post). A single full time step proceeds as:

```
GIVEN: S^n at time t^n, time step dt

═══════════════════════════════════════════════════════════════
STAGE 0  (nrk = 0):   slow_dt = dt/3
═══════════════════════════════════════════════════════════════

  S_data = S^n  (S_new was initialized from S_old)

  (a) slow_rhs_pre:  F_slow = f_fast(S^n)
      - cons_to_prim (already done before RK loop for nrk=0)
      - Compute gradp from S^n
      - Compute advection of rho, rho*theta from S^n
      - Compute advection of U, V, W from S^n
      - Add pressure gradient + sources

  (b) no_substep:    S_sum = S^n + (dt/3) * F_slow
      - Updates rho, rho*theta, U, V, W
      - apply_bcs on fast vars

  (c) slow_rhs_post: RHS_slow = f_slow(S^n, S_sum)
      - Advect slow scalars using momenta from S_sum
      - S_new[slow] = S^n + (dt/3) * RHS_slow
      - Copy fast vars from S_sum into S_new
      - apply_bcs on all vars

  RESULT: S_new = S* ≈ S^n + (dt/3) * f(S^n)

═══════════════════════════════════════════════════════════════
STAGE 1  (nrk = 1):   slow_dt = dt/2
═══════════════════════════════════════════════════════════════

  S_data = S* (S_new from previous stage)

  (a) slow_rhs_pre:  F_slow = f_fast(S*)
      - cons_to_prim(S*)  [nrk > 0, so recompute]
      - Compute gradp from S*
      - Compute advection of rho, rho*theta from S*
      - Compute advection of U, V, W from S*
      - Add pressure gradient + sources

  (b) no_substep:    S_sum = S^n + (dt/2) * F_slow
      - NOTE: always adds to S_old = S^n, not S*!

  (c) slow_rhs_post: RHS_slow = f_slow(S*, S_sum)
      - S_new[slow] = S^n + (dt/2) * RHS_slow

  RESULT: S_new = S** ≈ S^n + (dt/2) * f(S*)

═══════════════════════════════════════════════════════════════
STAGE 2  (nrk = 2):   slow_dt = dt
═══════════════════════════════════════════════════════════════

  S_data = S** (S_new from previous stage)

  (a) slow_rhs_pre:  F_slow = f_fast(S**)
      - cons_to_prim(S**)
      - Compute gradp from S**
      - Compute advection of rho, rho*theta from S**
      - Compute advection of U, V, W from S**
      - Add pressure gradient + sources

  (b) no_substep:    S_sum = S^n + dt * F_slow
      - NOTE: always adds to S_old = S^n, not S**!

  (c) slow_rhs_post: RHS_slow = f_slow(S**, S_sum)
      - S_new[slow] = S^n + dt * RHS_slow

  RESULT: S_new = S^{n+1} ≈ S^n + dt * f(S**)

═══════════════════════════════════════════════════════════════
```

**Important:** In each stage, the update formula is always `S^n + slow_dt * F_slow(S_stage)`, where `S^n` is the state at the **beginning** of the full time step (not the beginning of the RK stage). This is the defining property of a low-storage RK3 scheme.

---

## 8. Data Flow Diagram

```
                           ┌─────────┐
                           │  S^n    │ (state_old)
                           └────┬────┘
                                │
            ┌───────────────────┼───────────────────┐
            │                   │                   │
    ┌───────▼────────┐  ┌──────▼───────┐  ┌────────▼──────┐
    │ nrk=0          │  │ nrk=1        │  │ nrk=2         │
    │ eval f(S^n)    │  │ eval f(S*)   │  │ eval f(S**)   │
    │                │  │              │  │               │
    │ S^n+dt/3*f     │  │ S^n+dt/2*f   │  │ S^n+dt*f      │
    │  = S*          │  │  = S**       │  │  = S^{n+1}    │
    └───────┬────────┘  └──────┬───────┘  └────────┬──────┘
            │                   │                   │
            └─► S_new ─────────┘                   │
                (carries stage data between RK)    │
                                                    │
                                           ┌────────▼──────┐
                                           │   S^{n+1}     │
                                           └───────────────┘
```

---

## 9. Isentropic Vortex Simplification Summary

For the isentropic vortex, the full time step reduces to:

**Per RK stage, the effective update for each variable is:**

```
rho^{stage} = rho^n + slow_dt * [ -div(rho * u) ]

(rho*theta)^{stage} = (rho*theta)^n + slow_dt * [ -div(u * (rho*theta))_WS56 ]

U^{stage} = U^n + slow_dt * [ -div(u * U)_WS56  -  dp'/dx ]

V^{stage} = V^n + slow_dt * [ -div(u * V)_WS56  -  dp'/dy ]

W^{stage} = W^n + slow_dt * [ -div(u * W)_WS56  -  dp'/dz ]
```

where:
- `div()_WS56` denotes divergence computed with Wicker-Skamarock Blended 5th/6th order fluxes
- `p' = p_0 * (R_d * (rho*theta) / p_0)^Gamma - p0_hse` is the perturbation pressure
- `slow_dt` is `dt/3`, `dt/2`, or `dt` for stages 0, 1, 2 respectively
- All updates are relative to `S^n` (not the previous stage)

---

## 10. Source File Reference

| Step | Source File | Key Function |
|------|------------|--------------|
| Top-level advance | `Source/TimeIntegration/ERF_Advance.cpp` | `ERF::Advance()` |
| Dycore setup | `Source/TimeIntegration/ERF_AdvanceDycore.cpp` | `ERF::advance_dycore()` |
| RK3 integrator | `Source/TimeIntegration/ERF_MRI.H` | `MRISplitIntegrator::advance()` |
| Slow RHS pre (wrapper) | `Source/TimeIntegration/ERF_TI_slow_rhs_pre.H` | `slow_rhs_fun_pre` lambda |
| Slow RHS pre (core) | `Source/TimeIntegration/ERF_SlowRhsPre.cpp` | `erf_slow_rhs_pre()` |
| Slow RHS post (wrapper) | `Source/TimeIntegration/ERF_TI_slow_rhs_post.H` | `slow_rhs_fun_post` lambda |
| Slow RHS post (core) | `Source/TimeIntegration/ERF_SlowRhsPost.cpp` | `erf_slow_rhs_post()` |
| No-substep update | `Source/TimeIntegration/ERF_TI_no_substep_fun.H` | `no_substep_fun` lambda |
| Utility functions | `Source/TimeIntegration/ERF_TI_utils.H` | `cons_to_prim`, `apply_bcs` |
| Pressure gradient | `Source/SourceTerms/ERF_MakeGradP.cpp` | `make_gradp_pert()` |
| Buoyancy | `Source/SourceTerms/ERF_MakeBuoyancy.cpp` | `make_buoyancy()` |
| Buoyancy formulas | `Source/SourceTerms/ERF_BuoyancyUtils.H` | `buoyancy_rhopert()` etc. |
| EOS functions | `Source/Utils/ERF_EOS.H` | `getPgivenRTh()`, `getExnergivenRTh()` |
| Vel-to-Mom | `Source/Utils/ERF_VelocityToMomentum.cpp` | `VelocityToMomentum()` |
| Advection (scalars) | `Source/Advection/ERF_AdvectionSrcForState.cpp` | `AdvectionSrcForScalars()` |
| Advection (momentum) | `Source/Advection/ERF_AdvectionSrcForMom.cpp` | `AdvectionSrcForMom()` |
| Physical constants | `Source/ERF_Constants.H` | `R_d`, `Gamma`, `p_0`, etc. |
| Discretization notes | `docs/ERF_Discretization_Explained.md` | WS5/WS6 formulas |
| Governing equations | `docs/ERF_Governing_Equations.md` | Continuous equations |

---

*Cross-referenced from: `docs/ERF_Governing_Equations.md`, `docs/ERF_Discretization_Explained.md`, and all source files listed above.*
