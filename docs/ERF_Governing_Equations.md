# ERF (Energy Research and Forecasting) — Governing & Difference Equations

This document compiles the governing equations and their discretized (difference) forms used in the ERF atmospheric modeling code, as documented across the ERF readthedocs documentation and the published JAMES paper (Lattanzi et al., 2025).

---

## 1. Fully Compressible Governing Equations (Continuous Form)

ERF solves a system of four coupled PDEs for compressible atmospheric dynamics. These are written in conservation form as:

### 1.1 Continuity Equation (Conservation of Mass)

$$\frac{\partial \rho}{\partial t} = -\nabla \cdot (\rho \mathbf{u}) = R_\rho$$

**Variables:**
- **ρ** — Air density (kg/m³)
- **u** = (u, v, w) — Velocity vector with components in the x, y, and z directions (m/s)
- **R_ρ** — Shorthand for all right-hand-side terms contributing to the time evolution of density

**Physical meaning:** The rate of change of air density at a fixed point equals the negative divergence of mass flux. Mass is neither created nor destroyed; it is only transported by the flow.

---

### 1.2 Momentum Equation (Conservation of Momentum)

$$\frac{\partial (\rho \mathbf{u})}{\partial t} = -\nabla \cdot (\rho \mathbf{u} \mathbf{u} + p' \mathbf{I}) + \mathbf{F}_\mathbf{u} = \mathbf{R}_\mathbf{u}$$

**Variables:**
- **ρu** = (U, V, W) = (ρu, ρv, ρw) — Momentum vector (kg/m²s)
- **ρuu** — Momentum flux tensor (the outer product of momentum and velocity)
- **p'** — Perturbation pressure (Pa), i.e., deviation from a hydrostatically balanced base state
- **I** — Identity tensor
- **F_u** = (F_U, F_V, F_W) — External forcing terms per unit volume (including gravity, Coriolis, turbulent diffusion, etc.)
- **R_u** — Shorthand for all RHS terms

**Physical meaning:** Momentum changes due to advective transport of momentum, pressure gradient forces, and external body forces (gravity, Coriolis, diffusion).

---

### 1.3 Potential Temperature Equation (Thermodynamic Equation)

$$\frac{\partial \Theta}{\partial t} = -\nabla \cdot (\mathbf{u} \Theta) + \rho \alpha_T \nabla^2 \theta = R_\Theta$$

**Variables:**
- **Θ** = ρθ — Conserved potential temperature variable (kg·K/m³), the product of density and potential temperature
- **θ** — Potential temperature (K), defined as T(p₀/p)^κ where T is temperature, p₀ is a reference pressure (typically 1000 hPa), and κ = R_d/c_p
- **α_T** — Thermal diffusivity (m²/s)
- **R_Θ** — Shorthand for all RHS terms

**Physical meaning:** The conserved potential temperature evolves through advective transport and thermal diffusion. Potential temperature is conserved for adiabatic processes and serves as ERF's prognostic thermodynamic variable.

---

### 1.4 Scalar Transport Equation

$$\frac{\partial (\rho C)}{\partial t} = -\nabla \cdot (\rho \mathbf{u} C) + \rho \alpha_C \nabla^2 C = R_{\rho C}$$

**Variables:**
- **C** — Scalar mixing ratio (dimensionless or kg/kg), e.g., water vapor, cloud water, rain, or any passive tracer
- **ρC** — Scalar density (kg/m³)
- **α_C** — Scalar diffusivity (m²/s)
- **R_ρC** — Shorthand for all RHS terms

**Physical meaning:** Scalars (moisture species, tracers, etc.) are advected by the flow and subject to diffusive mixing. This equation is solved for each scalar quantity tracked by the model.

---

## 2. Equation of State

$$p = \left( \frac{\rho R_d \theta}{p_0^{R_d / c_p}} \right)^\gamma$$

**Variables:**
- **p** — Total pressure (Pa)
- **R_d** — Dry air gas constant (≈ 287 J/(kg·K))
- **c_p** — Specific heat at constant pressure (≈ 1004 J/(kg·K))
- **γ** = c_p/c_v — Ratio of specific heats (≈ 1.4 for dry air)
- **p₀** — Reference pressure (typically 100000 Pa = 1000 hPa)

**Physical meaning:** This is a diagnostic relationship that closes the system. Given density and potential temperature, pressure can be computed. It replaces the ideal gas law (p = ρRT) using potential temperature instead of absolute temperature.

---

## 3. Pressure–Theta Relation (Used in Acoustic Substepping)

$$\nabla p = \gamma R_d \pi \nabla \Theta$$

where the **Exner function** is:

$$\pi = \left(\frac{p}{p_0}\right)^\kappa, \quad \kappa = \frac{R_d}{c_p}$$

**Variables:**
- **π** — Exner function (dimensionless), a normalized pressure variable
- **κ** — Ratio R_d/c_p (≈ 0.286)

**Physical meaning:** This relation allows the pressure gradient to be expressed in terms of the gradient of potential temperature, which is a prognostic variable. This is crucial for the acoustic substepping algorithm.

---

## 4. Momentum Equations in Perturbational Form (for Acoustic Substepping)

Using the pressure–theta relation, the momentum equations are rewritten in component form:

$$\frac{\partial U}{\partial t} = -\nabla \cdot (\mathbf{u} U) - \gamma R_d \pi \frac{\partial \Theta'}{\partial x} + F_u$$

$$\frac{\partial V}{\partial t} = -\nabla \cdot (\mathbf{u} V) - \gamma R_d \pi \frac{\partial \Theta'}{\partial y} + F_v$$

$$\frac{\partial W}{\partial t} = -\nabla \cdot (\mathbf{u} W) - \gamma R_d \pi \frac{\partial \Theta'}{\partial z} - g\left(\bar{\rho} \frac{\pi'}{\bar{\pi}} - \rho'\right) + F_w$$

**Variables:**
- **U, V, W** = ρu, ρv, ρw — Momentum components (kg/m²s)
- **Θ'** — Perturbation of ρθ from the base state
- **π'** — Perturbation Exner function (deviation from the hydrostatically balanced base state)
- **π̄** — Base-state Exner function
- **ρ̄** — Base-state density
- **ρ'** — Perturbation density
- **g** — Gravitational acceleration (≈ 9.81 m/s²)
- **F_u, F_v, F_w** — Forcing terms in each direction

**Physical meaning:** The horizontal momentum responds to advection, pressure gradient forces (expressed via the potential temperature perturbation gradient), and external forcing. The vertical momentum equation additionally includes the buoyancy term involving both the Exner function perturbation and the density perturbation.

---

## 5. Time Integration: 3rd-Order Runge-Kutta (RK3) Scheme

For the fully compressible equations, ERF uses a low-storage third-order Runge-Kutta method. Denoting the state vector as **S** and the RHS operator as f(**S**):

$$\mathbf{S}^{(1)} = \mathbf{S}^n + \frac{\Delta t}{3} f(\mathbf{S}^n)$$

$$\mathbf{S}^{(2)} = \mathbf{S}^n + \frac{\Delta t}{2} f(\mathbf{S}^{(1)})$$

$$\mathbf{S}^{n+1} = \mathbf{S}^n + \Delta t \, f(\mathbf{S}^{(2)})$$

**Variables:**
- **S^n** — State vector at time level n (includes ρ, ρu, Θ, ρC)
- **Δt** — Time step (s)
- **f(S)** — Right-hand side operator evaluating all spatial terms
- **S^(1), S^(2)** — Intermediate Runge-Kutta stages
- Superscript **n** — Previous (known) time step

**Physical meaning:** This is a three-stage explicit time-stepping method. Each stage evaluates the spatial operators at an intermediate state, progressively refining the solution to achieve third-order accuracy in time.

---

## 6. Time Integration: 2nd-Order Runge-Kutta (for Anelastic Equations)

When solving the **anelastic** (rather than fully compressible) equations, ERF uses a simpler two-stage method:

$$\mathbf{S}^{*} = \mathbf{S}^n + \Delta t \, f(\mathbf{S}^n)$$

$$\mathbf{S}^{n+1} = \mathbf{S}^n + \frac{\Delta t}{2} \left( f(\mathbf{S}^{n}) + f(\mathbf{S}^{*}) \right)$$

**Variables:**
- **S*** — Predictor state (first stage)
- All other variables as defined above

**Physical meaning:** The predictor step takes a full Euler step forward. The corrector step averages the tendencies at the old and predicted states, yielding second-order accuracy. No acoustic substepping is needed because the anelastic approximation filters out sound waves.

---

## 7. Acoustic Substepping Equations

Within each RK3 stage, the fast-propagating acoustic and gravity wave modes are substep with a smaller time step **δτ**. Perturbational quantities are defined as:

$$\mathbf{U}'' = \mathbf{U} - \mathbf{U}^t$$

where **U^t** is the momentum at the most recent "large" time step (either t^n or an intermediate RK stage).

The acoustic substepping then evolves these perturbational quantities at a finer temporal resolution, solving for **U''**, **V''**, **W''**, **ρ''**, and **Θ''** at each sub-timestep τ + δτ.

**Procedure:**
1. Evolve U'' and V'' (horizontal momentum perturbations) explicitly using the known Θ''
2. Solve a tridiagonal system for W'' (vertical momentum perturbation) — this is the vertically implicit step
3. Update ρ'' and Θ'' using the newly computed W''

**Variables:**
- **δτ** — Acoustic sub-timestep (s), much smaller than Δt
- **U'', V'', W''** — Perturbational momentum components relative to the current RK stage
- **ρ''** — Perturbational density
- **Θ''** — Perturbational conserved potential temperature
- **U^t** — Reference momentum at the current large time step

**Physical meaning:** The fast acoustic modes (sound waves) propagate much faster than the meteorologically important advective motions. By substepping only the acoustic terms, ERF can use a large main time step determined by the advective CFL condition while maintaining stability for the fast waves.

---

## 8. Divergence Damping

To control horizontally propagating sound waves, ERF applies divergence damping via a forward-projected pressure:

$$p''^{\tau*} = p''^{\tau} + \beta_d \left( p''^{\tau} - p''^{\tau - \delta\tau} \right)$$

**Variables:**
- **p''^τ** — Perturbational pressure at acoustic sub-timestep τ
- **p''^(τ−δτ)** — Perturbational pressure at the previous acoustic sub-timestep
- **p''^τ*** — Forward-projected pressure value used in the RHS of acoustic substepping equations for horizontal momentum
- **β_d** — Divergence damping coefficient (dimensionless); typical WRF value is 0.1

**Physical meaning:** This forward projection of pressure acts like a temporal filter that damps spurious acoustic oscillations in the horizontal direction. It is equivalent to adding a horizontal diffusion term in the continuity equation.

---

## 9. Acoustic Off-Centering

The implicit vertical acoustic step uses off-centering controlled by a parameter **β_s** (default value 0.1):

This biases the time average of the vertically propagating acoustic/gravity wave terms toward the future sub-timestep, providing damping of both horizontally and vertically propagating sound waves.

**Variable:**
- **β_s** — Off-centering parameter (dimensionless, controllable via `erf.beta_s`)

---

## 10. Sponge Layer Damping

ERF applies sponge source terms near domain boundaries to prevent spurious wave reflections:

$$\frac{dQ}{dt} = \text{RHS} - A \xi^n (Q - Q_{\text{target}})$$

**Variables:**
- **Q** — Any prognostic variable (ρ, ρu, Θ, ρC, etc.)
- **RHS** — All other right-hand side terms in that variable's evolution equation
- **A** — Sponge amplitude (user-specified)
- **n** — Sponge strength exponent (user-specified)
- **ξ** — Linear coordinate ranging from 0 (start of sponge zone) to 1 (end/boundary)
- **Q_target** — Target solution value in the sponge zone (user-specified)

**Physical meaning:** Within the sponge zone, variables are nudged toward a prescribed target state. The nudging strength increases from zero at the inner edge of the sponge to maximum at the domain boundary. This prevents outgoing waves from reflecting back into the domain interior.

---

## 11. LES Turbulence Models

### 11.1 Smagorinsky Model

The turbulent (eddy) viscosity is computed as:

$$\mu_t = (C_s \Delta)^2 \sqrt{2 \tilde{S}_{ij} \tilde{S}_{ij}} \; \bar{\rho}$$

And the subgrid-scale stress tensor is:

$$\tau_{ij} = -2 \mu_t \tilde{\sigma}_{ij} = -K \tilde{\sigma}_{ij}$$

**Variables:**
- **μ_t** — Turbulent (eddy) viscosity (Pa·s)
- **C_s** — Smagorinsky constant (dimensionless, typically ~0.18)
- **Δ** — Filter width, taken as the cube root of cell volume (m)
- **S̃_ij** — Filtered strain-rate tensor
- **ρ̄** — Filtered (resolved) density
- **τ_ij** — Subgrid-scale stress tensor
- **σ̃_ij** — Filtered deviatoric strain-rate tensor
- **K** = 2μ_t — Eddy diffusivity coefficient

**Physical meaning:** The Smagorinsky model parameterizes unresolved turbulent motions by relating subgrid stresses to the resolved strain rate. The eddy viscosity scales with both grid spacing and the magnitude of resolved strain, providing more diffusion where there is more shear.

### 11.2 Filtered Equation of State (LES Approximation)

$$\overline{p} \approx \left( \frac{\bar{\rho} R_d \tilde{\theta}}{p_0^{R_d / c_p}} \right)^\gamma$$

**Physical meaning:** The nonlinearity in the equation of state is weak for γ = 1.4, so the subfilter contribution to pressure is neglected. The filtered pressure is approximated using the filtered density and Favre-averaged potential temperature.

---

## 12. Diffusion Coefficients

In the evolution equations, the diffusion terms use the following coefficients:

- **μ** = ρα_M — Dynamic viscosity in the momentum equation, where α_M is a kinematic viscosity (m²/s)
- **ρα_T** — Thermal diffusion coefficient in the potential temperature equation (kg/(m·s))
- **ρα_C** — Scalar diffusion coefficient in the scalar transport equation (kg/(m·s))

In DNS mode, these are the molecular transport coefficients. In LES mode, they are replaced by the turbulent equivalents computed from the chosen LES model (Smagorinsky or Deardorff).

---

## Summary of Key Notation

| Symbol | Meaning | Units |
|--------|---------|-------|
| ρ | Air density | kg/m³ |
| u = (u, v, w) | Velocity vector | m/s |
| U = (U, V, W) = ρu | Momentum vector | kg/(m²·s) |
| θ | Potential temperature | K |
| Θ = ρθ | Conserved potential temperature | kg·K/m³ |
| p | Pressure | Pa |
| p' | Perturbation pressure | Pa |
| π | Exner function (p/p₀)^κ | dimensionless |
| C | Scalar mixing ratio | kg/kg |
| R_d | Dry air gas constant | 287 J/(kg·K) |
| c_p | Specific heat at constant pressure | 1004 J/(kg·K) |
| γ | Ratio of specific heats c_p/c_v | ~1.4 |
| κ | R_d/c_p | ~0.286 |
| p₀ | Reference pressure | 100000 Pa |
| g | Gravitational acceleration | 9.81 m/s² |
| μ | Dynamic viscosity | Pa·s |
| α_T | Thermal diffusivity | m²/s |
| α_C | Scalar diffusivity | m²/s |
| Δt | Large (main) time step | s |
| δτ | Acoustic sub-timestep | s |
| β_d | Divergence damping coefficient | dimensionless |
| β_s | Acoustic off-centering parameter | dimensionless |

---

*Sources: ERF readthedocs documentation (erf.readthedocs.io), Lattanzi et al. (2025) — "ERF: Energy Research and Forecasting Model", JAMES, 17(11), e2024MS004884.*
