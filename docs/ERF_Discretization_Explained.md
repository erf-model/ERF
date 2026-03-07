# ERF Spatial Discretization — Step-by-Step Walkthrough

## 1. The Arakawa C-Grid: What Lives Where

ERF uses the classic **Arakawa C-grid** staggering. The key idea is that not everything is stored at the same location in a grid cell:

```
              y
              ^
              |
         +---------+-- v(i,j+1,k) --+---------+
         |         |                 |         |
         |         |                 |         |
    u(i,j,k)    [ρ, Θ, C, p]   u(i+1,j,k)
         |       (i,j,k)            |         |
         |     CELL CENTER           |         |
         +---------+-- v(i,j,k) ----+---------+  --> x
                   |
                   w(i,j,k) is on bottom z-face
                   w(i,j,k+1) is on top z-face
```

### Cell-centered quantities (defined at the volumetric center of each cell):
- **ρ** — density
- **Θ = ρθ** — conserved potential temperature
- **ρC** — scalar densities (moisture, tracers, etc.)
- **p** — pressure (diagnosed from equation of state)

### Face-centered quantities (defined at cell faces, normal to each respective direction):
- **u** — x-velocity, stored on the **x-faces** (i±½, j, k)
- **v** — y-velocity, stored on the **y-faces** (i, j±½, k)
- **w** — z-velocity, stored on the **z-faces** (i, j, k±½)

### Why this matters:
The staggering means that when you compute a derivative, you naturally get it at the *other* location. For example, ∂u/∂x evaluated between u(i+½) and u(i−½) lands exactly at the cell center (i) — where the scalars live. Conversely, ∂ρ/∂x evaluated between ρ(i) and ρ(i+1) lands at the face (i+½) — where u lives. This doubles the effective resolution compared to an unstaggered grid.

---

## 2. Every Equation Has the Same Flux-Divergence Structure

You're exactly right — every prognostic equation in ERF can be written in a single unifying structure:

$$\frac{\partial \phi}{\partial t} = -\frac{\partial (u \phi)}{\partial x} - \frac{\partial (v \phi)}{\partial y} - \frac{\partial (w \phi)}{\partial z} + \text{other terms}$$

where the advective fluxes are always of the form **velocity × transported quantity**:

| Equation | φ (transported) | x-flux | y-flux | z-flux |
|----------|----------------|--------|--------|--------|
| Continuity | ρ | ρu | ρv | ρw |
| x-Momentum | U = ρu | u·(ρu) | v·(ρu) | w·(ρu) |
| y-Momentum | V = ρv | u·(ρv) | v·(ρv) | w·(ρv) |
| z-Momentum | W = ρw | u·(ρw) | v·(ρw) | w·(ρw) |
| Pot. Temp. | Θ = ρθ | u·(ρθ) | v·(ρθ) | w·(ρθ) |
| Scalar | ρC | u·(ρC) | v·(ρC) | w·(ρC) |

Every flux is: **(advecting velocity) × (advected quantity)**. The advecting velocity is always at a face, and the advected quantity must be interpolated to that same face.

In semi-discrete form (continuous in time, discrete in space), for any cell-centered quantity φ at cell (i,j,k):

$$\frac{\partial \phi_{i,j,k}}{\partial t} = -\frac{F^x_{i+\frac{1}{2}} - F^x_{i-\frac{1}{2}}}{\Delta x} - \frac{F^y_{j+\frac{1}{2}} - F^y_{j-\frac{1}{2}}}{\Delta y} - \frac{F^z_{k+\frac{1}{2}} - F^z_{k-\frac{1}{2}}}{\Delta z} + \ldots$$

where F^x, F^y, F^z are the numerical fluxes at the cell faces.

---

## 3. The Critical Question: How Do You Compute the Flux at a Face?

This is where the order of the advection scheme comes in. The advecting velocity (say u) is *already* defined at the x-face — that's the C-grid staggering doing its job. But the advected quantity (say ρθ) is defined at cell centers. You need to **interpolate** ρθ to the face to form the flux F = u · (ρθ)_face.

The choice of interpolation stencil determines the order and character of the scheme.

---

## 4. The Wicker–Skamarock Flux Formulas (WS5 / WS6)

ERF (following WRF) uses the advection schemes from Wicker & Skamarock (2002). The key insight is:

> **The 5th-order upwind flux = 6th-order centered flux + a numerical dissipation term**

### 6th-order centered flux (WS6):

$$F^{6}_{i-\frac{1}{2}} = \frac{u_{i-\frac{1}{2}}}{60} \Big[ 37(\psi_i + \psi_{i-1}) - 8(\psi_{i+1} + \psi_{i-2}) + (\psi_{i+2} + \psi_{i-3}) \Big]$$

This is a **symmetric** interpolation — it uses 3 points on each side (6-point stencil). It's purely non-dissipative (no numerical diffusion), which means it preserves energy but can produce spurious oscillations.

### 5th-order upwind flux (WS5):

$$F^{5}_{i-\frac{1}{2}} = F^{6}_{i-\frac{1}{2}} - \frac{|u_{i-\frac{1}{2}}|}{60} \Big[ 10(\psi_i - \psi_{i-1}) - 5(\psi_{i+1} - \psi_{i-2}) + (\psi_{i+2} - \psi_{i-3}) \Big]$$

The second term is the **numerical dissipation**. Notice:
- It uses **|u|** (absolute value) — this ensures it always damps, regardless of flow direction
- It's constructed from **differences** (antisymmetric), not sums — it acts like a high-order diffusion
- It selectively damps the shortest, worst-resolved waves (near 2Δx–8Δx) while barely affecting long waves

### The "Blended_5th6th" option in ERF:

ERF provides a `upw_frac` parameter (between 0 and 1) that lets you blend between WS6 and WS5:

$$F^{\text{blended}}_{i-\frac{1}{2}} = F^{6}_{i-\frac{1}{2}} - \texttt{upw\_frac} \cdot \frac{|u_{i-\frac{1}{2}}|}{60} \Big[ 10(\psi_i - \psi_{i-1}) - 5(\psi_{i+1} - \psi_{i-2}) + (\psi_{i+2} - \psi_{i-3}) \Big]$$

- `upw_frac = 0` → pure 6th-order centered (no dissipation, energy-conserving)
- `upw_frac = 1` → full 5th-order upwind (stable, some numerical diffusion)
- `upw_frac = 0.5` → 50/50 blend (reduced dissipation, still stabilized)

This is the control knob for trading off **stability** (upwind) vs. **energy conservation** (centered).

---

## 5. Step-by-Step: Discretizing the Potential Temperature Equation

Let's walk through the Θ = ρθ equation as a concrete example:

$$\frac{\partial \Theta}{\partial t} = -\nabla \cdot (\mathbf{u} \Theta) + \rho \alpha_T \nabla^2 \theta$$

### Step 5a: Identify where Θ lives → **cell center** (i,j,k)

### Step 5b: Compute the x-flux at x-faces

At face (i+½, j, k):
- **u** is already defined here (C-grid!) ✓
- **Θ** needs interpolation from cell centers to this face

Using WS5:
$$F^x_{i+\frac{1}{2}} = \frac{u_{i+\frac{1}{2}}}{60}\Big[37(\Theta_{i+1} + \Theta_i) - 8(\Theta_{i+2} + \Theta_{i-1}) + (\Theta_{i+3} + \Theta_{i-2})\Big]$$
$$\quad - \frac{|u_{i+\frac{1}{2}}|}{60}\Big[10(\Theta_{i+1} - \Theta_i) - 5(\Theta_{i+2} - \Theta_{i-1}) + (\Theta_{i+3} - \Theta_{i-2})\Big]$$

### Step 5c: Compute y-flux and z-flux analogously

At face (i, j+½, k): use **v** (already there) and interpolate Θ in the j-direction with the same stencil.

At face (i, j, k+½): use **w** (already there) and interpolate Θ in the k-direction.

### Step 5d: Form the divergence

$$\left(\nabla \cdot (\mathbf{u}\Theta)\right)_{i,j,k} = \frac{F^x_{i+\frac{1}{2}} - F^x_{i-\frac{1}{2}}}{\Delta x} + \frac{F^y_{j+\frac{1}{2}} - F^y_{j-\frac{1}{2}}}{\Delta y} + \frac{F^z_{k+\frac{1}{2}} - F^z_{k-\frac{1}{2}}}{\Delta z}$$

This lives at the cell center (i,j,k), exactly where Θ is stored. ✓

### Step 5e: Add the diffusion term

The diffusion ρα_T ∇²θ is computed using standard second-order centered differences of θ at cell centers, yielding a result also at cell centers.

### Step 5f: Advance in time with RK3

Feed the total tendency into the Runge-Kutta time stepper.

---

## 6. Step-by-Step: Discretizing the x-Momentum Equation

The momentum equations are trickier because **U = ρu lives on x-faces**, not cell centers.

$$\frac{\partial U}{\partial t} = -\nabla \cdot (\mathbf{u} \, U) - \gamma R_d \pi \frac{\partial \Theta'}{\partial x} + F_U$$

### Step 6a: Identify where U lives → **x-face** (i+½, j, k)

### Step 6b: Compute the x-flux of U (advection of U by u)

The flux u·U must be evaluated at cell centers (i,j,k) — because that's where the x-faces of the "momentum control volume" are:
- **u** must be interpolated from x-faces to cell centers: u_c = ½(u_{i-½} + u_{i+½})
- **U = ρu** must also be interpolated to cell centers via the WS5/6 stencil
- Flux: F^x = u_c · (ρu)_{interpolated}

### Step 6c: Compute the y-flux of U (advection of U by v)

This flux lives on y-faces of the momentum cell, which corresponds to cell edges/corners:
- **v** is already on y-faces but at (i, j+½, k) — it needs averaging to (i+½, j+½, k)
- **U** is interpolated in the y-direction using the WS5/6 stencil

### Step 6d: Compute the z-flux of U (advection of U by w)

Similarly, w is on z-faces at (i, j, k+½) and needs averaging to (i+½, j, k+½).

### Step 6e: Pressure gradient term

$$-\gamma R_d \pi \frac{\partial \Theta'}{\partial x} \bigg|_{i+\frac{1}{2}} \approx -\gamma R_d \, \pi_{i+\frac{1}{2}} \, \frac{\Theta'_{i+1} - \Theta'_i}{\Delta x}$$

This naturally lives on the x-face because we're differencing Θ' between two adjacent cell centers. The Exner function π is interpolated to the face.

### Step 6f: Buoyancy (vertical momentum only)

Only the W equation gets the buoyancy term involving g, ρ', and π'. This lives on z-faces.

---

## 7. Summary of the Staggering Logic

| What you're updating | Where it lives | Advective fluxes evaluated at |
|---------------------|---------------|------------------------------|
| ρ (density) | Cell center | x-face, y-face, z-face |
| Θ = ρθ | Cell center | x-face, y-face, z-face |
| ρC (scalars) | Cell center | x-face, y-face, z-face |
| U = ρu | x-face | Cell center (x), edge (y), edge (z) |
| V = ρv | y-face | edge (x), Cell center (y), edge (z) |
| W = ρw | z-face | edge (x), edge (y), Cell center (z) |

The pattern is: each momentum component's advective flux in its *own* direction lands at cell centers, while cross-direction fluxes land at cell edges/corners.

---

## 8. Near Boundaries: Order Degradation

The WS5/6 stencil needs 3 cells on each side (7-point stencil). Near boundaries, walls, or terrain, ERF progressively degrades the scheme:

- **Interior**: 5th/6th order (7-point stencil, 3 ghost cells)
- **1 cell from boundary**: 3rd/4th order (5-point stencil)
- **Right at boundary**: 1st order upwind (2-point stencil)

This degradation is controlled by bit-flag arrays that mark each cell's advection order independently in each direction.

---

## 9. Available Schemes in ERF (for reference)

| Scheme name | Type | Order | Stencil width |
|-------------|------|-------|---------------|
| Centered_2nd | Centered | 2 | 3 points |
| Upwind_3rd | Upwind-biased | 3 | 5 points |
| Centered_4th | Centered | 4 | 5 points |
| Upwind_5th | Upwind-biased | 5 | 7 points |
| Centered_6th | Centered | 6 | 7 points |
| Blended_3rd4th | Blend of above | 3–4 | 5 points |
| Blended_5th6th | Blend of above | 5–6 | 7 points |
| WENO3/5/7 | Nonlinear (WENO-JS) | 3/5/7 | varies |
| WENOZ3/5/7 | Nonlinear (WENO-Z) | 3/5/7 | varies |

The **Blended_5th6th** scheme is the hybrid you asked about: it uses the `upw_frac` parameter to continuously interpolate between the non-dissipative 6th-order centered flux and the stabilizing 5th-order upwind flux.

---

*Sources: ERF readthedocs (Inputs, Discretizations, ERFvsWRF, BoundaryConditions pages); Wicker & Skamarock (2002), Mon. Wea. Rev. 130, 2088–2097; PALM model documentation; Lattanzi et al. (2025), JAMES.*
