#ifndef ERF_PBLMODELS_H_
#define ERF_PBLMODELS_H_

#include "ERF_DataStruct.H"

/**
 * Compute eddy diffusivities of momentum (eddy viscosity) and heat using the
 * Mellor-Yamada-Janjić (MYJ) boundary layer scheme.
 * Reference: Janjić, Z. I., 1994: The Step-Mountain Eta Coordinate Model:
 * Further Developments of the Convection, Viscous Sublayer, and Turbulence
 * Closure Schemes. Monthly Weather Review, 122, 927–945.
 * https://doi.org/10.1175/1520-0493(1994)122<0927:TSMECM>2.0.CO;2
 *
 * @param[in] xvel Velocity in x-dir
 * @param[in] yvel Velocity in y-dir
 * @param[in] cons_in Cell center conserved quantities
 * @param[out] eddyViscosity Holds turbulent viscosity
 * @param[in] geom Problem geometry
 * @param[in] turbChoice Container with turbulence parameters
 * @param[in] SurfLayer Pointer to Monin-Obukhov class if instantiated
 * @param[in] use_moisture If we have microphysics enabled
 * @param[in] level Current level
 * @param[in] bc_ptr Pointer to array with boundary condition info
 * @param[in] vert_only Only compute vertical eddy diffusivities
 * @param[in] z_phys_nd Physical location of grid nodes, if terrain (or grid stretching) is enabled
 * @param[in] z_phys_cc Physical location of cell centers
 */
void
ComputeDiffusivityMYJ (double dt,
                       const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       amrex::MultiFab& cons_in,
                       amrex::MultiFab& eddyViscosity,
                       const amrex::Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool use_moisture,
                       int level,
                       const amrex::BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                       const std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                       const MoistureComponentIndices& moisture_indices);

/**
 * Compute eddy diffusivities of momentum (eddy viscosity) and heat using the
 * Mellor-Yamada-Nakanishi-Niino Level amrex::Real(2.5) PBL scheme
 *
 * @param[in] xvel Velocity in x-dir
 * @param[in] yvel Velocity in y-dir
 * @param[in] cons_in Cell center conserved quantities
 * @param[out] eddyViscosity Holds turbulent viscosity
 * @param[in] geom Problem geometry
 * @param[in] turbChoice Container with turbulence parameters
 * @param[in] SurfLayer Pointer to Monin-Obukhov class if instantiated
 * @param[in] use_moisture If we have microphysics enabled
 * @param[in] level Current level
 * @param[in] bc_ptr Pointer to array with boundary condition info
 * @param[in] vert_only Only compute vertical eddy diffusivities
 * @param[in] z_phys_nd Physical location of grid nodes, if terrain (or grid stretching) is enabled
 * @param[in] z_phys_cc Physical location of cell centers
 */
void
ComputeDiffusivityMYNN25 (const amrex::MultiFab& xvel,
                          const amrex::MultiFab& yvel,
                          const amrex::MultiFab& cons_in,
                          amrex::MultiFab& eddyViscosity,
                          const amrex::Geometry& geom,
                          const TurbChoice& turbChoice,
                          std::unique_ptr<SurfaceLayer>& SurfLayer,
                          bool use_terrain_fitted_coords,
                          bool use_moisture,
                          int level,
                          const amrex::BCRec* bc_ptr,
                          bool /*vert_only*/,
                          const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                          const std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                          const MoistureComponentIndices& moisture_indices);

/**
 * Compute eddy diffusivities of momentum (eddy viscosity) and heat using the
 * Mellor-Yamada-Nakanishi-Niino Level amrex::Real(2.5) PBL scheme
 *
 * @param[in] xvel Velocity in x-dir
 * @param[in] yvel Velocity in y-dir
 * @param[in] cons_in Cell center conserved quantities
 * @param[out] eddyViscosity Holds turbulent viscosity
 * @param[in] geom Problem geometry
 * @param[in] turbChoice Container with turbulence parameters
 * @param[in] SurfLayer Pointer to Monin-Obukhov class if instantiated
 * @param[in] use_moisture If we have microphysics enabled
 * @param[in] level Current level
 * @param[in] bc_ptr Pointer to array with boundary condition info
 * @param[in] vert_only Only compute vertical eddy diffusivities
 * @param[in] z_phys_nd Physical location of grid nodes, if terrain (or grid stretching) is enabled
 * @param[in] z_phys_cc Physical location of cell centers
 */
void
ComputeDiffusivityMYNNEDMF (const amrex::MultiFab& xvel,
                            const amrex::MultiFab& yvel,
                            const amrex::MultiFab& cons_in,
                            amrex::MultiFab& eddyViscosity,
                            const amrex::Geometry& geom,
                            const TurbChoice& turbChoice,
                            std::unique_ptr<SurfaceLayer>& SurfLayer,
                            bool use_terrain_fitted_coords,
                            bool use_moisture,
                            int level,
                            const amrex::BCRec* bc_ptr,
                            bool /*vert_only*/,
                            const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                            const std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                            const MoistureComponentIndices& moisture_indices);

/**
 * Compute eddy diffusivities of momentum (eddy viscosity) and heat using the
 * Yonsei University PBL scheme
 *
 * @param[in] xvel Velocity in x-dir
 * @param[in] yvel Velocity in y-dir
 * @param[in] cons_in Cell center conserved quantities
 * @param[out] eddyViscosity Holds turbulent viscosity
 * @param[in] geom Problem geometry
 * @param[in] turbChoice Container with turbulence parameters
 * @param[in] SurfLayer Pointer to Monin-Obukhov class if instantiated
 * @param[in] use_moisture If we have microphysics enabled
 * @param[in] level Current level
 * @param[in] bc_ptr Pointer to array with boundary condition info
 * @param[in] vert_only Only compute vertical eddy diffusivities
 * @param[in] z_phys_nd Physical location of grid nodes, if terrain (or grid stretching) is enabled
 * @param[in] z_phys_cc Physical location of cell centers
 */
void
ComputeDiffusivityYSU (const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& cons_in,
                       amrex::MultiFab& eddyViscosity,
                       const amrex::Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool use_moisture,
                       int level,
                       const amrex::BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                       const std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                       const MoistureComponentIndices& moisture_indices);
/**
 * Compute eddy diffusivities of momentum (eddy viscosity) and heat using the
 * Medium Range Forecast (MRF) boundary layer scheme based on Hong and Pan (1996).
 * Reference: Hong, S. Y., and H.-L. Pan, 1996: Nonlocal Boundary Layer Vertical
 * Diffusion in a Medium-Range Forecast Model. Monthly Weather Review, 124,
 * 2322-2339. https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
 *
 * @param[in] xvel Velocity in x-dir
 * @param[in] yvel Velocity in y-dir
 * @param[in] cons_in Cell center conserved quantities
 * @param[out] eddyViscosity Holds turbulent viscosity
 * @param[in] geom Problem geometry
 * @param[in] turbChoice Container with turbulence parameters
 * @param[in] SurfLayer Pointer to Monin-Obukhov class if instantiated
 * @param[in] use_moisture If we have microphysics enabled
 * @param[in] level Current level
 * @param[in] bc_ptr Pointer to array with boundary condition info
 * @param[in] vert_only Only compute vertical eddy diffusivities
 * @param[in] z_phys_nd Physical location of grid nodes, if terrain (or grid stretching) is enabled
 * @param[in] z_phys_cc Physical location of cell centers
 */
void
ComputeDiffusivityMRF (const amrex::MultiFab& xvel,
                       const amrex::MultiFab& yvel,
                       const amrex::MultiFab& cons_in,
                       amrex::MultiFab& eddyViscosity,
                       const amrex::Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool use_moisture,
                       int level,
                       const amrex::BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                       const std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                       const MoistureComponentIndices& moisture_indices);

/**
 * Compute eddy diffusivities of momentum (eddy viscosity) and heat using the
 * WRF Yonsei University (YSU) PBL scheme.
 *
 * Based on:
 *   Hong, Noh &amp; Dudhia (2006), MWR, https://doi.org/10.1175/MWR3250.1
 *   Hong (2010), QJRMS, https://doi.org/10.1002/qj.665
 *
 * Features: three-pass bulk-Richardson PBLH diagnosis, nonlocal
 * countergradient fluxes, explicit entrainment at PBL top, and
 * grid-adaptive Richardson-number mixing in free atmosphere.
 *
 * @param[in]  xvel   x-direction velocity
 * @param[in]  yvel   y-direction velocity
 * @param[in]  cons_in  cell-center conserved quantities
 * @param[out] eddyViscosity  turbulent eddy viscosity
 * @param[in]  geom   problem geometry
 * @param[in]  turbChoice  turbulence parameter container
 * @param[in]  SurfLayer  Monin-Obukhov surface layer object
 * @param[in]  use_terrain_fitted_coords  flag for terrain coordinates
 * @param[in]  use_moisture  flag for moisture active
 * @param[in]  level  AMR level
 * @param[in]  bc_ptr  boundary condition array
 * @param[in]  vert_only  compute vertical diffusivity only
 * @param[in]  z_phys_nd  physical node heights (terrain)
 * @param[in]  z_phys_cc  physical cell-center heights
 * @param[in]  moisture_indices  component index mapping for moisture variables
 * @param[in]  qheating_rates  radiation heating rates (SW, LW components); optional
 */
void
ComputeDiffusivityYSUNew (const amrex::MultiFab& xvel,
                          const amrex::MultiFab& yvel,
                          const amrex::MultiFab& cons_in,
                          amrex::MultiFab& eddyViscosity,
                          const amrex::Geometry& geom,
                          const TurbChoice& turbChoice,
                          std::unique_ptr<SurfaceLayer>& SurfLayer,
                          bool use_terrain_fitted_coords,
                          bool use_moisture,
                          int level,
                          const amrex::BCRec* bc_ptr,
                          bool /*vert_only*/,
                          const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                          const std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                          const MoistureComponentIndices& moisture_indices,
                          const amrex::MultiFab* qheating_rates = nullptr);

/**
 * Function for computing vertical derivatives for use in PBL model
 *
 * @param[in] u velocity in x-dir
 * @param[in] v velocity in y-dir
 * @param[in] cell_data conserved cell center vars
 */
struct PBLDerivativeDzInv
{
    // Each value is the inverse physical distance between the two samples
    // used by its corresponding derivative stencil.
    amrex::Real interior;
    amrex::Real dirichlet_lo;
    amrex::Real dirichlet_hi;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ComputeVerticalDerivativesPBL (int i, int j, int k,
                               const amrex::Array4<const amrex::Real>& uvel,
                               const amrex::Array4<const amrex::Real>& vvel,
                               const amrex::Array4<const amrex::Real>& cell_data,
                               const int izmin,
                               const int izmax,
                               const PBLDerivativeDzInv& dz_inv,
                               const bool c_ext_dir_on_zlo,
                               const bool c_ext_dir_on_zhi,
                               const bool u_ext_dir_on_zlo,
                               const bool u_ext_dir_on_zhi,
                               const bool v_ext_dir_on_zlo,
                               const bool v_ext_dir_on_zhi,
                               amrex::Real& dthetadz,
                               amrex::Real& dudz,
                               amrex::Real& dvdz,
                               const MoistureComponentIndices& moisture_indices)
{
    // ext_dir ghost values hold the value at the physical boundary, which is
    // half a ghost-center spacing from the adjacent cell center. Use a
    // two-point, first-order boundary-to-center difference with inverse
    // distance 2/dz.
    if ( k==izmax && c_ext_dir_on_zhi ) {
        dthetadz = ( GetThetav(i,j,k+1,cell_data,moisture_indices)
                    -GetThetav(i,j,k  ,cell_data,moisture_indices) )*dz_inv.dirichlet_hi;
    } else if ( k==izmin && c_ext_dir_on_zlo ) {
        dthetadz = ( GetThetav(i,j,k  ,cell_data,moisture_indices)
                    -GetThetav(i,j,k-1,cell_data,moisture_indices) )*dz_inv.dirichlet_lo;
    } else {
        dthetadz = ( GetThetav(i,j,k+1,cell_data,moisture_indices)
                    -GetThetav(i,j,k-1,cell_data,moisture_indices) )*dz_inv.interior;
    }

    if ( k==izmax && u_ext_dir_on_zhi ) {
        dudz = myhalf*( uvel(i  ,j,k+1) - uvel(i  ,j,k)
                       +uvel(i+1,j,k+1) - uvel(i+1,j,k) )*dz_inv.dirichlet_hi;
    } else if ( k==izmin && u_ext_dir_on_zlo ) {
        dudz = myhalf*( uvel(i  ,j,k) - uvel(i  ,j,k-1)
                       +uvel(i+1,j,k) - uvel(i+1,j,k-1) )*dz_inv.dirichlet_lo;
    } else {
        dudz = myhalf*( uvel(i  ,j,k+1) - uvel(i  ,j,k-1)
                       +uvel(i+1,j,k+1) - uvel(i+1,j,k-1) )*dz_inv.interior;
    }

    if ( k==izmax && v_ext_dir_on_zhi ) {
        dvdz = myhalf*( vvel(i,j  ,k+1) - vvel(i,j  ,k)
                       +vvel(i,j+1,k+1) - vvel(i,j+1,k) )*dz_inv.dirichlet_hi;
    } else if ( k==izmin && v_ext_dir_on_zlo ) {
        dvdz = myhalf*( vvel(i,j  ,k) - vvel(i,j  ,k-1)
                       +vvel(i,j+1,k) - vvel(i,j+1,k-1) )*dz_inv.dirichlet_lo;
    } else {
        dvdz = myhalf*( vvel(i,j  ,k+1) - vvel(i,j  ,k-1)
                       +vvel(i,j+1,k+1) - vvel(i,j+1,k-1) )*dz_inv.interior;
    }
}

/**
 * Function for computing the QKE source terms (NN09, Eqn. 5).
 *
 * @param[in] u velocity in x-dir
 * @param[in] v velocity in y-dir
 * @param[in] cell_data conserved cell center vars
 * @param[in] cell_prim primitive cell center vars
 * @param[in] K_turb turbulent viscosity
 * @param[in] domain box of the whole domain
 * @param[in] dz_inv inverse spacings for interior and Dirichlet-boundary differences
 * @param[in] pbl_mynn_B1_l a parameter
 * @param[in] theta_mean average theta
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeQKESourceTerms (int i, int j, int k,
                       const amrex::Array4<const amrex::Real>& uvel,
                       const amrex::Array4<const amrex::Real>& vvel,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& K_turb,
                       const amrex::Box& domain,
                       const PBLDerivativeDzInv dz_inv,
                       amrex::Real pbl_mynn_B1_l,
                       const amrex::Real theta_mean,
                       const MoistureComponentIndices& moisture_indices,
                       bool c_ext_dir_on_zlo,
                       bool c_ext_dir_on_zhi,
                       bool u_ext_dir_on_zlo,
                       bool u_ext_dir_on_zhi,
                       bool v_ext_dir_on_zlo,
                       bool v_ext_dir_on_zhi)
{
    // Compute some relevant derivatives
    amrex::Real dthetadz, dudz, dvdz;
    amrex::Real source_term = zero;

    int izmin = domain.smallEnd(2);
    int izmax = domain.bigEnd(2);

    ComputeVerticalDerivativesPBL(i, j, k,
                                  uvel, vvel, cell_data, izmin, izmax, dz_inv,
                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                  u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                  v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                  dthetadz, dudz, dvdz,
                                  moisture_indices);

    // Notes:
    // - We transport TKE = myhalf*QKE rather than QKE, so the RHS terms do not
    //   have a factor of two.
    // - Transport terms due to turbulence and pressure are included when
    //   DiffusionSrcForState_* is called from ERF_slow_rhs_post.
    // - Eddy diffusivities are updated at the beginning of each time step only.

    // Second-order turbulent fluxes, e.g.:
    // -<uw> = L q SM dU/dz  (NN09, Eqn. 18)
    //       = Kmv/rho dU/dz

    // Shear Production
    source_term += K_turb(i,j,k,EddyDiff::Mom_v) * (dudz*dudz + dvdz*dvdz);

    // Buoyancy Production
    source_term -= (CONST_GRAV/theta_mean)*K_turb(i,j,k,EddyDiff::Theta_v)*dthetadz;

    // Dissipation (NN09, Eqn. 12)
    amrex::Real qke = two * cell_prim(i,j,k,PrimKE_comp);
    if (std::abs(qke) > zero) {
        source_term -= cell_data(i,j,k,Rho_comp) * std::pow(qke,amrex::Real(1.5)) /
                       (pbl_mynn_B1_l * K_turb(i,j,k,EddyDiff::Turb_lengthscale));
    }

    return source_term;
}

struct PBLDerivativeDzInv_N
{
    amrex::Real value;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    PBLDerivativeDzInv operator() (int, int, int) const noexcept
    {
        return {myhalf*value, two*value, two*value};
    }
};

struct PBLDerivativeDzInv_S
{
    const amrex::Real* dz;
    int klo;
    int khi;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    PBLDerivativeDzInv operator() (int, int, int k) const noexcept
    {
        // Adjacent cell centers are separated by half the sum of their widths.
        const amrex::Real dz_lo = dz[(k > klo) ? k-1 : k];
        const amrex::Real dz_hi = dz[(k < khi) ? k+1 : k];
        return {two / (dz_lo + two*dz[k] + dz_hi),
                two / dz[k],
                two / dz[k]};
    }
};

struct PBLDerivativeDzInv_T
{
    amrex::Array4<const amrex::Real> z_cc;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    PBLDerivativeDzInv operator() (int i, int j, int k) const noexcept
    {
        return {one / (z_cc(i,j,k+1) - z_cc(i,j,k-1)),
                two / (z_cc(i,j,k  ) - z_cc(i,j,k-1)),
                two / (z_cc(i,j,k+1) - z_cc(i,j,k  ))};
    }
};
#endif

/**
 * @brief Apply spatial smoothing to PBLH field using 5-point stencil.
 *
 * Applies a simple 5-point spatial averaging to reduce grid-to-grid noise in
 * the diagnosed PBLH field. This is standard post-processing practice for
 * discrete Rib-crossing detection as described in Seibert et al. (2000):
 * "Review and intercomparison of operational methods for the determination of
 * the mixing height", Atmospheric Environment, 34, 1001–1027.
 *
 * Stencil: PBLH_smooth(i,j) = w * PBLH(i,j) + (1-w)/4 * [PBLH(i±1,j) + PBLH(i,j±1)]
 *
 * @param[in,out] pblh_fab FArrayBox containing 2D PBLH field (component 0)
 * @param[in] xybx Valid box for smoothing (x,y indices only)
 * @param[in] weight Center cell weight in stencil (must be in [0,1])
 * @param[in] passes Number of smoothing iterations
 * @param[in] domain Domain box for reflective boundary conditions
 */
inline void
ApplyPBLHSmoothing (amrex::FArrayBox& pblh_fab,
                    const amrex::Box& xybx,
                    const amrex::Real weight,
                    const int passes,
                    const amrex::Box& domain)
{
    // Use a temporary to hold intermediate results
    amrex::FArrayBox pblh_temp(xybx, 1, amrex::The_Async_Arena());
    auto pblh = pblh_fab.array();
    auto pblh_tmp = pblh_temp.array();

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    const amrex::Real wt_side = (one - weight) / four; // Weight for each of 4 neighbors

    // Apply smoothing passes
    for (int pass = 0; pass < passes; ++pass) {
        amrex::ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            // Use reflective (Neumann) boundary conditions:
            // At domain boundaries, reuse the current cell's own value
            const int i_xlo = (i > dom_lo.x) ? (i - 1) : i;
            const int i_xhi = (i < dom_hi.x) ? (i + 1) : i;
            const int j_ylo = (j > dom_lo.y) ? (j - 1) : j;
            const int j_yhi = (j < dom_hi.y) ? (j + 1) : j;

            // 5-point stencil
            pblh_tmp(i, j, 0) = weight * pblh(i, j, 0)
                              + wt_side * (pblh(i_xlo, j, 0) + pblh(i_xhi, j, 0) +
                                          pblh(i, j_ylo, 0) + pblh(i, j_yhi, 0));
        });

        // Swap arrays for next iteration (or final result)
        if (pass < passes - 1) {
            pblh_fab.copy(pblh_temp, 0, 0, 1);
        } else {
            // Copy final result back
            pblh_fab.copy(pblh_temp, 0, 0, 1);
        }
    }
}