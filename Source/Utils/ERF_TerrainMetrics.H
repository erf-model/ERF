#ifndef ERF_TERRAIN_METRIC_H_
#define ERF_TERRAIN_METRIC_H_

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <ERF_IndexDefines.H>

/**
 * Routine to define default z_phys_nd and z_phys_cc
 */
void
init_default_zphys (int lev, const amrex::Geometry& geom,
                    amrex::MultiFab& z_phys_nd, amrex::MultiFab& z_phys_cc);

/**
 * Utility routines for constructing terrain metric terms
 */

// Declare functions for ERF.cpp
void
init_zlevels (amrex::Vector<amrex::Vector<amrex::Real>>& zlevels_stag,
              amrex::Vector<amrex::Vector<amrex::Real>>& stretched_dz_h,
              amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>>& stretched_dz_d,
              amrex::Vector<amrex::Geometry> const& geom,
              amrex::Vector<amrex::IntVect>  const& ref_ratio,
              const amrex::Real grid_stretching_ratio,
              const amrex::Real zsurf,
              const amrex::Real dz0);

void
make_terrain_fitted_coords (int lev, const amrex::Geometry& geom,
                            amrex::MultiFab& z_phys_nd,
                            amrex::Vector<amrex::Real> const& z_levels_h,
                            amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>& phys_bc_type);

void
init_which_terrain_grid (int lev, const amrex::Geometry& geom,
                         amrex::MultiFab& z_phys_nd,
                         amrex::Vector<amrex::Real> const& z_levels_h);

// Compute the min dz at cell center with terrain
amrex::Real
get_dzmin_terrain (amrex::MultiFab& z_phys_nd);

//*****************************************************************************************
// Compute terrain metric terms at cell-center
//*****************************************************************************************
// Metric is at cell center
AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
amrex::Real
Compute_h_zeta_AtCellCenter (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.25 * dzInv *
      ( z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1) + z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
       -z_nd(i+1,j,k  ) - z_nd(i+1,j+1,k  ) - z_nd(i,j,k  ) - z_nd(i,j+1,k  ) );
    return met_h_zeta;
}

// Metric is at cell center
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtCellCenter (const int &i, const int &j, const int &k,
                           const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                           const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi   = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j+1,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i  ,j,k) - z_nd(i  ,j+1,k) - z_nd(i  ,j,k+1) - z_nd(i  ,j+1,k+1) );
    return met_h_xi;
}

// Metric is at cell center
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtCellCenter (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta  = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i+1,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
       -z_nd(i,j  ,k) - z_nd(i+1,j  ,k) - z_nd(i,j  ,k+1) - z_nd(i+1,j  ,k+1) );
    return met_h_eta;
}


//*****************************************************************************************
// Compute terrain metric terms at face-centers
//*****************************************************************************************
// Metric coincides with U location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtIface (const int &i, const int &j, const int &k,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_zeta = 0.5 * cellSizeInv[2] * ( (z_nd(i,j  ,k+1) - z_nd(i,j  ,k))
                                                    + (z_nd(i,j+1,k+1) - z_nd(i,j+1,k)) );
    return met_h_zeta;
}

// Metric coincides with U location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtIface (const int &i, const int &j, const int &k,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_xi   = 0.125 * cellSizeInv[0] * (
      (z_nd(i+1,j  ,k  ) - z_nd(i-1,j  ,k  )) + (z_nd(i+1,j  ,k+1) - z_nd(i-1,j  ,k+1)) +
      (z_nd(i+1,j+1,k  ) - z_nd(i-1,j+1,k  )) + (z_nd(i+1,j+1,k+1) - z_nd(i-1,j+1,k+1)) );
    return met_h_xi;
}

// Metric coincides with U location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtIface (const int &i, const int &j, const int &k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_eta  = 0.5 * cellSizeInv[1] * ( (z_nd(i,j+1,k  ) - z_nd(i,j,k  ))
                                                    + (z_nd(i,j+1,k+1) - z_nd(i,j,k+1)) );
    return met_h_eta;
}

// Metric coincides with V location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtJface (const int &i, const int &j, const int &k,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_zeta = 0.5 * cellSizeInv[2] * ( (z_nd(i  ,j,k+1) - z_nd(i  ,j,k  ))
                                                    + (z_nd(i+1,j,k+1) - z_nd(i+1,j,k  )) );
    return met_h_zeta;
}

// Metric coincides with V location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtJface (const int &i, const int &j, const int &k,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_xi   = 0.5 * cellSizeInv[0] * ( (z_nd(i+1,j,k  ) - z_nd(i,j,k  ))
                                                    + (z_nd(i+1,j,k+1) - z_nd(i,j,k+1)) );
    return met_h_xi;
}

// Metric coincides with V location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtJface (const int &i, const int &j, const int &k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_eta  = 0.125 * cellSizeInv[1] *
      ( (z_nd(i  ,j+1,k  ) - z_nd(i  ,j-1,k  )) + (z_nd(i  ,j+1,k+1) - z_nd(i  ,j-1,k+1)) +
        (z_nd(i+1,j+1,k  ) - z_nd(i+1,j-1,k  )) + (z_nd(i+1,j+1,k+1) - z_nd(i+1,j-1,k+1)) );
    return met_h_eta;
}

// Metric coincides with K location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtKface (const int &i, const int &j, const int &k,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_zeta = 0.125 * cellSizeInv[2] *
      ( (z_nd(i  ,j  ,k+1) - z_nd(i  ,j  ,k-1)) + (z_nd(i+1,j  ,k+1) - z_nd(i+1,j  ,k-1)) +
        (z_nd(i  ,j+1,k+1) - z_nd(  i,j+1,k-1)) + (z_nd(i+1,j+1,k+1) - z_nd(i+1,j+1,k-1)) );
    return met_h_zeta;
}

// Metric coincides with K location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtKface (const int &i, const int &j, const int &k,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_xi   = 0.5 * cellSizeInv[0] * ( (z_nd(i+1,j  ,k) - z_nd(i,j  ,k))
                                                    + (z_nd(i+1,j+1,k) - z_nd(i,j+1,k)) );
    return met_h_xi;
}

// Metric coincides with K location
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtKface (const int &i, const int &j, const int &k,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real met_h_eta  = 0.5 * cellSizeInv[1] * ( (z_nd(i  ,j+1,k) - z_nd(i  ,j,k))
                                                    + (z_nd(i+1,j+1,k) - z_nd(i+1,j,k)) );
    return met_h_eta;
}


//*****************************************************************************************
// Compute terrain metric terms at edge-centers
//*****************************************************************************************
// -- EdgeCenterK --

// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtEdgeCenterK (const int &i, const int &j, const int &k,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = dzInv * (z_nd(i,j,k+1) - z_nd(i,j,k));
    return met_h_zeta;
}

// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtEdgeCenterK (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi  = 0.25 * dxInv *
      ( z_nd(i+1,j,k) + z_nd(i+1,j,k+1)
       -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) );
    return met_h_xi;
}

// Metric is at edge and center Z (red pentagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtEdgeCenterK (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta = 0.25 * dyInv *
      ( z_nd(i,j+1,k) + z_nd(i,j+1,k+1)
       -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) );
    return met_h_eta;
}

// -- EdgeCenterJ --

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtEdgeCenterJ (const int &i, const int &j, const int &k,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1)
                                             -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
    return met_h_zeta;
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtEdgeCenterJ (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi = 0.25 * dxInv *
      ( z_nd(i+1,j+1,k) + z_nd(i+1,j  ,k)
       -z_nd(i-1,j+1,k) - z_nd(i-1,j  ,k) );
    return met_h_xi;
}

// Metric is at edge and center Y (magenta cross)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtEdgeCenterJ (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta = dyInv * ( z_nd(i,j+1,k) - z_nd(i,j,k) );
    return met_h_eta;
}

// -- EdgeCenterI --

// Metric is at edge and center X (purple hexagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_zeta_AtEdgeCenterI (const int &i, const int &j, const int &k,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                              const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dzInv = cellSizeInv[2];
    amrex::Real met_h_zeta = 0.25 * dzInv * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1)
                                             -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
    return met_h_zeta;
}

// Metric is at edge and center X (purple hexagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_xi_AtEdgeCenterI (const int &i, const int &j, const int &k,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dxInv = cellSizeInv[0];
    amrex::Real met_h_xi  = dxInv * ( z_nd(i+1,j,k) - z_nd(i,j,k) );
    return met_h_xi;
}

// Metric is at edge and center X (purple hexagon)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_h_eta_AtEdgeCenterI (const int &i, const int &j, const int &k,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& z_nd)
{
    amrex::Real dyInv = cellSizeInv[1];
    amrex::Real met_h_eta = 0.25 * dyInv *
      ( z_nd(i+1,j+1,k) + z_nd(i,j+1,k)
       -z_nd(i+1,j-1,k) - z_nd(i,j-1,k) );
    return met_h_eta;
}

// Relative height above terrain surface at cell center from z_nd (nodal absolute height)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_Z_AtCellCenter (const int &i, const int &j, const int &k,
                        const amrex::Array4<const amrex::Real>& z_nd)
{
    const amrex::Real z_cc = 0.125*(  z_nd(i  ,j  ,k  ) + z_nd(i  ,j  ,k+1) +
                                    + z_nd(i+1,j  ,k  ) + z_nd(i  ,j  ,k+1)
                                    + z_nd(i  ,j+1,k  ) + z_nd(i  ,j+1,k+1)
                                    + z_nd(i+1,j+1,k  ) + z_nd(i  ,j+1,k+1));

    return z_cc;
}

// Relative height above terrain surface at cell center from z_nd (nodal absolute height)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_Z_AtWFace (const int &i, const int &j, const int &k,
                   const amrex::Array4<const amrex::Real>& z_nd)
{
    const amrex::Real z_wf = 0.25*( z_nd(i  ,j  ,k  ) + z_nd(i+1,j  ,k  )
                                  + z_nd(i  ,j+1,k  ) + z_nd(i+1,j+1,k  ) );

    return z_wf;
}

// Relative height above terrain surface at cell center from z_nd (nodal absolute height)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_Zrel_AtCellCenter (const int &i, const int &j, const int &k,
                           const amrex::Array4<const amrex::Real>& z_nd)
{
    const amrex::Real z_cc = 0.125*(  z_nd(i  ,j  ,k  ) + z_nd(i  ,j  ,k+1) +
                                    + z_nd(i+1,j  ,k  ) + z_nd(i+1,j  ,k+1)
                                    + z_nd(i  ,j+1,k  ) + z_nd(i  ,j+1,k+1)
                                    + z_nd(i+1,j+1,k  ) + z_nd(i+1,j+1,k+1));

    // Note: we assume the z_nd array spans from the bottom to top of the domain
    // i.e. no domain decomposition across processors in vertical direction
    const amrex::Real z0_cc = 0.25*(  z_nd(i  ,j  ,0) + z_nd(i  ,j+1,0)
                                    + z_nd(i+1,j  ,0) + z_nd(i+1,j+1,0));

    return (z_cc - z0_cc);
}


//*****************************************************************************************
// Contravariant velocity computation
//*****************************************************************************************

// Define omega given u,v and w (4th order)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
OmegaFromW (int& i, int& j, int& k, amrex::Real w,
            const amrex::Array4<const amrex::Real>& u_arr,
            const amrex::Array4<const amrex::Real>& v_arr,
            const amrex::Array4<const amrex::Real>& mf_u,
            const amrex::Array4<const amrex::Real>& mf_v,
            const amrex::Array4<const amrex::Real>& z_nd,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // This is dh/dxi at hi and lo edges
    amrex::Real z_x_p2  = Compute_Z_AtWFace(i+2,j,k,z_nd);
    amrex::Real z_x_p1  = Compute_Z_AtWFace(i+1,j,k,z_nd);
    amrex::Real z_x_m1  = Compute_Z_AtWFace(i-1,j,k,z_nd);
    amrex::Real z_x_m2  = Compute_Z_AtWFace(i-2,j,k,z_nd);
    amrex::Real met_xi  = (amrex::Real(1.0/12.0)*(z_x_m2 - z_x_p2)
                         + amrex::Real(8.0/12.0)*(z_x_p1 - z_x_m1)) * dxInv[0];

    // This is dh/deta at hi and lo edges
    amrex::Real z_y_p2   = Compute_Z_AtWFace(i,j+2,k,z_nd);
    amrex::Real z_y_p1   = Compute_Z_AtWFace(i,j+1,k,z_nd);
    amrex::Real z_y_m1   = Compute_Z_AtWFace(i,j-1,k,z_nd);
    amrex::Real z_y_m2   = Compute_Z_AtWFace(i,j-2,k,z_nd);
    amrex::Real met_eta  = (amrex::Real(1.0/12.0)*(z_y_m2 - z_y_p2)
                          + amrex::Real(8.0/12.0)*(z_y_p1 - z_y_m1)) * dxInv[1];

    // Use extrapolation instead of interpolation if at the bottom boundary
    amrex::Real u_hi = (k == 0) ? 1.5 *   u_arr(i+1,j  ,k  ) - 0.5 * u_arr(i+1,j  ,k+1)  :
                                  0.5 * ( u_arr(i+1,j  ,k-1) +       u_arr(i+1,j  ,k  ) );
    amrex::Real u_lo = (k == 0) ? 1.5 *   u_arr(i  ,j  ,k  ) - 0.5 * u_arr(i  ,j  ,k+1)  :
                                  0.5 * ( u_arr(i  ,j  ,k-1) +       u_arr(i  ,j  ,k  ) );
    amrex::Real mf_u_hi = mf_u(i+1,j,0);
    amrex::Real mf_u_lo = mf_u(i  ,j,0);

    amrex::Real v_hi = (k == 0) ? 1.5 *   v_arr(i  ,j+1,k  ) - 0.5 * v_arr(i  ,j+1,k+1)  :
                                  0.5 * ( v_arr(i  ,j+1,k-1) +       v_arr(i  ,j+1,k  ) );
    amrex::Real v_lo = (k == 0) ? 1.5 *   v_arr(i  ,j  ,k  ) - 0.5 * v_arr(i  ,j  ,k+1)  :
                                  0.5 * ( v_arr(i  ,j  ,k-1) +       v_arr(i  ,j  ,k  ) );
    amrex::Real mf_v_hi = mf_v(i,j+1,0);
    amrex::Real mf_v_lo = mf_v(i,j  ,0);

    amrex::Real u_met = met_xi  * 0.5 * ( u_hi*mf_u_hi + u_lo*mf_u_lo );
    amrex::Real v_met = met_eta * 0.5 * ( v_hi*mf_v_hi + v_lo*mf_v_lo );

    amrex::Real omega =  w - u_met - v_met;
    return omega;
}

// Define w given u and v arrays and scalar omega (4th order)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
WFromOmega (int& i, int& j, int& k, amrex::Real omega,
            const amrex::Array4<const amrex::Real>& u_arr,
            const amrex::Array4<const amrex::Real>& v_arr,
            const amrex::Array4<const amrex::Real>& mf_u,
            const amrex::Array4<const amrex::Real>& mf_v,
            const amrex::Array4<const amrex::Real>& z_nd,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // This is dh/dxi at hi and lo edges
    amrex::Real z_x_p2  = Compute_Z_AtWFace(i+2,j,k,z_nd);
    amrex::Real z_x_p1  = Compute_Z_AtWFace(i+1,j,k,z_nd);
    amrex::Real z_x_m1  = Compute_Z_AtWFace(i-1,j,k,z_nd);
    amrex::Real z_x_m2  = Compute_Z_AtWFace(i-2,j,k,z_nd);
    amrex::Real met_xi  = (amrex::Real(1.0/12.0)*(z_x_m2 - z_x_p2)
                         + amrex::Real(8.0/12.0)*(z_x_p1 - z_x_m1)) * dxInv[0];

    // This is dh/deta at hi and lo edges
    amrex::Real z_y_p2   = Compute_Z_AtWFace(i,j+2,k,z_nd);
    amrex::Real z_y_p1   = Compute_Z_AtWFace(i,j+1,k,z_nd);
    amrex::Real z_y_m1   = Compute_Z_AtWFace(i,j-1,k,z_nd);
    amrex::Real z_y_m2   = Compute_Z_AtWFace(i,j-2,k,z_nd);
    amrex::Real met_eta  = (amrex::Real(1.0/12.0)*(z_y_m2 - z_y_p2)
                          + amrex::Real(8.0/12.0)*(z_y_p1 - z_y_m1)) * dxInv[1];

    // Use extrapolation instead of interpolation if at the bottom boundary
    amrex::Real u_hi = (k == 0) ? 1.5 *   u_arr(i+1,j  ,k  ) - 0.5 * u_arr(i+1,j  ,k+1)  :
                                  0.5 * ( u_arr(i+1,j  ,k-1) +       u_arr(i+1,j  ,k  ) );
    amrex::Real u_lo = (k == 0) ? 1.5 *   u_arr(i  ,j  ,k  ) - 0.5 * u_arr(i  ,j  ,k+1)  :
                                  0.5 * ( u_arr(i  ,j  ,k-1) +       u_arr(i  ,j  ,k  ) );
    amrex::Real mf_u_hi = mf_u(i+1,j,0);
    amrex::Real mf_u_lo = mf_u(i  ,j,0);

    amrex::Real v_hi = (k == 0) ? 1.5 *   v_arr(i  ,j+1,k  ) - 0.5 * v_arr(i  ,j+1,k+1)  :
                                  0.5 * ( v_arr(i  ,j+1,k-1) +       v_arr(i  ,j+1,k  ) );
    amrex::Real v_lo = (k == 0) ? 1.5 *   v_arr(i  ,j  ,k  ) - 0.5 * v_arr(i  ,j  ,k+1)  :
                                  0.5 * ( v_arr(i  ,j  ,k-1) +       v_arr(i  ,j  ,k  ) );
    amrex::Real mf_v_hi = mf_v(i,j+1,0);
    amrex::Real mf_v_lo = mf_v(i,j  ,0);

    amrex::Real u_met = met_xi  * 0.5 * ( u_hi*mf_u_hi + u_lo*mf_u_lo );
    amrex::Real v_met = met_eta * 0.5 * ( v_hi*mf_v_hi + v_lo*mf_v_lo );

    amrex::Real w = omega + u_met + v_met;
    return w;
}

//
// NOTE: 2nd order kept here for debugging purposes
//

# if 0
// Define omega given u,v and w (2nd order)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
OmegaFromW (int& i, int& j, int& k, amrex::Real w,
            const amrex::Array4<const amrex::Real>& u_arr,
            const amrex::Array4<const amrex::Real>& v_arr,
            const amrex::Array4<const amrex::Real>& mf_u,
            const amrex::Array4<const amrex::Real>& mf_v,
            const amrex::Array4<const amrex::Real>& z_nd,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // This is dh/dxi at hi and lo edges
    amrex::Real z_x_h  = Compute_Z_AtWFace(i+1,j,k,z_nd);
    amrex::Real z_x_l  = Compute_Z_AtWFace(i-1,j,k,z_nd);
    amrex::Real met_xi = 0.5 * (z_x_h - z_x_l) * dxInv[0];

    // This is dh/deta at hi and lo edges
    amrex::Real z_y_h   = Compute_Z_AtWFace(i,j+1,k,z_nd);
    amrex::Real z_y_l   = Compute_Z_AtWFace(i,j-1,k,z_nd);
    amrex::Real met_eta = 0.5 * (z_y_h - z_y_l) * dxInv[1];

    // Use extrapolation instead of interpolation if at the bottom boundary
    amrex::Real u_hi = (k == 0) ? 1.5 *   u_arr(i+1,j  ,k  ) - 0.5 * u_arr(i+1,j  ,k+1)  :
                                  0.5 * ( u_arr(i+1,j  ,k-1) +       u_arr(i+1,j  ,k  ) );
    amrex::Real u_lo = (k == 0) ? 1.5 *   u_arr(i  ,j  ,k  ) - 0.5 * u_arr(i  ,j  ,k+1)  :
                                  0.5 * ( u_arr(i  ,j  ,k-1) +       u_arr(i  ,j  ,k  ) );
    amrex::Real mf_u_hi = mf_u(i+1,j,0);
    amrex::Real mf_u_lo = mf_u(i  ,j,0);

    amrex::Real v_hi = (k == 0) ? 1.5 *   v_arr(i  ,j+1,k  ) - 0.5 * v_arr(i  ,j+1,k+1)  :
                                  0.5 * ( v_arr(i  ,j+1,k-1) +       v_arr(i  ,j+1,k  ) );
    amrex::Real v_lo = (k == 0) ? 1.5 *   v_arr(i  ,j  ,k  ) - 0.5 * v_arr(i  ,j  ,k+1)  :
                                  0.5 * ( v_arr(i  ,j  ,k-1) +       v_arr(i  ,j  ,k  ) );
    amrex::Real mf_v_hi = mf_v(i,j+1,0);
    amrex::Real mf_v_lo = mf_v(i,j  ,0);

    amrex::Real u_met = met_xi  * 0.5 * ( u_hi*mf_u_hi + u_lo*mf_u_lo );
    amrex::Real v_met = met_eta * 0.5 * ( v_hi*mf_v_hi + v_lo*mf_v_lo );

    amrex::Real omega =  w - u_met - v_met;
    return omega;
}

// Define w given u and v arrays and scalar omega (2nd order)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
WFromOmega (int& i, int& j, int& k, amrex::Real omega,
            const amrex::Array4<const amrex::Real>& u_arr,
            const amrex::Array4<const amrex::Real>& v_arr,
            const amrex::Array4<const amrex::Real>& mf_u,
            const amrex::Array4<const amrex::Real>& mf_v,
            const amrex::Array4<const amrex::Real>& z_nd,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // This is dh/dxi at hi and lo edges
    amrex::Real z_x_h  = Compute_Z_AtWFace(i+1,j,k,z_nd);
    amrex::Real z_x_l  = Compute_Z_AtWFace(i-1,j,k,z_nd);
    amrex::Real met_xi = 0.5 * (z_x_h - z_x_l) * dxInv[0];

    // This is dh/deta at hi and lo edges
    amrex::Real z_y_h   = Compute_Z_AtWFace(i,j+1,k,z_nd);
    amrex::Real z_y_l   = Compute_Z_AtWFace(i,j-1,k,z_nd);
    amrex::Real met_eta = 0.5 * (z_y_h - z_y_l) * dxInv[1];

    // Use extrapolation instead of interpolation if at the bottom boundary
    amrex::Real u_hi = (k == 0) ? 1.5 *   u_arr(i+1,j  ,k  ) - 0.5 * u_arr(i+1,j  ,k+1)  :
                                  0.5 * ( u_arr(i+1,j  ,k-1) +       u_arr(i+1,j  ,k  ) );
    amrex::Real u_lo = (k == 0) ? 1.5 *   u_arr(i  ,j  ,k  ) - 0.5 * u_arr(i  ,j  ,k+1)  :
                                  0.5 * ( u_arr(i  ,j  ,k-1) +       u_arr(i  ,j  ,k  ) );
    amrex::Real mf_u_hi = mf_u(i+1,j,0);
    amrex::Real mf_u_lo = mf_u(i  ,j,0);

    amrex::Real v_hi = (k == 0) ? 1.5 *   v_arr(i  ,j+1,k  ) - 0.5 * v_arr(i  ,j+1,k+1)  :
                                  0.5 * ( v_arr(i  ,j+1,k-1) +       v_arr(i  ,j+1,k  ) );
    amrex::Real v_lo = (k == 0) ? 1.5 *   v_arr(i  ,j  ,k  ) - 0.5 * v_arr(i  ,j  ,k+1)  :
                                  0.5 * ( v_arr(i  ,j  ,k-1) +       v_arr(i  ,j  ,k  ) );
    amrex::Real mf_v_hi = mf_v(i,j+1,0);
    amrex::Real mf_v_lo = mf_v(i,j  ,0);

    amrex::Real u_met = met_xi  * 0.5 * ( u_hi*mf_u_hi + u_lo*mf_u_lo );
    amrex::Real v_met = met_eta * 0.5 * ( v_hi*mf_v_hi + v_lo*mf_v_lo );

    amrex::Real w = omega + u_met + v_met;
    return w;
}
#endif


//*****************************************************************************************
// Rotate scalar flux vector and stress tensor for MOST
//*****************************************************************************************
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
rotate_scalar_flux (const int& i,
                    const int& j,
                    const int& klo,
                    const amrex::Real& flux,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                    const amrex::Array4<const amrex::Real>& zphys_arr,
                    const amrex::Array4<amrex::Real>& phi1_arr,
                    const amrex::Array4<amrex::Real>& phi2_arr,
                    const amrex::Array4<amrex::Real>& phi3_arr)
{
    amrex::Real h_xi  = Compute_h_xi_AtCellCenter(i, j, klo, dxInv, zphys_arr);
    amrex::Real h_eta = Compute_h_eta_AtCellCenter(i, j, klo, dxInv, zphys_arr);
    amrex::Real InvNorm  = 1.0 / std::sqrt(1.0 + h_xi*h_xi + h_eta*h_eta);
    phi1_arr(i,j,klo) = -h_xi  * flux * InvNorm;
    phi2_arr(i,j,klo) = -h_eta * flux * InvNorm;
    phi3_arr(i,j,klo) =          flux * InvNorm;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
rotate_stress_tensor (const int& i,
                      const int& j,
                      const int& klo,
                      const amrex::Real& flux,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                      const amrex::Array4<const amrex::Real>& zphys_arr,
                      const amrex::Array4<const amrex::Real>& u_arr,
                      const amrex::Array4<const amrex::Real>& v_arr,
                      const amrex::Array4<const amrex::Real>& w_arr,
                      const amrex::Array4<amrex::Real>& tau11_arr,
                      const amrex::Array4<amrex::Real>& tau22_arr,
                      const amrex::Array4<amrex::Real>& tau33_arr,
                      const amrex::Array4<amrex::Real>& tau12_arr,
                      const amrex::Array4<amrex::Real>& tau21_arr,
                      const amrex::Array4<amrex::Real>& tau13_arr,
                      const amrex::Array4<amrex::Real>& tau31_arr,
                      const amrex::Array4<amrex::Real>& tau23_arr,
                      const amrex::Array4<amrex::Real>& tau32_arr)
{
    // Unit-normal vector
    amrex::Array1D<amrex::Real,0,2> n_hat;

    // Unit-tangent vectors
    amrex::Array1D<amrex::Real,0,2> t_hat_1;
    amrex::Array1D<amrex::Real,0,2> t_hat_2;
    amrex::Array1D<amrex::Real,0,2> u_t_hat;

    // Final basis vector for right-hand coord sys
    amrex::Array1D<amrex::Real,0,2> a_hat;

    // Rotation matrix
    amrex::Array2D<amrex::Real,0,2,0,2> R_mat;

    // Metric data
    amrex::Real h_xi   = Compute_h_xi_AtIface(i, j, klo, dxInv, zphys_arr);
    amrex::Real h_eta  = Compute_h_eta_AtIface(i, j, klo, dxInv, zphys_arr);
    amrex::Real h_zeta = Compute_h_zeta_AtIface(i, j, klo, dxInv, zphys_arr);

    // Populate the normal vector
    amrex::Real Inormn = 1./std::sqrt(1.0 + h_xi*h_xi + h_eta*h_eta);
    n_hat(0) = Inormn*h_xi; n_hat(1) = Inormn*h_eta; n_hat(2) = -Inormn;

    // Populate the tangent vectors
    amrex::Real Inorm1 = 1./std::sqrt(1.0 + h_xi*h_xi);
    amrex::Real Inorm2 = 1./std::sqrt(1.0 + h_eta*h_eta);
    t_hat_1(0) = -Inorm1;      t_hat_2(1) = -Inorm2;
    t_hat_1(2) = -Inorm1*h_xi; t_hat_2(2) = -Inorm2*h_eta;

    // Populate the u_t vector
    amrex::Real Norm_u_t = 0.0;
    amrex::Real mag1 = (u_arr(i,j,klo) + h_xi *w_arr(i,j,klo))*Inorm1;
    amrex::Real mag2 = (v_arr(i,j,klo) + h_eta*w_arr(i,j,klo))*Inorm2;
    for (int icol(0); icol<3; ++icol) {
        u_t_hat(icol) = mag1*t_hat_1(icol) + mag2*t_hat_2(icol);
        Norm_u_t  += u_t_hat(icol)*u_t_hat(icol);
    }
    for (int icol(0); icol<3; ++icol) {
        u_t_hat(icol) /= std::sqrt(Norm_u_t);
    }

    // Populate the a_hat vector
    a_hat(0) =   n_hat(1)*u_t_hat(2) - n_hat(2)*u_t_hat(1);
    a_hat(1) = -(n_hat(0)*u_t_hat(2) - n_hat(2)*u_t_hat(0));
    a_hat(2) =   n_hat(0)*u_t_hat(1) - n_hat(1)*u_t_hat(0);

    // Copy column vectors into R_mat
    int jrow;
    jrow = 0;
    for (int icol(0); icol<3; ++icol) {
        R_mat(icol,jrow) = u_t_hat(icol);
    }
    jrow = 1;
    for (int icol(0); icol<3; ++icol) {
        R_mat(icol,jrow) = a_hat(icol);
    }
    jrow = 2;
    for (int icol(0); icol<3; ++icol) {
        R_mat(icol,jrow) = n_hat(icol);
    }

    // Body-fixed tau
    amrex::Real T11    = 2.0*R_mat(0,0)*R_mat(2,0)*flux;
    amrex::Real T22    = 2.0*R_mat(0,1)*R_mat(2,1)*flux;
    amrex::Real T33    = 2.0*R_mat(0,2)*R_mat(2,2)*flux;
    amrex::Real T12    = (R_mat(0,0)*R_mat(2,1) + R_mat(2,0)*R_mat(0,1))*flux;
    amrex::Real T13    = (R_mat(0,0)*R_mat(2,2) + R_mat(0,2)*R_mat(2,0))*flux;
    amrex::Real T23    = (R_mat(0,1)*R_mat(2,2) + R_mat(0,2)*R_mat(2,1))*flux;

    // Rotated stress fluxes
    tau11_arr(i,j,klo) = h_zeta*T11;
    tau22_arr(i,j,klo) = h_zeta*T22;
    tau33_arr(i,j,klo) = -h_xi*T13 - h_eta*T23 + T33;

    tau12_arr(i,j,klo) = h_zeta*T12;
    tau21_arr(i,j,klo) = tau21_arr(i,j,klo);

    tau13_arr(i,j,klo) = -h_xi*T11 - h_eta*T12 + T13;
    tau31_arr(i,j,klo) = h_zeta*T13;

    tau23_arr(i,j,klo) = -h_xi*T12 - h_eta*T22 + T23;
    tau32_arr(i,j,klo) = h_zeta*T23;
}
#endif
