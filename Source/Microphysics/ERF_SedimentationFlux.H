/*! @file ERF_SedimentationFlux.H
 *  \brief The sedimentation face fluxes a precipitating microphysics scheme records
 */

#ifndef ERF_SEDIMENTATIONFLUX_H
#define ERF_SEDIMENTATIONFLUX_H

#include <memory>

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

/**
 * The hydrometeors a microphysics scheme moved across each z-face by sedimentation during
 * its most recent Advance on a level.
 *
 * down(i,j,k,n) is the amount of the conserved state component cons_comp[n] (rho*q, so
 * kg m^-3 for a mixing ratio and m^-3 for a number concentration) that crossed z-face
 * (i,j,k) downward during that Advance, per unit horizontal area: kg m^-2 or m^-2. Face k
 * is the bottom face of cell k, so the face at the domain bottom holds what reached the
 * ground. down is defined on the level's grids, z-nodal, and is zeroed at the start of
 * every Advance.
 *
 * A scheme records the fluxes its cell updates actually used, so that (up to the clipping
 * of negative values) the change of rho*q*dz in a cell is down on its top face minus down on
 * its bottom face. Schemes whose sedimentation has no face fluxes of its own (the
 * semi-Lagrangian WSM6/WDM6 solvers) reconstruct them with sedimentation_add_column_sums below.
 */
/**
 * The inflow a fine patch receives through its top face, handed to a scheme for one Advance
 * (NullMoist::Set_Sed_Inflow). It is z-nodal on the level's grids, with the components of
 * the level's SedimentationFlux plus a flag, comp inflow_flag_comp(): where the flag is 1,
 * face (i,j,k) is the top of a fine column that stops short of the domain top, and
 * component n holds the amount of cons_comp[n] entering through it during the Advance, per
 * unit horizontal area (as SedimentationFlux::down). Elsewhere everything is zero.
 */
struct SedimentationFlux
{
    std::shared_ptr<amrex::MultiFab> down;
    amrex::Vector<int> cons_comp;
    amrex::Real dt = 0.0;

    void define (const amrex::BoxArray& cc_ba, const amrex::DistributionMapping& dm,
                 const amrex::Vector<int>& a_cons_comp)
    {
        cons_comp = a_cons_comp;
        down = std::make_shared<amrex::MultiFab>(amrex::convert(cc_ba, amrex::IntVect(0,0,1)), dm,
                                                 static_cast<int>(cons_comp.size()), 0);
        down->setVal(0.0);
        dt = 0.0;
    }

    //! The component of a sedimentation inflow (see above) that flags the faces it applies at
    [[nodiscard]] int inflow_flag_comp () const { return static_cast<int>(cons_comp.size()); }

    //! Zero the fluxes at the start of an Advance of length a_dt
    void reset (amrex::Real a_dt)
    {
        if (down) { down->setVal(0.0); }
        dt = a_dt;
    }
};

/**
 * Add sign * (sum over m = k..khi of w(m)) to sed(i,j,k,n) for every k in [klo, khi].
 *
 * Called with sign = +1 and w(m) = rho*q*dz before a sedimentation step and with sign = -1
 * and the same w after it, this adds to each face k the mass that crossed it downward: what
 * the cells above it lost. The face above khi gets nothing (no inflow).
 */
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
sedimentation_add_column_sums (const amrex::Array4<amrex::Real>& sed, int n,
                               int i, int j, int klo, int khi, amrex::Real sign, F&& w) noexcept
{
    amrex::Real sum = 0.0;
    for (int k = khi; k >= klo; --k) {
        sum += w(k);
        sed(i,j,k,n) += sign * sum;
    }
}

#endif
