/**
 * \file ERF_init_from_hse.cpp
 */

#include <ERF.H>
#include <TileNoZ.H>
#include <EOS.H>

#include "AMReX_Print.H"

using namespace amrex;

/**
 * Initialize the background flow to have the calculated HSE density and
 * rho*theta calculated from the HSE pressure. In general, the hydrostatically
 * balanced density and pressure (r_hse and p_hse from base_state) used here may
 * be calculated through a solver path such as:
 *
 *   ERF::initHSE(lev)
 *   - call prob->erf_init_dens_hse(...)
 *     - call Problem::init_isentropic_hse(...), to simultaneously calculate
 *       r_hse and p_hse with Newton iteration -- assuming constant theta
 *     - save r_hse
 *   - call ERF::enforce_hse(...), calculates p_hse from saved r_hse (redundant,
 *     but needed because p_hse is not necessarily calculated by the Problem
 *     implementation) and pi_hse -- note: this pressure does not exactly match
 *     the p_hse from before because what is calculated by init_isentropic_hse
 *     comes from the EOS whereas what is calculated here comes from the hydro-
 *     static equation
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_hse (int lev)
{
    auto& lev_new = vars_new[lev];

    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi)
    {
        const Box &gbx = mfi.growntilebox(1);
        const Array4<Real      >& cons_arr  = lev_new[Vars::cons].array(mfi);
        const Array4<Real const>& r_hse_arr = r_hse.const_array(mfi);
        const Array4<Real const>& p_hse_arr = p_hse.const_array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            cons_arr(i,j,k,Rho_comp)      = r_hse_arr(i,j,k);
            cons_arr(i,j,k,RhoTheta_comp) = getRhoThetagivenP(p_hse_arr(i,j,k));
        });
    } //mfi
}
