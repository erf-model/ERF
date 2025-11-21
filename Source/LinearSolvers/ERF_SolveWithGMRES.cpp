#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainPoisson.H"

#include <AMReX_GMRES.H>

using namespace amrex;

/**
 * Solve the Poisson equation using FFT-preconditioned GMRES
 */
void ERF::solve_with_gmres (int lev, const Box& subdomain, MultiFab& rhs, MultiFab& phi,
                            Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                            MultiFab& ax_sub, MultiFab& ay_sub, MultiFab& az_sub,
                            MultiFab& dJ_sub, MultiFab& znd_sub)
{
#ifdef ERF_USE_FFT
    BL_PROFILE("ERF::solve_with_gmres()");

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    auto const dom_lo = lbound(Geom(lev).Domain());
    auto const dom_hi = ubound(Geom(lev).Domain());

    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);

    auto dx    = Geom(lev).CellSizeArray();

    Geometry my_geom;

    Array<int,AMREX_SPACEDIM> is_per; is_per[0] = 0; is_per[1] = 0; is_per[2] = 0;
    if (Geom(lev).isPeriodic(0) && sub_lo.x == dom_lo.x && sub_hi.x == dom_hi.x) { is_per[0] = 1;}
    if (Geom(lev).isPeriodic(1) && sub_lo.y == dom_lo.y && sub_hi.y == dom_hi.y) { is_per[1] = 1;}

    int coord_sys = 0;

    // If subdomain == domain then we pass Geom(lev) to the FFT solver
    if (subdomain == Geom(lev).Domain()) {
        my_geom.define(Geom(lev).Domain(), Geom(lev).ProbDomain(), coord_sys, is_per);
    } else {
        // else we create a new geometry based only on the subdomain
        // The information in my_geom used by the FFT routines is:
        //   1) my_geom.Domain()
        //   2) my_geom.CellSize()
        //   3) my_geom.isAllPeriodic() / my_geom.periodicity()
        RealBox rb( sub_lo.x   *dx[0],  sub_lo.y   *dx[1],  sub_lo.z   *dx[2],
                   (sub_hi.x+1)*dx[0], (sub_hi.y+1)*dx[1], (sub_hi.z+1)*dx[2]);
        my_geom.define(subdomain, rb, coord_sys, is_per);
    }

    amrex::GMRES<MultiFab, TerrainPoisson> gmsolver;

    TerrainPoisson tp(my_geom, rhs.boxArray(), rhs.DistributionMap(), domain_bc_type,
                      stretched_dz_d[lev], ax_sub, ay_sub, az_sub, dJ_sub, &znd_sub,
                      solverChoice.use_real_bcs);

    gmsolver.define(tp);

    gmsolver.setVerbose(mg_verbose);

    gmsolver.setRestartLength(50);

    tp.usePrecond(true);

    gmsolver.solve(phi, rhs, reltol, abstol);

    tp.getFluxes(phi, fluxes);

    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        const Array4<Real      >& fx_ar = fluxes[0].array(mfi);
        const Array4<Real      >& fy_ar = fluxes[1].array(mfi);
        const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
        ParallelFor(xbx,ybx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_ar(i,j,k) *= mf_ux(i,j,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_ar(i,j,k) *= mf_vy(i,j,0);
        });
    } // mfi
#else
    amrex::ignore_unused(lev, rhs, phi, fluxes, ax_sub, ay_sub, az_sub, dJ_sub, znd_sub);
#endif

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
}
