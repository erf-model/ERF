/**
 * \file ERF_SolveWithGMRES.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainPoisson.H"

#include <AMReX_GMRES.H>

using namespace amrex;

/**
 * Solve the Poisson equation using FFT-preconditioned GMRES
 *
 * @param lev Level index for the solve
 * @param subdomain Box over which the solve is performed
 * @param rhs Right-hand side field for the Poisson solve
 * @param p Solution field to fill
 * @param fluxes Face-centered gradient fluxes to fill
 * @param ax_sub Terrain metric coefficient on x-faces
 * @param ay_sub Terrain metric coefficient on y-faces
 * @param az_sub Terrain metric coefficient on z-faces
 * @param znd_sub Node-centered physical height field
 *
 * The unnamed metric argument in the class declaration is the cell-centered
 * Jacobian determinant.
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

    // The FFT preconditioner is built on the region, so it needs the dz of the region's own
    // layers only (make_sub_dz), as solve_with_fft passes on a stretched mesh.  A region that
    // stops below the top or above the ground (a level refined over part of the depth) would
    // otherwise hand it the whole column.
    Gpu::DeviceVector<Real> dz_region = make_sub_dz(lev, subdomain);

    amrex::GMRES<MultiFab, TerrainPoisson> gmsolver;

    TerrainPoisson tp(my_geom, Geom(lev), rhs.boxArray(), rhs.DistributionMap(), domain_bc_type,
                      dz_region, ax_sub, ay_sub, az_sub, dJ_sub, &znd_sub,
                      solverChoice.use_real_bcs);

    gmsolver.define(tp);

    gmsolver.setVerbose(mg_verbose);

    gmsolver.setRestartLength(50);

    tp.usePrecond(true);

    gmsolver.solve(phi, rhs, reltol, abstol);

    // The caller (project_momenta) multiplies the lateral fluxes by the map factor
    tp.getFluxes(phi, fluxes);
#else
    amrex::ignore_unused(lev, rhs, phi, fluxes, ax_sub, ay_sub, az_sub, dJ_sub, znd_sub);
#endif

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
}

/**
 * The stretched dz of a solve region's own layers, k = subdomain lo to hi in z, on the device.
 *
 * @param lev       Level index
 * @param subdomain Region solved on
 */
Gpu::DeviceVector<Real> ERF::make_sub_dz (int lev, const Box& subdomain)
{
    const int klo = subdomain.smallEnd(2);
    Vector<Real> dz_h(subdomain.length(2));
    for (int k = klo; k <= subdomain.bigEnd(2); k++) {
        dz_h[k-klo] = stretched_dz_h[lev][k];
    }
    Gpu::DeviceVector<Real> dz_d(dz_h.size());
    Gpu::copy(Gpu::hostToDevice, dz_h.begin(), dz_h.end(), dz_d.begin());
    return dz_d;
}
