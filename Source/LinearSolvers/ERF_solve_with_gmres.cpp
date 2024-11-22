#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainPoisson.H"

#include <AMReX_GMRES.H>

using namespace amrex;

/**
 * Solve the Poisson equation using GMRES
 */
void ERF::solve_with_gmres (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi, Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes)
{
    BL_PROFILE("ERF::solve_with_gmres()");

    amrex::GMRES<MultiFab, TerrainPoisson> gmsolver;

    TerrainPoisson tp(geom[lev], rhs[lev].boxArray(), rhs[lev].DistributionMap(), z_phys_nd[lev].get());

    gmsolver.define(tp);

    gmsolver.setVerbose(mg_verbose);

//   tp.usePrecond(false);

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    gmsolver.solve(phi[0], rhs[0], reltol, abstol);

    amrex::Print() << "PHI " << phi[0][0] << std::endl;

    tp.getFluxes(phi[0], fluxes[0]);
}
