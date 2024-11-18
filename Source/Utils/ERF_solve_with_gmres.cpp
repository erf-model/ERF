#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
//#include <AMReX_MLTerrainPoisson.H>
#include <AMReX_GMRES.H>
#include <AMReX_GMRES_MLMG.H>

using namespace amrex;

/**
 * Solve the Poisson equation using GMRES
 */
void ERF::solve_with_gmres (int lev, Vector<MultiFab>& /*rhs*/, Vector<MultiFab>& /*phi*/, Vector<Array<MultiFab,AMREX_SPACEDIM>>& /*fluxes*/)
//void ERF::solve_with_gmres (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi, Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes)
{
    BL_PROFILE("ERF::solve_with_gmres()");

#if 0
    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(rhs[lev].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(rhs[lev].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);

    // amrex::Print() << "BCLO " << bclo[0] << " " << bclo[1] << " " << bclo[2] << std::endl;
    // amrex::Print() << "BCHI " << bchi[0] << " " << bchi[1] << " " << bchi[2] << std::endl;

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    MLTerrainPoisson terrpoisson(geom_tmp, ba_tmp, dm_tmp, info);
    terrpoisson.setDomainBC(bclo, bchi);
    terrpoisson.setMaxOrder(2);

    terrpoisson.setZPhysNd(lev, *z_phys_nd[lev]);

    if (lev > 0) {
        terrpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    terrpoisson.setLevelBC(lev, &phi[lev]);

    MLMG mlmg(terrpoisson);
    GMRESMLMG gmsolver(mlmg);
    gmsolver.usePrecond(false);
    gmsolver.setVerbose(mg_verbose);
    gmsolver.solve(phi[0], rhs[0], reltol, abstol);

    Vector<MultiFab*> phi_vec; phi_vec.resize(1);
    phi_vec[0] = &phi[0];
    terrpoisson.getFluxes(GetVecOfArrOfPtrs(fluxes), phi_vec);
#endif
}
