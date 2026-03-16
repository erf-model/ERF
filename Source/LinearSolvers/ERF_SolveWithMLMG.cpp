#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

/**
 * Solve the Poisson equation using MLMG
 * Note that the level may or may not be level 0.
 *
 * Important: we solve on the whole level even if there are disjoint regions
 *
 */
void
solve_with_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi,
                 Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                 const Geometry& geom, const Vector<IntVect>& ref_ratio,
                 Array<std::string,2*AMREX_SPACEDIM> domain_bc_type,
                 int mg_verbose, Real reltol, Real abstol)
{
    BL_PROFILE("ERF::solve_with_mlmg()");

    auto const dom_lo = lbound(geom.Domain());
    auto const dom_hi = ubound(geom.Domain());

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(rhs[0].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(rhs[0].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom);

    auto bclo = get_lo_projection_bc(geom,domain_bc_type);
    auto bchi = get_hi_projection_bc(geom,domain_bc_type);

    // amrex::Print() << "BCLO " << bclo[0] << " " << bclo[1] << " " << bclo[2] << std::endl;
    // amrex::Print() << "BCHI " << bchi[0] << " " << bchi[1] << " " << bchi[2] << std::endl;

    // ****************************************************************************
    // Multigrid solve
    // ****************************************************************************

    MLPoisson mlpoisson(geom_tmp, ba_tmp, dm_tmp, info);
    mlpoisson.setDomainBC(bclo, bchi);
    if (lev > 0) {
        mlpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    mlpoisson.setLevelBC(0, nullptr);

    // Use low order for outflow at physical boundaries
    mlpoisson.setMaxOrder(2);

    MLMG mlmg(mlpoisson);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve(GetVecOfPtrs(phi),
               GetVecOfConstPtrs(rhs),
               reltol, abstol);
    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    // ImposeBCsOnPhi(lev, phi[0], geom.Domain()));
}
