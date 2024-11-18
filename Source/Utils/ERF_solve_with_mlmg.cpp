#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 * if we want to enforce incompressibility of the initial conditions
 */

using BCType = LinOpBCType;

Array<LinOpBCType,AMREX_SPACEDIM>
ERF::get_projection_bc (Orientation::Side side) const noexcept
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = domain_bc_type[Orientation(dir,side)];
            if (bc_type == "Outflow") {
                r[dir] = LinOpBCType::Dirichlet;
            } else
            {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    return r;
}

/**
 * Solve the Poisson equation using MLMG
 * Note that the level may or may not be level 0.
 */
void ERF::solve_with_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi, Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes)
{
    BL_PROFILE("ERF::solve_with_mlmg()");

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

    // ****************************************************************************
    // Multigrid solve
    // ****************************************************************************

    MLPoisson mlpoisson(geom_tmp, ba_tmp, dm_tmp, info);
    mlpoisson.setDomainBC(bclo, bchi);
    if (lev > 0) {
        mlpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    mlpoisson.setLevelBC(0, nullptr);

    MLMG mlmg(mlpoisson);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve(GetVecOfPtrs(phi),
               GetVecOfConstPtrs(rhs),
               reltol, abstol);
    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));
}
