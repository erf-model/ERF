#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

using namespace amrex;

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 */

using BCType = LinOpBCType;

/**
 * Solve the Poisson equation using EB_enabled MLMG
 * Note that the level may or may not be level 0.
 */
void ERF::solve_with_EB_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi, Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes)
{
    BL_PROFILE("ERF::solve_with_EB_mlmg()");

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
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(rhs[0].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(rhs[0].DistributionMap());
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

    MLEBABecLap mleb (geom_tmp, ba_tmp, dm_tmp, info, {&EBFactory(lev)});

    mleb.setMaxOrder(2);
    mleb.setDomainBC(bclo, bchi);
    mleb.setLevelBC(0, nullptr);

    //
    // This sets A = 0, B = 1 so that
    // the operator A alpha - b del dot beta grad to b
    // becomes  - del dot beta grad
    //
    mleb.setScalars(0.0, 1.0);

    Array<MultiFab,AMREX_SPACEDIM> bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        bcoef[idim].define(convert(ba_tmp[0],IntVect::TheDimensionVector(idim)),
                            dm_tmp[0], 1, 0, MFInfo(), EBFactory(lev));
        bcoef[idim].setVal(-1.0);
    }
    mleb.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));

    MLMG mlmg(mleb);

    int max_iter = 100;
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), reltol, abstol);

    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    ImposeBCsOnPhi(lev,phi[0]);

    //
    // This arises because we solve MINUS del dot beta grad phi = div (rho u)
    //
    fluxes[0][0].mult(-1.);
    fluxes[0][1].mult(-1.);
    fluxes[0][2].mult(-1.);
}
