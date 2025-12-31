#include "ERF.H"
#include "ERF_EB.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

using namespace amrex;

void
FillZeroAreaFaceFluxes (MultiFab& phi,
                        Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                        const Geometry& geom,
                        EBFArrayBoxFactory const& ebfact,
                        eb_aux_ const& ebfact_u,
                        eb_aux_ const& ebfact_v,
                        eb_aux_ const& ebfact_w);

/**
 * Solve the Poisson equation using EB_enabled MLMG
 * Note that the level may or may not be level 0.
 *
 * Important: we solve on the whole level even if there are disjoint regions
 *
 */
void
solve_with_EB_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi,
                    Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                    EBFArrayBoxFactory const& ebfact,
                    eb_aux_ const& ebfact_u,
                    eb_aux_ const& ebfact_v,
                    eb_aux_ const& ebfact_w,
                    const Geometry& geom, const Vector<amrex::IntVect>& ref_ratio,
                    Array<std::string,2*AMREX_SPACEDIM> domain_bc_type,
                    int mg_verbose, Real reltol, Real abstol)

{
    BL_PROFILE("ERF::solve_with_EB_mlmg()");

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

    amrex::Print() << "BCLO " << bclo[0] << " " << bclo[1] << " " << bclo[2] << std::endl;
    amrex::Print() << "BCHI " << bchi[0] << " " << bchi[1] << " " << bchi[2] << std::endl;

    // ****************************************************************************
    // Multigrid solve
    // ****************************************************************************

    MLEBABecLap mleb (geom_tmp, ba_tmp, dm_tmp, info, {&ebfact});

    mleb.setMaxOrder(2);
    mleb.setDomainBC(bclo, bchi);
    if (lev > 0) {
        mleb.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
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
                            dm_tmp[0], 1, 0, MFInfo(), ebfact);
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

    // Add the flux values (gradient phi) to the face-centered cell with zero area fraction
    // (mlmg.getFluxes does not fill gradient phi at faces with zero area fraction)
    FillZeroAreaFaceFluxes(phi[0], fluxes[0], geom, ebfact, ebfact_u, ebfact_v, ebfact_w);

    // ImposeBCsOnPhi(lev,phi[0], geom[lev].Domain());

    //
    // This arises because we solve MINUS del dot beta grad phi = div (rho u)
    //
    fluxes[0][0].mult(-1.);
    fluxes[0][1].mult(-1.);
    fluxes[0][2].mult(-1.);
}
