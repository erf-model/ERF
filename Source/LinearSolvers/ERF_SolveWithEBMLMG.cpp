/**
 * \file ERF_SolveWithEBMLMG.cpp
 */
#include "ERF.H"
#include "ERF_EB.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

using namespace amrex;

/**
 * Compute phi gradients where the area of the face is zero.
 *
 * @tparam T EB factory or auxiliary EB data type for face-centered grids
 * @param[in] phi Cell-centered solution used to compute gradients
 * @param[out] fluxes Face-centered gradient fluxes to fill
 * @param[in] geom Geometry used for inverse cell spacing
 * @param[in] ebfact Cell-centered embedded-boundary factory
 * @param[in] ebfact_u Embedded-boundary data on x-faces
 * @param[in] ebfact_v Embedded-boundary data on y-faces
 * @param[in] ebfact_w Embedded-boundary data on z-faces
 */
template <typename T>
void
FillZeroAreaFaceFluxes (MultiFab& phi,
                        Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                        const Geometry& geom,
                        EBFArrayBoxFactory const& ebfact,
                        T const& ebfact_u,
                        T const& ebfact_v,
                        T const& ebfact_w);

/**
 * Solve the Poisson equation using EB_enabled MLMG
 * Note that the level may or may not be level zero
 *
 * Important: we solve on the whole level even if there are disjoint regions
 *
 * @tparam T EB factory or auxiliary EB data type for face-centered grids
 * @param lev Level index for the solve
 * @param rhs Right-hand side MultiFab vector
 * @param phi Solution MultiFab vector to fill
 * @param fluxes Face-centered gradient fluxes to fill
 * @param ebfact Cell-centered embedded-boundary factory
 * @param ebfact_u Embedded-boundary data on x-faces
 * @param ebfact_v Embedded-boundary data on y-faces
 * @param ebfact_w Embedded-boundary data on z-faces
 * @param geom Geometry for the solve level
 * @param ref_ratio Coarse-fine refinement ratios
 * @param domain_bc_type Domain boundary-condition names
 * @param mg_verbose MLMG verbosity level
 * @param reltol Relative solver tolerance
 * @param abstol Absolute solver tolerance
 */
template <typename T>
void
solve_with_EB_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi,
                    Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                    EBFArrayBoxFactory const& ebfact,
                    T const& ebfact_u,
                    T const& ebfact_v,
                    T const& ebfact_w,
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
    mleb.setScalars(zero, one);

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
    //
    // That routine extrapolates the gradient onto a zero-area face from three cells on the
    // uncovered side, and at a face on the edge of a box the stencil reaches three cells
    // past the box -- well beyond the single ghost cell phi carries, and beyond whatever
    // the solve happened to leave in it (issue #3699). Do the extrapolation out of a
    // three-ghost-cell scratch copy instead of widening phi itself, which would cost MLMG
    // its alias fast path for the solution. Ghost cells outside the domain in a
    // non-periodic direction are still unfilled; FillZeroAreaFaceFluxes detects that and
    // drops to a lower order one-sided form there.
    MultiFab phi_ext(phi[0].boxArray(), phi[0].DistributionMap(), 1, 3);
    MultiFab::Copy(phi_ext, phi[0], 0, 0, 1, 0);
    phi_ext.FillBoundary(geom.periodicity());

    FillZeroAreaFaceFluxes(phi_ext, fluxes[0], geom, ebfact, ebfact_u, ebfact_v, ebfact_w);

    // ImposeBCsOnPhi(lev,phi[0], geom[lev].Domain());

    //
    // This arises because we solve MINUS del dot beta grad phi = div (rho u)
    //
    fluxes[0][0].mult(-one);
    fluxes[0][1].mult(-one);
    fluxes[0][2].mult(-one);
}

// Explicit template instantiations for the types we use
// When USE_FC_FACTORY=1, all are EBFArrayBoxFactory
template void solve_with_EB_mlmg(int, Vector<MultiFab>&, Vector<MultiFab>&,
                                 Vector<Array<MultiFab,AMREX_SPACEDIM>>&,
                                 EBFArrayBoxFactory const&,
                                 EBFArrayBoxFactory const&,
                                 EBFArrayBoxFactory const&,
                                 EBFArrayBoxFactory const&,
                                 const Geometry&, const Vector<amrex::IntVect>&,
                                 Array<std::string,2*AMREX_SPACEDIM>, int, Real, Real);
// When USE_FC_FACTORY=0, u/v/w are eb_aux_
template void solve_with_EB_mlmg(int, Vector<MultiFab>&, Vector<MultiFab>&,
                                 Vector<Array<MultiFab,AMREX_SPACEDIM>>&,
                                 EBFArrayBoxFactory const&,
                                 eb_aux_ const&,
                                 eb_aux_ const&,
                                 eb_aux_ const&,
                                 const Geometry&, const Vector<amrex::IntVect>&,
                                 Array<std::string,2*AMREX_SPACEDIM>, int, Real, Real);
