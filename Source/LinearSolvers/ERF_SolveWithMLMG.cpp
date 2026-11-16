/**
 * \file ERF_SolveWithMLMG.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#include <limits>
#include <memory>

using namespace amrex;

namespace {

/**
 * MLPoisson whose MG levels coarsen each direction on its own (as MLTerrainPoisson does): a
 * direction whose cells are at least 1.5 times longer than the shortest coarsenable direction
 * is not coarsened on that level, so the MG levels approach isotropic cells.
 */
class MLPoissonAniso : public MLPoisson {
public:
    MLPoissonAniso (const Vector<Geometry>& a_geom, const Vector<BoxArray>& a_grids,
                    const Vector<DistributionMapping>& a_dmap, const LPInfo& a_info)
    {
        mg_independent_coarsening = true;
        define(a_geom, a_grids, a_dmap, a_info);
    }
};

/**
 * Number of MG coarsenings for MLPoissonAniso. AMReX keeps coarsening the other directions
 * once one can no longer be coarsened (e.g. a domain few cells wide), which makes the coarse
 * levels anisotropic the other way. Stop instead at the level where a direction that cannot be
 * coarsened (by the domain or the grids) is 1.5 times shorter than the shortest one that can,
 * and leave the rest to the bottom solver.
 */
int
aniso_max_coarsening_level (const Geometry& geom, const BoxArray& grids, const LPInfo& info)
{
    const int hidden = info.hidden_direction;
    const Real* dx = geom.CellSize();
    const Box& domain = geom.Domain();
    IntVect rr(1);
    int nlev = 0;
    for (int lev = 0; lev < info.max_coarsening_level; ++lev) {
        Real h_all = std::numeric_limits<Real>::max();
        Real h_crse = std::numeric_limits<Real>::max();
        IntVect can_coarsen(0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim == hidden) { continue; }
            const Real h = dx[idim]*Real(rr[idim]);
            IntVect rr2 = rr; rr2[idim] *= 2;
            IntVect min_width(0); min_width[idim] = 2;
            can_coarsen[idim] = (domain.coarsenable(rr2, min_width) &&
                                 grids.coarsenable(rr2, min_width)) ? 1 : 0;
            h_all = std::min(h_all, h);
            if (can_coarsen[idim]) { h_crse = std::min(h_crse, h); }
        }
        if (h_crse == std::numeric_limits<Real>::max() || h_all*Real(1.5) <= h_crse) { break; }
        IntVect rr_level(1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (can_coarsen[idim] && dx[idim]*Real(rr[idim]) < Real(1.5)*h_crse) { rr_level[idim] = 2; }
        }
        rr *= rr_level;
        ++nlev;
    }
    return nlev;
}

} // namespace

/**
 * Solve the Poisson equation using MLMG
 * Note that the level may or may not be level zero
 *
 * Important: we solve on the whole level even if there are disjoint regions
 *
 * @param lev Level index for the solve
 * @param rhs Right-hand side MultiFab vector
 * @param phi Solution MultiFab vector to fill
 * @param fluxes Face-centered gradient fluxes to fill
 * @param geom Geometry for the solve level
 * @param ref_ratio Coarse-fine refinement ratios
 * @param domain_bc_type Domain boundary-condition names
 * @param mg_verbose MLMG verbosity level
 * @param reltol Relative solver tolerance
 * @param abstol Absolute solver tolerance
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

    // MLPoisson coarsens every direction by 2 on each MG level, so cells longer in some
    // direction stay that anisotropic all the way down. Red-black Gauss-Seidel then smooths only
    // along the short directions, and with piecewise-constant prolongation the coarse-grid
    // correction overshoots the error left along the long ones: at a cell aspect ratio of 4 the
    // V-cycle diverges (#5064). For anisotropic cells, coarsen the directions on their own.
    bool aniso = false;
    {
        const Real* dx = geom.CellSize();
        Real hmin = std::numeric_limits<Real>::max();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim != info.hidden_direction) { hmin = std::min(hmin, dx[idim]); }
        }
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim != info.hidden_direction && dx[idim] >= Real(1.5)*hmin) { aniso = true; }
        }
    }
    if (aniso) {
        info.setMaxCoarseningLevel(aniso_max_coarsening_level(geom, rhs[0].boxArray(), info));
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

    std::unique_ptr<MLPoisson> mlpoisson_ptr = aniso
        ? std::make_unique<MLPoissonAniso>(geom_tmp, ba_tmp, dm_tmp, info)
        : std::make_unique<MLPoisson>(geom_tmp, ba_tmp, dm_tmp, info);
    MLPoisson& mlpoisson = *mlpoisson_ptr;
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
