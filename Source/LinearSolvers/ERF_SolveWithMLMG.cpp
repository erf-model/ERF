/**
 * \file ERF_SolveWithMLMG.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#include <limits>
#include <memory>
#include <string>

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
 * MLABecLaplacian whose MG levels coarsen each direction on its own, as MLPoissonAniso.
 */
class MLABecLaplacianAniso : public MLABecLaplacian {
public:
    MLABecLaplacianAniso (const Vector<Geometry>& a_geom, const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap, const LPInfo& a_info)
    {
        mg_independent_coarsening = true;
        define(a_geom, a_grids, a_dmap, a_info);
    }
};

/**
 * Whether the cells of geom are anisotropic enough (a direction at least 1.5 times longer
 * than the shortest) for the MG levels to coarsen each direction on its own (#5064)
 */
bool
cells_are_anisotropic (const Geometry& geom, const LPInfo& info)
{
    const Real* dx = geom.CellSize();
    Real hmin = std::numeric_limits<Real>::max();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim != info.hidden_direction) { hmin = std::min(hmin, dx[idim]); }
    }
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim != info.hidden_direction && dx[idim] >= Real(1.5)*hmin) { return true; }
    }
    return false;
}

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
    const bool aniso = cells_are_anisotropic(geom, info);
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

/**
 * Solve, on one region of a flat mesh (ConstantDz or StretchedDz), the equation that
 * solve_with_terrain_mlmg solves, with amrex::MLABecLaplacian: MLTerrainPoisson needs grids
 * that cover the region's bounding box, MLABecLaplacian does not.  Used for a region whose boxes
 * do not make a rectangle, where the terrain solver and the hybrid FFT solver cannot be used.
 *
 * On a flat mesh MLTerrainPoisson has no cross terms: it solves
 *     div(area grad phi) = detJ rhs
 * with the x- and y-faces weighted by ax_sub and ay_sub, and the z-gradient taken over the
 * physical distance between cell centers.  That is MLABecLaplacian with alpha = 0 and
 *     B_x = ax_sub,  B_y = ay_sub,  B_z = az_sub dz / (distance between the centers),
 * dz being the cell size of the level's Geometry.  The fluxes it returns are -B grad phi,
 * with the gradient over the Geometry's cell size; they are divided by the area here, so that,
 * as from solve_with_terrain_mlmg and solve_with_fft, they are minus the gradient.
 *
 * The level's whole Geometry is used, as in solve_with_mlmg: the faces of the region that are
 * not domain faces are coarse/fine faces with zero flux, as the faces of a region are in the
 * other solvers, and a domain face is Dirichlet if it is an Outflow or Open face.
 *
 * @param lev Level index for the solve
 * @param subdomain Bounding box of the region
 * @param rhs Right-hand side, as passed to solve_with_terrain_mlmg (not changed)
 * @param phi Solution field to fill
 * @param fluxes Face-centered fluxes (minus the gradient) to fill
 * @param ax_sub Area weight on x-faces, with the map factors (and beta) of the solve
 * @param ay_sub Area weight on y-faces, likewise
 * @param az_sub Area weight on z-faces, likewise
 * @param dJ_sub Cell-centered Jacobian determinant
 * @param znd_sub Node-centered physical height
 */
void
ERF::solve_with_flat_mlmg (int lev, const Box& subdomain, MultiFab& rhs, MultiFab& phi,
                           Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                           const MultiFab& ax_sub, const MultiFab& ay_sub, const MultiFab& az_sub,
                           const MultiFab& dJ_sub, const MultiFab& znd_sub)
{
    BL_PROFILE("ERF::solve_with_flat_mlmg()");

    AMREX_ALWAYS_ASSERT(solverChoice.mesh_type != MeshType::VariableDz);

    const Geometry& lev_geom = Geom(lev);
    const Box& domain = lev_geom.Domain();
    const BoxArray& ba = rhs.boxArray();
    const DistributionMapping& dm = rhs.DistributionMap();

    // A direction one cell wide is dropped only where it is periodic, as in
    // solve_with_terrain_mlmg
    LPInfo info;
    if (domain.length(0) == 1 && lev_geom.isPeriodic(0)) {
        info.setHiddenDirection(0);
    } else if (domain.length(1) == 1 && lev_geom.isPeriodic(1)) {
        info.setHiddenDirection(1);
    }
    const bool aniso = cells_are_anisotropic(lev_geom, info);
    if (aniso) {
        info.setMaxCoarseningLevel(aniso_max_coarsening_level(lev_geom, ba, info));
    }

    // The face coefficients
    Array<MultiFab,AMREX_SPACEDIM> bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        bcoef[idim].define(convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 0);
    }
    MultiFab::Copy(bcoef[0], ax_sub, 0, 0, 1, 0);
    MultiFab::Copy(bcoef[1], ay_sub, 0, 0, 1, 0);
    // MLABecLaplacian keeps the faces of a hidden direction in its stencil (MLPoisson and
    // MLTerrainPoisson drop them).  Their flux is zero (the one cell is its own periodic
    // neighbor), but they would still add to the diagonal of the smoother, and on the coarser MG
    // levels, which do not coarsen that direction, they would swamp it: multigrid diverges.
    if (info.hasHiddenDimension()) {
        bcoef[info.hidden_direction].setVal(Real(0.0));
    }
    const Real dz = lev_geom.CellSize(2);
    for (MFIter mfi(bcoef[2]); mfi.isValid(); ++mfi) {
        const Array4<Real      >& bz  = bcoef[2].array(mfi);
        const Array4<Real const>& az_ = az_sub.const_array(mfi);
        const Array4<Real const>& zp  = znd_sub.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // The distance between the centers of cells k-1 and k, as MLTerrainPoisson's
            // metric_zf takes it
            const Real dzc = Real(0.125) * ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                            -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
            bz(i,j,k) = az_(i,j,k) * dz / dzc;
        });
    }

    // -div(B grad phi) = -detJ rhs
    MultiFab rhs_mlabec(ba, dm, 1, 0);
    MultiFab::Copy(rhs_mlabec, rhs, 0, 0, 1, 0);
    MultiFab::Multiply(rhs_mlabec, dJ_sub, 0, 0, 1, 0);
    rhs_mlabec.mult(Real(-1.0));

    Vector<Geometry>            geom_tmp; geom_tmp.push_back(lev_geom);
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(ba);
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(dm);

    std::unique_ptr<MLABecLaplacian> mlabec_ptr = aniso
        ? std::make_unique<MLABecLaplacianAniso>(geom_tmp, ba_tmp, dm_tmp, info)
        : std::make_unique<MLABecLaplacian>(geom_tmp, ba_tmp, dm_tmp, info);
    MLABecLaplacian& mlabec = *mlabec_ptr;
    mlabec.setDomainBC(get_lo_projection_bc(lev_geom, domain_bc_type),
                       get_hi_projection_bc(lev_geom, domain_bc_type));
    if (lev > 0) {
        mlabec.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    mlabec.setLevelBC(0, nullptr);
    // A Dirichlet face is half a cell from the center of the cell next to it, as in the
    // other solvers
    mlabec.setMaxOrder(2);
    mlabec.setScalars(Real(0.0), Real(1.0));
    mlabec.setBCoeffs(0, GetArrOfConstPtrs(bcoef));

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(100);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve({&phi}, {&rhs_mlabec}, solverChoice.poisson_reltol, solverChoice.poisson_abstol);
    mlmg.getFluxes({GetArrOfPtrs(fluxes)});

    // Minus the gradient, as the other solvers return it (zero in a hidden direction): the
    // fluxes are -B grad phi with grad phi on the Geometry's cell size, so -area times the
    // gradient over the physical distance
    const Array<MultiFab const*,AMREX_SPACEDIM> area{&ax_sub, &ay_sub, &az_sub};
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim == info.hidden_direction) {
            fluxes[idim].setVal(Real(0.0));
        } else {
            MultiFab::Divide(fluxes[idim], *area[idim], 0, 0, 1, 0);
        }
    }

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
}

/**
 * Multigrid solve of div(sigma grad phi) = rhs (erf.if_implicit_projection) on a constant-dz
 * mesh with map factors of one, for a build without FFT or a level whose regions are not
 * rectangles (erf-model/ERF#4887). As solve_with_mlmg, it takes the boxes of one region and the
 * geometry of the level, with zero normal gradient on the coarse/fine boundary.
 *
 * @param lev    Level index for the solve
 * @param rhs    Right-hand side, div(sigma m*) (not changed)
 * @param phi    Solution field to fill
 * @param fluxes Face fluxes to fill: -grad phi, without sigma
 * @param sigma  Face coefficients (make_if_projection_sigma)
 */
void
ERF::solve_with_if_mlmg (int lev, Vector<MultiFab>& rhs, Vector<MultiFab>& phi,
                         Vector<Array<MultiFab,AMREX_SPACEDIM>>& fluxes,
                         Array<MultiFab const*,AMREX_SPACEDIM> const& sigma)
{
    BL_PROFILE("ERF::solve_with_if_mlmg()");

    const Geometry& lev_geom = Geom(lev);
    auto const dom_lo = lbound(lev_geom.Domain());
    auto const dom_hi = ubound(lev_geom.Domain());
    const BoxArray& ba = rhs[0].boxArray();
    const DistributionMapping& dm = rhs[0].DistributionMap();

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }

    // Anisotropic cells: coarsen the directions on their own, as solve_with_mlmg does
    const bool aniso = cells_are_anisotropic(lev_geom, info);
    if (aniso) {
        info.setMaxCoarseningLevel(aniso_max_coarsening_level(lev_geom, ba, info));
    }

    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(ba);
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(dm);
    Vector<Geometry>            geom_tmp; geom_tmp.push_back(lev_geom);

    std::unique_ptr<MLABecLaplacian> linop_ptr = aniso
        ? std::make_unique<MLABecLaplacianAniso>(geom_tmp, ba_tmp, dm_tmp, info)
        : std::make_unique<MLABecLaplacian>(geom_tmp, ba_tmp, dm_tmp, info);
    MLABecLaplacian& linop = *linop_ptr;
    linop.setMaxOrder(2);
    linop.setDomainBC(get_lo_projection_bc(lev_geom, domain_bc_type),
                      get_hi_projection_bc(lev_geom, domain_bc_type));
    if (lev > 0) {
        linop.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    linop.setLevelBC(0, nullptr);

    // (A alpha - B div(b grad)) phi = rhs with A = 0, B = 1 and b = sigma: -div(sigma grad phi) = -rhs
    linop.setScalars(Real(0.0), Real(1.0));
    Array<MultiFab,AMREX_SPACEDIM> bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        bcoef[idim].define(sigma[idim]->boxArray(), sigma[idim]->DistributionMap(), 1, 0);
        MultiFab::Copy(bcoef[idim], *sigma[idim], 0, 0, 1, 0);
    }
    // MLABecLaplacian keeps the faces of a hidden direction in its stencil; see solve_with_flat_mlmg
    if (info.hasHiddenDimension()) {
        bcoef[info.hidden_direction].setVal(Real(0.0));
    }
    linop.setBCoeffs(0, GetArrOfConstPtrs(bcoef));

    MultiFab rhs_mlabec(ba, dm, 1, 0);
    MultiFab::Copy(rhs_mlabec, rhs[0], 0, 0, 1, 0);
    rhs_mlabec.mult(Real(-1.0));

    MLMG mlmg(linop);
    mlmg.setMaxIter(100);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);
    // As in solve_with_EB_mlmg: stop the bottom solve well below the tolerance of the solve
    mlmg.setBottomToleranceAbs(Real(1.e-4) * solverChoice.poisson_abstol);

    // MLMG stops the run when it does not converge; name the option and the inputs to change
    mlmg.setThrowException(true);
    try {
        mlmg.solve({&phi[0]}, {&rhs_mlabec},
                   solverChoice.poisson_reltol, solverChoice.poisson_abstol);
    } catch (MLMG::error const& err) {
        amrex::Abort(std::string("erf.if_implicit_projection, multigrid at level ") + std::to_string(lev)
                     + ": " + err.what() + "  Loosen erf.poisson_reltol / erf.poisson_abstol.");
    }

    // Return -grad phi, as solve_with_if_gmres does: the gradient of the solution, not the flux
    // sigma grad phi, so that no face is divided by sigma (zero in a hidden direction)
    mlmg.getGradSolution({GetArrOfPtrs(fluxes[0])});
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim == info.hidden_direction) {
            fluxes[0][idim].setVal(Real(0.0));
        } else {
            fluxes[0][idim].mult(Real(-1.0));
        }
    }
}
