/**
 * \file ERF_SolveWithTerrainMLMG.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLTerrainPoisson.H>

#include <string>

using namespace amrex;

/**
 * Solve the Poisson equation on a terrain-fitted mesh with amrex::MLTerrainPoisson
 * (AMReX-Codes/amrex#6030): the stencil of TerrainPoisson with zero flux through every
 * Neumann face, solved by multigrid that coarsens in the horizontal only and solves the
 * columns exactly, alone (erf.terrain_poisson_solver = mlmg) or as the preconditioner of
 * GMRES (gmres_mlmg).  It needs no FFT.
 *
 * The level is solved on its own region, as solve_with_gmres does: the faces of the region
 * inside the domain (coarse/fine faces) are Neumann faces; a domain face is Dirichlet only if
 * it is an Outflow or Open face (a lateral one only without real boundary conditions), as in
 * get_fft_bc; a direction is periodic only if the region spans it.  A region that reaches one
 * end of a periodic direction only is given a Neumann face there (coarse cells lie across
 * it), where get_fft_bc says periodic and the FFT preconditioner stops; project_momenta stops
 * a level whose regions meet across a periodic boundary.
 *
 * @param lev Level index for the solve
 * @param subdomain Box over which the solve is performed
 * @param rhs Right-hand side field for the Poisson solve
 * @param phi Solution field to fill
 * @param fluxes Face-centered fluxes (minus the terrain-following gradient) to fill
 * @param ax_sub Terrain metric coefficient on x-faces (map-factor scaled)
 * @param ay_sub Terrain metric coefficient on y-faces (map-factor scaled)
 * @param az_sub Terrain metric coefficient on z-faces (map-factor scaled)
 * @param dJ_sub Cell-centered Jacobian determinant
 * @param znd_sub Node-centered physical height field
 * @param use_gmres True to run GMRES with the multigrid as its preconditioner
 */
void ERF::solve_with_terrain_mlmg (int lev, const Box& subdomain, MultiFab& rhs, MultiFab& phi,
                                   Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                   MultiFab& ax_sub, MultiFab& ay_sub, MultiFab& az_sub,
                                   MultiFab& dJ_sub, MultiFab& znd_sub, bool use_gmres)
{
    BL_PROFILE("ERF::solve_with_terrain_mlmg()");

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    const Geometry& lev_geom = Geom(lev);
    const Box& domain = lev_geom.Domain();
    auto const dom_lo = lbound(domain);
    auto const dom_hi = ubound(domain);
    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);
    auto dx = lev_geom.CellSizeArray();

    // The geometry of the region solved on, built as in solve_with_gmres
    Array<int,AMREX_SPACEDIM> is_per; is_per[0] = 0; is_per[1] = 0; is_per[2] = 0;
    if (lev_geom.isPeriodic(0) && sub_lo.x == dom_lo.x && sub_hi.x == dom_hi.x) { is_per[0] = 1; }
    if (lev_geom.isPeriodic(1) && sub_lo.y == dom_lo.y && sub_hi.y == dom_hi.y) { is_per[1] = 1; }
    Geometry my_geom;
    if (subdomain == domain) {
        my_geom.define(domain, lev_geom.ProbDomain(), 0, is_per);
    } else {
        RealBox rb( sub_lo.x   *dx[0],  sub_lo.y   *dx[1],  sub_lo.z   *dx[2],
                   (sub_hi.x+1)*dx[0], (sub_hi.y+1)*dx[1], (sub_hi.z+1)*dx[2]);
        my_geom.define(subdomain, rb, 0, is_per);
    }

    // Boundary conditions as get_fft_bc sets them for the FFT preconditioner
    Array<LinOpBCType,AMREX_SPACEDIM> bclo, bchi;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (my_geom.isPeriodic(dir)) {
            bclo[dir] = LinOpBCType::Periodic;
            bchi[dir] = LinOpBCType::Periodic;
            continue;
        }
        const std::string& type_lo = domain_bc_type[Orientation(dir,Orientation::low)];
        const std::string& type_hi = domain_bc_type[Orientation(dir,Orientation::high)];
        const bool on_lo = (subdomain.smallEnd(dir) == domain.smallEnd(dir));
        const bool on_hi = (subdomain.bigEnd(dir)   == domain.bigEnd(dir));
        const bool real_ok = (dir == 2 || !solverChoice.use_real_bcs);
        bclo[dir] = (on_lo && real_ok && (type_lo == "Outflow" || type_lo == "Open")) ? LinOpBCType::Dirichlet
                                                                                      : LinOpBCType::Neumann;
        bchi[dir] = (on_hi && real_ok && (type_hi == "Outflow" || type_hi == "Open")) ? LinOpBCType::Dirichlet
                                                                                      : LinOpBCType::Neumann;
    }

    // A region one cell wide in x or y is solved as two-dimensional, as solve_with_mlmg does,
    // but only where that direction is periodic: MLTerrainPoisson then drops the face fluxes
    // and the terrain cross terms in it, which are zero only when the one cell is its own
    // periodic neighbour (its two node columns are the same).  Across a wall, inflow or
    // outflow face they are not, so such a region keeps all three directions.
    LPInfo info;
    if (sub_lo.x == sub_hi.x && is_per[0]) {
        info.setHiddenDirection(0);
    } else if (sub_lo.y == sub_hi.y && is_per[1]) {
        info.setHiddenDirection(1);
    }

    Vector<Geometry>            geom_tmp; geom_tmp.push_back(my_geom);
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(rhs.boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(rhs.DistributionMap());

    MLTerrainPoisson linop(geom_tmp, ba_tmp, dm_tmp, info);
    linop.setDomainBC(bclo, bchi);
    linop.setLevelBC(0, nullptr);
    linop.setZPhys(0, znd_sub);
    linop.setAreas(0, {&ax_sub, &ay_sub, &az_sub});
    linop.setDetJ(0, dJ_sub);

    MLMG mlmg(linop);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);
    mlmg.setMaxIter(100);

    if (use_gmres) {
        GMRESMLMG gmres(mlmg);
        gmres.usePrecond(true);
        gmres.setVerbose(mg_verbose);
        gmres.setMaxIters(100);
        gmres.solve(phi, rhs, reltol, abstol);
    } else {
        // MLMG stops the run when it does not converge (GMRES carries on with its last
        // iterate); say which inputs to change instead of AMReX's bare message
        mlmg.setThrowException(true);
        try {
            mlmg.solve({&phi}, {&rhs}, reltol, abstol);
        } catch (MLMG::error const& err) {
            amrex::Abort(std::string("erf.terrain_poisson_solver = mlmg at level ") + std::to_string(lev)
                         + ": " + err.what() + "  Loosen erf.poisson_reltol / erf.poisson_abstol, or use "
                         + "erf.terrain_poisson_solver = gmres_mlmg or gmres_fft, which continue with "
                         + "the last iterate.");
        }
    }

    // Minus the terrain-following gradient, zero on the Neumann faces: the convention of
    // TerrainPoisson::getFluxes.  The solution is passed in, since GMRES leaves its own
    // answer in phi rather than in the multigrid's solution.
    // The caller (project_momenta) multiplies the lateral fluxes by the map factor.
    mlmg.getFluxes({{&fluxes[0], &fluxes[1], &fluxes[2]}}, {&phi});

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
}
