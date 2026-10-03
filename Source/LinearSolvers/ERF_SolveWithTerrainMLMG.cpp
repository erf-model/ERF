/**
 * \file ERF_SolveWithTerrainMLMG.cpp
 */
#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_SolverUtils.H"
#include "ERF_MLTerrainPoisson.H"

#include <AMReX_MLMG.H>
#include <AMReX_GMRES_MLMG.H>
#include <memory>

using namespace amrex;

/**
 * A terrain multigrid operator kept between projections.  The metrics of a static
 * terrain do not change, so the coarse levels, ghost source maps and column coefficients
 * are built once; the operator is keyed on the layout of the solve and rebuilt when the
 * level is remade.  The finest-level metrics are aliases of ERF's own arrays (including
 * the map-factor scaling the projection applies to them before every solve).
 */
struct TerrainMLMGCache
{
    BoxArray ba;                        //!< layout of the solve this operator was built for
    DistributionMapping dm;
    bool regrid = false;                //!< whether the solve runs on re-gridded boxes
    BoxArray ba_mg;                     //!< the operator's layout (ba when not re-gridded)
    DistributionMapping dm_mg;
    MultiFab znd_mg, ax_mg, ay_mg, az_mg, dJ_mg;   //!< metrics copied onto the re-gridded layout
    std::unique_ptr<MLTerrainPoisson> op;
    std::unique_ptr<MLMG> mlmg;
};

/**
 * Solve the terrain-fitted Poisson equation with geometric multigrid on the full
 * terrain stencil, either as the solver (erf.terrain_poisson_solver = mlmg) or as a
 * fixed number of V-cycles preconditioning GMRES on the same operator
 * (erf.terrain_poisson_solver = gmres_mlmg, amrex::GMRESMLMG).  The fluxes come from the same
 * stencil in both cases, so the corrected velocity is discretely divergence-free
 * for the terrain operator.
 *
 * The multigrid can only coarsen as far as every box of the level allows.  When the
 * level's boxes stop the coarsening short of what the domain allows (a box 150
 * cells wide allows one coarsening, the 600-cell domain it tiles allows three) the
 * solve is re-gridded onto boxes whose extents are multiples of the full
 * coarsening ratio, the metrics and right-hand side are copied over and the
 * solution and fluxes copied back; the operator does not depend on the layout, so
 * only the solver's own round-off changes.  The operator (coarse metrics, ghost maps,
 * column coefficients) is kept between projections of a static terrain.
 *
 * @param lev Level index for the solve
 * @param subdomain Box over which the solve is performed
 * @param rhs Right-hand side field for the Poisson solve
 * @param phi Solution field to fill
 * @param fluxes Face-centered gradient fluxes to fill
 * @param ax_sub Terrain metric coefficient on x-faces (map-factor scaled)
 * @param ay_sub Terrain metric coefficient on y-faces (map-factor scaled)
 * @param az_sub Terrain metric coefficient on z-faces (map-factor scaled)
 * @param dJ_sub Cell-centered Jacobian determinant
 * @param znd_sub Node-centered physical height field
 * @param use_gmres True to run GMRES with multigrid as the preconditioner
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

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }
    const int hd = info.hidden_direction;

    // ****************************************************************************
    // How deep the multigrid can go: once a direction (typically z, with few cells)
    // stops being coarsenable the others continue alone (semicoarsening), which the
    // column smoother is made for.  AMReX does not combine that with a hidden direction.
    // ****************************************************************************
    const IntVect depth_dom = MLTerrainPoisson::coarsening_depth(domain, hd);
    if (hd < 0) {
        info.setSemicoarsening(true);
        info.setMaxSemicoarseningLevel(depth_dom.max());
        info.setSemicoarseningDirection(-1);
    }

    // ****************************************************************************
    // The operator: taken from the cache when one was built for this layout, else built
    // and kept.  A moving terrain changes the metrics every step and is never cached.
    // ****************************************************************************
    const BoxArray&            ba_lev = rhs.boxArray();
    const DistributionMapping& dm_lev = rhs.DistributionMap();
    const bool cacheable = (solverChoice.terrain_type != TerrainType::MovingFittedMesh);

    std::shared_ptr<TerrainMLMGCache> cache;
    if (cacheable) {
        for (auto& c : terrain_mlmg_cache[lev]) {
            if (c->ba == ba_lev && c->dm == dm_lev) { cache = c; break; }
        }
    }

    const Periodicity& period = lev_geom.periodicity();
    const IntVect ng1(1);
    const IntVect ng0(0);

    if (!cache)
    {
        cache = std::make_shared<TerrainMLMGCache>();
        cache->ba = ba_lev;
        cache->dm = dm_lev;

        // The layout of the solve: the level's boxes, or multigrid-friendly ones when
        // the level's boxes would stop the coarsening early in some direction
        const bool covers_domain = (ba_lev.minimalBox() == domain) && (ba_lev.numPts() == domain.numPts());
        const IntVect depth_box = MLTerrainPoisson::coarsening_depth(ba_lev, hd);
        cache->regrid = covers_domain && !depth_box.allGE(depth_dom);

        // A DistributionMapping copy shares its data, so the re-gridded one is built fresh
        cache->ba_mg = ba_lev;
        cache->dm_mg = dm_lev;
        if (cache->regrid) {
            cache->ba_mg = MLTerrainPoisson::multigrid_grids(domain, hd, depth_dom, 64);
            cache->dm_mg = DistributionMapping(cache->ba_mg);
            if (mg_verbose > 0) {
                amrex::Print() << "Terrain multigrid: the level's boxes allow " << depth_box
                               << " coarsenings and the domain " << depth_dom
                               << "; solving on " << cache->ba_mg.size() << " re-gridded boxes" << std::endl;
            }
        }

        const MultiFab* p_znd = &znd_sub;
        const MultiFab* p_ax  = &ax_sub;
        const MultiFab* p_ay  = &ay_sub;
        const MultiFab* p_az  = &az_sub;
        const MultiFab* p_dJ  = &dJ_sub;
        if (cache->regrid) {
            const BoxArray& ba_mg = cache->ba_mg;
            const DistributionMapping& dm_mg = cache->dm_mg;
            cache->znd_mg.define(amrex::convert(ba_mg, IntVect(1)),     dm_mg, 1, 1);
            cache->ax_mg.define (amrex::convert(ba_mg, IntVect(1,0,0)), dm_mg, 1, 1);
            cache->ay_mg.define (amrex::convert(ba_mg, IntVect(0,1,0)), dm_mg, 1, 1);
            cache->az_mg.define (amrex::convert(ba_mg, IntVect(0,0,1)), dm_mg, 1, 1);
            cache->dJ_mg.define (ba_mg, dm_mg, 1, 1);
            // The ghost values of the metrics outside the domain are ERF's (clamped and
            // extrapolated nodes), copied along so the operator is the same to the bit
            cache->znd_mg.ParallelCopy(znd_sub, 0, 0, 1, ng1, ng1, period);
            cache->ax_mg.ParallelCopy (ax_sub,  0, 0, 1, ng1, ng1, period);
            cache->ay_mg.ParallelCopy (ay_sub,  0, 0, 1, ng1, ng1, period);
            cache->az_mg.ParallelCopy (az_sub,  0, 0, 1, ng1, ng1, period);
            cache->dJ_mg.ParallelCopy (dJ_sub,  0, 0, 1, ng1, ng1, period);
            p_znd = &cache->znd_mg; p_ax = &cache->ax_mg; p_ay = &cache->ay_mg;
            p_az = &cache->az_mg; p_dJ = &cache->dJ_mg;
        }

        Vector<Geometry>            geom_tmp; geom_tmp.push_back(lev_geom);
        Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(cache->ba_mg);
        Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(cache->dm_mg);

        cache->op = std::make_unique<MLTerrainPoisson>(geom_tmp, ba_tmp, dm_tmp, info,
                                                       *p_znd, *p_ax, *p_ay, *p_az, *p_dJ);
        MLTerrainPoisson& mlop = *cache->op;

        Array<LinOpBCType,AMREX_SPACEDIM> bclo;
        Array<LinOpBCType,AMREX_SPACEDIM> bchi;
        get_terrain_projection_bc(lev_geom, domain_bc_type, solverChoice.use_real_bcs, bclo, bchi);
        mlop.setDomainBC(bclo, bchi);

        if (lev > 0) {
            mlop.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
        }
        mlop.setLevelBC(0, nullptr);

        // Second-order Dirichlet fill, the same as the GMRES operator
        mlop.setMaxOrder(2);

        cache->mlmg = std::make_unique<MLMG>(mlop);
        cache->mlmg->setBottomVerbose(0);

        if (cacheable) {
            terrain_mlmg_cache[lev].push_back(cache);
        }
    }

    const bool regrid = cache->regrid;
    const BoxArray& ba_mg = cache->ba_mg;
    const DistributionMapping& dm_mg = cache->dm_mg;
    MLMG& mlmg = *cache->mlmg;

    // The vectors of this solve, on the operator's layout
    MultiFab rhs_mg, phi_mg;
    Array<MultiFab,AMREX_SPACEDIM> flux_mg;
    MultiFab* p_rhs = &rhs;
    MultiFab* p_phi = &phi;
    Array<MultiFab*,AMREX_SPACEDIM> p_flux = {&fluxes[0], &fluxes[1], &fluxes[2]};
    if (regrid) {
        rhs_mg.define(ba_mg, dm_mg, 1, 0);
        phi_mg.define(ba_mg, dm_mg, 1, 1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            flux_mg[idim].define(amrex::convert(ba_mg, IntVect::TheDimensionVector(idim)), dm_mg, 1, 0);
        }
        rhs_mg.ParallelCopy(rhs, 0, 0, 1, ng0, ng0, period);
        phi_mg.setVal(0.0);
        phi_mg.ParallelCopy(phi, 0, 0, 1, ng0, ng0, period);
        p_rhs = &rhs_mg; p_phi = &phi_mg;
        p_flux = {&flux_mg[0], &flux_mg[1], &flux_mg[2]};
    }

    if (mg_verbose > 0) {
        amrex::Print() << "Solving the terrain Poisson equation with "
                       << (use_gmres ? "GMRES preconditioned by MLTerrainPoisson V-cycles"
                                     : "MLTerrainPoisson multigrid") << std::endl;
    }

    mlmg.setMaxIter(solverChoice.terrain_mlmg_max_iter);
    mlmg.setPreSmooth(solverChoice.terrain_mlmg_smooth_sweeps);
    mlmg.setPostSmooth(solverChoice.terrain_mlmg_smooth_sweeps);

    if (!use_gmres)
    {
        mlmg.setVerbose(mg_verbose);
        mlmg.solve({p_phi}, {p_rhs}, reltol, abstol);
        mlmg.getFluxes({p_flux}, {p_phi}, MLMG::Location::FaceCenter);

        if (regrid) {
            phi.ParallelCopy(phi_mg, 0, 0, 1, ng0, ng0, period);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fluxes[idim].ParallelCopy(flux_mg[idim], 0, 0, 1, ng0, ng0, period);
            }
        }
    }
    else
    {
        // GMRES on the multigrid operator with its V-cycles as the preconditioner
        // (amrex::GMRESMLMG).  The operator is the same stencil as TerrainPoisson's on the
        // finest level, and on a refined level it carries the zero-flux coarse/fine
        // condition the union of boxes needs, so the system GMRES solves is consistent.
        // The wrapper makes the preconditioner a fixed linear operator by using a fixed
        // number of smoothing sweeps as the bottom solve: a bottom solve iterated to a
        // tolerance would make the residual GMRES tracks drift from the true one.
        mlmg.setVerbose(amrex::max(0, mg_verbose-1));
        GMRESMLMG gmres(mlmg);
        gmres.usePrecond(true);
        gmres.setPrecondNumIters(solverChoice.terrain_mlmg_precond_iters);
        gmres.setVerbose(mg_verbose);
        gmres.setMaxIters(solverChoice.terrain_mlmg_max_iter);
        gmres.solve(*p_phi, *p_rhs, reltol, abstol);
        if (mg_verbose > 0) {
            amrex::Print() << "GMRES (terrain multigrid preconditioner): " << gmres.getNumIters()
                           << " iterations, residual " << gmres.getResidualNorm() << std::endl;
        }
        mlmg.getFluxes({p_flux}, {p_phi}, MLMG::Location::FaceCenter);

        if (regrid) {
            phi.ParallelCopy(phi_mg, 0, 0, 1, ng0, ng0, period);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                fluxes[idim].ParallelCopy(flux_mg[idim], 0, 0, 1, ng0, ng0, period);
            }
        }
    }

    // The lateral fluxes carry the map factor, as in solve_with_gmres
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        const Array4<Real      >& fx_ar = fluxes[0].array(mfi);
        const Array4<Real      >& fy_ar = fluxes[1].array(mfi);
        const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
        ParallelFor(xbx,ybx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_ar(i,j,k) *= mf_ux(i,j,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_ar(i,j,k) *= mf_vy(i,j,0);
        });
    } // mfi

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
}
