#include "ERF.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include "Utils.H"

#ifdef ERF_USE_POISSON_SOLVE

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
            //if (bc_type == "Outflow") {
            //    r[dir] = LinOpBCType::Dirichlet;
            //} else
            {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    //r[2] = LinOpBCType::Neumann;
    //amrex::Print() << "BCs for Poisson solve " << r[0] << " " << r[1] << " " << r[2] << std::endl;
    return r;
}
bool ERF::projection_has_dirichlet (Array<LinOpBCType,AMREX_SPACEDIM> bcs) const
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (bcs[dir] == LinOpBCType::Dirichlet) return true;
    }
    return false;
}


/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::project_velocities (int lev, Real l_dt, Vector<MultiFab>& vmf, MultiFab& pmf)
{
    AMREX_ALWAYS_ASSERT(!solverChoice.use_terrain);
    Vector<Vector<MultiFab>> tmpu(1);
    for (auto& mf : vmf) {
        tmpu[0].emplace_back(mf, amrex::make_alias, 0, mf.nComp());
    }

    Vector       <MultiFab> tmpp;
    tmpp.push_back(MultiFab(pmf, amrex::make_alias, 0, 1));

    project_velocities(lev,lev,l_dt,tmpu,tmpp);
}

/**
 * Project the multi-level velocity field to enforce incompressibility
 * Both vars and p come in with levels lev_min to lev_max, regardless of the number of levels in the hierarchy overall
 */
void
ERF::project_velocities (int lev_min, int lev_max, Real l_dt,
                         Vector<Vector<MultiFab>>& vars, Vector<MultiFab>& p)
{
    AMREX_ALWAYS_ASSERT(!solverChoice.use_terrain);
    BL_PROFILE("ERF::project_velocities()");

    const int nlevs = lev_max - lev_min + 1;

    // Make sure the solver only sees the levels over which we are solving
    LPInfo info;
    Vector<BoxArray> ba_tmp;
    Vector<DistributionMapping> dm_tmp;
    Vector<Geometry> geom_tmp;
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        amrex::Print() << "BA FOR POISSON SOLVE " << vars[ilev][Vars::cons].boxArray() << std::endl;
        ba_tmp.push_back(vars[ilev][Vars::cons].boxArray());
        dm_tmp.push_back(vars[ilev][Vars::cons].DistributionMap());
        geom_tmp.push_back(geom[lev_min+ilev]);
    }
    MLABecLaplacian mlabec(geom_tmp, ba_tmp, dm_tmp, info);

    //
    // This will hold (1/rho) on faces
    //
    Vector<Array<MultiFab,AMREX_SPACEDIM> > inv_rho(lev_max-lev_min+1);

    //
    // The operator is (alpha A - beta del dot B grad) phi = RHS
    // Here we set alpha to 0 and beta to -1
    // Then b is (dt/rho)
    //
    mlabec.setScalars(0.0, -1.0);
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        inv_rho[ilev][0].define(vars[ilev][Vars::xvel].boxArray(),dm_tmp[ilev],1,0,MFInfo());
        inv_rho[ilev][1].define(vars[ilev][Vars::yvel].boxArray(),dm_tmp[ilev],1,0,MFInfo());
        inv_rho[ilev][2].define(vars[ilev][Vars::zvel].boxArray(),dm_tmp[ilev],1,0,MFInfo());

#if 0
        MultiFab density(vars[ilev][Vars::cons], make_alias, 1, 1);
        density.FillBoundary(geom[ilev].periodicity());
        amrex::average_cellcenter_to_face(GetArrOfPtrs(inv_rho[ilev-lev_min]), density, geom[ilev]);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            inv_rho[ilev-lev_min][idim].invert(l_dt, 0);
        }
#else
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            inv_rho[ilev][idim].setVal(l_dt);
        }
#endif

        mlabec.setBCoeffs(ilev, GetArrOfConstPtrs(inv_rho[ilev]));
    }

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    bool need_adjust_rhs = (projection_has_dirichlet(bclo) || projection_has_dirichlet(bchi)) ? false : true;
    mlabec.setDomainBC(bclo, bchi);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
       if (lev_min+ilev > 0) {
           // mlabec.setCoarseFineBC(nullptr, ref_ratio[lev_min+ilev-1], LinOpBCType::Neumann);
           mlabec.setCoarseFineBC(&pp_inc[lev_min+ilev-1], ref_ratio[lev_min+ilev-1]);
       }
       mlabec.setLevelBC(ilev, nullptr);
    }

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

    rhs.resize(nlevs);
    phi.resize(nlevs);
    fluxes.resize(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        rhs[ilev].define(grids[lev_min+ilev], dmap[lev_min+ilev], 1, 0);
        phi[ilev].define(grids[lev_min+ilev], dmap[lev_min+ilev], 1, 1);
        rhs[ilev].setVal(0.0);
        phi[ilev].setVal(0.0);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            fluxes[ilev][idim].define(
                convert(grids[lev_min+ilev], IntVect::TheDimensionVector(idim)),
                dmap[lev_min+ilev], 1, 0);
        }
    }

    // Define a single Array of MultiFabs to hold the velocity components
    // Then define the RHS to be the divergence of the current velocities
    Array<MultiFab const*, AMREX_SPACEDIM> u;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        u[0] = &(vars[ilev][Vars::xvel]);
        u[1] = &(vars[ilev][Vars::yvel]);
        u[2] = &(vars[ilev][Vars::zvel]);
        computeDivergence(rhs[ilev], u, geom[ilev]);

        // If all Neumann BCs, adjust RHS to make sure we can converge
        if (need_adjust_rhs) {
            Real offset = volWgtSumMF(lev_min+ilev, rhs[ilev], 0, *mapfac_m[lev_min+ilev], false, false);
            amrex::Print() << "Poisson solvability offset = " << offset << std::endl;
            rhs[ilev].plus(-offset, 0, 1);
        }
    }

    // Initialize phi to 0
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        phi[ilev].setVal(0.0);
    }

    MLMG mlmg(mlabec);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    mlmg.setVerbose(mg_verbose);
    //mlmg.setBottomVerbose(mg_verbose);

    mlmg.solve(GetVecOfPtrs(phi),
               GetVecOfConstPtrs(rhs),
               solverChoice.poisson_reltol,
               solverChoice.poisson_abstol);

    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        // Update pressure variable with phi -- note that phi is change in pressure, not the full pressure
        MultiFab::Saxpy(p[ilev], 1.0, phi[ilev],0,0,1,0);
        p[ilev].FillBoundary(geom[ilev].periodicity());

        // Subtract (dt/rho) grad(phi) from the velocity components
        MultiFab::Add(vars[ilev][Vars::xvel], fluxes[ilev][0], 0,0,1,0);
        MultiFab::Add(vars[ilev][Vars::yvel], fluxes[ilev][1], 0,0,1,0);
        MultiFab::Add(vars[ilev][Vars::zvel], fluxes[ilev][2], 0,0,1,0);
    }

    // Average down the velocity from finest to coarsest to ensure consistency across levels
    Array<MultiFab const*, AMREX_SPACEDIM> u_fine;
    Array<MultiFab      *, AMREX_SPACEDIM> u_crse;
    for (int ilev = nlevs-1; ilev > 0; --ilev)
    {
        IntVect rr  = geom[ilev].Domain().size() / geom[ilev-1].Domain().size();
        u_fine[0] = &(vars[ilev  ][Vars::xvel]);
        u_fine[1] = &(vars[ilev  ][Vars::yvel]);
        u_fine[2] = &(vars[ilev  ][Vars::zvel]);
        u_crse[0] = &(vars[ilev-1][Vars::xvel]);
        u_crse[1] = &(vars[ilev-1][Vars::yvel]);
        u_crse[2] = &(vars[ilev-1][Vars::zvel]);
        average_down_faces(u_fine, u_crse, rr, geom[ilev-1]);
    }

#if 1
    // Confirm that the velocity is now divergence free
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        u[0] = &(vars[ilev][Vars::xvel]);
        u[1] = &(vars[ilev][Vars::yvel]);
        u[2] = &(vars[ilev][Vars::zvel]);
        rhs[ilev].setVal(0.0);
        computeDivergence(rhs[ilev], u, geom[lev_min+ilev]);
        Print() << "Max norm of divergence after solve at level " << lev_min+ilev << " : " << rhs[ilev].norm0() << std::endl;
        if (lev_min == 1) amrex::Print() << "RES " << rhs[ilev][0] << std::endl;
        if (lev_min == 1) exit(0);
    }
#endif
}
#endif
