#include "ERF.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include "Utils.H"

#ifdef ERF_USE_POISSON_SOLVE

using namespace amrex;

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 * if we want to enforce incompressibility of the initial conditions
 */

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
 * Note this assumes that we are projecting the level 0 field using domain bcs
 */
void ERF::project_velocities (Vector<MultiFab>& vmf)
{
    Vector<Vector<MultiFab>> tmpmf(1);
    for (auto& mf : vmf) {
        tmpmf[0].emplace_back(mf, amrex::make_alias, 0, mf.nComp());
    }
    project_velocities(0,0,tmpmf);
}

/**
 * Project the multi-level velocity field to enforce incompressibility
 */
void
ERF::project_velocities (int lev_min, int lev_max, Vector<Vector<MultiFab>>& vars)
{
    BL_PROFILE("ERF::project_velocities()");

    const int nlevs = lev_max - lev_min + 1;

    // Use the default settings
    LPInfo info;
    MLPoisson mlpoisson(geom, grids, dmap, info);

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    bool need_adjust_rhs = (projection_has_dirichlet(bclo) || projection_has_dirichlet(bchi)) ? false : true;
    mlpoisson.setDomainBC(bclo, bchi);

    for (int ilev = lev_min; ilev <= lev_max; ++ilev) {
       mlpoisson.setLevelBC(0, nullptr);
    }

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

    rhs.resize(nlevs);
    phi.resize(nlevs);
    fluxes.resize(nlevs);

    for (int ilev = lev_min; ilev <= lev_max; ++ilev) {
        rhs[ilev-lev_min].define(grids[ilev], dmap[ilev], 1, 0);
        phi[ilev-lev_min].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev-lev_min].setVal(0.0);
        phi[ilev-lev_min].setVal(0.0);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            fluxes[ilev-lev_min][idim].define(
                convert(grids[ilev], IntVect::TheDimensionVector(idim)),
                dmap[ilev], 1, 0);
        }
    }

    // Define a single Array of MultiFabs to hold the velocity components
    // Then define the RHS to be the divergence of the current velocities
    Array<MultiFab const*, AMREX_SPACEDIM> u;
    for (int ilev = lev_min; ilev <= lev_max; ++ilev)
    {
        u[0] = &(vars[ilev][Vars::xvel]);
        u[1] = &(vars[ilev][Vars::yvel]);
        u[2] = &(vars[ilev][Vars::zvel]);
        computeDivergence(rhs[ilev-lev_min], u, geom[ilev]);

        // If all Neumann BCs, adjust RHS to make sure we can converge
        if (need_adjust_rhs) {
            Real offset = volWgtSumMF(ilev, rhs[ilev], 0, *mapfac_m[ilev], false, false);
            amrex::Print() << "Poisson solvability offset = " << offset << std::endl;
            rhs[ilev-lev_min].plus(-offset, 0, 1);
        }
    }

    // Initialize phi to 0
    for (int ilev = lev_min; ilev <= lev_max; ++ilev)
    {
        phi[ilev-lev_min].setVal(0.0);
    }

    MLMG mlmg(mlpoisson);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    int verbose = 1;
    mlmg.setVerbose(verbose);
    //mlmg.setBottomVerbose(verbose);

    mlmg.solve(GetVecOfPtrs(phi),
               GetVecOfConstPtrs(rhs),
               solverChoice.poisson_reltol,
               solverChoice.poisson_abstol);

    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    // Subtract grad(phi) from the velocity components
    Real beta = 1.0;
    for (int ilev = lev_min; ilev <= lev_max; ++ilev) {
        MultiFab::Saxpy(vars[ilev][Vars::xvel], beta, fluxes[ilev-lev_min][0], 0,0,1,0);
        MultiFab::Saxpy(vars[ilev][Vars::yvel], beta, fluxes[ilev-lev_min][1], 0,0,1,0);
        MultiFab::Saxpy(vars[ilev][Vars::zvel], beta, fluxes[ilev-lev_min][2], 0,0,1,0);
    }

    // Average down the velocity from finest to coarsest to ensure consistency across levels
    Array<MultiFab const*, AMREX_SPACEDIM> u_fine;
    Array<MultiFab      *, AMREX_SPACEDIM> u_crse;
    for (int ilev = lev_min+1; ilev <= lev_max; ++ilev)
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
    for (int ilev = lev_min; ilev <= lev_max; ++ilev)
    {
        u[0] = &(vars[ilev][Vars::xvel]);
        u[1] = &(vars[ilev][Vars::yvel]);
        u[2] = &(vars[ilev][Vars::zvel]);
        computeDivergence(rhs[ilev], u, geom[ilev]);
        Print() << "Max norm of divergence after solve at level " << ilev << " : " << rhs[ilev].norm0() << std::endl;
    }
#endif
}
#endif
