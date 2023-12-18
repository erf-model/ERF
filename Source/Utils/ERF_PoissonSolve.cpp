#include "ERF.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

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
    for (int dir = 0; dir < 2; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = domain_bc_type[Orientation(dir,side)];
            if (bc_type == "outflow") {
                r[dir] = LinOpBCType::Dirichlet;
            } else {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    r[2] = LinOpBCType::Neumann;
    amrex::Print() << "BCs for Poisson solve " << r[0] << " " << r[1] << " " << r[2] << std::endl;
    return r;
}


/**
 * Project the single-level velocity field to enforce incompressibility
 */
void ERF::project_velocities (Vector<MultiFab>& vmf)
{
    Vector<Vector<MultiFab>> tmpmf(1);
    for (auto& mf : vmf) {
        tmpmf[0].emplace_back(mf, amrex::make_alias, 0, mf.nComp());
    }
    project_velocities(tmpmf);
}

/**
 * Project the multi-level velocity field to enforce incompressibility
 */
void
ERF::project_velocities (Vector<Vector<MultiFab>>& vars)
{
    BL_PROFILE("ERF::project_velocities()");

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 1.e-10;

    const int nlevs = geom.size();

    // Use the default settings
    LPInfo info;
    MLPoisson mlpoisson(geom, grids, dmap, info);

    // This is a 3d problem with Dirichlet BC
    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    mlpoisson.setDomainBC(bclo, bchi);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
       mlpoisson.setLevelBC(0, nullptr);
    }

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

    rhs.resize(nlevs);
    phi.resize(nlevs);
    fluxes.resize(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev].setVal(0.0);
        phi[ilev].setVal(0.0);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            fluxes[ilev][idim].define(
                convert(grids[ilev], IntVect::TheDimensionVector(idim)),
                dmap[ilev], 1, 0);
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
    }

    // Initialize phi to 0
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        phi[ilev].setVal(0.0);
    }

    MLMG mlmg(mlpoisson);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    int verbose = 1;
    mlmg.setVerbose(verbose);
    // mlmg.setBottomVerbose(bottom_verbose);

    mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    // Subtract grad(phi) from the velocity components
    Real beta = 1.0;
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        MultiFab::Saxpy(vars[ilev][Vars::xvel], beta, fluxes[ilev][0], 0,0,1,0);
        MultiFab::Saxpy(vars[ilev][Vars::yvel], beta, fluxes[ilev][1], 0,0,1,0);
        MultiFab::Saxpy(vars[ilev][Vars::zvel], beta, fluxes[ilev][2], 0,0,1,0);
    }

    // Average down the velocity from finest to coarsest to ensure consistency across levels
    int finest_level = nlevs - 1;
    Array<MultiFab const*, AMREX_SPACEDIM> u_fine;
    Array<MultiFab      *, AMREX_SPACEDIM> u_crse;
    for (int ilev = finest_level; ilev > 0; --ilev)
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
        computeDivergence(rhs[ilev], u, geom[ilev]);
        Print() << "Max norm of divergence after solve at level " << ilev << " : " << rhs[ilev].norm0() << std::endl;
    }
#endif
}
#endif
