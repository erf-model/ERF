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
    BL_PROFILE("ERF::project_velocities()");
    AMREX_ALWAYS_ASSERT(!solverChoice.use_terrain);

    // Make sure the solver only sees the levels over which we are solving
    LPInfo info;
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(vmf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(vmf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);
    // amrex::Print() << "AT LEVEL " << lev << " BA FOR POISSON SOLVE " << vmf[Vars::cons].boxArray() << std::endl;

    MLABecLaplacian mlabec(geom_tmp, ba_tmp, dm_tmp, info);

    //
    // This will hold (1/rho) on faces
    //
    Array<MultiFab,AMREX_SPACEDIM> inv_rho;

    //
    // The operator is (alpha A - beta del dot B grad) phi = RHS
    // Here we set alpha to 0 and beta to -1
    // Then b is (dt/rho)
    //
    mlabec.setScalars(0.0, -1.0);
    inv_rho[0].define(vmf[Vars::xvel].boxArray(),dm_tmp[0],1,0,MFInfo());
    inv_rho[1].define(vmf[Vars::yvel].boxArray(),dm_tmp[0],1,0,MFInfo());
    inv_rho[2].define(vmf[Vars::zvel].boxArray(),dm_tmp[0],1,0,MFInfo());

#if 0
    MultiFab density(vmf[Vars::cons], make_alias, 1, 1);
    density.FillBoundary(geom_tmp[0].periodicity());
    amrex::average_cellcenter_to_face(GetArrOfPtrs(inv_rho), density, geom_tmp[0]);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        inv_rho[idim].invert(l_dt, 0);
    }
#else
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        inv_rho[idim].setVal(l_dt);
    }
#endif

    mlabec.setBCoeffs(0, GetArrOfConstPtrs(inv_rho));

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    bool need_adjust_rhs = (projection_has_dirichlet(bclo) || projection_has_dirichlet(bchi)) ? false : true;
    mlabec.setDomainBC(bclo, bchi);

    if (lev > 0) {
        mlabec.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }
    mlabec.setLevelBC(0, nullptr);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

    rhs.resize(1);
    phi.resize(1);
    fluxes.resize(1);

    rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    phi[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    rhs[0].setVal(0.0);
    phi[0].setVal(0.0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
    }

    // Define a single Array of MultiFabs to hold the velocity components
    // Then define the RHS to be the divergence of the current velocities
    Array<MultiFab const*, AMREX_SPACEDIM> u;
    u[0] = &(vmf[Vars::xvel]);
    u[1] = &(vmf[Vars::yvel]);
    u[2] = &(vmf[Vars::zvel]);
    computeDivergence(rhs[0], u, geom_tmp[0]);

    // If all Neumann BCs, adjust RHS to make sure we can converge
    if (need_adjust_rhs)
    {
        Real offset = volWgtSumMF(lev, rhs[0], 0, *mapfac_m[lev], false, false);
        // amrex::Print() << "Poisson solvability offset = " << offset << std::endl;
        rhs[0].plus(-offset, 0, 1);
    }

    // Initialize phi to 0
    phi[0].setVal(0.0);

    MLMG mlmg(mlabec);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve(GetVecOfPtrs(phi),
               GetVecOfConstPtrs(rhs),
               solverChoice.poisson_reltol,
               solverChoice.poisson_abstol);

    mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

    // Update pressure variable with phi -- note that phi is change in pressure, not the full pressure
    MultiFab::Saxpy(pmf, 1.0, phi[0],0,0,1,0);
    pmf.FillBoundary(geom[lev].periodicity());

    // Subtract (dt/rho) grad(phi) from the velocity components
    MultiFab::Add(vmf[Vars::xvel], fluxes[0][0], 0,0,1,0);
    MultiFab::Add(vmf[Vars::yvel], fluxes[0][1], 0,0,1,0);
    MultiFab::Add(vmf[Vars::zvel], fluxes[0][2], 0,0,1,0);

#if 0
    // Confirm that the velocity is now divergence free
    u[0] = &(vmf[Vars::xvel]);
    u[1] = &(vmf[Vars::yvel]);
    u[2] = &(vmf[Vars::zvel]);
    rhs[0].setVal(0.0);
    computeDivergence(rhs[0], u, geom_tmp[0]);
    Print() << "Max norm of divergence after solve at level " << lev << " : " << rhs[0].norm0() << std::endl;
#endif
}
#endif
