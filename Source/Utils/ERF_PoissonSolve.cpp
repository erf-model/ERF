#include "ERF.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include "ERF_Utils.H"

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
            if (bc_type == "Outflow") {
                r[dir] = LinOpBCType::Dirichlet;
            } else
            {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    // r[2] = LinOpBCType::Neumann;
    // amrex::Print() << "BCs for Poisson solve " << r[0] << " " << r[1] << " " << r[2] << std::endl;
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

    MultiFab density(vmf[Vars::cons], make_alias, Rho_comp, 1);
    density.FillBoundary(geom_tmp[0].periodicity());

    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const& rho_arr     = density.const_array(mfi);
        Array4<Real const> const& rho_0_arr   = r_hse.const_array(mfi);

        Box const& bxx = mfi.nodaltilebox(0);
        Array4<Real      > const& inv_rhox_arr = inv_rho[0].array(mfi);
        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho_edge = Real(0.5) * (rho_arr(i,j,k) + rho_arr(i-1,j,k));
            inv_rhox_arr(i,j,k) = l_dt * rho_0_arr(i,j,k) / rho_edge;
        });

        Box const& bxy = mfi.nodaltilebox(1);
        Array4<Real      > const& inv_rhoy_arr = inv_rho[1].array(mfi);
        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho_edge = Real(0.5) * (rho_arr(i,j,k) + rho_arr(i,j-1,k));
            inv_rhoy_arr(i,j,k) = l_dt * rho_0_arr(i,j,k) / rho_edge;
        });

        Box const& bxz = mfi.nodaltilebox(2);
        Array4<Real      > const& inv_rhoz_arr = inv_rho[2].array(mfi);
        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho_edge = Real(0.5) * (rho_arr(i,j,k) + rho_arr(i,j,k-1));
            Real rho_0_edge = Real(0.5) * (rho_0_arr(i,j,k) + rho_0_arr(i,j,k-1));
            inv_rhoz_arr(i,j,k) = l_dt * rho_0_edge / rho_edge;
        });
    } // mfi

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

    Array<MultiFab,AMREX_SPACEDIM> rho0_u;
    rho0_u[0].define(vmf[Vars::xvel].boxArray(),dm_tmp[0],1,0,MFInfo());
    rho0_u[1].define(vmf[Vars::yvel].boxArray(),dm_tmp[0],1,0,MFInfo());
    rho0_u[2].define(vmf[Vars::zvel].boxArray(),dm_tmp[0],1,0,MFInfo());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const& rho0_arr  = r_hse.const_array(mfi);

        Box const& bxx = mfi.nodaltilebox(0);
        Array4<Real const> const&      u_arr = vmf[Vars::xvel].const_array(mfi);
        Array4<Real      > const& rho0_u_arr = rho0_u[0].array(mfi);
        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho0_u_arr(i,j,k) = u_arr(i,j,k) * rho0_arr(i,j,k);
        });

        Box const& bxy = mfi.nodaltilebox(1);
        Array4<Real const> const&      v_arr = vmf[Vars::yvel].const_array(mfi);
        Array4<Real      > const& rho0_v_arr = rho0_u[1].array(mfi);
        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho0_v_arr(i,j,k) = v_arr(i,j,k) * rho0_arr(i,j,k);
        });

        Box const& bxz = mfi.nodaltilebox(2);
        Array4<Real const> const&      w_arr = vmf[Vars::zvel].const_array(mfi);
        Array4<Real      > const& rho0_w_arr = rho0_u[2].array(mfi);
        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho0_edge = Real(0.5) * (rho0_arr(i,j,k) + rho0_arr(i,j,k-1));
            rho0_w_arr(i,j,k) = w_arr(i,j,k) * rho0_edge;
        });
    } // mfi

    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &rho0_u[0];
    rho0_u_const[1] = &rho0_u[1];
    rho0_u_const[2] = &rho0_u[2];

    computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);
    Print() << "Max norm of divergence after  at level " << lev << " : " << rhs[0].norm0() << std::endl;

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

    // Subtract (dt rho0/rho) grad(phi) from the rho0-weighted velocity components
    // MultiFab::Add(vmf[Vars::xvel], fluxes[0][0], 0,0,1,0);
    // MultiFab::Add(vmf[Vars::yvel], fluxes[0][1], 0,0,1,0);
    // MultiFab::Add(vmf[Vars::zvel], fluxes[0][2], 0,0,1,0);
    MultiFab::Add(rho0_u[0], fluxes[0][0], 0,0,1,0);
    MultiFab::Add(rho0_u[1], fluxes[0][1], 0,0,1,0);
    MultiFab::Add(rho0_u[2], fluxes[0][2], 0,0,1,0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const& rho0_arr  = r_hse.const_array(mfi);

        Box const& bxx = mfi.nodaltilebox(0);
        Array4<Real      > const&      u_arr = vmf[Vars::xvel].array(mfi);
        Array4<Real const> const& rho0_u_arr = rho0_u[0].array(mfi);
        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u_arr(i,j,k) = rho0_u_arr(i,j,k) / rho0_arr(i,j,k);
        });

        Box const& bxy = mfi.nodaltilebox(1);
        Array4<Real      > const&      v_arr = vmf[Vars::yvel].array(mfi);
        Array4<Real const> const& rho0_v_arr = rho0_u[1].array(mfi);
        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            v_arr(i,j,k) = rho0_v_arr(i,j,k) / rho0_arr(i,j,k);
        });

        Box const& bxz = mfi.nodaltilebox(2);
        Array4<Real      > const&      w_arr = vmf[Vars::zvel].array(mfi);
        Array4<Real const> const& rho0_w_arr = rho0_u[2].array(mfi);
        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho0_edge = Real(0.5) * (rho0_arr(i,j,k) + rho0_arr(i,j,k-1));
            w_arr(i,j,k) = rho0_w_arr(i,j,k) / rho0_edge;
        });
    } // mfi

#if 0
    //
    // BELOW IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const& rho0_arr  = r_hse.const_array(mfi);

        Box const& bxx = mfi.nodaltilebox(0);
        Array4<Real const> const&      u_arr = vmf[Vars::xvel].const_array(mfi);
        Array4<Real      > const& rho0_u_arr = rho0_u[0].array(mfi);
        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho0_u_arr(i,j,k) = u_arr(i,j,k) * rho0_arr(i,j,k);
        });

        Box const& bxy = mfi.nodaltilebox(1);
        Array4<Real const> const&      v_arr = vmf[Vars::yvel].const_array(mfi);
        Array4<Real      > const& rho0_v_arr = rho0_u[1].array(mfi);
        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho0_v_arr(i,j,k) = v_arr(i,j,k) * rho0_arr(i,j,k);
        });

        Box const& bxz = mfi.nodaltilebox(2);
        Array4<Real const> const&      w_arr = vmf[Vars::zvel].const_array(mfi);
        Array4<Real      > const& rho0_w_arr = rho0_u[2].array(mfi);
        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho0_edge = Real(0.5) * (rho0_arr(i,j,k) + rho0_arr(i,j,k-1));
            rho0_w_arr(i,j,k) = w_arr(i,j,k) * rho0_edge;
        });
    } // mfi

    computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);
    Print() << "Max norm of divergence after solve at level " << lev << " : " << rhs[0].norm0() << std::endl;
#endif
}
#endif
