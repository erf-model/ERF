#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
//#include <AMReX_MLTerrainPoisson.H>
#include <AMReX_GMRES.H>
#include <AMReX_GMRES_MLMG.H>

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
    return r;
}

#ifdef ERF_USE_FFT
Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
ERF::get_fft_bc () const noexcept
{
    AMREX_ALWAYS_ASSERT(geom[0].isPeriodic(0) && geom[0].isPeriodic(1));

    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;

    for (int dir = 0; dir <= 1; dir++) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = std::make_pair(FFT::Boundary::periodic,FFT::Boundary::periodic);
        }
    } // dir

    for (OrientationIter ori; ori != nullptr; ++ori) {
        const int dir  = ori().coordDir();
        if (!geom[0].isPeriodic(dir) && ori().faceDir() == Orientation::low) {
            auto bc_type_lo = domain_bc_type[Orientation(dir,Orientation::low)];
            auto bc_type_hi = domain_bc_type[Orientation(dir,Orientation::high)];
            if (bc_type_lo == "Outflow" && bc_type_hi == "Outflow") {
                r[dir] = std::make_pair(FFT::Boundary::odd,FFT::Boundary::odd);
            } else if (bc_type_lo != "Outflow" && bc_type_hi == "Outflow") {
                r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::odd);
            } else if (bc_type_lo == "Outflow" && bc_type_hi != "Outflow") {
                r[dir] = std::make_pair(FFT::Boundary::odd,FFT::Boundary::even);
            } else {
                r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::even);
            }
        } // not periodic
    } // ori

    return r;
}
#endif

/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::project_velocities (int lev, Real l_dt, Vector<MultiFab>& mom_mf, MultiFab& pmf)
{
    BL_PROFILE("ERF::project_velocities()");

    bool l_use_terrain = SolverChoice::terrain_type != TerrainType::None;
    bool use_gmres     = (l_use_terrain && !SolverChoice::terrain_is_flat);

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    LPInfo info;
    // Allow a hidden direction if the domain is one cell wide in any lateral direction
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
    }

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);

    // amrex::Print() << "BCLO " << bclo[0] << " " << bclo[1] << " " << bclo[2] << std::endl;
    // amrex::Print() << "BCHI " << bchi[0] << " " << bchi[1] << " " << bchi[2] << std::endl;

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;

    rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);

    fluxes.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
    }

    auto dxInv = geom[lev].InvCellSizeArray();

    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &mom_mf[IntVars::xmom];
    rho0_u_const[1] = &mom_mf[IntVars::ymom];
    rho0_u_const[2] = &mom_mf[IntVars::zmom];

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    // ****************************************************************************
    // Compute divergence which will form RHS
    // Note that we replace "rho0w" with the contravariant momentum, Omega
    // ****************************************************************************
    if (l_use_terrain)
    {
        for ( MFIter mfi(rhs[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Array4<Real const>& rho0u_arr = mom_mf[IntVars::xmom].const_array(mfi);
            const Array4<Real const>& rho0v_arr = mom_mf[IntVars::ymom].const_array(mfi);
            const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);

            const Array4<Real const>&     z_nd = z_phys_nd[lev]->const_array(mfi);
            //
            // Define Omega from (rho0 W) but store it in the same array
            //
            Box tbz = mfi.nodaltilebox(2);
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (k > dom_lo.z && k <= dom_hi.z) {
                    Real rho0w = rho0w_arr(i,j,k);
                    rho0w_arr(i,j,k) = OmegaFromW(i,j,k,rho0w,rho0u_arr,rho0v_arr,z_nd,dxInv);
                } else {
                    rho0w_arr(i,j,k) = Real(0.0);
                }
            });
        } // mfi
    } // use_terrain

    computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);

    if (l_use_terrain)
    {
        MultiFab::Divide(rhs[0],*detJ_cc[0],0,0,1,0);
    }

    Real rhsnorm = rhs[0].norm0();

    if (mg_verbose > 0) {
        Print() << "Max norm of divergence before at level " << lev << " : " << rhsnorm << std::endl;
    }

    if (rhsnorm <= abstol) return;

    Real start_step = static_cast<Real>(ParallelDescriptor::second());

#ifdef ERF_USE_FFT
    // ****************************************************************************
    // FFT solve
    // ****************************************************************************
    if (use_fft)
    {
        AMREX_ALWAYS_ASSERT(lev == 0);
        //
        // No terrain or stretched grids
        // This calls the full 3D FFT solver with bc's set through bc_fft
        //
        if (!l_use_terrain)
        {
            if (mg_verbose > 0) {
                amrex::Print() << "Using the 3D FFT solver..." << std::endl;
            }
            if (!m_3D_poisson) {
                auto bc_fft = get_fft_bc();
                m_3D_poisson = std::make_unique<FFT::Poisson<MultiFab>>(Geom(0),bc_fft);
            }
            m_3D_poisson->solve(phi[lev], rhs[lev]);

        //
        // Stretched grids
        // This calls the hybrid 2D FFT solver + tridiagonal in z
        //
        // For right now we can only do this solve for periodic in the x- and y-directions
        // We assume Neumann at top and bottom z-boundaries
        // This will be generalized in future
        //
        //
        } else if (l_use_terrain && SolverChoice::terrain_is_flat)
        {
            if (mg_verbose > 0) { amrex::Print() << "Using the hybrid FFT solver..." << std::endl;
            }
            if (!m_2D_poisson) {
                m_2D_poisson = std::make_unique<FFT::PoissonHybrid<MultiFab>>(Geom(0));
            }
            Gpu::DeviceVector<Real> stretched_dz(dom_hi.z+1, geom[lev].CellSize(2));
            m_2D_poisson->solve(phi[lev], rhs[lev], stretched_dz);

        } else {
            amrex::Abort("FFT isn't appropriate for spatially varying terrain");
        }

        phi[lev].FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real const> const&  p_arr  = phi[lev].array(mfi);

            Box const& xbx = mfi.nodaltilebox(0);
            const Real dx_inv = dxInv[0];
            Array4<Real> const& fx_arr  = fluxes[lev][0].array(mfi);
            ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                fx_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i-1,j,k)) * dx_inv;
            });

            Box const& ybx = mfi.nodaltilebox(1);
            const Real dy_inv = dxInv[1];
            Array4<Real> const& fy_arr  = fluxes[lev][1].array(mfi);
            ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                fy_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dy_inv;
            });

            Box const& zbx = mfi.nodaltilebox(2);
            const Real dz_inv = dxInv[2];
            Array4<Real> const& fz_arr  = fluxes[lev][2].array(mfi);
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == dom_lo.z || k == dom_hi.z+1) {
                    fz_arr(i,j,k) = 0.0;
                } else {
                    fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) * dz_inv;
                }
            });
        } // mfi
     } else
#endif
    {
        // Initialize phi to 0
        phi[0].setVal(0.0);

        // ****************************************************************************
        // GMRES solve
        // ****************************************************************************
        if (use_gmres)
        {
            amrex::Abort("GMRES isn't ready yet");
#if 0
            MLTerrainPoisson terrpoisson(geom_tmp, ba_tmp, dm_tmp, info);
            terrpoisson.setDomainBC(bclo, bchi);
            terrpoisson.setMaxOrder(2);

            terrpoisson.setZPhysNd(lev, *z_phys_nd[lev]);

            if (lev > 0) {
                terrpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
            }
            terrpoisson.setLevelBC(lev, &phi[lev]);

            MLMG mlmg(terrpoisson);
            GMRESMLMG gmsolver(mlmg);
            gmsolver.usePrecond(false);
            gmsolver.setVerbose(mg_verbose);
            gmsolver.solve(phi[0], rhs[0], reltol, abstol);

            Vector<MultiFab*> phi_vec; phi_vec.resize(1);
            phi_vec[0] = &phi[0];
            terrpoisson.getFluxes(GetVecOfArrOfPtrs(fluxes), phi_vec);

#if 0
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Array4<Real const> const&  p_arr  = phi[lev].array(mfi);

                Box const& xbx = mfi.nodaltilebox(0);
                const Real dx_inv = dxInv[0];
                Array4<Real> const& fx_arr  = fluxes[lev][0].array(mfi);
                ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fx_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i-1,j,k)) * dx_inv;
                });

                Box const& ybx = mfi.nodaltilebox(1);
                const Real dy_inv = dxInv[1];
                Array4<Real> const& fy_arr  = fluxes[lev][1].array(mfi);
                ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fy_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dy_inv;
                });

                auto const dom_lo = lbound(geom[lev].Domain());
                auto const dom_hi = ubound(geom[lev].Domain());

                Box const& zbx = mfi.nodaltilebox(2);
                const Real dz_inv = dxInv[2];
                Array4<Real> const& fz_arr  = fluxes[lev][2].array(mfi);
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k == dom_lo.z || k == dom_hi.z+1) {
                        fz_arr(i,j,k) = 0.0;
                    } else {
                        fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) * dz_inv;
                    }
                });
            } // mfi
#endif
#endif

        // ****************************************************************************
        // Multigrid solve
        // ****************************************************************************
        } else { // use MLMG

            MLPoisson mlpoisson(geom_tmp, ba_tmp, dm_tmp, info);
            mlpoisson.setDomainBC(bclo, bchi);
            if (lev > 0) {
                mlpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
            }
            mlpoisson.setLevelBC(0, nullptr);

            MLMG mlmg(mlpoisson);
            int max_iter = 100;
            mlmg.setMaxIter(max_iter);

            mlmg.setVerbose(mg_verbose);
            mlmg.setBottomVerbose(0);

            mlmg.solve(GetVecOfPtrs(phi),
                       GetVecOfConstPtrs(rhs),
                       reltol, abstol);
            mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));
        }
    } // not using fft

    // ****************************************************************************
    // Subtract dt grad(phi) from the momenta (rho0u, rho0v, Omega)
    // ****************************************************************************
    MultiFab::Add(mom_mf[IntVars::xmom],fluxes[0][0],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::ymom],fluxes[0][1],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::zmom],fluxes[0][2],0,0,1,0);

    // ****************************************************************************
    // Print time in solve
    // ****************************************************************************
    Real end_step = static_cast<Real>(ParallelDescriptor::second());
    if (mg_verbose > 0) {
        amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
    }

    //
    // This call is only to verify the divergence after the solve
    // It is important we do this before computing the rho0w_arr from Omega back to rho0w
    //
    //
    // ****************************************************************************
    // THIS IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    // ****************************************************************************
    //
    if (mg_verbose > 0)
    {
        computeDivergence(rhs[0], rho0_u_const, geom_tmp[0]);
        if (l_use_terrain)
        {
            MultiFab::Divide(rhs[0],*detJ_cc[0],0,0,1,0);
        } // use_terrain

        amrex::Print() << "Max norm of divergence after solve at level " << lev << " : " << rhs[0].norm0() << std::endl;
    }

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab back to holding (rho0w) rather than Omega
    // ****************************************************************************
    //
    if (l_use_terrain)
    {
        for (MFIter mfi(mom_mf[Vars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
             Box tbz = mfi.nodaltilebox(2);
             const Array4<Real      >& rho0u_arr = mom_mf[IntVars::xmom].array(mfi);
             const Array4<Real      >& rho0v_arr = mom_mf[IntVars::ymom].array(mfi);
             const Array4<Real      >& rho0w_arr = mom_mf[IntVars::zmom].array(mfi);
             const Array4<Real const>&      z_nd = z_phys_nd[lev]->const_array(mfi);
             ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                 Real omega = rho0w_arr(i,j,k);
                 rho0w_arr(i,j,k) = WFromOmega(i,j,k,omega,rho0u_arr,rho0v_arr,z_nd,dxInv);
             });
        } // mfi
    }

    // ****************************************************************************
    // Update pressure variable with phi -- note that phi is dt * change in pressure
    // ****************************************************************************
    MultiFab::Saxpy(pmf, 1.0/l_dt, phi[0],0,0,1,0);
    pmf.FillBoundary(geom[lev].periodicity());

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mom_mf[Vars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& pp_arr  = pmf.array(mfi);
        Box const& bx    = mfi.tilebox();
        auto const bx_lo = lbound(bx);
        auto const bx_hi = ubound(bx);
        if (bx_lo.x == dom_lo.x) {
            if (bclo[0] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i-1,j,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,0,dom_lo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i-1,j,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_lo.y == dom_lo.y) {
            if (bclo[1] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,1,dom_lo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j-1,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_lo.z == dom_lo.z) {
            ParallelFor(makeSlab(bx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k-1) = pp_arr(i,j,k);
            });
        }
        if (bx_hi.x == dom_hi.x) {
            if (bchi[0] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i+1,j,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,0,dom_hi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i+1,j,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_hi.y == dom_hi.y) {
            if (bchi[1] == LinOpBCType::Dirichlet) {
                ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = -pp_arr(i,j,k);
                });
            } else {
                ParallelFor(makeSlab(bx,1,dom_hi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pp_arr(i,j+1,k) = pp_arr(i,j,k);
                });
            }
        }
        if (bx_hi.z == dom_hi.z) {
            ParallelFor(makeSlab(bx,2,dom_hi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k+1) = pp_arr(i,j,k);
            });
        }
    } // mfi

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    pmf.FillBoundary(geom[lev].periodicity());
}
