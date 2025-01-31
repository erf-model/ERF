#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Project the single-level velocity field to enforce incompressibility
 * Note that the level may or may not be level 0.
 */
void ERF::project_velocities (int lev, Real l_dt, Vector<MultiFab>& mom_mf, MultiFab& pmf)
{
    BL_PROFILE("ERF::project_velocities()");

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(mom_mf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(mom_mf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;

    if (solverChoice.terrain_type == TerrainType::EB)
    {
        rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0, MFInfo(), Factory(lev));
        phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1, MFInfo(), Factory(lev));
    } else {
        rhs.resize(1); rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
        phi.resize(1); phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);
    }

    auto dxInv = geom[lev].InvCellSizeArray();

    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab to hold Omega rather than rhow
    // ****************************************************************************
    //
    if (solverChoice.mesh_type == MeshType::VariableDz)
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
    }

    // ****************************************************************************
    // Compute divergence which will form RHS
    // Note that we replace "rho0w" with the contravariant momentum, Omega
    // ****************************************************************************
    Array<MultiFab const*, AMREX_SPACEDIM> rho0_u_const;
    rho0_u_const[0] = &mom_mf[IntVars::xmom];
    rho0_u_const[1] = &mom_mf[IntVars::ymom];
    rho0_u_const[2] = &mom_mf[IntVars::zmom];

    compute_divergence(lev, rhs[0], rho0_u_const, geom_tmp[0]);

    Real rhsnorm = rhs[0].norm0();

    if (mg_verbose > 0) {
        Print() << "Max/L2 norm of divergence before solve at level " << lev << " : " << rhsnorm << " " << rhs[0].norm2() << std::endl;
    }

    // ****************************************************************************
    //
    // No need to build the solver if RHS == 0
    //
    if (rhsnorm <= solverChoice.poisson_abstol) return;
    // ****************************************************************************

    // ****************************************************************************
    // Initialize phi to 0
    // (It is essential that we do this in order to fill the corners; these are never
    //  used but the Saxpy requires the values to be initialized.)
    // ****************************************************************************
    phi[0].setVal(0.0);

    Real start_step = static_cast<Real>(ParallelDescriptor::second());

    // ****************************************************************************
    // Allocate fluxes
    // ****************************************************************************
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;
    fluxes.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (solverChoice.terrain_type == TerrainType::EB) {
            fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0, MFInfo(), Factory(lev));
        } else {
            fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
        }
    }

    // ****************************************************************************
    // Choose the solver and solve
    // ****************************************************************************

    // ****************************************************************************
    // EB
    // ****************************************************************************
    if (solverChoice.terrain_type == TerrainType::EB) {
        solve_with_EB_mlmg(lev, rhs, phi, fluxes);
    } else {

#ifdef ERF_USE_FFT
        Box my_region(ba_tmp[0].minimalBox());
        bool boxes_make_rectangle = (my_region.numPts() == ba_tmp[0].numPts());
#endif

    // ****************************************************************************
    // No terrain or grid stretching
    // ****************************************************************************
    if (solverChoice.mesh_type == MeshType::ConstantDz) {
#ifdef ERF_USE_FFT
        if (use_fft) {
            if (boxes_make_rectangle) {
                solve_with_fft(lev, rhs[0], phi[0], fluxes[0]);
            } else {
                amrex::Warning("FFT won't work unless the boxArray covers the domain: defaulting to MLMG");
                solve_with_mlmg(lev, rhs, phi, fluxes);
            }
        } else {
            solve_with_mlmg(lev, rhs, phi, fluxes);
        }
#else
        if (use_fft) {
            amrex::Warning("You set use_fft=true but didn't build with USE_FFT = TRUE; defaulting to MLMG");
        }
        solve_with_mlmg(lev, rhs, phi, fluxes);
#endif
    } // No terrain or grid stretching

    // ****************************************************************************
    // Grid stretching (flat terrain)
    // ****************************************************************************
    else if (solverChoice.mesh_type == MeshType::StretchedDz) {
#ifndef ERF_USE_FFT
        amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT solver");
#else
        if (!boxes_make_rectangle) {
            amrex::Abort("FFT won't work unless the boxArray covers the domain");
        } else {
            if (!use_fft) {
                amrex::Warning("Using FFT even though you didn't set use_fft to true; it's the best choice");
            }
            solve_with_fft(lev, rhs[0], phi[0], fluxes[0]);
        }
#endif
    } // grid stretching

    // ****************************************************************************
    // General terrain
    // ****************************************************************************
    else if (solverChoice.mesh_type == MeshType::VariableDz) {
#ifdef ERF_USE_FFT
        if (use_fft)
        {
            amrex::Warning("FFT solver does not work for general terrain: switching to FFT-preconditioned GMRES");
        }
        if (!boxes_make_rectangle) {
            amrex::Abort("FFT preconditioner for GMRES won't work unless the boxArray covers the domain");
        } else {
            solve_with_gmres(lev, rhs, phi, fluxes);
        }
#else
        amrex::Abort("Rebuild with USE_FFT = TRUE so you can use the FFT preconditioner for GMRES");
#endif
    } // general terrain

    } // not EB

    // ****************************************************************************
    // Print time in solve
    // ****************************************************************************
    Real end_step = static_cast<Real>(ParallelDescriptor::second());
    if (mg_verbose > 0) {
        amrex::Print() << "Time in solve " << end_step - start_step << std::endl;
    }

    // ****************************************************************************
    // Subtract dt grad(phi) from the momenta (rho0u, rho0v, Omega)
    // ****************************************************************************
    MultiFab::Add(mom_mf[IntVars::xmom],fluxes[0][0],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::ymom],fluxes[0][1],0,0,1,0);
    MultiFab::Add(mom_mf[IntVars::zmom],fluxes[0][2],0,0,1,0);

    //
    // This call is only to verify the divergence after the solve
    // It is important we do this before computing the rho0w_arr from Omega back to rho0w
    //
    // ****************************************************************************
    // THIS IS SIMPLY VERIFYING THE DIVERGENCE AFTER THE SOLVE
    // ****************************************************************************
    //
    if (mg_verbose > 0)
    {
        compute_divergence(lev, rhs[0], rho0_u_const, geom_tmp[0]);

        Print() << "Max/L2 norm of divergence  after solve at level " << lev << " : " << rhs[0].norm0() << " " << rhs[0].norm2() << std::endl;

#if 0
         // FOR DEBUGGING ONLY
         for ( MFIter mfi(rhs[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Array4<Real const>& rhs_arr = rhs[0].const_array(mfi);
            Box bx = mfi.validbox();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (std::abs(rhs_arr(i,j,k)) > 1.e-10) {
                    amrex::Print() << "RHS AT " << IntVect(i,j,k) << " " << rhs_arr(i,j,k) << std::endl;
                }
            });
         } // mfi
#endif

    } // mg_verbose


    //
    // ****************************************************************************
    // Now convert the rho0w MultiFab back to holding (rho0w) rather than Omega
    // ****************************************************************************
    //
    if (solverChoice.mesh_type == MeshType::VariableDz)
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
    MultiFab::Saxpy(pmf, 1.0/l_dt, phi[0],0,0,1,1);
    pmf.FillBoundary(geom[lev].periodicity());
}
