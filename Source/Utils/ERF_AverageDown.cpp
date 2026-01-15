/**
 * \file ERF_AverageDown.cpp
 */

/**
 * Main class in ERF code, instantiated from main.cpp
*/

#include <ERF.H>
#include <ERF_Utils.H>

using namespace amrex;

// Set covered coarse cells to be the average of overlying fine cells for all levels
void
ERF::AverageDown ()
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);

    int src_comp, num_comp;
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        // If anelastic we don't average down rho because rho == rho0.
        if (solverChoice.anelastic[lev]) {
            src_comp = 1;
        } else {
            src_comp = 0;
        }
        num_comp = vars_new[0][Vars::cons].nComp() - src_comp;
        AverageDownTo(lev,src_comp,num_comp);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
ERF::AverageDownTo (int crse_lev, int scomp, int ncomp) // NOLINT
{
    if (solverChoice.anelastic[crse_lev]) {
        AMREX_ALWAYS_ASSERT(scomp == 1);
    } else {
        AMREX_ALWAYS_ASSERT(scomp == 0);
    }

    AMREX_ALWAYS_ASSERT(ncomp == vars_new[crse_lev][Vars::cons].nComp() - scomp);
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);

    // ******************************************************************************************
    // First do cell-centered quantities
    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
    // Here we multiply (rho S) by detJ and divide (rho S) by m^2 before average down
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
      for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real>   cons_arr = vars_new[lev][Vars::cons].array(mfi);
        const Array4<const Real> mfx_arr = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mfy_arr = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        if (SolverChoice::mesh_type != MeshType::ConstantDz) {
            const Array4<const Real>   detJ_arr = detJ_cc[lev]->const_array(mfi);
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) *= detJ_arr(i,j,k) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
            });
        } else {
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) /= (mfx_arr(i,j,0)*mfy_arr(i,j,0));
            });
        }
      } // mfi
    } // lev

    int fine_lev = crse_lev+1;

    if (interpolation_type == StateInterpType::Perturbational) {
        // Make the fine rho and (rho theta) be perturbational
        MultiFab::Divide(vars_new[fine_lev][Vars::cons],vars_new[fine_lev][Vars::cons],
                         Rho_comp,RhoTheta_comp,1,IntVect{0});
        MultiFab::Subtract(vars_new[fine_lev][Vars::cons],base_state[fine_lev],
                           BaseState::r0_comp,Rho_comp,1,IntVect{0});
        MultiFab::Subtract(vars_new[fine_lev][Vars::cons],base_state[fine_lev],
                           BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});

        // Make the crse rho and (rho theta) be perturbational
        MultiFab::Divide(vars_new[crse_lev][Vars::cons],vars_new[crse_lev][Vars::cons],
                         Rho_comp,RhoTheta_comp,1,IntVect{0});
        MultiFab::Subtract(vars_new[crse_lev][Vars::cons],base_state[crse_lev],
                           BaseState::r0_comp,Rho_comp,1,IntVect{0});
        MultiFab::Subtract(vars_new[crse_lev][Vars::cons],base_state[crse_lev],
                           BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
    }

    if (SolverChoice::terrain_type != TerrainType::EB) {
        average_down(vars_new[crse_lev+1][Vars::cons],vars_new[crse_lev  ][Vars::cons],
                    scomp, ncomp, refRatio(crse_lev));
    } else {
        // const auto dx = geom[fine_lev].CellSize();
        // Setting cell_vol to the exact value may cause round-off errors in volume average.
        // const Real cell_vol = dx[0]*dx[1]*dx[2];
        constexpr Real cell_vol = 1.0;
        const BoxArray& ba = vars_new[fine_lev][IntVars::cons].boxArray();
        const DistributionMapping& dm = vars_new[fine_lev][IntVars::cons].DistributionMap();
        MultiFab vol_fine(ba, dm, 1, 0);
        vol_fine.setVal(cell_vol);
        EB_average_down(vars_new[fine_lev][Vars::cons],vars_new[crse_lev][Vars::cons],
                    vol_fine, *detJ_cc[fine_lev],
                    scomp, ncomp, refRatio(crse_lev));
    }

    if (interpolation_type == StateInterpType::Perturbational) {
        // Restore the fine data to what it was
        MultiFab::Add(vars_new[fine_lev][Vars::cons],base_state[fine_lev],
                      BaseState::r0_comp,Rho_comp,1,IntVect{0});
        MultiFab::Add(vars_new[fine_lev][Vars::cons],base_state[fine_lev],
                      BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
        MultiFab::Multiply(vars_new[fine_lev][Vars::cons],vars_new[fine_lev][Vars::cons],
                           Rho_comp,RhoTheta_comp,1,IntVect{0});

        // Make the crse data be full values not perturbational
        MultiFab::Add(vars_new[crse_lev][Vars::cons],base_state[crse_lev],
                      BaseState::r0_comp,Rho_comp,1,IntVect{0});
        MultiFab::Add(vars_new[crse_lev][Vars::cons],base_state[crse_lev],
                      BaseState::th0_comp,RhoTheta_comp,1,IntVect{0});
        MultiFab::Multiply(vars_new[crse_lev][Vars::cons],vars_new[crse_lev][Vars::cons],
                           Rho_comp,RhoTheta_comp,1,IntVect{0});
    }

    vars_new[crse_lev][Vars::cons].FillBoundary(geom[crse_lev].periodicity());

    // ******************************************************************************************
    // Here we multiply (rho S) by m^2 and divide by detJ after average down
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
      for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real> cons_arr = vars_new[lev][Vars::cons].array(mfi);
        const Array4<const Real> mfx_arr = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mfy_arr = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        if (SolverChoice::mesh_type != MeshType::ConstantDz) {
            const Array4<const Real> detJ_arr = detJ_cc[lev]->const_array(mfi);
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) *= (mfx_arr(i,j,0)*mfy_arr(i,j,0)) / detJ_arr(i,j,k);
            });
        } else { // MeshType::ConstantDz
            ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                cons_arr(i,j,k,scomp+n) *= (mfx_arr(i,j,0)*mfy_arr(i,j,0));
            });
        }
      } // mfi
    } // lev

    // Fill EB covered cells by old values
    // (This won't be needed because EB_average_down copyies the covered value.)
    if (SolverChoice::terrain_type == TerrainType::EB) {
        for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real> cons_new = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real> cons_old = vars_old[lev][Vars::cons].array(mfi);
                const Array4<const Real> detJ_arr = detJ_cc[lev]->const_array(mfi);
                ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (detJ_arr(i,j,k) == 0.0) {
                        cons_new(i,j,k,scomp+n) = cons_old(i,j,k,scomp+n);
                    }
                });
            } // mfi
        } // lev
    }

    // ******************************************************************************************
    // Now average down momenta.
    // Note that vars_new holds velocities not momenta, but we want to do conservative
    //    averaging so we first convert to momentum, then average down, then convert
    //    back to velocities -- only on the valid region
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1; lev++)
    {
        // FillBoundary for density so we can go back and forth between velocity and momentum
        vars_new[lev][Vars::cons].FillBoundary(geom[lev].periodicity());

        const MultiFab* c_vfrac = nullptr;
        if (SolverChoice::terrain_type == TerrainType::EB) {
            c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
        }

        VelocityToMomentum(vars_new[lev][Vars::xvel], IntVect(0,0,0),
                        vars_new[lev][Vars::yvel], IntVect(0,0,0),
                        vars_new[lev][Vars::zvel], IntVect(0,0,0),
                        vars_new[lev][Vars::cons],
                            rU_new[lev],
                            rV_new[lev],
                            rW_new[lev],
                        Geom(lev).Domain(),
                        domain_bcs_type,
                        c_vfrac);
    }

    if (SolverChoice::terrain_type != TerrainType::EB) {
        average_down_faces(rU_new[crse_lev+1], rU_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
        average_down_faces(rV_new[crse_lev+1], rV_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
        average_down_faces(rW_new[crse_lev+1], rW_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
    } else {
        EB_average_down_faces({&rU_new[crse_lev+1], &rV_new[crse_lev+1], &rW_new[crse_lev+1]},
                            {&rU_new[crse_lev], &rV_new[crse_lev], &rW_new[crse_lev]},
                            refRatio(crse_lev), 0);
    }

    for (int lev = crse_lev; lev <= crse_lev+1; lev++) {

        const MultiFab* c_vfrac = nullptr;
        if (SolverChoice::terrain_type == TerrainType::EB) {
            c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
        }

        MomentumToVelocity(vars_new[lev][Vars::xvel],
                        vars_new[lev][Vars::yvel],
                        vars_new[lev][Vars::zvel],
                        vars_new[lev][Vars::cons],
                            rU_new[lev],
                            rV_new[lev],
                            rW_new[lev],
                        Geom(lev).Domain(),
                        domain_bcs_type,
                        c_vfrac);
    }
}
