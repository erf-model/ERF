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
        // If anelastic we don't average down rho because rho == rho0 (a pseudo-incompressible
        // level evolves rho, which is averaged down as on a compressible level)
        src_comp = (solverChoice.compressibility_type[lev] == CompressibilityType::Anelastic) ? 1 : 0;
        num_comp = vars_new[0][Vars::cons].nComp() - src_comp;
        AverageDownTo(lev,src_comp,num_comp);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
ERF::AverageDownTo (int crse_lev, int scomp, int ncomp, bool do_full_state_and_momenta) // NOLINT
{
    if (do_full_state_and_momenta) {
        if (solverChoice.compressibility_type[crse_lev] == CompressibilityType::Anelastic) {
            AMREX_ALWAYS_ASSERT(scomp == 1);
        } else {
            AMREX_ALWAYS_ASSERT(scomp == 0);
        }
        AMREX_ALWAYS_ASSERT(ncomp == vars_new[crse_lev][Vars::cons].nComp() - scomp);
    } else {
        AMREX_ALWAYS_ASSERT(scomp >= 0);
        AMREX_ALWAYS_ASSERT(ncomp > 0);
        AMREX_ALWAYS_ASSERT(scomp + ncomp <= vars_new[crse_lev][Vars::cons].nComp());
    }
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);

    // Immersed forcing with erf.if_average_down = fluid_mass: each level forces its own solid
    // and partial cells, so the fine cells are weighted by their fluid mass and coarse cells that
    // are solid on their own level keep their value (the plain average mixes the values of solid
    // fine cells into coarse fluid cells). The branch carries the cell volume detJ/(m_x m_y) in
    // its weight itself (see if_average_down), so the scaling of the state around the plain
    // average is skipped, and it takes the coarse density as it is when the density is not
    // averaged, so the rescaling of the plain anelastic branch is built in.
    const bool if_fluid_mass = (solverChoice.if_average_down == "fluid_mass") &&
                               (solverChoice.terrain_type   == TerrainType::ImmersedForcing ||
                                solverChoice.buildings_type == BuildingsType::ImmersedForcing);

    // erf.perturbational_average_down: the cell-centered state is averaged by
    // AverageDownConsPerturbational, which reads the fine level as it is and writes the covered
    // coarse cells only, so the scaling of the state around the plain average is skipped
    const bool perturbational = perturbational_avgdown && !if_fluid_mass &&
                                (SolverChoice::terrain_type != TerrainType::EB);

    // ******************************************************************************************
    // First do cell-centered quantities
    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
    // Here we multiply (rho S) by detJ and divide (rho S) by m^2 before average down
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1 && !if_fluid_mass && !perturbational; lev++) {
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

    if (if_fluid_mass) {
        // the fine cell volume, detJ/(m_x m_y), as in the scaling of the plain branch
        const MultiFab& fcons = vars_new[fine_lev][Vars::cons];
        MultiFab vol_fine(fcons.boxArray(), fcons.DistributionMap(), 1, 0);
        for (MFIter mfi(vol_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Array4<      Real> vol_arr = vol_fine.array(mfi);
            const Array4<const Real> mfx_arr = mapfac[fine_lev][MapFacType::m_x]->const_array(mfi);
            const Array4<const Real> mfy_arr = mapfac[fine_lev][MapFacType::m_y]->const_array(mfi);
            if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                const Array4<const Real> detJ_arr = detJ_cc[fine_lev]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vol_arr(i,j,k) = detJ_arr(i,j,k) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vol_arr(i,j,k) = Real(1.0) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                });
            }
        }
        if_average_down(fcons, vars_new[crse_lev][Vars::cons],
                        *terrain_blanking[fine_lev], *terrain_blanking[crse_lev], vol_fine,
                        scomp, ncomp, refRatio(crse_lev));
    } else if (perturbational) {
        AverageDownConsPerturbational(crse_lev, scomp, ncomp);
    } else if (SolverChoice::terrain_type != TerrainType::EB) {
        // Anelastic: rho is not averaged down (rho == rho0 on each level), so averaging
        //    w (rho S) with w = detJ/(m_x m_y) gives rho_c S for a uniform S only if the fine
        //    w (rho) average exactly to the coarse w (rho), which fails when m varies (#4904).
        //    Rescale the covered coarse cells by w_c rho_c / avg(w_f rho_f) so that
        //    (rho S)_c = rho_c * (the same-weighted mean of S_f).  A pseudo-incompressible level
        //    evolves rho and averages it down, as a compressible level does, so it is not rescaled.
        const bool rescale_by_rho = solverChoice.compressibility_type[crse_lev] == CompressibilityType::Anelastic && (scomp > Rho_comp);

        // A flat mesh whose map factors vary: a covered coarse cell's volume 1/(m_x m_y)
        //    (detJ/(m_x m_y) on a stretched mesh) is not the mean of the fine volumes over it,
        //    and over terrain average_down_detJ makes it so, but a flat mesh keeps its own
        //    volumes (see average_down_detJ).  Divide by the sum of the fine volumes instead,
        //    w_c / avg(w_f) times the plain average, so that a uniform state stays uniform.  The
        //    rescaling above already divides by the sum of the fine masses.
        const bool rescale_by_vol = !rescale_by_rho && fine_volume_average(crse_lev);
        MultiFab vol_fac;
        if (rescale_by_vol) {
            MultiFab vol[2];
            for (int l = 0; l < 2; ++l) {
                const int lev = crse_lev + l;
                const MultiFab& cons = vars_new[lev][Vars::cons];
                vol[l].define(cons.boxArray(), cons.DistributionMap(), 1, 0);
                for (MFIter mfi(vol[l], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const Array4<      Real> v_arr   = vol[l].array(mfi);
                    const Array4<const Real> mfx_arr = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                    const Array4<const Real> mfy_arr = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                        const Array4<const Real> detJ_arr = detJ_cc[lev]->const_array(mfi);
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            v_arr(i,j,k) = detJ_arr(i,j,k) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            v_arr(i,j,k) = one / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                        });
                    }
                } // mfi
            } // l

            // Uncovered cells keep vol[0] unchanged, so their factor is exactly one
            vol_fac.define(vol[0].boxArray(), vol[0].DistributionMap(), 1, 0);
            MultiFab::Copy(vol_fac, vol[0], 0, 0, 1, 0);
            average_down(vol[1], vol[0], 0, 1, refRatio(crse_lev));
            MultiFab::Divide(vol_fac, vol[0], 0, 0, 1, 0);
        }

        MultiFab rho_w_fac;
        if (rescale_by_rho) {
            MultiFab rho_w[2];
            for (int l = 0; l < 2; ++l) {
                const int lev = crse_lev + l;
                const MultiFab& cons = vars_new[lev][Vars::cons];
                rho_w[l].define(cons.boxArray(), cons.DistributionMap(), 1, 0);
                for (MFIter mfi(rho_w[l], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const Array4<      Real> rw_arr   = rho_w[l].array(mfi);
                    const Array4<const Real> cons_arr = cons.const_array(mfi);
                    const Array4<const Real> mfx_arr  = mapfac[lev][MapFacType::m_x]->const_array(mfi);
                    const Array4<const Real> mfy_arr  = mapfac[lev][MapFacType::m_y]->const_array(mfi);
                    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
                        const Array4<const Real> detJ_arr = detJ_cc[lev]->const_array(mfi);
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            rw_arr(i,j,k) = cons_arr(i,j,k,Rho_comp) * detJ_arr(i,j,k)
                                          / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                        });
                    } else {
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            rw_arr(i,j,k) = cons_arr(i,j,k,Rho_comp) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
                        });
                    }
                } // mfi
            } // l

            // Uncovered cells keep rho_w[0] unchanged, so their factor is exactly one
            rho_w_fac.define(rho_w[0].boxArray(), rho_w[0].DistributionMap(), 1, 0);
            MultiFab::Copy(rho_w_fac, rho_w[0], 0, 0, 1, 0);
            average_down(rho_w[1], rho_w[0], 0, 1, refRatio(crse_lev));
            MultiFab::Divide(rho_w_fac, rho_w[0], 0, 0, 1, 0);
        }

        average_down(vars_new[crse_lev+1][Vars::cons],vars_new[crse_lev  ][Vars::cons],
                    scomp, ncomp, refRatio(crse_lev));

        if (rescale_by_rho) {
            for (int n = 0; n < ncomp; ++n) {
                MultiFab::Multiply(vars_new[crse_lev][Vars::cons], rho_w_fac, 0, scomp+n, 1, 0);
            }
        }
        if (rescale_by_vol) {
            for (int n = 0; n < ncomp; ++n) {
                MultiFab::Multiply(vars_new[crse_lev][Vars::cons], vol_fac, 0, scomp+n, 1, 0);
            }
        }
    } else {
        // const auto dx = geom[fine_lev].CellSize();
        // Setting cell_vol to the exact value may cause round-off errors in volume average.
        // const Real cell_vol = dx[0]*dx[1]*dx[2];
        constexpr Real cell_vol = one;
        const BoxArray& ba = vars_new[fine_lev][IntVars::cons].boxArray();
        const DistributionMapping& dm = vars_new[fine_lev][IntVars::cons].DistributionMap();
        MultiFab vol_fine(ba, dm, 1, 0);
        vol_fine.setVal(cell_vol);
        EB_average_down(vars_new[fine_lev][Vars::cons],vars_new[crse_lev][Vars::cons],
                    vol_fine, *detJ_cc[fine_lev],
                    scomp, ncomp, refRatio(crse_lev));
    }

    vars_new[crse_lev][Vars::cons].FillBoundary(geom[crse_lev].periodicity());

    // ******************************************************************************************
    // Here we multiply (rho S) by m^2 and divide by detJ after average down
    // ******************************************************************************************
    for (int lev = crse_lev; lev <= crse_lev+1 && !if_fluid_mass && !perturbational; lev++) {
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
    // (This won't be needed because EB_average_down copies the covered value.)
    if (SolverChoice::terrain_type == TerrainType::EB) {
        for (int lev = crse_lev; lev <= crse_lev+1; lev++) {
            for (MFIter mfi(vars_new[lev][Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<      Real> cons_new = vars_new[lev][Vars::cons].array(mfi);
                const Array4<const Real> cons_old = vars_old[lev][Vars::cons].array(mfi);
                const Array4<const Real> detJ_arr = detJ_cc[lev]->const_array(mfi);
                ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (detJ_arr(i,j,k) == zero) {
                        cons_new(i,j,k,scomp+n) = cons_old(i,j,k,scomp+n);
                    }
                });
            } // mfi
        } // lev
    }

    if (!do_full_state_and_momenta) return;

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

    // The momenta of immersed forcing keep the plain face average with either cell average: it
    // is the fine mass flux through the coarse face, which keeps the anelastic constraint (a
    // fluid-weighted face average does not, and advecting theta with it heats the flow).
    // On a pseudo-incompressible level it keeps the coarse mass flux consistent with the
    // averaged-down rho, but satisfies div(beta rho u) = S on the coarse level only to within
    // the variation of beta along a coarse face, since the coarse beta is not the fine one
    if (SolverChoice::terrain_type != TerrainType::EB) {
        // Anelastic or pseudo-incompressible over terrain, or on a flat mesh with map factors other
        // than one: the momenta as mass fluxes (see average_down_mass_flux)
        if (solverChoice.compressibility_type[crse_lev] != CompressibilityType::FullyCompressible && cf_interp_mass_flux(crse_lev+1)) {
            average_down_mass_flux(crse_lev);
        } else {
            average_down_faces(rU_new[crse_lev+1], rU_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
            average_down_faces(rV_new[crse_lev+1], rV_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
            average_down_faces(rW_new[crse_lev+1], rW_new[crse_lev], refRatio(crse_lev), geom[crse_lev]);
        }
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

/**
 * Average the cell-centered state of level crse_lev+1 down onto the covered cells of level
 * crse_lev, for erf.perturbational_average_down.
 *
 * With w = detJ / (m_x m_y) the volume of a cell (1 / (m_x m_y) on a constant-dz mesh), a covered
 * coarse cell gets
 *
 *   q_c = q0_c + sum_f w_f (q_f - q0_f) / (N w_c),
 *
 * with N the number of fine cells in a coarse cell, as in the plain, mass-exact average (on a flat
 * mesh whose map factors vary, sum_f w_f in place of N w_c: see fine_volume_average).  q0 is
 * rho0, rho0 theta0 and (with moisture) rho0 qv0 in Rho_comp, RhoTheta_comp and RhoQ1_comp, and
 * zero otherwise.  With the base state removed, a fine level in its base state averages to the
 * coarse base state.  The fine level is read only.
 *
 * On an anelastic level rho is not averaged down (it is rho0 on each level; a pseudo-incompressible
 * level averages it down, as a compressible level does), so the specific
 * value S = q / rho is averaged with the fine masses w_f rho_f, as the plain anelastic branch of
 * AverageDownTo does (#4904), and the coarse rho multiplies it:
 *
 *   q_c = rho_c [ S0_c + sum_f w_f rho_f (S_f - S0_f) / sum_f w_f rho_f ],
 *
 * with S0 = theta0 and (with moisture) qv0 in RhoTheta_comp and RhoQ1_comp, and zero otherwise.
 *
 * @param[in] crse_lev coarse level
 * @param[in] scomp    first component
 * @param[in] ncomp    number of components
 */
void
ERF::AverageDownConsPerturbational (int crse_lev, int scomp, int ncomp)
{
    const int fine_lev = crse_lev + 1;
    const IntVect rr = refRatio(crse_lev);
    const int rx = rr[0];
    const int ry = rr[1];
    const int rz = rr[2];
    const Real n_fine = Real(rx*ry*rz);

    const MultiFab& fcons = vars_new[fine_lev][Vars::cons];
          MultiFab& ccons = vars_new[crse_lev][Vars::cons];

    const bool anelastic = (solverChoice.compressibility_type[crse_lev] == CompressibilityType::Anelastic);
    AMREX_ALWAYS_ASSERT(!anelastic || scomp > Rho_comp);
    const bool use_detJ = (SolverChoice::mesh_type != MeshType::ConstantDz);
    const bool l_moist  = (solverChoice.moisture_type != MoistureType::None);
    // Divide by the sum of the fine volumes rather than N w_c (see fine_volume_average)
    const bool fine_vol = fine_volume_average(crse_lev);

    // The components, counted from scomp, whose base state is taken out (-1 for none)
    const int ir  = (Rho_comp      >= scomp && Rho_comp      < scomp+ncomp) ? Rho_comp      - scomp : -1;
    const int irt = (RhoTheta_comp >= scomp && RhoTheta_comp < scomp+ncomp) ? RhoTheta_comp - scomp : -1;
    const int iqv = (l_moist &&
                     RhoQ1_comp    >= scomp && RhoQ1_comp    < scomp+ncomp) ? RhoQ1_comp    - scomp : -1;

    // On the fine grids coarsened: sum_f w_f (q_f - q0_f) in components 0 to ncomp-1, sum_f w_f
    //    in component ncomp and sum_f w_f rho_f in component ncomp+1.  On an anelastic level
    //    q0_f is rho_f S0_f, so that the sums are sum_f w_f rho_f (S_f - S0_f)
    BoxArray cba = fcons.boxArray();
    cba.coarsen(rr);
    MultiFab sums(cba, fcons.DistributionMap(), ncomp+2, 0);

    for (MFIter mfi(sums, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real> s_arr   = sums.array(mfi);
        const Array4<const Real> f_arr   = fcons.const_array(mfi);
        const Array4<const Real> b_arr   = base_state[fine_lev].const_array(mfi);
        const Array4<const Real> mfx_arr = mapfac[fine_lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mfy_arr = mapfac[fine_lev][MapFacType::m_y]->const_array(mfi);
        const Array4<const Real> dJ_arr  = (use_detJ) ? detJ_cc[fine_lev]->const_array(mfi)
                                                      : Array4<const Real>{};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int n = 0; n <= ncomp+1; ++n) { s_arr(i,j,k,n) = zero; }
            for (int kk = 0; kk < rz; ++kk) {
            for (int jj = 0; jj < ry; ++jj) {
            for (int ii = 0; ii < rx; ++ii) {
                const int fi = i*rx + ii;
                const int fj = j*ry + jj;
                const int fk = k*rz + kk;
                const Real w = ((use_detJ) ? dJ_arr(fi,fj,fk) : one)
                             / (mfx_arr(fi,fj,0)*mfy_arr(fi,fj,0));
                // The density the base state is scaled by: rho0 (q0 = rho0 S0), or on an
                //    anelastic level the density itself (q0 = rho S0)
                const Real r0 = (anelastic) ? f_arr(fi,fj,fk,Rho_comp)
                                            : b_arr(fi,fj,fk,BaseState::r0_comp);
                for (int n = 0; n < ncomp; ++n) {
                    Real q = f_arr(fi,fj,fk,scomp+n);
                    if (n == ir ) { q -= r0; }
                    if (n == irt) { q -= r0 * b_arr(fi,fj,fk,BaseState::th0_comp); }
                    if (n == iqv) { q -= r0 * b_arr(fi,fj,fk,BaseState::qv0_comp); }
                    s_arr(i,j,k,n) += w * q;
                }
                s_arr(i,j,k,ncomp  ) += w;
                s_arr(i,j,k,ncomp+1) += w * f_arr(fi,fj,fk,Rho_comp);
            }}}
        });
    }

    // On the coarse grids; cells no fine cell covers keep a zero sum of the fine volumes
    MultiFab csums(ccons.boxArray(), ccons.DistributionMap(), ncomp+2, 0);
    csums.setVal(zero);
    csums.ParallelCopy(sums, 0, 0, ncomp+2);

    for (MFIter mfi(ccons, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real> c_arr   = ccons.array(mfi);
        const Array4<const Real> s_arr   = csums.const_array(mfi);
        const Array4<const Real> b_arr   = base_state[crse_lev].const_array(mfi);
        const Array4<const Real> mfx_arr = mapfac[crse_lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mfy_arr = mapfac[crse_lev][MapFacType::m_y]->const_array(mfi);
        const Array4<const Real> dJ_arr  = (use_detJ) ? detJ_cc[crse_lev]->const_array(mfi)
                                                      : Array4<const Real>{};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real wsum = s_arr(i,j,k,ncomp);
            if (wsum <= zero) { return; }
            const Real W = (fine_vol) ? wsum
                                      : n_fine * ((use_detJ) ? dJ_arr(i,j,k) : one)
                                               / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
            if (anelastic) {
                // rho_c is the coarse level's own, and is kept
                const Real rho_c = c_arr(i,j,k,Rho_comp);
                const Real mass  = s_arr(i,j,k,ncomp+1);
                for (int n = 0; n < ncomp; ++n) {
                    Real S0 = zero;
                    if (n == irt) { S0 = b_arr(i,j,k,BaseState::th0_comp); }
                    if (n == iqv) { S0 = b_arr(i,j,k,BaseState::qv0_comp); }
                    c_arr(i,j,k,scomp+n) = rho_c * (S0 + s_arr(i,j,k,n) / mass);
                }
                return;
            }
            const Real r0 = b_arr(i,j,k,BaseState::r0_comp);
            for (int n = 0; n < ncomp; ++n) {
                Real q0 = zero;
                if (n == ir ) { q0 = r0; }
                if (n == irt) { q0 = r0 * b_arr(i,j,k,BaseState::th0_comp); }
                if (n == iqv) { q0 = r0 * b_arr(i,j,k,BaseState::qv0_comp); }
                c_arr(i,j,k,scomp+n) = q0 + s_arr(i,j,k,n) / W;
            }
        });
    }
}

// Conservation-preserving fine-to-coarse average of RhoTheta and the active
// moist conserved components. Uses the same detJ/mfac weighting as AverageDownTo,
// but skips the full-state component checks and the momentum average-down.
void
ERF::AverageDownMoistStateTo (int crse_lev)
{
    AMREX_ALWAYS_ASSERT(solverChoice.coupling_type == CouplingType::TwoWay);
    AMREX_ALWAYS_ASSERT(crse_lev >= 0 && crse_lev < finest_level);

    AverageDownTo(crse_lev, RhoTheta_comp, 1, false);

    const int n_moist = (micro) ? micro->Get_Qstate_Size() : 0;
    if (n_moist > 0) {
        AverageDownTo(crse_lev, RhoQ1_comp, n_moist, false);
    }
}

/**
 * Whether the average-down onto the covered cells of level crse_lev divides by the sum of the fine
 * volumes rather than by N times the coarse volume: on a flat mesh (ConstantDz or StretchedDz)
 * whose map factors are not all one.  The two differ there, by O(dm^2), and a uniform state would
 * not stay uniform.  Over terrain average_down_detJ makes them the same, and with all map factors
 * one they are the same, so neither case changes.
 *
 * @param[in] crse_lev coarse level
 */
bool
ERF::fine_volume_average (int crse_lev) const
{
    return SolverChoice::mesh_type != MeshType::VariableDz &&
           SolverChoice::terrain_type != TerrainType::EB &&
           has_nonunit_mapfac(crse_lev+1);
}

/**
 * Set detJ_cc of level crse_lev, in the cells that level crse_lev+1 covers, to
 *
 *   detJ_c = m_x,c m_y,c (1/N) sum_f detJ_f / (m_x,f m_y,f),
 *
 * with N the number of fine cells in a coarse cell, so that each covered coarse cell's volume
 * detJ/(m_x m_y) is the mean of the fine volumes over it.  AverageDownTo weights the state by that
 * volume before averaging and divides by it afterwards, so with this the average-down both
 * conserves and keeps a uniform value uniform where the map factors vary.  With m = 1 it is the
 * plain average of detJ, computed in the same order (amrex::average_down), so bit for bit the
 * same.  The cells the finer level does not cover, and the ghost cells, are left as they are.
 *
 * @param[in] crse_lev coarse level
 */
void
ERF::average_down_detJ (int crse_lev)
{
    const int fine_lev = crse_lev + 1;
    const MultiFab& dJ_f = *detJ_cc[fine_lev];
          MultiFab& dJ_c = *detJ_cc[crse_lev];

    // On a stretched mesh without terrain detJ is a function of k alone: the projection's
    // divergence uses the stretched dz there (compute_divergence) while its operator and the
    // advection use detJ, so a covered coarse cell whose detJ took the map factors would keep a
    // divergence, and upwind advection would carry it out of the covered cells.  The plain
    // average leaves it as it was.
    if (solverChoice.mesh_type == MeshType::StretchedDz) {
        average_down(dJ_f, dJ_c, 0, 1, refRatio(crse_lev));
        return;
    }

    // The fine volumes, detJ / (m_x m_y)
    MultiFab vol_f(dJ_f.boxArray(), dJ_f.DistributionMap(), 1, 0);
    for (MFIter mfi(vol_f, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real> v_arr   = vol_f.array(mfi);
        const Array4<const Real> dJ_arr  = dJ_f.const_array(mfi);
        const Array4<const Real> mfx_arr = mapfac[fine_lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mfy_arr = mapfac[fine_lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            v_arr(i,j,k) = dJ_arr(i,j,k) / (mfx_arr(i,j,0)*mfy_arr(i,j,0));
        });
    }

    // Their mean over each covered coarse cell, times the coarse m_x m_y
    MultiFab vol_c(dJ_c.boxArray(), dJ_c.DistributionMap(), 1, 0);
    vol_c.setVal(zero);
    average_down(vol_f, vol_c, 0, 1, refRatio(crse_lev));
    const iMultiFab covered = makeFineMask(dJ_c, dJ_f.boxArray(), refRatio(crse_lev));

    for (MFIter mfi(dJ_c, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real> dJ_arr  = dJ_c.array(mfi);
        const Array4<const Real> v_arr   = vol_c.const_array(mfi);
        const Array4<const int > m_arr   = covered.const_array(mfi);
        const Array4<const Real> mfx_arr = mapfac[crse_lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real> mfy_arr = mapfac[crse_lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (m_arr(i,j,k) == 1) {
                dJ_arr(i,j,k) = v_arr(i,j,k) * (mfx_arr(i,j,0)*mfy_arr(i,j,0));
            }
        });
    }
}
