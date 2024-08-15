
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>

#include <TI_slow_headers.H>
#include <EOS.H>
#include <Utils.H>

using namespace amrex;

void erf_make_tau_terms (int level, int nrk,
                         const Vector<BCRec>& domain_bcs_type_h,
                         std::unique_ptr<MultiFab>& z_phys_nd,
                         Vector<MultiFab>& S_data,
                         const MultiFab& xvel,
                         const MultiFab& yvel,
                         const MultiFab& zvel,
                         MultiFab& Omega,
                         MultiFab* Tau11, MultiFab* Tau22, MultiFab* Tau33,
                         MultiFab* Tau12, MultiFab* Tau13, MultiFab* Tau21,
                         MultiFab* Tau23, MultiFab* Tau31, MultiFab* Tau32,
                         MultiFab* SmnSmn,
                         MultiFab* eddyDiffs,
                         const Geometry geom,
                         const SolverChoice& solverChoice,
                         std::unique_ptr<ABLMost>& most,
                         std::unique_ptr<MultiFab>& detJ,
                         std::unique_ptr<MultiFab>& mapfac_m,
                         std::unique_ptr<MultiFab>& mapfac_u,
                         std::unique_ptr<MultiFab>& mapfac_v)
{
    BL_PROFILE_REGION("erf_make_tau_terms()");

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const bool    l_use_terrain    = solverChoice.use_terrain;
    const bool    l_moving_terrain = (solverChoice.terrain_type == TerrainType::Moving);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT (l_use_terrain);


    const bool l_use_diff       = ( (dc.molec_diff_type != MolecDiffType::None) ||
                                    (tc.les_type        !=       LESType::None) ||
                                    (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_constAlpha = ( dc.molec_diff_type == MolecDiffType::ConstantAlpha );
    const bool l_use_turb       = ( tc.les_type == LESType::Smagorinsky ||
                                    tc.les_type == LESType::Deardorff   ||
                                    tc.pbl_type == PBLType::MYNN25      ||
                                    tc.pbl_type == PBLType::YSU );

    const bool use_most     = (most != nullptr);
    const bool exp_most     = (solverChoice.use_explicit_most);
    const bool rot_most     = (solverChoice.use_rotate_most);

    const Box& domain = geom.Domain();
    const int domlo_z = domain.smallEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Pre-computed quantities
    // *****************************************************************************
    int nvars                     = S_data[IntVars::cons].nComp();
    const BoxArray& ba            = S_data[IntVars::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVars::cons].DistributionMap();

    std::unique_ptr<MultiFab> expr;
    std::unique_ptr<MultiFab> dflux_x;
    std::unique_ptr<MultiFab> dflux_y;
    std::unique_ptr<MultiFab> dflux_z;

    if (l_use_diff) {
        expr    = std::make_unique<MultiFab>(ba  , dm, 1, IntVect(1,1,0));
        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, nvars, 0);

        // if using constant alpha (mu = rho * alpha), then first divide by the
        // reference density -- mu_eff will be scaled by the instantaneous
        // local density later when ComputeStress*Visc_*() is called
        Real mu_eff = (l_use_constAlpha) ? 2.0 * dc.dynamicViscosity / dc.rho0_trans
                                         : 2.0 * dc.dynamicViscosity;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& valid_bx = mfi.validbox();

            // Velocities
            const Array4<const Real> & u = xvel.array(mfi);
            const Array4<const Real> & v = yvel.array(mfi);
            const Array4<const Real> & w = zvel.array(mfi);

            // Contravariant velocity
            const Array4<Real>& omega_arr = Omega.array(mfi);

            // Map factors
            const Array4<const Real>& mf_m   = mapfac_m->const_array(mfi);
            const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
            const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

            // Eddy viscosity
            const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};
            const Array4<Real const>& cell_data = l_use_constAlpha ? S_data[IntVars::cons].const_array(mfi) : Array4<const Real>{};

            // Terrain metrics
            const Array4<const Real>& z_nd     = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
            const Array4<const Real>& detJ_arr = detJ->const_array(mfi);

            //-------------------------------------------------------------------------------
            // NOTE: Tile boxes with terrain are not intuitive. The linear combination of
            //       stress terms requires care. Create a tile box that intersects the
            //       valid box, then grow the box in x/y. Compute the strain on the local
            //       FAB over this grown tile box. Compute the stress over the tile box,
            //       except tau_ii which still needs the halo cells. Finally, write from
            //       the local FAB to the Tau MF but only on the tile box.
            //-------------------------------------------------------------------------------

            //-------------------------------------------------------------------------------
            // TODO: Avoid recomputing strain on the first RK stage. One could populate
            //       the FABs with tau_ij, compute stress, and then write to tau_ij. The
            //       problem with this approach is you will over-write the needed halo layer
            //       needed by subsequent tile boxes (particularly S_ii becomes Tau_ii).
            //-------------------------------------------------------------------------------

            // Strain/Stress tile boxes
            Box bxcc  = mfi.tilebox();
            Box tbxxy = mfi.tilebox(IntVect(1,1,0));
            Box tbxxz = mfi.tilebox(IntVect(1,0,1));
            Box tbxyz = mfi.tilebox(IntVect(0,1,1));

            // We need a halo cell for terrain
             bxcc.grow(IntVect(1,1,0));
            tbxxy.grow(IntVect(1,1,0));
            tbxxz.grow(IntVect(1,1,0));
            tbxyz.grow(IntVect(1,1,0));

            if (bxcc.smallEnd(2) != domain.smallEnd(2)) {
                 bxcc.growLo(2,1);
                tbxxy.growLo(2,1);
                tbxxz.growLo(2,1);
                tbxyz.growLo(2,1);
            }

            if (bxcc.bigEnd(2) != domain.bigEnd(2)) {
                 bxcc.growHi(2,1);
                tbxxy.growHi(2,1);
                tbxxz.growHi(2,1);
                tbxyz.growHi(2,1);
            }

            // Expansion rate
            Array4<Real> er_arr = expr->array(mfi);

            // Temporary storage for tiling/OMP
            FArrayBox S11,S22,S33;
            FArrayBox S12,S13,S23;
            S11.resize( bxcc,1,The_Async_Arena());  S22.resize(bxcc,1,The_Async_Arena());  S33.resize(bxcc,1,The_Async_Arena());
            S12.resize(tbxxy,1,The_Async_Arena()); S13.resize(tbxxz,1,The_Async_Arena()); S23.resize(tbxyz,1,The_Async_Arena());
            Array4<Real> s11 = S11.array();  Array4<Real> s22 = S22.array();  Array4<Real> s33 = S33.array();
            Array4<Real> s12 = S12.array();  Array4<Real> s13 = S13.array();  Array4<Real> s23 = S23.array();

            // Symmetric strain/stresses
            Array4<Real> tau11 = Tau11->array(mfi); Array4<Real> tau22 = Tau22->array(mfi); Array4<Real> tau33 = Tau33->array(mfi);
            Array4<Real> tau12 = Tau12->array(mfi); Array4<Real> tau13 = Tau13->array(mfi); Array4<Real> tau23 = Tau23->array(mfi);

            // Strain magnitude
            Array4<Real> SmnSmn_a;

            if (l_use_terrain) {
                // Terrain non-symmetric terms
                FArrayBox S21,S31,S32;
                S21.resize(tbxxy,1,The_Async_Arena()); S31.resize(tbxxz,1,The_Async_Arena()); S32.resize(tbxyz,1,The_Async_Arena());
                Array4<Real> s21   = S21.array();       Array4<Real> s31   = S31.array();       Array4<Real> s32   = S32.array();
                Array4<Real> tau21 = Tau21->array(mfi); Array4<Real> tau31 = Tau31->array(mfi); Array4<Real> tau32 = Tau32->array(mfi);

                // *****************************************************************************
                // Expansion rate compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_er_T");
                // First create Omega using velocity (not momentum)
                Box gbxo = surroundingNodes(bxcc,2);
                ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    omega_arr(i,j,k) = (k == 0) ? 0. : OmegaFromW(i,j,k,w(i,j,k),u,v,z_nd,dxInv);
                });

                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {

                    Real met_u_h_zeta_hi = Compute_h_zeta_AtIface(i+1, j  , k, dxInv, z_nd);
                    Real met_u_h_zeta_lo = Compute_h_zeta_AtIface(i  , j  , k, dxInv, z_nd);

                    Real met_v_h_zeta_hi = Compute_h_zeta_AtJface(i  , j+1, k, dxInv, z_nd);
                    Real met_v_h_zeta_lo = Compute_h_zeta_AtJface(i  , j  , k, dxInv, z_nd);

                    Real Omega_hi = omega_arr(i,j,k+1);
                    Real Omega_lo = omega_arr(i,j,k  );

                    Real mfsq = mf_m(i,j,0)*mf_m(i,j,0);

                    Real expansionRate = (u(i+1,j  ,k)/mf_u(i+1,j,0)*met_u_h_zeta_hi - u(i,j,k)/mf_u(i,j,0)*met_u_h_zeta_lo)*dxInv[0]*mfsq +
                                         (v(i  ,j+1,k)/mf_v(i,j+1,0)*met_v_h_zeta_hi - v(i,j,k)/mf_v(i,j,0)*met_v_h_zeta_lo)*dxInv[1]*mfsq +
                                         (Omega_hi - Omega_lo)*dxInv[2];

                    er_arr(i,j,k) = expansionRate / detJ_arr(i,j,k);
                });
                } // end profile

                // *****************************************************************************
                // Strain tensor compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_strain_T");
                ComputeStrain_T(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                s11, s22, s33,
                                s12, s13,
                                s21, s23,
                                s31, s32,
                                z_nd, detJ_arr, bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
                } // profile

                // Populate SmnSmn if using Deardorff (used as diff src in post)
                // and in the first RK stage (TKE tendencies constant for nrk>0, following WRF)
                if ((nrk==0) && (tc.les_type == LESType::Deardorff)) {
                    SmnSmn_a = SmnSmn->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,s11,s22,s33,s12,s13,s23,domlo_z,use_most,exp_most);
                    });
                }

                // We've updated the strains at all locations including the
                // surface. This is required to get the correct strain-rate
                // magnitude. Now, update the stress everywhere but the surface
                // to retain the values set by MOST.
                if (use_most && exp_most) {
                    // Don't overwrite modeled total stress value at boundary
                    tbxxz.setSmall(2,1);
                    tbxyz.setSmall(2,1);
                    if (rot_most) {
                        bxcc.setSmall(2,1);
                        tbxxy.setSmall(2,1);
                    }
                }

                // *****************************************************************************
                // Stress tensor compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_stress_T");

                // Remove Halo cells just for tau_ij comps
                tbxxy.grow(IntVect(-1,-1,0));
                tbxxz.grow(IntVect(-1,-1,0));
                tbxyz.grow(IntVect(-1,-1,0));

                if (!l_use_turb) {
                    ComputeStressConsVisc_T(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                            cell_data,
                                            s11, s22, s33,
                                            s12, s13,
                                            s21, s23,
                                            s31, s32,
                                            er_arr, z_nd, detJ_arr, dxInv);
                } else {
                    ComputeStressVarVisc_T(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s13,
                                           s21, s23,
                                           s31, s32,
                                           er_arr, z_nd, detJ_arr, dxInv);
                }

                // Remove halo cells from tau_ii but extend across valid_box bdry
                bxcc.grow(IntVect(-1,-1,0));
                if (bxcc.smallEnd(0) == valid_bx.smallEnd(0)) bxcc.growLo(0, 1);
                if (bxcc.bigEnd(0)   == valid_bx.bigEnd(0))   bxcc.growHi(0, 1);
                if (bxcc.smallEnd(1) == valid_bx.smallEnd(1)) bxcc.growLo(1, 1);
                if (bxcc.bigEnd(1)   == valid_bx.bigEnd(1))   bxcc.growHi(1, 1);

                // Copy from temp FABs back to tau
                ParallelFor(bxcc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau11(i,j,k) = s11(i,j,k);
                    tau22(i,j,k) = s22(i,j,k);
                    tau33(i,j,k) = s33(i,j,k);
                });

                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                    tau21(i,j,k) = s21(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                    tau31(i,j,k) = s31(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                    tau32(i,j,k) = s32(i,j,k);
                });
                } // end profile

            } else {

                // *****************************************************************************
                // Expansion rate compute no terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_er_N");
                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    Real mfsq = mf_m(i,j,0)*mf_m(i,j,0);
                    er_arr(i,j,k) = (u(i+1, j  , k  )/mf_u(i+1,j,0) - u(i, j, k)/mf_u(i,j,0))*dxInv[0]*mfsq +
                                    (v(i  , j+1, k  )/mf_v(i,j+1,0) - v(i, j, k)/mf_v(i,j,0))*dxInv[1]*mfsq +
                                    (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
                });
                } // end profile


                // *****************************************************************************
                // Strain tensor compute no terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_strain_N");
                ComputeStrain_N(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                s11, s22, s33,
                                s12, s13, s23,
                                bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
                } // end profile

                // Populate SmnSmn if using Deardorff (used as diff src in post)
                // and in the first RK stage (TKE tendencies constant for nrk>0, following WRF)
                if ((nrk==0) && (tc.les_type == LESType::Deardorff)) {
                    SmnSmn_a = SmnSmn->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,s11,s22,s33,s12,s13,s23,domlo_z,use_most,exp_most);
                    });
                }

                // We've updated the strains at all locations including the
                // surface. This is required to get the correct strain-rate
                // magnitude. Now, update the stress everywhere but the surface
                // to retain the values set by MOST.
                if (use_most && exp_most) {
                    // Don't overwrite modeled total stress value at boundary
                    tbxxz.setSmall(2,1);
                    tbxyz.setSmall(2,1);
                }

                // *****************************************************************************
                // Stress tensor compute no terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_stress_N");

                // Remove Halo cells just for tau_ij comps
                tbxxy.grow(IntVect(-1,-1,0));
                tbxxz.grow(IntVect(-1,-1,0));
                tbxyz.grow(IntVect(-1,-1,0));

                if (!l_use_turb) {
                    ComputeStressConsVisc_N(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                            cell_data,
                                            s11, s22, s33,
                                            s12, s13, s23,
                                            er_arr);
                } else {
                    ComputeStressVarVisc_N(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s13, s23,
                                           er_arr);
                }

                // Remove halo cells from tau_ii but extend across valid_box bdry
                bxcc.grow(IntVect(-1,-1,0));
                if (bxcc.smallEnd(0) == valid_bx.smallEnd(0)) bxcc.growLo(0, 1);
                if (bxcc.bigEnd(0)   == valid_bx.bigEnd(0))   bxcc.growHi(0, 1);
                if (bxcc.smallEnd(1) == valid_bx.smallEnd(1)) bxcc.growLo(1, 1);
                if (bxcc.bigEnd(1)   == valid_bx.bigEnd(1))   bxcc.growHi(1, 1);

                // Copy from temp FABs back to tau
                ParallelFor(bxcc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau11(i,j,k) = s11(i,j,k);
                    tau22(i,j,k) = s22(i,j,k);
                    tau33(i,j,k) = s33(i,j,k);
                });
                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                });
                } // end profile
            } // l_use_terrain
        } // MFIter
    } // l_use_diff
}
