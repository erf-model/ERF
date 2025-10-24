
#include "AMReX_ArrayLim.H"
#include "AMReX_BCRec.H"
#include "AMReX_GpuContainers.H"

#include "ERF_TI_slow_headers.H"
#include "ERF_EOS.H"
#include "ERF_Utils.H"

using namespace amrex;

void erf_make_tau_terms (int level, int nrk,
                         const Vector<BCRec>& domain_bcs_type_h,
                         const MultiFab& z_phys_nd,
                         Vector<MultiFab>& S_data,
                         const MultiFab& xvel,
                         const MultiFab& yvel,
                         const MultiFab& zvel,
                         Vector<std::unique_ptr<MultiFab>>& Tau_lev,
                         MultiFab* SmnSmn,
                         MultiFab* eddyDiffs,
                         const Geometry geom,
                         const SolverChoice& solverChoice,
                         std::unique_ptr<SurfaceLayer>& /*SurfLayer*/,
                         Gpu::DeviceVector<Real>& stretched_dz_d,
                         const MultiFab& detJ,
                         Vector<std::unique_ptr<MultiFab>>& mapfac)
{
    BL_PROFILE_REGION("erf_make_tau_terms()");

    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const bool    l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);
    const bool    l_moving_terrain = (solverChoice.terrain_type == TerrainType::MovingFittedMesh);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT (l_use_terrain_fitted_coords);


    const bool l_use_diff       = ( (dc.molec_diff_type != MolecDiffType::None) || tc.use_kturb );
    const bool l_use_constAlpha = ( dc.molec_diff_type == MolecDiffType::ConstantAlpha );
    const bool l_use_turb       = ( tc.les_type  == LESType::Smagorinsky ||
                                    tc.les_type  == LESType::Deardorff   ||
                                    tc.rans_type == RANSType::kEqn       ||
                                    tc.pbl_type  == PBLType::MYJ         ||
                                    tc.pbl_type  == PBLType::MYNN25      ||
                                    tc.pbl_type  == PBLType::MYNNEDMF    ||
                                    tc.pbl_type  == PBLType::YSU  ||
                                    tc.pbl_type  == PBLType::MRF);

    const bool need_SmnSmn      = (tc.les_type  == LESType::Deardorff ||
                                   tc.rans_type == RANSType::kEqn);

    const Real l_vert_implicit_fac = solverChoice.vert_implicit_fac[nrk];
    const bool need_tau31_tau32 = (solverChoice.mesh_type == MeshType::StretchedDz ||
                                   l_use_terrain_fitted_coords ||
                                   l_vert_implicit_fac > 0);

    const Box& domain = geom.Domain();
    const int domlo_z = domain.smallEnd(2);
    const int domhi_z = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Pre-computed quantities
    // *****************************************************************************
    const BoxArray& ba            = S_data[IntVars::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVars::cons].DistributionMap();

    std::unique_ptr<MultiFab> expr;

    if (l_use_diff) {
        expr    = std::make_unique<MultiFab>(ba, dm, 1, IntVect(1,1,1));

        // if using constant alpha (mu = rho * alpha), then first divide by the
        // reference density -- mu_eff will be scaled by the instantaneous
        // local density later when ComputeStress*Visc_*() is called
        Real mu_eff = (l_use_constAlpha) ? 2.0 * dc.dynamic_viscosity / dc.rho0_trans
                                         : 2.0 * dc.dynamic_viscosity;

        auto dz_ptr = stretched_dz_d.data();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
        {
            const Box& valid_bx = mfi.validbox();

            // Velocities
            const Array4<const Real> & u = xvel.array(mfi);
            const Array4<const Real> & v = yvel.array(mfi);
            const Array4<const Real> & w = zvel.array(mfi);

            // Map factors
            const Array4<const Real>& mf_mx  = mapfac[MapFacType::m_x]->const_array(mfi);
            const Array4<const Real>& mf_ux  = mapfac[MapFacType::u_x]->const_array(mfi);
            const Array4<const Real>& mf_vx  = mapfac[MapFacType::v_x]->const_array(mfi);
            const Array4<const Real>& mf_my  = mapfac[MapFacType::m_y]->const_array(mfi);
            const Array4<const Real>& mf_uy  = mapfac[MapFacType::u_y]->const_array(mfi);
            const Array4<const Real>& mf_vy  = mapfac[MapFacType::v_y]->const_array(mfi);

            // Eddy viscosity
            const Array4<Real const>& mu_turb   = l_use_turb       ? eddyDiffs->const_array(mfi) :
                                                                     Array4<const Real>{};
            const Array4<Real const>& cell_data = l_use_constAlpha ? S_data[IntVars::cons].const_array(mfi) :
                                                                     Array4<const Real>{};

            // Terrain metrics
            const Array4<const Real>& z_nd     = z_phys_nd.const_array(mfi);
            const Array4<const Real>& detJ_arr = detJ.const_array(mfi);

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
            FArrayBox S21,S31,S32;

            // Symmetric strain/stresses
            S11.resize( bxcc,1,The_Async_Arena()); S22.resize( bxcc,1,The_Async_Arena()); S33.resize( bxcc,1,The_Async_Arena());
            S12.resize(tbxxy,1,The_Async_Arena()); S13.resize(tbxxz,1,The_Async_Arena()); S23.resize(tbxyz,1,The_Async_Arena());
            Array4<Real> s11 = S11.array();  Array4<Real> s22 = S22.array();  Array4<Real> s33 = S33.array();
            Array4<Real> s12 = S12.array();  Array4<Real> s13 = S13.array();  Array4<Real> s23 = S23.array();
            Array4<Real> tau11 = Tau_lev[TauType::tau11]->array(mfi); Array4<Real> tau22 = Tau_lev[TauType::tau22]->array(mfi);
            Array4<Real> tau33 = Tau_lev[TauType::tau33]->array(mfi); Array4<Real> tau12 = Tau_lev[TauType::tau12]->array(mfi);
            Array4<Real> tau13 = Tau_lev[TauType::tau13]->array(mfi); Array4<Real> tau23 = Tau_lev[TauType::tau23]->array(mfi);

            // Terrain or implicit non-symmetric terms
            Array4<Real> s21{}, s31{}, s32{};
            Array4<Real> tau21{}, tau31{}, tau32{};
            if (solverChoice.mesh_type == MeshType::StretchedDz ||
                l_use_terrain_fitted_coords)
            {
                S21.resize(tbxxy,1,The_Async_Arena());
                s21 = S21.array();
                tau21 = Tau_lev[TauType::tau21]->array(mfi);
            }
            if (need_tau31_tau32)
            {
                S31.resize(tbxxz,1,The_Async_Arena());
                S32.resize(tbxyz,1,The_Async_Arena());
                s31 = S31.array();
                s32 = S32.array();
                tau31 = Tau_lev[TauType::tau31]->array(mfi);
                tau32 = Tau_lev[TauType::tau32]->array(mfi);
            }

            // Calculate strain-rate magnitude SmnSmn if using Deardorff or
            // or k-eqn RANS (included in diffusion source in post) and in the
            // first RK stage (TKE tendencies constant for nrk>0, following WRF)
            Array4<Real> SmnSmn_a = ((nrk==0) && need_SmnSmn) ? SmnSmn->array(mfi) : Array4<Real>{};

            if (solverChoice.mesh_type == MeshType::StretchedDz) {
                // *****************************************************************************
                // Expansion rate compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_er_S");
                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real mfsq = mf_mx(i,j,0)*mf_my(i,j,0);
                    er_arr(i,j,k) = (u(i+1, j  , k  )/mf_uy(i+1,j,0) - u(i, j, k)/mf_uy(i,j,0))*dxInv[0]*mfsq +
                                    (v(i  , j+1, k  )/mf_vx(i,j+1,0) - v(i, j, k)/mf_vx(i,j,0))*dxInv[1]*mfsq +
                                    (w(i  , j  , k+1) - w(i, j, k))/dz_ptr[k];
                });
                } // end profile

                // *****************************************************************************
                // Strain tensor compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_strain_S");
                ComputeStrain_S(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                u, v, w,
                                s11, s22, s33,
                                s12, s21,
                                s13, s31,
                                s23, s32,
                                stretched_dz_d, dxInv,
                                mf_mx, mf_ux, mf_vx,
                                mf_my, mf_uy, mf_vy, bc_ptr_h,
                                SmnSmn_a,
                                l_vert_implicit_fac);
                } // end profile

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
                    ComputeStressConsVisc_S(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                            cell_data,
                                            s11, s22, s33,
                                            s12, s21,
                                            s13, s31,
                                            s23, s32,
                                            er_arr, stretched_dz_d, dxInv,
                                            mf_mx, mf_ux, mf_vx,
                                            mf_my, mf_uy, mf_vy);
                } else {
                    ComputeStressVarVisc_S(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s21,
                                           s13, s31,
                                           s23, s32,
                                           er_arr, stretched_dz_d, dxInv,
                                           mf_mx, mf_ux, mf_vx,
                                           mf_my, mf_uy, mf_vy);
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

            } else if (l_use_terrain_fitted_coords) {
                // *****************************************************************************
                // Expansion rate compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_er_T");
                Box gbxo = surroundingNodes(bxcc,2);

                // We make a temporary container for contravariant velocity Omega here
                //     -- it is only used to compute er_arr below
                FArrayBox Omega;
                Omega.resize(gbxo,1,The_Async_Arena());

                // First create Omega using velocity (not momentum)
                Array4<Real> omega_arr = Omega.array();
                ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    omega_arr(i,j,k) = (k == 0) ? 0. : OmegaFromW(i,j,k,w(i,j,k),u,v,
                                                                  mf_ux,mf_vy,z_nd,dxInv);
                });

                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {

                    Real met_u_h_zeta_hi = Compute_h_zeta_AtIface(i+1, j  , k, dxInv, z_nd);
                    Real met_u_h_zeta_lo = Compute_h_zeta_AtIface(i  , j  , k, dxInv, z_nd);

                    Real met_v_h_zeta_hi = Compute_h_zeta_AtJface(i  , j+1, k, dxInv, z_nd);
                    Real met_v_h_zeta_lo = Compute_h_zeta_AtJface(i  , j  , k, dxInv, z_nd);

                    Real Omega_hi = omega_arr(i,j,k+1);
                    Real Omega_lo = omega_arr(i,j,k  );

                    Real mfsq = mf_mx(i,j,0)*mf_my(i,j,0);

                    Real expansionRate = (u(i+1,j  ,k)/mf_uy(i+1,j,0)*met_u_h_zeta_hi - u(i,j,k)/mf_uy(i,j,0)*met_u_h_zeta_lo)*dxInv[0]*mfsq +
                                         (v(i  ,j+1,k)/mf_vx(i,j+1,0)*met_v_h_zeta_hi - v(i,j,k)/mf_vx(i,j,0)*met_v_h_zeta_lo)*dxInv[1]*mfsq +
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
                                s12, s21,
                                s13, s31,
                                s23, s32,
                                z_nd, detJ_arr, dxInv,
                                mf_mx, mf_ux, mf_vx,
                                mf_my, mf_uy, mf_vy, bc_ptr_h,
                                SmnSmn_a,
                                l_vert_implicit_fac);
                } // end profile

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
                                            s12, s21,
                                            s13, s31,
                                            s23, s32,
                                            er_arr, z_nd, detJ_arr, dxInv,
                                            mf_mx, mf_ux, mf_vx,
                                            mf_my, mf_uy, mf_vy);
                } else {
                    ComputeStressVarVisc_T(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s21,
                                           s13, s31,
                                           s23, s32,
                                           er_arr, z_nd, detJ_arr, dxInv,
                                           mf_mx, mf_ux, mf_vx,
                                           mf_my, mf_uy, mf_vy);
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
                    Real mfsq = mf_mx(i,j,0)*mf_my(i,j,0);
                    er_arr(i,j,k) = (u(i+1, j  , k  )/mf_uy(i+1,j,0) - u(i, j, k)/mf_uy(i,j,0))*dxInv[0]*mfsq +
                                    (v(i  , j+1, k  )/mf_vx(i,j+1,0) - v(i, j, k)/mf_vx(i,j,0))*dxInv[1]*mfsq +
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
                                s12, /*s21,*/
                                s13, s31,
                                s23, s32,
                                dxInv,
                                mf_mx, mf_ux, mf_vx,
                                mf_my, mf_uy, mf_vy, bc_ptr_h,
                                SmnSmn_a,
                                l_vert_implicit_fac);
                } // end profile

                // *****************************************************************************
                // Stress tensor compute no terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_stress_N");

                // Remove Halo cells just for tau_ij comps
                tbxxy.grow(IntVect(-1,-1,0));
                tbxxz.grow(IntVect(-1,-1,0));
                tbxyz.grow(IntVect(-1,-1,0));
                if (tbxxy.smallEnd(2) > domlo_z) {
                    tbxxy.growLo(2,-1);
                    tbxxz.growLo(2,-1);
                    tbxyz.growLo(2,-1);
                }
                if (tbxxy.bigEnd(2) < domhi_z) {
                    tbxxy.growHi(2,-1);
                    tbxxz.growHi(2,-1);
                    tbxyz.growHi(2,-1);
                }

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
            } // l_use_terrain_fitted_coords
        } // MFIter
    } // l_use_diff
}
