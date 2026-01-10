
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
                         Vector<std::unique_ptr<MultiFab>>& Tau_corr_lev,
                         MultiFab* SmnSmn,
                         MultiFab* eddyDiffs,
                         const Geometry geom,
                         const SolverChoice& solverChoice,
                         std::unique_ptr<SurfaceLayer>& /*SurfLayer*/,
                         Gpu::DeviceVector<Real>& stretched_dz_d,
                         const MultiFab& detJ,
                         Vector<std::unique_ptr<MultiFab>>& mapfac,
                         const MultiFab& ax, const MultiFab& ay, const MultiFab& az,
                         const eb_& ebfact)
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

    const bool do_implicit = (solverChoice.vert_implicit_fac[nrk] > 0) && solverChoice.implicit_momentum_diffusion;

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

            // EB
            Array4<const EBCellFlag> cflag{};
            Array4<const Real> vfrac{};
            Array4<const Real> apx{};
            Array4<const Real> apy{};
            Array4<const Real> apz{};
            if (solverChoice.terrain_type == TerrainType::EB) {
                EBCellFlagFab const& cflag_fab = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi];
                cflag  = cflag_fab.const_array();
                if (cflag_fab.getType(valid_bx) == FabType::singlevalued) {
                    vfrac = (ebfact.get_const_factory())->getVolFrac().const_array(mfi);
                    apx = (ebfact.get_const_factory())->getAreaFrac()[0]->const_array(mfi);
                    apy = (ebfact.get_const_factory())->getAreaFrac()[1]->const_array(mfi);
                    apz = (ebfact.get_const_factory())->getAreaFrac()[2]->const_array(mfi);
                } else {
                    vfrac = detJ.const_array(mfi);
                    apx = ax.const_array(mfi);
                    apy = ay.const_array(mfi);
                    apz = az.const_array(mfi);
                }
            }

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
            Box bx    = mfi.tilebox();
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

            // Symmetric strain/stresses
            S11.resize( bxcc,1,The_Async_Arena()); S22.resize( bxcc,1,The_Async_Arena()); S33.resize( bxcc,1,The_Async_Arena());
            S12.resize(tbxxy,1,The_Async_Arena()); S13.resize(tbxxz,1,The_Async_Arena()); S23.resize(tbxyz,1,The_Async_Arena());
            Array4<Real> s11 = S11.array();  Array4<Real> s22 = S22.array();  Array4<Real> s33 = S33.array();
            Array4<Real> s12 = S12.array();  Array4<Real> s13 = S13.array();  Array4<Real> s23 = S23.array();
            Array4<Real> tau11 = Tau_lev[TauType::tau11]->array(mfi); Array4<Real> tau22 = Tau_lev[TauType::tau22]->array(mfi);
            Array4<Real> tau33 = Tau_lev[TauType::tau33]->array(mfi); Array4<Real> tau12 = Tau_lev[TauType::tau12]->array(mfi);
            Array4<Real> tau13 = Tau_lev[TauType::tau13]->array(mfi); Array4<Real> tau23 = Tau_lev[TauType::tau23]->array(mfi);

            // We cannot simply scale the tau3* terms since our implicit
            // correction to vertical diffusion only applies to the
            // second-order derivatives in the vertical and we don't want to
            // touch the cross terms -- we save the terms here and
            // manipulate them later.
            FArrayBox S13_for_impl, S23_for_impl;
            S13_for_impl.resize(tbxxz,1,The_Async_Arena());
            S23_for_impl.resize(tbxyz,1,The_Async_Arena());
            Array4<Real> s13_corr = (do_implicit) ? S13_for_impl.array() : Array4<Real>{};
            Array4<Real> s23_corr = (do_implicit) ? S23_for_impl.array() : Array4<Real>{};
            Array4<Real> tau13_corr = (do_implicit) ? Tau_corr_lev[0]->array(mfi) : Array4<Real>{};
            Array4<Real> tau23_corr = (do_implicit) ? Tau_corr_lev[1]->array(mfi) : Array4<Real>{};
#ifdef ERF_IMPLICIT_W
            FArrayBox S33_for_impl;
            S33_for_impl.resize( bxcc,1,The_Async_Arena());
            Array4<Real> s33_corr = (do_implicit) ? S33_for_impl.array() : Array4<Real>{};
            Array4<Real> tau33_corr = (do_implicit) ? Tau_corr_lev[2]->array(mfi) : Array4<Real>{};
#else
            Array4<Real> s33_corr = Array4<Real>{};
            Array4<Real> tau33_corr = Array4<Real>{};
#endif

            // Calculate the magnitude of the strain-rate tensor squared if
            // using Deardorff or k-eqn RANS. This contributes to the production
            // term, included in diffusion source in slow RHS post and in the
            // first RK stage (TKE tendencies constant for nrk>0, following WRF)
            Array4<Real> SmnSmn_a = ((nrk==0) && need_SmnSmn) ? SmnSmn->array(mfi) : Array4<Real>{};

            // ****************************************************************
            //
            // These are the steps taken below...
            //
            // 1. Calculate expansion rate
            //    - will be added to the normal strain rates in ComputeStress
            //
            // 2. Call ComputeStrain
            //    - IMPLICIT path: s31_corr and s32_corr are modified in here
            //
            // 3. Call ComputeSmnSmn, if needed for turbulence model
            //
            // 4. Call ComputeStress
            //    - add expansion rates to terms on diagonal
            //    - multiply strain rates by diffusivities, with the total
            //      viscosity calculated as the sum of a constant viscosity (or
            //      constant alpha with mu = rho*alpha) and the eddy viscosity
            //      from the turbulence model
            //    - IMPLICIT path: s33_corr is modified in here
            //
            // 5. Copy temp Sij fabs into Tau_lev multifabs
            //    - stress tensor is symmetric if no terrain and no implicit diffusion
            //    - otherwise, stress tensor is asymmetric
            //
            // ****************************************************************

            if (solverChoice.mesh_type == MeshType::StretchedDz) {
                // Terrain non-symmetric terms
                FArrayBox S21,S31,S32;
                S21.resize(tbxxy,1,The_Async_Arena()); S31.resize(tbxxz,1,The_Async_Arena()); S32.resize(tbxyz,1,The_Async_Arena());
                Array4<Real> s21   = S21.array();       Array4<Real> s31   = S31.array();       Array4<Real> s32   = S32.array();
                Array4<Real> tau21 = Tau_lev[TauType::tau21]->array(mfi);
                Array4<Real> tau31 = Tau_lev[TauType::tau31]->array(mfi);
                Array4<Real> tau32 = Tau_lev[TauType::tau32]->array(mfi);

                // *****************************************************************************
                // Expansion rate compute terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_er_S");
                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real mfsq = mf_mx(i,j,0)*mf_my(i,j,0);
                    er_arr(i,j,k) = (u(i+1, j  , k  )/mf_uy(i+1,j,0) - u(i, j, k)/mf_uy(i,j,0))*dxInv[0] * mfsq +  // == du / (dη/dy) * (1/dξ) * (dξ/dx)*(dη/dy) = du/dx
                                    (v(i  , j+1, k  )/mf_vx(i,j+1,0) - v(i, j, k)/mf_vx(i,j,0))*dxInv[1] * mfsq +  // == dv / (dξ/dx) * (1/dη) * (dξ/dx)*(dη/dy) = dv/dy
                                    (w(i  , j  , k+1)                - w(i, j, k)             )/dz_ptr[k];
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
                                s13_corr, s23_corr);
                } // end profile

                if (SmnSmn_a) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,
                                                        s11,s22,s33,
                                                        s12,s13,s23);
                    });
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
                    ComputeStressConsVisc_S(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                            cell_data,
                                            s11, s22, s33,
                                            s12, s21,
                                            s13, s31,
                                            s23, s32,
                                            er_arr,
                                            mf_mx, mf_ux, mf_vx,
                                            mf_my, mf_uy, mf_vy,
                                            s13_corr, s23_corr, s33_corr);
                } else {
                    ComputeStressVarVisc_S(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s21,
                                           s13, s31,
                                           s23, s32,
                                           er_arr,
                                           mf_mx, mf_ux, mf_vx,
                                           mf_my, mf_uy, mf_vy,
                                           s13_corr, s23_corr, s33_corr);
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
                    if (tau33_corr) tau33_corr(i,j,k) = s33_corr(i,j,k);
                });

                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                    tau21(i,j,k) = s21(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                    tau31(i,j,k) = s31(i,j,k);
                    if (tau13_corr) tau13_corr(i,j,k) = s13_corr(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                    tau32(i,j,k) = s32(i,j,k);
                    if (tau23_corr) tau23_corr(i,j,k) = s23_corr(i,j,k);
                });
                } // end profile

            } else if (l_use_terrain_fitted_coords) {

                // Terrain non-symmetric terms
                FArrayBox S21,S31,S32;
                S21.resize(tbxxy,1,The_Async_Arena()); S31.resize(tbxxz,1,The_Async_Arena()); S32.resize(tbxyz,1,The_Async_Arena());
                Array4<Real> s21   = S21.array();       Array4<Real> s31   = S31.array();       Array4<Real> s32   = S32.array();
                Array4<Real> tau21 = Tau_lev[TauType::tau21]->array(mfi);
                Array4<Real> tau31 = Tau_lev[TauType::tau31]->array(mfi);
                Array4<Real> tau32 = Tau_lev[TauType::tau32]->array(mfi);


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

                    // Note:
                    //   expansionRate ~ du / (dη/dy) * (dz/dζ) * (1/dξ) * (dξ/dx)*(dη/dy)
                    //                 + dv / (dξ/dx) * (dz/dζ) * (1/dη) * (dξ/dx)*(dη/dy)
                    //                 + dΩ/dζ
                    //     ~ (du/dx)*(dz/dζ) + (dv/dy)*(dz/dζ) + dΩ/dζ
                    // Dividing by detJ==dz/dζ gives du/dx + dv/dy + dΩ/dz
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
                                s13_corr, s23_corr);
                } // end profile

                if (SmnSmn_a) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,
                                                        s11,s22,s33,
                                                        s12,s13,s23);
                    });
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
                                            s12, s21,
                                            s13, s31,
                                            s23, s32,
                                            er_arr, z_nd, detJ_arr, dxInv,
                                            mf_mx, mf_ux, mf_vx,
                                            mf_my, mf_uy, mf_vy,
                                            s13_corr, s23_corr, s33_corr);
                } else {
                    ComputeStressVarVisc_T(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s21,
                                           s13, s31,
                                           s23, s32,
                                           er_arr, z_nd, detJ_arr, dxInv,
                                           mf_mx, mf_ux, mf_vx,
                                           mf_my, mf_uy, mf_vy,
                                           s13_corr, s23_corr, s33_corr);
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
                    if (tau33_corr) tau33_corr(i,j,k) = s33_corr(i,j,k);
                });

                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                    tau21(i,j,k) = s21(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                    tau31(i,j,k) = s31(i,j,k);
                    if(tau13_corr) tau13_corr(i,j,k) = s13_corr(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                    tau32(i,j,k) = s32(i,j,k);
                    if(tau23_corr) tau23_corr(i,j,k) = s23_corr(i,j,k);
                });
                } // end profile

            } else {

                // *****************************************************************************
                // Expansion rate compute no terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_er_N");
                if (solverChoice.terrain_type != TerrainType::EB) {
                    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        Real mfsq = mf_mx(i,j,0)*mf_my(i,j,0);
                        er_arr(i,j,k) = (u(i+1, j  , k  )/mf_uy(i+1,j,0) - u(i, j, k)/mf_uy(i,j,0))*dxInv[0]*mfsq +
                                        (v(i  , j+1, k  )/mf_vx(i,j+1,0) - v(i, j, k)/mf_vx(i,j,0))*dxInv[1]*mfsq +
                                        (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
                    });
                } else {
                    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (cflag(i,j,k).isSingleValued()) {
                            er_arr(i,j,k) = (Real(1.0)/vfrac(i,j,k)) * (
                            dxInv[0] * ( apx(i+1,j,k)*u(i+1,j,k) - apx(i,j,k)*u(i,j,k) )
                            + dxInv[1] * ( apy(i,j+1,k)*v(i,j+1,k) - apy(i,j,k)*v(i,j,k) )
                            + dxInv[2] * ( apz(i,j,k+1)*w(i,j,k+1) - apz(i,j,k)*w(i,j,k) ) );
                        } else if (cflag(i,j,k).isRegular()) {
                            er_arr(i,j,k) = (u(i+1, j  , k  ) - u(i, j, k))*dxInv[0] +
                                            (v(i  , j+1, k  ) - v(i, j, k))*dxInv[1] +
                                            (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
                        } else {
                            er_arr(i,j,k) = 0.0;
                        }
                    });
                }
                } // end profile

                // *****************************************************************************
                // Strain tensor compute no terrain
                // *****************************************************************************
                {
                BL_PROFILE("slow_rhs_making_strain_N");
                if (solverChoice.terrain_type != TerrainType::EB) {
                    ComputeStrain_N(bxcc, tbxxy, tbxxz, tbxyz, domain,
                                    u, v, w,
                                    s11, s22, s33,
                                    s12, s13, s23,
                                    dxInv,
                                    mf_mx, mf_ux, mf_vx,
                                    mf_my, mf_uy, mf_vy, bc_ptr_h,
                                    s13_corr, s23_corr);
                } else {
                    ComputeStrain_EB(mfi, bxcc, tbxxy, tbxxz, tbxyz, domain,
                                    u, v, w,
                                    s11, s22, s33,
                                    s12, s13, s23,
                                    dxInv,
                                    bc_ptr_h,
                                    ebfact,
                                    s13_corr, s23_corr);
                }
                } // end profile

                if (SmnSmn_a) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,
                                                        s11,s22,s33,
                                                        s12,s13,s23);
                    });
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
                    if (solverChoice.terrain_type != TerrainType::EB) {
                        ComputeStressConsVisc_N(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                                cell_data,
                                                s11, s22, s33,
                                                s12, s13, s23,
                                                er_arr,
                                                s13_corr, s23_corr, s33_corr);
                    } else {
                        ComputeStressConsVisc_EB(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                                cell_data,
                                                s11, s22, s33,
                                                s12, s13, s23,
                                                er_arr,
                                                vfrac,
                                                s13_corr, s23_corr, s33_corr);
                    }
                } else {
                    ComputeStressVarVisc_N(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           cell_data,
                                           s11, s22, s33,
                                           s12, s13, s23,
                                           er_arr,
                                           s13_corr, s23_corr, s33_corr);
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
                    if (tau33_corr) tau33_corr(i,j,k) = s33_corr(i,j,k);
                });
                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                    if (tau13_corr) tau13_corr(i,j,k) = s13_corr(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                    if (tau23_corr) tau23_corr(i,j,k) = s23_corr(i,j,k);
                });
                } // end profile
            } // no terrain
        } // MFIter
    } // l_use_diff
}


void copy_surface_tau_for_implicit (
    Vector<std::unique_ptr<MultiFab>>& Tau_lev,
    Vector<std::unique_ptr<MultiFab>>& Tau_corr_lev)
{
    // This is only needed if we're using a surface layer, which overwrites the
    // shear stresses at klo -- at the moment, this is for testing

    for ( MFIter mfi(*Tau_lev[TauType::tau11],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Array4<Real> tau13 = Tau_lev[TauType::tau13]->array(mfi);
        Array4<Real> tau23 = Tau_lev[TauType::tau23]->array(mfi);

        Array4<Real> tau13_corr = Tau_corr_lev[0]->array(mfi);
        Array4<Real> tau23_corr = Tau_corr_lev[1]->array(mfi);

        const int klo{0};
        Box bx = mfi.tilebox();
        bx.makeSlab(2,klo);
        Box bxx = surroundingNodes(bx,0);
        Box bxy = surroundingNodes(bx,1);

        ParallelFor(bxx, bxy,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau13_corr(i,j,k) = tau13(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau23_corr(i,j,k) = tau23(i,j,k);
        });
    }
}
