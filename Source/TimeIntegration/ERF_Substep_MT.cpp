
#include <ERF_TI_fast_headers.H>

using namespace amrex;

/**
 * Function for computing the fast RHS with moving terrain
 *
 * @param[in]    step  which fast time step within each Runge-Kutta step
 * @param[in   ] nrk   which Runge-Kutta step
 * @param[in   ] level level of resolution
 * @param[in   ] finest_level finest level of resolution
 * @param[in   ] S_slow_rhs slow RHS computed in erf_slow_rhs_pre
 * @param[in   ] S_prev previous solution
 * @param[in   ] S_stg_data solution            at previous RK stage
 * @param[in   ] S_stg_prim primitive variables at previous RK stage
 * @param[in   ] pi_stage   Exner function      at previous RK stage
 * @param[in   ] fast_coeffs coefficients for the tridiagonal solve used in the fast integrator
 * @param[  out] S_data current solution
 * @param[inout] lagged_delta_rt
 * @param[inout] avg_xmom: time-averaged x-momentum to be used for updating slow variables
 * @param[inout] avg_ymom: time-averaged y-momentum to be used for updating slow variables
 * @param[inout] avg_zmom: time-averaged z-momentum to be used for updating slow variables
 * @param[in   ] cc_src source terms for conserved variables
 * @param[in   ] xmom_src source terms for x-momentum
 * @param[in   ] ymom_src source terms for y-momentum
 * @param[in   ] zmom_src source terms for z-momentum
 * @param[in   ] geom container for geometric information
 * @param[in   ] gravity Magnitude of gravity
 * @param[in   ] use_lagged_delta_rt define lagged_delta_rt for our next step
 * @param[in   ] z_t_rk rate of change of grid height -- only relevant for moving terrain
 * @param[in   ] z_t_pert rate of change of grid height -- interpolated between RK stages
 * @param[in   ] z_phys_nd_old height coordinate at nodes at old time
 * @param[in   ] z_phys_nd_new height coordinate at nodes at new time
 * @param[in   ] z_phys_nd_stg height coordinate at nodes at previous stage
 * @param[in   ] detJ_cc_old Jacobian of the metric transformation at old time
 * @param[in   ] detJ_cc_new Jacobian of the metric transformation at new time
 * @param[in   ] detJ_cc_stg Jacobian of the metric transformation at previous stage
 * @param[in   ] dtau fast time step
 * @param[in   ] beta_s  Coefficient which determines how implicit vs explicit the solve is
 * @param[in   ] facinv inverse factor for time-averaging the momenta
 * @param[in   ] mapfac vector of map factors
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 * @param[in   ]  l_use_moisture
 * @param[in   ]  l_reflux should we add fluxes to the FluxRegisters?
 */

void erf_substep_MT (int step, int /*nrk*/,
                     int level, int finest_level,
                     Vector<MultiFab>& S_slow_rhs,                  // the slow RHS already computed
                     const Vector<MultiFab>& S_prev,                // if step == 0, this is S_old, else the previous solution
                     Vector<MultiFab>& S_stg_data,                  // at last RK stg: S^n, S^* or S^**
                     const MultiFab& S_stg_prim,                    // Primitive version of S_stg_data[IntVars::cons]
                     const MultiFab& qt,                            // Total moisture
                     const MultiFab& pi_stage,                      // Exner function evaluated at last RK stg
                     const MultiFab& fast_coeffs,                   // Coeffs for tridiagonal solve
                     Vector<MultiFab>& S_data,                      // S_sum = state at end of this substep
                     MultiFab& lagged_delta_rt,
                     MultiFab& avg_xmom,
                     MultiFab& avg_ymom,
                     MultiFab& avg_zmom,
                     const MultiFab& cc_src,
                     const MultiFab& xmom_src,
                     const MultiFab& ymom_src,
                     const MultiFab& zmom_src,
                     const Geometry geom,
                     const Real gravity,
                     const bool use_lagged_delta_rt,
                     std::unique_ptr<MultiFab>& z_t_rk,             // evaluated from previous RK stg to next RK stg
                     const MultiFab* z_t_pert,                      // evaluated from tau to (tau + delta tau) - z_t_rk
                     std::unique_ptr<MultiFab>& z_phys_nd_old,      // at previous substep time (tau)
                     std::unique_ptr<MultiFab>& z_phys_nd_new,      // at      new substep time (tau + delta tau)
                     std::unique_ptr<MultiFab>& z_phys_nd_stg,      // at last RK stg
                     std::unique_ptr<MultiFab>& detJ_cc_old,        // at previous substep time (tau)
                     std::unique_ptr<MultiFab>& detJ_cc_new,        // at      new substep time (tau + delta tau)
                     std::unique_ptr<MultiFab>& detJ_cc_stg,        // at last RK stg
                     const Real dtau, const Real beta_s,
                     const Real facinv,
                     Vector<std::unique_ptr<MultiFab>>& mapfac,
                     YAFluxRegister* fr_as_crse,
                     YAFluxRegister* fr_as_fine,
                     bool l_use_moisture,
                     bool l_reflux)
{
    BL_PROFILE_REGION("erf_substep_MT()");

    Real beta_1 = 0.5 * (1.0 - beta_s);  // multiplies explicit terms
    Real beta_2 = 0.5 * (1.0 + beta_s);  // multiplies implicit terms

    // How much do we project forward the (rho theta) that is used in the horizontal momentum equations
    Real beta_d = 0.1;

    Real RvOverRd = R_v / R_d;

    const Real* dx = geom.CellSize();
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    Real dxi = dxInv[0];
    Real dyi = dxInv[1];
    Real dzi = dxInv[2];

    MultiFab     coeff_A_mf(fast_coeffs, make_alias, 0, 1);
    MultiFab inv_coeff_B_mf(fast_coeffs, make_alias, 1, 1);
    MultiFab     coeff_C_mf(fast_coeffs, make_alias, 2, 1);
    MultiFab     coeff_P_mf(fast_coeffs, make_alias, 3, 1);
    MultiFab     coeff_Q_mf(fast_coeffs, make_alias, 4, 1);

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    MultiFab extrap(S_data[IntVars::cons].boxArray(),S_data[IntVars::cons].DistributionMap(),1,1);

    MultiFab Omega(S_data[IntVars::zmom].boxArray(), S_data[IntVars::zmom].DistributionMap(), 1, 1);

    // *************************************************************************
    // Define updates in the current RK stg
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
    FArrayBox temp_rhs_fab;

    FArrayBox RHS_fab;
    FArrayBox soln_fab;

    std::array<FArrayBox,AMREX_SPACEDIM> flux;

    //  NOTE: we leave tiling off here for efficiency -- to make this loop work with tiling
    //        will require additional changes
    for ( MFIter mfi(S_stg_data[IntVars::cons],false); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();
        Box tbx = surroundingNodes(bx,0);
        Box tby = surroundingNodes(bx,1);
        Box tbz = surroundingNodes(bx,2);

        Box vbx = mfi.validbox();
        const auto& vbx_hi = ubound(vbx);

        const Array4<Real const>& xmom_src_arr   = xmom_src.const_array(mfi);
        const Array4<Real const>& ymom_src_arr   = ymom_src.const_array(mfi);
        const Array4<Real const>& zmom_src_arr   = zmom_src.const_array(mfi);
        const Array4<Real const>& cc_src_arr     = cc_src.const_array(mfi);

        const Array4<const Real> & stg_cons = S_stg_data[IntVars::cons].const_array(mfi);
        const Array4<const Real> & stg_xmom = S_stg_data[IntVars::xmom].const_array(mfi);
        const Array4<const Real> & stg_ymom = S_stg_data[IntVars::ymom].const_array(mfi);
        const Array4<const Real> & stg_zmom = S_stg_data[IntVars::zmom].const_array(mfi);
        const Array4<const Real> & prim     = S_stg_prim.const_array(mfi);
        const Array4<const Real> & qt_arr   = qt.const_array(mfi);

        const Array4<const Real>& slow_rhs_cons  = S_slow_rhs[IntVars::cons].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_u = S_slow_rhs[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_v = S_slow_rhs[IntVars::ymom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_w = S_slow_rhs[IntVars::zmom].const_array(mfi);

        const Array4<Real>& cur_cons = S_data[IntVars::cons].array(mfi);
        const Array4<Real>& cur_xmom = S_data[IntVars::xmom].array(mfi);
        const Array4<Real>& cur_ymom = S_data[IntVars::ymom].array(mfi);
        const Array4<Real>& cur_zmom = S_data[IntVars::zmom].array(mfi);

        const Array4<Real>& lagged = lagged_delta_rt.array(mfi);

        const Array4<const Real>& prev_cons = S_prev[IntVars::cons].const_array(mfi);
        const Array4<const Real>& prev_xmom = S_prev[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& prev_ymom = S_prev[IntVars::ymom].const_array(mfi);
        const Array4<const Real>& prev_zmom = S_prev[IntVars::zmom].const_array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<Real>& avg_xmom_arr = avg_xmom.array(mfi);
        const Array4<Real>& avg_ymom_arr = avg_ymom.array(mfi);
        const Array4<Real>& avg_zmom_arr = avg_zmom.array(mfi);

        const Array4<const Real>& z_nd_old = z_phys_nd_old->const_array(mfi);
        const Array4<const Real>& z_nd_new = z_phys_nd_new->const_array(mfi);
        const Array4<const Real>& z_nd_stg = z_phys_nd_stg->const_array(mfi);
        const Array4<const Real>& detJ_old =   detJ_cc_old->const_array(mfi);
        const Array4<const Real>& detJ_new =   detJ_cc_new->const_array(mfi);
        const Array4<const Real>& detJ_stg =   detJ_cc_stg->const_array(mfi);

        const Array4<const Real>& z_t_arr  = z_t_rk->const_array(mfi);
        const Array4<const Real>& zp_t_arr = z_t_pert->const_array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        const Array4<const Real>& pi_stage_ca = pi_stage.const_array(mfi);

        const Array4<Real>& theta_extrap = extrap.array(mfi);

        // Map factors
        const Array4<const Real>& mf_mx = mapfac[MapFacType::m_x]->const_array(mfi);
        const Array4<const Real>& mf_my = mapfac[MapFacType::m_y]->const_array(mfi);
        const Array4<const Real>& mf_ux = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy = mapfac[MapFacType::v_y]->const_array(mfi);

        // Note: it is important to grow the tilebox rather than use growntilebox because
        //       we need to fill the ghost cells of the tilebox so we can use them below
        Box gbx   = mfi.tilebox();  gbx.grow(1);
        Box gtbx  = mfi.nodaltilebox(0); gtbx.grow(1); gtbx.setSmall(2,0);
        Box gtby  = mfi.nodaltilebox(1); gtby.grow(1); gtby.setSmall(2,0);

        if (step == 0) {
            ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cur_cons(i,j,k,Rho_comp)      = prev_cons(i,j,k,Rho_comp);
                cur_cons(i,j,k,RhoTheta_comp) = prev_cons(i,j,k,RhoTheta_comp);

                Real delta_rt  =  cur_cons(i,j,k,RhoTheta_comp) - stg_cons(i,j,k,RhoTheta_comp);
                theta_extrap(i,j,k) = delta_rt;

                // NOTE: qv is not changing over the fast steps so we use the stage data
                Real qv = (l_use_moisture) ? prim(i,j,k,PrimQ1_comp) : 0.0;
                theta_extrap(i,j,k) *= (1.0 + RvOverRd*qv);

                // We define lagged_delta_rt for our next step as the current delta_rt
                lagged(i,j,k) = delta_rt;
            });
        } else if (use_lagged_delta_rt) {
            // This is the default for cases with no or static terrain
            ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real delta_rt = cur_cons(i,j,k,RhoTheta_comp) - stg_cons(i,j,k,RhoTheta_comp);
                theta_extrap(i,j,k) = delta_rt + beta_d * (delta_rt - lagged(i,j,k));

                // NOTE: qv is not changing over the fast steps so we use the stage data
                Real qv = (l_use_moisture) ? prim(i,j,k,PrimQ1_comp) : 0.0;
                theta_extrap(i,j,k) *= (1.0 + RvOverRd*qv);

                // We define lagged_delta_rt for our next step as the current delta_rt
                lagged(i,j,k) = delta_rt;
            });
        } else {
            // For the moving wave problem, this choice seems more robust
            ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                theta_extrap(i,j,k) = cur_cons(i,j,k,RhoTheta_comp) - stg_cons(i,j,k,RhoTheta_comp);

                // NOTE: qv is not changing over the fast steps so we use the stage data
                Real qv = (l_use_moisture) ? prim(i,j,k,PrimQ1_comp) : 0.0;
                theta_extrap(i,j,k) *= (1.0 + RvOverRd*qv);
            });
        } // if step

        RHS_fab.resize     (tbz,1, The_Async_Arena());
        soln_fab.resize    (tbz,1, The_Async_Arena());
        temp_rhs_fab.resize(tbz,2, The_Async_Arena());

        auto const& RHS_a        =      RHS_fab.array();
        auto const& soln_a       =     soln_fab.array();
        auto const& temp_rhs_arr = temp_rhs_fab.array();

        auto const&     coeffA_a =     coeff_A_mf.array(mfi);
        auto const& inv_coeffB_a = inv_coeff_B_mf.array(mfi);
        auto const&     coeffC_a =     coeff_C_mf.array(mfi);
        auto const&     coeffP_a =     coeff_P_mf.array(mfi);
        auto const&     coeffQ_a =     coeff_Q_mf.array(mfi);

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        {
        BL_PROFILE("substep_xymom_T");
        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real h_xi_old   = Compute_h_xi_AtIface(i, j, k, dxInv, z_nd_old);
                Real h_zeta_old = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd_old);
                Real gp_xi = (theta_extrap(i,j,k) - theta_extrap(i-1,j,k)) * dxi;
                Real gp_zeta_on_iface = (k == 0) ?
                   0.5  * dzi * ( theta_extrap(i-1,j,k+1) + theta_extrap(i,j,k+1)
                                 -theta_extrap(i-1,j,k  ) - theta_extrap(i,j,k  ) ) :
                   0.25 * dzi * ( theta_extrap(i-1,j,k+1) + theta_extrap(i,j,k+1)
                                 -theta_extrap(i-1,j,k-1) - theta_extrap(i,j,k-1) );
                Real gpx = h_zeta_old * gp_xi - h_xi_old * gp_zeta_on_iface;
                gpx *= mf_ux(i,j,0);

                Real q = (l_use_moisture) ? 0.5 * (qt_arr(i-1,j,k) + qt_arr(i,j,k)) : 0.0;

                Real pi_c =  0.5 * (pi_stage_ca(i-1,j,k,0) + pi_stage_ca(i  ,j,k,0));
                Real fast_rhs_rho_u = -Gamma * R_d * pi_c * gpx / (1.0 + q);

                // We have already scaled the source terms to have the extra factor of dJ
                cur_xmom(i,j,k) = h_zeta_old * prev_xmom(i,j,k) + dtau * fast_rhs_rho_u
                                                                + dtau * slow_rhs_rho_u(i,j,k)
                                                                + dtau * xmom_src_arr(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real h_eta_old  = Compute_h_eta_AtJface(i, j, k, dxInv, z_nd_old);
                Real h_zeta_old = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd_old);
                Real gp_eta = (theta_extrap(i,j,k) -theta_extrap(i,j-1,k)) * dyi;
                Real gp_zeta_on_jface = (k == 0) ?
                    0.5  * dzi * ( theta_extrap(i,j,k+1) + theta_extrap(i,j-1,k+1)
                                  -theta_extrap(i,j,k  ) - theta_extrap(i,j-1,k  ) ) :
                    0.25 * dzi * ( theta_extrap(i,j,k+1) + theta_extrap(i,j-1,k+1)
                                  -theta_extrap(i,j,k-1) - theta_extrap(i,j-1,k-1) );
                Real gpy = h_zeta_old * gp_eta - h_eta_old  * gp_zeta_on_jface;
                gpy *= mf_vy(i,j,0);

                Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j-1,k) + qt_arr(i,j,k)) : 0.0;

                Real pi_c =  0.5 * (pi_stage_ca(i,j-1,k,0) + pi_stage_ca(i,j  ,k,0));
                Real fast_rhs_rho_v = -Gamma * R_d * pi_c * gpy / (1.0 + q);

                // We have already scaled the source terms to have the extra factor of dJ
                cur_ymom(i, j, k) = h_zeta_old * prev_ymom(i,j,k) + dtau * fast_rhs_rho_v
                                                                  + dtau * slow_rhs_rho_v(i,j,k)
                                                                  + dtau * ymom_src_arr(i,j,k);
        });
        } // end profile

        // *************************************************************************
        // Define flux arrays for use in advection
        // *************************************************************************
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flux[dir].resize(surroundingNodes(bx,dir),2,The_Async_Arena());
            flux[dir].setVal<RunOn::Device>(0.);
        }
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>
            flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};

        // *********************************************************************
        {
        BL_PROFILE("fast_T_making_rho_rhs");
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real h_zeta_stg_xlo = Compute_h_zeta_AtIface(i,  j  , k, dxInv, z_nd_stg);
            Real h_zeta_stg_xhi = Compute_h_zeta_AtIface(i+1,j  , k, dxInv, z_nd_stg);
            Real xflux_lo = cur_xmom(i  ,j,k) - stg_xmom(i  ,j,k)*h_zeta_stg_xlo;
            Real xflux_hi = cur_xmom(i+1,j,k) - stg_xmom(i+1,j,k)*h_zeta_stg_xhi;

            Real h_zeta_stg_yhi = Compute_h_zeta_AtJface(i,  j+1, k, dxInv, z_nd_stg);
            Real h_zeta_stg_ylo = Compute_h_zeta_AtJface(i,  j  , k, dxInv, z_nd_stg);
            Real yflux_lo = cur_ymom(i,j  ,k) - stg_ymom(i,j  ,k)*h_zeta_stg_ylo;
            Real yflux_hi = cur_ymom(i,j+1,k) - stg_ymom(i,j+1,k)*h_zeta_stg_yhi;

            // NOTE: we are saving the (1/J) weighting for later when we add this to rho and theta
            temp_rhs_arr(i,j,k,0) = ( ( xflux_hi - xflux_lo ) * dxi + ( yflux_hi - yflux_lo ) * dyi );
            temp_rhs_arr(i,j,k,1) = ( (( xflux_hi * (prim(i,j,k,0) + prim(i+1,j,k,0)) -
                                         xflux_lo * (prim(i,j,k,0) + prim(i-1,j,k,0)) ) * dxi +
                                       ( yflux_hi * (prim(i,j,k,0) + prim(i,j+1,k,0)) -
                                         yflux_lo * (prim(i,j,k,0) + prim(i,j-1,k,0)) ) * dyi) * 0.5 );

            if (l_reflux) {
                (flx_arr[0])(i,j,k,0) = xflux_lo;
                (flx_arr[0])(i,j,k,1) = (flx_arr[0])(i  ,j,k,0) * 0.5 * (prim(i,j,k,0) + prim(i-1,j,k,0));

                (flx_arr[1])(i,j,k,0) = yflux_lo;
                (flx_arr[1])(i,j,k,1) = (flx_arr[1])(i,j  ,k,0) * 0.5 * (prim(i,j,k,0) + prim(i,j-1,k,0));

                if (i == vbx_hi.x) {
                    (flx_arr[0])(i+1,j,k,0) = xflux_hi;
                    (flx_arr[0])(i+1,j,k,1) = (flx_arr[0])(i+1,j,k,0) * 0.5 * (prim(i,j,k,0) + prim(i+1,j,k,0));
                }
                if (j == vbx_hi.y) {
                    (flx_arr[1])(i,j+1,k,0) = yflux_hi;
                    (flx_arr[1])(i,j+1,k,1) = (flx_arr[1])(i,j+1,k,0) * 0.5 * (prim(i,j,k,0) + prim(i,j+1,k,0));
                }
            }
        });
        } // end profile

        // *********************************************************************
        // This must be done before we set cur_xmom and cur_ymom, since those
        //      in fact point to the same array as prev_xmom and prev_ymom
        // *********************************************************************
        Box gbxo = mfi.nodaltilebox(2);
        {
        BL_PROFILE("fast_MT_making_omega");
        Box gbxo_lo = gbxo; gbxo_lo.setBig(2,0);
        ParallelFor(gbxo_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            omega_arr(i,j,k) = 0.;
        });
        Box gbxo_hi = gbxo; gbxo_hi.setSmall(2,gbxo.bigEnd(2));
        ParallelFor(gbxo_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            omega_arr(i,j,k) = prev_zmom(i,j,k) - stg_zmom(i,j,k) - zp_t_arr(i,j,k);
        });
        Box gbxo_mid = gbxo; gbxo_mid.setSmall(2,1); gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
        ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            omega_arr(i,j,k) =
                ( OmegaFromW(i,j,k,prev_zmom(i,j,k),prev_xmom,prev_ymom,mf_ux,mf_vy,z_nd_old,dxInv)
                 -OmegaFromW(i,j,k, stg_zmom(i,j,k), stg_xmom, stg_ymom,mf_ux,mf_vy,z_nd_old,dxInv) )
                - zp_t_arr(i,j,k);
        });
        } // end profile
        // *********************************************************************

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real h_zeta_new = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd_new);
            cur_xmom(i, j, k) /= h_zeta_new;
            avg_xmom_arr(i,j,k) += facinv*(cur_xmom(i,j,k) - stg_xmom(i,j,k));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real h_zeta_new = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd_new);
            cur_ymom(i, j, k) /= h_zeta_new;
            avg_ymom_arr(i,j,k) += facinv*(cur_ymom(i,j,k) - stg_ymom(i,j,k));
        });

        Box bx_shrunk_in_k = bx;
        int klo = tbz.smallEnd(2);
        int khi = tbz.bigEnd(2);
        bx_shrunk_in_k.setSmall(2,klo+1);
        bx_shrunk_in_k.setBig(2,khi-1);

        // Note that the notes use "g" to mean the magnitude of gravity, so it is positive
        // We set grav_gpu[2] to be the vector component which is negative
        // We define halfg to match the notes (which is why we take the absolute value)
        Real halfg = std::abs(0.5 * grav_gpu[2]);

        {
        BL_PROFILE("fast_loop_on_shrunk_t");
        //Note we don't act on the bottom or top boundaries of the domain
        ParallelFor(bx_shrunk_in_k, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k-1) + qt_arr(i,j,k)) : 0.0;

            Real coeff_P = coeffP_a(i,j,k) / (1.0 + q);
            Real coeff_Q = coeffQ_a(i,j,k) / (1.0 + q);

            Real theta_t_lo  = 0.5 * ( prim(i,j,k-2,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
            Real theta_t_mid = 0.5 * ( prim(i,j,k-1,PrimTheta_comp) + prim(i,j,k  ,PrimTheta_comp) );
            Real theta_t_hi  = 0.5 * ( prim(i,j,k  ,PrimTheta_comp) + prim(i,j,k+1,PrimTheta_comp) );

            // line 2 last two terms (order dtau)
            Real R0_tmp  = coeff_P * cur_cons(i,j,k  ,RhoTheta_comp)
                         + coeff_Q * cur_cons(i,j,k-1,RhoTheta_comp)
                         - coeff_P * stg_cons(i,j,k  ,RhoTheta_comp)
                         - coeff_Q * stg_cons(i,j,k-1,RhoTheta_comp)
                       - halfg   * ( cur_cons(i,j,k,Rho_comp) + cur_cons(i,j,k-1,Rho_comp) )
                       + halfg   * ( stg_cons(i,j,k,Rho_comp) + stg_cons(i,j,k-1,Rho_comp) );

            // line 3 residuals (order dtau^2) 1.0 <-> beta_2
            Real R1_tmp = - halfg * ( slow_rhs_cons(i,j,k  ,Rho_comp)
                                    + slow_rhs_cons(i,j,k-1,Rho_comp) )
                          + coeff_P * slow_rhs_cons(i,j,k  ,RhoTheta_comp)
                          + coeff_Q * slow_rhs_cons(i,j,k-1,RhoTheta_comp);

            Real Omega_kp1 = omega_arr(i,j,k+1);
            Real Omega_k   = omega_arr(i,j,k  );
            Real Omega_km1 = omega_arr(i,j,k-1);

            Real detJdiff = (detJ_old(i,j,k) - detJ_old(i,j,k-1)) / (detJ_old(i,j,k)*detJ_old(i,j,k-1));

            // consolidate lines 4&5 (order dtau^2)
            R1_tmp += halfg * ( beta_1 * dzi * (Omega_kp1/detJ_old(i,j,k) + detJdiff*Omega_k - Omega_km1/detJ_old(i,j,k-1))
                                + temp_rhs_arr(i,j,k,Rho_comp)/detJ_old(i,j,k) + temp_rhs_arr(i,j,k-1,Rho_comp)/detJ_old(i,j,k-1) );

            // consolidate lines 6&7 (order dtau^2)
            R1_tmp += -( coeff_P/detJ_stg(i,j,k  ) * ( beta_1 * dzi * (Omega_kp1*theta_t_hi - Omega_k*theta_t_mid) +
                        temp_rhs_arr(i,j,k  ,RhoTheta_comp) )
                       + coeff_Q/detJ_stg(i,j,k-1) * ( beta_1 * dzi * (Omega_k*theta_t_mid - Omega_km1*theta_t_lo) +
                        temp_rhs_arr(i,j,k-1,RhoTheta_comp) ) );

            // line 1
            Real detJ_half = 0.5 * (detJ_stg(i,j,k) + detJ_stg(i,j,k-1)); // TODO: THIS MAY NOT BE RIGHT
            RHS_a(i,j,k) = prev_zmom(i,j,k) - stg_zmom(i,j,k)
                + dtau * (slow_rhs_rho_w(i,j,k) + zmom_src_arr(i,j,k)) / detJ_half
                + dtau * (R0_tmp + dtau*beta_2*R1_tmp);

            // We cannot use omega_arr here since that was built with old_rho_u and old_rho_v ...
            Real UppVpp = OmegaFromW(i,j,k,0.,cur_xmom,cur_ymom,mf_ux,mf_vy,z_nd_new,dxInv)
                        - OmegaFromW(i,j,k,0.,stg_xmom,stg_ymom,mf_ux,mf_vy,z_nd_stg,dxInv);
            RHS_a(i,j,k) += UppVpp;
        });
        } // end profile

        Box b2d = tbz; // Copy constructor
        b2d.setRange(2,0);

        auto const lo = lbound(bx);
        auto const hi = ubound(bx);

        {
        BL_PROFILE("substep_b2d_loop_t");

#ifdef AMREX_USE_GPU
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Moving terrain
            Real rho_on_bdy = 0.5 * ( prev_cons(i,j,lo.z) + prev_cons(i,j,lo.z-1) );
            RHS_a(i,j,lo.z) = rho_on_bdy * zp_t_arr(i,j,0);

            soln_a(i,j,lo.z) = RHS_a(i,j,lo.z) * inv_coeffB_a(i,j,lo.z);

            // w_khi = 0
            RHS_a(i,j,hi.z+1)     =  0.0;

            for (int k = lo.z+1; k <= hi.z+1; k++) {
                soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
            }

            for (int k = hi.z; k >= lo.z; k--) {
                soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) * soln_a(i,j,k+1);
            }

           // We assume that Omega == w at the top boundary and that changes in J there are irrelevant
            cur_zmom(i,j,hi.z+1) = stg_zmom(i,j,hi.z+1) + soln_a(i,j,hi.z+1);
        });
#else
        for (int j = lo.y; j <= hi.y; ++j) {
             AMREX_PRAGMA_SIMD
             for (int i = lo.x; i <= hi.x; ++i) {

                 Real rho_on_bdy = 0.5 * ( prev_cons(i,j,lo.z) + prev_cons(i,j,lo.z-1) );
                 RHS_a(i,j,lo.z) = rho_on_bdy * zp_t_arr(i,j,lo.z);

                 soln_a(i,j,lo.z) = RHS_a(i,j,lo.z) * inv_coeffB_a(i,j,lo.z);
             }
        }

        for (int j = lo.y; j <= hi.y; ++j) {
             AMREX_PRAGMA_SIMD
             for (int i = lo.x; i <= hi.x; ++i) {
                 RHS_a (i,j,hi.z+1) =  0.0;
             }
        }
        for (int k = lo.z+1; k <= hi.z+1; ++k) {
             for (int j = lo.y; j <= hi.y; ++j) {
                 AMREX_PRAGMA_SIMD
                 for (int i = lo.x; i <= hi.x; ++i) {
                     soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
                 }
           }
        }
        for (int k = hi.z; k >= lo.z; --k) {
             for (int j = lo.y; j <= hi.y; ++j) {
                 AMREX_PRAGMA_SIMD
                 for (int i = lo.x; i <= hi.x; ++i) {
                     soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) * soln_a(i,j,k+1);
                 }
             }
        }

        // We assume that Omega == w at the top boundary and that changes in J there are irrelevant
        for (int j = lo.y; j <= hi.y; ++j) {
             AMREX_PRAGMA_SIMD
             for (int i = lo.x; i <= hi.x; ++i) {
                cur_zmom(i,j,hi.z+1) = stg_zmom(i,j,hi.z+1) + soln_a(i,j,hi.z+1);
            }
        }
#endif
        } // end profile

        {
        BL_PROFILE("substep_new_drhow");
        tbz.setBig(2,hi.z);
        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
             Real rho_on_face = 0.5 * (cur_cons(i,j,k,Rho_comp) + cur_cons(i,j,k-1,Rho_comp));

             if (k == lo.z) {
                 cur_zmom(i,j,k) = WFromOmega(i,j,k,rho_on_face*(z_t_arr(i,j,k)+zp_t_arr(i,j,k)),
                                              cur_xmom,cur_ymom,mf_ux,mf_vy,z_nd_new,dxInv);

                 // We need to set this here because it is used to define zflux_lo below
                 soln_a(i,j,k) = 0.;

             } else {

                 Real UppVpp = WFromOmega(i,j,k,0.0,cur_xmom,cur_ymom,mf_ux,mf_vy,z_nd_new,dxInv)
                             - WFromOmega(i,j,k,0.0,stg_xmom,stg_ymom,mf_ux,mf_vy,z_nd_stg,dxInv);
                 Real wpp = soln_a(i,j,k) + UppVpp;
                 Real dJ_old_kface = 0.5 * (detJ_old(i,j,k) + detJ_old(i,j,k-1));
                 Real dJ_new_kface = 0.5 * (detJ_new(i,j,k) + detJ_new(i,j,k-1));

                 cur_zmom(i,j,k) = dJ_old_kface * (stg_zmom(i,j,k) + wpp);
                 cur_zmom(i,j,k) /= dJ_new_kface;

                 soln_a(i,j,k) = OmegaFromW(i,j,k,cur_zmom(i,j,k),cur_xmom,cur_ymom,mf_ux,mf_vy,z_nd_new,dxInv)
                               - OmegaFromW(i,j,k,stg_zmom(i,j,k),stg_xmom,stg_ymom,mf_ux,mf_vy,z_nd_stg,dxInv);
                 soln_a(i,j,k) -= rho_on_face * zp_t_arr(i,j,k);
             }
        });
        } // end profile

        // **************************************************************************
        // Define updates in the RHS of rho and (rho theta)
        // **************************************************************************
        {
        BL_PROFILE("fast_rho_final_update");
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
              Real zflux_lo = beta_2 * soln_a(i,j,k  ) + beta_1 * omega_arr(i,j,k);
              Real zflux_hi = beta_2 * soln_a(i,j,k+1) + beta_1 * omega_arr(i,j,k+1);

              // Note that in the solve we effectively impose new_drho_w(i,j,vbx_hi.z+1)=0
              // so we don't update avg_zmom at k=vbx_hi.z+1
              avg_zmom_arr(i,j,k)      += facinv*zflux_lo / (mf_mx(i,j,0) * mf_my(i,j,0));
              if (l_reflux) {
                  (flx_arr[2])(i,j,k,0) =        zflux_lo / (mf_mx(i,j,0) * mf_my(i,j,0));
              }

              // Note that the factor of (1/J) in the fast source term is canceled
              // when we multiply old and new by detJ_old and detJ_new , respectively
              // We have already scaled the slow source term to have the extra factor of dJ
              Real fast_rhs_rho = -(temp_rhs_arr(i,j,k,0) + ( zflux_hi - zflux_lo ) * dzi);
              Real temp_rho = detJ_old(i,j,k) * cur_cons(i,j,k,0) +
                              dtau * ( slow_rhs_cons(i,j,k,0)*detJ_stg(i,j,k) + fast_rhs_rho );
              cur_cons(i,j,k,0) = temp_rho / detJ_new(i,j,k);

              // Note that the factor of (1/J) in the fast source term is canceled
              // when we multiply old and new by detJ_old and detJ_new , respectively
              // We have already scaled the slow source term to have the extra factor of dJ
              Real fast_rhs_rhotheta = -( temp_rhs_arr(i,j,k,1) + 0.5 *
                                        ( zflux_hi * (prim(i,j,k) + prim(i,j,k+1))
                                        - zflux_lo * (prim(i,j,k) + prim(i,j,k-1)) ) * dzi );
              Real temp_rth = detJ_old(i,j,k) * cur_cons(i,j,k,1) +
                              dtau * ( slow_rhs_cons(i,j,k,1)*detJ_stg(i,j,k) + fast_rhs_rhotheta );
              cur_cons(i,j,k,1) = temp_rth / detJ_new(i,j,k);
              if (l_reflux) {
                  (flx_arr[2])(i,j,k,1) = (flx_arr[2])(i,j,k,0) * 0.5 * (prim(i,j,k) + prim(i,j,k-1));
              }

              if (k == vbx_hi.z) {
                  avg_zmom_arr(i,j,k+1)      += facinv * zflux_hi / (mf_mx(i,j,0) * mf_my(i,j,0));
                  if (l_reflux) {
                      (flx_arr[2])(i,j,k+1,0) =          zflux_hi / (mf_mx(i,j,0) * mf_my(i,j,0));
                      (flx_arr[2])(i,j,k+1,1) = (flx_arr[2])(i,j,k+1,0) * 0.5 * (prim(i,j,k) + prim(i,j,k+1));
                  }
              }

              // add in source terms for cell-centered conserved variables
              cur_cons(i,j,k,Rho_comp)      += dtau * cc_src_arr(i,j,k,Rho_comp);
              cur_cons(i,j,k,RhoTheta_comp) += dtau * cc_src_arr(i,j,k,RhoTheta_comp);
        });
        } // end profile

        // We only add to the flux registers in the final RK step
        if (l_reflux) {
            int strt_comp_reflux = 0;
            // For now we don't reflux (rho theta) because it seems to create issues at c/f boundaries
            int  num_comp_reflux = 1;
            if (level < finest_level) {
                fr_as_crse->CrseAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dtau, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, RunOn::Device);
            }
            if (level > 0) {
                fr_as_fine->FineAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dtau, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, RunOn::Device);
            }

            // This is necessary here so we don't go on to the next FArrayBox without
            // having finished copying the fluxes into the FluxRegisters (since the fluxes
            // are stored in temporary FArrayBox's)
            Gpu::streamSynchronize();

        } // two-way coupling

    } // mfi
    } // OMP
}
