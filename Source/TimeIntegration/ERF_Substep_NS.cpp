
#include <ERF_TI_fast_headers.H>

using namespace amrex;

/**
 * Function for computing the fast RHS with no terrain and variable vertical spacing

 *
 * @param[in   ]  step  which fast time step within each Runge-Kutta step
 * @param[in   ]  nrk   which Runge-Kutta step
 * @param[in   ]  level level of resolution
 * @param[in   ]  finest_level finest level of resolution
 * @param[in   ]  S_slow_rhs slow RHS computed in erf_slow_rhs_pre
 * @param[in   ]  S_prev                           if step == 0, this is S_old, else the previous fast solution
 * @param[in   ]  S_stage_data solution            at previous RK stage
 * @param[in   ]  S_stage_prim primitive variables at previous RK stage
 * @param[in   ]  pi_stage   Exner function      at previous RK stage
 * @param[in   ]  fast_coeffs coefficients for the tridiagonal solve used in the fast integrator
 * @param[  out]  S_data    current solution
 * @param[inout]  lagged_delta_rt
 * @param[inout] avg_xmom: time-averaged x-momentum to be used for updating slow variables
 * @param[inout] avg_ymom: time-averaged y-momentum to be used for updating slow variables
 * @param[inout] avg_zmom: time-averaged z-momentum to be used for updating slow variables
 * @param[in   ]  cc_src source terms for conserved variables
 * @param[in   ]  xmom_src source terms for x-momentum
 * @param[in   ]  ymom_src source terms for y-momentum
 * @param[in   ]  zmom_src source terms for z-momentum
 * @param[in   ]  geom container for geometric information
 * @param[in   ]  gravity magnitude of gravity
 * @param[in   ]  stretched_dz_d
 * @param[in   ]  dtau fast time step
 * @param[in   ]  beta_s  Coefficient which determines how implicit vs explicit the solve is
 * @param[in   ]  facinv inverse factor for time-averaging the momenta
 * @param[in   ]  mapfac   vector of map factors
 * @param[inout]  fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout]  fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 * @param[in   ]  l_use_moisture
 * @param[in   ]  l_reflux should we add fluxes to the FluxRegisters?
 */

void erf_substep_NS (int step, int nrk,
                     int level, int finest_level,
                     Vector<MultiFab>& S_slow_rhs,                   // the slow RHS already computed
                     const Vector<MultiFab>& S_prev,                 // if step == 0, this is S_old, else the previous solution
                     Vector<MultiFab>& S_stage_data,                 // S_stage = S^n, S^* or S^**
                     const  MultiFab & S_stage_prim,                 // Primitive version of S_stage_data[IntVars::cons]
                     const  MultiFab & qt,                           // Total moisture
                     const  MultiFab & pi_stage,                     // Exner function evaluated at last stage
                     const  MultiFab &fast_coeffs,                   // Coeffs for tridiagonal solve
                     Vector<MultiFab>& S_data,                       // S_sum = most recent full solution
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
                     amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                     const Real dtau, const Real beta_s,
                     const Real facinv,
                     Vector<std::unique_ptr<MultiFab>>& mapfac,
                     YAFluxRegister* fr_as_crse,
                     YAFluxRegister* fr_as_fine,
                     bool l_use_moisture,
                     bool l_reflux,
                     const amrex::Real* sinesq_stag_d,
                     const Real l_damp_coef)
{
    //
    // NOTE: for step > 0, S_data and S_prev point to the same MultiFab data!!
    //

    BL_PROFILE_REGION("erf_substep_S()");

    Real beta_1 = 0.5 * (1.0 - beta_s);  // multiplies explicit terms
    Real beta_2 = 0.5 * (1.0 + beta_s);  // multiplies implicit terms

    // How much do we project forward the (rho theta) that is used in the horizontal momentum equations
    Real beta_d = 0.1;

    Real RvOverRd = R_v / R_d;

    bool l_rayleigh_impl_for_w = (sinesq_stag_d != nullptr);

    const Real* dx = geom.CellSize();
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    Real dxi = dxInv[0];
    Real dyi = dxInv[1];

    auto dz_ptr = stretched_dz_d.data();

    const auto& ba = S_stage_data[IntVars::cons].boxArray();
    const auto& dm = S_stage_data[IntVars::cons].DistributionMap();

    MultiFab Delta_rho_theta(        ba                , dm, 1, 1);
    MultiFab Delta_rho_w    (convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,0));

    MultiFab     coeff_A_mf(fast_coeffs, make_alias, 0, 1);
    MultiFab inv_coeff_B_mf(fast_coeffs, make_alias, 1, 1);
    MultiFab     coeff_C_mf(fast_coeffs, make_alias, 2, 1);
    MultiFab     coeff_P_mf(fast_coeffs, make_alias, 3, 1);
    MultiFab     coeff_Q_mf(fast_coeffs, make_alias, 4, 1);

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // This will hold theta extrapolated forward in time
    MultiFab extrap(S_data[IntVars::cons].boxArray(),S_data[IntVars::cons].DistributionMap(),1,1);

    // This will hold the update for (rho) and (rho theta)
    MultiFab temp_rhs(S_stage_data[IntVars::zmom].boxArray(),S_stage_data[IntVars::zmom].DistributionMap(),2,0);

    // This will hold the new x- and y-momenta temporarily (so that we don't overwrite values we need when tiling)
    MultiFab temp_cur_xmom(S_stage_data[IntVars::xmom].boxArray(),S_stage_data[IntVars::xmom].DistributionMap(),1,0);
    MultiFab temp_cur_ymom(S_stage_data[IntVars::ymom].boxArray(),S_stage_data[IntVars::ymom].DistributionMap(),1,0);

    // We assume that in the first step (nrk == 0) we are only doing one substep.
    AMREX_ALWAYS_ASSERT(nrk > 0 || step == 0);

    // *************************************************************************
    // First set up some arrays we'll need
    // *************************************************************************

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_stage_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<const Real>& prev_cons  = S_prev[IntVars::cons].const_array(mfi);
        const Array4<const Real>& prev_zmom =  S_prev[IntVars::zmom].const_array(mfi);

        const Array4<const Real>& stage_cons = S_stage_data[IntVars::cons].const_array(mfi);
        const Array4<const Real>& stage_zmom = S_stage_data[IntVars::zmom].const_array(mfi);

        const Array4<Real>& prev_drho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& prev_drho_theta = Delta_rho_theta.array(mfi);
        const Array4<Real>& lagged_arr      = lagged_delta_rt.array(mfi);
        const Array4<Real>& theta_extrap    = extrap.array(mfi);
        const Array4<const Real>& prim      = S_stage_prim.const_array(mfi);

        Box gbx = mfi.growntilebox(1);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            prev_drho_theta(i,j,k) = prev_cons(i,j,k,RhoTheta_comp) - stage_cons(i,j,k,RhoTheta_comp);

            if (step == 0) {
                theta_extrap(i,j,k) = prev_drho_theta(i,j,k);
            } else {
                theta_extrap(i,j,k) = prev_drho_theta(i,j,k) + beta_d *
                  ( prev_drho_theta(i,j,k) - lagged_arr(i,j,k) );
            }

            // NOTE: qv is not changing over the fast steps so we use the stage data
            Real qv = (l_use_moisture) ? prim(i,j,k,PrimQ1_comp) : 0.0;
            theta_extrap(i,j,k) *= (1.0 + RvOverRd*qv);

            // We define lagged_delta_rt for our next step as the current delta_rt
            // (after using it above to extrapolate theta for this step)
            lagged_arr(i,j,k) = prev_drho_theta(i,j,k);
        });

        // NOTE: We must do this here because for step > 0, prev_zmom and cur_zmom both point to the same data,
        //       so by the time we would use prev_zmom to define zflux, it would have already been over-written.
        Box gtbz = mfi.nodaltilebox(2);
        gtbz.grow(IntVect(1,1,0));
        ParallelFor(gtbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            prev_drho_w(i,j,k) = prev_zmom(i,j,k) - stage_zmom(i,j,k);
        });
    } // mfi

    // *************************************************************************
    // Define updates in the current RK stage
    // *************************************************************************

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_stage_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        const Array4<Real const>& xmom_src_arr   = xmom_src.const_array(mfi);
        const Array4<Real const>& ymom_src_arr   = ymom_src.const_array(mfi);

        const Array4<const Real> & stage_xmom = S_stage_data[IntVars::xmom].const_array(mfi);
        const Array4<const Real> & stage_ymom = S_stage_data[IntVars::ymom].const_array(mfi);
        const Array4<const Real> & qt_arr     = qt.const_array(mfi);

        const Array4<const Real>& slow_rhs_rho_u = S_slow_rhs[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_v = S_slow_rhs[IntVars::ymom].const_array(mfi);

        const Array4<Real>& temp_cur_xmom_arr  = temp_cur_xmom.array(mfi);
        const Array4<Real>& temp_cur_ymom_arr  = temp_cur_ymom.array(mfi);

        const Array4<const Real>& prev_xmom = S_prev[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& prev_ymom = S_prev[IntVars::ymom].const_array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<      Real>& avg_xmom_arr = avg_xmom.array(mfi);
        const Array4<      Real>& avg_ymom_arr = avg_ymom.array(mfi);

        const Array4<const Real>& pi_stage_ca = pi_stage.const_array(mfi);

        const Array4<Real>& theta_extrap = extrap.array(mfi);

        // Map factors
        const Array4<const Real>& mf_ux = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy = mapfac[MapFacType::v_y]->const_array(mfi);

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        if (nrk == 0 and step == 0) { // prev == stage
            ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real new_drho_u = dtau * slow_rhs_rho_u(i,j,k) + dtau * xmom_src_arr(i,j,k);;
                avg_xmom_arr(i,j,k) += facinv*new_drho_u;
                temp_cur_xmom_arr(i,j,k) = stage_xmom(i,j,k) + new_drho_u;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real new_drho_v = dtau * slow_rhs_rho_v(i,j,k) + dtau * ymom_src_arr(i,j,k);
                avg_ymom_arr(i,j,k) += facinv*new_drho_v;
                temp_cur_ymom_arr(i,j,k) = stage_ymom(i,j,k) + new_drho_v;
            });
        } else {
            ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real gpx = (theta_extrap(i,j,k) - theta_extrap(i-1,j,k))*dxi;
                gpx *= mf_ux(i,j,0);

                Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i-1,j,k)) : 0.0;

                Real pi_c =  0.5 * (pi_stage_ca(i-1,j,k,0) + pi_stage_ca(i,j,k,0));
                Real fast_rhs_rho_u = -Gamma * R_d * pi_c * gpx / (1.0 + q);

                Real new_drho_u = prev_xmom(i,j,k) - stage_xmom(i,j,k)
                                + dtau * fast_rhs_rho_u + dtau * slow_rhs_rho_u(i,j,k)
                                + dtau * xmom_src_arr(i,j,k);

                avg_xmom_arr(i,j,k) += facinv*new_drho_u;

                temp_cur_xmom_arr(i,j,k) = stage_xmom(i,j,k) + new_drho_u;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real gpy = (theta_extrap(i,j,k) - theta_extrap(i,j-1,k))*dyi;
                gpy *= mf_vy(i,j,0);

                Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i,j-1,k)) : 0.0;

                Real pi_c =  0.5 * (pi_stage_ca(i,j-1,k,0) + pi_stage_ca(i,j,k,0));
                Real fast_rhs_rho_v = -Gamma * R_d * pi_c * gpy / (1.0 + q);

                Real new_drho_v = prev_ymom(i,j,k) - stage_ymom(i,j,k)
                                + dtau * fast_rhs_rho_v + dtau * slow_rhs_rho_v(i,j,k)
                                + dtau * ymom_src_arr(i,j,k);

                avg_ymom_arr(i,j,k) += facinv*new_drho_v;

                temp_cur_ymom_arr(i,j,k) = stage_ymom(i,j,k) + new_drho_v;
            });
        } // nrk > 0 and/or step > 0
    } //mfi

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
    std::array<FArrayBox,AMREX_SPACEDIM> flux;
    for ( MFIter mfi(S_stage_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();
        Box tbz = surroundingNodes(bx,2);

        Box vbx = mfi.validbox();
        const auto& vbx_hi = ubound(vbx);

        const Array4<Real const>& zmom_src_arr   = zmom_src.const_array(mfi);

        const Array4<const Real>& stage_xmom = S_stage_data[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& stage_ymom = S_stage_data[IntVars::ymom].const_array(mfi);
        const Array4<const Real>& stage_zmom = S_stage_data[IntVars::zmom].const_array(mfi);
        const Array4<const Real> & prim       = S_stage_prim.const_array(mfi);
        const Array4<const Real> & qt_arr     = qt.const_array(mfi);

        const Array4<const Real>& prev_drho_theta = Delta_rho_theta.array(mfi);

        const Array4<const Real>& prev_cons  = S_prev[IntVars::cons].const_array(mfi);
        const Array4<const Real>& stage_cons = S_stage_data[IntVars::cons].const_array(mfi);

        const Array4<const Real>& slow_rhs_cons  = S_slow_rhs[IntVars::cons].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_w = S_slow_rhs[IntVars::zmom].const_array(mfi);

        const Array4<const Real>& prev_zmom = S_prev[IntVars::zmom].const_array(mfi);
        const Array4<      Real>&  cur_zmom = S_data[IntVars::zmom].array(mfi);

        const Array4<Real>& temp_cur_xmom_arr  = temp_cur_xmom.array(mfi);
        const Array4<Real>& temp_cur_ymom_arr  = temp_cur_ymom.array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<      Real>& avg_zmom_arr = avg_zmom.array(mfi);

        // Map factors
        const Array4<const Real>& mf_mx = mapfac[MapFacType::m_x]->const_array(mfi);
        const Array4<const Real>& mf_my = mapfac[MapFacType::m_y]->const_array(mfi);
        const Array4<const Real>& mf_uy = mapfac[MapFacType::u_y]->const_array(mfi);
        const Array4<const Real>& mf_vx = mapfac[MapFacType::v_x]->const_array(mfi);

        FArrayBox RHS_fab;
        RHS_fab.resize(tbz,1, The_Async_Arena());

        FArrayBox soln_fab;
        soln_fab.resize(tbz,1, The_Async_Arena());

        auto const& RHS_a  = RHS_fab.array();
        auto const& soln_a = soln_fab.array();

        auto const& temp_rhs_arr = temp_rhs.array(mfi);

        auto const&     coeffA_a =     coeff_A_mf.array(mfi);
        auto const& inv_coeffB_a = inv_coeff_B_mf.array(mfi);
        auto const&     coeffC_a =     coeff_C_mf.array(mfi);
        auto const&     coeffP_a =     coeff_P_mf.array(mfi);
        auto const&     coeffQ_a =     coeff_Q_mf.array(mfi);

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
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real xflux_lo = (temp_cur_xmom_arr(i  ,j,k) - stage_xmom(i  ,j,k)) / mf_uy(i  ,j,0);
            Real xflux_hi = (temp_cur_xmom_arr(i+1,j,k) - stage_xmom(i+1,j,k)) / mf_uy(i+1,j,0);
            Real yflux_lo = (temp_cur_ymom_arr(i,j  ,k) - stage_ymom(i,j  ,k)) / mf_vx(i,j  ,0);
            Real yflux_hi = (temp_cur_ymom_arr(i,j+1,k) - stage_ymom(i,j+1,k)) / mf_vx(i,j+1,0);

            Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

            temp_rhs_arr(i,j,k,Rho_comp     ) =  ( xflux_hi - xflux_lo ) * dxi * mfsq
                                               + ( yflux_hi - yflux_lo ) * dyi * mfsq;
            temp_rhs_arr(i,j,k,RhoTheta_comp) = (( xflux_hi * (prim(i,j,k,0) + prim(i+1,j,k,0)) -
                                                   xflux_lo * (prim(i,j,k,0) + prim(i-1,j,k,0)) ) * dxi * mfsq +
                                                 ( yflux_hi * (prim(i,j,k,0) + prim(i,j+1,k,0)) -
                                                   yflux_lo * (prim(i,j,k,0) + prim(i,j-1,k,0)) ) * dyi * mfsq) * 0.5;

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

        Box bx_shrunk_in_k = bx;
        int klo = tbz.smallEnd(2);
        int khi = tbz.bigEnd(2);
        bx_shrunk_in_k.setSmall(2,klo+1);
        bx_shrunk_in_k.setBig(2,khi-1);

        // Note that the notes use "g" to mean the magnitude of gravity, so it is positive
        // We set grav_gpu[2] to be the vector component which is negative
        // We define halfg to match the notes (which is why we take the absolute value)
        Real halfg = std::abs(0.5 * grav_gpu[2]);

        // *********************************************************************
        // fast_loop_on_shrunk
        // *********************************************************************
        //Note we don't act on the bottom or top boundaries of the domain
        ParallelFor(bx_shrunk_in_k, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i,j,k-1)) : 0.0;

            Real coeff_P = coeffP_a(i,j,k) / (1.0 + q);
            Real coeff_Q = coeffQ_a(i,j,k) / (1.0 + q);

            Real theta_t_lo  = 0.5 * ( prim(i,j,k-2,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
            Real theta_t_mid = 0.5 * ( prim(i,j,k-1,PrimTheta_comp) + prim(i,j,k  ,PrimTheta_comp) );
            Real theta_t_hi  = 0.5 * ( prim(i,j,k  ,PrimTheta_comp) + prim(i,j,k+1,PrimTheta_comp) );

            Real Omega_kp1 = prev_zmom(i,j,k+1) - stage_zmom(i,j,k+1);
            Real Omega_k   = prev_zmom(i,j,k  ) - stage_zmom(i,j,k  );
            Real Omega_km1 = prev_zmom(i,j,k-1) - stage_zmom(i,j,k-1);

            // line 2 last two terms (order dtau)
            Real old_drho_k   = prev_cons(i,j,k  ,Rho_comp) - stage_cons(i,j,k  ,Rho_comp);
            Real old_drho_km1 = prev_cons(i,j,k-1,Rho_comp) - stage_cons(i,j,k-1,Rho_comp);
            Real R0_tmp = coeff_P * prev_drho_theta(i,j,k) + coeff_Q * prev_drho_theta(i,j,k-1)
                         - halfg * ( old_drho_k + old_drho_km1 );

            // lines 3-5 residuals (order dtau^2) 1.0 <-> beta_2
            Real R1_tmp =  halfg * (-slow_rhs_cons(i,j,k  ,Rho_comp)
                                    -slow_rhs_cons(i,j,k-1,Rho_comp)
                                    +temp_rhs_arr(i,j,k,0) + temp_rhs_arr(i,j,k-1) )
                + ( coeff_P * (slow_rhs_cons(i,j,k  ,RhoTheta_comp) - temp_rhs_arr(i,j,k  ,RhoTheta_comp)) +
                    coeff_Q * (slow_rhs_cons(i,j,k-1,RhoTheta_comp) - temp_rhs_arr(i,j,k-1,RhoTheta_comp)) );

            // lines 6&7 consolidated (reuse Omega & metrics) (order dtau^2)
            Real dz_inv = 1.0 / dz_ptr[k];
            R1_tmp +=  beta_1 * dz_inv * ( (Omega_kp1 - Omega_km1)                         * halfg
                                          -(Omega_kp1*theta_t_hi  - Omega_k  *theta_t_mid) * coeff_P
                                          -(Omega_k  *theta_t_mid - Omega_km1*theta_t_lo ) * coeff_Q );

            // line 1
            RHS_a(i,j,k) = Omega_k + dtau * (slow_rhs_rho_w(i,j,k) + R0_tmp + dtau * beta_2 * R1_tmp + zmom_src_arr(i,j,k));

        }); // bx_shrunk_in_k

        Box b2d = tbz; // Copy constructor
        b2d.setRange(2,0);

        auto const lo = lbound(bx);
        auto const hi = ubound(bx);

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // w at bottom boundary of grid is 0 if at domain boundary, otherwise w = w_old + dtau * slow_rhs
            RHS_a (i,j,lo.z) = prev_zmom(i,j,lo.z) - stage_zmom(i,j,lo.z)
                             + dtau * slow_rhs_rho_w(i,j,lo.z)
                             + dtau * zmom_src_arr(i,j,lo.z);

            // w at top boundary of grid is 0 if at domain boundary, otherwise w = w_old + dtau * slow_rhs
            RHS_a (i,j,hi.z+1) = prev_zmom(i,j,hi.z+1) - stage_zmom(i,j,hi.z+1)
                               + dtau * slow_rhs_rho_w(i,j,hi.z+1)
                               + dtau * zmom_src_arr(i,j,hi.z+1);
        }); // b2d

#ifdef AMREX_USE_GPU
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
          // w = specified Dirichlet value at k = lo.z
            soln_a(i,j,lo.z) = RHS_a(i,j,lo.z) * inv_coeffB_a(i,j,lo.z);
          cur_zmom(i,j,lo.z) = stage_zmom(i,j,lo.z) + soln_a(i,j,lo.z);

          for (int k = lo.z+1; k <= hi.z+1; k++) {
              soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
          }

          cur_zmom(i,j,hi.z+1) = stage_zmom(i,j,hi.z+1) + soln_a(i,j,hi.z+1);

          for (int k = hi.z; k >= lo.z; k--) {
              soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) *soln_a(i,j,k+1);
              cur_zmom(i,j,k) = stage_zmom(i,j,k) + soln_a(i,j,k);
          }
        }); // b2d
#else
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                soln_a(i,j,lo.z) = RHS_a(i,j,lo.z) * inv_coeffB_a(i,j,lo.z);
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
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                cur_zmom(i,j,hi.z+1) = stage_zmom(i,j,hi.z+1) + soln_a(i,j,hi.z+1);
            }
        }
        for (int k = hi.z; k >= lo.z; --k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) * soln_a(i,j,k+1);
                    cur_zmom(i,j,k) = stage_zmom(i,j,k) + soln_a(i,j,k);
                }
            }
        }
#endif
        if (l_rayleigh_impl_for_w) {
            ParallelFor(bx_shrunk_in_k, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
              Real damping_coeff = l_damp_coef * dtau * sinesq_stag_d[k];
              cur_zmom(i,j,k) /= (1.0 + damping_coeff);
            });
        }

        // **************************************************************************
        // Define updates in the RHS of rho and (rho theta)
        // **************************************************************************
        const Array4<Real>&  prev_drho_w = Delta_rho_w.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real zflux_lo = beta_2 * soln_a(i,j,k  ) + beta_1 * prev_drho_w(i,j,k  );
            Real zflux_hi = beta_2 * soln_a(i,j,k+1) + beta_1 * prev_drho_w(i,j,k+1);

            avg_zmom_arr(i,j,k)      += facinv*zflux_lo / (mf_mx(i,j,0) * mf_my(i,j,0));
            if (l_reflux) {
                (flx_arr[2])(i,j,k,0) =        zflux_lo / (mf_mx(i,j,0) * mf_my(i,j,0));
                (flx_arr[2])(i,j,k,1) = (flx_arr[2])(i,j,k,0) * 0.5 * (prim(i,j,k) + prim(i,j,k-1));
            }

            if (k == vbx_hi.z) {
                avg_zmom_arr(i,j,k+1)      += facinv * zflux_hi / (mf_mx(i,j,0) * mf_my(i,j,0));
                if (l_reflux) {
                    (flx_arr[2])(i,j,k+1,0) =          zflux_hi / (mf_mx(i,j,0) * mf_my(i,j,0));
                    (flx_arr[2])(i,j,k+1,1) = (flx_arr[2])(i,j,k+1,0) * 0.5 * (prim(i,j,k) + prim(i,j,k+1));
                }
            }

            Real dz_inv = 1.0 / dz_ptr[k];
            temp_rhs_arr(i,j,k,Rho_comp     ) += dz_inv * ( zflux_hi - zflux_lo );
            temp_rhs_arr(i,j,k,RhoTheta_comp) += 0.5 * dz_inv * ( zflux_hi * (prim(i,j,k) + prim(i,j,k+1))
                                                                - zflux_lo * (prim(i,j,k) + prim(i,j,k-1)) );
        });

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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_stage_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const Array4<      Real>&  cur_cons = S_data[IntVars::cons].array(mfi);
        const Array4<const Real>& prev_cons = S_prev[IntVars::cons].const_array(mfi);
        auto const& temp_rhs_arr     = temp_rhs.const_array(mfi);
        auto const& slow_rhs_cons    = S_slow_rhs[IntVars::cons].const_array(mfi);
        const Array4<Real const>& cc_src_arr   = cc_src.const_array(mfi);

        if (step == 0) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cur_cons(i,j,k,Rho_comp)      = prev_cons(i,j,k,Rho_comp) +
                                                dtau * (slow_rhs_cons(i,j,k,Rho_comp) - temp_rhs_arr(i,j,k,Rho_comp));
                cur_cons(i,j,k,RhoTheta_comp) = prev_cons(i,j,k,RhoTheta_comp) +
                                                dtau * (slow_rhs_cons(i,j,k,RhoTheta_comp) - temp_rhs_arr(i,j,k,RhoTheta_comp));

                // add in source terms for cell-centered conserved variables
                cur_cons(i,j,k,Rho_comp)      += dtau * cc_src_arr(i,j,k,Rho_comp);
                cur_cons(i,j,k,RhoTheta_comp) += dtau * cc_src_arr(i,j,k,RhoTheta_comp);
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                //
                // We didn't need to set cur_cons = prev_cons above because they point to the same data for step > 0
                //
                cur_cons(i,j,k,Rho_comp)      += dtau * (slow_rhs_cons(i,j,k,Rho_comp) - temp_rhs_arr(i,j,k,Rho_comp));
                cur_cons(i,j,k,RhoTheta_comp) += dtau * (slow_rhs_cons(i,j,k,RhoTheta_comp) - temp_rhs_arr(i,j,k,RhoTheta_comp));

                // add in source terms for cell-centered conserved variables
                cur_cons(i,j,k,Rho_comp)      += dtau * cc_src_arr(i,j,k,Rho_comp);
                cur_cons(i,j,k,RhoTheta_comp) += dtau * cc_src_arr(i,j,k,RhoTheta_comp);
            });
        } // step = 0

        const Array4<Real>& cur_xmom = S_data[IntVars::xmom].array(mfi);
        const Array4<Real>& cur_ymom = S_data[IntVars::ymom].array(mfi);

        const Array4<Real const>& temp_cur_xmom_arr = temp_cur_xmom.const_array(mfi);
        const Array4<Real const>& temp_cur_ymom_arr = temp_cur_ymom.const_array(mfi);

        Box tbx = surroundingNodes(bx,0);
        Box tby = surroundingNodes(bx,1);

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cur_xmom(i,j,k) = temp_cur_xmom_arr(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cur_ymom(i,j,k) = temp_cur_ymom_arr(i,j,k);
        });

    } // mfi
}
