
#include <ERF_TI_fast_headers.H>

using namespace amrex;

/**
 * Function for computing the fast RHS with fixed-in-time terrain
 *
 * @param[in   ] step  which fast time step within each Runge-Kutta step
 * @param[in   ] nrk   which Runge-Kutta step
 * @param[in   ] level level of resolution
 * @param[in   ] finest_level finest level of resolution
 * @param[in   ] S_slow_rhs slow RHS computed in erf_slow_rhs_pre
 * @param[in   ] S_prev previous solution
 * @param[in   ] S_stage_data solution            at previous RK stage
 * @param[in   ] S_stage_prim primitive variables at previous RK stage
 * @param[in   ] pi_stage     Exner function      at previous RK stage
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
 * @param[in   ] gravity magnitude of gravity
 * @param[in   ] z_phys_nd height coordinate at nodes
 * @param[in   ] detJ_cc Jacobian of the metric transformation
 * @param[in   ] dtau fast time step
 * @param[in   ] beta_s  Coefficient which determines how implicit vs explicit the solve is
 * @param[in   ] facinv inverse factor for time-averaging the momenta
 * @param[in   ] mapfac vector of map factors
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 * @param[in   ] l_use_moisture
 * @param[in   ] l_reflux should we add fluxes to the FluxRegisters?
 * @param[in   ] l_damp_coef
 */

void erf_substep_T (int step, int /*nrk*/,
                    int level, int finest_level,
                    Vector<MultiFab>& S_slow_rhs,                   // the slow RHS already computed
                    const Vector<MultiFab>& S_prev,                 // if step == 0, this is S_old, else the previous solution
                    Vector<MultiFab>& S_stage_data,                 // S_stage = S^n, S^* or S^**
                    const MultiFab& S_stage_prim,                   // Primitive version of S_stage_data[IntVars::cons]
                    const MultiFab& qt,                             // Total moisture
                    const MultiFab& pi_stage,                       // Exner function evaluated at last stage
                    const MultiFab& fast_coeffs,                    // Coeffs for tridiagonal solve
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
                    std::unique_ptr<MultiFab>& z_phys_nd,
                    std::unique_ptr<MultiFab>& detJ_cc,
                    const Real dtau, const Real beta_s,
                    const Real facinv,
                    Vector<std::unique_ptr<MultiFab>>& mapfac,
                    YAFluxRegister* fr_as_crse,
                    YAFluxRegister* fr_as_fine,
                    bool l_use_moisture,
                    bool l_reflux,
                    const Real* sinesq_stag_d,
                    const Real l_damp_coef)
{
    BL_PROFILE_REGION("erf_substep_T()");

    const Box& domain = geom.Domain();
    auto const domlo = lbound(domain);
    auto const domhi = ubound(domain);

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
    Real dzi = dxInv[2];
    const auto& ba = S_stage_data[IntVars::cons].boxArray();
    const auto& dm = S_stage_data[IntVars::cons].DistributionMap();

    MultiFab Delta_rho_u(    convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab Delta_rho_v(    convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab Delta_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,0));
    MultiFab Delta_rho  (            ba                , dm, 1, 1);
    MultiFab Delta_rho_theta(        ba                , dm, 1, 1);

    MultiFab New_rho_u(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab New_rho_v(convert(ba,IntVect(0,1,0)), dm, 1, 1);

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

    // *************************************************************************
    // First set up some arrays we'll need
    // *************************************************************************

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_stage_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<Real>       & cur_cons  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real>& prev_cons  = S_prev[IntVars::cons].const_array(mfi);
        const Array4<const Real>& stage_cons = S_stage_data[IntVars::cons].const_array(mfi);
        const Array4<Real>& lagged_arr       = lagged_delta_rt.array(mfi);

        const Array4<Real>& old_drho       = Delta_rho.array(mfi);
        const Array4<Real>& old_drho_u     = Delta_rho_u.array(mfi);
        const Array4<Real>& old_drho_v     = Delta_rho_v.array(mfi);
        const Array4<Real>& old_drho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& old_drho_theta = Delta_rho_theta.array(mfi);

        const Array4<const Real>&  prev_xmom = S_prev[IntVars::xmom].const_array(mfi);
        const Array4<const Real>&  prev_ymom = S_prev[IntVars::ymom].const_array(mfi);
        const Array4<const Real>&  prev_zmom = S_prev[IntVars::zmom].const_array(mfi);

        const Array4<const Real>& stage_xmom = S_stage_data[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& stage_ymom = S_stage_data[IntVars::ymom].const_array(mfi);
        const Array4<const Real>& stage_zmom = S_stage_data[IntVars::zmom].const_array(mfi);

        Box  bx = mfi.validbox();
        Box gbx = mfi.tilebox(); gbx.grow(1);

        if (step == 0) {
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cur_cons(i,j,k,Rho_comp)      = prev_cons(i,j,k,Rho_comp);
                cur_cons(i,j,k,RhoTheta_comp) = prev_cons(i,j,k,RhoTheta_comp);
            });
        } // step = 0

        Box gtbx = mfi.nodaltilebox(0); gtbx.grow(IntVect(1,1,0));
        Box gtby = mfi.nodaltilebox(1); gtby.grow(IntVect(1,1,0));
        Box gtbz = mfi.nodaltilebox(2); gtbz.grow(IntVect(1,1,0));

        const auto& bx_lo = lbound(bx);
        const auto& bx_hi = ubound(bx);

        ParallelFor(gtbx, gtby, gtbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_u(i,j,k) = prev_xmom(i,j,k) - stage_xmom(i,j,k);
            if (k == bx_lo.z && k != domlo.z) {
                old_drho_u(i,j,k-1) = old_drho_u(i,j,k);
            } else if (k == bx_hi.z) {
                old_drho_u(i,j,k+1) = old_drho_u(i,j,k);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_v(i,j,k) = prev_ymom(i,j,k) - stage_ymom(i,j,k);
            if (k == bx_lo.z && k != domlo.z) {
                old_drho_v(i,j,k-1) = old_drho_v(i,j,k);
            } else if (k == bx_hi.z) {
                old_drho_v(i,j,k+1) = old_drho_v(i,j,k);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_w(i,j,k) = prev_zmom(i,j,k) - stage_zmom(i,j,k);
        });

        const Array4<Real>& theta_extrap = extrap.array(mfi);
        const Array4<const Real>& prim   = S_stage_prim.const_array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho(i,j,k)       = cur_cons(i,j,k,Rho_comp)      - stage_cons(i,j,k,Rho_comp);
            old_drho_theta(i,j,k) = cur_cons(i,j,k,RhoTheta_comp) - stage_cons(i,j,k,RhoTheta_comp);
            if (step == 0) {
                theta_extrap(i,j,k) = old_drho_theta(i,j,k);
            } else {
                theta_extrap(i,j,k) = old_drho_theta(i,j,k) + beta_d *
                  ( old_drho_theta(i,j,k) - lagged_arr(i,j,k) );
            }

            // NOTE: qv is not changing over the fast steps so we use the stage data
            Real qv = (l_use_moisture) ? prim(i,j,k,PrimQ1_comp) : 0.0;
            theta_extrap(i,j,k) *= (1.0 + RvOverRd*qv);
        });
    } // mfi

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_stage_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // We define lagged_delta_rt for our next step as the current delta_rt
        Box gbx = mfi.tilebox(); gbx.grow(1);
        const Array4<Real>& old_drho_theta = Delta_rho_theta.array(mfi);
        const Array4<Real>& lagged_arr     = lagged_delta_rt.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            lagged_arr(i,j,k) = old_drho_theta(i,j,k);
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
        Box  bx = mfi.validbox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        const Array4<Real const>& xmom_src_arr   = xmom_src.const_array(mfi);
        const Array4<Real const>& ymom_src_arr   = ymom_src.const_array(mfi);

        const Array4<const Real> & stage_xmom = S_stage_data[IntVars::xmom].const_array(mfi);
        const Array4<const Real> & stage_ymom = S_stage_data[IntVars::ymom].const_array(mfi);
        const Array4<const Real> & qt_arr     = qt.const_array(mfi);

        const Array4<Real>& old_drho_u     = Delta_rho_u.array(mfi);
        const Array4<Real>& old_drho_v     = Delta_rho_v.array(mfi);

        const Array4<const Real>& slow_rhs_rho_u = S_slow_rhs[IntVars::xmom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_v = S_slow_rhs[IntVars::ymom].const_array(mfi);

        const Array4<Real>& new_drho_u = New_rho_u.array(mfi);
        const Array4<Real>& new_drho_v = New_rho_v.array(mfi);

        const Array4<Real>& cur_xmom = S_data[IntVars::xmom].array(mfi);
        const Array4<Real>& cur_ymom = S_data[IntVars::ymom].array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<Real>& avg_xmom_arr = avg_xmom.array(mfi);
        const Array4<Real>& avg_ymom_arr = avg_ymom.array(mfi);

        const Array4<const Real>& z_nd   = z_phys_nd->const_array(mfi);

        const Array4<const Real>& pi_stage_ca = pi_stage.const_array(mfi);

        const Array4<Real>& theta_extrap = extrap.array(mfi);

        // Map factors
        const Array4<const Real>& mf_ux = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy = mapfac[MapFacType::v_y]->const_array(mfi);

        // Create old_drho_u/v/w/theta  = U'', V'', W'', Theta'' in the docs
        // Note that we do the Copy and Subtract including one ghost cell
        //    so that we don't have to fill ghost cells of the new MultiFabs
        // Initialize New_rho_u/v/w to Delta_rho_u/v/w so that
        // the ghost cells in New_rho_u/v/w will match old_drho_u/v/w

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        {
        BL_PROFILE("substep_xymom_T");

        const auto& bx_lo = lbound(bx);
        const auto& bx_hi = ubound(bx);

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real met_h_xi   = Compute_h_xi_AtIface  (i, j, k, dxInv, z_nd);
                Real met_h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);
                Real gp_xi = (theta_extrap(i,j,k) - theta_extrap(i-1,j,k)) * dxi;
                Real gp_zeta_on_iface = (k == 0) ?
                   0.5  * dzi * ( theta_extrap(i-1,j,k+1) + theta_extrap(i,j,k+1)
                                 -theta_extrap(i-1,j,k  ) - theta_extrap(i,j,k  ) ) :
                   0.25 * dzi * ( theta_extrap(i-1,j,k+1) + theta_extrap(i,j,k+1)
                                 -theta_extrap(i-1,j,k-1) - theta_extrap(i,j,k-1) );
                Real gpx = gp_xi - (met_h_xi / met_h_zeta) * gp_zeta_on_iface;
                gpx *= mf_ux(i,j,0);

                Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i-1,j,k)) : 0.0;

                Real pi_c =  0.5 * (pi_stage_ca(i-1,j,k,0) + pi_stage_ca(i  ,j,k,0));
                Real fast_rhs_rho_u = -Gamma * R_d * pi_c * gpx / (1.0 + q);

                new_drho_u(i, j, k) = old_drho_u(i,j,k) + dtau * fast_rhs_rho_u
                                                        + dtau * slow_rhs_rho_u(i,j,k)
                                                        + dtau * xmom_src_arr(i,j,k);
                if (k == bx_lo.z && k != domlo.z) {
                    new_drho_u(i,j,k-1) = new_drho_u(i,j,k);
                } else if (k == bx_hi.z) {
                    new_drho_u(i,j,k+1) = new_drho_u(i,j,k);
                }

                avg_xmom_arr(i,j,k) += facinv*new_drho_u(i,j,k);

                cur_xmom(i,j,k) = stage_xmom(i,j,k) + new_drho_u(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real met_h_eta  = Compute_h_eta_AtJface(i, j, k, dxInv, z_nd);
                Real met_h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
                Real gp_eta = (theta_extrap(i,j,k) -theta_extrap(i,j-1,k)) * dyi;
                Real gp_zeta_on_jface = (k == 0) ?
                    0.5  * dzi * ( theta_extrap(i,j,k+1) + theta_extrap(i,j-1,k+1)
                                  -theta_extrap(i,j,k  ) - theta_extrap(i,j-1,k  ) ) :
                    0.25 * dzi * ( theta_extrap(i,j,k+1) + theta_extrap(i,j-1,k+1)
                                  -theta_extrap(i,j,k-1) - theta_extrap(i,j-1,k-1) );
                Real gpy = gp_eta - (met_h_eta / met_h_zeta) * gp_zeta_on_jface;
                gpy *= mf_vy(i,j,0);

                Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i,j-1,k)) : 0.0;

                Real pi_c =  0.5 * (pi_stage_ca(i,j-1,k,0) + pi_stage_ca(i,j  ,k,0));
                Real fast_rhs_rho_v = -Gamma * R_d * pi_c * gpy / (1.0 + q);

                new_drho_v(i, j, k) = old_drho_v(i,j,k) + dtau * fast_rhs_rho_v
                                                        + dtau * slow_rhs_rho_v(i,j,k)
                                                        + dtau * ymom_src_arr(i,j,k);

                if (k == bx_lo.z && k != domlo.z) {
                    new_drho_v(i,j,k-1) = new_drho_v(i,j,k);
                } else if (k == bx_hi.z) {
                    new_drho_v(i,j,k+1) = new_drho_v(i,j,k);
                }

                avg_ymom_arr(i,j,k) += facinv*new_drho_v(i,j,k);

                cur_ymom(i,j,k) = stage_ymom(i,j,k) + new_drho_v(i,j,k);
            });
        } // end profile
    }

    MultiFab Omega(S_data[IntVars::zmom].boxArray(), dm, 1, 1);

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
        const Array4<Real const>& cc_src_arr     = cc_src.const_array(mfi);

        const Array4<const Real> & stage_zmom = S_stage_data[IntVars::zmom].const_array(mfi);
        const Array4<const Real> & prim       = S_stage_prim.const_array(mfi);
        const Array4<const Real> & qt_arr     = qt.const_array(mfi);

        const Array4<Real>& old_drho_u     = Delta_rho_u.array(mfi);
        const Array4<Real>& old_drho_v     = Delta_rho_v.array(mfi);
        const Array4<Real>& old_drho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& old_drho       = Delta_rho.array(mfi);
        const Array4<Real>& old_drho_theta = Delta_rho_theta.array(mfi);

        const Array4<const Real>& slow_rhs_cons  = S_slow_rhs[IntVars::cons].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_w = S_slow_rhs[IntVars::zmom].const_array(mfi);

        const Array4<Real>& new_drho_u = New_rho_u.array(mfi);
        const Array4<Real>& new_drho_v = New_rho_v.array(mfi);

        const Array4<Real>& cur_cons = S_data[IntVars::cons].array(mfi);
        const Array4<Real>& cur_zmom = S_data[IntVars::zmom].array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<Real>& avg_zmom_arr = avg_zmom.array(mfi);

        const Array4<const Real>& z_nd   = z_phys_nd->const_array(mfi);
        const Array4<const Real>& detJ   = detJ_cc->const_array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        // Map factors
        const Array4<const Real>& mf_mx = mapfac[MapFacType::m_x]->const_array(mfi);
        const Array4<const Real>& mf_my = mapfac[MapFacType::m_y]->const_array(mfi);
        const Array4<const Real>& mf_ux = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_uy = mapfac[MapFacType::u_y]->const_array(mfi);
        const Array4<const Real>& mf_vx = mapfac[MapFacType::v_x]->const_array(mfi);
        const Array4<const Real>& mf_vy = mapfac[MapFacType::v_y]->const_array(mfi);

        // Create old_drho_u/v/w/theta  = U'', V'', W'', Theta'' in the docs
        // Note that we do the Copy and Subtract including one ghost cell
        //    so that we don't have to fill ghost cells of the new MultiFabs
        // Initialize New_rho_u/v/w to Delta_rho_u/v/w so that
        // the ghost cells in New_rho_u/v/w will match old_drho_u/v/w

        FArrayBox temp_rhs_fab;
        FArrayBox RHS_fab;
        FArrayBox soln_fab;

        RHS_fab.resize     (tbz,1,The_Async_Arena());
        soln_fab.resize    (tbz,1,The_Async_Arena());
        temp_rhs_fab.resize(tbz,2,The_Async_Arena());

        auto const& RHS_a        =      RHS_fab.array();
        auto const& soln_a       =     soln_fab.array();
        auto const& temp_rhs_arr = temp_rhs_fab.array();

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
        {
        BL_PROFILE("fast_T_making_rho_rhs");
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real h_zeta_cc_xface_hi = 0.5 * dzi *
              (  z_nd(i+1,j  ,k+1) + z_nd(i+1,j+1,k+1)
                -z_nd(i+1,j  ,k  ) - z_nd(i+1,j+1,k  ) );

            Real h_zeta_cc_xface_lo = 0.5 * dzi *
              (  z_nd(i  ,j  ,k+1) + z_nd(i  ,j+1,k+1)
                -z_nd(i  ,j  ,k  ) - z_nd(i  ,j+1,k  ) );

            Real h_zeta_cc_yface_hi = 0.5 * dzi *
              (  z_nd(i  ,j+1,k+1) + z_nd(i+1,j+1,k+1)
                -z_nd(i  ,j+1,k  ) - z_nd(i+1,j+1,k  ) );

            Real h_zeta_cc_yface_lo = 0.5 * dzi *
              (  z_nd(i  ,j  ,k+1) + z_nd(i+1,j  ,k+1)
                -z_nd(i  ,j  ,k  ) - z_nd(i+1,j  ,k  ) );

            Real xflux_lo = new_drho_u(i  ,j,k)*h_zeta_cc_xface_lo / mf_uy(i  ,j,0);
            Real xflux_hi = new_drho_u(i+1,j,k)*h_zeta_cc_xface_hi / mf_uy(i+1,j,0);
            Real yflux_lo = new_drho_v(i,j  ,k)*h_zeta_cc_yface_lo / mf_vx(i,j  ,0);
            Real yflux_hi = new_drho_v(i,j+1,k)*h_zeta_cc_yface_hi / mf_vx(i,j+1,0);

            Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

            // NOTE: we are saving the (1/J) weighting for later when we add this to rho and theta
            temp_rhs_arr(i,j,k,0) =  ( xflux_hi - xflux_lo ) * dxi * mfsq +
                                     ( yflux_hi - yflux_lo ) * dyi * mfsq;
            temp_rhs_arr(i,j,k,1) = (( xflux_hi * (prim(i,j,k,0) + prim(i+1,j,k,0)) -
                                       xflux_lo * (prim(i,j,k,0) + prim(i-1,j,k,0)) ) * dxi * mfsq+
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
        } // end profile

        // *********************************************************************
        {
        Box gbxo = mfi.nodaltilebox(2);
        Box gbxo_mid = gbxo;

        if (gbxo.smallEnd(2) == domlo.z) {
            Box gbxo_lo = gbxo; gbxo_lo.setBig(2,gbxo.smallEnd(2));
            gbxo_mid.setSmall(2,gbxo.smallEnd(2)+1);
            ParallelFor(gbxo_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                omega_arr(i,j,k) = 0.;
            });
        }
        if (gbxo.bigEnd(2) == domhi.z+1) {
            Box gbxo_hi = gbxo; gbxo_hi.setSmall(2,gbxo.bigEnd(2));
            gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
            ParallelFor(gbxo_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                omega_arr(i,j,k) = old_drho_w(i,j,k);
            });
        }
        ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            omega_arr(i,j,k) = OmegaFromW(i,j,k,old_drho_w(i,j,k),
                                          old_drho_u,old_drho_v,
                                          mf_ux,mf_vy,z_nd,dxInv);
        });
        } // end profile
        // *********************************************************************

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
            Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i,j,k-1)) : 0.0;

            Real coeff_P = coeffP_a(i,j,k) / (1.0 + q);
            Real coeff_Q = coeffQ_a(i,j,k) / (1.0 + q);

            Real theta_t_lo  = 0.5 * ( prim(i,j,k-2,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
            Real theta_t_mid = 0.5 * ( prim(i,j,k-1,PrimTheta_comp) + prim(i,j,k  ,PrimTheta_comp) );
            Real theta_t_hi  = 0.5 * ( prim(i,j,k  ,PrimTheta_comp) + prim(i,j,k+1,PrimTheta_comp) );

            // line 2 last two terms (order dtau)
            Real R0_tmp  =  -halfg * old_drho(i,j,k  ) + coeff_P * old_drho_theta(i,j,k  )
                            -halfg * old_drho(i,j,k-1) + coeff_Q * old_drho_theta(i,j,k-1);

            // line 3 residuals (order dtau^2) 1.0 <-> beta_2
            Real R1_tmp =  -halfg * (  slow_rhs_cons(i,j,k  ,Rho_comp) + slow_rhs_cons(i,j,k-1,Rho_comp) )
                           + coeff_P * slow_rhs_cons(i,j,k  ,RhoTheta_comp)
                           + coeff_Q * slow_rhs_cons(i,j,k-1,RhoTheta_comp);

            Real Omega_kp1 = omega_arr(i,j,k+1);
            Real Omega_k   = omega_arr(i,j,k  );
            Real Omega_km1 = omega_arr(i,j,k-1);

            Real detJdiff = (detJ(i,j,k) - detJ(i,j,k-1)) / (detJ(i,j,k)*detJ(i,j,k-1));

            // consolidate lines 4&5 (order dtau^2)
            R1_tmp += halfg * ( beta_1 * dzi * (Omega_kp1/detJ(i,j,k) + detJdiff*Omega_k - Omega_km1/detJ(i,j,k-1))
                              + temp_rhs_arr(i,j,k,Rho_comp)/detJ(i,j,k) + temp_rhs_arr(i,j,k-1,Rho_comp)/detJ(i,j,k-1) );

            // consolidate lines 6&7 (order dtau^2)
            R1_tmp += -( coeff_P/detJ(i,j,k  ) * ( beta_1 * dzi * (Omega_kp1*theta_t_hi - Omega_k*theta_t_mid) + temp_rhs_arr(i,j,k  ,RhoTheta_comp) )
                       + coeff_Q/detJ(i,j,k-1) * ( beta_1 * dzi * (Omega_k*theta_t_mid - Omega_km1*theta_t_lo) + temp_rhs_arr(i,j,k-1,RhoTheta_comp) ) );

            // line 1
            RHS_a(i,j,k) = old_drho_w(i,j,k)
                         + dtau * (slow_rhs_rho_w(i,j,k) + zmom_src_arr(i,j,k) + R0_tmp + dtau*beta_2*R1_tmp);

            // We cannot use omega_arr here since that was built with old_rho_u and old_rho_v ...
            RHS_a(i,j,k) += OmegaFromW(i,j,k,0.,
                                       new_drho_u,new_drho_v,
                                       mf_ux,mf_vy,z_nd,dxInv);
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
            // w_klo, w_khi given by specified Dirichlet values
            RHS_a(i,j,lo.z  ) = dtau * (slow_rhs_rho_w(i,j,lo.z) + zmom_src_arr(i,j,lo.z));
            RHS_a(i,j,hi.z+1) = dtau * (slow_rhs_rho_w(i,j,hi.z+1) + zmom_src_arr(i,j,hi.z+1));

            // w = specified Dirichlet value at k = lo.z
            soln_a(i,j,lo.z) = RHS_a(i,j,lo.z) * inv_coeffB_a(i,j,lo.z);

            for (int k = lo.z+1; k <= hi.z+1; k++) {
                soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
            }

            cur_zmom(i,j,lo.z  ) = stage_zmom(i,j,lo.z  ) + soln_a(i,j,lo.z  );
            cur_zmom(i,j,hi.z+1) = stage_zmom(i,j,hi.z+1) + soln_a(i,j,hi.z+1);

            for (int k = hi.z; k >= lo.z; k--) {
                soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) *soln_a(i,j,k+1);
            }
        });
#else
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i)
            {
                RHS_a(i,j,lo.z) = dtau * (slow_rhs_rho_w(i,j,lo.z) + zmom_src_arr(i,j,lo.z));
               soln_a(i,j,lo.z) = RHS_a(i,j,lo.z) * inv_coeffB_a(i,j,lo.z);
            }

            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i)
            {
                RHS_a(i,j,hi.z+1) = dtau * (slow_rhs_rho_w(i,j,hi.z+1) + zmom_src_arr(i,j,hi.z+1));
               soln_a(i,j,hi.z+1) = RHS_a(i,j,hi.z+1) * inv_coeffB_a(i,j,hi.z+1);
            }
        }

        for (int k = lo.z+1; k <= hi.z; ++k) {
             for (int j = lo.y; j <= hi.y; ++j) {
                 AMREX_PRAGMA_SIMD
                 for (int i = lo.x; i <= hi.x; ++i) {
                     soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
                 }
           }
        }
        for (int k = hi.z; k > lo.z; --k) {
             for (int j = lo.y; j <= hi.y; ++j) {
                 AMREX_PRAGMA_SIMD
                 for (int i = lo.x; i <= hi.x; ++i) {
                     soln_a(i,j,k) -= (coeffC_a(i,j,k) * inv_coeffB_a(i,j,k)) * soln_a(i,j,k+1);
                 }
             }
        }
        if (hi.z == domhi.z) {
            for (int j = lo.y; j <= hi.y; ++j) {
                 AMREX_PRAGMA_SIMD
                 for (int i = lo.x; i <= hi.x; ++i) {
                    cur_zmom(i,j,hi.z+1) = stage_zmom(i,j,hi.z+1) + soln_a(i,j,hi.z+1);
                }
            }
        }
#endif
        } // end profile

        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cur_zmom(i,j,k) = stage_zmom(i,j,k);
        });

        if (lo.z == domlo.z) {
            tbz.setSmall(2,domlo.z+1);
        }
        if (hi.z == domhi.z) {
            tbz.setBig(2,domhi.z);
        }
        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real wpp = WFromOmega(i,j,k,soln_a(i,j,k),
                                  new_drho_u,new_drho_v,
                                  mf_ux,mf_vy,z_nd,dxInv);

            cur_zmom(i,j,k) += wpp;

            if (l_rayleigh_impl_for_w) {
              Real damping_coeff = l_damp_coef * dtau * sinesq_stag_d[k];
              cur_zmom(i,j,k) /= (1.0 + damping_coeff);
            }
        });

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
                  (flx_arr[2])(i,j,k,0) =    zflux_lo / (mf_mx(i,j,0) * mf_my(i,j,0));
              }

              if (k == vbx_hi.z) {
                  avg_zmom_arr(i,j,k+1)      += facinv * zflux_hi / (mf_mx(i,j,0) * mf_my(i,j,0));
                  if (l_reflux) {
                      (flx_arr[2])(i,j,k+1,0) =      zflux_hi / (mf_mx(i,j,0) * mf_my(i,j,0));
                      (flx_arr[2])(i,j,k+1,1) = (flx_arr[2])(i,j,k+1,0) * 0.5 * (prim(i,j,k) + prim(i,j,k+1));
                  }
              }

              Real fast_rhs_rho = -(temp_rhs_arr(i,j,k,0) + ( zflux_hi - zflux_lo ) * dzi) / detJ(i,j,k);
              cur_cons(i,j,k,0) += dtau * (slow_rhs_cons(i,j,k,0) + fast_rhs_rho);

              Real fast_rhs_rhotheta = -( temp_rhs_arr(i,j,k,1) + 0.5 *
                ( zflux_hi * (prim(i,j,k) + prim(i,j,k+1)) -
                  zflux_lo * (prim(i,j,k) + prim(i,j,k-1)) ) * dzi ) / detJ(i,j,k);

              cur_cons(i,j,k,1) += dtau * (slow_rhs_cons(i,j,k,1) + fast_rhs_rhotheta);

              if (l_reflux) {
                  (flx_arr[2])(i,j,k,1) = (flx_arr[2])(i,j,k,0) * 0.5 * (prim(i,j,k) + prim(i,j,k-1));
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
