#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <IndexDefines.H>
#include <SpatialStencils.H>
#include <TimeIntegration.H>
#include <prob_common.H>

using namespace amrex;

void erf_fast_rhs_N (int step, int level, const Real /*time*/,
                     Vector<MultiFab>& S_slow_rhs,                   // the slow RHS already computed
                     const Vector<MultiFab>& S_prev,                 // if step == 0, this is S_old, else the previous solution
                     Vector<MultiFab>& S_stage_data,                 // S_bar = S^n, S^* or S^**
                     const MultiFab& S_stage_prim,
                     const MultiFab& pi_stage,                       // Exner function evaluted at least stage
                     const MultiFab& fast_coeffs,
                     Vector<MultiFab>& S_data,                       // S_sum = most recent full solution
                     Vector<MultiFab>& S_scratch,                    // S_sum_old at most recent fast timestep for (rho theta)
                     const amrex::Geometry geom,
                     amrex::InterpFaceRegister* ifr,
                     const SolverChoice& solverChoice,
                     const MultiFab* /*z_t_pert*/,
                     std::unique_ptr<MultiFab>& /*z_phys_nd*/,
                     std::unique_ptr<MultiFab>& /*detJ_cc*/,
                     const MultiFab* r0, const MultiFab* pi0,
                     const amrex::Real dtau, const amrex::Real facinv)
{
    BL_PROFILE_REGION("erf_fast_rhs_N()");

    bool l_use_terrain  = solverChoice.use_terrain;
    AMREX_ALWAYS_ASSERT(!l_use_terrain);

    // Per p2902 of Klemp-Skamarock-Dudhia-2007
    // beta_s = -1.0 : fully explicit
    // beta_s =  1.0 : fully implicit
    Real beta_s = 0.1;
    Real beta_1 = 0.5 * (1.0 - beta_s);  // multiplies explicit terms
    Real beta_2 = 0.5 * (1.0 + beta_s);  // multiplies implicit terms

    // How much do we project forward the (rho theta) that is used in the horizontal momentum equations
    Real beta_d = 0.1;

    Real c_v = c_p - R_d;

    const Box domain(geom.Domain());
    const int domhi_z = domain.bigEnd()[2];

    const GpuArray<Real, AMREX_SPACEDIM> dx    = geom.CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    Real dxi = dxInv[0];
    Real dyi = dxInv[1];
    Real dzi = dxInv[2];
    const auto& ba = S_stage_data[IntVar::cons].boxArray();
    const auto& dm = S_stage_data[IntVar::cons].DistributionMap();

    MultiFab Delta_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,0));
    MultiFab Delta_rho  (            ba                , dm, 1, 1);
    MultiFab Delta_rho_theta(        ba                , dm, 1, 1);

    MultiFab New_rho_u(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab New_rho_v(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab New_rho_w(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    MultiFab     coeff_A(fast_coeffs, amrex::make_alias, 0, 1);
    MultiFab inv_coeff_B(fast_coeffs, amrex::make_alias, 1, 1);
    MultiFab     coeff_C(fast_coeffs, amrex::make_alias, 2, 1);
    MultiFab     coeff_P(fast_coeffs, amrex::make_alias, 3, 1);
    MultiFab     coeff_Q(fast_coeffs, amrex::make_alias, 4, 1);

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const iMultiFab *mlo_mf_x, *mhi_mf_x;
    const iMultiFab *mlo_mf_y, *mhi_mf_y;

    if (level > 0)
    {
        mlo_mf_x = &(ifr->mask(Orientation(0,Orientation::low)));
        mhi_mf_x = &(ifr->mask(Orientation(0,Orientation::high)));
        mlo_mf_y = &(ifr->mask(Orientation(1,Orientation::low)));
        mhi_mf_y = &(ifr->mask(Orientation(1,Orientation::high)));
    }

    MultiFab extrap(S_data[IntVar::cons].boxArray(),S_data[IntVar::cons].DistributionMap(),1,1);

    // *************************************************************************
    // Define updates in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {

    FArrayBox temp_rhs_fab;
    FArrayBox RHS_fab;

    for ( MFIter mfi(S_stage_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& valid_bx = mfi.validbox();

        const Box& bx = mfi.tilebox();

        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        if (level > 0) {
            int vlo_x = valid_bx.smallEnd(0);
            int vhi_x = valid_bx.bigEnd(0);
            int vlo_y = valid_bx.smallEnd(1);
            int vhi_y = valid_bx.bigEnd(1);
            int vlo_z = valid_bx.smallEnd(2);
            int vhi_z = valid_bx.bigEnd(2);

            auto mlo_x = (level > 0) ? mlo_mf_x->const_array(mfi) : Array4<const int>{};
            auto mhi_x = (level > 0) ? mhi_mf_x->const_array(mfi) : Array4<const int>{};
            auto mlo_y = (level > 0) ? mlo_mf_y->const_array(mfi) : Array4<const int>{};
            auto mhi_y = (level > 0) ? mhi_mf_y->const_array(mfi) : Array4<const int>{};

            // ******************************************************************
            // This assumes that refined regions are always rectangular
            // ******************************************************************
            bool left_edge_dirichlet = ( level > 0 && mlo_x(vlo_x  ,vlo_y  ,vlo_z) );
            bool rght_edge_dirichlet = ( level > 0 && mhi_x(vhi_x+1,vhi_y  ,vhi_z) );
            bool  bot_edge_dirichlet = ( level > 0 && mlo_y(vlo_x  ,vlo_y  ,vlo_z) );
            bool  top_edge_dirichlet = ( level > 0 && mhi_y(vhi_x  ,vhi_y+1,vhi_z) );

            if (left_edge_dirichlet) tbx.growLo(0,-1);
            if (rght_edge_dirichlet) tbx.growHi(0,-1);
            if ( bot_edge_dirichlet) tby.growLo(1,-1);
            if ( top_edge_dirichlet) tby.growHi(1,-1);
        }  // level > 0

        const Array4<const Real> & stage_cons = S_stage_data[IntVar::cons].const_array(mfi);
        const Array4<const Real> & stage_xmom = S_stage_data[IntVar::xmom].const_array(mfi);
        const Array4<const Real> & stage_ymom = S_stage_data[IntVar::ymom].const_array(mfi);
        const Array4<const Real> & stage_zmom = S_stage_data[IntVar::zmom].const_array(mfi);
        const Array4<const Real> & prim       = S_stage_prim.const_array(mfi);

        const Array4<Real>& old_drho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& old_drho       = Delta_rho.array(mfi);
        const Array4<Real>& old_drho_theta = Delta_rho_theta.array(mfi);

        const Array4<const Real>& slow_rhs_cons  = S_slow_rhs[IntVar::cons].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_u = S_slow_rhs[IntVar::xmom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_v = S_slow_rhs[IntVar::ymom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_w = S_slow_rhs[IntVar::zmom].const_array(mfi);

        const Array4<      Real>& new_drho_u = New_rho_u.array(mfi);
        const Array4<      Real>& new_drho_v = New_rho_v.array(mfi);
        const Array4<      Real>& new_drho_w = New_rho_w.array(mfi);

        const Array4<Real>& cur_cons       = S_data[IntVar::cons].array(mfi);
        const Array4<Real>& cur_xmom       = S_data[IntVar::xmom].array(mfi);
        const Array4<Real>& cur_ymom       = S_data[IntVar::ymom].array(mfi);
        const Array4<Real>& cur_zmom       = S_data[IntVar::zmom].array(mfi);

        const Array4<Real>& scratch_rtheta = S_scratch[IntVar::cons].array(mfi);

        const Array4<const Real>& prev_cons = S_prev[IntVar::cons].const_array(mfi);
        const Array4<const Real>& prev_xmom = S_prev[IntVar::xmom].const_array(mfi);
        const Array4<const Real>& prev_ymom = S_prev[IntVar::ymom].const_array(mfi);
        const Array4<const Real>& prev_zmom = S_prev[IntVar::zmom].const_array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<      Real>& avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        const Array4<      Real>& avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        const Array4<      Real>& avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real>& r0_ca       = r0->const_array(mfi);
        const Array4<const Real>& pi0_ca      = pi0->const_array(mfi);
        const Array4<const Real>& pi_stage_ca = pi_stage.const_array(mfi);

        const Array4<Real>& extrap_arr = extrap.array(mfi);

        // Create old_drho_u/v/w/theta  = U'', V'', W'', Theta'' in the docs
        // Note that we do the Copy and Subtract including one ghost cell
        //    so that we don't have to fill ghost cells of the new MultiFabs
        // Initialize New_rho_u/v/w to Delta_rho_u/v/w so that
        // the ghost cells in New_rho_u/v/w will match old_drho_u/v/w
        Box gbx   = mfi.growntilebox(1);
        Box gtbx  = mfi.nodaltilebox(0).grow(1); gtbx.setSmall(2,0);
        Box gtby  = mfi.nodaltilebox(1).grow(1); gtby.setSmall(2,0);
        Box gtbz  = mfi.nodaltilebox(2).grow(IntVect(1,1,0));

        {
        BL_PROFILE("fast_rhs_copies_0");
        if (step == 0) {
            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cur_cons(i,j,k,Rho_comp)            = prev_cons(i,j,k,Rho_comp);
                cur_cons(i,j,k,RhoTheta_comp)       = prev_cons(i,j,k,RhoTheta_comp);
                scratch_rtheta(i,j,k,RhoTheta_comp) = prev_cons(i,j,k,RhoTheta_comp);
            });
        }
        } // end profile

        {
        BL_PROFILE("fast_rhs_copies_1");
        amrex::ParallelFor(gtbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_w(i,j,k) = prev_zmom(i,j,k) - stage_zmom(i,j,k);
        });
        } // end profile

        {
        BL_PROFILE("fast_rhs_copies_2");
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho(i,j,k)       = cur_cons(i,j,k,Rho_comp)      - stage_cons(i,j,k,Rho_comp);
            old_drho_theta(i,j,k) = cur_cons(i,j,k,RhoTheta_comp) - stage_cons(i,j,k,RhoTheta_comp);
            extrap_arr(i,j,k)     = old_drho_theta(i,j,k) + beta_d * (
              (cur_cons(i,j  ,k,RhoTheta_comp) - scratch_rtheta(i,j  ,k,RhoTheta_comp)));
        });
        } // end profile

        RHS_fab.resize(tbz,1);
        auto const& RHS_a  = RHS_fab.array();
        Elixir rCeli       = RHS_fab.elixir();

        temp_rhs_fab.resize(tbz,2);
        auto const& temp_rhs_arr = temp_rhs_fab.array();
        Elixir temp_rhs_eli      = temp_rhs_fab.elixir();

        auto const&     coeffA_a =     coeff_A.array(mfi);
        auto const& inv_coeffB_a = inv_coeff_B.array(mfi);
        auto const&     coeffC_a =     coeff_C.array(mfi);
        auto const&     coeffP_a =     coeff_P.array(mfi);
        auto const&     coeffQ_a =     coeff_Q.array(mfi);

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        {
        BL_PROFILE("fast_rhs_xymom");
        amrex::ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
            Real pi_c =  0.5 * (pi_stage_ca(i-1,j,k,0) + pi_stage_ca(i,j,k,0));

            Real drho_theta_hi = extrap_arr(i  ,j,k);
            Real drho_theta_lo = extrap_arr(i-1,j,k);

            Real gpx = (drho_theta_hi - drho_theta_lo)*dxi;

#ifdef ERF_USE_MOISTURE
            Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i-1,j,k,PrimQv_comp)
                            +prim(i,j,k,PrimQc_comp) + prim(i-1,j,k,PrimQc_comp) );
            gpx /= (1.0 + q);
#endif
            Real fast_rhs_rho_u = -Gamma * R_d * pi_c * gpx;

            new_drho_u(i, j, k) = prev_xmom(i,j,k) - stage_xmom(i,j,k)
                + dtau * fast_rhs_rho_u + dtau * slow_rhs_rho_u(i,j,k);

            avg_xmom(i,j,k) += facinv*new_drho_u(i,j,k);

            cur_xmom(i,j,k) = stage_xmom(i,j,k) + new_drho_u(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
            Real pi_c =  0.5 * (pi_stage_ca(i,j-1,k,0) + pi_stage_ca(i,j,k,0));

            Real drho_theta_hi = extrap_arr(i,j,k);
            Real drho_theta_lo = extrap_arr(i,j-1,k);

            Real gpy = (drho_theta_hi - drho_theta_lo)*dyi;

#ifdef ERF_USE_MOISTURE
            Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i,j-1,k,PrimQv_comp)
                            +prim(i,j,k,PrimQc_comp) + prim(i,j-1,k,PrimQc_comp) );
            gpy /= (1.0 + q);
#endif
            Real fast_rhs_rho_v = -Gamma * R_d * pi_c * gpy;

            new_drho_v(i, j, k) = prev_ymom(i,j,k) - stage_ymom(i,j,k)
                 + dtau * fast_rhs_rho_v + dtau * slow_rhs_rho_v(i,j,k);

            avg_ymom(i,j,k) += facinv*new_drho_v(i,j,k);

            cur_ymom(i,j,k) = stage_ymom(i,j,k) + new_drho_v(i,j,k);
        });
        } // end profile

        // *********************************************************************
        {
        BL_PROFILE("making_rho_rhs");
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_lo = new_drho_u(i  ,j,k);
            Real xflux_hi = new_drho_u(i+1,j,k);
            Real yflux_lo = new_drho_v(i,j  ,k);
            Real yflux_hi = new_drho_v(i,j+1,k);

            temp_rhs_arr(i,j,k,0) =  ( xflux_hi - xflux_lo ) * dxi + ( yflux_hi - yflux_lo ) * dyi;
            temp_rhs_arr(i,j,k,1) = (( xflux_hi * (prim(i,j,k) + prim(i+1,j,k)) -
                                       xflux_lo * (prim(i,j,k) + prim(i-1,j,k)) ) * dxi +
                                     ( yflux_hi * (prim(i,j,k) + prim(i,j+1,k)) -
                                       yflux_lo * (prim(i,j,k) + prim(i,j-1,k)) ) * dyi) * 0.5;
        });
        } // end profile


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
        {
        BL_PROFILE("fast_loop_on_shrunk");
        //Note we don't act on the bottom or top boundaries of the domain
        ParallelFor(bx_shrunk_in_k, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
#if 0
            Real rhobar_lo, rhobar_hi, pibar_lo, pibar_hi;
            rhobar_lo =  r0_ca(i,j,k-1);
            rhobar_hi =  r0_ca(i,j,k  );
             pibar_lo = pi0_ca(i,j,k-1);
             pibar_hi = pi0_ca(i,j,k  );

             Real pi_lo = pi_stage_ca(i,j,k-1,0);
             Real pi_hi = pi_stage_ca(i,j,k  ,0);
             Real pi_c =  0.5 * (pi_lo + pi_hi);

             Real coeff_P = -Gamma * R_d * pi_c * dzi
                          +  halfg * R_d * rhobar_hi * pi_hi  /
                          (  c_v * pibar_hi * stage_cons(i,j,k,RhoTheta_comp) );

             Real coeff_Q = Gamma * R_d * pi_c * dzi
                          + halfg * R_d * rhobar_lo * pi_lo  /
                          ( c_v  * pibar_lo * stage_cons(i,j,k-1,RhoTheta_comp) );
#else
             Real coeff_P = coeffP_a(i,j,k);
             Real coeff_Q = coeffQ_a(i,j,k);
#endif

#ifdef ERF_USE_MOISTURE
            Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i,j,k-1,PrimQv_comp)
                            +prim(i,j,k,PrimQc_comp) + prim(i,j,k-1,PrimQc_comp) );
            coeff_P /= (1.0 + q);
            coeff_Q /= (1.0 + q);
#endif

            Real theta_t_lo  = 0.5 * ( prim(i,j,k-2,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
            Real theta_t_mid = 0.5 * ( prim(i,j,k-1,PrimTheta_comp) + prim(i,j,k  ,PrimTheta_comp) );
            Real theta_t_hi  = 0.5 * ( prim(i,j,k  ,PrimTheta_comp) + prim(i,j,k+1,PrimTheta_comp) );

            // line 2 last two terms (order dtau)
            Real R_tmp = coeff_P * old_drho_theta(i,j,k) + coeff_Q * old_drho_theta(i,j,k-1)
                         - halfg * ( old_drho(i,j,k) + old_drho(i,j,k-1) );

            // line 3 residuals (order dtau^2) 1.0 <-> beta_2
            R_tmp += -dtau * beta_2 * halfg * ( slow_rhs_cons(i,j,k  ,Rho_comp) +
                                                slow_rhs_cons(i,j,k-1,Rho_comp) )
                   +  dtau * beta_2 * ( coeff_P * slow_rhs_cons(i,j,k  ,RhoTheta_comp) +
                                        coeff_Q * slow_rhs_cons(i,j,k-1,RhoTheta_comp) );

            // line 4 (order dtau^2)
            Real Omega_kp1 = prev_zmom(i,j,k+1) - stage_zmom(i,j,k+1);
            Real Omega_k   = prev_zmom(i,j,k  ) - stage_zmom(i,j,k  );
            Real Omega_km1 = prev_zmom(i,j,k-1) - stage_zmom(i,j,k-1);
            R_tmp += ( dtau * beta_2 * halfg ) *
                     ( beta_1 * dzi * (Omega_kp1 - Omega_k) + temp_rhs_arr(i,j,k) );

            // line 6 (reuse Omega & metrics) (order dtau^2)
            R_tmp += -( dtau * beta_2 * coeff_P ) *
                      ( beta_1 * dzi * (Omega_kp1*theta_t_hi - Omega_k*theta_t_mid) + temp_rhs_arr(i,j,k,1) );

            // line 5 (order dtau^2)
            R_tmp += ( dtau * beta_2 * halfg ) *
                     ( beta_1 * dzi * (Omega_k - Omega_km1) + temp_rhs_arr(i,j,k-1) );

            // line 7 (reuse Omega & metrics) (order dtau^2)
            R_tmp += -( dtau * beta_2 * coeff_Q ) *
                      ( beta_1 * dzi * (Omega_k*theta_t_mid - Omega_km1*theta_t_lo) + temp_rhs_arr(i,j,k-1,1) );

            // line 1
            RHS_a(i,j,k) = Omega_k + dtau * (slow_rhs_rho_w(i,j,k) + R_tmp);
        });
        } // end profile

        amrex::Box b2d = tbz; // Copy constructor
        b2d.setRange(2,0);

        {
        BL_PROFILE("fast_rhs_b2d_loop");
#ifdef AMREX_USE_GPU
        auto const hi = amrex::ubound(bx);
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
          // w_0 = 0
          RHS_a   (i,j,0) =  0.0;

          // w_khi = 0
          // Note that if we ever change this, we will need to include it in avg_zmom at the top
          RHS_a   (i,j,hi.z+1) =  0.0;

          // w = 0 at k = 0
          new_drho_w(i,j,0) = RHS_a(i,j,0) * inv_coeffB_a(i,j,0);
          cur_zmom(i,j,0) = stage_zmom(i,j,0) + new_drho_w(i,j,0);

          for (int k = 1; k <= hi.z+1; k++) {
              new_drho_w(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*new_drho_w(i,j,k-1)) * inv_coeffB_a(i,j,k);
          }
          cur_zmom(i,j,hi.z+1) = stage_zmom(i,j,hi.z+1) + new_drho_w(i,j,hi.z+1);
          for (int k = hi.z; k >= 0; k--) {
              new_drho_w(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) *new_drho_w(i,j,k+1);
              cur_zmom(i,j,k) = stage_zmom(i,j,k) + new_drho_w(i,j,k);
          }
        }); // b2d
#else
        auto const lo = amrex::lbound(bx);
        auto const hi = amrex::ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                RHS_a   (i,j,0) =  0.0;
                new_drho_w(i,j,0) = RHS_a(i,j,0) * inv_coeffB_a(i,j,0);
            }
        }
        // Note that if we ever change this, we will need to include it in avg_zmom at the top
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                RHS_a   (i,j,hi.z+1) =  0.0;
            }
        }
        for (int k = lo.z+1; k <= hi.z+1; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    new_drho_w(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*new_drho_w(i,j,k-1)) * inv_coeffB_a(i,j,k);
                }
            }
        }
        int k = hi.z+1;
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    cur_zmom(i,j,k) = stage_zmom(i,j,k) + new_drho_w(i,j,k);
                }
            }
        for (int k = hi.z; k >= lo.z; --k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    new_drho_w(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) * new_drho_w(i,j,k+1);
                    cur_zmom(i,j,k) = stage_zmom(i,j,k) + new_drho_w(i,j,k);
                }
            }
        }
#endif
        } // end profile

        // **************************************************************************
        // Define updates in the RHS of rho and (rho theta)
        // **************************************************************************

        // We note that valid_bx is the actual grid, while bx may be a tile within that grid
        const auto& vbx_hi = amrex::ubound(valid_bx);

        {
        BL_PROFILE("fast_rho_final_update");
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real zflux_lo = beta_2 * new_drho_w(i,j,k  ) + beta_1 * old_drho_w(i,j,k  );
            Real zflux_hi = beta_2 * new_drho_w(i,j,k+1) + beta_1 * old_drho_w(i,j,k+1);

            avg_zmom(i,j,k  ) += facinv*zflux_lo;
            // Note that in the solve we effectively impose new_drho_w(i,j,vbx_hi.z+1)=0
            // so we don't update avg_zmom at k=vbx_hi.z+1

            cur_cons(i,j,k,0) += dtau * (slow_rhs_cons(i,j,k,0) - temp_rhs_arr(i,j,k,0) - ( zflux_hi - zflux_lo ) * dzi );

            cur_cons(i,j,k,1) += dtau * (slow_rhs_cons(i,j,k,1) - temp_rhs_arr(i,j,k,1) - 0.5 * (
              ( zflux_hi * (prim(i,j,k) + prim(i,j,k+1)) - zflux_lo * (prim(i,j,k) + prim(i,j,k-1)) ) * dzi ) );
        });
        } // end profile
    } // mfi
    }
}
