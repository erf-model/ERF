#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <ERF_Constants.H>
#include <SpatialStencils.H>
#include <TimeIntegration.H>
#include <EOS.H>

using namespace amrex;

void erf_implicit_fast_rhs (int level,
                   Vector<MultiFab >& S_rhs,                        // the fast RHS we will return
                   Vector<MultiFab >& S_slow_rhs,                   // the slow RHS already computed
                   Vector<MultiFab >& S_stage_data,                 // S_bar = S^n, S^* or S^**
                   const MultiFab& S_stage_prim,
                   const Vector<MultiFab >& S_data,                 // S_sum = most recent full solution
                   const Vector<MultiFab >& S_data_old,             // S_sum_old at most recent fast timestep
                   std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
                   const amrex::Geometry geom,
                         amrex::InterpFaceRegister* ifr,
                   const SolverChoice& solverChoice,
#ifdef ERF_USE_TERRAIN
                   const MultiFab& z_phys_nd,
                   const MultiFab& detJ_cc,
                   const MultiFab& r0,
                   const MultiFab& p0,
#else
                   const amrex::Real* dptr_dens_hse, const amrex::Real* dptr_pres_hse,
#endif
                   const amrex::Real /*time*/, const amrex::Real dtau)
{
    BL_PROFILE_VAR("erf_fast_rhs()",erf_fast_rhs);

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

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    Real dxi = dxInv[0];
    Real dyi = dxInv[1];
    Real dzi = dxInv[2];
    const auto& ba = S_stage_data[IntVar::cons].boxArray();
    const auto& dm = S_stage_data[IntVar::cons].DistributionMap();

    MultiFab Delta_rho_u(    convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab Delta_rho_v(    convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab Delta_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, 1);
    MultiFab Delta_rho  (            ba                , dm, 1, 1);
    MultiFab Delta_rho_theta(        ba                , dm, 1, 1);

    // Create old_drho_u/v/w/theta  = U'', V'', W'', Theta'' in the docs
    MultiFab::Copy(Delta_rho_u    , S_data[IntVar::xmom], 0, 0, 1, 1);
    MultiFab::Copy(Delta_rho_v    , S_data[IntVar::ymom], 0, 0, 1, 1);
    MultiFab::Copy(Delta_rho_w    , S_data[IntVar::zmom], 0, 0, 1, 1);
    MultiFab::Copy(Delta_rho      , S_data[IntVar::cons], Rho_comp     , 0, 1,1);
    MultiFab::Copy(Delta_rho_theta, S_data[IntVar::cons], RhoTheta_comp, 0, 1,1);

    MultiFab::Subtract(Delta_rho_u    , S_stage_data[IntVar::xmom], 0, 0, 1, 1);
    MultiFab::Subtract(Delta_rho_v    , S_stage_data[IntVar::ymom], 0, 0, 1, 1);
    MultiFab::Subtract(Delta_rho_w    , S_stage_data[IntVar::zmom], 0, 0, 1, 1);
    MultiFab::Subtract(Delta_rho      , S_stage_data[IntVar::cons], Rho_comp     , 0, 1, 1);
    MultiFab::Subtract(Delta_rho_theta, S_stage_data[IntVar::cons], RhoTheta_comp, 0, 1, 1);

    // Not sure if we need these
    Delta_rho_u.FillBoundary(geom.periodicity());
    Delta_rho_v.FillBoundary(geom.periodicity());
    Delta_rho_w.FillBoundary(geom.periodicity());
    Delta_rho.FillBoundary(geom.periodicity());
    Delta_rho_theta.FillBoundary(geom.periodicity());

    amrex::MultiFab New_rho_u(    convert(ba,IntVect(1,0,0)), dm, 1, 1);
    amrex::MultiFab New_rho_v(    convert(ba,IntVect(0,1,0)), dm, 1, 1);
    amrex::MultiFab New_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, 1);

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const iMultiFab *mlo_mf_x, *mhi_mf_x;
    const iMultiFab *mlo_mf_y, *mhi_mf_y;
    const iMultiFab *mlo_mf_z, *mhi_mf_z;

    if (level > 0)
    {
        mlo_mf_x = &(ifr->mask(Orientation(0,Orientation::low)));
        mhi_mf_x = &(ifr->mask(Orientation(0,Orientation::high)));
        mlo_mf_y = &(ifr->mask(Orientation(1,Orientation::low)));
        mhi_mf_y = &(ifr->mask(Orientation(1,Orientation::high)));
        mlo_mf_z = &(ifr->mask(Orientation(2,Orientation::low)));
        mhi_mf_z = &(ifr->mask(Orientation(2,Orientation::high)));
    }

    // *************************************************************************
    // Define updates in the current RK stage, fluxes are computed here itself
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_stage_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        Box const& valid_bx = mfi.validbox();
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
        auto mlo_z = (level > 0) ? mlo_mf_z->const_array(mfi) : Array4<const int>{};
        auto mhi_z = (level > 0) ? mhi_mf_z->const_array(mfi) : Array4<const int>{};

        const Array4<      Real> & fast_rhs_cell = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_stage    = S_stage_data[IntVar::cons].const_array(mfi);
        const Array4<const Real> & theta         = S_stage_prim.const_array(mfi);

        const Array4<Real>& old_drho_u     = Delta_rho_u.array(mfi);
        const Array4<Real>& old_drho_v     = Delta_rho_v.array(mfi);
        const Array4<Real>& old_drho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& old_drho       = Delta_rho.array(mfi);
        const Array4<Real>& old_drho_theta = Delta_rho_theta.array(mfi);

        const Array4<Real>& fast_rhs_rho_u = S_rhs[IntVar::xmom].array(mfi);
        const Array4<Real>& fast_rhs_rho_v = S_rhs[IntVar::ymom].array(mfi);
        const Array4<Real>& fast_rhs_rho_w = S_rhs[IntVar::zmom].array(mfi);

        const Array4<const Real>& slow_rhs_cons     = S_slow_rhs[IntVar::cons].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_u    = S_slow_rhs[IntVar::xmom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_v    = S_slow_rhs[IntVar::ymom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_w    = S_slow_rhs[IntVar::zmom].const_array(mfi);

        const Array4<Real>& new_drho_u = New_rho_u.array(mfi);
        const Array4<Real>& new_drho_v = New_rho_v.array(mfi);
        const Array4<Real>& new_drho_w = New_rho_w.array(mfi);

        const Array4<Real>& xflux_rhs = S_rhs[IntVar::xflux].array(mfi);
        const Array4<Real>& yflux_rhs = S_rhs[IntVar::yflux].array(mfi);
        const Array4<Real>& zflux_rhs = S_rhs[IntVar::zflux].array(mfi);

        // These are temporaries we use to add to the S_rhs for the fluxes
        const Array4<Real>& advflux_x = advflux[0].array(mfi);
        const Array4<Real>& advflux_y = advflux[1].array(mfi);
        const Array4<Real>& advflux_z = advflux[2].array(mfi);

        const Array4<const Real>& cur_data = S_data[IntVar::cons].const_array(mfi);
        const Array4<const Real>& old_data = S_data_old[IntVar::cons].const_array(mfi);

#ifdef ERF_USE_TERRAIN
        const Array4<const Real>& z_nd   = z_phys_nd.const_array(mfi);
        const Array4<const Real>& detJ   = detJ_cc.const_array(mfi);
        const Array4<const Real>& r0_arr = r0.const_array(mfi);
        const Array4<const Real>& p0_arr = p0.const_array(mfi);
#endif

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // x-momentum equation

            fast_rhs_rho_u(i, j, k) = 0.0; // Initialize the updated x-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real pi_l = getExnergivenRTh(cell_stage(i-1,j,k,RhoTheta_comp));
                Real pi_r = getExnergivenRTh(cell_stage(i  ,j,k,RhoTheta_comp));
                Real pi_c =  0.5 * (pi_l + pi_r);

                Real drho_theta_hi = old_drho_theta(i,j,k) + beta_d * (
                  (cur_data(i,j  ,k,RhoTheta_comp) - old_data(i,j  ,k,RhoTheta_comp)));
                Real drho_theta_lo = old_drho_theta(i-1,j,k) + beta_d * (
                  (cur_data(i-1,j,k,RhoTheta_comp) - old_data(i-1,j,k,RhoTheta_comp)));

#ifdef ERF_USE_TERRAIN
                Real gp_xi = (drho_theta_hi - drho_theta_lo) * dxi;
                Real h_xi_on_iface = 0.125 * (
                    z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
                   -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k-1) );
                Real h_zeta_on_iface = 0.5 * (
                    z_nd(i,j,k+1) + z_nd(i,j+1,k+1) - z_nd(i,j,k) + z_nd(i,j+1,k) );
                Real gp_zeta_on_iface = 0.25 * dzi * (
                  (old_drho_theta(i  ,j,k+1) + beta_d*((cur_data(i  ,j,k+1,RhoTheta_comp)-old_data(i  ,j,k+1,RhoTheta_comp))))
                 +(old_drho_theta(i-1,j,k+1) + beta_d*((cur_data(i-1,j,k+1,RhoTheta_comp)-old_data(i-1,j,k+1,RhoTheta_comp))))
                 -(old_drho_theta(i  ,j,k-1) + beta_d*((cur_data(i  ,j,k-1,RhoTheta_comp)-old_data(i  ,j,k-1,RhoTheta_comp))))
                 -(old_drho_theta(i-1,j,k-1) + beta_d*((cur_data(i-1,j,k-1,RhoTheta_comp)-old_data(i-1,j,k-1,RhoTheta_comp)))) );
                Real gpx = gp_xi - (h_xi_on_iface / h_zeta_on_iface) * gp_zeta_on_iface;
#else
                Real gpx = (drho_theta_hi - drho_theta_lo)*dxi;
#endif
                fast_rhs_rho_u(i, j, k) = -Gamma * R_d * pi_c * gpx;

                new_drho_u(i, j, k) = old_drho_u(i,j,k) + dtau * fast_rhs_rho_u(i,j,k)
                                                        + dtau * slow_rhs_rho_u(i,j,k);
            } // not on coarse-fine boundary
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // y-momentum equation

            fast_rhs_rho_v(i, j, k) = 0.0; // Initialize the updated y-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (j == vlo_y && mlo_y(i,j,k)) || (j == vhi_y+1 && mhi_y(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real pi_l = getExnergivenRTh(cell_stage(i,j-1,k,RhoTheta_comp));
                Real pi_r = getExnergivenRTh(cell_stage(i,j  ,k,RhoTheta_comp));
                Real pi_c =  0.5 * (pi_l + pi_r);

                Real drho_theta_hi = old_drho_theta(i,j,k) + beta_d * (
                  (cur_data(i,j  ,k,RhoTheta_comp) - old_data(i,j  ,k,RhoTheta_comp)));
                Real drho_theta_lo = old_drho_theta(i,j-1,k) + beta_d * (
                  (cur_data(i,j-1,k,RhoTheta_comp) - old_data(i,j-1,k,RhoTheta_comp)));

#ifdef ERF_USE_TERRAIN
                Real gp_eta = (drho_theta_hi - drho_theta_lo) * dyi;
                Real h_eta_on_jface = 0.125 * (
                    z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
                   -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k-1) );
                Real h_zeta_on_jface = 0.5 * (
                    z_nd(i,j,k+1) + z_nd(i+1,j,k+1) - z_nd(i,j,k) + z_nd(i+1,j,k) );
                Real gp_zeta_on_jface = 0.25 * dzi * (
                  (old_drho_theta(i,j  ,k+1) + beta_d*((cur_data(i,j  ,k+1,RhoTheta_comp)-old_data(i,j  ,k+1,RhoTheta_comp))))
                 +(old_drho_theta(i,j-1,k+1) + beta_d*((cur_data(i,j-1,k+1,RhoTheta_comp)-old_data(i,j-1,k+1,RhoTheta_comp))))
                 -(old_drho_theta(i,j  ,k-1) + beta_d*((cur_data(i,j  ,k-1,RhoTheta_comp)-old_data(i,j  ,k-1,RhoTheta_comp))))
                 -(old_drho_theta(i,j-1,k-1) + beta_d*((cur_data(i,j-1,k-1,RhoTheta_comp)-old_data(i,j-1,k-1,RhoTheta_comp)))) );
                Real gpy = gp_eta - (h_eta_on_jface / h_zeta_on_jface) * gp_zeta_on_jface;
#else
                Real gpy = (drho_theta_hi - drho_theta_lo)*dyi;
#endif
                fast_rhs_rho_v(i, j, k) = -Gamma * R_d * pi_c * gpy;

                new_drho_v(i, j, k) = old_drho_v(i,j,k) + dtau * fast_rhs_rho_v(i,j,k)
                                                        + dtau * slow_rhs_rho_v(i,j,k);
            } // not on coarse-fine boundary
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation

            fast_rhs_rho_w(i, j, k) = 0.0; // Initialize the updated z-mom eqn term to zero

        });

        // *********************************************************************
        // *********************************************************************
        // *********************************************************************

        int klen = bx.bigEnd(2) - bx.smallEnd(2) + 1;
        if (klen > 256) amrex::Abort("Tridiagonal solver is not big enough!");

        amrex::Box b2d = tbz; // Copy constructor
        int klo = tbz.smallEnd(2);
        int khi = tbz.bigEnd(2);
        b2d.setRange(2,0);
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            Array1D<Real,0,255> coeff_A;
            Array1D<Real,0,255> coeff_B;
            Array1D<Real,0,255> coeff_C;
            Array1D<Real,0,255> RHS;
            Array1D<Real,0,255> soln;
            Array1D<Real,0,255> gam;

            //Note we don't act on the bottom or top boundaries of the domain
            for (int k = klo+1; k < khi; ++k)
            {
#ifdef ERF_USE_TERRAIN
                Real rhobar_lo = r0_arr(i,j,k-1);
                Real rhobar_hi = r0_arr(i,j,k  );
                Real  pibar_lo = getExnergivenRTh(p0_arr(i,j,k-1));
                Real  pibar_hi = getExnergivenRTh(p0_arr(i,j,k  ));
#else
                Real rhobar_lo = dptr_dens_hse[k-1];
                Real rhobar_hi = dptr_dens_hse[k];
                Real  pibar_lo = getExnergivenRTh(dptr_pres_hse[k-1]);
                Real  pibar_hi = getExnergivenRTh(dptr_pres_hse[k  ]);
#endif

                // Note that the notes use "g" to mean the magnitude of gravity, so it is positive
                // We set grav_gpu[2] to be the vector component which is negative
                // We define halfg to match the notes (which is why we take the absolute value)
                Real halfg = std::abs(0.5 * grav_gpu[2]);

                Real pi_lo = getExnergivenRTh(cell_stage(i,j,k-1,RhoTheta_comp));
                Real pi_hi = getExnergivenRTh(cell_stage(i,j,k  ,RhoTheta_comp));
                Real pi_c =  0.5 * (pi_lo + pi_hi);

                Real coeff_P = 0.5   * beta_2 * R_d * ( -Gamma * dzi* pi_c )
                             + halfg * beta_2 * R_d  * rhobar_hi * pi_hi  /
                              (c_v * pibar_hi * cell_stage(i,j,k,RhoTheta_comp));

                Real coeff_Q = 0.5 * beta_2 * R_d * (  Gamma * dzi* pi_c )
                             + halfg * beta_2 * R_d * rhobar_lo * pi_lo  /
                              (c_v * pibar_lo * cell_stage(i,j,k-1,RhoTheta_comp));

                Real theta_t_lo  = 0.5 * ( theta(i,j,k-2) + theta(i,j,k-1) );
                Real theta_t_mid = 0.5 * ( theta(i,j,k-1) + theta(i,j,k  ) );
                Real theta_t_hi  = 0.5 * ( theta(i,j,k  ) + theta(i,j,k+1) );

                coeff_A(k) = (dtau * dtau * dzi) * (-beta_2 * coeff_Q *theta_t_lo + halfg*beta_2*beta_2);

                coeff_C(k) = (dtau * dtau * dzi) * ( beta_2 * coeff_P *theta_t_hi - halfg*beta_2*beta_2);

                coeff_B(k) = 1.0 + (beta_2 * dtau * dtau * dzi) * (coeff_Q - coeff_P) * theta_t_mid;

                amrex::Real R_tmp = 0.;

                // first and last parts of lines 2-3
                R_tmp += coeff_P * (old_drho_theta(i,j,k  ) + dtau * slow_rhs_cons(i,j,k  ,RhoTheta_comp))
                        +coeff_Q * (old_drho_theta(i,j,k-1) + dtau * slow_rhs_cons(i,j,k-1,RhoTheta_comp));

                // first part of second half of lines 2-3
                R_tmp += dtau * (
                  -coeff_P * dxi * ( new_drho_u(i+1,j,k  ) * 0.5 * ( theta(i+1,j  ,k  ) + theta(i,j,k) )
                                    -new_drho_u(i  ,j,k  ) * 0.5 * ( theta(i-1,j  ,k  ) + theta(i,j,k) ) )
                  -coeff_P * dyi * ( new_drho_v(i,j+1,k  ) * 0.5 * ( theta(i  ,j+1,k  ) + theta(i,j,k) )
                                    -new_drho_v(i,j  ,k  ) * 0.5 * ( theta(i  ,j-1,k  ) + theta(i,j,k) ) )
                  -coeff_Q * dxi * ( new_drho_u(i+1,j,k-1) * 0.5 * ( theta(i+1,j  ,k-1) + theta(i,j,k-1) )
                                    -new_drho_u(i  ,j,k-1) * 0.5 * ( theta(i-1,j  ,k-1) + theta(i,j,k-1) ) )
                  -coeff_Q * dyi * ( new_drho_v(i,j+1,k-1) * 0.5 * ( theta(i  ,j+1,k-1) + theta(i,j,k-1) )
                                    -new_drho_v(i,j  ,k-1) * 0.5 * ( theta(i  ,j-1,k-1) + theta(i,j,k-1) ) ) );

                // second part of second half of lines 2-3
                R_tmp += -dtau * beta_1 * dzi * (
                   coeff_P * ( old_drho_w(i,j,k+1) * theta_t_hi  - old_drho_w(i,j,k  ) * theta_t_mid )
                  +coeff_Q * ( old_drho_w(i,j,k  ) * theta_t_mid - old_drho_w(i,j,k-1) * theta_t_lo ) );

                // line 4
                R_tmp +=  -beta_1 * Gamma * R_d * pi_c *
                          (old_drho_theta(i,j,k) - old_drho_theta(i,j,k-1))*dzi;

                // line 5
                R_tmp +=  beta_1 * halfg * (R_d / c_v) * (
                    pi_lo * rhobar_lo / pibar_lo * (old_drho_theta(i,j,k-1) / cell_stage(i,j,k-1,RhoTheta_comp))
                   +pi_hi * rhobar_hi / pibar_hi * (old_drho_theta(i,j,k  ) / cell_stage(i,j,k  ,RhoTheta_comp)) );

                // line 6
                R_tmp += -(beta_2+beta_1) * halfg * (
                         ( old_drho(i,j,k  ) + old_drho(i,j,k-1) )
                         + dtau * (slow_rhs_cons(i,j,k,Rho_comp) + slow_rhs_cons(i,j,k-1,Rho_comp)) );

                // lines 7-8
                R_tmp += dtau * beta_2 * halfg * (
                      dxi * (new_drho_u(i+1,j,k  ) - new_drho_u(i,j,k  ))
                    + dxi * (new_drho_u(i+1,j,k-1) - new_drho_u(i,j,k-1))
                    + dyi * (new_drho_v(i,j+1,k  ) - new_drho_v(i,j,k  ))
                    + dyi * (new_drho_v(i,j+1,k-1) - new_drho_v(i,j,k-1))
                    + dzi * (old_drho_w(i,j  ,k+1) - old_drho_w(i,j,k  )) * beta_1
                    + dzi * (old_drho_w(i,j  ,k  ) - old_drho_w(i,j,k-1)) * beta_1 );

                // line 1
                RHS(k) = old_drho_w(i,j,k) + dtau * (slow_rhs_rho_w(i,j,k) + R_tmp);
          } // k

          // w_0 = 0
          coeff_A(0) =  0.0;
          coeff_B(0) =  1.0;
          coeff_C(0) =  0.0;
          RHS(0)     =  0.0;

          // w_khi = 0
          coeff_A(klen) =  0.0;
          coeff_B(klen) =  1.0;
          coeff_C(klen) =  0.0;
          RHS(klen)     =  0.0;

          // w = 0 at k = 0
          Real bet = coeff_B(0);
          soln(0) = RHS(0) / bet;

          for (int k = 1; k <= klen; k++) {
              gam(k) = coeff_C(k-1) / bet;
              bet = coeff_B(k) - coeff_A(k)*gam(k);
              if (bet == 0) amrex::Abort(">>>TRIDIAG FAILED");
              soln(k) = (RHS(k)-coeff_A(k)*soln(k-1)) / bet;
          }
          for (int k = klen-1; k >= 0; k--) {
              soln(k) = soln(k) - gam(k+1)*soln(k+1);
          }

          for (int k = 0; k <= klen; k++) {
              new_drho_w(i,j,k) = soln(k);
              fast_rhs_rho_w(i,j,k) = ( new_drho_w(i,j,k) - old_drho_w(i,j,k) - dtau * slow_rhs_rho_w(i,j,k)) / dtau;
          }
        });

        // **************************************************************************
        // Define updates in the RHS of rho and (rho theta)
        // **************************************************************************

        const int l_spatial_order = 2;
        amrex::ParallelFor(bx, S_stage_data[IntVar::cons].nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {

            if (n == Rho_comp)
            {
                // Note that we pass theta but it won't be used here
                fast_rhs_cell(i, j, k, n) = -AdvectionContributionForState(i, j, k, new_drho_u, new_drho_v, new_drho_w,
                                                                           theta, n, advflux_x, advflux_y, advflux_z,
#ifdef ERF_USE_TERRAIN
                                                                           z_nd, detJ,
#endif
                                                                           dxInv, l_spatial_order);
            } else if (n == RhoTheta_comp) {
                fast_rhs_cell(i, j, k, n) = -AdvectionContributionForState(i, j, k, new_drho_u, new_drho_v, new_drho_w,
                                                                           theta, n, advflux_x, advflux_y, advflux_z,
#ifdef ERF_USE_TERRAIN
                                                                           z_nd, detJ,
#endif
                                                                           dxInv, l_spatial_order);
            } else {
                fast_rhs_cell(i, j, k, n) = 0.0;
            }
        });

        // Compute the RHS for the flux terms from this stage --
        //     we do it this way so we don't double count
        amrex::ParallelFor(tbx, S_stage_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == Rho_comp && n == RhoTheta_comp)
                xflux_rhs(i,j,k,n) = advflux_x(i,j,k,n);
            else
                xflux_rhs(i,j,k,n) = 0.0;
        });
        amrex::ParallelFor(tby, S_stage_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == Rho_comp && n == RhoTheta_comp)
                yflux_rhs(i,j,k,n) = advflux_y(i,j,k,n);
            else
                yflux_rhs(i,j,k,n) = 0.0;
        });
        amrex::ParallelFor(tbz, S_stage_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == Rho_comp && n == RhoTheta_comp)
                zflux_rhs(i,j,k,n) = advflux_z(i,j,k,n);
            else
                zflux_rhs(i,j,k,n) = 0.0;
        });

    } // mfi
}
