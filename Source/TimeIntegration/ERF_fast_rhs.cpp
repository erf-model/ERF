#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <IndexDefines.H>
#include <SpatialStencils.H>
#include <TerrainMetrics.H>
#include <TimeIntegration.H>

using namespace amrex;

void erf_fast_rhs (int level,
                   Vector<MultiFab>& S_rhs,                        // the fast RHS we will return
                   Vector<MultiFab>& S_slow_rhs,                   // the slow RHS already computed
                   Vector<MultiFab>& S_stage_data,                 // S_bar = S^n, S^* or S^**
                   const MultiFab& S_stage_prim,
                   const Vector<MultiFab>& S_data,                 // S_sum = most recent full solution
                         Vector<MultiFab>& S_scratch,              // S_sum_old at most recent fast timestep for (rho theta)
                   const amrex::Geometry geom,
                   amrex::InterpFaceRegister* ifr,
                   const SolverChoice& solverChoice,
                   std::unique_ptr<MultiFab>& z_phys_nd,
                   std::unique_ptr<MultiFab>& detJ_cc,
                   const MultiFab* r0,
                   const MultiFab* p0,
                   const amrex::Real dtau, const amrex::Real facinv)
{
    BL_PROFILE_VAR("erf_fast_rhs()",erf_fast_rhs);

    bool l_use_terrain = solverChoice.use_terrain;

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
    MultiFab Delta_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,0));
    MultiFab Delta_rho  (            ba                , dm, 1, 1);
    MultiFab Delta_rho_theta(        ba                , dm, 1, 1);

    MultiFab New_rho_u(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab New_rho_v(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab New_rho_w(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    MultiFab Exner_mf(convert(ba,IntVect(1,1,1)), dm, 1, 1);

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
    // Define updates in the current RK stage, fluxes are computed here itself
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {

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

        auto mlo_x = (level > 0) ? mlo_mf_x->const_array(mfi) : Array4<const int>{};
        auto mhi_x = (level > 0) ? mhi_mf_x->const_array(mfi) : Array4<const int>{};
        auto mlo_y = (level > 0) ? mlo_mf_y->const_array(mfi) : Array4<const int>{};
        auto mhi_y = (level > 0) ? mhi_mf_y->const_array(mfi) : Array4<const int>{};

        const Array4<      Real> & fast_rhs_cell = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_stage    = S_stage_data[IntVar::cons].const_array(mfi);
        const Array4<const Real> & cell_stage_xmom  = S_stage_data[IntVar::xmom].const_array(mfi);
        const Array4<const Real> & cell_stage_ymom  = S_stage_data[IntVar::ymom].const_array(mfi);
        const Array4<const Real> & cell_stage_zmom  = S_stage_data[IntVar::zmom].const_array(mfi);
        const Array4<const Real> & prim          = S_stage_prim.const_array(mfi);

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

        const Array4<      Real>& new_drho_u = New_rho_u.array(mfi);
        const Array4<      Real>& new_drho_v = New_rho_v.array(mfi);
        const Array4<      Real>& new_drho_w = New_rho_w.array(mfi);

        const Array4<      Real>& Exner_arr  = Exner_mf.array(mfi);

        const Array4<      Real>& xflux_rhs = S_rhs[IntVar::xflux].array(mfi);
        const Array4<      Real>& yflux_rhs = S_rhs[IntVar::yflux].array(mfi);
        const Array4<      Real>& zflux_rhs = S_rhs[IntVar::zflux].array(mfi);

        const Array4<const Real>& cur_data = S_data[IntVar::cons].const_array(mfi);
        const Array4<const Real>& cur_data_xmom = S_data[IntVar::xmom].const_array(mfi);
        const Array4<const Real>& cur_data_ymom = S_data[IntVar::ymom].const_array(mfi);
        const Array4<const Real>& cur_data_zmom = S_data[IntVar::zmom].const_array(mfi);
        const Array4<const Real>& old_data = S_scratch[IntVar::cons].const_array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<      Real>& avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        const Array4<      Real>& avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        const Array4<      Real>& avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real>& z_nd   = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ   = l_use_terrain ?   detJ_cc->const_array(mfi) : Array4<const Real>{};

        const Array4<const Real>& r0_arr = r0->const_array(mfi);
        const Array4<const Real>& p0_arr = p0->const_array(mfi);

        const Array4<Real>& extrap_arr = extrap.array(mfi);

        // Create old_drho_u/v/w/theta  = U'', V'', W'', Theta'' in the docs
        // Note that we do the Copy and Subtract including one ghost cell
        //    so that we don't have to fill ghost cells of the new MultiFabs
        // Initialize New_rho_u/v/w to Delta_rho_u/v/w so that
        // the ghost cells in New_rho_u/v/w will match old_drho_u/v/w
        const Box gbx   = mfi.growntilebox(1);
        const Box gtbx  = mfi.nodaltilebox(0).grow(1);
        const Box gtby  = mfi.nodaltilebox(1).grow(1);
        const Box gtbz  = mfi.nodaltilebox(2).grow(IntVect(1,1,0));
        const Box gntbx = mfi.nodaltilebox().grow(1);

        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho(i,j,k) = cur_data(i,j,k,Rho_comp) - cell_stage(i,j,k,Rho_comp);
            old_drho_theta(i,j,k) = cur_data(i,j,k,RhoTheta_comp) - cell_stage(i,j,k,RhoTheta_comp);
            extrap_arr(i,j,k) = old_drho_theta(i,j,k) + beta_d * (
              (cur_data(i,j  ,k,RhoTheta_comp) - old_data(i,j  ,k,RhoTheta_comp)));
        });

        amrex::ParallelFor(gtbx, gtby, gtbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_u(i,j,k) = cur_data_xmom(i,j,k) - cell_stage_xmom(i,j,k);
            new_drho_u(i,j,k) = old_drho_u(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_v(i,j,k) = cur_data_ymom(i,j,k) - cell_stage_ymom(i,j,k);
            new_drho_v(i,j,k) = old_drho_v(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho_w(i,j,k) = cur_data_zmom(i,j,k) - cell_stage_zmom(i,j,k);
            new_drho_w(i,j,k) = old_drho_w(i,j,k);
        });

        amrex::ParallelFor(gntbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Exner_arr(i,j,k) = getExnergivenRTh(cell_stage(i,j,k,RhoTheta_comp));
        });

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
                Real pi_l = Exner_arr(i-1,j,k);
                Real pi_r = Exner_arr(i  ,j,k);
                Real pi_c =  0.5 * (pi_l + pi_r);

                Real drho_theta_hi = extrap_arr(i  ,j,k);
                Real drho_theta_lo = extrap_arr(i-1,j,k);

                Real gpx;
                if (l_use_terrain) {
                    Real met_h_xi,met_h_eta,met_h_zeta;
                    ComputeMetricAtIface(i,j,k,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::h_xi_zeta);
                    Real gp_xi = (drho_theta_hi - drho_theta_lo) * dxi;
                    Real gp_zeta_on_iface;
                    if(k==0) {
                        gp_zeta_on_iface = 0.5 * dzi * (
                                                       extrap_arr(i-1,j,k+1) + extrap_arr(i,j,k+1)
                                                     - extrap_arr(i-1,j,k  ) - extrap_arr(i,j,k  ) );
                    } else if(k==domhi_z) {
                        gp_zeta_on_iface = 0.5 * dzi * (
                                                       extrap_arr(i-1,j,k  ) + extrap_arr(i,j,k  )
                                                     - extrap_arr(i-1,j,k-1) - extrap_arr(i,j,k-1) );
                    } else {
                        gp_zeta_on_iface = 0.25 * dzi * (
                                                        extrap_arr(i-1,j,k+1) + extrap_arr(i,j,k+1)
                                                      - extrap_arr(i-1,j,k-1) - extrap_arr(i,j,k-1) );
                    }
                    gpx = gp_xi - (met_h_xi / met_h_zeta) * gp_zeta_on_iface;
                } else {
                    gpx = (drho_theta_hi - drho_theta_lo)*dxi;
                }

#ifdef ERF_USE_MOISTURE
                Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i-1,j,k,PrimQv_comp)
                                +prim(i,j,k,PrimQc_comp) + prim(i-1,j,k,PrimQc_comp) );
                gpx /= (1.0 + q);
#endif
                fast_rhs_rho_u(i, j, k) = -Gamma * R_d * pi_c * gpx;

                new_drho_u(i, j, k) = old_drho_u(i,j,k) + dtau * fast_rhs_rho_u(i,j,k)
                                                        + dtau * slow_rhs_rho_u(i,j,k);

                if (k == domhi_z) new_drho_u(i,j,k+1) = new_drho_u(i,j,k);

                avg_xmom(i,j,k) += facinv * new_drho_u(i,j,k);

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
                Real pi_l = Exner_arr(i,j-1,k);
                Real pi_r = Exner_arr(i,j  ,k);
                Real pi_c =  0.5 * (pi_l + pi_r);

                Real drho_theta_hi = extrap_arr(i,j,k);
                Real drho_theta_lo = extrap_arr(i,j-1,k);

                Real gpy;
                if (l_use_terrain) {
                    Real met_h_xi,met_h_eta,met_h_zeta;
                    ComputeMetricAtJface(i,j,k,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::h_eta_zeta);
                    Real gp_eta = (drho_theta_hi - drho_theta_lo) * dyi;
                    Real gp_zeta_on_jface;
                    if(k==0) {
                        gp_zeta_on_jface = 0.5 * dzi * (
                                                       extrap_arr(i,j,k+1) + extrap_arr(i,j-1,k+1)
                                                     - extrap_arr(i,j,k  ) - extrap_arr(i,j-1,k  ) );
                    } else if(k==domhi_z) {
                        gp_zeta_on_jface = 0.5 * dzi * (
                                                       extrap_arr(i,j,k  ) + extrap_arr(i,j-1,k  )
                                                     - extrap_arr(i,j,k-1) - extrap_arr(i,j-1,k-1) );
                    } else {
                        gp_zeta_on_jface = 0.25 * dzi * (
                                                        extrap_arr(i,j,k+1) + extrap_arr(i,j-1,k+1)
                                                      - extrap_arr(i,j,k-1) - extrap_arr(i,j-1,k-1) );
                    }
                    gpy = gp_eta - (met_h_eta / met_h_zeta) * gp_zeta_on_jface;
               } else {
                    gpy = (drho_theta_hi - drho_theta_lo)*dyi;
               }

#ifdef ERF_USE_MOISTURE
                Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i,j-1,k,PrimQv_comp)
                                +prim(i,j,k,PrimQc_comp) + prim(i,j-1,k,PrimQc_comp) );
                gpy /= (1.0 + q);
#endif
                fast_rhs_rho_v(i, j, k) = -Gamma * R_d * pi_c * gpy;

                new_drho_v(i, j, k) = old_drho_v(i,j,k) + dtau * fast_rhs_rho_v(i,j,k)
                                                        + dtau * slow_rhs_rho_v(i,j,k);
                if (k == domhi_z) new_drho_v(i,j,k+1) = new_drho_v(i,j,k);

                avg_ymom(i,j,k) += facinv * new_drho_v(i,j,k);

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

            Array1D<Real,0,256> coeff_A;
            Array1D<Real,0,256> coeff_B;
            Array1D<Real,0,256> coeff_C;
            Array1D<Real,0,256> RHS;
            Array1D<Real,0,256> soln;
            Array1D<Real,0,256> gam;

            //Note we don't act on the bottom or top boundaries of the domain
            for (int k = klo+1; k < khi; ++k)
            {
                Real rhobar_lo, rhobar_hi, pibar_lo, pibar_hi;
                rhobar_lo = (k == 0) ?  r0_arr(i,j,k) : r0_arr(i,j,k-1);
                rhobar_hi = r0_arr(i,j,k  );
                 pibar_lo = getExnergivenRTh(p0_arr(i,j,k-1));
                 pibar_hi = getExnergivenRTh(p0_arr(i,j,k  ));

                // Note that the notes use "g" to mean the magnitude of gravity, so it is positive
                // We set grav_gpu[2] to be the vector component which is negative
                // We define halfg to match the notes (which is why we take the absolute value)
                Real halfg = std::abs(0.5 * grav_gpu[2]);

                Real pi_lo = Exner_arr(i,j,k-1);
                Real pi_hi = Exner_arr(i,j,k  );
                Real pi_c =  0.5 * (pi_lo + pi_hi);

                Real detJ_on_kface      = 1.0;
                Real h_zeta_on_kface    = 1.0;
                Real h_zeta_cc_xface_hi = 1.0;
                Real h_zeta_cc_xface_lo = 1.0;
                Real h_zeta_cc_yface_hi = 1.0;
                Real h_zeta_cc_yface_lo = 1.0;

                if (l_use_terrain) {
                    h_zeta_on_kface = 0.125 * dzi * (
                        z_nd(i,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k+1)
                       -z_nd(i,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j,k-1) - z_nd(i+1,j-1,k-1) );

                    detJ_on_kface = 0.5 * (detJ(i,j,k) + detJ(i,j,k-1));
                }

                Real coeff_P = -Gamma * R_d * pi_c * dzi / h_zeta_on_kface
                             +  halfg * R_d * rhobar_hi * pi_hi  /
                             (  c_v * pibar_hi * cell_stage(i,j,k,RhoTheta_comp) );

                Real coeff_Q = Gamma * R_d * pi_c * dzi / h_zeta_on_kface
                             + halfg * R_d * rhobar_lo * pi_lo  /
                             ( c_v  * pibar_lo * cell_stage(i,j,k-1,RhoTheta_comp) );

#ifdef ERF_USE_MOISTURE
                Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i,j,k-1,PrimQv_comp)
                                +prim(i,j,k,PrimQc_comp) + prim(i,j,k-1,PrimQc_comp) );
                coeff_P /= (1.0 + q);
                coeff_Q /= (1.0 + q);
#endif

                Real theta_t_lo  = 0.5 * ( prim(i,j,k-2,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
                Real theta_t_mid = 0.5 * ( prim(i,j,k-1,PrimTheta_comp) + prim(i,j,k  ,PrimTheta_comp) );
                Real theta_t_hi  = 0.5 * ( prim(i,j,k  ,PrimTheta_comp) + prim(i,j,k+1,PrimTheta_comp) );

                // LHS for tri-diagonal system
                Real D = dtau * dtau * beta_2 * beta_2 * dzi / detJ_on_kface;
                coeff_A(k) = D * ( halfg - coeff_Q * theta_t_lo );
                coeff_C(k) = D * (-halfg + coeff_P * theta_t_hi );
                coeff_B(k) = 1.0 + D * (coeff_Q - coeff_P) * theta_t_mid;

                amrex::Real R_tmp = 0.;

                // line 2 last two terms (order dtau)
                R_tmp += coeff_P * old_drho_theta(i,j,k) + coeff_Q * old_drho_theta(i,j,k-1)
                       - halfg * ( old_drho(i,j,k) + old_drho(i,j,k-1) );

                // line 3 residuals (order dtau^2) 1.0 <-> beta_2
                R_tmp += -dtau * beta_2 * halfg * ( slow_rhs_cons(i,j,k  ,Rho_comp) +
                                                    slow_rhs_cons(i,j,k-1,Rho_comp) )
                       +  dtau * beta_2 * ( coeff_P * slow_rhs_cons(i,j,k  ,RhoTheta_comp) +
                                            coeff_Q * slow_rhs_cons(i,j,k-1,RhoTheta_comp) );

                // line 4 (order dtau^2)
                Real Omega_hi;
                Real Omega_lo;
                if (l_use_terrain) {
                    Omega_hi = OmegaFromW(i,j,k+1,old_drho_w(i,j,k+1),old_drho_u,old_drho_v,z_nd,dxInv);
                    Omega_lo = OmegaFromW(i,j,k  ,old_drho_w(i,j,k  ),old_drho_u,old_drho_v,z_nd,dxInv);
                    h_zeta_cc_xface_hi = 0.5 * dzi *
                      (  z_nd(i+1,j  ,k+1) + z_nd(i+1,j+1,k+1)
                        -z_nd(i+1,j  ,k  ) - z_nd(i+1,j+1,k  ) );

                    h_zeta_cc_xface_lo = 0.5 * dzi *
                      (  z_nd(i  ,j  ,k+1) + z_nd(i  ,j+1,k+1)
                        -z_nd(i  ,j  ,k  ) - z_nd(i  ,j+1,k  ) );

                    h_zeta_cc_yface_hi = 0.5 * dzi *
                      (  z_nd(i  ,j+1,k+1) + z_nd(i+1,j+1,k+1)
                        -z_nd(i  ,j+1,k  ) - z_nd(i+1,j+1,k  ) );

                    h_zeta_cc_yface_lo = 0.5 * dzi *
                      (  z_nd(i  ,j  ,k+1) + z_nd(i+1,j  ,k+1)
                        -z_nd(i  ,j  ,k  ) - z_nd(i+1,j  ,k  ) );
                } else {
                    Omega_hi = old_drho_w(i,j,k+1);
                    Omega_lo = old_drho_w(i,j,k  );
                }
                R_tmp += ( dtau * beta_2 * halfg / detJ_on_kface ) *
                         ( beta_1 * dzi * (Omega_hi - Omega_lo)    +
                                    dxi * (new_drho_u(i+1,j,k)*h_zeta_cc_xface_hi  -
                                           new_drho_u(i  ,j,k)*h_zeta_cc_xface_lo) +
                                    dyi * (new_drho_v(i,j+1,k)*h_zeta_cc_yface_hi  -
                                           new_drho_v(i,j  ,k)*h_zeta_cc_yface_lo) );

                // line 6 (reuse Omega & metrics) (order dtau^2)
                Real Theta_x_hi = 0.5 * ( prim(i+1,j  ,k,PrimTheta_comp) + prim(i,j,k,PrimTheta_comp) );
                Real Theta_x_lo = 0.5 * ( prim(i-1,j  ,k,PrimTheta_comp) + prim(i,j,k,PrimTheta_comp) );
                Real Theta_y_hi = 0.5 * ( prim(i  ,j+1,k,PrimTheta_comp) + prim(i,j,k,PrimTheta_comp) );
                Real Theta_y_lo = 0.5 * ( prim(i  ,j-1,k,PrimTheta_comp) + prim(i,j,k,PrimTheta_comp) );
                R_tmp += -( dtau * beta_2 * coeff_P / detJ_on_kface ) *
                          ( beta_1 * dzi * (Omega_hi*theta_t_hi - Omega_lo*theta_t_mid) +
                                     dxi * (new_drho_u(i+1,j,k)*Theta_x_hi*h_zeta_cc_xface_hi  -
                                            new_drho_u(i  ,j,k)*Theta_x_lo*h_zeta_cc_xface_lo) +
                                     dyi * (new_drho_v(i,j+1,k)*Theta_y_hi*h_zeta_cc_yface_hi  -
                                            new_drho_v(i,j  ,k)*Theta_y_lo*h_zeta_cc_yface_lo) );

                // line 5 (order dtau^2)
                if (l_use_terrain) {
                    Omega_hi = OmegaFromW(i,j,k  ,old_drho_w(i,j,k  ),old_drho_u,old_drho_v,z_nd,dxInv);
                    Omega_lo = OmegaFromW(i,j,k-1,old_drho_w(i,j,k-1),old_drho_u,old_drho_v,z_nd,dxInv);
                    h_zeta_cc_xface_hi = 0.5 * dzi *
                      (  z_nd(i+1,j  ,k  ) + z_nd(i+1,j+1,k  )
                        -z_nd(i+1,j  ,k-1) - z_nd(i+1,j+1,k-1) );

                    h_zeta_cc_xface_lo = 0.5 * dzi *
                      (  z_nd(i  ,j  ,k  ) + z_nd(i  ,j+1,k  )
                        -z_nd(i  ,j  ,k-1) - z_nd(i  ,j+1,k-1) );

                    h_zeta_cc_yface_hi = 0.5 * dzi *
                      (  z_nd(i  ,j+1,k  ) + z_nd(i+1,j+1,k  )
                        -z_nd(i  ,j+1,k-1) - z_nd(i+1,j+1,k-1) );

                    h_zeta_cc_yface_lo = 0.5 * dzi *
                      (  z_nd(i  ,j  ,k  ) + z_nd(i+1,j  ,k  )
                        -z_nd(i  ,j  ,k-1) - z_nd(i+1,j  ,k-1) );
                } else {
                    Omega_hi = old_drho_w(i,j,k  );
                    Omega_lo = old_drho_w(i,j,k-1);
                }
                R_tmp += ( dtau * beta_2 * halfg / detJ_on_kface ) *
                         ( beta_1 * dzi * (Omega_hi - Omega_lo) +
                                    dxi * (new_drho_u(i+1,j,k-1)*h_zeta_cc_xface_hi  -
                                           new_drho_u(i  ,j,k-1)*h_zeta_cc_xface_lo) +
                                    dyi * (new_drho_v(i,j+1,k-1)*h_zeta_cc_yface_hi  -
                                           new_drho_v(i,j  ,k-1)*h_zeta_cc_yface_lo) );

                // line 7 (reuse Omega & metrics) (order dtau^2)
                Theta_x_hi = 0.5 * ( prim(i+1,j  ,k-1,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
                Theta_x_lo = 0.5 * ( prim(i-1,j  ,k-1,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
                Theta_y_hi = 0.5 * ( prim(i  ,j+1,k-1,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
                Theta_y_lo = 0.5 * ( prim(i  ,j-1,k-1,PrimTheta_comp) + prim(i,j,k-1,PrimTheta_comp) );
                R_tmp += -( dtau * beta_2 * coeff_Q / detJ_on_kface ) *
                          ( beta_1 * dzi * (Omega_hi*theta_t_mid - Omega_lo*theta_t_lo) +
                                     dxi * (new_drho_u(i+1,j,k-1)*Theta_x_hi*h_zeta_cc_xface_hi  -
                                            new_drho_u(i  ,j,k-1)*Theta_x_lo*h_zeta_cc_xface_lo) +
                                     dyi * (new_drho_v(i,j+1,k-1)*Theta_y_hi*h_zeta_cc_yface_hi  -
                                            new_drho_v(i,j  ,k-1)*Theta_y_lo*h_zeta_cc_yface_lo) );

                // line 1
                RHS(k) = old_drho_w(i,j,k) + dtau * (slow_rhs_rho_w(i,j,k) + R_tmp);
                if (l_use_terrain) {
                    RHS(k) += OmegaFromW(i,j,k,0.,new_drho_u,new_drho_v,z_nd,dxInv);
                }
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
              if (l_use_terrain) {
                  new_drho_w(i,j,k) = WFromOmega(i,j,k,soln(k),new_drho_u,new_drho_v,z_nd,dxInv);
              } else {
                  new_drho_w(i,j,k) = soln(k);
              }
              fast_rhs_rho_w(i,j,k) = ( new_drho_w(i,j,k) - old_drho_w(i,j,k) - dtau * slow_rhs_rho_w(i,j,k)) / dtau;

              avg_zmom(i,j,k) += facinv * new_drho_w(i,j,k);
          }
        });

        // **************************************************************************
        // Define updates in the RHS of rho and (rho theta)
        // **************************************************************************

        const int l_spatial_order = 2;

        // We only update (rho) and (rho theta) here
        int start_comp = 0;
        int ncomp      = 2;
        AdvectionSrcForState(bx, start_comp, ncomp, new_drho_u, new_drho_v, new_drho_w,
                             prim, fast_rhs_cell, xflux_rhs, yflux_rhs, zflux_rhs,
                             z_nd, detJ, dxInv, l_spatial_order, l_use_terrain, false, false);
    } // mfi
    }
}
