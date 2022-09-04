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
#include <prob_common.H>

using namespace amrex;

void erf_fast_rhs_T (int step, int level, const Real time,
                     Vector<MultiFab>& S_slow_rhs,                   // the slow RHS already computed
                     const Vector<MultiFab>& S_old,
                     Vector<MultiFab>& S_stage_data,                 // S_bar = S^n, S^* or S^**
                     const MultiFab& S_stage_prim,
                     Vector<MultiFab>& S_data,                       // S_sum = most recent full solution
                     Vector<MultiFab>& S_scratch,                    // S_sum_old at most recent fast timestep for (rho theta)
                     const amrex::Geometry geom,
                     amrex::InterpFaceRegister* ifr,
                     const SolverChoice& solverChoice,
                     const MultiFab* z_t_pert,
                     std::unique_ptr<MultiFab>& z_phys_nd,
                     std::unique_ptr<MultiFab>& detJ_cc,
                     const MultiFab* r0, const MultiFab* pi0,
                     const amrex::Real dtau, const amrex::Real facinv)
{
    BL_PROFILE_REGION("erf_fast_rhs_T()");

    bool l_use_terrain  = solverChoice.use_terrain;
    bool l_move_terrain = (solverChoice.terrain_type == 1);

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

    MultiFab Delta_rho_u(    convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab Delta_rho_v(    convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab Delta_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, IntVect(1,1,0));
    MultiFab Delta_rho  (            ba                , dm, 1, 1);
    MultiFab Delta_rho_theta(        ba                , dm, 1, 1);

    MultiFab New_rho_u(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab New_rho_v(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab New_rho_w(convert(ba,IntVect(0,0,1)), dm, 1, 1);

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

    FArrayBox pifab;
    FArrayBox srcfab;

    FArrayBox coeff_A_fab;
    FArrayBox coeff_B_fab;
    FArrayBox coeff_C_fab;
    FArrayBox RHS_fab;
    FArrayBox soln_fab;
    FArrayBox gam_fab;

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

        const Array4<Real>& old_drho_u     = Delta_rho_u.array(mfi);
        const Array4<Real>& old_drho_v     = Delta_rho_v.array(mfi);
        const Array4<Real>& old_drho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& old_drho       = Delta_rho.array(mfi);
        const Array4<Real>& old_drho_theta = Delta_rho_theta.array(mfi);

        const Array4<const Real>& slow_rhs_cons  = S_slow_rhs[IntVar::cons].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_u = S_slow_rhs[IntVar::xmom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_v = S_slow_rhs[IntVar::ymom].const_array(mfi);
        const Array4<const Real>& slow_rhs_rho_w = S_slow_rhs[IntVar::zmom].const_array(mfi);

        const Array4<Real>& new_drho_u = New_rho_u.array(mfi);
        const Array4<Real>& new_drho_v = New_rho_v.array(mfi);
        const Array4<Real>& new_drho_w = New_rho_w.array(mfi);

        const Array4<Real>& cur_cons = S_data[IntVar::cons].array(mfi);
        const Array4<Real>& cur_xmom = S_data[IntVar::xmom].array(mfi);
        const Array4<Real>& cur_ymom = S_data[IntVar::ymom].array(mfi);
        const Array4<Real>& cur_zmom = S_data[IntVar::zmom].array(mfi);

        const Array4<Real>& scratch_rtheta = S_scratch[IntVar::cons].array(mfi);

        const Array4<const Real>& old_cons = S_old[IntVar::cons].const_array(mfi);
        const Array4<const Real>& old_xmom = S_old[IntVar::xmom].const_array(mfi);
        const Array4<const Real>& old_ymom = S_old[IntVar::ymom].const_array(mfi);
        const Array4<const Real>& old_zmom = S_old[IntVar::zmom].const_array(mfi);

        // These store the advection momenta which we will use to update the slow variables
        const Array4<Real>& avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        const Array4<Real>& avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        const Array4<Real>& avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real>& z_nd   = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ   = l_use_terrain ?   detJ_cc->const_array(mfi) : Array4<const Real>{};

        const Array4<const Real>& zp_t_arr = l_move_terrain ? z_t_pert->const_array(mfi) : Array4<const Real>{};

        const Array4<const Real>& r0_arr  = r0->const_array(mfi);
        const Array4<const Real>& pi0_arr = pi0->const_array(mfi);

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

        if (step == 0) {
            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cur_cons(i,j,k,Rho_comp)            = old_cons(i,j,k,Rho_comp);
                cur_cons(i,j,k,RhoTheta_comp)       = old_cons(i,j,k,RhoTheta_comp);
                scratch_rtheta(i,j,k,RhoTheta_comp) = old_cons(i,j,k,RhoTheta_comp);
            });

            amrex::ParallelFor(gtbx, gtby, gtbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                old_drho_u(i,j,k) = old_xmom(i,j,k) - stage_xmom(i,j,k);
            }, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                old_drho_v(i,j,k) = old_ymom(i,j,k) - stage_ymom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                old_drho_w(i,j,k) = old_zmom(i,j,k) - stage_zmom(i,j,k);
            });
        } else {
            amrex::ParallelFor(gtbx, gtby, gtbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                old_drho_u(i,j,k) = cur_xmom(i,j,k) - stage_xmom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                old_drho_v(i,j,k) = cur_ymom(i,j,k) - stage_ymom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                old_drho_w(i,j,k) = cur_zmom(i,j,k) - stage_zmom(i,j,k);
            });
        }

        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            old_drho(i,j,k)       = cur_cons(i,j,k,Rho_comp)      - stage_cons(i,j,k,Rho_comp);
            old_drho_theta(i,j,k) = cur_cons(i,j,k,RhoTheta_comp) - stage_cons(i,j,k,RhoTheta_comp);
            extrap_arr(i,j,k)     = old_drho_theta(i,j,k) + beta_d * (
              (cur_cons(i,j  ,k,RhoTheta_comp) - scratch_rtheta(i,j  ,k,RhoTheta_comp)));
        });

        Box tmpbox = bx;
        tmpbox.grow(Direction::x,1).grow(Direction::y,1);
        pifab.resize(tmpbox,1);
        auto const& pi_a = pifab.array();
        auto const& pi_ca = pifab.const_array();
        Elixir pieli = pifab.elixir();

        coeff_A_fab.resize(tbz,1);
        coeff_B_fab.resize(tbz,1);
        coeff_C_fab.resize(tbz,1);
        soln_fab.resize(tbz,1);
        gam_fab.resize(tbz,1);
        RHS_fab.resize(tbz,1);

        Elixir cAeli = coeff_A_fab.elixir();
        Elixir cBeli = coeff_B_fab.elixir();
        Elixir cCeli = coeff_C_fab.elixir();
        Elixir gCeli =     gam_fab.elixir();
        Elixir rCeli =     RHS_fab.elixir();
        Elixir sCeli =    soln_fab.elixir();

        auto const& coeffA_a  = coeff_A_fab.array();
        auto const& coeffB_a  = coeff_B_fab.array();
        auto const& coeffC_a  = coeff_C_fab.array();
        auto const& RHS_a     = RHS_fab.array();
        auto const& gam_a     = gam_fab.array();
        auto const& soln_a    = soln_fab.array();

        {
        BL_PROFILE("fast_rhs_eos_t");
            amrex::ParallelFor(tmpbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                pi_a(i,j,k,0) = getExnergivenRTh(stage_cons(i  ,j,k,RhoTheta_comp));
            });
        } // end profile

        if (l_use_terrain)
        {
        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        amrex::ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real pi_l = pi_ca(i-1,j,k,0);
                Real pi_r = pi_ca(i  ,j,k,0);
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
                Real fast_rhs_rho_u = -Gamma * R_d * pi_c * gpx;

                new_drho_u(i, j, k) = old_drho_u(i,j,k) + dtau * fast_rhs_rho_u
                                                        + dtau * slow_rhs_rho_u(i,j,k);

                if (k == domhi_z) new_drho_u(i,j,k+1) = new_drho_u(i,j,k);

                cur_xmom(i,j,k) = stage_xmom(i,j,k) + new_drho_u(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real pi_l = pi_ca(i,j-1,k,0);
                Real pi_r = pi_ca(i,j  ,k,0);
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
                Real fast_rhs_rho_v = -Gamma * R_d * pi_c * gpy;

                new_drho_v(i, j, k) = old_drho_v(i,j,k) + dtau * fast_rhs_rho_v
                                                        + dtau * slow_rhs_rho_v(i,j,k);
                if (k == domhi_z) new_drho_v(i,j,k+1) = new_drho_v(i,j,k);

                cur_ymom(i,j,k) = stage_ymom(i,j,k) + new_drho_v(i,j,k);
            });
    } else { // no terrain
            amrex::ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real pi_l = pi_ca(i-1,j,k,0);
                Real pi_r = pi_ca(i  ,j,k,0);
                Real pi_c =  0.5 * (pi_l + pi_r);
                Real drho_theta_hi = extrap_arr(i  ,j,k);
                Real drho_theta_lo = extrap_arr(i-1,j,k);

                Real gpx = (drho_theta_hi - drho_theta_lo)*dxi;

#ifdef ERF_USE_MOISTURE
                Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i-1,j,k,PrimQv_comp)
                                +prim(i,j,k,PrimQc_comp) + prim(i-1,j,k,PrimQc_comp) );
                gpx /= (1.0 + q);
#endif
                Real fast_rhs_rho_u = -Gamma * R_d * pi_c * gpx;

                new_drho_u(i, j, k) = old_drho_u(i,j,k) + dtau * fast_rhs_rho_u
                                                        + dtau * slow_rhs_rho_u(i,j,k);

                if (k == domhi_z) new_drho_u(i,j,k+1) = new_drho_u(i,j,k);

                cur_xmom(i,j,k) = stage_xmom(i,j,k) + new_drho_u(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real pi_l = pi_ca(i,j-1,k,0);
                Real pi_r = pi_ca(i,j  ,k,0);
                Real pi_c =  0.5 * (pi_l + pi_r);

                Real drho_theta_hi = extrap_arr(i,j,k);
                Real drho_theta_lo = extrap_arr(i,j-1,k);

                Real gpy = (drho_theta_hi - drho_theta_lo)*dyi;

#ifdef ERF_USE_MOISTURE
                Real q = 0.5 * ( prim(i,j,k,PrimQv_comp) + prim(i,j-1,k,PrimQv_comp)
                                +prim(i,j,k,PrimQc_comp) + prim(i,j-1,k,PrimQc_comp) );
                gpy /= (1.0 + q);
#endif
                Real fast_rhs_rho_v = -Gamma * R_d * pi_c * gpy;

                new_drho_v(i, j, k) = old_drho_v(i,j,k) + dtau * fast_rhs_rho_v
                                                        + dtau * slow_rhs_rho_v(i,j,k);
                if (k == domhi_z) new_drho_v(i,j,k+1) = new_drho_v(i,j,k);

                cur_ymom(i,j,k) = stage_ymom(i,j,k) + new_drho_v(i,j,k);
            });
    } // end of no terrain

        // *********************************************************************
        // *********************************************************************
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
            Real rhobar_lo, rhobar_hi, pibar_lo, pibar_hi;
            rhobar_lo =  r0_arr(i,j,k-1);
            rhobar_hi =  r0_arr(i,j,k  );
             pibar_lo = pi0_arr(i,j,k-1);
             pibar_hi = pi0_arr(i,j,k  );

            Real pi_lo = pi_ca(i,j,k-1,0);
            Real pi_hi = pi_ca(i,j,k  ,0);
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
                         (  c_v * pibar_hi * stage_cons(i,j,k,RhoTheta_comp) );

            Real coeff_Q = Gamma * R_d * pi_c * dzi / h_zeta_on_kface
                         + halfg * R_d * rhobar_lo * pi_lo  /
                         ( c_v  * pibar_lo * stage_cons(i,j,k-1,RhoTheta_comp) );

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
            coeffA_a(i,j,k) = D * ( halfg - coeff_Q * theta_t_lo );
            coeffC_a(i,j,k) = D * (-halfg + coeff_P * theta_t_hi );
            coeffB_a(i,j,k) = 1.0 + D * (coeff_Q - coeff_P) * theta_t_mid;

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
            RHS_a(i,j,k) = old_drho_w(i,j,k) + dtau * (slow_rhs_rho_w(i,j,k) + R_tmp);

            // Terrain solves for Omega
            if (l_use_terrain)
                RHS_a(i,j,k) += OmegaFromW(i,j,k,0.,new_drho_u,new_drho_v,z_nd,dxInv);
      });
      } // end profile

      int klen = bx.bigEnd(2) - bx.smallEnd(2) + 1;

      amrex::Box b2d = tbz; // Copy constructor
      b2d.setRange(2,0);

      {
      BL_PROFILE("fast_rhs_b2d_loop_t");
#ifdef AMREX_USE_GPU
      ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
      {
          // w_0 = 0
          coeffA_a(i,j,0) =  0.0;
          coeffB_a(i,j,0) =  1.0;
          coeffC_a(i,j,0) =  0.0;
          RHS_a   (i,j,0) =  0.0;

          // Moving terrain
          if (l_use_terrain && l_move_terrain) {
              Real time_mt  = time - 0.5*dtau;
              Real omega_bc = dhdt(i,j,dx,time_mt,dtau) - stage_zmom(i,j,0) / stage_cons(i,j,0); // w'' = w(t) - w^{RK stage}
              RHS_a(i,j,0) = cur_cons(i,j,0) * omega_bc;
          }

          // w_khi = 0
          coeffA_a(i,j,klen) =  0.0;
          coeffB_a(i,j,klen) =  1.0;
          coeffC_a(i,j,klen) =  0.0;
          RHS_a(i,j,klen)     =  0.0;

          // w = 0 at k = 0
          Real bet = coeffB_a(i,j,0);
          soln_a(i,j,0) = RHS_a(i,j,0) / bet;

          for (int k = 1; k <= klen; k++) {
              gam_a(i,j,k) = coeffC_a(i,j,k-1) / bet;
              bet = coeffB_a(i,j,k) - coeffA_a(i,j,k)*gam_a(i,j,k);
              soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) / bet;
          }
          for (int k = klen-1; k >= 0; k--) {
              soln_a(i,j,k) -= gam_a(i,j,k+1)*soln_a(i,j,k+1);
          }
      });
#else
        auto const lo = amrex::lbound(bx);
        auto const hi = amrex::ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                coeffA_a(i,j,0) =  0.0;
                coeffB_a(i,j,0) =  1.0;
                coeffC_a(i,j,0) =  0.0;
                RHS_a   (i,j,0) =  0.0;
                soln_a(i,j,0) = RHS_a(i,j,0) / coeffB_a(i,j,0);
            }
        }
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                coeffA_a(i,j,hi.z+1) =  0.0;
                coeffB_a(i,j,hi.z+1) =  1.0;
                coeffC_a(i,j,hi.z+1) =  0.0;
                RHS_a   (i,j,hi.z+1) =  0.0;
            }
        }
        for (int k = lo.z+1; k <= hi.z+1; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    gam_a(i,j,k) = coeffC_a(i,j,k-1) / coeffB_a(i,j,k-1);
                    Real bet = coeffB_a(i,j,k) - coeffA_a(i,j,k)*gam_a(i,j,k);
                    coeffB_a(i,j,k) = bet;
                    soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) / bet;
                }
            }
        }
        for (int k = hi.z; k >= lo.z; --k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    soln_a(i,j,k) -= gam_a(i,j,k+1)*soln_a(i,j,k+1);
                }
            }
        }
#endif
      } // end profile

      {
      BL_PROFILE("fast_rhs_new_drhow_t");
      ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
            if (l_use_terrain) {
                new_drho_w(i,j,k) = WFromOmega(i,j,k,soln_a(i,j,k),new_drho_u,new_drho_v,z_nd,dxInv);
            } else {
                new_drho_w(i,j,k) = soln_a(i,j,k);
            }
            cur_zmom(i,j,k) = stage_zmom(i,j,k) + new_drho_w(i,j,k);

            // Sum implicit and explicit W for AdvSrc
            new_drho_w(i,j,k) *= beta_2;
            new_drho_w(i,j,k) += beta_1 * old_drho_w(i,j,k);
      });
      } // end profile

      // **************************************************************************
      // Define updates in the RHS of rho and (rho theta)
      // **************************************************************************

      const int l_spatial_order = 2;

      srcfab.resize(bx,2);
      auto const& src_a = srcfab.array();
      Elixir src_eli = srcfab.elixir();

      AdvectionSrcForRhoAndTheta(bx, valid_bx, src_a,
                                 new_drho_u, new_drho_v, new_drho_w,      // these are being used to build the fluxes
                                 facinv, avg_xmom, avg_ymom, avg_zmom,    // these are being defined from the rho fluxes
                                 zp_t_arr,  // z_t
                                 prim, z_nd, detJ, dxInv, l_spatial_order, l_use_terrain);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          cur_cons(i,j,k,0) += dtau * (slow_rhs_cons(i,j,k,0) + src_a(i,j,k,0));
          cur_cons(i,j,k,1) += dtau * (slow_rhs_cons(i,j,k,1) + src_a(i,j,k,1));
      });
    } // mfi
    }
}
