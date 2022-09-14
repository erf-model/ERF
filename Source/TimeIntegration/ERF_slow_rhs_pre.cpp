
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <ERF_Constants.H>
//#include <ABLMost.H>
#include <Advection.H>
#include <Diffusion.H>
#include <TimeIntegration.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>

using namespace amrex;

void erf_slow_rhs_pre (int level,
                   Vector<MultiFab>& S_rhs,
                   Vector<MultiFab>& S_data,
                   const MultiFab& S_prim,
                         Vector<MultiFab>& S_scratch,
                   const MultiFab& xvel,
                   const MultiFab& yvel,
                   const MultiFab& zvel,
                   std::unique_ptr<MultiFab>& z_t_mf,
                         MultiFab& Omega,
                   const MultiFab& source,
                   const MultiFab& eddyDiffs,
                   std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
                   const amrex::Geometry geom,
                         amrex::InterpFaceRegister* ifr,
                   const SolverChoice& solverChoice,
                   std::unique_ptr<ABLMost>& most,
                   const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                   std::unique_ptr<MultiFab>& z_phys_nd, std::unique_ptr<MultiFab>& dJ,
                   const MultiFab* r0, const MultiFab* p0,
                   const amrex::Real* dptr_rayleigh_tau, const amrex::Real* dptr_rayleigh_ubar,
                   const amrex::Real* dptr_rayleigh_vbar, const amrex::Real* dptr_rayleigh_thetabar,
                   bool ingested_bcs)
{
    BL_PROFILE_REGION("erf_slow_rhs_pre()");

    amrex::Real theta_mean;
    if (most) theta_mean = most->theta_mean;

    int start_comp = 0;
    int   num_comp = 2;

    const int l_spatial_order = solverChoice.spatial_order;
    const int l_use_terrain   = solverChoice.use_terrain;

    const amrex::BCRec* bc_ptr = domain_bcs_type_d.data();

    const Box& domain = geom.Domain();
    const int domhi_z = domain.bigEnd()[2];

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *************************************************************************
    // Combine external forcing terms
    // *************************************************************************
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const GpuArray<Real, AMREX_SPACEDIM> ext_forcing = {
       -solverChoice.abl_pressure_grad[0] + solverChoice.abl_geo_forcing[0],
       -solverChoice.abl_pressure_grad[1] + solverChoice.abl_geo_forcing[1],
       -solverChoice.abl_pressure_grad[2] + solverChoice.abl_geo_forcing[2]};

    const iMultiFab *mlo_mf_x, *mhi_mf_x;
    const iMultiFab *mlo_mf_y, *mhi_mf_y;

    if (level > 0)
    {
        mlo_mf_x = &(ifr->mask(Orientation(0,Orientation::low)));
        mhi_mf_x = &(ifr->mask(Orientation(0,Orientation::high)));
        mlo_mf_y = &(ifr->mask(Orientation(1,Orientation::low)));
        mhi_mf_y = &(ifr->mask(Orientation(1,Orientation::high)));
    }

    MultiFab pprime(S_data[IntVar::cons].boxArray(), S_data[IntVar::cons].DistributionMap(), 1, 1);
    MultiFab expr(S_data[IntVar::cons].boxArray(), S_data[IntVar::cons].DistributionMap(), 1, IntVect(1,1,0));

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        bool left_edge_dirichlet = false;
        bool rght_edge_dirichlet = false;
        bool  bot_edge_dirichlet = false;
        bool  top_edge_dirichlet = false;

        const Box& valid_bx = mfi.validbox();

        if (level == 0 && ingested_bcs) {
            left_edge_dirichlet = (bx.smallEnd(0) == domain.smallEnd(0));
            rght_edge_dirichlet = (bx.bigEnd(1)   == domain.bigEnd(0));
            bot_edge_dirichlet  = (bx.smallEnd(0) == domain.smallEnd(1));
            top_edge_dirichlet  = (bx.bigEnd(1)   == domain.bigEnd(1));
        } else if (level > 0) {
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
            left_edge_dirichlet = mlo_x(vlo_x  ,vlo_y  ,vlo_z);
            rght_edge_dirichlet = mhi_x(vhi_x+1,vhi_y  ,vhi_z);
            bot_edge_dirichlet  = mlo_y(vlo_x  ,vlo_y  ,vlo_z);
            top_edge_dirichlet  = mhi_y(vhi_x  ,vhi_y+1,vhi_z);
        } // level > 0

        if (left_edge_dirichlet) tbx.growLo(0,-1);
        if (rght_edge_dirichlet) tbx.growHi(0,-1);
        if ( bot_edge_dirichlet) tby.growLo(1,-1);
        if ( top_edge_dirichlet) tby.growHi(1,-1);

        // We don't compute a source term for z-momentum on the bottom or top boundary
        tbz.growLo(2,-1);
        tbz.growHi(2,-1);

        const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real> & cell_rhs   = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & source_fab = source.const_array(mfi);

        // We must initialize these to zero each RK step
        S_scratch[IntVar::xmom][mfi].template setVal<RunOn::Device>(0.);
        S_scratch[IntVar::ymom][mfi].template setVal<RunOn::Device>(0.);
        S_scratch[IntVar::zmom][mfi].template setVal<RunOn::Device>(0.);

        Array4<Real> avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<const Real>& rho_u = S_data[IntVar::xmom].array(mfi);
        const Array4<const Real>& rho_v = S_data[IntVar::ymom].array(mfi);
        const Array4<const Real>& rho_w = S_data[IntVar::zmom].array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        Array4<const Real> z_t;
        if (z_t_mf)
            z_t = z_t_mf->array(mfi);
        else
            z_t = Array4<const Real>{};

        const Array4<Real>& rho_u_rhs = S_rhs[IntVar::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVar::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVar::zmom].array(mfi);

        // These are temporaries we use to add to the S_rhs for the fluxes
        const Array4<Real>& diffflux_x = diffflux[0].array(mfi);
        const Array4<Real>& diffflux_y = diffflux[1].array(mfi);
        const Array4<Real>& diffflux_z = diffflux[2].array(mfi);

        const Array4<const Real>& K_turb = eddyDiffs.const_array(mfi);

        const Array4<const Real>& z_nd = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ = l_use_terrain ?        dJ->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& r0_arr = r0->const_array(mfi);
        const Array4<const Real>& p0_arr = p0->const_array(mfi);

        const Box& gbx = mfi.growntilebox(1);
        const Array4<Real> & pp_arr  = pprime.array(mfi);
        {
        BL_PROFILE("slow_rhs_pre_pprime");
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            //if (cell_data(i,j,k,RhoTheta_comp) < 0.) printf("BAD THETA AT %d %d %d %e %e \n",
            //    i,j,k,cell_data(i,j,k,RhoTheta_comp),cell_data(i,j,k+1,RhoTheta_comp));
            AMREX_ASSERT(cell_data(i,j,k,RhoTheta_comp) > 0.);
            pp_arr(i,j,k) = getPprimegivenRTh(cell_data(i,j,k,RhoTheta_comp),p0_arr(i,j,k));
        });
        } // end profile

        const Array4<Real> & er_arr = expr.array(mfi);

        {
        BL_PROFILE("slow_rhs_making_er");
        if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
             (solverChoice.les_type        !=       LESType::None) ||
             (solverChoice.pbl_type        !=       PBLType::None) )
        {
            const Box& gbx2 = mfi.growntilebox(IntVect(1,1,0));

            if (l_use_terrain) {
                // First create Omega using velocity (not momentum)
                Box gbxo = mfi.nodaltilebox(2);gbxo.grow(IntVect(1,1,0));
                amrex::ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = (k == 0) ? 0. : OmegaFromW(i,j,k,w(i,j,k),u,v,z_nd,dxInv);
                });

                amrex::ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

                    Real met_u_h_zeta_hi = Compute_h_zeta_AtIface(i+1, j  , k, dxInv, z_nd);
                    Real met_u_h_zeta_lo = Compute_h_zeta_AtIface(i  , j  , k, dxInv, z_nd);

                    Real met_v_h_zeta_hi = Compute_h_zeta_AtJface(i  , j+1, k, dxInv, z_nd);
                    Real met_v_h_zeta_lo = Compute_h_zeta_AtJface(i  , j  , k, dxInv, z_nd);

                    Real Omega_hi = omega_arr(i,j,k+1);
                    Real Omega_lo = omega_arr(i,j,k  );

                    Real expansionRate = (u(i+1,j  ,k)*met_u_h_zeta_hi - u(i,j,k)*met_u_h_zeta_lo)*dxInv[0] +
                                         (v(i  ,j+1,k)*met_v_h_zeta_hi - v(i,j,k)*met_v_h_zeta_lo)*dxInv[1] +
                                         (Omega_hi - Omega_lo)*dxInv[2];

                    er_arr(i,j,k) = expansionRate / detJ(i,j,k);
                });

            } else {
                amrex::ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    er_arr(i,j,k) = (u(i+1, j  , k  ) - u(i, j, k))*dxInv[0] +
                                    (v(i  , j+1, k  ) - v(i, j, k))*dxInv[1] +
                                    (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
                });
            }
        }
        } // end profile

        {
        BL_PROFILE("slow_rhs_making_omega");
            Box gbxo = mfi.nodaltilebox(2);gbxo.grow(IntVect(1,1,0));
            // Now create Omega with momentum (not velocity) with z_t subtracted if moving terrain
            if (l_use_terrain) {
                if (z_t) {
                    amrex::ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        Real rho_at_face;
                        if (w(i,j,k) != 0.)
                            rho_at_face = rho_w(i,j,k) / w(i,j,k);
                        else
                            rho_at_face = 0.;
                        omega_arr(i,j,k) = (k == 0) ? 0. : OmegaFromW(i,j,k,rho_w(i,j,k),rho_u,rho_v,z_nd,dxInv) -
                            rho_at_face * z_t(i,j,k);
                    });
                } else {
                    amrex::ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = (k == 0) ? 0. : OmegaFromW(i,j,k,rho_w(i,j,k),rho_u,rho_v,z_nd,dxInv);
                    });
                }
            } else {
                amrex::ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = rho_w(i,j,k);
                });
            }
        } // end profile

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************

        Real fac = 1.0;
        AdvectionSrcForRhoAndTheta(bx, valid_bx, cell_rhs,       // these are being used to build the fluxes
                                   rho_u, rho_v, omega_arr, fac,
                                   avg_xmom, avg_ymom, avg_zmom, // these are being defined from the rho fluxes
                                   cell_prim, z_nd, detJ,
                                   dxInv, l_spatial_order, l_use_terrain);

        // NOTE: No diffusion for continuity, so n starts at 1.
        //       KE calls moved inside DiffSrcForState.
        int n_start = amrex::max(start_comp,RhoTheta_comp);
        int n_end   = start_comp + num_comp - 1;
        DiffusionSrcForState(bx, domain, n_start, n_end, u, v, w,
                             cell_data, cell_prim, source_fab, cell_rhs,
                             diffflux_x, diffflux_y, diffflux_z,
                             dxInv, K_turb, solverChoice, theta_mean, grav_gpu, bc_ptr);

        // Add Rayleigh damping
        if (solverChoice.use_rayleigh_damping) {
            int n  = RhoTheta_comp;
            int nr = Rho_comp;
            int np = PrimTheta_comp;
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real theta = cell_prim(i,j,k,np);
                cell_rhs(i, j, k, n) -= dptr_rayleigh_tau[k] * (theta - dptr_rayleigh_thetabar[k]) * cell_data(i,j,k,nr);
            });
        }

        AdvectionSrcForMom(tbx, tby, tbz,
                           rho_u_rhs, rho_v_rhs, rho_w_rhs, u, v, w,
                           rho_u    , rho_v    , omega_arr,
                           z_nd, detJ, dxInv, l_spatial_order, l_use_terrain, domhi_z);

        if (l_use_terrain) {
            DiffusionSrcForMom_T(tbx, tby, tbz, domain, rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                 u, v, w, K_turb, cell_data, er_arr,
                                 solverChoice, bc_ptr, z_nd, detJ, dxInv);
        } else {
            DiffusionSrcForMom_N(tbx, tby, tbz, domain,  rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                 u, v, w, K_turb, cell_data, er_arr,
                                 solverChoice, bc_ptr, dxInv);
        }

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        {
        BL_PROFILE("slow_rhs_pre_xmom");
        amrex::ParallelFor(tbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // x-momentum equation
            // Add pressure gradient
            amrex::Real gpx;
            if (l_use_terrain) {

                Real met_h_xi   = Compute_h_xi_AtIface  (i, j, k, dxInv, z_nd);
                Real met_h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);

                Real gp_xi = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));
                Real gp_zeta_on_iface;
                if(k==0) {
                    gp_zeta_on_iface = 0.5 * dxInv[2] * (
                                                        pp_arr(i-1,j,k+1) + pp_arr(i,j,k+1)
                                                      - pp_arr(i-1,j,k  ) - pp_arr(i,j,k  ) );
                } else if (k==domhi_z) {
                    gp_zeta_on_iface = 0.5 * dxInv[2] * (
                                                        pp_arr(i-1,j,k  ) + pp_arr(i,j,k  )
                                                      - pp_arr(i-1,j,k-1) - pp_arr(i,j,k-1) );
                } else {
                    gp_zeta_on_iface = 0.25 * dxInv[2] * (
                                                         pp_arr(i-1,j,k+1) + pp_arr(i,j,k+1)
                                                       - pp_arr(i-1,j,k-1) - pp_arr(i,j,k-1) );
                }
                gpx = gp_xi - (met_h_xi/ met_h_zeta) * gp_zeta_on_iface;
            } else {
                gpx = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));
            }

#ifdef ERF_USE_MOISUTRE
            Real q = 0.5 * ( cell_prim(i,j,k,PrimQv_comp) + cell_prim(i-1,j,k,PrimQv_comp)
                            +cell_prim(i,j,k,PrimQc_comp) + cell_prim(i-1,j,k,PrimQc_comp) );
            rho_u_rhs(i, j, k) -= gpx / (1.0 + q);
#else
            rho_u_rhs(i, j, k) -= gpx;
#endif

            // Add external drivers
            rho_u_rhs(i, j, k) += ext_forcing[0];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                rho_u_rhs(i, j, k) += solverChoice.coriolis_factor *
                        (rho_v_loc * solverChoice.sinphi - rho_w_loc * solverChoice.cosphi);
            }

            // Add Rayleigh damping
            if (solverChoice.use_rayleigh_damping)
            {
                Real uu = rho_u(i,j,k) / cell_data(i,j,k,Rho_comp);
                rho_u_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (uu - dptr_rayleigh_ubar[k]) * cell_data(i,j,k,Rho_comp);
            }
        });
        } // end profile
        {
        BL_PROFILE("slow_rhs_pre_ymom");
        amrex::ParallelFor(tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // y-momentum equation
                // Add pressure gradient
                amrex::Real gpy;
                if (l_use_terrain) {

                    Real met_h_eta  = Compute_h_eta_AtJface (i, j, k, dxInv, z_nd);
                    Real met_h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);

                    Real gp_eta = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));
                    Real gp_zeta_on_jface;
                    if(k==0) {
                        gp_zeta_on_jface = 0.5 * dxInv[2] * (
                                                            pp_arr(i,j,k+1) + pp_arr(i,j-1,k+1)
                                                          - pp_arr(i,j,k  ) - pp_arr(i,j-1,k  ) );
                    } else if (k==domhi_z) {
                        gp_zeta_on_jface = 0.5 * dxInv[2] * (
                                                            pp_arr(i,j,k  ) + pp_arr(i,j-1,k  )
                                                          - pp_arr(i,j,k-1) - pp_arr(i,j-1,k-1) );
                    } else {
                        gp_zeta_on_jface = 0.25 * dxInv[2] * (
                                                             pp_arr(i,j,k+1) + pp_arr(i,j-1,k+1)
                                                           - pp_arr(i,j,k-1) - pp_arr(i,j-1,k-1) );
                    }
                    gpy = gp_eta - (met_h_eta / met_h_zeta) * gp_zeta_on_jface;
                } else {
                    gpy = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));
                }

#ifdef ERF_USE_MOISUTRE
                Real q = 0.5 * ( cell_prim(i,j,k,PrimQv_comp) + cell_prim(i,j-1,k,PrimQv_comp)
                                +cell_prim(i,j,k,PrimQc_comp) + cell_prim(i,j-1,k,PrimQc_comp) );
                rho_v_rhs(i, j, k) -= gpy / (1.0_rt + q);
#else
                rho_v_rhs(i, j, k) -= gpy;
#endif

                // Add external drivers
                rho_v_rhs(i, j, k) += ext_forcing[1];

                // Add Coriolis forcing (that assumes east is +x, north is +y)
                if (solverChoice.use_coriolis)
                {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                    rho_v_rhs(i, j, k) += -solverChoice.coriolis_factor * rho_u_loc * solverChoice.sinphi;
                }

                // Add Rayleigh damping
                if (solverChoice.use_rayleigh_damping)
                {
                    Real vv = rho_v(i,j,k) / cell_data(i,j,k,Rho_comp);
                    rho_v_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (vv - dptr_rayleigh_vbar[k]) * cell_data(i,j,k,Rho_comp);
                }
        });
        } // end profile
        {
        BL_PROFILE("slow_rhs_pre_zmom_2d");
        amrex::Box b2d = tbz;
        b2d.setSmall(2,0);
        b2d.setBig(2,0);
        // Enforce no forcing term at top and bottom boundaries
        amrex::ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rho_w_rhs(i,j,        0) = 0.;
            rho_w_rhs(i,j,domhi_z+1) = 0.;
        });
        } // end profile

        {
        BL_PROFILE("slow_rhs_pre_zmom");
        amrex::ParallelFor(tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation
                // Add pressure gradient
                amrex::Real gpz;
                if (l_use_terrain) {
                    Real met_h_zeta = Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd);
                    gpz = dxInv[2] * (pp_arr(i,j,k) - pp_arr(i,j,k-1)) / met_h_zeta;
                } else {
                    gpz = dxInv[2] * (pp_arr(i,j,k) - pp_arr(i,j,k-1));
                }

#ifdef ERF_USE_MOISUTRE
                Real q = 0.5 * ( cell_prim(i,j,k,PrimQv_comp) + cell_prim(i,j,k-1,PrimQv_comp)
                                +cell_prim(i,j,k,PrimQc_comp) + cell_prim(i,j,k-1,PrimQc_comp) );
                rho_w_rhs(i, j, k) -= gpz / (1.0_rt + q);
#else
                rho_w_rhs(i, j, k) -= gpz;
#endif

                // Add gravity term
                if (solverChoice.use_gravity)
                {
                    Real rho_prime = 0.5 * (cell_data(i,j,k) + cell_data(i,j,k-1) - r0_arr(i,j,k) - r0_arr(i,j,k-1));
                    rho_w_rhs(i, j, k) += grav_gpu[2] * rho_prime;
                }

                // Add external drivers
                rho_w_rhs(i, j, k) += ext_forcing[2];

                // Add Coriolis forcing (that assumes east is +x, north is +y)
                if (solverChoice.use_coriolis)
                {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                    rho_w_rhs(i, j, k) += solverChoice.coriolis_factor * rho_u_loc * solverChoice.cosphi;
                }

                // Add Rayleigh damping
                if (solverChoice.use_rayleigh_damping)
                {
                    rho_w_rhs(i, j, k) -= dptr_rayleigh_tau[k] * rho_w(i,j,k);
                }
        });
        } // end profile
    } // mfi
}
