#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <ERF_Constants.H>
#include <EddyViscosity.H>
#include <ABLMost.H>
#include <PBLModels.H>
#include <SpatialStencils.H>
#include <TimeIntegration.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>
#include <Interpolation.H>

using namespace amrex;

void erf_slow_rhs (int level,
                   Vector<MultiFab>& S_rhs,
                   const Vector<MultiFab>& S_data,
                   const MultiFab& S_prim,
                         Vector<MultiFab>& S_scratch,
                   const MultiFab& xvel,
                   const MultiFab& yvel,
                   const MultiFab& zvel,
                   MultiFab& source,
                   std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
                   const amrex::Geometry geom,
                         amrex::InterpFaceRegister* ifr,
                   const SolverChoice& solverChoice,
                   std::unique_ptr<ABLMost>& most,
                   const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                   std::unique_ptr<MultiFab>& z0, std::unique_ptr<MultiFab>& dJ,
                   const MultiFab* r0, const MultiFab* p0,
                   const amrex::Real* dptr_rayleigh_tau, const amrex::Real* dptr_rayleigh_ubar,
                   const amrex::Real* dptr_rayleigh_vbar, const amrex::Real* dptr_rayleigh_thetabar,
                   const int rhs_vars)
{
    BL_PROFILE_VAR("erf_slow_rhs()",erf_slow_rhs);

    amrex::Real theta_mean;
    if (most) theta_mean = most->theta_mean;
      
    int start_comp;
    int   num_comp;
    if (rhs_vars == RHSVar::fast) {
        start_comp = 0;
          num_comp = 2; // Rho_comp and RhoTheta_comp
    } else if (rhs_vars == RHSVar::slow) {
        start_comp = 2;
          num_comp = S_data[IntVar::cons].nComp() - 2;
    } else if (rhs_vars == RHSVar::all) {
        start_comp = 0;
          num_comp = S_data[IntVar::cons].nComp();
    }

    const int l_spatial_order = solverChoice.spatial_order;
    const int l_use_terrain   = solverChoice.use_terrain;

    const amrex::BCRec* bc_ptr = domain_bcs_type_d.data();

    const Box& domain = geom.Domain();
    const int domhi_z = domain.bigEnd()[2];

    const GpuArray<Real, AMREX_SPACEDIM> dx    = geom.CellSizeArray();
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *************************************************************************
    // Calculate cell-centered eddy viscosity & diffusivities
    //
    // Notes -- we fill all the data in ghost cells before calling this so
    //    that we can fill the eddy viscosity in the ghost regions and
    //    not have to call a boundary filler on this data itself
    //
    // LES - updates both horizontal and vertical eddy viscosity components
    // PBL - only updates vertical eddy viscosity components so horizontal
    //       components come from the LES model or are left as zero.
    // *************************************************************************
    MultiFab eddyDiffs(S_data[IntVar::cons].boxArray(),S_data[IntVar::cons].DistributionMap(),EddyDiff::NumDiffs,1);
    ComputeTurbulentViscosity(xvel, yvel, zvel, S_data[IntVar::cons],
                              eddyDiffs, geom, solverChoice, most, domain_bcs_type_d);

    const iMultiFab *mlo_mf_x, *mhi_mf_x;
    const iMultiFab *mlo_mf_y, *mhi_mf_y;
    const iMultiFab *mlo_mf_z, *mhi_mf_z;

    bool l_use_QKE       = solverChoice.use_QKE && solverChoice.advect_QKE;
    bool l_use_deardorff = (solverChoice.les_type == LESType::Deardorff);

    if (level > 0)
    {
        mlo_mf_x = &(ifr->mask(Orientation(0,Orientation::low)));
        mhi_mf_x = &(ifr->mask(Orientation(0,Orientation::high)));
        mlo_mf_y = &(ifr->mask(Orientation(1,Orientation::low)));
        mhi_mf_y = &(ifr->mask(Orientation(1,Orientation::high)));
        mlo_mf_z = &(ifr->mask(Orientation(2,Orientation::low)));
        mhi_mf_z = &(ifr->mask(Orientation(2,Orientation::high)));
    }

    MultiFab pprime(S_data[IntVar::cons].boxArray(), S_data[IntVar::cons].DistributionMap(), 1, 1);
    MultiFab expr(S_data[IntVar::cons].boxArray(), S_data[IntVar::cons].DistributionMap(), 1, IntVect(1,1,0));

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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

        const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real> & cell_rhs   = S_rhs[IntVar::cons].array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        Array4<Real> avg_xmom;
        Array4<Real> avg_ymom;
        Array4<Real> avg_zmom;
        if (S_scratch.size() > 0) {
            avg_xmom = S_scratch[IntVar::xmom].array(mfi);
            avg_ymom = S_scratch[IntVar::ymom].array(mfi);
            avg_zmom = S_scratch[IntVar::zmom].array(mfi);
        }

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<const Real>& rho_u = S_data[IntVar::xmom].array(mfi);
        const Array4<const Real>& rho_v = S_data[IntVar::ymom].array(mfi);
        const Array4<const Real>& rho_w = S_data[IntVar::zmom].array(mfi);

        const Array4<Real>& rho_u_rhs = S_rhs[IntVar::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVar::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVar::zmom].array(mfi);

        const Array4<Real>& xflux_rhs = S_rhs[IntVar::xflux].array(mfi);
        const Array4<Real>& yflux_rhs = S_rhs[IntVar::yflux].array(mfi);
        const Array4<Real>& zflux_rhs = S_rhs[IntVar::zflux].array(mfi);

        // These are temporaries we use to add to the S_rhs for the fluxes
        const Array4<Real>& diffflux_x = diffflux[0].array(mfi);
        const Array4<Real>& diffflux_y = diffflux[1].array(mfi);
        const Array4<Real>& diffflux_z = diffflux[2].array(mfi);

        const Array4<Real>& K_turb = eddyDiffs.array(mfi);

        const Array4<const Real>& z_nd = l_use_terrain ? z0->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ = l_use_terrain ? dJ->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& r0_arr = r0->const_array(mfi);
        const Array4<const Real>& p0_arr = p0->const_array(mfi);

        const Box& gbx = mfi.growntilebox(1);
        const Array4<Real> & pp_arr  = pprime.array(mfi);
        if (rhs_vars != RHSVar::slow) {
            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (cell_data(i,j,k,RhoTheta_comp) < 0.) printf("BAD THETA AT %d %d %d %e %e \n",
                    i,j,k,cell_data(i,j,k,RhoTheta_comp),cell_data(i,j,k+1,RhoTheta_comp));
                AMREX_ALWAYS_ASSERT(cell_data(i,j,k,RhoTheta_comp) > 0.);
                pp_arr(i,j,k) = getPprimegivenRTh(cell_data(i,j,k,RhoTheta_comp),p0_arr(i,j,k));
            });
        }

        const Box& gbx2 = mfi.growntilebox(IntVect(1,1,0));
        const Array4<Real> & er_arr  = expr.array(mfi);
        amrex::ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                er_arr(i,j,k) = (u(i+1, j  , k  ) - u(i, j, k))*dxInv[0] +
                                (v(i  , j+1, k  ) - v(i, j, k))*dxInv[1] +
                                (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
        });

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************

        // NOTE: given how we fill the fluxes, we must call AdvectionSrcForState before
        //       we call DiffusionSrcForState
        if (rhs_vars == RHSVar::fast || rhs_vars == RHSVar::all) {
            AdvectionSrcForState(bx, start_comp, num_comp, rho_u, rho_v, rho_w, cell_prim,
                                 cell_rhs, xflux_rhs, yflux_rhs, zflux_rhs, z_nd, detJ,
                                 dxInv, l_spatial_order, l_use_terrain, l_use_deardorff, l_use_QKE);
        } else if (rhs_vars == RHSVar::slow) {
            AdvectionSrcForState(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom, cell_prim,
                                 cell_rhs, xflux_rhs, yflux_rhs, zflux_rhs, z_nd, detJ,
                                 dxInv, l_spatial_order, l_use_terrain, l_use_deardorff, l_use_QKE);
        }

        // NOTE: No diffusion for continuity, so n starts at 1.
        //       KE calls moved inside DiffSrcForState.
        int n_start = amrex::max(start_comp,RhoTheta_comp);
        int n_end   = start_comp + num_comp;
        for (int n(n_start); n<n_end; ++n) {
            DiffusionSrcForState(bx, domain, n, u, v, w,
                                 cell_data, cell_prim, source_fab, cell_rhs,
                                 diffflux_x, diffflux_y, diffflux_z,
                                 dxInv, K_turb, solverChoice, theta_mean, grav_gpu,
                                 bc_ptr, dptr_rayleigh_tau, dptr_rayleigh_thetabar);
        }

        // Compute the RHS for the flux terms from this stage -- we do it this way so we don't double count
        //         fluxes at fine-fine interfaces
        int rc = Rho_comp;
        int update_mom = (S_scratch.size() > 0);
        amrex::ParallelFor(tbx, num_comp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n_in) noexcept
        {
             int n = start_comp + n_in;
             if (update_mom && n == rc) // Only update if we are computing source for density
                 avg_xmom(i,j,k) = xflux_rhs(i,j,k,0);

             xflux_rhs(i,j,k,n) += diffflux_x(i,j,k,n);
        });
        amrex::ParallelFor(tby, num_comp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n_in) noexcept
        {
             int n = start_comp + n_in;
             if (update_mom && n == rc) // Only update if we are computing source for density
                 avg_ymom(i,j,k) = yflux_rhs(i,j,k,0);

             yflux_rhs(i,j,k,n) += diffflux_y(i,j,k,n);
        });
        amrex::ParallelFor(tbz, num_comp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n_in) noexcept
        {
             int n = start_comp + n_in;
             if (update_mom && n == rc) // Only update if we are computing source for density
                 avg_zmom(i,j,k) = zflux_rhs(i,j,k,0);

             zflux_rhs(i,j,k,n) += diffflux_z(i,j,k,n);
        });

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        if (rhs_vars != RHSVar::slow) {
        amrex::ParallelFor(tbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // x-momentum equation

            rho_u_rhs(i, j, k) = 0.0; // Initialize the updated x-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

            // Add advective terms
            rho_u_rhs(i, j, k) += -AdvectionSrcForXMom(i, j, k, rho_u, rho_v, rho_w, u, z_nd, detJ,
                                                       dxInv, l_spatial_order, l_use_terrain);

            // Add diffusive terms
            amrex::Real diff_update;
            if (k < domhi_z && l_use_terrain) {
                diff_update = DiffusionSrcForMomWithTerrain(i, j, k, u, v, w, cell_data,
                                                            MomentumEqn::x, dxInv, K_turb, solverChoice,
                                                            z_nd, detJ, domain, bc_ptr);
            } else {
                diff_update = DiffusionSrcForMom(i, j, k, u, v, w, cell_data,
                                                 MomentumEqn::x, dxInv, K_turb, solverChoice,
                                                 domain, bc_ptr, er_arr);
            }
            rho_u_rhs(i, j, k) += diff_update;

            // Add pressure gradient
            amrex::Real gpx;
            if (l_use_terrain) {
                Real met_h_xi, met_h_eta, met_h_zeta;

                ComputeMetricAtIface(i,j,k,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::h_xi_zeta);

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

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_u_rhs(i, j, k) += -solverChoice.abl_pressure_grad[0];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                rho_u_rhs(i, j, k) += solverChoice.coriolis_factor *
                        (rho_v_loc * solverChoice.sinphi - rho_w_loc * solverChoice.cosphi);
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_u_rhs(i, j, k) += solverChoice.abl_geo_forcing[0];

            // Add Rayleigh damping
            if (solverChoice.use_rayleigh_damping)
            {
                Real uu = rho_u(i,j,k) / cell_data(i,j,k,Rho_comp);
                rho_u_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (uu - dptr_rayleigh_ubar[k]) * cell_data(i,j,k,Rho_comp);
            }

            } // not on coarse-fine boundary
        });
        amrex::ParallelFor(tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // y-momentum equation
            rho_v_rhs(i, j, k) = 0.0; // Initialize the updated y-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (j == vlo_y && mlo_y(i,j,k)) || (j == vhi_y+1 && mhi_y(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

            // Add advective terms
            rho_v_rhs(i, j, k) += -AdvectionSrcForYMom(i, j, k, rho_u, rho_v, rho_w, v, z_nd, detJ,
                                                       dxInv, l_spatial_order, l_use_terrain);

            // Add diffusive terms
            amrex::Real diff_update;
            if (k < domhi_z && l_use_terrain) {
                diff_update = DiffusionSrcForMomWithTerrain(i, j, k, u, v, w, cell_data,
                                                            MomentumEqn::y, dxInv, K_turb, solverChoice,
                                                            z_nd, detJ, domain, bc_ptr);
            } else {
                diff_update = DiffusionSrcForMom(i, j, k, u, v, w, cell_data,
                                                 MomentumEqn::y, dxInv, K_turb, solverChoice,
                                                 domain, bc_ptr, er_arr);
            }
            rho_v_rhs(i, j, k) += diff_update;

            // Add pressure gradient
            amrex::Real gpy;
            if (l_use_terrain) {
                Real met_h_xi,met_h_eta,met_h_zeta;

                ComputeMetricAtJface(i,j,k,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::h_eta_zeta);

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

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_v_rhs(i, j, k) += -solverChoice.abl_pressure_grad[1];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                rho_v_rhs(i, j, k) += -solverChoice.coriolis_factor * rho_u_loc * solverChoice.sinphi;
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_v_rhs(i, j, k) += solverChoice.abl_geo_forcing[1];

            // Add Rayleigh damping
            if (solverChoice.use_rayleigh_damping)
            {
                Real vv = rho_v(i,j,k) / cell_data(i,j,k,Rho_comp);
                rho_v_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (vv - dptr_rayleigh_vbar[k]) * cell_data(i,j,k,Rho_comp);
            }

            } // not on coarse-fine boundary
        });
        amrex::ParallelFor(tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation
            rho_w_rhs(i, j, k) = 0.0; // Initialize the updated z-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (k == vlo_z && mlo_z(i,j,k)) || (k == vhi_z+1 && mhi_z(i,j,k)) );
            }

            if (!on_coarse_fine_boundary && k > 0)
            {
            // Add advective terms
            rho_w_rhs(i, j, k) += -AdvectionSrcForZMom(i, j, k, rho_u, rho_v, rho_w, w,
                                       z_nd, detJ, dxInv, l_spatial_order, l_use_terrain, domhi_z);

            // Add diffusive terms
            int k_diff = (k == domhi_z+1) ? domhi_z : k;
            amrex::Real diff_update;
            if (k < domhi_z && l_use_terrain) {
                diff_update = DiffusionSrcForMomWithTerrain(i, j, k_diff, u, v, w, cell_data,
                                       MomentumEqn::z, dxInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr);
            } else {
                diff_update = DiffusionSrcForMom(i, j, k_diff, u, v, w, cell_data,
                                                 MomentumEqn::z, dxInv, K_turb, solverChoice,
                                                 domain, bc_ptr, er_arr);
            }
            rho_w_rhs(i, j, k) += diff_update;

            // Add pressure gradient
            amrex::Real gpz;
            if (l_use_terrain) {
                Real met_h_xi,met_h_eta,met_h_zeta;
                ComputeMetricAtKface(i,j,k,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::h_zeta);
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
                int local_spatial_order = 2;
                rho_w_rhs(i, j, k) += grav_gpu[2] *
                     InterpolateDensityPertFromCellToFace(i, j, k, cell_data, rho_w(i,j,k),
                                                          Coord::z, local_spatial_order, r0_arr);
            }

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_w_rhs(i, j, k) += -solverChoice.abl_pressure_grad[2];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                rho_w_rhs(i, j, k) += solverChoice.coriolis_factor * rho_u_loc * solverChoice.cosphi;
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_w_rhs(i, j, k) += solverChoice.abl_geo_forcing[2];

            // Add Rayleigh damping
            if (solverChoice.use_rayleigh_damping)
            {
                rho_w_rhs(i, j, k) -= dptr_rayleigh_tau[k] * rho_w(i,j,k);
            }

            // Enforce no forcing term at bottom boundary
            if (k == 0) {
                rho_w_rhs(i,j,k) = 0.;
            } else if (k == domhi_z+1) {
                rho_w_rhs(i, j, k) = 0.;
            }
            } // not on coarse-fine boundary
        });
        } // not (rhs_vars == RHSVar::slow)
    }
}
