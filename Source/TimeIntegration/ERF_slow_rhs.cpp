#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <ERF_Constants.H>
#include <EddyViscosity.H>
#include <PBLModels.H>
#include <SpatialStencils.H>
#include <TimeIntegration.H>
#include <EOS.H>
#include <ERF.H>

using namespace amrex;

void erf_rhs (int level,
              Vector<MultiFab>& S_rhs,
              const Vector<MultiFab>& S_data,
              const MultiFab& S_prim,
              const MultiFab& xvel,
              const MultiFab& yvel,
              const MultiFab& zvel,
              MultiFab& source,
              std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
              std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
              const amrex::Geometry geom,
                    amrex::InterpFaceRegister* ifr,
              const SolverChoice& solverChoice,
              const ABLMost& most,
              const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
#ifdef ERF_USE_TERRAIN
              const MultiFab& z_phys_nd,
              const MultiFab& detJ_cc,
              const MultiFab& r0,
              const MultiFab& p0,
#else
              const amrex::Real* dptr_dens_hse, const amrex::Real* dptr_pres_hse,
#endif
              const amrex::Real* dptr_rayleigh_tau, const amrex::Real* dptr_rayleigh_ubar,
              const amrex::Real* dptr_rayleigh_vbar, const amrex::Real* dptr_rayleigh_thetabar)
{
    BL_PROFILE_VAR("erf_slow_rhs()",erf_slow_rhs);

    const int l_spatial_order = solverChoice.spatial_order;

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
    eddyDiffs.setVal(0.0);
    if (solverChoice.les_type == LESType::Smagorinsky ||
        solverChoice.les_type == LESType::Deardorff) {
        ComputeTurbulentViscosity(xvel, yvel, zvel, S_data[IntVar::cons],
                                  eddyDiffs, geom, solverChoice, domain_bcs_type_d);
    }
    if (solverChoice.pbl_type != PBLType::None) {
        ComputeTurbulentViscosityPBL(xvel, yvel, zvel, S_data[IntVar::cons],
                                     eddyDiffs, geom, solverChoice, most);
    }

    const iMultiFab *mlo_mf_x, *mhi_mf_x;
    const iMultiFab *mlo_mf_y, *mhi_mf_y;
    const iMultiFab *mlo_mf_z, *mhi_mf_z;

    bool l_use_QKE       = solverChoice.use_QKE;
    bool l_use_deardorff = (solverChoice.les_type == LESType::Deardorff);
    Real l_Delta         = std::pow(dx[0] * dx[1] * dx[2],1./3.);
    Real l_C_e           = solverChoice.Ce;

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
        const Array4<Real>& advflux_x = advflux[0].array(mfi);
        const Array4<Real>& advflux_y = advflux[1].array(mfi);
        const Array4<Real>& advflux_z = advflux[2].array(mfi);

        // These are temporaries we use to add to the S_rhs for the fluxes
        const Array4<Real>& diffflux_x = diffflux[0].array(mfi);
        const Array4<Real>& diffflux_y = diffflux[1].array(mfi);
        const Array4<Real>& diffflux_z = diffflux[2].array(mfi);

#ifdef ERF_USE_TERRAIN
        // These are metric terms for terrain-fitted coordiantes
        const Array4<const Real>& z_nd = z_phys_nd.const_array(mfi);
        const Array4<const Real>& detJ = detJ_cc.const_array(mfi);
#endif

        const Array4<Real>& K_turb = eddyDiffs.array(mfi);

#ifdef ERF_USE_TERRAIN
        const Array4<const Real>& r0_arr = r0.const_array(mfi);
        const Array4<const Real>& p0_arr = p0.const_array(mfi);
#endif
        const Box& gbx = mfi.growntilebox(1);
        const Array4<Real> & pp_arr  = pprime.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
#ifdef ERF_USE_TERRAIN
           pp_arr(i,j,k) = getPprimegivenRTh(cell_data(i,j,k,RhoTheta_comp),p0_arr(i,j,k));
#else
           pp_arr(i,j,k) = getPprimegivenRTh(cell_data(i,j,k,RhoTheta_comp),dptr_pres_hse[k]);
#endif
        });

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        int ncomp = S_data[IntVar::cons].nComp();

        int start_comp = 0;
        AdvectionSrcForState(bx, start_comp, ncomp, rho_u, rho_v, rho_w, cell_prim,
                             cell_rhs, advflux_x, advflux_y, advflux_z,
#ifdef ERF_USE_TERRAIN
                             z_nd, detJ,
#endif
                             dxInv, l_spatial_order, l_use_deardorff, l_use_QKE);

        amrex::ParallelFor(bx, ncomp,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {

            // Add diffusive terms.
            if (n == RhoTheta_comp)
                cell_rhs(i, j, k, n) += DiffusionSrcForState(i, j, k, cell_data, cell_prim, RhoTheta_comp,
                                        diffflux_x, diffflux_y, diffflux_z, dxInv, K_turb, solverChoice);
            if (n == RhoScalar_comp)
                cell_rhs(i, j, k, n) += DiffusionSrcForState(i, j, k, cell_data, cell_prim, RhoScalar_comp,
                                        diffflux_x, diffflux_y, diffflux_z, dxInv, K_turb, solverChoice);
            if (l_use_deardorff && n == RhoKE_comp)
                cell_rhs(i, j, k, n) += DiffusionSrcForState(i, j, k, cell_data, cell_prim, RhoKE_comp,
                                        diffflux_x, diffflux_y, diffflux_z, dxInv, K_turb, solverChoice);
            if (l_use_QKE && n == RhoQKE_comp)
                cell_rhs(i, j, k, n) += DiffusionSrcForState(i, j, k, cell_data, cell_prim, RhoQKE_comp,
                                        diffflux_x, diffflux_y, diffflux_z, dxInv, K_turb, solverChoice);

            // Add Rayleigh damping
            if (solverChoice.use_rayleigh_damping && n == RhoTheta_comp)
            {
                Real theta = cell_prim(i,j,k,PrimTheta_comp);
                cell_rhs(i, j, k, n) -= dptr_rayleigh_tau[k] * (theta - dptr_rayleigh_thetabar[k]) * cell_data(i,j,k,Rho_comp);
            }

            if (l_use_deardorff && n == RhoKE_comp)
            {
                // Add Buoyancy Source
                Real theta     = cell_prim(i,j,k,PrimTheta_comp);
                Real dtheta_dz = 0.5*(cell_prim(i,j,k+1,PrimTheta_comp)-cell_prim(i,j,k-1,PrimTheta_comp))*dxInv[2];
                Real E         = cell_prim(i,j,k,PrimKE_comp);
                Real length;
                if (dtheta_dz <= 0.) {
                   length = l_Delta;
                } else {
                   length = 0.76*std::sqrt(E)*(grav_gpu[2]/theta)*dtheta_dz;
                }
                Real KH   = 0.1 * (1.+2.*length/l_Delta) * std::sqrt(E);
                cell_rhs(i, j, k, n) += cell_data(i,j,k,Rho_comp) * grav_gpu[2] * KH * dtheta_dz;

                // Add TKE production
                cell_rhs(i, j, k, n) += ComputeTKEProduction(i,j,k,u,v,w,K_turb,dxInv,domain,bc_ptr);

                // Add dissipation
                if (std::abs(E) > 0.) {
                    cell_rhs(i, j, k, n) += cell_data(i,j,k,Rho_comp) * l_C_e *
                        std::pow(E,1.5) / length;
                }
            }

            // QKE : similar terms to TKE
            if (l_use_QKE && n == RhoQKE_comp) {
                cell_rhs(i,j,k,n) += ComputeQKESourceTerms(i,j,k,u,v,cell_data,cell_prim,K_turb,dxInv,domain,solverChoice,most);
            }

            // Add source terms. TODO: Put this under an if condition when we implement source term
            cell_rhs(i, j, k, n) += source_fab(i, j, k, n);
        });

        // Compute the RHS for the flux terms from this stage -- we do it this way so we don't double count
        //         fluxes at fine-fine interfaces
        amrex::ParallelFor(tbx, S_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
             xflux_rhs(i,j,k,n) = advflux_x(i,j,k,n) + diffflux_x(i,j,k,n);
        });
        amrex::ParallelFor(tby, S_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
             yflux_rhs(i,j,k,n) = advflux_y(i,j,k,n) + diffflux_y(i,j,k,n);
        });
        amrex::ParallelFor(tbz, S_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
             zflux_rhs(i,j,k,n) = advflux_z(i,j,k,n) + diffflux_z(i,j,k,n);
        });

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
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
            rho_u_rhs(i, j, k) += -AdvectionSrcForXMom(i, j, k, rho_u, rho_v, rho_w, u,
#ifdef ERF_USE_TERRAIN
                                                                z_nd, detJ,
#endif
                                                                dxInv, l_spatial_order);

            // Add diffusive terms
            rho_u_rhs(i, j, k) += DiffusionSrcForMom(i, j, k, u, v, w, cell_data,
                                                              MomentumEqn::x, dxInv, K_turb, solverChoice,
                                                              domain, bc_ptr);

            // Add pressure gradient
#ifdef ERF_USE_TERRAIN
            Real gp_xi = dxInv[0] *
                (pp_arr(i,j,k) - pp_arr(i-1,j,k));
            amrex::Real h_xi_on_iface = 0.125 * dxInv[0] * (
                z_nd(i+1,j,k) + z_nd(i+1,j,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
               -z_nd(i-1,j,k) - z_nd(i-1,j,k+1) - z_nd(i-1,j+1,k) - z_nd(i-1,j+1,k+1) );
            amrex::Real h_zeta_on_iface = 0.5 * dxInv[2] * (
                z_nd(i,j,k+1) + z_nd(i,j+1,k+1) - z_nd(i,j,k) - z_nd(i,j+1,k) );

            Real gp_zeta_on_iface = (k == 0) ?
                0.5 * dxInv[2] * (
                pp_arr(i,j,k+1) + pp_arr(i-1,j,k+1) - pp_arr(i,j,k) - pp_arr(i-1,j,k)):
                0.25 * dxInv[2] * (
                  pp_arr(i,j,k+1) + pp_arr(i-1,j,k+1) - pp_arr(i,j,k-1) - pp_arr(i-1,j,k-1));
            amrex::Real gpx = gp_xi - (h_xi_on_iface / h_zeta_on_iface) * gp_zeta_on_iface;
#else
            amrex::Real gpx = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));
#endif
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
            rho_v_rhs(i, j, k) += -AdvectionSrcForYMom(i, j, k, rho_u, rho_v, rho_w, v,
#ifdef ERF_USE_TERRAIN
                                                                z_nd, detJ,
#endif
                                                                dxInv, l_spatial_order);

            // Add diffusive terms
            rho_v_rhs(i, j, k) += DiffusionSrcForMom(i, j, k, u, v, w, cell_data,
                                                              MomentumEqn::y, dxInv, K_turb, solverChoice,
                                                              domain, bc_ptr);

            // Add pressure gradient
#ifdef ERF_USE_TERRAIN
            Real gp_eta = dxInv[1] *
                (pp_arr(i,j,k) - pp_arr(i,j-1,k));
            amrex::Real h_eta_on_jface = 0.125 * dxInv[1] * (
                z_nd(i,j+1,k) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k) + z_nd(i+1,j+1,k+1)
               -z_nd(i,j-1,k) - z_nd(i,j-1,k+1) - z_nd(i+1,j-1,k) - z_nd(i+1,j-1,k+1) );
            amrex::Real h_zeta_on_jface = 0.5 * dxInv[2] * (
                z_nd(i,j,k+1) + z_nd(i+1,j,k+1) - z_nd(i,j,k) - z_nd(i+1,j,k) );

            Real gp_zeta_on_jface = (k == 0) ?
                0.5 * dxInv[2] * (
                  pp_arr(i,j,k+1) + pp_arr(i,j-1,k+1) - pp_arr(i,j,k) - pp_arr(i,j-1,k)):
                0.25 * dxInv[2] * (
                  pp_arr(i,j,k+1) + pp_arr(i,j-1,k+1) - pp_arr(i,j,k-1) - pp_arr(i,j-1,k-1));
            amrex::Real gpy = gp_eta - (h_eta_on_jface / h_zeta_on_jface) * gp_zeta_on_jface;
#else
            amrex::Real gpy = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));
#endif
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

            if (!on_coarse_fine_boundary)
            {

            // Add advective terms
            rho_w_rhs(i, j, k) += -AdvectionSrcForZMom(i, j, k, rho_u, rho_v, rho_w, w,
#ifdef ERF_USE_TERRAIN
                                                                z_nd, detJ,
#endif
                                                                dxInv, l_spatial_order);

            // Add diffusive terms
            rho_w_rhs(i, j, k) += DiffusionSrcForMom(i, j, k, u, v, w, cell_data,
                                                              MomentumEqn::z, dxInv, K_turb, solverChoice,
                                                              domain, bc_ptr);

            // Add pressure gradient
#ifdef ERF_USE_TERRAIN
            amrex::Real h_zeta = 0.125 * dxInv[2] * (
                z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
               -z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
            amrex::Real gpz = dxInv[2] * (pp_arr(i,j,k) - pp_arr(i,j,k-1)) / h_zeta;
#else
            amrex::Real gpz = dxInv[2] * (pp_arr(i,j,k) - pp_arr(i,j,k-1));
#endif
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
#ifdef ERF_USE_TERRAIN
                     InterpolateDensityPertFromCellToFace(i, j, k, cell_data, rho_w(i,j,k),
                                                          Coord::z, local_spatial_order, r0_arr);
#else
                     InterpolateDensityPertFromCellToFace(i, j, k, cell_data, rho_w(i,j,k),
                                                          Coord::z, local_spatial_order, dptr_dens_hse);
#endif
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
                //rho_w_rhs(i, j, k) = rho_w_rhs(i,j,k-1);
                rho_w_rhs(i, j, k) = 0.;
            }

            } // not on coarse-fine boundary
        });
    }
}
