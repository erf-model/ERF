#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>

#include <TI_slow_headers.H>
#include <Utils.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in]  level level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[out]  S_rhs RHS computed here
 * @param[in]  S_old  old-time solution
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in]  S_scratch scratch space
 * @param[in]  xvel x-component of velocity
 * @param[in]  yvel y-component of velocity
 * @param[in]  zvel z-component of velocity
 * @param[in]  qv   water vapor
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] cc_src source terms for conserved variables
 * @param[in] xmom_src source terms for x-momentum
 * @param[in] ymom_src source terms for y-momentum
 * @param[in] zmom_src source terms for z-momentum
 * @param[in] Tau11 tau_11 component of stress tensor
 * @param[in] Tau22 tau_22 component of stress tensor
 * @param[in] Tau33 tau_33 component of stress tensor
 * @param[in] Tau12 tau_12 component of stress tensor
 * @param[in] Tau12 tau_13 component of stress tensor
 * @param[in] Tau21 tau_21 component of stress tensor
 * @param[in] Tau23 tau_23 component of stress tensor
 * @param[in] Tau31 tau_31 component of stress tensor
 * @param[in] Tau32 tau_32 component of stress tensor
 * @param[in] SmnSmn strain rate magnitude
 * @param[in] eddyDiffs diffusion coefficients for LES turbulence models
 * @param[in] Hfx3 heat flux in z-dir
 * @param[in] Diss dissipation of turbulent kinetic energy
 * @param[in]  geom   Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in]  most  Pointer to MOST class for Monin-Obukhov Similarity Theory boundary condition
 * @param[in]  domain_bcs_type_d device vector for domain boundary conditions
 * @param[in]  domain_bcs_type_h   host vector for domain boundary conditions
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] ax area fractions on x-faces
 * @param[in] ay area fractions on y-faces
 * @param[in] az area fractions on z-faces
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in]  p0     Reference (hydrostatically stratified) pressure
 * @param[in] mapfac_m map factor at cell centers
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 * @param[in] pp_inc   perturbational pressure for incompressible flows
 */

void erf_slow_rhs_inc (int level, int nrk,
                       Real dt,
                       Vector<MultiFab>& S_rhs,
                       Vector<MultiFab>& S_old,
                       Vector<MultiFab>& S_data,
                       const MultiFab& S_prim,
                       Vector<MultiFab>& S_scratch,
                       const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& zvel,
                       MultiFab& Omega,
                       const MultiFab& cc_src,
                       const MultiFab& xmom_src,
                       const MultiFab& ymom_src,
                       const MultiFab& zmom_src,
                       MultiFab* Tau11, MultiFab* Tau22, MultiFab* Tau33,
                       MultiFab* Tau12, MultiFab* Tau13, MultiFab* Tau21,
                       MultiFab* Tau23, MultiFab* Tau31, MultiFab* Tau32,
                       MultiFab* SmnSmn,
                       MultiFab* eddyDiffs,
                       MultiFab* Hfx3, MultiFab* Diss,
                       const Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<ABLMost>& most,
                       const Gpu::DeviceVector<BCRec>& domain_bcs_type_d,
                       const Vector<BCRec>& domain_bcs_type_h,
                       std::unique_ptr<MultiFab>& z_phys_nd,
                       std::unique_ptr<MultiFab>& ax,
                       std::unique_ptr<MultiFab>& ay,
                       std::unique_ptr<MultiFab>& az,
                       std::unique_ptr<MultiFab>& detJ,
                       const MultiFab* /*p0*/,
                       std::unique_ptr<MultiFab>& mapfac_m,
                       std::unique_ptr<MultiFab>& mapfac_u,
                       std::unique_ptr<MultiFab>& mapfac_v,
                       const MultiFab& pp_inc)
{
    BL_PROFILE_REGION("erf_slow_rhs_pre_inc()");

    const BCRec* bc_ptr_d = domain_bcs_type_d.data();
    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const MultiFab* t_mean_mf = nullptr;
    if (most) t_mean_mf = most->get_mac_avg(0,2);

    int start_comp = 0;
    int   num_comp = 2;
    int   end_comp = start_comp + num_comp - 1;

    const bool    l_const_rho      = solverChoice.constant_density;
    const AdvType l_horiz_adv_type = solverChoice.advChoice.dycore_horiz_adv_type;
    const AdvType l_vert_adv_type  = solverChoice.advChoice.dycore_vert_adv_type;
    const Real    l_horiz_upw_frac = solverChoice.advChoice.dycore_horiz_upw_frac;
    const Real    l_vert_upw_frac  = solverChoice.advChoice.dycore_vert_upw_frac;
    const bool    l_use_terrain    = solverChoice.use_terrain;

    AMREX_ALWAYS_ASSERT (!l_use_terrain);

    const bool l_use_diff       = ( (dc.molec_diff_type != MolecDiffType::None) ||
                                    (tc.les_type        !=       LESType::None) ||
                                    (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_constAlpha = ( dc.molec_diff_type == MolecDiffType::ConstantAlpha );
    const bool l_use_turb       = ( tc.les_type == LESType::Smagorinsky ||
                                    tc.les_type == LESType::Deardorff   ||
                                    tc.pbl_type == PBLType::MYNN25      ||
                                    tc.pbl_type == PBLType::YSU );

    const bool use_most     = (most != nullptr);
    const bool exp_most     = (solverChoice.use_explicit_most);

    const Box& domain = geom.Domain();
    const int domhi_z = domain.bigEnd(2);
    const int domlo_z = domain.smallEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Combine external forcing terms
    // *****************************************************************************
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *****************************************************************************
    // Pre-computed quantities
    // *****************************************************************************
    int nvars                     = S_data[IntVars::cons].nComp();
    const BoxArray& ba            = S_data[IntVars::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVars::cons].DistributionMap();

    std::unique_ptr<MultiFab> expr;
    std::unique_ptr<MultiFab> dflux_x;
    std::unique_ptr<MultiFab> dflux_y;
    std::unique_ptr<MultiFab> dflux_z;

    if (l_use_diff) {
        erf_make_tau_terms(level,nrk,domain_bcs_type_d,domain_bcs_type_h,z_phys_nd,
                           S_data,xvel,yvel,zvel,Omega,
                           Tau11,Tau22,Tau33,Tau12,Tau13,Tau21,Tau23,Tau31,Tau32,
                           SmnSmn,eddyDiffs,Hfx3,Diss,geom,solverChoice,most,
                           detJ,mapfac_m,mapfac_u,mapfac_v);

        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, nvars, 0);
    } // l_use_diff

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************

    // Open bc will be imposed upon all vars (we only access cons here for simplicity)
    const bool xlo_open = (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::open);
    const bool xhi_open = (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::open);
    const bool ylo_open = (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::open);
    const bool yhi_open = (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::open);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
    std::array<FArrayBox,AMREX_SPACEDIM> flux;

    for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        // If we are imposing open bc's then don't add rhs terms at the boundary locations
        if ( xlo_open && (tbx.smallEnd(0) == domain.smallEnd(0)) ) {tbx.growLo(0,-1);}
        if ( xhi_open && (tbx.bigEnd(0)   == domain.bigEnd(0)+1) ) {tbx.growHi(0,-1);}
        if ( ylo_open && (tby.smallEnd(1) == domain.smallEnd(1)) ) {tby.growLo(1,-1);}
        if ( yhi_open && (tby.bigEnd(1)   == domain.bigEnd(1)+1) ) {tby.growHi(1,-1);}

        // We don't compute a source term for z-momentum on the bottom or top boundary
        tbz.growLo(2,-1);
        tbz.growHi(2,-1);

        const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real> &       cell_rhs   = S_rhs[IntVars::cons].array(mfi);

        const Array4<const Real> & cell_data_old  = S_old[IntVars::cons].array(mfi);

        // We must initialize these to zero each RK step
        S_scratch[IntVars::xmom][mfi].template setVal<RunOn::Device>(0.,tbx);
        S_scratch[IntVars::ymom][mfi].template setVal<RunOn::Device>(0.,tby);
        S_scratch[IntVars::zmom][mfi].template setVal<RunOn::Device>(0.,tbz);

        Array4<Real> avg_xmom = S_scratch[IntVars::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVars::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVars::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<const Real>& rho_u = S_data[IntVars::xmom].array(mfi);
        const Array4<const Real>& rho_v = S_data[IntVars::ymom].array(mfi);
        const Array4<const Real>& rho_w = S_data[IntVars::zmom].array(mfi);

        const Array4<const Real>& rho_u_old = S_old[IntVars::xmom].array(mfi);
        const Array4<const Real>& rho_v_old = S_old[IntVars::ymom].array(mfi);
        const Array4<const Real>& rho_w_old = S_old[IntVars::zmom].array(mfi);

        const Array4<Real const>& xmom_src_arr   = xmom_src.const_array(mfi);
        const Array4<Real const>& ymom_src_arr   = ymom_src.const_array(mfi);
        const Array4<Real const>& zmom_src_arr   = zmom_src.const_array(mfi);

        // Perturbational pressure
        const Array4<const Real>& pp_arr = pp_inc.const_array(mfi);

        // Map factors
        const Array4<const Real>& mf_m   = mapfac_m->const_array(mfi);
        const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
        const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        const Array4<Real>& rho_u_rhs = S_rhs[IntVars::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVars::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVars::zmom].array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        // Terrain metrics
        const Array4<const Real>& z_nd     = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ_arr = detJ->const_array(mfi);

        // Base state
        // const Array4<const Real>& p0_arr = p0->const_array(mfi);

        // *****************************************************************************
        // *****************************************************************************
        // Define flux arrays for use in advection
        // *****************************************************************************
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flux[dir].resize(surroundingNodes(bx,dir),2);
            flux[dir].setVal<RunOn::Device>(0.);
        }
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>
            flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};

        // *****************************************************************************
        // Contravariant flux field
        // *****************************************************************************
        {
        BL_PROFILE("slow_rhs_making_omega");
            Box gbxo = surroundingNodes(bx,2); gbxo.grow(IntVect(1,1,0));
            // Now create Omega with momentum (not velocity)
            ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                omega_arr(i,j,k) = rho_w(i,j,k);
            });
        } // end profile


        // *****************************************************************************
        // Diffusive terms (pre-computed above)
        // *****************************************************************************
        // Expansion
        Array4<Real> er_arr;
        if (expr) {
            er_arr = expr->array(mfi);
        } else {
            er_arr = Array4<Real>{};
        }

        // No terrain diffusion
        Array4<Real> tau11,tau22,tau33;
        Array4<Real> tau12,tau13,tau23;
        if (Tau11) {
            tau11 = Tau11->array(mfi); tau22 = Tau22->array(mfi); tau33 = Tau33->array(mfi);
            tau12 = Tau12->array(mfi); tau13 = Tau13->array(mfi); tau23 = Tau23->array(mfi);
        } else {
            tau11 = Array4<Real>{}; tau22 = Array4<Real>{}; tau33 = Array4<Real>{};
            tau12 = Array4<Real>{}; tau13 = Array4<Real>{}; tau23 = Array4<Real>{};
        }
        // Terrain diffusion
        Array4<Real> tau21,tau31,tau32;
        if (Tau21) {
            tau21 = Tau21->array(mfi); tau31 = Tau31->array(mfi); tau32 = Tau32->array(mfi);
        } else {
            tau21 = Array4<Real>{}; tau31 = Array4<Real>{}; tau32 = Array4<Real>{};
        }

        // Strain magnitude
        Array4<Real> SmnSmn_a;
        if (tc.les_type == LESType::Deardorff) {
            SmnSmn_a = SmnSmn->array(mfi);
        } else {
            SmnSmn_a = Array4<Real>{};
        }

        // **************************************************************************
        // Define updates in the RHS of continuity and potential temperature equations
        // **************************************************************************
        auto const& ax_arr = ax->const_array(mfi);
        auto const& ay_arr = ay->const_array(mfi);
        auto const& az_arr = az->const_array(mfi);

        AdvectionSrcForRho(bx, cell_rhs,
                           rho_u, rho_v, omega_arr,      // these are being used to build the fluxes
                           avg_xmom, avg_ymom, avg_zmom, // these are being defined from the fluxes
                           ax_arr, ay_arr, az_arr, detJ_arr,
                           dxInv, mf_m, mf_u, mf_v,
                           flx_arr, l_const_rho);

        int icomp = RhoTheta_comp; int ncomp = 1;
        AdvectionSrcForScalars(bx, icomp, ncomp,
                               avg_xmom, avg_ymom, avg_zmom,
                               cell_prim, cell_rhs, detJ_arr,
                               dxInv, mf_m,
                               l_horiz_adv_type, l_vert_adv_type,
                               l_horiz_upw_frac, l_vert_upw_frac,
                               flx_arr, domain, bc_ptr_h);

        if (l_use_diff) {
            Array4<Real> diffflux_x = dflux_x->array(mfi);
            Array4<Real> diffflux_y = dflux_y->array(mfi);
            Array4<Real> diffflux_z = dflux_z->array(mfi);

            Array4<Real> hfx_z = Hfx3->array(mfi);
            Array4<Real> diss  = Diss->array(mfi);

            const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

            // NOTE: No diffusion for continuity, so n starts at 1.
            int n_start = amrex::max(start_comp,RhoTheta_comp);
            int n_comp  = end_comp - n_start + 1;

            DiffusionSrcForState_N(bx, domain, n_start, n_comp, exp_most, u, v,
                                   cell_data, cell_prim, cell_rhs,
                                   diffflux_x, diffflux_y, diffflux_z,
                                   dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                   hfx_z, diss,
                                   mu_turb, dc, tc,
                                   tm_arr, grav_gpu, bc_ptr_d, use_most);
        }

        // Add source terms for (rho theta)
        {
            auto const& src_arr = cc_src.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_rhs(i,j,k,RhoTheta_comp) += src_arr(i,j,k,RhoTheta_comp);
            });
        }

        // If in second RK stage, take average of old-time and new-time source
        if (nrk == 1)
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_rhs(i,j,k,     Rho_comp) *= 0.5;
                cell_rhs(i,j,k,RhoTheta_comp) *= 0.5;

                cell_rhs(i,j,k,     Rho_comp) += 0.5 / dt * (cell_data(i,j,k,     Rho_comp) - cell_data_old(i,j,k,     Rho_comp));
                cell_rhs(i,j,k,RhoTheta_comp) += 0.5 / dt * (cell_data(i,j,k,RhoTheta_comp) - cell_data_old(i,j,k,RhoTheta_comp));
            });
        }

        // *****************************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *****************************************************************************
        int lo_z_face;
        int hi_z_face;
        if (level == 0) {
            lo_z_face = domain.smallEnd(2);
            hi_z_face = domain.bigEnd(2)+1;
        } else {
            lo_z_face = mfi.validbox().smallEnd(2);
            hi_z_face = mfi.validbox().bigEnd(2)+1;
        }
        AdvectionSrcForMom(bx, tbx, tby, tbz,
                           rho_u_rhs, rho_v_rhs, rho_w_rhs,
                           cell_data, u, v, w,
                           rho_u, rho_v, omega_arr,
                           z_nd, ax_arr, ay_arr, az_arr, detJ_arr,
                           dxInv, mf_m, mf_u, mf_v,
                           l_horiz_adv_type, l_vert_adv_type,
                           l_horiz_upw_frac, l_vert_upw_frac,
                           l_use_terrain, lo_z_face, hi_z_face,
                           domain, bc_ptr_h);

        if (l_use_diff) {
            DiffusionSrcForMom_N(tbx, tby, tbz,
                                 rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                 tau11, tau22, tau33,
                                 tau12, tau13, tau23,
                                 dxInv,
                                 mf_m, mf_u, mf_v);
        }

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // x-momentum equation

            Real gpx = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));

            rho_u_rhs(i, j, k) += xmom_src_arr(i,j,k) - gpx - solverChoice.abl_pressure_grad[0];

            if (nrk == 1) {
              rho_u_rhs(i,j,k) *= 0.5;
              rho_u_rhs(i,j,k) += 0.5 / dt * (rho_u(i,j,k) - rho_u_old(i,j,k));
            }
        });

        ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // y-momentum equation

            Real gpy = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));

            rho_v_rhs(i, j, k) += ymom_src_arr(i,j,k) - gpy - solverChoice.abl_pressure_grad[1];

            if (nrk == 1) {
              rho_v_rhs(i,j,k) *= 0.5;
              rho_v_rhs(i,j,k) += 0.5 / dt * (rho_v(i,j,k) - rho_v_old(i,j,k));
            }
        });

        // *****************************************************************************
        // Zero out source terms for x- and y- momenta if at walls or inflow
        // We need to do this -- even though we call the boundary conditions later --
        // because the slow source is used to update the state in the fast interpolater.
        // *****************************************************************************
        if ( (bx.smallEnd(0) == domain.smallEnd(0)) &&
             (bc_ptr_h[BCVars::xvel_bc].lo(0) == ERFBCType::ext_dir) ) {
            Box lo_x_dom_face(bx); lo_x_dom_face.setBig(0,bx.smallEnd(0));
            ParallelFor(lo_x_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                rho_u_rhs(i,j,k) = 0.;
            });
        }
        if ( (bx.bigEnd(0) == domain.bigEnd(0)) &&
             (bc_ptr_h[BCVars::xvel_bc].hi(0) == ERFBCType::ext_dir) ) {
            Box hi_x_dom_face(bx); hi_x_dom_face.setSmall(0,bx.bigEnd(0)+1); hi_x_dom_face.setBig(0,bx.bigEnd(0)+1);
            ParallelFor(hi_x_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                rho_u_rhs(i,j,k) = 0.;
            });
        }
        if ( (bx.smallEnd(1) == domain.smallEnd(1)) &&
             (bc_ptr_h[BCVars::yvel_bc].lo(1) == ERFBCType::ext_dir) ) {
            Box lo_y_dom_face(bx); lo_y_dom_face.setBig(1,bx.smallEnd(1));
            ParallelFor(lo_y_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                rho_v_rhs(i,j,k) = 0.;
            });
        }
        if ( (bx.bigEnd(1) == domain.bigEnd(1)) &&
             (bc_ptr_h[BCVars::yvel_bc].hi(1) == ERFBCType::ext_dir) ) {
            Box hi_y_dom_face(bx); hi_y_dom_face.setSmall(1,bx.bigEnd(1)+1); hi_y_dom_face.setBig(1,bx.bigEnd(1)+1);;
            ParallelFor(hi_y_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                rho_v_rhs(i,j,k) = 0.;
            });
        }

        amrex::Box b2d = tbz;
        b2d.setSmall(2,0);
        b2d.setBig(2,0);
        // Enforce no forcing term at top and bottom boundaries
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rho_w_rhs(i,j,        0) = 0.;
            rho_w_rhs(i,j,domhi_z+1) = 0.;
        });

        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // z-momentum equation

            Real gpz = dxInv[2] * (pp_arr(i,j,k) - pp_arr(i,j,k-1));

            rho_w_rhs(i, j, k) += zmom_src_arr(i,j,k) - gpz - solverChoice.abl_pressure_grad[2];

            if (nrk == 1) {
                rho_w_rhs(i,j,k) *= 0.5;
                rho_w_rhs(i,j,k) += 0.5 / dt * (rho_w(i,j,k) - rho_w_old(i,j,k));
            }
        });
    } // mfi
    } // omp
}
