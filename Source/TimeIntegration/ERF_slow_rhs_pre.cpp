#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>

#include <TI_slow_headers.H>
#include <EOS.H>
#include <Utils.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in]  level level of resolution
 * @param[in]  finest_level finest level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[out]  S_rhs RHS computed here
 * @param[in]  S_old  old-time solution -- used only for incompressible
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in]  S_scratch scratch space
 * @param[in]  xvel x-component of velocity
 * @param[in]  yvel y-component of velocity
 * @param[in]  zvel z-component of velocity
 * @param[in]  qv   water vapor
 * @param[in]  z_t_ mf rate of change of grid height -- only relevant for moving terrain
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
 * @param[in] pp_inc  Perturbational pressure only used for incompressible flow
 * @param[in] mapfac_m map factor at cell centers
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 */

void erf_slow_rhs_pre (int level, int finest_level,
                       int nrk,
                       Real dt,
                       Vector<MultiFab>& S_rhs,
                       Vector<MultiFab>& S_old,
                       Vector<MultiFab>& S_data,
                       const MultiFab& S_prim,
                       Vector<MultiFab>& S_scratch,
                       const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& zvel,
                       std::unique_ptr<MultiFab>& z_t_mf,
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
                       const MultiFab* p0,
#ifdef ERF_USE_POISSON_SOLVE
                       const MultiFab& pp_inc,
#endif
                       std::unique_ptr<MultiFab>& mapfac_m,
                       std::unique_ptr<MultiFab>& mapfac_u,
                       std::unique_ptr<MultiFab>& mapfac_v,
#ifdef ERF_USE_EB
                       EBFArrayBoxFactory const& ebfact,
#endif
                       YAFluxRegister* fr_as_crse,
                       YAFluxRegister* fr_as_fine)
{
    BL_PROFILE_REGION("erf_slow_rhs_pre()");

#ifdef ERF_USE_EB
    amrex::ignore_unused(ax,ay,az,detJ);
#endif

    const BCRec* bc_ptr_d = domain_bcs_type_d.data();
    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const MultiFab* t_mean_mf = nullptr;
    if (most) t_mean_mf = most->get_mac_avg(0,2);

    int start_comp = 0;
    int   num_comp = 2;
    int   end_comp = start_comp + num_comp - 1;

    const AdvType l_horiz_adv_type = solverChoice.advChoice.dycore_horiz_adv_type;
    const AdvType l_vert_adv_type  = solverChoice.advChoice.dycore_vert_adv_type;
    const Real    l_horiz_upw_frac = solverChoice.advChoice.dycore_horiz_upw_frac;
    const Real    l_vert_upw_frac  = solverChoice.advChoice.dycore_vert_upw_frac;
    const bool    l_use_terrain    = solverChoice.use_terrain;
    const bool    l_moving_terrain = (solverChoice.terrain_type == TerrainType::Moving);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT (l_use_terrain);

    const bool l_reflux = (solverChoice.coupling_type == CouplingType::TwoWay);

    const bool l_use_diff       = ( (dc.molec_diff_type != MolecDiffType::None) ||
                                    (tc.les_type        !=       LESType::None) ||
                                    (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_turb       = ( tc.les_type == LESType::Smagorinsky ||
                                    tc.les_type == LESType::Deardorff   ||
                                    tc.pbl_type == PBLType::MYNN25      ||
                                    tc.pbl_type == PBLType::YSU );

    const bool l_use_moisture = (solverChoice.moisture_type != MoistureType::None);
    const bool l_use_most     = (most != nullptr);
    const bool l_exp_most     = (solverChoice.use_explicit_most);

#ifdef ERF_USE_POISSON_SOLVE
    const bool l_incompressible = solverChoice.incompressible[level];
    const bool l_const_rho      = solverChoice.constant_density;

    // We cannot use incompressible with terrain or with moisture
    AMREX_ALWAYS_ASSERT(!l_use_terrain  || !l_incompressible);
    AMREX_ALWAYS_ASSERT(!l_use_moisture || !l_incompressible);
#else
    const bool l_incompressible = false;
    const bool l_const_rho      = false;
#endif

    const Box& domain = geom.Domain();
    const int domhi_z = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const Real* dx = geom.CellSize();

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
        erf_make_tau_terms(level,nrk,domain_bcs_type_h,z_phys_nd,
                           S_data,xvel,yvel,zvel,Omega,
                           Tau11,Tau22,Tau33,Tau12,Tau13,Tau21,Tau23,Tau31,Tau32,
                           SmnSmn,eddyDiffs,geom,solverChoice,most,
                           detJ,mapfac_m,mapfac_u,mapfac_v);

        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, nvars, 0);
    } // l_use_diff

    // *****************************************************************************
    // Define updates and fluxes in the current RK stage
    // *****************************************************************************

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

        // We don't compute a source term for z-momentum on the bottom or top boundary
        tbz.growLo(2,-1);
        tbz.growHi(2,-1);

        const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real>       & cell_rhs   = S_rhs[IntVars::cons].array(mfi);

        const Array4<const Real> & cell_old   = S_old[IntVars::cons].array(mfi);

        const Array4<Real const>& xmom_src_arr   = xmom_src.const_array(mfi);
        const Array4<Real const>& ymom_src_arr   = ymom_src.const_array(mfi);
        const Array4<Real const>& zmom_src_arr   = zmom_src.const_array(mfi);

        const Array4<Real>& rho_u_old = S_old[IntVars::xmom].array(mfi);
        const Array4<Real>& rho_v_old = S_old[IntVars::ymom].array(mfi);

        if (l_incompressible) {
            // When incompressible we must reset these to 0 each RK step
            S_scratch[IntVars::xmom][mfi].template setVal<RunOn::Device>(0.0,tbx);
            S_scratch[IntVars::ymom][mfi].template setVal<RunOn::Device>(0.0,tby);
            S_scratch[IntVars::zmom][mfi].template setVal<RunOn::Device>(0.0,tbz);
        }

        Array4<Real> avg_xmom = S_scratch[IntVars::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVars::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVars::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<const Real>& rho_u = S_data[IntVars::xmom].array(mfi);
        const Array4<const Real>& rho_v = S_data[IntVars::ymom].array(mfi);
        const Array4<const Real>& rho_w = S_data[IntVars::zmom].array(mfi);

        // Map factors
        const Array4<const Real>& mf_m   = mapfac_m->const_array(mfi);
        const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
        const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        Array4<const Real> z_t;
        if (z_t_mf)
            z_t = z_t_mf->array(mfi);
        else
            z_t = Array4<const Real>{};

        const Array4<Real>& rho_u_rhs = S_rhs[IntVars::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVars::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVars::zmom].array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        // Terrain metrics
        const Array4<const Real>& z_nd     = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};

        // Base state
        const Array4<const Real>& p0_arr = p0->const_array(mfi);

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
        // Perturbational pressure field
        // *****************************************************************************
        FArrayBox pprime;
        if (!l_incompressible) {
            Box gbx = mfi.tilebox(); gbx.grow(IntVect(1,1,1));
            if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);
            pprime.resize(gbx,1,The_Async_Arena());
            const Array4<Real>& pptemp_arr = pprime.array();
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                //if (cell_data(i,j,k,RhoTheta_comp) < 0.) printf("BAD THETA AT %d %d %d %e %e \n",
                //    i,j,k,cell_data(i,j,k,RhoTheta_comp),cell_data(i,j,k+1,RhoTheta_comp));
                AMREX_ASSERT(cell_data(i,j,k,RhoTheta_comp) > 0.);
                Real qv_for_p = (l_use_moisture) ? cell_data(i,j,k,RhoQ1_comp)/cell_data(i,j,k,Rho_comp) : 0.0;
                pptemp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p) - p0_arr(i,j,k);
            });
        }
#ifdef ERF_USE_POISSON_SOLVE
        const Array4<const Real>& pp_arr = (l_incompressible) ? pp_inc.const_array(mfi) : pprime.const_array();
#else
        const Array4<const Real>& pp_arr = pprime.const_array();
#endif

        // *****************************************************************************
        // Contravariant flux field
        // *****************************************************************************
        {
        BL_PROFILE("slow_rhs_making_omega");
            Box gbxo = surroundingNodes(bx,2); gbxo.grow(IntVect(1,1,0));
            // Now create Omega with momentum (not velocity) with z_t subtracted if moving terrain
            if (l_use_terrain) {

                Box gbxo_lo = gbxo; gbxo_lo.setBig(2,0);
                if (gbxo_lo.smallEnd(2) <= 0) {
                    ParallelFor(gbxo_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = 0.;
                    });
                }
                Box gbxo_hi = gbxo; gbxo_hi.setSmall(2,gbxo.bigEnd(2));
                int hi_z_face = domain.bigEnd(2)+1;
                if (gbxo_hi.bigEnd(2) >= hi_z_face) {
                    ParallelFor(gbxo_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = rho_w(i,j,k);
                    });
                }

                if (z_t) {
                    Box gbxo_mid = gbxo; gbxo_mid.setSmall(2,1); gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
                    ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        // We define rho on the z-face the same way as in MomentumToVelocity/VelocityToMomentum
                        Real rho_at_face = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp));
                        omega_arr(i,j,k) = OmegaFromW(i,j,k,rho_w(i,j,k),rho_u,rho_v,z_nd,dxInv) -
                            rho_at_face * z_t(i,j,k);
                    });
                } else {
                    Box gbxo_mid = gbxo;
                    if (gbxo_mid.smallEnd(2) <= 0) {
                        gbxo_mid.setSmall(2,1);
                    }
                    if (gbxo_mid.bigEnd(2) >= domain.bigEnd(2)+1) {
                        gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
                    }
                    ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = OmegaFromW(i,j,k,rho_w(i,j,k),rho_u,rho_v,z_nd,dxInv);
                    });
                }
            } else {
                ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = rho_w(i,j,k);
                });
            }
        } // end profile


        // *****************************************************************************
        // Diffusive terms (pre-computed above)
        // *****************************************************************************
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

        // *****************************************************************************
        // Define updates in the RHS of continuity and potential temperature equations
        // *****************************************************************************
#ifdef ERF_USE_EB
        auto const& ax_arr   = ebfact.getAreaFrac()[0]->const_array(mfi);
        auto const& ay_arr   = ebfact.getAreaFrac()[1]->const_array(mfi);
        auto const& az_arr   = ebfact.getAreaFrac()[2]->const_array(mfi);
        const auto& detJ_arr = ebfact.getVolFrac().const_array(mfi);
#else
        auto const& ax_arr   = ax->const_array(mfi);
        auto const& ay_arr   = ay->const_array(mfi);
        auto const& az_arr   = az->const_array(mfi);
        auto const& detJ_arr = detJ->const_array(mfi);
#endif

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

            if (l_use_terrain) {
                DiffusionSrcForState_T(bx, domain, n_start, n_comp, l_exp_most, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       z_nd, ax_arr, ay_arr, az_arr, detJ_arr,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       hfx_z, diss, mu_turb, dc, tc,
                                       tm_arr, grav_gpu, bc_ptr_d, l_use_most);
            } else {
                DiffusionSrcForState_N(bx, domain, n_start, n_comp, l_exp_most, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       hfx_z, diss,
                                       mu_turb, dc, tc,
                                       tm_arr, grav_gpu, bc_ptr_d, l_use_most);
            }
        }

        const Array4<Real const>& source_arr   = cc_src.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cell_rhs(i,j,k,Rho_comp)      += source_arr(i,j,k,Rho_comp);
            cell_rhs(i,j,k,RhoTheta_comp) += source_arr(i,j,k,RhoTheta_comp);
        });

        // Multiply the slow RHS for rho and rhotheta by detJ here so we don't have to later
        if (l_use_terrain && l_moving_terrain) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_rhs(i,j,k,Rho_comp)      *= detJ_arr(i,j,k);
                cell_rhs(i,j,k,RhoTheta_comp) *= detJ_arr(i,j,k);
            });
        }

        // If incompressible and in second RK stage, take average of old-time and new-time source
        if ( l_incompressible && (nrk == 1) )
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_rhs(i,j,k,     Rho_comp) *= 0.5;
                cell_rhs(i,j,k,RhoTheta_comp) *= 0.5;

                cell_rhs(i,j,k,     Rho_comp) += 0.5 / dt * (cell_data(i,j,k,     Rho_comp) - cell_old(i,j,k,     Rho_comp));
                cell_rhs(i,j,k,RhoTheta_comp) += 0.5 / dt * (cell_data(i,j,k,RhoTheta_comp) - cell_old(i,j,k,RhoTheta_comp));
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
            // Note: tau** were calculated with calls to
            // ComputeStress[Cons|Var]Visc_[N|T] in which ConsVisc ("constant
            // viscosity") means that there is no contribution from a
            // turbulence model. However, whether this field truly is constant
            // depends on whether MolecDiffType is Constant or ConstantAlpha.
            if (l_use_terrain) {
                DiffusionSrcForMom_T(tbx, tby, tbz,
                                     rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                     tau11, tau22, tau33,
                                     tau12, tau13,
                                     tau21, tau23,
                                     tau31, tau32,
                                     detJ_arr, dxInv,
                                     mf_m, mf_u, mf_v);
            } else {
                DiffusionSrcForMom_N(tbx, tby, tbz,
                                     rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                     tau11, tau22, tau33,
                                     tau12, tau13, tau23,
                                     dxInv,
                                     mf_m, mf_u, mf_v);
            }
        }

        auto abl_pressure_grad    = solverChoice.abl_pressure_grad;

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // x-momentum equation

            //Note : mx/my == 1, so no map factor needed here
            Real gp_xi = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));
            Real gpx = gp_xi;

            if (l_use_terrain) {
                Real met_h_xi   = Compute_h_xi_AtIface  (i, j, k, dxInv, z_nd);
                Real met_h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);

                Real gp_zeta_on_iface;
                if (k==0) {
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
                gpx -= (met_h_xi/ met_h_zeta) * gp_zeta_on_iface;
            }

            gpx *= mf_u(i,j,0);

            Real q = 0.0;
            if (l_use_moisture) {
                q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i-1,j,k,PrimQ1_comp)
                           +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i-1,j,k,PrimQ2_comp) );
            }

            rho_u_rhs(i, j, k) += (-gpx - abl_pressure_grad[0]) / (1.0 + q)
                                  + xmom_src_arr(i,j,k);

            if (l_moving_terrain) {
                Real h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);
                rho_u_rhs(i, j, k) *= h_zeta;
            }

            if ( l_incompressible && (nrk == 1) ) {
              rho_u_rhs(i,j,k) *= 0.5;
              rho_u_rhs(i,j,k) += 0.5 / dt * (rho_u(i,j,k) - rho_u_old(i,j,k));
            }
        });

        ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // y-momentum equation

            //Note : mx/my == 1, so no map factor needed here
            Real gp_eta = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));
            Real gpy = gp_eta;

            if (l_use_terrain) {
                Real met_h_eta  = Compute_h_eta_AtJface (i, j, k, dxInv, z_nd);
                Real met_h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
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
                gpy -= (met_h_eta / met_h_zeta) * gp_zeta_on_jface;
            } // l_use_terrain

            gpy *= mf_v(i,j,0);

            Real q = 0.0;
            if (l_use_moisture) {
                q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i,j-1,k,PrimQ1_comp)
                           +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i,j-1,k,PrimQ2_comp) );
            }

            rho_v_rhs(i, j, k) += (-gpy - abl_pressure_grad[1]) / (1.0_rt + q) + ymom_src_arr(i,j,k);

            if (l_moving_terrain) {
                Real h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
                rho_v_rhs(i, j, k) *= h_zeta;
            }

            if ( l_incompressible && (nrk == 1) ) {
              rho_v_rhs(i,j,k) *= 0.5;
              rho_v_rhs(i,j,k) += 0.5 / dt * (rho_v(i,j,k) - rho_v_old(i,j,k));
            }
        });

        // *****************************************************************************
        // Zero out source terms for x- and y- momenta if at walls or inflow
        // We need to do this -- even though we call the boundary conditions later --
        // because the slow source is used to update the state in the fast interpolater.
        // *****************************************************************************
        {
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
        }

        // *****************************************************************************
        // Zero out source term for z-momentum at top and bottom of grid
        // *****************************************************************************
        Box b2d = tbz;
        b2d.setSmall(2,lo_z_face);
        b2d.setBig(2,lo_z_face);
        // Enforce no forcing term at top and bottom boundaries of this grid
        // We do this even when not at top or bottom of the domain because
        //    z-vel on the coarse/fine boundary is given by the coarse value
        //    (suitably interpolated tangentially and in time)
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rho_w_rhs(i,j,lo_z_face) = 0.;
            rho_w_rhs(i,j,hi_z_face) = 0.;
        });

        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // z-momentum equation

            Real met_h_zeta = (l_use_terrain) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd) : 1;
            Real gpz = dxInv[2] * ( pp_arr(i,j,k)-pp_arr(i,j,k-1) )  / met_h_zeta;

            Real q = 0.0;
            if (l_use_moisture) {
                q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i,j,k-1,PrimQ1_comp)
                           +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i,j,k-1,PrimQ2_comp) );
            }
            rho_w_rhs(i, j, k) += (zmom_src_arr(i,j,k) - gpz - abl_pressure_grad[2]) / (1.0_rt + q);

            if (l_use_terrain && l_moving_terrain) {
                 rho_w_rhs(i, j, k) *= 0.5 * (detJ_arr(i,j,k) + detJ_arr(i,j,k-1));
            }
        });

        {
        BL_PROFILE("slow_rhs_pre_fluxreg");
        // We only add to the flux registers in the final RK step
        // NOTE: for now we are only refluxing density not (rho theta) since the latter seems to introduce
        //       a problem at top and bottom boundaries
        if (l_reflux && nrk == 2) {
            int strt_comp_reflux = (l_const_rho) ? 1 : 0;
            int  num_comp_reflux = 1;
            if (level < finest_level) {
                fr_as_crse->CrseAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dt, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, RunOn::Device);
            }
            if (level > 0) {
                fr_as_fine->FineAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dt, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, RunOn::Device);
            }
        } // two-way coupling
        } // end profile
    } // mfi
    } // OMP
}
