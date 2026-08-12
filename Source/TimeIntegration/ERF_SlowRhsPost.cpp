#include <AMReX.H>
#include <ERF_SrcHeaders.H>
#include <ERF_TI_slow_headers.H>
#include <ERF_ShocDriver.H>
#include <ERF_EBAdvection.H>
#include <ERF_EBRedistribute.H>
#include "ERF_ResolvedWallFlux.H"
#include "Prob/ERF_CloudChamberBudget.H"

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the scalars other than density or potential temperature
 *
 * @param[in   ] evel level of resolution
 * @param[in   ] finest_level finest level of resolution
 * @param[in   ] nrk   which RK stage
 * @param[in   ] dt    slow time step
 * @param[  out] S_rhs RHS computed here
 * @param[in   ] S_old solution at start of time step
 * @param[in   ] S_new solution at end of current RK stage
 * @param[in   ] S_data current solution
 * @param[in   ] S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in   ] avg_xmom
 * @param[in   ] avg_ymom
 * @param[in   ] avg_zmom
 * @param[in   ] xvel x-component of velocity
 * @param[in   ] yvel y-component of velocity
 * @param[in   ] zvel z-component of velocity
 * @param[in   ] source source terms for conserved variables
 * @param[in   ] terrain_blank immersed forcing mask
 * @param[in   ] SmnSmn strain rate magnitude
 * @param[in   ] eddyDiffs diffusion coefficients for LES turbulence models
 * @param[in   ] Hfx3 heat flux in z-dir
 * @param[in   ] Diss dissipation of turbulent kinetic energy
 * @param[in   ] geom   Container for geometric information
 * @param[in   ] solverChoice  Container for solver parameters
 * @param[in   ] SurfLayer  Pointer to SurfaceLayer class for Monin-Obukhov Similarity Theory boundary condition
 * @param[in   ] domain_bcs_type_d device vector for domain boundary conditions
 * @param[in   ] z_phys_nd height coordinate at nodes
 * @param[in   ] ax area fractions on x-faces
 * @param[in   ] ay area fractions on y-faces
 * @param[in   ] az area fractions on z-faces
 * @param[in   ] detJ     Jacobian of the metric transformation at start of time step (= 1 if use_terrain is false)
 * @param[in   ] detJ_new Jacobian of the metric transformation at new RK stage time (= 1 if use_terrain is false)
 * @param[in   ] mapfac map factors
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 */

void erf_slow_rhs_post (int level, int finest_level,
                        int nrk,
                        double dt_d,
                        int n_qstate,
                        Vector<MultiFab>& S_rhs,
                        Vector<MultiFab>& S_old,
                        Vector<MultiFab>& S_new,
                        Vector<MultiFab>& S_data,
                        const MultiFab& S_prim,
                              MultiFab& avg_xmom,
                              MultiFab& avg_ymom,
                              MultiFab& avg_zmom,
                        const MultiFab& xvel,
                        const MultiFab& yvel,
                        const MultiFab& /*zvel*/,
                        const MultiFab& source,
                              MultiFab* terrain_blank,
                              MultiFab* terrain_blank_xface,
                              MultiFab* terrain_blank_yface,
                              MultiFab* terrain_blank_zface,
                        const MultiFab* SmnSmn,
                        const MultiFab* eddyDiffs,
                        MultiFab* Hfx1, MultiFab* Hfx2, MultiFab* Hfx3,
                        MultiFab* Q1fx1, MultiFab* Q1fx2,
                        MultiFab* Q1fx3, MultiFab* Q2fx3,
                        MultiFab* Diss,
                        const Geometry geom,
                        const SolverChoice& solverChoice,
                        std::unique_ptr<SurfaceLayer>& SurfLayer,
                        const Gpu::DeviceVector<BCRec>& domain_bcs_type_d,
                        const Vector<BCRec>& domain_bcs_type_h,
                        std::unique_ptr<MultiFab>& z_phys_nd,
                        std::unique_ptr<MultiFab>& z_phys_cc,
                        std::unique_ptr<MultiFab>& ax,
                        std::unique_ptr<MultiFab>& ay,
                        std::unique_ptr<MultiFab>& az,
                        std::unique_ptr<MultiFab>& detJ,
                        MultiFab* detJ_new,
                        Gpu::DeviceVector<Real>& stretched_dz_d,
                        Vector<std::unique_ptr<MultiFab>>& mapfac,
                        amrex::EBFArrayBoxFactory const& ebfact,
#ifdef ERF_USE_EAMXX_SHOC
                        SHOCInterface* eamxx_shoc_lev,
#endif
                        ShocDriver* native_shoc_lev,
                        YAFluxRegister* fr_as_crse,
                        YAFluxRegister* fr_as_fine,
                        std::unique_ptr<ReadBndryPlanes>& m_r2d,
                        const MultiFab* cloud_chamber_base_state,
                        const erf_cloud_chamber::Config* cloud_chamber_config,
                        CloudChamberBudget* cloud_budget)
{
    BL_PROFILE_REGION("erf_slow_rhs_post()");

    Real dt = static_cast<Real>(dt_d);

    const BCRec* bc_ptr_d = domain_bcs_type_d.data();
    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    AdvChoice  ac = solverChoice.advChoice;
    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const MultiFab*  t_mean_mf = nullptr;
    if (SurfLayer) { t_mean_mf = SurfLayer->get_mac_avg(level,2); }

    const bool l_use_terrain      = (solverChoice.mesh_type != MeshType::ConstantDz);
    const bool l_moving_terrain   = (solverChoice.terrain_type == TerrainType::MovingFittedMesh);
    const bool l_reflux = ( (solverChoice.coupling_type == CouplingType::TwoWay) && (nrk == 2) && (finest_level > 0) );
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT(l_use_terrain);

    const bool l_anelastic   = solverChoice.anelastic[level];

    const bool l_use_KE         = ( tc.use_tke );
    const bool l_need_SmnSmn    = ( tc.les_type  == LESType::Deardorff ||
                                    tc.rans_type == RANSType::kEqn );
    const bool l_advect_KE      = ( tc.use_tke && tc.advect_tke );
    const bool l_use_diff       = ((dc.molec_diff_type != MolecDiffType::None) ||
                                   (tc.les_type        !=       LESType::None) ||
                                   (tc.rans_type       !=      RANSType::None) ||
                                   (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_turb       = tc.use_kturb;
    const bool l_rotate         = (solverChoice.use_rotate_surface_flux);
    const bool l_do_scalar      = (solverChoice.transport_scalar);
    const bool l_use_eb         = (solverChoice.terrain_type == TerrainType::EB);

    amrex::ignore_unused(m_r2d);

    const Box& domain = geom.Domain();

    bool l_apply_surface_layer_fluxes_in_diffusion = (SurfLayer != nullptr);
#ifdef ERF_USE_EAMXX_SHOC
    if (tc.uses_eamxx_shoc()) {
        l_apply_surface_layer_fluxes_in_diffusion = false;
    }
#endif
    if (tc.uses_native_shoc()) {
        AMREX_ALWAYS_ASSERT(native_shoc_lev != nullptr);
        l_apply_surface_layer_fluxes_in_diffusion =
            l_apply_surface_layer_fluxes_in_diffusion &&
            native_shoc_lev->uses_host_diffusion();
    }

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const Real* dx = geom.CellSize();

    // *************************************************************************
    // Set gravity as a vector
    // *************************************************************************
    const    Array<Real,AMREX_SPACEDIM> grav{zero, zero, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *************************************************************************
    // Pre-computed quantities
    // *************************************************************************
    int nvars                     = S_data[IntVars::cons].nComp();
    const BoxArray& ba            = S_data[IntVars::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVars::cons].DistributionMap();
    const bool use_physical_chamber_wall_flux =
        cloud_chamber_config != nullptr && cloud_chamber_base_state != nullptr &&
        cloud_chamber_config->physical_initialization;
    const erf_wall_thermodynamics::Boundary chamber_walls =
        use_physical_chamber_wall_flux ? cloud_chamber_config->wall_boundary() :
                                         erf_wall_thermodynamics::Boundary{};

    std::unique_ptr<MultiFab> dflux_x;
    std::unique_ptr<MultiFab> dflux_y;
    std::unique_ptr<MultiFab> dflux_z;

    if (l_use_diff) {
        IntVect ng(0,0,1);
        // The physical chamber needs one persistent component for each moist
        // state even when budgets are disabled: qv and qc are corrected and
        // retained independently.  All other configurations retain ERF's
        // established one-component reusable diffusion storage.
        const int n_flux_components = use_physical_chamber_wall_flux ?
            std::max(1, n_qstate) : 1;
        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, n_flux_components, ng);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, n_flux_components, ng);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, n_flux_components, 0);
        // Every physical wall override reads the old face flux before
        // replacing it; make that read deterministic regardless of budget
        // diagnostics.
        if (use_physical_chamber_wall_flux) {
            dflux_x->setVal(0.0);
            dflux_y->setVal(0.0);
            dflux_z->setVal(0.0);
        }
    } else {
        dflux_x = nullptr;
        dflux_y = nullptr;
        dflux_z = nullptr;
    }

    // Valid vars
    Vector<int> is_valid_slow_var; is_valid_slow_var.resize(RhoQ1_comp+1,0);
    if (l_use_KE)    { is_valid_slow_var[    RhoKE_comp] = 1; }
    if (l_do_scalar) { is_valid_slow_var[RhoScalar_comp] = 1; }
    if (solverChoice.moisture_type != MoistureType::None) {
         is_valid_slow_var[RhoQ1_comp] = 1;
    }

    // *************************************************************************
    // Calculate cell-centered eddy viscosity & diffusivities
    //
    // Notes -- we fill all the data in ghost cells before calling this so
    //    that we can fill the eddy viscosity in the ghost regions and
    //    not have to call a boundary filler on this data itself
    //
    // LES - updates both horizontal and vertical eddy viscosityS_tmp components
    // PBL - only updates vertical eddy viscosity components so horizontal
    //       components come from the LES model or are left as zero.
    // *************************************************************************

    // EB Anelastic: Copy projected momentum with ghost-cell synchronization
    if (l_anelastic && l_use_eb) {
        avg_xmom.ParallelCopy(S_data[IntVars::xmom], 0, 0, 1, 0, 1, geom.periodicity());
        avg_ymom.ParallelCopy(S_data[IntVars::ymom], 0, 0, 1, 0, 1, geom.periodicity());
        avg_zmom.ParallelCopy(S_data[IntVars::zmom], 0, 0, 1, 0, 1, geom.periodicity());
    }

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      std::array<FArrayBox,AMREX_SPACEDIM> flux;

      int start_comp;
      int   num_comp;

      // Cell-centered masks for EB (used for flux interpolation)
      iMultiFab physbnd_mask;
      bool already_on_centroids = false;
      if (l_use_eb) {
          physbnd_mask.define(S_data[IntVars::cons].boxArray(), S_data[IntVars::cons].DistributionMap(), 1, 1);
          physbnd_mask.BuildMask(geom.Domain(), geom.periodicity(), 1, 1, 0, 1);
      }

      for (MFIter mfi(S_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box tbx  = mfi.tilebox();

        // *************************************************************************
        // Define flux arrays for use in advection
        // *************************************************************************
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (!l_use_eb) {
                flux[dir].resize(surroundingNodes(tbx,dir),nvars,The_Async_Arena());
            } else {
                flux[dir].resize(surroundingNodes(tbx,dir).grow(1),nvars,The_Async_Arena());
            }
            flux[dir].setVal<RunOn::Device>(0);
        }
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>
            flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};

        // *************************************************************************
        // Define Array4's
        // *************************************************************************
        const Array4<const Real> & old_cons   = S_old[IntVars::cons].array(mfi);
        const Array4<      Real> & cell_rhs   = S_rhs[IntVars::cons].array(mfi);

        const Array4<      Real> & new_cons  = S_new[IntVars::cons].array(mfi);
        const Array4<      Real> & new_xmom  = S_new[IntVars::xmom].array(mfi);
        const Array4<      Real> & new_ymom  = S_new[IntVars::ymom].array(mfi);
        const Array4<      Real> & new_zmom  = S_new[IntVars::zmom].array(mfi);

        const Array4<      Real> & cur_cons  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real> & cur_prim  = S_prim.array(mfi);
        const Array4<      Real> & cur_xmom  = S_data[IntVars::xmom].array(mfi);
        const Array4<      Real> & cur_ymom  = S_data[IntVars::ymom].array(mfi);
        const Array4<      Real> & cur_zmom  = S_data[IntVars::zmom].array(mfi);

        Array4<Real> avg_xmom_arr = avg_xmom.array(mfi);
        Array4<Real> avg_ymom_arr = avg_ymom.array(mfi);
        Array4<Real> avg_zmom_arr = avg_zmom.array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);

        const Array4<const Real>& z_nd         = z_phys_nd->const_array(mfi);
        const Array4<const Real>& z_cc         = z_phys_cc->const_array(mfi);
        const Array4<const Real>& detJ_new_arr = l_moving_terrain ? detJ_new->const_array(mfi)    : Array4<const Real>{};

        const Array4<const Real>& t_blank_arr = (terrain_blank) ? terrain_blank->const_array(mfi) :
                                                                Array4<const Real>{};
        const Array4<const Real>& t_blank_xface_arr = (terrain_blank_xface) ? terrain_blank_xface->const_array(mfi) :
                                                                              Array4<const Real>{};
        const Array4<const Real>& t_blank_yface_arr = (terrain_blank_yface) ? terrain_blank_yface->const_array(mfi) :
                                                                              Array4<const Real>{};
        const Array4<const Real>& t_blank_zface_arr = (terrain_blank_zface) ? terrain_blank_zface->const_array(mfi) :
                                                                              Array4<const Real>{};

        // Map factors
        const Array4<const Real>& mf_mx = mapfac[MapFacType::m_x]->const_array(mfi);
        const Array4<const Real>& mf_ux = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vx = mapfac[MapFacType::v_x]->const_array(mfi);
        const Array4<const Real>& mf_my = mapfac[MapFacType::m_y]->const_array(mfi);
        const Array4<const Real>& mf_uy = mapfac[MapFacType::u_y]->const_array(mfi);
        const Array4<const Real>& mf_vy = mapfac[MapFacType::v_y]->const_array(mfi);

        // SmnSmn for KE src with Deardorff or k-eqn RANS
        const Array4<const Real>& SmnSmn_a = l_need_SmnSmn ? SmnSmn->const_array(mfi) : Array4<const Real>{};

        // **************************************************************************
        // Here we fill the "current" data with "new" data because that is the result of the previous RK stage
        // **************************************************************************
        int nsv = S_old[IntVars::cons].nComp() - 2;
        const GpuArray<int, IntVars::NumTypes> scomp_slow = {  2,0,0,0};
        const GpuArray<int, IntVars::NumTypes> ncomp_slow = {nsv,0,0,0};

        // **************************************************************************
        // Note that here we do copy only the "slow" variables, not (rho) or (rho theta)
        // **************************************************************************
        ParallelFor(tbx, ncomp_slow[IntVars::cons],
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
            const int n = scomp_slow[IntVars::cons] + nn;
            cur_cons(i,j,k,n) = new_cons(i,j,k,n);
        });

        // Non-EB Anelastic: Per-tile copy of projected momentum (EB done above)
        if (l_anelastic && !l_use_eb) {
            Box tbx_inc = mfi.nodaltilebox(0);
            Box tby_inc = mfi.nodaltilebox(1);
            Box tbz_inc = mfi.nodaltilebox(2);

            ParallelFor(tbx_inc, tby_inc, tbz_inc,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                avg_xmom_arr(i,j,k) = cur_xmom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                avg_ymom_arr(i,j,k) = cur_ymom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                avg_zmom_arr(i,j,k) = cur_zmom(i,j,k);
            });
        }

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        bool l_eb_terrain_cc = false; // EB terrain on cell-centered grid
        Array4<const int> mask_arr{};
        Array4<const EBCellFlag> cfg_arr{};
        Array4<const Real> ax_arr{};
        Array4<const Real> ay_arr{};
        Array4<const Real> az_arr{};
        Array4<const Real> fcx_arr{};
        Array4<const Real> fcy_arr{};
        Array4<const Real> fcz_arr{};
        Array4<const Real> detJ_arr{};
        Array4<const Real> barea_arr{};
        Array4<const Real> bcent_arr{};

        if (l_use_eb) {
            EBCellFlagFab const& cfg = ebfact.getMultiEBCellFlagFab()[mfi];
            cfg_arr  = cfg.const_array();
            if (cfg.getType(tbx) == FabType::singlevalued) {
                l_eb_terrain_cc = true;
                ax_arr   = ebfact.getAreaFrac()[0]->const_array(mfi);
                ay_arr   = ebfact.getAreaFrac()[1]->const_array(mfi);
                az_arr   = ebfact.getAreaFrac()[2]->const_array(mfi);
                fcx_arr  = ebfact.getFaceCent()[0]->const_array(mfi);
                fcy_arr  = ebfact.getFaceCent()[1]->const_array(mfi);
                fcz_arr  = ebfact.getFaceCent()[2]->const_array(mfi);
                detJ_arr = ebfact.getVolFrac().const_array(mfi);
                mask_arr = physbnd_mask.const_array(mfi);
                barea_arr = ebfact.getBndryArea().const_array(mfi);
                bcent_arr = ebfact.getBndryCent().const_array(mfi);
            } else {
                ax_arr   = ax->const_array(mfi);
                ay_arr   = ay->const_array(mfi);
                az_arr   = az->const_array(mfi);
                detJ_arr = detJ->const_array(mfi);
            }
        } else {
            ax_arr   = ax->const_array(mfi);
            ay_arr   = ay->const_array(mfi);
            az_arr   = az->const_array(mfi);
            detJ_arr = detJ->const_array(mfi);
        }

        AdvType horiz_adv_type, vert_adv_type;
        Real    horiz_upw_frac, vert_upw_frac;

        Array4<Real> diffflux_x, diffflux_y, diffflux_z;
        Array4<Real> hfx_x, hfx_y, hfx_z, diss;
        Array4<Real> q1fx_x, q1fx_y, q1fx_z, q2fx_z;
        Array4<Real> hfx_EB{};

        if (l_use_diff) {
            diffflux_x = dflux_x->array(mfi);
            diffflux_y = dflux_y->array(mfi);
            diffflux_z = dflux_z->array(mfi);

            hfx_x = Hfx1->array(mfi);
            hfx_y = Hfx2->array(mfi);
            hfx_z = Hfx3->array(mfi);
            diss  = Diss->array(mfi);

            if (Q1fx1) q1fx_x = Q1fx1->array(mfi);
            if (Q1fx2) q1fx_y = Q1fx2->array(mfi);
            if (Q1fx3) q1fx_z = Q1fx3->array(mfi);
            if (Q2fx3) q2fx_z = Q2fx3->array(mfi);
        }

        if (l_use_diff && l_use_turb) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                eddyDiffs != nullptr,
                "erf_slow_rhs_post: active turbulence requires non-null eddyDiffs");
        }
        const Array4<const Real>& mu_turb =
            l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        //
        // Note that we either advect and diffuse all or none of the moisture variables
        //
        for (int ivar(RhoKE_comp); ivar<= RhoQ1_comp; ++ivar)
        {

            if (is_valid_slow_var[ivar])
            {
                start_comp = ivar;
                num_comp = 1;

                if (ivar == RhoQ1_comp) {
                    horiz_adv_type = ac.moistscal_horiz_adv_type;
                     vert_adv_type = ac.moistscal_vert_adv_type;
                    horiz_upw_frac = ac.moistscal_horiz_upw_frac;
                     vert_upw_frac = ac.moistscal_vert_upw_frac;

                    if (ac.use_efficient_advection){
                         horiz_adv_type = EfficientAdvType(nrk,ac.moistscal_horiz_adv_type);
                          vert_adv_type = EfficientAdvType(nrk,ac.moistscal_vert_adv_type);
                    }

                    num_comp = n_qstate;

                } else {
                    horiz_adv_type = ac.dryscal_horiz_adv_type;
                     vert_adv_type = ac.dryscal_vert_adv_type;
                    horiz_upw_frac = ac.dryscal_horiz_upw_frac;
                     vert_upw_frac = ac.dryscal_vert_upw_frac;

                    if (ac.use_efficient_advection){
                         horiz_adv_type = EfficientAdvType(nrk,ac.dryscal_horiz_adv_type);
                          vert_adv_type = EfficientAdvType(nrk,ac.dryscal_vert_adv_type);
                    }

                    if (ivar == RhoScalar_comp) {
                        num_comp = NSCALARS;
                    }
                }

                if (( ivar != RhoKE_comp                 ) ||
                    ((ivar == RhoKE_comp) && l_advect_KE))
                {
                    if (!l_eb_terrain_cc){
                        AdvectionSrcForScalars(tbx, start_comp, num_comp,
                                               avg_xmom_arr, avg_ymom_arr, avg_zmom_arr,
                                               cur_prim, cell_rhs,
                                               detJ_arr, dxInv, mf_mx, mf_my,
                                               horiz_adv_type, vert_adv_type,
                                               horiz_upw_frac, vert_upw_frac,
                                               flx_arr, domain, bc_ptr_h);
                    } else {
                        EBAdvectionSrcForScalars(tbx, start_comp, num_comp,
                                                 avg_xmom_arr, avg_ymom_arr, avg_zmom_arr,
                                                 cur_prim, cell_rhs,
                                                 mask_arr, cfg_arr, ax_arr, ay_arr, az_arr,
                                                 fcx_arr, fcy_arr, fcz_arr,
                                                 detJ_arr, dxInv, mf_mx, mf_my,
                                                 horiz_adv_type, vert_adv_type,
                                                 horiz_upw_frac, vert_upw_frac,
                                                 flx_arr, domain, bc_ptr_h,
                                                 already_on_centroids);
                    }
                }

                if (l_use_diff)
                {
                    // Allow for implicit moisture diffusion
                    Real l_vert_implicit_fac = zero;
                    if ( (ivar == RhoKE_comp && solverChoice.implicit_ke_diffusion      ) ||
                         (ivar == RhoQ1_comp && solverChoice.implicit_moisture_diffusion) ) {
                        l_vert_implicit_fac = solverChoice.vert_implicit_fac[level][nrk];
                    }

                    const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

                    // Only the physical chamber needs separate qv/qc calls:
                    // its wall correction must be applied to distinct flux
                    // components.  Generic moisture models retain the
                    // established multi-component diffusion call.
                    const bool componentwise_moisture =
                        use_physical_chamber_wall_flux && ivar == RhoQ1_comp;
                    const int n_diff_calls = componentwise_moisture ? n_qstate : 1;
                    for (int qstate = 0; qstate < n_diff_calls; ++qstate) {
                        const int state_comp = componentwise_moisture ?
                            RhoQ1_comp + qstate : start_comp;
                        const int diffusion_start = state_comp;
                        const int diffusion_num = componentwise_moisture ? 1 : num_comp;
                        const int flux_comp = componentwise_moisture ? qstate : 0;
                        AMREX_ALWAYS_ASSERT(state_comp >= 0 && state_comp < nvars);
                        AMREX_ALWAYS_ASSERT(flux_comp < dflux_x->nComp());
                        AMREX_ALWAYS_ASSERT(flux_comp < dflux_y->nComp());
                        AMREX_ALWAYS_ASSERT(flux_comp < dflux_z->nComp());
                        const Array4<Real> diffusion_x = dflux_x->array(mfi, flux_comp);
                        const Array4<Real> diffusion_y = dflux_y->array(mfi, flux_comp);
                        const Array4<Real> diffusion_z = dflux_z->array(mfi, flux_comp);

                    if (solverChoice.mesh_type == MeshType::StretchedDz) {
                        DiffusionSrcForState_S(tbx, domain, diffusion_start, diffusion_num, u, v,
                                               new_cons, cur_prim, cell_rhs,
                                               diffusion_x, diffusion_y, diffusion_z,
                                               stretched_dz_d, dxInv, SmnSmn_a,
                                               mf_mx, mf_ux, mf_vx,
                                               mf_my, mf_uy, mf_vy,
                                               hfx_z, q1fx_z, q2fx_z, diss,
                                               mu_turb, solverChoice, level,
                                               tm_arr, grav_gpu, bc_ptr_d, l_apply_surface_layer_fluxes_in_diffusion, l_vert_implicit_fac);
                    } else if (l_use_terrain) {
                        DiffusionSrcForState_T(tbx, domain, diffusion_start, diffusion_num, l_rotate, u, v,
                                               new_cons, cur_prim, cell_rhs,
                                               diffusion_x, diffusion_y, diffusion_z,
                                               z_nd, z_cc, ax_arr, ay_arr, az_arr,
                                               detJ_arr, dxInv, SmnSmn_a,
                                               mf_mx, mf_ux, mf_vx,
                                               mf_my, mf_uy, mf_vy,
                                               hfx_x, hfx_y, hfx_z, q1fx_x, q1fx_y, q1fx_z,q2fx_z, diss,
                                               mu_turb, solverChoice, level,
                                               tm_arr, grav_gpu, bc_ptr_d, l_apply_surface_layer_fluxes_in_diffusion, l_vert_implicit_fac);
                    } else if (l_use_eb) {
                        DiffusionSrcForState_EB(tbx, domain, diffusion_start, diffusion_num, u, v,
                                                new_cons, cur_prim, cell_rhs,
                                                diffusion_x, diffusion_y, diffusion_z,
                                                cfg_arr, ax_arr, ay_arr, az_arr, detJ_arr,
                                                barea_arr, bcent_arr,
                                                dx, dxInv,
                                                hfx_z, q1fx_z, q2fx_z, hfx_EB,
                                                mu_turb, solverChoice, level,
                                                bc_ptr_d, l_apply_surface_layer_fluxes_in_diffusion);
                    } else {
                        DiffusionSrcForState_N(tbx, domain, diffusion_start, diffusion_num, u, v,
                                               new_cons, cur_prim, cell_rhs,
                                               diffusion_x, diffusion_y, diffusion_z, dxInv, SmnSmn_a,
                                               mf_mx, mf_ux, mf_vx,
                                               mf_my, mf_uy, mf_vy,
                                               hfx_z, q1fx_z, q2fx_z, diss,
                                               mu_turb, solverChoice, level,
                                               tm_arr, grav_gpu, bc_ptr_d, l_apply_surface_layer_fluxes_in_diffusion, l_vert_implicit_fac);
                    }
                    if (use_physical_chamber_wall_flux) {
                        // Apply the physical wall correction immediately to
                        // the flux component just computed.  This keeps the
                        // q-state diffusion path identical with budgets on
                        // and off and guarantees no stale flux is consumed.
                        // The diffusion views are component-shifted; the
                        // wall helper receives the unshifted views and the
                        // explicit flux component index.
                        erf_resolved_wall_flux::apply(
                            tbx, domain, state_comp, flux_comp, new_cons, cur_prim,
                            cloud_chamber_base_state->const_array(mfi), cell_rhs,
                            diffflux_x, diffflux_y, diffflux_z, dxInv,
                            chamber_walls, dc.alpha_T, dc.alpha_C,
                            solverChoice.rdOcp);
                    }
                    }
                } // use_diff


            } // valid slow var
        } // loop ivar

#ifdef ERF_USE_EAMXX_SHOC
        if (tc.uses_eamxx_shoc() && eamxx_shoc_lev) {
            eamxx_shoc_lev->add_slow_tend(mfi,tbx,cell_rhs);
        }
#endif
        // This updates just the "slow" conserved variables
        {
        BL_PROFILE("rhs_post_8");

        const Real eps = std::numeric_limits<Real>::epsilon();

        auto const& src_arr = source.const_array(mfi);

        for (int ivar(RhoKE_comp); ivar<= RhoQ1_comp; ++ivar)
        {
            if (is_valid_slow_var[ivar])
            {
                start_comp = ivar;
                num_comp = 1;
                if (ivar == RhoQ1_comp) {
                    num_comp = nvars - RhoQ1_comp;
                } else if (ivar == RhoScalar_comp) {
                    num_comp = NSCALARS;
                }

               if (l_moving_terrain)
               {
                    ParallelFor(tbx, num_comp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = start_comp + nn;
                        cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                        Real temp_val = detJ_arr(i,j,k) * old_cons(i,j,k,n) + dt * detJ_arr(i,j,k) * cell_rhs(i,j,k,n);
                        cur_cons(i,j,k,n) = temp_val / detJ_new_arr(i,j,k);
                        if (ivar == RhoKE_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), eps);
                        }
                    });

                } else if (l_anelastic && (nrk == 1)) { // not moving and ( (anelastic) and second RK stage) )

                    ParallelFor(tbx, num_comp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = start_comp + nn;
                        cell_rhs(i,j,k,n) += src_arr(i,j,k,n);

                        // Re-construct the cell_rhs used in the first RK stage
                        Real dt_times_old_cell_rhs = cur_cons(i,j,k,n) - old_cons(i,j,k,n);

                        // Add the time-averaged RHS to the old state
                        cur_cons(i,j,k,n) = old_cons(i,j,k,n) + myhalf * (dt_times_old_cell_rhs + dt * cell_rhs(i,j,k,n));

                        if (ivar == RhoKE_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), eps);
                        } else if (ivar >= RhoQ1_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), amrex::Real(0));
                        }
                    });

                } else { // not moving and ( (not anelastic) or (first RK stage) )

                    ParallelFor(tbx, num_comp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = start_comp + nn;
                        cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                        cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
                        if (ivar == RhoKE_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), eps);
                        } else if (ivar >= RhoQ1_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), amrex::Real(0));
                        }
                    });

                } // moving, anelastic or neither?

            } // is_valid
        } // ivar
        } // profile

        {
        BL_PROFILE("rhs_post_9");
        // This updates all the conserved variables (not just the "slow" ones)
        int   num_comp_all = S_data[IntVars::cons].nComp();
        ParallelFor(tbx, num_comp_all,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            new_cons(i,j,k,n)  = cur_cons(i,j,k,n);
        });
        } // end profile

        Box xtbx = mfi.nodaltilebox(0);
        Box ytbx = mfi.nodaltilebox(1);
        Box ztbx = mfi.nodaltilebox(2);

        {
        BL_PROFILE("rhs_post_10()");
        if (l_anelastic && terrain_blank) { // explicitly set fully immersed cells to have 0 velocities for anelastic (unstable for fully compressible).
            ParallelFor(xtbx, ytbx, ztbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // Use face-centered terrain_blanking if available, otherwise average from cell centers
                Real t_blank = (t_blank_xface_arr) ? t_blank_xface_arr(i, j, k) :
                               myhalf * (t_blank_arr(i, j, k) + t_blank_arr(i-1, j, k));
                if (t_blank == one) {
                    new_xmom(i,j,k) = zero;
                } else {
                    new_xmom(i,j,k) = cur_xmom(i,j,k);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // Use face-centered terrain_blanking if available, otherwise average from cell centers
                Real t_blank = (t_blank_yface_arr) ? t_blank_yface_arr(i, j, k) :
                               myhalf * (t_blank_arr(i, j, k) + t_blank_arr(i, j-1, k));
                if (t_blank == one) {
                    new_ymom(i,j,k) = zero;
                } else {
                    new_ymom(i,j,k) = cur_ymom(i,j,k);
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // Use face-centered terrain_blanking if available, otherwise average from cell centers
                Real t_blank = (t_blank_zface_arr) ? t_blank_zface_arr(i, j, k) :
                               myhalf * (t_blank_arr(i, j, k) + t_blank_arr(i, j, k-1));
                if (t_blank == one) {
                    new_zmom(i,j,k) = zero;
                } else {
                    new_zmom(i,j,k) = cur_zmom(i,j,k);
                }
            });
            ParallelFor(tbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real t_blank       = t_blank_arr(i, j, k  );
                if (t_blank == one) { // don't update rho and theta for fully immersed cells
                    new_cons(i, j, k, Rho_comp)      = old_cons(i, j, k, Rho_comp);
                    new_cons(i, j, k, RhoTheta_comp) = old_cons(i, j, k, RhoTheta_comp);
                }
            });
        } else {
            ParallelFor(xtbx, ytbx, ztbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                new_xmom(i,j,k) = cur_xmom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                new_ymom(i,j,k) = cur_ymom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                new_zmom(i,j,k) = cur_zmom(i,j,k);
            });
        }
        } // end profile

        {
        BL_PROFILE("rhs_post_10");
        // We only add to the flux registers in the final RK step
        if (l_reflux) {
            int strt_comp_reflux = RhoTheta_comp + 1;
            int  num_comp_reflux = nvars - strt_comp_reflux;
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

            // This is necessary here so we don't go on to the next FArrayBox without
            // having finished copying the fluxes into the FluxRegisters (since the fluxes
            // are stored in temporary FArrayBox's)
            Gpu::streamSynchronize();

        } // two-way coupling
        } // end profile
      } // mfi
    } // OMP
    if (cloud_budget && l_use_diff && n_qstate > 0) {
        for (int qstate = 0; qstate < n_qstate; ++qstate) {
            MultiFab qflux_x(*dflux_x, make_alias, qstate, 1);
            MultiFab qflux_y(*dflux_y, make_alias, qstate, 1);
            MultiFab qflux_z(*dflux_z, make_alias, qstate, 1);
            cloud_budget->capture_stage(
                qstate == 0 ? CloudChamberBudget::RhoQv : CloudChamberBudget::RhoQc,
                nrk, static_cast<Real>(dt_d), qflux_x, qflux_y, qflux_z, geom);
        }
    }
}
