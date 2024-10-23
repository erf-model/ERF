#include <AMReX.H>
#include <ERF_Src_headers.H>
#include <ERF_TI_slow_headers.H>

#if defined(ERF_USE_NETCDF)
// #include <ERF_Utils.H>
#endif

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the scalars other than density or potential temperature
 *
 * @param[in]  level level of resolution
 * @param[in]  finest_level finest level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[out]  S_rhs RHS computed here
 * @param[in]  S_old solution at start of time step
 * @param[in]  S_new solution at end of current RK stage
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in]  S_scratch scratch space
 * @param[in]  xvel x-component of velocity
 * @param[in]  yvel y-component of velocity
 * @param[in]  zvel z-component of velocity
 * @param[in] source source terms for conserved variables
 * @param[in] SmnSmn strain rate magnitude
 * @param[in] eddyDiffs diffusion coefficients for LES turbulence models
 * @param[in] Hfx3 heat flux in z-dir
 * @param[in] Diss dissipation of turbulent kinetic energy
 * @param[in]  geom   Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in]  most  Pointer to MOST class for Monin-Obukhov Similarity Theory boundary condition
 * @param[in]  domain_bcs_type_d device vector for domain boundary conditions
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] ax area fractions on x-faces
 * @param[in] ay area fractions on y-faces
 * @param[in] az area fractions on z-faces
 * @param[in] detJ     Jacobian of the metric transformation at start of time step (= 1 if use_terrain is false)
 * @param[in] detJ_new Jacobian of the metric transformation at new RK stage time (= 1 if use_terrain is false)
 * @param[in] mapfac_m map factor at cell centers
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 */

void erf_slow_rhs_post (int level, int finest_level,
                        int nrk,
                        Real dt,
                        int n_qstate,
                        Vector<MultiFab>& S_rhs,
                        Vector<MultiFab>& S_old,
                        Vector<MultiFab>& S_new,
                        Vector<MultiFab>& S_data,
                        const MultiFab& S_prim,
                        Vector<MultiFab>& S_scratch,
                        const MultiFab& xvel,
                        const MultiFab& yvel,
                        const MultiFab& /*zvel*/,
                        const MultiFab& source,
                        const MultiFab* SmnSmn,
                        const MultiFab* eddyDiffs,
                        MultiFab* Hfx1,
                        MultiFab* Hfx2,
                        MultiFab* Hfx3,
                        MultiFab* Q1fx1,
                        MultiFab* Q1fx2,
                        MultiFab* Q1fx3,
                        MultiFab* Q2fx3,
                        MultiFab* Diss,
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
                        std::unique_ptr<MultiFab>& detJ_new,
                        std::unique_ptr<MultiFab>& mapfac_m,
                        std::unique_ptr<MultiFab>& mapfac_u,
                        std::unique_ptr<MultiFab>& mapfac_v,
#ifdef ERF_USE_EB
                        amrex::EBFArrayBoxFactory const& ebfact,
#endif
#if defined(ERF_USE_NETCDF)
                        const bool& moist_set_rhs_bool,
                        const Real& bdy_time_interval,
                        const Real& start_bdy_time,
                        const Real& new_stage_time,
                        int  width,
                        int  set_width,
                        Vector<Vector<FArrayBox>>& bdy_data_xlo,
                        Vector<Vector<FArrayBox>>& bdy_data_xhi,
                        Vector<Vector<FArrayBox>>& bdy_data_ylo,
                        Vector<Vector<FArrayBox>>& bdy_data_yhi,
#endif
                        YAFluxRegister* fr_as_crse,
                        YAFluxRegister* fr_as_fine)
{
    BL_PROFILE_REGION("erf_slow_rhs_post()");

    const BCRec* bc_ptr_d = domain_bcs_type_d.data();
    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    AdvChoice  ac = solverChoice.advChoice;
    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const MultiFab* t_mean_mf = nullptr;
    if (most) t_mean_mf = most->get_mac_avg(0,2);

    const bool l_use_terrain      = solverChoice.use_terrain;
    const bool l_reflux = (solverChoice.coupling_type != CouplingType::OneWay);
    const bool l_moving_terrain   = (solverChoice.terrain_type == TerrainType::Moving);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT(l_use_terrain);

    const bool l_use_mono_adv   = solverChoice.use_mono_adv;
    const bool l_use_KE         = ( (tc.les_type == LESType::Deardorff) ||
                                    (tc.pbl_type == PBLType::MYNN) );
    const bool l_use_diff       = ((dc.molec_diff_type != MolecDiffType::None) ||
                                   (tc.les_type        !=       LESType::None) ||
                                   (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_turb       = ( tc.les_type == LESType::Smagorinsky ||
                                    tc.les_type == LESType::Deardorff   ||
                                    tc.pbl_type == PBLType::MYNN25      ||
                                    tc.pbl_type == PBLType::YSU );
    const bool exp_most         = (solverChoice.use_explicit_most);
    const bool rot_most         = (solverChoice.use_rotate_most);

    const Box& domain = geom.Domain();

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const Real* dx = geom.CellSize();

    // *************************************************************************
    // Set gravity as a vector
    // *************************************************************************
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *************************************************************************
    // Pre-computed quantities
    // *************************************************************************
    int nvars                     = S_data[IntVars::cons].nComp();
    const BoxArray& ba            = S_data[IntVars::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVars::cons].DistributionMap();

    std::unique_ptr<MultiFab> dflux_x;
    std::unique_ptr<MultiFab> dflux_y;
    std::unique_ptr<MultiFab> dflux_z;

    if (l_use_diff) {
        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, nvars, 0);
    } else {
        dflux_x = nullptr;
        dflux_y = nullptr;
        dflux_z = nullptr;
    }

    // Valid vars
    Vector<int> is_valid_slow_var; is_valid_slow_var.resize(RhoQ1_comp+1,0);
    if (l_use_KE) {is_valid_slow_var[    RhoKE_comp] = 1;}
                   is_valid_slow_var[RhoScalar_comp] = 1;
    if (solverChoice.moisture_type != MoistureType::None) {
         is_valid_slow_var[RhoQ1_comp] = 1;
    }

    // *****************************************************************************
    // Monotonic advection for scalars
    // *****************************************************************************
    int nvar = S_new[IntVars::cons].nComp();
    Vector<Real> max_scal(nvar, 1.0e34); Gpu::DeviceVector<Real> max_scal_d(nvar);
    Vector<Real> min_scal(nvar,-1.0e34); Gpu::DeviceVector<Real> min_scal_d(nvar);
    if (l_use_mono_adv) {
        auto const& ma_s_arr = S_new[IntVars::cons].const_arrays();
        for (int ivar(RhoKE_comp); ivar<nvar; ++ivar) {
            GpuTuple<Real,Real> mm = ParReduce(TypeList<ReduceOpMax,ReduceOpMin>{},
                                               TypeList<Real, Real>{},
                                               S_new[IntVars::cons], IntVect(0),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                -> GpuTuple<Real,Real>
                {
                    return { ma_s_arr[box_no](i,j,k,ivar), ma_s_arr[box_no](i,j,k,ivar) };
                });
            max_scal[ivar] = get<0>(mm);
            min_scal[ivar] = get<1>(mm);
        }
    }
    Gpu::copy(Gpu::hostToDevice, max_scal.begin(), max_scal.end(), max_scal_d.begin());
    Gpu::copy(Gpu::hostToDevice, min_scal.begin(), min_scal.end(), min_scal_d.begin());
    Real* max_s_ptr = max_scal_d.data();
    Real* min_s_ptr = min_scal_d.data();

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

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      std::array<FArrayBox,AMREX_SPACEDIM> flux;
      std::array<FArrayBox,AMREX_SPACEDIM> flux_tmp;

      int start_comp;
      int   num_comp;

      for ( MFIter mfi(S_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box tbx  = mfi.tilebox();

        // *************************************************************************
        // Define flux arrays for use in advection
        // *************************************************************************
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flux[dir].resize(surroundingNodes(tbx,dir),nvars);
            flux[dir].setVal<RunOn::Device>(0.);
            if (l_use_mono_adv) {
                flux_tmp[dir].resize(surroundingNodes(tbx,dir),nvars);
                flux_tmp[dir].setVal<RunOn::Device>(0.);
            }
        }
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>
            flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};
        Array4<Real> tmpx = (l_use_mono_adv) ? flux_tmp[0].array() : Array4<Real>{};
        Array4<Real> tmpy = (l_use_mono_adv) ? flux_tmp[1].array() : Array4<Real>{};
        Array4<Real> tmpz = (l_use_mono_adv) ? flux_tmp[2].array() : Array4<Real>{};
        const GpuArray<Array4<Real>, AMREX_SPACEDIM> flx_tmp_arr{{AMREX_D_DECL(tmpx,tmpy,tmpz)}};

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

        Array4<Real> avg_xmom = S_scratch[IntVars::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVars::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVars::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        const Array4<const Real>& z_nd         = l_use_terrain    ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ_new_arr = l_moving_terrain ? detJ_new->const_array(mfi)    : Array4<const Real>{};

        // Map factors
        const Array4<const Real>& mf_m = mapfac_m->const_array(mfi);
        const Array4<const Real>& mf_u = mapfac_u->const_array(mfi);
        const Array4<const Real>& mf_v = mapfac_v->const_array(mfi);

        // SmnSmn for KE src with Deardorff
        const Array4<const Real>& SmnSmn_a = l_use_KE ? SmnSmn->const_array(mfi) : Array4<const Real>{};

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

        // We have projected the velocities stored in S_data but we will use
        //    the velocities stored in S_scratch to update the scalars, so
        //    we need to copy from S_data (projected) into S_scratch
        if (solverChoice.anelastic[level]) {
            Box tbx_inc = mfi.nodaltilebox(0);
            Box tby_inc = mfi.nodaltilebox(1);
            Box tbz_inc = mfi.nodaltilebox(2);

            ParallelFor(tbx_inc, tby_inc, tbz_inc,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                avg_xmom(i,j,k) = cur_xmom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                avg_ymom(i,j,k) = cur_ymom(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                avg_zmom(i,j,k) = cur_zmom(i,j,k);
            });
        }

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        // Metric terms
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

        AdvType horiz_adv_type, vert_adv_type;
        Real    horiz_upw_frac, vert_upw_frac;

        Array4<Real> diffflux_x, diffflux_y, diffflux_z;
        Array4<Real> hfx_x, hfx_y, hfx_z, diss;
        Array4<Real> q1fx_x, q1fx_y, q1fx_z, q2fx_z;
        const bool use_most = (most != nullptr);

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

        //
        // Note that we either advect and diffuse all or none of the moisture variables
        //
        for (int ivar(RhoKE_comp); ivar<= RhoQ1_comp; ++ivar)
        {
            if (is_valid_slow_var[ivar])
            {
                start_comp = ivar;

                if (ivar >= RhoQ1_comp) {
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
                    num_comp = 1;
                }

                if (l_use_diff) {

                    const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

                    if (l_use_terrain) {
                        DiffusionSrcForState_T(tbx, domain, start_comp, num_comp, exp_most, rot_most, u, v,
                                               new_cons, cur_prim, cell_rhs,
                                               diffflux_x, diffflux_y, diffflux_z,
                                               z_nd, ax_arr, ay_arr, az_arr, detJ_arr,
                                               dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                               hfx_x, hfx_y, hfx_z, q1fx_x, q1fx_y, q1fx_z,q2fx_z, diss,
                                               mu_turb, solverChoice, level,
                                               tm_arr, grav_gpu, bc_ptr_d, use_most);
                    } else {
                        DiffusionSrcForState_N(tbx, domain, start_comp, num_comp, exp_most, u, v,
                                               new_cons, cur_prim, cell_rhs,
                                               diffflux_x, diffflux_y, diffflux_z,
                                               dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                               hfx_z, q1fx_z, q2fx_z, diss,
                                               mu_turb, solverChoice, level,
                                               tm_arr, grav_gpu, bc_ptr_d, use_most);
                    }
                } // use_diff
            } // valid slow var
        } // loop ivar

#if defined(ERF_USE_NETCDF)
        if (moist_set_rhs_bool)
        {
            Box gtbx_moist  = mfi.tilebox(IntVect(0),IntVect(2,2,0));
            const Array4<const Real> & old_cons_const = S_old[IntVars::cons].const_array(mfi);
            const Array4<const Real> & new_cons_const = S_new[IntVars::cons].const_array(mfi);
            moist_set_rhs(tbx, gtbx_moist, old_cons_const, new_cons_const, cell_rhs,
                          bdy_time_interval, start_bdy_time, new_stage_time, dt, width, set_width, domain,
                          bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi);
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

                if (ivar >= RhoQ1_comp) {
                    num_comp = nvars - RhoQ1_comp;
                } else {
                    num_comp = 1;
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

                } else {

                    ParallelFor(tbx, num_comp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                        const int n = start_comp + nn;
                        cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                        cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
                        if (ivar == RhoKE_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), eps);
                        } else if (ivar >= RhoQ1_comp) {
                            cur_cons(i,j,k,n) = amrex::max(cur_cons(i,j,k,n), 0.0);
                        }
                    });

                } // moving or not?

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
        } // end profile

        {
        BL_PROFILE("rhs_post_10");
        // We only add to the flux registers in the final RK step
        if (l_reflux && nrk == 2) {
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
}
