#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <Advection.H>
#include <Diffusion.H>
#include <NumericalDiffusion.H>
#include <TI_headers.H>
#include <TileNoZ.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the scalars other than density or potential temperature
 *
 * @param[in]  level level of resolution
 * @param[in]  dt    slow time step
 * @param[in]  grids_to_evolve the region in the domain excluding the relaxation and specified zones
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
 * @param[in]  geom   Container for geometric informaiton
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in]  most  Pointer to MOST class for Monin-Obukhov Similarity Theory boundary condition
 * @param[in]  domain_bcs_type_d device vector for domain boundary conditions
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] detJ     Jacobian of the metric transformation at start of time step (= 1 if use_terrain is false)
 * @param[in] detJ_new Jacobian of the metric transformation at new RK stage time (= 1 if use_terrain is false)
 * @param[in] mapfac_m map factor at cell centers
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 */

void erf_slow_rhs_post (int /*level*/, Real dt,
                        BoxArray& grids_to_evolve,
                        Vector<MultiFab>& S_rhs,
                        Vector<MultiFab>& S_old,
                        Vector<MultiFab>& S_new,
                        Vector<MultiFab>& S_data,
                        const MultiFab& S_prim,
                        Vector<MultiFab>& S_scratch,
                        const MultiFab& xvel,
                        const MultiFab& yvel,
                        const MultiFab& zvel,
                        const MultiFab& source,
                        const MultiFab* SmnSmn,
                        const MultiFab* eddyDiffs,
                        MultiFab* Hfx3, MultiFab* Diss,
                        const amrex::Geometry geom,
                        const SolverChoice& solverChoice,
                        std::unique_ptr<ABLMost>& most,
                        const Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                        std::unique_ptr<MultiFab>& z_phys_nd,
                        std::unique_ptr<MultiFab>& detJ,
                        std::unique_ptr<MultiFab>& detJ_new,
                        std::unique_ptr<MultiFab>& mapfac_m,
                        std::unique_ptr<MultiFab>& mapfac_u,
                        std::unique_ptr<MultiFab>& mapfac_v)
{
    BL_PROFILE_REGION("erf_slow_rhs_post()");

    const MultiFab* t_mean_mf = nullptr;
    if (most) t_mean_mf = most->get_mac_avg(0,2);

    const int  l_horiz_spatial_order = solverChoice.horiz_spatial_order;
    const int  l_vert_spatial_order  = solverChoice.vert_spatial_order;
    const bool l_use_terrain    = solverChoice.use_terrain;
    const bool l_moving_terrain = (solverChoice.terrain_type == 1);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT(l_use_terrain);

    const bool l_use_ndiff      = solverChoice.use_NumDiff;
    const bool l_use_QKE        = solverChoice.use_QKE && solverChoice.advect_QKE;
    const bool l_use_deardorff  = (solverChoice.les_type == LESType::Deardorff);
    const bool l_use_diff       = ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
                                    (solverChoice.les_type        !=       LESType::None) ||
                                    (solverChoice.pbl_type        !=       PBLType::None) );
    const bool l_use_turb       = ( solverChoice.les_type == LESType::Smagorinsky ||
                                    solverChoice.les_type == LESType::Deardorff   ||
                                    solverChoice.pbl_type == PBLType::MYNN25 );
    const bool l_all_WENO       = solverChoice.all_use_WENO;
    const bool l_moist_WENO     = solverChoice.moist_use_WENO;
    const bool l_all_WENO_Z     = solverChoice.all_use_WENO_Z;
    const bool l_moist_WENO_Z   = solverChoice.moist_use_WENO_Z;
    const int  l_spatial_order_WENO = solverChoice.spatial_order_WENO;

    const amrex::BCRec* bc_ptr = domain_bcs_type_d.data();

    const Box& domain = geom.Domain();

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *************************************************************************
    // Pre-computed quantities
    // *************************************************************************
    int nvars                     = S_data[IntVar::cons].nComp();
    const BoxArray& ba            = S_data[IntVar::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVar::cons].DistributionMap();

    MultiFab* dflux_x = nullptr;
    MultiFab* dflux_y = nullptr;
    MultiFab* dflux_z = nullptr;

    if (l_use_diff) {
        dflux_x = new MultiFab(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = new MultiFab(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = new MultiFab(convert(ba,IntVect(0,0,1)), dm, nvars, 0);
    }

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

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& valid_bx = grids_to_evolve[mfi.index()];

        const Box& tbx = mfi.tilebox();

        // Construct intersection of current tilebox and valid region for updating
        Box bx = tbx & valid_bx;

        const Array4<      Real> & old_cons   = S_old[IntVar::cons].array(mfi);
        const Array4<      Real> & cell_rhs   = S_rhs[IntVar::cons].array(mfi);

        const Array4<      Real> & new_cons   = S_new[IntVar::cons].array(mfi);
        const Array4<      Real> & new_xmom  = S_new[IntVar::xmom].array(mfi);
        const Array4<      Real> & new_ymom  = S_new[IntVar::ymom].array(mfi);
        const Array4<      Real> & new_zmom  = S_new[IntVar::zmom].array(mfi);

        const Array4<      Real> & cur_cons  = S_data[IntVar::cons].array(mfi);
        const Array4<const Real> & cur_prim  = S_prim.array(mfi);
        const Array4<      Real> & cur_xmom  = S_data[IntVar::xmom].array(mfi);
        const Array4<      Real> & cur_ymom  = S_data[IntVar::ymom].array(mfi);
        const Array4<      Real> & cur_zmom  = S_data[IntVar::zmom].array(mfi);

        Array4<Real> avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        // Metric terms
        const Array4<const Real>& z_nd         = l_use_terrain    ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ_arr     = l_use_terrain    ? detJ->const_array(mfi)        : Array4<const Real>{};
        const Array4<const Real>& detJ_new_arr = l_moving_terrain ? detJ_new->const_array(mfi)    : Array4<const Real>{};

        // Map factors
        const Array4<const Real>& mf_m = mapfac_m->const_array(mfi);
        const Array4<const Real>& mf_u = mapfac_u->const_array(mfi);
        const Array4<const Real>& mf_v = mapfac_v->const_array(mfi);

        // SmnSmn for KE src with Deardorff
        const Array4<const Real>& SmnSmn_a = l_use_deardorff ? SmnSmn->const_array(mfi) : Array4<const Real>{};

        // **************************************************************************
        // Here we fill the "current" data with "new" data because that is the result of the previous RK stage
        // **************************************************************************
        int nsv = Cons::NumVars-2;
        const amrex::GpuArray<int, IntVar::NumVars> scomp_slow = {  2,0,0,0};
        const amrex::GpuArray<int, IntVar::NumVars> ncomp_slow = {nsv,0,0,0};

        {
        BL_PROFILE("rhs_post_7");
        ParallelFor(bx, ncomp_slow[IntVar::cons],
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
            const int n = scomp_slow[IntVar::cons] + nn;
            cur_cons(i,j,k,n) = new_cons(i,j,k,n);
        });
        } // end profile

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        int start_comp;
        int   num_comp;
        if (l_use_deardorff) {
            start_comp = RhoKE_comp;
              num_comp = 1;
            AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                                   cur_prim, cell_rhs, detJ_arr,
                                   dxInv, mf_m, l_all_WENO, l_moist_WENO, l_all_WENO_Z, l_moist_WENO_Z,
                                   l_spatial_order_WENO,l_horiz_spatial_order, l_vert_spatial_order,
                                   l_use_terrain);
        }
        if (l_use_QKE) {
            start_comp = RhoQKE_comp;
              num_comp = 1;
            AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                                   cur_prim, cell_rhs, detJ_arr,
                                   dxInv, mf_m, l_all_WENO, l_moist_WENO, l_all_WENO_Z, l_moist_WENO_Z,
                                   l_spatial_order_WENO,l_horiz_spatial_order, l_vert_spatial_order,
                                   l_use_terrain);
        }
        start_comp = RhoScalar_comp;
          num_comp = S_data[IntVar::cons].nComp() - start_comp;
        AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                               cur_prim, cell_rhs, detJ_arr,
                               dxInv, mf_m, l_all_WENO, l_moist_WENO, l_all_WENO_Z, l_moist_WENO_Z,
                               l_spatial_order_WENO,l_horiz_spatial_order, l_vert_spatial_order,
                               l_use_terrain);

        if (l_use_diff) {
            Array4<Real> diffflux_x = dflux_x->array(mfi);
            Array4<Real> diffflux_y = dflux_y->array(mfi);
            Array4<Real> diffflux_z = dflux_z->array(mfi);

            Array4<Real> hfx_z = Hfx3->array(mfi);
            Array4<Real> diss  = Diss->array(mfi);

            const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

            if (l_use_deardorff) {
                start_comp = RhoKE_comp;
                  num_comp = 1;
                if (l_use_terrain) {
                    DiffusionSrcForState_T(bx, domain, start_comp, num_comp, u, v,
                                           cur_cons, cur_prim, cell_rhs,
                                           diffflux_x, diffflux_y, diffflux_z, z_nd, detJ_arr,
                                           dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                           hfx_z, diss,
                                           mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
                } else {
                    DiffusionSrcForState_N(bx, domain, start_comp, num_comp, u, v,
                                           cur_cons, cur_prim, cell_rhs,
                                           diffflux_x, diffflux_y, diffflux_z,
                                           dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                           hfx_z, diss,
                                           mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
                }
                if (l_use_ndiff) {
                    NumericalDiffusion(bx, start_comp, num_comp, dt, solverChoice,
                                       new_cons, cell_rhs, mf_u, mf_v, false, false);
                }
            }
            if (l_use_QKE) {
                start_comp = RhoQKE_comp;
                  num_comp = 1;
                if (l_use_terrain) {
                    DiffusionSrcForState_T(bx, domain, start_comp, num_comp, u, v,
                                           cur_cons, cur_prim, cell_rhs,
                                           diffflux_x, diffflux_y, diffflux_z, z_nd, detJ_arr,
                                           dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                           hfx_z, diss,
                                           mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
                } else {
                    DiffusionSrcForState_N(bx, domain, start_comp, num_comp, u, v,
                                           cur_cons, cur_prim, cell_rhs,
                                           diffflux_x, diffflux_y, diffflux_z,
                                           dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                           hfx_z, diss,
                                           mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
                }
                if (l_use_ndiff) {
                    NumericalDiffusion(bx, start_comp, num_comp, dt, solverChoice,
                                       new_cons, cell_rhs, mf_u, mf_v, false, false);
                }
            }
            start_comp = RhoScalar_comp;
              num_comp = S_data[IntVar::cons].nComp() - start_comp;
            if (l_use_terrain) {
                DiffusionSrcForState_T(bx, domain, start_comp, num_comp, u, v,
                                       cur_cons, cur_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z, z_nd, detJ_arr,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       hfx_z, diss,
                                       mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
            } else {
                DiffusionSrcForState_N(bx, domain, start_comp, num_comp, u, v,
                                       cur_cons, cur_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       hfx_z, diss,
                                       mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
            }
            if (l_use_ndiff) {
                NumericalDiffusion(bx, start_comp, num_comp, dt, solverChoice,
                                   new_cons, cell_rhs, mf_u, mf_v, false, false);
            }
        }

        // NOTE: Computing the RHS is done over bx (union w/ grids to evolve).
        //       However, the update is over tbx (no union). The interior ghost
        //       cells have their RHS populated already.

        // This updates just the "slow" conserved variables
        {
        BL_PROFILE("rhs_post_8");

        if (l_moving_terrain)
        {
            auto const& src_arr = source.const_array(mfi);
            num_comp = S_data[IntVar::cons].nComp() - start_comp;
            ParallelFor(tbx, num_comp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                // NOTE: we don't include additional source terms when terrain is moving
                Real temp_val = detJ_arr(i,j,k) * old_cons(i,j,k,n) + dt * detJ_arr(i,j,k) * cell_rhs(i,j,k,n);
                cur_cons(i,j,k,n) = temp_val / detJ_new_arr(i,j,k);
            });

            if (l_use_deardorff) {
              start_comp = RhoKE_comp;
              num_comp   = 1;
              ParallelFor(tbx, num_comp,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                // NOTE: we don't include additional source terms when terrain is moving
                Real temp_val = detJ_arr(i,j,k) * old_cons(i,j,k,n) + dt * detJ_arr(i,j,k) * cell_rhs(i,j,k,n);
                cur_cons(i,j,k,n) = temp_val / detJ_new_arr(i,j,k);
              });
            }
            if (l_use_QKE) {
              start_comp = RhoQKE_comp;
              num_comp   = 1;
              ParallelFor(tbx, num_comp,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                // NOTE: we don't include additional source terms when terrain is moving
                Real temp_val = detJ_arr(i,j,k) * old_cons(i,j,k,n) + dt * detJ_arr(i,j,k) * cell_rhs(i,j,k,n);
                cur_cons(i,j,k,n) = temp_val / detJ_new_arr(i,j,k);
              });
            }
        } else {
            auto const& src_arr = source.const_array(mfi);
            num_comp = S_data[IntVar::cons].nComp() - start_comp;
            ParallelFor(tbx, num_comp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);

            });

            if (l_use_deardorff) {
              start_comp = RhoKE_comp;
              num_comp = 1;
              amrex::Real eps = std::numeric_limits<Real>::epsilon();
              ParallelFor(tbx, num_comp,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
                // make sure rho*e is positive
                if (cur_cons(i,j,k,n) < eps) cur_cons(i,j,k,n) = eps;
              });
            }
            if (l_use_QKE) {
              start_comp = RhoQKE_comp;
              num_comp   = 1;
              ParallelFor(tbx, num_comp,
              [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
              });
            }
        }
        } // end profile

        {
        BL_PROFILE("rhs_post_9");
        // This updates all the conserved variables (not just the "slow" ones)
        int   num_comp_all = S_data[IntVar::cons].nComp();
        ParallelFor(tbx, num_comp_all,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            new_cons(i,j,k,n)  = cur_cons(i,j,k,n);
        });
        } // end profile

        Box gtbx = mfi.tilebox(IntVect(1,0,0),S_old[IntVar::xmom].nGrowVect());
        Box gtby = mfi.tilebox(IntVect(0,1,0),S_old[IntVar::ymom].nGrowVect());
        Box gtbz = mfi.tilebox(IntVect(0,0,1),S_old[IntVar::zmom].nGrowVect());

        {
        BL_PROFILE("rhs_post_10()");
        ParallelFor(gtbx, gtby, gtbz,
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
    } // mfi

    if (l_use_diff) {
        delete dflux_x;
        delete dflux_y;
        delete dflux_z;
    }
}
