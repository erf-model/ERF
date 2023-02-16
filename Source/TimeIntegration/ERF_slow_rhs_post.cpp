#include <AMReX.H>
#include <AMReX_MultiFab.H>
//#include <AMReX_ArrayLim.H>
//#include <AMReX_BCRec.H>
//#include <ERF_Constants.H>
//#include <ABLMost.H>
#include <Advection.H>
#include <Diffusion.H>
#include <TimeIntegration.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>

using namespace amrex;

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
                        const amrex::Geometry geom,
                        const SolverChoice& solverChoice,
                        std::unique_ptr<ABLMost>& most,
                        const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                        std::unique_ptr<MultiFab>& z_phys_nd,
                        std::unique_ptr<MultiFab>& dJ,
                        std::unique_ptr<MultiFab>& dJ_new,
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
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        // Metric terms
        const Array4<const Real>& z_nd     = l_use_terrain    ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ     = l_use_terrain    ? dJ->const_array(mfi)        : Array4<const Real>{};
        const Array4<const Real>& detJ_new = l_moving_terrain ? dJ_new->const_array(mfi)    : Array4<const Real>{};

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
        ParallelFor(valid_bx, ncomp_slow[IntVar::cons],
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) {
            const int n = scomp_slow[IntVar::cons] + nn;
            cur_cons(i,j,k,n) = new_cons(i,j,k,n);
        });
        } // end profile

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        if (l_use_deardorff) {
            int start_comp = RhoKE_comp;
            int   num_comp = 1;
            AdvectionSrcForScalars(valid_bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                                   cur_prim, cell_rhs, detJ,
                                   dxInv, mf_m, l_all_WENO, l_moist_WENO, l_spatial_order_WENO,
                                   l_horiz_spatial_order, l_vert_spatial_order, l_use_terrain);
        }
        if (l_use_QKE) {
            int start_comp = RhoQKE_comp;
            int   num_comp = 1;
            AdvectionSrcForScalars(valid_bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                                   cur_prim, cell_rhs, detJ,
                                   dxInv, mf_m, l_all_WENO, l_moist_WENO, l_spatial_order_WENO,
                                   l_horiz_spatial_order, l_vert_spatial_order, l_use_terrain);
        }
        int start_comp = RhoScalar_comp;
        int   num_comp = S_data[IntVar::cons].nComp() - start_comp;
        AdvectionSrcForScalars(valid_bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                               cur_prim, cell_rhs, detJ,
                               dxInv, mf_m, l_all_WENO, l_moist_WENO, l_spatial_order_WENO,
                               l_horiz_spatial_order, l_vert_spatial_order, l_use_terrain);

        if (l_use_diff) {
            Array4<Real> diffflux_x = dflux_x->array(mfi);
            Array4<Real> diffflux_y = dflux_y->array(mfi);
            Array4<Real> diffflux_z = dflux_z->array(mfi);

            const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

            // NOTE: No diffusion for continuity, so n starts at 1.
            //       KE calls moved inside DiffSrcForState.
            int n_start = amrex::max(start_comp,RhoTheta_comp);
            int n_end   = start_comp + num_comp - 1;

            if (l_use_terrain) {
                DiffusionSrcForState_T(valid_bx, domain, n_start, n_end, u, v, w,
                                       cur_cons, cur_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z, z_nd, detJ,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
            } else {
                DiffusionSrcForState_N(valid_bx, domain, n_start, n_end, u, v, w,
                                       cur_cons, cur_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       mu_turb, solverChoice, tm_arr, grav_gpu, bc_ptr);
            }
        }

        // This updates just the "slow" conserved variables
        {
        BL_PROFILE("rhs_post_8");

        if (l_moving_terrain)
        {
            auto const& src_arr = source.const_array(mfi);
            ParallelFor(valid_bx, num_comp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                // NOTE: we don't include additional source terms when terrain is moving
                Real temp_val = detJ(i,j,k) * old_cons(i,j,k,n) + dt * detJ(i,j,k) * cell_rhs(i,j,k,n);
                cur_cons(i,j,k,n) = temp_val / detJ_new(i,j,k);
            });
        } else {
            auto const& src_arr = source.const_array(mfi);
            ParallelFor(valid_bx, num_comp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                cell_rhs(i,j,k,n) += src_arr(i,j,k,n);
                cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
            });
        }
        } // end profile

        {
        BL_PROFILE("rhs_post_9");
        // This updates all the conserved variables (not just the "slow" ones)
        int   num_comp_all = S_data[IntVar::cons].nComp();
        ParallelFor(valid_bx, num_comp_all,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            new_cons(i,j,k,n)  = cur_cons(i,j,k,n);
        });
        } // end profile

        // const Box gtbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xmom].nGrowVect());
        // const Box gtby = mfi.nodaltilebox(1).grow(S_old[IntVar::ymom].nGrowVect());
        // const Box gtbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zmom].nGrowVect());
        Box gtbx = surroundingNodes(valid_bx,0);  gtbx.grow(S_old[IntVar::xmom].nGrowVect());
        Box gtby = surroundingNodes(valid_bx,1);  gtby.grow(S_old[IntVar::ymom].nGrowVect());
        Box gtbz = surroundingNodes(valid_bx,2);  gtbz.grow(S_old[IntVar::zmom].nGrowVect());

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
