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
                        Vector<MultiFab>& S_rhs,
                        Vector<MultiFab>& S_old,
                        Vector<MultiFab>& S_new,
                        Vector<MultiFab>& S_data,
                        const MultiFab& S_prim,
                              Vector<MultiFab>& S_scratch,
                        const MultiFab& xvel,
                        const MultiFab& yvel,
                        const MultiFab& zvel,
                        const MultiFab& source, const MultiFab& eddyDiffs,
                        const amrex::Geometry geom,
                        const SolverChoice& solverChoice,
                        std::unique_ptr<ABLMost>& most,
                        const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                        std::unique_ptr<MultiFab>& dJ,
                        std::unique_ptr<MultiFab>& dJ_new)
{
    BL_PROFILE_REGION("erf_slow_rhs_post()");

    amrex::Real theta_mean;
    if (most) theta_mean = most->theta_mean;

    const int l_spatial_order   = solverChoice.spatial_order;
    const bool l_use_terrain    = solverChoice.use_terrain;
    const bool l_moving_terrain = (solverChoice.terrain_type == 1);

    const bool l_use_QKE        = solverChoice.use_QKE && solverChoice.advect_QKE;
    const bool l_use_deardorff  = (solverChoice.les_type == LESType::Deardorff);
    const bool l_use_diff       = ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
                                    (solverChoice.les_type        !=       LESType::None) ||
                                    (solverChoice.pbl_type        !=       PBLType::None) );

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

        const Box& bx = mfi.tilebox();

        const Array4<      Real> & old_cons   = S_old[IntVar::cons].array(mfi);
        const Array4<      Real> & cell_rhs   = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & source_fab = source.const_array(mfi);

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

        const Array4<Real const>& K_turb = eddyDiffs.const_array(mfi);

        const Array4<const Real>& detJ     = l_use_terrain ? dJ->const_array(mfi)     : Array4<const Real>{};
        const Array4<const Real>& detJ_new = l_use_terrain ? dJ_new->const_array(mfi) : Array4<const Real>{};

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
        if (l_use_deardorff) {
            int start_comp = RhoKE_comp;
            int   num_comp = 1;
            AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                                   cur_prim, cell_rhs, detJ,
                                   dxInv, l_spatial_order, l_use_terrain);
        }
        if (l_use_QKE) {
            int start_comp = RhoQKE_comp;
            int   num_comp = 1;
            AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                                   cur_prim, cell_rhs, detJ,
                                   dxInv, l_spatial_order, l_use_terrain);
        }
        int start_comp = RhoScalar_comp;
        int   num_comp = S_data[IntVar::cons].nComp() - start_comp;
        AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                               cur_prim, cell_rhs, detJ,
                               dxInv, l_spatial_order, l_use_terrain);

        if (l_use_diff && !l_use_terrain) {
            Array4<Real> diffflux_x = dflux_x->array(mfi);
            Array4<Real> diffflux_y = dflux_y->array(mfi);
            Array4<Real> diffflux_z = dflux_z->array(mfi);

            // NOTE: No diffusion for continuity, so n starts at 1.
            //       KE calls moved inside DiffSrcForState.
            int n_start = amrex::max(start_comp,RhoTheta_comp);
            int n_end   = start_comp + num_comp - 1;
            DiffusionSrcForState_N(bx, domain, n_start, n_end, u, v, w,
                                   cur_cons, cur_prim, source_fab, cell_rhs,
                                   diffflux_x, diffflux_y, diffflux_z,
                                   dxInv, K_turb, solverChoice, theta_mean, grav_gpu, bc_ptr);
        }

        // This updates just the "slow" conserved variables
        {
        BL_PROFILE("rhs_post_8");

        if ( solverChoice.use_terrain && solverChoice.terrain_type == 1 )
        {
            ParallelFor(bx, num_comp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                cur_cons(i,j,k,n) =  old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
                cur_cons(i,j,k,n) *= detJ(i,j,k) / detJ_new(i,j,k);
            });
        } else {
            ParallelFor(bx, num_comp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int nn) noexcept {
                const int n = start_comp + nn;
                cur_cons(i,j,k,n) = old_cons(i,j,k,n) + dt * cell_rhs(i,j,k,n);
            });
        }
        } // end profile

        {
        BL_PROFILE("rhs_post_9");
        // This updates all the conserved variables (not just the "slow" ones)
        int   num_comp_all = S_data[IntVar::cons].nComp();
        ParallelFor(bx, num_comp_all,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            new_cons(i,j,k,n)  = cur_cons(i,j,k,n);
        });
        } // end profile

        const Box gtbx = mfi.nodaltilebox(0).grow(S_old[IntVar::xmom].nGrowVect());
        const Box gtby = mfi.nodaltilebox(1).grow(S_old[IntVar::ymom].nGrowVect());
        const Box gtbz = mfi.nodaltilebox(2).grow(S_old[IntVar::zmom].nGrowVect());

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
}
