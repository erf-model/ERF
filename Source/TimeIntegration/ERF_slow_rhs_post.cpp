#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <ERF_Constants.H>
#include <ABLMost.H>
#include <SpatialStencils.H>
#include <TimeIntegration.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>

using namespace amrex;

void erf_slow_rhs_post (int level,
                        Vector<MultiFab>& S_rhs,
                        const Vector<MultiFab>& S_data,
                        const MultiFab& S_prim,
                              Vector<MultiFab>& S_scratch,
                        const MultiFab& xvel,
                        const MultiFab& yvel,
                        const MultiFab& zvel,
                        std::unique_ptr<MultiFab>& z_t_mf,
                        const MultiFab& source, const MultiFab& eddyDiffs,
                        std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
                        const amrex::Geometry geom,
                        const SolverChoice& solverChoice,
                        std::unique_ptr<ABLMost>& most,
                        const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                        std::unique_ptr<MultiFab>& z0, std::unique_ptr<MultiFab>& dJ)
{
    BL_PROFILE_REGION("erf_slow_rhs_post()");

    amrex::Real theta_mean;
    if (most) theta_mean = most->theta_mean;

    const int l_spatial_order = solverChoice.spatial_order;
    const int l_use_terrain   = solverChoice.use_terrain;

    const amrex::BCRec* bc_ptr = domain_bcs_type_d.data();

    const Box& domain = geom.Domain();
    const int domhi_z = domain.bigEnd()[2];

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

    bool l_use_QKE       = solverChoice.use_QKE && solverChoice.advect_QKE;
    bool l_use_deardorff = (solverChoice.les_type == LESType::Deardorff);

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(S_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<      Real> & cell_rhs   = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & source_fab = source.const_array(mfi);

        Array4<Real> avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        Array4<const Real> z_t;
        if (z_t_mf)
            z_t = z_t_mf->array(mfi);
        else
            z_t = Array4<const Real>{};

        // These are temporaries we use to add to the S_rhs for the fluxes
        const Array4<Real>& diffflux_x = diffflux[0].array(mfi);
        const Array4<Real>& diffflux_y = diffflux[1].array(mfi);
        const Array4<Real>& diffflux_z = diffflux[2].array(mfi);

        const Array4<Real const>& K_turb = eddyDiffs.const_array(mfi);

        const Array4<const Real>& detJ = l_use_terrain ? dJ->const_array(mfi) : Array4<const Real>{};

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************

        int start_comp = 2;
        int   num_comp = S_data[IntVar::cons].nComp() - 2;
        AdvectionSrcForScalars(bx, start_comp, num_comp, avg_xmom, avg_ymom, avg_zmom,
                               cell_prim, cell_rhs, detJ,
                               dxInv, l_spatial_order, l_use_terrain, l_use_deardorff, l_use_QKE);

        // NOTE: No diffusion for continuity, so n starts at 1.
        //       KE calls moved inside DiffSrcForState.
        int n_start = amrex::max(start_comp,RhoTheta_comp);
        int n_end   = start_comp + num_comp - 1;
        DiffusionSrcForState(bx, domain, n_start, n_end, u, v, w,
                             cell_data, cell_prim, source_fab, cell_rhs,
                             diffflux_x, diffflux_y, diffflux_z,
                             dxInv, K_turb, solverChoice, theta_mean, grav_gpu, bc_ptr);
    } // mfi
}
