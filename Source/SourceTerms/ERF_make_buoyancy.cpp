#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_GpuContainers.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_PlaneAverage.H>
#include <ERF_Src_headers.H>
#include <ERF_buoyancy_utils.H>

using namespace amrex;

/**
 * Function for computing the buoyancy term to be used in the evolution
 * equation for the z-component of momentum in the slow integrator.  There
 * are three options for how buoyancy is computed (two are the same in the absence of moisture).
 *
 * @param[in]  S_data        current solution
 * @param[in]  S_prim        primitive variables (i.e. conserved variables divided by density)
 * @param[out] buoyancy      buoyancy term computed here
 * @param[in]  qmoist        moisture variables (in order: qv, qc, qi, ...)
 * @param[in]  qv_d          lateral average of cloud vapor
 * @param[in]  qc_d          lateral average of cloud vapor
 * @param[in]  qd_d          lateral average of cloud vapor
 * @param[in]  geom          Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in]  r0            Reference (hydrostatically stratified) density
 * @param[in]  n_qstate      Number of moist variables used by the current model
 */

void make_buoyancy (Vector<MultiFab>& S_data,
                    const MultiFab& S_prim,
                          MultiFab& buoyancy,
                    const amrex::Geometry geom,
                    const SolverChoice& solverChoice,
                    const MultiFab& base_state,
                    const int n_qstate,
                    const int anelastic)
{
    BL_PROFILE("make_buoyancy()");

    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const int klo = geom.Domain().smallEnd()[2];
    const int khi = geom.Domain().bigEnd()[2] + 1;

    Real rd_over_cp = solverChoice.rdOcp;
    Real rv_over_rd = R_v/R_d;

    MultiFab r0 (base_state, make_alias, BaseState::r0_comp , 1);
    MultiFab p0 (base_state, make_alias, BaseState::p0_comp , 1);
    MultiFab th0(base_state, make_alias, BaseState::th0_comp, 1);

    if (anelastic == 1) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbz = mfi.tilebox();

            // We don't compute a source term for z-momentum on the bottom or top boundary
            if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
            if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

            const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
            const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

            // Base state density and pressure
            const Array4<const Real>&  r0_arr =  r0.const_array(mfi);
            const Array4<const Real>& th0_arr = th0.const_array(mfi);

            if (solverChoice.moisture_type == MoistureType::None) {
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    //
                    // Return -rho0 g (thetaprime / theta0)
                    //
                    buoyancy_fab(i, j, k) = buoyancy_dry_anelastic(i,j,k,
                                                                   grav_gpu[2],
                                                                   r0_arr,th0_arr,cell_data);
                });
            } else {
                // NOTE: For decomposition in the vertical direction, klo may not
                //       reside in the valid box and this call will yield an out
                //       of bounds error since it depends upon the surface theta_l
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    //
                    // Return -rho0 g (thetaprime / theta0)
                    //
                    buoyancy_fab(i, j, k) = buoyancy_moist_anelastic(i,j,k,
                                                                     grav_gpu[2],rv_over_rd,
                                                                     r0_arr,th0_arr,cell_data);
                });
            }
        } // mfi
    }
    else
    {
        // ******************************************************************************************
        // Dry versions of buoyancy expressions (type 1 and type 2/3 -- types 2 and 3 are equivalent)
        // ******************************************************************************************
        if (solverChoice.moisture_type == MoistureType::None)
        {
            int n_q_dry = 0;
            if (solverChoice.buoyancy_type == 1) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Box tbz = mfi.tilebox();

                    // We don't compute a source term for z-momentum on the bottom or top domain boundary
                    if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                    if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                    const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                    const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                    // Base state density
                    const Array4<const Real>& r0_arr = r0.const_array(mfi);

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_rhopert(i,j,k,n_q_dry,grav_gpu[2],r0_arr,cell_data);
                    });
                } // mfi

            }
            else // (buoyancy_type != 1)
            {
                // We now use the base state rather than planar average because
                //     1) we don't want to average over the limited region of the fine level if doing multilevel.
                //     2) it's cheaper to use the base state than to compute the horizontal averages
                //     3) when running in a smallish domain, the horizontal average may evolve over time,
                //        which is not necessarily the intended behavior
                //
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Box tbz = mfi.tilebox();

                    // We don't compute a source term for z-momentum on the bottom or top boundary
                    if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                    if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                    // Base state density and pressure
                    const Array4<const Real>&  r0_arr = r0.const_array(mfi);
                    const Array4<const Real>&  p0_arr = p0.const_array(mfi);
                    const Array4<const Real>& th0_arr = th0.const_array(mfi);

                    const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                    const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_dry_default(i,j,k,
                                                                     grav_gpu[2],rd_over_cp,
                                                                     r0_arr,p0_arr,th0_arr,cell_data);
                    });
                } // mfi
            } // buoyancy_type
        } // moisture type
        else
        {
        // ******************************************************************************************
        // Moist versions of buoyancy expressions
        // ******************************************************************************************

          if ( (solverChoice.moisture_type == MoistureType::Kessler_NoRain) ||
               (solverChoice.moisture_type == MoistureType::SAM)            ||
               (solverChoice.moisture_type == MoistureType::SAM_NoPrecip_NoIce) )
          {
              AMREX_ALWAYS_ASSERT(solverChoice.buoyancy_type == 1);
          }

          if (solverChoice.buoyancy_type == 1) {

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box tbz = mfi.tilebox();

                // We don't compute a source term for z-momentum on the bottom or top domain boundary
                if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                // Base state density
                const Array4<const Real>& r0_arr = r0.const_array(mfi);

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    buoyancy_fab(i, j, k) = buoyancy_rhopert(i,j,k,n_qstate,
                                                           grav_gpu[2],r0_arr,cell_data);
                });
            } // mfi

          } else {

            PlaneAverage state_ave(&(S_data[IntVars::cons]), geom, solverChoice.ave_plane);
            PlaneAverage  prim_ave(&S_prim                 , geom, solverChoice.ave_plane);

            // Compute horizontal averages of all components of each field
            state_ave.compute_averages(ZDir(), state_ave.field());
            prim_ave.compute_averages(ZDir(), prim_ave.field());

            int ncell = state_ave.ncell_line();

            Gpu::HostVector  <Real> rho_h(ncell), theta_h(ncell);
            Gpu::DeviceVector<Real> rho_d(ncell), theta_d(ncell);

            state_ave.line_average(Rho_comp, rho_h);
            Gpu::copyAsync(Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());

            prim_ave.line_average(PrimTheta_comp, theta_h);
            Gpu::copyAsync(Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());

            Real*   rho_d_ptr =   rho_d.data();
            Real* theta_d_ptr = theta_d.data();

            // Average valid moisture vars
            Gpu::HostVector  <Real> qv_h(ncell)    , qc_h(ncell)    , qp_h(ncell);
            Gpu::DeviceVector<Real> qv_d(ncell,0.0), qc_d(ncell,0.0), qp_d(ncell,0.0);
            if (n_qstate >=1) {
                prim_ave.line_average(PrimQ1_comp, qv_h);
               Gpu::copyAsync(Gpu::hostToDevice,  qv_h.begin(), qv_h.end(), qv_d.begin());
            }
            if (n_qstate >=2) {
                prim_ave.line_average(PrimQ2_comp, qc_h);
                Gpu::copyAsync(Gpu::hostToDevice,  qc_h.begin(), qc_h.end(), qc_d.begin());
            }
            if (n_qstate >=3) {
                prim_ave.line_average(PrimQ3_comp, qp_h);
                Gpu::copyAsync(Gpu::hostToDevice,  qp_h.begin(), qp_h.end(), qp_d.begin());
            }
            Real* qv_d_ptr = qv_d.data();
            Real* qc_d_ptr = qc_d.data();
            Real* qp_d_ptr = qp_d.data();

            if (solverChoice.buoyancy_type == 2 || solverChoice.buoyancy_type == 4 ) {

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Box tbz = mfi.tilebox();

                    // We don't compute a source term for z-momentum on the bottom or top domain boundary
                    if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                    if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                    const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                    const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                    const Array4<const Real> & cell_prim  = S_prim.array(mfi);

                    // TODO: ice has not been dealt with (q1=qv, q2=qv, q3=qp)
                    if (solverChoice.buoyancy_type == 2) {
                        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            buoyancy_fab(i, j, k) = buoyancy_type2(i,j,k,n_qstate,grav_gpu[2],
                                                                   rho_d_ptr,theta_d_ptr,
                                                                   qv_d_ptr,qc_d_ptr,qp_d_ptr,
                                                                   cell_prim,cell_data);
                        });
                    } else {
                        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            buoyancy_fab(i, j, k) = buoyancy_type4(i,j,k,n_qstate,grav_gpu[2],
                                                                   rho_d_ptr,theta_d_ptr,
                                                                   qv_d_ptr,qc_d_ptr,qp_d_ptr,
                                                                   cell_prim,cell_data);
                        });
                    }
                } // mfi

            } else if (solverChoice.buoyancy_type == 3) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Box tbz = mfi.tilebox();

                    // We don't compute a source term for z-momentum on the bottom or top domain boundary
                    if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                    if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                    const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                    const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                    const Array4<const Real> & cell_prim  = S_prim.array(mfi);

                    // TODO: ice has not been dealt with (q1=qv, q2=qv, q3=qp)

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_type3(i,j,k,n_qstate,grav_gpu[2],
                                                               rho_d_ptr,theta_d_ptr,qv_d_ptr,
                                                               cell_prim,cell_data);
                    });
                } // mfi
            }  // buoyancy_type
          } // not buoyancy_type == 1
        } // has moisture
    } // anelastic?
}
