#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_GpuContainers.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_SrcHeaders.H>
#include <ERF_BuoyancyUtils.H>
#include <ERF_EB.H>

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

void make_buoyancy (int lev,
                    const Vector<MultiFab>& S_data,
                    const MultiFab& S_prim,
                    const MultiFab& qt,
                          MultiFab& buoyancy,
                    const Geometry geom,
                    const SolverChoice& solverChoice,
                    const MultiFab& base_state,
                    const int n_qstate,
                    const eb_& ebfact,
                    const int anelastic)
{
    BL_PROFILE("make_buoyancy()");

    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const int klo = geom.Domain().smallEnd()[2];
    const int khi = geom.Domain().bigEnd()[2] + 1;

    Real rd_over_cp = solverChoice.rdOcp;
    //Real rv_over_rd = R_v/R_d;

    MultiFab r0 (base_state, make_alias, BaseState::r0_comp , 1);
    MultiFab p0 (base_state, make_alias, BaseState::p0_comp , 1);
    MultiFab th0(base_state, make_alias, BaseState::th0_comp, 1);
    MultiFab qv0(base_state, make_alias, BaseState::qv0_comp, 1);

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
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<const Real> &    qt_arr  = qt.array(mfi);
        const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

        // Base state density and pressure
        const Array4<const Real>&  r0_arr =  r0.const_array(mfi);
        const Array4<const Real>&  p0_arr =  p0.const_array(mfi);
        const Array4<const Real>& th0_arr = th0.const_array(mfi);
        const Array4<const Real>& qv0_arr = qv0.const_array(mfi);

        if (solverChoice.terrain_type != TerrainType::EB) {

            if ( anelastic && (solverChoice.moisture_type == MoistureType::None) )
            {
                // ******************************************************************************************
                // Dry anelastic
                // ******************************************************************************************
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    //
                    // Return -rho0 g (thetaprime / theta0)
                    //
                    buoyancy_fab(i, j, k) = buoyancy_dry_anelastic(i,j,k,grav_gpu[2],
                                                                r0_arr,th0_arr,cell_data);
                });
            }
            else if ( anelastic && (solverChoice.moisture_type != MoistureType::None) )
            {
                // ******************************************************************************************
                // Moist anelastic
                // ******************************************************************************************
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    //
                    // Return -rho0 g (thetaprime / theta0)
                    //
                    //buoyancy_fab(i, j, k) = buoyancy_moist_anelastic(i,j,k,grav_gpu[2],rv_over_rd,
                    //                                                 r0_arr,th0_arr,qv0_arr,cell_data,qt_arr);

                    // NOTE: Using the type 4, which we formally derived.
                    //       The above has errors and needs rederiving.
                    buoyancy_fab(i, j, k) = buoyancy_moist_Thpert(i,j,k,n_qstate,grav_gpu[2],
                                                                    r0_arr,th0_arr,qv0_arr,cell_prim,qt_arr);
                });
            }
            else if ( !anelastic && (solverChoice.moisture_type == MoistureType::None) )
            {
                // ******************************************************************************************
                // Dry compressible
                // ******************************************************************************************
                if (solverChoice.buoyancy_type[lev] == 1) {

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        //
                        // Return -rho0 g (thetaprime / theta0)
                        //
                        buoyancy_fab(i, j, k) = buoyancy_rhopert(i,j,k,grav_gpu[2],
                                                                r0_arr,cell_data,qt_arr);
                    });
                }
                else if (solverChoice.buoyancy_type[lev] == 2 || solverChoice.buoyancy_type[lev] == 3)
                {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        //
                        // Return -rho0 g (Tprime / T0)
                        //
                        buoyancy_fab(i, j, k) = buoyancy_dry_Tpert(i,j,k,grav_gpu[2],rd_over_cp,
                                                                r0_arr,p0_arr,th0_arr,cell_data);
                    });
                }
                else if (solverChoice.buoyancy_type[lev] == 4)
                {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        //
                        // Return -rho0 g (Theta_prime / Theta_0)
                        //
                        buoyancy_fab(i, j, k) = buoyancy_dry_Thpert(i,j,k,grav_gpu[2],
                                                                    r0_arr,th0_arr,cell_data);
                    });
                } // buoyancy_type for dry compressible
            }
            else // if ( !anelastic && (solverChoice.moisture_type != MoistureType::None) )
            {
                // ******************************************************************************************
                // Moist compressible
                // ******************************************************************************************

                if ( (solverChoice.moisture_type == MoistureType::Kessler_NoRain) ||
                    (solverChoice.moisture_type == MoistureType::SAM)            ||
                    (solverChoice.moisture_type == MoistureType::SAM_NoPrecip_NoIce) )
                {
                    AMREX_ALWAYS_ASSERT(solverChoice.buoyancy_type[lev] == 1);
                }

                if (solverChoice.buoyancy_type[lev] == 1)
                {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_rhopert(i,j,k,grav_gpu[2],
                                                                r0_arr,cell_data,qt_arr);
                    });
                }
                else if (solverChoice.buoyancy_type[lev] == 2 || solverChoice.buoyancy_type[lev] == 3)
                {

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_moist_Tpert(i,j,k,n_qstate,grav_gpu[2],rd_over_cp,
                                                                    r0_arr,th0_arr,qv0_arr,p0_arr,
                                                                    cell_prim,cell_data,qt_arr);
                    });
                }
                else if (solverChoice.buoyancy_type[lev] == 4)
                {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_moist_Thpert(i,j,k,n_qstate,grav_gpu[2],
                                                                    r0_arr,th0_arr,qv0_arr,cell_prim,qt_arr);
                    });
                }
            } // moist compressible

        } else {

            if ( anelastic && (solverChoice.moisture_type == MoistureType::None) ) {

                if (grav_gpu[2]==0) {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = 0.0;
                    });
                } else {

                    Array4<const EBCellFlag> cellflg = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_dry_anelastic_eb(i,j,k,grav_gpu[2],
                                                                    r0_arr,th0_arr,cell_data,cellflg);
                    });
                }
            }
            else
            {
                if (grav_gpu[2]==0) {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = 0.0;
                    });
                } else {
                    // Currently, only dry compressible is supported
                    AMREX_ASSERT( !anelastic && (solverChoice.moisture_type == MoistureType::None) && solverChoice.buoyancy_type[lev] == 1 );

                    Array4<const EBCellFlag> cellflg = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();

                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        buoyancy_fab(i, j, k) = buoyancy_rhopert_eb(i,j,k,grav_gpu[2],r0_arr,cell_data,qt_arr,cellflg);
                    });
                }
            }
        } // TerrainType::EB
    } // mfi
}
