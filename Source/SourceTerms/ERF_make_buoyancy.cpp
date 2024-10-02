#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_GpuContainers.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_PlaneAverage.H>
#include <ERF_Src_headers.H>
//#include <ERF_TerrainMetrics.H>

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
                    const MultiFab* r0,
                    const MultiFab* p0,
                    const int n_qstate,
                    const int anelastic)
{
    BL_PROFILE("make_buoyancy()");

    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const int klo = geom.Domain().smallEnd()[2];
    const int khi = geom.Domain().bigEnd()[2] + 1;

    Real rd_over_cp = solverChoice.rdOcp;

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

            // Base state density, pressure
            const Array4<const Real>& r0_arr = r0->const_array(mfi);
            const Array4<const Real>& p0_arr = p0->const_array(mfi);

            ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rt0_hi = getRhoThetagivenP(p0_arr(i,j,k));
                Real  t0_hi = getTgivenPandTh(p0_arr(i,j,k), rt0_hi/r0_arr(i,j,k), rd_over_cp);
                Real   t_hi = getTgivenPandTh(p0_arr(i,j,k), cell_data(i,j,k,RhoTheta_comp)/r0_arr(i,j,k), rd_over_cp);
                Real qplus   = (t_hi-t0_hi)/t0_hi;

                Real rt0_lo = getRhoThetagivenP(p0_arr(i,j,k-1));
                Real  t0_lo = getTgivenPandTh(p0_arr(i,j,k-1), rt0_lo/r0_arr(i,j,k-1), rd_over_cp);
                Real   t_lo = getTgivenPandTh(p0_arr(i,j,k-1), cell_data(i,j,k-1,RhoTheta_comp)/r0_arr(i,j,k-1), rd_over_cp);
                Real qminus  = (t_lo-t0_lo)/t0_lo;

                Real r0_q_avg = Real(0.5) * (r0_arr(i,j,k) * qplus + r0_arr(i,j,k-1) * qminus);
                buoyancy_fab(i, j, k) = -r0_q_avg * grav_gpu[2];
            });
        } // mfi
    }
    else
    {
        // ******************************************************************************************
        // Dry versions of buoyancy expressions (type 1 and type 2/3 -- types 2 and 3 are equivalent)
        // ******************************************************************************************
        if (solverChoice.moisture_type == MoistureType::None)
        {

          if (solverChoice.buoyancy_type == 1) {
            for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box tbz = mfi.tilebox();

                // We don't compute a source term for z-momentum on the bottom or top domain boundary
                if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                // Base state density
                const Array4<const Real>& r0_arr = r0->const_array(mfi);

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    buoyancy_fab(i, j, k) = grav_gpu[2] * 0.5 * ( (cell_data(i,j,k  ) - r0_arr(i,j,k  ))
                                                                 +(cell_data(i,j,k-1) - r0_arr(i,j,k-1)) );
                });
            } // mfi

          }
          else if (solverChoice.buoyancy_type == 2 || solverChoice.buoyancy_type == 3)
          {
            PlaneAverage state_ave(&(S_data[IntVars::cons]), geom, solverChoice.ave_plane);
            PlaneAverage prim_ave(&S_prim, geom, solverChoice.ave_plane);

            int ncell = state_ave.ncell_line();

            state_ave.compute_averages(ZDir(), state_ave.field());
            prim_ave.compute_averages(ZDir(), prim_ave.field());

            Gpu::HostVector<Real> rho_h(ncell), theta_h(ncell);
            state_ave.line_average(Rho_comp, rho_h);
            prim_ave.line_average(PrimTheta_comp, theta_h);

            Gpu::DeviceVector<Real>   rho_d(ncell);
            Gpu::DeviceVector<Real> theta_d(ncell);

            Gpu::copyAsync(Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
            Gpu::copyAsync(Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());

            Real*   rho_d_ptr =   rho_d.data();
            Real* theta_d_ptr = theta_d.data();

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

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ]);
                    Real tempp3d = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));
                    Real qplus   = (tempp3d-tempp1d)/tempp1d;

                    Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1]);
                    Real tempm3d = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));
                    Real qminus  = (tempm3d-tempm1d)/tempm1d;

                    Real r0_q_avg = Real(0.5) * (rho_d_ptr[k]*qplus + rho_d_ptr[k-1]*qminus);

                    buoyancy_fab(i, j, k) = -r0_q_avg * grav_gpu[2];
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

            for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box tbz = mfi.tilebox();

                // We don't compute a source term for z-momentum on the bottom or top domain boundary
                if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
                if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

                const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
                const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

                // Base state density
                const Array4<const Real>& r0_arr = r0->const_array(mfi);

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rhop_hi = cell_data(i,j,k  ,Rho_comp) + cell_data(i,j,k  ,RhoQ1_comp)
                                 + cell_data(i,j,k  ,RhoQ2_comp) - r0_arr(i,j,k  );
                    Real rhop_lo = cell_data(i,j,k-1,Rho_comp) + cell_data(i,j,k-1,RhoQ1_comp)
                                 + cell_data(i,j,k-1,RhoQ2_comp) - r0_arr(i,j,k-1);

                    for (int q_offset(2); q_offset<n_qstate; ++q_offset) {
                        rhop_hi += cell_data(i,j,k  ,RhoQ1_comp+q_offset);
                        rhop_lo += cell_data(i,j,k-1,RhoQ1_comp+q_offset);
                    }
                    buoyancy_fab(i, j, k) = grav_gpu[2] * 0.5 * ( rhop_hi + rhop_lo );
                });
            } // mfi

          } else {

            PlaneAverage state_ave(&(S_data[IntVars::cons]), geom, solverChoice.ave_plane);
            PlaneAverage  prim_ave(&S_prim                , geom, solverChoice.ave_plane);

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
                    ParallelFor(tbz, [=, buoyancy_type=solverChoice.buoyancy_type] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ], qv_d_ptr[k  ]);
                        Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1], qv_d_ptr[k-1]);

                        Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp),
                                                         cell_data(i,j,k  ,RhoTheta_comp),
                                                         cell_data(i,j,k  ,RhoQ1_comp)/cell_data(i,j,k  ,Rho_comp));
                        Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp),
                                                         cell_data(i,j,k-1,RhoTheta_comp),
                                                         cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp));

                        Real qplus, qminus;

                        Real qv_plus  = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
                        Real qv_minus = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

                        Real qc_plus  = (n_qstate >= 2) ? cell_prim(i,j,k  ,PrimQ2_comp) : 0.0;
                        Real qc_minus = (n_qstate >= 2) ? cell_prim(i,j,k-1,PrimQ2_comp) : 0.0;

                        Real qp_plus  = (n_qstate >= 3) ? cell_prim(i,j,k  ,PrimQ3_comp) : 0.0;
                        Real qp_minus = (n_qstate >= 3) ? cell_prim(i,j,k-1,PrimQ3_comp) : 0.0;

                        if (buoyancy_type == 2) {
                            qplus  = 0.61 * ( qv_plus - qv_d_ptr[k] ) -
                                            ( qc_plus - qc_d_ptr[k]   +
                                              qp_plus - qp_d_ptr[k] )
                                   + (tempp3d-tempp1d)/tempp1d*(Real(1.0) + Real(0.61)*qv_d_ptr[k]-qc_d_ptr[k]-qp_d_ptr[k]);

                            qminus = 0.61 * ( qv_minus - qv_d_ptr[k-1] ) -
                                            ( qc_minus - qc_d_ptr[k-1]   +
                                              qp_minus - qp_d_ptr[k-1] )
                                   + (tempm3d-tempm1d)/tempm1d*(Real(1.0) + Real(0.61)*qv_d_ptr[k-1]-qc_d_ptr[k-1]-qp_d_ptr[k-1]);

                        } else { // (buoyancy_type == 4)
                            qplus  = 0.61 * ( qv_plus - qv_d_ptr[k] ) -
                                            ( qc_plus - qc_d_ptr[k]   +
                                              qp_plus - qp_d_ptr[k] )
                                   + (cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp) - theta_d_ptr[k  ])/theta_d_ptr[k  ];

                            qminus = 0.61 * ( qv_minus - qv_d_ptr[k-1] ) -
                                            ( qc_minus - qc_d_ptr[k-1]   +
                                              qp_minus - qp_d_ptr[k-1] )
                                   + (cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp) - theta_d_ptr[k-1])/theta_d_ptr[k-1];
                        }

                        Real qavg  = Real(0.5) * (qplus + qminus);
                        Real r0avg = Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

                        buoyancy_fab(i, j, k) = -qavg * r0avg * grav_gpu[2];
                    });
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
                        Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ], qv_d_ptr[k  ]);
                        Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1], qv_d_ptr[k-1]);

                        Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp),
                                                         cell_data(i,j,k  ,RhoTheta_comp),
                                                         cell_data(i,j,k  ,RhoQ1_comp)/cell_data(i,j,k  ,Rho_comp));
                        Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp),
                                                         cell_data(i,j,k-1,RhoTheta_comp),
                                                         cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp));

                        Real qv_plus  = (n_qstate >= 1) ? cell_prim(i,j,k  ,PrimQ1_comp) : 0.0;
                        Real qv_minus = (n_qstate >= 1) ? cell_prim(i,j,k-1,PrimQ1_comp) : 0.0;

                        Real qc_plus  = (n_qstate >= 2) ? cell_prim(i,j,k  ,PrimQ2_comp) : 0.0;
                        Real qc_minus = (n_qstate >= 2) ? cell_prim(i,j,k-1,PrimQ2_comp) : 0.0;

                        Real qp_plus  = (n_qstate >= 3) ? cell_prim(i,j,k  ,PrimQ3_comp) : 0.0;
                        Real qp_minus = (n_qstate >= 3) ? cell_prim(i,j,k-1,PrimQ3_comp) : 0.0;

                        Real qplus  = 0.61 * qv_plus  - (qc_plus  + qp_plus)  + (tempp3d-tempp1d)/tempp1d;

                        Real qminus = 0.61 * qv_minus - (qc_minus + qp_minus) + (tempm3d-tempm1d)/tempm1d;

                        Real qavg  = Real(0.5) * (qplus + qminus);
                        Real r0avg = Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

                        buoyancy_fab(i, j, k) = -qavg * r0avg * grav_gpu[2];
                    });
                } // mfi
            }  // buoyancy_type
          } // not buoyancy_type == 1
        } // has moisture
    } // anelastic?
}
