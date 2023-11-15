#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>
#include <ERF_Constants.H>
#include <Advection.H>
#include <Diffusion.H>
#include <TI_headers.H>
#include <TileNoZ.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>
#include <PlaneAverage.H>

using namespace amrex;

/**
 * Function for computing the buoyancy term to be used in the evolution
 * equation for the z-component of momentum in the slow integrator.  There
 * are three options for how buoyancy is computed (two are the same in the absence of moisture).
 *
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[out] buoyancy the buoyancy term computed here
 * @param[in]  qvapor water vapor
 * @param[in]  qcloud cloud water
 * @param[in]  qice   cloud ice
 * @param[in]  qv_d   lateral average of cloud vapor
 * @param[in]  qc_d   lateral average of cloud vapor
 * @param[in]  qd_d   lateral average of cloud vapor
 * @param[in]  geom   Container for geometric informaiton
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in]  r0     Reference (hydrostatically stratified) density
 */

void make_buoyancy (Vector<MultiFab>& S_data,
                    const MultiFab& S_prim,
                          MultiFab& buoyancy,
#if defined(ERF_USE_MOISTURE)
                    const MultiFab& qmoist,
                    Gpu::DeviceVector<Real> qv_d,
                    Gpu::DeviceVector<Real> qc_d,
                    Gpu::DeviceVector<Real> qi_d,
#endif
                    const amrex::Geometry geom,
                    const SolverChoice& solverChoice,
                    const MultiFab* r0)
{
    BL_PROFILE_REGION("make_buoyancy()");

    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const int klo = 0;
    const int khi = geom.Domain().bigEnd()[2] + 1;

#if defined(ERF_USE_MOISTURE)
    MultiFab qvapor(qmoist, make_alias, 0, 1);
    MultiFab qcloud(qmoist, make_alias, 1, 1);
    MultiFab qice  (qmoist, make_alias, 2, 1);
#endif

    // ******************************************************************************************
    // Dry versions of buoyancy expressions (type 1 and type 2/3 -- types 2 and 3 are equivalent)
    // ******************************************************************************************
#if !defined(ERF_USE_MOISTURE) && !defined(ERF_USE_WARM_NO_PRECIP)

    if (solverChoice.buoyancy_type == 1) {
        for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbz = mfi.tilebox();

            // We don't compute a source term for z-momentum on the bottom or top domain boundary
            if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
            if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

            // Base state density
            const Array4<const Real>& r0_arr = r0->const_array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                buoyancy_fab(i, j, k) = grav_gpu[2] * 0.5 * ( (cell_data(i,j,k  ) - r0_arr(i,j,k  ))
                                                             +(cell_data(i,j,k-1) - r0_arr(i,j,k-1)) );
            });
        } // mfi

    } else if (solverChoice.buoyancy_type == 2 || solverChoice.buoyancy_type == 3) {

        PlaneAverage state_ave(&(S_data[IntVar::cons]), geom, solverChoice.ave_plane);
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

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ]);
                Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1]);

                Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));
                Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));

                Real qplus  = (tempp3d-tempp1d)/tempp1d;
                Real qminus = (tempm3d-tempm1d)/tempm1d;

                Real qavg  = Real(0.5) * (qplus + qminus);
                Real r0avg = Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

                buoyancy_fab(i, j, k) = -qavg * r0avg * grav_gpu[2];
            });
        } // mfi
    } // buoyancy_type
#endif

    // ******************************************************************************************
    // Moist versions of buoyancy expressions
    // ******************************************************************************************
#if defined(ERF_USE_MOISTURE)
    if (solverChoice.buoyancy_type == 1) {

        for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbz = mfi.tilebox();

            // We don't compute a source term for z-momentum on the bottom or top domain boundary
            if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
            if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

            // Base state density
            const Array4<const Real>& r0_arr = r0->const_array(mfi);

            const Array4<const Real> & qv_data    = qvapor.array(mfi);
            const Array4<const Real> & qc_data    = qcloud.array(mfi);
            const Array4<const Real> & qi_data    = qice.array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rhop_hi = cell_data(i,j,k  ,Rho_comp) * (1.0 + qv_data(i,j,k  ) + qc_data(i,j,k  )
                                                                  + qi_data(i,j,k  ) + cell_data(i,j,k,RhoQp_comp)/cell_data(i,j,k,Rho_comp)) - r0_arr(i,j,k  );
                Real rhop_lo = cell_data(i,j,k-1,Rho_comp) * (1.0 + qv_data(i,j,k-1) + qc_data(i,j,k-1)
                                                                  + qi_data(i,j,k-1) + + cell_data(i,j,k-1,RhoQp_comp)/cell_data(i,j,k-1,Rho_comp)) - r0_arr(i,j,k-1);
                buoyancy_fab(i, j, k) = grav_gpu[2] * 0.5 * ( rhop_hi + rhop_lo );
            });
        } // mfi

    } else {

    PlaneAverage state_ave(&(S_data[IntVar::cons]), geom, solverChoice.ave_plane);
    PlaneAverage  prim_ave(&S_prim                , geom, solverChoice.ave_plane);

    // Compute horizontal averages of all components of each field
    state_ave.compute_averages(ZDir(), state_ave.field());
     prim_ave.compute_averages(ZDir(), prim_ave.field());

    int ncell = state_ave.ncell_line();

    Gpu::HostVector  <Real> rho_h(ncell), theta_h(ncell), qp_h(ncell);
    Gpu::DeviceVector<Real> rho_d(ncell), theta_d(ncell);

    state_ave.line_average(Rho_comp, rho_h);
    Gpu::copyAsync(Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());

    prim_ave.line_average(PrimTheta_comp, theta_h);
    Gpu::copyAsync(Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());

    Real*   rho_d_ptr =   rho_d.data();
    Real* theta_d_ptr = theta_d.data();

    if (solverChoice.buoyancy_type == 2) {

        prim_ave.line_average(PrimQp_comp, qp_h);

        Gpu::DeviceVector<Real>    qp_d(ncell);

        // Copy data to device
        Gpu::copyAsync(Gpu::hostToDevice, qp_h.begin(), qp_h.end(), qp_d.begin());
        Gpu::streamSynchronize();

        Real*    qp_d_ptr =    qp_d.data();
        Real*    qv_d_ptr =    qv_d.data();
        Real*    qc_d_ptr =    qc_d.data();
        Real*    qi_d_ptr =    qi_d.data();

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

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<const Real> & cell_prim  = S_prim.array(mfi);
            const Array4<const Real> & qv_data    = qvapor.array(mfi);
            const Array4<const Real> & qc_data    = qcloud.array(mfi);
            const Array4<const Real> & qi_data    = qice.array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ]);
                Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1]);

                Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));
                Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));

                Real qplus = 0.61* ( qv_data(i,j,k)-qv_d_ptr[k]) -
                                    (qc_data(i,j,k)-qc_d_ptr[k]+
                                     qi_data(i,j,k)-qi_d_ptr[k]+
                                     cell_prim(i,j,k,PrimQp_comp)-qp_d_ptr[k])
                           + (tempp3d-tempp1d)/tempp1d*(Real(1.0) + Real(0.61)*qv_d_ptr[k]-qc_d_ptr[k]-qi_d_ptr[k]-qp_d_ptr[k]);

                Real qminus = 0.61 *( qv_data(i,j,k-1)-qv_d_ptr[k-1]) -
                                     (qc_data(i,j,k-1)-qc_d_ptr[k-1]+
                                      qi_data(i,j,k-1)-qi_d_ptr[k-1]+
                                      cell_prim(i,j,k-1,PrimQp_comp)-qp_d_ptr[k-1])
                           + (tempm3d-tempm1d)/tempm1d*(Real(1.0) + Real(0.61)*qv_d_ptr[k-1]-qi_d_ptr[k-1]-qc_d_ptr[k-1]-qp_d_ptr[k-1]);

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

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<const Real> & cell_prim  = S_prim.array(mfi);
            const Array4<const Real> & qv_data    = qvapor.array(mfi);
            const Array4<const Real> & qc_data    = qcloud.array(mfi);
            const Array4<const Real> & qi_data    = qice.array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real tempp1d = getTgivenRandRTh(rho_d_ptr[k  ], rho_d_ptr[k  ]*theta_d_ptr[k  ]);
                Real tempm1d = getTgivenRandRTh(rho_d_ptr[k-1], rho_d_ptr[k-1]*theta_d_ptr[k-1]);

                Real tempp3d  = getTgivenRandRTh(cell_data(i,j,k  ,Rho_comp), cell_data(i,j,k  ,RhoTheta_comp));
                Real tempm3d  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp), cell_data(i,j,k-1,RhoTheta_comp));

                Real qplus = 0.61 * qv_data(i,j,k) - (qc_data(i,j,k)+ qi_data(i,j,k)+ cell_prim(i,j,k,PrimQp_comp))
                           + (tempp3d-tempp1d)/tempp1d;

                Real qminus = 0.61 *qv_data(i,j,k-1) - (qc_data(i,j,k-1)+ qi_data(i,j,k-1)+ cell_prim(i,j,k-1,PrimQp_comp))
                           + (tempm3d-tempm1d)/tempm1d;

                Real qavg  = Real(0.5) * (qplus + qminus);
                Real r0avg = Real(0.5) * (rho_d_ptr[k] + rho_d_ptr[k-1]);

                buoyancy_fab(i, j, k) = -qavg * r0avg * grav_gpu[2];
            });
        } // mfi
    } // buoyancy_type
    } // not buoyancy_type == 1
#endif

#if defined(ERF_USE_WARM_NO_PRECIP)

    AMREX_ALWAYS_ASSERT(solverChoice.buoyancy_type == 1);

        for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbz = mfi.tilebox();

            // We don't compute a source term for z-momentum on the bottom or top domain boundary
            if (tbz.smallEnd(2) == klo) tbz.growLo(2,-1);
            if (tbz.bigEnd(2)   == khi) tbz.growHi(2,-1);

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

            // Base state density
            const Array4<const Real>& r0_arr = r0->const_array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rhop_hi = cell_data(i,j,k  ,Rho_comp)   + cell_data(i,j,k  ,RhoQv_comp)
                             + cell_data(i,j,k  ,RhoQc_comp) - r0_arr(i,j,k  );
                Real rhop_lo = cell_data(i,j,k-1,Rho_comp)   + cell_data(i,j,k-1,RhoQv_comp)
                             + cell_data(i,j,k-1,RhoQc_comp) - r0_arr(i,j,k-1);
                buoyancy_fab(i, j, k) = grav_gpu[2] * 0.5 * ( rhop_hi + rhop_lo );
            });
        } // mfi
#endif
}
