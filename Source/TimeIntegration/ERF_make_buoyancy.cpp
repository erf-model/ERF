#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>
#include <ERF_Constants.H>
#include <Advection.H>
#include <Diffusion.H>
#include <TimeIntegration.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>
#include <PlaneAverage.H>

using namespace amrex;

void make_buoyancy (BoxArray& grids_to_evolve,
                    Vector<MultiFab>& S_data,
                    const MultiFab& S_prim,
                          MultiFab& buoyancy,
#ifdef ERF_USE_MOISTURE
                    const MultiFab& qvapor,
                    const MultiFab& qcloud,
                    const MultiFab& qice,
#endif
                    const amrex::Geometry geom,
                    const SolverChoice& solverChoice,
                    const MultiFab* r0)
{
    BL_PROFILE_REGION("make_buoyancy()");

    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // ******************************************************************************************
    // Dry versions of buoyancy expressions (type 1 and type 2/3 -- types 2 and 3 are equivalent)
    // ******************************************************************************************
#ifndef ERF_USE_MOISTURE

    if (solverChoice.buoyancy_type == 1) {
        for ( MFIter mfi(buoyancy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& valid_bx = surroundingNodes(grids_to_evolve[mfi.index()],2);

            // Construct intersection of current tilebox and valid region for updating
            Box tbz = mfi.tilebox() & valid_bx;

            // We don't compute a source term for z-momentum on the bottom or top boundary
            tbz.growLo(2,-1);
            tbz.growHi(2,-1);

            const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
            const Array4<      Real> & buoyancy_fab = buoyancy.array(mfi);

            // Base state density
            const Array4<const Real>& r0_arr = r0->const_array(mfi);

            amrex::ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                buoyancy_fab(i, j, k) = grav_gpu[2] * 0.5 * ( cell_data(i,j,k) + cell_data(i,j,k-1)
                                                               - r0_arr(i,j,k) -    r0_arr(i,j,k-1) );
            });
        } // mfi

    } else if (solverChoice.buoyancy_type == 2 || solverChoice.buoyancy_type == 3) {

        PlaneAverage state_ave(&(S_data[IntVar::cons]), geom, 2);
        PlaneAverage prim_ave(&S_prim, geom, 2);

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
            const Box& valid_bx = surroundingNodes(grids_to_evolve[mfi.index()],2);

            // Construct intersection of current tilebox and valid region for updating
            Box tbz = mfi.tilebox() & valid_bx;

            // We don't compute a source term for z-momentum on the bottom or top boundary
            tbz.growLo(2,-1);
            tbz.growHi(2,-1);

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
    // Moist versions of buoyancy expressions -- only types 2 and 3 allowed
    // ******************************************************************************************
#ifdef ERF_USE_MOISTURE
    PlaneAverage state_ave(&(S_data[IntVar::cons]), geom, 2);
    PlaneAverage prim_ave(&S_prim, geom, 2);

    state_ave.compute_averages(ZDir(), state_ave.field());
    prim_ave.compute_averages(ZDir(), prim_ave.field());

    int ncell = state_ave.ncell_line();

    Gpu::HostVector<Real> rho_h(ncell), theta_h(ncell),
                          qp_h(ncell), qv_h(ncell),
                          qi_h(ncell), qc_h(ncell);


    state_ave.line_average(Rho_comp, rho_h);
    Gpu::DeviceVector<Real>   rho_d(ncell);
    Gpu::copyAsync(Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());

    prim_ave.line_average(PrimTheta_comp, theta_h);
    Gpu::DeviceVector<Real> theta_d(ncell);
    Gpu::copyAsync(Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());

    Real*   rho_d_ptr =   rho_d.data();
    Real* theta_d_ptr = theta_d.data();

    if (solverChoice.buoyancy_type == 2) {

        PlaneAverage qv_ave(&qvapor, geom, 2);
        PlaneAverage qc_ave(&qcloud, geom, 2);
        PlaneAverage qi_ave(&qice, geom, 2);

        qv_ave.compute_averages(ZDir(), qv_ave.field());
        qc_ave.compute_averages(ZDir(), qc_ave.field());
        qi_ave.compute_averages(ZDir(), qi_ave.field());

        prim_ave.line_average(PrimQp_comp, qp_h);
        qv_ave.line_average(0, qv_h);
        qi_ave.line_average(0, qi_h);
        qc_ave.line_average(0, qc_h);

        // copy data to device
        Gpu::DeviceVector<Real>    qp_d(ncell);
        Gpu::DeviceVector<Real>    qv_d(ncell);
        Gpu::DeviceVector<Real>    qc_d(ncell);
        Gpu::DeviceVector<Real>    qi_d(ncell);

        Gpu::copyAsync(Gpu::hostToDevice, qp_h.begin(), qp_h.end(), qp_d.begin());
        Gpu::copyAsync(Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
        Gpu::copyAsync(Gpu::hostToDevice, qi_h.begin(), qi_h.end(), qi_d.begin());
        Gpu::copyAsync(Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
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
            const Box& valid_bx = surroundingNodes(grids_to_evolve[mfi.index()],2);

            // Construct intersection of current tilebox and valid region for updating
            Box tbz = mfi.tilebox() & valid_bx;

            // We don't compute a source term for z-momentum on the bottom or top boundary
            tbz.growLo(2,-1);
            tbz.growHi(2,-1);

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
            const Box& valid_bx = surroundingNodes(grids_to_evolve[mfi.index()],2);

            // Construct intersection of current tilebox and valid region for updating
            Box tbz = mfi.tilebox() & valid_bx;

            // We don't compute a source term for z-momentum on the bottom or top boundary
            tbz.growLo(2,-1);
            tbz.growHi(2,-1);

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
#endif
}
