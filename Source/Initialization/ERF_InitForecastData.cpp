/**
 * \file ERF_InitForecastData.cpp
 */
#include <ERF.H>
#include "ERF_ReadCustomBinaryIC.H"
#include "ERF_Interpolation_Bilinear.H"

using namespace amrex;

void fill_weather_data_multifab(MultiFab& mf,
                                const Geometry& geom,
                                const int nx,
                                const int ny,
                                const int nz,
                                const Vector<Real>& zvec_h,
                                const Vector<Real>& rho_h,
                                const Vector<Real>& uvel_h,
                                const Vector<Real>& vvel_h,
                                const Vector<Real>& wvel_h,
                                const Vector<Real>& theta_h,
                                const Vector<Real>& qv_h,
                                const Vector<Real>& qc_h,
                                const Vector<Real>& qr_h)
{


    const int nx_d = nx;
    const int ny_d = ny;
    const int nz_d = nz;
    const int tot_size = nx*ny*nz;

    amrex::Gpu::DeviceVector<Real> zvec_d(nz), rho_d(tot_size), uvel_d(tot_size),
                                   vvel_d(tot_size), wvel_d(tot_size), theta_d(tot_size),
                                   qv_d(tot_size), qc_d(tot_size), qr_d(tot_size);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, uvel_h.begin(), uvel_h.end(), uvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, vvel_h.begin(), vvel_h.end(), vvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, wvel_h.begin(), wvel_h.end(), wvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qr_h.begin(), qr_h.end(), qr_d.begin());

    Real* zvec_d_ptr  = zvec_d.data();
    Real* rho_d_ptr   = rho_d.data();
    Real* uvel_d_ptr  = uvel_d.data();
    Real* vvel_d_ptr  = vvel_d.data();
    Real* wvel_d_ptr  = wvel_d.data();
    Real* theta_d_ptr = theta_d.data();
    Real* qv_d_ptr = qv_d.data();
    Real* qc_d_ptr = qc_d.data();
    Real* qr_d_ptr = qr_d.data();

    const int ncomp = mf.nComp();

    const auto prob_lo  = geom.ProbLoArray();
    const auto dx       = geom.CellSizeArray();

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.nodaltilebox();
        Array4<Real> const& arr = mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real x = prob_lo[0] + i * dx[0];
            const Real y = prob_lo[1] + j * dx[1];
            const Real z = prob_lo[2] + k * dx[2];
            int kloc = -1;
            for (int kk = 0; kk < nz_d; kk++) {
                if (zvec_d_ptr[kk] >= z) {
                    kloc = kk;
                    break;
                }
            }
            // Replace this with logic using your coarse input data
            int idx = get_single_index(i,j,kloc,nx_d,ny_d);
            arr(i,j,k,0) = rho_d_ptr[idx];
            arr(i,j,k,1) = uvel_d_ptr[idx];
            arr(i,j,k,2) = vvel_d_ptr[idx];
            arr(i,j,k,3) = wvel_d_ptr[idx];
            arr(i,j,k,4) = theta_d_ptr[idx];
            arr(i,j,k,5) = qv_d_ptr[idx];
            arr(i,j,k,6) = qc_d_ptr[idx];
            arr(i,j,k,7) = qr_d_ptr[idx];
        });
    }
}


void
ERF::init_coarse_weather_data()
{
    Vector<Real> xvec_h, yvec_h, zvec_h, rho_h, uvel_h, vvel_h, wvel_h, theta_h, qv_h, qc_h, qr_h;

    std::string filename;
    ParmParse pp("erf");
    pp.query("IC_file", filename);

    if (filename.empty()) {
        amrex::Abort("Error: IC_file is not specified in the input file.");
    }

    ReadCustomBinaryIC(filename, xvec_h, yvec_h, zvec_h, rho_h,
                                uvel_h, vvel_h, wvel_h,
                                theta_h, qv_h, qc_h, qr_h);

    // Number of cells
    int nx_cells = xvec_h.size()-1;
    int ny_cells = yvec_h.size()-1;
    int nz_cells = zvec_h.size()-1;

    IntVect dom_lo(0, 0, 0);
    IntVect dom_hi(nx_cells-1, ny_cells-1, nz_cells-1);
    Box domain(dom_lo, dom_hi);

    // Define the extents of the physical domain box
    RealBox real_box({xvec_h[0], yvec_h[0], zvec_h[0]}, {xvec_h[nx_cells], yvec_h[ny_cells], zvec_h[nz_cells]});

    int coord = 0; // Cartesian
    Array<int, AMREX_SPACEDIM> is_periodic{0, 0, 0}; // non-periodic

    Geometry geom(domain, real_box, coord, is_periodic);

    BoxArray ba(domain);
    ba.maxSize(64);
    BoxArray nba = amrex::convert(ba, IntVect::TheNodeVector()); // nodal in all directions

    // Create DistributionMapping
    DistributionMapping dm(nba);

    int ncomp = 8;
    int ngrow = 0;

    int n_time = 1;      // or however many time slices you want
    weather_forecast_data.resize(n_time);
    MultiFab& weather_mf = weather_forecast_data[0];
    weather_mf.define(nba, dm, ncomp, ngrow);

    fill_weather_data_multifab(weather_mf, geom, nx_cells+1, ny_cells+1, nz_cells+1,
                               zvec_h, rho_h,uvel_h, vvel_h, wvel_h,
                               theta_h, qv_h, qc_h, qr_h);

    Vector<std::string> varnames = {
    "rho", "uvel", "vvel", "wvel", "theta", "qv", "qc", "qr"
    }; // Customize variable names

    const std::string plotfilename = "plt_weather"; // or any name you want

    const Real time = 0.0;

    // Assume weather_mf is nodal in all directions
    BoxArray cba = amrex::convert(weather_mf.boxArray(), IntVect::TheCellVector());

    MultiFab cell_centered_mf(cba, weather_mf.DistributionMap(),
                          weather_mf.nComp(), 0);

    amrex::average_node_to_cellcenter(cell_centered_mf, 0, weather_mf, 0, weather_mf.nComp());


    WriteSingleLevelPlotfile(
        plotfilename,
        cell_centered_mf,
        varnames,
        geom,
        time,
        0 // level
    );

}
