/**
 * \file ERF_ReadCustomBinaryForecastData.cpp
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
                                const Vector<Real>& xvec_h,
                                const Vector<Real>& yvec_h,
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

    const Real xmin = xvec_h[0];
    const Real ymin = yvec_h[0];
    const Real zmin = zvec_h[0];

    amrex::Gpu::DeviceVector<Real> rho_d(nx*ny*nz), uvel_d(nx*ny*nz), vvel_d(nx*ny*nz), wvel_d(nx*ny*nz),
                                   theta_d(nx*ny*nz), qv_d(nx*ny*nz), qc_d(nx*ny*nz), qr_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, uvel_h.begin(), uvel_h.end(), uvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, vvel_h.begin(), vvel_h.end(), vvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, wvel_h.begin(), wvel_h.end(), wvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qr_h.begin(), qr_h.end(), qr_d.begin());

    Real* rho_d_ptr   = rho_d.data();
    Real* uvel_d_ptr  = uvel_d.data();
    Real* vvel_d_ptr  = vvel_d.data();
    Real* wvel_d_ptr  = wvel_d.data();
    Real* theta_d_ptr = theta_d.data();
    Real* qv_d_ptr = qv_d.data();
    Real* qc_d_ptr = qc_d.data();
    Real* qr_d_ptr = qr_d.data();

    const int ncomp = mf.nComp();

    const auto prob_lo  = geom.ProbLo();
    const auto dx       = geom.CellSize();

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

    const Box& bx = mfi.nodaltilebox();
    Array4<Real> const& arr = mf.array(mfi);

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real x        = prob_lo[0] + i * dx[0];
        const Real y        = prob_lo[1] + j * dx[1];
        const Real z        = prob_lo[2] + k * dx[2];

        int iloc=-1, jloc=-1, kloc=-1;
        iloc = std::floor((x-xmin)/dx[0]);
        jloc = std::floor((y-ymin)/dx[1]);
        kloc = std::floor((z-zmin)/dx[2]);

        int idx = get_single_index(iloc,jloc,kloc,nx_d,ny_d);
        if(n==0) {
            arr(i,j,k,n) = rho_d_ptr[idx];
        } else if (n==1) {
            arr(i,j,k,n) = uvel_d_ptr[idx];
        } else if (n==2) {
            arr(i,j,k,n) = vvel_d_ptr[idx];
        } else if (n==3) {
            arr(i,j,k,n) = wvel_d_ptr[idx];
        } else if (n==4) {
            arr(i,j,k,n) = theta_d_ptr[idx];
        } else if (n==5) {
            arr(i,j,k,n) = qv_d_ptr[idx];
        } else if (n==6) {
            arr(i,j,k,n) = qc_d_ptr[idx];
        } else if (n==7) {
            arr(i,j,k,n) = qr_d_ptr[idx];
        }
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
    int nx = xvec_h.size()-1;
    int ny = yvec_h.size()-1;
    int nz = zvec_h.size()-1;

    IntVect dom_lo(0, 0, 0);
    IntVect dom_hi(nx-1, ny-1, nz-1); // 64x64x32
    Box domain(dom_lo, dom_hi);

    // Define the extents of the physical domain box
    RealBox real_box({xvec_h[0], yvec_h[0], zvec_h[0]}, {xvec_h[nx], yvec_h[ny], zvec_h[nz]});

    int coord = 0; // Cartesian
    Array<int, AMREX_SPACEDIM> is_periodic{0, 0, 0}; // non-periodic

    Geometry geom(domain, real_box, coord, is_periodic);

    // -------------------------
    // Make BoxArray and MultiFab
    // -------------------------

    BoxArray ba(domain);
    ba.maxSize(64);
    BoxArray nba = amrex::convert(ba, IntVect::TheNodeVector()); // nodal in all directions

    // Create DistributionMapping
    DistributionMapping dm(nba);

    int ncomp = 8;
    int ngrow = 0;

    MultiFab weather_mf(nba, dm, ncomp, ngrow);

    // -------------------------
    // Fill the MultiFab
    // -------------------------
    fill_weather_data_multifab(weather_mf, geom, nx+1, ny+1, nz+1,
                               xvec_h, yvec_h, zvec_h,
                               rho_h,uvel_h, vvel_h, wvel_h,
                               theta_h, qv_h, qc_h, qr_h);

    /*Vector<std::string> varnames = {
    "var0", "var1", "var2", "var3", "var4", "var5", "var6", "var7"
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
    );*/

}
