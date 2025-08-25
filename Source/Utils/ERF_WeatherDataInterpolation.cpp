#ifndef ERF_WEATHERDATAINTERPOLATION_H_
#define ERF_WEATHERDATAINTERPOLATION_H_
/**
 * Trilinear interpolation of weather forecast data onto the simulation mesh
 * The coarse weather forecast data is interpolated in time first to get the forecast
 * at the current time, and then spatially interpolated onto the simulation mesh
 */

#include <filesystem>
#include <stdexcept>
#include "ERF.H"
#include "ERF_ReadCustomBinaryIC.H"
#include "ERF_Interpolation_Bilinear.H"

using namespace amrex;
namespace fs = std::filesystem;

void fill_weather_data_multifab(MultiFab& mf,
     const Geometry& geom_weather,
     const int nx,
     const int ny,
     const int nz,
     const Vector<Real>& latvec_h,
     const Vector<Real>& lonvec_h,
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

    amrex::Gpu::DeviceVector<Real> latvec_d(nx*ny), lonvec_d(nx*ny), zvec_d(nz);
    amrex::Gpu::DeviceVector<Real> rho_d(tot_size), uvel_d(tot_size),
                                   vvel_d(tot_size), wvel_d(tot_size), theta_d(tot_size),
                                   qv_d(tot_size), qc_d(tot_size), qr_d(tot_size);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, latvec_h.begin(), latvec_h.end(), latvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, lonvec_h.begin(), lonvec_h.end(), lonvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, uvel_h.begin(), uvel_h.end(), uvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, vvel_h.begin(), vvel_h.end(), vvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, wvel_h.begin(), wvel_h.end(), wvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qr_h.begin(), qr_h.end(), qr_d.begin());

    Real* latvec_d_ptr  = latvec_d.data();
    Real* lonvec_d_ptr  = lonvec_d.data();
    Real* zvec_d_ptr  = zvec_d.data();
    Real* rho_d_ptr   = rho_d.data();
    Real* uvel_d_ptr  = uvel_d.data();
    Real* vvel_d_ptr  = vvel_d.data();
    Real* wvel_d_ptr  = wvel_d.data();
    Real* theta_d_ptr = theta_d.data();
    Real* qv_d_ptr = qv_d.data();
    Real* qc_d_ptr = qc_d.data();
    Real* qr_d_ptr = qr_d.data();

    const auto prob_lo  = geom_weather.ProbLoArray();
    const auto dx       = geom_weather.CellSizeArray();

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.nodaltilebox();
        Array4<Real> const& arr = mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real z = prob_lo[2] + k * dx[2];
            int kloc = -1;
            for (int kk = 0; kk < nz_d; kk++) {
                 if (zvec_d_ptr[kk] > z) {
                     kloc = kk;
                     break;
                 }
            }
            // Replace this with logic using your coarse input data
            int idx1 = get_single_index(i,j,kloc-1,nx_d,ny_d);
            int idx2 = get_single_index(i,j,kloc,nx_d,ny_d);
            Real dz = zvec_d_ptr[kloc] - zvec_d_ptr[kloc-1];
            Real fac = (z-zvec_d_ptr[kloc-1])/dz;
            arr(i,j,k,0) = rho_d_ptr[idx1]   + fac*(rho_d_ptr[idx2]-rho_d_ptr[idx1]);
            arr(i,j,k,1) = uvel_d_ptr[idx1]  + fac*(uvel_d_ptr[idx2]-uvel_d_ptr[idx1]);
            arr(i,j,k,2) = vvel_d_ptr[idx1]  + fac*(vvel_d_ptr[idx2]-vvel_d_ptr[idx1]);
            arr(i,j,k,3) = wvel_d_ptr[idx1]  + fac*(wvel_d_ptr[idx2]-wvel_d_ptr[idx1]);
            arr(i,j,k,4) = theta_d_ptr[idx1] + fac*(theta_d_ptr[idx2]-theta_d_ptr[idx1]);
            arr(i,j,k,5) = qv_d_ptr[idx1]    + fac*(qv_d_ptr[idx2]-qv_d_ptr[idx1]);
            arr(i,j,k,6) = qc_d_ptr[idx1]    + fac*(qc_d_ptr[idx2]-qc_d_ptr[idx1]);
            arr(i,j,k,7) = qr_d_ptr[idx1]    + fac*(qr_d_ptr[idx2]-qr_d_ptr[idx1]);
            idx1 = get_single_index(i,j,0,nx_d,ny_d);
            arr(i,j,k,8) = latvec_d_ptr[idx1];
            arr(i,j,k,9) = lonvec_d_ptr[idx1];
        });
    }
}

enum class BoundType { Lo, Hi };
enum class MultiFabType { CC, NC };


void PlotMultiFab(const MultiFab& mf,
                  const Geometry& geom_mf,
                  const std::string plotfilename,
                  MultiFabType mftype)
{

    Vector<std::string> varnames = {
    "rho", "uvel", "vvel", "wvel", "theta", "qv", "qc", "qr", "latitude", "longitude"
    }; // Customize variable names

    const Real time = 0.0;


    // Assume weather_mf is nodal in all directions
    if(mftype == MultiFabType::NC) {
        BoxArray cba = mf.boxArray();
        cba = amrex::convert(mf.boxArray(), IntVect::TheCellVector());

        MultiFab cc_mf(cba, mf.DistributionMap(),
               mf.nComp(), 0);

        amrex::average_node_to_cellcenter(cc_mf, 0, mf, 0, mf.nComp());

        WriteSingleLevelPlotfile(
            plotfilename,
            cc_mf,
            varnames,
            geom_mf,
            time,
            0 // level
        );
    } else {
        WriteSingleLevelPlotfile(
            plotfilename,
            mf,
            varnames,
            geom_mf,
            time,
            0 // level
        );
    }
}

void
ERF::CreateWeatherDataGeomBoxArrayDistMap(const std::string& filename,
                                          Geometry& geom_weather,
                                          BoxArray& nba,
                                          DistributionMapping& dm)
{
    Vector<Real> latvec_h, lonvec_h, xvec_h, yvec_h, zvec_h;
    Vector<Real> rho_h, uvel_h, vvel_h, wvel_h, theta_h, qv_h, qc_h, qr_h;

    ReadCustomBinaryIC(filename, latvec_h, lonvec_h,
                       xvec_h, yvec_h, zvec_h, rho_h,
                       uvel_h, vvel_h, wvel_h,
                       theta_h, qv_h, qc_h, qr_h);

    const auto prob_lo_erf  = geom[0].ProbLoArray();
    const auto prob_hi_erf  = geom[0].ProbHiArray();
    const auto dx_erf       = geom[0].CellSizeArray();

    if(prob_lo_erf[0] < xvec_h.front() + 4*dx_erf[0]){
        amrex::Abort("The xlo value of the domain has to be greater than " + std::to_string(xvec_h.front() + 4*dx_erf[0]));
    }
    if(prob_hi_erf[0] > xvec_h.back() - 4*dx_erf[0]){
        amrex::Abort("The xhi value of the domain has to be less than " + std::to_string(xvec_h.back() - 4*dx_erf[0]));
    }
    if(prob_lo_erf[1] < yvec_h.front() + 4*dx_erf[1]){
        amrex::Abort("The ylo value of the domain has to be greater than " + std::to_string(yvec_h.front() + 4*dx_erf[1]));
    }
    if(prob_hi_erf[1] > yvec_h.back() - 4*dx_erf[1]){
        amrex::Abort("The yhi value of the domain has to be less than " + std::to_string(yvec_h.back() - 4*dx_erf[1]));
    }

    // Number of cells
    int nx_cells = xvec_h.size()-1;
    int ny_cells = yvec_h.size()-1;

    const amrex::Geometry& geom0 = geom[0]; // or whatever your Geometry vector is called
    const amrex::Box& domainBox = geom0.Domain();
    const amrex::IntVect& domainSize = domainBox.size(); // Number of cells in each direction
    int nz_cells = domainSize[2];

    IntVect dom_lo(0, 0, 0);
    IntVect dom_hi(nx_cells-1, ny_cells-1, nz_cells-1);
    Box domain(dom_lo, dom_hi);

    const amrex::Real* prob_hi = geom0.ProbHi();

    // Define the extents of the physical domain box
    RealBox real_box({xvec_h[0], yvec_h[0], zvec_h[0]}, {xvec_h[nx_cells], yvec_h[ny_cells], prob_hi[2]});

    int coord = 0; // Cartesian
    Array<int, AMREX_SPACEDIM> is_periodic{0, 0, 0}; // non-periodic

    geom_weather.define(domain, real_box, coord, is_periodic);

    BoxArray ba(domain);
    ba.maxSize(64);
    nba = amrex::convert(ba, IntVect::TheNodeVector()); // nodal in all directions

    // Create DistributionMapping
    dm = DistributionMapping(nba);
 }

IntVect
find_bound_idx(const Real& x, const Real& y, const Real& z,
               const BoxList& bl_weather, const Geometry& geom_weather,
               BoundType bound_type)
{
    const auto prob_lo_weather  = geom_weather.ProbLoArray();
    const auto dx_weather       = geom_weather.CellSizeArray();

    int i, j, k;

    if (bound_type == BoundType::Lo) {
        i = static_cast<int>(std::floor((x - prob_lo_weather[0]) / dx_weather[0]));
        j = static_cast<int>(std::floor((y - prob_lo_weather[1]) / dx_weather[1]));
        k = static_cast<int>(std::floor((z - prob_lo_weather[2]) / dx_weather[2]));
    } else { // BoundType::Hi
        i = static_cast<int>(std::ceil((x - prob_lo_weather[0]) / dx_weather[0]));
        j = static_cast<int>(std::ceil((y - prob_lo_weather[1]) / dx_weather[1]));
        k = static_cast<int>(std::ceil((z - prob_lo_weather[2]) / dx_weather[2]));
    }

    IntVect idx(i, j, k);

    for (const auto& b : bl_weather) {
        if (b.contains(idx)) {
            return idx;
        }
    }

    amrex::Abort("Bound index not found in any box in BoxList!");
    return IntVect::TheZeroVector(); // unreachable if Abort
}

void
ERF::InterpWeatherDataOntoMesh (const Geometry& geom_weather,
                                    MultiFab& weather_forecast_interp)
{

    ParmParse pp_erf("erf");
    bool is_lateral_sponges_hurricanes = false;
    if (pp_erf.query("is_lateral_sponges_hurricanes", is_lateral_sponges_hurricanes)) {
        initial_state.resize(max_level+1);
        for (int lev = 0; lev < max_level+1; ++lev) {
            initial_state[lev].resize(vars_new[lev].size()+1);
            for (int comp = 0; comp < vars_new[lev].size(); ++comp) {
                const MultiFab& src = vars_new[lev][comp];
                initial_state[lev][comp].define(src.boxArray(), src.DistributionMap(),
                                        src.nComp(), src.nGrow());
            }
            int comp = vars_new[lev].size();
            const MultiFab& src = vars_new[lev][0];
            initial_state[lev][comp].define(src.boxArray(), src.DistributionMap(),
                                        2, src.nGrow());
        }
    }

    MultiFab& weather_mf    = weather_forecast_interp;
    MultiFab& erf_mf_cons   = initial_state[0][Vars::cons];
    MultiFab& erf_mf_xvel   = initial_state[0][Vars::xvel];
    MultiFab& erf_mf_yvel   = initial_state[0][Vars::yvel];
    MultiFab& erf_mf_zvel   = initial_state[0][Vars::zvel];
    MultiFab& erf_mf_latlon = initial_state[0][4];

    erf_mf_cons.setVal(0.0);
    erf_mf_xvel.setVal(0.0);
    erf_mf_yvel.setVal(0.0);
    erf_mf_zvel.setVal(0.0);
    erf_mf_latlon.setVal(0.0);

    BoxList bl_erf     = erf_mf_cons.boxArray().boxList();
    BoxList bl_weather = weather_mf.boxArray().boxList();

    const auto prob_lo_erf  = geom[0].ProbLoArray();
    const auto dx_erf       = geom[0].CellSizeArray();

    for (auto& b : bl_erf) {
        // You look at the lo corner of b, and find out the lowest cell in
        // coarse weather data you need for the interpolation. That gives
        // you the lo corner of the new b. Similarly, you can find out the
        // hi corner of the new b. For cells outside the coarse_weath_data's
        // bounding data, it's up to you. You probably want to use a biased
        // interpolation stencil.

        // Get the cell indices of the bottom corner and top corner
        const IntVect& lo_erf = b.smallEnd();  // Lower corner (inclusive)
        const IntVect& hi_erf = b.bigEnd();    // Upper corner (inclusive)

        Real x = prob_lo_erf[0] + lo_erf[0] * dx_erf[0];
        Real y = prob_lo_erf[1] + lo_erf[1] * dx_erf[1];
        Real z = prob_lo_erf[2] + lo_erf[2] * dx_erf[2];

        auto idx_lo = find_bound_idx(x, y, z, bl_weather, geom_weather, BoundType::Lo);

        x = prob_lo_erf[0] + (hi_erf[0]+1) * dx_erf[0];
        y = prob_lo_erf[1] + (hi_erf[1]+1) * dx_erf[1];
        z = prob_lo_erf[2] + (hi_erf[2]+1) * dx_erf[2];

        auto idx_hi = find_bound_idx(x, y, z, bl_weather, geom_weather, BoundType::Hi);

        b.setSmall(idx_lo);
        b.setBig(idx_hi);
    }

    BoxArray cba(std::move(bl_erf));
    cba.convert(IndexType::TheNodeType());  // <-- Make it nodal in all directions
    MultiFab tmp_coarse_data(cba, erf_mf_cons.DistributionMap(), weather_mf.nComp(), 0);
    tmp_coarse_data.ParallelCopy(weather_mf);

    PlotMultiFab(weather_mf, geom_weather, "plt_coarse_weather_par_copy",MultiFabType::NC);

    const auto prob_lo_weather  = geom_weather.ProbLoArray();
    const auto dx_weather       = geom_weather.CellSizeArray();

    for (MFIter mfi(erf_mf_cons); mfi.isValid(); ++mfi) {
        const Array4<Real> &fine_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &fine_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &fine_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &fine_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &fine_latlon_arr = erf_mf_latlon.array(mfi);

        const Array4<Real> &crse_arr = tmp_coarse_data.array(mfi);

        const Box& gbx = mfi.growntilebox(); // tilebox + ghost cells

        const Box &gtbx = mfi.tilebox(IntVect(1,0,0));
        const Box &gtby = mfi.tilebox(IntVect(0,1,0));
        const Box &gtbz = mfi.tilebox(IntVect(0,0,1));

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+0.5) * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+0.5) * dx_erf[1];
            Real z = prob_lo_erf[2] + (k+0.5) * dx_erf[2];

            Real rho    = interpolate_from_coarse(crse_arr, 0, x, y, z, prob_lo_weather.data(), dx_weather.data());
            Real theta  = interpolate_from_coarse(crse_arr, 4, x, y, z, prob_lo_weather.data(), dx_weather.data());
            Real qv     = interpolate_from_coarse(crse_arr, 5, x, y, z, prob_lo_weather.data(), dx_weather.data());
            Real qc     = interpolate_from_coarse(crse_arr, 6, x, y, z, prob_lo_weather.data(), dx_weather.data());
            Real qr     = interpolate_from_coarse(crse_arr, 7, x, y, z, prob_lo_weather.data(), dx_weather.data());
            Real lat    = interpolate_from_coarse(crse_arr, 8, x, y, z, prob_lo_weather.data(), dx_weather.data());
            Real lon    = interpolate_from_coarse(crse_arr, 9, x, y, z, prob_lo_weather.data(), dx_weather.data());

            fine_cons_arr(i,j,k,Rho_comp) = rho;
            fine_cons_arr(i,j,k,RhoTheta_comp) = rho*theta;
            fine_cons_arr(i,j,k,RhoQ1_comp) = rho*qv;
            fine_cons_arr(i,j,k,RhoQ2_comp) = rho*qc;
            fine_cons_arr(i,j,k,RhoQ3_comp) = rho*qr;

            fine_latlon_arr(i,j,k,0) = lat;
            fine_latlon_arr(i,j,k,1) = lon;
        });

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + i * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+0.5) * dx_erf[1];
            Real z = prob_lo_erf[2] + (k+0.5) * dx_erf[2];
            fine_xvel_arr(i, j, k, 0) = interpolate_from_coarse(crse_arr, 1, x, y, z, prob_lo_weather.data(), dx_weather.data());
        });

        ParallelFor(gtby, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+0.5) * dx_erf[0];
            Real y = prob_lo_erf[1] + j       * dx_erf[1];
            Real z = prob_lo_erf[2] + (k+0.5) * dx_erf[2];
            fine_yvel_arr(i, j, k, 0) = interpolate_from_coarse(crse_arr, 2, x, y, z, prob_lo_weather.data(), dx_weather.data());
        });

        ParallelFor(gtbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+0.5) * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+0.5) * dx_erf[1];
            Real z = prob_lo_erf[2] + k       * dx_erf[2];
            fine_zvel_arr(i, j, k, 0) = interpolate_from_coarse(crse_arr, 3, x, y, z, prob_lo_weather.data(), dx_weather.data());
        });
    }

    Vector<std::string> varnames = {
    "rho", "uvel", "vvel", "wvel", "theta", "qv", "qc", "qr"
    }; // Customize variable names

     Vector<std::string> varnames_cons = {
    "rho", "rhotheta", "ke", "sc", "rhoqv", "rhoqc", "rhoqr"
    }; // Customize variable names

    Vector<std::string> varnames_plot_mf = {
    "rho", "rhotheta", "rhoqv", "rhoqc", "rhoqr", "xvel", "yvel", "zvel", "latitude", "longitude"
    }; // Customize variable names


    const Real time = 0.0;

    std::string pltname = "plt_interp";

    MultiFab plot_mf(erf_mf_cons.boxArray(), erf_mf_cons.DistributionMap(),
                     10, 0);

    plot_mf.setVal(0.0);

    for (MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
        const Array4<Real> &plot_mf_arr = plot_mf.array(mfi);
        const Array4<Real> &erf_mf_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &erf_mf_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &erf_mf_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &erf_mf_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &erf_mf_latlon_arr = erf_mf_latlon.array(mfi);

        const Box& bx = mfi.validbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plot_mf_arr(i,j,k,0) = erf_mf_cons_arr(i,j,k,Rho_comp);
            plot_mf_arr(i,j,k,1) = erf_mf_cons_arr(i,j,k,RhoTheta_comp);
            plot_mf_arr(i,j,k,2) = erf_mf_cons_arr(i,j,k,RhoQ1_comp);
            plot_mf_arr(i,j,k,3) = erf_mf_cons_arr(i,j,k,RhoQ2_comp);
            plot_mf_arr(i,j,k,4) = erf_mf_cons_arr(i,j,k,RhoQ3_comp);

            plot_mf_arr(i,j,k,5) = (erf_mf_xvel_arr(i,j,k,0) + erf_mf_xvel_arr(i+1,j,k,0))/2.0;
            plot_mf_arr(i,j,k,6) = (erf_mf_yvel_arr(i,j,k,0) + erf_mf_yvel_arr(i,j+1,k,0))/2.0;
            plot_mf_arr(i,j,k,7) = (erf_mf_zvel_arr(i,j,k,0) + erf_mf_zvel_arr(i,j,k+1,0))/2.0;

            plot_mf_arr(i,j,k,8) = erf_mf_latlon_arr(i,j,k,0);
            plot_mf_arr(i,j,k,9) = erf_mf_latlon_arr(i,j,k,1);
        });
    }


    WriteSingleLevelPlotfile(
            pltname,
            plot_mf,
            varnames_plot_mf,
            geom[0],
            time,
            0 // level
        );
}

void
ERF::FillWeatherDataMultiFab(const std::string& filename,
                             const Geometry& geom_weather,
                             const BoxArray& nba,
                             const DistributionMapping& dm,
                             Vector<MultiFab>& weather_forecast_data)
{

    Vector<Real> latvec_h, lonvec_h, xvec_h, yvec_h, zvec_h;
    Vector<Real> rho_h, uvel_h, vvel_h, wvel_h, theta_h, qv_h, qc_h, qr_h;

    ReadCustomBinaryIC(filename, latvec_h, lonvec_h,
                       xvec_h, yvec_h, zvec_h, rho_h,
                       uvel_h, vvel_h, wvel_h,
                       theta_h, qv_h, qc_h, qr_h);

    const auto prob_lo_erf  = geom[0].ProbLoArray();
    const auto prob_hi_erf  = geom[0].ProbHiArray();
    const auto dx_erf       = geom[0].CellSizeArray();

    if(prob_lo_erf[0] < xvec_h.front() + 4*dx_erf[0]){
        amrex::Abort("The xlo value of the domain has to be greater than " + std::to_string(xvec_h.front() + 4*dx_erf[0]));
    }
    if(prob_hi_erf[0] > xvec_h.back() - 4*dx_erf[0]){
        amrex::Abort("The xhi value of the domain has to be less than " + std::to_string(xvec_h.back() - 4*dx_erf[0]));
    }
    if(prob_lo_erf[1] < yvec_h.front() + 4*dx_erf[1]){
        amrex::Abort("The ylo value of the domain has to be greater than " + std::to_string(yvec_h.front() + 4*dx_erf[1]));
    }
    if(prob_hi_erf[1] > yvec_h.back() - 4*dx_erf[1]){
        amrex::Abort("The yhi value of the domain has to be less than " + std::to_string(yvec_h.back() - 4*dx_erf[1]));
    }

    // Number of cells
    int nx_cells = xvec_h.size()-1;
    int ny_cells = yvec_h.size()-1;

    const amrex::Geometry& geom0 = geom[0]; // or whatever your Geometry vector is called
    const amrex::Box& domainBox = geom0.Domain();
    const amrex::IntVect& domainSize = domainBox.size(); // Number of cells in each direction
    int nz_cells = domainSize[2];


    int ncomp = 10;
    int ngrow = 0;

    int n_time = 1;      // or however many time slices you want
    weather_forecast_data.resize(n_time);
    MultiFab& weather_mf = weather_forecast_data[0];
    weather_mf.define(nba, dm, ncomp, ngrow);

    fill_weather_data_multifab(weather_mf, geom_weather, nx_cells+1, ny_cells+1, nz_cells+1,
                               latvec_h, lonvec_h, zvec_h,
                               rho_h,uvel_h, vvel_h, wvel_h,
                               theta_h, qv_h, qc_h, qr_h);

    PlotMultiFab(weather_mf, geom_weather, "plt_coarse_weather", MultiFabType::NC);
}


void
ERF::WeatherDataInterpolation(const Real time)
{
    static Real next_read_forecast_time = -1.0;

    if (next_read_forecast_time < 0.0) {
        int next_multiple = static_cast<int>(time / 10800.0);
        next_read_forecast_time = next_multiple * 10800.0;
    }

    if (time >= next_read_forecast_time) {

        std::string folder = "WeatherData";

        // Check if folder exists and is a directory
        if (!fs::exists(folder) || !fs::is_directory(folder)) {
            throw std::runtime_error("Error: Folder '" + folder + "' does not exist or is not a directory.");
        }

        std::vector<std::string> bin_files;
        for (const auto& entry : fs::directory_iterator(folder)) {
            if (entry.is_regular_file() && entry.path().extension() == ".bin") {
                bin_files.push_back(entry.path().string());
            }
        }

        // 2. Sort lexicographically
        std::sort(bin_files.begin(), bin_files.end());


        for (const auto& entry : fs::directory_iterator(folder)) {
            if (entry.is_regular_file() && entry.path().extension() == ".bin") {
                bin_files.push_back(entry.path().string());
            }
        }

    // Check if no .bin files were found
        if (bin_files.empty()) {
            throw std::runtime_error("Error: No .bin files found in folder '" + folder + "'.");
        }

        std::string filename1, filename2;
        Vector<MultiFab> weather_forecast_data_1, weather_forecast_data_2;
        amrex::Geometry geom_weather;
        BoxArray nba;
        DistributionMapping dm;

        int idx1 = static_cast<int>(time / 10800.0);
        int idx2 = static_cast<int>(time / 10800.0)+1;

        if (idx2 >= static_cast<int>(bin_files.size())) {
            throw std::runtime_error("Error: Not enough .bin files to cover time " + std::to_string(time));
        }

        filename1 = bin_files[idx1];
        filename2 = bin_files[idx2];

        //Read in weather_forecast_1
        CreateWeatherDataGeomBoxArrayDistMap(filename1,
                                             geom_weather,
                                             nba,
                                             dm);

        FillWeatherDataMultiFab(filename1,
                                geom_weather,
                                nba,
                                dm,
                                weather_forecast_data_1);

        FillWeatherDataMultiFab(filename2,
                                geom_weather,
                                nba,
                                dm,
                                weather_forecast_data_2);

        //Interpolate in time to get the weather_forecast_interp
        int ncomp = weather_forecast_data_1[0].nComp();
        MultiFab weather_forecast_interp(nba, dm, ncomp, 0);
        Real alpha1 = 1.0 - (time - next_read_forecast_time)/10800.0;
        Real alpha2 = 1.0 - alpha1;

        MultiFab::LinComb(weather_forecast_interp,
                          alpha1, weather_forecast_data_1[0], 0,
                          alpha2, weather_forecast_data_2[0], 0,
                             0, ncomp, 0);

        //Interpolate in space to get the erf_forecast_interp
        InterpWeatherDataOntoMesh(geom_weather, weather_forecast_data_1[0]);
        next_read_forecast_time += 10800.0;
    }
}

#endif
