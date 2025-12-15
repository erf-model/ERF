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
ERF::FillForecastStateMultiFabs(const int lev,
                                const std::string& filename,
                                const std::unique_ptr<MultiFab>& a_z_phys_nd,
                                Vector<Vector<MultiFab>>& forecast_state)
{

    Vector<Real> latvec_h, lonvec_h, xvec_h, yvec_h, zvec_h;
    Vector<Real> rho_h, uvel_h, vvel_h, wvel_h, theta_h, qv_h, qc_h, qr_h;

    ReadCustomBinaryIC(filename, latvec_h, lonvec_h,
                       xvec_h, yvec_h, zvec_h, rho_h,
                       uvel_h, vvel_h, wvel_h,
                       theta_h, qv_h, qc_h, qr_h);

    Real zmax = *std::max_element(zvec_h.begin(), zvec_h.end());

    const auto prob_lo_erf  = geom[lev].ProbLoArray();
    const auto prob_hi_erf  = geom[lev].ProbHiArray();
    const auto dx_erf       = geom[lev].CellSizeArray();

    if (prob_hi_erf[2] >= zmax) {
        Abort("ERROR: the maximum z of the domain (" + std::to_string(prob_hi_erf[2]) +
        ") should be less than the maximum z in the forecast data (" + std::to_string(zmax) +
        "). Change geometry.prob_hi[2] in the inputs to be less than " + std::to_string(zmax) + "."
        );
    }

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


    int nx = xvec_h.size();
    int ny = yvec_h.size();
    int nz = zvec_h.size();

    amrex::Real dxvec = (xvec_h[nx-1]-xvec_h[0])/(nx-1);
    amrex::Real dyvec = (yvec_h[ny-1]-yvec_h[0])/(ny-1);

    amrex::Gpu::DeviceVector<Real> latvec_d(nx*ny), lonvec_d(nx*ny), zvec_d(nz);
    amrex::Gpu::DeviceVector<Real> xvec_d(nx*ny*nz), yvec_d(nx*ny*nz);
    amrex::Gpu::DeviceVector<Real> rho_d(nx*ny*nz), uvel_d(nx*ny*nz), vvel_d(nx*ny*nz), wvel_d(nx*ny*nz),
                                   theta_d(nx*ny*nz), qv_d(nx*ny*nz), qc_d(nx*ny*nz), qr_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, latvec_h.begin(), latvec_h.end(), latvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, lonvec_h.begin(), lonvec_h.end(), lonvec_d.begin());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xvec_h.begin(), xvec_h.end(), xvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yvec_h.begin(), yvec_h.end(), yvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, uvel_h.begin(), uvel_h.end(), uvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, vvel_h.begin(), vvel_h.end(), vvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, wvel_h.begin(), wvel_h.end(), wvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qr_h.begin(), qr_h.end(), qr_d.begin());

    amrex::Gpu::streamSynchronize();

    Real* latvec_d_ptr = latvec_d.data();
    Real* lonvec_d_ptr = lonvec_d.data();
    Real* xvec_d_ptr = xvec_d.data();
    Real* yvec_d_ptr = yvec_d.data();
    Real* zvec_d_ptr = zvec_d.data();
    Real* rho_d_ptr   = rho_d.data();
    Real* uvel_d_ptr  = uvel_d.data();
    Real* vvel_d_ptr  = vvel_d.data();
    Real* wvel_d_ptr  = wvel_d.data();
    Real* theta_d_ptr = theta_d.data();
    Real* qv_d_ptr = qv_d.data();
    Real* qc_d_ptr = qc_d.data();
    Real* qr_d_ptr = qr_d.data();

    MultiFab& erf_mf_cons   = forecast_state[lev][Vars::cons];
    MultiFab& erf_mf_xvel   = forecast_state[lev][Vars::xvel];
    MultiFab& erf_mf_yvel   = forecast_state[lev][Vars::yvel];
    MultiFab& erf_mf_zvel   = forecast_state[lev][Vars::zvel];
    MultiFab& erf_mf_latlon = forecast_state[lev][4];

    erf_mf_cons.setVal(0.0);
    erf_mf_xvel.setVal(0.0);
    erf_mf_yvel.setVal(0.0);
    erf_mf_zvel.setVal(0.0);
    erf_mf_latlon.setVal(0.0);

    // Interpolate the data on to the ERF mesh

     for (MFIter mfi(erf_mf_cons); mfi.isValid(); ++mfi) {
        const auto z_arr    = (a_z_phys_nd) ? a_z_phys_nd->const_array(mfi) :
                                            Array4<const Real> {};
        const Array4<Real> &fine_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &fine_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &fine_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &fine_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &fine_latlon_arr = erf_mf_latlon.array(mfi);


        const Box& gbx = mfi.growntilebox(); // tilebox + ghost cells

        const Box &gtbx = mfi.tilebox(IntVect(1,0,0));
        const Box &gtby = mfi.tilebox(IntVect(0,1,0));
        const Box &gtbz = mfi.tilebox(IntVect(0,0,1));
        const auto prob_lo  = geom[lev].ProbLoArray();
        const auto dx       = geom[lev].CellSizeArray();
       //const Box &gtbz = mfi.tilebox(IntVect(0,0,1));

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Geometry (note we must include these here to get the data on device)
            const Real x        = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y        = prob_lo[1] + (j + 0.5) * dx[1];
            //const Real z        = prob_lo[2] + (k + 0.5) * dx[2];
            const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/2.0;

            // First interpolate where the weather data is available from
            Real tmp_rho, tmp_theta, tmp_qv, tmp_qc, tmp_qr, tmp_lat, tmp_lon;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   rho_d_ptr, tmp_rho);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   theta_d_ptr, tmp_theta);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   qv_d_ptr, tmp_qv);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   qc_d_ptr, tmp_qc);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   qr_d_ptr, tmp_qr);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, 1,
                                   x, y, 0.0,
                                   latvec_d_ptr, tmp_lat);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, 1,
                                   x, y, 0.0,
                                   lonvec_d_ptr, tmp_lon);

            fine_cons_arr(i,j,k,Rho_comp) = tmp_rho;
            fine_latlon_arr(i,j,k,0) = tmp_lat;
            fine_latlon_arr(i,j,k,1) = tmp_lon;
        });

        ParallelFor(gtbx, gtby, gtbz,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + i       * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+0.5) * dx_erf[1];
            //Real z = prob_lo_erf[2] + (k+0.5) * dx_erf[2];
            const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/2.0;

            Real tmp_uvel;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   uvel_d_ptr, tmp_uvel);

            fine_xvel_arr(i, j, k, 0) = tmp_uvel;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+0.5) * dx_erf[0];
            Real y = prob_lo_erf[1] + j       * dx_erf[1];
            //Real z = prob_lo_erf[2] + (k+0.5) * dx_erf[2];
            const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/2.0;

            Real tmp_vvel;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   vvel_d_ptr, tmp_vvel);

            fine_yvel_arr(i, j, k, 0) = tmp_vvel;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+0.5) * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+0.5) * dx_erf[1];
            Real z = prob_lo_erf[2] + k       * dx_erf[2];
            //const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/2.0;

            Real tmp_wvel;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   wvel_d_ptr, tmp_wvel);

            fine_zvel_arr(i, j, k, 0) = tmp_wvel;
        });
    }

    /*Vector<std::string> varnames = {
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
        );*/
}

void
ERF::WeatherDataInterpolation(const int lev,
                              const Real time,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd,
                              bool regrid_forces_file_read)
{

    static amrex::Vector<Real> next_read_forecast_time;
    static amrex::Vector<Real> last_read_forecast_time;

    const int nlevs = a_z_phys_nd.size();

    Real hindcast_data_interval = solverChoice.hindcast_data_interval_in_hrs*3600.0;

    // Initialize static vectors once
    if (next_read_forecast_time.empty()) {
        next_read_forecast_time.resize(nlevs, -1.0);
        last_read_forecast_time.resize(nlevs, -1.0);
        Print() << "Initializing the time vector values here by " << lev << std::endl;
    }

    if (next_read_forecast_time[lev] < 0.0) {
        int next_multiple = static_cast<int>(time / hindcast_data_interval);
        next_read_forecast_time[lev] = next_multiple * hindcast_data_interval;
        last_read_forecast_time[lev] = next_read_forecast_time[lev];
    }

    if (time >= next_read_forecast_time[lev] or regrid_forces_file_read) {

        Print() << "Data reading happening at level " << lev << std::endl;

        std::string folder = solverChoice.hindcast_boundary_data_dir;

        // Check if folder exists and is a directory
        if (!fs::exists(folder) || !fs::is_directory(folder)) {
            throw std::runtime_error("Error: Folder '" + folder + "' does not exist or is not a directory.");
        }

        std::vector<std::string> bin_files;

        for (const auto& entry : fs::directory_iterator(folder)) {
            if (!entry.is_regular_file()) continue;

            std::string fname = entry.path().filename().string();
            if (fname.size() >= 4 && fname.substr(fname.size() - 4) == ".bin") {
                bin_files.push_back(entry.path().string());
            }
        }
        std::sort(bin_files.begin(), bin_files.end());

    // Check if no .bin files were found
        if (bin_files.empty()) {
            throw std::runtime_error("Error: No .bin files found in folder '" + folder + "'.");
        }

        std::string filename1, filename2;

        int idx1 = static_cast<int>(time / hindcast_data_interval);
        int idx2 = static_cast<int>(time / hindcast_data_interval)+1;
        std::cout << "Reading weather data " << time << " " << idx1 << " " << idx2 <<" " << bin_files.size() << std::endl;

        if (idx2 >= static_cast<int>(bin_files.size())) {
            throw std::runtime_error("Error: Not enough .bin files to cover time " + std::to_string(time));
        }

        filename1 = bin_files[idx1];
        filename2 = bin_files[idx2];

        FillForecastStateMultiFabs(lev, filename1, a_z_phys_nd[lev], forecast_state_1);
        FillForecastStateMultiFabs(lev, filename2, a_z_phys_nd[lev], forecast_state_2);

         // Create the time-interpolated forecast state
        //CreateForecastStateMultiFabs(forecast_state_interp);
        if(!regrid_forces_file_read){
            last_read_forecast_time[lev] = next_read_forecast_time[lev];
            next_read_forecast_time[lev] += hindcast_data_interval;
            Print() << "Next forecast time getting updated here " << std::endl;
        }
    }

    Real prev_read_time = last_read_forecast_time[lev];
    Real alpha1 = 1.0 - (time - prev_read_time)/hindcast_data_interval;
    Real alpha2 = 1.0 - alpha1;

    amrex::Print()<< "The values of alpha1 and alpha2 are " << alpha1 << " "<< alpha2 <<std::endl;

    if (alpha1 < 0.0 || alpha1 > 1.0 ||
    alpha2 < 0.0 || alpha2 > 1.0)
    {
        std::stringstream ss;
        ss << "Interpolation weights for hindcast files are incorrect: "
        << "alpha1 = " << alpha1 << ", alpha2 = " << alpha2;
       Abort(ss.str());
    }

    MultiFab& erf_mf_cons   = forecast_state_interp[lev][Vars::cons];
    MultiFab& erf_mf_xvel   = forecast_state_interp[lev][Vars::xvel];
    MultiFab& erf_mf_yvel   = forecast_state_interp[lev][Vars::yvel];
    //MultiFab& erf_mf_zvel   = forecast_state_interp[0][Vars::zvel];
    MultiFab& erf_mf_latlon = forecast_state_interp[lev][4];

    // Fill the time-interpolated forecast states
    MultiFab::LinComb(forecast_state_interp[lev][Vars::cons],
                      alpha1, forecast_state_1[lev][Vars::cons], 0,
                      alpha2, forecast_state_2[lev][Vars::cons], 0,
                      0, erf_mf_cons.nComp(), forecast_state_interp[lev][Vars::cons].nGrow());
    MultiFab::LinComb(forecast_state_interp[lev][Vars::xvel],
                      alpha1, forecast_state_1[lev][Vars::xvel], 0,
                      alpha2, forecast_state_2[lev][Vars::xvel], 0,
                      0, erf_mf_xvel.nComp(), forecast_state_interp[lev][Vars::xvel].nGrow());
    MultiFab::LinComb(forecast_state_interp[lev][Vars::yvel],
                      alpha1, forecast_state_1[lev][Vars::yvel], 0,
                      alpha2, forecast_state_2[lev][Vars::yvel], 0,
                      0, erf_mf_yvel.nComp(), forecast_state_interp[lev][Vars::yvel].nGrow());
    MultiFab::LinComb(forecast_state_interp[lev][4],
                      alpha1, forecast_state_1[lev][4], 0,
                      alpha2, forecast_state_2[lev][4], 0,
                      0, erf_mf_latlon.nComp(), forecast_state_interp[lev][4].nGrow());

    /*Vector<std::string> varnames_plot_mf = {
    "rho", "rhotheta", "rhoqv", "rhoqc", "rhoqr", "xvel", "yvel", "zvel", "latitude", "longitude"
    }; // Customize variable names

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
        );*/
}

#endif
