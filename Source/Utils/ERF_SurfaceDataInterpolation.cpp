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

void
ERF::FillSurfaceStateMultiFabs(const int lev,
                               const std::string& filename,
                               Vector<MultiFab>& surface_state)
{
       // Open the binary file in input mode
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }
    Vector<Real> xvec_h, yvec_h, zvec_h;
    Vector<Real> sst_h, q_star_h, t_star_h, u_star_h, ls_mask_h;

    int nx, ny, nz, ndata;
    float value;

    // Read the four integers
    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ndata), sizeof(int));

    amrex::Gpu::DeviceVector<Real> xvec_d(nx*ny*nz), yvec_d(nx*ny*nz), zvec_d(nz);
    for(int i=0; i<nx; i++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        xvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xvec_h.begin(), xvec_h.end(), xvec_d.begin());

    for(int j=0; j<ny; j++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        yvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yvec_h.begin(), yvec_h.end(), yvec_d.begin());

    for(int k=0; k<nz; k++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        zvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());

         // Vector to store the data

    Vector<Real>* data_h = nullptr; // Declare pointer outside the loop

    Real* xvec_d_ptr = xvec_d.data();
    Real* yvec_d_ptr = yvec_d.data();

    Real dxvec = (xvec_h[nx-1]-xvec_h[0])/(nx-1);
    Real dyvec = (yvec_h[ny-1]-yvec_h[0])/(ny-1);

    // Read the file
    for(int idx=0; idx<ndata; idx++){
        if(idx == 0){
            data_h = &sst_h;
        } else if (idx==1) {
            data_h = &q_star_h;
        } else if (idx==2) {
            data_h = &t_star_h;
        } else if (idx==3) {
            data_h = &u_star_h;
        } else if(idx==4) {
            data_h = &ls_mask_h;
        }
        for(int k=0; k<nz; k++) {
            for(int j=0; j<ny; j++) {
                for(int i=0; i<nx; i++) {
                    infile.read(reinterpret_cast<char*>(&value), sizeof(float));
                    //if(idx == 3) {
                        //printf("theta is %0.15g, %0.15g, %0.15g %0.15g\n", xvec_h[i], yvec_h[j], zvec_h[k], value);
                    //}
                    data_h->emplace_back(value);
                }
            }
        }
    }

    infile.close();

    amrex::Gpu::DeviceVector<Real> ls_mask_d(nx*ny*nz), sst_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ls_mask_h.begin(), ls_mask_h.end(), ls_mask_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, sst_h.begin(), sst_h.end(), sst_d.begin());

    Real* ls_mask_d_ptr = ls_mask_d.data();
    Real* sst_d_ptr   = sst_d.data();

    const auto prob_lo  = geom[lev].ProbLo();
    const auto dx       = geom[lev].CellSize();

    for (amrex::MFIter mfi(surface_state[lev]); mfi.isValid(); ++mfi) {
        const Box gbx = mfi.growntilebox();
        const Array4<Real>& surf_arr = surface_state[lev].array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            if(k == 0) {
                const Real x        = prob_lo[0] + (i + myhalf) * dx[0];
                const Real y        = prob_lo[1] + (j + myhalf) * dx[1];

                // First interpolate where the weather data is available from
                Real tmp_ls_mask, tmp_sst;

                bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                          dxvec, dyvec,
                                          nx, ny,
                                          x, y,
                                          ls_mask_d_ptr, tmp_ls_mask);

                bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                          dxvec, dyvec,
                                          nx, ny,
                                          x, y,
                                          sst_d_ptr, tmp_sst);

                surf_arr(i, j, 0) = std::min(tmp_ls_mask, one);
                surf_arr(i, j, 1) = tmp_sst;
            }
        });
    }

}

void
ERF::SurfaceDataInterpolation(const int lev,
                              const Real time,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd,
                              bool regrid_forces_file_read)
{

    static amrex::Vector<Real> next_read_forecast_time;
    static amrex::Vector<Real> last_read_forecast_time;

    const int nlevs = a_z_phys_nd.size();

    Real hindcast_data_interval = solverChoice.hindcast_data_interval_in_hrs*Real(3600.0);

    // Initialize static vectors once
    if (next_read_forecast_time.empty()) {
        next_read_forecast_time.resize(nlevs, -one);
        last_read_forecast_time.resize(nlevs, -one);
        Print() << "Initializing the time vector values here by " << lev << std::endl;
    }

    if (next_read_forecast_time[lev] < zero) {
        int next_multiple = static_cast<int>(time / hindcast_data_interval);
        next_read_forecast_time[lev] = next_multiple * hindcast_data_interval;
        last_read_forecast_time[lev] = next_read_forecast_time[lev];
    }

    if (time >= next_read_forecast_time[lev] or regrid_forces_file_read) {

        Print() << "Data reading happening at level " << lev << std::endl;

        std::string folder = solverChoice.hindcast_surface_data_dir;

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
        Print() << "Reading surface data " << time << " " << idx1 << " " << idx2 <<" " << bin_files.size() << std::endl;

        if (idx2 >= static_cast<int>(bin_files.size())) {
            throw std::runtime_error("Error: Not enough .bin files to cover time " + std::to_string(time));
        }

        filename1 = bin_files[idx1];
        filename2 = bin_files[idx2];

        FillSurfaceStateMultiFabs(lev, filename1, surface_state_1);
        FillSurfaceStateMultiFabs(lev, filename2, surface_state_2);

         // Create the time-interpolated forecast state
        //CreateForecastStateMultiFabs(forecast_state_interp);
        if(!regrid_forces_file_read){
            last_read_forecast_time[lev] = next_read_forecast_time[lev];
            next_read_forecast_time[lev] += hindcast_data_interval;
            Print() << "Next forecast time getting updated here " << std::endl;
        }
    }

    Real prev_read_time = last_read_forecast_time[lev];
    Real alpha1 = one - (time - prev_read_time)/hindcast_data_interval;
    Real alpha2 = one - alpha1;

    amrex::Print()<< "The values of alpha1 and alpha2 are " << alpha1 << " "<< alpha2 <<std::endl;

    if (alpha1 < zero || alpha1 > one ||
    alpha2 < zero || alpha2 > one)
    {
        std::stringstream ss;
        ss << "Interpolation weights for hindcast files are incorrect: "
        << "alpha1 = " << alpha1 << ", alpha2 = " << alpha2;
       Abort(ss.str());
    }

    /*MultiFab& mf_surf_interp   = surface_state_interp[lev];

    // Fill the time-interpolated forecast states
    MultiFab::LinComb(surface_state_interp[lev],
                      alpha1, surface_state_1[lev], 0,
                      alpha2, surface_state_2[lev], 0,
                      0, mf_surf_interp.nComp(), mf_surf_interp.nGrow());

    std::string pltname = "plt_interp_surface";
    Vector<std::string> varnames_plot_mf = {"ls_mask", "SST"};

    const MultiFab& src = vars_new[0][0];

    MultiFab plot_mf(src.boxArray(),
                     src.DistributionMap(),
                     2, 0);

    plot_mf.setVal(0.0);

    for (MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
        const Array4<Real> &plot_mf_arr = plot_mf.array(mfi);
        const Array4<Real> &surf_mf_arr = surface_state_1[0].array(mfi);

        const Box& bx = mfi.validbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plot_mf_arr(i,j,k,0) = surf_mf_arr(i,j,0);
            plot_mf_arr(i,j,k,1) = surf_mf_arr(i,j,1);
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
