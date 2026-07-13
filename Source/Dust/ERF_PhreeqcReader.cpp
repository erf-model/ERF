#ifdef ERF_USE_DUST

#include <ERF_PhreeqcReader.H>
#include <ERF_DustGrid.H>
#include <ERF_DustParams.H>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

using namespace amrex;

bool read_phreeqc_csv(MultiFab& mf, const DustGrid& dg,
                      const std::string& filename,
                      const std::string& col_name,
                      Real nodata_fill)
{
    // Return false if filename is empty
    if (filename.empty()) {
        return false;
    }

    int nx = 0, ny = 0, col_idx = -1;
    std::vector<Real> data;

    // Rank 0 reads the file
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            amrex::Print() << "[DUST] WARNING: Could not open file: " << filename << "\n";
            return false;
        }

        std::string line;

        // Read header row
        if (!std::getline(file, line)) {
            amrex::Print() << "[DUST] read_phreeqc_csv: could not read header from " << filename << "\n";
            file.close();
            return false;
        }

        // Parse header to find column index
        std::istringstream header_stream(line);
        std::string col;
        int idx = 0;
        while (std::getline(header_stream, col, ',')) {
            // Trim whitespace
            col.erase(0, col.find_first_not_of(" \t\r\n"));
            col.erase(col.find_last_not_of(" \t\r\n") + 1);
            if (col == col_name) {
                col_idx = idx;
                break;
            }
            ++idx;
        }

        if (col_idx < 0) {
            amrex::Print() << "[DUST] read_phreeqc_csv: column '" << col_name
                           << "' not found in " << filename << "\n";
            file.close();
            return false;
        }

        // Get dust grid dimensions
        Box dust_domain = dg.ba.minimalBox();
        nx = dust_domain.length(0);
        ny = dust_domain.length(1);

        // Read data rows
        // CSV row order: row 0 = southernmost (j=0). No reversal needed, unlike ESRI ASCII.
        // Each CSV data row contains: i, j, and then the various property columns.
        // We extract the column specified by col_name.
        data.resize(nx * ny);
        int row_count = 0;
        while (std::getline(file, line) && row_count < nx * ny) {
            std::istringstream row_stream(line);
            std::string cell;
            int cell_idx = 0;
            bool found = false;
            while (std::getline(row_stream, cell, ',')) {
                // Trim whitespace
                cell.erase(0, cell.find_first_not_of(" \t\r\n"));
                cell.erase(cell.find_last_not_of(" \t\r\n") + 1);
                
                if (cell_idx == col_idx) {
                    try {
                        Real val = std::stod(cell);
                        if (std::abs(val - PhreeqcDustConst::NODATA_CSV) < 1.0e-6) {
                            val = nodata_fill;
                        }
                        data[row_count] = val;
                        found = true;
                    } catch (...) {
                        amrex::Print() << "[DUST] read_phreeqc_csv: could not parse cell at row "
                                      << row_count << ", col " << col_idx << "\n";
                        file.close();
                        return false;
                    }
                    break;
                }
                ++cell_idx;
            }
            if (!found) {
                amrex::Print() << "[DUST] read_phreeqc_csv: column " << col_idx
                              << " not found in row " << row_count << "\n";
                file.close();
                return false;
            }
            ++row_count;
        }

        if (row_count < nx * ny) {
            amrex::Print() << "[DUST] read_phreeqc_csv: expected " << (nx*ny) << " rows but got "
                          << row_count << "\n";
            file.close();
            return false;
        }

        file.close();
    }

    // Broadcast dimensions from rank 0
    ParallelDescriptor::Bcast(&nx, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&ny, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&col_idx, 1, ParallelDescriptor::IOProcessorNumber());

    // Resize data on non-IO ranks
    if (!ParallelDescriptor::IOProcessor()) {
        data.resize(nx * ny);
    }

    // Broadcast data
    ParallelDescriptor::Bcast(data.data(), nx * ny, 
                               ParallelDescriptor::IOProcessorNumber());

    // If col_idx < 0, the column was not found on rank 0
    if (col_idx < 0) {
        return false;
    }

    // Copy data to device
    Gpu::DeviceVector<Real> d_data(data.size());
    Gpu::copy(Gpu::hostToDevice, data.begin(), data.end(), d_data.begin());

    // GPU ParallelFor over MultiFab tiles
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto mf_arr = mf.array(mfi);
        Real* d_data_ptr = d_data.data();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Clamp to domain bounds
            int ii = amrex::min(i, nx - 1);
            int jj = amrex::min(j, ny - 1);
            mf_arr(i, j, k) = d_data_ptr[jj * nx + ii];
        });
    }

    return true;
}

bool read_phreeqc_netcdf(MultiFab& mf, const DustGrid& dg,
                         const std::string& filename,
                         const std::string& varname, Real nodata_fill)
{
#ifdef ERF_USE_NETCDF
    // NetCDF-C API: https://docs.unidata.ucar.edu/netcdf-c/current/
    // Full implementation deferred. Use CSV format for now.
    amrex::Abort("[DUST] read_phreeqc_netcdf: implementation not yet complete. "
                 "Use CSV format for PHREEQC output.");
    return false;
#else
    amrex::Abort("[DUST] read_phreeqc_netcdf requires ERF_ENABLE_NETCDF=ON");
    return false;
#endif
}

void update_ustar_t_from_chemistry(MultiFab& ustar_t, const MultiFab& ustar_base,
                                   const MultiFab& crust, const MultiFab& efflor)
{
    // u*_t reduction from mineral crust and salt efflorescence.
    // Marticorena & Bergametti (1995), https://doi.org/10.1029/95JD00690
    const Real alpha_c    = PhreeqcDustConst::ALPHA_CRUST;
    const Real alpha_e    = PhreeqcDustConst::ALPHA_EFFLOR;
    const Real ustar_tmin = PhreeqcDustConst::USTAR_T_MIN;

    for (MFIter mfi(ustar_t, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx       = mfi.tilebox();
        auto ut_arr         = ustar_t.array(mfi);
        auto ut_base_arr    = ustar_base.const_array(mfi);
        auto crust_arr      = crust.const_array(mfi);
        auto efflor_arr     = efflor.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real ut = ut_base_arr(i,j,k)
                    * (Real(1.0) - alpha_c * crust_arr(i,j,k))
                    * (Real(1.0) - alpha_e * efflor_arr(i,j,k));
            ut_arr(i,j,k) = amrex::max(ut, ustar_tmin);
        });
    }
}

void update_dust_from_phreeqc(MultiFab&       dust_ustar_t,
                              const MultiFab& dust_ustar_base,
                              MultiFab&       dust_crust_index,
                              MultiFab&       dust_silt_frac,
                              MultiFab&       dust_efflor,
                              MultiFab&       dust_suppression,
                              MultiFab&       dust_emission,
                              const DustGrid& dg,
                              const DustParams& params)
{
    // Return immediately if phreeqc_output_file is empty
    if (params.phreeqc_output_file.empty()) {
        return;
    }

    // Helper lambda to select reader by file extension
    auto read_field = [&](MultiFab& mf, const std::string& varname) -> bool {
        const std::string& filename = params.phreeqc_output_file;
        if (filename.size() > 3) {
            std::string ext = filename.substr(filename.size() - 3);
            if (ext == ".nc") {
                return read_phreeqc_netcdf(mf, dg, filename, varname, 0.0);
            }
        }
        // Default to CSV
        return read_phreeqc_csv(mf, dg, filename, varname, 0.0);
    };

    // Read crust_index if variable name is set
    if (!params.phreeqc_crust_var.empty()) {
        read_field(dust_crust_index, params.phreeqc_crust_var);
    }

    // Read efflorescence if variable name is set
    if (!params.phreeqc_efflor_var.empty()) {
        read_field(dust_efflor, params.phreeqc_efflor_var);
    }

    // Update u*_t from chemistry after reading crust and efflorescence
    update_ustar_t_from_chemistry(dust_ustar_t, dust_ustar_base,
                                  dust_crust_index, dust_efflor);

    // Read silt_fraction if variable name is set
    if (!params.phreeqc_silt_var.empty()) {
        read_field(dust_silt_frac, params.phreeqc_silt_var);
    }

    // Read suppression modifier if variable name is set
    if (!params.phreeqc_supp_var.empty()) {
        read_field(dust_suppression, params.phreeqc_supp_var);
    }

    // Read metal mass fraction (As) for bin 0 if variable name is set
    if (!params.phreeqc_metal_var.empty()) {
        // Create a temporary MultiFab for reading
        MultiFab temp_mf(dg.ba, dg.dm, 1, IntVect(1,1,0));
        if (read_field(temp_mf, params.phreeqc_metal_var)) {
            // Copy to component 0 of dust_emission
            MultiFab::Copy(dust_emission, temp_mf, 0, 0, 1, IntVect(1,1,0));
        }
    }

    // Print status message
    amrex::Print() << "[DUST] PHREEQC update from file: " 
                   << params.phreeqc_output_file << "\n";
}

#endif // ERF_USE_DUST
