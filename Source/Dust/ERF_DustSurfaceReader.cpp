#ifdef ERF_USE_DUST

#include <ERF_DustSurfaceReader.H>
#include <ERF_DustGrid.H>
#include <ERF_DustParams.H>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

using namespace amrex;

bool read_ascii_surface_map(MultiFab& mf, const DustGrid& dg,
                            const std::string& filename, Real nodata_fill)
{
    // Return false if filename is empty
    if (filename.empty()) {
        return false;
    }

    int ncols = 0, nrows = 0;
    Real xllcorner = 0.0, yllcorner = 0.0, cellsize = 0.0, nodata_value = -9999.0;
    std::vector<Real> data;

    // Rank 0 reads the file
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            amrex::Print() << "[DUST] WARNING: Could not open file: " << filename << "\n";
            return false;
        }

        std::string line;
        
        // Helper lambda to parse header line (case-insensitive key)
        auto parse_header = [](const std::string& line, const std::string& key_lower)
            -> std::pair<bool, Real> {
            std::string line_lower = line;
            std::transform(line_lower.begin(), line_lower.end(), 
                         line_lower.begin(), 
                         [](unsigned char c) { return std::tolower(c); });
            
            std::istringstream iss(line_lower);
            std::string key;
            Real value;
            if (iss >> key >> value) {
                if (key == key_lower) {
                    return {true, value};
                }
            }
            return {false, 0.0};
        };

        // Read 6 header lines
        // Line 1: ncols
        if (std::getline(file, line)) {
            auto [ok, val] = parse_header(line, "ncols");
            if (ok) ncols = static_cast<int>(val);
        }
        // Line 2: nrows
        if (std::getline(file, line)) {
            auto [ok, val] = parse_header(line, "nrows");
            if (ok) nrows = static_cast<int>(val);
        }
        // Line 3: xllcorner
        if (std::getline(file, line)) {
            auto [ok, val] = parse_header(line, "xllcorner");
            if (ok) xllcorner = val;
        }
        // Line 4: yllcorner
        if (std::getline(file, line)) {
            auto [ok, val] = parse_header(line, "yllcorner");
            if (ok) yllcorner = val;
        }
        // Line 5: cellsize
        if (std::getline(file, line)) {
            auto [ok, val] = parse_header(line, "cellsize");
            if (ok) cellsize = val;
        }
        // Line 6: nodata_value
        if (std::getline(file, line)) {
            auto [ok, val] = parse_header(line, "nodata_value");
            if (ok) nodata_value = val;
        }

        // Read data rows
        data.resize(ncols * nrows);
        for (int j = 0; j < nrows; ++j) {
            for (int i = 0; i < ncols; ++i) {
                if (!(file >> data[j * ncols + i])) {
                    amrex::Print() << "[DUST] ERROR: Could not read data at row " << j << ", col " << i << "\n";
                    file.close();
                    return false;
                }
            }
        }
        file.close();

        // Row reversal: file row 0 is northernmost; domain row 0 is southernmost. 
        // Matches ERF_FuelMap.H convention.
        std::vector<Real> data_reversed(ncols * nrows);
        for (int j = 0; j < nrows; ++j) {
            for (int i = 0; i < ncols; ++i) {
                Real val = data[j * ncols + i];
                if (std::abs(val - nodata_value) < 1.0e-10) {
                    val = nodata_fill;
                }
                data_reversed[(nrows - 1 - j) * ncols + i] = val;
            }
        }
        data = data_reversed;
    }

    // Broadcast dimensions from rank 0
    ParallelDescriptor::Bcast(&ncols, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&nrows, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&xllcorner, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&yllcorner, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&cellsize, 1, ParallelDescriptor::IOProcessorNumber());

    // Resize data on non-IO ranks
    if (!ParallelDescriptor::IOProcessor()) {
        data.resize(ncols * nrows);
    }

    // Broadcast data
    ParallelDescriptor::Bcast(data.data(), ncols * nrows, 
                               ParallelDescriptor::IOProcessorNumber());

    // Copy data to device
    Gpu::DeviceVector<Real> d_data(data.size());
    Gpu::copy(Gpu::hostToDevice, data.begin(), data.end(), d_data.begin());

    // Build uniform source-grid coordinate vectors on device
    Gpu::DeviceVector<Real> x_src(ncols);
    Gpu::DeviceVector<Real> y_src(nrows);
    
    // Host vectors to copy from
    std::vector<Real> x_src_host(ncols), y_src_host(nrows);
    for (int i = 0; i < ncols; ++i) {
        x_src_host[i] = xllcorner + (i + 0.5) * cellsize;
    }
    for (int j = 0; j < nrows; ++j) {
        y_src_host[j] = yllcorner + (j + 0.5) * cellsize;
    }
    Gpu::copy(Gpu::hostToDevice, x_src_host.begin(), x_src_host.end(), x_src.begin());
    Gpu::copy(Gpu::hostToDevice, y_src_host.begin(), y_src_host.end(), y_src.begin());

    // Get dust grid data
    const Geometry& geom = dg.geom;
    // GpuArray by value: CellSize()/ProbLo() return pointers into host Geometry
    // members, and the interpolation kernel below dereferences these.
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    // GPU kernel: bilinear interpolation over MultiFab tiles
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto arr = mf.array(mfi);

        // Capture device pointers
        const Real* d_data_ptr = d_data.data();
        const Real* x_src_ptr = x_src.data();
        const Real* y_src_ptr = y_src.data();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Compute dust-grid cell center
            Real x_dust = prob_lo[0] + (i + 0.5) * dx[0];
            Real y_dust = prob_lo[1] + (j + 0.5) * dx[1];

            // Clamp to source domain (for ghost cells)
            Real x_src_min = x_src_ptr[0];
            Real x_src_max = x_src_ptr[ncols - 1];
            Real y_src_min = y_src_ptr[0];
            Real y_src_max = y_src_ptr[nrows - 1];

            x_dust = amrex::max(amrex::min(x_dust, x_src_max), x_src_min);
            y_dust = amrex::max(amrex::min(y_dust, y_src_max), y_src_min);

            // Find bounding source cells
            int i_left = 0, j_bottom = 0;
            
            // Linear search for x
            for (int ii = 0; ii < ncols - 1; ++ii) {
                if (x_dust >= x_src_ptr[ii] && x_dust <= x_src_ptr[ii + 1]) {
                    i_left = ii;
                    break;
                }
            }
            
            // Linear search for y
            for (int jj = 0; jj < nrows - 1; ++jj) {
                if (y_dust >= y_src_ptr[jj] && y_dust <= y_src_ptr[jj + 1]) {
                    j_bottom = jj;
                    break;
                }
            }

            // Get four corner values
            Real v00 = d_data_ptr[j_bottom * ncols + i_left];
            Real v10 = d_data_ptr[j_bottom * ncols + i_left + 1];
            Real v01 = d_data_ptr[(j_bottom + 1) * ncols + i_left];
            Real v11 = d_data_ptr[(j_bottom + 1) * ncols + i_left + 1];

            // Bilinear interpolation weights
            Real dx_cell = x_src_ptr[i_left + 1] - x_src_ptr[i_left];
            Real dy_cell = y_src_ptr[j_bottom + 1] - y_src_ptr[j_bottom];
            
            Real wx = (x_dust - x_src_ptr[i_left]) / dx_cell;
            Real wy = (y_dust - y_src_ptr[j_bottom]) / dy_cell;

            Real v0 = v00 * (1.0 - wx) + v10 * wx;
            Real v1 = v01 * (1.0 - wx) + v11 * wx;
            Real v = v0 * (1.0 - wy) + v1 * wy;

            arr(i, j, k, 0) = v;
        });
    }

    // ParallelFor is asynchronous; d_data, x_src and y_src free their device
    // allocations when this function returns. Let the kernels finish first.
    Gpu::streamSynchronize();

    return true;
}

bool read_netcdf_surface_map(MultiFab& mf, const DustGrid& dg,
                             const std::string& filename,
                             const std::string& varname,
                             Real nodata_fill)
{
    // Both branches below abort, so nothing here consumes the arguments.
    amrex::ignore_unused(mf, dg, filename, varname, nodata_fill);
#ifdef ERF_USE_NETCDF
    // NetCDF-C API: https://docs.unidata.ucar.edu/netcdf-c/current/
    // Full implementation deferred. Use ESRI ASCII format for now.
    amrex::Abort("[DUST] read_netcdf_surface_map: implementation not yet complete. "
                 "Use ESRI ASCII (.asc) format.");
    return false;
#else
    amrex::Abort("[DUST] read_netcdf_surface_map requires ERF_ENABLE_NETCDF=ON");
    return false;
#endif
}

void populate_dust_surface_maps(MultiFab& soil, MultiFab& silt,
                                MultiFab& crust, MultiFab& moist,
                                MultiFab& supp,
                                const DustGrid& dg, const DustParams& p)
{
    // Select reader by file extension: .nc -> NetCDF, else ESRI ASCII.
    auto read_map = [&](MultiFab& mf, const std::string& fname,
                        Real fill, const std::string& vname) {
        if (fname.empty()) return;
        bool nc = fname.size() >= 3 &&
                  fname.substr(fname.size()-3) == ".nc";
        bool ok = nc ? read_netcdf_surface_map(mf, dg, fname, vname, fill)
                     : read_ascii_surface_map(mf, dg, fname, fill);
        if (ok) amrex::Print() << "[DUST] Loaded " << vname
                               << " from: " << fname << "\n";
    };

    read_map(soil,  p.soil_type_file,     0.0,           "soil_type");
    read_map(silt,  p.silt_fraction_file, p.silt_fraction,"silt_fraction");
    read_map(crust, p.crust_index_file,   p.crust_index,  "crust_index");
    read_map(moist, p.moisture_flag_file, 0.0,            "moisture_flag");
    read_map(supp,  p.suppression_file,   0.0,            "suppression");
}

#endif // ERF_USE_DUST
