#ifndef ERF_FUEL_MAP_H
#define ERF_FUEL_MAP_H

#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>

/**
 * @file ERF_FuelMap.H
 * @brief Spatial fuel map readers (ASCII and binary LCP formats)
 *
 * Ported from wildfire_levelset/src/:
 *   - landscape_file.H (ASCII ESRI grid reader)
 *   - lcp_binary_reader.H (FARSITE LCP reader)
 *
 * Key design patterns:
 * - File I/O on rank 0 only, broadcast via MPI
 * - Device copy after broadcast for GPU kernel access
 * - Proper row reordering: file row 0 = northernmost (highest y)
 *   → reverse to match domain ordering
 */

/**
 * @brief Read an ESRI ASCII grid fuel map file.
 *
 * Format:
 *   ncols <N>
 *   nrows <N>
 *   xllcorner <x>
 *   yllcorner <y>
 *   cellsize <dx>
 *   nodata_value <val>
 *   <rows of integer codes, row 0 = northernmost>
 *
 * Only rank 0 reads the file. The vector is broadcast to all ranks via MPI_Bcast:
 * first the count, then the data.
 *
 * Cells with nodata_value are assigned code 0.
 *
 * Rows are reordered so that fuel_codes[j * ncols + i] corresponds to
 * domain row j at prob_lo_y (file row 0 = northernmost = highest y, so reverse).
 *
 * @param[in]  filename     Path to ASCII fuel map file
 * @param[in]  fire_nx      Fire grid resolution in x
 * @param[in]  fire_ny      Fire grid resolution in y
 * @param[out] fuel_codes   Vector of integer fuel codes (size fire_nx * fire_ny)
 * @param[out] nodata_val   The nodata_value from the header
 * @return true if file read successfully, false if file cannot be opened
 */
inline bool read_ascii_fuel_map(const std::string& filename,
                                int fire_nx, int fire_ny,
                                std::vector<int>& fuel_codes,
                                int& nodata_val)
{
    int ncols = 0, nrows = 0;
    double xll = 0.0, yll = 0.0, cellsize = 0.0;
    nodata_val = -9999;
    fuel_codes.clear();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename);
        if (!file.good()) {
            return false;  // File cannot be opened
        }

        // Read header
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);
        std::string key;
        iss >> key >> ncols;

        std::getline(file, line);
        iss.clear();
        iss.str(line);
        iss >> key >> nrows;

        std::getline(file, line);
        iss.clear();
        iss.str(line);
        iss >> key >> xll;

        std::getline(file, line);
        iss.clear();
        iss.str(line);
        iss >> key >> yll;

        std::getline(file, line);
        iss.clear();
        iss.str(line);
        iss >> key >> cellsize;

        std::getline(file, line);
        iss.clear();
        iss.str(line);
        iss >> key >> nodata_val;

        // Read data rows (row 0 = northernmost)
        // Reverse row order for domain consistency
        std::vector<int> temp_data(static_cast<size_t>(ncols * nrows));
        for (int row = nrows - 1; row >= 0; --row) {
            for (int col = 0; col < ncols; ++col) {
                int val;
                file >> val;
                if (val == nodata_val) {
                    val = 0;  // Assign code 0 to nodata cells
                }
                temp_data[static_cast<size_t>((nrows - 1 - row) * ncols + col)] = val;
            }
        }

        file.close();
        fuel_codes = temp_data;
    }

    // Broadcast count first
    int count = static_cast<int>(fuel_codes.size());
    amrex::ParallelDescriptor::Bcast(&count, 1);
    amrex::ParallelDescriptor::Bcast(&nodata_val, 1);

    // Resize on non-IO processors
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        fuel_codes.resize(count);
    }

    // Broadcast data
    if (count > 0) {
        amrex::ParallelDescriptor::Bcast(fuel_codes.data(), count);
    }

    return true;
}

// Include LCP reader implementation
#include <ERF_LcpReader.H>

/**
 * @brief Fill a MultiFab with fuel model codes from a device-side array.
 *
 * GPU kernel: for each cell (i,j), set fire_fuel_model(i,j,0) to the
 * integer code at d_fuel_codes[j * fire_nx + i], converted to Real.
 *
 * @param[out] fire_fuel_model MultiFab (1 component) to fill
 * @param[in]  d_fuel_codes    Device pointer to integer fuel codes
 * @param[in]  geom            Fire grid geometry
 * @param[in]  fire_nx         Fire grid resolution in x
 */
inline void fill_fuel_model_mf(amrex::MultiFab& fire_fuel_model,
                                const int* d_fuel_codes,
                                const amrex::Geometry& geom,
                                int fire_nx)
{
    // ParallelFor over the fire grid
    for (amrex::MFIter mfi(fire_fuel_model); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto& fuel_arr = fire_fuel_model.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            int idx = j * fire_nx + i;
            fuel_arr(i, j, k) = static_cast<amrex::Real>(d_fuel_codes[idx]);
        });
    }
}

#endif  // ERF_FUEL_MAP_H
