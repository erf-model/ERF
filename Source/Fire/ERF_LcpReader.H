#ifndef ERF_LCP_READER_H
#define ERF_LCP_READER_H

#include <AMReX_ParallelDescriptor.H>
#include <string>
#include <vector>
#include <fstream>
#include <cstring>

/**
 * @file ERF_LcpReader.H
 * @brief FARSITE LCP binary format reader
 *
 * Ported from wildfire_levelset/src/lcp_binary_reader.H
 *
 * LCP (Landscape) file format:
 * - 256-byte header
 * - Stacked 16-bit integer layers in the standard FARSITE layer order:
 *   0: elevation [m]
 *   1: slope [%]
 *   2: aspect [degrees, 0-360, 0=N]
 *   3: canopy cover [%]
 *   4: fuel model [code, typically 1-13]
 *   ...and more
 *
 * The fuel model is at layer index 4, stored as 16-bit integers.
 * File rows: row 0 = northernmost (highest y). Must reverse to match domain.
 *
 * Only rank 0 reads the file; MPI_Bcast distributes to all ranks.
 */

/**
 * @brief Read a FARSITE LCP binary fuel map file.
 *
 * Reads the 256-byte header to determine grid dimensions and number of layers.
 * Skips to layer 4 (fuel model) and reads 16-bit integer codes.
 * Reverses row order for domain consistency.
 *
 * Only rank 0 reads; broadcasts via MPI_Bcast (count first, then data).
 *
 * @param[in]  filename     Path to LCP file
 * @param[in]  fire_nx      Fire grid resolution in x (must match file width)
 * @param[in]  fire_ny      Fire grid resolution in y (must match file height)
 * @param[out] fuel_codes   Vector of integer fuel codes (size fire_nx * fire_ny)
 * @return true if file read successfully, false otherwise
 */
inline bool read_lcp_fuel_map(const std::string& filename,
                              int fire_nx, int fire_ny,
                              std::vector<int>& fuel_codes)
{
    fuel_codes.clear();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.good()) {
            return false;
        }

        // Read 256-byte header
        char header[256];
        file.read(header, 256);
        if (!file.good()) {
            return false;
        }

        // Parse header to get number of layers (at offset 52, 2-byte integer)
        // Standard LCP header structure (from FARSITE documentation):
        //   Offset 0-3:   "LOLA" or similar magic
        //   Offset 4-7:   Version
        //   Offset 8-11:  Latitude
        //   Offset 12-15: Longitude
        //   ...
        //   Offset 52-53: Number of layers
        int16_t n_layers = 0;
        std::memcpy(&n_layers, header + 52, sizeof(int16_t));

        // Ensure at least 5 layers (0-4 for fuel model at index 4)
        if (n_layers < 5) {
            return false;
        }

        // Grid dimensions are typically at offset 10-11 (nrows) and 12-13 (ncols)
        // But infer from fire_nx, fire_ny for safety
        int expected_size = fire_nx * fire_ny;

        // Each layer is ncols * nrows 16-bit integers
        int layer_size_bytes = expected_size * static_cast<int>(sizeof(int16_t));

        // Skip to layer 4 (fuel model layer)
        // Offset: 256 (header) + 4 * layer_size_bytes
        file.seekg(256 + 4 * layer_size_bytes, std::ios::beg);
        if (!file.good()) {
            return false;
        }

        // Read fuel model data (16-bit integers, row 0 = northernmost)
        std::vector<int16_t> temp_i16(expected_size);
        file.read(reinterpret_cast<char*>(temp_i16.data()), layer_size_bytes);
        if (!file.good()) {
            return false;
        }

        file.close();

        // Convert to int and reverse row order
        std::vector<int> temp_data(expected_size);
        for (int j = 0; j < fire_ny; ++j) {
            for (int i = 0; i < fire_nx; ++i) {
                int file_row = fire_ny - 1 - j;  // Reverse row order
                int src_idx = file_row * fire_nx + i;
                int dst_idx = j * fire_nx + i;
                temp_data[dst_idx] = static_cast<int>(temp_i16[src_idx]);
            }
        }

        fuel_codes = temp_data;
    }

    // Broadcast count and data via MPI
    int count = static_cast<int>(fuel_codes.size());
    amrex::ParallelDescriptor::Bcast(&count, 1);

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        fuel_codes.resize(count);
    }

    if (count > 0) {
        amrex::ParallelDescriptor::Bcast(fuel_codes.data(), count);
    }

    return true;
}

#endif  // ERF_LCP_READER_H
