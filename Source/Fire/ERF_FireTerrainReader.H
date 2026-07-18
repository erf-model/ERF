#ifndef ERF_FIRE_TERRAIN_READER_H
#define ERF_FIRE_TERRAIN_READER_H

#include <AMReX_MultiFab.H>
#include <ERF_FireGrid.H>
#include <string>

/**
 * @file ERF_FireTerrainReader.H
 * @brief Self-contained terrain reader for fire grid
 *
 * Reads terrain data from a standard ERF text file format and interpolates
 * onto the fire grid using bilinear interpolation.
 *
 * File format (line-based):
 *   Line 1: nx (integer)
 *   Line 2: ny (integer)
 *   Next nx lines: x-coordinates [m]
 *   Next ny lines: y-coordinates [m]
 *   Next nx*ny lines: z-values stored contiguous in y, i.e., z[ix*ny + iy]
 */

/**
 * @brief Read terrain data onto fire grid with bilinear interpolation
 *
 * Reads ERF standard terrain text file format:
 * 1. IO rank reads nx, ny from first two lines
 * 2. Reads nx x-coordinates, ny y-coordinates, nx*ny z-values
 * 3. Broadcasts nx, ny, and arrays to all MPI ranks
 * 4. Copies arrays to GPU with Gpu::DeviceVector
 * 5. Bilinearly interpolates terrain heights at fire-grid nodes
 * 6. Fills z_fire_nd MultiFab with interpolated values
 *
 * @param[out] z_fire_nd   Terrain height at fire-grid nodes (1 component, 1 ghost cell)
 *                         On fire grid, indexed by fire grid cell indices
 * @param[in]  fg          Fire grid structure
 * @param[in]  fname       File path (empty string means no file, returns false)
 *
 * @return true if file was successfully read and interpolation completed,
 *         false if fname is empty or file cannot be opened
 *
 * Ghost cells are handled by clamping indices to valid domain bounds before
 * computing physical positions.
 *
 * @note This is a host-side I/O operation followed by GPU-side interpolation.
 *       The function broadcasts to all MPI ranks so all can compute in parallel.
 */
bool read_terrain_onto_fire_grid(
    amrex::MultiFab&   z_fire_nd,
    const FireGrid&    fg,
    const std::string& fname);

#endif
