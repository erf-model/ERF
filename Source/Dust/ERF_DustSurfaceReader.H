#ifndef ERF_DUST_SURFACE_READER_H
#define ERF_DUST_SURFACE_READER_H

#ifdef ERF_USE_DUST

/**
 * @file ERF_DustSurfaceReader.H
 * @brief Surface property map reader for the ERF-Dust 2D emission layer.
 *
 * Reads raster files and populates MultiFabs on the 2D dust grid via
 * bilinear interpolation. I/O pattern follows ERF_FireTerrainReader.cpp
 * and ERF_FuelMap.H: rank-0 reads, MPI broadcast, GPU interpolation kernel.
 *
 * Supported formats:
 *   ESRI ASCII grid (.asc): 6-line header, data rows, row 0 = northernmost.
 *     https://desktop.arcgis.com/en/arcmap/latest/manage-data/
 *     raster-and-images/esri-ascii-raster-format.htm
 *   NetCDF (.nc): requires ERF_ENABLE_NETCDF=ON.
 *     https://docs.unidata.ucar.edu/netcdf-c/current/
 *
 * References:
 *   Marticorena & Bergametti (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao & Lu (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 *   U.S. EPA AP-42, Chapter 13.2.2 (Unpaved Roads).
 *     https://www.epa.gov/air-emissions-factors-and-quantification/
 *     ap-42-compilation-air-emission-factors
 */

#include <AMReX_MultiFab.H>
#include <string>

// Forward declarations
struct DustParams;
struct DustGrid;

/**
 * @brief Read ESRI ASCII raster and fill MultiFab via bilinear interpolation.
 *
 * Follows ERF_FireTerrainReader.cpp pattern:
 *   rank-0 reads 6-line header + data -> Bcast(ncols,nrows) -> Bcast(data)
 *   -> Gpu::copy(hostToDevice) -> ParallelFor bilinear interp over tiles.
 * Row 0 in file = northernmost; reversed on read to match domain convention.
 * Nodata cells set to nodata_fill.
 *
 * @param[out] mf           MultiFab to fill (1 component, on DustGrid).
 * @param[in]  dg           DustGrid geometry for interpolation.
 * @param[in]  filename     Path to .asc file. Returns false if empty.
 * @param[in]  nodata_fill  Value for nodata cells.
 * @return true on success.
 */
bool read_ascii_surface_map(amrex::MultiFab&   mf,
                            const DustGrid&    dg,
                            const std::string& filename,
                            amrex::Real        nodata_fill);

/**
 * @brief Read NetCDF 2D variable and fill MultiFab via bilinear interpolation.
 *
 * Requires ERF_ENABLE_NETCDF=ON. Calls amrex::Abort if built without NetCDF.
 * Uses ERF_USE_NETCDF macro (set by CMake), matching ERF.cpp convention.
 */
bool read_netcdf_surface_map(amrex::MultiFab&   mf,
                             const DustGrid&    dg,
                             const std::string& filename,
                             const std::string& varname,
                             amrex::Real        nodata_fill);

/**
 * @brief Populate all dust surface MultiFabs from file paths in DustParams.
 *
 * Called from DustLayer::initialize after MultiFab allocation. Selects
 * reader by extension (.nc -> NetCDF, else ASCII). If a path is empty
 * the MultiFab retains the default from DustLayer::initialize.
 */
void populate_dust_surface_maps(amrex::MultiFab&  dust_soil_type,
                                amrex::MultiFab&  dust_silt_fraction,
                                amrex::MultiFab&  dust_crust_index,
                                amrex::MultiFab&  dust_moisture_flag,
                                amrex::MultiFab&  dust_suppression,
                                const DustGrid&   dg,
                                const DustParams& params);

#endif // ERF_USE_DUST

#endif // ERF_DUST_SURFACE_READER_H
