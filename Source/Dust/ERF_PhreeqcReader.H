/**
 * @file ERF_PhreeqcReader.H
 * @brief Reader for PHREEQC geochemical output files for the ERF-Dust layer.
 *
 * PHREEQC (Parkhurst & Appelo 2013) computes geochemical speciation and
 * mineral equilibria offline. This reader ingests PHREEQC output files at
 * prescribed intervals and updates dust grid MultiFabs. No runtime PHREEQC
 * calls are made; coupling is one-way and file-based.
 *
 * Timescale justification: geochemical processes (crust formation,
 * efflorescence, salt precipitation) evolve over days to weeks. The
 * atmospheric timestep is seconds. File-based periodic updates are sufficient.
 *
 * Supported file formats:
 *   CSV (.csv): comma-separated, first row = header with column names,
 *     subsequent rows = one row per dust grid cell in row-major order
 *     (j=0..ny-1, i=0..nx-1 within each row).
 *   NetCDF (.nc): requires ERF_ENABLE_NETCDF=ON.
 *     Variable name specified by caller. Coordinate variables used for
 *     bilinear interpolation if grid does not match dust grid exactly.
 *
 * PHREEQC output fields read by this module:
 *   crust_index      [0,1]     — mineral crust strength index.
 *   silt_fraction    [0,1]     — updated silt mass fraction.
 *   efflorescence    [0,1]     — salt efflorescence fraction.
 *   suppression_mod  [0,1]     — suppression agent binding modifier
 *                                from PHREEQC salt hygroscopicity.
 *   metal_as_bin0..N [kg/kg]   — toxic metal (As, Pb, Cd, Se) mass
 *                                fraction per size bin per cell.
 *
 * Threshold friction velocity update:
 *   u*_t_new = u*_t_base * (1 + alpha_c * crust_index)
 *                        * (1 + alpha_e * efflorescence)
 *   where alpha_c and alpha_e are crust and efflorescence enhancement
 *   coefficients from DustParams. Higher crust increases u*_t (harder to emit).
 *   Reference: Marticorena & Bergametti (1995), J. Geophys. Res., 100,
 *     16415-16430. https://doi.org/10.1029/95JD00690
 *
 * I/O pattern follows ERF_FuelMap.H and ERF_FireTerrainReader.cpp:
 *   rank-0 reads -> ParallelDescriptor::Bcast -> Gpu::copy -> ParallelFor.
 *
 * References:
 *   Parkhurst, D.L., & Appelo, C.A.J. (2013). PHREEQC version 3.
 *     USGS Techniques and Methods, book 6, chap. A43.
 *     https://pubs.usgs.gov/tm/06/a43/
 *   Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao, Y., & Lu, H. (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 *   Bagnold, R.A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 */

#ifndef ERF_PHREEQC_READER_H
#define ERF_PHREEQC_READER_H

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>
#include <string>

// Forward declarations
struct DustParams;
struct DustGrid;

/**
 * @namespace PhreeqcDustConst
 * @brief Physical constants and coefficients for PHREEQC geochemical coupling.
 *
 * Modeled on FuelMoistureConst in ERF_FuelMoisture.H.
 * Reference: Marticorena & Bergametti (1995), J. Geophys. Res., 100, 16415.
 *   https://doi.org/10.1029/95JD00690
 */
namespace PhreeqcDustConst {
    /// Crust strength enhancement coefficient for u*_t [dimensionless].
    /// Higher crust_index (fully crusted, protected soil) increases u*_t (harder to emit).
    /// Lower crust_index (bare, uncrusted soil) decreases u*_t (easier to emit).
    /// Formula: f_chem = (1 + ALPHA_CRUST * crust_index)
    /// Reference: Marticorena & Bergametti (1995), https://doi.org/10.1029/95JD00690
    constexpr amrex::Real ALPHA_CRUST       = 0.5;

    /// Efflorescence enhancement coefficient for u*_t [dimensionless].
    /// Higher efflorescence increases u*_t (protective effect).
    /// Formula: f_efflor = (1 + ALPHA_EFFLOR * efflorescence)
    constexpr amrex::Real ALPHA_EFFLOR      = 0.3;

    /// Minimum allowed u*_t after crust/efflorescence modulation [m/s].
    constexpr amrex::Real USTAR_T_MIN       = 0.05;

    /// Suppression agent binding factor reduction per unit salt hygroscopicity.
    constexpr amrex::Real SUPP_HYGRO_FACTOR = 0.2;

    /// Nodata fill value for CSV files when a cell entry is missing [-9999].
    constexpr amrex::Real NODATA_CSV        = -9999.0;
}

/**
 * @brief Read a PHREEQC CSV output file and scatter values to a MultiFab.
 *
 * CSV format: first row is a header. Subsequent rows contain one value per
 * dust grid cell in row-major order: j=0 (southernmost row), i=0..nx-1,
 * then j=1, etc. The column named col_name is extracted. Cells with value
 * PhreeqcDustConst::NODATA_CSV are set to nodata_fill.
 *
 * I/O pattern follows ERF_FuelMap.H: rank-0 reads, Bcast, GPU fill.
 *
 * @param[out] mf           MultiFab to fill (1 component, on DustGrid).
 * @param[in]  dg           DustGrid.
 * @param[in]  filename     Path to CSV file. Returns false if empty.
 * @param[in]  col_name     Column name to extract from CSV header.
 * @param[in]  nodata_fill  Value for nodata cells.
 * @return true on success.
 */
bool read_phreeqc_csv(amrex::MultiFab&   mf,
                      const DustGrid&    dg,
                      const std::string& filename,
                      const std::string& col_name,
                      amrex::Real        nodata_fill);

/**
 * @brief Read a PHREEQC NetCDF output file and scatter values to a MultiFab.
 *
 * Requires ERF_ENABLE_NETCDF=ON. Calls amrex::Abort if built without NetCDF.
 * Uses ERF_USE_NETCDF macro (set by CMake), matching ERF.cpp convention.
 */
bool read_phreeqc_netcdf(amrex::MultiFab&   mf,
                         const DustGrid&    dg,
                         const std::string& filename,
                         const std::string& varname,
                         amrex::Real        nodata_fill);

/**
 * @brief Update u*_t from crust index and efflorescence using the
 * Marticorena & Bergametti (1995) enhancement formula:
 *
 *   u*_t_new = u*_t_base * (1 + alpha_c * crust_index)
 *                        * (1 + alpha_e * efflorescence)
 *
 * clamped to [PhreeqcDustConst::USTAR_T_MIN, ...].
 * GPU kernel over dust grid cells. Called after crust_index and
 * efflorescence MultiFabs have been updated.
 * Higher crust and efflorescence both increase u*_t (harder to emit).
 *
 * @param[in,out] dust_ustar_t    Threshold friction velocity [m/s].
 * @param[in]     dust_ustar_base Base u*_t before chemistry [m/s].
 * @param[in]     dust_crust      Crust index MultiFab [0,1].
 * @param[in]     dust_efflor     Efflorescence fraction MultiFab [0,1].
 */
void update_ustar_t_from_chemistry(amrex::MultiFab&       dust_ustar_t,
                                   const amrex::MultiFab& dust_ustar_base,
                                   const amrex::MultiFab& dust_crust,
                                   const amrex::MultiFab& dust_efflor);

/**
 * @brief Update all dust MultiFabs from PHREEQC output files.
 *
 * Selects reader by file extension (.nc -> NetCDF, else CSV).
 * Called from DustLayer::advance when phreeqc_update_interval_s has elapsed.
 * If phreeqc_output_file is empty, returns immediately without modification.
 *
 * @param[in,out] dust_ustar_t     Threshold friction velocity [m/s].
 * @param[in]     dust_ustar_base  Base u*_t MultiFab (retained from init).
 * @param[in,out] dust_crust_index Crust strength index [0,1].
 * @param[in,out] dust_silt_frac   Silt mass fraction [-].
 * @param[in,out] dust_efflor      Efflorescence fraction [0,1].
 * @param[in,out] dust_suppression Suppression coverage [-].
 * @param[in,out] dust_emission    Emission flux MultiFab (ncomp=n_size_bins).
 *                                 Component 0 updated with As mass fraction.
 * @param[in]     dg               DustGrid.
 * @param[in]     params           DustParams.
 */
void update_dust_from_phreeqc(amrex::MultiFab&       dust_ustar_t,
                              const amrex::MultiFab& dust_ustar_base,
                              amrex::MultiFab&       dust_crust_index,
                              amrex::MultiFab&       dust_silt_frac,
                              amrex::MultiFab&       dust_efflor,
                              amrex::MultiFab&       dust_suppression,
                              amrex::MultiFab&       dust_emission,
                              const DustGrid&        dg,
                              const DustParams&      params);

#endif // ERF_USE_DUST

#endif // ERF_PHREEQC_READER_H
