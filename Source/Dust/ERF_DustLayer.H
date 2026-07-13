/**
 * @file ERF_DustLayer.H
 * @brief Container class for the ERF-Dust 2D surface emission layer.
 *
 * DustLayer owns all dust-state MultiFabs on the 2D dust grid and
 * provides the initialize() and advance() interface to the ERF solver.
 *
 * In Phase 2, initialize() allocates MultiFabs and sets them to their
 * default values. advance() is a stub that returns without computation.
 * Subsequent phases replace the stub advance() with physics.
 *
 * The dust grid is a 2D slab (k-extent = 1) at k=0 of the atmospheric
 * domain, optionally refined by DustParams::grid_ratio in x and y.
 *
 * MultiFab layout (all on DustGrid, grid_ratio*nx x grid_ratio*ny x 1):
 *
 * | Name                 | ncomp           | Description                             |
 * |----------------------|-----------------|-----------------------------------------|
 * | dust_ustar_t         | 1               | Threshold friction velocity [m/s]       |
 * | dust_soil_type       | 1               | Soil type index [integer stored as Real]|
 * | dust_silt_fraction   | 1               | Surface silt mass fraction [-]          |
 * | dust_crust_index     | 1               | Mineral crust strength index [0,1]      |
 * | dust_moisture_flag   | 1               | Surface moisture inhibition flag [0,1]  |
 * | dust_suppression     | 1               | Suppression agent coverage fraction [-] |
 * | dust_emission_flux   | n_size_bins     | Emission flux per size bin [kg/m^2/s]  |
 * | dust_efflor          | 1               | Efflorescence fraction [0,1]            |
 * | dust_ustar_base      | 1               | Base u*_t before chemistry [m/s]        |
 * | dust_ustar_in        | 1               | Friction velocity u* [m/s]              |
 *
 * dust_soil_type stores an integer code as amrex::Real. Codes follow the
 * STATSGO classification for natural soils (1-16) and reserve codes
 * >= 100 for mine-specific surface types (defined in Phase 3).
 *
 * References:
 *   Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao, Y., & Lu, H. (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 *   Bagnold, R. A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 */

#ifndef ERF_DUST_LAYER_H
#define ERF_DUST_LAYER_H

#include <ERF_DustGrid.H>
#include <ERF_DustParams.H>
#include <memory>
#include <AMReX_MultiFab.H>

#ifdef ERF_USE_DUST
#include <ERF_DustBlastSchedule.H>
#endif

// Forward declarations
class SurfaceLayer;
class ERF;

class DustLayer
{
public:
    DustLayer()  = default;
    ~DustLayer() = default;

    /**
     * @brief Initialize the dust layer after SurfaceLayer is ready.
     *
     * Performs prerequisite verification, creates the 2D dust grid,
     * allocates all MultiFabs listed in the class-level table, and
     * fills them with default values from dust_params.
     *
     * Called once from ERF::InitData_post() when erf.dust.enable = true.
     *
     * @param[in] erf              ERF instance (public AmrCore accessors only).
     * @param[in] surface_layer    Pointer to the initialized SurfaceLayer.
     * @param[in] dust_params      Dust configuration parameters.
     */
    void initialize(const ERF&          erf,
                    const SurfaceLayer* surface_layer,
                    const DustParams&   dust_params);

    /**
     * @brief Advance the dust layer by one atmospheric timestep.
     *
     * In Phase 2 this function is a stub and returns without computation.
     * Phases 5-13 replace it with emission, coupling, and deposition physics.
     *
     * @param[in] dt                 Atmospheric timestep [s].
     * @param[in] dust_params        Dust configuration parameters.
     * @param[in] surface_layer      Pointer to SurfaceLayer (nullptr = placeholder fallback).
     * @param[in] xvel_mf            3D u velocity field (nullptr = placeholder fallback).
     * @param[in] yvel_mf            3D v velocity field (nullptr = placeholder fallback).
     * @param[in] z_phys_cc_mf       3D terrain height field (nullptr = placeholder fallback).
     * @param[in] nz                 Number of vertical levels (0 = placeholder fallback).
     */
    void advance(amrex::Real             dt,
                 const DustParams&       dust_params,
                 const SurfaceLayer*     surface_layer = nullptr,
                 const amrex::MultiFab*  xvel_mf = nullptr,
                 const amrex::MultiFab*  yvel_mf = nullptr,
                 const amrex::MultiFab*  z_phys_cc_mf = nullptr,
                 int                     nz = 0);

    // Const accessors for plotfile and diagnostic output
    const amrex::MultiFab* get_ustar_t()         const { return dust_ustar_t.get(); }
    const amrex::MultiFab* get_soil_type()        const { return dust_soil_type.get(); }
    const amrex::MultiFab* get_silt_fraction()    const { return dust_silt_fraction.get(); }
    const amrex::MultiFab* get_crust_index()      const { return dust_crust_index.get(); }
    const amrex::MultiFab* get_moisture_flag()    const { return dust_moisture_flag.get(); }
    const amrex::MultiFab* get_suppression()      const { return dust_suppression.get(); }
    const amrex::MultiFab* get_efflor()           const { return dust_efflor.get(); }
    const amrex::MultiFab* get_emission_flux()    const { return dust_emission_flux.get(); }
    const amrex::MultiFab* get_surf_moist()       const { return dust_surf_moist.get(); }
    const amrex::MultiFab* get_ustar_in()         const { return dust_ustar_in.get(); }
    const amrex::MultiFab* get_retreat_flag()     const { return dust_retreat_flag.get(); }
    const amrex::MultiFab* get_wind_ref()         const { return dust_wind_ref.get(); }
    const amrex::MultiFab* get_pblh()             const { return dust_pblh.get(); }
    const amrex::MultiFab* get_tsfc()             const { return dust_tsfc.get(); }
    const DustGrid&        get_dust_grid()        const { return m_dg; }
    const amrex::Geometry& get_dust_geom()        const { return m_dg.geom; }
    const DustParams&      get_params()           const { return m_params; }

private:
    DustGrid                              m_dg;

    std::unique_ptr<amrex::MultiFab>      dust_ustar_t;         ///< Threshold u* [m/s]
    std::unique_ptr<amrex::MultiFab>      dust_soil_type;       ///< Soil type index (as Real)
    std::unique_ptr<amrex::MultiFab>      dust_silt_fraction;   ///< Silt mass fraction [-]
    std::unique_ptr<amrex::MultiFab>      dust_crust_index;     ///< Crust strength index [0,1]
    std::unique_ptr<amrex::MultiFab>      dust_moisture_flag;   ///< Moisture inhibition flag [0,1]
    std::unique_ptr<amrex::MultiFab>      dust_suppression;     ///< Suppression coverage fraction [-]
    std::unique_ptr<amrex::MultiFab>      dust_emission_flux;   ///< Emission flux per bin [kg/m^2/s]
    std::unique_ptr<amrex::MultiFab>      dust_efflor;          ///< Efflorescence fraction [0,1]. Updated by PHREEQC reader.
    std::unique_ptr<amrex::MultiFab>      dust_ustar_base;      ///< Base threshold friction velocity before chemistry [m/s]. Set at init, not modified.
    std::unique_ptr<amrex::MultiFab>      dust_surf_moist;      ///< Surface wetness [0,1]. Set to 0 at init. Driven by ERF surface layer in Phase 9.
    std::unique_ptr<amrex::MultiFab>      dust_ustar_in;        ///< Friction velocity u* [m/s].
                                                                  ///<   Phase 6: uniform value from DustParams.
                                                                  ///<   Phase 9: extracted from MRF surface layer.
    std::unique_ptr<amrex::MultiFab>      dust_retreat_flag;    ///< Re-treatment flag [0=ok, 1=retreat needed].
                                                                  ///<   Set each timestep by advance_dust_suppression.
                                                                  ///<   Phase 17 exports this for operational use.

    // Phase 9: atmospheric fields extracted each timestep
    std::unique_ptr<amrex::MultiFab>      dust_wind_ref;        ///< u,v at zref [m/s] (ncomp=2)
    std::unique_ptr<amrex::MultiFab>      dust_pblh;            ///< MRF PBL height [m]
    std::unique_ptr<amrex::MultiFab>      dust_tsfc;            ///< Surface temperature [K]

    DustParams                            m_params;
    int                                   m_step = 0;           ///< Timestep counter
    amrex::Real                           m_time = 0.0;         ///< Simulation time [s]
    amrex::Real                           m_last_phreeqc_update = -1.0e30;  ///< Time of last PHREEQC update [s]
#ifdef ERF_USE_DUST
    BlastSchedule                         m_blast_schedule;     ///< Loaded at init if blast_schedule_file is non-empty.
    bool                                  m_has_blast_schedule = false; ///< True if blast schedule was loaded.
#endif
};

#endif // ERF_DUST_LAYER_H
