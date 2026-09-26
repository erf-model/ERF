#ifndef ERF_SLM_H
#define ERF_SLM_H

#include <string>
#include <vector>
#include <memory>
#include <ctime>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_TableData.H>
#include <AMReX_ParmParse.H>

#include <ERF_NullSurf.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_LSM_Types.H>

#include <ERF_NCInterface.H>

namespace LsmVar_SLM {
   enum {
      // independent variables
      tsurf = 0,           // surface temperature
      ustar,               // ustar for MOST
      tstar,               // tstar for MOST
      qstar,               // qstar for MOST
      // inputs from ERF
      tref,            // temperature at reference level
      qref,            // specific humidity at ref
      pref,            // pressure at reference level
      dref,            // density at reference level
      uref,            // zonal wind speed at ref
      vref,            // meridional wind speed at ref
      precipref,       // precipitation rate at ref
      swdsvisxyref,    // downwelling SW radiation, direct visible flux
      swdsnirxyref,    // downwelling SW radiation, direct near-IR flux
      swdsvisdxyref,   // downwelling SW radiation, diffuse visible flux
      swdsnirdxyref,   // downwelling SW radiation, diffuse near-IR flux
      lwref,           // downwelling LW radiation
      coszrsxy,        // zenith angle
      // model variables
      tv,              // vegetation temperature
      mv,              // vegetation moisture storage
      soilt,           // soil temperature
      soilw,           // soil wetness
      // soil parameters
      sand,            // sand percentage
      clay,            // clay percentage
      s_depth,         // thickness of each soil layer in meters
      node_z,          // depth of the center of each soil layer
      interface_z,     // depth of the soil layer interface
      soilw_nudge,     // nudging of soilw
      soilt_nudge,     // nudging of soilt
      sw_over_sat,     // soilw over saturation
      soil_relax_hgt,  // depth dependent soil relaxation function
      soilt_obs,
      soilw_obs,
       // vegetation
       soil_transp_frac,// normalized fraction of transpiration extracted from each soil layer
       w_s_FC,          // wetness at field capacity
      w_s_WP,          // wetness at wilting point
      // soil parameters calculated from sand/clay
      sst_capa,        // heat capacity for soil solids
      sst_cond,        // heat conductivity of soil solids
      poro_soil,       // saturation moisture content, porosity
      theta_FC,        // volumetric moisture content at field capacity
      theta_WP,        // volumetric moisture content at wilting point
      m_pot_sat,       // moisture potential at saturation
      Bconst,          // B constant
      ks,              // hydraulic conductivity at saturation [m/s]
      // 2D inputs populated from WRFInput
      lai,             // leaf area index (2D)
      vegtype,         // vegetation type (2D)
      soiltype,        // soil type (2D)
      veg_frac,        // vegetation fraction
      veg_frac_min,    // yearly min vegetation fraction
      veg_frac_max,    // yearly max vegetation fraction
      // model outputs
      flbu,            // surface u-momentum
      flbv,            // surface v-momentum
      flbq,            // surface vapor flux (latent heat flux)
      flbt,            // surface heat flux
      prsfc,           // precipitation rate reaching soil
      emis_sfc,        // IR emissivity of surface (veg+soil)
      alb_nir_sfc,     // NIR albedo of surface (veg+soil)
      alb_vis_sfc,     // VIS albedo of surface (veg+soil)
      alb_nir_sfc_diff,// NIR albedo of surface (veg+soil), diffuse beam
      alb_vis_sfc_diff,// VIS albedo of surface (veg+soil), diffuse beam
      NumVars
  };
}

namespace LsmFlux_SLM {
    enum {
      // 2D SurfaceLayer flux vars
      t_flux = 0,
      q_flux,
      tau13,
      tau23,
      olen,         // Obukhov length for MOST/PBL
      NumVars
    };
}

namespace SLM_NetRad {
   enum {
     net_swup1 = 0,
     net_swup2,
     net_swdn1,
     net_swdn2,
     net_sw1,
     net_sw2,
     tir1,
     tir2,
     net_lwup1,
     net_lwup2,
     net_lwdn1,
     net_lwdn2,
     net_lw1,
     net_lw2,
     net_rad1,
     net_rad2,
     NumVars
   };
}

namespace SLM_Diag {
   enum {
     precip = 0,
     precip_in,
     precip_sfc, // saved in prsfc too, but 3D
     run_off_sfc,
     drain,
     drain_flux,
     evapo_wet,
     grflux,
     NumVars
    };
}

// refers to device variables used in soil temperature and moisture update
namespace SLM_DST {
    enum {
        st_cond,
        st_capa,
        st_eff_cond,
        alpha,
        beta,
        NumVars
    };
}
namespace SLM_DSW {
    enum {
        sw_wgt,
        sh_eff_cond,
        sh_eff_vel,
        alpha,
        beta,
        sdepth_mm,
        NumVars
    };
}

using SLMParameterTable = std::unordered_map<std::string, amrex::Vector<std::pair<std::string, amrex::Vector<amrex::Real>>>>;

/**
 * The Simple Land-surface Model (SLM) supplies lower-boundary
 * sensible-heat, latent-heat, and momentum fluxes for land cells.
 *
 * The original formulation and evaluation of SLM are documented in Lee and
 * Khairoutdinov (2015), "A Simplified Land Model (SLM) for use in
 * cloud-resolving models: Formulation and evaluation," Journal of Advances
 * in Modeling Earth Systems, 7, 1368--1392,
 * https://doi.org/10.1002/2014MS000419.
 *
 * SLM includes selected formulations adapted from NoahMP, updates inherited
 * from gSAM-SLM, and SLM-specific corrections to the surface energy and water
 * budgets. It remains a simplified land model and is not a complete
 * implementation of NoahMP.
 */
class SLM : public NullSurf {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    using NullSurf::Init;

    // Constructor
    SLM () {}

    // Destructor
    virtual ~SLM () = default;

    // Set thermo and grid properties
    void
    Define (SolverChoice& sc) override
    {
        // SLM exchanges water vapor with the atmosphere unconditionally: Copy_State_to_Lsm
        // reads RhoQ1_comp out of the state and Copy_Lsm_to_State writes it back.  Without a
        // moisture model that component is not part of the state, so those accesses run off
        // the end of the conserved-variable FAB.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            sc.moisture_type != MoistureType::None,
            "The SLM land surface model requires a moisture model: it exchanges water vapor "
            "with the atmospheric state, which has no RhoQ1 component when "
            "erf.moisture_model = None.  Set erf.moisture_model, "
            "or choose a different erf.land_surface_model.");

        use_terrain = sc.terrain_type != TerrainType::None;

        // Check if initializing from WRFInput
        if (sc.init_type == InitType::WRFInput) {
            amrex::Print() << " SLM is initializing from WRFInput!" << std::endl;
            use_wrfinput = true;
        } else {
            use_wrfinput = false;
        }
    }

    // Initialize data structures
    void
    Init (const int& /*lev*/,
          const amrex::MultiFab& cons_in,
          const amrex::MultiFab& u_in,
          const amrex::MultiFab& v_in,
          const amrex::Geometry& geom,
          const amrex::Geometry& geom0,
          amrex::Vector<amrex::BCRec>& domain_bcs_type,
          amrex::IntVect& refRatio,
          const amrex::Real& dt,
          std::unique_ptr<amrex::MultiFab>& z_phys_nd_in,
          amrex::Vector<amrex::Vector<std::string>>& nc_init_file) override;

    /**
     * Rebuilds SLM arrays on a remade atmospheric grid and transfers state.
     */
    void
    Lsm_Regrid_Level (const int& lev,
                      const amrex::MultiFab& cons_in,
                      const amrex::MultiFab& u_in,
                      const amrex::MultiFab& v_in,
                      const amrex::Geometry& geom,
                      const amrex::Geometry& geom0,
                      amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      amrex::IntVect& refRatio,
                      const amrex::Real& dt,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd_in,
                      amrex::Vector<amrex::Vector<std::string>>& nc_init_file,
                      const amrex::Geometry& source_geom,
                      const amrex::Geometry& source_lsm_geom,
                      const amrex::IntVect& source_ref_ratio,
                      const bool has_source) override;

    // Wrapper to do all the updating
    void
    Advance (const double& dt) override
    {
        Advance(static_cast<amrex::Real>(dt), time + static_cast<amrex::Real>(dt), start_time);
    }

    void
    Advance (const amrex::Real& dt,
             const amrex::Real& t_new,
             const amrex::Real& t_start) override
    {
        m_updated = false;
        m_dt = dt;
        time = t_new;
        start_time = t_start;

        if (first_step)
        {
            first_step = false;

            if (use_wrfinput && !wrfinput_initialized) {
                this->init_wrfinput_vars();
            }

            initialize_zrefxy();
        }

        this->AdvanceSLM();
        m_updated = true;

    }

    // Advance the lsm state vars
    void
    AdvanceSLM ();

    // Get state vars from lsm class
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_data_size && varIdx >= 0);
        int lsmIdx = LsmDataMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < LsmVar_SLM::NumVars && lsmIdx>=0);
        return lsm_fab_vars[lsmIdx].get();
    }

    // Get flux vars from lsm class
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_flux_size && varIdx >= 0);
        int lsmIdx = LsmFluxMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < LsmFlux_SLM::NumVars && lsmIdx>=0);
        return lsm_fab_flux[lsmIdx].get();
    }

    // Get lsm geometry
    amrex::Geometry
    Lsm_Geom ( ) override { return m_lsm_geom; }

    /**
     * Returns the atmospheric geometry used by SLM coupling fields.
     */
    amrex::Geometry
    Lsm_Atmos_Geom () override { return m_geom; }

    /**
     * Registers an exposed source data field.
     */
    void
    Lsm_Set_Source_Data_Ptr (const int& varIdx,
                             amrex::MultiFab* source_data_ptr) override
    {
        AMREX_ALWAYS_ASSERT(varIdx >= 0 && varIdx < m_lsm_data_size);
        if (static_cast<int>(source_data.size()) != m_lsm_data_size) {
            source_data.resize(m_lsm_data_size, nullptr);
        }
        source_data[varIdx] = source_data_ptr;
    }

    /**
     * Registers an exposed source flux field.
     */
    void
    Lsm_Set_Source_Flux_Ptr (const int& varIdx,
                             amrex::MultiFab* source_flux_ptr) override
    {
        AMREX_ALWAYS_ASSERT(varIdx >= 0 && varIdx < m_lsm_flux_size);
        if (static_cast<int>(source_flux.size()) != m_lsm_flux_size) {
            source_flux.resize(m_lsm_flux_size, nullptr);
        }
        source_flux[varIdx] = source_flux_ptr;
    }

    /**
     * Returns the number of private SLM state fields.
     */
    int
    Lsm_Unmapped_Size () override { return static_cast<int>(unmapped_fields.size()); }

    /**
     * Returns a private SLM state field.
     */
    amrex::MultiFab*
    Lsm_Unmapped_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx >= 0 && varIdx < static_cast<int>(unmapped_fields.size()));
        return lsm_fab_vars[unmapped_fields[varIdx]].get();
    }

    /**
     * Registers a private source SLM state field.
     */
    void
    Lsm_Set_Source_Unmapped_Ptr (const int& varIdx,
                                  amrex::MultiFab* source_data_ptr) override
    {
        AMREX_ALWAYS_ASSERT(varIdx >= 0 && varIdx < static_cast<int>(unmapped_fields.size()));
        if (source_unmapped.size() != unmapped_fields.size()) {
            source_unmapped.resize(unmapped_fields.size(), nullptr);
        }
        source_unmapped[varIdx] = source_data_ptr;
    }

    /**
     * Transfers exposed state from a source AMR level.
     */
    void
    Lsm_Interpolate_From_Source (const amrex::Geometry& source_geom,
                                 const amrex::Geometry& source_lsm_geom,
                                 const amrex::IntVect& ref_ratio,
                                 LSMTransferMode mode) override;

    // Get number of vars lsm class contains
    int
    Lsm_Data_Size () override { return SLM::m_lsm_data_size; }

    // Get number of vars lsm class contains
    int
    Lsm_Flux_Size () override { return SLM::m_lsm_flux_size; }

    // Get variable names
    std::string
    Lsm_DataName (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_data_size && varIdx >= 0);
        return LsmDataName[varIdx];
    }

    int
    Lsm_DataIndex (const std::string &varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<SLM::m_lsm_data_size; ++idx) {
            if (lc_varname == amrex::toLower(Lsm_DataName(idx))) {
                varIdx = idx;
            }
        }
        AMREX_ASSERT(varIdx != -1); // no name found
        return varIdx;
    }

    // Get flux variable names
    std::string
    Lsm_FluxName (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_flux_size && varIdx >= 0);
        return LsmFluxName[varIdx];
    }

    int
    Lsm_FluxIndex (const std::string &varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<SLM::m_lsm_flux_size; ++idx) {
            if (lc_varname == amrex::toLower(Lsm_FluxName(idx))) {
                varIdx = idx;
            }
        }
        AMREX_ALWAYS_ASSERT(varIdx != -1); // no name found
        return varIdx;
    }

    std::unordered_map<std::string, std::string>&
    Lsm_WRFInputNames() override {
        return SLM::wrfinput_map;
    }

    void
    Copy_State_to_Lsm (const amrex::MultiFab& cons_in, const amrex::MultiFab& u_in, const amrex::MultiFab& v_in) override;

    void
    Copy_Lsm_to_State (amrex::MultiFab& cons_in) override;

    void
    Update_Lsm_Vars (amrex::MultiFab& cons_in, amrex::MultiFab& u_in, amrex::MultiFab& v_in) override
    {
        this->Copy_State_to_Lsm(cons_in, u_in, v_in);
        if (first_step) {
            if (use_wrfinput && !wrfinput_initialized) {
                this->init_wrfinput_vars();
            } else if (!use_wrfinput) {
                this->init_slm_vars();
            }
        }
    }

    bool Get_Update_Status () override { return m_updated; }

    void
    Update_State_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_Lsm_to_State(cons_in);
    }

    void
    set_precip_input(const amrex::MultiFab* precip_in) override;

    void
    set_terrain_inputs(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& sst_in,
                       const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& lmask_in) override;

    int Lsm_Type() override
    {
        return LSM_Type::SLM;
    }

    void writeSLM_Data(const PlotFileType plotfile_type, const amrex::Real time, const std::string plot_prefix, const int level_step, const int lev, const int finest_lev, amrex::MultiFab &fab, amrex::Geometry &geom, amrex::Vector<std::string> &varnames);

    void init_from_inputs();

    void slm_init();

    void init_layer_depths();

    void init_soil_vars();

    void init_landtype();

    void init_soil_tw();

    void init_slm_vars();

    /**
     * Initializes local SLM fields from WRFInput data.
     */
    void init_wrfinput_vars();

    /**
     * Initializes diagnostics and radiation coupling fields from the SLM state.
     */
    void initialize_surface_outputs();

    /**
     * Rebuilds restart fields that can be derived from the checkpoint state.
     */
    void rebuild_restart_fields();

    SLMParameterTable ReadParameterFile(const std::string &filename, amrex::Vector<std::string> &veg_categories, amrex::Vector<std::string> &soil_categories);

    /**
     * Validates the selected Noah-MP parameter blocks before copying them to the device.
     */
    void validate_parameter_tables();

    void init_from_params();

    /**
     * Loads built-in Noah-MP radiation lookup values for idealized SLM cases.
     */
    void load_builtin_radiation_parameters();

    void UpdateLAI(const amrex::MFIter &mfi);

    void UpdateLAIParameters(const amrex::MFIter &mfi);

    // NOAHMP radiation functions
    void radiation_noahmp(const amrex::MFIter &mfi);

    void albedo_noahmp(const amrex::MFIter &mfi, int i, int j,
                       amrex::Real dt, amrex::Real cosz, amrex::Real snowh, amrex::Real fsno,
                       amrex::Real fwet, amrex::Real sneqvo, amrex::Real sneqv, amrex::Real qsnow,
                       amrex::Real tg, amrex::Real tv, amrex::Real elai, amrex::Real esai,
                       const amrex::Real* smc, int nsoil, amrex::Real fveg,
                       amrex::Real& albold, amrex::Real& tauss,
                       amrex::Real* albgrd, amrex::Real* albgri, amrex::Real* albd, amrex::Real* albi,
                       amrex::Real* fabd, amrex::Real* fabi, amrex::Real* ftdd, amrex::Real* ftid, amrex::Real* ftii,
                       amrex::Real& fsun, amrex::Real* frevd, amrex::Real* frevi,
                       amrex::Real* fregd, amrex::Real* fregi, amrex::Real& bgap, amrex::Real& wgap,
                       amrex::Real* albsnd, amrex::Real* albsni);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void surrad_noahmp(amrex::Real mpe, amrex::Real fsun, amrex::Real fsha, amrex::Real elai, amrex::Real vai,
                              amrex::Real laisun, amrex::Real laisha, const amrex::Real* solad, const amrex::Real* solai,
                              const amrex::Real* fabd, const amrex::Real* fabi, const amrex::Real* ftdd,
                              const amrex::Real* ftid, const amrex::Real* ftii, const amrex::Real* albgrd,
                              const amrex::Real* albgri, const amrex::Real* albd, const amrex::Real* albi,
                              const amrex::Real* frevd, const amrex::Real* frevi, const amrex::Real* fregd, const amrex::Real* fregi,
                              amrex::Real& parsun, amrex::Real& parsha, amrex::Real& sav, amrex::Real& sag,
                              amrex::Real& fsa, amrex::Real& fsr, amrex::Real& fsrv, amrex::Real& fsrg);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void snow_age_noahmp(amrex::Real dt, amrex::Real tg, amrex::Real sneqvo, amrex::Real sneqv,
                                amrex::Real tau0, amrex::Real grain_growth, amrex::Real extra_growth,
                                amrex::Real dirt_soot, amrex::Real swemx,
                                amrex::Real& tauss, amrex::Real& fage);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void snowalb_bats_noahmp(int nband, amrex::Real fsno, amrex::Real cosz, amrex::Real fage,
                                    amrex::Real bats_cosz, amrex::Real bats_vis_new, amrex::Real bats_nir_new,
                                    amrex::Real bats_vis_age, amrex::Real bats_nir_age,
                                    amrex::Real bats_vis_dir, amrex::Real bats_nir_dir,
                                    amrex::Real* albsnd, amrex::Real* albsni);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void groundalb_noahmp(int nsoil, int nband, int ice, int ist, amrex::Real fsno,
                                 const amrex::Real* smc, const amrex::Real* albsnd, const amrex::Real* albsni,
                                 amrex::Real cosz, amrex::Real tg,
                                 const amrex::Real* albsat, const amrex::Real* albdry, const amrex::Real* alblak,
                                 amrex::Real* albgrd, amrex::Real* albgri);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static void twostream_noahmp(int ib, int ic, int vegtyp, amrex::Real cosz, amrex::Real vai,
                                 amrex::Real fwet, amrex::Real t, const amrex::Real* albgrd, const amrex::Real* albgri,
                                 const amrex::Real* rho, const amrex::Real* tau, amrex::Real fveg, int ist,
                                 amrex::Real xl, amrex::Real omegas_param, amrex::Real betads, amrex::Real betais,
                                 int opt_rad, amrex::Real rc, amrex::Real hvt, amrex::Real hvb, amrex::Real den,
                                 amrex::Real* fab, amrex::Real* fre, amrex::Real* ftd, amrex::Real* fti,
                                 amrex::Real& gdir, amrex::Real* frev, amrex::Real* freg,
                                 amrex::Real& bgap, amrex::Real& wgap,
                                 amrex::Real& xl_out, amrex::Real& chil_out, amrex::Real& phi1_out, amrex::Real& phi2_out);

    void transfer_coeff(const amrex::MFIter &mfi);

    void resistances(const amrex::MFIter &mfi);

    void fluxes_canopy(const amrex::MFIter &mfi);

    void solve_ground_skin_temperature(const amrex::MFIter &mfi);

    void soil_water(const amrex::MFIter &mfi);

    void soil_temperature(const amrex::MFIter &mfi);

    void soil_nudging(const amrex::MFIter &mfi3d);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real fh_calc(const amrex::Real &t, const amrex::Real &mps,
                               const amrex::Real &sw, const amrex::Real &B);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real linear_interp(const amrex::Real t0, const amrex::Real t1,
                                     const amrex::Real t, const amrex::Real x,
                                     const amrex::Real y);

#ifdef ERF_USE_NETCDF
    void writeSLM_NetCDF(const amrex::MultiFab& mf, const amrex::Vector<std::string>& varnames,
                         const amrex::Real time, const std::string plot_prefix, const int level_step);

    void writeNCHeader(ncutils::NCFile &nc_file, const amrex::Geometry &geom);

    void writeMFtoNC(ncutils::NCFile &nc_file, const amrex::MultiFab* mf,
                     const std::string name, amrex::Real time, bool write_ghost = false);

    void writeMFVecToNC(ncutils::NCFile &nc_file,
                        const amrex::Vector<const amrex::MultiFab*> &mf_vec,
                        const amrex::Vector<std::string> &mf_names,
                        const amrex::Geometry& geom) const;
#endif

    void set_dt(const amrex::Real &dt) {
        m_dt = dt;
    }

    const amrex::Vector<const amrex::MultiFab*> &export_to_RRTMGP()
    {
        return slm_to_rad_vars;
    }

    std::unordered_map<std::string, std::string>& get_rad_output_map()
    {
        return rad_output_map;
    }

    // utility to skip to next line in Header
    //  -- taken from ERF_Checkpoint
    static void GotoNextLine (std::istream& is);

    void WriteCheckpoint(const int &lev, const std::string &checkpointname) const override;

    void ReadCheckpoint(const int &lev, const std::string &checkpointname) override;

    /**
     * Initializes the canopy-adjusted reference height on the current grid.
     */
    void initialize_zrefxy();

    /**
     * Rebuilds SLM-local fields after transferring an evolved state.
     */
    void initialize_processed_state(const bool preserve_checkpoint_fields = false);

private:
    // number of lsm variables (theta)
    int m_lsm_data_size = 44;

    // number of flux variables (qflux, tflux, uflux, vflux)
    int m_lsm_flux_size = 5;

    // LsmVar map (state indices -> LsmVar enum)
    amrex::Vector<int> LsmDataMap;
    amrex::Vector<int> LsmFluxMap;

    // Lsm varnames
    amrex::Vector<std::string> LsmDataName;
    amrex::Vector<std::string> LsmFluxName;

    // geometry for atmosphere
    amrex::Geometry m_geom;

    // geometry for atmospheric level 0
    amrex::Geometry m_geom0;

    // geometry for lsm
    amrex::Geometry m_lsm_geom;

    // geometry for the level-0 soil grid
    amrex::Geometry m_lsm_geom0;

    // domain boundary conditions used by interpolation
    amrex::Vector<amrex::BCRec> m_domain_bcs_type;

    amrex::BoxArray ba_lsm_2d;

    // timestep
    amrex::Real m_dt;

    // domain klo-1 or lsm khi
    int khi_lsm;

    int klo_lsm;

    // independent variables
    amrex::Array<FabPtr, LsmVar_SLM::NumVars> lsm_fab_vars;

    // flux array for conjugate transfer
    amrex::Array<FabPtr, LsmFlux_SLM::NumVars> lsm_fab_flux;

    // Source fields used to transfer SLM state between AMR levels.
    amrex::Vector<amrex::MultiFab*> source_data;
    amrex::Vector<amrex::MultiFab*> source_flux;
    amrex::Vector<amrex::MultiFab*> source_unmapped;

    // Number of grid points in z
    int m_nz_lsm = 7;

    // Size of grid spacing in z
    std::vector<amrex::Real> m_dz_lsm;
    std::string rad_input_file;

    // Input variables
    int landtype0 = 16; // landtype over whole domain

    amrex::Real LAI0 = zero; // default leaf area index
    amrex::Vector<amrex::Real> clay0; // uniform clay content (%)
    amrex::Vector<amrex::Real> sand0; // uniform sand content (%)
    amrex::Vector<amrex::Real> sw0; // uniform soil wetness (fraction)
    amrex::Vector<amrex::Real> st0; // uniform soil temperature (K)
    amrex::Vector<amrex::Real> relax_hgt; // depth dependent soil relaxation function

    amrex::Real tabs_s = zero; // surface temperature (K)
    amrex::Real t00 = 300.0; // constant offset for sstxy
    static constexpr amrex::Real t_canop_max = 343.; // maximum plant temperature

    bool use_terrain = false;

    bool dosoilwnudging = false;
    bool dosoiltnudging = false;
    amrex::Real tausoil = 86400.0;

    amrex::Real z0_soil = 0.01; // baresoil roughness length (m)
    static constexpr amrex::Real tfriz = 273.15;
    amrex::Real mws_mx0 = 50.0; // maximum puddle water storage (mm)

    static constexpr amrex::Real epsv = 0.61;

    static constexpr amrex::Real cp = 1004.0; // specific heat of air [J/kg/K]
    static constexpr amrex::Real cp_water = 4.1796e6; // heat capacity of water [J/m^3K]

    static constexpr amrex::Real sigma = 5.67e-8; // Stefen-Boltzmann constant [W/m^2/K^4]

    // parameters for stomatal resistance
    amrex::Real Rc_max = 5000.0; // maximum stomatal resistance
    amrex::Real T_opt = 298.0; // optimum temperature for transpiration
    amrex::Real zref = myhalf; // height of reference level (m), used in transfer coef
    amrex::MultiFab zrefxy;

    static constexpr amrex::Real IR_emis_urban = 0.9;

    // Model variables
    //   NOTE: these MultiFabs are defined over the grid only in X and Y (no Z component)
    amrex::MultiFab t_canop; // canopy leaf temperature
    amrex::MultiFab mw; // water storage on leaves, kg/m^2 = mm of water
    amrex::MultiFab mws; // water storage in puddle, kg/m^2 = mm of water
    amrex::MultiFab t_skin; // skin temperature, surface level viewed from atmosphere
    amrex::MultiFab t_ground_skin; // massless ground-atmosphere interface temperature
    amrex::MultiFab t_sfc; // 'surface' level temperature. for vegetated, t_sfc = t_cas, for non-vegetated = t_sfc=top soil t
    amrex::MultiFab q_sfc; // 'surface' level specific humidity, for vegetated, q_sfc = q_cas, for non-vegetated = top soil specific humidity
    amrex::MultiFab t_cas; // canopy air space temperature
    amrex::MultiFab q_cas; // canopy air space specific humidity

    amrex::MultiFab q_gr; // soil pore space specific humidity
    amrex::MultiFab sdew; // soil dew condition

    // TODO: placeholder landmask array - should probably use the one defined by ERF
    amrex::iMultiFab landmask;
    amrex::iMultiFab landtype; // IGBP landtype (1 to 16)
    amrex::iMultiFab nroot; // number of rooted soil layers used by the Noah BTR=1 stress factor
    amrex::MultiFab LAI; // leaf area index
    amrex::MultiFab SAI; // stem area index
    amrex::MultiFab sstxy;

    // Vegetation parameter lists
    amrex::iMultiFab vegetype; // vegetation type (0 - baresoil, non-zero for vegetated)
    amrex::MultiFab vege_YES; // 0-baresoil, 1-vegetated
    amrex::MultiFab cp_vege; // dry heat capacity of vegetation [J/m2/K]
    amrex::MultiFab cbiom; // NoahMP canopy biomass heat capacity parameter [m]
    amrex::MultiFab dleaf; // NoahMP characteristic leaf dimension [m]
    amrex::MultiFab z0_sfc; // surface level roughness: for vegetated z0_sfc=z0_vege, for non-vegetated, z0_sfc=z0_soil
    amrex::MultiFab Khai_L; // leaf distribution factor (called XL in NOAHMP radiation)
    amrex::MultiFab phi_1;
    amrex::MultiFab phi_2;
    amrex::MultiFab IR_emis_vege; // vegetation emissivity for infrared radiation
    amrex::MultiFab IR_emis_soil; // soil emissivity for infrared radiation
    amrex::MultiFab IR_emis_grnd; // ground emissivity for infrared radiation
    amrex::MultiFab ztop; // height of vegetation (called HVT in NOAHMP radiation)
    amrex::MultiFab disp_hgt; // displacement height
    amrex::MultiFab Rgl; // parameter for stomatal resistance
    amrex::MultiFab Rc_min; // minimum stomatal resistance
    amrex::MultiFab hs_rc; // parameter for stomatal resistance
    amrex::MultiFab precip_extinc; // precipitation extinction coef
    amrex::MultiFab precip_sfc; // precipitation reaching to soil top through canopy
    amrex::MultiFab mw_mx; // max of mw
    amrex::MultiFab mws_mx; // max of mws
    amrex::MultiFab BAI; // Basal area index, sq.feet/acre
    amrex::MultiFab IMPERV;

    // Packed MF containing 1D arrays of variables used in soil temperature and soil water
    //  soilt: st_cond, st_capa, st_eff_cond, alpha, beta
    //  soilw: st_wgt, sh_eff_cond, sh_eff_vel, alpha, beta, sdepth_mm
    amrex::MultiFab soilt_vars;
    amrex::MultiFab soilw_vars;

    amrex::MultiFab mw_inc;

    // Soil hydrology variables
    amrex::MultiFab evapo_dry;
    amrex::MultiFab evapo_wet;

    // Sensible heat fluxes
    amrex::MultiFab shf_canop; // flux transfer between canopy and air space [W/m^2]
    amrex::MultiFab shf_soil; // flux transfer between soil surface and air space [W/m^2]
    amrex::MultiFab shf_air; // flux transfer between air space and the reference level [W/m^2]

    // Latent heat fluxes
    amrex::MultiFab lhf_canop;
    amrex::MultiFab lhf_soil;
    amrex::MultiFab lhf_air;
    amrex::MultiFab evp_canop;
    amrex::MultiFab evp_soil;
    amrex::MultiFab evp_air;

    // Albedos
    amrex::MultiFab albedovis_v, albedonir_v; // visible and NIR albedos for vegetation
    amrex::MultiFab albedovis_s, albedonir_s; // visible and NIR albedos for soil


    // Resistances
    amrex::MultiFab r_a; // aerodynamic resistance between surface and ref level
    amrex::MultiFab r_b; // leaf boundary layer resistance
    amrex::MultiFab r_c; // leaf-level stomatal resistance
    amrex::MultiFab r_d; // undercanopy aerodynamic resistance
    amrex::MultiFab r_soil;
    static constexpr amrex::Real r_litter = zero; // litter resistance [s/m]


    amrex::MultiFab net_rad;
    amrex::MultiFab wet_canop; // wetted portion of canopy

    amrex::Real time = zero;

    amrex::Real start_time = -1.0;
    bool first_step = true;
    bool m_updated = false;

    std::string parameter_file = "NoahmpTable.TBL";
    std::string veg_dataset = "modis"; // modis or usgs
    std::string soil_dataset = "stas"; // stas or stas_ruc
    int num_landtypes = -1; // Number of land types in sai/lai, parameter tables
    int num_veg_params = -1;
    int num_soil_params = -1;

    SLMParameterTable full_params;
    amrex::Vector<std::string> param_veg_categories;
    amrex::Vector<std::string> param_soil_categories;

    std::unordered_map<std::string, std::unique_ptr<amrex::Gpu::DeviceVector<amrex::Real>>> d_soil_params;
    std::unordered_map<std::string, std::unique_ptr<amrex::Gpu::DeviceVector<amrex::Real>>> d_veg_params;
    std::unordered_map<std::string, std::unique_ptr<amrex::Gpu::DeviceVector<amrex::Real>>> d_rad_params;

    std::vector<std::vector<amrex::Real>> sai_table; // table containing monthly SAI values for each land type
    std::vector<std::vector<amrex::Real>> lai_table; // table containing monthly LAI values for each land type

    amrex::Gpu::DeviceVector<amrex::Real> d_lai_curr; // GPU LAI values for current month only
    amrex::Gpu::DeviceVector<amrex::Real> d_lai_next; // GPU LAI values for next month
    amrex::Gpu::DeviceVector<amrex::Real> d_sai_curr; // GPU SAI values for current month only
    amrex::Gpu::DeviceVector<amrex::Real> d_sai_next; // GPU SAI values for next month

    bool use_param_file = false;  // whether to initialize using parameters from file
    bool interpolate_lai = false; // whether to interpolate LAI+SAI between monthly values for each landtype (Noah-MP style, dynamic veg=4)
    bool use_wrf_lai = false;     // whether to use WRF-Noah style LAI (interpolate between min and max lai)
    bool params_updated = false;  // whether soil/vegetation parameters have been updated from file to recompute

    // NOAHMP radiation parameters
    int opt_rad = 1;  // Radiation option: 1=modified two-stream (default), 2=two-stream, 3=gap=f(veg)
    int opt_alb = 1;  // Snow albedo option: 1=BATS, 2=CLASS

    // NOAHMP radiation parameters from NoahmpTable.TBL
    // ALBSAT_VIS: saturated soil albedos VIS (20 soil color classes)
    // ALBSAT_NIR: saturated soil albedos NIR (20 soil color classes)
    // ALBDRY_VIS: dry soil albedos VIS (20 soil color classes)
    // ALBDRY_NIR: dry soil albedos NIR (20 soil color classes)
    // ALBLAK: albedo frozen lakes: 1=vis, 2=nir
    // OMEGAS: two-stream parameter omega for snow: 1=vis, 2=nir
    amrex::Real betads_rad = myhalf;  // two-stream parameter betad for snow
    amrex::Real betais_rad = myhalf;  // two-stream parameter betai for snow

    // Snow albedo parameters (from &noahmp_global_parameters)
    amrex::Real tau0_rad = 1.0e6;         // tau0 from yang97 eqn. 10a
    amrex::Real grain_growth_rad = 5000.; // growth from vapor diffusion yang97 eqn. 10b
    amrex::Real extra_growth_rad = 10.;   // extra growth near freezing yang97 eqn. 10c
    amrex::Real dirt_soot_rad = 0.3;      // dirt and soot term yang97 eqn. 10d
    amrex::Real rsurf_exp = 5.0;          // NoahMP option-1 soil resistance exponent
    amrex::Real bats_cosz_rad = two;      // zenith angle snow albedo adjustment; b in yang97 eqn. 15
    amrex::Real bats_vis_new_rad = 0.95;  // new snow visible albedo
    amrex::Real bats_nir_new_rad = 0.65;  // new snow nir albedo
    amrex::Real bats_vis_age_rad = 0.2;   // age factor for diffuse visible snow albedo yang97 eqn. 17
    amrex::Real bats_nir_age_rad = myhalf;   // age factor for diffuse nir snow albedo yang97 eqn. 18
    amrex::Real bats_vis_dir_rad = 0.4;   // cosz factor for direct visible snow albedo yang97 eqn. 15
    amrex::Real bats_nir_dir_rad = 0.4;   // cosz factor for direct nir snow albedo yang97 eqn. 16
    amrex::Real swemx_rad = one;          // new snow mass to fully cover old snow (mm)

    // NOAHMP emissivity parameters for longwave radiation from NoahmpTable.TBL
    amrex::Real snow_emis_rad = 0.95;     // snow emissivity
    amrex::Real eg_soil_rad = 0.97;       // ground emissivity for soil (IST=1)
    amrex::Real eg_lake_rad = 0.98;       // ground emissivity for lake (IST=2)

    // NOAHMP vegetation radiation parameters are stored in d_veg_params:
    // RHOL_VIS: leaf reflectance: visible
    // RHOL_NIR: leaf reflectance: NIR
    // RHOS_VIS: stem reflectance: visible
    // RHOS_NIR: stem reflectance: NIR
    // TAUL_VIS: leaf transmittance: visible
    // TAUL_NIR: leaf transmittance: NIR
    // TAUS_VIS: stem transmittance: visible
    // TAUS_NIR: stem transmittance: NIR
    // RC: crown radius (m)
    // HVB: bottom of canopy (m) - NOAHMP name: HVB
    // DEN: tree density (no. of trunks per m2) - NOAHMP name: DEN
    // CWPVT: canopy wind extinction parameter - NOAHMP name: CWPVT
    // NOTE: NOAHMP parameter "xl" (leaf/stem orientation index) uses existing Khai_L from init_landtype()
    // NOTE: NOAHMP parameter "hvt" (top of canopy height) uses existing ztop from init_landtype()

    // MultiFabs to store NOAHMP radiation state variables
    amrex::MultiFab albold_noahmp;  // snow albedo at last time step (for CLASS option)
    amrex::MultiFab tauss_noahmp;   // non-dimensional snow age

    std::vector<std::vector<amrex::Real>> param_table;
    amrex::TableData<amrex::Real, 2> d_param_table;

    amrex::Real m_calday = -1.0;
    int m_orbital_year = -1;
    int m_orbital_mon = -1;
    int m_orbital_day = -1;
    int m_orbital_sec = -1;

    amrex::MultiFab *z_phys_nd; // input z_phys_nd from ERF

    // Extra MFs for output:
    amrex::MultiFab slm_diag;

    amrex::Vector<const amrex::MultiFab*> slm_to_rad_vars;

    // Defines the WRFInput names and corresponding SLM name used to
    // fill variables at initialization
    inline static std::unordered_map<std::string, std::string> wrfinput_map = {
        {"DZS", "soil_thickness"},
        {"ZS", "node_z"},
        {"TSLB", "tsoil"},
        {"SMOIS", "wsoil"},
        {"LAI", "lai"},
        {"IVGTYP", "vegtype"},
        {"ISLTYP", "soiltype"},
        {"TSK", "tsurf"},
        {"VEGFRA", "veg_frac"},
        {"SHDMIN", "veg_frac_min"},
        {"SHDMAX", "veg_frac_max"},
        //{"TMN", "soilt_obs"}
    };

    // map RRTMGP names to SLM names for radiation data
    inline static std::unordered_map<std::string, std::string> rad_output_map = {
        {"cos_zenith_angle",   "cos_zenith"},
        {"sw_flux_dn", ""},
        {"sw_flux_dn_dir_vis", "SW_dw_dir_vis"},
        {"sw_flux_dn_dir_nir", "SW_dw_dir_nir"},
        {"sw_flux_dn_dif_vis", "SW_dw_dif_vis"},
        {"sw_flux_dn_dif_nir", "SW_dw_dif_nir"},
        {"lw_flux_dn", "LW_dw"}
    };

    // whether to initialize SLM from inputs file or from wrfinput
    bool use_wrfinput = false;
    bool wrfinput_initialized = false;

    // list of extra MFs in lsm_fab_vars that are not checkpointed by ERF
    amrex::Vector<int> unmapped_fields;

    // List of 2D output variables that are declared as 3D
    // These will only be saved at k = 0 (surface level) to save space
    // TODO: these variables should be redefined to be 2D
    inline static const std::unordered_set<int> surf_vars = {
      LsmVar_SLM::tsurf,
      LsmVar_SLM::ustar,
      LsmVar_SLM::tstar,
      LsmVar_SLM::qstar,
      LsmVar_SLM::tref,
      LsmVar_SLM::qref,
      LsmVar_SLM::pref,
      LsmVar_SLM::dref,
      LsmVar_SLM::uref,
      LsmVar_SLM::vref,
      LsmVar_SLM::precipref,
      LsmVar_SLM::swdsvisxyref,
      LsmVar_SLM::swdsnirxyref,
      LsmVar_SLM::swdsvisdxyref,
      LsmVar_SLM::swdsnirdxyref,
      LsmVar_SLM::lwref,
      LsmVar_SLM::coszrsxy,
      LsmVar_SLM::tv,
      LsmVar_SLM::mv,
      LsmVar_SLM::flbu,
      LsmVar_SLM::flbv,
      LsmVar_SLM::flbq,
      LsmVar_SLM::flbt,
      LsmVar_SLM::prsfc,
      LsmVar_SLM::lai,
      LsmVar_SLM::vegtype,
      LsmVar_SLM::soiltype,
      LsmVar_SLM::veg_frac,
      LsmVar_SLM::veg_frac_min,
      LsmVar_SLM::veg_frac_max,
      LsmVar_SLM::emis_sfc,
      LsmVar_SLM::alb_nir_sfc,
      LsmVar_SLM::alb_vis_sfc,
      LsmVar_SLM::alb_nir_sfc_diff,
      LsmVar_SLM::alb_vis_sfc_diff,
    };

    static constexpr const char* LsmVarName_Full[] = {
      "tsurf", "ustar", "tstar", "qstar", "ref_t", "ref_q", "ref_p", "ref_d", "ref_u", "ref_v", "ref_precip",
      "SW_dw_dir_vis", "SW_dw_dir_nir", "SW_dw_dif_vis", "SW_dw_dif_nir", "LW_dw", "zenith",
      "tveg", "mv", "tsoil", "wsoil", "sand", "clay", "soil_thickness", "node_z", "interface_z",
      "soilw_nudge", "soilt_nudge", "sw_over_sat", "soil_relax_hgt", "soilt_obs", "soilw_obs", "soil_transp_frac", "w_s_FC", "w_s_WP",
      "sst_capa", "sst_cond", "poro_soil", "theta_FC", "theta_WP", "m_pot_sat", "Bconst", "ks","lai","vegtype","soiltype",
      "veg_frac", "veg_frac_min", "veg_frac_max", "surface_u", "surface_v", "surface_vapor", "surface_heat", "precip_soil", "emis_sfc", "alb_nir_sfc", "alb_vis_sfc", "alb_nir_sfc_diff", "alb_vis_sfc_diff"
    };

    static_assert(sizeof(LsmVarName_Full)/sizeof(*LsmVarName_Full) == LsmVar_SLM::NumVars);

    inline static const std::vector<std::string> diag_names = {
        "precip",
        "precip_in",
        "precip_sfc",
        "runoff_sfc",
        "drain",
        "drain_flux",
        "evapo_wet",
        "grflux"
    };

#ifdef ERF_USE_NETCDF
    // List of output variables that are constant in time
    inline static const std::unordered_set<int> const_vars = {
      LsmVar_SLM::sand,
      LsmVar_SLM::clay,
      LsmVar_SLM::s_depth,
      LsmVar_SLM::node_z,
      LsmVar_SLM::interface_z,
      LsmVar_SLM::sw_over_sat,
      LsmVar_SLM::soil_relax_hgt,
      LsmVar_SLM::soilt_obs,
      LsmVar_SLM::soilw_obs,
      LsmVar_SLM::w_s_FC,
      LsmVar_SLM::w_s_WP,
      LsmVar_SLM::sst_capa,
      LsmVar_SLM::sst_cond,
      LsmVar_SLM::poro_soil,
      LsmVar_SLM::theta_FC,
      LsmVar_SLM::theta_WP,
      LsmVar_SLM::m_pot_sat,
      LsmVar_SLM::Bconst,
      LsmVar_SLM::ks
    };

    inline static const std::vector<std::string> soilt_var_names = {
        "st_cond",
        "st_capa",
        "st_eff_cond",
        "st_alpha",
        "st_beta"
    };

    inline static const std::vector<std::string> soilw_var_names = {
        "sw_wgt",
        "sw_sh_eff_cond",
        "sw_sh_eff_vel",
        "sw_alpha",
        "sw_beta",
        "sw_sdepth_mm"
    };

    amrex::MultiFab rad_input_data;
    int num_rad_times = 0;
    amrex::Vector<amrex::Real> rad_times;
    inline static const std::vector<std::string> rad_names = {
      "S_SWVIS", "S_SWNIR", "S_SWVISD", "S_SWNIRD", "S_COSZRS", "S_LWDS"};

#endif
};
#endif
