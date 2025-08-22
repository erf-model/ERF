#ifndef ERF_SLM_H
#define ERF_SLM_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_TableData.H>
#include <AMReX_ParmParse.H>

#include <ERF_NullSurf.H>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>

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
      //tsurf,           // surface temperature
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
      rootF,           // fraction of total root length of a vegetation in each soil layer
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
      // model outputs
      flbu,            // surface u-momentum
      flbv,            // surface v-momentum
      flbq,            // surface vapor flux (latent heat flux)
      flbt,            // surface heat flux
      prsfc,           // precipitation rate reaching soil
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

/* See Chen & Dudhia (2001) https://doi.org/10.1175/1520-0493(2001)129<0569:CAALSH>2.0.CO;2 */
class SLM : public NullSurf {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // Constructor
    SLM () {}

    // Destructor
    virtual ~SLM () = default;

    // Set thermo and grid properties
    void
    Define (SolverChoice& sc) override
    {
        // NOTE: We should parse things from sc here,
        //       but they are hard coded because this
        //       is a demonstration for now.

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
          const amrex::Real& dt,
          std::unique_ptr<amrex::MultiFab>& z_phys_cc_in) override;

    // Wrapper to do all the updating
    void
    Advance (const amrex::Real& dt) override
    {
        m_dt = dt;
        time += m_dt;

        // override ERF inputs from reference data if testing flag is set
        this->time_interp_from_ref();

        if (first_step)
        {
            // initialize values using first reference data
            //this->init_slm_vars();
            first_step = false;

            if (use_wrfinput) {
                // TODO: fix - WRFinput sets 3D MFs in lsm_data - these need to be copied into SLM local 2D MFs
                int d_khi_lsm = khi_lsm;

                for ( amrex::MFIter mfi(landtype,TileNoZ()); mfi.isValid(); ++mfi) {
                    amrex::Box bx2d = mfi.tilebox();

                    amrex::Box bx3d_mf = bx2d;
                    bx3d_mf.setBig(2, khi_lsm);
                    bx3d_mf.setSmall(2, khi_lsm);

                    // Copies from 3D MF at k=-1 filled from WRFInput to 2D MF defined at k=0
                    //landtype[mfi].template copy<amrex::RunOn::Device>((*lsm_fab_vars[LsmVar_SLM::vegtype])[mfi], bx2d);
                    //LAI[mfi].template copy<amrex::RunOn::Device>((*lsm_fab_vars[LsmVar_SLM::lai])[mfi], bx2d);

                    auto vegtype_arr = lsm_fab_vars[LsmVar_SLM::vegtype]->const_array(mfi);
                    auto in_soiltype_arr = lsm_fab_vars[LsmVar_SLM::soiltype]->const_array(mfi);
                    auto in_lai_arr = lsm_fab_vars[LsmVar_SLM::lai]->const_array(mfi);
                    auto in_tsk_arr = lsm_fab_vars[LsmVar_SLM::tsurf]->array(mfi);

                    auto landmask_arr = landmask.array(mfi);
                    auto landtype_arr = landtype.array(mfi);
                    auto lai_arr = LAI.array(mfi);

                    auto sstxy_arr = sstxy.array(mfi);

                    amrex::ParallelFor(bx2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                        landtype_arr(i,j,0) = vegtype_arr(i,j,d_khi_lsm);
                        //soiltype_arr(i,j,0) = in_soiltype_arr(i,j,d_khi_lsm);
                        lai_arr(i,j,0) = in_lai_arr(i,j,d_khi_lsm);

                        sstxy_arr(i,j,0) = in_tsk_arr(i,j,d_khi_lsm);

                        in_tsk_arr(i,j,0) = in_tsk_arr(i,j,d_khi_lsm);

                        //amrex::Print() << " i = " << i << " j = " << j << ": landtype = " << landtype_arr(i,j,0) << " soiltype = " << in_soiltype_arr(i,j,d_khi_lsm) << " LAI = " << lai_arr(i, j, 0) << std::endl;

                        if (landtype_arr(i,j,0) == 17) {
                            amrex::Print() << "     -- water cell - setting landtype=0 landmask=0" << std::endl;
                            landmask_arr(i,j,0) = 0;
                            landtype_arr(i,j,0) = 0;
                        } else if (landtype_arr(i,j,0) > 16) {
                            amrex::Print() << "     -- urban cell - setting landtype=13" << std::endl;
                            landtype_arr(i,j,0) = 13;
                        }
                    });
                }

                // If using WRFInput, SLM init needs to be delayed here to use correct values set from WRFInput
                slm_init();

                init_slm_vars();
            }

            // update zref using updated z_phys_cc
            // TODO: fix this - use_terrain flag is deprecated
            if (!set_from_file) {
                amrex::ParmParse pp_erf("erf");
                pp_erf.query("use_terrain", use_terrain);

                amrex::Real zlo      = m_geom.ProbLo(2);
                amrex::Real dz       = m_geom.CellSize(2);
                for ( amrex::MFIter mfi(landtype,TileNoZ()); mfi.isValid(); ++mfi) {
                    const amrex::Box& xybx      = mfi.growntilebox(0);
                    const amrex::Array4<const amrex::Real>& z_cc_arr = (use_terrain) ? z_phys_cc->const_array(mfi) : amrex::Array4<amrex::Real>{};
                    auto ztop_arr = ztop.array(mfi);
                    auto zrefxy_arr = zrefxy.array(mfi);

                    ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                        amrex::Real zcc = (z_cc_arr) ? z_cc_arr(i,j,0) : zlo + 0.5*dz;
                        zrefxy_arr(i,j,0) = ztop_arr(i,j,0) + zcc;
                        //amrex::Print() << " i = " << i << " j = " << j << ": zrefxy = " << zrefxy_arr(i,j,0) << std::endl;
                    });
                }
            }
        }

        this->AdvanceSLM();

    }

    // Advance the lsm state vars
    void
    AdvanceSLM ();

    // Get state vars from lsm class
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_size && varIdx >= 0);
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < LsmVar_SLM::NumVars && lsmIdx>=0);
        return lsm_fab_vars[lsmIdx].get();
    }

    // Get flux vars from lsm class
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_size && varIdx >= 0);
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < LsmVar_SLM::NumVars && lsmIdx>=0);
        return lsm_fab_flux[lsmIdx].get();
    }

    // Get lsm geometry
    amrex::Geometry
    Lsm_Geom ( ) override { return m_lsm_geom; }

    // Get number of vars lsm class contains
    int
    Lsm_Data_Size () override { return SLM::m_lsm_size; }

    // Get variable names
    std::string
    Lsm_VarName (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < SLM::m_lsm_size && varIdx >= 0);
        return LsmVarName[varIdx];
    }

    int
    Lsm_VarIndex (std::string varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<SLM::m_lsm_size; ++idx) {
            if (lc_varname == amrex::toLower(Lsm_VarName(idx))) {
                varIdx = idx;
            }
        }
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
        if (first_step && time == 0.0) {
            this->init_slm_vars();
            time += m_dt; // update time here so that first timestep time is correct
        }
        //this->Diagnose();
        //this->Compute_Coefficients();
    }

    void
    Update_State_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_Lsm_to_State(cons_in);
    }

    void
    set_flux_inputs(const amrex::MultiFab* sw_lw_fluxes_in,
                    const amrex::MultiFab* zenith_in) override;

    void
    set_precip_input(const amrex::MultiFab* precip_in) override;

    void
    set_terrain_inputs(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& sst_in,
                       const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& lmask_in) override;

    static std::vector<std::vector<amrex::Real>>
    read_cols(const std::string& fname, const int skip_nlines = 1);

    int Lsm_Type() override
    {
        return LSM_Type::SLM;
    }

    void writeSLM_Data(const PlotFileType plotfile_type, const amrex::Real time, const std::string plot_prefix, const int level_step);

    void init_from_file();

    void time_interp_from_ref();

    void slm_init();

    void init_landtype();

    void init_soil_tw();

    void vege_root_init();

    void init_slm_vars();

    void radiative_fluxes(const amrex::MFIter &mfi);

    void transfer_coeff(const amrex::MFIter &mfi);

    void resistances(const amrex::MFIter &mfi);

    void vapor_fluxes(const amrex::MFIter &mfi);

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

private:
    // number of lsm variables (theta)
    int m_lsm_size = 34;

    // LsmVar map (state indices -> LsmVar enum)
    amrex::Vector<int> LsmVarMap;

    // Lsm varnames
    amrex::Vector<std::string> LsmVarName;

    // geometry for atmosphere
    amrex::Geometry m_geom;

    // geometry for lsm
    amrex::Geometry m_lsm_geom;

    amrex::BoxArray ba_lsm_2d;

    // timestep
    amrex::Real m_dt;

    // domain klo-1 or lsm khi
    int khi_lsm;

    int klo_lsm;

    // constants
    // amrex::Real m_fac_cond;
    // amrex::Real m_fac_fus;
    // amrex::Real m_fac_sub;
    // amrex::Real m_gOcp;

    // independent variables
    amrex::Array<FabPtr, LsmVar_SLM::NumVars> lsm_fab_vars;

    // flux array for conjugate transfer
    amrex::Array<FabPtr, LsmVar_SLM::NumVars> lsm_fab_flux;

    // Vars that should be parsed
    // ==========================
    // Number of grid points in z
    int m_nz_lsm = 7;

    // Size of grid spacing in z
    //amrex::Real m_dz_lsm = 0.1; // 3 [m] below surface
    std::vector<amrex::Real> m_dz_lsm;

    // Specific heat
    amrex::Real m_cp_soil =  1.26e6; // [J/m^3 K]

    // Conductivity
    amrex::Real m_k_soil = 0.2; // dry soil [W/m K]

    // Theta Dirichlet value at lowest point below surface
    amrex::Real m_theta_dir = 400.0; // [K]

    // Thermal diffusivity
    amrex::Real m_d_soil = m_k_soil / m_cp_soil; // [m^2/s]

    std::string soil_file; // file for reading in initial soil values at each later
    std::string rad_input_file;

    // Input variables
    // TODO: move these to inputs parsing
    int landtype0 = 16; // landtype over whole domain

    amrex::Real LAI0 = 0.0; // default leaf area index
    amrex::Vector<amrex::Real> clay0; // uniform clay content (%)
    amrex::Vector<amrex::Real> sand0; // uniform sand content (%)
    amrex::Vector<amrex::Real> sw0; // uniform soil wetness (fraction)
    amrex::Vector<amrex::Real> st0; // uniform soil temperature (K)
    amrex::Vector<amrex::Real> relax_hgt; // depth dependent soil relaxation function

    amrex::Real tabs_s = 0.0; // surface temperature (K)
    amrex::Real t00 = 300.0; // constant offset for sstxy
    
    bool use_terrain = false;

    //   TODO: implement reading these input files
    bool readlandtype = false; // read landtype map from file
    bool readLAI = false; // read LAI file
    std::string landtypefile;
    std::string LAIfile;

    bool dosoilwnudging = false;
    bool dosoiltnudging = false;
    amrex::Real tausoil = 86400.0;

    amrex::Real z0_soil = 0.0387; // baresoil roughness length (m)
    static constexpr amrex::Real tfriz = 273.15;
    amrex::Real mws_mx0 = 50.0; // maximum puddle water storage (mm)

    static constexpr amrex::Real epsv = 0.61;

    amrex::Real pres0 = p_0; // reference surface pressure [Pa]

    static constexpr amrex::Real cp = 1004.0; // specific heat of air [J/kg/K]
    static constexpr amrex::Real cp_water = 4.1796e6; // heat capacity of water [J/m^3K]

    static constexpr amrex::Real sigma = 5.67e-8; // Stefen-Boltzmann constant [W/m^2/K^4]

    // parameters for stomatal resistance
    amrex::Real Rc_max = 5000.0; // maximum stomatal resistance
    amrex::Real T_opt = 298.0; // optimum temperature for transpiration

    amrex::Real zref = 0.5; // height of reference level (m), used in transer coef
    amrex::MultiFab zrefxy;

    // Model variables
    //   NOTE: these MultiFabs are defined over the grid only in X and Y (no Z component)
    amrex::MultiFab t_canop; // canopy leaf temperature
    amrex::MultiFab mw; // water storage on leaves, kg/m^2 = mm of water
    amrex::MultiFab mws; // water storage in puddle, kg/m^2 = mm of water
    amrex::MultiFab t_skin; // skin temperature, surface level viewed from atmosphere

    amrex::MultiFab t_cas; // canopy air space temperature
    amrex::MultiFab q_cas; // canopy air space specific humidity

    amrex::TableData<amrex::Real, 1> soilw_inc;

    // TODO: placeholder landmask array - should probably use the one defined by ERF
    amrex::iMultiFab landmask;
    amrex::iMultiFab landtype; // IGBP landtype (1 to 16)
    amrex::MultiFab LAI; // leaf area index

    amrex::MultiFab sstxy;

    // Vegetation parameter lists
    amrex::iMultiFab vegetype; // vegetation type (0 - baresoil, non-zero for vegetated)
    amrex::MultiFab vege_YES; // 0-baresoil, 1-vegetated
    amrex::MultiFab cp_vege; // heat capacity of vegetation
    amrex::MultiFab z0_sfc; // surface level roughness: for vegetated z0_sfc=z0_vege, for non-vegetated, z0_sfc=z0_soil
    amrex::MultiFab Khai_L; // leaf distribution factor
    amrex::MultiFab phi_1;
    amrex::MultiFab phi_2;
    amrex::MultiFab IR_emis_vege; // vegetation emissivity for infrared radiation
    amrex::MultiFab IR_emis_soil; // soil emissivity for infrared radiation
    amrex::MultiFab ztop; // height of vegetation
    amrex::MultiFab disp_hgt; // displacement height
    amrex::MultiFab Rgl; // parameter for stomatal resistance
    amrex::MultiFab Rc_min; // minimum stomatal resistance
    amrex::MultiFab hs_rc; // parameter for stomatal resistance
    amrex::MultiFab rootL; // root length
    amrex::MultiFab root_a; // coef a for root density
    amrex::MultiFab root_b; // coef b for root density
    amrex::MultiFab precip_extinc; // precipitation extinction coef
    amrex::MultiFab mw_mx; // max of mw
    amrex::MultiFab mws_mx; // max of mws
    amrex::MultiFab BAI; // Basal area index, sq.feet/acre

    amrex::TableData<amrex::Real, 1> alpha;
    amrex::TableData<amrex::Real, 1> beta;

    // Packed MF containing 1D arrays of variables used in soil temperature and soil water
    //  soilt: st_cond, st_capa, st_eff_cond, alpha, beta
    //  soilw: st_wgt, sh_eff_cond, sh_eff_vel, alpha, beta, sdepth_mm
    amrex::MultiFab soilt_vars;
    amrex::MultiFab soilw_vars;

    amrex::MultiFab mw_inc;

    //amrex::MultiFab ustar, tstar;

    //amrex::Real mom_trans_coef; // surface momentum transfer coef.
    //amrex::Real heat_trans_coef; // surface heat transfer coef.
    amrex::Real xsi; // M-O stability parameter
    amrex::Real RiB; // bulk Richardson number

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

    // Albedos
    amrex::MultiFab albedovis_v, albedonir_v; // visible and NIR albedos for vegetation
    amrex::MultiFab albedovis_s, albedonir_s; // visible and NIR albedos for soil


    // Resistances
    //   TODO: in the fortran version, these are scalar quantities but are updated over i,j?
    amrex::MultiFab r_a; // aerodynamic resistance between surface and ref level
    amrex::MultiFab r_b; // leaf boundary layer resistance
    amrex::MultiFab r_c; // stomatal resistance
    amrex::MultiFab r_d; // undercanopy aerodynamic resistance
    amrex::MultiFab r_soil;
    static constexpr amrex::Real r_litter = 0.0; // not used currently


    amrex::MultiFab net_rad;
    amrex::MultiFab wet_canop; // wetted portion of canopy

    // Variables related to setting initial conditions (for testing purposes only)
    bool set_from_file = false;
    int num_ref_inputs = 1;
    std::string ref_sounding_file, ref_flux_file, ref_sst_file;
    std::vector<std::vector<amrex::Real>> sounding;
    std::vector<std::vector<amrex::Real>> fluxes;
    std::vector<std::vector<amrex::Real>> sst;
    int time_index = 0;
    int start_time_index = 0;
    amrex::Real time = 0.0;

    amrex::Real start_time = -1.0;
    amrex::Real time_unit = 1.0;
    bool first_step = true;

    amrex::MultiFab *z_phys_cc; // input z_phys_cc from ERF

    // Extra MFs for output:
    amrex::MultiFab slm_diag;

    amrex::Vector<const amrex::MultiFab*> slm_to_rad_vars;

    // Defines the WRFInput names and corresponding SLM name used to
    // fill variables at initialization
    inline static std::unordered_map<std::string, std::string> wrfinput_map = {
        {"ZS", "soil_thickness"},
        {"DZS", "node_z"},
        {"TSLB", "tsoil"},
        {"SMOIS", "wsoil"},
        {"LAI", "lai"},
        {"IVGTYP", "vegtype"}, // TODO: these two should be replaced by ERF input directly
        {"ISLTYP", "soiltype"},
        {"TSK", "tsurf"},
    };

    // whether to initialize SLM from inputs file or from wrfinput
    bool use_wrfinput = false;

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
    };

    static constexpr const char* LsmVarName_Full[] = {
      "tsurf", "ustar", "tstar", "qstar", "ref_t", "ref_q", "ref_p", "ref_d", "ref_u", "ref_v", "ref_precip",
      "SW_dw_dir_vis", "SW_dw_dir_nir", "SW_dw_dif_vis", "SW_dw_dif_nir", "LW_dw", "zenith",
      "tveg", "mv", "tsoil", "wsoil", "sand", "clay", "soil_thickness", "node_z", "interface_z",
      "soilw_nudge", "soilt_nudge", "sw_over_sat", "soil_relax_hgt", "soilt_obs", "soilw_obs", "rootF", "w_s_FC", "w_s_WP",
      "sst_capa", "sst_cond", "poro_soil", "theta_FC", "theta_WP", "m_pot_sat", "Bconst", "ks","lai","vegtype","soiltype",
      "surface_u", "surface_v", "surface_vapor", "surface_heat", "precip_soil"
    };

    static_assert(sizeof(LsmVarName_Full)/sizeof(*LsmVarName_Full) == LsmVar_SLM::NumVars);

#ifdef ERF_USE_NETCDF
    int ncfid = -1;

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
      LsmVar_SLM::rootF,
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
