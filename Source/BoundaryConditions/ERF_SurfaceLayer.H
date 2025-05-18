#ifndef ERF_SURFACELAYER_H
#define ERF_SURFACELAYER_H

#include "AMReX_Geometry.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_Interpolater.H"

#include "ERF_IndexDefines.H"
#include "ERF_Constants.H"
#include "ERF_MOSTAverage.H"
#include "ERF_MOSTStress.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_PBLHeight.H"

/** Abstraction layer for different surface layer schemes (e.g. MOST, Cd)
 *
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */
class SurfaceLayer
{

public:

    // Constructor
    explicit SurfaceLayer (const amrex::Vector<amrex::Geometry>& geom,
                           bool& use_rot_surface_flux,
                           std::string a_pp_prefix,
                           amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qv_prim,
                           amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd,
                           const TerrainType& a_terrain_type,
                           amrex::Real start_bdy_time = 0.0,
                           amrex::Real bdy_time_interval = 0.0)
    : m_geom(geom),
      m_rotate(use_rot_surface_flux),
      m_start_bdy_time(start_bdy_time),
      m_bdy_time_interval(bdy_time_interval),
      m_ma(geom,(z_phys_nd[0]!=nullptr),a_pp_prefix,a_terrain_type)
    {
        // We have a moisture model if Qv_prim is a valid pointer
        use_moisture = (Qv_prim[0].get());

        // Get roughness
        amrex::ParmParse pp("erf");
        pp.query("most.z0", z0_const);

        // Specify how to compute the flux
        if (use_rot_surface_flux) {
            flux_type = FluxCalcType::ROTATE;
        } else {
            std::string flux_string_in;
            std::string flux_string{"moeng"};
            auto read_flux = pp.query("surface_layer.flux_type", flux_string_in);
            if (read_flux) { flux_string = amrex::toLower(flux_string_in); }
            if (flux_string == "donelan") {
                flux_type = FluxCalcType::DONELAN;
            } else if (flux_string == "moeng") {
                flux_type = FluxCalcType::MOENG;
            } else if (flux_string == "custom") {
                flux_type = FluxCalcType::CUSTOM;
            } else {
                amrex::Abort("Undefined MOST flux type!");
            }
        }

        // Include w* to handle free convection (Beljaars 1995, QJRMS)
        pp.query("most.include_wstar", m_include_wstar);

        std::string pblh_string_in;
        std::string pblh_string{"none"};
        auto read_pblh = pp.query("most.pblh_calc", pblh_string_in);
        if (read_pblh) { pblh_string = amrex::toLower(pblh_string_in); }
        if (pblh_string == "none") {
            pblh_type = PBLHeightCalcType::None;
        } else if (pblh_string == "mynn25") {
            pblh_type = PBLHeightCalcType::MYNN25;
        } else if (pblh_string == "mynnedmf") {
            pblh_type = PBLHeightCalcType::MYNN25;
        } else if (pblh_string == "ysu") {
            pblh_type = PBLHeightCalcType::YSU;
        } else {
            amrex::Abort("Undefined PBLH calc type!");
        }

        // Get surface temperature
        auto erf_st = pp.query("most.surf_temp", surf_temp);
        if (erf_st) { default_land_surf_temp = surf_temp; }

        // Custom type user must specify the fluxes
        if (flux_type == FluxCalcType::CUSTOM) {
            theta_type = ThetaCalcType::HEAT_FLUX;
            pp.query("most.ustar", custom_ustar);
            pp.query("most.tstar", custom_tstar);
            pp.query("most.qstar", custom_qstar);
            if (custom_qstar != 0) {
                AMREX_ASSERT_WITH_MESSAGE(use_moisture, "Specified custom MOST qv flux without moisture model!");
            }
            amrex::Print() << "Using specified ustar, tstar, qstar for MOST = "
                << custom_ustar << " "
                << custom_tstar << " "
                << custom_qstar << std::endl;

        // Specify surface temperature or surface flux
        } else {
            if (erf_st) {
                theta_type = ThetaCalcType::SURFACE_TEMPERATURE;
                pp.query("most.surf_heating_rate", surf_heating_rate); // [K/h]
                surf_heating_rate = surf_heating_rate / 3600.0; // [K/s]
                if (pp.query("most.surf_temp_flux", surf_temp_flux)) {
                    amrex::Abort("Can only specify one of surf_temp_flux or surf_heating_rate");
                }
            } else {
                pp.query("most.surf_temp_flux", surf_temp_flux);
                if (pp.query("most.surf_heating_rate", surf_heating_rate)) {
                    amrex::Abort("Can only specify one of surf_temp_flux or surf_heating_rate");
                }
                if (std::abs(surf_temp_flux) > std::numeric_limits<amrex::Real>::epsilon()) {
                    theta_type = ThetaCalcType::HEAT_FLUX;
                } else {
                    theta_type = ThetaCalcType::ADIABATIC;
                }
            }
        }

        // Make sure the inputs file doesn't try to use most.roughness_type
        std::string bogus_input;
        if (pp.query("most.roughness_type", bogus_input) > 0) {
            amrex::Abort("most.roughness_type is deprecated; use most.roughness_type_land and/or most.roughness_type_sea");
        }

        // Specify how to compute the surface flux over land (if there is any)
        std::string rough_land_string_in;
        std::string rough_land_string{"constant"};
        auto read_rough_land = pp.query("most.roughness_type_land", rough_land_string_in);
        if (read_rough_land) { rough_land_string = amrex::toLower(rough_land_string_in); }
        if (rough_land_string == "constant") {
            rough_type_land = RoughCalcType::CONSTANT;
        } else {
            amrex::Abort("Undefined MOST roughness type for land!");
        }

        // Specify how to compute the surface flux over sea (if there is any)
        std::string rough_sea_string_in;
        std::string rough_sea_string{"charnock"};
        auto read_rough_sea = pp.query("most.roughness_type_sea", rough_sea_string_in);
        if (read_rough_sea) { rough_sea_string = amrex::toLower(rough_sea_string_in); }
        if (rough_sea_string == "charnock") {
            rough_type_sea = RoughCalcType::CHARNOCK;
            pp.query("most.charnock_constant",cnk_a);
            pp.query("most.charnock_viscosity",cnk_visc);
            if (cnk_a > 0) {
                amrex::Print() << "If there is water, Charnock relation with C_a=" << cnk_a
                    << (cnk_visc? " and viscosity" : "")
                    << " will be used"
                    << std::endl;
            } else {
                amrex::Print() << "If there is water, Charnock relation with variable Charnock parameter (COARE3.0)"
                    << (cnk_visc? " and viscosity" : "")
                    << " will be used"
                    << std::endl;
            }
        } else if (rough_sea_string == "coare3.0") {
            rough_type_sea = RoughCalcType::CHARNOCK;
            amrex::Print() << "If there is water, Charnock relation with variable Charnock parameter (COARE3.0)"
                << (cnk_visc? " and viscosity" : "")
                << " will be used"
                << std::endl;
            cnk_a = -1;
        } else if (rough_sea_string == "donelan") {
            rough_type_sea = RoughCalcType::DONELAN;
        } else if (rough_sea_string == "modified_charnock") {
            rough_type_sea = RoughCalcType::MODIFIED_CHARNOCK;
            pp.query("most.modified_charnock_depth",depth);
        } else if (rough_sea_string == "wave_coupled") {
            rough_type_sea = RoughCalcType::WAVE_COUPLED;
        } else if (rough_sea_string == "constant") {
            rough_type_sea = RoughCalcType::CONSTANT;
        } else {
            amrex::Abort("Undefined MOST roughness type for sea!");
        }
    } // constructor

    void
    make_SurfaceLayer_at_level (const int& lev,
                                int nlevs,
                                const amrex::Vector<amrex::MultiFab*>& mfv,
                                std::unique_ptr<amrex::MultiFab>& Theta_prim,
                                std::unique_ptr<amrex::MultiFab>& Qv_prim,
                                std::unique_ptr<amrex::MultiFab>& Qr_prim,
                                std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                                amrex::MultiFab* Hwave,
                                amrex::MultiFab* Lwave,
                                amrex::MultiFab* eddyDiffs,
                                amrex::Vector<amrex::MultiFab*> lsm_data,
                                amrex::Vector<amrex::MultiFab*> lsm_flux,
                                amrex::Vector<std::unique_ptr<amrex::MultiFab>>& sst_lev,
                                amrex::Vector<std::unique_ptr<amrex::MultiFab>>& tsk_lev,
                                amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& lmask_lev)
    {
        // Update MOST Average
        m_ma.make_MOSTAverage_at_level(lev, mfv, Theta_prim, Qv_prim, Qr_prim, z_phys_nd);

        // Get CC vars
        amrex::MultiFab& mf = *(mfv[0]);

        amrex::ParmParse pp("erf");

        // Do we have a time-varying surface roughness that needs to be saved?
        if (lev == 0)
        {
            const int nghost = 0; // ghost cells not included
            int lmask_min = lmask_min_reduce(*lmask_lev[0].get(), nghost);
            amrex::ParallelDescriptor::ReduceIntMin(lmask_min);

            m_var_z0 = (lmask_min < 1) & (rough_type_sea != RoughCalcType::CONSTANT);
            if (m_var_z0) {
                std::string rough_sea_string{"charnock"};
                pp.query("most.roughness_type_sea", rough_sea_string);
                amrex::Print() << "Variable sea roughness (type " << rough_sea_string << ")" << std::endl;
            }
        }

        if ( (lev == 0) || (lev > nlevs-1) ) {
            m_Hwave_lev.resize(nlevs);
            m_Lwave_lev.resize(nlevs);
            m_eddyDiffs_lev.resize(nlevs);

            m_lsm_data_lev.resize(nlevs);
            m_lsm_flux_lev.resize(nlevs);

            m_sst_lev.resize(nlevs);
            m_tsk_lev.resize(nlevs);
            m_lmask_lev.resize(nlevs);

            // Size the MOST params for all levels
            z_0.resize(nlevs);
            u_star.resize(nlevs);
            w_star.resize(nlevs);
            t_star.resize(nlevs);
            q_star.resize(nlevs);
            t_surf.resize(nlevs);
            olen.resize(nlevs);
            pblh.resize(nlevs);
        }


        // Get pointers to SST,TSK and LANDMASK data
        int nt_tot_sst = sst_lev.size();
        m_sst_lev[lev].resize(nt_tot_sst);
        for (int nt(0); nt<nt_tot_sst; ++nt) {
            m_sst_lev[lev][nt] = sst_lev[nt].get();
        }
        int nt_tot_tsk = tsk_lev.size();
        m_tsk_lev[lev].resize(nt_tot_tsk);
        for (int nt(0); nt<nt_tot_tsk; ++nt) {
            m_tsk_lev[lev][nt] = tsk_lev[nt].get();
        }
        int nt_tot_lmask = lmask_lev.size();
        m_lmask_lev[lev].resize(nt_tot_lmask);
        for (int nt(0); nt<nt_tot_lmask; ++nt) {
            m_lmask_lev[lev][nt] = lmask_lev[nt].get();
        }

        // Get pointers to wave data
        m_Hwave_lev[lev] = Hwave;
        m_Lwave_lev[lev] = Lwave;
        m_eddyDiffs_lev[lev] = eddyDiffs;

        // Get pointers to LSM data and Fluxes
        int nvar = lsm_data.size();
        m_lsm_data_lev[lev].resize(nvar);
        m_lsm_flux_lev[lev].resize(nvar);
        for (int n(0); n<nvar; ++n) {
            m_lsm_data_lev[lev][n] = lsm_data[n];
            m_lsm_flux_lev[lev][n] = lsm_flux[n];
        }

        // Check if there is a user-specified roughness file to be read
        std::string fname;
        bool read_z0 = false;
        if ( (flux_type == FluxCalcType::MOENG) ||
             (flux_type == FluxCalcType::ROTATE) ) {
            read_z0 = pp.query("most.roughness_file_name", fname);
        }

        // Attributes for MFs and FABs
        //--------------------------------------------------------
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
            b.setRange(2,0);
        }
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

        // Z0 heights FAB
        //--------------------------------------------------------
        amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
        Arena_Used = amrex::The_Pinned_Arena();
#endif
        amrex::Box bx = amrex::grow(m_geom[lev].Domain(),ng);
        bx.setSmall(2,0);
        bx.setBig(2,0);
        z_0[lev].resize(bx,1,Arena_Used);
        z_0[lev].setVal<amrex::RunOn::Device>(z0_const);
        if (read_z0) read_custom_roughness(lev,fname);

        // 2D MFs for U*, T*, T_surf
        //--------------------------------------------------------
        u_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        u_star[lev]->setVal(1.E34);

        w_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        w_star[lev]->setVal(1.E34);

        t_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        t_star[lev]->setVal(0.0); // default to neutral

        q_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        q_star[lev]->setVal(0.0); // default to dry

        olen[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        olen[lev]->setVal(1.E34);

        pblh[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        pblh[lev]->setVal(1.E34);

        t_surf[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
        t_surf[lev]->setVal(default_land_surf_temp);

        // TODO: Do we want lsm_data to have theta at 0 index always?
        //       Do we want an external enum struct for indexing?
        if (m_sst_lev[lev][0] || m_tsk_lev[lev][0] || m_lsm_data_lev[lev][0]) {
            // Valid SST, TSK or LSM data; t_surf set before computing fluxes (avoids extended lambda capture)
            // Note that land temp will be set from m_tsk_lev
            //      while sea temp will be set from m_sst_lev
            theta_type = ThetaCalcType::SURFACE_TEMPERATURE;
        }
    }

    void
    update_fluxes (const int& lev,
                   const amrex::Real& time,
                   int max_iters = 500);

    template <typename FluxIter>
    void
    compute_fluxes (const int& lev,
                    const int& max_iters,
                    const FluxIter& most_flux,
                    bool is_land);

    void
    impose_SurfaceLayer_bcs (const int& lev,
                             amrex::Vector<const amrex::MultiFab*> mfs,
                             amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                             amrex::MultiFab* xheat_flux,
                             amrex::MultiFab* yheat_flux,
                             amrex::MultiFab* zheat_flux,
                             amrex::MultiFab* xqv_flux,
                             amrex::MultiFab* yqv_flux,
                             amrex::MultiFab* zqv_flux,
                             const amrex::MultiFab* z_phys);

    template <typename FluxCalc>
    void
    compute_SurfaceLayer_bcs (const int& lev,
                              amrex::Vector<const amrex::MultiFab*> mfs,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                              amrex::MultiFab* xheat_flux,
                              amrex::MultiFab* yheat_flux,
                              amrex::MultiFab* zheat_flux,
                              amrex::MultiFab* xqv_flux,
                              amrex::MultiFab* yqv_flux,
                              amrex::MultiFab* zqv_flux,
                              const amrex::MultiFab* z_phys,
                              const FluxCalc& flux_comp);

    void
    fill_tsurf_with_sst_and_tsk (const int& lev,
                                 const amrex::Real& time);

    void
    get_lsm_tsurf (const int& lev);

    /* wrapper around compute_pblh */
    void
    update_pblh (const int& lev,
                 amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars,
                 amrex::MultiFab* z_phys_cc,
                 const int RhoQv_comp,
                 const int RhoQc_comp,
                 const int RhoQr_comp);

    template <typename PBLHeightEstimator>
    void
    compute_pblh (const int& lev,
                  amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars,
                  amrex::MultiFab* z_phys_cc,
                  const PBLHeightEstimator& est,
                  const int RhoQv_comp,
                  const int RhoQc_comp,
                  const int RhoQr_comp);

    void
    read_custom_roughness (const int& lev,
                           const std::string& fname);

    void
    update_surf_temp (const amrex::Real& time)
    {
        if (surf_heating_rate != 0) {
            int nlevs = m_geom.size();
            for (int lev = 0; lev < nlevs; lev++)
            {
               t_surf[lev]->setVal(surf_temp + surf_heating_rate*time);
               amrex::Print() << "Surface temp at t=" << time
                   << ": "
                   << surf_temp + surf_heating_rate*time
                   << std::endl;
            }
        }
    }

    void
    update_mac_ptrs (const int& lev,
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qv_prim,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qr_prim)
    { m_ma.update_field_ptrs(lev,vars_old,Theta_prim,Qv_prim,Qr_prim); }

    amrex::MultiFab*
    get_u_star (const int& lev) { return u_star[lev].get(); }

    amrex::MultiFab*
    get_w_star (const int& lev) { return w_star[lev].get(); }

    amrex::MultiFab*
    get_t_star (const int& lev) { return t_star[lev].get(); }

    amrex::MultiFab*
    get_q_star (const int& lev) { return q_star[lev].get(); }

    amrex::MultiFab*
    get_olen (const int& lev) { return olen[lev].get(); }

    amrex::MultiFab*
    get_pblh (const int& lev) { return pblh[lev].get(); }

    const amrex::MultiFab*
    get_mac_avg (const int& lev, int comp) { return m_ma.get_average(lev,comp); }

    amrex::MultiFab*
    get_t_surf (const int& lev) { return t_surf[lev].get(); }

    amrex::Real
    get_zref () { return m_ma.get_zref(); }

    amrex::FArrayBox*
    get_z0 (const int& lev) { return &z_0[lev]; }

    bool have_variable_sea_roughness() { return m_var_z0; }

    amrex::iMultiFab*
    get_lmask(const int& lev) { return m_lmask_lev[lev][0]; }

    int
    lmask_min_reduce (amrex::iMultiFab& lmask, const int& nghost)
    {
        int lmask_min = amrex::ReduceMin(lmask, nghost,
                                         [=] AMREX_GPU_HOST_DEVICE (amrex::Box const& bx, amrex::Array4<int const> const& lm_arr) -> int
        {
            int locmin = std::numeric_limits<int>::max();
            const auto lo = lbound(bx);
            const auto hi = ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    locmin = std::min(locmin, lm_arr(i,j,0));
                }
            }
            return locmin;
        });

        return lmask_min;
    }

    enum struct FluxCalcType {
        MOENG = 0,      ///< Moeng functional form
        DONELAN,        ///< Donelan functional form
        CUSTOM,         ///< Custom constant flux functional form
        ROTATE          ///< Terrain rotation flux functional form
    };

    enum struct ThetaCalcType {
        ADIABATIC = 0,
        HEAT_FLUX,          ///< Heat-flux specified
        SURFACE_TEMPERATURE ///< Surface temperature specified
    };

    enum struct RoughCalcType {
        CONSTANT = 0,      ///< Constant z0
        CHARNOCK,
        MODIFIED_CHARNOCK,
        DONELAN,
        WAVE_COUPLED
    };

    enum struct PBLHeightCalcType {
        None, MYNN25, YSU
    };

    FluxCalcType   flux_type{FluxCalcType::MOENG};
    ThetaCalcType theta_type{ThetaCalcType::ADIABATIC};
    RoughCalcType rough_type_land{RoughCalcType::CONSTANT};
    RoughCalcType rough_type_sea{RoughCalcType::CHARNOCK};
    PBLHeightCalcType pblh_type{PBLHeightCalcType::None};

private:
    // Set in constructor
    amrex::Vector<amrex::Geometry>  m_geom;
    bool m_rotate   = false;
    amrex::Real m_start_bdy_time;
    amrex::Real m_bdy_time_interval;

    bool m_include_wstar = false;
    amrex::Real z0_const{0.1};
    amrex::Real default_land_surf_temp = 300.;
    amrex::Real surf_temp;
    amrex::Real surf_heating_rate{0};
    amrex::Real surf_temp_flux{0};
    amrex::Real custom_ustar{0};
    amrex::Real custom_tstar{0};
    amrex::Real custom_qstar{0};
    amrex::Real cnk_a{0.0185};
    bool        cnk_visc{false};
    amrex::Real depth{30.0};
    amrex::Vector<amrex::FArrayBox> z_0;
    bool m_var_z0{false};

    bool use_moisture;

    MOSTAverage m_ma;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> u_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> w_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> t_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> q_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> olen;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> pblh;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> t_surf;

    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_sst_lev;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_tsk_lev;
    amrex::Vector<amrex::Vector<amrex::iMultiFab*>> m_lmask_lev;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_lsm_data_lev;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_lsm_flux_lev;
    amrex::Vector<amrex::MultiFab*>  m_Hwave_lev;
    amrex::Vector<amrex::MultiFab*>  m_Lwave_lev;
    amrex::Vector<amrex::MultiFab*>  m_eddyDiffs_lev;
};

#endif /* SURFACELAYER_H */
