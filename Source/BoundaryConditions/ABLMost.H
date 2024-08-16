#ifndef ABLMOST_H
#define ABLMOST_H

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <IndexDefines.H>
#include <ERF_Constants.H>
#include <MOSTAverage.H>
#include <MOSTStress.H>
#include <TerrainMetrics.H>

/** Monin-Obukhov surface layer profile
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
class ABLMost
{

public:

    // Constructor
    explicit ABLMost (const amrex::Vector<amrex::Geometry>& geom,
                      bool& use_exp_most,
                      bool& use_rot_most,
                      amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qv_prim,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qr_prim,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd,
                      amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>& sst_lev,
                      amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>>& lmask_lev,
                      amrex::Vector<amrex::Vector<amrex::MultiFab*>> lsm_data,
                      amrex::Vector<amrex::Vector<amrex::MultiFab*>> lsm_flux,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Hwave,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Lwave,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& eddyDiffs,
                      amrex::Real start_bdy_time = 0.0,
                      amrex::Real bdy_time_interval = 0.0)
    : m_exp_most(use_exp_most),
      m_rotate(use_rot_most),
      m_start_bdy_time(start_bdy_time),
      m_bdy_time_interval(bdy_time_interval),
      m_geom(geom),
      m_ma(geom,vars_old,Theta_prim,Qv_prim,Qr_prim,z_phys_nd)
    {
        // We have a moisture model if Qv_prim is a valid pointer
        use_moisture = (Qv_prim[0].get());

        // Get roughness
        amrex::ParmParse pp("erf");
        pp.query("most.z0", z0_const);

        // Specify how to compute the flux
        if (use_rot_most) {
            flux_type = FluxCalcType::ROTATE;
        } else {
            std::string flux_string{"moeng"};
            pp.query("most.flux_type", flux_string);
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

        // Get surface temperature
        auto erf_st = pp.query("most.surf_temp", surf_temp);

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

        // Specify how to compute the surface flux over land (if there is any)
        std::string rough_land_string{"constant"};
        pp.query("most.roughness_type_land", rough_land_string);
        if (rough_land_string == "constant") {
            rough_type_land = RoughCalcType::CONSTANT;
        } else {
            amrex::Abort("Undefined MOST roughness type for land!");
        }

        // Specify how to compute the surface flux over sea (if there is any)
        std::string rough_sea_string{"charnock"};
        pp.query("most.roughness_type_sea", rough_sea_string);
        if (rough_sea_string == "charnock") {
            rough_type_sea = RoughCalcType::CHARNOCK;
            pp.query("most.charnock_constant",cnk_a);
        } else if (rough_sea_string == "modified_charnock") {
            rough_type_sea = RoughCalcType::MODIFIED_CHARNOCK;
            pp.query("most.modified_charnock_depth",depth);
        } else if (rough_sea_string == "wave_coupled") {
            rough_type_sea = RoughCalcType::WAVE_COUPLED;
        } else {
            amrex::Abort("Undefined MOST roughness type for sea!");
        }

        // Size the MOST params for all levels
        int nlevs = m_geom.size();
        z_0.resize(nlevs);
        u_star.resize(nlevs);
        t_star.resize(nlevs);
        q_star.resize(nlevs);
        t_surf.resize(nlevs);
        olen.resize(nlevs);

        // Get pointers to SST and LANDMASK data
        m_sst_lev.resize(nlevs);
        m_lmask_lev.resize(nlevs);

        for (int lev(0); lev<nlevs; ++lev) {
            int nt_tot_sst = sst_lev[lev].size();
            m_sst_lev[lev].resize(nt_tot_sst);
            for (int nt(0); nt<nt_tot_sst; ++nt) {
                m_sst_lev[lev][nt]   = sst_lev[lev][nt].get();
            }
            int nt_tot_lmask = lmask_lev[lev].size();
            m_lmask_lev[lev].resize(nt_tot_lmask);
            for (int nt(0); nt<nt_tot_lmask; ++nt) {
                m_lmask_lev[lev][nt] = lmask_lev[lev][nt].get();
            }
        } // lev

        // Get pointers to LSM data and Fluxes
        m_lsm_data_lev.resize(nlevs);
        m_lsm_flux_lev.resize(nlevs);
        for (int lev(0); lev<nlevs; ++lev) {
            int nvar = lsm_data[lev].size();
            m_lsm_data_lev[lev].resize(nvar);
            m_lsm_flux_lev[lev].resize(nvar);
            for (int n(0); n<nvar; ++n) {
                m_lsm_data_lev[lev][n] = lsm_data[lev][n];
                m_lsm_flux_lev[lev][n] = lsm_flux[lev][n];
            }
        } // lev

        // Get pointers to wave data
        m_Hwave_lev.resize(nlevs);
        m_Lwave_lev.resize(nlevs);
        m_eddyDiffs_lev.resize(nlevs);
        for (int lev(0); lev<nlevs; ++lev) {
            m_Hwave_lev[lev] = Hwave[lev].get();
            m_Lwave_lev[lev] = Lwave[lev].get();
            m_eddyDiffs_lev[lev] = eddyDiffs[lev].get();
        }

        for (int lev = 0; lev < nlevs; lev++) {
            // Attributes for MFs and FABs
            //--------------------------------------------------------
            auto& mf = vars_old[lev][Vars::cons];
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
            amrex::Box bx = amrex::grow(m_geom[lev].Domain(),ng);
            bx.setSmall(2,0);
            bx.setBig(2,0);
            z_0[lev].resize(bx,1);
            z_0[lev].setVal<amrex::RunOn::Device>(z0_const);

            // 2D MFs for U*, T*, T_surf
            //--------------------------------------------------------
            u_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
            u_star[lev]->setVal(1.E34);

            t_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
            t_star[lev]->setVal(1.E34);

            q_star[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
            q_star[lev]->setVal(0.0); // default to dry

            olen[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);
            olen[lev]->setVal(1.E34);

            t_surf[lev] = std::make_unique<amrex::MultiFab>(ba2d,dm,ncomp,ng);

            // TODO: Do we want lsm_data to have theta at 0 index always?
            //       Do we want an external enum struct for indexing?
            if (m_sst_lev[lev][0] || m_lsm_data_lev[lev][0]) {
                // Valid SST or LSM data; t_surf set before computing fluxes (avoids extended lambda capture)
                theta_type = ThetaCalcType::SURFACE_TEMPERATURE;
            } else if (erf_st) {
                t_surf[lev]->setVal(surf_temp);
            } else {
                t_surf[lev]->setVal(0.0);
            }
        }// lev
    }

    void
    update_fluxes (const int& lev,
                   const amrex::Real& time,
                   int max_iters = 25);

    template <typename FluxIter>
    void
    compute_fluxes (const int& lev,
                    const int& max_iters,
                    const FluxIter& most_flux,
                    bool is_land);

    void
    impose_most_bcs (const int& lev,
                     const amrex::Vector<amrex::MultiFab*>& mfs,
                     amrex::MultiFab* xxmom_flux,
                     amrex::MultiFab* yymom_flux,
                     amrex::MultiFab* zzmom_flux,
                     amrex::MultiFab* xymom_flux, amrex::MultiFab* yxmom_flux,
                     amrex::MultiFab* xzmom_flux, amrex::MultiFab* zxmom_flux,
                     amrex::MultiFab* yzmom_flux, amrex::MultiFab* zymom_flux,
                     amrex::MultiFab* xheat_flux,
                     amrex::MultiFab* yheat_flux,
                     amrex::MultiFab* zheat_flux,
                     amrex::MultiFab* xqv_flux,
                     amrex::MultiFab* yqv_flux,
                     amrex::MultiFab* zqv_flux,
                     amrex::MultiFab* z_phys);

    template<typename FluxCalc>
    void
    compute_most_bcs (const int& lev,
                      const amrex::Vector<amrex::MultiFab*>& mfs,
                      amrex::MultiFab* xxmom_flux,
                      amrex::MultiFab* yymom_flux,
                      amrex::MultiFab* zzmom_flux,
                      amrex::MultiFab* xymom_flux, amrex::MultiFab* yxmom_flux,
                      amrex::MultiFab* xzmom_flux, amrex::MultiFab* zxmom_flux,
                      amrex::MultiFab* yzmom_flux, amrex::MultiFab* zymom_flux,
                      amrex::MultiFab* xheat_flux,
                      amrex::MultiFab* yheat_flux,
                      amrex::MultiFab* zheat_flux,
                      amrex::MultiFab* xqv_flux,
                      amrex::MultiFab* yqv_flux,
                      amrex::MultiFab* zqv_flux,
                      amrex::MultiFab* z_phys,
                      const FluxCalc& flux_comp);

    void
    time_interp_sst (const int& lev,
                     const amrex::Real& time);

    void
    get_lsm_tsurf (const int& lev);

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

    const amrex::MultiFab*
    get_u_star (const int& lev) { return u_star[lev].get(); }

    const amrex::MultiFab*
    get_t_star (const int& lev) { return t_star[lev].get(); }

    const amrex::MultiFab*
    get_q_star (const int& lev) { return q_star[lev].get(); }

    const amrex::MultiFab*
    get_olen (const int& lev) { return olen[lev].get(); }

    const amrex::MultiFab*
    get_mac_avg (const int& lev, int comp) { return m_ma.get_average(lev,comp); }

    const amrex::MultiFab*
    get_t_surf (const int& lev) { return t_surf[lev].get(); }

    amrex::Real
    get_zref () {return m_ma.get_zref();}

    const amrex::FArrayBox*
    get_z0 (const int& lev) {return &z_0[lev];}

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
        WAVE_COUPLED
    };

    FluxCalcType   flux_type{FluxCalcType::MOENG};
    ThetaCalcType theta_type{ThetaCalcType::ADIABATIC};
    RoughCalcType rough_type_land{RoughCalcType::CONSTANT};
    RoughCalcType rough_type_sea{RoughCalcType::CHARNOCK};

private:
    bool use_moisture;
    bool m_exp_most = false;
    bool m_rotate   = false;
    amrex::Real z0_const{0.1};
    amrex::Real surf_temp;
    amrex::Real surf_heating_rate{0};
    amrex::Real surf_temp_flux{0};
    amrex::Real custom_ustar{0};
    amrex::Real custom_tstar{0};
    amrex::Real custom_qstar{0};
    amrex::Real cnk_a{0.0185};
    amrex::Real depth{30.0};
    amrex::Real m_start_bdy_time;
    amrex::Real m_bdy_time_interval;
    amrex::Vector<amrex::Geometry>  m_geom;
    amrex::Vector<amrex::FArrayBox> z_0;

    MOSTAverage m_ma;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> u_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> t_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> q_star;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> olen;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> t_surf;

    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_sst_lev;
    amrex::Vector<amrex::Vector<amrex::iMultiFab*>> m_lmask_lev;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_lsm_data_lev;
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_lsm_flux_lev;
    amrex::Vector<amrex::MultiFab*>  m_Hwave_lev;
    amrex::Vector<amrex::MultiFab*>  m_Lwave_lev;
    amrex::Vector<amrex::MultiFab*>  m_eddyDiffs_lev;
};

#endif /* ABLMOST_H */
