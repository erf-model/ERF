#ifndef ERF_URBAN_H
#define ERF_URBAN_H

#include <memory>
#include <string>

#include <ERF_NullUrban.H>

/**
 * @brief Wrapper class for managing urban models across AMR levels.
 */
class Urban {

public:

    Urban () = default;

    ~Urban () = default;

    /**
     * @brief Resizes urban-model storage and records active AMR levels.
     */
    void ReSize (const int& nlev, const amrex::Vector<int>& active_levs)
    {
        AMREX_ALWAYS_ASSERT(active_levs.size() == nlev);
        m_urban_model.resize(nlev);
        m_active_levs = active_levs;
    }

    /**
     * @brief Instantiates an urban model on each active AMR level.
     */
    template<class NewUrbanModel>
    void SetModel ()
    {
        m_first_active_lev = -1;
        for (int lev = 0; lev < m_urban_model.size(); ++lev) {
            if (m_active_levs[lev] == 1) {
                m_urban_model[lev] = std::make_unique<NewUrbanModel>();
                if (m_first_active_lev == -1) {
                    m_first_active_lev = lev;
                }
            }
        }
    }

    /**
     * @brief Defines the urban model on an active level.
     */
    void Define (const int& lev, SolverChoice& sc)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->Define(sc);
        }
    }

    /**
     * @brief Initializes the urban model on an active level.
     */
    void Init (const int& lev, const amrex::MultiFab& cons_in,
               const amrex::MultiFab& u_in, const amrex::MultiFab& v_in,
               const amrex::Geometry& geom, const amrex::Real& dt_advance,
               const amrex::MultiFab& z_phys_nd, const amrex::iMultiFab& land_type,
               const amrex::MultiFab& urb_frac)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->Init(cons_in, u_in, v_in, geom, dt_advance,
                                      z_phys_nd, land_type, urb_frac);
        }
    }

    /**
     * @brief Advances the urban model on an active level.
     */
    void Advance (const int& lev, const amrex::Real& dt_advance,
                  const amrex::Real& time, const amrex::Real& start_time,
                  const amrex::Geometry& geom, const amrex::MultiFab* z_phys_nd,
                  amrex::MultiFab& eddyDiffs_lev, const amrex::Real& calday)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->Advance(dt_advance, time, start_time, geom,
                                        z_phys_nd, eddyDiffs_lev, calday);
        }
    }

    /**
     * @brief Updates coupled urban fields on an active level.
     */
    void Update_Urban_Vars_Lev (const int& lev, amrex::MultiFab& cons_in,
                                amrex::MultiFab& u_in, amrex::MultiFab& v_in)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->Update_Urban_Vars(cons_in, u_in, v_in);
        }
    }

    /**
     * @brief Updates atmospheric state fields on an active level.
     */
    void Update_State_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->Update_State_Vars(cons_in);
        }
    }

    /**
     * @brief Returns a model data field, or null for an inactive level.
     */
    amrex::MultiFab* Get_Data_Ptr (const int& lev, const int& varIdx)
    {
        return m_active_levs[lev] ? m_urban_model[lev]->Urban_Data_Ptr(varIdx) : nullptr;
    }

    /**
     * @brief Returns a model flux field, or null for an inactive level.
     */
    amrex::MultiFab* Get_Flux_Ptr (const int& lev, const int& varIdx)
    {
        return m_active_levs[lev] ? m_urban_model[lev]->Urban_Flux_Ptr(varIdx) : nullptr;
    }

    /**
     * @brief Returns the urban-model geometry for a level.
     */
    amrex::Geometry Get_Urban_Geom (const int& lev)
    {
        return m_active_levs[lev] ? m_urban_model[lev]->Urban_Geom() : amrex::Geometry();
    }

    /**
     * @brief Returns the number of model data fields.
     */
    int Get_Data_Size ()
    {
        return (m_first_active_lev != -1)
            ? m_urban_model[m_first_active_lev]->Urban_Data_Size() : 0;
    }

    /**
     * @brief Returns a model data-field name.
     */
    std::string Get_DataName (const int& varIdx)
    {
        return (m_first_active_lev != -1)
            ? m_urban_model[m_first_active_lev]->Urban_DataName(varIdx) : "";
    }

    /**
     * @brief Returns a model data-field index.
     */
    int Get_DataIdx (const int& lev, const std::string& varname)
    {
        return m_active_levs[lev] ? m_urban_model[lev]->Urban_DataIndex(varname) : -1;
    }

    /**
     * @brief Sets flux inputs on an active level.
     */
    void set_urban_flux_inputs (const int& lev, const amrex::Real& declin)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->set_flux_inputs(declin);
        }
    }

    /**
     * @brief Sets terrain inputs on an active level.
     */
    void set_urban_terrain_inputs (
        const int& lev,
        const amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>& tsk_in,
        const amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>>& lmask_in,
        const amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>>& land_type_in,
        const amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>& urb_frac_in,
        amrex::MultiFab* lat_in, amrex::MultiFab* lon_in)
    {
        if (m_active_levs[lev] == 1) {
            m_urban_model[lev]->set_terrain_inputs(
                tsk_in[lev], lmask_in[lev], land_type_in[lev], urb_frac_in[lev],
                lat_in, lon_in);
        }
    }

private:
    amrex::Vector<std::unique_ptr<NullUrban>> m_urban_model;
    amrex::Vector<int> m_active_levs;
    int m_first_active_lev = -1;
};

#endif
