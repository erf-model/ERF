#ifndef ERF_LANDSURFACE_H
#define ERF_LANDSURFACE_H

#include <typeinfo>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_NullSurf.H>
#include <ERF_SLM.H>
#include <ERF_OceanSurf.H>

#if ERF_USE_NOAHMP
#include <ERF_NOAHMP.H>
#endif

/**
 * @brief Wrapper class for managing land surface models across AMR levels.
 */
class LandSurface {

public:

    LandSurface () { }

    ~LandSurface () = default;

    void
    ReSize (const int& nlev) { m_lsm_model.resize(nlev); }

    template<class NewSurfModel>
    void
    SetModel ()
    {
        for (int lev(0); lev<m_lsm_model.size(); ++lev) {
            m_lsm_model[lev] = std::make_unique<NewSurfModel>();
        }
    }

    void
    Define (const int& lev,
            SolverChoice& sc)
    {
        m_lsm_model[lev]->Define(sc);
    }

    void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Geometry& geom0,
          amrex::Vector<amrex::BCRec>& domain_bcs_type,
          amrex::IntVect& refRatio,
          const amrex::Real& dt_advance,
          amrex::Vector<amrex::Vector<std::string>>& nc_init_file)
    {
        m_lsm_model[lev]->Init(lev, cons_in, geom, geom0,
                               domain_bcs_type, refRatio,
                               dt_advance, nc_init_file);
    }

    void
    Advance (const int& lev,
             amrex::MultiFab& cons_in,
             amrex::MultiFab& xvel_in,
             amrex::MultiFab& yvel_in,
             amrex::MultiFab* hfx3_out,
             amrex::MultiFab* qfx3_out,
             const SurfacePrecipAccumulationSources& precip_sources,
             const double& time_d,
             const double& dt_advance_d,
             const int& nstep,
             const bool updated_lev0)
    {
        amrex::Real time = static_cast<amrex::Real>(time_d);
        amrex::Real dt_advance = static_cast<amrex::Real>(dt_advance_d);
        m_lsm_model[lev]->Advance_With_State(lev, cons_in, xvel_in, yvel_in,
                                             hfx3_out, qfx3_out,
                                             precip_sources,
                                             time, dt_advance, nstep, updated_lev0);
    }

    void
    Advance (const int& lev, const double& dt_advance_d)
    {
        amrex::Real dt_advance = static_cast<amrex::Real>(dt_advance_d);
        m_lsm_model[lev]->Advance(dt_advance);
    }

    void
    Plot(const int& lev, const int& nstep)
    {
        m_lsm_model[lev]->Plot_Landfile(nstep);
    }

    void
    Update_Micro_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        m_lsm_model[lev]->Update_Micro_Vars(cons_in);
    }

    void
    Update_State_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        m_lsm_model[lev]->Update_State_Vars(cons_in);
    }

    amrex::MultiFab*
    Get_Data_Ptr (const int& lev, const int& varIdx) { return m_lsm_model[lev]->Lsm_Data_Ptr(varIdx); }

    amrex::MultiFab*
    Get_Data_Ptr (const int& lev, const std::string& varname)
    {
        std::string name = varname;
        const int index = m_lsm_model[lev]->Lsm_DataIndex(name);
        return (index >= 0) ? m_lsm_model[lev]->Lsm_Data_Ptr(index) : nullptr;
    }

    bool
    Has_Model () const { return !m_lsm_model.empty() && m_lsm_model[0] != nullptr; }

    amrex::MultiFab*
    Get_Flux_Ptr (const int& lev, const int& varIdx) { return m_lsm_model[lev]->Lsm_Flux_Ptr(varIdx); }

    amrex::Geometry
    Get_Lsm_Geom (const int& lev ) { return m_lsm_model[lev]->Lsm_Geom(); }

    int
    Get_Data_Size () { return m_lsm_model[0]->Lsm_Data_Size(); }

    int
    Get_Flux_Size () { return m_lsm_model[0]->Lsm_Flux_Size(); }

    std::string
    Get_DataName (const int& varIdx) { return m_lsm_model[0]->Lsm_DataName(varIdx); }

    amrex::Vector<std::string>
    Get_DataNames ()
    {
        amrex::Vector<std::string> names;
        if (m_lsm_model.empty() || m_lsm_model[0] == nullptr) {
            return names;
        }
        const int nvars = Get_Data_Size();
        names.reserve(nvars);
        for (int var = 0; var < nvars; ++var) {
            names.push_back(Get_DataName(var));
        }
        return names;
    }

    int
    Get_DataIdx (const int& lev, std::string& varname) { return m_lsm_model[lev]->Lsm_DataIndex(varname); }

    void
    Set_Lev0_Data_Ptr (const int& lev) {
        // Ensure the same type of model for copy in place of flux data
        if ( (lev > 0) && typeid(m_lsm_model[lev]) == typeid(m_lsm_model[0]) ) {
            int ndata = this->Get_Data_Size();
            for (int varIdx(0); varIdx<ndata; ++varIdx) {
                amrex::MultiFab* Lev0_Data_Ptr = this->Get_Data_Ptr(0, varIdx);
                m_lsm_model[lev]->Lsm_Set_Lev0_Data_Ptr(varIdx, Lev0_Data_Ptr);
            }
        } else {
            amrex::Abort("Differing types of land surface models across levels is not supported!");
        }
    }

    std::string
    Get_FluxName (const int& varIdx) { return m_lsm_model[0]->Lsm_FluxName(varIdx); }

    int
    Get_FluxIdx (const int& lev, std::string& varname) { return m_lsm_model[lev]->Lsm_FluxIndex(varname); }

    void
    Set_Lev0_Flux_Ptr (const int& lev) {
        // Ensure the same type of model for copy in place of flux data
        if ( (lev > 0) && typeid(m_lsm_model[lev]) == typeid(m_lsm_model[0]) ) {
            int nflux = this->Get_Flux_Size();
            for (int varIdx(0); varIdx<nflux; ++varIdx) {
                amrex::MultiFab* Lev0_Flux_Ptr = this->Get_Flux_Ptr(0, varIdx);
                m_lsm_model[lev]->Lsm_Set_Lev0_Flux_Ptr(varIdx, Lev0_Flux_Ptr);
            }
        } else {
            amrex::Abort("Differing types of land surface models across levels is not supported!");
        }
    }

    bool
    Get_LSM_Update_Status (const int& lev) { return m_lsm_model[lev]->Get_Update_Status(); }

    std::unordered_map<std::string,std::string>&
    Get_WRFInputNames() { return m_lsm_model[0]->Lsm_WRFInputNames(); }

    // Scalar LSM step counter (e.g., NoahMP itimestep) for checkpoint/restart.
    int
    Get_LSM_Step (const int& lev) const { return m_lsm_model[lev]->Get_LSM_Step(); }

    void
    Set_LSM_Step (const int& lev, int step) { m_lsm_model[lev]->Set_LSM_Step(step); }

    // Full LSM prognostic-state restart I/O (e.g., NoahMP soil/snow/canopy).
    // Write is const (called from the const WriteCheckpointFile); it only
    // serializes state, it does not mutate the ERF-side model wrapper.
    void
    Write_Lsm_Restart (const int& lev, const std::string& dir) const { m_lsm_model[lev]->Write_Lsm_Restart(dir); }

    void
    Read_Lsm_Restart (const int& lev, const std::string& dir) { m_lsm_model[lev]->Read_Lsm_Restart(dir); }

    void
    Plot_Lsm_Data (double time,
                   int finest_level,
                   const amrex::Vector<int> &level_steps,
                   const amrex::Vector<amrex::IntVect> &ref_ratio)
    {
        // Note we must use finest_level here rather than m_lsm_model.size(), which is
        //    (max_level+1) -- the models at levels that have not (yet) been created
        //    have been allocated but not Define'd or Init'ed, so their data pointers
        //    are still null
        int nlev = finest_level + 1;
        int nvar = this->Get_Data_Size();
        amrex::Vector<std::string> varnames;

        // Only write if we have valid pointers
        if (this->Get_Data_Ptr(0,0)) {
            varnames.resize(nvar);
            m_lsm_geom_lev.resize(nlev);
            m_lsm_data_lev.resize(nlev);
            std::string plotfilename = amrex::Concatenate(plot_file_lsm, level_steps[0], 5);
            for (int lev(0); lev<nlev; ++lev) {
                m_lsm_geom_lev[lev] = this->Get_Lsm_Geom(lev);

                amrex::MultiFab* mf_lsm = this->Get_Data_Ptr(lev,0);
                amrex::IntVect ng(0,0,0);
                amrex::BoxArray ba = mf_lsm->boxArray();
                amrex::DistributionMapping dm = mf_lsm->DistributionMap();
                m_lsm_data_lev[lev].define(ba, dm, nvar, ng);
                for (int n(0); n<nvar; ++n) {
                    mf_lsm = this->Get_Data_Ptr(lev,n);
                    amrex::MultiFab::Copy(m_lsm_data_lev[lev],*(mf_lsm),0,n,1,0);
                    if (lev==0) varnames[n] = this->Get_DataName(n);
                }
            }
            WriteMultiLevelPlotfile (plotfilename, nlev, GetVecOfConstPtrs(m_lsm_data_lev),
                                     varnames, m_lsm_geom_lev, static_cast<amrex::Real>(time),
                                     level_steps, ref_ratio);
            m_lsm_geom_lev.clear();
            m_lsm_data_lev.clear();
        }
    }

    template<class SurfModelType>
    SurfModelType* get_model_lev(const int &lev) {
        return static_cast<SurfModelType*>(m_lsm_model[lev].get());
    }

private:
    // lsm model at each level
    amrex::Vector<std::unique_ptr<NullSurf>> m_lsm_model;

    // plotfile prefix
    std::string plot_file_lsm {"plt_lsm_"};

    // Vector of geometry
    amrex::Vector<amrex::Geometry> m_lsm_geom_lev;

    // Vector of data pointers
    amrex::Vector<amrex::MultiFab> m_lsm_data_lev;
};
#endif
