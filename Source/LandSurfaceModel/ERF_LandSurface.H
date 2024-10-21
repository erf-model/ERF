#ifndef ERF_LANDSURFACE_H
#define ERF_LANDSURFACE_H

#include <AMReX_PlotFileUtil.H>

#include <ERF_NullSurf.H>
#include <ERF_SLM.H>
#include <ERF_MM5.H>

#if ERF_USE_NOAH
#include <NOAH.H>
#endif

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
          const amrex::Real& dt_advance)
    {
        m_lsm_model[lev]->Init(cons_in, geom, dt_advance);
    }

    void
    Advance (const int& lev, const amrex::Real& dt_advance)
    {
        m_lsm_model[lev]->Advance(dt_advance);
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
    Get_Flux_Ptr (const int& lev, const int& varIdx) { return m_lsm_model[lev]->Lsm_Flux_Ptr(varIdx); }

    amrex::Geometry
    Get_Lsm_Geom (const int& lev ) { return m_lsm_model[lev]->Lsm_Geom(); }

    int
    Get_Data_Size () { return m_lsm_model[0]->Lsm_Data_Size(); }

    std::string
    Get_VarName (const int& varIdx) { return m_lsm_model[0]->Lsm_VarName(varIdx); }

    void
    Plot_Lsm_Data (amrex::Real time,
                   const amrex::Vector<int> &level_steps,
                   const amrex::Vector<amrex::IntVect> &ref_ratio)
    {
        int nlev = m_lsm_model.size();
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
                    if (lev==0) varnames[n] = this->Get_VarName(n);
                }
            }
            WriteMultiLevelPlotfile (plotfilename, nlev, GetVecOfConstPtrs(m_lsm_data_lev),
                                     varnames, m_lsm_geom_lev, time, level_steps, ref_ratio);
            m_lsm_geom_lev.clear();
            m_lsm_data_lev.clear();
        }
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
