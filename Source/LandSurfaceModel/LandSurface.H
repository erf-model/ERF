#ifndef LANDSURFACE_H
#define LANDSURFACE_H

#include <NullSurf.H>
#include <SLM.H>
#include <MM5.H>

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

    int
    Get_Data_Size () { return m_lsm_model[0]->Lsm_Data_Size(); }

private:
    amrex::Vector<std::unique_ptr<NullSurf>> m_lsm_model;
};
#endif
