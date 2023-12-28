#ifndef MICROPHYSICS_H
#define MICROPHYSICS_H

#include <NullMoist.H>
#include <SAM.H>
#include <Kessler.H>
#include <FastEddy.H>

class Microphysics {

public:

    Microphysics () { }

    ~Microphysics () = default;

    void
    ReSize (const int& nlev) { m_moist_model.resize(nlev); }

    template<class NewMoistModel>
    void
    SetModel ()
    {
        for (int lev(0); lev<m_moist_model.size(); ++lev) {
            m_moist_model[lev] = std::make_unique<NewMoistModel>();
        }
    }

    void
    Define (const int& lev,
            SolverChoice& sc)
    {
        m_moist_model[lev]->Define(sc);
    }

    void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance)
    {
        m_moist_model[lev]->Init(cons_in, grids, geom, dt_advance);
    }

    void
    Advance (const int& lev, const amrex::Real& dt_advance)
    {
        m_moist_model[lev]->Advance(dt_advance);
    }

    void
    Diagnose (const int& lev)
    {
        m_moist_model[lev]->Diagnose();
    }

    void
    Update (const int& lev,
            amrex::MultiFab& cons_in,
            amrex::MultiFab& qmoist)
    {
        m_moist_model[lev]->Update(cons_in, qmoist);
    }

    void
    Update_Micro_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        m_moist_model[lev]->Update_Micro_Vars(cons_in);
    }

    void
    Update_State_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        m_moist_model[lev]->Update_State_Vars(cons_in);
    }

    amrex::MultiFab*
    Get_Qmoist_Ptr (const int& lev, const int& varIdx) { return m_moist_model[lev]->Qmoist_Ptr(varIdx); }

    int
    Get_Qmoist_Size () { return m_moist_model[0]->Qmoist_Size(); }

private:
    amrex::Vector<std::unique_ptr<NullMoist>> m_moist_model;
};
#endif
