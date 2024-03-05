#ifndef EULERIANMICROPHYSICS_H
#define EULERIANMICROPHYSICS_H

#include "NullMoist.H"
#include "SAM.H"
#include "Kessler.H"
#include "FastEddy.H"
#include "Microphysics.H"

class EulerianMicrophysics : public Microphysics {

public:

    EulerianMicrophysics () { }

    ~EulerianMicrophysics () = default;

    EulerianMicrophysics ( const int& nlev,
                           const MoistureType& a_model_type )
    {
        AMREX_ASSERT( Microphysics::modelType(a_model_type) == MoistureModelType::Eulerian );
        m_moist_model.resize(nlev);
        if (a_model_type == MoistureType::SAM) {
            SetModel<SAM>();
            amrex::Print() << "SAM moisture model!\n";
        } else if (a_model_type == MoistureType::Kessler or
                   a_model_type == MoistureType::Kessler_NoRain) {
            SetModel<Kessler>();
            amrex::Print() << "Kessler moisture model!\n";
        } else if (a_model_type == MoistureType::FastEddy) {
            SetModel<FastEddy>();
            amrex::Print() << "FastEddy moisture model!\n";
        } else if (a_model_type == MoistureType::None) {
            SetModel<NullMoist>();
            amrex::Print() << "WARNING: Compiled with moisture but using NullMoist model!\n";
        } else {
            amrex::Abort("EulerianMicrophysics: Dont know this moisture_type!") ;
        }
    }

    void
    Define (const int& lev,
            SolverChoice& sc) override
    {
        m_moist_model[lev]->Define(sc);
    }

    void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance) override
    {
        m_moist_model[lev]->Init(cons_in, grids, geom, dt_advance);
    }

    void
    Advance (   const int& lev,
                const amrex::Real& dt_advance,
                const SolverChoice &solverChoice,
                amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ) override
    {
        m_moist_model[lev]->Advance(dt_advance, solverChoice);
    }

    void
    Diagnose (const int& lev) override
    {
        m_moist_model[lev]->Diagnose();
    }

    void
    Update_Micro_Vars_Lev (const int& lev, amrex::MultiFab& cons_in) override
    {
        m_moist_model[lev]->Update_Micro_Vars(cons_in);
    }

    void
    Update_State_Vars_Lev (const int& lev, amrex::MultiFab& cons_in) override
    {
        m_moist_model[lev]->Update_State_Vars(cons_in);
    }

    amrex::MultiFab*
    Get_Qmoist_Ptr (const int& lev, const int& varIdx) override
    {
        return m_moist_model[lev]->Qmoist_Ptr(varIdx);
    }

    int
    Get_Qmoist_Size (const int& /* lev */) override
    {
        return m_moist_model[0]->Qmoist_Size();
    }

    int
    Get_Qstate_Size () override
    {
        return m_moist_model[0]->Qstate_Size();
    }

protected:

    template<class NewMoistModel>
    void
    SetModel ()
    {
        for (int lev(0); lev<m_moist_model.size(); ++lev) {
            m_moist_model[lev] = std::make_unique<NewMoistModel>();
        }
    }

private:
    amrex::Vector<std::unique_ptr<NullMoist>> m_moist_model;
};
#endif
