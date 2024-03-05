#ifndef LAGRANGIANMICROPHYSICS_H
#define LAGRANGIANMICROPHYSICS_H

#ifdef ERF_USE_PARTICLES

#include <utility>
#include <string>

#include "NullMoistLagrangian.H"
#include "Microphysics.H"

/* forward declaration */
class ERFPC;

class LagrangianMicrophysics : public Microphysics {

public:

    LagrangianMicrophysics () { }

    ~LagrangianMicrophysics () = default;

    LagrangianMicrophysics (  const int& /* nlev */,
                              const MoistureType& a_model_type )
    {
        AMREX_ASSERT( Microphysics::modelType(a_model_type) == MoistureModelType::Lagrangian );
        amrex::Abort("No Lagrangian moisture model implemented yet!") ;
    }

    void Define (const int& lev,
                 SolverChoice& sc) override
    {
        if (lev > 0) return;
        m_moist_model->Define(sc);
    }

    void Init (const int& lev,
               const amrex::MultiFab& cons_in,
               const amrex::BoxArray& grids,
               const amrex::Geometry& geom,
               const amrex::Real& dt_advance) override
    {
        if (lev > 0) return;
        m_moist_model->Init(cons_in, grids, geom, dt_advance);
    }

    void Advance ( const int& lev,
                   const amrex::Real& dt_advance,
                   const SolverChoice&,
                   amrex::Vector<amrex::Vector<amrex::MultiFab>>& a_vars,
                   const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z) override
    {
        if (lev > 0) return;
        m_moist_model->Advance(dt_advance, a_vars, a_z);
    }

    void Diagnose (const int& lev) override
    {
        if (lev > 0) return;
        m_moist_model->Diagnose();
    }

    void Update_Micro_Vars_Lev (const int& lev, amrex::MultiFab& cons_in) override
    {
        if (lev > 0) return;
        m_moist_model->Update_Micro_Vars(cons_in);
    }

    void Update_State_Vars_Lev (const int& lev, amrex::MultiFab& cons_in) override
    {
        if (lev > 0) return;
        m_moist_model->Update_State_Vars(cons_in);
    }

    amrex::MultiFab* Get_Qmoist_Ptr (const int& lev, const int& varIdx) override
    {
        return (lev > 0 ? nullptr : m_moist_model->Qmoist_Ptr(varIdx));
    }

    int Get_Qmoist_Size (const int& a_lev) override
    {
        return (a_lev > 0 ? 0 : m_moist_model->Qmoist_Size());
    }

    int Get_Qstate_Size () override
    {
        return m_moist_model->Qstate_Size();
    }

    inline ERFPC* getParticleContainer() const
    {
        return m_moist_model->getParticleContainer();
    }

    inline const std::string& getName() const
    {
        return m_moist_model->getName();
    }

protected:

    template<class NewMoistModel>
    void
    SetModel ()
    {
        m_moist_model = std::make_unique<NewMoistModel>();
    }

    std::unique_ptr<NullMoistLagrangian> m_moist_model;
};

#endif
#endif
