#ifndef NULLMOIST_H
#define NULLMOIST_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>
#include <DataStruct.H>

class NullMoist {

public:
    NullMoist () {}

    virtual ~NullMoist () = default;

    virtual
    void
    Define (SolverChoice& /*sc*/) { }

    virtual
    void Init (const amrex::MultiFab& /*cons_in*/,
               const amrex::BoxArray& /*grids*/,
               const amrex::Geometry& /*geom*/,
               const amrex::Real& /*dt_advance*/) { }

    virtual
    void
    Advance (const amrex::Real& /*dt_advance*/) { }

    virtual
    void
    Update_Micro_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Update_State_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Diagnose () { }

    virtual
    void
    Copy_State_to_Micro (const amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_Micro_to_State (amrex::MultiFab& /*cons_in*/) { }

    virtual
    amrex::MultiFab*
    Qmoist_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    int
    Qmoist_Size () { return NullMoist::m_qmoist_size; }

    virtual
    int
    Qstate_Size () { return NullMoist::m_qstate_size; }

private:
    int m_qmoist_size = 1;
    int m_qstate_size = 0;
};

#ifdef ERF_USE_PARTICLES

/* forward declaration */
class ERFPC;

class NullMoistLagrangian : public NullMoist {

public:
    NullMoistLagrangian() {}

    virtual ~NullMoistLagrangian() = default;

    virtual ERFPC* getParticleContainer()
    {
        return nullptr;
    }

    virtual const std::string& getName() const
    {
        return m_name;
    }

    using NullMoist::Advance;

    virtual void
    Advance (const amrex::Real&,
             amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ) { }

protected:

private:
    const std::string m_name = "null";
};

#endif
#endif
