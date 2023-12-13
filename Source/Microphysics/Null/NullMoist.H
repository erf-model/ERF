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
    define (SolverChoice& /*sc*/) { }

    virtual
    void Init (const amrex::MultiFab& /*cons_in*/,
                     amrex::MultiFab& /*qmoist*/,
               const amrex::BoxArray& /*grids*/,
               const amrex::Geometry& /*geom*/,
               const amrex::Real& /*dt_advance*/) { }

    virtual
    void
    Advance ( ) { }

    virtual
    void
    Update (amrex::MultiFab& /*cons_in*/,
            amrex::MultiFab& /*qmoist*/) { }

    virtual
    void
    Update_Micro_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Update_State_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Diagnose ( ) { }

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
    Qmoist_Size ( ) { return NullMoist::m_qmoist_size; }

 private:
    int m_qmoist_size = 1;
};
#endif
