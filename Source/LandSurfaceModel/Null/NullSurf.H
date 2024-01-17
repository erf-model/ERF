#ifndef NULLSURF_H
#define NULLSURF_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>
#include <DataStruct.H>

class NullSurf {

 public:
    NullSurf () {}

    virtual ~NullSurf () = default;

    virtual
    void
    Define (SolverChoice& /*sc*/) { }

    virtual
    void Init (const amrex::MultiFab& /*cons_in*/,
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
    Copy_State_to_Micro (const amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_Micro_to_State (amrex::MultiFab& /*cons_in*/) { }

    virtual
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    int
    Lsm_Data_Size () { return NullSurf::m_lsm_size; }

 private:
    int m_lsm_size = 1;
};
#endif
