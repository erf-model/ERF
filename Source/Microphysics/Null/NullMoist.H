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
    Diagnose ( ) { }

    virtual
    void
    Update (amrex::MultiFab& /*cons_in*/,
            amrex::MultiFab& /*qmoist*/) { }

    virtual
    int
    Qmoist_Size ( ) { return NullMoist::m_qmoist_size; }

 private:
    int m_qmoist_size = 1;
};
#endif
