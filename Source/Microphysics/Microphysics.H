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

    template<class NewMoistModel>
    void
    SetModel ()
    {
        m_moist_model = std::make_unique<NewMoistModel>();
    }

    void
    define (SolverChoice& sc)
    {
        m_moist_model->define(sc);
    }

    void
    Init (const amrex::MultiFab& cons_in,
                amrex::MultiFab& qmoist,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance)
    {
        m_moist_model->Init(cons_in, qmoist, grids, geom, dt_advance);
    }

    void
    Advance ( )
    {
        m_moist_model->Advance();
    }

    void
    Diagnose ( )
    {
        m_moist_model->Diagnose();
    }

    void
    Update (amrex::MultiFab& cons_in,
            amrex::MultiFab& qmoist)
    {
        m_moist_model->Update(cons_in, qmoist);
    }

    int
    Get_Qmoist_Size ( ) { return m_moist_model->Qmoist_Size(); }


private:
    std::unique_ptr<NullMoist> m_moist_model;
};
#endif
