#ifndef MICROPHYSICS_H
#define MICROPHYSICS_H

#include <NullMoist.H>
#include <SAM.H>

template<class MoistModel>
class Microphysics : public MoistModel {

public:

    Microphysics () { }

   ~Microphysics () = default;

    template<class NewMoistModel>
    void
    SetModel (NewMoistModel& T )
    {
        m_moist_model.reset(new NewMoistModel());
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
    Update (amrex::MultiFab& cons_in,
            amrex::MultiFab& qmoist)
    {
        m_moist_model->Update(cons_in, qmoist);
    }


private:
    std::unique_ptr<MoistModel> m_moist_model;
};
#endif


/*
#ifndef MICROPHYSICS_H
#define MICROPHYSICS_H

#include <NullMoist.H>
#include <SAM.H>

// Available models
namespace MicMod {
   enum {
      Null = 0,
      SAM,
      Kessler,
      NumVars
  };
}

class Microphysics {

public:

    Microphysics ( ) { }

   ~Microphysics ( ) = default;

    void
    SetModel (std::string moisture_model)
    {
        if ( moisture_model == "SAM" ) {
            m_MicMod = MicMod::SAM;
            m_sam_model = std::make_unique<SAM>();
        } else if ( moisture_model == "SAM" ) {
            // TODO: Code block
        }else {
            amrex::Print() << "WARNING: Compiled with moisture support but no valid moisture model found\n";
            m_MicMod = MicMod::Null;
            m_null_model = std::make_unique<NullMoist>();
        }
    }

    void
    define (SolverChoice& sc)
    {
        switch(m_MicMod) {
        case MicMod::SAM:
            m_sam_model->define(sc);
            break;
        case MicMod::Kessler:
            // TODO: code block
            break;
        default:
            m_null_model->define(sc);
        }
    }

    void
    Init (const amrex::MultiFab& cons_in,
                amrex::MultiFab& qmoist,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance)
    {
        switch(m_MicMod) {
        case MicMod::SAM:
            m_sam_model->Init(cons_in, qmoist, grids, geom, dt_advance);
            break;
        case MicMod::Kessler:
            // TODO: code block
            break;
        default:
            m_null_model->Init(cons_in, qmoist, grids, geom, dt_advance);
        }
    }

    void
    Advance ( )
    {
        switch(m_MicMod) {
        case MicMod::SAM:
            m_sam_model->Advance();
            break;
        case MicMod::Kessler:
            // TODO: code block
            break;
        default:
            m_null_model->Advance();
        }
    }

    void
    Update (amrex::MultiFab& cons_in,
            amrex::MultiFab& qmoist)
    {
        switch(m_MicMod) {
        case MicMod::SAM:
            m_sam_model->Update(cons_in, qmoist);
            break;
        case MicMod::Kessler:
            // TODO: code block
            break;
        default:
            m_null_model->Update(cons_in, qmoist);
        }
    }

private:
    int m_MicMod;
    std::unique_ptr<NullMoist> m_null_model = nullptr;
    std::unique_ptr<SAM> m_sam_model = nullptr;
    //std::unique_ptr<Kessler> m_kessler_model = nullptr;
};
#endif
*/
