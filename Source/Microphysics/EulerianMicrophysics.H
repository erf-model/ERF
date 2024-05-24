/*! @file EulerianMicrophysics.H
 *  \brief Containes the Eulerian microphysics class
 */

#ifndef EULERIANMICROPHYSICS_H
#define EULERIANMICROPHYSICS_H

#include "NullMoist.H"
#include "SAM.H"
#include "Kessler.H"
#include "Microphysics.H"

/*! \brief Eulerian microphysics interface */
class EulerianMicrophysics : public Microphysics {

public:

    /*! \brief Null constructor */
    EulerianMicrophysics () { }

    /*! \brief default destructor */
    ~EulerianMicrophysics () = default;

    /*! \brief Constructor: create the moisture model */
    EulerianMicrophysics (const int& nlev, /*!< Number of AMR levels */
                          const MoistureType& a_model_type /*!< moisture model */)
    {
        AMREX_ASSERT( Microphysics::modelType(a_model_type) == MoistureModelType::Eulerian );
        m_moist_model.resize(nlev);
        if (a_model_type == MoistureType::SAM or
		    a_model_type == MoistureType::SAM_NoPrecip_NoIce) {
            SetModel<SAM>();
            amrex::Print() << "SAM moisture model!\n";
        } else if (a_model_type == MoistureType::Kessler or
                   a_model_type == MoistureType::Kessler_NoRain) {
            SetModel<Kessler>();
            amrex::Print() << "Kessler moisture model!\n";
        } else if (a_model_type == MoistureType::None) {
            SetModel<NullMoist>();
            amrex::Print() << "No moisture model!\n";
        } else {
            amrex::Abort("EulerianMicrophysics: Dont know this moisture_type!") ;
        }
    }

    /*! \brief Define the moisture model */
    void Define (const int& lev, /*!< AMR level */
                 SolverChoice& sc /*!< Solver choice object */) override
    {
        m_moist_model[lev]->Define(sc);
    }

    /*! \brief Initialize the moisture model */
    void Init (const int& lev, /*!< AMR level */
               const amrex::MultiFab& cons_in, /*!< Conserved state variables */
               const amrex::BoxArray& grids,   /*!< Grids */
               const amrex::Geometry& geom,    /*!< Geometry */
               const amrex::Real& dt_advance,   /*!< Time step */
               std::unique_ptr<amrex::MultiFab>& z_phys_nd, /*< Nodal z heights */
               std::unique_ptr<amrex::MultiFab>& detJ_cc    /*< CC Jacobian determinants */) override
    {
        m_moist_model[lev]->Init(cons_in, grids, geom, dt_advance,
                                 z_phys_nd, detJ_cc);
    }

    /*! \brief Advance the moisture model for one time step */
    void Advance (const int& lev, /*!< AMR level */
                  const amrex::Real& dt_advance, /*!< Time step */
                  const int&, /*!< iteration number */
                  const amrex::Real&, /*!< current time */
                  const SolverChoice &solverChoice, /*!< Solver choice object */
                  amrex::Vector<amrex::Vector<amrex::MultiFab>>&, /*!< Dycore state variables */
                  const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& /*!< terrain */) override
    {
        m_moist_model[lev]->Advance(dt_advance, solverChoice);
    }

    /*! \brief update microphysics variables from ERF state variables */
    void Update_Micro_Vars_Lev (const int& lev, /*! AMR level */
                                amrex::MultiFab& cons_in /*!< Conserved state variables */) override
    {
        m_moist_model[lev]->Update_Micro_Vars(cons_in);
    }

    /*! \brief update ERF state variables from microphysics variables */
    void Update_State_Vars_Lev (const int& lev, /*!< AMR level */
                                amrex::MultiFab& cons_in /*!< Conserved state variables */) override
    {
        m_moist_model[lev]->Update_State_Vars(cons_in);
    }

    /*! \brief get pointer to a moisture variable */
    amrex::MultiFab* Get_Qmoist_Ptr (const int& lev, /*!< AMR level */
                                     const int& varIdx /*!< moisture variable index */) override
    {
        return m_moist_model[lev]->Qmoist_Ptr(varIdx);
    }

    /*! \brief get the number of moisture model variables */
    int Get_Qmoist_Size (const int& /* lev */) override
    {
        return m_moist_model[0]->Qmoist_Size();
    }

    /*! \brief get the number of moisture-model-related conserved state variables */
    int Get_Qstate_Size () override
    {
        return m_moist_model[0]->Qstate_Size();
    }

protected:

    /*! \brief Create and set the specified moisture model */
    template<class NewMoistModel>
    void SetModel ()
    {
        for (int lev(0); lev<m_moist_model.size(); ++lev) {
            m_moist_model[lev] = std::make_unique<NewMoistModel>();
        }
    }

private:
    amrex::Vector<std::unique_ptr<NullMoist>> m_moist_model; /*!< moisture model */
};
#endif
