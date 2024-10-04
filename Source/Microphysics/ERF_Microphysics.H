/*! @file ERF_Microphysics.H
 *  \brief Contains the base class for microphysics
 */

#ifndef ERF_MICROPHYSICS_H
#define ERF_MICROPHYSICS_H

#include <vector>
#include <string>
#include "ERF_DataStruct.H"

/*! \brief Base class for microphysics interface */
class Microphysics {

public:

    /*! \brief Null constructor */
    Microphysics () { }

    /*! \brief default destructor */
    virtual ~Microphysics () = default;

    /*! \brief define the microphysics object */
    virtual void Define (const int&, SolverChoice&) = 0;

    /*! \brief initialize the microphysics object */
    virtual void Init (const int&,
                       const amrex::MultiFab&,
                       const amrex::BoxArray&,
                       const amrex::Geometry&,
                       const amrex::Real&,
                       std::unique_ptr<amrex::MultiFab>&,
                       std::unique_ptr<amrex::MultiFab>&) = 0;

    /*! \brief advance microphysics for one time step */
    virtual void Advance (const int&,
                          const amrex::Real&,
                          const int&,
                          const amrex::Real&,
                          const SolverChoice&,
                          amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ) = 0;

    /*! \brief update microphysics variables from ERF state variables */
    virtual void Update_Micro_Vars_Lev (const int&, amrex::MultiFab&) = 0;

    /*! \brief update ERF state variables from microphysics variables */
    virtual void Update_State_Vars_Lev (const int&, amrex::MultiFab&) = 0;

    /*! \brief get pointer to a moisture variable */
    virtual amrex::MultiFab* Get_Qmoist_Ptr (const int&, const int&) = 0;

    /*! \brief get the number of moisture model variables */
    virtual int Get_Qmoist_Size (const int&) = 0;

    /*! \brief get the number of moisture-model-related conserved state variables */
    virtual int Get_Qstate_Size () = 0;

    /*! \brief get the indices and names of moisture model variables for restart
               at a given level */
    virtual void Get_Qmoist_Restart_Vars ( int, std::vector<int>&, std::vector<std::string>& ) const = 0;

    /*! \brief query if a specified moisture model is Eulerian or Lagrangian */
    static MoistureModelType modelType (const MoistureType a_moisture_type)
    {
        if (    (a_moisture_type == MoistureType::SAM)
             || (a_moisture_type == MoistureType::SAM_NoIce)
             || (a_moisture_type == MoistureType::SAM_NoPrecip_NoIce)
             || (a_moisture_type == MoistureType::Kessler)
             || (a_moisture_type == MoistureType::Kessler_NoRain)
             || (a_moisture_type == MoistureType::None) ) {
            return MoistureModelType::Eulerian;
        } else {
            amrex::Abort("Dont know this moisture_type!") ;
            return MoistureModelType::Undefined;
        }
    }

private:

};
#endif
