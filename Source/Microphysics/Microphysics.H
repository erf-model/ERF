#ifndef MICROPHYSICS_H
#define MICROPHYSICS_H

#include <DataStruct.H>

class Microphysics {

public:

    Microphysics () { }

    virtual ~Microphysics () = default;

    virtual void Define (const int&, SolverChoice&) = 0;

    virtual void Init (const int&,
                       const amrex::MultiFab&,
                       const amrex::BoxArray&,
                       const amrex::Geometry&,
                       const amrex::Real&) = 0;

    virtual void Advance (  const int&,
                            const amrex::Real&,
                            amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                            const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ) = 0;

    virtual void Diagnose (const int&) = 0;

    virtual void Update_Micro_Vars_Lev (const int&, amrex::MultiFab&) = 0;

    virtual void Update_State_Vars_Lev (const int&, amrex::MultiFab&) = 0;

    virtual amrex::MultiFab* Get_Qmoist_Ptr (const int&, const int&) = 0;

    virtual int Get_Qmoist_Size (const int&) = 0;

    virtual int Get_Qstate_Size () = 0;

    static MoistureModelType modelType(const MoistureType a_moisture_type)
    {
        if (    (a_moisture_type == MoistureType::SAM)
             || (a_moisture_type == MoistureType::Kessler)
             || (a_moisture_type == MoistureType::FastEddy)
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
