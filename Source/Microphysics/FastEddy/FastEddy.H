/*
 */
#ifndef FastEddy_H
#define FastEddy_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Constants.H"
#include "Microphysics_Utils.H"
#include "IndexDefines.H"
#include "DataStruct.H"

namespace MicVar_FE {
   enum {
      // independent variables
      qv = 0,
      qc,
      qt,
      rho,   // density
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      pres,  // pressure
      NumVars
  };
}

//
// use MultiFab for 3D data, but table for 1D data
//
class FastEddy : public NullMoist {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // constructor
    FastEddy () {}

    // destructor
    virtual ~FastEddy () = default;

    // cloud physics
    void AdvanceFE ();

    // diagnose
    void Diagnose () override;

    // Set up for first time
    void
    Define (SolverChoice& sc) override
    {
        docloud = sc.do_cloud;
        doprecip = sc.do_precip;
        m_fac_cond = lcond / sc.c_p;
        m_fac_fus = lfus / sc.c_p;
        m_fac_sub = lsub / sc.c_p;
        m_gOcp = CONST_GRAV / sc.c_p;
        m_axis = sc.ave_plane;
    }

    // init
    void
    Init (const amrex::MultiFab& cons_in,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance) override;

    // Copy state into micro vars
    void
    Copy_State_to_Micro (const amrex::MultiFab& cons_in) override;

    // Copy micro into state vars
    void
    Copy_Micro_to_State (amrex::MultiFab& cons_in) override;

    void
    Update_Micro_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_State_to_Micro(cons_in);
        this->Diagnose();
    }

    void
    Update_State_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_Micro_to_State(cons_in);
    }

    // wrapper to do all the updating
    void
    Advance (const amrex::Real& dt_advance) override
    {
        dt = dt_advance;

        this->AdvanceFE();
        this->Diagnose();
    }

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    int
    Qmoist_Size () override { return FastEddy::m_qmoist_size; }

    int
    Qstate_Size () { return FastEddy::m_qstate_size; }

private:
    // Number of qmoist variables (qt, qv, qc)
    int m_qmoist_size = 3;

    // Number of qstate variables
    int m_qstate_size = 2;

    // MicVar map (Qmoist indices -> MicVar enum)
    amrex::Vector<int> MicVarMap;

    // geometry
    amrex::Geometry m_geom;

    // valid boxes on which to evolve the solution
    amrex::BoxArray m_gtoe;

    // timestep
    amrex::Real dt;

    // number of vertical levels
    int nlev, zlo, zhi;

    // plane average axis
    int m_axis;

    // model options
    bool docloud, doprecip;

    // constants
    amrex::Real m_fac_cond;
    amrex::Real m_fac_fus;
    amrex::Real m_fac_sub;
    amrex::Real m_gOcp;

    // independent variables
    amrex::Array<FabPtr, MicVar_FE::NumVars> mic_fab_vars;
};
#endif
