/*
 * Implementation 1-moment microphysics model
 * NOTE: this model is based on the Kessler code, and the Klemp's paper
 * 1): Joseph, Klemp, the simulation of three-dimensional convective storm dynamics,
 * Journal of the atmospheric sciences, vol35, p1070
 * 2): Marat Khairoutdinov and David Randall, cloud resolving modeling of the ARM summer 1997 IOP:
 * model formulation, results, unvertainties, and sensitivities, Journal of the atmospheric sciences, vol60, p607
 */
#ifndef ERF_Kessler_H
#define ERF_Kessler_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Constants.H"
#include "ERF_Microphysics_Utils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"

namespace MicVar_Kess {
   enum {
      // independent variables
      rho=0, // density
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      pres,  // pressure
      // non-precipitating vars
      qt,    // total cloud
      qv,    // cloud vapor
      qcl,   // cloud water
      // precipitating vars
      qp,    // total precip
      // derived vars
      rain_accum,
      NumVars
  };
}

class Kessler : public NullMoist {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // constructor
    Kessler () {}

    // destructor
    virtual ~Kessler () = default;

    // cloud physics
    void AdvanceKessler (const SolverChoice &solverChoice);

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
          const amrex::Real& dt_advance,
          std::unique_ptr<amrex::MultiFab>& z_phys_nd,
          std::unique_ptr<amrex::MultiFab>& detJ_cc) override;

    // Copy state into micro vars
    void
    Copy_State_to_Micro (const amrex::MultiFab& cons_in) override;

    // Copy state into micro vars
    void
    Copy_Micro_to_State (amrex::MultiFab& cons_in) override;

    // update micro vars
    void
    Update_Micro_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_State_to_Micro(cons_in);
    }

    // update state vars
    void
    Update_State_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_Micro_to_State(cons_in);
    }

    // wrapper to advance micro vars
    void
    Advance (const amrex::Real& dt_advance,
             const SolverChoice& solverChoice) override
    {
        dt = dt_advance;

        this->AdvanceKessler(solverChoice);
    }

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    int
    Qmoist_Size () override { return Kessler::m_qmoist_size; }

    int
    Qstate_Size () override { return Kessler::m_qstate_size; }

    void
    Qmoist_Restart_Vars ( std::vector<int>& a_idx,
                          std::vector<std::string>& a_names) const override
    {
        a_idx.clear();
        a_names.clear();
        a_idx.push_back(4); a_names.push_back("RainAccum");
    }

private:
    // Number of qmoist variables (qt, qv, qcl, qp)
    int m_qmoist_size = 5;

    // Number of qstate variables
    int m_qstate_size = 3;

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

    // Pointer to terrain data
    amrex::MultiFab* m_z_phys_nd;
    amrex::MultiFab* m_detJ_cc;

    // independent variables
    amrex::Array<FabPtr, MicVar_Kess::NumVars> mic_fab_vars;
};
#endif
