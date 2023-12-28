/*
 * Implementation 1-moment microphysics model
 * NOTE: this model is based on the Kessler code, and the Klemp's paper
 * 1): Joseph, Klemp, the simulation of three-dimensional convective storm dynamics,
 * Journal of the atmospheric sciences, vol35, p1070
 * 2): Marat Khairoutdinov and David Randall, cloud resolving modeling of the ARM summer 1997 IOP:
 * model formulation, results, unvertainties, and sensitivities, Journal of the atmospheric sciences, vol60, p607
 */
#ifndef Kessler_H
#define Kessler_H

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

namespace MicVar_Kess {
   enum {
      // independent variables
      qt = 0,
      qp,
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      rho,   // density
      pres,  // pressure
      // derived variables
      qr,   // rain water
      qv,   // water vapor
      qn,   // cloud condensate (liquid+ice), initial to zero
      qci,  // cloud ice
      qcl,  // cloud water
      qpl,  // precip rain
      qpi,  // precip ice
      qg,   // graupel
      // temporary variable
      omega,
      NumVars
  };
}

//
// use MultiFab for 3D data, but table for 1D data
//
class Kessler : public NullMoist {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // constructor
    Kessler () {}

    // destructor
    virtual ~Kessler () = default;

    // cloud physics
    void AdvanceKessler ();

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

    // Copy state into micro vars
    void
    Copy_Micro_to_State (amrex::MultiFab& cons_in) override;

    // update ERF variables
    void
    Update (amrex::MultiFab& cons_in,
            amrex::MultiFab& qmoist) override;

    // update micro vars
    void
    Update_Micro_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_State_to_Micro(cons_in);
        this->Diagnose();
        this->Compute_Tabs_and_Pres();
    }

    // update state vars
    void
    Update_State_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_Micro_to_State(cons_in);
    }

    // wrapper to advance micro vars
    void
    Advance (const amrex::Real& dt_advance) override
    {
        dt = dt_advance;

        this->AdvanceKessler();
        this->Diagnose();
        this->Compute_Tabs_and_Pres();
    }

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    void
    Compute_Tabs_and_Pres ();

    int
    Qmoist_Size () override { return Kessler::m_qmoist_size; }

private:
    // Number of qmoist variables
    int m_qmoist_size = 6;

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
    amrex::Array<FabPtr, MicVar_Kess::NumVars> mic_fab_vars;

#if 0
    // microphysics parameters/coefficients
    amrex::TableData<amrex::Real, 1> accrrc;
    amrex::TableData<amrex::Real, 1> accrsi;
    amrex::TableData<amrex::Real, 1> accrsc;
    amrex::TableData<amrex::Real, 1> coefice;
    amrex::TableData<amrex::Real, 1> evaps1;
    amrex::TableData<amrex::Real, 1> evaps2;
    amrex::TableData<amrex::Real, 1> accrgi;
    amrex::TableData<amrex::Real, 1> accrgc;
    amrex::TableData<amrex::Real, 1> evapg1;
    amrex::TableData<amrex::Real, 1> evapg2;
    amrex::TableData<amrex::Real, 1> evapr1;
    amrex::TableData<amrex::Real, 1> evapr2;

    // vertical plane average data
    amrex::TableData<amrex::Real, 1> rho1d;
    amrex::TableData<amrex::Real, 1> pres1d;
    amrex::TableData<amrex::Real, 1> tabs1d;
    amrex::TableData<amrex::Real, 1> qt1d;
    amrex::TableData<amrex::Real, 1> qv1d;
    amrex::TableData<amrex::Real, 1> qn1d;

    // independent variables
    amrex::Array<FabPtr, MicVar_Kess::NumVars> mic_fab_vars;

    amrex::TableData<amrex::Real, 1> gamaz;
    amrex::TableData<amrex::Real, 1> zmid; // mid value of vertical coordinate in physical domain

    // data (output)
    amrex::TableData<amrex::Real, 1> qifall;
    amrex::TableData<amrex::Real, 1> tlatqi;
#endif

};
#endif
