/*
 * Implementation of a 2-moment microphysics model
 * This model is based on the Morrison model as implemented in WRF (in Fortran)
 */
#ifndef ERF_Morrison_H
#define ERF_Morrison_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"
#include <ERF_Morrison_Fortran_Interface.H>

namespace MicVar_Morr {
   enum {
      // independent variables
      rho=0, // density
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      pres,  // pressure
      // non-precipitating vars
      qt,    // total cloud
      qn,    // cloud condensate (liquid + ice)
      qv,    // cloud vapor
      qcl,   // cloud water
      qci,   // cloud ice
      // precipitating vars
      qp,    // total precip
      qpr,   // precip rain
      qps,   // precip ice
      qpg,   // graupel
      // The "numbers" and the "accums" are written/read at checkpoint/restart
      // number concentrations
      nc,    // cloud droplet number
      nr,    // rain number
      ni,    // cloud ice number
      ns,    // snow number
      ng,    // graupel number
      // derived vars
      rain_accum,
      snow_accum,
      graup_accum,
      omega,
      NumVars
  };
}

class Morrison : public NullMoist {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // constructor
    Morrison () {}

    // destructor
    virtual ~Morrison () = default;

    // Set up for first time
    void
    Define (SolverChoice& sc) override
    {
        m_axis     = sc.ave_plane;
        m_rdOcp    = sc.rdOcp;
        m_do_cond  = (!sc.use_shoc);
    }

    // init
    void
    Init (const amrex::MultiFab& cons_in,
          const amrex::BoxArray& grids,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance,
          std::unique_ptr<amrex::MultiFab>& z_phys_nd,
          std::unique_ptr<amrex::MultiFab>& detJ_cc) override;

    // Import minimum dz at this level
    void
    Set_dzmin (const amrex::Real dz_min) override
    {
        m_dzmin = dz_min;
    }

    // Copy state into micro vars
    void
    Copy_State_to_Micro (const amrex::MultiFab& cons_in) override;

    // Copy state into micro vars
    void
    Copy_Micro_to_State (amrex::MultiFab& cons_in) override;

    void
    Update_Micro_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_State_to_Micro(cons_in);
        this->Compute_Coefficients();
    }

    void
    Update_State_Vars (amrex::MultiFab& cons_in,
                       const amrex::MultiFab& /*z_phys_nd*/) override
    {
        this->Copy_Micro_to_State(cons_in);
    }

    // wrapper to do all the updating
    void
    Advance (const amrex::Real& dt_advance,
             const SolverChoice& sc) override;

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    void
    Compute_Coefficients ();

    int
    Qmoist_Size () override { return Morrison::m_qmoist_size; }

    int
    Qstate_Moist_Size () override { return Morrison::n_qstate_moist_size; }

    void
    Set_RealWidth (const int real_width) override { m_real_width = real_width; }

    void
    Qmoist_Restart_Vars (const SolverChoice& /*a_sc*/,
                         std::vector<int>& a_idx,
                         std::vector<std::string>& a_names) const override
    {
        a_idx.clear();
        a_names.clear();

        // The ordering here needs to match that in
        // MicVarMap = {MicVar_Morr::nc, MicVar_Morr::nr, MicVar_Morr::ni, MicVar_Morr::ns, MicVar_Morr::ng,
        //              MicVar_Morr::rain_accum, MicVar_Morr::snow_accum, MicVar_Morr::graup_accum};
        //
        a_idx.push_back(0); a_names.push_back("Nc");
        a_idx.push_back(1); a_names.push_back("Nr");
        a_idx.push_back(2); a_names.push_back("Ni");
        a_idx.push_back(3); a_names.push_back("Ns");
        a_idx.push_back(4); a_names.push_back("Ng");
        a_idx.push_back(5); a_names.push_back("RainAccum");
        a_idx.push_back(6); a_names.push_back("SnowAccum");
        a_idx.push_back(7); a_names.push_back("GraupAccum");
    }

    void
    GetPlotVarNames (amrex::Vector<std::string>& a_vec) const override;

    void
    GetPlotVar (const std::string& a_name,
                amrex::MultiFab& a_mf) const override;

    void
    GetPlotVar (const std::string& a_name,
                amrex::MultiFab& a_mf,
                const int /*a_lev*/) const override
    {
        GetPlotVar(a_name, a_mf);
    }

private:
    // Number of qmoist variables (Ncloud, Nrain, Nice, Nsnow, Ngraup, rain_accum, snow_accum, graup_accum)
    int m_qmoist_size = 8;

    // Number of qstate variables
    int n_qstate_moist_size = 6;

    // CFL MAX for vertical advection
    static constexpr amrex::Real CFL_MAX = myhalf;

    // MicVar map (Qmoist indices -> MicVar enum)
    amrex::Vector<int> MicVarMap;

    // geometry
    amrex::Geometry m_geom;

    // Nudging + set width
    int m_real_width{0};

    // number of vertical levels
    int nlev, zlo, zhi;

    // plane average axis
    int m_axis;

    // constants
    amrex::Real m_rdOcp;
    bool m_do_cond;

    // Minimum dz at this level
    amrex::Real m_dzmin;

    // Pointer to terrain data
    amrex::MultiFab* m_z_phys_nd;
    amrex::MultiFab* m_detJ_cc;

    // independent variables
    amrex::Array<FabPtr, MicVar_Morr::NumVars> mic_fab_vars;

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
};
#endif
