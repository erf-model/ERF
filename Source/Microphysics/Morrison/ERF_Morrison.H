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
#include "ERF_SurfacePrecipitation.H"
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

struct MorrisonThermoState {
    amrex::Real pressure_pa;
    amrex::Real temperature;
};

// Morrison's C++ and Fortran paths use Pa internally. Keep this small adapter
// beside the scheme so the production copy-in wiring and its unit coverage use
// exactly the same stored-pi0 diagnosis.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MorrisonThermoState morrison_copy_in_thermo (
    const amrex::Real rho,
    const amrex::Real rho_theta,
    const amrex::Real qv,
    const bool use_anelastic_reference_pressure,
    const amrex::Real p0,
    const amrex::Real pi0) noexcept
{
    const MicrophysicsThermoState thermo = diagnose_microphysics_thermo_state(
        rho, rho_theta, qv, use_anelastic_reference_pressure, p0, pi0);
    return {thermo.pressure_pa, thermo.temperature};
}

// Array-level copy-in used by production and by the wiring test.  The number
// concentration arrays remain in the caller because they are independent of
// the thermodynamic/base-state contract.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void morrison_copy_state_to_micro_cell (
    const amrex::Array4<const amrex::Real>& states,
    const amrex::Array4<const amrex::Real>& base,
    const amrex::Array4<amrex::Real>& rho,
    const amrex::Array4<amrex::Real>& theta,
    const amrex::Array4<amrex::Real>& qv,
    const amrex::Array4<amrex::Real>& qc,
    const amrex::Array4<amrex::Real>& qi,
    const amrex::Array4<amrex::Real>& qn,
    const amrex::Array4<amrex::Real>& qt,
    const amrex::Array4<amrex::Real>& qpr,
    const amrex::Array4<amrex::Real>& qps,
    const amrex::Array4<amrex::Real>& qpg,
    const amrex::Array4<amrex::Real>& qp,
    const amrex::Array4<amrex::Real>& tabs,
    const amrex::Array4<amrex::Real>& pres,
    const bool use_anelastic_reference_pressure,
    const int i, const int j, const int k) noexcept
{
    const amrex::Real rho_value = states(i,j,k,Rho_comp);
    const amrex::Real rho_theta = states(i,j,k,RhoTheta_comp);
    const amrex::Real qv_value =
        amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ1_comp) / rho_value);
    rho(i,j,k) = rho_value;
    theta(i,j,k) = rho_theta / rho_value;
    qv(i,j,k) = qv_value;
    qc(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ2_comp) / rho_value);
    qi(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ3_comp) / rho_value);
    qn(i,j,k) = qc(i,j,k) + qi(i,j,k);
    qt(i,j,k) = qv(i,j,k) + qn(i,j,k);
    qpr(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ4_comp) / rho_value);
    qps(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ5_comp) / rho_value);
    qpg(i,j,k) = amrex::max(amrex::Real(0.0), states(i,j,k,RhoQ6_comp) / rho_value);
    qp(i,j,k) = qpr(i,j,k) + qps(i,j,k) + qpg(i,j,k);

    const amrex::Real p0 = use_anelastic_reference_pressure
        ? base(i,j,k,BaseState::p0_comp) : amrex::Real(0.0);
    const amrex::Real pi0 = use_anelastic_reference_pressure
        ? base(i,j,k,BaseState::pi0_comp) : amrex::Real(0.0);
    const MorrisonThermoState thermo = morrison_copy_in_thermo(
        rho_value, rho_theta, qv_value,
        use_anelastic_reference_pressure, p0, pi0);
    tabs(i,j,k) = thermo.temperature;
    pres(i,j,k) = thermo.pressure_pa;
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
        m_moisture_type = sc.moisture_type;
        m_rdOcp    = sc.rdOcp;
        m_do_cond  = (!sc.uses_shoc_family());
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

    void
    Copy_State_to_Micro (const amrex::MultiFab& cons_in,
                         const amrex::MultiFab* base_state);

    void
    Update_Micro_Vars (amrex::MultiFab& cons_in,
                       const amrex::MultiFab* base_state) override;

    // Copy state into micro vars
    void
    Copy_Micro_to_State (amrex::MultiFab& cons_in) override;

    void
    Update_Micro_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_State_to_Micro(cons_in);
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

    int
    Qmoist_Size () override { return Morrison::m_qmoist_size; }

    int
    Qstate_Moist_Size () override { return Morrison::n_qstate_moist_size; }

    int
    Qstate_Moist_NumConc_Size () override { return Morrison::n_qstate_moist_numconc_size; }

    void
    Qmoist_Restart_Vars (const SolverChoice& /*a_sc*/,
                         std::vector<int>& a_idx,
                         std::vector<std::string>& a_names) const override
    {
        a_idx.clear();
        a_names.clear();

        // The ordering here needs to match that in
        // MicVarMap = {MicVar_Morr::rain_accum, MicVar_Morr::snow_accum, MicVar_Morr::graup_accum};
        //
        a_idx.push_back(0); a_names.push_back("RainAccum");
        a_idx.push_back(1); a_names.push_back("SnowAccum");
        a_idx.push_back(2); a_names.push_back("GraupAccum");
    }

    SurfacePrecipAccumulationSources
    Get_Surface_Precip_Accumulation_Ptrs (const int&) const override
    {
        SurfacePrecipAccumulationSources sources;
        sources.total = {mic_fab_vars[MicVar_Morr::rain_accum].get(), rhor / amrex::Real(1000.0)};

        if (m_moisture_type == MoistureType::Morrison) {
            sources.snow    = {mic_fab_vars[MicVar_Morr::snow_accum].get(), rhor / amrex::Real(1000.0)};
            sources.graupel = {mic_fab_vars[MicVar_Morr::graup_accum].get(), rhor / amrex::Real(1000.0)};
        }

        return sources;
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
    // Number of qmoist variables (rain_accum, snow_accum, graup_accum)
    int m_qmoist_size = 3;

    // Number of qstate variables
    int n_qstate_moist_size = 11;

    // Number of qstate variables that are number concentrations
    int n_qstate_moist_numconc_size = 5;

    // MicVar map (Qmoist indices -> MicVar enum)
    amrex::Vector<int> MicVarMap;

    // geometry
    amrex::Geometry m_geom;

    // constants
    amrex::Real m_rdOcp;
    bool m_do_cond;
    MoistureType m_moisture_type = MoistureType::None;

    // Minimum dz at this level
    amrex::Real m_dzmin;

    // Pointer to terrain data
    amrex::MultiFab* m_z_phys_nd;
    amrex::MultiFab* m_detJ_cc;

    // independent variables
    amrex::Array<FabPtr, MicVar_Morr::NumVars> mic_fab_vars;
};
#endif
