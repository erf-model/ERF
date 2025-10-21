#ifndef ERF_SATADJ_H
#define ERF_SATADJ_H

/*
 * Implementation saturation adjustment microphysics model
 * The model transports qv and qc and does Newton iterations
 * to complete condensation/evaporation.
 */

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_GpuContainers.H>

#include "ERF_EOS.H"
#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"
#include "ERF_TileNoZ.H"

namespace MicVar_SatAdj {
   enum {
      // independent variables
      rho=0, // density
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      pres,  // pressure
      // non-precipitating vars
      qv,    // cloud vapor
      qc,   // cloud water
      NumVars
  };
}

class SatAdj : public NullMoist {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // constructor
    SatAdj () {}

    // destructor
    virtual ~SatAdj () = default;

    // cloud physics
    void AdvanceSatAdj (const SolverChoice& /*solverChoice*/);

    // Set up for first time
    void
    Define (SolverChoice& sc) override
    {
        m_fac_cond = lcond / sc.c_p;
        m_rdOcp    = sc.rdOcp;
        m_do_cond  = (!sc.use_shoc);
    }

    // init
    void
    Init (const amrex::MultiFab& cons_in,
          const amrex::BoxArray& /*grids*/,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance,
          std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
          std::unique_ptr<amrex::MultiFab>& /*detJ_cc*/) override;

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

        this->AdvanceSatAdj(solverChoice);
    }

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return nullptr;
    }

    int
    Qmoist_Size () override { return SatAdj::m_qmoist_size; }

    int
    Qstate_Moist_Size () override { return SatAdj::n_qstate_moist_size; }

    void
    Qmoist_Restart_Vars (const SolverChoice& /*a_sc*/,
                         std::vector<int>& a_idx,
                         std::vector<std::string>& a_names) const override
    {
        a_idx.clear();
        a_names.clear();
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    NewtonIterSat (int& i, int& j, int& k,
                   const amrex::Real& fac_cond,
                   const amrex::Array4<amrex::Real>& tabs_array,
                   const amrex::Array4<amrex::Real>& pres_array,
                   const amrex::Array4<amrex::Real>& qv_array,
                   const amrex::Array4<amrex::Real>& qc_array)
    {
        // Solution tolerance
        amrex::Real tol = 1.0e-8;

        // Saturation moisture fractions
        amrex::Real qsat, dqsat;

        // Newton iteration vars
        int  niter;
        amrex::Real fff, dfff;
        amrex::Real dtabs, delta_qv;

        // Initial guess for temperature & pressure
        amrex::Real tabs = tabs_array(i,j,k);
        amrex::Real pres = pres_array(i,j,k);

        niter = 0;
        dtabs = 1;

        //==================================================
        // Newton iteration to qv=qsat (cloud phase only)
        //==================================================
        do {
            // Saturation moisture fractions
            erf_qsatw(tabs, pres, qsat);
            erf_dtqsatw(tabs, pres, dqsat);

            // Function for root finding:
            // 0 = -T_new + T_old + L_eff/C_p * (qv - qsat)
            fff   = -tabs + tabs_array(i,j,k) + fac_cond*(qv_array(i,j,k) - qsat);

            // Derivative of function (T_new iterated on)
            dfff  = -1.0 - fac_cond*dqsat;

            // Update the temperature
            dtabs = -fff/dfff;
            tabs += dtabs;

            // Update iteration
            niter = niter+1;
        } while(std::abs(dtabs) > tol && niter < 20);

        // Update qsat from last iteration (dq = dq/dt * dt)
        qsat += dqsat*dtabs;

        // Changes in each component
        delta_qv = qv_array(i,j,k) - qsat;

        // Partition the change in non-precipitating q
        qv_array(i,j,k)  = qsat;
        qc_array(i,j,k) += delta_qv;

        // Return to temperature
        return tabs;
    }

private:
    // Number of qmoist variables (no precipitating comps to accumulate)
    int m_qmoist_size = 0;

    // Number of qstate variables
    int n_qstate_moist_size = 2;

    // geometry
    amrex::Geometry m_geom;

    // timestep
    amrex::Real dt;

    // constants
    amrex::Real m_fac_cond;
    amrex::Real m_rdOcp;
    bool m_do_cond;

    // independent variables
    amrex::Array<FabPtr, MicVar_SatAdj::NumVars> mic_fab_vars;
};
#endif
