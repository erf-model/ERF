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
    Update_State_Vars (amrex::MultiFab& cons_in,
                       const amrex::MultiFab& /*z_phys_nd*/) override
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

    void
    Set_RealWidth (const int real_width) override { m_real_width = real_width; }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    NewtonSolveSatTemperature (const amrex::Real fac_cond,
                               const amrex::Real tabs_old,
                               const amrex::Real pres_mbar,
                               const amrex::Real qv_old,
                               amrex::Real& qsat)
    {
        constexpr amrex::Real tol = amrex::Real(1.0e-8);
        constexpr int max_niter = 20;

        amrex::Real dqsat;

        int  niter = 0;
        amrex::Real fff, dfff;
        amrex::Real dtabs = one;

        amrex::Real tabs = tabs_old;

        //==================================================
        // Newton iteration to qv=qsat (cloud phase only)
        //==================================================
        do {
            // Saturation moisture fractions
            erf_qsatw(tabs, pres_mbar, qsat);
            erf_dtqsatw(tabs, pres_mbar, dqsat);

            // Function for root finding:
            // 0 = -T_new + T_old + L_eff/C_p * (qv - qsat)
            fff   = -tabs + tabs_old + fac_cond*(qv_old - qsat);

            // Derivative of function (T_new iterated on)
            dfff  = -one - fac_cond*dqsat;

            // Update the temperature
            dtabs = -fff/dfff;
            tabs += dtabs;

            // Update iteration
            niter = niter+1;
        } while(std::abs(dtabs) > tol && niter < max_niter);

        // Update qsat from last iteration (dq = dq/dt * dt)
        qsat += dqsat*dtabs;

        return tabs;
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static void
    AdjustSatAdjCell (const amrex::Real fac_cond,
                      const amrex::Real rdOcp,
                      amrex::Real& tabs,
                      const amrex::Real pres_mbar,
                      amrex::Real& theta,
                      amrex::Real& qv,
                      amrex::Real& qc)
    {
        amrex::Real qsat;
        erf_qsatw(tabs, pres_mbar, qsat);

        if ((qv + qc) > qsat) {
#if defined(AMREX_DEBUG)
            const amrex::Real qvprev = qv;
            const amrex::Real qcprev = qc;
#endif

            if (qc < 0) {
                qv += qc;
                qc = zero;
            }

            tabs = NewtonSolveSatTemperature(fac_cond, tabs, pres_mbar, qv, qsat);

            const amrex::Real delta_qv = qv - qsat;
            qv = qsat;
            qc += delta_qv;

#if defined(AMREX_DEBUG)
            amrex::Real qsatnew;
            erf_qsatw(tabs, pres_mbar, qsatnew);
            AMREX_ASSERT(std::abs(qv-qsatnew) < 1e-12);
            AMREX_ASSERT(std::abs(qv+qc-qvprev-qcprev) < 1e-14);
#endif

            theta = getThgivenTandP(tabs, amrex::Real(100.0)*pres_mbar, rdOcp);
        } else {
            const amrex::Real delta_qc = qc;

            qv += qc;
            qc = zero;

            tabs -= fac_cond * delta_qc;
            theta = getThgivenTandP(tabs, amrex::Real(100.0)*pres_mbar, rdOcp);

            erf_qsatw(tabs, pres_mbar, qsat);
            if (qv > qsat) {
#if defined(AMREX_DEBUG)
                const amrex::Real qvprev = qv;
                const amrex::Real qcprev = qc;
                const amrex::Real Tprev = tabs;
#endif

                tabs = NewtonSolveSatTemperature(fac_cond, tabs, pres_mbar, qv, qsat);

                const amrex::Real delta_qv = qv - qsat;
                qv = qsat;
                qc += delta_qv;

#if defined(AMREX_DEBUG)
                amrex::Real qsatnew;
                erf_qsatw(tabs, pres_mbar, qsatnew);
                AMREX_ASSERT(qv < qvprev);
                AMREX_ASSERT(qc > qcprev);
                AMREX_ASSERT(tabs > Tprev);
                AMREX_ASSERT(std::abs(qv-qsatnew) < 1e-14);
                AMREX_ASSERT(std::abs(qv+qc-qvprev-qcprev) < 1e-14);
#endif

                theta = getThgivenTandP(tabs, amrex::Real(100.0)*pres_mbar, rdOcp);
            }
        }
    }

private:
    // Number of qmoist variables (no precipitating comps to accumulate)
    int m_qmoist_size = 0;

    // Number of qstate variables
    int n_qstate_moist_size = 2;

    // geometry
    amrex::Geometry m_geom;

    // Nudging + set width
    int m_real_width{0};

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
