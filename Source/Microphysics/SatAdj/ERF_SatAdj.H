#ifndef ERF_SATADJ_H
#define ERF_SATADJ_H

/*
 * SatAdj is a cell-local, fixed-pressure saturation-adjustment source step.
 * It is not a flux-form finite-volume update.
 *
 * The scheme reads rho, rhoTheta, rhoQv, and rhoQc from the conserved state,
 * adjusts theta/qv/qc locally, and writes rhoTheta/rhoQv/rhoQc back.
 * Density is copied for diagnostics and conservative reconstruction, but rho is
 * not a SatAdj prognostic variable and is not modified by this module.
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
    rho=0, // density copied from conserved state; not modified by SatAdj
    theta, // dry potential temperature
    tabs,  // absolute temperature [K]
    pres,  // pressure [mbar/hPa] for saturation helpers
      // non-precipitating vars
    qv,    // water-vapor mixing ratio
    qc,    // cloud-water mixing ratio
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

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    // Newton solve for the saturated final temperature at fixed pressure.
    // The invariant is T + (L/cp) * qv for the local cell. The solve finds
    // T_new such that qv_new = qsat(T_new, p). The pressure argument is held
    // fixed in mbar/hPa throughout the solve.
    //
    // Solves the fixed-pressure saturation equation
    //   0 = -T_new + T_old + (L/cp) * (qv_old - qsat(T_new, p))
    // for T_new. On return, qsat is evaluated exactly at the returned
    // temperature. The caller uses that qsat to set qv and qc while
    // conserving qv + qc. Density does not enter this solve.
    // If the Newton iteration reaches max_niter before the tolerance, the
    // routine returns the last iterate and recomputes qsat at that temperature.
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

        // Return qsat evaluated at the final temperature, not a tangent extrapolation.
        erf_qsatw(tabs, pres_mbar, qsat);

        return tabs;
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static void
    // Cell-local fixed-pressure saturation adjustment.
    // Inputs/outputs are specific quantities for one cell. The density is not
    // used or modified. The update conserves qv + qc and, at the pressure
    // passed in, conserves T + (L/cp) * qv. The final state should satisfy
    // qc >= 0 and either qc == 0 or qv == qsat(T,p).
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
            // Total water exceeds saturation at the initial temperature, so the final state
            // is cloudy and saturated. Solve directly for the final saturated temperature.
            // This is algebraically equivalent to evaporating existing cloud water first
            // and then recondensing to the fixed-pressure equilibrium state.
#if defined(AMREX_DEBUG)
            const amrex::Real qvprev = qv;
            const amrex::Real qcprev = qc;
#endif

            if (qc < 0) {
                // Repair negative cloud water by transferring the deficit back to vapor
                // before solving the saturated equilibrium and clamping qc to zero.
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
            // Total water does not exceed initial saturation. Evaporate all cloud water,
            // cooling the cell by latent heat absorption. The cooled cell can become
            // supersaturated, so recondense if needed.
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
