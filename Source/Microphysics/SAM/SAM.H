/*
 * Implementation 1-moment microphysics model
 * NOTE: this model is based on the SAM code, and the Klemp's paper
 * 1): Joseph, Klemp, the simulation of three-dimensional convective storm dynamics,
 * Journal of the atmospheric sciences, vol35, p1070
 * 2): Marat Khairoutdinov and David Randall, cloud resolving modeling of the ARM summer 1997 IOP:
 * model formulation, results, unvertainties, and sensitivities, Journal of the atmospheric sciences, vol60, p607
 */
#ifndef SAM_H
#define SAM_H

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
#include "NullMoist.H"

namespace MicVar {
   enum {
      // independent variables
      rho=0, // density
      theta, // liquid/ice water potential temperature
      tabs,  // temperature
      pres,  // pressure
      // non-precipitating vars
      qt,    // total cloud
      qn,    // cloud condesnate (liquid + ice)
      qv,    // cloud vapor
      qcl,   // cloud water
      qci,   // cloud ice
      // precipitating vars
      qp,    // total precip
      qpr,   // precip rain
      qps,   // precip ice
      qpg,   // graupel
      // derived vars
      rain_accum,
      snow_accum,
      graup_accum,
      omega,
      NumVars
  };
}

class SAM : public NullMoist {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // constructor
    SAM () {}

    // destructor
    virtual ~SAM () = default;

    // cloud physics
    void Cloud (const SolverChoice& sc);

    // ice physics
    void IceFall ();

    // precip
    void Precip ();

    // precip fall
    void PrecipFall ();

    // Set up for first time
    void
    Define (SolverChoice& sc) override
    {
        docloud    = sc.do_cloud;
        doprecip   = sc.do_precip;
        m_fac_cond = lcond / sc.c_p;
        m_fac_fus  = lfus / sc.c_p;
        m_fac_sub  = lsub / sc.c_p;
        m_gOcp     = CONST_GRAV / sc.c_p;
        m_axis     = sc.ave_plane;
        m_rdOcp    = sc.rdOcp;
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

    void
    Update_Micro_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_State_to_Micro(cons_in);
        this->Compute_Coefficients();
    }

    void
    Update_State_Vars (amrex::MultiFab& cons_in) override
    {
        this->Copy_Micro_to_State(cons_in);
    }

    // wrapper to do all the updating
    void
    Advance (const amrex::Real& dt_advance,
             const SolverChoice& sc) override
    {
        dt = dt_advance;

        this->Cloud(sc);
        if(sc.moisture_type == MoistureType::SAM) {
            this->IceFall();
            this->Precip();
            this->PrecipFall();
        }
    }

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    void
    Compute_Coefficients ();

    int
    Qmoist_Size () override { return SAM::m_qmoist_size; }

    int
    Qstate_Size () override { return SAM::m_qstate_size; }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    NewtonIterSat (int& i, int& j, int& k,
                   const int& SAM_moisture_type,
                   const amrex::Real& fac_cond,
                   const amrex::Real& fac_fus,
                   const amrex::Real& fac_sub,
                   const amrex::Real& an,
                   const amrex::Real& bn,
                   const amrex::Array4<amrex::Real>& tabs_array,
                   const amrex::Array4<amrex::Real>& pres_array,
                   const amrex::Array4<amrex::Real>& qv_array,
                   const amrex::Array4<amrex::Real>& qc_array,
                   const amrex::Array4<amrex::Real>& qi_array,
                   const amrex::Array4<amrex::Real>& qn_array,
                   const amrex::Array4<amrex::Real>& qt_array)
    {
        // Solution tolerance
        amrex::Real tol = 1.0e-4;

        // Saturation moisture fractions
        amrex::Real omn, domn;
        amrex::Real qsat, dqsat;
        amrex::Real qsatw, dqsatw;
        amrex::Real qsati, dqsati;

        // Newton iteration vars
        int  niter;
        amrex::Real fff, dfff, dtabs;
        amrex::Real lstar, dlstar;
        amrex::Real lstarw, lstari;
        amrex::Real delta_qv, delta_qc, delta_qi;

        // Initial guess for temperature & pressure
        amrex::Real tabs = tabs_array(i,j,k);
        amrex::Real pres = pres_array(i,j,k);

        niter = 0;
        dtabs = 1;
        //==================================================
        // Newton iteration to qv=qsat (cloud phase only)
        //==================================================
        do {
            // Latent heats and their derivatives wrt to T
            lstarw  = fac_cond;
            lstari  = fac_fus;
            domn    = 0.0;

            // Saturation moisture fractions
            erf_qsatw(tabs, pres, qsatw);
            erf_qsati(tabs, pres, qsati);
            erf_dtqsatw(tabs, pres, dqsatw);
            erf_dtqsati(tabs, pres, dqsati);

            if (SAM_moisture_type == 1) {
                // Cloud ice not permitted (condensation & fusion)
                if(tabs >= tbgmax) {
                    omn   = 1.0;
                }
                // Cloud water not permitted (sublimation & fusion)
                else if(tabs <= tbgmin) {
                    omn    = 0.0;
                    lstarw = fac_sub;
                }
                // Mixed cloud phase (condensation & fusion)
                else {
                    omn   = an*tabs-bn;
                    domn  = an;
                }
            } else if (SAM_moisture_type == 2) {
                omn  = 1.0;
                domn = 0.0;
            }

            // Linear combination of each component
            qsat   =  omn * qsatw  + (1.0-omn) * qsati;
            dqsat  =  omn * dqsatw + (1.0-omn) * dqsati
                   + domn *  qsatw -     domn  *  qsati;
            lstar  =  omn * lstarw + (1.0-omn) * lstari;
            dlstar = domn * lstarw -     domn  * lstari;

            // Function for root finding:
            // 0 = -T_new + T_old + L_eff/C_p * (qv - qsat)
            fff   = -tabs + tabs_array(i,j,k) +  lstar*(qv_array(i,j,k) - qsat);

            // Derivative of function (T_new iterated on)
            dfff  = -1.0 + dlstar*(qv_array(i,j,k) - qsat) - lstar*dqsat;

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
        delta_qc = std::max(-qc_array(i,j,k), delta_qv * omn);
        delta_qi = std::max(-qi_array(i,j,k), delta_qv * (1.0-omn));

        // Partition the change in non-precipitating q
        qv_array(i,j,k)  = qsat;
        qc_array(i,j,k) += delta_qc;
        qi_array(i,j,k) += delta_qi;
        qn_array(i,j,k)  = qc_array(i,j,k) + qi_array(i,j,k);
        qt_array(i,j,k)  = qv_array(i,j,k) + qn_array(i,j,k);

        // Return to temperature
        return tabs;
    }

private:
    // Number of qmoist variables (qt, qv, qcl, qci, qp, qpr, qps, qpg)
    int m_qmoist_size = 11;

    // Number of qmoist variables
    int m_qstate_size = 6;

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
    amrex::Real m_rdOcp;

    // Pointer to terrain data
    amrex::MultiFab* m_z_phys_nd;
    amrex::MultiFab* m_detJ_cc;

    // independent variables
    amrex::Array<FabPtr, MicVar::NumVars> mic_fab_vars;

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

    amrex::TableData<amrex::Real, 1> gamaz;
    amrex::TableData<amrex::Real, 1> zmid; // mid value of vertical coordinate in physical domain
};
#endif
