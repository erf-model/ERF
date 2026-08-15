/*
 * Implementation 1-moment microphysics model
 * NOTE: this model is based on the SAM code, and the Klemp's paper
 * 1): Joseph, Klemp, the simulation of three-dimensional convective storm dynamics,
 * Journal of the atmospheric sciences, vol35, p1070
 * 2): Marat Khairoutdinov and David Randall, cloud resolving modeling of the ARM summer 1997 IOP:
 * model formulation, results, uncertainties, and sensitivities, Journal of the atmospheric sciences, vol60, p607
 */
#ifndef ERF_SAM_H
#define ERF_SAM_H

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
#include "ERF_SAMUtils.H"
#include "ERF_SurfacePrecipitation.H"

namespace MicVar {
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
      // derived vars
      rain_accum,
      snow_accum,
      graup_accum,
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
    void IceFall (const SolverChoice& sc);

    // precip
    void Precip (const SolverChoice& sc);

    // precip fall
    void PrecipFall (const SolverChoice& sc);

    // Set up for first time
    void
    Define (SolverChoice& sc) override
    {
        m_moisture_type = sc.moisture_type;
        m_fac_cond = lcond / sc.c_p;
        m_fac_fus  = lfus / sc.c_p;
        m_fac_sub  = lsub / sc.c_p;
        m_axis     = sc.ave_plane;
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

    // Copy state into micro vars
    void
    Copy_Micro_to_State (amrex::MultiFab& cons_in) override;

    using NullMoist::Update_Micro_Vars;

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
             const SolverChoice& sc) override
    {
        dt = dt_advance;

        this->Cloud(sc);
        this->IceFall(sc);
        this->Precip(sc);
        this->PrecipFall(sc);
    }

    amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    const amrex::MultiFab*
    Qmoist_Ptr (const int& varIdx) const
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    void
    Compute_Coefficients ();

    SAMCoefficientRow
    CoefficientRowAt (int k) const;

    int
    Qmoist_Size () override { return SAM::m_qmoist_size; }

    int
    Qstate_Moist_Size () override { return SAM::n_qstate_moist_size; }

    int
    Qstate_Moist_NumConc_Size () override { return SAM::n_qstate_moist_numconc_size; }

    void
    Qmoist_Restart_Vars (const SolverChoice& /*a_sc*/,
                         std::vector<int>& a_idx,
                         std::vector<std::string>& a_names) const override
    {
        a_idx.clear();
        a_names.clear();

        // NOTE: These are the indices to access into qmoist (not mic_fab_vars)
        a_idx.push_back(0); a_names.push_back("RainAccum");
        a_idx.push_back(1); a_names.push_back("SnowAccum");
        a_idx.push_back(2); a_names.push_back("GraupAccum");
    }

    SurfacePrecipAccumulationSources
    Get_Surface_Precip_Accumulation_Ptrs (const int&) const override
    {
        if (m_moisture_type == MoistureType::SAM_NoPrecip_NoIce) {
            return {};
        }

        SurfacePrecipAccumulationSources sources;
        sources.rain = {mic_fab_vars[MicVar::rain_accum].get(), rhor / amrex::Real(1000.0)};

        if (m_moisture_type == MoistureType::SAM) {
            sources.snow    = {mic_fab_vars[MicVar::snow_accum].get(), rhos / amrex::Real(1000.0)};
            sources.graupel = {mic_fab_vars[MicVar::graup_accum].get(), rhog / amrex::Real(1000.0)};
        }

        return sources;
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static amrex::Real
    NewtonIterSat (int& i, int& j, int& k,
                   const int& SAM_moisture_type,
                   const amrex::Real& fac_cond,
                   const amrex::Real& fac_fus,
                   const amrex::Real& /*fac_sub*/,
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
#ifdef AMREX_USE_FLOAT
        constexpr amrex::Real tol = amrex::Real(1.e-4);
#else
        constexpr amrex::Real tol = amrex::Real(1.e-8);
#endif
        constexpr int max_iter = 20;

        // Saturation moisture fractions
        amrex::Real omn, domn;
        amrex::Real qsatm, dqsatm;
        amrex::Real qsatw, dqsatw;
        amrex::Real qsati, dqsati;
        amrex::Real qn, qif, dqif;

        // Newton iteration vars
        int niter    = 0;
        amrex::Real F, dFdT;

        //=======================================================================
        // Mixed-phase Newton iteration to qv=qsat (cloud phase only)
        //=======================================================================

        // Initial guess for temperature & pressure
        amrex::Real tabs = tabs_array(i,j,k);
        amrex::Real pres = pres_array(i,j,k);

        // Mixed saturation moisture fractions
        erf_qsatw(tabs, pres, qsatw);
        erf_qsati(tabs, pres, qsati);
        omn   = sam_cloud_liquid_fraction(SAM_moisture_type, tabs, an, bn);
        qsatm = sam_mixed_qsat(omn, qsatw, qsati);

        // Saturation derivative WRT temperature
        domn = zero;
        if ( (SAM_moisture_type == 1) &&
             (tabs < tbgmax) && (tabs > tbgmin) ) {
            domn = an;
        }
        erf_dtqsatw(tabs, pres, dqsatw);
        erf_dtqsati(tabs, pres, dqsati);
        dqsatm = sam_mixed_dqsat_dT(omn, domn, qsatw, qsati, dqsatw, dqsati);

        // Enforce omega splitting
        qn   = (qt_array(i,j,k) - qsatm);
        qif  = qn * (one - omn);
        dqif = sam_mixed_dqif_dT(omn, domn, qn, dqsatm);

        // Function for root finding (conserve latent proxy):
        // 0 = -T_new + T_old + L_v/C_p * (qv,0 - qs,m) - L_f/C_p * (qi,0 - qi,f)
        F = sam_newton_residual(tabs, tabs_array(i,j,k),
                                fac_cond, qv_array(i,j,k), qsatm,
                                fac_fus , qi_array(i,j,k), qif  );
        while (std::abs(F) > tol && niter < max_iter) {
            // Derivative of function (T_new iterated on)
            dFdT  = sam_newton_residual_derivative(fac_cond, fac_fus, dqsatm, dqif);

            // Update the temperature
            tabs -= F/dFdT;

            // Mixed saturation moisture fractions
            erf_qsatw(tabs, pres, qsatw);
            erf_qsati(tabs, pres, qsati);
            omn   = sam_cloud_liquid_fraction(SAM_moisture_type, tabs, an, bn);
            qsatm = sam_mixed_qsat(omn, qsatw, qsati);

            // Saturation derivative WRT temperature
            domn = zero;
            if ( (SAM_moisture_type == 1) &&
                 (tabs < tbgmax) && (tabs > tbgmin) ) {
                domn = an;
            }
            erf_dtqsatw(tabs, pres, dqsatw);
            erf_dtqsati(tabs, pres, dqsati);
            dqsatm = sam_mixed_dqsat_dT(omn, domn, qsatw, qsati, dqsatw, dqsati);

            // Enforce omega splitting
            qn   = (qt_array(i,j,k) - qsatm);
            qif  = qn * (one - omn);
            dqif = sam_mixed_dqif_dT(omn, domn, qn, dqsatm);

            // Function for root finding (conserve latent proxy):
            // 0 = -T_new + T_old + L_v/C_p * (qv,0 - qs,m) - L_f/C_p * (qi,0 - qi,f)
            F = sam_newton_residual(tabs, tabs_array(i,j,k),
                                    fac_cond, qv_array(i,j,k), qsatm,
                                    fac_fus , qi_array(i,j,k), qif  );

            // Update iteration
            ++niter;
        }

        // NOTE: the while loop already enforces the iteration limit, so the
        //       residual alone determines whether we found a valid root; testing
        //       niter < max_iter here would spuriously fail a solve that
        //       converged on exactly the last allowed iteration.
        bool valid_solution = (std::abs(F) <= tol);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(valid_solution, "SAM mixed phase saturation adjustment failed to converge.");

        // Partition the change in non-precipitating q
        qv_array(i,j,k) = qsatm;
        qc_array(i,j,k) = qn * omn;
        qi_array(i,j,k) = qn * (one - omn);
        qn_array(i,j,k) = qn;
        qt_array(i,j,k) = qv_array(i,j,k) + qn_array(i,j,k);

        // Return to temperature
        return tabs;
    }

private:
    // Number of qmoist variables (rain_accum, snow_accum, graup_accum)
    int m_qmoist_size = 3;

    // Number of qstate variables
    int n_qstate_moist_size = 6;

    // Number of qstate variables that are number concentrations
    int n_qstate_moist_numconc_size = 0;

    // CFL MAX for vertical advection
    static constexpr amrex::Real CFL_MAX = myhalf;

    // MicVar map (Qmoist indices -> MicVar enum)
    amrex::Vector<int> MicVarMap;

    // geometry
    amrex::Geometry m_geom;

    // timestep
    amrex::Real dt;

    // number of vertical levels
    int nlev, zlo, zhi;

    // plane average axis
    int m_axis;

    // constants
    amrex::Real m_fac_cond;
    amrex::Real m_fac_fus;
    amrex::Real m_fac_sub;
    amrex::Real m_rdOcp;
    bool m_do_cond;
    MoistureType m_moisture_type = MoistureType::None;

    // Minimum dz at this level
    amrex::Real m_dzmin;

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
};
#endif
