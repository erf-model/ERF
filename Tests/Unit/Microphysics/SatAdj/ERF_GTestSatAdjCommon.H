#ifndef ERF_GTEST_SATADJ_COMMON_H_
#define ERF_GTEST_SATADJ_COMMON_H_

#include <algorithm>
#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_MicrophysicsUtils.H>
#include <ERF_SatAdj.H>

namespace satadj_test {

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kThermoTolFactor = amrex::Real(2.0e-5);
constexpr amrex::Real kStateTolFactor = amrex::Real(5.0e-5);
constexpr amrex::Real kAsymptoticRelTol = amrex::Real(2.0e-1);
#else
constexpr amrex::Real kThermoTolFactor = amrex::Real(1.0e-12);
constexpr amrex::Real kStateTolFactor = amrex::Real(1.0e-12);
constexpr amrex::Real kAsymptoticRelTol = amrex::Real(2.0e-2);
#endif

constexpr amrex::Real kFacCond = lcond / Cp_d;
constexpr amrex::Real kRdOcp = R_d / Cp_d;

struct CellState {
    amrex::Real tabs;
    amrex::Real pres_mbar;
    amrex::Real theta;
    amrex::Real qv;
    amrex::Real qc;
};

struct ConservedState {
    amrex::Real rho;
    amrex::Real rhotheta;
    amrex::Real rhoqv;
    amrex::Real rhoqc;
};

inline amrex::Real scaled_tol (const amrex::Real reference, const amrex::Real factor)
{
    return factor * std::max(amrex::Real(1), std::abs(reference));
}

inline amrex::Real qsat (const amrex::Real tabs, const amrex::Real pres_mbar)
{
    amrex::Real qsat_local;
    erf_qsatw(tabs, pres_mbar, qsat_local);
    return qsat_local;
}

inline amrex::Real dqsat_dt (const amrex::Real tabs, const amrex::Real pres_mbar)
{
    amrex::Real dqsat_local;
    erf_dtqsatw(tabs, pres_mbar, dqsat_local);
    return dqsat_local;
}

inline CellState make_cell_state (const amrex::Real tabs,
                                  const amrex::Real pres_mbar,
                                  const amrex::Real qv,
                                  const amrex::Real qc)
{
    return CellState{tabs,
                     pres_mbar,
                     getThgivenTandP(tabs, amrex::Real(100.0) * pres_mbar, kRdOcp),
                     qv,
                     qc};
}

inline void adjust (CellState& state)
{
    SatAdj::AdjustSatAdjCell(kFacCond,
                             kRdOcp,
                             state.tabs,
                             state.pres_mbar,
                             state.theta,
                             state.qv,
                             state.qc);
}

inline ConservedState make_conserved_state (const CellState& state)
{
    const amrex::Real pres_pa = amrex::Real(100.0) * state.pres_mbar;
    const amrex::Real rho = getRhogivenTandPress(state.tabs, pres_pa, state.qv);
    const amrex::Real theta = getThgivenTandP(state.tabs, pres_pa, kRdOcp);

    return ConservedState{rho, rho * theta, rho * state.qv, rho * state.qc};
}

inline CellState make_state_from_conserved (const amrex::Real rho,
                                            const amrex::Real rhotheta,
                                            const amrex::Real rhoqv,
                                            const amrex::Real rhoqc)
{
    const amrex::Real qv = rhoqv / rho;
    return CellState{getTgivenRandRTh(rho, rhotheta, qv),
                     getPgivenRTh(rhotheta, qv) * amrex::Real(0.01),
                     rhotheta / rho,
                     qv,
                     rhoqc / rho};
}

inline SolverChoice make_solver_choice (const bool use_shoc)
{
    SolverChoice sc;
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.use_shoc = use_shoc;
    return sc;
}

inline bool find_evaporation_then_recondensation_state (CellState& state)
{
    const amrex::Real tabs = amrex::Real(285.0);
    const amrex::Real pres_mbar = amrex::Real(900.0);
    const amrex::Real qsat_initial = qsat(tabs, pres_mbar);
    const amrex::Real qv = amrex::Real(0.8) * qsat_initial;

    for (int iqc = 1; iqc <= 80; ++iqc) {
        const amrex::Real qc = amrex::Real(iqc) * amrex::Real(1.0e-4);
        if ((qv + qc) >= qsat_initial) {
            break;
        }

        const amrex::Real tabs_all_evap = tabs - kFacCond * qc;
        if ((qv + qc) > qsat(tabs_all_evap, pres_mbar)) {
            state = make_cell_state(tabs, pres_mbar, qv, qc);
            return true;
        }
    }

    return false;
}

inline amrex::Geometry make_geometry (const int nx, const int ny, const int nz)
{
    amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                      amrex::IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1)));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(1, 1, 1)};
    return amrex::Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodicity.data());
}

} // namespace satadj_test

#endif