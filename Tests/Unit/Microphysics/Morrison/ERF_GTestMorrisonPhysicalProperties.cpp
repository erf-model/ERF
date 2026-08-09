#include <memory>
#include <limits>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include "../../Utils/Microphysics/ERF_GTestMicrophysicsCommon.H"
#include <ERF_EOS.H>
#include "ERF_Morrison.H"

using namespace microphysics_test;

namespace {

constexpr amrex::Real kRdOcp = R_d / Cp_d;

amrex::Real exact_zero_or_near_zero_tol ()
{
#ifdef AMREX_USE_FLOAT
    return amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon();
#else
    return amrex::Real(32.0) * std::numeric_limits<amrex::Real>::epsilon();
#endif
}

class ScopedParmParseString
{
public:
    ScopedParmParseString (const char* prefix,
                           const char* name,
                           const std::string& value)
        : m_pp(prefix),
          m_name(name)
    {
        m_had_previous = m_pp.query(m_name, m_previous);
        m_pp.remove(m_name);
        m_pp.add(m_name, value);
    }

    ~ScopedParmParseString ()
    {
        m_pp.remove(m_name);
        if (m_had_previous) {
            m_pp.add(m_name, m_previous);
        }
    }

private:
    amrex::ParmParse m_pp;
    std::string m_name;
    std::string m_previous;
    bool m_had_previous = false;
};

struct MorrisonCellState {
    amrex::Real rho;
    amrex::Real theta;
    amrex::Real tabs;
    amrex::Real pres_pa;
    amrex::Real qv;
    amrex::Real qc;
    amrex::Real qi;
    amrex::Real qpr;
    amrex::Real qps;
    amrex::Real qpg;
    amrex::Real nc;
    amrex::Real ni;
    amrex::Real nr;
    amrex::Real ns;
    amrex::Real ng;
};

amrex::Geometry make_geometry ()
{
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(1.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return amrex::Geometry(domain, &real_box, 0, is_periodic.data());
}

MorrisonCellState make_morrison_cell_state (const amrex::Real tabs,
                                            const amrex::Real pres_pa,
                                            const amrex::Real qv,
                                            const amrex::Real qc = amrex::Real(0.0),
                                            const amrex::Real qi = amrex::Real(0.0),
                                            const amrex::Real qpr = amrex::Real(0.0),
                                            const amrex::Real qps = amrex::Real(0.0),
                                            const amrex::Real qpg = amrex::Real(0.0),
                                            const amrex::Real nc = amrex::Real(1.0),
                                            const amrex::Real ni = amrex::Real(1.0),
                                            const amrex::Real nr = amrex::Real(1.0),
                                            const amrex::Real ns = amrex::Real(1.0),
                                            const amrex::Real ng = amrex::Real(1.0))
{
    MorrisonCellState state{};
    state.tabs = tabs;
    state.pres_pa = pres_pa;
    state.qv = qv;
    state.qc = qc;
    state.qi = qi;
    state.qpr = qpr;
    state.qps = qps;
    state.qpg = qpg;
    state.nc = nc;
    state.ni = ni;
    state.nr = nr;
    state.ns = ns;
    state.ng = ng;
    state.theta = getThgivenTandP(tabs, pres_pa, kRdOcp);
    state.rho = getRhogivenTandPress(tabs, pres_pa, qv);
    return state;
}

void fill_single_cell_from_morrison_state (amrex::MultiFab& cons,
                                           const MorrisonCellState& state)
{
    cons.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.fabbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k,Rho_comp) = state.rho;
            arr(i,j,k,RhoTheta_comp) = state.rho * state.theta;
            arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
            arr(i,j,k,RhoQ1_comp) = state.rho * state.qv;
            arr(i,j,k,RhoQ2_comp) = state.rho * state.qc;
            arr(i,j,k,RhoQ3_comp) = state.rho * state.qi;
            arr(i,j,k,RhoQ4_comp) = state.rho * state.qpr;
            arr(i,j,k,RhoQ5_comp) = state.rho * state.qps;
            arr(i,j,k,RhoQ6_comp) = state.rho * state.qpg;
            arr(i,j,k,RhoQ7_comp) = state.rho * state.nc;
            arr(i,j,k,RhoQ8_comp) = state.rho * state.ni;
            arr(i,j,k,RhoQ9_comp) = state.rho * state.nr;
            arr(i,j,k,RhoQ10_comp) = state.rho * state.ns;
            arr(i,j,k,RhoQ11_comp) = state.rho * state.ng;
        });
    }

    amrex::Gpu::streamSynchronize();
}

SolverChoice make_morrison_solver_choice ()
{
    SolverChoice sc{};
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.moisture_type = MoistureType::Morrison;
    sc.use_eamxx_shoc = false;
    sc.use_native_shoc = true;
    sc.turbChoice.resize(1);
    sc.turbChoice[0].pbl_type = PBLType::NATIVE_SHOC;
    return sc;
}

} // namespace

// Motivation: The C++ Morrison saturation-adjustment block must respect the
// SHOC-family condensation-suppression contract. With SHOC active, the public
// C++ path should leave a warm supersaturated cloud-free state unchanged.
TEST(MorrisonPhysicalProperties, NativeShocSuppressesCppSaturationAdjustment)
{
    [[maybe_unused]] ScopedParmParseString use_cpp("erf", "use_morr_cpp_answer", "true");

    const amrex::Geometry geom = make_geometry();
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ11_comp + 1, 3);
    amrex::MultiFab before(ba, dm, RhoQ11_comp + 1, 3);

    const amrex::Real tabs = amrex::Real(280.0);
    const amrex::Real pres_pa = amrex::Real(90000.0);
    amrex::Real qsatw = amrex::Real(0.0);
    erf_qsatw(tabs, pres_pa * amrex::Real(0.01), qsatw);

    const MorrisonCellState state = make_morrison_cell_state(
        tabs, pres_pa, qsatw + amrex::Real(1.0e-4));

    fill_single_cell_from_morrison_state(cons, state);
    amrex::MultiFab::Copy(before, cons, 0, 0, cons.nComp(), cons.nGrowVect());

    Morrison morrison;
    SolverChoice sc = make_morrison_solver_choice();
    morrison.Define(sc);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    morrison.Set_dzmin(geom.CellSize(2));
    morrison.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
    morrison.Copy_State_to_Micro(cons);
    morrison.Advance(amrex::Real(1.0), sc);
    morrison.Copy_Micro_to_State(cons);
    amrex::Gpu::streamSynchronize();

    for (int comp = Rho_comp; comp <= RhoQ6_comp; ++comp) {
        SCOPED_TRACE("comp=" + std::to_string(comp));
        EXPECT_NEAR(cons.max(comp), before.max(comp), exact_zero_or_near_zero_tol());
    }
}
