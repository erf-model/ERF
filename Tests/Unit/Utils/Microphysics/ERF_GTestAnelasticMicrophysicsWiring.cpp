#include <cmath>
#include <memory>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>

#include <gtest/gtest.h>

#include <ERF_DataStruct.H>
#include <ERF_IndexDefines.H>
#include <ERF_Kessler.H>
#include <ERF_Morrison.H>
#include <ERF_SAM.H>
#include <ERF_WDM6.H>
#include <ERF_WSM6.H>

#include "ERF_GTestMicrophysicsCommon.H"

namespace {

using amrex::Array4;
using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::Real;
using namespace microphysics_test;

constexpr Real kRdOcp = RdoCp;

Geometry make_test_geometry ()
{
    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                  {AMREX_D_DECL(2.0, 1.0, 1.0)});
    amrex::Array<int, AMREX_SPACEDIM> periodicity{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, &real_box, amrex::CoordSys::cartesian, periodicity.data());
}

void initialize_test_state (MultiFab& states, MultiFab& base_state)
{
    states.setVal(Real(0.0));
    base_state.setVal(Real(-12345.0));

    for (MFIter mfi(states); mfi.isValid(); ++mfi) {
        const Box box = mfi.growntilebox();
        const auto state = states.array(mfi);
        const auto base = base_state.array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int cell = amrex::max(0, amrex::min(1, i));
            const Real rho = cell == 0 ? Real(1.10) : Real(0.90);
            const Real theta = cell == 0 ? Real(303.0) : Real(287.0);
            const Real p0 = cell == 0 ? Real(90000.0) : Real(82000.0);
            const Real qv = cell == 0 ? Real(0.0030) : Real(0.0070);
            const Real qc = cell == 0 ? Real(0.0010) : Real(0.0020);
            const Real qi = cell == 0 ? Real(0.0004) : Real(0.0008);
            const Real qr = cell == 0 ? Real(0.0007) : Real(0.0011);
            const Real qs = cell == 0 ? Real(0.0002) : Real(0.0005);
            const Real qg = cell == 0 ? Real(0.0003) : Real(0.0006);

            state(i,j,k,Rho_comp) = rho;
            state(i,j,k,RhoTheta_comp) = rho * theta;
            state(i,j,k,RhoQ1_comp) = rho * qv;
            state(i,j,k,RhoQ2_comp) = rho * qc;
            state(i,j,k,RhoQ3_comp) = rho * qi;
            state(i,j,k,RhoQ4_comp) = rho * qr;
            state(i,j,k,RhoQ5_comp) = rho * qs;
            state(i,j,k,RhoQ6_comp) = rho * qg;
            state(i,j,k,RhoQ7_comp) = rho * Real(2.0);
            state(i,j,k,RhoQ8_comp) = Real(0.0);
            state(i,j,k,RhoQ9_comp) = rho * Real(3.0);
            state(i,j,k,RhoQ10_comp) = rho * Real(4.0);
            state(i,j,k,RhoQ11_comp) = rho * Real(5.0);

            // Only p0 is a physical input to the anelastic microphysics path.
            // The other base-state fields deliberately remain sentinels.
            base(i,j,k,BaseState::p0_comp) = p0;
        });
    }

    amrex::Gpu::streamSynchronize();
}

struct TestState {
    BoxArray boxes;
    DistributionMapping dm;
    MultiFab states;
    MultiFab base_state;

    TestState ()
        : boxes(Box(IntVect(0, 0, 0), IntVect(1, 0, 0))),
          dm(boxes),
          states(boxes, dm, RhoQ11_comp + 1, 1),
          base_state(boxes, dm, BaseState::num_comps, 1)
    {
        initialize_test_state(states, base_state);
    }
};

struct WorkingArrays {
    MultiFab rho;
    MultiFab theta;
    MultiFab qv;
    MultiFab qc;
    MultiFab qi;
    MultiFab qn;
    MultiFab qt;
    MultiFab qpr;
    MultiFab qps;
    MultiFab qpg;
    MultiFab qp;
    MultiFab tabs;
    MultiFab pres;

    explicit WorkingArrays (const TestState& input)
        : rho(input.boxes, input.dm, 1, 1),
          theta(input.boxes, input.dm, 1, 1),
          qv(input.boxes, input.dm, 1, 1),
          qc(input.boxes, input.dm, 1, 1),
          qi(input.boxes, input.dm, 1, 1),
          qn(input.boxes, input.dm, 1, 1),
          qt(input.boxes, input.dm, 1, 1),
          qpr(input.boxes, input.dm, 1, 1),
          qps(input.boxes, input.dm, 1, 1),
          qpg(input.boxes, input.dm, 1, 1),
          qp(input.boxes, input.dm, 1, 1),
          tabs(input.boxes, input.dm, 1, 1),
          pres(input.boxes, input.dm, 1, 1)
    {
        set_all(Real(-999.0));
    }

    void set_all (const Real value)
    {
        rho.setVal(value);
        theta.setVal(value);
        qv.setVal(value);
        qc.setVal(value);
        qi.setVal(value);
        qn.setVal(value);
        qt.setVal(value);
        qpr.setVal(value);
        qps.setVal(value);
        qpg.setVal(value);
        qp.setVal(value);
        tabs.setVal(value);
        pres.setVal(value);
    }
};

void copy_kessler (const TestState& input, WorkingArrays& work)
{
    for (MFIter mfi(input.states); mfi.isValid(); ++mfi) {
        const Box box = mfi.growntilebox();
        const auto states = input.states.const_array(mfi);
        const auto base = input.base_state.const_array(mfi);
        const auto rho = work.rho.array(mfi);
        const auto theta = work.theta.array(mfi);
        const auto qv = work.qv.array(mfi);
        const auto qc = work.qc.array(mfi);
        const auto qp = work.qp.array(mfi);
        const auto qt = work.qt.array(mfi);
        const auto tabs = work.tabs.array(mfi);
        const auto pres = work.pres.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            kessler_copy_state_to_micro_cell(
                states, base, rho, theta, qv, qc, qp, qt, tabs, pres,
                kRdOcp, true, i, j, k);
        });
    }
    amrex::Gpu::streamSynchronize();
}

void copy_sam (const TestState& input, WorkingArrays& work)
{
    for (MFIter mfi(input.states); mfi.isValid(); ++mfi) {
        const Box box = mfi.growntilebox();
        const auto states = input.states.const_array(mfi);
        const auto base = input.base_state.const_array(mfi);
        const auto rho = work.rho.array(mfi);
        const auto theta = work.theta.array(mfi);
        const auto qv = work.qv.array(mfi);
        const auto qc = work.qc.array(mfi);
        const auto qi = work.qi.array(mfi);
        const auto qn = work.qn.array(mfi);
        const auto qt = work.qt.array(mfi);
        const auto qpr = work.qpr.array(mfi);
        const auto qps = work.qps.array(mfi);
        const auto qpg = work.qpg.array(mfi);
        const auto qp = work.qp.array(mfi);
        const auto tabs = work.tabs.array(mfi);
        const auto pres = work.pres.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            sam_copy_state_to_micro_cell(
                states, base, rho, theta, qv, qc, qi, qn, qt, qpr, qps,
                qpg, qp, tabs, pres, kRdOcp, true, i, j, k);
        });
    }
    amrex::Gpu::streamSynchronize();
}

void copy_morrison (const TestState& input, WorkingArrays& work)
{
    for (MFIter mfi(input.states); mfi.isValid(); ++mfi) {
        const Box box = mfi.growntilebox();
        const auto states = input.states.const_array(mfi);
        const auto base = input.base_state.const_array(mfi);
        const auto rho = work.rho.array(mfi);
        const auto theta = work.theta.array(mfi);
        const auto qv = work.qv.array(mfi);
        const auto qc = work.qc.array(mfi);
        const auto qi = work.qi.array(mfi);
        const auto qn = work.qn.array(mfi);
        const auto qt = work.qt.array(mfi);
        const auto qpr = work.qpr.array(mfi);
        const auto qps = work.qps.array(mfi);
        const auto qpg = work.qpg.array(mfi);
        const auto qp = work.qp.array(mfi);
        const auto tabs = work.tabs.array(mfi);
        const auto pres = work.pres.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            morrison_copy_state_to_micro_cell(
                states, base, rho, theta, qv, qc, qi, qn, qt, qpr, qps,
                qpg, qp, tabs, pres, kRdOcp, true, i, j, k);
        });
    }
    amrex::Gpu::streamSynchronize();
}

void copy_wsm6 (const TestState& input, WorkingArrays& work)
{
    for (MFIter mfi(input.states); mfi.isValid(); ++mfi) {
        const Box box = mfi.growntilebox();
        const auto states = input.states.const_array(mfi);
        const auto base = input.base_state.const_array(mfi);
        const auto rho = work.rho.array(mfi);
        const auto theta = work.theta.array(mfi);
        const auto tabs = work.tabs.array(mfi);
        const auto pres = work.pres.array(mfi);
        const auto qv = work.qv.array(mfi);
        const auto qc = work.qc.array(mfi);
        const auto qi = work.qi.array(mfi);
        const auto qr = work.qpr.array(mfi);
        const auto qs = work.qps.array(mfi);
        const auto qg = work.qpg.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            wsm6_copy_state_to_micro_cell(
                states, base, rho, theta, tabs, pres, qv, qc, qi, qr, qs,
                qg, kRdOcp, true, i, j, k);
        });
    }
    amrex::Gpu::streamSynchronize();
}

void copy_wdm6 (const TestState& input, WorkingArrays& work)
{
    for (MFIter mfi(input.states); mfi.isValid(); ++mfi) {
        const Box box = mfi.growntilebox();
        const auto states = input.states.const_array(mfi);
        const auto base = input.base_state.const_array(mfi);
        const auto rho = work.rho.array(mfi);
        const auto theta = work.theta.array(mfi);
        const auto tabs = work.tabs.array(mfi);
        const auto pres = work.pres.array(mfi);
        const auto qv = work.qv.array(mfi);
        const auto qc = work.qc.array(mfi);
        const auto qi = work.qi.array(mfi);
        const auto qr = work.qpr.array(mfi);
        const auto qs = work.qps.array(mfi);
        const auto qg = work.qpg.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            wdm6_copy_state_to_micro_cell(
                states, base, rho, theta, tabs, pres, qv, qc, qi, qr, qs,
                qg, kRdOcp, true, i, j, k);
        });
    }
    amrex::Gpu::streamSynchronize();
}

void expect_reference_pressure_diagnosis (const TestState& input,
                                         const WorkingArrays& work,
                                         const Real pressure_scale)
{
    MultiFab errors(input.boxes, input.dm, 4, 0);
    errors.setVal(Real(0.0));

    for (MFIter mfi(input.states); mfi.isValid(); ++mfi) {
        const Box box = mfi.validbox();
        const auto states = input.states.const_array(mfi);
        const auto base = input.base_state.const_array(mfi);
        const auto rho = work.rho.const_array(mfi);
        const auto theta = work.theta.const_array(mfi);
        const auto tabs = work.tabs.const_array(mfi);
        const auto pres = work.pres.const_array(mfi);
        const auto error = errors.array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real expected_rho = states(i,j,k,Rho_comp);
            const Real expected_theta = states(i,j,k,RhoTheta_comp) / expected_rho;
            const Real expected_p0 = base(i,j,k,BaseState::p0_comp);
            const Real expected_tabs = getTgivenPandTh(expected_p0, expected_theta, kRdOcp);
            error(i,j,k,0) = normalized_error(rho(i,j,k), expected_rho, kValueRelTol);
            error(i,j,k,1) = normalized_error(theta(i,j,k), expected_theta, kValueRelTol);
            error(i,j,k,2) = normalized_error(tabs(i,j,k), expected_tabs, kValueRelTol);
            error(i,j,k,3) = normalized_error(pres(i,j,k), expected_p0 * pressure_scale,
                                               kValueRelTol);
        });
    }

    amrex::Gpu::streamSynchronize();
    for (int comp = 0; comp < errors.nComp(); ++comp) {
        EXPECT_LE(errors.max(comp), Real(20.0));
    }
    EXPECT_LT(work.tabs.min(0), work.tabs.max(0));
    EXPECT_LT(work.pres.min(0), work.pres.max(0));
}

SolverChoice make_scheme_choice (const MoistureType moisture_type,
                                 const bool anelastic = false)
{
    SolverChoice sc{};
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.moisture_type = moisture_type;
    sc.anelastic = {anelastic ? 1 : 0};
    sc.use_eamxx_shoc = false;
    sc.use_native_shoc = false;
    sc.ave_plane = 2;
    return sc;
}

void expect_coefficients_near (const SAMCoefficientRow& actual,
                               const SAMCoefficientRow& expected)
{
    constexpr Real coefficient_tol = Real(1.0e-10);
    EXPECT_NEAR(actual.accrrc, expected.accrrc,
                scaled_tol(actual.accrrc, expected.accrrc, coefficient_tol));
    EXPECT_NEAR(actual.accrsi, expected.accrsi,
                scaled_tol(actual.accrsi, expected.accrsi, coefficient_tol));
    EXPECT_NEAR(actual.accrsc, expected.accrsc,
                scaled_tol(actual.accrsc, expected.accrsc, coefficient_tol));
    EXPECT_NEAR(actual.coefice, expected.coefice,
                scaled_tol(actual.coefice, expected.coefice, coefficient_tol));
    EXPECT_NEAR(actual.evaps1, expected.evaps1,
                scaled_tol(actual.evaps1, expected.evaps1, coefficient_tol));
    EXPECT_NEAR(actual.evaps2, expected.evaps2,
                scaled_tol(actual.evaps2, expected.evaps2, coefficient_tol));
    EXPECT_NEAR(actual.accrgi, expected.accrgi,
                scaled_tol(actual.accrgi, expected.accrgi, coefficient_tol));
    EXPECT_NEAR(actual.accrgc, expected.accrgc,
                scaled_tol(actual.accrgc, expected.accrgc, coefficient_tol));
    EXPECT_NEAR(actual.evapg1, expected.evapg1,
                scaled_tol(actual.evapg1, expected.evapg1, coefficient_tol));
    EXPECT_NEAR(actual.evapg2, expected.evapg2,
                scaled_tol(actual.evapg2, expected.evapg2, coefficient_tol));
    EXPECT_NEAR(actual.evapr1, expected.evapr1,
                scaled_tol(actual.evapr1, expected.evapr1, coefficient_tol));
    EXPECT_NEAR(actual.evapr2, expected.evapr2,
                scaled_tol(actual.evapr2, expected.evapr2, coefficient_tol));
}

} // namespace

// Motivation: Each supported Eulerian scheme must receive the anelastic
// reference pressure instead of rebuilding pressure from rho and rho*theta.
TEST(AnelasticMicrophysicsWiring, KesslerCopyInUsesReferencePressure)
{
    TestState input;
    WorkingArrays work(input);
    copy_kessler(input, work);
    expect_reference_pressure_diagnosis(input, work, Real(0.01));
}

// Motivation: SAM stores pressure in mbar, but its anelastic temperature must
// still be diagnosed from the BaseState pressure in Pa.
TEST(AnelasticMicrophysicsWiring, SAMCopyInUsesReferencePressure)
{
    TestState input;
    WorkingArrays work(input);
    copy_sam(input, work);
    expect_reference_pressure_diagnosis(input, work, Real(0.01));
}

// Motivation: Morrison's source kernels consume pressure in Pa and must see
// the same BaseState pressure used to diagnose their absolute temperature.
TEST(AnelasticMicrophysicsWiring, MorrisonCopyInUsesReferencePressure)
{
    TestState input;
    WorkingArrays work(input);
    copy_morrison(input, work);
    expect_reference_pressure_diagnosis(input, work, Real(1.0));
}

// Motivation: WSM6 copy-in and copy-out must use the same held pressure so an
// anelastic source update cannot project theta through a rho-based EOS.
TEST(AnelasticMicrophysicsWiring, WSM6CopyInUsesReferencePressure)
{
    TestState input;
    WorkingArrays work(input);
    copy_wsm6(input, work);
    expect_reference_pressure_diagnosis(input, work, Real(1.0));

    SolverChoice sc = make_scheme_choice(MoistureType::WSM6, true);
    WSM6 wsm6;
    wsm6.SetCurrentLevel(0);
    wsm6.Define(sc);
    std::unique_ptr<MultiFab> z_phys_nd;
    std::unique_ptr<MultiFab> detJ_cc;
    const Geometry geom = make_test_geometry();
    wsm6.Init(input.states, input.boxes, geom, Real(1.0), z_phys_nd, detJ_cc);
    wsm6.Copy_State_to_Micro(input.states, &input.base_state);
    wsm6.Copy_Micro_to_State(input.states);
    amrex::Gpu::streamSynchronize();

    EXPECT_NEAR(input.states.max(RhoTheta_comp), Real(1.10) * Real(303.0),
                scaled_tol(Real(1.10) * Real(303.0), Real(1.0), Real(1.0e-11)));
}

// Motivation: WDM6 must receive p0 during copy-in while retaining the scheme's
// existing pressure units and moisture-component mapping.
TEST(AnelasticMicrophysicsWiring, WDM6CopyInUsesReferencePressure)
{
    TestState input;
    WorkingArrays work(input);
    copy_wdm6(input, work);
    expect_reference_pressure_diagnosis(input, work, Real(1.0));

    const Real rho_theta_before = input.states.max(RhoTheta_comp);
    SolverChoice sc = make_scheme_choice(MoistureType::WDM6, true);
    WDM6 wdm6;
    wdm6.SetCurrentLevel(0);
    wdm6.Define(sc);
    std::unique_ptr<MultiFab> z_phys_nd;
    std::unique_ptr<MultiFab> detJ_cc;
    const Geometry geom = make_test_geometry();
    wdm6.Init(input.states, input.boxes, geom, Real(1.0), z_phys_nd, detJ_cc);
    wdm6.Copy_State_to_Micro(input.states, &input.base_state);
    wdm6.Copy_Micro_to_State(input.states);
    amrex::Gpu::streamSynchronize();

    // With no Advance call, copy-in/copy-out must not create a thermodynamic
    // tendency merely because the two cells use different reference pressures.
    EXPECT_NEAR(input.states.max(RhoTheta_comp), rho_theta_before,
                scaled_tol(rho_theta_before, rho_theta_before, Real(1.0e-11)));
}

// Motivation: WSM6's compressible copy-out must retain the historical RdoCp
// exponent, while anelastic copy-out must use the configured SolverChoice
// exponent with the held reference pressure.
TEST(AnelasticMicrophysicsWiring, WSM6CustomCpKeepsModeSpecificExponent)
{
    constexpr Real custom_cp = Real(900.0);

    {
        TestState input;
        SolverChoice sc = make_scheme_choice(MoistureType::WSM6, false);
        sc.c_p = custom_cp;
        sc.rdOcp = R_d / custom_cp;

        WSM6 wsm6;
        wsm6.SetCurrentLevel(0);
        wsm6.Define(sc);
        std::unique_ptr<MultiFab> z_phys_nd;
        std::unique_ptr<MultiFab> detJ_cc;
        const Geometry geom = make_test_geometry();
        wsm6.Init(input.states, input.boxes, geom, Real(1.0), z_phys_nd, detJ_cc);
        wsm6.Copy_State_to_Micro(input.states);
        wsm6.Copy_Micro_to_State(input.states);
        amrex::Gpu::streamSynchronize();

        EXPECT_NEAR(input.states.min(RhoTheta_comp), Real(0.90) * Real(287.0),
                    scaled_tol(Real(0.90) * Real(287.0), Real(1.0), Real(1.0e-11)));
        EXPECT_NEAR(input.states.max(RhoTheta_comp), Real(1.10) * Real(303.0),
                    scaled_tol(Real(1.10) * Real(303.0), Real(1.0), Real(1.0e-11)));

        const Real temperature = Real(279.0);
        const Real pressure = Real(87000.0);
        const Real expected = getThgivenRandT(
            Real(1.0), temperature, RdoCp, Real(0.004));
        const Real configured = getThgivenRandT(
            Real(1.0), temperature, sc.rdOcp, Real(0.004));
        const Real actual = wsm6_theta_from_temperature_and_pressure(
            Real(1.0), temperature, pressure, Real(0.004), sc.rdOcp, false);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, Real(1.0e-12)));
        EXPECT_GT(std::abs(expected - configured), Real(1.0e-2));
    }

    {
        TestState input;
        SolverChoice sc = make_scheme_choice(MoistureType::WSM6, true);
        sc.c_p = custom_cp;
        sc.rdOcp = R_d / custom_cp;

        WSM6 wsm6;
        wsm6.SetCurrentLevel(0);
        wsm6.Define(sc);
        std::unique_ptr<MultiFab> z_phys_nd;
        std::unique_ptr<MultiFab> detJ_cc;
        const Geometry geom = make_test_geometry();
        wsm6.Init(input.states, input.boxes, geom, Real(1.0), z_phys_nd, detJ_cc);
        wsm6.Copy_State_to_Micro(input.states, &input.base_state);
        wsm6.Copy_Micro_to_State(input.states);
        amrex::Gpu::streamSynchronize();

        EXPECT_NEAR(input.states.min(RhoTheta_comp), Real(0.90) * Real(287.0),
                    scaled_tol(Real(0.90) * Real(287.0), Real(1.0), Real(1.0e-11)));
        EXPECT_NEAR(input.states.max(RhoTheta_comp), Real(1.10) * Real(303.0),
                    scaled_tol(Real(1.10) * Real(303.0), Real(1.0), Real(1.0e-11)));

        const Real temperature = Real(279.0);
        const Real pressure = Real(87000.0);
        const Real expected = getThgivenTandP(temperature, pressure, sc.rdOcp);
        const Real actual = wsm6_theta_from_temperature_and_pressure(
            Real(1.0), temperature, pressure, Real(0.004), sc.rdOcp, true);
        EXPECT_NEAR(actual, expected, scaled_tol(actual, expected, Real(1.0e-12)));
    }
}

// Motivation: Compressible calls must keep the established rho/rho*theta EOS
// diagnosis and must not accidentally consume a caller's base-state pressure.
TEST(AnelasticMicrophysicsWiring, CompressibleDiagnosisIsPreserved)
{
    const Real rho = Real(1.2);
    const Real theta = Real(300.0);
    const Real rho_theta = rho * theta;
    const Real qv = Real(0.004);
    const Real expected_pressure = getPgivenRTh(rho_theta, qv);
    const Real expected_temperature = getTgivenRandRTh(rho, rho_theta, qv);

    const MicrophysicsThermoState thermo = diagnose_microphysics_thermo_state(
        rho, rho_theta, qv, kRdOcp, false, Real(82000.0));
    EXPECT_NEAR(thermo.pressure_pa, expected_pressure,
                scaled_tol(expected_pressure, expected_pressure, Real(1.0e-12)));
    EXPECT_NEAR(thermo.temperature, expected_temperature,
                scaled_tol(expected_temperature, expected_temperature, Real(1.0e-12)));

    const SAMPrimitiveCell sam_state = sam_cons_to_primitive(
        rho, rho_theta, rho * qv, Real(0.0), Real(0.0), Real(0.0), Real(0.0), Real(0.0));
    EXPECT_NEAR(sam_state.pres_mbar, Real(0.01) * expected_pressure,
                scaled_tol(expected_pressure, expected_pressure, Real(1.0e-12)));
    EXPECT_NEAR(sam_state.tabs, expected_temperature,
                scaled_tol(expected_temperature, expected_temperature, Real(1.0e-12)));
}

// Motivation: SAM's compressible coefficient path must average rho, theta, and
// qv before applying the nonlinear EOS, rather than averaging diagnosed T. A
// non-null base-state pointer must not override the configured compressible mode.
TEST(AnelasticMicrophysicsWiring, SAMHeterogeneousCompressiblePlaneIgnoresBasePointer)
{
    TestState input;
    SolverChoice sc = make_scheme_choice(MoistureType::SAM);
    SAM sam;
    sam.SetCurrentLevel(0);
    sam.Define(sc);
    std::unique_ptr<MultiFab> z_phys_nd;
    std::unique_ptr<MultiFab> detJ_cc;
    const Geometry geom = make_test_geometry();
    sam.Init(input.states, input.boxes, geom, Real(1.0), z_phys_nd, detJ_cc);
    sam.Update_Micro_Vars(input.states, &input.base_state);

    const Real rho_bar = Real(1.0);
    const Real theta_bar = Real(295.0);
    const Real qv_bar = Real(0.005);
    const Real rho_theta_bar = rho_bar * theta_bar;
    const Real averaged_state_temperature = getTgivenRandRTh(
        rho_bar, rho_theta_bar, qv_bar);
    const SAMCoefficientRow expected = sam_compute_coefficient_row(
        rho_bar, amrex::min(averaged_state_temperature, Real(273.16)),
        erf_gammafff(three + b_rain),
        erf_gammafff((Real(5.0) + b_rain) / two),
        erf_gammafff(three + b_snow),
        erf_gammafff((Real(5.0) + b_snow) / two),
        erf_gammafff(three + b_grau),
        erf_gammafff((Real(5.0) + b_grau) / two));

    expect_coefficients_near(sam.CoefficientRowAt(0), expected);

    const Real cell_temperature_0 = getTgivenRandRTh(
        Real(1.10), Real(1.10) * Real(303.0), Real(0.0030));
    const Real cell_temperature_1 = getTgivenRandRTh(
        Real(0.90), Real(0.90) * Real(287.0), Real(0.0070));
    EXPECT_GT(std::abs(averaged_state_temperature -
                       Real(0.5) * (cell_temperature_0 + cell_temperature_1)),
              Real(1.0e-2));
}

// Motivation: WDM6 uses the same held pressure in both conversion paths but
// retains RdoCp for compressible runs and uses configured rdOcp for anelastic.
TEST(AnelasticMicrophysicsWiring, WDM6CustomCpConversionIsModeSpecific)
{
    constexpr Real custom_cp = Real(900.0);
    const Real configured_rdOcp = R_d / custom_cp;
    const Real temperature = Real(279.0);
    const Real pressure = Real(87000.0);
    const Real expected_compressible = getThgivenTandP(
        temperature, pressure, RdoCp);
    const Real expected_anelastic = getThgivenTandP(
        temperature, pressure, configured_rdOcp);

    const Real actual_compressible = wdm6_theta_from_temperature_and_pressure(
        temperature, pressure, configured_rdOcp, false);
    const Real actual_anelastic = wdm6_theta_from_temperature_and_pressure(
        temperature, pressure, configured_rdOcp, true);

    EXPECT_NEAR(actual_compressible, expected_compressible,
                scaled_tol(actual_compressible, expected_compressible, Real(1.0e-12)));
    EXPECT_NEAR(actual_anelastic, expected_anelastic,
                scaled_tol(actual_anelastic, expected_anelastic, Real(1.0e-12)));
    EXPECT_GT(std::abs(expected_compressible - expected_anelastic), Real(1.0e-2));
}

// Motivation: SAM's heterogeneous anelastic plane must average the diagnosed
// pressure/temperature fields, not average theta and then rediagnose a local
// compressible EOS.
TEST(AnelasticMicrophysicsWiring, SAMHeterogeneousPlaneUsesDiagnosedTemperature)
{
    TestState input;
    SolverChoice sc = make_scheme_choice(MoistureType::SAM, true);
    SAM sam;
    sam.SetCurrentLevel(0);
    sam.Define(sc);
    std::unique_ptr<MultiFab> z_phys_nd;
    std::unique_ptr<MultiFab> detJ_cc;
    const Geometry geom = make_test_geometry();
    sam.Init(input.states, input.boxes, geom, Real(1.0), z_phys_nd, detJ_cc);
    sam.Update_Micro_Vars(input.states, &input.base_state);

    const Real temp0 = getTgivenPandTh(Real(90000.0), Real(303.0), kRdOcp);
    const Real temp1 = getTgivenPandTh(Real(82000.0), Real(287.0), kRdOcp);
    const Real tabs = amrex::min(Real(0.5) * (temp0 + temp1), Real(273.16));
    const SAMCoefficientRow expected = sam_compute_coefficient_row(
        Real(1.0), tabs,
        erf_gammafff(three + b_rain),
        erf_gammafff((Real(5.0) + b_rain) / two),
        erf_gammafff(three + b_snow),
        erf_gammafff((Real(5.0) + b_snow) / two),
        erf_gammafff(three + b_grau),
        erf_gammafff((Real(5.0) + b_grau) / two));

    expect_coefficients_near(sam.CoefficientRowAt(0), expected);
}

// Motivation: The invalid SuperDroplets/anelastic combination must be
// recognized before the simulation reaches a compressible-only source path.
TEST(AnelasticMicrophysicsWiring, SuperDropletsAnelasticConfigurationIsInvalid)
{
    EXPECT_TRUE(anelastic_superdroplets_configuration_invalid(
        MoistureType::SuperDroplets, true));
    EXPECT_FALSE(anelastic_superdroplets_configuration_invalid(
        MoistureType::SuperDroplets, false));
    EXPECT_FALSE(anelastic_superdroplets_configuration_invalid(
        MoistureType::SAM, true));
}
