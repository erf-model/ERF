#include <vector>

#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestKesslerCommon.H"

// These tests exercise the AMReX-portable Kessler helper path. Setup and
// error computation use ParallelFor so the same test code runs in CPU and GPU
// builds. Host-side GTest assertions inspect reduced errors after
// synchronization.

using namespace kessler_test;

namespace {

struct KernelOutputs {
    amrex::Real sat_cond;
    amrex::Real sat_evap;
    amrex::Real source_cond;
    amrex::Real source_evap;
    amrex::Real source_cloud_to_rain;
    amrex::Real source_rain_to_vapor;
    amrex::Real velocity;
    amrex::Real flux;
    amrex::Real face_rho;
    amrex::Real face_qp;
    amrex::Real sedimentation;
    int is_small;
    int substeps;
};

// Keep device lambdas in namespace-scope helpers. CUDA can reject extended
// device lambdas enclosed by GTest's private TestBody() member function.
void launch_helper_kernel (const int ncases,
                           const KernelCase* cases_ptr,
                           KernelOutputs* outputs_ptr)
{
    amrex::ParallelFor(ncases, [=] AMREX_GPU_DEVICE (int idx) noexcept {
        const KernelCase test_case = cases_ptr[idx];
        const KesslerSaturationAdjustment sat = kessler_saturation_adjustment(
            test_case.qv, test_case.qc, test_case.qsat, test_case.dtqsat,
            test_case.do_cond, kSatAdjLatentOverCp);
        const KesslerSourceTerms sources = kessler_warm_rain_sources(
            test_case.qv, test_case.qc, test_case.qp, test_case.rho, test_case.pressure_mbar,
            test_case.qsat, test_case.dtqsat, test_case.dt, test_case.do_cond,
            kSatAdjLatentOverCp);
        const amrex::Real velocity = kessler_terminal_velocity(test_case.rho, test_case.qp);
        const amrex::Real flux = kessler_precip_flux(test_case.rho, velocity, test_case.qp);
        const KesslerFaceState face = kessler_face_state(
            test_case.face_k, test_case.face_k_hi,
            test_case.face_rho_km1, test_case.face_rho_k,
            test_case.face_qp_km1, test_case.face_qp_k);
        const amrex::Real sedimentation = kessler_sedimentation_tendency(
            test_case.fz_hi, test_case.fz_lo, test_case.rho, test_case.dJinv, test_case.coef);
        const bool is_small = kessler_is_small_sedimentation_value(test_case.reduced_value);
        const int substeps = kessler_num_sedimentation_substeps(
            test_case.reduced_value, test_case.substep_dt, test_case.substep_dz);

        outputs_ptr[idx] = KernelOutputs{
            sat.dq_vapor_to_cloud,
            sat.dq_cloud_to_vapor,
            sources.dq_vapor_to_cloud,
            sources.dq_cloud_to_vapor,
            sources.dq_cloud_to_rain,
            sources.dq_rain_to_vapor,
            velocity,
            flux,
            face.rho,
            face.qp,
            sedimentation,
            static_cast<int>(is_small),
            substeps};
    });

    kessler_test::sync();
}

std::vector<KernelCase> make_branch_coverage_cases ()
{
    std::vector<KernelCase> cases = make_kernel_cases();
    cases.push_back(KernelCase{
        amrex::Real(1.0e-2), amrex::Real(1.5e-3), amrex::Real(1.0e-6), amrex::Real(1.05), amrex::Real(900.0),
        amrex::Real(1.1e-2), amrex::Real(0.02), amrex::Real(1.0), true,
        1, 0, 2, amrex::Real(1.02), amrex::Real(1.05), amrex::Real(5.0e-7), amrex::Real(1.0e-6),
        amrex::Real(9.0e-15), amrex::Real(8.0e-15), amrex::Real(1.0), amrex::Real(0.5),
        amrex::Real(9.0e-15), amrex::Real(1.0), amrex::Real(1.0)});
    cases.push_back(KernelCase{
        amrex::Real(1.3e-2), amrex::Real(2.0e-3), amrex::Real(1.0e-2), amrex::Real(1.35), amrex::Real(950.0),
        amrex::Real(1.0e-2), amrex::Real(0.03), amrex::Real(0.75), true,
        2, 0, 2, amrex::Real(1.30), amrex::Real(1.35), amrex::Real(5.0e-3), amrex::Real(1.0e-2),
        amrex::Real(4.0e-2), amrex::Real(5.0e-3), amrex::Real(1.1), amrex::Real(0.375),
        amrex::Real(1.0e-14), amrex::Real(2.0), amrex::Real(0.25)});
    return cases;
}

KernelOutputs host_reference (const KernelCase& test_case)
{
    const KesslerSaturationAdjustment sat = kessler_saturation_adjustment(
        test_case.qv, test_case.qc, test_case.qsat, test_case.dtqsat,
        test_case.do_cond, kSatAdjLatentOverCp);
    const KesslerSourceTerms sources = kessler_warm_rain_sources(
        test_case.qv, test_case.qc, test_case.qp, test_case.rho, test_case.pressure_mbar,
        test_case.qsat, test_case.dtqsat, test_case.dt, test_case.do_cond,
        kSatAdjLatentOverCp);
    const amrex::Real velocity = kessler_terminal_velocity(test_case.rho, test_case.qp);
    const amrex::Real flux = kessler_precip_flux(test_case.rho, velocity, test_case.qp);
    const KesslerFaceState face = kessler_face_state(
        test_case.face_k, test_case.face_k_hi,
        test_case.face_rho_km1, test_case.face_rho_k,
        test_case.face_qp_km1, test_case.face_qp_k);

    return KernelOutputs{
        sat.dq_vapor_to_cloud,
        sat.dq_cloud_to_vapor,
        sources.dq_vapor_to_cloud,
        sources.dq_cloud_to_vapor,
        sources.dq_cloud_to_rain,
        sources.dq_rain_to_vapor,
        velocity,
        flux,
        face.rho,
        face.qp,
        kessler_sedimentation_tendency(test_case.fz_hi, test_case.fz_lo,
                                       test_case.rho, test_case.dJinv, test_case.coef),
        static_cast<int>(kessler_is_small_sedimentation_value(test_case.reduced_value)),
        kessler_num_sedimentation_substeps(test_case.reduced_value,
                                           test_case.substep_dt,
                                           test_case.substep_dz)};
}

} // namespace

// Motivation: Production Kessler calls these helpers from AMReX kernels. This
// host/device equivalence sweep covers saturation no-op, condensation, cloud
// evaporation, rain evaporation, autoconversion off/on, face branches,
// threshold branches, and integer substep outputs across deterministic cases.
TEST(KesslerKernel, HostDeviceHelperEquivalenceCoversBranches)
{
    const std::vector<KernelCase> cases = make_branch_coverage_cases();
    std::vector<KernelOutputs> host_outputs(cases.size());

    for (int idx = 0; idx < static_cast<int>(cases.size()); ++idx) {
        host_outputs[idx] = host_reference(cases[idx]);
    }

    amrex::Gpu::DeviceVector<KernelCase> d_cases(cases.size());
    amrex::Gpu::DeviceVector<KernelOutputs> d_outputs(cases.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), d_cases.begin());

    launch_helper_kernel(static_cast<int>(cases.size()), d_cases.data(), d_outputs.data());

    std::vector<KernelOutputs> device_outputs(cases.size());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_outputs.begin(), d_outputs.end(), device_outputs.begin());

    for (int idx = 0; idx < static_cast<int>(cases.size()); ++idx) {
        SCOPED_TRACE(case_label(cases[idx]));

        EXPECT_NEAR(device_outputs[idx].sat_cond, host_outputs[idx].sat_cond,
                    formula_abs_tol(host_outputs[idx].sat_cond));
        EXPECT_NEAR(device_outputs[idx].sat_evap, host_outputs[idx].sat_evap,
                    formula_abs_tol(host_outputs[idx].sat_evap));
        EXPECT_NEAR(device_outputs[idx].source_cond, host_outputs[idx].source_cond,
                    formula_abs_tol(host_outputs[idx].source_cond));
        EXPECT_NEAR(device_outputs[idx].source_evap, host_outputs[idx].source_evap,
                    formula_abs_tol(host_outputs[idx].source_evap));
        EXPECT_NEAR(device_outputs[idx].source_cloud_to_rain, host_outputs[idx].source_cloud_to_rain,
                    backend_math_abs_tol(host_outputs[idx].source_cloud_to_rain));
        EXPECT_NEAR(device_outputs[idx].source_rain_to_vapor, host_outputs[idx].source_rain_to_vapor,
                    backend_math_abs_tol(host_outputs[idx].source_rain_to_vapor));
        EXPECT_NEAR(device_outputs[idx].velocity, host_outputs[idx].velocity,
                    backend_math_abs_tol(host_outputs[idx].velocity));
        EXPECT_NEAR(device_outputs[idx].flux, host_outputs[idx].flux,
                    backend_math_abs_tol(host_outputs[idx].flux));
        EXPECT_NEAR(device_outputs[idx].face_rho, host_outputs[idx].face_rho,
                    formula_abs_tol(host_outputs[idx].face_rho));
        EXPECT_NEAR(device_outputs[idx].face_qp, host_outputs[idx].face_qp,
                    formula_abs_tol(host_outputs[idx].face_qp));
        EXPECT_NEAR(device_outputs[idx].sedimentation, host_outputs[idx].sedimentation,
                    formula_abs_tol(host_outputs[idx].sedimentation));
        EXPECT_EQ(device_outputs[idx].is_small, host_outputs[idx].is_small);
        EXPECT_EQ(device_outputs[idx].substeps, host_outputs[idx].substeps);
    }
}