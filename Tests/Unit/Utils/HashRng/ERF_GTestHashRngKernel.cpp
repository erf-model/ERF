#include <cstdint>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <gtest/gtest.h>

#include "ERF_GTestHashRngCommon.H"

using namespace hash_rng_test;

namespace {

// Keep the extended device lambda outside GoogleTest's private TestBody so
// CUDA extended-lambda diagnostics do not reject the test.
void launch_hash_rng_cases (amrex::Gpu::DeviceVector<std::uint64_t>& results)
{
    std::uint64_t* result_ptr = results.data();
    amrex::ParallelFor(
        kHostDeviceCaseCount,
        [=] AMREX_GPU_DEVICE (int case_id) noexcept {
            result_ptr[case_id] = hash_rng_case_word(case_id);
        });
}

} // namespace

TEST(HashRngKernel, HostDeviceWordsAgreeExactly)
{
    amrex::Gpu::DeviceVector<std::uint64_t> device_results(kHostDeviceCaseCount);
    amrex::Gpu::HostVector<std::uint64_t> host_results(kHostDeviceCaseCount);

    launch_hash_rng_cases(device_results);
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     device_results.begin(), device_results.end(), host_results.begin());

    for (int case_id = 0; case_id < kHostDeviceCaseCount; ++case_id) {
        const std::uint64_t expected = hash_rng_case_word(case_id);
        // Unlike floating-point physics kernels, this hash uses only integer
        // operations and exactly representable power-of-two mappings. Exact
        // uint64_t equality is therefore the CPU/GPU portability contract.
        if (host_results[case_id] != expected) {
            ADD_FAILURE() << "case_id=" << case_id
                          << " host=" << expected
                          << " device=" << host_results[case_id];
            return;
        }
    }
}
