// Test-only deterministic mass profiles in copied Morrison or Kessler
// checkpoints.
// No model evolution or production initialization is changed by this utility.
#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include "ERF_DataStruct.H"

using namespace amrex;

int
main (int argc, char** argv)
{
    const bool kessler = argc == 3 && std::string(argv[1]) == "--kessler";
    if (!(argc == 2 || kessler)) {
        return 2;
    }
    const char* checkpoint = kessler ? argv[2] : argv[1];
    int amrex_argc = 1;
    char* amrex_argv_storage[] = {argv[0], nullptr};
    char** amrex_argv = amrex_argv_storage;
    Initialize(amrex_argc, amrex_argv, false);
    {
        const auto mi = MoistureComponentIndices::from_moisture_model(
            kessler ? MoistureType::Kessler : MoistureType::Morrison);
        const GpuArray<int, 5> components{{mi.qc, mi.qi, mi.qr, mi.qs, mi.qg}};
        for (int lev = 0;; ++lev) {
            const std::string prefix = std::string(checkpoint) + "/Level_" +
                                       std::to_string(lev) + "/Cell";
            if (!FileExists(prefix + "_H")) {
                break;
            }
            MultiFab state;
            VisMF::Read(state, prefix);
            if (kessler) {
                AMREX_ALWAYS_ASSERT(mi.qc >= 0 && mi.qc < state.nComp());
                AMREX_ALWAYS_ASSERT(mi.qr >= 0 && mi.qr < state.nComp());
            } else {
                for (int comp : components) {
                    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < state.nComp());
                }
            }
            const Box domain = state.boxArray().minimalBox();
            const int middle = domain.smallEnd(2) + domain.length(2) / 2;
            for (MFIter mfi(state); mfi.isValid(); ++mfi) {
                const auto a = state.array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(
                                                int i, int j, int k) noexcept {
                    if (kessler) {
                        // Distinct, positive conserved liquid and rain masses
                        // make both non-default optical overrides measurable.
                        a(i, j, k, components[0]) = Real(1.e-5);
                        a(i, j, k, components[2]) = Real(2.e-5);
                    } else {
                        // All five species are positive and distinct. Most mass
                        // is above the future partial fine level, making its
                        // omission visible.
                        const Real mass = k < middle ? Real(1.e-6) : Real(1.e-5);
                        for (int n = 0; n < 5; ++n) {
                            a(i, j, k, components[n]) = Real(n + 1) * mass;
                        }
                    }
                });
            }
            Gpu::streamSynchronize();
            VisMF::Write(state, prefix);
            if (kessler) {
                Print() << "Synthetic Kessler qc/qr masses written at level "
                        << lev << '\n';
            } else {
                Print() << "Synthetic Morrison mass profile written at level "
                        << lev << ", upper half starts at k=" << middle << '\n';
            }
        }
    }
    Finalize();
    return 0;
}
