// Test-only deterministic mass profiles in a copied Morrison checkpoint.
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
    if (argc != 2) {
        return 2;
    }
    Initialize(argc, argv, false);
    {
        const auto mi = MoistureComponentIndices::from_moisture_model(
            MoistureType::Morrison);
        const GpuArray<int, 5> components{{mi.qc, mi.qi, mi.qr, mi.qs, mi.qg}};
        for (int lev = 0;; ++lev) {
            const std::string prefix = std::string(argv[1]) + "/Level_" +
                                       std::to_string(lev) + "/Cell";
            if (!FileExists(prefix + "_H")) {
                break;
            }
            MultiFab state;
            VisMF::Read(state, prefix);
            for (int comp : components) {
                AMREX_ALWAYS_ASSERT(comp >= 0 && comp < state.nComp());
            }
            const Box domain = state.boxArray().minimalBox();
            const int middle = domain.smallEnd(2) + domain.length(2) / 2;
            for (MFIter mfi(state); mfi.isValid(); ++mfi) {
                const auto a = state.array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(
                                                int i, int j, int k) noexcept {
                    // All five species are positive and distinct. Most mass is
                    // above the future partial fine level, making its omission
                    // visible.
                    const Real mass = k < middle ? Real(1.e-6) : Real(1.e-5);
                    for (int n = 0; n < 5; ++n) {
                        a(i, j, k, components[n]) = Real(n + 1) * mass;
                    }
                });
            }
            Gpu::streamSynchronize();
            VisMF::Write(state, prefix);
            Print() << "Synthetic Morrison mass profile written at level "
                    << lev << ", upper half starts at k=" << middle << '\n';
        }
    }
    Finalize();
    return 0;
}
