#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_VisMF.H>
#include <ERF_IndexDefines.H>

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void require (const bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void compare_component (const amrex::MultiFab& lhs,
                        const amrex::MultiFab& rhs,
                        const int component,
                        const std::string& label)
{
    using namespace amrex;
    require(lhs.boxArray() == rhs.boxArray(), label + ": checkpoint BoxArrays differ");
    require(lhs.DistributionMap() == rhs.DistributionMap(),
            label + ": checkpoint DistributionMaps differ");
    require(lhs.nComp() == rhs.nComp(), label + ": checkpoint component counts differ");
    require(component >= 0 && component < lhs.nComp(),
            label + ": requested component is outside the checkpoint state");

    ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    const Real eps = std::numeric_limits<Real>::epsilon();
    const Real nonfinite_sentinel = std::numeric_limits<Real>::max();
    for (MFIter mfi(lhs); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto lhs_arr = lhs.const_array(mfi);
        const auto rhs_arr = rhs.const_array(mfi);
        const int comp = component;
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
                -> GpuTuple<Real, Real> {
                const Real a = lhs_arr(i, j, k, comp);
                const Real b = rhs_arr(i, j, k, comp);
                if (!Math::isfinite(a) || !Math::isfinite(b)) {
                    return {nonfinite_sentinel, nonfinite_sentinel};
                }
                const Real error = Math::abs(a - b);
                const Real scale = max(Real(1.0), max(Math::abs(a), Math::abs(b)));
                const Real tolerance = Real(512.0) * eps * scale;
                return {error, error / tolerance};
            });
    }

    const auto reduced = reduce_data.value(reduce_op);
    Real max_error = get<0>(reduced);
    Real max_tolerance_ratio = get<1>(reduced);
    ParallelDescriptor::ReduceRealMax(&max_error, 1);
    ParallelDescriptor::ReduceRealMax(&max_tolerance_ratio, 1);
    Print() << label << ": max_abs_difference=" << max_error
            << " max_error_over_tolerance=" << max_tolerance_ratio << '\n';
    require(max_tolerance_ratio <= Real(1.0),
            label + ": full valid-cell data exceed 512-epsilon consistency tolerance");
}

void compare_all_components (const amrex::MultiFab& lhs,
                             const amrex::MultiFab& rhs,
                             const std::string& label)
{
    for (int component = 0; component < lhs.nComp(); ++component) {
        compare_component(lhs, rhs, component,
                          label + " component " + std::to_string(component));
    }
}

//
// The restarted checkpoint as it would read on the continuous one's grids. A restart that
// re-made the level-0 grids (issue 4225) writes its checkpoint on the new ones, so with
// allow_different_grids the valid cells are copied across first; both must cover the same
// region. Without it the grids must match, as they always did.
//
amrex::MultiFab on_grids_of (const amrex::MultiFab& lhs,
                             amrex::MultiFab&& rhs,
                             const bool allow_different_grids,
                             const std::string& label)
{
    if (lhs.boxArray() == rhs.boxArray() || !allow_different_grids) {
        return std::move(rhs);
    }
    require(lhs.boxArray().minimalBox() == rhs.boxArray().minimalBox() &&
            lhs.boxArray().numPts() == rhs.boxArray().numPts(),
            label + ": checkpoints cover different regions");
    require(lhs.nComp() == rhs.nComp(), label + ": checkpoint component counts differ");
    amrex::MultiFab aligned(lhs.boxArray(), lhs.DistributionMap(), lhs.nComp(), 0);
    aligned.ParallelCopy(rhs, 0, 0, rhs.nComp());
    return aligned;
}

void compare_checkpoint (const std::string& lhs_path, const std::string& rhs_path,
                         const bool allow_different_grids)
{
    amrex::MultiFab lhs_spectrum;
    amrex::MultiFab rhs_spectrum;
    amrex::VisMF::Read(lhs_spectrum, lhs_path + "/Level_0/SBMSpectrum");
    amrex::VisMF::Read(rhs_spectrum, rhs_path + "/Level_0/SBMSpectrum");
    compare_all_components(lhs_spectrum,
                           on_grids_of(lhs_spectrum, std::move(rhs_spectrum),
                                       allow_different_grids, "SBMSpectrum"),
                           "SBMSpectrum");

    amrex::MultiFab lhs_cell;
    amrex::MultiFab rhs_cell;
    amrex::VisMF::Read(lhs_cell, lhs_path + "/Level_0/Cell");
    amrex::VisMF::Read(rhs_cell, rhs_path + "/Level_0/Cell");
    const amrex::MultiFab rhs_cell_aligned =
        on_grids_of(lhs_cell, std::move(rhs_cell), allow_different_grids, "Cell");
    compare_component(lhs_cell, rhs_cell_aligned, RhoQ2_comp, "projected qc");
    compare_component(lhs_cell, rhs_cell_aligned, RhoQ3_comp, "projected qr");
}

} // namespace

int main (int argc, char* argv[])
{
    // The two arguments are checkpoint directories, not ParmParse input. With
    // the command line parsed, AMReX takes argv[1] (no "=" in it) for an inputs
    // file and opens it: a directory opens, its tellg() is LLONG_MAX, and the
    // read buffer resize throws std::length_error before main's try block, so
    // the comparison dies with "terminate called ... vector::_M_default_append"
    // instead of running. An empty ParmParse is still built.
    amrex::Initialize(argc, argv, false);
    int result = 0;
    try {
        const bool allow_different_grids =
            (argc == 4 && std::string(argv[3]) == "--allow-different-grids");
        require(argc == 3 || allow_different_grids,
                "usage: erf_sbm_checkpoint_compare <continuous-checkpoint> <restarted-checkpoint>"
                " [--allow-different-grids]");
        compare_checkpoint(argv[1], argv[2], allow_different_grids);
    } catch (const std::exception& error) {
        amrex::Print() << "SBM checkpoint comparison failed: " << error.what() << '\n';
        result = 1;
    }
    amrex::Finalize();
    return result;
}
