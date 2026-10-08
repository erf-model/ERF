#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_VisMF.H>
#include <ERF_IndexDefines.H>

#include <limits>
#include <stdexcept>
#include <string>

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

void compare_checkpoint (const std::string& lhs_path, const std::string& rhs_path)
{
    amrex::MultiFab lhs_spectrum;
    amrex::MultiFab rhs_spectrum;
    amrex::VisMF::Read(lhs_spectrum, lhs_path + "/Level_0/SBMSpectrum");
    amrex::VisMF::Read(rhs_spectrum, rhs_path + "/Level_0/SBMSpectrum");
    compare_all_components(lhs_spectrum, rhs_spectrum, "SBMSpectrum");

    amrex::MultiFab lhs_cell;
    amrex::MultiFab rhs_cell;
    amrex::VisMF::Read(lhs_cell, lhs_path + "/Level_0/Cell");
    amrex::VisMF::Read(rhs_cell, rhs_path + "/Level_0/Cell");
    compare_component(lhs_cell, rhs_cell, RhoQ2_comp, "projected qc");
    compare_component(lhs_cell, rhs_cell, RhoQ3_comp, "projected qr");
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
        require(argc == 3,
                "usage: erf_sbm_checkpoint_compare <continuous-checkpoint> <restarted-checkpoint>");
        compare_checkpoint(argv[1], argv[2]);
    } catch (const std::exception& error) {
        amrex::Print() << "SBM checkpoint comparison failed: " << error.what() << '\n';
        result = 1;
    }
    amrex::Finalize();
    return result;
}
