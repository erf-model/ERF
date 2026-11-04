// Inspect actual native Plot2D data, independently of the diagnostic helper.
#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <limits>

using namespace amrex;

int
main (int argc, char** argv)
{
    if (argc != 2 && argc != 3) {
        return 2;
    }
    Initialize(argc, argv, false);
    int result = 0;
    {
        PlotFileData plot(argv[1]);
        if (argc == 3) {
            // Compare every common field, including single-name requests or
            // all existing diagnostics in a base-commit regression plotfile.
            PlotFileData reference(argv[2]);
            int compared = 0;
            for (const auto& name : plot.varNames()) {
                if (std::find(reference.varNames().begin(),
                              reference.varNames().end(),
                              name) == reference.varNames().end()) {
                    continue;
                }
                ++compared;
                auto data = plot.get(0, name);
                auto ref = reference.get(0, name);
                MultiFab aligned(data.boxArray(), data.DistributionMap(), 1, 0);
                aligned.ParallelCopy(ref);
                MultiFab::Subtract(data, aligned, 0, 0, 1, 0);
                const Real error = data.norm0(0);
                const Real tolerance = 512. *
                                       std::numeric_limits<Real>::epsilon() *
                                       amrex::max(Real(1.), ref.norm0(0));
                Print() << name << " max_abs_difference=" << error << '\n';
                if (!data.is_finite() || !Math::isfinite(error) ||
                    error > tolerance) {
                    result = 1;
                }
            }
            if (compared == 0) {
                result = 1;
            }
        } else {
            auto tau = plot.get(0, "cloud_optical_depth_sw");
            auto alb = plot.get(0, "cloud_albedo_est");
            auto qc = plot.get(0, "integrated_qc");
            auto qr = plot.get(0, "integrated_qr");
            MultiFab error(tau.boxArray(), tau.DistributionMap(), 2, 0);
            for (MFIter mfi(error); mfi.isValid(); ++mfi) {
                const auto t = tau.const_array(mfi);
                const auto a = alb.const_array(mfi);
                const auto c = qc.const_array(mfi);
                const auto r = qr.const_array(mfi);
                const auto out = error.array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(
                                                int i, int j, int k) noexcept {
                    const Real expected_tau =
                        150. * c(i, j, k) + 3. * r(i, j, k);
                    const Real expected_albedo =
                        t(i, j, k) / (t(i, j, k) + 7.7);
                    out(i, j, k, 0) =
                        Math::abs(t(i, j, k) - expected_tau) /
                        amrex::max(Real(1.), Math::abs(expected_tau));
                    out(i, j, k, 1) = Math::abs(a(i, j, k) - expected_albedo);
                });
            }
            const Real tau_error = error.norm0(0);
            const Real alb_error = error.norm0(1);
            const Real tolerance = 512. * std::numeric_limits<Real>::epsilon();
            Print() << "tau_range=[" << tau.min(0) << ',' << tau.max(0)
                    << "] albedo_range=[" << alb.min(0) << ',' << alb.max(0)
                    << "] tau_relative_error=" << tau_error
                    << " albedo_absolute_error=" << alb_error << '\n';
            if (!tau.is_finite() || !alb.is_finite() || tau.min(0) < 0. ||
                tau.max(0) <= 0. || alb.min(0) < 0. || alb.max(0) > 1. ||
                tau_error > tolerance || alb_error > tolerance) {
                result = 1;
            }
        }
    }
    Finalize();
    return result;
}
