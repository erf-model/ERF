#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TwoStream column-physics regression checker.
//
// Reads a plotfile written after at least one slow step of a TwoStream
// radiation run and verifies the vertical structure of the per-level
// heating rates that the two-stream sweep wrote into qsrc_sw / qsrc_lw:
//
//   1. both fields exist and are finite everywhere;
//   2. shortwave heating is non-negative in every cell (the direct beam
//      only loses energy on the way down) and strictly positive somewhere;
//   3. shortwave heating in the top layer exceeds that in the bottom layer:
//      the beam is strongest at the top and the air densest at the bottom,
//      so a sweep that entered the column at k = 0 would invert this;
//   4. the top layer is the strongest longwave cooling layer in the column
//      (cooling to space) and it cools, i.e. qsrc_lw < 0 there;
//   5. the column-integrated longwave tendency is a net cooling.
//
// Checks 3 and 4 fail if the sweep runs upside down; check 4 also fails if
// the LW heating-rate sign is flipped. The domain is expected to be a
// single box in the vertical (amrex.max_grid_size_z >= n_cell_z).

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

int fail (const std::string& message)
{
    std::cerr << "TwoStreamRadiationCheck: FAIL: " << message << "\n";
    return 1;
}

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

// Horizontal mean of a field at each k of the level-0 domain.
std::vector<Real> horizontal_mean_profile (const MultiFab& mf, const amrex::Box& domain)
{
    const int kmin = domain.smallEnd(2);
    const int kmax = domain.bigEnd(2);
    const int nz = kmax - kmin + 1;
    std::vector<Real> sum(nz, Real(0.0));
    std::vector<Real> count(nz, Real(0.0));

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto& arr = mf.const_array(mfi);
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    sum[k - kmin] += arr(i, j, k);
                    count[k - kmin] += Real(1.0);
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(sum.data(), nz);
    amrex::ParallelDescriptor::ReduceRealSum(count.data(), nz);
    for (int k = 0; k < nz; ++k) {
        sum[k] = (count[k] > Real(0.0)) ? sum[k] / count[k] : Real(0.0);
    }
    return sum;
}

bool all_finite (const MultiFab& mf)
{
    return !mf.contains_nan(0, 1, 0) && !mf.contains_inf(0, 1, 0);
}

int run_checks (const std::string& plotfile_path)
{
    PlotFileData plotfile(plotfile_path);
    if (!has_variable(plotfile, "qsrc_sw")) {
        return fail("plotfile has no qsrc_sw (TwoStream heating not written)");
    }
    if (!has_variable(plotfile, "qsrc_lw")) {
        return fail("plotfile has no qsrc_lw (TwoStream heating not written)");
    }

    MultiFab qsrc_sw = plotfile.get(0, "qsrc_sw");
    MultiFab qsrc_lw = plotfile.get(0, "qsrc_lw");
    if (!all_finite(qsrc_sw)) { return fail("qsrc_sw contains NaN or Inf"); }
    if (!all_finite(qsrc_lw)) { return fail("qsrc_lw contains NaN or Inf"); }

    const amrex::Box domain = plotfile.probDomain(0);
    const std::vector<Real> sw = horizontal_mean_profile(qsrc_sw, domain);
    const std::vector<Real> lw = horizontal_mean_profile(qsrc_lw, domain);
    const int nz = static_cast<int>(sw.size());
    if (nz < 3) { return fail("need at least 3 vertical levels"); }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cout << "k  <qsrc_sw> [K/s]  <qsrc_lw> [K/s]\n";
        for (int k = 0; k < nz; ++k) {
            std::cout << k << "  " << sw[k] << "  " << lw[k] << "\n";
        }
    }

    // 2. SW heating non-negative everywhere, positive somewhere.
    const Real sw_min = qsrc_sw.min(0);
    const Real sw_max = qsrc_sw.max(0);
    if (sw_min < Real(0.0)) {
        return fail("qsrc_sw has a negative value: " + std::to_string(sw_min));
    }
    if (!(sw_max > Real(0.0))) {
        return fail("qsrc_sw is zero everywhere; radiation heating was not applied");
    }

    // 3. Orientation: SW heating at the top layer exceeds the bottom layer.
    if (!(sw[nz - 1] > sw[0])) {
        return fail("qsrc_sw at the top layer (" + std::to_string(sw[nz - 1]) +
                    ") does not exceed the bottom layer (" + std::to_string(sw[0]) +
                    "); the sweep appears to run upside down");
    }

    // 4. Top layer is the strongest LW cooling layer and it cools.
    const int k_min_lw = static_cast<int>(std::min_element(lw.begin(), lw.end()) - lw.begin());
    if (!(lw[nz - 1] < Real(0.0))) {
        return fail("qsrc_lw at the top layer is not negative (" + std::to_string(lw[nz - 1]) +
                    "); no cooling to space, check orientation and heating sign");
    }
    if (k_min_lw != nz - 1) {
        return fail("strongest LW cooling is at k = " + std::to_string(k_min_lw) +
                    " rather than the top layer k = " + std::to_string(nz - 1));
    }

    // 5. Column-integrated LW tendency is a net cooling.
    Real lw_sum = Real(0.0);
    for (Real v : lw) { lw_sum += v; }
    if (!(lw_sum < Real(0.0))) {
        return fail("column-integrated qsrc_lw is not a net cooling: " + std::to_string(lw_sum));
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cout << "TwoStreamRadiationCheck: PASS\n";
    }
    return 0;
}

} // namespace

int main (int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <plotfile>\n";
        return 2;
    }
    amrex::Initialize(argc, argv);
    int result = 0;
    {
        result = run_checks(argv[1]);
    }
    amrex::Finalize();
    return result;
}
