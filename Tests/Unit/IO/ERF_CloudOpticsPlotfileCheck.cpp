// Writer integration oracle: written water paths, optical fields, and station.
// Independent mass integration and edge cases belong to the CloudOptics gtests.
#include <AMReX.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
using namespace amrex;

namespace {
void require (bool ok, const std::string& message)
{
    if (!ok) { throw std::runtime_error(message); }
}
Real check (Real actual, Real expected, const std::string& field,
            const IntVect& cell, Real allowance = 0.)
{
    const Real tolerance = allowance + 512 * std::numeric_limits<Real>::epsilon() *
                                       std::max(Real(1.), std::abs(expected));
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance) {
        std::ostringstream error;
        error.precision(std::numeric_limits<Real>::max_digits10);
        error << field << " level=0 cell=" << cell << " expected=" << expected
              << " actual=" << actual << " tolerance=" << tolerance;
        require(false, error.str());
    }
    return std::abs(actual - expected);
}

void check_writer (PlotFileData& plot, PlotFileData& other, const char* station)
{
    require(ParallelDescriptor::NProcs() == 1, "checker requires one process");
    require(plot.finestLevel() == 0 && other.finestLevel() == 0,
            "writer fixture must be single-level");
    require(plot.time() > 0. && plot.time() == other.time(),
            "output streams must match at a positive simulation time");
    const char* names[] = {"integrated_qc", "integrated_qr",
                          "cloud_optical_depth_sw", "cloud_albedo_est"};
    Vector<MultiFab> fields;
    fields.reserve(5);
    for (const auto* name : names) { fields.push_back(plot.get(0, name)); }
    auto second = other.get(0, "cloud_albedo_est");
    fields.emplace_back(fields[0].boxArray(), fields[0].DistributionMap(), 1, 0);
    fields.back().setVal(std::numeric_limits<Real>::quiet_NaN());
    fields.back().ParallelCopy(second);

    // Same interior station coordinates as the CMake fixture, deliberately
    // between cell centers to exercise bilinear interpolation.
    const auto dx = plot.cellSize(0), lo = plot.probLo();
    const auto domain = plot.probDomain(0);
    const Real fx = (8250 - lo[0]) / dx[0] + domain.smallEnd(0) - Real(.5);
    const Real fy = (175 - lo[1]) / dx[1] + domain.smallEnd(1) - Real(.5);
    const int si = int(std::floor(fx)), sj = int(std::floor(fy));
    const Real wx = fx - si, wy = fy - sj;
    require(si >= domain.smallEnd(0) && si + 1 <= domain.bigEnd(0) &&
            sj >= domain.smallEnd(1) && sj + 1 <= domain.bigEnd(1),
            "station interpolation stencil must be interior");
    Real station_tau = 0., station_albedo = 0., max_qc = 0., max_qr = 0.;
    Real tau_error = 0., albedo_error = 0., default_delta = 0., max_tau = 0.;
    Real min_path = std::numeric_limits<Real>::max();
    int stencil_cells = 0;
    for (MFIter mfi(fields[0]); mfi.isValid(); ++mfi) {
        // Copy to host explicitly so this small serial reader also works when
        // AMReX uses CUDA/HIP device allocations. No device test lambdas.
        Vector<FArrayBox> host;
        host.reserve(fields.size());
        for (auto& field : fields) {
            const auto& fab = field[mfi];
            host.emplace_back(fab.box(), 1, The_Pinned_Arena());
            Gpu::copy(Gpu::deviceToHost, fab.dataPtr(), fab.dataPtr() + fab.size(),
                      host.back().dataPtr());
        }
        for (BoxIterator bit(mfi.validbox()); bit.ok(); ++bit) {
            const auto& iv = bit();
            const Real qc = host[0](iv), qr = host[1](iv);
            const Real tau = host[2](iv), albedo = host[3](iv);
            const Real path_tol = 512 * std::numeric_limits<Real>::epsilon();
            require(std::isfinite(qc) && std::isfinite(qr) &&
                    qc >= -path_tol && qr >= -path_tol,
                    "signed water paths outside nonnegative fixture regime at " +
                    std::to_string(iv[0]) + "," + std::to_string(iv[1]));
            require(tau >= 0. && albedo >= 0. && albedo <= 1.,
                    "optical range/sentinel failure at " + std::to_string(iv[0]) +
                    "," + std::to_string(iv[1]));
            // Do not clip the written paths to conceal a negative-mass mismatch.
            tau_error = std::max(tau_error, check(tau, 50 * qc + 2 * qr, "tau", iv));
            albedo_error = std::max(albedo_error,
                                   check(albedo, tau / (tau + 3), "albedo", iv));
            check(host[4](iv), albedo, "stream2 albedo", iv);
            max_tau = std::max(max_tau, tau);
            max_qc = std::max(max_qc, qc); max_qr = std::max(max_qr, qr);
            min_path = std::min(min_path, std::min(qc, qr));
            default_delta = std::max(default_delta, std::abs(tau - (150 * qc + 3 * qr)));
            if (iv[0] >= si && iv[0] <= si + 1 && iv[1] >= sj && iv[1] <= sj + 1) {
                const Real weight = (iv[0] == si ? 1 - wx : wx) *
                                    (iv[1] == sj ? 1 - wy : wy);
                station_tau += weight * tau; station_albedo += weight * albedo;
                ++stencil_cells;
            }
        }
    }
    require(max_tau > 0. && max_qc > 0. && default_delta > 512 * std::numeric_limits<Real>::epsilon(),
            "fixture must have nonzero liquid and distinguish overridden coefficients");
    require(stencil_cells == 4, "incomplete station interpolation stencil");
    std::ifstream input(station);
    require(input.good(), "cannot read station output");
    std::string line;
    bool matched = false;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') { continue; }
        Real time, albedo, tau;
        std::istringstream row(line);
        require(bool(row >> time >> albedo >> tau), "malformed station row");
        require(std::isfinite(time), "nonfinite station time");
        if (std::abs(time - plot.time()) > 1.e-9) { continue; }
        // Station files retain 10 significant decimal digits.
        const IntVect cell(si, sj, 0);
        check(tau, station_tau, "station tau", cell,
              Real(2.e-9) * std::max(Real(1.), std::abs(station_tau)));
        check(albedo, station_albedo, "station albedo", cell, Real(2.e-9));
        Print() << "station time=" << time << " tau=" << tau << " expected="
                << station_tau << " albedo=" << albedo << " expected=" << station_albedo << '\n';
        matched = true;
    }
    require(matched, "no station row matches Plot2D time");
    Print() << "writer time=" << plot.time() << " min_water_path=" << min_path
            << " max_qc_contribution=" << 50 * max_qc << " max_qr_contribution=" << 2 * max_qr
            << " max_tau=" << max_tau << " tau_max_error=" << tau_error << " albedo_max_error=" << albedo_error
            << " tau_default_delta=" << default_delta << '\n';
}
} // namespace

int main (int argc, char** argv)
{
    if (argc != 4) { return 2; }
    Initialize(argc, argv, false);
    OutStream().precision(std::numeric_limits<Real>::max_digits10);
    int status = 0;
    try {
        PlotFileData plot(argv[1]), other(argv[2]);
        check_writer(plot, other, argv[3]);
    } catch (const std::exception& error) {
        Print() << "Cloud optics writer check failed: " << error.what() << '\n';
        status = 1;
    }
    Finalize();
    return status;
}
