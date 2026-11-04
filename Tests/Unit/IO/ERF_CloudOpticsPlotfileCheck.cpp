// Independent checks of native Plot2D/3D data and the station text export.
#include <AMReX.H>
#include <AMReX_Math.H>
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
const char* paths[] = {"integrated_qc", "integrated_qi", "integrated_qr",
                       "integrated_qs", "integrated_qg"};
// Public native conserved names for this explicitly Morrison test fixture.
const char* masses[] = {"rhoQ2", "rhoQ3", "rhoQ4", "rhoQ5", "rhoQ6"};
const Real kext[] = {Real(150.), Real(65.4), Real(3.), Real(30.), Real(7.5)};

void
require (bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool
has (const PlotFileData& plot, const char* name)
{
    return std::find(plot.varNames().begin(), plot.varNames().end(), name) !=
           plot.varNames().end();
}

void
compare (MultiFab& data, const MultiFab& reference, const std::string& name)
{
    MultiFab aligned(data.boxArray(), data.DistributionMap(), 1, 0);
    aligned.setVal(std::numeric_limits<Real>::quiet_NaN());
    aligned.ParallelCopy(reference);
    MultiFab::Subtract(data, aligned, 0, 0, 1, 0);
    const Real error = data.norm0(0);
    const Real tol = 512. * std::numeric_limits<Real>::epsilon() *
                     amrex::max(Real(1.), reference.norm0(0));
    Print() << name << " max_abs_difference=" << error << '\n';
    require(data.is_finite() && Math::isfinite(error) && error <= tol,
            name + " differs");
}

void
compare_plots (PlotFileData& plot, PlotFileData& reference)
{
    require(plot.finestLevel() == reference.finestLevel(),
            "plot level counts differ");
    int count = 0;
    for (int lev = 0; lev <= plot.finestLevel(); ++lev) {
        for (const auto& name : plot.varNames()) {
            if (!has(reference, name.c_str())) {
                continue;
            }
            auto data = plot.get(lev, name);
            auto ref = reference.get(lev, name);
            compare(data, ref, "level " + std::to_string(lev) + " " + name);
            ++count;
        }
    }
    require(count > 0, "no common fields compared");
}

void
check_formula (PlotFileData& plot, const Vector<MultiFab>* native_tau = nullptr)
{
    for (int lev = 0; lev <= plot.finestLevel(); ++lev) {
        auto tau = plot.get(lev, "cloud_optical_depth_sw");
        auto alb = plot.get(lev, "cloud_albedo_est");
        MultiFab expected(tau.boxArray(), tau.DistributionMap(), 1, 0);
        expected.setVal(0.);
        if (native_tau) {
            expected.ParallelCopy((*native_tau)[lev]);
        } else {
            for (int n = 0; n < 5; ++n) {
                if (!has(plot, paths[n])) {
                    continue;
                }
                auto path = plot.get(lev, paths[n]);
                MultiFab::Saxpy(expected, kext[n], path, 0, 0, 1, 0);
            }
        }
        Print() << "level " << lev << " tau_range=[" << tau.min(0) << ','
                << tau.max(0) << "] albedo_range=[" << alb.min(0) << ','
                << alb.max(0) << "]\n";
        require(tau.is_finite() && alb.is_finite() && tau.min(0) >= 0. &&
                    tau.max(0) > 0. && alb.min(0) >= 0. && alb.max(0) <= 1.,
                "invalid optical range");
        auto tau_error = plot.get(lev, "cloud_optical_depth_sw");
        compare(tau_error, expected,
                "level " + std::to_string(lev) +
                    (native_tau ? " clipped native 3D mass"
                                : " weighted water paths"));
        for (MFIter mfi(expected); mfi.isValid(); ++mfi) {
            const auto t = tau.const_array(mfi);
            const auto out = expected.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j,
                                                             int k) noexcept {
                out(i, j, k) = t(i, j, k) / (t(i, j, k) + Real(7.7));
            });
        }
        compare(alb, expected,
                "level " + std::to_string(lev) + " reflectance formula");
    }
}

// Independently integrate native 3D conserved values using a serial k loop
// per (i,j), then combine partial boxes with AMReX communication. Neither the
// production water-path helper nor the cloud-optics helper is called here.
MultiFab
mass_path (const MultiFab& mass, Real dz, bool clip)
{
    BoxList bl = mass.boxArray().boxList();
    for (auto& b : bl) {
        b.setRange(2, 0);
    }
    BoxArray slab(std::move(bl));
    MultiFab partial(slab, mass.DistributionMap(), 1, 0);
    for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
        const auto in = mass.const_array(mfi);
        const auto out = partial.array(mfi);
        const int lo = mfi.validbox().smallEnd(2),
                  hi = mfi.validbox().bigEnd(2);
        ParallelFor(slab[mfi.index()], [=] AMREX_GPU_DEVICE(int i, int j,
                                                            int k) noexcept {
            Real sum = 0.;
            for (int z = lo; z <= hi; ++z) {
                sum += clip ? amrex::max(Real(0.), in(i, j, z)) : in(i, j, z);
            }
            out(i, j, k) = sum * dz;
        });
    }
    BoxArray unique(slab.boxList());
    unique.removeOverlap();
    MultiFab sum(unique, DistributionMapping(unique), 1, 0);
    sum.setVal(0.);
    sum.ParallelAdd(partial);
    return sum;
}

void
check_native (PlotFileData& plot, PlotFileData& volume, bool partial)
{
    require(plot.finestLevel() == volume.finestLevel(), "2D/3D levels differ");
    if (partial) {
        require(plot.finestLevel() == 1,
                "partial fixture needs two AMR levels");
    }
    Vector<MultiFab> native_tau(plot.finestLevel() + 1);
    for (int lev = 0; lev <= plot.finestLevel(); ++lev) {
        native_tau[lev].define(plot.boxArray(lev), plot.DistributionMap(lev), 1,
                               0);
        native_tau[lev].setVal(0.);
        for (int n = 0; n < 5; ++n) {
            auto mass = volume.get(lev, masses[n]);
            Print() << "level " << lev << ' ' << masses[n]
                    << " native minimum=" << mass.min(0) << '\n';
            // Existing integrated_q* fields preserve signed mass. Cloud optics
            // instead clips negative numerical noise, as its public definition
            // requires.
            auto expected = mass_path(mass, volume.cellSize(lev)[2], false);
            auto clipped = mass_path(mass, volume.cellSize(lev)[2], true);
            MultiFab aligned(plot.boxArray(lev), plot.DistributionMap(lev), 1,
                             0);
            aligned.ParallelCopy(clipped);
            MultiFab::Saxpy(native_tau[lev], kext[n], aligned, 0, 0, 1, 0);
            auto actual = plot.get(lev, paths[n]);
            if (!partial) {
                require(actual.max(0) > 0.,
                        std::string(paths[n]) + " was not exercised");
            }
            Print() << "level " << lev << ' ' << paths[n] << " range=["
                    << actual.min(0) << ',' << actual.max(0) << "]\n";
            compare(actual, expected,
                    "level " + std::to_string(lev) + " " + paths[n] +
                        " from native 3D");
        }
    }
    check_formula(plot, &native_tau);
    if (partial) {
        const Box coverage = volume.boxArray(1).minimalBox();
        const Box domain = volume.probDomain(1);
        require(coverage.smallEnd(2) == domain.smallEnd(2) &&
                    coverage.bigEnd(2) < domain.bigEnd(2),
                "fine level must reach ground and stop below the domain top");
        auto coarse = plot.get(0, "cloud_optical_depth_sw");
        auto fine = plot.get(1, "cloud_optical_depth_sw");
        require(coarse.min(0) > fine.max(0),
                "fixture must show omitted upper-cloud mass");
        const Real top =
            volume.probLo()[2] + (coverage.bigEnd(2) - domain.smallEnd(2) + 1) *
                                     volume.cellSize(1)[2];
        Print() << "PARTIAL-LEVEL PROXIES, NOT COMPOSITE ALBEDO: fine_top="
                << top << " domain_top=" << volume.probHi()[2]
                << " coarse_tau_min=" << coarse.min(0)
                << " partial_fine_tau_max=" << fine.max(0) << '\n';
    }
}

Real
point_value (const MultiFab& field, int i, int j)
{
    BoxArray ba(Box(IntVect(i, j, 0), IntVect(i, j, 0)));
    MultiFab point(ba, DistributionMapping(ba), 1, 0);
    point.setVal(std::numeric_limits<Real>::quiet_NaN());
    point.ParallelCopy(field);
    require(point.is_finite(), "station reference point is not covered");
    return point.sum(0);
}

void
check_station (PlotFileData& plot, const char* file, Real x, Real y)
{
    // These fixtures cover the station horizontally on every level and have
    // ground sources, so the existing 2D-only station rule selects finest.
    const int lev = plot.finestLevel();
    const auto dx = plot.cellSize(lev), lo = plot.probLo();
    const Box domain = plot.probDomain(lev);
    const Real fx = (x - lo[0]) / dx[0] + domain.smallEnd(0) - Real(.5);
    const Real fy = (y - lo[1]) / dx[1] + domain.smallEnd(1) - Real(.5);
    const int i = static_cast<int>(std::floor(fx)),
              j = static_cast<int>(std::floor(fy));
    const Real wx = fx - i, wy = fy - j;
    require(i >= domain.smallEnd(0) && i + 1 <= domain.bigEnd(0) &&
                j >= domain.smallEnd(1) && j + 1 <= domain.bigEnd(1),
            "fixture station must be interior");
    auto expected = [&] (const char* name) {
        auto f = plot.get(lev, name);
        return (1 - wy) * ((1 - wx) * point_value(f, i, j) +
                           wx * point_value(f, i + 1, j)) +
               wy * ((1 - wx) * point_value(f, i, j + 1) +
                     wx * point_value(f, i + 1, j + 1));
    };
    const Real tau = expected("cloud_optical_depth_sw"),
               alb = expected("cloud_albedo_est");
    std::ifstream in(file);
    require(in.good(), "station file is unreadable");
    bool found = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        std::istringstream row(line);
        Real time, observed_alb, observed_tau;
        require(bool(row >> time >> observed_alb >> observed_tau),
                "invalid station row");
        if (std::abs(time - plot.time()) > 1.e-9) {
            continue;
        }
        // Station output has 10 significant decimal digits; allow its text
        // rounding plus a precision-scaled arithmetic tolerance.
        const Real tol =
            Real(2.e-9) + 512 * std::numeric_limits<Real>::epsilon();
        require(std::isfinite(observed_tau) && std::isfinite(observed_alb) &&
                    std::abs(observed_tau - tau) <=
                        tol * amrex::max(Real(1.), std::abs(tau)) &&
                    std::abs(observed_alb - alb) <= tol,
                "station differs from Plot2D column");
        Print() << "station time=" << time << " level=" << lev
                << " tau=" << observed_tau << " Plot2D_tau=" << tau
                << " albedo=" << observed_alb << " Plot2D_albedo=" << alb
                << '\n';
        found = true;
    }
    require(found, "no station row matches plotfile time");
}
} // namespace

int
main (int argc, char** argv)
{
    if (argc < 2) {
        return 2;
    }
    Initialize(argc, argv, false);
    int result = 0;
    try {
        if (argc == 6 && std::string(argv[1]) == "--station") {
            PlotFileData plot(argv[2]);
            check_station(plot, argv[3], Real(std::stod(argv[4])),
                          Real(std::stod(argv[5])));
        } else if (argc == 4 && (std::string(argv[1]) == "--mixed" ||
                                 std::string(argv[1]) == "--partial")) {
            PlotFileData plot(argv[2]), volume(argv[3]);
            check_native(plot, volume, std::string(argv[1]) == "--partial");
        } else if (argc == 2) {
            PlotFileData plot(argv[1]);
            check_formula(plot);
        } else if (argc == 3) {
            PlotFileData plot(argv[1]), reference(argv[2]);
            compare_plots(plot, reference);
        } else {
            require(false, "invalid checker arguments");
        }
    } catch (const std::exception& error) {
        Print() << "Cloud optics native check failed: " << error.what() << '\n';
        result = 1;
    }
    Finalize();
    return result;
}
