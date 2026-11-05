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

bool
has_vertical_stack (const BoxArray& boxes)
{
    for (int i = 0; i < boxes.size(); ++i) {
        for (int j = i + 1; j < boxes.size(); ++j) {
            const Box& a = boxes[i];
            const Box& b = boxes[j];
            const bool same_horizontal_footprint =
                a.smallEnd(0) == b.smallEnd(0) &&
                a.bigEnd(0) == b.bigEnd(0) &&
                a.smallEnd(1) == b.smallEnd(1) &&
                a.bigEnd(1) == b.bigEnd(1);
            const bool adjacent_in_z = a.bigEnd(2) + 1 == b.smallEnd(2) ||
                                       b.bigEnd(2) + 1 == a.smallEnd(2);
            if (same_horizontal_footprint && adjacent_in_z) {
                return true;
            }
        }
    }
    return false;
}

void
check_runtime_override (PlotFileData& plot, PlotFileData& volume,
                        Real kext_qc, Real kext_qr, Real scale)
{
    require(plot.finestLevel() == volume.finestLevel(),
            "override Plot2D/3D level counts differ");
    require(has(volume, "rhoQ2") && has(volume, "rhoQ3"),
            "Kessler native qc/rain fields are required");
    bool z_stacked_boxes = false;
    for (int lev = 0; lev <= volume.finestLevel(); ++lev) {
        z_stacked_boxes = z_stacked_boxes ||
                          has_vertical_stack(volume.boxArray(lev));
    }
    require(z_stacked_boxes,
            "override fixture must contain vertically stacked 3D boxes");

    for (int lev = 0; lev <= plot.finestLevel(); ++lev) {
        const auto dz = volume.cellSize(lev)[2];
        auto qc_mass = volume.get(lev, "rhoQ2");
        auto qr_mass = volume.get(lev, "rhoQ3");
        auto qc_path = mass_path(qc_mass, dz, true);
        auto qr_path = mass_path(qr_mass, dz, true);

        const Real qc_contribution = kext_qc * qc_path.norm0(0);
        const Real qr_contribution = kext_qr * qr_path.norm0(0);
        const Real contribution_tolerance =
            512. * std::numeric_limits<Real>::epsilon() *
            amrex::max(Real(1.), amrex::max(qc_contribution, qr_contribution));
        require(qc_contribution > 0. && qr_contribution > 0. &&
                    std::abs(qc_contribution - qr_contribution) >
                        contribution_tolerance,
                "override fixture needs nonzero, distinguishable qc and rain contributions");
        Print() << "level " << lev << " override weighted_paths qc="
                << qc_path.norm0(0) << " qr=" << qr_path.norm0(0)
                << " tau_contributions qc=" << qc_contribution
                << " qr=" << qr_contribution << '\n';

        auto tau = plot.get(lev, "cloud_optical_depth_sw");
        auto albedo = plot.get(lev, "cloud_albedo_est");
        MultiFab qc_aligned(tau.boxArray(), tau.DistributionMap(), 1, 0);
        MultiFab qr_aligned(tau.boxArray(), tau.DistributionMap(), 1, 0);
        qc_aligned.setVal(std::numeric_limits<Real>::quiet_NaN());
        qr_aligned.setVal(std::numeric_limits<Real>::quiet_NaN());
        qc_aligned.ParallelCopy(qc_path);
        qr_aligned.ParallelCopy(qr_path);

        // This oracle integrates native rhoQ fields independently of both
        // fill_cloud_optics() and the production water-path diagnostics.
        MultiFab expected_tau(tau.boxArray(), tau.DistributionMap(), 1, 0);
        expected_tau.setVal(0.);
        MultiFab::Saxpy(expected_tau, kext_qc, qc_aligned, 0, 0, 1, 0);
        MultiFab::Saxpy(expected_tau, kext_qr, qr_aligned, 0, 0, 1, 0);
        MultiFab expected_albedo(tau.boxArray(), tau.DistributionMap(), 1, 0);
        for (MFIter mfi(expected_tau); mfi.isValid(); ++mfi) {
            const auto t = expected_tau.const_array(mfi);
            const auto a = expected_albedo.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j,
                                                             int k) noexcept {
                a(i, j, k) = t(i, j, k) / (t(i, j, k) + scale);
            });
        }

        require(expected_tau.is_finite() && expected_tau.min(0) >= 0. &&
                    expected_tau.max(0) > 0.,
                "independent override tau is not positive and finite");
        compare(tau, expected_tau,
                "level " + std::to_string(lev) +
                    " runtime-override tau from native 3D rhoQ2/rhoQ3");
        compare(albedo, expected_albedo,
                "level " + std::to_string(lev) +
                    " runtime-override reflectance from native 3D rhoQ2/rhoQ3");

        MultiFab default_tau(tau.boxArray(), tau.DistributionMap(), 1, 0);
        default_tau.setVal(0.);
        MultiFab::Saxpy(default_tau, kext[0], qc_aligned, 0, 0, 1, 0);
        MultiFab::Saxpy(default_tau, kext[2], qr_aligned, 0, 0, 1, 0);
        MultiFab default_albedo(tau.boxArray(), tau.DistributionMap(), 1, 0);
        for (MFIter mfi(default_tau); mfi.isValid(); ++mfi) {
            const auto t = default_tau.const_array(mfi);
            const auto a = default_albedo.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j,
                                                             int k) noexcept {
                a(i, j, k) = t(i, j, k) / (t(i, j, k) + Real(7.7));
            });
        }
        MultiFab tau_difference(tau.boxArray(), tau.DistributionMap(), 1, 0);
        tau_difference.ParallelCopy(default_tau);
        MultiFab::Subtract(tau_difference, expected_tau, 0, 0, 1, 0);
        MultiFab albedo_difference(tau.boxArray(), tau.DistributionMap(), 1, 0);
        albedo_difference.ParallelCopy(default_albedo);
        MultiFab::Subtract(albedo_difference, expected_albedo, 0, 0, 1, 0);
        const Real tau_tolerance =
            512. * std::numeric_limits<Real>::epsilon() *
            amrex::max(Real(1.), default_tau.norm0(0));
        const Real albedo_tolerance =
            512. * std::numeric_limits<Real>::epsilon() *
            amrex::max(Real(1.), default_albedo.norm0(0));
        require(tau_difference.norm0(0) > tau_tolerance &&
                    albedo_difference.norm0(0) > albedo_tolerance,
                "override values must differ from the production defaults");
        Print() << "level " << lev << " override-default deltas tau="
                << tau_difference.norm0(0) << " albedo="
                << albedo_difference.norm0(0) << '\n';
    }
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
    // Keep the native integration errors useful in reproducible test logs.
    OutStream().precision(std::numeric_limits<Real>::max_digits10);
    int result = 0;
    try {
        if (argc == 6 && std::string(argv[1]) == "--station") {
            PlotFileData plot(argv[2]);
            check_station(plot, argv[3], Real(std::stod(argv[4])),
                          Real(std::stod(argv[5])));
        } else if (argc == 7 && std::string(argv[1]) == "--override") {
            PlotFileData plot(argv[2]), volume(argv[3]);
            check_runtime_override(plot, volume, Real(std::stod(argv[4])),
                                   Real(std::stod(argv[5])),
                                   Real(std::stod(argv[6])));
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
