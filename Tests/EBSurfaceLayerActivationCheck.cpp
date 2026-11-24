// The physical zlo BC must not own an elevated EB wall's MOST state or drag.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_VisMF.H>
#include <ERF_Constants.H>

#include <cmath>
#include <limits>
#include <string>

using namespace amrex;

namespace {

void require (bool condition, const std::string& message)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(condition, message);
}

void check_state (const std::string& checkpoint, int finest, bool initial)
{
    for (int lev = 0; lev <= finest; ++lev) {
        const std::string prefix = checkpoint + "/Level_" + std::to_string(lev) + "/";
        MultiFab z0, ustar;
        VisMF::Read(z0, prefix + "Z0");
        VisMF::Read(ustar, prefix + "Ustar");
        require(z0.is_finite() && ustar.is_finite(), "Non-finite EB MOST state");
        require(std::abs(z0.min(0) - Real(0.1)) < Real(1.e-7) &&
                std::abs(z0.max(0) - Real(0.1)) < Real(1.e-7),
                "EB roughness was not initialized on the level");
        require(z0.nGrowVect()[2] == 1, "EB MOST must retain its 3D layout");

        // The deck's plane at z = 3.25 cuts k = 1 at dx = 2 and k = 3 at dx = 1.
        const int kcut = lev == 0 ? 1 : 3;
        ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpSum> op;
        ReduceData<Real, Real, Long> data(op);
        using Tuple = typename decltype(data)::Type;
        for (MFIter mfi(ustar); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            if (bx.smallEnd(2) > kcut || bx.bigEnd(2) < kcut) { continue; }
            bx.setRange(2, kcut);
            const auto src = ustar.const_array(mfi);
            op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Tuple {
                return {src(i,j,k), src(i,j,k), Long(1)};
            });
        }
        auto result = data.value();
        Real lo = get<0>(result), hi = get<1>(result);
        Long count = get<2>(result);
        ParallelDescriptor::ReduceRealMin(lo);
        ParallelDescriptor::ReduceRealMax(hi);
        ParallelDescriptor::ReduceLongSum(count);
        Print() << "Level " << lev << " cut-cell u*: min=" << lo
                << " max=" << hi << " count=" << count << "\n";
        require(count > 0 && lo > 0 && hi < Real(2),
                "The elevated EB cut cells did not receive finite positive friction velocity");
        if (initial) {
            const Real expected = KAPPA * Real(10) / std::log(Real(10) / Real(0.1));
            const Real tol = Real(256) * std::numeric_limits<Real>::epsilon();
            require(std::abs(lo - expected) < tol && std::abs(hi - expected) < tol,
                    "Initial EB friction velocity does not match the neutral MOST relation");
        }
    }
}

void compare_plotfiles (const std::string& first, const std::string& second)
{
    PlotFileData a(first), b(second);
    require(a.finestLevel() == b.finestLevel() && a.varNames() == b.varNames(),
            "SurfaceLayer parity outputs have different inventories");
    const Real tol = Real(2048) * std::numeric_limits<Real>::epsilon();
    for (int lev = 0; lev <= a.finestLevel(); ++lev) {
        for (const auto& name : a.varNames()) {
            MultiFab x = a.get(lev, name), y = b.get(lev, name);
            require(x.is_finite() && y.is_finite(), "Non-finite output for " + name);
            const Real scale = amrex::max(Real(1), x.norm0());
            MultiFab aligned(x.boxArray(), x.DistributionMap(), 1, 0);
            aligned.ParallelCopy(y);
            MultiFab::Subtract(x, aligned, 0, 0, 1, 0);
            require(x.norm0() <= tol * scale, "SurfaceLayer parity failed for " + name);
        }
    }
}

void check_drag (const std::string& active, const std::string& inactive)
{
    PlotFileData a(active), b(inactive);
    MultiFab u = a.get(0, "x_velocity"), uoff = b.get(0, "x_velocity");
    require(u.is_finite() && uoff.is_finite(), "Non-finite velocities in EB drag check");
    MultiFab::Subtract(uoff, u, 0, 0, 1, 0);
    // Covered cells are zero in both outputs. A positive difference in a fluid
    // cell establishes applied drag, beyond merely constructing a MOST object.
    require(uoff.max(0) > Real(1.e-5), "EB SurfaceLayer did not reduce tangential wind");
    Print() << "Maximum velocity reduction from EB MOST: " << uoff.max(0) << "\n";
}

void check_roughness (const std::string& checkpoint, int finest,
                      bool initial, bool interpolated)
{
    const Real tol = Real(512) * std::numeric_limits<Real>::epsilon();
    for (int lev = 0; lev <= finest; ++lev) {
        const std::string prefix = checkpoint + "/Level_" + std::to_string(lev) + "/";
        MultiFab z0, ustar;
        VisMF::Read(z0, prefix + "Z0");
        VisMF::Read(ustar, prefix + "Ustar");
        require(z0.is_finite(0, 1, z0.nGrowVect()) && ustar.is_finite(),
                "Non-finite varying-roughness MOST state");
        require(z0.nGrowVect()[2] == 1, "EB roughness must have vertical ghosts");
        // Fine-file transition is intentionally shifted by one fine cell to
        // distinguish reading that file from interpolating the coarse map.
        const int nx = 16 << lev;
        const int split = lev == 0 ? 8 : 17;
        const int kcut = lev == 0 ? 1 : 3;
        const bool blended = interpolated && lev > 0;
        ReduceOps<ReduceOpMax, ReduceOpMin, ReduceOpMax> op;
        ReduceData<Real, Real, Real> data(op);
        using Tuple = typename decltype(data)::Type;
        for (MFIter mfi(z0); mfi.isValid(); ++mfi) {
            // Include vertical ghosts and every vertically split FAB. Horizontal
            // ghosts may be periodic or clamped depending on restart/regrid.
            Box bx = amrex::grow(mfi.validbox(), IntVect(0, 0, 1));
            const int ka = mfi.validbox().smallEnd(2);
            const auto a = z0.const_array(mfi);
            op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Tuple {
                const Real value = a(i,j,k);
                Real error = amrex::Math::abs(value - a(i,j,ka));
                if (!blended || i < nx/2-2 || i > nx/2+1) {
                    const Real expected = i < (blended ? nx/2 : split) ? Real(0.01) : Real(0.1);
                    error = amrex::max(error, amrex::Math::abs(value - expected));
                }
                return {error, value, value};
            });
        }
        auto result = data.value();
        Real error = get<0>(result), lo = get<1>(result), hi = get<2>(result);
        ParallelDescriptor::ReduceRealMax(error);
        ParallelDescriptor::ReduceRealMin(lo);
        ParallelDescriptor::ReduceRealMax(hi);
        require(error <= tol, "Roughness map or vertical extrusion changed");
        require(std::abs(lo-Real(0.01)) <= tol && std::abs(hi-Real(0.1)) <= tol,
                "Both prescribed roughness regions must survive without overshoot");

        Real region_min[2], region_max[2];
        for (int region = 0; region < 2; ++region) {
            ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpMax, ReduceOpSum> stress_op;
            ReduceData<Real, Real, Real, Long> stress_data(stress_op);
            using StressTuple = typename decltype(stress_data)::Type;
            const Real expected = KAPPA * Real(10) /
                                  std::log(Real(10) / (region == 0 ? Real(0.01) : Real(0.1)));
            for (MFIter mfi(ustar); mfi.isValid(); ++mfi) {
                Box bx = mfi.validbox();
                if (bx.smallEnd(2) > kcut || bx.bigEnd(2) < kcut) { continue; }
                bx.setRange(2, kcut);
                // Stay away from the transition when interpolation is requested.
                bx.setSmall(0, amrex::max(bx.smallEnd(0), region == 0 ? 0 : nx/2+2));
                bx.setBig(0, amrex::min(bx.bigEnd(0), region == 0 ? nx/2-3 : nx-1));
                if (!bx.ok()) { continue; }
                const auto a = ustar.const_array(mfi);
                stress_op.eval(bx, stress_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> StressTuple {
                        return {a(i,j,k), a(i,j,k), amrex::Math::abs(a(i,j,k)-expected), Long(1)};
                    });
            }
            auto stress = stress_data.value();
            Real umin = get<0>(stress), umax = get<1>(stress), uerror = get<2>(stress);
            Long count = get<3>(stress);
            ParallelDescriptor::ReduceRealMin(umin);
            ParallelDescriptor::ReduceRealMax(umax);
            ParallelDescriptor::ReduceRealMax(uerror);
            ParallelDescriptor::ReduceLongSum(count);
            require(count > 0 && umin > zero && umax < Real(2),
                    "Each roughness region must have positive finite cut-cell friction velocity");
            if (initial) {
                require(uerror <= tol, "Varying roughness does not match neutral MOST");
            }
            region_min[region] = umin;
            region_max[region] = umax;
            Print() << "Level " << lev << " region " << region << ": u* = ["
                    << umin << ", " << umax << "], samples = " << count << "\n";
        }
        require(region_min[1] > region_max[0],
                "The rougher region must produce stronger MOST friction velocity");
        Print() << "Level " << lev << ": z0 = [" << lo << ", " << hi
                << "], maximum map/column error = " << error << "\n";
    }
}

void check_regional_drag (const std::string& active, const std::string& inactive)
{
    PlotFileData a(active), b(inactive);
    MultiFab u = a.get(0, "x_velocity"), uoff = b.get(0, "x_velocity");
    require(u.is_finite() && uoff.is_finite(), "Non-finite varying-roughness velocities");
    MultiFab difference(u.boxArray(), u.DistributionMap(), 1, 0);
    difference.ParallelCopy(uoff);
    MultiFab::Subtract(difference, u, 0, 0, 1, 0);
    Real reduction[2];
    for (int region = 0; region < 2; ++region) {
        // Interior cut cells avoid periodic seams and the material interface.
        const Box sample(IntVect(region == 0 ? 2 : 10, 2, 1),
                         IntVect(region == 0 ? 5 : 13, 13, 1));
        ReduceOps<ReduceOpMin, ReduceOpSum, ReduceOpSum> op;
        ReduceData<Real, Real, Long> data(op);
        using Tuple = typename decltype(data)::Type;
        for (MFIter mfi(difference); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox() & sample;
            if (!bx.ok()) { continue; }
            const auto d = difference.const_array(mfi);
            op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Tuple {
                return {d(i,j,k), d(i,j,k), Long(1)};
            });
        }
        auto result = data.value();
        Real min_drag = get<0>(result), sum = get<1>(result);
        Long count = get<2>(result);
        ParallelDescriptor::ReduceRealMin(min_drag);
        ParallelDescriptor::ReduceRealSum(sum);
        ParallelDescriptor::ReduceLongSum(count);
        require(count == sample.numPts() && min_drag > Real(1.e-5),
                "EB MOST must slow the wind in both material regions");
        reduction[region] = sum / Real(count);
    }
    require(reduction[1] > reduction[0],
            "The rougher region must experience stronger applied momentum drag");
    Print() << "Mean wind reduction: smooth=" << reduction[0] << ", rough="
            << reduction[1] << ", ratio=" << reduction[1]/reduction[0] << "\n";
}

} // namespace

int main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);
    {
        require(argc >= 4, "Usage: checker state chk finest [initial] | parity plt plt | drag plt plt | "
                          "roughness chk finest [initial|interpolated] | regional_drag plt plt | "
                          "roughness_remade old_chk new_chk [interpolated]");
        const std::string mode(argv[1]);
        if (mode == "state") {
            check_state(argv[2], std::stoi(argv[3]), argc > 4);
        } else if (mode == "parity") {
            compare_plotfiles(argv[2], argv[3]);
        } else if (mode == "drag") {
            check_drag(argv[2], argv[3]);
        } else if (mode == "roughness") {
            check_roughness(argv[2], std::stoi(argv[3]),
                            argc > 4 && std::string(argv[4]) == "initial",
                            argc > 4 && std::string(argv[4]) == "interpolated");
        } else if (mode == "roughness_remade") {
            MultiFab old_z0, new_z0;
            VisMF::Read(old_z0, std::string(argv[2]) + "/Level_1/Z0");
            VisMF::Read(new_z0, std::string(argv[3]) + "/Level_1/Z0");
            require(old_z0.boxArray() != new_z0.boxArray(),
                    "The roughness AMR leg must remake the existing fine-level grids");
            check_roughness(argv[3], 1, false, argc > 4 && std::string(argv[4]) == "interpolated");
        } else if (mode == "regional_drag") {
            check_regional_drag(argv[2], argv[3]);
        } else {
            Abort("Unknown EB SurfaceLayer check mode");
        }
        Print() << "EB SurfaceLayer " << mode << " check passed\n";
    }
    amrex::Finalize();
}
