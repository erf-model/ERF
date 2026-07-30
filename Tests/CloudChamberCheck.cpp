#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "../Source/Utils/ERF_EOS.H"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::PlotFileData;
using amrex::Real;
using amrex::TilingIfNotGPU;

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

int fail (const std::string& message)
{
    std::cerr << "CloudChamberCheck: " << message << "\n";
    return 1;
}

} // namespace

int main (int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: checker mode initial_plotfile final_plotfile\n";
        return 2;
    }

    amrex::Initialize(argc, argv, false);
    const std::string mode(argv[1]);
    PlotFileData initial(argv[2]);
    PlotFileData final(argv[3]);
    const bool cloudy = (mode == "cloudy");
    if (!cloudy && mode != "dry") {
        amrex::Finalize();
        return fail("mode must be dry or cloudy");
    }

    for (const std::string& name : {"density", "theta", "temp", "x_velocity",
                                    "y_velocity", "z_velocity"}) {
        if (!has_variable(initial, name) || !has_variable(final, name)) {
            amrex::Finalize();
            return fail("missing required variable " + name);
        }
    }
    if (cloudy) {
        for (const std::string& name : {"qv", "qc"}) {
            if (!has_variable(initial, name) || !has_variable(final, name)) {
                amrex::Finalize();
                return fail("missing cloudy variable " + name);
            }
        }
    } else if (has_variable(initial, "qv") || has_variable(initial, "qc")) {
        amrex::Finalize();
        return fail("dry case unexpectedly contains moisture variables");
    }

    MultiFab density = initial.get(0, "density");
    MultiFab theta = initial.get(0, "theta");
    MultiFab temp = initial.get(0, "temp");
    MultiFab initial_qv;
    if (cloudy) {
        initial_qv = initial.get(0, "qv");
    }
    MultiFab x_velocity = initial.get(0, "x_velocity");
    MultiFab y_velocity = initial.get(0, "y_velocity");
    MultiFab z_velocity = initial.get(0, "z_velocity");
    if (!density.is_finite() || !theta.is_finite() || !temp.is_finite() ||
        !x_velocity.is_finite() || !y_velocity.is_finite() ||
        !z_velocity.is_finite()) {
        amrex::Finalize();
        return fail("initial required field is non-finite");
    }

    const auto domain = initial.probDomain(0);
    const auto prob_lo = initial.probLo();
    const auto dx = initial.cellSize(0);
    const Real lx = dx[0] * Real(domain.length(0));
    const Real ly = dx[1] * Real(domain.length(1));
    const Real lz = dx[2] * Real(domain.length(2));
    const Real theta_bottom = Real(299.0);
    const Real theta_top = Real(280.0);
    const Real amplitude = Real(0.05);
    const Real qv_bottom = Real(0.030);
    const Real qv_top = Real(0.030);
    MultiFab theta_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab temperature_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab qv_error;
    if (cloudy) {
        qv_error.define(theta.boxArray(), theta.DistributionMap(), 1, 0);
    }
    MultiFab velocity_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(theta_error, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto bx = mfi.tilebox();
        const auto rho = density.const_array(mfi);
        const auto th = theta.const_array(mfi);
        const auto t = temp.const_array(mfi);
        const auto qv = cloudy ? initial_qv.const_array(mfi) : amrex::Array4<const Real>{};
        const auto u = x_velocity.const_array(mfi);
        const auto v = y_velocity.const_array(mfi);
        const auto w = z_velocity.const_array(mfi);
        const auto therr = theta_error.array(mfi);
        const auto terr = temperature_error.array(mfi);
        const auto qverr = cloudy ? qv_error.array(mfi) : amrex::Array4<Real>{};
        const auto velerr = velocity_error.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real x = prob_lo[0] + (Real(i-domain.smallEnd(0)) + Real(0.5))*dx[0];
            const Real y = prob_lo[1] + (Real(j-domain.smallEnd(1)) + Real(0.5))*dx[1];
            const Real z = prob_lo[2] + (Real(k-domain.smallEnd(2)) + Real(0.5))*dx[2];
            const Real expected = theta_bottom + (theta_top-theta_bottom)*(z-prob_lo[2])/lz +
                amplitude * std::sin(Real(2.0)*Real(3.14159265358979323846)*(x-prob_lo[0])/lx) *
                std::sin(Real(2.0)*Real(3.14159265358979323846)*(y-prob_lo[1])/ly) *
                std::sin(Real(3.14159265358979323846)*(z-prob_lo[2])/lz);
            const Real expected_qv = qv_bottom + (qv_top-qv_bottom)*(z-prob_lo[2])/lz;
            const Real expected_temp = getTgivenRandRTh(
                rho(i,j,k), rho(i,j,k)*expected, cloudy ? expected_qv : Real(0.0));
            therr(i,j,k) = std::abs(th(i,j,k)-expected);
            terr(i,j,k) = std::abs(t(i,j,k)-expected_temp);
            if (cloudy) {
                qverr(i,j,k) = std::abs(qv(i,j,k)-expected_qv);
            }
            velerr(i,j,k) = std::max(std::abs(u(i,j,k)),
                                     std::max(std::abs(v(i,j,k)),std::abs(w(i,j,k))));
        });
    }
    amrex::Gpu::streamSynchronize();

    const Real theta_error_max = theta_error.norm0(0, 0, false);
    const Real temperature_error_max = temperature_error.norm0(0, 0, false);
    const Real qv_error_max = cloudy ? qv_error.norm0(0, 0, false) : Real(0.0);
    const Real initial_velocity_max = velocity_error.norm0(0, 0, false);
    const Real theta_tolerance = Real(2.0e-10) * std::max(Real(1.0), theta_bottom);
    if (theta_error_max > theta_tolerance) {
        amrex::Finalize();
        return fail("initial theta profile mismatch: max error=" +
                    std::to_string(static_cast<double>(theta_error_max)));
    }
    const Real temperature_tolerance = Real(2.0e-10) * Real(300.0);
    if (temperature_error_max > temperature_tolerance) {
        amrex::Finalize();
        return fail("initial temperature mismatch: max error=" +
                    std::to_string(static_cast<double>(temperature_error_max)));
    }
    if (cloudy && qv_error_max > Real(2.0e-12)) {
        amrex::Finalize();
        return fail("initial qv profile mismatch: max error=" +
                    std::to_string(static_cast<double>(qv_error_max)));
    }
    if (initial_velocity_max > Real(2.0e-12)) {
        amrex::Finalize();
        return fail("initial velocity is not zero: max=" +
                    std::to_string(static_cast<double>(initial_velocity_max)));
    }

    for (const std::string& name : {"density", "theta", "temp", "x_velocity",
                                    "y_velocity", "z_velocity"}) {
        if (!final.get(0, name).is_finite()) {
            amrex::Finalize();
            return fail("final field is non-finite: " + name);
        }
    }

    if (cloudy) {
        MultiFab qv = final.get(0, "qv");
        MultiFab qc = final.get(0, "qc");
        if (!qv.is_finite() || !qc.is_finite()) {
            amrex::Finalize();
            return fail("cloudy scalar is non-finite");
        }
        if (qv.min(0) < Real(-1.0e-12) || qc.min(0) < Real(-1.0e-12)) {
            amrex::Finalize();
            return fail("cloudy scalar became negative");
        }
        if (qc.max(0) <= Real(1.0e-12)) {
            amrex::Finalize();
            return fail("SatAdj did not form cloud water above roundoff");
        }
    }

    const MultiFab final_u = final.get(0, "x_velocity");
    const MultiFab final_v = final.get(0, "y_velocity");
    const MultiFab final_w = final.get(0, "z_velocity");
    const Real evolved_velocity = std::max(final_u.norm0(0,0,false),
        std::max(final_v.norm0(0,0,false), final_w.norm0(0,0,false)));
    std::cout << "mode=" << mode << " initial_theta_error=" << theta_error_max
              << " initial_temperature_error=" << temperature_error_max
              << " evolved_velocity_max=" << evolved_velocity;
    if (cloudy) {
        std::cout << " final_qv_min=" << final.get(0,"qv").min(0)
                  << " final_qc_max=" << final.get(0,"qc").max(0);
    }
    std::cout << "\n";

    amrex::Finalize();
    return evolved_velocity > Real(1.0e-12) ? 0 : fail("thermal perturbation produced no buoyant response");
}
