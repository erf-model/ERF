/**
 * \file ERF_Plotfile2DCloudOptics.cpp
 */
#include "ERF_Plotfile2DCloudOptics.H"
#include "ERF_PlanarBoundary.H"
#include "ERF_TileNoZ.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <cmath>
#include <limits>

using namespace amrex;

namespace plotfile2d {
namespace {
constexpr const char* coefficient_keys[] = {"kext_qc", "kext_qi", "kext_qr",
                                            "kext_qs", "kext_qg"};

GpuArray<int, 5>
mass_components (const MoistureComponentIndices& mi) noexcept
{
    return {{mi.qc, mi.qi, mi.qr, mi.qs, mi.qg}};
}
} // namespace

std::string
cloud_optics_parameter_error (const CloudOpticsParameters& parameters)
{
    for (int n = 0; n < 5; ++n) {
        if (!std::isfinite(parameters.kext[n]) || parameters.kext[n] < 0.) {
            return std::string("erf.plot2d.cloud_optics.") +
                   coefficient_keys[n] + " must be finite and nonnegative";
        }
    }
    if (!std::isfinite(parameters.reflectance_scale) ||
        parameters.reflectance_scale <= 0.) {
        return "erf.plot2d.cloud_optics.reflectance_scale must be finite and "
               "positive";
    }
    return {};
}

CloudOpticsParameters
parse_cloud_optics_parameters ()
{
    CloudOpticsParameters parameters;
    ParmParse pp("erf.plot2d.cloud_optics");
    for (int n = 0; n < 5; ++n) {
        pp.query(coefficient_keys[n], parameters.kext[n]);
    }
    pp.query("reflectance_scale", parameters.reflectance_scale);
    const auto error = cloud_optics_parameter_error(parameters);
    if (!error.empty()) {
        Abort(error);
    }
    // The writer caches this result, so runtime settings are logged once,
    // on first requested use, including when only a station needs the fields.
    amrex::Print output;
    output << "Plot2D cloud optics: illustrative kext [m^2/kg]";
    for (int n = 0; n < 5; ++n) {
        output << " " << coefficient_keys[n] << "=" << parameters.kext[n];
    }
    output << " reflectance_scale=" << parameters.reflectance_scale << "\n";
    return parameters;
}

bool
cloud_optics_available (const MoistureComponentIndices& indices) noexcept
{
    for (int comp : mass_components(indices)) {
        if (comp >= 0) {
            return true;
        }
    }
    return false;
}

void
fill_cloud_optics (MultiFab& dst,
                   int tau_component,
                   int albedo_component,
                   const MultiFab& cons,
                   const MoistureComponentIndices& indices,
                   const CloudOpticsParameters& parameters,
                   const Geometry& geom,
                   const MultiFab& detJ)
{
    if (tau_component < 0 && albedo_component < 0) {
        return;
    }
    const auto error = cloud_optics_parameter_error(parameters);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error.empty(), error);
    AMREX_ALWAYS_ASSERT(tau_component < dst.nComp() &&
                        albedo_component < dst.nComp());
    AMREX_ALWAYS_ASSERT(tau_component < 0 || tau_component != albedo_component);
    const auto components = mass_components(indices);
    for (int comp : components) {
        AMREX_ALWAYS_ASSERT(comp < cons.nComp());
    }
    if (!cloud_optics_available(indices)) {
        if (tau_component >= 0) {
            dst.setVal(-999., tau_component, 1, 0);
        }
        if (albedo_component >= 0) {
            dst.setVal(-999., albedo_component, 1, 0);
        }
        return;
    }

    // Preserve Plot2D's per-level water-path convention: cons may cover only
    // part of the vertical domain. Such a bottom-reaching level reports its
    // partial column, without reconstructing or masking missing upper cells.
    MultiFab tau(dst.boxArray(), dst.DistributionMap(), 1, 0);
    tau.setVal(0.);
    const bool use_metric = SolverChoice::mesh_type != MeshType::ConstantDz;
    const auto coefficients = parameters.kext;
    const Real invalid = std::numeric_limits<Real>::quiet_NaN();

    // Match water paths: columns must stay inside one host tile. GPU atomics
    // accumulate per-cell contributions into each box's partial column.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cons, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto src = cons.const_array(mfi);
        const auto out = tau.array(mfi);
        const auto metric =
            use_metric ? detJ.const_array(mfi) : Array4<const Real>{};
        ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        Real weighted_mass = 0.;
                        for (int n = 0; n < 5; ++n) {
                            if (components[n] >= 0) {
                                const Real mass = src(i, j, k, components[n]);
                                // Check before clipping: max(0, NaN) could
                                // silently hide bad state.
                                if (!Math::isfinite(mass)) {
                                    weighted_mass = invalid;
                                    break;
                                }
                                weighted_mass += coefficients[n] *
                                                 amrex::max(Real(0.), mass);
                            }
                        }
                        if (use_metric) {
                            const Real jac = metric(i, j, k);
                            weighted_mass = Math::isfinite(jac) && jac > 0.
                                                ? weighted_mass * jac
                                                : invalid;
                        }
                        HostDevice::Atomic::Add(&out(i, j, 0), weighted_mass);
                    });
    }
    tau.mult(geom.CellSize(2), 0, 1, 0);
    // Never transform partial vertical sums: reflectance is nonlinear.
    sum_over_stacked_boxes(tau, 0, 1);
    const Real scale = parameters.reflectance_scale;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto in = tau.const_array(mfi);
        const auto out = dst.array(mfi);
        ParallelFor(mfi.tilebox(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const Real value = in(i, j, k);
                        const bool valid = Math::isfinite(value) && value >= 0.;
                        if (tau_component >= 0) {
                            out(i, j, k, tau_component) = valid ? value : -999.;
                        }
                        if (albedo_component >= 0) {
                            out(i, j, k, albedo_component) =
                                valid ? cloud_reflectance(value, scale) : -999.;
                        }
                    });
    }
}
} // namespace plotfile2d
