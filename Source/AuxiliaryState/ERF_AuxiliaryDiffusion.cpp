#include "ERF_AuxiliaryDiffusion.H"

#include <AMReX_Math.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Reduce.H>

#include <cmath>

namespace erf_auxiliary {
namespace {

using amrex::Array4;
using amrex::Box;
using amrex::MultiFab;
using amrex::Real;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
on_physical_face (const int index, const int lo, const int hi, const bool periodic) noexcept
{
    return !periodic && (index == lo || index == hi + 1);
}

// Counts nonfinite or negative coefficients, and nonpositive densities scaled
// by a kinematic coefficient, in every cell a non-physical face reads.
int count_invalid_coefficient_reads (const amrex::Geometry& geom,
                                     const MappedDiffusionCoefficients& coefficients,
                                     const MultiFab& density,
                                     const int density_component,
                                     const MultiFab* eddy_diffusivity)
{
    const Box& domain = geom.Domain();
    const bool use_eddy = coefficients.mode.add_turbulence;
    const bool use_density = coefficients.mode.molecular_is_kinematic;
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<int> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(density, false); mfi.isValid(); ++mfi) {
        const auto rho = density.const_array(mfi);
        const auto mu = use_eddy ? eddy_diffusivity->const_array(mfi)
                                 : Array4<const Real>{};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            Box read = amrex::grow(mfi.validbox(), dir, 1);
            if (!geom.isPeriodic(dir)) {
                read.setSmall(dir, amrex::max(read.smallEnd(dir), domain.smallEnd(dir)));
                read.setBig(dir, amrex::min(read.bigEnd(dir), domain.bigEnd(dir)));
            }
            const int eddy_comp = dir < 2 ? coefficients.horizontal.eddy_h_comp
                                          : coefficients.vertical.eddy_v_comp;
            const int rho_comp = density_component;
            reduce_op.eval(read, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> ReduceTuple {
                    int bad = 0;
                    if (use_eddy) {
                        const Real value = mu(i, j, k, eddy_comp);
                        if (!(amrex::Math::isfinite(value) && value >= Real(0.0))) { bad = 1; }
                    }
                    if (use_density) {
                        const Real value = rho(i, j, k, rho_comp);
                        if (!(amrex::Math::isfinite(value) && value > Real(0.0))) { bad = 1; }
                    }
                    return {bad};
                });
        }
    }
    int count = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelAllReduce::Sum(count, amrex::ParallelContext::CommunicatorSub());
    return count;
}

template <bool MultiplyMolecularByDensity, bool AddTurbulence>
void fill_transmissibility (MappedFaceTransmissibility& transmissibility,
                            const amrex::Geometry& geom,
                            const MappedDiffusionMetrics& metrics,
                            const MappedDiffusionCoefficients& coefficients,
                            const MultiFab& density,
                            const int density_component,
                            const MultiFab* eddy_diffusivity)
{
    const Box& domain = geom.Domain();
    const auto dom_lo = amrex::lbound(domain);
    const auto dom_hi = amrex::ubound(domain);
    const auto dx_inv = geom.InvCellSizeArray();
    const bool periodic_x = geom.isPeriodic(0);
    const bool periodic_y = geom.isPeriodic(1);
    const bool periodic_z = geom.isPeriodic(2);
    const bool stretched = metrics.mesh == DiffusionMesh::StretchedDz;
    const Real* dz = metrics.stretched_dz;
    const auto horizontal = coefficients.horizontal;
    const auto vertical = coefficients.vertical;
    const int rho_comp = density_component;

    auto& tx_face = transmissibility.dir(0);
    for (amrex::MFIter mfi(tx_face, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto rho = density.const_array(mfi);
        const auto mu = AddTurbulence ? eddy_diffusivity->const_array(mfi)
                                      : Array4<const Real>{};
        const auto map_num = metrics.mf_ux->const_array(mfi);
        const auto map_den = metrics.mf_uy->const_array(mfi);
        const auto area = stretched ? metrics.ax->const_array(mfi) : Array4<const Real>{};
        const auto tx = tx_face.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (on_physical_face(i, dom_lo.x, dom_hi.x, periodic_x)) {
                tx(i, j, k, 0) = Real(0.0);
                return;
            }
            const Real coeff = ScalarDiffusionFaceCoefficient<
                MultiplyMolecularByDensity, AddTurbulence>(
                rho, rho_comp, mu, horizontal, i, j, k, 1, 0, 0, horizontal.eddy_h_comp);
            tx(i, j, k, 0) = ScalarDiffusionMappedHorizontalTransfer(
                coeff * dx_inv[0] * map_num(i, j, 0) / map_den(i, j, 0),
                stretched ? area(i, j, k) : Real(1.0));
        });
    }

    auto& ty_face = transmissibility.dir(1);
    for (amrex::MFIter mfi(ty_face, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto rho = density.const_array(mfi);
        const auto mu = AddTurbulence ? eddy_diffusivity->const_array(mfi)
                                      : Array4<const Real>{};
        const auto map_num = metrics.mf_vy->const_array(mfi);
        const auto map_den = metrics.mf_vx->const_array(mfi);
        const auto area = stretched ? metrics.ay->const_array(mfi) : Array4<const Real>{};
        const auto ty = ty_face.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (on_physical_face(j, dom_lo.y, dom_hi.y, periodic_y)) {
                ty(i, j, k, 0) = Real(0.0);
                return;
            }
            const Real coeff = ScalarDiffusionFaceCoefficient<
                MultiplyMolecularByDensity, AddTurbulence>(
                rho, rho_comp, mu, horizontal, i, j, k, 0, 1, 0, horizontal.eddy_h_comp);
            ty(i, j, k, 0) = ScalarDiffusionMappedHorizontalTransfer(
                coeff * dx_inv[1] * map_num(i, j, 0) / map_den(i, j, 0),
                stretched ? area(i, j, k) : Real(1.0));
        });
    }

    auto& tz_face = transmissibility.dir(2);
    for (amrex::MFIter mfi(tz_face, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto rho = density.const_array(mfi);
        const auto mu = AddTurbulence ? eddy_diffusivity->const_array(mfi)
                                      : Array4<const Real>{};
        const auto mx = metrics.mf_mx->const_array(mfi);
        const auto my = metrics.mf_my->const_array(mfi);
        const auto tz = tz_face.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (on_physical_face(k, dom_lo.z, dom_hi.z, periodic_z)) {
                tz(i, j, k, 0) = Real(0.0);
                return;
            }
            const Real coeff = ScalarDiffusionFaceCoefficient<
                MultiplyMolecularByDensity, AddTurbulence>(
                rho, rho_comp, mu, vertical, i, j, k, 0, 0, 1, vertical.eddy_v_comp);
            // The host S stencil spans the two cell centres of an interior face.
            const Real inv_spacing = stretched ? Real(2.0) / (dz[k] + dz[k - 1])
                                               : dx_inv[2];
            tz(i, j, k, 0) = ScalarDiffusionMappedVerticalTransfer(
                coeff * inv_spacing, mx(i, j, 0), my(i, j, 0));
        });
    }
}

bool same_face_layout (const MultiFab& face, const MultiFab* field)
{
    return field != nullptr && field->boxArray() == face.boxArray() &&
           field->DistributionMap() == face.DistributionMap();
}

} // namespace

bool BuildMappedDiffusionTransmissibility (
    MappedFaceTransmissibility& transmissibility,
    const amrex::Geometry& geom,
    const MappedDiffusionMetrics& metrics,
    const MappedDiffusionCoefficients& coefficients,
    const MultiFab& density,
    const int density_component,
    const MultiFab* eddy_diffusivity,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (AMREX_SPACEDIM != 3 || !transmissibility.is_defined() ||
        transmissibility.nComp() != 1 ||
        !MappedFaceLayoutMatchesCellLayout(transmissibility, density)) {
        diagnostic = "diffusion transmissibility requires three dimensions, one "
                     "component, and faces matching the density layout";
        return false;
    }
    if (density_component < 0 || density_component >= density.nComp() ||
        !density.nGrowVect().allGE(amrex::IntVect(1))) {
        diagnostic = "diffusion transmissibility requires a valid density component "
                     "with at least one ghost cell";
        return false;
    }
    for (const Real value : {coefficients.horizontal.molecular_coeff,
                             coefficients.vertical.molecular_coeff}) {
        if (!(std::isfinite(value) && value >= Real(0.0))) {
            diagnostic = "diffusion molecular coefficients must be finite and nonnegative";
            return false;
        }
    }
    if (coefficients.mode.add_turbulence) {
        if (eddy_diffusivity == nullptr || !SameCellLayout(*eddy_diffusivity, density) ||
            !eddy_diffusivity->nGrowVect().allGE(amrex::IntVect(1))) {
            diagnostic = "diffusion with turbulent coefficients requires the host eddy "
                         "diffusivity on the density layout with at least one ghost cell";
            return false;
        }
        for (const int comp : {coefficients.horizontal.eddy_h_comp,
                               coefficients.vertical.eddy_v_comp}) {
            if (comp < 0 || comp >= eddy_diffusivity->nComp()) {
                diagnostic = "diffusion eddy coefficient component is out of range";
                return false;
            }
        }
    }
    const auto& tx = transmissibility.dir(0);
    const auto& ty = transmissibility.dir(1);
    const auto& tz = transmissibility.dir(2);
    if (metrics.mf_ux == nullptr || metrics.mf_uy == nullptr ||
        metrics.mf_vx == nullptr || metrics.mf_vy == nullptr ||
        metrics.mf_mx == nullptr || metrics.mf_my == nullptr ||
        !SameHorizontalLayout(tx, *metrics.mf_ux) || !SameHorizontalLayout(tx, *metrics.mf_uy) ||
        !SameHorizontalLayout(ty, *metrics.mf_vx) || !SameHorizontalLayout(ty, *metrics.mf_vy) ||
        !SameHorizontalLayout(tz, *metrics.mf_mx) || !SameHorizontalLayout(tz, *metrics.mf_my)) {
        diagnostic = "diffusion transmissibility requires ERF's staggered k=0 map factors";
        return false;
    }
    if (metrics.mesh == DiffusionMesh::StretchedDz) {
        if (!same_face_layout(tx, metrics.ax) || !same_face_layout(ty, metrics.ay) ||
            metrics.stretched_dz == nullptr) {
            diagnostic = "stretched diffusion transmissibility requires ax, ay and the "
                         "stretched vertical spacing";
            return false;
        }
        if (geom.isPeriodic(2)) {
            diagnostic = "stretched diffusion transmissibility has no cell-centre spacing "
                         "across a periodic z face";
            return false;
        }
    }

    const int invalid = count_invalid_coefficient_reads(
        geom, coefficients, density, density_component, eddy_diffusivity);
    if (invalid != 0) {
        diagnostic = "diffusion coefficients must be finite and nonnegative, and densities "
                     "finite and positive, in every cell a face reads (" +
                     std::to_string(invalid) + " invalid reads)";
        return false;
    }

    const auto& mode = coefficients.mode;
    if (mode.molecular_is_kinematic) {
        if (mode.add_turbulence) {
            fill_transmissibility<true, true>(transmissibility, geom, metrics, coefficients,
                                              density, density_component, eddy_diffusivity);
        } else {
            fill_transmissibility<true, false>(transmissibility, geom, metrics, coefficients,
                                               density, density_component, eddy_diffusivity);
        }
    } else if (mode.add_turbulence) {
        fill_transmissibility<false, true>(transmissibility, geom, metrics, coefficients,
                                           density, density_component, eddy_diffusivity);
    } else {
        fill_transmissibility<false, false>(transmissibility, geom, metrics, coefficients,
                                            density, density_component, eddy_diffusivity);
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& face = transmissibility.dir(dir);
        if (!face.is_finite(0, 1, 0) || face.min(0) < Real(0.0)) {
            diagnostic = "diffusion transmissibility is nonfinite or negative; check the "
                         "map factors, face areas and vertical spacing";
            return false;
        }
    }
    return true;
}

bool BuildLowOrderDiffusionRate (
    MappedFaceFluxRate& rate,
    const MappedFaceTransmissibility& transmissibility,
    const amrex::Geometry& geom,
    const MultiFab& intensive,
    const std::vector<int>& components,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (!rate.is_defined() || !transmissibility.is_defined() ||
        transmissibility.nComp() != 1 ||
        static_cast<int>(components.size()) > rate.nComp() ||
        !SameMappedFaceLayout(rate, transmissibility) ||
        !MappedFaceLayoutMatchesCellLayout(transmissibility, intensive) ||
        !intensive.nGrowVect().allGE(amrex::IntVect(1))) {
        diagnostic = "low-order diffusion rates require matching face layouts, enough "
                     "rate slots, and an intensive field with at least one ghost cell";
        return false;
    }
    for (const int comp : components) {
        if (comp < 0 || comp >= intensive.nComp()) {
            diagnostic = "low-order diffusion component is out of range";
            return false;
        }
    }

    const Box& domain = geom.Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const int lo = domain.smallEnd(dir);
        const int hi = domain.bigEnd(dir);
        const bool periodic = geom.isPeriodic(dir);
        const amrex::IntVect shift = amrex::IntVect::TheDimensionVector(dir);
        auto& face = rate.dir(dir);
        for (amrex::MFIter mfi(face, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box bx = mfi.tilebox();
            const auto t = transmissibility.dir(dir).const_array(mfi);
            const auto z = intensive.const_array(mfi);
            const auto out = face.array(mfi);
            for (std::size_t local = 0; local < components.size(); ++local) {
                const int slot = static_cast<int>(local);
                const int comp = components[local];
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::IntVect cell(AMREX_D_DECL(i, j, k));
                    if (on_physical_face(cell[dir], lo, hi, periodic)) {
                        out(i, j, k, slot) = Real(0.0);
                        return;
                    }
                    out(i, j, k, slot) =
                        -t(i, j, k, 0) * (z(cell, comp) - z(cell - shift, comp));
                });
            }
        }
    }
    const int nslots = static_cast<int>(components.size());
    for (int dir = 0; dir < AMREX_SPACEDIM && nslots > 0; ++dir) {
        if (!rate.dir(dir).is_finite(0, nslots, 0)) {
            diagnostic = "low-order diffusion rates are nonfinite; check the intensive field, "
                         "including the ghost cells read across periodic faces";
            return false;
        }
    }
    return true;
}

bool ComputeMaxMappedCombinedOutgoingRate (
    const MappedFaceFluxRate& carrier,
    const MappedFaceTransmissibility& transmissibility,
    const MultiFab& measure,
    const MultiFab& conserved,
    const int density_component,
    const amrex::GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    Real& max_rate,
    std::string& diagnostic)
{
    diagnostic.clear();
    max_rate = Real(0.0);
    if (!carrier.is_defined() || carrier.nComp() != 1 ||
        !transmissibility.is_defined() || transmissibility.nComp() != 1 ||
        !MappedFaceLayoutMatchesCellLayout(carrier, measure) ||
        !SameMappedFaceLayout(carrier, transmissibility) ||
        !SameCellLayout(measure, conserved) || measure.nComp() < 1) {
        diagnostic = "combined donor-rate reduction received incompatible cell or face layouts";
        return false;
    }
    if (!ValidatePositiveFiniteComponent(measure, 0, diagnostic)) { return false; }
    if (!ValidatePositiveFiniteComponent(conserved, density_component, diagnostic)) {
        return false;
    }

    MultiFab cell_rate(measure.boxArray(), measure.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(cell_rate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto fx = carrier.dir(0).const_array(mfi);
        const auto fy = carrier.dir(1).const_array(mfi);
        const auto fz = carrier.dir(2).const_array(mfi);
        const auto tx = transmissibility.dir(0).const_array(mfi);
        const auto ty = transmissibility.dir(1).const_array(mfi);
        const auto tz = transmissibility.dir(2).const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto rho = conserved.const_array(mfi);
        const auto out = cell_rate.array(mfi);
        const int rho_comp = density_component;
        const auto inv = dx_inv;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = MappedCombinedOutgoingDemandRate(
                fx(i + 1, j, k, 0), fx(i, j, k, 0),
                fy(i, j + 1, k, 0), fy(i, j, k, 0),
                fz(i, j, k + 1, 0), fz(i, j, k, 0),
                tx(i + 1, j, k, 0), tx(i, j, k, 0),
                ty(i, j + 1, k, 0), ty(i, j, k, 0),
                tz(i, j, k + 1, 0), tz(i, j, k, 0),
                omega(i, j, k, 0), rho(i, j, k, rho_comp), inv);
        });
    }
    max_rate = cell_rate.max(0);
    if (!std::isfinite(max_rate) || max_rate < Real(0.0)) {
        diagnostic = "combined donor outgoing-rate maximum is nonfinite or negative";
        return false;
    }
    return true;
}

} // namespace erf_auxiliary
