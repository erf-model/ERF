#include "ERF_AuxiliaryMappedTransport.H"

#include <AMReX_Math.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <limits>

namespace erf_auxiliary {
namespace {

bool same_horizontal_layout (const amrex::MultiFab& cell_field,
                           const amrex::MultiFab& map_field)
{
    if (cell_field.DistributionMap() != map_field.DistributionMap() ||
        cell_field.boxArray().size() != map_field.boxArray().size()) {
        return false;
    }
    for (int index = 0; index < static_cast<int>(cell_field.boxArray().size()); ++index) {
        const auto& cell_box = cell_field.boxArray()[index];
        const auto& map_box = map_field.boxArray()[index];
        if (cell_box.smallEnd(0) != map_box.smallEnd(0) ||
            cell_box.bigEnd(0) != map_box.bigEnd(0) ||
            cell_box.smallEnd(1) != map_box.smallEnd(1) ||
            cell_box.bigEnd(1) != map_box.bigEnd(1) ||
            map_box.smallEnd(2) != 0 || map_box.bigEnd(2) != 0) {
            return false;
        }
    }
    return true;
}

} // namespace

bool SameCellLayout (const amrex::MultiFab& lhs, const amrex::MultiFab& rhs)
{
    return lhs.boxArray() == rhs.boxArray() &&
           lhs.DistributionMap() == rhs.DistributionMap();
}

bool BuildMappedCellMeasure (amrex::MultiFab& omega,
                             const amrex::MultiFab& detJ,
                             const amrex::MultiFab& mx,
                             const amrex::MultiFab& my,
                             std::string& diagnostic)
{
    diagnostic.clear();
    if (omega.nComp() != 1 || detJ.nComp() < 1 || mx.nComp() < 1 || my.nComp() < 1 ||
        !SameCellLayout(omega, detJ) || !same_horizontal_layout(omega, mx) ||
        !same_horizontal_layout(omega, my) ||
        mx.boxArray() != my.boxArray() || mx.DistributionMap() != my.DistributionMap()) {
        diagnostic = "mapped measure requires a cell-centered Jacobian, matching horizontal map-factor layouts, and one output component";
        return false;
    }

    amrex::MultiFab invalid(omega.boxArray(), omega.DistributionMap(), 1, 0);
    invalid.setVal(amrex::Real(0.0));
    for (amrex::MFIter mfi(omega, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto det = detJ.const_array(mfi);
        const auto map_x = mx.const_array(mfi);
        const auto map_y = my.const_array(mfi);
        const auto out = omega.array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real d = det(i, j, k, 0);
            // ERF's horizontal map factors are vertically invariant and all
            // native mapped divergence paths read them from k=0.
            const amrex::Real x = map_x(i, j, 0, 0);
            const amrex::Real y = map_y(i, j, 0, 0);
            const amrex::Real measure = d / (x * y);
            const bool valid = amrex::Math::isfinite(d) && d > amrex::Real(0.0) &&
                               amrex::Math::isfinite(x) && x > amrex::Real(0.0) &&
                               amrex::Math::isfinite(y) && y > amrex::Real(0.0) &&
                               amrex::Math::isfinite(measure) && measure > amrex::Real(0.0);
            out(i, j, k, 0) = valid ? measure : amrex::Real(0.0);
            bad(i, j, k, 0) = valid ? amrex::Real(0.0) : amrex::Real(1.0);
        });
    }

    const amrex::Real invalid_count = invalid.sum(0);
    if (invalid_count != amrex::Real(0.0)) {
        diagnostic = "detJ, horizontal map factors, and mapped measure must be finite and positive";
        return false;
    }
    return true;
}

bool BuildMappedDryAirCarrierFluxRate (
    MappedFaceFluxRate& rate,
    const amrex::MultiFab& rho_u,
    const amrex::MultiFab& rho_v,
    const amrex::MultiFab& omega,
    const amrex::MultiFab& ax,
    const amrex::MultiFab& ay,
    const amrex::MultiFab& az,
    const amrex::MultiFab& mf_uy,
    const amrex::MultiFab& mf_vx,
    const amrex::MultiFab& mx,
    const amrex::MultiFab& my,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (!rate.is_defined() || rate.nComp() != 1 || rho_u.nComp() < 1 ||
        rho_v.nComp() < 1 || omega.nComp() < 1 || ax.nComp() < 1 ||
        ay.nComp() < 1 || az.nComp() < 1 || mf_uy.nComp() < 1 ||
        mf_vx.nComp() < 1 || mx.nComp() < 1 || my.nComp() < 1 ||
        rho_u.boxArray() != rate.dir(0).boxArray() ||
        rho_u.DistributionMap() != rate.dir(0).DistributionMap() ||
        ax.boxArray() != rate.dir(0).boxArray() ||
        ax.DistributionMap() != rate.dir(0).DistributionMap() ||
        rho_v.boxArray() != rate.dir(1).boxArray() ||
        rho_v.DistributionMap() != rate.dir(1).DistributionMap() ||
        ay.boxArray() != rate.dir(1).boxArray() ||
        ay.DistributionMap() != rate.dir(1).DistributionMap() ||
        omega.boxArray() != rate.dir(2).boxArray() ||
        omega.DistributionMap() != rate.dir(2).DistributionMap() ||
        az.boxArray() != rate.dir(2).boxArray() ||
        az.DistributionMap() != rate.dir(2).DistributionMap() ||
        !same_horizontal_layout(rate.dir(0), mf_uy) ||
        !same_horizontal_layout(rate.dir(1), mf_vx) ||
        mx.boxArray() != my.boxArray() ||
        mx.DistributionMap() != my.DistributionMap() ||
        !same_horizontal_layout(rate.dir(2), mx) ||
        !same_horizontal_layout(rate.dir(2), my)) {
        diagnostic = "mapped dry-air carrier rates require matching staggered momentum, area, map-factor, and output layouts";
        return false;
    }

    for (amrex::MFIter mfi(rate.dir(0), amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto momentum = rho_u.const_array(mfi);
        const auto area = ax.const_array(mfi);
        const auto map = mf_uy.const_array(mfi);
        const auto out = rate.dir(0).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = area(i, j, k, 0) * momentum(i, j, k, 0) /
                              map(i, j, 0, 0);
        });
    }
    for (amrex::MFIter mfi(rate.dir(1), amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto momentum = rho_v.const_array(mfi);
        const auto area = ay.const_array(mfi);
        const auto map = mf_vx.const_array(mfi);
        const auto out = rate.dir(1).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = area(i, j, k, 0) * momentum(i, j, k, 0) /
                              map(i, j, 0, 0);
        });
    }
    for (amrex::MFIter mfi(rate.dir(2), amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto contravariant_momentum = omega.const_array(mfi);
        const auto area = az.const_array(mfi);
        const auto map_x = mx.const_array(mfi);
        const auto map_y = my.const_array(mfi);
        const auto out = rate.dir(2).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = area(i, j, k, 0) *
                              contravariant_momentum(i, j, k, 0) /
                              (map_x(i, j, 0, 0) * map_y(i, j, 0, 0));
        });
    }

    if (!rate.dir(0).is_finite(0, 1, 0) ||
        !rate.dir(1).is_finite(0, 1, 0) ||
        !rate.dir(2).is_finite(0, 1, 0)) {
        diagnostic = "mapped dry-air carrier face rates are nonfinite";
        return false;
    }
    return true;
}

bool CopyNativeMappedDryAirCarrierFluxRate (
    MappedFaceFluxRate& rate,
    const amrex::MultiFab& rho_u,
    const amrex::MultiFab& rho_v,
    const amrex::MultiFab& rho_w,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (AMREX_SPACEDIM != 3 || !rate.is_defined() || rate.nComp() != 1 ||
        rho_u.nComp() < 1 || rho_v.nComp() < 1 || rho_w.nComp() < 1) {
        diagnostic = "native mapped dry-air carrier requires three dimensions, "
                     "one output component, and one component per input face";
        return false;
    }

    const amrex::MultiFab* input_faces[3] = {&rho_u, &rho_v, &rho_w};
    for (int dir = 0; dir < 3; ++dir) {
        const auto& output = rate.dir(dir);
        const auto& input = *input_faces[dir];
        if (output.boxArray() != input.boxArray() ||
            output.DistributionMap() != input.DistributionMap()) {
            diagnostic = "native mapped dry-air carrier input must match the "
                         "output's x/y/z face staggering and distribution map";
            return false;
        }
    }

    for (int dir = 0; dir < 3; ++dir) {
        amrex::MultiFab::Copy(rate.dir(dir), *input_faces[dir], 0, 0, 1, 0);
    }
    if (!rate.dir(0).is_finite(0, 1, 0) ||
        !rate.dir(1).is_finite(0, 1, 0) ||
        !rate.dir(2).is_finite(0, 1, 0)) {
        diagnostic = "native mapped dry-air carrier face values are nonfinite";
        return false;
    }
    return true;
}

bool ComputeMaxMappedOutgoingRate (
    const MappedFaceFluxRate& rate,
    const amrex::MultiFab& measure,
    const amrex::MultiFab& conserved,
    const int density_component,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
    amrex::Real& max_rate,
    std::string& diagnostic)
{
    diagnostic.clear();
    max_rate = amrex::Real(0.0);
    if (!rate.is_defined() || rate.nComp() != 1 ||
        !MappedFaceLayoutMatchesCellLayout(rate, measure) ||
        !SameCellLayout(measure, conserved) || measure.nComp() < 1) {
        diagnostic = "mapped donor-rate reduction received incompatible cell or face layouts";
        return false;
    }
    if (!ValidatePositiveFiniteComponent(measure, 0, diagnostic)) return false;
    if (!ValidatePositiveFiniteComponent(conserved, density_component,
                                         diagnostic)) return false;

    amrex::MultiFab cell_rate(measure.boxArray(), measure.DistributionMap(), 1,
                              0);
    for (amrex::MFIter mfi(cell_rate, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto fx = rate.dir(0).const_array(mfi);
        const auto fy = rate.dir(1).const_array(mfi);
        const auto fz = rate.dir(2).const_array(mfi);
        const auto omega = measure.const_array(mfi);
        const auto rho = conserved.const_array(mfi);
        const auto out = cell_rate.array(mfi);
        const int rho_comp = density_component;
        const auto inv = dx_inv;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out(i, j, k, 0) = MappedOutgoingDemandRate(
                fx(i + 1, j, k, 0), fx(i, j, k, 0),
                fy(i, j + 1, k, 0), fy(i, j, k, 0),
                fz(i, j, k + 1, 0), fz(i, j, k, 0),
                omega(i, j, k, 0), rho(i, j, k, rho_comp), inv);
        });
    }
    max_rate = cell_rate.max(0);
    if (!std::isfinite(max_rate) || max_rate < amrex::Real(0.0)) {
        diagnostic = "mapped donor outgoing-rate maximum is nonfinite or negative";
        return false;
    }
    return true;
}

bool FixedDtExceedsMappedDonorLimit (const double fixed_dt,
                                     const amrex::Real max_outgoing_rate,
                                     double& hard_limit) noexcept
{
    hard_limit = std::numeric_limits<double>::infinity();
    if (std::isfinite(max_outgoing_rate) &&
        max_outgoing_rate > amrex::Real(0.0)) {
        hard_limit = 1.0 / static_cast<double>(max_outgoing_rate);
    }
    if (max_outgoing_rate == amrex::Real(0.0) || fixed_dt <= 0.0)
        return false;
    if (!std::isfinite(fixed_dt) ||
        !std::isfinite(max_outgoing_rate) ||
        max_outgoing_rate < amrex::Real(0.0)) {
        hard_limit = 0.0;
        return true;
    }

    if (!std::isfinite(hard_limit)) return false;
    const double scale = std::max(std::abs(fixed_dt), std::abs(hard_limit));
    const double tolerance = 64.0 *
        static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()) *
        scale;
    return fixed_dt - hard_limit > tolerance;
}

bool ValidatePositiveFiniteComponent (const amrex::MultiFab& field,
                                      const int component,
                                      std::string& diagnostic)
{
    diagnostic.clear();
    if (component < 0 || component >= field.nComp()) {
        diagnostic = "positive finite field validation received an invalid component";
        return false;
    }
    amrex::MultiFab invalid(field.boxArray(), field.DistributionMap(), 1, 0);
    invalid.setVal(amrex::Real(0.0));
    for (amrex::MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto values = field.const_array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real value = values(i, j, k, component);
            const bool valid = amrex::Math::isfinite(value) && value > amrex::Real(0.0);
            bad(i, j, k, 0) = valid ? amrex::Real(0.0) : amrex::Real(1.0);
        });
    }
    if (invalid.sum(0) != amrex::Real(0.0)) {
        diagnostic = "selected field component must be finite and positive in every valid cell";
        return false;
    }
    return true;
}

bool ValidateFiniteComponent (const amrex::MultiFab& field,
                              const int component,
                              std::string& diagnostic)
{
    diagnostic.clear();
    if (component < 0 || component >= field.nComp()) {
        diagnostic = "finite field validation received an invalid component";
        return false;
    }
    amrex::MultiFab invalid(field.boxArray(), field.DistributionMap(), 1, 0);
    invalid.setVal(amrex::Real(0.0));
    for (amrex::MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto values = field.const_array(mfi);
        const auto bad = invalid.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            bad(i, j, k, 0) = amrex::Math::isfinite(values(i, j, k, component)) ?
                amrex::Real(0.0) : amrex::Real(1.0);
        });
    }
    if (invalid.sum(0) != amrex::Real(0.0)) {
        diagnostic = "selected field component must be finite in every valid cell";
        return false;
    }
    return true;
}

amrex::Real MaxFaceFieldDifference (const MappedFaceFluxRate& lhs,
                                  const int lhs_comp,
                                  const MappedFaceFluxRate& rhs,
                                  const int rhs_comp)
{
    AMREX_ALWAYS_ASSERT(lhs.is_defined() && rhs.is_defined());
    AMREX_ALWAYS_ASSERT(lhs_comp >= 0 && lhs_comp < lhs.nComp());
    AMREX_ALWAYS_ASSERT(rhs_comp >= 0 && rhs_comp < rhs.nComp());
    AMREX_ALWAYS_ASSERT(SameMappedFaceLayout(lhs, rhs));
    amrex::Real maximum = amrex::Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& a = lhs.dir(dir);
        const auto& b = rhs.dir(dir);
        amrex::MultiFab difference(a.boxArray(), a.DistributionMap(), 1, 0);
        for (amrex::MFIter mfi(difference, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto av = a.const_array(mfi);
            const auto bv = b.const_array(mfi);
            const auto out = difference.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                out(i, j, k, 0) = amrex::Math::abs(av(i, j, k, lhs_comp) -
                                                    bv(i, j, k, rhs_comp));
            });
        }
        maximum = amrex::max(maximum, difference.norm0(0));
    }
    return maximum;
}

void AccumulateIntegratedFaceFlux (IntegratedMappedFaceFlux& ledger,
                                 const MappedFaceFluxRate& rate,
                                 const amrex::Real weight)
{
    AMREX_ALWAYS_ASSERT(ledger.is_defined() && rate.is_defined());
    AMREX_ALWAYS_ASSERT(ledger.nComp() == rate.nComp());
    AMREX_ALWAYS_ASSERT(SameMappedFaceLayout(ledger, rate));
    AMREX_ALWAYS_ASSERT(std::isfinite(weight) && weight >= amrex::Real(0.0));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Saxpy(ledger.dir(dir), weight, rate.dir(dir),
                               0, 0, ledger.nComp(), 0);
    }
}

} // namespace erf_auxiliary
