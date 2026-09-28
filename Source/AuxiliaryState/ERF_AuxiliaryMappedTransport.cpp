#include "ERF_AuxiliaryMappedTransport.H"

#include <AMReX_Math.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>

namespace erf_auxiliary {
namespace {

bool same_cell_layout(const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    return a.boxArray() == b.boxArray() && a.DistributionMap() == b.DistributionMap();
}

bool same_horizontal_layout(const amrex::MultiFab& cell_field,
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

bool BuildMappedCellMeasure(amrex::MultiFab& omega,
                            const amrex::MultiFab& detJ,
                            const amrex::MultiFab& mx,
                            const amrex::MultiFab& my,
                            std::string& diagnostic)
{
    diagnostic.clear();
    if (omega.nComp() != 1 || detJ.nComp() < 1 || mx.nComp() < 1 || my.nComp() < 1 ||
        !same_cell_layout(omega, detJ) || !same_horizontal_layout(omega, mx) ||
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

bool ValidatePositiveFiniteComponent(const amrex::MultiFab& field,
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

bool ValidateFiniteComponent(const amrex::MultiFab& field,
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

amrex::Real MaxFaceFieldDifference(const MappedFaceFluxRate& lhs,
                                  const int lhs_comp,
                                  const MappedFaceFluxRate& rhs,
                                  const int rhs_comp)
{
    AMREX_ALWAYS_ASSERT(lhs.is_defined() && rhs.is_defined());
    AMREX_ALWAYS_ASSERT(lhs_comp >= 0 && lhs_comp < lhs.nComp());
    AMREX_ALWAYS_ASSERT(rhs_comp >= 0 && rhs_comp < rhs.nComp());
    amrex::Real maximum = amrex::Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& a = lhs.dir(dir);
        const auto& b = rhs.dir(dir);
        AMREX_ALWAYS_ASSERT(a.boxArray() == b.boxArray());
        AMREX_ALWAYS_ASSERT(a.DistributionMap() == b.DistributionMap());
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

void AccumulateIntegratedFaceFlux(IntegratedMappedFaceFlux& ledger,
                                 const MappedFaceFluxRate& rate,
                                 const amrex::Real weight)
{
    AMREX_ALWAYS_ASSERT(ledger.is_defined() && rate.is_defined());
    AMREX_ALWAYS_ASSERT(ledger.nComp() == rate.nComp());
    AMREX_ALWAYS_ASSERT(std::isfinite(weight) && weight >= amrex::Real(0.0));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Saxpy(ledger.dir(dir), weight, rate.dir(dir),
                               0, 0, ledger.nComp(), 0);
    }
}

} // namespace erf_auxiliary
