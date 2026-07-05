#include "ERF_Plotfile2DInterpolator.H"
#include "ERF_Plotfile2DFill.H"

#include <cmath>
#include <sstream>

using namespace amrex;

namespace plotfile2d
{

namespace
{

// A target is bracketed when it lies between adjacent coordinate values. This
// test works for increasing height and decreasing pressure.
SampledBracket find_sampled_bracket (SampledCoordinate coordinate,
                                     Real target,
                                     int klo,
                                     int khi,
                                     int i, int j,
                                     const Array4<const Real>& cons_arr,
                                     const Array4<const Real>& z_phys_cc_arr,
                                     const Array4<const Real>& z_phys_nd_arr,
                                     bool have_z_phys_cc,
                                     const MoistureComponentIndices& moisture_indices) noexcept
{
    SampledBracket bracket;
    if (khi < klo) {
        return bracket;
    }

    auto coordinate_value = [=] AMREX_GPU_HOST_DEVICE (int kk) noexcept -> Real {
        if (coordinate == SampledCoordinate::ModelIndex) {
            return static_cast<Real>(kk);
        }

        const auto field_id = (coordinate == SampledCoordinate::HeightMSL)
            ? SampledFieldID::HeightMSL
            : (coordinate == SampledCoordinate::HeightAGL)
                ? SampledFieldID::HeightAGL
                : SampledFieldID::Pressure;
        return sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                   have_z_phys_cc, i, j, kk, moisture_indices);
    };

    const Real c0 = coordinate_value(klo);
    if (target == c0) {
        bracket.klo = klo;
        bracket.khi = klo;
        bracket.found = true;
        return bracket;
    }

    Real lo = c0;
    for (int k = klo; k < khi; ++k) {
        const Real hi = coordinate_value(k + 1);
        if (sampled_target_is_bracketed(target, lo, hi)) {
            bracket.klo = k;
            bracket.khi = k + 1;
            bracket.found = true;
            return bracket;
        }
        lo = hi;
    }

    return bracket;
}

} // namespace

void
fill_sampled_level_component (MultiFab& dst,
                              int dst_comp,
                              const Plotfile2DOutputDescriptor& descriptor,
                              const MultiFab& cons,
                              const MultiFab* z_phys_cc,
                              const MultiFab& z_phys_nd,
                              bool have_z_phys_cc,
                              const MoistureComponentIndices& moisture_indices,
                              int klo,
                              int khi)
{
    const auto* sampled = descriptor.sampled_level ? &*descriptor.sampled_level : nullptr;
    if (sampled == nullptr) {
        fill_component_with_value(dst, dst_comp, descriptor.missing_value);
        return;
    }

    const auto* field_descriptor = find_sampled_field(sampled->source_field);
    if (field_descriptor == nullptr) {
        Abort("Unknown sampled-level source field '" + sampled->source_field + "'");
    }

    const auto coordinate = sampled_coordinate_from_string(sampled->vertical_coordinate.type);
    const auto field_id = field_descriptor->id;
    const Real target = sampled->vertical_coordinate.canonical_value;
    const int target_k = static_cast<int>(std::nearbyint(target));
    const Real missing_value = descriptor.missing_value;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& dst_arr = dst.array(mfi);
        const auto& cons_arr = cons.const_array(mfi);
        const auto& z_phys_nd_arr = z_phys_nd.const_array(mfi);
        const auto& z_phys_cc_arr = have_z_phys_cc && z_phys_cc != nullptr
            ? z_phys_cc->const_array(mfi)
            : cons.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (coordinate == SampledCoordinate::ModelIndex) {
                if (target_k < klo || target_k > khi) {
                    dst_arr(i, j, k, dst_comp) = missing_value;
                    return;
                }

                dst_arr(i, j, k, dst_comp) =
                    sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                        have_z_phys_cc, i, j, target_k, moisture_indices);
                return;
            }

            const auto bracket = find_sampled_bracket(coordinate, target, klo, khi, i, j,
                                                      cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                                      have_z_phys_cc, moisture_indices);
            if (!bracket.found) {
                dst_arr(i, j, k, dst_comp) = missing_value;
                return;
            }

            const Real field_lo =
                sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                    have_z_phys_cc, i, j, bracket.klo, moisture_indices);
            if (bracket.klo == bracket.khi) {
                dst_arr(i, j, k, dst_comp) = field_lo;
                return;
            }

            const Real field_hi =
                sampled_field_value(field_id, cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                    have_z_phys_cc, i, j, bracket.khi, moisture_indices);
            const Real coord_lo =
                sampled_field_value((coordinate == SampledCoordinate::HeightMSL)
                                        ? SampledFieldID::HeightMSL
                                        : (coordinate == SampledCoordinate::HeightAGL)
                                            ? SampledFieldID::HeightAGL
                                            : SampledFieldID::Pressure,
                                    cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                    have_z_phys_cc, i, j, bracket.klo, moisture_indices);
            const Real coord_hi =
                sampled_field_value((coordinate == SampledCoordinate::HeightMSL)
                                        ? SampledFieldID::HeightMSL
                                        : (coordinate == SampledCoordinate::HeightAGL)
                                            ? SampledFieldID::HeightAGL
                                            : SampledFieldID::Pressure,
                                    cons_arr, z_phys_cc_arr, z_phys_nd_arr,
                                    have_z_phys_cc, i, j, bracket.khi, moisture_indices);
            dst_arr(i, j, k, dst_comp) =
                linear_interpolate(field_lo, field_hi, coord_lo, coord_hi, target);
        });
    }
}

} // namespace plotfile2d
