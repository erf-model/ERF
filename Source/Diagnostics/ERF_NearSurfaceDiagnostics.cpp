#include "Diagnostics/ERF_NearSurfaceDiagnostics.H"

#include <AMReX_Gpu.H>

#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_MOSTStress.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_IndexDefines.H"

using namespace amrex;

namespace near_surface_diagnostics
{

namespace
{

Array4<const Real> array_or_empty (const MultiFab* mf, const MFIter& mfi)
{
    return mf ? mf->const_array(mfi) : Array4<const Real>{};
}

Array4<const int> int_array_or_empty (const iMultiFab* mf, const MFIter& mfi)
{
    return mf ? mf->const_array(mfi) : Array4<const int>{};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool source_value_is_valid (Real value) noexcept
{
    if (!valid_real(value) || value < Real(1.0) || value > Real(6.0)) {
        return false;
    }
    const int code = static_cast<int>(value);
    return value == static_cast<Real>(code);
}

} // namespace

void fill (MultiFab& dst,
           int temperature_comp,
           int mixing_ratio_comp,
           int source_comp,
           const Sources& sources,
           Real missing_value)
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto dst_arr = dst.array(mfi);

        const auto tv = array_or_empty(sources.native_temperature_vegetated, mfi);
        const auto tb = array_or_empty(sources.native_temperature_bare, mfi);
        const auto qv = array_or_empty(sources.native_mixing_ratio_vegetated, mfi);
        const auto qb = array_or_empty(sources.native_mixing_ratio_bare, mfi);
        const auto fv = array_or_empty(sources.native_vegetation_fraction, mfi);
        const auto ts = array_or_empty(sources.theta_surface, mfi);
        const auto tstar = array_or_empty(sources.theta_star, mfi);
        const auto qs = array_or_empty(sources.mixing_ratio_surface, mfi);
        const auto qstar = array_or_empty(sources.mixing_ratio_star, mfi);
        const auto z0 = array_or_empty(sources.roughness_height, mfi);
        const auto olen = array_or_empty(sources.obukhov_length, mfi);
        const auto source_mask = array_or_empty(sources.source_mask, mfi);
        const auto lmask = int_array_or_empty(sources.land_mask, mfi);
        const auto cons = array_or_empty(sources.cons, mfi);
        const auto znd = array_or_empty(sources.z_phys_nd, mfi);

        const bool request_temperature = temperature_comp >= 0;
        const bool request_mixing_ratio = mixing_ratio_comp >= 0;
        const bool request_source = source_comp >= 0;
        // Provenance follows the temperature path only for a source-only
        // request. When paired with humidity, it describes the humidity
        // source without adding an unrequested temperature prerequisite.
        const bool source_only = request_source && !request_temperature &&
                                 !request_mixing_ratio;
        const bool need_temperature = request_temperature || source_only;
        const bool need_mixing_ratio = request_mixing_ratio;
        const bool have_native_temperature = sources.native_temperature_vegetated &&
                                              sources.native_temperature_bare &&
                                              sources.native_vegetation_fraction;
        const bool have_native_mixing_ratio = sources.native_mixing_ratio_vegetated &&
                                               sources.native_mixing_ratio_bare &&
                                               sources.native_vegetation_fraction;
        const bool have_most_temperature = sources.theta_surface && sources.theta_star &&
                                           sources.roughness_height && sources.obukhov_length;
        const bool have_most_mixing_ratio = sources.mixing_ratio_surface &&
                                            sources.mixing_ratio_star &&
                                            sources.roughness_height && sources.obukhov_length;
        const bool have_atmosphere = sources.cons != nullptr;
        const bool have_source_mask = sources.source_mask != nullptr;
        const bool have_land_mask = sources.land_mask != nullptr;
        const int klo = sources.klo;
        const Real dz = sources.dz;
        const bool moist = sources.moist;
        const bool has_lsm = sources.has_lsm;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real temperature = missing_value;
            Real mixing_ratio = missing_value;
            Real source_value = surface_diagnostics::to_plot_value(
                surface_diagnostics::SurfaceDiagnosticSource::Missing);

            const NativeBundle native_bundle{
                have_native_temperature ? tv(i,j,0) : missing_value,
                have_native_temperature ? tb(i,j,0) : missing_value,
                have_native_mixing_ratio ? qv(i,j,0) : missing_value,
                have_native_mixing_ratio ? qb(i,j,0) : missing_value,
                (have_native_temperature || have_native_mixing_ratio) ? fv(i,j,0) : missing_value
            };
            const bool native_temperature_valid = have_native_temperature &&
                native_bundle_temperature_is_valid(native_bundle);
            const bool native_mixing_ratio_valid = moist && have_native_mixing_ratio &&
                native_bundle_mixing_ratio_is_valid(native_bundle);
            const bool native_selected = (!need_temperature || native_temperature_valid) &&
                                         (!need_mixing_ratio || native_mixing_ratio_valid);

            bool selected = false;
            if (native_selected && (need_temperature || need_mixing_ratio)) {
                const NativeAggregate aggregate = aggregate_native_bundle(native_bundle);
                if (need_temperature) {
                    temperature = aggregate.temperature;
                }
                if (need_mixing_ratio) {
                    mixing_ratio = aggregate.mixing_ratio;
                }
                source_value = surface_diagnostics::to_plot_value(
                    surface_diagnostics::SurfaceDiagnosticSource::LSMLand);
                selected = true;
            }

            if (!selected && (need_temperature || need_mixing_ratio)) {
                const MostState state{
                    have_most_temperature ? ts(i,j,0) : missing_value,
                    have_most_temperature ? tstar(i,j,0) : missing_value,
                    have_most_mixing_ratio ? qs(i,j,0) : missing_value,
                    have_most_mixing_ratio ? qstar(i,j,0) : missing_value,
                    (have_most_temperature || have_most_mixing_ratio) ? z0(i,j,0) : missing_value,
                    (have_most_temperature || have_most_mixing_ratio) ? olen(i,j,0) : missing_value
                };
                Real factor = zero;
                const bool common_most_valid = evaluate_most_factor(state, Real(2.0), factor);
                bool most_temperature_valid = false;
                bool most_mixing_ratio_valid = false;
                Real theta_2m = missing_value;
                Real q_2m = missing_value;
                if (common_most_valid && need_temperature && have_most_temperature) {
                    theta_2m = state.theta_surface + state.theta_star / KAPPA * factor;
                    if (valid_real(theta_2m) && theta_2m > zero && have_atmosphere) {
                        const Real rho = cons(i,j,klo,Rho_comp);
                        const Real rhotheta = cons(i,j,klo,RhoTheta_comp);
                        const Real q_cell = moist ? cons(i,j,klo,RhoQ1_comp) / rho : zero;
                        if (valid_real(rho) && rho > zero && valid_real(rhotheta) &&
                            valid_real(q_cell) && q_cell >= zero) {
                            const Real p_cc = getPgivenRTh(rhotheta, q_cell);
                            const Real z_cc_agl = znd ? Compute_Zrel_AtCellCenter(i,j,klo,znd)
                                                       : Real(0.5) * dz;
                            const Real p_2m = p_cc + rho * CONST_GRAV *
                                              (z_cc_agl - Real(2.0));
                            const Real t_2m = (valid_real(p_2m) && p_2m > zero)
                                ? getTgivenPandTh(p_2m, theta_2m, R_d/Cp_d)
                                : missing_value;
                            most_temperature_valid = valid_real(t_2m) && t_2m > zero;
                            if (most_temperature_valid) { temperature = t_2m; }
                        }
                    }
                }
                if (common_most_valid && need_mixing_ratio && moist &&
                    have_most_mixing_ratio && valid_real(state.mixing_ratio_surface) &&
                    valid_real(state.mixing_ratio_star)) {
                    q_2m = state.mixing_ratio_surface + state.mixing_ratio_star / KAPPA * factor;
                    most_mixing_ratio_valid = valid_real(q_2m) && q_2m >= zero;
                    if (most_mixing_ratio_valid) {
                        mixing_ratio = q_2m;
                    }
                }
                const bool most_selected = (!need_temperature || most_temperature_valid) &&
                                           (!need_mixing_ratio || most_mixing_ratio_valid);
                // Do not expose a valid scalar from an incomplete MOST bundle
                // when both requested fields must come from one source.
                if (!most_selected) {
                    temperature = missing_value;
                    mixing_ratio = missing_value;
                }
                if (most_selected) {
                    if (have_source_mask && source_value_is_valid(source_mask(i,j,0))) {
                        source_value = source_mask(i,j,0);
                        if (static_cast<int>(source_value) ==
                            surface_diagnostics::to_int(surface_diagnostics::SurfaceDiagnosticSource::LSMLand)) {
                            source_value = surface_diagnostics::to_plot_value(
                                surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerFallback);
                        }
                    } else if (have_land_mask && lmask(i,j,0) == 0) {
                        source_value = surface_diagnostics::to_plot_value(
                            surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerSea);
                    } else {
                        source_value = surface_diagnostics::to_plot_value(
                            has_lsm ? surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerFallback
                                    : surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerLand);
                    }
                    selected = true;
                }
            }

            if (temperature_comp >= 0) {
                dst_arr(i,j,k,temperature_comp) = temperature;
            }
            if (mixing_ratio_comp >= 0) {
                dst_arr(i,j,k,mixing_ratio_comp) = (selected && moist)
                    ? mixing_ratio : missing_value;
            }
            if (source_comp >= 0) {
                dst_arr(i,j,k,source_comp) = source_value;
            }
        });
    }
}

} // namespace near_surface_diagnostics
