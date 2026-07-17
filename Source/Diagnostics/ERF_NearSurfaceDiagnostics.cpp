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
    return valid_real(value) && value >= Real(1.0) && value <= Real(6.0);
}

} // namespace

void fill (MultiFab& dst,
           int temperature_comp,
           int mixing_ratio_comp,
           int source_comp,
           const Sources& sources,
           Real missing_value)
{
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

        const bool have_native = sources.native_temperature_vegetated &&
                                 sources.native_temperature_bare &&
                                 sources.native_mixing_ratio_vegetated &&
                                 sources.native_mixing_ratio_bare &&
                                 sources.native_vegetation_fraction;
        const bool have_most = sources.theta_surface && sources.theta_star &&
                               sources.mixing_ratio_surface && sources.mixing_ratio_star &&
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
            Real source_value = missing_value;
            bool temperature_valid = false;

            if (have_native) {
                const NativeBundle bundle{
                    tv(i,j,0), tb(i,j,0), qv(i,j,0), qb(i,j,0), fv(i,j,0)
                };
                if (native_bundle_is_valid(bundle)) {
                    const NativeAggregate aggregate = aggregate_native_bundle(bundle);
                    temperature = aggregate.temperature;
                    mixing_ratio = aggregate.mixing_ratio;
                    temperature_valid = true;
                    source_value = static_cast<Real>(surface_diagnostics::to_int(
                        surface_diagnostics::SurfaceDiagnosticSource::LSMLand));
                }
            }

            if (!temperature_valid && have_most) {
                const MostState state{ts(i,j,0), tstar(i,j,0), qs(i,j,0),
                                      qstar(i,j,0), z0(i,j,0), olen(i,j,0)};
                MostProfile profile{};
                if (evaluate_most_profile(state, Real(2.0), profile) && have_atmosphere) {
                    const Real rho = cons(i,j,klo,Rho_comp);
                    const Real rhotheta = cons(i,j,klo,RhoTheta_comp);
                    const Real q_cell = moist ? cons(i,j,klo,RhoQ1_comp) / rho : zero;
                    if (valid_real(rho) && rho > zero && valid_real(rhotheta) &&
                        valid_real(q_cell) && q_cell >= zero) {
                        const Real p_cc = getPgivenRTh(rhotheta, q_cell);
                        const Real z_cc_agl = znd ? Compute_Zrel_AtCellCenter(i,j,klo,znd)
                                                   : Real(0.5) * dz;
                        const Real p_2m = p_cc + rho * CONST_GRAV * (z_cc_agl - Real(2.0));
                        const Real temperature_2m = (p_2m > zero && valid_real(p_2m))
                            ? getTgivenPandTh(p_2m, profile.potential_temperature, R_d/Cp_d)
                            : missing_value;
                        if (valid_real(temperature_2m) && temperature_2m > zero) {
                        temperature = temperature_2m;
                        temperature_valid = true;
                        mixing_ratio = (moist && profile.mixing_ratio >= zero)
                            ? profile.mixing_ratio : missing_value;

                        Real selected_source = missing_value;
                        if (have_source_mask && source_value_is_valid(source_mask(i,j,0))) {
                            selected_source = source_mask(i,j,0);
                            if (static_cast<int>(selected_source) ==
                                surface_diagnostics::to_int(surface_diagnostics::SurfaceDiagnosticSource::LSMLand)) {
                                selected_source = static_cast<Real>(surface_diagnostics::to_int(
                                    surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerFallback));
                            }
                        } else if (have_land_mask && lmask(i,j,0) == 0) {
                            selected_source = static_cast<Real>(surface_diagnostics::to_int(
                                surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerSea));
                        } else {
                            selected_source = static_cast<Real>(surface_diagnostics::to_int(
                                has_lsm ? surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerFallback
                                        : surface_diagnostics::SurfaceDiagnosticSource::SurfaceLayerLand));
                        }
                        source_value = selected_source;
                        }
                    }
                }
            }

            if (!temperature_valid) {
                source_value = missing_value;
            }
            if (temperature_comp >= 0) {
                dst_arr(i,j,k,temperature_comp) = temperature;
            }
            if (mixing_ratio_comp >= 0) {
                dst_arr(i,j,k,mixing_ratio_comp) = (temperature_valid && moist)
                    ? mixing_ratio : missing_value;
            }
            if (source_comp >= 0) {
                dst_arr(i,j,k,source_comp) = source_value;
            }
        });
    }
}

} // namespace near_surface_diagnostics
