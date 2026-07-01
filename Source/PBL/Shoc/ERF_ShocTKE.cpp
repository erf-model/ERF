#include "ERF_ShocTKE.H"

#include "ERF_Constants.H"
#include "ERF_ShocGpuUtils.H"

#include <algorithm>
#include <cmath>

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_base_temp () noexcept { return 300.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_min_tke () noexcept { return 4.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_max_tke () noexcept { return 50.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_max_iso () noexcept { return 2.0e4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_tabs_crit () noexcept { return 182.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pbl_trans () noexcept { return 200.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_stable_ckh () noexcept { return 0.1_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_stable_ckm () noexcept { return 0.1_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_shear_ck () noexcept { return 0.1_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_tke_cs () noexcept { return 0.15_rt; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real shoc_tke_ce () noexcept
    {
        return (shoc_shear_ck() * shoc_shear_ck() * shoc_shear_ck()) /
               ((shoc_tke_cs() * shoc_tke_cs()) * (shoc_tke_cs() * shoc_tke_cs()));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real shoc_tke_cee () noexcept
    {
        return (shoc_tke_ce() / 0.7_rt) * 0.19_rt + (shoc_tke_ce() / 0.7_rt) * 0.51_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_trop_pres () noexcept { return 8.0e4_rt; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real weighted_linear_interp (Real x0, Real x1, Real y0, Real y1, Real x)
    {
        const Real denom = x1 - x0;
        if (amrex::Math::abs(denom) <= 1.0e-12_rt) {
            return 0.5_rt * (y0 + y1);
        }
        return y0 + (y1 - y0) * (x - x0) / denom;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real interface_spacing (const Array4<const Real>& zt,
                            const Array4<const Real>& zi,
                            const ShocColumnLayout& layout,
                            int ic, int k)
    {
        if (k <= 0) {
            return amrex::max(zt(ic,0,0) - zi(ic,0,0), 1.0e-12_rt);
        }
        if (k >= layout.nlev) {
            return amrex::max(zi(ic,layout.nlev,0) - zt(ic,layout.nlev-1,0), 1.0e-12_rt);
        }
        return amrex::max(zt(ic,k,0) - zt(ic,k-1,0), 1.0e-12_rt);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real top_taper_factor (const Array4<const Real>& zt,
                           const Array4<const Real>& zi,
                           const ShocColumnLayout& layout,
                           const ShocRuntimeOptions& opts,
                           int ic, int k)
    {
        if (opts.top_taper_depth <= 0.0) {
            return 1.0;
        }

        const Real ztop = zi(ic,layout.nlev,0);
        const Real distance_from_top = amrex::max(0.0_rt, ztop - zt(ic,k,0));
        const Real x = shoc_clamp(distance_from_top / opts.top_taper_depth, 0.0_rt, 1.0_rt);
        const Real smooth = x * x * (3.0_rt - 2.0_rt * x);
        return opts.top_taper_min_factor + (1.0_rt - opts.top_taper_min_factor) * smooth;
    }

    void interpolate_iface_to_cc (const ShocColumnData& col,
                                  const FArrayBox& iface_in,
                                  FArrayBox& cc_out,
                                  Real min_thresh)
    {
        const auto zt = col.zt.const_array();
        const auto zi = col.zi.const_array();
        const auto in = iface_in.const_array();
        auto out = cc_out.array();
        const auto layout = col.layout;
        const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));

        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            for (int k = 0; k < layout.nlev; ++k) {
                const Real interp = weighted_linear_interp(zi(ic,k,0), zi(ic,k+1,0),
                                                           in(ic,k,0), in(ic,k+1,0),
                                                           zt(ic,k,0));
                out(ic,k,0) = amrex::max(min_thresh, interp);
            }
        });
    }
}

void
ShocTKE::compute_shear_production (const ShocColumnData& col,
                                   FArrayBox& sterm_iface)
{
    const Box iface_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, col.layout.nlev, 0));
    sterm_iface.resize(iface_box, 1, The_Async_Arena());
    shoc::set_fab_val(sterm_iface, 0.0);

    auto sterm = sterm_iface.array();
    const auto u = col.u.const_array();
    const auto v = col.v.const_array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));

    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        sterm(ic,0,0) = 0.0_rt;
        sterm(ic,layout.nlev,0) = 0.0_rt;
        for (int k = 1; k < layout.nlev; ++k) {
            const Real dz_zi = interface_spacing(zt, zi, layout, ic, k);
            const Real du_dz = (u(ic,k-1,0) - u(ic,k,0)) / dz_zi;
            const Real dv_dz = (v(ic,k-1,0) - v(ic,k,0)) / dz_zi;
            sterm(ic,k,0) = shoc_shear_ck() * (du_dz * du_dz + dv_dz * dv_dz);
        }
    });
}

void
ShocTKE::integrate_column_stability (const ShocColumnData& col,
                                     Vector<Real>& brunt_int)
{
    brunt_int.assign(col.layout.ncell, 0.0);
    const auto dz = col.dz.const_array();
    const auto p_mid = col.p_mid.const_array();
    const auto brunt = col.brunt.const_array();

    for (int ic = 0; ic < col.layout.ncell; ++ic) {
        for (int k = 0; k < col.layout.nlev; ++k) {
            if (p_mid(ic,k,0) > shoc_trop_pres()) {
                brunt_int[ic] += dz(ic,k,0) * brunt(ic,k,0);
            }
        }
    }
}

void
ShocTKE::diagnose_tke_and_diffusivities (ShocColumnData& col,
                                         const ShocRuntimeOptions& opts,
                                         Real dt)
{
    AMREX_ALWAYS_ASSERT(dt > 0.0);

    const Box cell_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, col.layout.nlev - 1, 0));
    FArrayBox sterm_iface;
    FArrayBox sterm_zt(cell_box, 1, The_Async_Arena());
    FArrayBox a_diss(cell_box, 1, The_Async_Arena());

    compute_shear_production(col, sterm_iface);
    interpolate_iface_to_cc(col, sterm_iface, sterm_zt, 0.0);

    auto tke = col.tke.array();
    auto isotropy = col.isotropy.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    auto tke_tend = col.tke_tend.array();
    auto shear_prod_arr = col.shear_prod.array();
    auto buoy_prod_arr = col.buoy_prod.array();
    auto diss_arr = col.diss_tke.array();

    const auto wthv_sec = col.wthv_sec.const_array();
    const auto shoc_mix = col.shoc_mix.const_array();
    const auto brunt = col.brunt.const_array();
    const auto zt = col.zt.const_array();
    const auto tabs = col.tabs.const_array();
    const auto pblh = col.pblh.const_array();
    const auto sterm = sterm_zt.const_array();
    auto a_diss_arr = a_diss.array();
    const auto dz = col.dz.const_array();
    const auto p_mid = col.p_mid.const_array();
    const auto zi = col.zi.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));

    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        const Real z_sfc = zi(ic,0,0);
        Real brunt_int = 0.0_rt;
        for (int k = 0; k < layout.nlev; ++k) {
            if (p_mid(ic,k,0) > shoc_trop_pres()) {
                brunt_int += dz(ic,k,0) * brunt(ic,k,0);
            }

            const Real old_tke = amrex::max(0.0_rt, tke(ic,k,0));
            const Real buoy_prod = opts.shoc_1p5tke
                ? -old_tke * brunt(ic,k,0)
                : (CONST_GRAV / shoc_base_temp()) * wthv_sec(ic,k,0);
            const Real shear_prod = tk(ic,k,0) * sterm(ic,k,0);
            const Real mix = amrex::max(shoc_mix(ic,k,0), 1.0e-12_rt);
            const Real diss = shoc_tke_cee() / mix * old_tke * std::sqrt(old_tke);
            const Real net_prod = opts.signed_tke_production
                ? (shear_prod + buoy_prod)
                : amrex::max(0.0_rt, shear_prod + buoy_prod);
            const Real raw_new_tke = shoc_clamp(old_tke + dt * (net_prod - diss),
                                                shoc_min_tke(), shoc_max_tke());
            const Real top_factor = top_taper_factor(zt, zi, layout, opts, ic, k);
            const Real new_tke = shoc_clamp(shoc_min_tke() +
                                            top_factor * (raw_new_tke - shoc_min_tke()),
                                            shoc_min_tke(), shoc_max_tke());

            a_diss_arr(ic,k,0) = top_factor * diss;
            shear_prod_arr(ic,k,0) = top_factor * shear_prod;
            buoy_prod_arr(ic,k,0) = top_factor * buoy_prod;
            diss_arr(ic,k,0) = top_factor * diss;
            tke_tend(ic,k,0) = new_tke - old_tke;
            tke(ic,k,0) = new_tke;
        }

        Real lambda = opts.lambda_low + ((brunt_int / CONST_GRAV) - opts.lambda_thresh) * opts.lambda_slope;
        lambda = shoc_clamp(lambda, opts.lambda_low, opts.lambda_high);

        for (int k = 0; k < layout.nlev; ++k) {
            const Real diss = amrex::max(a_diss_arr(ic,k,0), 1.0e-12_rt);
            const Real tscale = 2.0_rt * tke(ic,k,0) / diss;
            const Real local_lambda = (brunt(ic,k,0) <= 0.0_rt) ? 0.0_rt : lambda;
            isotropy(ic,k,0) = amrex::min(shoc_max_iso(),
                                          tscale / (1.0_rt + local_lambda * brunt(ic,k,0) * tscale * tscale));

            const Real zt_agl = shoc::height_agl(zt(ic,k,0), z_sfc);
            const bool use_stable_mix = (zt_agl < pblh(ic,0,0) + shoc_pbl_trans()) &&
                                        (tabs(ic,0,0) < shoc_tabs_crit());
            if (use_stable_mix) {
                const Real sterm_sqrt = std::sqrt(amrex::max(sterm(ic,k,0), 0.0_rt));
                tkh(ic,k,0) = shoc_stable_ckh() * shoc_mix(ic,k,0) * shoc_mix(ic,k,0) * sterm_sqrt;
                tk(ic,k,0) = shoc_stable_ckm() * shoc_mix(ic,k,0) * shoc_mix(ic,k,0) * sterm_sqrt;
            } else {
                tkh(ic,k,0) = opts.coeff_kh * isotropy(ic,k,0) * tke(ic,k,0);
                tk(ic,k,0) = opts.coeff_km * isotropy(ic,k,0) * tke(ic,k,0);
            }

            const Real top_factor = top_taper_factor(zt, zi, layout, opts, ic, k);
            isotropy(ic,k,0) *= top_factor;
            tkh(ic,k,0) = top_factor * amrex::max(0.0_rt, tkh(ic,k,0));
            tk(ic,k,0) = top_factor * amrex::max(0.0_rt, tk(ic,k,0));
        }
    });
}
