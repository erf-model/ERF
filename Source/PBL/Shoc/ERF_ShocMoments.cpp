#include "ERF_ShocMoments.H"

#include "ERF_Constants.H"

#include <algorithm>
#include <cmath>

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_min_tke () noexcept { return 4.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_large_neg () noexcept { return -99999999.99_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_w3clip () noexcept { return 1.2_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_w3clipdef () noexcept { return 0.02_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_base_temp () noexcept { return 300.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_ufmin () noexcept { return 0.01_rt; }

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
    Real signed_denominator (Real value, Real eps = 1.0e-12)
    {
        if (amrex::Math::abs(value) >= eps) {
            return value;
        }
        return std::copysign(eps, value == 0.0_rt ? 1.0_rt : value);
    }

    void interpolate_cc_to_iface (const ShocColumnData& col,
                                  const FArrayBox& cc_in,
                                  FArrayBox& iface_out,
                                  Real min_thresh)
    {
        const auto zt = col.zt.const_array();
        const auto zi = col.zi.const_array();
        const auto in = cc_in.const_array();
        auto out = iface_out.array();

        const auto layout = col.layout;
        const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            if (layout.nlev == 1) {
                out(ic,0,0) = amrex::max(min_thresh, in(ic,0,0));
                out(ic,1,0) = amrex::max(min_thresh, in(ic,0,0));
                return;
            }

            out(ic,0,0) = amrex::max(min_thresh,
                                     weighted_linear_interp(zt(ic,0,0), zt(ic,1,0),
                                                            in(ic,0,0), in(ic,1,0),
                                                            zi(ic,0,0)));
            for (int k = 1; k < layout.nlev; ++k) {
                const Real interp = weighted_linear_interp(zt(ic,k-1,0), zt(ic,k,0),
                                                           in(ic,k-1,0), in(ic,k,0),
                                                           zi(ic,k,0));
                out(ic,k,0) = amrex::max(min_thresh, interp);
            }
            out(ic,layout.nlev,0) = amrex::max(min_thresh,
                                               weighted_linear_interp(zt(ic,layout.nlev-2,0),
                                                                      zt(ic,layout.nlev-1,0),
                                                                      in(ic,layout.nlev-2,0),
                                                                      in(ic,layout.nlev-1,0),
                                                                      zi(ic,layout.nlev,0)));
        });
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real cell_spacing_at_iface (const Array4<const Real>& zt,
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
    Real top_taper_from_distance (const ShocRuntimeOptions& opts,
                                  Real distance_from_top)
    {
        if (opts.top_taper_depth <= 0.0) {
            return 1.0;
        }

        const Real x = shoc_clamp(amrex::max(0.0_rt, distance_from_top) / opts.top_taper_depth,
                                  0.0_rt, 1.0_rt);
        const Real smooth = x * x * (3.0_rt - 2.0_rt * x);
        return opts.top_taper_min_factor + (1.0_rt - opts.top_taper_min_factor) * smooth;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real top_taper_cell (const Array4<const Real>& zt,
                         const Array4<const Real>& zi,
                         const ShocColumnLayout& layout,
                         const ShocRuntimeOptions& opts,
                         int ic,
                         int k)
    {
        return top_taper_from_distance(opts, zi(ic,layout.nlev,0) - zt(ic,k,0));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real top_taper_iface (const Array4<const Real>& zi,
                          const ShocColumnLayout& layout,
                          const ShocRuntimeOptions& opts,
                          int ic,
                          int k)
    {
        return top_taper_from_distance(opts, zi(ic,layout.nlev,0) - zi(ic,k,0));
    }

    void apply_second_moment_boundary_conditions (ShocColumnData& col)
    {
        auto thl_sec = col.thl_sec.array();
        auto qw_sec = col.qw_sec.array();
        auto qwthl_sec = col.qwthl_sec.array();
        auto wthl_sec = col.wthl_sec.array();
        auto wqw_sec = col.wqw_sec.array();
        auto uw_sec = col.uw_sec.array();
        auto vw_sec = col.vw_sec.array();
        auto wtke_sec = col.wtke_sec.array();

        const auto sflux = col.surf_sens_flux.const_array();
        const auto lflux = col.surf_lat_flux.const_array();
        const auto tauu = col.surf_tau_u.const_array();
        const auto tauv = col.surf_tau_v.const_array();

        const auto layout = col.layout;
        const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            const Real wthl_sfc = sflux(ic,0,0);
            const Real wqw_sfc = lflux(ic,0,0);
            const Real uw_sfc = tauu(ic,0,0);
            const Real vw_sfc = tauv(ic,0,0);
            const Real ustar2 = std::sqrt(uw_sfc * uw_sfc + vw_sfc * vw_sfc);
            const Real wstar = (wthl_sfc >= 0.0_rt)
                ? std::cbrt((CONST_GRAV / shoc_base_temp()) * wthl_sfc)
                : 0.0_rt;
            const Real uf = amrex::max(shoc_ufmin(),
                                       std::sqrt(ustar2 + 0.3_rt * wstar * wstar));
            const Real wthl_ratio = wthl_sfc / uf;
            const Real wqw_ratio = wqw_sfc / uf;
            const Real wtke_scale = amrex::max(std::sqrt(ustar2), shoc_ufmin());

            thl_sec(ic,0,0) = 0.72_rt * wthl_ratio * wthl_ratio;
            qw_sec(ic,0,0) = 0.72_rt * wqw_ratio * wqw_ratio;
            qwthl_sec(ic,0,0) = 0.36_rt * wthl_ratio * wqw_ratio;
            wthl_sec(ic,0,0) = wthl_sfc;
            wqw_sec(ic,0,0) = wqw_sfc;
            uw_sec(ic,0,0) = uw_sfc;
            vw_sec(ic,0,0) = vw_sfc;
            wtke_sec(ic,0,0) = wtke_scale * wtke_scale * wtke_scale;

            const int ktop = layout.nlev;
            thl_sec(ic,ktop,0) = 0.0_rt;
            qw_sec(ic,ktop,0) = 0.0_rt;
            qwthl_sec(ic,ktop,0) = 0.0_rt;
            wthl_sec(ic,layout.nlev,0) = 0.0_rt;
            wqw_sec(ic,layout.nlev,0) = 0.0_rt;
            uw_sec(ic,layout.nlev,0) = 0.0_rt;
            vw_sec(ic,layout.nlev,0) = 0.0_rt;
            wtke_sec(ic,layout.nlev,0) = 0.0_rt;
        });
    }

    void apply_top_taper_to_second_moments (ShocColumnData& col,
                                            const ShocRuntimeOptions& opts)
    {
        if (opts.top_taper_depth <= 0.0) {
            return;
        }

        auto thl_sec = col.thl_sec.array();
        auto qw_sec = col.qw_sec.array();
        auto qwthl_sec = col.qwthl_sec.array();
        auto wthl_sec = col.wthl_sec.array();
        auto wqw_sec = col.wqw_sec.array();
        auto uw_sec = col.uw_sec.array();
        auto vw_sec = col.vw_sec.array();
        auto wtke_sec = col.wtke_sec.array();
        auto w_sec = col.w_sec.array();

        const auto zt = col.zt.const_array();
        const auto zi = col.zi.const_array();
        const auto layout = col.layout;
        const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            for (int k = 0; k < layout.nlev; ++k) {
                w_sec(ic,k,0) *= top_taper_cell(zt, zi, layout, opts, ic, k);
            }
            for (int k = 0; k <= layout.nlev; ++k) {
                const Real factor = top_taper_iface(zi, layout, opts, ic, k);
                thl_sec(ic,k,0) *= factor;
                qw_sec(ic,k,0) *= factor;
                qwthl_sec(ic,k,0) *= factor;
                wthl_sec(ic,k,0) *= factor;
                wqw_sec(ic,k,0) *= factor;
                uw_sec(ic,k,0) *= factor;
                vw_sec(ic,k,0) *= factor;
                wtke_sec(ic,k,0) *= factor;
            }
        });
    }

    void apply_top_taper_to_third_moments (ShocColumnData& col,
                                           const ShocRuntimeOptions& opts)
    {
        if (opts.top_taper_depth <= 0.0) {
            return;
        }

        auto w3 = col.w3.array();
        const auto zi = col.zi.const_array();
        const auto layout = col.layout;
        const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            for (int k = 0; k <= layout.nlev; ++k) {
                w3(ic,k,0) *= top_taper_iface(zi, layout, opts, ic, k);
            }
        });
    }
}

void
ShocMoments::calc_var_or_covar (const ShocColumnData& col,
                                Real tunefac,
                                const FArrayBox& isotropy_zi,
                                const FArrayBox& tkh_zi,
                                const FArrayBox& invar1,
                                const FArrayBox& invar2,
                                FArrayBox& outvar)
{
    auto out = outvar.array();
    const auto iso = isotropy_zi.const_array();
    const auto tkh = tkh_zi.const_array();
    const auto v1 = invar1.const_array();
    const auto v2 = invar2.const_array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        for (int k = 1; k < layout.nlev; ++k) {
            const Real dz_zi = cell_spacing_at_iface(zt, zi, layout, ic, k);
            const Real grid_dz2 = 1.0_rt / (dz_zi * dz_zi);
            out(ic,k,0) = tunefac * (iso(ic,k,0) * tkh(ic,k,0)) * grid_dz2 *
                          (v1(ic,k-1,0) - v1(ic,k,0)) * (v2(ic,k-1,0) - v2(ic,k,0));
        }
    });
}

void
ShocMoments::calc_vertflux (const ShocColumnData& col,
                            const FArrayBox& tkh_zi,
                            const FArrayBox& invar,
                            FArrayBox& vertflux)
{
    auto out = vertflux.array();
    const auto tkh = tkh_zi.const_array();
    const auto v = invar.const_array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        for (int k = 1; k < layout.nlev; ++k) {
            const Real dz_zi = cell_spacing_at_iface(zt, zi, layout, ic, k);
            // E3SM uses top-down vertical indexing. ERF columns are bottom-up,
            // so the interface gradient must use upper-minus-lower to preserve
            // the same physical downgradient flux sign.
            out(ic,k,0) = -(tkh(ic,k,0) / dz_zi) * (v(ic,k,0) - v(ic,k-1,0));
        }
    });
}

void
ShocMoments::diagnose_second_moments (ShocColumnData& col,
                                      const ShocRuntimeOptions& opts)
{
    FArrayBox isotropy_zi, tkh_zi, tk_zi;
    const Box iface_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, col.layout.nlev, 0));
    isotropy_zi.resize(iface_box, 1, The_Async_Arena());
    tkh_zi.resize(iface_box, 1, The_Async_Arena());
    tk_zi.resize(iface_box, 1, The_Async_Arena());

    interpolate_cc_to_iface(col, col.isotropy, isotropy_zi, 0.0);
    interpolate_cc_to_iface(col, col.tkh, tkh_zi, 0.0);
    interpolate_cc_to_iface(col, col.tk, tk_zi, 0.0);

    auto thl_sec = col.thl_sec.array();
    auto qw_sec = col.qw_sec.array();
    auto qwthl_sec = col.qwthl_sec.array();
    auto w_sec = col.w_sec.array();

    const auto tke = col.tke.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        for (int k = 0; k < layout.nlev; ++k) {
            w_sec(ic,k,0) = opts.shoc_1p5tke ? 0.0_rt : opts.w2tune * (2.0_rt / 3.0_rt) * tke(ic,k,0);
        }
    });

    if (opts.shoc_1p5tke) {
        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            for (int k = 0; k <= layout.nlev; ++k) {
                thl_sec(ic,k,0) = 0.0_rt;
                qw_sec(ic,k,0) = 0.0_rt;
                qwthl_sec(ic,k,0) = 0.0_rt;
            }
        });
    } else {
        calc_var_or_covar(col, opts.thl2tune, isotropy_zi, tkh_zi, col.thetal, col.thetal, col.thl_sec);
        calc_var_or_covar(col, opts.qw2tune, isotropy_zi, tkh_zi, col.qw, col.qw, col.qw_sec);
        calc_var_or_covar(col, opts.qwthl2tune, isotropy_zi, tkh_zi, col.thetal, col.qw, col.qwthl_sec);
    }

    calc_vertflux(col, tkh_zi, col.thetal, col.wthl_sec);
    calc_vertflux(col, tkh_zi, col.qw, col.wqw_sec);
    calc_vertflux(col, tkh_zi, col.tke, col.wtke_sec);
    calc_vertflux(col, tk_zi, col.u, col.uw_sec);
    calc_vertflux(col, tk_zi, col.v, col.vw_sec);

    apply_second_moment_boundary_conditions(col);
    apply_top_taper_to_second_moments(col, opts);
}

void
ShocMoments::clip_third_moments (const ShocColumnData& col,
                                 const FArrayBox& w_sec_zi,
                                 FArrayBox& w3_fab)
{
    auto w3 = w3_fab.array();
    const auto wsec = w_sec_zi.const_array();

    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        for (int k = 0; k <= layout.nlev; ++k) {
            const Real wsec3 = wsec(ic,k,0) * wsec(ic,k,0) * wsec(ic,k,0);
            const Real clip_cond = shoc_w3clip() * std::sqrt(amrex::max(0.0_rt, 2.0_rt * wsec3));
            if (amrex::Math::abs(w3(ic,k,0)) > clip_cond) {
                w3(ic,k,0) = shoc_w3clipdef();
            }
        }
    });
}

void
ShocMoments::diagnose_third_moments (ShocColumnData& col,
                                     const ShocRuntimeOptions& opts)
{
    FArrayBox isotropy_zi, brunt_zi, w_sec_zi, thetal_zi;
    const Box iface_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, col.layout.nlev, 0));
    isotropy_zi.resize(iface_box, 1, The_Async_Arena());
    brunt_zi.resize(iface_box, 1, The_Async_Arena());
    w_sec_zi.resize(iface_box, 1, The_Async_Arena());
    thetal_zi.resize(iface_box, 1, The_Async_Arena());

    interpolate_cc_to_iface(col, col.isotropy, isotropy_zi, 0.0);
    interpolate_cc_to_iface(col, col.brunt, brunt_zi, shoc_large_neg());
    interpolate_cc_to_iface(col, col.w_sec, w_sec_zi, (2.0_rt / 3.0_rt) * shoc_min_tke());
    interpolate_cc_to_iface(col, col.thetal, thetal_zi, 0.0);

    auto w3 = col.w3.array();
    const auto w_sec = col.w_sec.const_array();
    const auto thl_sec = col.thl_sec.const_array();
    const auto wthl_sec = col.wthl_sec.const_array();
    const auto tke = col.tke.const_array();
    const auto dz = col.dz.const_array();
    const auto isotropy_i = isotropy_zi.const_array();
    const auto brunt_i = brunt_zi.const_array();
    const auto w_sec_i = w_sec_zi.const_array();
    const auto thetal_i = thetal_zi.const_array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        w3(ic,0,0) = 0.0_rt;
        w3(ic,layout.nlev,0) = 0.0_rt;

        for (int k = 1; k < layout.nlev; ++k) {
            if (opts.shoc_1p5tke) {
                w3(ic,k,0) = 0.0_rt;
                continue;
            }

            const Real dz_zt_k = amrex::max(dz(ic,k,0), 1.0e-12_rt);
            const Real dz_zt_km1 = amrex::max(dz(ic,k-1,0), 1.0e-12_rt);
            const Real dz_zi = cell_spacing_at_iface(zt, zi, layout, ic, k);
            const Real thedz = 1.0_rt / dz_zi;
            const Real thedz2 = 1.0_rt / (dz_zt_k + dz_zt_km1);
            const Real iso = isotropy_i(ic,k,0);
            const Real isosq = iso * iso;
            const Real buoy_sgs2 = isosq * brunt_i(ic,k,0);
            const Real bet2 = CONST_GRAV / amrex::max(thetal_i(ic,k,0), 1.0e-12_rt);

            // E3SM's top-down indexing forms above-minus-below centered
            // differences here. In ERF's bottom-up ordering, the upper
            // neighbor is k+1 and the lower neighbor is k-1, while the cell
            // just above interface k is k and the one below is k-1.
            const Real thl_sec_diff = thl_sec(ic,amrex::min(k+1, layout.nlev),0) - thl_sec(ic,k-1,0);
            const Real wthl_sec_diff = wthl_sec(ic,amrex::min(k+1, layout.nlev),0) - wthl_sec(ic,k-1,0);
            const Real wsec_diff = w_sec(ic,k,0) - w_sec(ic,k-1,0);
            const Real tke_diff = tke(ic,k,0) - tke(ic,k-1,0);

            const Real c = opts.c_diag_3rd_mom;
            const Real a0 = (0.52_rt / (c * c)) / amrex::max(c - 2.0_rt, 1.0e-12_rt);
            const Real a1 = 0.87_rt / (c * c);
            const Real a2 = 0.5_rt / c;
            const Real a3 = 0.6_rt / (c * amrex::max(c - 2.0_rt, 1.0e-12_rt));
            const Real a4 = 2.4_rt / (3.0_rt * c + 5.0_rt);
            const Real a5 = 0.6_rt / (c * (3.0_rt + 5.0_rt * c));
            const Real bet2_sq = bet2 * bet2;
            const Real bet2_cu = bet2_sq * bet2;
            const Real iso_cu = isosq * iso;

            const Real f0 = thedz2 * bet2_cu * isosq * isosq *
                            wthl_sec(ic,k,0) * thl_sec_diff;
            const Real f1 = thedz2 * bet2_sq * iso_cu *
                            (wthl_sec(ic,k,0) * wthl_sec_diff +
                            0.5_rt * w_sec_i(ic,k,0) * thl_sec_diff);
            const Real f2 = thedz * bet2 * isosq * wthl_sec(ic,k,0) * wsec_diff +
                            2.0_rt * thedz2 * bet2 * isosq * w_sec_i(ic,k,0) * wthl_sec_diff;
            const Real f3 = thedz2 * bet2 * isosq * w_sec_i(ic,k,0) * wthl_sec_diff +
                            thedz * bet2 * isosq * (wthl_sec(ic,k,0) * tke_diff);
            const Real f4 = thedz * iso * w_sec_i(ic,k,0) * (wsec_diff + tke_diff);
            const Real f5 = thedz * iso * w_sec_i(ic,k,0) * wsec_diff;

            const Real denom0 = signed_denominator(1.0_rt - (a1 + a3) * buoy_sgs2);
            const Real omega0 = a4 / signed_denominator(1.0_rt - a5 * buoy_sgs2);
            const Real omega1 = omega0 / (2.0_rt * c);
            const Real omega2 = omega1 * f3 + 1.25_rt * omega0 * f4;
            const Real x0 = (a2 * buoy_sgs2 * (1.0_rt - a3 * buoy_sgs2)) / denom0;
            const Real y0 = (2.0_rt * a2 * buoy_sgs2 * x0) / signed_denominator(1.0_rt - a3 * buoy_sgs2);
            const Real x1 = (a0 * f0 + a1 * f1 + a2 * (1.0_rt - a3 * buoy_sgs2) * f2) / denom0;
            const Real y1 = (2.0_rt * a2 * (buoy_sgs2 * x1 + (a0 / amrex::max(a1, 1.0e-12_rt)) * f0 + f1)) /
                            signed_denominator(1.0_rt - a3 * buoy_sgs2);
            const Real aa0 = omega0 * x0 + omega1 * y0;
            const Real aa1 = omega0 * x1 + omega1 * y1 + omega2;

            const Real denom = signed_denominator(c - 1.2_rt * x0 + aa0);
            Real w3_val = (aa1 - 1.2_rt * x1 - 1.5_rt * f5) / denom;
            w3(ic,k,0) = w3_val;
        }
    });

    clip_third_moments(col, w_sec_zi, col.w3);
    apply_top_taper_to_third_moments(col, opts);
}

void
ShocMoments::diagnose_moments (ShocColumnData& col,
                               const ShocRuntimeOptions& opts)
{
    diagnose_second_moments(col, opts);
    diagnose_third_moments(col, opts);
}
