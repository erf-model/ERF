#include "ERF_ShocStructure.H"

#include "ERF_Constants.H"
#include "ERF_ShocThermoUtils.H"

#include <algorithm>
#include <cmath>

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_zvir () noexcept { return 0.61_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_u_star_min () noexcept { return 0.01_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_min_tke () noexcept { return 4.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_min_len () noexcept { return 20.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_max_len () noexcept { return 2.0e4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pblmaxp () noexcept { return 4.0e4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pbl_ricr () noexcept { return 0.3_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pbl_fac () noexcept { return 100.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pbl_fak () noexcept { return 8.5_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pbl_betam () noexcept { return 15.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_pbl_sffrac () noexcept { return 0.1_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_kbfs_eps () noexcept { return 1.0e-10_rt; }

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
    int diagnose_npbl (const Array4<const Real>& p_mid,
                       const ShocColumnLayout& layout,
                       int ic)
    {
        int npbl = 1;
        for (int k = 0; k < layout.nlev; ++k) {
            if (p_mid(ic,k,0) >= shoc_pblmaxp()) {
                npbl = k + 1;
            }
        }
        return npbl;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real theta_from_shoc_state (Real thetal, Real qc, Real qi, Real exner)
    {
        return shoc::theta_from_thetal(thetal, qc, qi, exner);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real virtual_theta_from_shoc_state (Real thetal, Real qc, Real qi, Real qv, Real exner)
    {
        const Real theta = theta_from_shoc_state(thetal, qc, qi, exner);
        const Real ql = qc + qi;
        return theta * (1.0_rt + shoc_zvir() * qv - ql);
    }
}

void
ShocStructure::diagnose_surface_layer (ShocColumnData& col)
{
    auto ustar = col.ustar.array();
    auto obklen = col.obklen.array();
    auto wthv_sec = col.wthv_sec.array();
    auto pblh = col.pblh.array();
    const auto zi = col.zi.const_array();
    const auto thetal = col.thetal.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();
    const auto qv = col.qv.const_array();
    const auto exner = col.exner.const_array();
    const auto sflux = col.surf_sens_flux.const_array();
    const auto lflux = col.surf_lat_flux.const_array();
    const auto tauu = col.surf_tau_u.const_array();
    const auto tauv = col.surf_tau_v.const_array();
    const auto zt = col.zt.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));

    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        const Real cldliq_sfc = qc(ic,0,0) + qi(ic,0,0);
        const Real th_sfc = theta_from_shoc_state(thetal(ic,0,0), qc(ic,0,0), qi(ic,0,0), exner(ic,0,0));
        const Real thv_sfc = th_sfc * (1.0_rt + shoc_zvir() * qv(ic,0,0) - cldliq_sfc);
        const Real stress_mag = std::sqrt(tauu(ic,0,0) * tauu(ic,0,0) +
                                          tauv(ic,0,0) * tauv(ic,0,0));
        const Real ustar_val = std::sqrt(stress_mag);
        const Real kbfs = sflux(ic,0,0) + shoc_zvir() * th_sfc * lflux(ic,0,0);
        const Real sign_val = (kbfs >= 0.0_rt) ? shoc_kbfs_eps() : -shoc_kbfs_eps();

        ustar(ic,0,0) = amrex::max(shoc_u_star_min(), ustar_val);
        const Real ustar_cu = ustar(ic,0,0) * ustar(ic,0,0) * ustar(ic,0,0);
        obklen(ic,0,0) = -thv_sfc * ustar_cu /
                         (CONST_GRAV * KAPPA * (kbfs + sign_val));
        pblh(ic,0,0) = shoc::height_agl(zt(ic,0,0), zi(ic,0,0));

        // E3SM SHOC advances TKE using the carried buoyancy-flux profile from
        // the previous SHOC call. ERF overwrites only the surface entry here
        // with the current lower-boundary forcing and preserves the interior
        // profile that was diagnosed later in the previous SHOC call.
        wthv_sec(ic,0,0) = kbfs;
    });
}

void
ShocStructure::diagnose_pblh (ShocColumnData& col)
{
    auto pblh = col.pblh.array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto u = col.u.const_array();
    const auto v = col.v.const_array();
    const auto ustar = col.ustar.const_array();
    const auto obklen = col.obklen.const_array();
    const auto wthv_sec = col.wthv_sec.const_array();
    const auto thetal = col.thetal.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();
    const auto qv = col.qv.const_array();
    const auto exner = col.exner.const_array();
    const auto p_mid = col.p_mid.const_array();
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));

    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        const int npbl = diagnose_npbl(p_mid, layout, ic);
        const Real ustar_loc = ustar(ic,0,0);
        const Real z_sfc = zi(ic,0,0);
        const Real zt0_agl = shoc::height_agl(zt(ic,0,0), z_sfc);
        const Real thv0 = virtual_theta_from_shoc_state(thetal(ic,0,0), qc(ic,0,0), qi(ic,0,0),
                                                        qv(ic,0,0), exner(ic,0,0));
        Real pblh_loc = shoc::height_agl(zt(ic,npbl-1,0), z_sfc);
        Real prev_rino = 0.0_rt;
        bool found_pblh = false;

        for (int k = 1; k < npbl; ++k) {
            const Real ztk_agl = shoc::height_agl(zt(ic,k,0), z_sfc);
            const Real thvk = virtual_theta_from_shoc_state(thetal(ic,k,0), qc(ic,k,0), qi(ic,k,0),
                                                            qv(ic,k,0), exner(ic,k,0));
            const Real du = u(ic,k,0) - u(ic,0,0);
            const Real dv = v(ic,k,0) - v(ic,0,0);
            const Real vvk = amrex::max(1.0e-36_rt,
                                        du * du +
                                        dv * dv +
                                        shoc_pbl_fac() * ustar_loc * ustar_loc);
            const Real rino = CONST_GRAV * (thvk - thv0) * (ztk_agl - zt0_agl) /
                              (amrex::max(thv0, 1.0e-12_rt) * vvk);
            if (rino >= shoc_pbl_ricr()) {
                if (k == 1 || amrex::Math::abs(rino - prev_rino) <= 1.0e-12_rt) {
                    pblh_loc = ztk_agl;
                } else {
                    const Real ztkm1_agl = shoc::height_agl(zt(ic,k-1,0), z_sfc);
                    pblh_loc = ztkm1_agl + (shoc_pbl_ricr() - prev_rino) *
                               (ztk_agl - ztkm1_agl) / (rino - prev_rino);
                }
                found_pblh = true;
                break;
            }
            prev_rino = rino;
        }

        if (!found_pblh) {
            pblh_loc = shoc::height_agl(zt(ic,npbl-1,0), z_sfc);
        }

        if (wthv_sec(ic,0,0) > 0.0_rt) {
            const Real obk_abs = amrex::max(amrex::Math::abs(obklen(ic,0,0)), 1.0e-6_rt);
            const Real obk = std::copysign(obk_abs, obklen(ic,0,0));
            const Real binm = shoc_pbl_betam() * shoc_pbl_sffrac();
            const Real phiminv = std::cbrt(amrex::max(1.0e-12_rt, 1.0_rt - binm * pblh_loc / obk));
            const Real tlv = thv0 + wthv_sec(ic,0,0) * shoc_pbl_fak() /
                                   (amrex::max(ustar_loc, shoc_u_star_min()) * phiminv);
            prev_rino = 0.0_rt;
            for (int k = 1; k < npbl; ++k) {
                const Real ztk_agl = shoc::height_agl(zt(ic,k,0), z_sfc);
                const Real thvk = virtual_theta_from_shoc_state(thetal(ic,k,0), qc(ic,k,0), qi(ic,k,0),
                                                                qv(ic,k,0), exner(ic,k,0));
                const Real du = u(ic,k,0) - u(ic,0,0);
                const Real dv = v(ic,k,0) - v(ic,0,0);
                const Real vvk = amrex::max(1.0e-36_rt,
                                            du * du +
                                            dv * dv +
                                            shoc_pbl_fac() * ustar_loc * ustar_loc);
                const Real rino = CONST_GRAV * (thvk - tlv) * (ztk_agl - zt0_agl) /
                                  (amrex::max(thv0, 1.0e-12_rt) * vvk);
                if (rino >= shoc_pbl_ricr()) {
                    if (k == 1 || amrex::Math::abs(rino - prev_rino) <= 1.0e-12_rt) {
                        pblh_loc = ztk_agl;
                    } else {
                        const Real ztkm1_agl = shoc::height_agl(zt(ic,k-1,0), z_sfc);
                        pblh_loc = ztkm1_agl + (shoc_pbl_ricr() - prev_rino) *
                                   (ztk_agl - ztkm1_agl) / (rino - prev_rino);
                    }
                    break;
                }
                prev_rino = rino;
            }
        }

        pblh_loc = amrex::max(pblh_loc, 700.0_rt * ustar_loc);
        if (qc(ic,0,0) + qi(ic,0,0) > 0.0_rt && layout.nlev > 1) {
            pblh_loc = amrex::max(pblh_loc, shoc::height_agl(zi(ic,1,0), z_sfc) + 50.0_rt);
        }
        pblh(ic,0,0) = pblh_loc;
    });
}

void
ShocStructure::diagnose_length_and_brunt (ShocColumnData& col,
                                          const ShocRuntimeOptions& opts,
                                          Real dx,
                                          Real dy)
{
    auto brunt = col.brunt.array();
    auto shoc_mix = col.shoc_mix.array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto dz = col.dz.const_array();
    const auto thetal = col.thetal.const_array();
    const auto qc = col.qc.const_array();
    const auto qi = col.qi.const_array();
    const auto qv = col.qv.const_array();
    const auto exner = col.exner.const_array();
    const auto tke = col.tke.const_array();
    const Real max_horiz_len = std::sqrt(dx * dy);
    const auto layout = col.layout;
    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));

    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        const Real z_sfc = zi(ic,0,0);
        Real numerator = 0.0_rt;
        Real denom = 0.0_rt;
        for (int k = 0; k < layout.nlev; ++k) {
            const Real tke_sqrt = std::sqrt(amrex::max(tke(ic,k,0), shoc_min_tke()));
            const Real zt_agl = shoc::height_agl(zt(ic,k,0), z_sfc);
            numerator += tke_sqrt * zt_agl * dz(ic,k,0);
            denom += tke_sqrt * dz(ic,k,0);
        }
        const Real l_inf = (denom > 0.0_rt) ? 0.1_rt * numerator / denom : shoc_min_len();

        for (int k = 0; k < layout.nlev; ++k) {
            const Real theta_v_k = virtual_theta_from_shoc_state(thetal(ic,k,0), qc(ic,k,0), qi(ic,k,0),
                                                                 qv(ic,k,0), exner(ic,k,0));
            Real theta_v_lo = theta_v_k;
            Real theta_v_hi = theta_v_k;
            if (layout.nlev > 1) {
                if (k == 0) {
                    const Real theta_v_kp1 = virtual_theta_from_shoc_state(thetal(ic,1,0), qc(ic,1,0), qi(ic,1,0),
                                                                           qv(ic,1,0), exner(ic,1,0));
                    theta_v_lo = weighted_linear_interp(zt(ic,0,0), zt(ic,1,0),
                                                        theta_v_k, theta_v_kp1, zi(ic,0,0));
                } else {
                    const Real theta_v_km1 = virtual_theta_from_shoc_state(thetal(ic,k-1,0), qc(ic,k-1,0), qi(ic,k-1,0),
                                                                           qv(ic,k-1,0), exner(ic,k-1,0));
                    theta_v_lo = weighted_linear_interp(zt(ic,k-1,0), zt(ic,k,0),
                                                        theta_v_km1, theta_v_k, zi(ic,k,0));
                }

                if (k == layout.nlev - 1) {
                    const Real theta_v_km1 = virtual_theta_from_shoc_state(thetal(ic,layout.nlev-2,0),
                                                                           qc(ic,layout.nlev-2,0), qi(ic,layout.nlev-2,0),
                                                                           qv(ic,layout.nlev-2,0), exner(ic,layout.nlev-2,0));
                    theta_v_hi = weighted_linear_interp(zt(ic,layout.nlev-2,0), zt(ic,layout.nlev-1,0),
                                                        theta_v_km1, theta_v_k, zi(ic,layout.nlev,0));
                } else {
                    const Real theta_v_kp1 = virtual_theta_from_shoc_state(thetal(ic,k+1,0), qc(ic,k+1,0), qi(ic,k+1,0),
                                                                           qv(ic,k+1,0), exner(ic,k+1,0));
                    theta_v_hi = weighted_linear_interp(zt(ic,k,0), zt(ic,k+1,0),
                                                        theta_v_k, theta_v_kp1, zi(ic,k+1,0));
                }
            }

            // ERF columns use bottom-up indexing, so the upper interface is
            // k+1 and the lower interface is k. Stable stratification must
            // therefore yield positive Brunt-Vaisala frequency.
            brunt(ic,k,0) = (CONST_GRAV / amrex::max(theta_v_k, 1.0e-12_rt)) *
                            (theta_v_hi - theta_v_lo) / amrex::max(dz(ic,k,0), 1.0e-12_rt);
            const Real tkes = std::sqrt(amrex::max(tke(ic,k,0), shoc_min_tke()));
            const Real brunt_pos = amrex::max(brunt(ic,k,0), 0.0_rt);
            const Real zt_agl = shoc::height_agl(zt(ic,k,0), z_sfc);
            const Real inv_term = (1.0_rt / amrex::max(400.0_rt * tkes * KAPPA * amrex::max(zt_agl, 1.0_rt), 1.0e-12_rt)) +
                                  (1.0_rt / amrex::max(400.0_rt * tkes * amrex::max(l_inf, shoc_min_len()), 1.0e-12_rt)) +
                                  0.01_rt * brunt_pos / amrex::max(tke(ic,k,0), shoc_min_tke());
            Real mix = amrex::min(shoc_max_len(),
                                  2.8284_rt * std::sqrt(1.0_rt / amrex::max(inv_term, 1.0e-12_rt)) /
                                  amrex::max(opts.length_fac, 1.0e-12_rt));
            mix = amrex::min(max_horiz_len, amrex::max(shoc_min_len(), mix));
            shoc_mix(ic,k,0) = mix;
        }
    });
}
