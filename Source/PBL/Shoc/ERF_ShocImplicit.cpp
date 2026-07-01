#include "ERF_ShocImplicit.H"

#include "ERF_Constants.H"
#include "ERF_ShocEnergyFixer.H"
#include "ERF_ShocThermoUtils.H"
#include "ERF_SolveTridiag.H"

#include <AMReX_BLProfiler.H>

#include <algorithm>
#include <cmath>

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_zvir () noexcept { return 0.61_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_u_ws_min () noexcept { return 1.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_ksrf_min () noexcept { return 1.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_ustar_min () noexcept { return 0.01_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_min_tke () noexcept { return 4.0e-4_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_max_tke () noexcept { return 50.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real shoc_min_qw () noexcept { return 0.0_rt; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real weighted_linear_interp (Real x0, Real x1, Real y0, Real y1, Real x)
    {
        const Real denom = x1 - x0;
        if (amrex::Math::abs(denom) <= 1.0e-12_rt) {
            return 0.5_rt * (y0 + y1);
        }
        return y0 + (y1 - y0) * (x - x0) / denom;
    }

    void interpolate_cc_to_iface (const ShocColumnData& col,
                                  const FArrayBox& cc_in,
                                  FArrayBox& iface_out)
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
                out(ic,0,0) = in(ic,0,0);
                out(ic,1,0) = in(ic,0,0);
                return;
            }

            out(ic,0,0) = weighted_linear_interp(zt(ic,0,0), zt(ic,1,0),
                                                 in(ic,0,0), in(ic,1,0),
                                                 zi(ic,0,0));
            for (int k = 1; k < layout.nlev; ++k) {
                out(ic,k,0) = weighted_linear_interp(zt(ic,k-1,0), zt(ic,k,0),
                                                     in(ic,k-1,0), in(ic,k,0),
                                                     zi(ic,k,0));
            }
            out(ic,layout.nlev,0) =
                weighted_linear_interp(zt(ic,layout.nlev-2,0), zt(ic,layout.nlev-1,0),
                                       in(ic,layout.nlev-2,0), in(ic,layout.nlev-1,0),
                                       zi(ic,layout.nlev,0));
        });
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
}

void
ShocImplicit::cache_baseline_state (ShocColumnData& col)
{
    auto thetal = col.thetal.const_array();
    auto theta = col.theta.const_array();
    auto qv = col.qv.const_array();
    auto qc = col.qc.const_array();
    auto qi = col.qi.const_array();
    auto u = col.u.const_array();
    auto v = col.v.const_array();
    auto tke = col.tke.const_array();

    auto thetal_base = col.thetal_base.array();
    auto theta_base = col.theta_base.array();
    auto qv_base = col.qv_base.array();
    auto qc_base = col.qc_base.array();
    auto qi_base = col.qi_base.array();
    auto u_base = col.u_base.array();
    auto v_base = col.v_base.array();
    auto tke_base = col.tke_base_state.array();

    const Box cell_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, col.layout.nlev - 1, 0));
    ParallelFor(cell_box, [=] AMREX_GPU_DEVICE (int ic, int k, int) noexcept
    {
        thetal_base(ic,k,0) = thetal(ic,k,0);
        theta_base(ic,k,0) = theta(ic,k,0);
        qv_base(ic,k,0) = qv(ic,k,0);
        qc_base(ic,k,0) = qc(ic,k,0);
        qi_base(ic,k,0) = qi(ic,k,0);
        u_base(ic,k,0) = u(ic,k,0);
        v_base(ic,k,0) = v(ic,k,0);
        tke_base(ic,k,0) = tke(ic,k,0);
    });
}

void
ShocImplicit::compute_tmpi (const ShocColumnData& col,
                            int ic,
                            Real dt,
                            const Vector<Real>& rho_zi,
    Vector<Real>& tmpi)
{
    tmpi.assign(col.layout.nlev + 1, 0.0_rt);
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto layout = col.layout;
    for (int k = 0; k <= col.layout.nlev; ++k) {
        tmpi[k] = dt * CONST_GRAV * amrex::max(rho_zi[k], 1.0e-12_rt) /
                  interface_spacing(zt, zi, layout, ic, k);
    }
}

void
ShocImplicit::compute_dp_inverse (const ShocColumnData& col,
                                  int ic,
                                  Vector<Real>& rdp_zt)
{
    rdp_zt.assign(col.layout.nlev, 0.0);
    const auto rho = col.rho.const_array();
    const auto dz = col.dz.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        rdp_zt[k] = 1.0_rt / amrex::max(CONST_GRAV * rho(ic,k,0) * dz(ic,k,0), 1.0e-12_rt);
    }
}

AMREX_GPU_HOST_DEVICE
Real
ShocImplicit::compute_temperature (Real thetal,
                                   Real ql,
                                   Real exner)
{
    return shoc::temperature_from_thetal(thetal, ql, 0.0_rt, exner);
}

AMREX_GPU_HOST_DEVICE
Real
ShocImplicit::compute_vapor (Real qw,
                             Real ql)
{
    return amrex::max(0.0_rt, qw - amrex::max(0.0_rt, ql));
}

void
ShocImplicit::advance_implicit_state (ShocColumnData& col,
                                      const ShocRuntimeOptions& opts,
                                      Real dt)
{
    AMREX_ALWAYS_ASSERT(dt > 0.0);
    static_cast<void>(opts);

    const int nlev = col.layout.nlev;
    if (nlev == 0) {
        return;
    }

    const Box iface_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, nlev, 0));
    const Box tri_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, 0, nlev - 1));
    FArrayBox rho_zi_fab(iface_box, 1, The_Async_Arena());
    FArrayBox tk_zi_fab(iface_box, 1, The_Async_Arena());
    FArrayBox tkh_zi_fab(iface_box, 1, The_Async_Arena());
    FArrayBox tmpi_fab(iface_box, 1, The_Async_Arena());
    FArrayBox rdp_zt_fab(Box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, nlev - 1, 0)), 1, The_Async_Arena());
    FArrayBox rhs_fab(tri_box, 1, The_Async_Arena());
    FArrayBox soln_fab(tri_box, 1, The_Async_Arena());
    FArrayBox coeffA_fab(tri_box, 1, The_Async_Arena());
    FArrayBox coeffB_fab(tri_box, 1, The_Async_Arena());
    FArrayBox inv_coeffB_fab(tri_box, 1, The_Async_Arena());
    FArrayBox coeffC_fab(tri_box, 1, The_Async_Arena());

    interpolate_cc_to_iface(col, col.rho, rho_zi_fab);
    interpolate_cc_to_iface(col, col.tk, tk_zi_fab);
    interpolate_cc_to_iface(col, col.tkh, tkh_zi_fab);

    auto thetal = col.thetal.array();
    auto qw = col.qw.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();

    const auto surf_tau_u = col.surf_tau_u.const_array();
    const auto surf_tau_v = col.surf_tau_v.const_array();
    const auto surf_sens_flux = col.surf_sens_flux.const_array();
    const auto surf_lat_flux = col.surf_lat_flux.const_array();
    const auto rho = col.rho.const_array();
    const auto dz = col.dz.const_array();
    const auto zt = col.zt.const_array();
    const auto zi = col.zi.const_array();
    const auto rho_zi = rho_zi_fab.const_array();
    const auto tk_zi = tk_zi_fab.const_array();
    const auto tkh_zi = tkh_zi_fab.const_array();
    const auto layout = col.layout;
    auto tmpi = tmpi_fab.array();
    auto rdp_zt = rdp_zt_fab.array();
    auto rhs = rhs_fab.array();
    auto soln = soln_fab.array();
    auto coeffA = coeffA_fab.array();
    auto coeffB = coeffB_fab.array();
    auto inv_coeffB = inv_coeffB_fab.array();
    auto coeffC = coeffC_fab.array();

    const Box col_box(IntVect(0,0,0), IntVect(layout.ncell - 1, 0, 0));
    ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
    {
        for (int k = 0; k <= nlev; ++k) {
            tmpi(ic,k,0) = dt * CONST_GRAV * amrex::max(rho_zi(ic,k,0), 1.0e-12_rt) /
                           interface_spacing(zt, zi, layout, ic, k);
            if (k < nlev) {
                rdp_zt(ic,k,0) = 1.0_rt / amrex::max(CONST_GRAV * rho(ic,k,0) * dz(ic,k,0), 1.0e-12_rt);
            }
        }

        // E3SM: cmnfac = dtime * g * rho_zi(nlevi-1) * rdp_zt(nlev-1)
        // The tmpi array includes a 1/dz_zi interface-spacing factor that
        // belongs only in the tridiagonal diffusion matrix, NOT in the
        // explicit surface-flux injection.
        const Real cmnfac = dt * CONST_GRAV * amrex::max(rho_zi(ic,0,0), 1.0e-12_rt) * rdp_zt(ic,0,0);

        const Real uw_sfc = surf_tau_u(ic,0,0);
        const Real vw_sfc = surf_tau_v(ic,0,0);
        const Real stress_mag = std::sqrt(uw_sfc * uw_sfc + vw_sfc * vw_sfc);
        Real ksrf = 0.0_rt;
        Real wtke_sfc = 0.0_rt;
        if (stress_mag > 1.0e-12_rt) {
            const Real ws = amrex::max(std::sqrt(u(ic,0,0) * u(ic,0,0) + v(ic,0,0) * v(ic,0,0)), shoc_u_ws_min());
            const Real rho_zi_sfc = amrex::max(rho_zi(ic,0,0), 1.0e-12_rt);
            const Real tau_u = rho_zi_sfc * uw_sfc;
            const Real tau_v = rho_zi_sfc * vw_sfc;
            const Real tau = std::sqrt(tau_u * tau_u + tau_v * tau_v);
            ksrf = amrex::max(tau / ws, shoc_ksrf_min());
            const Real ustar = amrex::max(std::sqrt(stress_mag), shoc_ustar_min());
            wtke_sfc = ustar * ustar * ustar;
        }

        for (int k = 0; k < nlev; ++k) {
            rhs(ic,0,k) = u(ic,k,0);
            coeffA(ic,0,k) = (k > 0) ? -tk_zi(ic,k,0) * tmpi(ic,k,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffC(ic,0,k) = (k < nlev - 1) ? -tk_zi(ic,k+1,0) * tmpi(ic,k+1,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffB(ic,0,k) = 1.0_rt - coeffA(ic,0,k) - coeffC(ic,0,k);
        }
        coeffB(ic,0,0) += ksrf * dt * CONST_GRAV * rdp_zt(ic,0,0);
        SolveTridiag(ic, 0, 0, nlev - 1, soln, coeffA, coeffB, inv_coeffB, coeffC, rhs);
        for (int k = 0; k < nlev; ++k) {
            u(ic,k,0) = soln(ic,0,k);
        }

        for (int k = 0; k < nlev; ++k) {
            rhs(ic,0,k) = v(ic,k,0);
            coeffA(ic,0,k) = (k > 0) ? -tk_zi(ic,k,0) * tmpi(ic,k,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffC(ic,0,k) = (k < nlev - 1) ? -tk_zi(ic,k+1,0) * tmpi(ic,k+1,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffB(ic,0,k) = 1.0_rt - coeffA(ic,0,k) - coeffC(ic,0,k);
        }
        coeffB(ic,0,0) += ksrf * dt * CONST_GRAV * rdp_zt(ic,0,0);
        SolveTridiag(ic, 0, 0, nlev - 1, soln, coeffA, coeffB, inv_coeffB, coeffC, rhs);
        for (int k = 0; k < nlev; ++k) {
            v(ic,k,0) = soln(ic,0,k);
        }

        for (int k = 0; k < nlev; ++k) {
            rhs(ic,0,k) = thetal(ic,k,0);
            coeffA(ic,0,k) = (k > 0) ? -tkh_zi(ic,k,0) * tmpi(ic,k,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffC(ic,0,k) = (k < nlev - 1) ? -tkh_zi(ic,k+1,0) * tmpi(ic,k+1,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffB(ic,0,k) = 1.0_rt - coeffA(ic,0,k) - coeffC(ic,0,k);
        }
        rhs(ic,0,0) += cmnfac * surf_sens_flux(ic,0,0);
        SolveTridiag(ic, 0, 0, nlev - 1, soln, coeffA, coeffB, inv_coeffB, coeffC, rhs);
        for (int k = 0; k < nlev; ++k) {
            thetal(ic,k,0) = soln(ic,0,k);
        }

        for (int k = 0; k < nlev; ++k) {
            rhs(ic,0,k) = qw(ic,k,0);
            coeffA(ic,0,k) = (k > 0) ? -tkh_zi(ic,k,0) * tmpi(ic,k,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffC(ic,0,k) = (k < nlev - 1) ? -tkh_zi(ic,k+1,0) * tmpi(ic,k+1,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffB(ic,0,k) = 1.0_rt - coeffA(ic,0,k) - coeffC(ic,0,k);
        }
        rhs(ic,0,0) += cmnfac * surf_lat_flux(ic,0,0);
        SolveTridiag(ic, 0, 0, nlev - 1, soln, coeffA, coeffB, inv_coeffB, coeffC, rhs);
        for (int k = 0; k < nlev; ++k) {
            qw(ic,k,0) = amrex::max(soln(ic,0,k), shoc_min_qw());
        }

        for (int k = 0; k < nlev; ++k) {
            rhs(ic,0,k) = shoc_clamp(tke(ic,k,0), shoc_min_tke(), shoc_max_tke());
            coeffA(ic,0,k) = (k > 0) ? -tkh_zi(ic,k,0) * tmpi(ic,k,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffC(ic,0,k) = (k < nlev - 1) ? -tkh_zi(ic,k+1,0) * tmpi(ic,k+1,0) * rdp_zt(ic,k,0) : 0.0_rt;
            coeffB(ic,0,k) = 1.0_rt - coeffA(ic,0,k) - coeffC(ic,0,k);
        }
        rhs(ic,0,0) += cmnfac * wtke_sfc;
        SolveTridiag(ic, 0, 0, nlev - 1, soln, coeffA, coeffB, inv_coeffB, coeffC, rhs);
        for (int k = 0; k < nlev; ++k) {
            tke(ic,k,0) = shoc_clamp(soln(ic,0,k), shoc_min_tke(), shoc_max_tke());
        }
    });
}

void
ShocImplicit::finalize_from_pdf (ShocColumnData& col,
                                 const ShocRuntimeOptions& opts,
                                 Real dt)
{
    AMREX_ALWAYS_ASSERT(dt > 0.0);

    const int nlev = col.layout.nlev;
    if (nlev == 0) {
        return;
    }

    auto thetal = col.thetal.array();
    auto theta = col.theta.array();
    auto theta_v = col.theta_v.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto tabs = col.tabs.array();
    auto host_dse = col.host_dse.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();
    auto shoc_ql = col.shoc_ql.array();
    auto theta_tend = col.theta_tend.array();
    auto qv_tend = col.qv_tend.array();
    auto qc_tend = col.qc_tend.array();
    auto qi_tend = col.qi_tend.array();
    auto u_tend = col.u_tend.array();
    auto v_tend = col.v_tend.array();
    auto tke_tend = col.tke_tend.array();

    const auto zt = col.zt.const_array();
    const auto exner = col.exner.const_array();
    const auto shoc_ql_pdf = col.shoc_ql.const_array();
    const auto theta_base = col.theta_base.const_array();
    const auto qv_base = col.qv_base.const_array();
    const auto qc_base = col.qc_base.const_array();
    const auto qi_base = col.qi_base.const_array();
    const auto u_base = col.u_base.const_array();
    const auto v_base = col.v_base.const_array();
    const auto tke_base = col.tke_base_state.const_array();
    const auto energy_view = shoc::make_energy_fixer_view(col);
    const Box col_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, 0, 0));
    const Box cell_box(IntVect(0,0,0), IntVect(col.layout.ncell - 1, nlev - 1, 0));

    {
        BL_PROFILE("SHOC::advance::finalize::energy_fix");
        ParallelFor(col_box, [=] AMREX_GPU_DEVICE (int ic, int, int) noexcept
        {
            shoc::apply_energy_fix_column(energy_view, ic, dt);
        });
    }

    {
        BL_PROFILE("SHOC::advance::finalize::reconstruct");
        ParallelFor(cell_box, [=] AMREX_GPU_DEVICE (int ic, int k, int) noexcept
        {
            const Real qw_new_loc = amrex::max(qw(ic,k,0), shoc_min_qw());
            const Real ql_base = amrex::max(0.0_rt, qc_base(ic,k,0) + qi_base(ic,k,0));
            const Real pdf_ql_raw = shoc_clamp(shoc_ql_pdf(ic,k,0), 0.0_rt, qw_new_loc);
            const Real pdf_ql = opts.debug_disable_pdf_cloud_increment
                ? shoc_clamp(ql_base, 0.0_rt, qw_new_loc)
                : pdf_ql_raw;
            Real tabs_new = 0.0_rt;
            Real qv_new_loc = 0.0_rt;
            Real qc_new_loc = 0.0_rt;
            Real qi_new_loc = 0.0_rt;
            shoc::reconstruct_pdf_state(thetal(ic,k,0), qw_new_loc,
                                        exner(ic,k,0), qi_base(ic,k,0), pdf_ql,
                                        tabs_new, qv_new_loc, qc_new_loc, qi_new_loc);

            const Real ql_total = qc_new_loc + qi_new_loc;
            tabs(ic,k,0) = tabs_new;
            theta(ic,k,0) = tabs_new / amrex::max(exner(ic,k,0), 1.0e-12_rt);
            qw(ic,k,0) = qw_new_loc;
            qv(ic,k,0) = qv_new_loc;
            qc(ic,k,0) = qc_new_loc;
            qi(ic,k,0) = qi_new_loc;
            shoc_ql(ic,k,0) = ql_total;
            theta_v(ic,k,0) = theta(ic,k,0) * (1.0_rt + shoc_zvir() * qv_new_loc - ql_total);
            host_dse(ic,k,0) = Cp_d * tabs_new + CONST_GRAV * zt(ic,k,0);

            theta_tend(ic,k,0) = (theta(ic,k,0) - theta_base(ic,k,0)) / dt;
            qv_tend(ic,k,0) = (qv_new_loc - qv_base(ic,k,0)) / dt;
            qc_tend(ic,k,0) = (qc_new_loc - qc_base(ic,k,0)) / dt;
            qi_tend(ic,k,0) = (qi_new_loc - qi_base(ic,k,0)) / dt;
            u_tend(ic,k,0) = (u(ic,k,0) - u_base(ic,k,0)) / dt;
            v_tend(ic,k,0) = (v(ic,k,0) - v_base(ic,k,0)) / dt;
            tke_tend(ic,k,0) = (tke(ic,k,0) - tke_base(ic,k,0)) / dt;
        });
    }
}

void
ShocImplicit::update_prognostics (ShocColumnData& col,
                                  const ShocRuntimeOptions& opts,
                                  Real dt)
{
    cache_baseline_state(col);
    advance_implicit_state(col, opts, dt);
    finalize_from_pdf(col, opts, dt);
}
