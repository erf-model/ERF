#include "ERF_ShocPreprocess.H"

#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_MoistUtils.H"
#include "ERF_ShocThermoUtils.H"

using namespace amrex;

namespace
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real load_q (const Array4<const Real>& cons_arr,
             int i, int j, int k,
             int comp, Real rho, int ncomp)
{
    return shoc_valid_comp(comp, ncomp) ? cons_arr(i,j,k,comp) / rho : 0.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real four_node_z_avg (const Array4<const Real>& z, int i, int j, int k) noexcept
{
    return 0.25_rt * (z(i  ,j  ,k) +
                      z(i+1,j  ,k) +
                      z(i  ,j+1,k) +
                      z(i+1,j+1,k));
}
}

void
ShocPreprocess::fill_columns (ShocColumnData& col,
                              const MFIter& mfi,
                              const MultiFab& cons,
                              const MultiFab& xvel,
                              const MultiFab& yvel,
                              const MultiFab& zvel,
                              const MultiFab* hfx3,
                              const MultiFab* qfx3,
                              const MultiFab* tau13,
                              const MultiFab* tau23,
                              const MultiFab& z_phys_nd,
                              const Geometry& geom,
                              const MoistureComponentIndices& moisture_indices)
{
    const auto& cons_arr = cons.const_array(mfi);
    const auto& u_arr = xvel.const_array(mfi);
    const auto& v_arr = yvel.const_array(mfi);
    const auto& w_arr = zvel.const_array(mfi);
    const auto& z_arr = z_phys_nd.const_array(mfi);

    const Array4<const Real> hfx_arr = hfx3 ? hfx3->const_array(mfi) : Array4<const Real>{};
    const Array4<const Real> qfx_arr = qfx3 ? qfx3->const_array(mfi) : Array4<const Real>{};
    const Array4<const Real> t13_arr = tau13 ? tau13->const_array(mfi) : Array4<const Real>{};
    const Array4<const Real> t23_arr = tau23 ? tau23->const_array(mfi) : Array4<const Real>{};

    auto zt_arr = col.zt.array();
    auto zi_arr = col.zi.array();
    auto dz_arr = col.dz.array();
    auto p_mid_arr = col.p_mid.array();
    auto p_int_arr = col.p_int.array();
    auto rho_arr = col.rho.array();
    auto theta_arr = col.theta.array();
    auto exner_arr = col.exner.array();
    auto theta_v_arr = col.theta_v.array();
    auto thetal_arr = col.thetal.array();
    auto qv_arr = col.qv.array();
    auto qc_arr = col.qc.array();
    auto qi_arr = col.qi.array();
    auto qw_arr = col.qw.array();
    auto shoc_ql_arr = col.shoc_ql.array();
    auto tabs_arr = col.tabs.array();
    auto tke_arr = col.tke.array();
    auto ucol_arr = col.u.array();
    auto vcol_arr = col.v.array();
    auto wcol_arr = col.w.array();
    auto dse_arr = col.host_dse.array();
    auto sflux_arr = col.surf_sens_flux.array();
    auto lflux_arr = col.surf_lat_flux.array();
    auto tauu_arr = col.surf_tau_u.array();
    auto tauv_arr = col.surf_tau_v.array();

    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const int klo = col.layout.kmin;
    const int khi = col.layout.kmax;
    const int ncomp = cons.nComp();
    const auto layout = col.layout;
    amrex::ignore_unused(problo, dx);
    const Box xy_box = amrex::makeSlab(mfi.validbox(), 2, klo);

    ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        const int ic = shoc_column_index(layout, i, j);
        const Real rho_sfc = amrex::max(cons_arr(i,j,klo,Rho_comp), 1.0e-12_rt);

        // ERF host fluxes/stresses are stored as density-weighted lower
        // boundary fluxes. SHOC consumes kinematic surface fluxes.
        if (hfx3) sflux_arr(ic,0,0) = hfx_arr(i,j,klo) / rho_sfc;
        if (qfx3) lflux_arr(ic,0,0) = qfx_arr(i,j,klo) / rho_sfc;
        if (tau13) {
            tauu_arr(ic,0,0) = 0.5_rt * (t13_arr(i,j,klo) + t13_arr(i+1,j,klo)) / rho_sfc;
        }
        if (tau23) {
            tauv_arr(ic,0,0) = 0.5_rt * (t23_arr(i,j,klo) + t23_arr(i,j+1,klo)) / rho_sfc;
        }

        for (int k = klo; k <= khi; ++k) {
            const int kk = k - klo;
            // z_phys_nd is nodal in x, y, and z. Pack SHOC column interface
            // heights as the four-node horizontal average over the ERF cell
            // footprint, matching the EAMxx SHOC coupling.
            const Real zlo = four_node_z_avg(z_arr, i, j, k);
            const Real zhi = four_node_z_avg(z_arr, i, j, k+1);
            const Real zc = 0.5_rt * (zlo + zhi);
            const Real dz = zhi - zlo;
            const Real rho = cons_arr(i,j,k,Rho_comp);
            const Real theta = cons_arr(i,j,k,RhoTheta_comp) / rho;
            const Real qv = load_q(cons_arr, i, j, k, moisture_indices.qv, rho, ncomp);
            const Real qc = load_q(cons_arr, i, j, k, moisture_indices.qc, rho, ncomp);
            const Real qi = load_q(cons_arr, i, j, k, moisture_indices.qi, rho, ncomp);
            const Real p = getPgivenRTh(cons_arr(i,j,k,RhoTheta_comp), qv);
            const Real tabs = getTgivenRandRTh(rho, cons_arr(i,j,k,RhoTheta_comp), qv);
            const Real ql_np = qc + qi;
            const Real exner = tabs / amrex::max(theta, 1.0e-12_rt);
            // SHOC carries a liquid-water potential temperature-like variable.
            // For native SHOC coupling, include phase-aware condensate
            // heating so theta = theta_l + (Lv*qc + Ls*qi)/(Cp*exner).
            const Real thetal = shoc::thetal_from_theta(theta, qc, qi, exner);
            const Real theta_v = theta * (1.0_rt + 0.61_rt * qv - ql_np);
            const Real qke = cons_arr(i,j,k,RhoKE_comp) / rho;

            zt_arr(ic,kk,0) = zc;
            zi_arr(ic,kk,0) = zlo;
            zi_arr(ic,kk+1,0) = zhi;
            dz_arr(ic,kk,0) = dz;
            p_mid_arr(ic,kk,0) = p;
            rho_arr(ic,kk,0) = rho;
            theta_arr(ic,kk,0) = theta;
            exner_arr(ic,kk,0) = exner;
            theta_v_arr(ic,kk,0) = theta_v;
            thetal_arr(ic,kk,0) = thetal;
            qv_arr(ic,kk,0) = qv;
            qc_arr(ic,kk,0) = qc;
            qi_arr(ic,kk,0) = qi;
            qw_arr(ic,kk,0) = qv + qc + qi;
            shoc_ql_arr(ic,kk,0) = ql_np;
            tabs_arr(ic,kk,0) = tabs;
            tke_arr(ic,kk,0) = qke;
            dse_arr(ic,kk,0) = Cp_d * tabs + CONST_GRAV * zc;

            ucol_arr(ic,kk,0) = 0.5_rt * (u_arr(i,j,k) + u_arr(i+1,j,k));
            vcol_arr(ic,kk,0) = 0.5_rt * (v_arr(i,j,k) + v_arr(i,j+1,k));
            wcol_arr(ic,kk,0) = 0.5_rt * (w_arr(i,j,k) + w_arr(i,j,k+1));
        }

        p_int_arr(ic,0,0) = p_mid_arr(ic,0,0);
        for (int k = klo+1; k <= khi; ++k) {
            const int kk = k - klo;
            p_int_arr(ic,kk,0) = 0.5_rt * (p_mid_arr(ic,kk-1,0) + p_mid_arr(ic,kk,0));
        }
        p_int_arr(ic,layout.nlev,0) = p_mid_arr(ic,layout.nlev-1,0);
    });
}
