#include "ERF_ShocEnergyFixer.H"

using namespace amrex;

int
ShocEnergyFixer::diagnose_active_top (const Vector<Real>& tke)
{
    int shoc_top = -1;
    for (int k = static_cast<int>(tke.size()) - 1; k >= 0; --k) {
        if (tke[k] > shoc::constants::min_tke()) {
            shoc_top = k;
            break;
        }
    }
    return shoc_top;
}

void
ShocEnergyFixer::apply_column (const ShocColumnData& col,
                               int ic,
                               Real dt,
                               const Vector<Real>& thl_old,
                               const Vector<Real>& qv_old,
                               const Vector<Real>& qc_old,
                               const Vector<Real>& qi_old,
                               const Vector<Real>& u_old,
                               const Vector<Real>& v_old,
                               const Vector<Real>& tke_old,
                               Vector<Real>& thl_new,
                               const Vector<Real>& qv_new,
                               const Vector<Real>& qc_new,
                               const Vector<Real>& qi_new,
                               const Vector<Real>& u_new,
                               const Vector<Real>& v_new,
                               const Vector<Real>& tke_new)
{
    ShocColumnData tmp;
    define_shoc_column_data(tmp, col.layout, amrex::The_Arena(), shoc::InitRunOn::Host);

    auto rho = tmp.rho.array();
    auto dz = tmp.dz.array();
    auto zt = tmp.zt.array();
    auto exner = tmp.exner.array();
    auto surf_sens_flux = tmp.surf_sens_flux.array();
    auto surf_lat_flux = tmp.surf_lat_flux.array();
    auto thetal_base = tmp.thetal_base.array();
    auto qv_base = tmp.qv_base.array();
    auto qc_base = tmp.qc_base.array();
    auto qi_base = tmp.qi_base.array();
    auto u_base = tmp.u_base.array();
    auto v_base = tmp.v_base.array();
    auto tke_base = tmp.tke_base_state.array();
    auto thetal = tmp.thetal.array();
    auto qw = tmp.qw.array();
    auto u = tmp.u.array();
    auto v = tmp.v.array();
    auto tke = tmp.tke.array();
    auto shoc_ql = tmp.shoc_ql.array();

    const auto rho_in = col.rho.const_array();
    const auto dz_in = col.dz.const_array();
    const auto zt_in = col.zt.const_array();
    const auto exner_in = col.exner.const_array();
    const auto surf_sens_flux_in = col.surf_sens_flux.const_array();
    const auto surf_lat_flux_in = col.surf_lat_flux.const_array();

    for (int k = 0; k < col.layout.nlev; ++k) {
        rho(ic,k,0) = rho_in(ic,k,0);
        dz(ic,k,0) = dz_in(ic,k,0);
        zt(ic,k,0) = zt_in(ic,k,0);
        exner(ic,k,0) = exner_in(ic,k,0);
        thetal_base(ic,k,0) = thl_old[k];
        qv_base(ic,k,0) = qv_old[k];
        qc_base(ic,k,0) = qc_old[k];
        qi_base(ic,k,0) = qi_old[k];
        u_base(ic,k,0) = u_old[k];
        v_base(ic,k,0) = v_old[k];
        tke_base(ic,k,0) = tke_old[k];
        thetal(ic,k,0) = thl_new[k];
        qw(ic,k,0) = qv_new[k] + qc_new[k] + qi_new[k];
        u(ic,k,0) = u_new[k];
        v(ic,k,0) = v_new[k];
        tke(ic,k,0) = tke_new[k];
        shoc_ql(ic,k,0) = qc_new[k] + qi_new[k];
    }

    surf_sens_flux(ic,0,0) = surf_sens_flux_in(ic,0,0);
    surf_lat_flux(ic,0,0) = surf_lat_flux_in(ic,0,0);

    shoc::apply_energy_fix_column(shoc::make_energy_fixer_view(tmp), ic, dt);

    const auto thetal_out = tmp.thetal.const_array();
    for (int k = 0; k < col.layout.nlev; ++k) {
        thl_new[k] = thetal_out(ic,k,0);
    }
}
