#ifndef ERF_SHOC_ENERGY_FIXER_H_
#define ERF_SHOC_ENERGY_FIXER_H_

#include "ERF_ShocColumnData.H"
#include "ERF_ShocConstants.H"
#include "ERF_ShocThermoUtils.H"

#include "ERF_Constants.H"

namespace shoc {

struct ShocEnergyFixerView
{
    amrex::Array4<const amrex::Real> rho;
    amrex::Array4<const amrex::Real> dz;
    amrex::Array4<const amrex::Real> zt;
    amrex::Array4<const amrex::Real> exner;
    amrex::Array4<const amrex::Real> surf_sens_flux;
    amrex::Array4<const amrex::Real> surf_lat_flux;

    amrex::Array4<const amrex::Real> thetal_base;
    amrex::Array4<const amrex::Real> qv_base;
    amrex::Array4<const amrex::Real> qc_base;
    amrex::Array4<const amrex::Real> qi_base;
    amrex::Array4<const amrex::Real> u_base;
    amrex::Array4<const amrex::Real> v_base;
    amrex::Array4<const amrex::Real> tke_base;

    amrex::Array4<amrex::Real> thetal;
    amrex::Array4<const amrex::Real> qw;
    amrex::Array4<const amrex::Real> u;
    amrex::Array4<const amrex::Real> v;
    amrex::Array4<const amrex::Real> tke;
    amrex::Array4<const amrex::Real> shoc_ql;

    int nlev = 0;
};

inline ShocEnergyFixerView
make_energy_fixer_view (ShocColumnData& col)
{
    return {
        col.rho.const_array(),
        col.dz.const_array(),
        col.zt.const_array(),
        col.exner.const_array(),
        col.surf_sens_flux.const_array(),
        col.surf_lat_flux.const_array(),
        col.thetal_base.const_array(),
        col.qv_base.const_array(),
        col.qc_base.const_array(),
        col.qi_base.const_array(),
        col.u_base.const_array(),
        col.v_base.const_array(),
        col.tke_base_state.const_array(),
        col.thetal.array(),
        col.qw.const_array(),
        col.u.const_array(),
        col.v.const_array(),
        col.tke.const_array(),
        col.shoc_ql.const_array(),
        col.layout.nlev
    };
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
diagnose_active_top (const amrex::Array4<const amrex::Real>& tke,
                     int ic,
                     int nlev) noexcept
{
    int shoc_top = -1;
    for (int k = nlev - 1; k >= 0; --k) {
        if (tke(ic,k,0) > constants::min_tke()) {
            shoc_top = k;
            break;
        }
    }
    return shoc_top;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
apply_energy_fix_column (const ShocEnergyFixerView& v,
                         int ic,
                         amrex::Real dt) noexcept
{
    const int shoc_top = diagnose_active_top(v.tke, ic, v.nlev);
    if (shoc_top < 0) {
        return;
    }

    amrex::Real energy_before = amrex::Real(0.0);
    amrex::Real energy_after = amrex::Real(0.0);
    amrex::Real air_mass = amrex::Real(0.0);

    for (int k = 0; k <= shoc_top; ++k) {
        const amrex::Real mass = v.rho(ic,k,0) * v.dz(ic,k,0);
        const amrex::Real tabs_old = shoc::temperature_from_thetal(v.thetal_base(ic,k,0),
                                                                    v.qc_base(ic,k,0),
                                                                    v.qi_base(ic,k,0),
                                                                    v.exner(ic,k,0));

        amrex::Real tabs_new = amrex::Real(0.0);
        amrex::Real qv_new = amrex::Real(0.0);
        amrex::Real qc_new = amrex::Real(0.0);
        amrex::Real qi_new = amrex::Real(0.0);
        shoc::reconstruct_pdf_state(v.thetal(ic,k,0), v.qw(ic,k,0), v.exner(ic,k,0),
                                    v.qi_base(ic,k,0), v.shoc_ql(ic,k,0),
                                    tabs_new, qv_new, qc_new, qi_new);

        energy_before += mass * moist_energy(tabs_old, v.zt(ic,k,0),
                                             v.qv_base(ic,k,0), v.qc_base(ic,k,0), v.qi_base(ic,k,0),
                                             v.u_base(ic,k,0), v.v_base(ic,k,0), v.tke_base(ic,k,0));
        energy_after += mass * moist_energy(tabs_new, v.zt(ic,k,0),
                                            qv_new, qc_new, qi_new,
                                            v.u(ic,k,0), v.v(ic,k,0), v.tke(ic,k,0));
        air_mass += mass;
    }

    if (air_mass <= amrex::Real(0.0)) {
        return;
    }

    const amrex::Real latent_flux_coeff = L_v + constants::latent_ice();
    const amrex::Real energy_target = energy_before
        + dt * v.rho(ic,0,0)
        * (Cp_d * v.exner(ic,0,0) * v.surf_sens_flux(ic,0,0)
           + latent_flux_coeff * v.surf_lat_flux(ic,0,0));
    const amrex::Real delta_tabs = (energy_target - energy_after) / (Cp_d * air_mass);

    for (int k = 0; k <= shoc_top; ++k) {
        v.thetal(ic,k,0) += delta_tabs / amrex::max(v.exner(ic,k,0), amrex::Real(1.0e-12));
    }
}

} // namespace shoc

class ShocEnergyFixer
{
public:
    static int diagnose_active_top (const amrex::Vector<amrex::Real>& tke);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static int diagnose_active_top (const amrex::Array4<const amrex::Real>& tke,
                                    int ic,
                                    int nlev) noexcept
    {
        return shoc::diagnose_active_top(tke, ic, nlev);
    }

    static void apply_column (const ShocColumnData& col,
                              int ic,
                              amrex::Real dt,
                              const amrex::Vector<amrex::Real>& thl_old,
                              const amrex::Vector<amrex::Real>& qv_old,
                              const amrex::Vector<amrex::Real>& qc_old,
                              const amrex::Vector<amrex::Real>& qi_old,
                              const amrex::Vector<amrex::Real>& u_old,
                              const amrex::Vector<amrex::Real>& v_old,
                              const amrex::Vector<amrex::Real>& tke_old,
                              amrex::Vector<amrex::Real>& thl_new,
                              const amrex::Vector<amrex::Real>& qv_new,
                              const amrex::Vector<amrex::Real>& qc_new,
                              const amrex::Vector<amrex::Real>& qi_new,
                              const amrex::Vector<amrex::Real>& u_new,
                              const amrex::Vector<amrex::Real>& v_new,
                              const amrex::Vector<amrex::Real>& tke_new);
};

#endif
