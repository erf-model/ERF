#ifndef ERF_SHOC_COLUMN_DATA_H_
#define ERF_SHOC_COLUMN_DATA_H_

#include "ERF_ShocTypes.H"
#include "ERF_ShocGpuUtils.H"

#include <AMReX_Arena.H>
#include <AMReX_Box.H>

inline ShocColumnLayout
make_shoc_layout (const amrex::Box& bx, const amrex::Geometry& geom)
{
    ShocColumnLayout layout;
    layout.nx = bx.length(0);
    layout.ny = bx.length(1);
    layout.ncell = layout.nx * layout.ny;
    layout.nlev = geom.Domain().length(2);
    layout.imin = bx.smallEnd(0);
    layout.jmin = bx.smallEnd(1);
    layout.kmin = geom.Domain().smallEnd(2);
    layout.kmax = geom.Domain().bigEnd(2);
    return layout;
}

template <typename F>
inline void
for_each_shoc_fab (ShocColumnData& col, F&& f)
{
    for (amrex::FArrayBox* fab : {&col.zt,&col.zi,&col.dz,&col.p_mid,&col.p_int,&col.rho,&col.theta,
                                  &col.exner,&col.theta_v,&col.thetal,&col.qv,&col.qc,&col.qi,&col.qw,&col.tabs,&col.tke,
                                  &col.u,&col.v,&col.w,&col.host_dse,&col.pblh,&col.obklen,&col.ustar,
                                  &col.shoc_mix,&col.brunt,&col.isotropy,&col.tk,&col.tkh,&col.shear_prod,&col.buoy_prod,&col.diss_tke,&col.wthv_sec,
                                  &col.shoc_cldfrac,&col.shoc_ql,&col.shoc_ql2,&col.wqls_sec,&col.shoc_cond,&col.shoc_evap,
                                  &col.thl_sec,&col.qw_sec,&col.qwthl_sec,&col.wthl_sec,
                                  &col.wqw_sec,&col.uw_sec,&col.vw_sec,&col.wtke_sec,&col.w3,&col.w_sec,
                                  &col.thetal_base,&col.theta_base,&col.qv_base,&col.qc_base,&col.qi_base,
                                  &col.u_base,&col.v_base,&col.tke_base_state,
                                  &col.theta_tend,&col.qv_tend,&col.qc_tend,&col.qi_tend,&col.u_tend,&col.v_tend,&col.tke_tend,
                                  &col.surf_sens_flux,&col.surf_lat_flux,&col.surf_tau_u,&col.surf_tau_v}) {
        f(*fab);
    }
}

inline void
reset_shoc_column_data (ShocColumnData& col,
                        shoc::InitRunOn init_run_on = shoc::default_init_run_on())
{
    for_each_shoc_fab(col, [&] (amrex::FArrayBox& fab) {
        shoc::set_fab_val(fab, amrex::Real(0.0), init_run_on);
    });
}

inline void
define_shoc_column_data (ShocColumnData& col,
                         const ShocColumnLayout& layout,
                         amrex::Arena* arena = amrex::The_Async_Arena(),
                         shoc::InitRunOn init_run_on = shoc::default_init_run_on())
{
    col.layout = layout;

    const amrex::Box cell_box(amrex::IntVect(0,0,0),
                              amrex::IntVect(layout.ncell-1, layout.nlev-1, 0));
    const amrex::Box iface_box(amrex::IntVect(0,0,0),
                               amrex::IntVect(layout.ncell-1, layout.nlev, 0));
    const amrex::Box scal_box(amrex::IntVect(0,0,0),
                              amrex::IntVect(layout.ncell-1, 0, 0));

    auto make2 = [&] (amrex::FArrayBox& fab) { fab.resize(cell_box, 1, arena); };
    auto make2i = [&] (amrex::FArrayBox& fab) { fab.resize(iface_box, 1, arena); };
    auto make1 = [&] (amrex::FArrayBox& fab) { fab.resize(scal_box, 1, arena); };

    make2(col.zt); make2i(col.zi); make2(col.dz); make2(col.p_mid); make2i(col.p_int);
    make2(col.rho); make2(col.theta); make2(col.exner); make2(col.theta_v); make2(col.thetal);
    make2(col.qv); make2(col.qc); make2(col.qi); make2(col.qw); make2(col.tabs); make2(col.tke);
    make2(col.u); make2(col.v); make2(col.w); make2(col.host_dse);
    make1(col.pblh); make1(col.obklen); make1(col.ustar);
    make2(col.shoc_mix); make2(col.brunt); make2(col.isotropy); make2(col.tk); make2(col.tkh);
    make2(col.shear_prod); make2(col.buoy_prod); make2(col.diss_tke);
    make2(col.wthv_sec); make2(col.shoc_cldfrac); make2(col.shoc_ql); make2(col.shoc_ql2);
    make2(col.wqls_sec); make2(col.shoc_cond); make2(col.shoc_evap); make2i(col.thl_sec); make2i(col.qw_sec);
    make2i(col.qwthl_sec); make2i(col.wthl_sec); make2i(col.wqw_sec); make2i(col.uw_sec);
    make2i(col.vw_sec); make2i(col.wtke_sec); make2i(col.w3); make2(col.w_sec);
    make2(col.thetal_base); make2(col.theta_base); make2(col.qv_base); make2(col.qc_base); make2(col.qi_base);
    make2(col.u_base); make2(col.v_base); make2(col.tke_base_state);
    make2(col.theta_tend); make2(col.qv_tend); make2(col.qc_tend); make2(col.qi_tend); make2(col.u_tend);
    make2(col.v_tend); make2(col.tke_tend);
    make1(col.surf_sens_flux); make1(col.surf_lat_flux); make1(col.surf_tau_u); make1(col.surf_tau_v);

    reset_shoc_column_data(col, init_run_on);
}

inline bool
shoc_layout_fits (const ShocColumnLayout& active,
                  const ShocColumnLayout& capacity)
{
    return capacity.ncell >= active.ncell &&
           capacity.nlev >= active.nlev &&
           capacity.kmin == active.kmin &&
           capacity.kmax == active.kmax;
}

class ShocColumnWorkspace
{
public:
    ShocColumnData col;

    void ensure_capacity (const ShocColumnLayout& layout,
                          amrex::Arena* arena = amrex::The_Async_Arena(),
                          shoc::InitRunOn init_run_on = shoc::default_init_run_on())
    {
        if (!m_defined || arena != m_arena || !shoc_layout_fits(layout, m_capacity_layout)) {
            m_capacity_layout = layout;
            define_shoc_column_data(col, m_capacity_layout, arena, init_run_on);
            m_defined = true;
            m_arena = arena;
        } else {
            col.layout = layout;
            reset_shoc_column_data(col, init_run_on);
        }
    }

private:
    ShocColumnLayout m_capacity_layout{};
    amrex::Arena* m_arena = nullptr;
    bool m_defined = false;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
shoc_column_index (const ShocColumnLayout& layout, int i, int j)
{
    return (j - layout.jmin) * layout.nx + (i - layout.imin);
}

#endif
