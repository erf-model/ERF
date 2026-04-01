#include "ERF_WSM6.H"
#include "ERF_WSM6_Fortran_Interface.H"

using namespace amrex;

void
WSM6::Advance(const Real& dt_advance,
              const SolverChoice&)
{
    dt = dt_advance;

#ifdef ERF_USE_WSM6_FORT
    static bool wsm6_inited = false;

    // Minimal phase-1 initialization for single-moment WSM6.
    if (!wsm6_inited) {
        constexpr double den0 = 1.28;                 // Standard dry-air density (kg/m^3)
        constexpr double denr = static_cast<double>(rhoh2o);
        constexpr double dens = static_cast<double>(rhos);
        constexpr double cl = static_cast<double>(Cp_l);
        constexpr double cpv = static_cast<double>(Cp_v);
        constexpr int hail_opt = 0;                   // Graupel mode
        mp_wsm6_init_c(den0, denr, dens, cl, cpv, hail_opt);
        wsm6_inited = true;
    }

    const auto domain = m_geom.Domain();
    const int i_lo = domain.smallEnd(0);
    const int i_hi = domain.bigEnd(0);
    const int j_lo = domain.smallEnd(1);
    const int j_hi = domain.bigEnd(1);
    const int k_lo = domain.smallEnd(2);
    const int k_hi = domain.bigEnd(2);

    constexpr double g = static_cast<double>(CONST_GRAV);
    constexpr double cpd = static_cast<double>(Cp_d);
    constexpr double cpv = static_cast<double>(Cp_v);
    constexpr double rd = static_cast<double>(R_d);
    constexpr double rv = static_cast<double>(R_v);
    constexpr double t0c = 273.15;
    constexpr double ep1 = static_cast<double>(R_v / R_d - one);
    constexpr double ep2 = static_cast<double>(R_d / R_v);
    constexpr double qmin = 1.0e-12;
    constexpr double xls = static_cast<double>(lsub);
    constexpr double xlv0 = static_cast<double>(lat_vap);
    constexpr double xlf0 = static_cast<double>(lat_ice);
    constexpr double den0 = 1.28;
    constexpr double denr = static_cast<double>(rhoh2o);
    constexpr double cliq = static_cast<double>(Cp_l);
    constexpr double cice = 2106.0;
    constexpr double psat = 610.78;

    for (MFIter mfi(*mic_fab_vars[MicVar_WSM6::qv], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();

        auto const& t_arr = mic_fab_vars[MicVar_WSM6::tabs]->array(mfi);
        auto const& qv_arr = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto const& qc_arr = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto const& qi_arr = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto const& qr_arr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto const& qs_arr = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto const& qg_arr = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);
        auto const& den_arr = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto const& p_arr = mic_fab_vars[MicVar_WSM6::pres]->array(mfi);
        auto const& rain_arr = mic_fab_vars[MicVar_WSM6::rain_accum]->array(mfi);
        auto const& snow_arr = mic_fab_vars[MicVar_WSM6::snow_accum]->array(mfi);
        auto const& graup_arr = mic_fab_vars[MicVar_WSM6::graup_accum]->array(mfi);

        const int ilo = box.smallEnd(0);
        const int ihi = box.bigEnd(0);
        const int jlo = box.smallEnd(1);
        const int jhi = box.bigEnd(1);
        const int klo = box.smallEnd(2);
        const int khi = box.bigEnd(2);

        FArrayBox delz_fab(box, 1);
        auto const& delz_arr = delz_fab.array();

        const Real dz_val = m_geom.CellSize(m_axis);
        const Array4<const Real> z_arr = (m_z_phys_nd) ? m_z_phys_nd->const_array(mfi) : Array4<const Real> {};
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = (z_arr) ? Real(0.25) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                     + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                     + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                     + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz_val;
        });

        Box box2d(box);
        box2d.makeSlab(2, 0);
        FArrayBox rainncv_fab(box2d, 1);
        FArrayBox sr_fab(box2d, 1);
        FArrayBox snowncv_fab(box2d, 1);
        FArrayBox graupelncv_fab(box2d, 1);

        auto const& rainncv_arr = rainncv_fab.array();
        auto const& sr_arr = sr_fab.array();
        auto const& snowncv_arr = snowncv_fab.array();
        auto const& graupelncv_arr = graupelncv_fab.array();
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rainncv_arr(i,j,0) = Real(0.0);
            sr_arr(i,j,0) = Real(0.0);
            snowncv_arr(i,j,0) = Real(0.0);
            graupelncv_arr(i,j,0) = Real(0.0);
        });

        mp_wsm6_run_c(
            &(t_arr(ilo,jlo,klo)),
            &(qv_arr(ilo,jlo,klo)), &(qc_arr(ilo,jlo,klo)), &(qi_arr(ilo,jlo,klo)),
            &(qr_arr(ilo,jlo,klo)), &(qs_arr(ilo,jlo,klo)), &(qg_arr(ilo,jlo,klo)),
            &(den_arr(ilo,jlo,klo)), &(p_arr(ilo,jlo,klo)), &(delz_arr(ilo,jlo,klo)),
            static_cast<double>(dt), g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin,
            xls, xlv0, xlf0, den0, denr, cliq, cice, psat,
            &(rain_arr(ilo,jlo,0)), &(rainncv_arr(ilo,jlo,0)), &(sr_arr(ilo,jlo,0)),
            &(snow_arr(ilo,jlo,0)), &(snowncv_arr(ilo,jlo,0)),
            &(graup_arr(ilo,jlo,0)), &(graupelncv_arr(ilo,jlo,0)),
            i_lo, i_hi, j_lo, j_hi, k_lo, k_hi,
            i_lo, i_hi, j_lo, j_hi, k_lo, k_hi,
            ilo, ihi, jlo, jhi, klo, khi);
    }
#else
    amrex::Abort("WSM6 Fortran bridge requested but ERF was not built with ERF_USE_WSM6_FORT");
#endif
}
