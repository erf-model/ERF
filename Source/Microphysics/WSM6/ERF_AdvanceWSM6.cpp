#include "ERF_WSM6.H"
#include "ERF_WSM6_Fortran_Interface.H"
#include <cmath>

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
        const Box fab_box = mfi.fabbox();

        auto const& t_arr = mic_fab_vars[MicVar_WSM6::tabs]->array(mfi);
        auto const& qv_arr = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto const& qc_arr = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto const& qi_arr = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto const& qr_arr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto const& qs_arr = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto const& qg_arr = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);
        auto const& den_arr = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto const& p_arr = mic_fab_vars[MicVar_WSM6::pres]->array(mfi);
        auto rain_arr = mic_fab_vars[MicVar_WSM6::rain_accum]->array(mfi);
        auto snow_arr = mic_fab_vars[MicVar_WSM6::snow_accum]->array(mfi);
        auto graup_arr = mic_fab_vars[MicVar_WSM6::graup_accum]->array(mfi);

        const int ilo = box.smallEnd(0);
        const int ihi = box.bigEnd(0);
        const int jlo = box.smallEnd(1);
        const int jhi = box.bigEnd(1);
        const int klo = box.smallEnd(2);
        const int khi = box.bigEnd(2);

        const int imlo = fab_box.smallEnd(0);
        const int imhi = fab_box.bigEnd(0);
        const int jmlo = fab_box.smallEnd(1);
        const int jmhi = fab_box.bigEnd(1);
        const int kmlo = fab_box.smallEnd(2);
        const int kmhi = fab_box.bigEnd(2);

        const Real dz_val = m_geom.CellSize(2);
        FArrayBox delz_fab(fab_box, 1);
        auto const& delz_arr = delz_fab.array();
        delz_fab.setVal(dz_val);

        const Array4<const Real> z_arr = (m_z_phys_nd) ? m_z_phys_nd->const_array(mfi) : Array4<const Real> {};
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = (z_arr) ? Real(0.25) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                     + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                     + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                     + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz_val;
        });

        Box box2d(fab_box);
        box2d.makeSlab(2, 0);
        FArrayBox rainncv_fab(box2d, 1);
        FArrayBox sr_fab(box2d, 1);
        FArrayBox snowncv_fab(box2d, 1);
        FArrayBox graupelncv_fab(box2d, 1);
        FArrayBox rainacc_fab(box2d, 1);
        FArrayBox snowacc_fab(box2d, 1);
        FArrayBox graupacc_fab(box2d, 1);

        auto const& rainncv_arr = rainncv_fab.array();
        auto const& sr_arr = sr_fab.array();
        auto const& snowncv_arr = snowncv_fab.array();
        auto const& graupelncv_arr = graupelncv_fab.array();
        auto const& rainacc_arr = rainacc_fab.array();
        auto const& snowacc_arr = snowacc_fab.array();
        auto const& graupacc_arr = graupacc_fab.array();
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rainacc_arr(i,j,0) = rain_arr(i,j,klo);
            snowacc_arr(i,j,0) = snow_arr(i,j,klo);
            graupacc_arr(i,j,0) = graup_arr(i,j,klo);
            rainncv_arr(i,j,0) = Real(0.0);
            sr_arr(i,j,0) = Real(0.0);
            snowncv_arr(i,j,0) = Real(0.0);
            graupelncv_arr(i,j,0) = Real(0.0);
        });

        mp_wsm6_run_c(
            t_arr.dataPtr(),
            qv_arr.dataPtr(), qc_arr.dataPtr(), qi_arr.dataPtr(),
            qr_arr.dataPtr(), qs_arr.dataPtr(), qg_arr.dataPtr(),
            den_arr.dataPtr(), p_arr.dataPtr(), delz_arr.dataPtr(),
            static_cast<double>(dt), g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin,
            xls, xlv0, xlf0, den0, denr, cliq, cice, psat,
            rainacc_arr.dataPtr(), rainncv_arr.dataPtr(), sr_arr.dataPtr(),
            snowacc_arr.dataPtr(), snowncv_arr.dataPtr(),
            graupacc_arr.dataPtr(), graupelncv_arr.dataPtr(),
            imlo, imhi, jmlo, jmhi, kmlo, kmhi,
            ilo, ihi, jlo, jhi, klo, khi);

        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rain_arr(i,j,klo) = rainacc_arr(i,j,0);
            snow_arr(i,j,klo) = snowacc_arr(i,j,0);
            graup_arr(i,j,klo) = graupacc_arr(i,j,0);
        });
    }
#else
    amrex::Abort("WSM6 Fortran bridge requested but ERF was not built with ERF_USE_WSM6_FORT");
#endif
}
