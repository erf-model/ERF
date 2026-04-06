#include "ERF_WSM6.H"
#include "ERF_WSM6_Fortran_Interface.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace amrex;

#ifdef ERF_USE_WSM6_FORT
namespace {

struct MPDbgVarStats {
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    double l1 = 0.0;
    long long n = 0;
};

struct MPDbgStats {
    MPDbgVarStats qrain;
    MPDbgVarStats qc;
    MPDbgVarStats rhoQ2;
    MPDbgVarStats qsnow;
    MPDbgVarStats qgraup;
    MPDbgVarStats zvel;
};

inline void
update_stat (MPDbgVarStats& s, double v)
{
    s.min = std::min(s.min, v);
    s.max = std::max(s.max, v);
    s.l1 += std::abs(v);
    s.n += 1;
}

inline void
reduce_stat (MPDbgVarStats& s)
{
    ParallelDescriptor::ReduceRealMin(s.min);
    ParallelDescriptor::ReduceRealMax(s.max);
    ParallelDescriptor::ReduceRealSum(s.l1);
    amrex::Long n = static_cast<amrex::Long>(s.n);
    ParallelDescriptor::ReduceLongSum(n);
    s.n = static_cast<long long>(n);
}

inline void
reduce_stats (MPDbgStats& s)
{
    reduce_stat(s.qrain);
    reduce_stat(s.qc);
    reduce_stat(s.rhoQ2);
    reduce_stat(s.qsnow);
    reduce_stat(s.qgraup);
    reduce_stat(s.zvel);
}

inline void
print_stats (const char* tag, const char* model, int step, const MPDbgStats& s)
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::printf(
        "%s %s step=%04d "
        "qrain_min=%.12e qrain_max=%.12e qrain_L1=%.12e "
        "qc_min=%.12e qc_max=%.12e qc_L1=%.12e "
        "rhoQ2_min=%.12e rhoQ2_max=%.12e rhoQ2_L1=%.12e "
        "qsnow_min=%.12e qsnow_max=%.12e qsnow_L1=%.12e "
        "qgraup_min=%.12e qgraup_max=%.12e qgraup_L1=%.12e "
        "z_velocity_min=%.12e z_velocity_max=%.12e z_velocity_L1=%.12e\n",
        tag, model, step,
        s.qrain.min, s.qrain.max, s.qrain.l1,
        s.qc.min, s.qc.max, s.qc.l1,
        s.rhoQ2.min, s.rhoQ2.max, s.rhoQ2.l1,
        s.qsnow.min, s.qsnow.max, s.qsnow.l1,
        s.qgraup.min, s.qgraup.max, s.qgraup.l1,
        s.zvel.min, s.zvel.max, s.zvel.l1);
    std::fflush(stdout);
}

inline void
mark_zvel_nan (MPDbgVarStats& s)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    s.min = nan;
    s.max = nan;
    s.l1 = nan;
    s.n = 1;
}

} // namespace
#endif

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

    static int mpdbg_step = 0;
    ParmParse pp("erf");
    int microphysics_debug = 0;
    pp.query("microphysics_debug", microphysics_debug);
    microphysics_debug = std::max(0, std::min(2, microphysics_debug));
    int microphysics_debug_seam = 0;
    pp.query("microphysics_debug_seam", microphysics_debug_seam);
    int seam_i_lo = 94;
    int seam_i_hi = 97;
    int seam_j_lo = 0;
    int seam_j_hi = 0;
    int seam_k_lo = 0;
    int seam_k_hi = 38;
    pp.query("microphysics_debug_seam_i_lo", seam_i_lo);
    pp.query("microphysics_debug_seam_i_hi", seam_i_hi);
    pp.query("microphysics_debug_seam_j_lo", seam_j_lo);
    pp.query("microphysics_debug_seam_j_hi", seam_j_hi);
    pp.query("microphysics_debug_seam_k_lo", seam_k_lo);
    pp.query("microphysics_debug_seam_k_hi", seam_k_hi);

    MPDbgStats pre_stats;
    MPDbgStats post_stats;

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

        const Real dz_val = m_geom.CellSize(m_axis);
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

        if (microphysics_debug > 0) {
            for (int k = klo; k <= khi; ++k) {
                for (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        const double qr = qr_arr(i,j,k);
                        const double qc = qc_arr(i,j,k);
                        const double rhoqc = den_arr(i,j,k) * qc;
                        const double qs = qs_arr(i,j,k);
                        const double qg = qg_arr(i,j,k);
                        update_stat(pre_stats.qrain, qr);
                        update_stat(pre_stats.qc, qc);
                        update_stat(pre_stats.rhoQ2, rhoqc);
                        update_stat(pre_stats.qsnow, qs);
                        update_stat(pre_stats.qgraup, qg);
                    }
                }
            }
        }

        if (microphysics_debug > 1) {
            const int rank = ParallelDescriptor::MyProc();
            for (int k = klo; k <= khi; ++k) {
                for (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        std::printf(
                            "MPDBG PRE WSM6 rank=%d step=%04d i=%d j=%d k=%d "
                            "qrain=%.12e qc=%.12e rhoQ2=%.12e qsnow=%.12e qgraup=%.12e\n",
                            rank, mpdbg_step, i, j, k,
                            static_cast<double>(qr_arr(i,j,k)),
                            static_cast<double>(qc_arr(i,j,k)),
                            static_cast<double>(den_arr(i,j,k) * qc_arr(i,j,k)),
                            static_cast<double>(qs_arr(i,j,k)),
                            static_cast<double>(qg_arr(i,j,k)));
                    }
                }
            }
            std::fflush(stdout);
        }

        if (microphysics_debug_seam > 0) {
            const int rank = ParallelDescriptor::MyProc();
            const int si_lo = std::max(imlo, seam_i_lo);
            const int si_hi = std::min(imhi, seam_i_hi);
            const int sj_lo = std::max(jmlo, seam_j_lo);
            const int sj_hi = std::min(jmhi, seam_j_hi);
            const int sk_lo = std::max(kmlo, seam_k_lo);
            const int sk_hi = std::min(kmhi, seam_k_hi);
            for (int k = sk_lo; k <= sk_hi; ++k) {
                for (int j = sj_lo; j <= sj_hi; ++j) {
                    for (int i = si_lo; i <= si_hi; ++i) {
                        const int in_valid = (i >= ilo && i <= ihi &&
                                              j >= jlo && j <= jhi &&
                                              k >= klo && k <= khi) ? 1 : 0;
                        std::printf(
                            "MPDBG SEAM PRE WSM6 rank=%d step=%04d i=%d j=%d k=%d valid=%d "
                            "qv=%.12e qc=%.12e qi=%.12e qrain=%.12e qsnow=%.12e qgraup=%.12e\n",
                            rank, mpdbg_step, i, j, k, in_valid,
                            static_cast<double>(qv_arr(i,j,k)),
                            static_cast<double>(qc_arr(i,j,k)),
                            static_cast<double>(qi_arr(i,j,k)),
                            static_cast<double>(qr_arr(i,j,k)),
                            static_cast<double>(qs_arr(i,j,k)),
                            static_cast<double>(qg_arr(i,j,k)));
                    }
                }
            }
            std::fflush(stdout);
        }

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

        if (microphysics_debug > 0) {
            for (int k = klo; k <= khi; ++k) {
                for (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        const double qr = qr_arr(i,j,k);
                        const double qc = qc_arr(i,j,k);
                        const double rhoqc = den_arr(i,j,k) * qc;
                        const double qs = qs_arr(i,j,k);
                        const double qg = qg_arr(i,j,k);
                        update_stat(post_stats.qrain, qr);
                        update_stat(post_stats.qc, qc);
                        update_stat(post_stats.rhoQ2, rhoqc);
                        update_stat(post_stats.qsnow, qs);
                        update_stat(post_stats.qgraup, qg);
                    }
                }
            }
        }

        if (microphysics_debug > 1) {
            const int rank = ParallelDescriptor::MyProc();
            for (int k = klo; k <= khi; ++k) {
                for (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {
                        std::printf(
                            "MPDBG POST WSM6 rank=%d step=%04d i=%d j=%d k=%d "
                            "qrain=%.12e qc=%.12e rhoQ2=%.12e qsnow=%.12e qgraup=%.12e\n",
                            rank, mpdbg_step, i, j, k,
                            static_cast<double>(qr_arr(i,j,k)),
                            static_cast<double>(qc_arr(i,j,k)),
                            static_cast<double>(den_arr(i,j,k) * qc_arr(i,j,k)),
                            static_cast<double>(qs_arr(i,j,k)),
                            static_cast<double>(qg_arr(i,j,k)));
                    }
                }
            }
            std::fflush(stdout);
        }

        if (microphysics_debug_seam > 0) {
            const int rank = ParallelDescriptor::MyProc();
            const int si_lo = std::max(imlo, seam_i_lo);
            const int si_hi = std::min(imhi, seam_i_hi);
            const int sj_lo = std::max(jmlo, seam_j_lo);
            const int sj_hi = std::min(jmhi, seam_j_hi);
            const int sk_lo = std::max(kmlo, seam_k_lo);
            const int sk_hi = std::min(kmhi, seam_k_hi);
            for (int k = sk_lo; k <= sk_hi; ++k) {
                for (int j = sj_lo; j <= sj_hi; ++j) {
                    for (int i = si_lo; i <= si_hi; ++i) {
                        const int in_valid = (i >= ilo && i <= ihi &&
                                              j >= jlo && j <= jhi &&
                                              k >= klo && k <= khi) ? 1 : 0;
                        std::printf(
                            "MPDBG SEAM POST WSM6 rank=%d step=%04d i=%d j=%d k=%d valid=%d "
                            "qv=%.12e qc=%.12e qi=%.12e qrain=%.12e qsnow=%.12e qgraup=%.12e\n",
                            rank, mpdbg_step, i, j, k, in_valid,
                            static_cast<double>(qv_arr(i,j,k)),
                            static_cast<double>(qc_arr(i,j,k)),
                            static_cast<double>(qi_arr(i,j,k)),
                            static_cast<double>(qr_arr(i,j,k)),
                            static_cast<double>(qs_arr(i,j,k)),
                            static_cast<double>(qg_arr(i,j,k)));
                    }
                }
            }
            std::fflush(stdout);
        }
    }

    if (microphysics_debug > 0) {
        mark_zvel_nan(pre_stats.zvel);
        mark_zvel_nan(post_stats.zvel);
        reduce_stats(pre_stats);
        reduce_stats(post_stats);
        print_stats("MPDBG PRE", "WSM6", mpdbg_step, pre_stats);
        print_stats("MPDBG POST", "WSM6", mpdbg_step, post_stats);
        print_stats("MPDBG STEP", "WSM6", mpdbg_step, post_stats);
    }
    ++mpdbg_step;
#else
    amrex::Abort("WSM6 Fortran bridge requested but ERF was not built with ERF_USE_WSM6_FORT");
#endif
}
