#include "ERF_WSM6.H"
#include "ERF_WSM6_Fortran_Interface.H"
#include <cmath>

using namespace amrex;

// ---------------------------------------------------------------
// WSM6 device-callable free functions (Rule 16, Rule 18)
// Statement functions from mp_wsm6_run declaration section
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_cpmcal (Real x, Real qmin_arg, Real cpd_arg, Real cpv_arg) {
    return cpd_arg*(Real(1.0)-amrex::max(x,qmin_arg))
          +amrex::max(x,qmin_arg)*cpv_arg;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_xlcal (Real x, Real xlv0_arg, Real xlv1_arg, Real t0c_arg) {
    return xlv0_arg - xlv1_arg*(x - t0c_arg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_diffus (Real x, Real y) {
    return Real(8.794e-5)*std::exp(std::log(x)*Real(1.81))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_viscos (Real x, Real y) {
    return Real(1.496e-6)*(x*std::sqrt(x))/(x+Real(120.0))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_xka (Real x, Real y) {
    return Real(1.414e3)*wsm6_viscos(x,y)*y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_diffac (Real a, Real b, Real c, Real d, Real e,
                   Real rv_arg) {
    return d*a*a/(wsm6_xka(c,d)*rv_arg*c*c)
          +Real(1.0)/(e*wsm6_diffus(c,b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_venfac (Real a, Real b, Real c, Real den0_arg) {
    return std::exp(std::log(wsm6_viscos(b,c)/wsm6_diffus(b,a))
                   *Real(0.3333333))
          /std::sqrt(wsm6_viscos(b,c))
          *std::sqrt(std::sqrt(den0_arg/c));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_conden (Real a, Real b, Real c, Real d, Real e,
                   Real qmin_arg, Real rv_arg) {
    return (amrex::max(b,qmin_arg)-c)
          /(Real(1.0)+d*d/(rv_arg*e)*c/(a*a));
}

// ---------------------------------------------------------------
// Slope parameter lambda functions (statement functions from
// slope_wsm6, slope_rain, slope_snow, slope_graup)
// pidn0r, pidn0s, pidn0g are class constexpr members from ERF_WSM6.H
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_lamdar (Real x, Real y, Real pidn0r_arg) {
    return std::sqrt(std::sqrt(pidn0r_arg/(x*y)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_lamdas (Real x, Real y, Real z, Real pidn0s_arg) {
    return std::sqrt(std::sqrt(pidn0s_arg*z/(x*y)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wsm6_lamdag (Real x, Real y, Real pidn0g_arg) {
    return std::sqrt(std::sqrt(pidn0g_arg/(x*y)));
}

// ---------------------------------------------------------------
// Full slope subroutine device functions (Rule 18)
// Each takes single-cell scalar inputs, returns slope params
// by reference — loop over (i,j,k) is provided by ParallelFor
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_slope_rain_cell (Real qr, Real den, Real denfac,
                            Real pidn0r_arg,
                            Real qcrmin_arg, Real rslopermax_arg,
                            Real rsloperbmax_arg, Real rsloper2max_arg,
                            Real rsloper3max_arg, Real bvtr_arg,
                            Real pvtr_arg,
                            Real& rslope, Real& rslopeb,
                            Real& rslope2, Real& rslope3, Real& vt)
{
    if (qr <= qcrmin_arg) {
        rslope  = rslopermax_arg;
        rslopeb = rsloperbmax_arg;
        rslope2 = rsloper2max_arg;
        rslope3 = rsloper3max_arg;
    } else {
        rslope  = Real(1.0)/wsm6_lamdar(qr,den,pidn0r_arg);
        rslopeb = std::pow(rslope,bvtr_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvtr_arg*rslopeb*denfac;
    if (qr <= Real(0.0)) vt = Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_slope_snow_cell (Real qs, Real den, Real denfac, Real t,
                            Real pidn0s_arg, Real alpha_arg,
                            Real n0smax_arg, Real n0s_arg,
                            Real t0c_arg, Real qcrmin_arg,
                            Real rslopesmax_arg, Real rslopesbmax_arg,
                            Real rslopes2max_arg, Real rslopes3max_arg,
                            Real bvts_arg, Real pvts_arg,
                            Real& rslope, Real& rslopeb,
                            Real& rslope2, Real& rslope3, Real& vt,
                            Real& n0sfac)
{
    Real supcol = t0c_arg - t;
    n0sfac = amrex::max(amrex::min(std::exp(alpha_arg*supcol),
                                    n0smax_arg/n0s_arg), Real(1.0));
    if (qs <= qcrmin_arg) {
        rslope  = rslopesmax_arg;
        rslopeb = rslopesbmax_arg;
        rslope2 = rslopes2max_arg;
        rslope3 = rslopes3max_arg;
    } else {
        rslope  = Real(1.0)/wsm6_lamdas(qs,den,n0sfac,pidn0s_arg);
        rslopeb = std::pow(rslope,bvts_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvts_arg*rslopeb*denfac;
    if (qs <= Real(0.0)) vt = Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wsm6_slope_graup_cell (Real qg, Real den, Real denfac,
                             Real pidn0g_arg, Real qcrmin_arg,
                             Real rslopegmax_arg, Real rslopegbmax_arg,
                             Real rslopeg2max_arg, Real rslopeg3max_arg,
                             Real bvtg_arg, Real pvtg_arg,
                             Real& rslope, Real& rslopeb,
                             Real& rslope2, Real& rslope3, Real& vt)
{
    if (qg <= qcrmin_arg) {
        rslope  = rslopegmax_arg;
        rslopeb = rslopegbmax_arg;
        rslope2 = rslopeg2max_arg;
        rslope3 = rslopeg3max_arg;
    } else {
        rslope  = Real(1.0)/wsm6_lamdag(qg,den,pidn0g_arg);
        rslopeb = std::pow(rslope,bvtg_arg);
        rslope2 = rslope*rslope;
        rslope3 = rslope2*rslope;
    }
    vt = pvtg_arg*rslopeb*denfac;
    if (qg <= Real(0.0)) vt = Real(0.0);
}

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
#endif

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

#ifdef ERF_USE_WSM6_FORT
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
#else
        // --- Phase 4 native C++ kernel ---

        // box2d for 1D per-column arrays (already defined above)
        // delqrs1/2/3, delqi: surface precipitation flux accumulators
        FArrayBox delqrs1_fab(box2d,1); delqrs1_fab.setVal(Real(0.0));
        FArrayBox delqrs2_fab(box2d,1); delqrs2_fab.setVal(Real(0.0));
        FArrayBox delqrs3_fab(box2d,1); delqrs3_fab.setVal(Real(0.0));
        FArrayBox delqi_fab(box2d,1);   delqi_fab.setVal(Real(0.0));
        FArrayBox tstepsnow_fab(box2d,1);  tstepsnow_fab.setVal(Real(0.0));
        FArrayBox tstepgraup_fab(box2d,1); tstepgraup_fab.setVal(Real(0.0));
        auto const& delqrs1_arr   = delqrs1_fab.array();
        auto const& delqrs2_arr   = delqrs2_fab.array();
        auto const& delqrs3_arr   = delqrs3_fab.array();
        auto const& delqi_arr     = delqi_fab.array();
        auto const& tstepsnow_arr = tstepsnow_fab.array();
        auto const& tstepgraup_arr= tstepgraup_fab.array();

        // 3D working FABs
        FArrayBox denfac_fab(fab_box,1);  FArrayBox xni_fab(fab_box,1);
        FArrayBox cpm_fab(fab_box,1);     FArrayBox xl_fab(fab_box,1);
        FArrayBox qsatw_fab(fab_box,1);   FArrayBox qsati_fab(fab_box,1);
        FArrayBox rhw_fab(fab_box,1);     FArrayBox rhi_fab(fab_box,1);
        FArrayBox den_tmp_fab(fab_box,1); FArrayBox delz_tmp_fab(fab_box,1);
        FArrayBox n0sfac_fab(fab_box,1);
        FArrayBox qrs_tmp_r_fab(fab_box,1); FArrayBox qrs_tmp_s_fab(fab_box,1);
        FArrayBox qrs_tmp_g_fab(fab_box,1);
        FArrayBox rslope_r_fab(fab_box,1);  FArrayBox rslope_s_fab(fab_box,1);
        FArrayBox rslope_g_fab(fab_box,1);
        FArrayBox rslope2_r_fab(fab_box,1); FArrayBox rslope2_s_fab(fab_box,1);
        FArrayBox rslope2_g_fab(fab_box,1);
        FArrayBox rslope3_r_fab(fab_box,1); FArrayBox rslope3_s_fab(fab_box,1);
        FArrayBox rslope3_g_fab(fab_box,1);
        FArrayBox rslopeb_r_fab(fab_box,1); FArrayBox rslopeb_s_fab(fab_box,1);
        FArrayBox rslopeb_g_fab(fab_box,1);
        FArrayBox work1_r_fab(fab_box,1);   FArrayBox work1_s_fab(fab_box,1);
        FArrayBox work1_g_fab(fab_box,1);
        FArrayBox work2_fab(fab_box,1);     FArrayBox workdiffw_fab(fab_box,1);
        FArrayBox workdiffi_fab(fab_box,1);
        FArrayBox workr_fab(fab_box,1);     FArrayBox worka_fab(fab_box,1);
        FArrayBox work1c_fab(fab_box,1);
        FArrayBox denqrs1_fab(fab_box,1);   FArrayBox denqrs2_fab(fab_box,1);
        FArrayBox denqrs3_fab(fab_box,1);   FArrayBox denqci_fab(fab_box,1);
        FArrayBox fall_r_fab(fab_box,1);    FArrayBox fall_s_fab(fab_box,1);
        FArrayBox fall_g_fab(fab_box,1);    FArrayBox fallc_fab(fab_box,1);
        FArrayBox qsum_fab(fab_box,1);
        // process rates
        FArrayBox praut_fab(fab_box,1); FArrayBox pracw_fab(fab_box,1);
        FArrayBox prevp_fab(fab_box,1); FArrayBox psdep_fab(fab_box,1);
        FArrayBox pgdep_fab(fab_box,1); FArrayBox psaut_fab(fab_box,1);
        FArrayBox pgaut_fab(fab_box,1); FArrayBox praci_fab(fab_box,1);
        FArrayBox piacr_fab(fab_box,1); FArrayBox psaci_fab(fab_box,1);
        FArrayBox psacw_fab(fab_box,1); FArrayBox pgacw_fab(fab_box,1);
        FArrayBox pgaci_fab(fab_box,1); FArrayBox paacw_fab(fab_box,1);
        FArrayBox pracs_fab(fab_box,1); FArrayBox psacr_fab(fab_box,1);
        FArrayBox pgacr_fab(fab_box,1); FArrayBox pgacs_fab(fab_box,1);
        FArrayBox pigen_fab(fab_box,1); FArrayBox pidep_fab(fab_box,1);
        FArrayBox pcond_fab(fab_box,1); FArrayBox psmlt_fab(fab_box,1);
        FArrayBox pgmlt_fab(fab_box,1); FArrayBox pseml_fab(fab_box,1);
        FArrayBox pgeml_fab(fab_box,1); FArrayBox psevp_fab(fab_box,1);
        FArrayBox pgevp_fab(fab_box,1);

        auto const& denfac_arr    = denfac_fab.array();
        auto const& xni_arr       = xni_fab.array();
        auto const& cpm_arr       = cpm_fab.array();
        auto const& xl_arr        = xl_fab.array();
        auto const& qsatw_arr     = qsatw_fab.array();
        auto const& qsati_arr     = qsati_fab.array();
        auto const& rhw_arr       = rhw_fab.array();
        auto const& rhi_arr       = rhi_fab.array();
        auto const& den_tmp_arr   = den_tmp_fab.array();
        auto const& delz_tmp_arr  = delz_tmp_fab.array();
        auto const& n0sfac_arr    = n0sfac_fab.array();
        auto const& qrs_tmp_r_arr = qrs_tmp_r_fab.array();
        auto const& qrs_tmp_s_arr = qrs_tmp_s_fab.array();
        auto const& qrs_tmp_g_arr = qrs_tmp_g_fab.array();
        auto const& rslope_r_arr  = rslope_r_fab.array();
        auto const& rslope_s_arr  = rslope_s_fab.array();
        auto const& rslope_g_arr  = rslope_g_fab.array();
        auto const& rslope2_r_arr = rslope2_r_fab.array();
        auto const& rslope2_s_arr = rslope2_s_fab.array();
        auto const& rslope2_g_arr = rslope2_g_fab.array();
        auto const& rslope3_r_arr = rslope3_r_fab.array();
        auto const& rslope3_s_arr = rslope3_s_fab.array();
        auto const& rslope3_g_arr = rslope3_g_fab.array();
        auto const& rslopeb_r_arr = rslopeb_r_fab.array();
        auto const& rslopeb_s_arr = rslopeb_s_fab.array();
        auto const& rslopeb_g_arr = rslopeb_g_fab.array();
        auto const& work1_r_arr   = work1_r_fab.array();
        auto const& work1_s_arr   = work1_s_fab.array();
        auto const& work1_g_arr   = work1_g_fab.array();
        auto const& work2_arr     = work2_fab.array();
        auto const& workdiffw_arr = workdiffw_fab.array();
        auto const& workdiffi_arr = workdiffi_fab.array();
        auto const& workr_arr     = workr_fab.array();
        auto const& worka_arr     = worka_fab.array();
        auto const& work1c_arr    = work1c_fab.array();
        auto const& denqrs1_arr   = denqrs1_fab.array();
        auto const& denqrs2_arr   = denqrs2_fab.array();
        auto const& denqrs3_arr   = denqrs3_fab.array();
        auto const& denqci_arr    = denqci_fab.array();
        auto const& fall_r_arr    = fall_r_fab.array();
        auto const& fall_s_arr    = fall_s_fab.array();
        auto const& fall_g_arr    = fall_g_fab.array();
        auto const& fallc_arr     = fallc_fab.array();
        auto const& qsum_arr      = qsum_fab.array();
        auto const& praut_arr     = praut_fab.array();
        auto const& pracw_arr     = pracw_fab.array();
        auto const& prevp_arr     = prevp_fab.array();
        auto const& psdep_arr     = psdep_fab.array();
        auto const& pgdep_arr     = pgdep_fab.array();
        auto const& psaut_arr     = psaut_fab.array();
        auto const& pgaut_arr     = pgaut_fab.array();
        auto const& praci_arr     = praci_fab.array();
        auto const& piacr_arr     = piacr_fab.array();
        auto const& psaci_arr     = psaci_fab.array();
        auto const& psacw_arr     = psacw_fab.array();
        auto const& pgacw_arr     = pgacw_fab.array();
        auto const& pgaci_arr     = pgaci_fab.array();
        auto const& paacw_arr     = paacw_fab.array();
        auto const& pracs_arr     = pracs_fab.array();
        auto const& psacr_arr     = psacr_fab.array();
        auto const& pgacr_arr     = pgacr_fab.array();
        auto const& pgacs_arr     = pgacs_fab.array();
        auto const& pigen_arr     = pigen_fab.array();
        auto const& pidep_arr     = pidep_fab.array();
        auto const& pcond_arr     = pcond_fab.array();
        auto const& psmlt_arr     = psmlt_fab.array();
        auto const& pgmlt_arr     = pgmlt_fab.array();
        auto const& pseml_arr     = pseml_fab.array();
        auto const& pgeml_arr     = pgeml_fab.array();
        auto const& psevp_arr     = psevp_fab.array();
        auto const& pgevp_arr     = pgevp_fab.array();

        // Groups A-E: pre-loop setup
        // Clamp negative values (Group A)
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k), Real(0.0));
            qr_arr(i,j,k) = amrex::max(qr_arr(i,j,k), Real(0.0));
            qi_arr(i,j,k) = amrex::max(qi_arr(i,j,k), Real(0.0));
            qs_arr(i,j,k) = amrex::max(qs_arr(i,j,k), Real(0.0));
            qg_arr(i,j,k) = amrex::max(qg_arr(i,j,k), Real(0.0));
            den_tmp_arr(i,j,k)  = den_arr(i,j,k);
            delz_tmp_arr(i,j,k) = delz_arr(i,j,k);
        });

        // Group B: cpm, xl — computed once from initial state [lines 455-460]
        const Real xlv1_loc = m_xlv1;
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            cpm_arr(i,j,k) = wsm6_cpmcal(qv_arr(i,j,k), Real(qmin), Real(cpd), Real(cpv));
            xl_arr(i,j,k)  = wsm6_xlcal(t_arr(i,j,k), Real(xlv0), xlv1_loc, Real(t0c));
        });

        // Outer minor timestep loop (Rule 29)
        const int wsm6_loops = std::max(
            static_cast<int>(std::round(dt / dtcldcr)), 1);
        const Real dtcld = dt / static_cast<Real>(wsm6_loops);

        for (int loop = 0; loop < wsm6_loops; ++loop) {
            // G1b: denfac = sqrt(den0/den)  [lines 503-515]
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                denfac_arr(i,j,k) = std::sqrt(Real(den0)/den_arr(i,j,k));
            });

            // G1c: qsatw, qsati, rhw, rhi  [lines 517-549]
            {
                const Real ttp  = Real(t0c) + Real(0.01);
                const Real dldt = Real(cpv) - Real(cliq);
                const Real xa   = -dldt / Real(rv);
                const Real xb   =  xa + Real(xlv0) / (Real(rv)*ttp);
                const Real dldti= Real(cpv) - Real(cice);
                const Real xai  = -dldti / Real(rv);
                const Real xbi  =  xai + Real(xls) / (Real(rv)*ttp);
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real tr = ttp / t_arr(i,j,k);
                    Real qsw = Real(psat)*std::exp(std::log(tr)*xa)*std::exp(xb*(Real(1.0)-tr));
                    qsw = amrex::min(qsw, Real(0.99)*p_arr(i,j,k));
                    qsw = Real(ep2)*qsw / (p_arr(i,j,k) - qsw);
                    qsw = amrex::max(qsw, Real(qmin));
                    qsatw_arr(i,j,k) = qsw;
                    rhw_arr(i,j,k)   = amrex::max(qv_arr(i,j,k)/qsw, Real(qmin));
                    Real qsi = (t_arr(i,j,k) < ttp)
                        ? Real(psat)*std::exp(std::log(tr)*xai)*std::exp(xbi*(Real(1.0)-tr))
                        : Real(psat)*std::exp(std::log(tr)*xa )*std::exp(xb *(Real(1.0)-tr));
                    qsi = amrex::min(qsi, Real(0.99)*p_arr(i,j,k));
                    qsi = Real(ep2)*qsi / (p_arr(i,j,k) - qsi);
                    qsi = amrex::max(qsi, Real(qmin));
                    qsati_arr(i,j,k) = qsi;
                    rhi_arr(i,j,k)   = amrex::max(qv_arr(i,j,k)/qsi, Real(qmin));
                });
            }

            // G2: zero all process rates each sub-step  [lines 555-594]
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                prevp_arr(i,j,k) = Real(0.0); psdep_arr(i,j,k) = Real(0.0);
                pgdep_arr(i,j,k) = Real(0.0); praut_arr(i,j,k) = Real(0.0);
                psaut_arr(i,j,k) = Real(0.0); pgaut_arr(i,j,k) = Real(0.0);
                pracw_arr(i,j,k) = Real(0.0); praci_arr(i,j,k) = Real(0.0);
                piacr_arr(i,j,k) = Real(0.0); psaci_arr(i,j,k) = Real(0.0);
                psacw_arr(i,j,k) = Real(0.0); pracs_arr(i,j,k) = Real(0.0);
                psacr_arr(i,j,k) = Real(0.0); pgacw_arr(i,j,k) = Real(0.0);
                paacw_arr(i,j,k) = Real(0.0); pgaci_arr(i,j,k) = Real(0.0);
                pgacr_arr(i,j,k) = Real(0.0); pgacs_arr(i,j,k) = Real(0.0);
                pigen_arr(i,j,k) = Real(0.0); pidep_arr(i,j,k) = Real(0.0);
                pcond_arr(i,j,k) = Real(0.0); psmlt_arr(i,j,k) = Real(0.0);
                pgmlt_arr(i,j,k) = Real(0.0); pseml_arr(i,j,k) = Real(0.0);
                pgeml_arr(i,j,k) = Real(0.0); psevp_arr(i,j,k) = Real(0.0);
                pgevp_arr(i,j,k) = Real(0.0);
                fall_r_arr(i,j,k) = Real(0.0); fall_s_arr(i,j,k) = Real(0.0);
                fall_g_arr(i,j,k) = Real(0.0); fallc_arr(i,j,k) = Real(0.0);
            });

            // G3: xni ice crystal number concentration  [lines 598-604]
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                const Real tmp = den_arr(i,j,k)*amrex::max(qi_arr(i,j,k), Real(qmin));
                xni_arr(i,j,k) = amrex::min(
                    amrex::max(Real(5.38e7)*std::sqrt(std::sqrt(tmp*tmp*tmp)), Real(1.e3)),
                    Real(1.e6));
            });

            // G4: pack qrs_tmp, first slope_wsm6 [lines 610-618]
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qrs_tmp_r_arr(i,j,k) = qr_arr(i,j,k);
                qrs_tmp_s_arr(i,j,k) = qs_arr(i,j,k);
                qrs_tmp_g_arr(i,j,k) = qg_arr(i,j,k);
                Real dummy_n0sfac;
                wsm6_slope_rain_cell(
                    qrs_tmp_r_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    m_pidn0r, Real(qcrmin), Real(rslopermax), Real(rsloperbmax),
                    Real(rsloper2max), Real(rsloper3max), Real(bvtr), Real(pvtr),
                    rslope_r_arr(i,j,k), rslopeb_r_arr(i,j,k),
                    rslope2_r_arr(i,j,k), rslope3_r_arr(i,j,k),
                    work1_r_arr(i,j,k));
                wsm6_slope_snow_cell(
                    qrs_tmp_s_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    t_arr(i,j,k), m_pidn0s, Real(alpha_wsm6),
                    Real(n0smax), Real(n0s), Real(t0c), Real(qcrmin),
                    Real(rslopesmax), Real(rslopesbmax),
                    Real(rslopes2max), Real(rslopes3max),
                    Real(bvts), Real(pvts),
                    rslope_s_arr(i,j,k), rslopeb_s_arr(i,j,k),
                    rslope2_s_arr(i,j,k), rslope3_s_arr(i,j,k),
                    work1_s_arr(i,j,k), dummy_n0sfac);
                wsm6_slope_graup_cell(
                    qrs_tmp_g_arr(i,j,k), den_arr(i,j,k), denfac_arr(i,j,k),
                    m_pidn0g, Real(qcrmin),
                    Real(rslopegmax), Real(rslopegbmax),
                    Real(rslopeg2max), Real(rslopeg3max),
                    Real(bvtg), Real(pvtg),
                    rslope_g_arr(i,j,k), rslopeb_g_arr(i,j,k),
                    rslope2_g_arr(i,j,k), rslope3_g_arr(i,j,k),
                    work1_g_arr(i,j,k));
                n0sfac_arr(i,j,k) = dummy_n0sfac;
            });

            amrex::ignore_unused(dtcld,
                work2_arr, workdiffw_arr, workdiffi_arr,
                workr_arr, worka_arr, work1c_arr,
                denqrs1_arr, denqrs2_arr, denqrs3_arr, denqci_arr,
                qsum_arr, delqrs1_arr, delqrs2_arr, delqrs3_arr, delqi_arr,
                tstepsnow_arr, tstepgraup_arr,
                ep1, xls, xlf0, denr);
        }
#endif
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rain_arr(i,j,klo) = rainacc_arr(i,j,0);
            snow_arr(i,j,klo) = snowacc_arr(i,j,0);
            graup_arr(i,j,klo) = graupacc_arr(i,j,0);
        });
    }
}
