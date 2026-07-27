#include "ERF_WDM6.H"
#include <AMReX_Reduce.H>
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#ifdef ERF_USE_WDM6_FORT
#include "ERF_WDM6_Fortran_Interface.H"
#endif

using namespace amrex;

// ---------------------------------------------------------------
// WDM6 device-callable helper functions
// Statement functions from WRF WDM6
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_cpmcal (Real x, Real qmin_arg, Real cpd_arg, Real cpv_arg) {
    return cpd_arg*(Real(1.0)-amrex::max(x,qmin_arg))
          +amrex::max(x,qmin_arg)*cpv_arg;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_xlcal (Real x, Real xlv0_arg, Real xlv1_arg, Real t0c_arg) {
    return xlv0_arg - xlv1_arg*(x - t0c_arg);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_diffus (Real x, Real y) {
    return Real(8.794e-5)*std::exp(std::log(x)*Real(1.81))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_viscos (Real x, Real y) {
    return Real(1.496e-6)*(x*std::sqrt(x))/(x+Real(120.0))/y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_xka (Real x, Real y) {
    return Real(1.414e3)*wdm6_viscos(x,y)*y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_diffac (Real a, Real b, Real c, Real d, Real e,
                   Real rv_arg) {
    return d*a*a/(wdm6_xka(c,d)*rv_arg*c*c)
          +Real(1.0)/(e*wdm6_diffus(c,b));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_venfac (Real a, Real b, Real c, Real den0_arg) {
    return std::exp(std::log(wdm6_viscos(b,c)/wdm6_diffus(b,a))
                   *Real(0.3333333))
          /std::sqrt(wdm6_viscos(b,c))
          *std::sqrt(std::sqrt(den0_arg/c));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_conden (Real a, Real b, Real c, Real d, Real e,
                   Real qmin_arg, Real rv_arg) {
    return (amrex::max(b,qmin_arg)-c)
          /(Real(1.0)+d*d/(rv_arg*e)*c/(a*a));
}

// ---------------------------------------------------------------
// Main Advance routine
// ---------------------------------------------------------------

void WDM6::Advance(const Real& dt_advance,
                   const SolverChoice& solverChoice)
{
    // ---------------------------------------------------------------
    // Dual-mode implementation following WSM6 pattern:
    // - With ERF_USE_WDM6_FORT: Call Fortran bridge (CPU-only)
    // - Without: Use C++ GPU kernels (not yet implemented)
    // ---------------------------------------------------------------

#ifdef ERF_USE_WDM6_FORT
    // Fortran bridge mode - initialize once
    static bool wdm6_inited = false;
    if (!wdm6_inited) {
        constexpr double den0 = 1.28;                 // Standard dry-air density (kg/m^3)
        constexpr double denr = static_cast<double>(rhoh2o);
        constexpr double dens = static_cast<double>(rhos);
        constexpr double cl = static_cast<double>(Cp_l);
        constexpr double cpv = static_cast<double>(Cp_v);
        const double ccn0 = static_cast<double>(m_ccn0);
        constexpr int hail_opt = 0;                   // Graupel mode
        mp_wdm6_init_c(den0, denr, dens, cl, cpv, ccn0, hail_opt);
        wdm6_inited = true;
        amrex::Print() << "WDM6 Fortran bridge initialized\n";
    }
#endif

    // Physical constants
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
    const double ccn0 = static_cast<double>(m_ccn0);

    for (MFIter mfi(*mic_fab_vars[MicVar_WDM6::qv], TileNoZ()); mfi.isValid(); ++mfi) {
        const Box box = mfi.tilebox();
        const Box fab_box = mfi.fabbox();

        // Get array pointers for all WDM6 variables
        auto const& t_arr = mic_fab_vars[MicVar_WDM6::tabs]->array(mfi);
        auto const& qv_arr = mic_fab_vars[MicVar_WDM6::qv]->array(mfi);
        auto const& qc_arr = mic_fab_vars[MicVar_WDM6::qc]->array(mfi);
        auto const& qi_arr = mic_fab_vars[MicVar_WDM6::qi]->array(mfi);
        auto const& qr_arr = mic_fab_vars[MicVar_WDM6::qr]->array(mfi);
        auto const& qs_arr = mic_fab_vars[MicVar_WDM6::qs]->array(mfi);
        auto const& qg_arr = mic_fab_vars[MicVar_WDM6::qg]->array(mfi);
        auto const& nn_arr = mic_fab_vars[MicVar_WDM6::nn]->array(mfi);  // Aerosol number
        auto const& nc_arr = mic_fab_vars[MicVar_WDM6::nc]->array(mfi);  // Cloud droplet number
        auto const& nr_arr = mic_fab_vars[MicVar_WDM6::nr]->array(mfi);  // Rain drop number
        auto const& den_arr = mic_fab_vars[MicVar_WDM6::rho]->array(mfi);
        auto const& p_arr = mic_fab_vars[MicVar_WDM6::pres]->array(mfi);
        auto rain_arr = mic_fab_vars[MicVar_WDM6::rain_accum]->array(mfi);
        auto snow_arr = mic_fab_vars[MicVar_WDM6::snow_accum]->array(mfi);
        auto graup_arr = mic_fab_vars[MicVar_WDM6::graup_accum]->array(mfi);

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

#ifdef ERF_USE_WDM6_FORT
        // Fortran bridge path

        // Create delz array (cell thickness)
        const Real dz_val = m_geom.CellSize(2);
        FArrayBox delz_fab(fab_box, 1, The_Pinned_Arena());
        auto const& delz_arr = delz_fab.array();
        ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            delz_arr(i,j,k) = dz_val;
        });

        // Create landmask array (xland: 0=water, 1=land)
        // TODO: Get from ERF's lmask_lev when available
        // For now, default to land (continental CCN)
        FArrayBox xland_fab(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), 1, The_Pinned_Arena());
        auto const& xland_arr = xland_fab.array();
        ParallelFor(Box(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            xland_arr(i,j,k) = Real(1.0);  // Default to land
        });

        // Create 2D accumulation arrays
        Box box2d(IntVect(imlo,jmlo,0), IntVect(imhi,jmhi,0));
        FArrayBox rainacc_fab(box2d, 1, The_Pinned_Arena());
        FArrayBox rainncv_fab(box2d, 1, The_Pinned_Arena());
        FArrayBox sr_fab(box2d, 1, The_Pinned_Arena());
        FArrayBox snowacc_fab(box2d, 1, The_Pinned_Arena());
        FArrayBox snowncv_fab(box2d, 1, The_Pinned_Arena());
        FArrayBox graupacc_fab(box2d, 1, The_Pinned_Arena());
        FArrayBox graupelncv_fab(box2d, 1, The_Pinned_Arena());

        auto const& rainacc_arr = rainacc_fab.array();
        auto const& rainncv_arr = rainncv_fab.array();
        auto const& sr_arr = sr_fab.array();
        auto const& snowacc_arr = snowacc_fab.array();
        auto const& snowncv_arr = snowncv_fab.array();
        auto const& graupacc_arr = graupacc_fab.array();
        auto const& graupelncv_arr = graupelncv_fab.array();

        // Initialize 2D arrays to zero
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainacc_arr(i,j,k) = Real(0.0);
            rainncv_arr(i,j,k) = Real(0.0);
            sr_arr(i,j,k) = Real(0.0);
            snowacc_arr(i,j,k) = Real(0.0);
            snowncv_arr(i,j,k) = Real(0.0);
            graupacc_arr(i,j,k) = Real(0.0);
            graupelncv_arr(i,j,k) = Real(0.0);
        });

        // Call Fortran WDM6
        mp_wdm6_run_c(
            t_arr.dataPtr(),
            qv_arr.dataPtr(), qc_arr.dataPtr(), qi_arr.dataPtr(),
            qr_arr.dataPtr(), qs_arr.dataPtr(), qg_arr.dataPtr(),
            nn_arr.dataPtr(), nc_arr.dataPtr(), nr_arr.dataPtr(),  // WDM6 number concentrations
            den_arr.dataPtr(), p_arr.dataPtr(), delz_arr.dataPtr(),
            static_cast<double>(dt_advance), g, cpd, cpv, rd, rv, t0c, ep1, ep2, qmin,
            xls, xlv0, xlf0, den0, denr, cliq, cice, psat,
            ccn0, xland_arr.dataPtr(),  // WDM6 parameters
            rainacc_arr.dataPtr(), rainncv_arr.dataPtr(), sr_arr.dataPtr(),
            snowacc_arr.dataPtr(), snowncv_arr.dataPtr(),
            graupacc_arr.dataPtr(), graupelncv_arr.dataPtr(),
            imlo, imhi, jmlo, jmhi, kmlo, kmhi,
            ilo, ihi, jlo, jhi, klo, khi,
            0, ilo, jlo);  // microphysics_debug=0, diag_i=ilo, diag_j=jlo

        // Accumulate precipitation
        ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rain_arr(i,j,k) += rainacc_arr(i,j,k);
            snow_arr(i,j,k) += snowacc_arr(i,j,k);
            graup_arr(i,j,k) += graupacc_arr(i,j,k);
        });

#else
        // C++ GPU kernel path (not yet implemented)

        // TODO: Port full WDM6 physics to C++ GPU kernels
        // - CCN activation (activ_conc)
        // - Cloud droplet nucleation
        // - Autoconversion (mass and number)
        // - Accretion processes
        // - Sedimentation (coupled mass/number for rain)
        // - Phase changes
        // - Process rate calculations

        amrex::Print() << "WDM6::Advance() - C++ GPU kernels not yet implemented\n";
        amrex::Print() << "  CCN0 = " << m_ccn0 << " m^-3\n";
        amrex::Print() << "  dt = " << dt_advance << " s\n";
        amrex::Print() << "  Build with -DERF_ENABLE_WDM6_FORT=ON to use Fortran bridge\n";

        // For now, just enforce minimum number concentrations
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(1e1));   // ncmin
            nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(1e-2));  // nrmin
            nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k), Real(0.0));
        });
#endif

    } // MFIter loop
}
