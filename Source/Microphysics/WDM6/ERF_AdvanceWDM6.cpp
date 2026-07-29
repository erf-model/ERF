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
// WDM6 double-moment slope parameter functions
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdar (Real qr, Real den, Real nr, Real pidnr_arg) {
    // Rain slope parameter using number concentration nr
    // qr = rain mixing ratio (kg/kg)
    // den = air density (kg/m^3)
    // nr = rain number concentration (#/kg)
    // Returns 1/lambda (m)
    if (qr <= Real(1.e-9) || nr <= Real(1.e-2)) {
        return Real(1.0) / Real(5.0e4);  // lamdarmax
    }
    return std::pow((pidnr_arg * nr * den) / (den * qr), Real(1.0)/Real(3.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_lamdac (Real qc, Real den, Real nc, Real pidnc_arg) {
    // Cloud droplet slope parameter using number concentration nc
    // qc = cloud mixing ratio (kg/kg)
    // den = air density (kg/m^3)
    // nc = cloud droplet number concentration (#/kg)
    // Returns 1/lambda (m)
    if (qc <= Real(1.e-9) || nc <= Real(1.e1)) {
        return Real(1.0) / Real(5.0e5);  // lamdacmax
    }
    return std::pow((pidnc_arg * nc * den) / (den * qc), Real(1.0)/Real(3.0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real wdm6_mean_droplet_diameter (Real qc, Real nc, Real den, Real pidnc_arg) {
    // Volume-weighted mean diameter of cloud droplets
    // Returns diameter in meters
    if (nc < Real(1.e1) || qc < Real(1.e-9)) return Real(0.0);

    // Mean mass = rho * qc / nc
    // Mean volume = mean mass / rho_water
    // Mean diameter = (6 * volume / pi)^(1/3)
    Real mean_mass = den * qc / (nc * den);  // kg per droplet
    Real mean_volume = mean_mass / Real(1000.0);  // m^3 per droplet
    Real diameter = std::pow(Real(6.0) * mean_volume / Real(3.14159265359), Real(1.0)/Real(3.0));

    return diameter;
}

// ---------------------------------------------------------------
// WDM6 CCN activation (simplified Abdul-Razzak & Ghan 2000)
// ---------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wdm6_ccn_activation (
    Real& nc,           // cloud droplet number (#/kg) - in/out
    Real& nn,           // aerosol number (#/kg) - in/out
    const Real qv,      // water vapor mixing ratio (kg/kg)
    const Real qc,      // cloud water mixing ratio (kg/kg)
    const Real qvs,     // saturation mixing ratio (kg/kg)
    const Real temp,    // temperature (K)
    const Real w,       // vertical velocity (m/s)
    const Real dt,      // timestep (s)
    const Real ccn0,    // background CCN (#/m^3)
    const Real den      // air density (kg/m^3)
) {
    // Only activate if supersaturated and cloud exists
    if (qv <= qvs || qc < Real(1.e-8)) return;

    // Supersaturation
    Real supersaturation = qv / qvs - Real(1.0);
    supersaturation = amrex::min(supersaturation, Real(0.0048));  // satmax - 1

    // Critical supersaturation based on updraft strength
    // From WRF: stronger updrafts -> more activation
    Real s_crit = Real(0.6) * std::pow(amrex::max(w, Real(0.01)), Real(0.5));

    if (supersaturation > s_crit && nn > Real(0.0)) {
        // Activation rate: convert aerosols to droplets
        Real nc_activate = amrex::min(nn, Real(0.1) * (supersaturation - s_crit) / s_crit * nn);

        // Update concentrations
        nc += nc_activate;
        nn -= nc_activate;

        // Enforce physical bounds
        nc = amrex::max(nc, Real(1.e1));
        nn = amrex::max(nn, Real(0.0));
    }
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
        // ===================================================================
        // WDM6 C++ GPU kernel path (adapted from WSM6 with double-moment)
        // ===================================================================

        // Working FABs (similar to WSM6 but with WDM6-specific additions)
#if defined(AMREX_USE_GPU)
        Arena* Arena_Used = The_Async_Arena();
#else
        Arena* Arena_Used = The_Async_Arena();
#endif

        // 3D working arrays
        FArrayBox denfac_fab(fab_box,1, Arena_Used);
        FArrayBox xni_fab(fab_box,1, Arena_Used);
        FArrayBox cpm_fab(fab_box,1, Arena_Used);
        FArrayBox xl_fab(fab_box,1, Arena_Used);
        FArrayBox qsatw_fab(fab_box,1, Arena_Used);
        FArrayBox qsati_fab(fab_box,1, Arena_Used);
        FArrayBox rhw_fab(fab_box,1, Arena_Used);
        FArrayBox rhi_fab(fab_box,1, Arena_Used);

        // Process rate arrays
        FArrayBox praut_fab(fab_box,1, Arena_Used);
        FArrayBox pracw_fab(fab_box,1, Arena_Used);
        FArrayBox prevp_fab(fab_box,1, Arena_Used);
        FArrayBox pcond_fab(fab_box,1, Arena_Used);

        // Number concentration process rates (WDM6-specific)
        FArrayBox ncauto_fab(fab_box,1, Arena_Used);  // nc lost to autoconversion
        FArrayBox ncaccr_fab(fab_box,1, Arena_Used);  // nc lost to accretion by rain
        FArrayBox nrauto_fab(fab_box,1, Arena_Used);  // nr gained from autoconversion
        FArrayBox nraccr_fab(fab_box,1, Arena_Used);  // nr gained from accretion
        FArrayBox nrevp_fab(fab_box,1, Arena_Used);   // nr lost to evaporation

        auto const& denfac_arr = denfac_fab.array();
        auto const& xni_arr = xni_fab.array();
        auto const& cpm_arr = cpm_fab.array();
        auto const& xl_arr = xl_fab.array();
        auto const& qsatw_arr = qsatw_fab.array();
        auto const& qsati_arr = qsati_fab.array();
        auto const& rhw_arr = rhw_fab.array();
        auto const& rhi_arr = rhi_fab.array();
        auto const& praut_arr = praut_fab.array();
        auto const& pracw_arr = pracw_fab.array();
        auto const& prevp_arr = prevp_fab.array();
        auto const& pcond_arr = pcond_fab.array();
        auto const& ncauto_arr = ncauto_fab.array();
        auto const& ncaccr_arr = ncaccr_fab.array();
        auto const& nrauto_arr = nrauto_fab.array();
        auto const& nraccr_arr = nraccr_fab.array();
        auto const& nrevp_arr = nrevp_fab.array();

        // Clamp negative values and enforce minimums
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k), Real(0.0));
            qr_arr(i,j,k) = amrex::max(qr_arr(i,j,k), Real(0.0));
            qi_arr(i,j,k) = amrex::max(qi_arr(i,j,k), Real(0.0));
            qs_arr(i,j,k) = amrex::max(qs_arr(i,j,k), Real(0.0));
            qg_arr(i,j,k) = amrex::max(qg_arr(i,j,k), Real(0.0));

            // WDM6: Enforce minimum number concentrations
            nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(1.e1));   // ncmin
            nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(1.e-2));  // nrmin
            nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k), Real(0.0));
        });

        // Compute cpm and xl once from initial state
        const Real xlv1_loc = m_xlv1;
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            cpm_arr(i,j,k) = wdm6_cpmcal(qv_arr(i,j,k), Real(qmin), Real(cpd), Real(cpv));
            xl_arr(i,j,k) = wdm6_xlcal(t_arr(i,j,k), Real(xlv0), xlv1_loc, Real(t0c));
        });

        // Minor timestep loop (match WSM6 structure)
        const int wdm6_loops = std::max(
            static_cast<int>(std::round(dt_advance / Real(120.0))), 1);  // dtcldcr = 120s
        const Real dtcld = dt_advance / static_cast<Real>(wdm6_loops);

        // Extract WDM6 coefficients
        const Real qc0_loc = m_qc0;
        const Real qc1_loc = m_qc1;
        const Real qck1_loc = m_qck1;
        const Real pidnc_loc = m_pidnc;
        const Real pidnr_loc = m_pidnr;

        for (int loop = 0; loop < wdm6_loops; ++loop) {
            // ============================================================
            // Step 1: Density factor
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                denfac_arr(i,j,k) = std::sqrt(Real(den0) / den_arr(i,j,k));
            });

            // ============================================================
            // Step 2: Saturation calculations
            // ============================================================
            {
                const Real ttp = Real(t0c) + Real(0.01);
                const Real dldt = Real(cpv) - Real(cliq);
                const Real xa = -dldt / Real(rv);
                const Real xb = xa + Real(xlv0) / (Real(rv) * ttp);
                const Real dldti = Real(cpv) - Real(cice);
                const Real xai = -dldti / Real(rv);
                const Real xbi = xai + Real(xls) / (Real(rv) * ttp);

                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real tr = ttp / t_arr(i,j,k);

                    // Saturation over water
                    Real qsw = Real(psat) * std::exp(std::log(tr) * xa) * std::exp(xb * (Real(1.0) - tr));
                    qsw = amrex::min(qsw, Real(0.99) * p_arr(i,j,k));
                    qsw = Real(ep2) * qsw / (p_arr(i,j,k) - qsw);
                    qsw = amrex::max(qsw, Real(qmin));
                    qsatw_arr(i,j,k) = qsw;
                    rhw_arr(i,j,k) = amrex::max(qv_arr(i,j,k) / qsw, Real(qmin));

                    // Saturation over ice
                    Real qsi = (t_arr(i,j,k) < ttp)
                        ? Real(psat) * std::exp(std::log(tr) * xai) * std::exp(xbi * (Real(1.0) - tr))
                        : Real(psat) * std::exp(std::log(tr) * xa) * std::exp(xb * (Real(1.0) - tr));
                    qsi = amrex::min(qsi, Real(0.99) * p_arr(i,j,k));
                    qsi = Real(ep2) * qsi / (p_arr(i,j,k) - qsi);
                    qsi = amrex::max(qsi, Real(qmin));
                    qsati_arr(i,j,k) = qsi;
                    rhi_arr(i,j,k) = amrex::max(qv_arr(i,j,k) / qsi, Real(qmin));
                });
            }

            // ============================================================
            // Step 3: Zero process rates
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                praut_arr(i,j,k) = Real(0.0);
                pracw_arr(i,j,k) = Real(0.0);
                prevp_arr(i,j,k) = Real(0.0);
                pcond_arr(i,j,k) = Real(0.0);
                ncauto_arr(i,j,k) = Real(0.0);
                ncaccr_arr(i,j,k) = Real(0.0);
                nrauto_arr(i,j,k) = Real(0.0);
                nraccr_arr(i,j,k) = Real(0.0);
                nrevp_arr(i,j,k) = Real(0.0);
            });

            // ============================================================
            // Step 4: WDM6 CCN Activation
            // ============================================================
            // TODO: Need vertical velocity from ERF state
            // For now, use a placeholder (could extract from momentum)
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Placeholder: assume small updraft of 0.1 m/s
                // In production, extract from w-velocity component
                Real w_velocity = Real(0.1);

                wdm6_ccn_activation(
                    nc_arr(i,j,k), nn_arr(i,j,k),
                    qv_arr(i,j,k), qc_arr(i,j,k),
                    qsatw_arr(i,j,k), t_arr(i,j,k),
                    w_velocity, dtcld, ccn0, den_arr(i,j,k)
                );
            });

            // ============================================================
            // Step 5: Condensation/Evaporation (adapted from WSM6)
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (t_arr(i,j,k) > Real(t0c)) {
                    // Warm rain: condensation/evaporation
                    Real qcond = wdm6_conden(
                        t_arr(i,j,k), qv_arr(i,j,k), qsatw_arr(i,j,k),
                        xl_arr(i,j,k), cpm_arr(i,j,k),
                        Real(qmin), Real(rv)
                    );

                    pcond_arr(i,j,k) = qcond / dtcld;
                    qv_arr(i,j,k) -= qcond;
                    qc_arr(i,j,k) += qcond;
                    t_arr(i,j,k) += qcond * xl_arr(i,j,k) / cpm_arr(i,j,k);
                }
            });

            // ============================================================
            // Step 6: WDM6 Double-Moment Autoconversion
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qc_arr(i,j,k) > qc1_loc && nc_arr(i,j,k) > Real(1.e1)) {
                    // Calculate mean droplet diameter
                    Real mean_dia = wdm6_mean_droplet_diameter(
                        qc_arr(i,j,k), nc_arr(i,j,k), den_arr(i,j,k), pidnc_loc
                    );

                    // Autoconversion only if droplets exceed threshold size
                    if (mean_dia > Real(15.e-6)) {  // di15
                        // Mass autoconversion rate (Berry & Reinhardt 1974)
                        Real auto_qc = qck1_loc * qc_arr(i,j,k) * qc_arr(i,j,k) * nc_arr(i,j,k);
                        auto_qc = amrex::min(auto_qc * dtcld, qc_arr(i,j,k));

                        // Number autoconversion: Long's collection kernel
                        // Simplified: convert proportional to mass, with efficiency < 1
                        Real auto_nc = auto_qc * nc_arr(i,j,k) / qc_arr(i,j,k) * Real(0.5);
                        auto_nc = amrex::min(auto_nc, nc_arr(i,j,k));

                        // Number of new raindrops created (much less than droplets consumed)
                        Real auto_nr = auto_nc * Real(0.01);  // ~1% of droplets become rain drops

                        // Apply autoconversion
                        praut_arr(i,j,k) = auto_qc / dtcld;
                        qc_arr(i,j,k) -= auto_qc;
                        qr_arr(i,j,k) += auto_qc;

                        ncauto_arr(i,j,k) = auto_nc / dtcld;
                        nc_arr(i,j,k) -= auto_nc;

                        nrauto_arr(i,j,k) = auto_nr / dtcld;
                        nr_arr(i,j,k) += auto_nr;
                    }
                }
            });

            // ============================================================
            // Step 7: Accretion (cloud water collected by rain)
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qr_arr(i,j,k) > Real(1.e-9) && qc_arr(i,j,k) > Real(1.e-9)) {
                    // Accretion rate (continuous collection)
                    Real accr_qc = Real(6.0) * qc_arr(i,j,k) * qr_arr(i,j,k);
                    accr_qc = amrex::min(accr_qc * dtcld, qc_arr(i,j,k));

                    // Cloud droplets collected by rain
                    Real accr_nc = accr_qc * nc_arr(i,j,k) / qc_arr(i,j,k);
                    accr_nc = amrex::min(accr_nc, nc_arr(i,j,k));

                    // Apply accretion
                    pracw_arr(i,j,k) = accr_qc / dtcld;
                    qc_arr(i,j,k) -= accr_qc;
                    qr_arr(i,j,k) += accr_qc;

                    ncaccr_arr(i,j,k) = accr_nc / dtcld;
                    nc_arr(i,j,k) -= accr_nc;
                    // Note: nr unchanged (rain drops grow but don't multiply)
                }
            });

            // ============================================================
            // Step 8: Rain evaporation
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (qr_arr(i,j,k) > Real(1.e-9) && rhw_arr(i,j,k) < Real(1.0)) {
                    // Simplified evaporation (full version in WRF)
                    Real qrevp = Real(0.001) * qr_arr(i,j,k) * (Real(1.0) - rhw_arr(i,j,k));
                    qrevp = amrex::min(qrevp * dtcld, qr_arr(i,j,k));

                    // Rain number lost to evaporation
                    Real nrevp = qrevp * nr_arr(i,j,k) / qr_arr(i,j,k);
                    nrevp = amrex::min(nrevp, nr_arr(i,j,k));

                    prevp_arr(i,j,k) = qrevp / dtcld;
                    qr_arr(i,j,k) -= qrevp;
                    qv_arr(i,j,k) += qrevp;
                    t_arr(i,j,k) -= qrevp * xl_arr(i,j,k) / cpm_arr(i,j,k);

                    nrevp_arr(i,j,k) = nrevp / dtcld;
                    nr_arr(i,j,k) -= nrevp;
                }
            });

            // ============================================================
            // Step 9: Enforce bounds
            // ============================================================
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                qc_arr(i,j,k) = amrex::max(qc_arr(i,j,k), Real(0.0));
                qr_arr(i,j,k) = amrex::max(qr_arr(i,j,k), Real(0.0));
                nc_arr(i,j,k) = amrex::max(nc_arr(i,j,k), Real(1.e1));
                nr_arr(i,j,k) = amrex::max(nr_arr(i,j,k), Real(1.e-2));
                nn_arr(i,j,k) = amrex::max(nn_arr(i,j,k), Real(0.0));
            });

        } // End minor timestep loop

        // ============================================================
        // Step 10: Simplified sedimentation (placeholder)
        // ============================================================
        // TODO: Implement full PLM sedimentation from WRF
        // For now, use simplified top-down fall
        const Real dz_sedi = m_geom.CellSize(2);
        ParallelFor(amrex::makeSlab(box,2,klo), [=] AMREX_GPU_DEVICE (int i, int j, int) {
            Real precip_rain = Real(0.0);

            // Top-down sedimentation
            for (int kk = khi; kk >= klo; --kk) {
                if (qr_arr(i,j,kk) > Real(1.e-9) && nr_arr(i,j,kk) > Real(1.e-2)) {
                    // Terminal velocity using nr
                    Real lamdar = std::pow(
                        (pidnr_loc * nr_arr(i,j,kk) * den_arr(i,j,kk)) / (den_arr(i,j,kk) * qr_arr(i,j,kk)),
                        Real(1.0)/Real(3.0)
                    );
                    Real vt = Real(841.9) * std::pow(Real(1.0) / lamdar, Real(0.8));

                    // Flux out of cell
                    Real flux_qr = qr_arr(i,j,kk) * vt * dt_advance / dz_sedi;
                    Real flux_nr = nr_arr(i,j,kk) * vt * dt_advance / dz_sedi;

                    flux_qr = amrex::min(flux_qr, qr_arr(i,j,kk));
                    flux_nr = amrex::min(flux_nr, nr_arr(i,j,kk));

                    qr_arr(i,j,kk) -= flux_qr;
                    nr_arr(i,j,kk) -= flux_nr;

                    // Add to cell below or accumulate at surface
                    if (kk > klo) {
                        qr_arr(i,j,kk-1) += flux_qr;
                        nr_arr(i,j,kk-1) += flux_nr;
                    } else {
                        precip_rain += flux_qr * den_arr(i,j,kk) * dz_sedi;
                    }
                }
            }

            // Accumulate precipitation (mm)
            rain_arr(i,j,klo) += precip_rain;
        });

        amrex::Print() << "WDM6::Advance() - C++ GPU kernels executed (simplified implementation)\n";
        amrex::Print() << "  CCN0 = " << m_ccn0 << " m^-3\n";
        amrex::Print() << "  dt = " << dt_advance << " s\n";
        amrex::Print() << "  Minor timesteps = " << wdm6_loops << "\n";
        amrex::Print() << "  NOTE: Simplified sedimentation - full PLM not yet implemented\n";
        amrex::Print() << "  NOTE: Ice processes (qi, qs, qg) not yet implemented\n";
        amrex::Print() << "  For production, build with -DERF_ENABLE_WDM6_FORT=ON\n";

#endif

    } // MFIter loop
}
