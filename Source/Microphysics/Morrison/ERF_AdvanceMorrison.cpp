#include <string>
#include <vector>
#include <memory>
#include <complex>
#include <cmath>

#include <AMReX_Math.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF.H"
#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"
#include "ERF_Morrison.H"
#ifdef ERF_USE_MORR_FORT
#include <ERF_Morrison_Fortran_Interface.H>
#endif

using namespace amrex;

namespace MORRInd {
    enum  {
           qrcuten_arr = 0, // cumulus tendencies
           qscuten_arr,
           qicuten_arr,
           lamc,            // Slope parameters and PSD variables
           lami,
           lams,
           lamr,
           lamg,
           cdist1,
           n0i,
           n0s,
           n0r,
           n0g,
           pgam,
           qc3dten, // CLOUD WATER MIXING RATIO TENDENCY
           qi3dten, // CLOUD ICE MIXING RATIO TENDENCY
           qni3dten, // SNOW MIXING RATIO TENDENCY
           qr3dten, // RAIN MIXING RATIO TENDENCY
           ni3dten, // CLOUD ICE NUMBER CONCENTRATION
           ns3dten, // SNOW NUMBER CONCENTRATION
           nr3dten, // RAIN NUMBER CONCENTRATION
           qc3d, // CLOUD WATER MIXING RATIO
           qi3d, // CLOUD ICE MIXING RATIO
           qni3d, // SNOW MIXING RATIO
           qr3d, // RAIN MIXING RATIO
           ni3d, // CLOUD ICE NUMBER CONCENTRATION
           ns3d, // SNOW NUMBER CONCENTRATION
           nr3d, // RAIN NUMBER CONCENTRATION
           t3dten, // TEMPERATURE TENDENCY
           qv3dten, // WATER VAPOR MIXING RATIO TENDENCY
           t3d, // TEMPERATURE
           qv3d, // WATER VAPOR MIXING RATIO
           pres, // ATMOSPHERIC PRESSURE
           dzq, // DIFFERENCE IN HEIGHT ACROSS LEVEL
           w3d, // GRID-SCALE VERTICAL VELOCITY
           nc3d, // CLOUD DROPLET NUMBER CONCENTRATION
           nc3dten, // CLOUD DROPLET NUMBER CONCENTRATION TENDENCY
           qg3dten, // GRAUPEL MIX RATIO TENDENCY
           ng3dten, // GRAUPEL NUMB CONC TENDENCY
           qg3d, // GRAUPEL MIX RATIO
           ng3d, // GRAUPEL NUMBER CONC
           qgsten, // GRAUPEL SED TEND
           qrsten, // RAIN SED TEND
           qisten, // CLOUD ICE SED TEND
           qnisten, // SNOW SED TEND
           qcsten, // CLOUD WAT SED TEND
           qrcu1d, // RAIN FROM CUMULUS PARAMETERIZATION
           qscu1d, // SNOW FROM CUMULUS PARAMETERIZATION
           qicu1d, // ICE FROM CUMULUS PARAMETERIZATION
           precrt, // TOTAL PRECIP PER TIME STEP
           snowrt, // SNOW PER TIME STEP
           snowprt, // TOTAL CLOUD ICE PLUS SNOW PER TIME STEP
           grplprt, // TOTAL GRAUPEL PER TIME STEP
           effc,  // DROPLET EFFECTIVE RADIUS
           effi, // CLOUD ICE EFFECTIVE RADIUS
           effs, // SNOW EFFECTIVE RADIUS
           effr, // RAIN EFFECTIVE RADIUS
           effg, // GRAUPEL EFFECTIVE RADIUS
           rho,
           mu,
           ain,
           arn,
           asn,
           acn,
           agn,
           dumi, // MASSIVE LIST STARTS
           dumr,
           dumfni,
           dumg,
           dumfng,
           uni,
           umi,
           umr,
           fr,
           fi,
           fni,
           fg,
           fng,
           rgvm,
           faloutr,
           falouti,
           faloutni,
           faltndr,
           faltndi,
           faltndni,
           dumqs,
           dumfns,
           ums,
           uns,
           fs,
           fns,
           falouts,
           faloutns,
           faloutg,
           faloutng,
           faltnds,
           faltndns,
           unr,
           faltndg,
           faltndng,
           dumc,
           dumfnc,
           unc,
           umc,
           ung,
           umg,
           fc,
           faloutc,
           faloutnc,
           faltndc,
           faltndnc,
           fnc,
           dumfnr,
           faloutnr,
           faltndnr,
           fnr,
           dlams,
           dlamr,
           dlami,
           dlamc,
           dlamg,
           xxls, // Latent heat of sublimation
           xxlv, // Latent heat of vaporization
           cpm, // Specific heat at constant pressure for moist air
           xlf, // Latent heat of freezing
           NumInds
    };
}

constexpr Real xxx = 0.9189385332046727417803297;
/*
!------------------------------------------------------------------------------

      REAL(C_DOUBLE) FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL(C_DOUBLE) ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL(C_DOUBLE), SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wrf_gamma (amrex::Real x)
{
    // Debug: using printf since it's GPU compatible
//    printf("wrf_gamma: Input value x = %g\n", x);

    // Local variables
    int i, n;
    bool parity = false;
    amrex::Real fact, half, one, res, sum, twelve, two, xbig, xden, xinf, xminin;
    amrex::Real xnum, y, y1, ysq, z, zero;
    amrex::Real c[7];
    amrex::Real p[8];
    amrex::Real q[8];

    // Mathematical constants
    one = 1.0;
    half = 0.5;
    twelve = 12.0;
    two = 2.0;
    zero = 0.0;

    // Machine dependent parameters
    xbig = 35.040;
    xminin = 1.18e-38;
    amrex::Real eps = 1.19e-7;
    xinf = 3.4e38;

    // Numerator and denominator coefficients for rational minimax approximation over (1,2)
    p[0] = -1.71618513886549492533811e+0;
    p[1] = 2.47656508055759199108314e+1;
    p[2] = -3.79804256470945635097577e+2;
    p[3] = 6.29331155312818442661052e+2;
    p[4] = 8.66966202790413211295064e+2;
    p[5] = -3.14512729688483675254357e+4;
    p[6] = -3.61444134186911729807069e+4;
    p[7] = 6.64561438202405440627855e+4;

    q[0] = -3.08402300119738975254353e+1;
    q[1] = 3.15350626979604161529144e+2;
    q[2] = -1.01515636749021914166146e+3;
    q[3] = -3.10777167157231109440444e+3;
    q[4] = 2.25381184209801510330112e+4;
    q[5] = 4.75584627752788110767815e+3;
    q[6] = -1.34659959864969306392456e+5;
    q[7] = -1.15132259675553483497211e+5;

    // Coefficients for minimax approximation over (12, inf)
    c[0] = -1.910444077728e-03;
    c[1] = 8.4171387781295e-04;
    c[2] = -5.952379913043012e-04;
    c[3] = 7.93650793500350248e-04;
    c[4] = -2.777777777777681622553e-03;
    c[5] = 8.333333333333333331554247e-02;
    c[6] = 5.7083835261e-03;

    // Initialize variables
    parity = false;
    fact = one;
    n = 0;
    y = x;

//    printf("wrf_gamma: Initial y = %g\n", y);

    if (y <= zero) {
        // Argument is negative
//        printf("wrf_gamma: Handling negative argument\n");
        y = -x;
        y1 = std::floor(y);
        res = y - y1;
        if (res != zero) {
            if (y1 != std::floor(y1 * half) * two)
                parity = true;
            Real pi=amrex::Math::pi<Real>();
            fact = -pi / std::sin(pi * res);
            y = y + one;
//            printf("wrf_gamma: After reflection formula: y = %g, fact = %g, parity = %d\n",
//                   y, fact, parity);
        }
        else {
//            printf("wrf_gamma: Singularity detected, returning xinf = %g\n", xinf);
            res = xinf;
            return res;
        }
    }

    // Argument is positive
    if (y < eps) {
        // Argument < eps
//        printf("wrf_gamma: Small argument branch (y < eps)\n");
        if (y >= xminin) {
            res = one / y;
//            printf("wrf_gamma: Small argument result: res = %g\n", res);
        }
        else {
//            printf("wrf_gamma: Argument too small, returning xinf = %g\n", xinf);
            res = xinf;
            return res;
        }
    }
    else if (y < twelve) {
        // Medium range argument
//        printf("wrf_gamma: Medium range branch (eps <= y < 12)\n");
        y1 = y;
        if (y < one) {
            // 0.0 < argument < 1.0
//            printf("wrf_gamma: Sub-branch: 0 < y < 1\n");
            z = y;
            y = y + one;
        }
        else {
            // 1.0 < argument < 12.0, reduce argument if necessary
            n = static_cast<int>(y) - 1;
//            printf("wrf_gamma: Sub-branch: 1 <= y < 12, n = %d\n", n);
            y = y - static_cast<amrex::Real>(n);
            z = y - one;
        }

        // Evaluate approximation
//        printf("wrf_gamma: Before approximation: z = %g, y = %g\n", z, y);
        xnum = zero;
        xden = one;
        for (i = 0; i < 8; i++) {
            xnum = (xnum + p[i]) * z;
            xden = xden * z + q[i];
        }
        res = xnum / xden + one;
//        printf("wrf_gamma: After approximation: res = %g\n", res);

        if (y1 < y) {
            // Adjust result for case 0.0 < argument < 1.0
            res = res / y1;
//            printf("wrf_gamma: Adjusted for y < 1: res = %g\n", res);
        }
        else if (y1 > y) {
            // Adjust for 2.0 < argument < 12.0
//            printf("wrf_gamma: Adjusting for y > 2 with %d multiplications\n", n);
            for (i = 0; i < n; i++) {
                res = res * y;
                y = y + one;
//                printf("wrf_gamma: Multiplication %d: res = %g, y = %g\n", i+1, res, y);
            }
        }
    }
    else {
        // Large argument
//        printf("wrf_gamma: Large argument branch (y >= 12)\n");
        if (y <= xbig) {
            ysq = y * y;
            sum = c[6];
            for (i = 0; i < 6; i++) {
                sum = sum / ysq + c[i];
//                printf("wrf_gamma: Sum step %d: sum = %g\n", i+1, sum);
            }
            sum = sum / y - y + xxx;
            sum = sum + (y - half) * std::log(y);
//            printf("wrf_gamma: Before exp: sum = %g\n", sum);
            res = std::exp(sum);
//            printf("wrf_gamma: After exp: res = %g\n", res);
        }
        else {
//            printf("wrf_gamma: Argument too large, returning xinf = %g\n", xinf);
            res = xinf;
            return res;
        }
    }

    // Final adjustments
    if (parity) {
        res = -res;
//        printf("wrf_gamma: Applied parity adjustment: res = %g\n", res);
    }
    if (fact != one) {
        res = fact / res;
//        printf("wrf_gamma: Applied reflection adjustment: res = %g\n", res);
    }

//    printf("wrf_gamma: Final result = %g\n", res);
    return res;
}

// Gamma function using the custom wrf implementation of the gamma function
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real gamma_function(Real x) {
  return wrf_gamma(x);
}
  /**
   * Helper function to calculate saturation vapor pressure for water or ice.
   * This corresponds to the POLYSVP function in the Fortran code (line ~5580).
   *
   * @param[in] T Temperature in Kelvin
   * @param[in] type 0 for liquid water, 1 for ice
   * @return Saturation vapor pressure in Pascals
   */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real
  calc_saturation_vapor_pressure (const amrex::Real T, const int type)
  {
    amrex::Real polysvp = 0.0;
    amrex::Real del_T = T - 273.15;  // Convert to Celsius

    if (type == 1) {  // Ice (lines ~5631-5644)
        if (T >= 195.8) {
            // Flatau et al. formula for ice
            const amrex::Real a0i = 6.11147274;
            const amrex::Real a1i = 0.503160820;
            const amrex::Real a2i = 0.188439774e-1;
            const amrex::Real a3i = 0.420895665e-3;
            const amrex::Real a4i = 0.615021634e-5;
            const amrex::Real a5i = 0.602588177e-7;
            const amrex::Real a6i = 0.385852041e-9;
            const amrex::Real a7i = 0.146898966e-11;
            const amrex::Real a8i = 0.252751365e-14;

            polysvp = a0i + del_T*(a1i + del_T*(a2i + del_T*(a3i + del_T*(a4i + del_T*(a5i + del_T*(a6i + del_T*(a7i + a8i*del_T)))))));
            polysvp *= 100.0;  // Convert from hPa to Pa
        } else {
            // Goff-Gratch formula for ice at cold temperatures
            polysvp = std::pow(10.0, (-9.09718*(273.16/T-1.0) - 3.56654*std::log10(273.16/T) +
                                      0.876793*(1.0-T/273.16) + std::log10(6.1071))) * 100.0;
            polysvp = 0.0;
        } // T
    } else {  // Water (lines ~5648-5665)
      if (T >= 202.0) {
        // Flatau et al. formula for liquid water
        const amrex::Real a0 = 6.11239921;
        const amrex::Real a1 = 0.443987641;
        const amrex::Real a2 = 0.142986287e-1;
        const amrex::Real a3 = 0.264847430e-3;
        const amrex::Real a4 = 0.302950461e-5;
        const amrex::Real a5 = 0.206739458e-7;
        const amrex::Real a6 = 0.640689451e-10;
        const amrex::Real a7 = -0.952447341e-13;
        const amrex::Real a8 = -0.976195544e-15;

        polysvp = a0 + del_T*(a1 + del_T*(a2 + del_T*(a3 + del_T*(a4 + del_T*(a5 + del_T*(a6 + del_T*(a7 + a8*del_T)))))));
        polysvp *= 100.0;  // Convert from hPa to Pa
      } else {
        // Goff-Gratch formula for water at cold temperatures
        polysvp = std::pow(10.0, (-7.90298*(373.16/T-1.0) + 5.02808*std::log10(373.16/T) -
                                  1.3816e-7*(std::pow(10.0, (11.344*(1.0-T/373.16)))-1.0) +
                                  8.1328e-3*(std::pow(10.0, (-3.49149*(373.16/T-1.0)))-1.0) +
                                  std::log10(1013.246))) * 100.0;
      }
    }

    return polysvp;
  }

    // wrapper to do all the updating
    void
    Morrison::Advance (const amrex::Real& dt_advance,
                       const SolverChoice& sc)
    {
        // Expose for GPU
        bool do_cond = m_do_cond;

        // Store timestep
        amrex::Real dt = dt_advance;

        // Loop through the grids
        for (amrex::MFIter mfi(*mic_fab_vars[MicVar_Morr::qcl],TileNoZ()); mfi.isValid(); ++mfi)
        {
          const amrex::Box& box = mfi.tilebox();

          // Get array data from class member variables
          auto const& theta_arr = mic_fab_vars[MicVar_Morr::theta]->array(mfi);
          auto const& qv_arr = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
          auto const& qcl_arr = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
          auto const& qpr_arr = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
          auto const& qci_arr = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
          auto const& qps_arr = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
          auto const& qpg_arr = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);
          auto const& ni_arr = mic_fab_vars[MicVar_Morr::ni]->array(mfi);
          [[maybe_unused]] auto const& nc_arr = mic_fab_vars[MicVar_Morr::nc]->array(mfi);
          auto const& ns_arr = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
          auto const& nr_arr = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
          auto const& ng_arr = mic_fab_vars[MicVar_Morr::ng]->array(mfi);
          [[maybe_unused]] auto const& rho_arr = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
          auto const& pres_arr = mic_fab_vars[MicVar_Morr::pres]->array(mfi);
          [[maybe_unused]] auto const& tabs_arr = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
          auto const& rain_accum_arr = mic_fab_vars[MicVar_Morr::rain_accum]->array(mfi);
          auto const& snow_accum_arr = mic_fab_vars[MicVar_Morr::snow_accum]->array(mfi);
          auto const& graup_accum_arr = mic_fab_vars[MicVar_Morr::graup_accum]->array(mfi);
          auto const& w_arr = mic_fab_vars[MicVar_Morr::omega]->array(mfi);

          // Get radar reflectivity array if radar diagnostics enabled
          //        auto const& refl_arr = m_do_radar_ref ? m_radar->array(mfi) : nullptr;
          //        auto const& refl_arr = m_radar->array(mfi);

          // Extract box dimensions
          const int ilo = box.loVect()[0];
          const int ihi = box.hiVect()[0];
          const int jlo = box.loVect()[1];
          const int jhi = box.hiVect()[1];
          const int klo = box.loVect()[2];
          const int khi = box.hiVect()[2];

          amrex::Box grown_box(box); grown_box.grow(3);
#ifdef ERF_USE_MORR_FORT
          const int ilom = grown_box.loVect()[0];
          const int ihim = grown_box.hiVect()[0];
          const int jlom = grown_box.loVect()[1];
          const int jhim = grown_box.hiVect()[1];
          const int klom = grown_box.loVect()[2];
          const int khim = grown_box.hiVect()[2];
#endif
          // Calculate Exner function (PII) to convert potential temperature to temperature
          // PII = (P/P0)^(R/cp)
          FArrayBox pii_fab(grown_box, 1);
          auto const& pii_arr = pii_fab.array();

          const amrex::Real p0 = 100000.0; // Reference pressure (Pa)

          const amrex::Real rdcp = m_rdOcp; // R/cp ratio

          // Calculate Exner function
          ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // NOTE: the Morrison Fortran version uses Pa not hPa so we didn't divide p by 100
            //       so we don't need to multiply by 100 here
            pii_arr(i,j,k) = std::pow((pres_arr(i,j,k)) / p0, rdcp);
          });

          // Create arrays for height differences (dz)
          FArrayBox dz_fab(grown_box, 1);
          auto const& dz_arr = dz_fab.array();

          // Calculate height differences
          const amrex::Real dz_val = m_geom.CellSize(m_axis);
          ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            dz_arr(i,j,k) = dz_val;
          });
          amrex::Box grown_boxD(grown_box); grown_boxD.makeSlab(2,0);

          // Arrays to store precipitation rates
          FArrayBox    rainncv_fab(grown_boxD, 1);
          FArrayBox         sr_fab(grown_boxD, 1);     // Ratio of snow to total precipitation
          FArrayBox    snowncv_fab(grown_boxD, 1);
          FArrayBox graupelncv_fab(grown_boxD, 1);

          auto const& rainncv_arr = rainncv_fab.array();
          auto const& sr_arr = sr_fab.array();
          auto const& snowncv_arr = snowncv_fab.array();
          auto const& graupelncv_arr = graupelncv_fab.array();

          // Initialize precipitation rate arrays to zero
          ParallelFor(grown_boxD, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainncv_arr(i,j,k) = 0.0;
            sr_arr(i,j,k) = 0.0;
            snowncv_arr(i,j,k) = 0.0;
            graupelncv_arr(i,j,k) = 0.0;
          });

          // Create terrain height array (not actually used by Morrison scheme)
          FArrayBox ht_fab(amrex::Box(amrex::IntVect(ilo, jlo, 0), amrex::IntVect(ihi, jhi, 0)), 1);
          [[maybe_unused]] auto const& ht_arr = ht_fab.array();
          ParallelFor(amrex::Box(amrex::IntVect(ilo, jlo, 0), amrex::IntVect(ihi, jhi, 0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            ht_arr(i,j,k) = 0.0;  // Not used by Morrison scheme
          });

#ifdef ERF_USE_MORR_FORT
          // Create dummy arrays for cumulus tendencies (if needed)
          FArrayBox qrcuten_fab(grown_box, 1);
          FArrayBox qscuten_fab(grown_box, 1);
          FArrayBox qicuten_fab(grown_box, 1);
          auto const& qrcuten_arr = qrcuten_fab.array();
          auto const& qscuten_arr = qscuten_fab.array();
          auto const& qicuten_arr = qicuten_fab.array();

          // Initialize tendencies to zero (no cumulus parameterization in this example)
          ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            qrcuten_arr(i,j,k) = 0.0;
            qscuten_arr(i,j,k) = 0.0;
            qicuten_arr(i,j,k) = 0.0;
          });

          // WRF-Chem related variables (optional)
          bool flag_qndrop = false;  // Flag to indicate droplet number prediction

          // Now create arrays for other optional variables
          FArrayBox rainprod_fab(grown_box, 1);
          FArrayBox evapprod_fab(grown_box, 1);
          FArrayBox qlsink_fab(grown_box, 1);
          FArrayBox precr_fab(grown_box, 1);
          FArrayBox preci_fab(grown_box, 1);
          FArrayBox precs_fab(grown_box, 1);
          FArrayBox precg_fab(grown_box, 1);

          auto const& rainprod_arr = rainprod_fab.array();
          auto const& evapprod_arr = evapprod_fab.array();
          auto const& qlsink_arr = qlsink_fab.array();
          auto const& precr_arr = precr_fab.array();
          auto const& preci_arr = preci_fab.array();
          auto const& precs_arr = precs_fab.array();
          auto const& precg_arr = precg_fab.array();

          // Initialize WRF-Chem arrays to zero
          ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainprod_arr(i,j,k) = 0.0;
            evapprod_arr(i,j,k) = 0.0;
            qlsink_arr(i,j,k) = 0.0;
            precr_arr(i,j,k) = 0.0;
            preci_arr(i,j,k) = 0.0;
            precs_arr(i,j,k) = 0.0;
            precg_arr(i,j,k) = 0.0;
          });
#endif

#ifdef ERF_USE_MORR_FORT
          // Prepare data pointers for Fortran call
          // These would be passed directly to the Fortran interface
          double dummy_reflectivity = 0.0;
          double* dummy_reflectivity_ptr = &dummy_reflectivity;
#endif
          // Example call (pseudo-code - actual interface would depend on your Fortran interop setup)

          // Microphysics options/switches
          int m_iact = 2;    // CCN activation option (1: power-law, 2: lognormal aerosol)
          int m_inum = 1;    // Droplet number option (0: predict, 1: constant)

          int m_iliq = 0;    // Liquid-only option (0: include ice, 1: liquid only)
          int m_inuc = 0;    // Ice nucleation option (0: mid-latitude, 1: arctic)
          [[maybe_unused]] int m_ibase = 2;   // Cloud base activation option
          [[maybe_unused]] int m_isub = 0;    // Sub-grid vertical velocity option
          int m_igraup = 0;  // Graupel option (0: include graupel, 1: no graupel)
          int m_ihail = 0;   // Graupel/hail option (0: graupel, 1: hail)

          if(sc.moisture_type == MoistureType::Morrison_NoIce) {
            m_iliq = 1;    // Liquid-only option (0: include ice, 1: liquid only)
            m_inuc = 0;    // Ice nucleation option (0: mid-latitude, 1: arctic)
            m_ibase = 2;   // Cloud base activation option
            m_isub = 0;    // Sub-grid vertical velocity option
            m_igraup = 1;  // Graupel option (0: include graupel, 1: no graupel)
            m_ihail = 0;   // Graupel/hail option (0: graupel, 1: hail)
          }
          [[maybe_unused]] bool m_do_radar_ref = false;  // Radar reflectivity calculation flag

          // Physical constants
          amrex::Real m_pi;          // Pi constant
          amrex::Real m_R;           // Gas constant for dry air (J/kg/K)
          amrex::Real m_Rd;           // Gas constant for dry air (J/kg/K)
          amrex::Real m_Rv;          // Gas constant for water vapor (J/kg/K)
          [[maybe_unused]] amrex::Real m_cp;          // Specific heat at constant pressure (J/kg/K)
          amrex::Real m_g;           // Gravitational acceleration (m/s^2)
          amrex::Real m_ep_2;        // Molecular weight ratio (Rd/Rv)

          // Reference density and species densities
          amrex::Real m_rhosu;       // Standard air density at 850 mb (kg/m^3)
          amrex::Real m_rhow;        // Density of liquid water (kg/m^3)
          amrex::Real m_rhoi;        // Bulk density of cloud ice (kg/m^3)
          amrex::Real m_rhosn;       // Bulk density of snow (kg/m^3)
          amrex::Real m_rhog;        // Bulk density of graupel/hail (kg/m^3)

          // Fall speed parameters (V=AD^B)
          amrex::Real m_ai, m_bi;    // Cloud ice fall speed parameters
          [[maybe_unused]] amrex::Real m_ac, m_bc;    // Cloud droplet fall speed parameters
          amrex::Real m_as, m_bs;    // Snow fall speed parameters
          amrex::Real m_ar, m_br;    // Rain fall speed parameters
          amrex::Real m_ag, m_bg;    // Graupel/hail fall speed parameters

          // Microphysical parameters
          amrex::Real m_aimm;        // Parameter in Bigg immersion freezing
          amrex::Real m_bimm;        // Parameter in Bigg immersion freezing
          amrex::Real m_ecr;         // Collection efficiency between droplets/rain and snow/rain
          amrex::Real m_dcs;         // Threshold size for cloud ice autoconversion (m)
          amrex::Real m_mi0;         // Initial mass of nucleated ice crystal (kg)
          amrex::Real m_mg0;         // Mass of embryo graupel (kg)
          amrex::Real m_f1s;         // Ventilation parameter for snow
          amrex::Real m_f2s;         // Ventilation parameter for snow
          amrex::Real m_f1r;         // Ventilation parameter for rain
          amrex::Real m_f2r;         // Ventilation parameter for rain
          amrex::Real m_qsmall;      // Smallest allowed hydrometeor mixing ratio
          amrex::Real m_eii;         // Collection efficiency, ice-ice collisions
          amrex::Real m_eci;         // Collection efficiency, ice-droplet collisions
          amrex::Real m_cpw;         // Specific heat of liquid water (J/kg/K)
          amrex::Real m_rin;         // Radius of contact nuclei (m)
          amrex::Real m_mmult;       // Mass of splintered ice particle (kg)

          // Size distribution parameters
          amrex::Real m_ci, m_di;    // Cloud ice size distribution parameters
          amrex::Real m_cs, m_ds;    // Snow size distribution parameters
          amrex::Real m_cg, m_dg;    // Graupel size distribution parameters

          // Lambda limits for size distributions
          amrex::Real m_lammaxi, m_lammini;    // Cloud ice lambda limits
          amrex::Real m_lammaxr, m_lamminr;    // Rain lambda limits
          amrex::Real m_lammaxs, m_lammins;    // Snow lambda limits
          amrex::Real m_lammaxg, m_lamming;    // Graupel lambda limits

          // Constant droplet concentration (if INUM = 1)
          amrex::Real m_ndcnst = 250.0;  // Droplet number concentration (cm^-3)

          // CCN spectra parameters (for IACT = 1)
          [[maybe_unused]] amrex::Real m_k1;          // Exponent in CCN activation formula
          [[maybe_unused]] amrex::Real m_c1;          // Coefficient in CCN activation formula (cm^-3)

          // Aerosol activation parameters (for IACT = 2)
          [[maybe_unused]] amrex::Real m_mw;          // Molecular weight water (kg/mol)
          [[maybe_unused]] amrex::Real m_osm;         // Osmotic coefficient
          [[maybe_unused]] amrex::Real m_vi;          // Number of ions dissociated in solution
          [[maybe_unused]] amrex::Real m_epsm;        // Aerosol soluble fraction
          [[maybe_unused]] amrex::Real m_rhoa;        // Aerosol bulk density (kg/m^3)
          [[maybe_unused]] amrex::Real m_map;         // Molecular weight aerosol (kg/mol)
          [[maybe_unused]] amrex::Real m_ma;          // Molecular weight of air (kg/mol)
          [[maybe_unused]] amrex::Real m_rr;          // Universal gas constant (J/mol/K)
          [[maybe_unused]] amrex::Real m_bact;        // Activation parameter
          [[maybe_unused]] amrex::Real m_rm1;         // Geometric mean radius, mode 1 (m)
          [[maybe_unused]] amrex::Real m_rm2;         // Geometric mean radius, mode 2 (m)
          amrex::Real m_nanew1;      // Total aerosol concentration, mode 1 (m^-3)
          amrex::Real m_nanew2;      // Total aerosol concentration, mode 2 (m^-3)
          [[maybe_unused]] amrex::Real m_sig1;        // Standard deviation of aerosol dist, mode 1
          [[maybe_unused]] amrex::Real m_sig2;        // Standard deviation of aerosol dist, mode 2
          [[maybe_unused]] amrex::Real m_f11;         // Correction factor for activation, mode 1
          [[maybe_unused]] amrex::Real m_f12;         // Correction factor for activation, mode 1
          [[maybe_unused]] amrex::Real m_f21;         // Correction factor for activation, mode 2
          [[maybe_unused]] amrex::Real m_f22;         // Correction factor for activation, mode 2

          // Precomputed constants for efficiency
          amrex::Real m_cons1, m_cons2, m_cons3, m_cons4, m_cons5;
          amrex::Real m_cons6, m_cons7, m_cons8, m_cons9, m_cons10;
          amrex::Real m_cons11, m_cons12, m_cons13, m_cons14, m_cons15;
          amrex::Real m_cons16, m_cons17, m_cons18, m_cons19, m_cons20;
          amrex::Real m_cons21, m_cons22, m_cons23, m_cons24, m_cons25;
          amrex::Real m_cons26, m_cons27, m_cons28, m_cons29; [[maybe_unused]] amrex::Real m_cons30;
          amrex::Real m_cons31, m_cons32, m_cons34, m_cons35; [[maybe_unused]] amrex::Real m_cons33;
          amrex::Real m_cons36, m_cons37, m_cons38, m_cons39, m_cons40;
          amrex::Real m_cons41;
          // Set microphysics control parameters
          m_inum = 1;           // Use constant droplet number concentration
          m_ndcnst = 250.0;     // Droplet number concentration (cm^-3)
          // Mathematical constants
          m_pi = 3.1415926535897932384626434;

          m_R = 287.0;         // Gas constant for dry air (J/kg/K)
          m_Rd = 287.0;         // Gas constant for dry air (J/kg/K)
          m_Rv = 461.6;        // Gas constant for water vapor (J/kg/K)
          m_cp = 7.0*287.0/2.0;        // Specific heat at constant pressure (J/kg/K)
          m_g = 9.81;           // Gravitational acceleration (m/s^2)
          m_ep_2 = m_Rd / m_Rv;     // Molecular weight ratio (Rd/Rv)

          // Reference density
          m_rhosu = 85000.0/(287.15*273.15);  // Standard air density at 850 mb (kg/m^3)

          // Densities for different hydrometeor species
          m_rhow = 997.0;     // Density of liquid water (kg/m^3)
          m_rhoi = 500.0;     // Bulk density of cloud ice (kg/m^3)
          m_rhosn = 100.0;    // Bulk density of snow (kg/m^3)

          // Set density for graupel or hail based on configuration
          if (m_ihail == 0) {
            m_rhog = 400.0; // Bulk density of graupel (kg/m^3)
          } else {
            m_rhog = 900.0; // Bulk density of hail (kg/m^3)
          }

          // Fall speed parameters (V=AD^B) for different hydrometeors
          // Cloud ice
          m_ai = 700.0;
          m_bi = 1.0;

          // Cloud droplets
          m_ac = 3.0E7;
          m_bc = 2.0;

          // Snow
          m_as = 11.72;
          m_bs = 0.41;

          // Rain
          m_ar = 841.99667;
          m_br = 0.8;

          // Graupel/hail (dependent on configuration)
          if (m_ihail == 0) {
            // Graupel parameters
            m_ag = 19.3;
            m_bg = 0.37;
          } else {
            // Hail parameters (Matsun and Huggins 1980)
            m_ag = 114.5;
            m_bg = 0.5;
          }

          // Microphysical parameters
          m_aimm = 0.66;       // Parameter in Bigg immersion freezing
          m_bimm = 100.0;      // Parameter in Bigg immersion freezing
          m_ecr = 1.0;         // Collection efficiency between rain and snow/graupel
          m_dcs = 125.0E-6;    // Threshold size for cloud ice autoconversion (m)
          m_mi0 = 4.0/3.0*m_pi*m_rhoi*std::pow(10.0E-6, 3);  // Initial mass of nucleated ice crystal (kg)
          m_mg0 = 1.6E-10;     // Mass of embryo graupel (kg)

          // Ventilation parameters
          m_f1s = 0.86;        // Ventilation parameter for snow
          m_f2s = 0.28;        // Ventilation parameter for snow
          m_f1r = 0.78;        // Ventilation parameter for rain
          m_f2r = 0.308;       // Ventilation parameter for rain

          // Smallest allowed hydrometeor mixing ratio
          m_qsmall = 1.0E-14;

          // Collection efficiencies
          m_eii = 0.1;         // Ice-ice collision efficiency
          m_eci = 0.7;         // Ice-droplet collision efficiency

          // Specific heat of liquid water (J/kg/K)
          m_cpw = 4187.0;

          // Size distribution parameters
          m_ci = m_rhoi * m_pi / 6.0;
          m_di = 3.0;
          m_cs = m_rhosn * m_pi / 6.0;
          m_ds = 3.0;
          m_cg = m_rhog * m_pi / 6.0;
          m_dg = 3.0;

          // Radius of contact nuclei (m)
          m_rin = 0.1E-6;

          // Mass of splintered ice particle (kg)
          m_mmult = 4.0/3.0*m_pi*m_rhoi*std::pow(5.0E-6, 3);

          // Set lambda limits for size distributions
          // Maximum and minimum values for lambda parameter in size distributions
          m_lammaxi = 1.0/1.0E-6;
          m_lammini = 1.0/(2.0*m_dcs + 100.0E-6);
          m_lammaxr = 1.0/20.0E-6;
          m_lamminr = 1.0/2800.0E-6;
          m_lammaxs = 1.0/10.0E-6;
          m_lammins = 1.0/2000.0E-6;
          m_lammaxg = 1.0/20.0E-6;
          m_lamming = 1.0/2000.0E-6;

          // Set CCN parameters for different environments
          if (m_iact == 1) {
            // Maritime CCN spectrum parameters (modified from Rasmussen et al. 2002)
            // NCCN = C*S^K, where S is supersaturation in %
            m_k1 = 0.4;        // Exponent in CCN activation formula
            m_c1 = 120.0;      // Coefficient in CCN activation formula (cm^-3)
          }

          // Initialize aerosol activation parameters for lognormal distribution
          if (m_iact == 2) {
            // Parameters for ammonium sulfate
            m_mw = 0.018;      // Molecular weight of water (kg/mol)
            m_osm = 1.0;       // Osmotic coefficient
            m_vi = 3.0;        // Number of ions dissociated in solution
            m_epsm = 0.7;      // Aerosol soluble fraction
            m_rhoa = 1777.0;   // Aerosol bulk density (kg/m^3)
            m_map = 0.132;     // Molecular weight of aerosol (kg/mol)
            m_ma = 0.0284;     // Molecular weight of air (kg/mol)
            m_rr = 8.3145;     // Universal gas constant (J/mol/K)
            m_bact = m_vi * m_osm * m_epsm * m_mw * m_rhoa / (m_map * m_rhow);
            //            m_a_w = 2.0 * m_mw * 0.0761 / (m_rhow * m_r_v * 293.15);  // "A" parameter

            // Aerosol size distribution parameters for MPACE (Morrison et al. 2007, JGR)
            // Mode 1
            m_rm1 = 0.052E-6;  // Geometric mean radius, mode 1 (m)
            m_sig1 = 2.04;     // Standard deviation of aerosol size distribution, mode 1
            m_nanew1 = 72.2E6; // Total aerosol concentration, mode 1 (m^-3)
            m_f11 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig1), 2));
            m_f21 = 1.0 + 0.25 * std::log(m_sig1);

            // Mode 2
            m_rm2 = 1.3E-6;    // Geometric mean radius, mode 2 (m)
            m_sig2 = 2.5;      // Standard deviation of aerosol size distribution, mode 2
            m_nanew2 = 1.8E6;  // Total aerosol concentration, mode 2 (m^-3)
            m_f12 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig2), 2));
            m_f22 = 1.0 + 0.25 * std::log(m_sig2);
          }

          // Precompute constants for efficiency
          m_cons1 = gamma_function(1.0 + m_ds) * m_cs;
          m_cons2 = gamma_function(1.0 + m_dg) * m_cg;
          m_cons3 = gamma_function(4.0 + m_bs) / 6.0;
          m_cons4 = gamma_function(4.0 + m_br) / 6.0;
          m_cons5 = gamma_function(1.0 + m_bs);
          m_cons6 = gamma_function(1.0 + m_br);
          m_cons7 = gamma_function(4.0 + m_bg) / 6.0;
          m_cons8 = gamma_function(1.0 + m_bg);
          m_cons9 = gamma_function(5.0/2.0 + m_br/2.0);
          m_cons10 = gamma_function(5.0/2.0 + m_bs/2.0);
          m_cons11 = gamma_function(5.0/2.0 + m_bg/2.0);
          m_cons12 = gamma_function(1.0 + m_di) * m_ci;
          m_cons13 = gamma_function(m_bs + 3.0) * m_pi / 4.0 * m_eci;
          m_cons14 = gamma_function(m_bg + 3.0) * m_pi / 4.0 * m_eci;
          m_cons15 = -1108.0 * m_eii * std::pow(m_pi, (1.0-m_bs)/3.0) *
            std::pow(m_rhosn, (-2.0-m_bs)/3.0) / (4.0*720.0);
          m_cons16 = gamma_function(m_bi + 3.0) * m_pi / 4.0 * m_eci;
          m_cons17 = 4.0 * 2.0 * 3.0 * m_rhosu * m_pi * m_eci * m_eci *
            gamma_function(2.0*m_bs + 2.0) / (8.0*(m_rhog-m_rhosn));
          m_cons18 = m_rhosn * m_rhosn;
          m_cons19 = m_rhow * m_rhow;
          m_cons20 = 20.0 * m_pi * m_pi * m_rhow * m_bimm;
          m_cons21 = 4.0 / (m_dcs * m_rhoi);
          m_cons22 = m_pi * m_rhoi * std::pow(m_dcs, 3) / 6.0;
          m_cons23 = m_pi / 4.0 * m_eii * gamma_function(m_bs + 3.0);
          m_cons24 = m_pi / 4.0 * m_ecr * gamma_function(m_br + 3.0);
          m_cons25 = m_pi * m_pi / 24.0 * m_rhow * m_ecr * gamma_function(m_br + 6.0);
          m_cons26 = m_pi / 6.0 * m_rhow;
          m_cons27 = gamma_function(1.0 + m_bi);
          m_cons28 = gamma_function(4.0 + m_bi) / 6.0;
          m_cons29 = 4.0/3.0 * m_pi * m_rhow * std::pow(25.0E-6, 3);
          m_cons30 = 4.0/3.0 * m_pi * m_rhow;
          m_cons31 = m_pi * m_pi * m_ecr * m_rhosn;
          m_cons32 = m_pi / 2.0 * m_ecr;
          m_cons33 = m_pi * m_pi * m_ecr * m_rhog;
          m_cons34 = 5.0/2.0 + m_br/2.0;
          m_cons35 = 5.0/2.0 + m_bs/2.0;
          m_cons36 = 5.0/2.0 + m_bg/2.0;
          m_cons37 = 4.0 * m_pi * 1.38E-23 / (6.0 * m_pi * m_rin);
          m_cons38 = m_pi * m_pi / 3.0 * m_rhow;
          m_cons39 = m_pi * m_pi / 36.0 * m_rhow * m_bimm;
          m_cons40 = m_pi / 6.0 * m_bimm;
          m_cons41 = m_pi * m_pi * m_ecr * m_rhow;

          // Set CCN parameters for different environments
          if (m_iact == 1) {
            // Maritime CCN spectrum parameters (modified from Rasmussen et al. 2002)
            // NCCN = C*S^K, where S is supersaturation in %
            m_k1 = 0.4;        // Exponent in CCN activation formula
            m_c1 = 120.0;      // Coefficient in CCN activation formula (cm^-3)
          }

          // Initialize aerosol activation parameters for IACT=2
          if (m_iact == 2) {
            // Parameters for ammonium sulfate
            m_mw = 0.018;      // Molecular weight of water (kg/mol)
            m_osm = 1.0;       // Osmotic coefficient
            m_vi = 3.0;        // Number of ions dissociated in solution
            m_epsm = 0.7;      // Aerosol soluble fraction
            m_rhoa = 1777.0;   // Aerosol bulk density (kg/m^3)
            m_map = 0.132;     // Molecular weight of aerosol (kg/mol)
            m_ma = 0.0284;     // Molecular weight of air (kg/mol)
            m_rr = 8.3145;     // Universal gas constant (J/mol/K)
            m_bact = m_vi * m_osm * m_epsm * m_mw * m_rhoa / (m_map * m_rhow);

            // Aerosol size distribution parameters for MPACE (Morrison et al. 2007, JGR)
            // Mode 1
            m_rm1 = 0.052E-6;  // Geometric mean radius, mode 1 (m)
            m_sig1 = 2.04;     // Standard deviation of aerosol size distribution, mode 1
            m_nanew1 = 72.2E6; // Total aerosol concentration, mode 1 (m^-3)
            m_f11 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig1), 2));
            m_f21 = 1.0 + 0.25 * std::log(m_sig1);

            // Mode 2
            m_rm2 = 1.3E-6;    // Geometric mean radius, mode 2 (m)
            m_sig2 = 2.5;      // Standard deviation of aerosol size distribution, mode 2
            m_nanew2 = 1.8E6;  // Total aerosol concentration, mode 2 (m^-3)
            m_f12 = 0.5 * std::exp(2.5 * std::pow(std::log(m_sig2), 2));
            m_f22 = 1.0 + 0.25 * std::log(m_sig2);
          }
          // Set microphysics control parameters
          m_iact = 2;  // Lognormal aerosol activation
          m_inuc = 0;      // Mid-latitude ice nucleation (Cooper)
          if (sc.moisture_type == MoistureType::Morrison_NoIce) {
              m_iliq = 1;           // Include ice processes
              m_igraup = 1;         // Include graupel processes
          } else {
              m_iliq = 0;           // Include ice processes
              m_igraup = 0;         // Include graupel processes
          }
          m_ihail = 0;          // Use graupel (0) instead of hail (1)
          m_isub = 0;           // Sub-grid vertical velocity option
          m_do_radar_ref = false; // Disable radar reflectivity by default
          amrex::Box boxD(box); boxD.makeSlab(2,0);

          ParmParse pp("erf");

          // Allow user to override constant droplet concentration from inputs file
          pp.query("morrison_ndcnst", m_ndcnst);

#ifdef ERF_USE_MORR_FORT
          // If using Fortran version, update the Fortran module variable as well
          set_morrison_ndcnst_c(m_ndcnst);
#endif

          bool run_morr_cpp = true;

          bool use_morr_cpp_answer = true;
          pp.query("use_morr_cpp_answer", use_morr_cpp_answer);
          Print() << "use_morr_cpp_answer" << use_morr_cpp_answer <<std::endl;

          bool run_morr_fort = !use_morr_cpp_answer;

          std::string filename = std::string("output_cpp") + std::to_string(use_morr_cpp_answer) + ".txt";

          if(run_morr_cpp) {

            // One FAB to rule them all
            FArrayBox morr_fab(grown_box, MORRInd::NumInds);
            morr_fab.template setVal<RunOn::Device>(0.0);
            auto const& morr_arr = morr_fab.array();

          ////////////////////////////////////////////////////////////
          // ParallelFor for testing partial C++ implementation
          // NOTE: Currently all Array4 values are copied to locals
          //       This means we're not updating or outputting anything
          ////////////////////////////////////////////////////////////
          ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
            // Tendencies and mixing ratios
            morr_arr(i,j,k,MORRInd::qc3d) = qcl_arr(i,j,k);   // CLOUD WATER MIXING RATIO
            morr_arr(i,j,k,MORRInd::qi3d) = qci_arr(i,j,k);   // CLOUD ICE MIXING RATIO
            morr_arr(i,j,k,MORRInd::qni3d) = qps_arr(i,j,k);  // SNOW MIXING RATIO
            morr_arr(i,j,k,MORRInd::qr3d) = qpr_arr(i,j,k);   // RAIN MIXING RATIO
            morr_arr(i,j,k,MORRInd::ni3d) = ni_arr(i,j,k);    // CLOUD ICE NUMBER CONCENTRATION
            morr_arr(i,j,k,MORRInd::ns3d) = ns_arr(i,j,k);    // SNOW NUMBER CONCENTRATION
            morr_arr(i,j,k,MORRInd::nr3d) = nr_arr(i,j,k);    // RAIN NUMBER CONCENTRATION
            morr_arr(i,j,k,MORRInd::nc3d) = nc_arr(i,j,k);    // RAIN NUMBER CONCENTRATION

            morr_arr(i,j,k,MORRInd::t3d) = theta_arr(i,j,k) * pii_arr(i,j,k);  // TEMPERATURE
            morr_arr(i,j,k,MORRInd::qv3d) = qv_arr(i,j,k);                     // WATER VAPOR MIXING RATIO
            morr_arr(i,j,k,MORRInd::pres) = pres_arr(i,j,k);                   // ATMOSPHERIC PRESSURE
            morr_arr(i,j,k,MORRInd::dzq) = dz_arr(i,j,k);                      // DIFFERENCE IN HEIGHT ACROSS LEVEL
            morr_arr(i,j,k,MORRInd::w3d) = w_arr(i,j,k);                       // GRID-SCALE VERTICAL VELOCITY
            morr_arr(i,j,k,MORRInd::qg3d) = qpg_arr(i,j,k);                    // GRAUPEL MIX RATIO
            morr_arr(i,j,k,MORRInd::ng3d) = ng_arr(i,j,k);                     // GRAUPEL NUMBER CONC
            morr_arr(i,j,k,MORRInd::qrcu1d) = morr_arr(i,j,k,MORRInd::qrcuten_arr);              // RAIN FROM CUMULUS PARAMETERIZATION
            morr_arr(i,j,k,MORRInd::qscu1d) = morr_arr(i,j,k,MORRInd::qscuten_arr);              // SNOW FROM CUMULUS PARAMETERIZATION
            morr_arr(i,j,k,MORRInd::qicu1d) = morr_arr(i,j,k,MORRInd::qicuten_arr);              // ICE FROM CUMULUS PARAMETERIZATION
         });
          ParallelFor( boxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
         {
           int ltrue=0;                      // LTRUE: SWITCH = 0: NO HYDROMETEORS IN COLUMN, = 1: HYDROMETEORS IN COLUMN
           int nstep;                        // NSTEP: Timestep counter
           int iinum=m_inum;                      // iinum: Integer control variable

           for(int k=klo; k<=khi; k++) {
            // Model input parameters
            //amrex::Real dt;                 // DT: MODEL TIME STEP (SEC)
            //amrex::Real morr_arr(i,j,k,MORRInd::lami);               // LAMI: Slope parameter for cloud ice (m^-1)

            // Microphysical processes
            [[maybe_unused]] amrex::Real nsubc;              // NSUBC: Loss of NC during evaporation
            amrex::Real nsubi;              // NSUBI: Loss of NI during sublimation
            amrex::Real nsubs;              // NSUBS: Loss of NS during sublimation
            amrex::Real nsubr;              // NSUBR: Loss of NR during evaporation
            amrex::Real prd;                // PRD: Deposition cloud ice
            amrex::Real pre;                // PRE: Evaporation of rain
            amrex::Real prds;               // PRDS: Deposition snow
            amrex::Real nnuccc;             // NNUCCC: Change N due to contact freezing droplets
            amrex::Real mnuccc;             // MNUCCC: Change Q due to contact freezing droplets
            amrex::Real pra;                // PRA: Accretion droplets by rain
            amrex::Real prc;                // PRC: Autoconversion droplets
            amrex::Real pcc;                // PCC: Condensation/evaporation droplets
            amrex::Real nnuccd;             // NNUCCD: Change N freezing aerosol (primary ice nucleation)
            amrex::Real mnuccd;             // MNUCCD: Change Q freezing aerosol (primary ice nucleation)
            amrex::Real mnuccr;             // MNUCCR: Change Q due to contact freezing rain
            amrex::Real nnuccr;             // NNUCCR: Change N due to contact freezing rain
            amrex::Real npra;               // NPRA: Change N due to droplet accretion by rain
            amrex::Real nragg;              // NRAGG: Self-collection/breakup of rain
            amrex::Real nsagg;              // NSAGG: Self-collection of snow
            amrex::Real nprc;               // NPRC: Change NC autoconversion droplets
            amrex::Real nprc1;              // NPRC1: Change NR autoconversion droplets
            amrex::Real prai;               // PRAI: Change Q accretion cloud ice by snow
            amrex::Real prci;               // PRCI: Change Q autoconversion cloud ice to snow
            amrex::Real psacws;             // PSACWS: Change Q droplet accretion by snow
            amrex::Real npsacws;            // NPSACWS: Change N droplet accretion by snow
            amrex::Real psacwi;             // PSACWI: Change Q droplet accretion by cloud ice
            amrex::Real npsacwi;            // NPSACWI: Change N droplet accretion by cloud ice
            amrex::Real nprci;              // NPRCI: Change N autoconversion cloud ice by snow
            amrex::Real nprai;              // NPRAI: Change N accretion cloud ice
            amrex::Real nmults;             // NMULTS: Ice multiplication due to riming droplets by snow
            amrex::Real nmultr;             // NMULTR: Ice multiplication due to riming rain by snow
            amrex::Real qmults;             // QMULTS: Change Q due to ice multiplication droplets/snow
            amrex::Real qmultr;             // QMULTR: Change Q due to ice multiplication rain/snow
            amrex::Real pracs;              // PRACS: Change Q rain-snow collection
            amrex::Real npracs;             // NPRACS: Change N rain-snow collection
            [[maybe_unused]] amrex::Real pccn;               // PCCN: Change Q droplet activation
            amrex::Real psmlt;              // PSMLT: Change Q melting snow to rain
            amrex::Real evpms;              // EVPMS: Change Q melting snow evaporating
            amrex::Real nsmlts;             // NSMLTS: Change N melting snow
            amrex::Real nsmltr;             // NSMLTR: Change N melting snow to rain
            amrex::Real piacr;              // PIACR: Change QR, ice-rain collection
            amrex::Real niacr;              // NIACR: Change N, ice-rain collection
            amrex::Real praci;              // PRACI: Change QI, ice-rain collection
            amrex::Real piacrs;             // PIACRS: Change QR, ice rain collision, added to snow
            amrex::Real niacrs;             // NIACRS: Change N, ice rain collision, added to snow
            amrex::Real pracis;             // PRACIS: Change QI, ice rain collision, added to snow
            amrex::Real eprd;               // EPRD: Sublimation cloud ice
            amrex::Real eprds;              // EPRDS: Sublimation snow

            // Graupel processes
            amrex::Real pracg;              // PRACG: Change in Q collection rain by graupel
            amrex::Real psacwg;             // PSACWG: Change in Q collection droplets by graupel
            amrex::Real pgsacw;             // PGSACW: Conversion Q to graupel due to collection droplets by snow
            amrex::Real pgracs;             // PGRACS: Conversion Q to graupel due to collection rain by snow
            amrex::Real prdg;               // PRDG: Deposition of graupel
            amrex::Real eprdg;              // EPRDG: Sublimation of graupel
            amrex::Real evpmg;              // EVPMG: Change Q melting of graupel and evaporation
            amrex::Real pgmlt;              // PGMLT: Change Q melting of graupel
            amrex::Real npracg;             // NPRACG: Change N collection rain by graupel
            amrex::Real npsacwg;            // NPSACWG: Change N collection droplets by graupel
            amrex::Real nscng;              // NSCNG: Change N conversion to graupel due to collection droplets by snow
            amrex::Real ngracs;             // NGRACS: Change N conversion to graupel due to collection rain by snow
            amrex::Real ngmltg;             // NGMLTG: Change N melting graupel
            amrex::Real ngmltr;             // NGMLTR: Change N melting graupel to rain
            amrex::Real nsubg;              // NSUBG: Change N sublimation/deposition of graupel
            amrex::Real psacr;              // PSACR: Conversion due to collection of snow by rain
            amrex::Real nmultg;             // NMULTG: Ice multiplication due to accretion droplets by graupel
            amrex::Real nmultrg;            // NMULTRG: Ice multiplication due to accretion rain by graupel
            amrex::Real qmultg;             // QMULTG: Change Q due to ice multiplication droplets/graupel
            amrex::Real qmultrg;            // QMULTRG: Change Q due to ice multiplication rain/graupel

            // Time-varying atmospheric parameters
            amrex::Real kap;                // KAP: Thermal conductivity of air
            amrex::Real evs;                // EVS: Saturation vapor pressure
            amrex::Real eis;                // EIS: Ice saturation vapor pressure
            amrex::Real qvs;                // QVS: Saturation mixing ratio
            amrex::Real qvi;                // QVI: Ice saturation mixing ratio
            amrex::Real qvqvs;              // QVQVS: Saturation ratio
            amrex::Real qvqvsi;             // QVQVSI: Ice saturation ratio
            amrex::Real dv;                 // DV: Diffusivity of water vapor in air
            amrex::Real sc_schmidt;         // SC: Schmidt number
            amrex::Real ab;                 // AB: Correction to condensation rate due to latent heating
            amrex::Real abi;                // ABI: Correction to deposition rate due to latent heating

            // Dummy variables
            amrex::Real dum;                // DUM: General dummy variable
            amrex::Real dum1;               // DUM1: General dummy variable
            [[maybe_unused]] amrex::Real dum2;               // DUM2: General dummy variable
            amrex::Real dumt;               // DUMT: Dummy variable for temperature
            amrex::Real dumqv;              // DUMQV: Dummy variable for water vapor
            amrex::Real dumqss;             // DUMQSS: Dummy saturation mixing ratio
            [[maybe_unused]] amrex::Real dumqsi;             // DUMQSI: Dummy ice saturation mixing ratio
            amrex::Real dums;               // DUMS: General dummy variable

            // Prognostic supersaturation
            amrex::Real dqsdt;              // DQSDT: Change of saturation mixing ratio with temperature
            amrex::Real dqsidt;             // DQSIDT: Change in ice saturation mixing ratio with temperature

            amrex::Real epsi;               // EPSI: 1/phase relaxation time (see M2005), ice
            amrex::Real epss;               // EPSS: 1/phase relaxation time (see M2005), snow
            amrex::Real epsr;               // EPSR: 1/phase relaxation time (see M2005), rain
            amrex::Real epsg;               // EPSG: 1/phase relaxation time (see M2005), graupel
            amrex::Real kc2;                // KC2: Total ice nucleation rate
            amrex::Real di0;                // DC0: Characteristic diameter for ice
            [[maybe_unused]] amrex::Real dc0;                // DC0: Characteristic diameter for cloud droplets
            amrex::Real ds0;                // DS0: Characteristic diameter for snow
            amrex::Real dg0;                // DG0: Characteristic diameter for graupel
            amrex::Real dumqc;              // DUMQC: Dummy variable for cloud water mixing ratio
            [[maybe_unused]] amrex::Real dumqr;              // DUMQR: Dummy variable for rain mixing ratio
            amrex::Real ratio;              // RATIO: General ratio variable
            amrex::Real sum_dep;            // SUM_DEP: Sum of deposition/sublimation
            amrex::Real fudgef;             // FUDGEF: Adjustment factor
            // For WRF-CHEM
            [[maybe_unused]] amrex::Real c2prec;             // C2PREC: Cloud to precipitation conversion
            [[maybe_unused]] amrex::Real csed;               // CSED: Cloud sedimentation
            [[maybe_unused]] amrex::Real ised;               // ISED: Ice sedimentation
            [[maybe_unused]] amrex::Real ssed;               // SSED: Snow sedimentation
            [[maybe_unused]] amrex::Real gsed;               // GSED: Graupel sedimentation
            [[maybe_unused]] amrex::Real rsed;               // RSED: Rain sedimentation
            [[maybe_unused]] amrex::Real tqimelt;            // tqimelt: Melting of cloud ice (tendency)

            // NC3DTEN LOCAL ARRAY INITIALIZED
            morr_arr(i,j,k,MORRInd::nc3dten) = 0.0;

            // INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO
            c2prec = 0.0;
            csed = 0.0;
            ised = 0.0;
            ssed = 0.0;
            gsed = 0.0;
            rsed = 0.0;

            // LATENT HEAT OF VAPORIZATION
            morr_arr(i,j,k,MORRInd::xxlv) = 3.1484E6 - 2370.0 * morr_arr(i,j,k,MORRInd::t3d);
            // LATENT HEAT OF SUBLIMATION
            morr_arr(i,j,k,MORRInd::xxls) = 3.15E6 - 2370.0 * morr_arr(i,j,k,MORRInd::t3d) + 0.3337E6;

            // Assuming CP is a constant defined elsewhere (specific heat of dry air at constant pressure)
            const amrex::Real CP = 1004.5; // J/kg/K
            morr_arr(i,j,k,MORRInd::cpm) = CP * (1.0 + 0.887 * morr_arr(i,j,k,MORRInd::qv3d));

            // SATURATION VAPOR PRESSURE AND MIXING RATIO
            // hm, add fix for low pressure, 5/12/10
            // Assuming POLYSVP is defined elsewhere
            evs = std::min(0.99 * morr_arr(i,j,k,MORRInd::pres), calc_saturation_vapor_pressure(morr_arr(i,j,k,MORRInd::t3d), 0));  // PA
            eis = std::min(0.99 * morr_arr(i,j,k,MORRInd::pres), calc_saturation_vapor_pressure(morr_arr(i,j,k,MORRInd::t3d), 1));  // PA
            // MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING
            if (eis > evs) {
              eis = evs; // temporary update: adjust ice saturation pressure
            }

            // SATURATION MIXING RATIOS
            qvs = m_ep_2 * evs / (morr_arr(i,j,k,MORRInd::pres) - evs); // budget equation: calculate water saturation mixing ratio
            qvi = m_ep_2 * eis / (morr_arr(i,j,k,MORRInd::pres) - eis); // budget equation: calculate ice saturation mixing ratio

            // SATURATION RATIOS
            qvqvs = morr_arr(i,j,k,MORRInd::qv3d) / qvs; // budget equation: calculate water saturation ratio
            qvqvsi = morr_arr(i,j,k,MORRInd::qv3d) / qvi; // budget equation: calculate ice saturation ratio

            // AIR DENSITY
            morr_arr(i,j,k,MORRInd::rho) = morr_arr(i,j,k,MORRInd::pres) / (m_R * morr_arr(i,j,k,MORRInd::t3d)); // budget equation: calculate air density

            ds0 = 3.0;       // Size distribution parameter for snow
            di0 = 3.0;       // Size distribution parameter for cloud ice
            dg0 = 3.0;       // Size distribution parameter for graupel
            const double CI = 800.0;     // Mass-diameter relationship parameter for cloud ice
            // ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
            // ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
            // ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
            // FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON
            if (morr_arr(i,j,k,MORRInd::qrcu1d) >= 1.0e-10) {
              dum = 1.8e5 * std::pow(morr_arr(i,j,k,MORRInd::qrcu1d) * dt / (m_pi * m_rhow * std::pow(morr_arr(i,j,k,MORRInd::rho), 3)), 0.25); // rate equation: calculate rain number concentration from cumulus
              morr_arr(i,j,k,MORRInd::nr3d) += dum; // budget equation: update rain number concentration
            }
            if (morr_arr(i,j,k,MORRInd::qscu1d) >= 1.0e-10) {
              dum = 3.e5 * std::pow(morr_arr(i,j,k,MORRInd::qscu1d) * dt / (m_cons1 * std::pow(morr_arr(i,j,k,MORRInd::rho), 3)), 1.0 / (ds0 + 1.0)); // rate equation: calculate snow number concentration from cumulus
              morr_arr(i,j,k,MORRInd::ns3d) += dum; // budget equation: update snow number concentration
            }
            if (morr_arr(i,j,k,MORRInd::qicu1d) >= 1.0e-10) {
              dum = morr_arr(i,j,k,MORRInd::qicu1d) * dt / (CI * std::pow(80.0e-6, di0)); // rate equation: calculate cloud ice number concentration from cumulus
              morr_arr(i,j,k,MORRInd::ni3d) += dum; // budget equation: update cloud ice number concentration
            }

            // AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
            // hm modify 7/0/09 change limit to 1.e-8
            if (qvqvs < 0.9) {
              if (morr_arr(i,j,k,MORRInd::qr3d) < 1.0e-8) {
                morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qr3d); // budget equation: transfer rain to vapor
                morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qr3d) * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: adjust temperature
                morr_arr(i,j,k,MORRInd::qr3d) = 0.0; // temporary update: set rain to zero
              }
              if (morr_arr(i,j,k,MORRInd::qc3d) < 1.0e-8) {
                morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qc3d); // budget equation: transfer cloud water to vapor
                morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: adjust temperature
                morr_arr(i,j,k,MORRInd::qc3d) = 0.0; // temporary update: set cloud water to zero
              }
            }
            if (qvqvsi < 0.9) {
              if (morr_arr(i,j,k,MORRInd::qi3d) < 1.0e-8) {
                morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qi3d); // budget equation: transfer cloud ice to vapor
                morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qi3d) * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: adjust temperature
                morr_arr(i,j,k,MORRInd::qi3d) = 0.0; // temporary update: set cloud ice to zero
              }
              if (morr_arr(i,j,k,MORRInd::qni3d) < 1.0e-8) {
                morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qni3d); // budget equation: transfer snow to vapor
                morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qni3d) * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: adjust temperature
                morr_arr(i,j,k,MORRInd::qni3d) = 0.0; // temporary update: set snow to zero
              }
              if (morr_arr(i,j,k,MORRInd::qg3d) < 1.0e-8) {
                morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qg3d); // budget equation: transfer graupel to vapor
                morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qg3d) * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: adjust temperature
                morr_arr(i,j,k,MORRInd::qg3d) = 0.0; // temporary update: set graupel to zero
              }
            }
            // HEAT OF FUSION
            morr_arr(i,j,k,MORRInd::xlf) = morr_arr(i,j,k,MORRInd::xxls) - morr_arr(i,j,k,MORRInd::xxlv);

            // IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO
            // Note: QSMALL is not defined in the variable list, so I'll define it
            const amrex::Real QSMALL = m_qsmall;

            if (morr_arr(i,j,k,MORRInd::qc3d) < QSMALL) {
              morr_arr(i,j,k,MORRInd::qc3d) = 0.0;
              morr_arr(i,j,k,MORRInd::nc3d) = 0.0;
              morr_arr(i,j,k,MORRInd::effc) = 0.0;
            }
            if (morr_arr(i,j,k,MORRInd::qr3d) < QSMALL) {
              morr_arr(i,j,k,MORRInd::qr3d) = 0.0;
              morr_arr(i,j,k,MORRInd::nr3d) = 0.0;
              morr_arr(i,j,k,MORRInd::effr) = 0.0;
            }
            if (morr_arr(i,j,k,MORRInd::qi3d) < QSMALL) {
              morr_arr(i,j,k,MORRInd::qi3d) = 0.0;
              morr_arr(i,j,k,MORRInd::ni3d) = 0.0;
              morr_arr(i,j,k,MORRInd::effi) = 0.0;
            }
            if (morr_arr(i,j,k,MORRInd::qni3d) < QSMALL) {
              morr_arr(i,j,k,MORRInd::qni3d) = 0.0;
              morr_arr(i,j,k,MORRInd::ns3d) = 0.0;
              morr_arr(i,j,k,MORRInd::effs) = 0.0;
            }
            if (morr_arr(i,j,k,MORRInd::qg3d) < QSMALL) {
              morr_arr(i,j,k,MORRInd::qg3d) = 0.0;
              morr_arr(i,j,k,MORRInd::ng3d) = 0.0;
              morr_arr(i,j,k,MORRInd::effg) = 0.0;
            }
            // INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO
            morr_arr(i,j,k,MORRInd::qrsten) = 0.0;  // temporary update: initialize QRSTEN
            morr_arr(i,j,k,MORRInd::qisten) = 0.0;  // temporary update: initialize QISTEN
            morr_arr(i,j,k,MORRInd::qnisten) = 0.0; // temporary update: initialize QNISTEN
            morr_arr(i,j,k,MORRInd::qcsten) = 0.0;  // temporary update: initialize QCSTEN
            morr_arr(i,j,k,MORRInd::qgsten) = 0.0;  // temporary update: initialize QGSTEN

            // MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT
            morr_arr(i,j,k,MORRInd::mu) = 1.496e-6 * std::pow(morr_arr(i,j,k,MORRInd::t3d), 1.5) / (morr_arr(i,j,k,MORRInd::t3d) + 120.0); // budget equation: calculate air viscosity

            // Fall speed with density correction (Heymsfield and Benssemer 2006)
            dum = std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.54); // temporary update: calculate density correction factor

            // AA revision 4/1/11: Ikawa and Saito 1991 air-density correction
            morr_arr(i,j,k,MORRInd::ain) = std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.35) * m_ai; // budget equation: calculate ice fall speed parameter
            morr_arr(i,j,k,MORRInd::arn) = dum * m_ar; // budget equation: calculate rain fall speed parameter
            morr_arr(i,j,k,MORRInd::asn) = dum * m_as; // budget equation: calculate snow fall speed parameter

            // AA revision 4/1/11: temperature-dependent Stokes fall speed
            morr_arr(i,j,k,MORRInd::acn) = m_g * m_rhow / (18.0 * morr_arr(i,j,k,MORRInd::mu)); // budget equation: calculate cloud droplet fall speed parameter

            // HM ADD GRAUPEL 8/28/06
            morr_arr(i,j,k,MORRInd::agn) = dum * m_ag; // budget equation: calculate graupel fall speed parameter
            // hm 4/7/09 bug fix, initialize morr_arr(i,j,k,MORRInd::lami) to prevent later division by zero
            morr_arr(i,j,k,MORRInd::lami) = 0.0; // temporary update: initialize LAMI

            // If there is no cloud/precip water, and if subsaturated, then skip microphysics for this level
            bool skipMicrophysics = false;
            bool skipConcentrations = false;
            if (morr_arr(i,j,k,MORRInd::qc3d) < QSMALL && morr_arr(i,j,k,MORRInd::qi3d) < QSMALL && morr_arr(i,j,k,MORRInd::qni3d) < QSMALL && morr_arr(i,j,k,MORRInd::qr3d) < QSMALL && morr_arr(i,j,k,MORRInd::qg3d) < QSMALL) {
              if ((morr_arr(i,j,k,MORRInd::t3d) < 273.15 && qvqvsi < 0.999) || (morr_arr(i,j,k,MORRInd::t3d) >= 273.15 && qvqvs < 0.999)) {
                skipMicrophysics = true;//                goto label_200;
              }
            }

            if(!skipMicrophysics) {

            // Thermal conductivity for air
              kap = 1.414e3 * morr_arr(i,j,k,MORRInd::mu); // budget equation: calculate thermal conductivity

            // Diffusivity of water vapor
            dv = 8.794e-5 * std::pow(morr_arr(i,j,k,MORRInd::t3d), 1.81) / morr_arr(i,j,k,MORRInd::pres); // budget equation: calculate vapor diffusivity

            // Schmidt number
            sc_schmidt = morr_arr(i,j,k,MORRInd::mu) / (morr_arr(i,j,k,MORRInd::rho) * dv); // budget equation: calculate Schmidt number

            // Psychometric corrections
            // Rate of change sat. mix. ratio with temperature
            dum = (m_Rv * std::pow(morr_arr(i,j,k,MORRInd::t3d),2)); // temporary update: calculate temperature factor
            dqsdt = morr_arr(i,j,k,MORRInd::xxlv) * qvs / dum; // budget equation: calculate DQSDT
            dqsidt = morr_arr(i,j,k,MORRInd::xxls) * qvi / dum; // budget equation: calculate DQSIDT
            abi = 1.0 + dqsidt * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: calculate ABI
            ab = 1.0 + dqsdt * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm); // budget equation: calculate AB

            // CASE FOR TEMPERATURE ABOVE FREEZING
            if (morr_arr(i,j,k,MORRInd::t3d) >= 273.15) {
              //......................................................................
              // ALLOW FOR CONSTANT DROPLET NUMBER
              // INUM = 0, PREDICT DROPLET NUMBER
              // INUM = 1, SET CONSTANT DROPLET NUMBER

              if (m_inum == 1) {
                // CONVERT NDCNST FROM CM-3 TO KG-1
                // Note: NDCNST constant would need to be defined elsewhere
                morr_arr(i,j,k,MORRInd::nc3d) = m_ndcnst * 1.0e6 / morr_arr(i,j,k,MORRInd::rho); // Set cloud droplet number concentration
              }

              // GET SIZE DISTRIBUTION PARAMETERS
              // MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
              if (morr_arr(i,j,k,MORRInd::qni3d) < 1.0e-6) {
                morr_arr(i,j,k,MORRInd::qr3d) = morr_arr(i,j,k,MORRInd::qr3d) + morr_arr(i,j,k,MORRInd::qni3d);         // Transfer snow to rain
                morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::nr3d) + morr_arr(i,j,k,MORRInd::ns3d);          // Transfer snow number to rain
                morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) - morr_arr(i,j,k,MORRInd::qni3d) * morr_arr(i,j,k,MORRInd::xlf) / morr_arr(i,j,k,MORRInd::cpm); // Adjust temperature
                morr_arr(i,j,k,MORRInd::qni3d) = 0.0;                 // Set snow to zero
                morr_arr(i,j,k,MORRInd::ns3d) = 0.0;                  // Set snow number to zero
              }

              if (morr_arr(i,j,k,MORRInd::qg3d) < 1.0e-6) {
                morr_arr(i,j,k,MORRInd::qr3d) = morr_arr(i,j,k,MORRInd::qr3d) + morr_arr(i,j,k,MORRInd::qg3d);          // Transfer graupel to rain
                morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::nr3d) + morr_arr(i,j,k,MORRInd::ng3d);          // Transfer graupel number to rain
                morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) - morr_arr(i,j,k,MORRInd::qg3d) * morr_arr(i,j,k,MORRInd::xlf) / morr_arr(i,j,k,MORRInd::cpm);  // Adjust temperature
                morr_arr(i,j,k,MORRInd::qg3d) = 0.0;                  // Set graupel to zero
                morr_arr(i,j,k,MORRInd::ng3d) = 0.0;                  // Set graupel number to zero
              }
              // Skip to label 300 if concentrations are below thresholds
              if (morr_arr(i,j,k,MORRInd::qc3d) < m_qsmall && morr_arr(i,j,k,MORRInd::qni3d) < 1.0e-8 && morr_arr(i,j,k,MORRInd::qr3d) < m_qsmall && morr_arr(i,j,k,MORRInd::qg3d) < 1.0e-8) {
                skipConcentrations=true;//                goto label_300;
              }
              if(!skipConcentrations) {
                morr_arr(i,j,k,MORRInd::ns3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::ns3d));
                morr_arr(i,j,k,MORRInd::nc3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::nc3d));
                morr_arr(i,j,k,MORRInd::nr3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::nr3d));
                morr_arr(i,j,k,MORRInd::ng3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::ng3d));

                // ========================================================================
                // USING WRF APPROACH FOR SIZE DISTRIBUTION PARAMETERS
                // ========================================================================
                // Rain
                if (morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  // Calculate lambda parameter using cons26 (pi*rhow/6)
                  morr_arr(i,j,k,MORRInd::lamr) = pow(m_pi * m_rhow * morr_arr(i,j,k,MORRInd::nr3d) / morr_arr(i,j,k,MORRInd::qr3d), 1.0/3.0);
                  morr_arr(i,j,k,MORRInd::n0r) = morr_arr(i,j,k,MORRInd::nr3d)*morr_arr(i,j,k,MORRInd::lamr);

                  // Check for slope and adjust vars
                  if (morr_arr(i,j,k,MORRInd::lamr) < m_lamminr) {
                    morr_arr(i,j,k,MORRInd::lamr) = m_lamminr;
                    morr_arr(i,j,k,MORRInd::n0r) = pow(morr_arr(i,j,k,MORRInd::lamr), 4.0) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                    morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::n0r) / morr_arr(i,j,k,MORRInd::lamr);  // Update number concentration
                  } else if (morr_arr(i,j,k,MORRInd::lamr) > m_lammaxr) {
                    morr_arr(i,j,k,MORRInd::lamr) = m_lammaxr;
                    morr_arr(i,j,k,MORRInd::n0r) = pow(morr_arr(i,j,k,MORRInd::lamr), 4.0) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                    morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::n0r) / morr_arr(i,j,k,MORRInd::lamr);  // Update number concentration
                  }
                }

                // Cloud droplets
                if (morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  // Calculate air density factor (moist air density)
                  dum = morr_arr(i,j,k,MORRInd::pres)/(287.15*morr_arr(i,j,k,MORRInd::t3d));

                  // MARTIN ET AL. (1994) FORMULA FOR PGAM (WRF implementation)
                  morr_arr(i,j,k,MORRInd::pgam) = 0.0005714*(morr_arr(i,j,k,MORRInd::nc3d)/1.0e6*dum) + 0.2714;
                  morr_arr(i,j,k,MORRInd::pgam) = 1.0/(morr_arr(i,j,k,MORRInd::pgam)*morr_arr(i,j,k,MORRInd::pgam)) - 1.0;
                  morr_arr(i,j,k,MORRInd::pgam) = amrex::max(morr_arr(i,j,k,MORRInd::pgam), 2.0);
                  morr_arr(i,j,k,MORRInd::pgam) = amrex::min(morr_arr(i,j,k,MORRInd::pgam), 10.0);

                  // Calculate gamma function values
                  amrex::Real gamma_pgam_plus_1 = gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.0);
                  amrex::Real gamma_pgam_plus_4 = gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0);

                  // Calculate lambda parameter
                  morr_arr(i,j,k,MORRInd::lamc) = pow((m_cons26 * morr_arr(i,j,k,MORRInd::nc3d) * gamma_pgam_plus_4) / (morr_arr(i,j,k,MORRInd::qc3d) * gamma_pgam_plus_1), 1.0/3.0);

                  // Lambda bounds from WRF - 60 micron max diameter, 1 micron min diameter
                  amrex::Real lambda_min = (morr_arr(i,j,k,MORRInd::pgam) + 1.0)/60.0e-6;
                  amrex::Real lambda_max = (morr_arr(i,j,k,MORRInd::pgam) + 1.0)/1.0e-6;

                  // Check bounds and update number concentration if needed
                  if (morr_arr(i,j,k,MORRInd::lamc) < lambda_min) {
                    morr_arr(i,j,k,MORRInd::lamc) = lambda_min;
                    // Update cloud droplet number using the same formula as in WRF
                    morr_arr(i,j,k,MORRInd::nc3d) = exp(3.0*log(morr_arr(i,j,k,MORRInd::lamc)) + log(morr_arr(i,j,k,MORRInd::qc3d)) +
                               log(gamma_pgam_plus_1) - log(gamma_pgam_plus_4))/ m_cons26;
                  } else if (morr_arr(i,j,k,MORRInd::lamc) > lambda_max) {
                    morr_arr(i,j,k,MORRInd::lamc) = lambda_max;
                    // Update cloud droplet number using the same formula as in WRF
                    morr_arr(i,j,k,MORRInd::nc3d) = exp(3.0*log(morr_arr(i,j,k,MORRInd::lamc)) + log(morr_arr(i,j,k,MORRInd::qc3d)) +
                               log(gamma_pgam_plus_1) - log(gamma_pgam_plus_4))/ m_cons26;
                  }

                  // Calculate intercept parameter
                  morr_arr(i,j,k,MORRInd::cdist1) = morr_arr(i,j,k,MORRInd::nc3d) * pow(morr_arr(i,j,k,MORRInd::lamc), morr_arr(i,j,k,MORRInd::pgam)+1) / gamma_pgam_plus_1;
                }

                // Snow
                if (morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                  // Calculate lambda parameter
                  morr_arr(i,j,k,MORRInd::lams) = pow(m_cons1 * morr_arr(i,j,k,MORRInd::ns3d) / morr_arr(i,j,k,MORRInd::qni3d), 1.0/ds0);

                  // Calculate intercept parameter
                  morr_arr(i,j,k,MORRInd::n0s) = morr_arr(i,j,k,MORRInd::ns3d) * morr_arr(i,j,k,MORRInd::lams);

                  // Check for slope and adjust vars
                  if (morr_arr(i,j,k,MORRInd::lams) < m_lammins) {
                    morr_arr(i,j,k,MORRInd::lams) = m_lammins;
                    morr_arr(i,j,k,MORRInd::n0s) = pow(morr_arr(i,j,k,MORRInd::lams), 4.0) * morr_arr(i,j,k,MORRInd::qni3d) / m_cons1;
                    morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::n0s) / morr_arr(i,j,k,MORRInd::lams);  // Update number concentration
                  } else if (morr_arr(i,j,k,MORRInd::lams) > m_lammaxs) {
                    morr_arr(i,j,k,MORRInd::lams) = m_lammaxs;
                    morr_arr(i,j,k,MORRInd::n0s) = pow(morr_arr(i,j,k,MORRInd::lams), 4.0) * morr_arr(i,j,k,MORRInd::qni3d) / m_cons1;
                    morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::n0s) / morr_arr(i,j,k,MORRInd::lams);  // Update number concentration
                  }
                }

                // Graupel
                if (morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                  // Calculate lambda parameter
                  morr_arr(i,j,k,MORRInd::lamg) = pow(m_cons2 * morr_arr(i,j,k,MORRInd::ng3d) / morr_arr(i,j,k,MORRInd::qg3d), 1.0/dg0);

                  // Calculate intercept parameter
                  morr_arr(i,j,k,MORRInd::n0g) = morr_arr(i,j,k,MORRInd::ng3d) * morr_arr(i,j,k,MORRInd::lamg);

                  // Check for slope and adjust vars
                  if (morr_arr(i,j,k,MORRInd::lamg) < m_lamming) {
                    morr_arr(i,j,k,MORRInd::lamg) = m_lamming;
                    morr_arr(i,j,k,MORRInd::n0g) = pow(morr_arr(i,j,k,MORRInd::lamg), 4.0) * morr_arr(i,j,k,MORRInd::qg3d) / m_cons2;
                    morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::n0g) / morr_arr(i,j,k,MORRInd::lamg);  // Update number concentration
                  } else if (morr_arr(i,j,k,MORRInd::lamg) > m_lammaxg) {
                    morr_arr(i,j,k,MORRInd::lamg) = m_lammaxg;
                    morr_arr(i,j,k,MORRInd::n0g) = pow(morr_arr(i,j,k,MORRInd::lamg), 4.0) * morr_arr(i,j,k,MORRInd::qg3d) / m_cons2;
                    morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::n0g) / morr_arr(i,j,k,MORRInd::lamg);  // Update number concentration
                  }
                }
                ////////////////////// First instance of ZERO OUT PROCESS RATES
                // Zero out process rates
                prc = 0.0;         // Cloud water to rain conversion rate (PRC)
                nprc = 0.0;        // Change in cloud droplet number due to autoconversion (NPRC)
                nprc1 = 0.0;       // Change in rain number due to autoconversion (NPRC1)
                pra = 0.0;         // Accretion of cloud water by rain (PRA)
                npra = 0.0;        // Change in cloud droplet number due to accretion by rain (NPRA)
                nragg = 0.0;       // Self-collection/breakup of rain (NRAGG)
                nsmlts = 0.0;      // Loss of snow number during melting (NSMLTS)
                nsmltr = 0.0;      // Change in rain number due to snow melting (NSMLTR)
                evpms = 0.0;       // Melting snow evaporation rate (EVPMS)
                pcc = 0.0;         // Condensation/evaporation of cloud water (PCC)
                pre = 0.0;         // Evaporation of rain (PRE)
                nsubc = 0.0;       // Loss of cloud droplet number during evaporation (NSUBC)
                nsubr = 0.0;       // Loss of rain number during evaporation (NSUBR)
                pracg = 0.0;       // Collection of rain by graupel (PRACG)
                npracg = 0.0;      // Change in number due to collection of rain by graupel (NPRACG)
                psmlt = 0.0;       // Melting of snow (PSMLT)
                pgmlt = 0.0;       // Melting of graupel (PGMLT)
                evpmg = 0.0;       // Evaporation of melting graupel (EVPMG)
                pracs = 0.0;       // Collection of snow by rain (PRACS)
                npracs = 0.0;      // Change in number due to collection of snow by rain (NPRACS)
                ngmltg = 0.0;      // Loss of graupel number during melting (NGMLTG)
                ngmltr = 0.0;      // Change in rain number due to graupel melting (NGMLTR)

                // CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

                // AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
                // FORMULA FROM BEHENG (1994)
                // USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
                // AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
                // AS A GAMMA DISTRIBUTION

                // USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

                if (morr_arr(i,j,k,MORRInd::qc3d) >= 1.0e-6) {
                  // HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
                  // FROM KHAIROUTDINOV AND KOGAN 2000, MWR
                  prc = 1350.0 * std::pow(morr_arr(i,j,k,MORRInd::qc3d), 2.47) *
                    std::pow((morr_arr(i,j,k,MORRInd::nc3d)/1.0e6*morr_arr(i,j,k,MORRInd::rho)), -1.79);

                  // note: nprc1 is change in Nr,
                  // nprc is change in Nc
                  nprc1 = prc / m_cons29;
                  nprc = prc / (morr_arr(i,j,k,MORRInd::qc3d) / morr_arr(i,j,k,MORRInd::nc3d));

                  // hm bug fix 3/20/12
                  nprc = std::min(nprc, morr_arr(i,j,k,MORRInd::nc3d) / dt);
                  nprc1 = std::min(nprc1, nprc);
                }

                // HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
                // FORMULA FROM IKAWA AND SAITO (1991)

                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qni3d) >= 1.0e-8) {
                  amrex::Real ums_local = morr_arr(i,j,k,MORRInd::asn) * m_cons3 / std::pow(morr_arr(i,j,k,MORRInd::lams), m_bs);
                  amrex::Real umr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons4 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);
                  amrex::Real uns_local = morr_arr(i,j,k,MORRInd::asn) * m_cons5 / std::pow(morr_arr(i,j,k,MORRInd::lams), m_bs);
                  amrex::Real unr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons6 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);

                  // SET REALISTIC LIMITS ON FALLSPEEDS
                  // bug fix, 10/08/09
                  dum = std::pow(m_rhosu/morr_arr(i,j,k,MORRInd::rho), 0.54);
                  ums_local = std::min(ums_local, 1.2*dum);
                  uns_local = std::min(uns_local, 1.2*dum);
                  umr_local = std::min(umr_local, 9.1*dum);
                  unr_local = std::min(unr_local, 9.1*dum);


                  // hm fix, 2/12/13
                  // for above freezing conditions to get accelerated melting of snow,
                  // we need collection of rain by snow (following Lin et al. 1983)
                  ////////////////////////Might need pow expanding
                  pracs = m_cons41 * (std::sqrt(std::pow(1.2*umr_local-0.95*ums_local, 2) +
                                                0.08*ums_local*umr_local) * morr_arr(i,j,k,MORRInd::rho) *
                                      morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0s) / std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) *
                                      (5.0/(std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lams)) +
                                       2.0/(std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lams), 2)) +
                                       0.5/(morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lams), 3))));
                }
                // ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
                // ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
                // ASSUME SHED DROPS ARE 1 MM IN SIZE

                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qg3d) >= 1.0e-8) {

                  amrex::Real umg_local = morr_arr(i,j,k,MORRInd::agn) * m_cons7 / std::pow(morr_arr(i,j,k,MORRInd::lamg), m_bg);
                  amrex::Real umr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons4 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);
                  amrex::Real ung_local = morr_arr(i,j,k,MORRInd::agn) * m_cons8 / std::pow(morr_arr(i,j,k,MORRInd::lamg), m_bg);
                  amrex::Real unr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons6 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);

                  // SET REALISTIC LIMITS ON FALLSPEEDS
                  // bug fix, 10/08/09
                  dum = std::pow(m_rhosu/morr_arr(i,j,k,MORRInd::rho), 0.54);
                  umg_local = std::min(umg_local, 20.0*dum);
                  ung_local = std::min(ung_local, 20.0*dum);
                  umr_local = std::min(umr_local, 9.1*dum);
                  unr_local = std::min(unr_local, 9.1*dum);

                  // PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
                  pracg = m_cons41 * (std::sqrt(std::pow(1.2*umr_local-0.95*umg_local, 2) +
                                                0.08*umg_local*umr_local) * morr_arr(i,j,k,MORRInd::rho) *
                                      morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0g) / std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) *
                                      (5.0/(std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lamg)) +
                                       2.0/(std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 2)) +
                                       0.5/(morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 3))));

                  // ASSUME 1 MM DROPS ARE SHED, GET NUMBER SHED PER SEC
                  dum = pracg/5.2e-7;

                  npracg = m_cons32 * morr_arr(i,j,k,MORRInd::rho) * (std::sqrt(1.7*std::pow(unr_local-ung_local, 2) +
                                                              0.3*unr_local*ung_local) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0g) *
                                                    (1.0/(std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lamg)) +
                                                     1.0/(std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 2)) +
                                                     1.0/(morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 3))));
                  // hm 7/15/13, remove limit so that the number of collected drops can smaller than
                  // number of shed drops
                  npracg = npracg - dum;
                }
                // ACCRETION OF CLOUD LIQUID WATER BY RAIN
                // CONTINUOUS COLLECTION EQUATION WITH
                // GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qc3d) >= 1.0e-8) {
                  // 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
                  // KHAIROUTDINOV AND KOGAN 2000, MWR
                  dum = morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::qr3d);
                  pra = 67.0 * std::pow(dum, 1.15);
                  npra = pra / (morr_arr(i,j,k,MORRInd::qc3d) / morr_arr(i,j,k,MORRInd::nc3d));
                }

                // SELF-COLLECTION OF RAIN DROPS
                // FROM BEHENG(1994)
                // FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
                // AS DESCRIBED ABOVE FOR AUTOCONVERSION

                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8) {
                  // include breakup add 10/09/09
                  dum1 = 300.0e-6;
                  if (1.0/morr_arr(i,j,k,MORRInd::lamr) < dum1) {
                    dum = 1.0;
                  } else {
                    dum = 2.0 - std::exp(2300.0 * (1.0/morr_arr(i,j,k,MORRInd::lamr) - dum1));
                  }
                  nragg = -5.78 * dum * morr_arr(i,j,k,MORRInd::nr3d) * morr_arr(i,j,k,MORRInd::qr3d) * morr_arr(i,j,k,MORRInd::rho);
                }
                // CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)
                if (morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  epsr = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::rho) * dv *
                    (m_f1r/(morr_arr(i,j,k,MORRInd::lamr)*morr_arr(i,j,k,MORRInd::lamr)) +
                     m_f2r * std::sqrt(morr_arr(i,j,k,MORRInd::arn)*morr_arr(i,j,k,MORRInd::rho)/morr_arr(i,j,k,MORRInd::mu)) *
                     std::pow(sc_schmidt, 1.0/3.0) * m_cons9 /
                     std::pow(morr_arr(i,j,k,MORRInd::lamr), m_cons34));
                } else {
                  epsr = 0.0;
                }
                // NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED
                if (morr_arr(i,j,k,MORRInd::qv3d) < qvs) {
                  pre = epsr * (morr_arr(i,j,k,MORRInd::qv3d) - qvs) / ab;
                  pre = std::min(pre, 0.0);
                } else {
                  pre = 0.0;
                }
                // MELTING OF SNOW
                // SNOW MAY PERSIST ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
                // IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

                if (morr_arr(i,j,k,MORRInd::qni3d) >= 1.0e-8) {
                  // fix 053011
                  // HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
                  dum = -m_cpw/morr_arr(i,j,k,MORRInd::xlf) * (morr_arr(i,j,k,MORRInd::t3d) - 273.15) * pracs;

                  // hm fix 1/20/15
                  psmlt = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0s) * kap * (273.15 - morr_arr(i,j,k,MORRInd::t3d)) /
                    morr_arr(i,j,k,MORRInd::xlf) * (m_f1s/(morr_arr(i,j,k,MORRInd::lams)*morr_arr(i,j,k,MORRInd::lams)) +
                           m_f2s * std::sqrt(morr_arr(i,j,k,MORRInd::asn)*morr_arr(i,j,k,MORRInd::rho)/morr_arr(i,j,k,MORRInd::mu)) *
                           std::pow(sc_schmidt, 1.0/3.0) * m_cons10 /
                           std::pow(morr_arr(i,j,k,MORRInd::lams), m_cons35)) + dum;

                  // IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES
                  if (qvqvs < 1.0) {
                    epss = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0s) * morr_arr(i,j,k,MORRInd::rho) * dv *
                      (m_f1s/(morr_arr(i,j,k,MORRInd::lams)*morr_arr(i,j,k,MORRInd::lams)) +
                       m_f2s * std::sqrt(morr_arr(i,j,k,MORRInd::asn)*morr_arr(i,j,k,MORRInd::rho)/morr_arr(i,j,k,MORRInd::mu)) *
                       std::pow(sc_schmidt, 1.0/3.0) * m_cons10 /
                       std::pow(morr_arr(i,j,k,MORRInd::lams), m_cons35));

                    // hm fix 8/4/08
                    evpms = (morr_arr(i,j,k,MORRInd::qv3d) - qvs) * epss / ab;
                    evpms = std::max(evpms, psmlt);
                    psmlt = psmlt - evpms;
                  }
                }
                // MELTING OF GRAUPEL
                // GRAUPEL MAY PERSIST ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
                // IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

                if (morr_arr(i,j,k,MORRInd::qg3d) >= 1.0e-8) {
                  // fix 053011
                  // HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN

                  dum = -m_cpw/morr_arr(i,j,k,MORRInd::xlf) * (morr_arr(i,j,k,MORRInd::t3d) - 273.15) * pracg;

                  // hm fix 1/20/15
                  pgmlt = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0g) * kap * (273.15 - morr_arr(i,j,k,MORRInd::t3d)) /
                    morr_arr(i,j,k,MORRInd::xlf) * (m_f1s/(morr_arr(i,j,k,MORRInd::lamg)*morr_arr(i,j,k,MORRInd::lamg)) +
                           m_f2s * std::sqrt(morr_arr(i,j,k,MORRInd::agn)*morr_arr(i,j,k,MORRInd::rho)/morr_arr(i,j,k,MORRInd::mu)) *
                           std::pow(sc_schmidt, 1.0/3.0) * m_cons11 /
                           std::pow(morr_arr(i,j,k,MORRInd::lamg), m_cons36)) + dum;

                  // IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES
                  if (qvqvs < 1.0) {
                    epsg = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0g) * morr_arr(i,j,k,MORRInd::rho) * dv *
                      (m_f1s/(morr_arr(i,j,k,MORRInd::lamg)*morr_arr(i,j,k,MORRInd::lamg)) +
                       m_f2s * std::sqrt(morr_arr(i,j,k,MORRInd::agn)*morr_arr(i,j,k,MORRInd::rho)/morr_arr(i,j,k,MORRInd::mu)) *
                       std::pow(sc_schmidt, 1.0/3.0) * m_cons11 /
                       std::pow(morr_arr(i,j,k,MORRInd::lamg), m_cons36));

                    // hm fix 8/4/08
                    evpmg = (morr_arr(i,j,k,MORRInd::qv3d) - qvs) * epsg / ab;
                    evpmg = std::max(evpmg, pgmlt);
                    pgmlt = pgmlt - evpmg;
                  }
                }
                // HM, V3.2
                // RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
                // TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
                // ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

                pracg = 0.0;
                pracs = 0.0;
                // CONSERVATION OF QC
                dum = (prc + pra) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qc3d) && morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  ratio = morr_arr(i,j,k,MORRInd::qc3d) / dum;
                  prc = prc * ratio;
                  pra = pra * ratio;
                }

                // CONSERVATION OF SNOW
                dum = (-psmlt - evpms + pracs) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qni3d) && morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                  // NO SOURCE TERMS FOR SNOW AT T > FREEZING
                  ratio = morr_arr(i,j,k,MORRInd::qni3d) / dum;
                  psmlt = psmlt * ratio;
                  evpms = evpms * ratio;
                  pracs = pracs * ratio;
                }

                // CONSERVATION OF GRAUPEL
                dum = (-pgmlt - evpmg + pracg) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qg3d) && morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                  // NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
                  ratio = morr_arr(i,j,k,MORRInd::qg3d) / dum;
                  pgmlt = pgmlt * ratio;
                  evpmg = evpmg * ratio;
                  pracg = pracg * ratio;
                }

                // CONSERVATION OF QR
                // HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

                dum = (-pracs - pracg - pre - pra - prc + psmlt + pgmlt) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qr3d) && morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  ratio = (morr_arr(i,j,k,MORRInd::qr3d)/dt + pracs + pracg + pra + prc - psmlt - pgmlt) / (-pre);
                  pre = pre * ratio;
                }
                // Update tendencies
                morr_arr(i,j,k,MORRInd::qv3dten) = morr_arr(i,j,k,MORRInd::qv3dten) + (-pre - evpms - evpmg);

                morr_arr(i,j,k,MORRInd::t3dten) = morr_arr(i,j,k,MORRInd::t3dten) + (pre * morr_arr(i,j,k,MORRInd::xxlv) +
                                                 (evpms + evpmg) * morr_arr(i,j,k,MORRInd::xxls) +
                                                 (psmlt + pgmlt - pracs - pracg) * morr_arr(i,j,k,MORRInd::xlf)) / morr_arr(i,j,k,MORRInd::cpm);

                morr_arr(i,j,k,MORRInd::qc3dten) = morr_arr(i,j,k,MORRInd::qc3dten) + (-pra - prc);
                morr_arr(i,j,k,MORRInd::qr3dten) = morr_arr(i,j,k,MORRInd::qr3dten) + (pre + pra + prc - psmlt - pgmlt + pracs + pracg);
                morr_arr(i,j,k,MORRInd::qni3dten) = morr_arr(i,j,k,MORRInd::qni3dten) + (psmlt + evpms - pracs);
                morr_arr(i,j,k,MORRInd::qg3dten) = morr_arr(i,j,k,MORRInd::qg3dten) + (pgmlt + evpmg - pracg);

                // fix 053011
                // HM, bug fix 5/12/08, npracg is subtracted from nr not ng
                morr_arr(i,j,k,MORRInd::nc3dten) = morr_arr(i,j,k,MORRInd::nc3dten) + (-npra - nprc);
                morr_arr(i,j,k,MORRInd::nr3dten) = morr_arr(i,j,k,MORRInd::nr3dten) + (nprc1 + nragg - npracg);

                // HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC
                c2prec = pra + prc;

                if (pre < 0.0) {
                  dum = pre * dt / morr_arr(i,j,k,MORRInd::qr3d);
                  dum = std::max(-1.0, dum);
                  nsubr = dum * morr_arr(i,j,k,MORRInd::nr3d) / dt;
                }

                if (evpms + psmlt < 0.0) {
                  dum = (evpms + psmlt) * dt / morr_arr(i,j,k,MORRInd::qni3d);
                  dum = std::max(-1.0, dum);
                  nsmlts = dum * morr_arr(i,j,k,MORRInd::ns3d) / dt;
                }

                if (psmlt < 0.0) {
                  dum = psmlt * dt / morr_arr(i,j,k,MORRInd::qni3d);
                  dum = std::max(-1.0, dum);
                  nsmltr = dum * morr_arr(i,j,k,MORRInd::ns3d) / dt;
                }

                if (evpmg + pgmlt < 0.0) {
                  dum = (evpmg + pgmlt) * dt / morr_arr(i,j,k,MORRInd::qg3d);
                  dum = std::max(-1.0, dum);
                  ngmltg = dum * morr_arr(i,j,k,MORRInd::ng3d) / dt;
                }

                if (pgmlt < 0.0) {
                  dum = pgmlt * dt / morr_arr(i,j,k,MORRInd::qg3d);
                  dum = std::max(-1.0, dum);
                  ngmltr = dum * morr_arr(i,j,k,MORRInd::ng3d) / dt;
                }

                morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) + nsmlts;
                morr_arr(i,j,k,MORRInd::ng3dten) = morr_arr(i,j,k,MORRInd::ng3dten) + ngmltg;
                morr_arr(i,j,k,MORRInd::nr3dten) = morr_arr(i,j,k,MORRInd::nr3dten) + (nsubr - nsmltr - ngmltr);

              }
              //Right after 300 CONTINUE
//            label_300:
              // Calculate saturation adjustment to condense extra vapor above water saturation
              dumt = morr_arr(i,j,k,MORRInd::t3d) + dt * morr_arr(i,j,k,MORRInd::t3dten);
              dumqv = morr_arr(i,j,k,MORRInd::qv3d) + dt * morr_arr(i,j,k,MORRInd::qv3dten);

              // Fix for low pressure (added 5/12/10)
              dum = std::min(0.99 * morr_arr(i,j,k,MORRInd::pres), calc_saturation_vapor_pressure(dumt, 0));
              dumqss = m_ep_2 * dum / (morr_arr(i,j,k,MORRInd::pres) - dum);
              dumqc = morr_arr(i,j,k,MORRInd::qc3d) + dt * morr_arr(i,j,k,MORRInd::qc3dten);
              dumqc = std::max(dumqc, 0.0);

              // Saturation adjustment for liquid
              dums = dumqv - dumqss;
              pcc = dums / (1.0 + std::pow(morr_arr(i,j,k,MORRInd::xxlv), 2) * dumqss / (morr_arr(i,j,k,MORRInd::cpm) * m_Rv * std::pow(dumt, 2))) / dt;
              if (pcc * dt + dumqc < 0.0) {
                pcc = -dumqc / dt;
              }

              if (!do_cond) { pcc = 0.0; }

              // Update tendencies
              morr_arr(i,j,k,MORRInd::qv3dten) -= pcc;
              morr_arr(i,j,k,MORRInd::t3dten)  += pcc * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm);
              morr_arr(i,j,k,MORRInd::qc3dten) += pcc;
            } else { //cold
              //......................................................................
              // ALLOW FOR CONSTANT DROPLET NUMBER
              // INUM = 0, PREDICT DROPLET NUMBER
              // INUM = 1, SET CONSTANT DROPLET NUMBER

              if (m_inum == 1) {
                // CONVERT NDCNST FROM CM-3 TO KG-1
                // Note: NDCNST constant would need to be defined elsewhere
                morr_arr(i,j,k,MORRInd::nc3d) = m_ndcnst * 1.0e6 / morr_arr(i,j,k,MORRInd::rho); // Set cloud droplet number concentration
              }

              morr_arr(i,j,k,MORRInd::ni3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::ni3d));
              morr_arr(i,j,k,MORRInd::ns3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::ns3d));
              morr_arr(i,j,k,MORRInd::nc3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::nc3d));
              morr_arr(i,j,k,MORRInd::nr3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::nr3d));
              morr_arr(i,j,k,MORRInd::ng3d) = amrex::max(0.0,morr_arr(i,j,k,MORRInd::ng3d));

              // ========================================================================
              // USING WRF APPROACH FOR SIZE DISTRIBUTION PARAMETERS
              // ========================================================================
              // Rain
              if (morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                // Calculate lambda parameter using cons26 (pi*rhow/6)
                morr_arr(i,j,k,MORRInd::lamr) = pow(m_pi * m_rhow * morr_arr(i,j,k,MORRInd::nr3d) / morr_arr(i,j,k,MORRInd::qr3d), 1.0/3.0);

                // Check for slope and adjust vars
                if (morr_arr(i,j,k,MORRInd::lamr) < m_lamminr) {
                  morr_arr(i,j,k,MORRInd::lamr) = m_lamminr;
                  morr_arr(i,j,k,MORRInd::n0r) = pow(morr_arr(i,j,k,MORRInd::lamr), 4.0) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                  morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::n0r) / morr_arr(i,j,k,MORRInd::lamr);  // Update number concentration
                } else if (morr_arr(i,j,k,MORRInd::lamr) > m_lammaxr) {
                  morr_arr(i,j,k,MORRInd::lamr) = m_lammaxr;
                  morr_arr(i,j,k,MORRInd::n0r) = pow(morr_arr(i,j,k,MORRInd::lamr), 4.0) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                  morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::n0r) / morr_arr(i,j,k,MORRInd::lamr);  // Update number concentration
                } else {
                  // Calculate intercept parameter using WRF formula
                  morr_arr(i,j,k,MORRInd::n0r) = pow(morr_arr(i,j,k,MORRInd::lamr), 4.0) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                }
              }


              // Cloud droplets
              if (morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                // Calculate air density factor (moist air density)
                dum = morr_arr(i,j,k,MORRInd::pres)/(287.15*morr_arr(i,j,k,MORRInd::t3d));

                // MARTIN ET AL. (1994) FORMULA FOR PGAM (WRF implementation)
                morr_arr(i,j,k,MORRInd::pgam) = 0.0005714*(morr_arr(i,j,k,MORRInd::nc3d)/1.0e6*dum) + 0.2714;
                morr_arr(i,j,k,MORRInd::pgam) = 1.0/(std::pow(morr_arr(i,j,k,MORRInd::pgam), 2)) - 1.0;
                morr_arr(i,j,k,MORRInd::pgam) = amrex::max(morr_arr(i,j,k,MORRInd::pgam), 2.0);
                morr_arr(i,j,k,MORRInd::pgam) = amrex::min(morr_arr(i,j,k,MORRInd::pgam), 10.0);

                // Calculate gamma function values
                amrex::Real gamma_pgam_plus_1 = gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.0);
                amrex::Real gamma_pgam_plus_4 = gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0);

                // Calculate lambda parameter
                morr_arr(i,j,k,MORRInd::lamc) = pow((m_cons26 * morr_arr(i,j,k,MORRInd::nc3d) * gamma_pgam_plus_4) / (morr_arr(i,j,k,MORRInd::qc3d) * gamma_pgam_plus_1), 1.0/3.0);

                // Lambda bounds from WRF - 60 micron max diameter, 1 micron min diameter
                amrex::Real lambda_min = (morr_arr(i,j,k,MORRInd::pgam) + 1.0)/60.0e-6;
                amrex::Real lambda_max = (morr_arr(i,j,k,MORRInd::pgam) + 1.0)/1.0e-6;

                // Check bounds and update number concentration if needed
                if (morr_arr(i,j,k,MORRInd::lamc) < lambda_min) {
                  morr_arr(i,j,k,MORRInd::lamc) = lambda_min;
                  // Update cloud droplet number using the same formula as in WRF
                  morr_arr(i,j,k,MORRInd::nc3d) = exp(3.0*log(morr_arr(i,j,k,MORRInd::lamc)) + log(morr_arr(i,j,k,MORRInd::qc3d)) +
                                    log(gamma_pgam_plus_1) - log(gamma_pgam_plus_4))/ m_cons26;
                } else if (morr_arr(i,j,k,MORRInd::lamc) > lambda_max) {
                  morr_arr(i,j,k,MORRInd::lamc) = lambda_max;
                  // Update cloud droplet number using the same formula as in WRF
                  morr_arr(i,j,k,MORRInd::nc3d) = exp(3.0*log(morr_arr(i,j,k,MORRInd::lamc)) + log(morr_arr(i,j,k,MORRInd::qc3d)) +
                                    log(gamma_pgam_plus_1) - log(gamma_pgam_plus_4))/ m_cons26;
                }

                // Calculate intercept parameter
                morr_arr(i,j,k,MORRInd::cdist1) = morr_arr(i,j,k,MORRInd::nc3d) / gamma_pgam_plus_1;
              }

              // Snow
              if (morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                // Calculate lambda parameter
                morr_arr(i,j,k,MORRInd::lams) = pow(m_cons1 * morr_arr(i,j,k,MORRInd::ns3d) / morr_arr(i,j,k,MORRInd::qni3d), 1.0/ds0);

                // Calculate intercept parameter
                morr_arr(i,j,k,MORRInd::n0s) = morr_arr(i,j,k,MORRInd::ns3d) * morr_arr(i,j,k,MORRInd::lams);

                // Check for slope and adjust vars
                if (morr_arr(i,j,k,MORRInd::lams) < m_lammins) {
                  morr_arr(i,j,k,MORRInd::lams) = m_lammins;
                  morr_arr(i,j,k,MORRInd::n0s) = pow(morr_arr(i,j,k,MORRInd::lams), 4.0) * morr_arr(i,j,k,MORRInd::qni3d) / m_cons1;
                  morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::n0s) / morr_arr(i,j,k,MORRInd::lams);  // Update number concentration
                } else if (morr_arr(i,j,k,MORRInd::lams) > m_lammaxs) {
                  morr_arr(i,j,k,MORRInd::lams) = m_lammaxs;
                  morr_arr(i,j,k,MORRInd::n0s) = pow(morr_arr(i,j,k,MORRInd::lams), 4.0) * morr_arr(i,j,k,MORRInd::qni3d) / m_cons1;
                  morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::n0s) / morr_arr(i,j,k,MORRInd::lams);  // Update number concentration
                }
              }

              // Cloud ice
              if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                // Calculate lambda parameter
                morr_arr(i,j,k,MORRInd::lami) = pow(m_cons12 * morr_arr(i,j,k,MORRInd::ni3d) / morr_arr(i,j,k,MORRInd::qi3d), 1.0/3.0);

                // Calculate intercept parameter (initial calculation)
                morr_arr(i,j,k,MORRInd::n0i) = morr_arr(i,j,k,MORRInd::ni3d) * morr_arr(i,j,k,MORRInd::lami);

                // Check for slope (apply bounds)
                if (morr_arr(i,j,k,MORRInd::lami) < m_lammini) {
                  morr_arr(i,j,k,MORRInd::lami) = m_lammini;
                  // Recalculate morr_arr(i,j,k,MORRInd::n0i) when lambda is adjusted
                  morr_arr(i,j,k,MORRInd::n0i) = pow(morr_arr(i,j,k,MORRInd::lami), 4.0) * morr_arr(i,j,k,MORRInd::qi3d) / m_cons12;
                  // Update ni3d when lambda is adjusted
                  morr_arr(i,j,k,MORRInd::ni3d) = morr_arr(i,j,k,MORRInd::n0i) / morr_arr(i,j,k,MORRInd::lami);
                } else if (morr_arr(i,j,k,MORRInd::lami) > m_lammaxi) {
                  morr_arr(i,j,k,MORRInd::lami) = m_lammaxi;
                  // Recalculate morr_arr(i,j,k,MORRInd::n0i) when lambda is adjusted
                  morr_arr(i,j,k,MORRInd::n0i) = pow(morr_arr(i,j,k,MORRInd::lami), 4.0) * morr_arr(i,j,k,MORRInd::qi3d) / m_cons12;
                  // Update ni3d when lambda is adjusted
                  morr_arr(i,j,k,MORRInd::ni3d) = morr_arr(i,j,k,MORRInd::n0i) / morr_arr(i,j,k,MORRInd::lami);
                }
              }
              // Graupel
              if (morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                // Calculate lambda parameter
                morr_arr(i,j,k,MORRInd::lamg) = pow(m_cons2 * morr_arr(i,j,k,MORRInd::ng3d) / morr_arr(i,j,k,MORRInd::qg3d), 1.0/dg0);

                // Calculate intercept parameter
                morr_arr(i,j,k,MORRInd::n0g) = morr_arr(i,j,k,MORRInd::ng3d) * morr_arr(i,j,k,MORRInd::lamg);

                // Check for slope and adjust vars
                if (morr_arr(i,j,k,MORRInd::lamg) < m_lamming) {
                  morr_arr(i,j,k,MORRInd::lamg) = m_lamming;
                  morr_arr(i,j,k,MORRInd::n0g) = pow(morr_arr(i,j,k,MORRInd::lamg), 4.0) * morr_arr(i,j,k,MORRInd::qg3d) / m_cons2;
                  morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::n0g) / morr_arr(i,j,k,MORRInd::lamg);  // Update number concentration
                } else if (morr_arr(i,j,k,MORRInd::lamg) > m_lammaxg) {
                  morr_arr(i,j,k,MORRInd::lamg) = m_lammaxg;
                  morr_arr(i,j,k,MORRInd::n0g) = pow(morr_arr(i,j,k,MORRInd::lamg), 4.0) * morr_arr(i,j,k,MORRInd::qg3d) / m_cons2;
                  morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::n0g) / morr_arr(i,j,k,MORRInd::lamg);  // Update number concentration
                }
              }
                ////////////////////// Second instance of ZERO OUT PROCESS RATES
                // Zero out process rates
                mnuccc = 0.0;      // Change Q due to contact freezing droplets (MNUCCC)
                nnuccc = 0.0;      // Change N due to contact freezing droplets (NNUCCC)
                prc = 0.0;         // Autoconversion droplets (PRC)
                nprc = 0.0;        // Change NC autoconversion droplets (NPRC)
                nprc1 = 0.0;       // Change NR autoconversion droplets (NPRC1)
                nsagg = 0.0;       // Self-collection of snow (NSAGG)
                psacws = 0.0;      // Change Q droplet accretion by snow (PSACWS)
                npsacws = 0.0;     // Change N droplet accretion by snow (NPSACWS)
                psacwi = 0.0;      // Change Q droplet accretion by cloud ice (PSACWI)
                npsacwi = 0.0;     // Change N droplet accretion by cloud ice (NPSACWI)
                pracs = 0.0;       // Change Q rain-snow collection (PRACS)
                npracs = 0.0;      // Change N rain-snow collection (NPRACS)
                nmults = 0.0;      // Ice multiplication due to riming droplets by snow (NMULTS)
                qmults = 0.0;      // Change Q due to ice multiplication droplets/snow (QMULTS)
                nmultr = 0.0;      // Ice multiplication due to riming rain by snow (NMULTR)
                qmultr = 0.0;      // Change Q due to ice multiplication rain/snow (QMULTR)
                nmultg = 0.0;      // Ice multiplication due to accretion droplets by graupel (NMULTG)
                qmultg = 0.0;      // Change Q due to ice multiplication droplets/graupel (QMULTG)
                nmultrg = 0.0;     // Ice multiplication due to accretion rain by graupel (NMULTRG)
                qmultrg = 0.0;     // Change Q due to ice multiplication rain/graupel (QMULTRG)
                mnuccr = 0.0;      // Change Q due to contact freezing rain (MNUCCR)
                nnuccr = 0.0;      // Change N due to contact freezing rain (NNUCCR)
                pra = 0.0;         // Accretion droplets by rain (PRA)
                npra = 0.0;        // Change N due to droplet accretion by rain (NPRA)
                nragg = 0.0;       // Self-collection/breakup of rain (NRAGG)
                prci = 0.0;        // Change Q autoconversion cloud ice to snow (PRCI)
                nprci = 0.0;       // Change N autoconversion cloud ice by snow (NPRCI)
                prai = 0.0;        // Change Q accretion cloud ice by snow (PRAI)
                nprai = 0.0;       // Change N accretion cloud ice (NPRAI)
                nnuccd = 0.0;      // Change N freezing aerosol (primary ice nucleation) (NNUCCD)
                mnuccd = 0.0;      // Change Q freezing aerosol (primary ice nucleation) (MNUCCD)
                pcc = 0.0;         // Condensation/evaporation droplets (PCC)
                pre = 0.0;         // Evaporation of rain (PRE)
                prd = 0.0;         // Deposition cloud ice (PRD)
                prds = 0.0;        // Deposition snow (PRDS)
                eprd = 0.0;        // Sublimation cloud ice (EPRD)
                eprds = 0.0;       // Sublimation snow (EPRDS)
                nsubc = 0.0;       // Loss of NC during evaporation (NSUBC)
                nsubi = 0.0;       // Loss of NI during sublimation (NSUBI)
                nsubs = 0.0;       // Loss of NS during sublimation (NSUBS)
                nsubr = 0.0;       // Loss of NR during evaporation (NSUBR)
                piacr = 0.0;       // Change QR, ice-rain collection (PIACR)
                niacr = 0.0;       // Change N, ice-rain collection (NIACR)
                praci = 0.0;       // Change QI, ice-rain collection (PRACI)
                piacrs = 0.0;      // Change QR, ice rain collision, added to snow (PIACRS)
                niacrs = 0.0;      // Change N, ice rain collision, added to snow (NIACRS)
                pracis = 0.0;      // Change QI, ice rain collision, added to snow (PRACIS)

                // Graupel processes
                pracg = 0.0;       // Change in Q collection rain by graupel (PRACG)
                psacr = 0.0;       // Conversion due to collection of snow by rain (PSACR)
                psacwg = 0.0;      // Change in Q collection droplets by graupel (PSACWG)
                pgsacw = 0.0;      // Conversion Q to graupel due to collection droplets by snow (PGSACW)
                pgracs = 0.0;      // Conversion Q to graupel due to collection rain by snow (PGRACS)
                prdg = 0.0;        // Deposition of graupel (PRDG)
                eprdg = 0.0;       // Sublimation of graupel (EPRDG)
                npracg = 0.0;      // Change N collection rain by graupel (NPRACG)
                npsacwg = 0.0;     // Change N collection droplets by graupel (NPSACWG)
                nscng = 0.0;       // Change N conversion to graupel due to collection droplets by snow (NSCNG)
                ngracs = 0.0;      // Change N conversion to graupel due to collection rain by snow (NGRACS)
                nsubg = 0.0;       // Change N sublimation/deposition of graupel (NSUBG)

                ////////////////////// CALCULATION OF MICROPHYSICAL PROCESS RATES
                // FREEZING OF CLOUD DROPLETS - ONLY ALLOWED BELOW -4C
                if (morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall && morr_arr(i,j,k,MORRInd::t3d) < 269.15) {
                  // NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
                  // FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3
                  // MEYERS CURVE
                  amrex::Real nacnt = std::exp(-2.80 + 0.262 * (273.15 - morr_arr(i,j,k,MORRInd::t3d))) * 1000.0;

                  // MEAN FREE PATH
                  dum = 7.37 * morr_arr(i,j,k,MORRInd::t3d) / (288.0 * 10.0 * morr_arr(i,j,k,MORRInd::pres)) / 100.0;

                  // EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
                  // BASED ON BROWNIAN DIFFUSION
                  amrex::Real dap = m_cons37 * morr_arr(i,j,k,MORRInd::t3d) * (1.0 + dum / m_rin) / morr_arr(i,j,k,MORRInd::mu);

                  // CONTACT FREEZING
                  mnuccc = m_cons38 * dap * nacnt * std::exp(std::log(morr_arr(i,j,k,MORRInd::cdist1)) +
                                                             std::log(gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 5.0)) - 4.0 * std::log(morr_arr(i,j,k,MORRInd::lamc)));
                  nnuccc = 2.0 * m_pi * dap * nacnt * morr_arr(i,j,k,MORRInd::cdist1) *
                    gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 2.0) / morr_arr(i,j,k,MORRInd::lamc);

                  // IMMERSION FREEZING (BIGG 1953)
                  // hm 7/15/13 fix for consistency w/ original formula
                  mnuccc = mnuccc + m_cons39 *
                    std::exp(std::log(morr_arr(i,j,k,MORRInd::cdist1)) + std::log(gamma_function(7.0 + morr_arr(i,j,k,MORRInd::pgam))) - 6.0 * std::log(morr_arr(i,j,k,MORRInd::lamc))) *
                    (std::exp(m_aimm * (273.15 - morr_arr(i,j,k,MORRInd::t3d))) - 1.0);

                  nnuccc = nnuccc +
                    m_cons40 * std::exp(std::log(morr_arr(i,j,k,MORRInd::cdist1)) + std::log(gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0)) - 3.0 * std::log(morr_arr(i,j,k,MORRInd::lamc))) *
                    (std::exp(m_aimm * (273.15 - morr_arr(i,j,k,MORRInd::t3d))) - 1.0);

                  // PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
                  // MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC
                  nnuccc = std::min(nnuccc, morr_arr(i,j,k,MORRInd::nc3d) / dt);
                }

                // AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
                // FORMULA FROM BEHENG (1994)
                // USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
                // AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
                // AS A GAMMA DISTRIBUTION

                // USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR
                if (morr_arr(i,j,k,MORRInd::qc3d) >= 1.0e-6) {
                  // hm add 12/13/06, replace with newer formula
                  // from khairoutdinov and kogan 2000, mwr
                  prc = 1350.0 * std::pow(morr_arr(i,j,k,MORRInd::qc3d), 2.47) *
                    std::pow((morr_arr(i,j,k,MORRInd::nc3d) / 1.0e6 * morr_arr(i,j,k,MORRInd::rho)), -1.79);

                  // note: nprc1 is change in nr,
                  // nprc is change in nc
                  nprc1 = prc / m_cons29;
                  nprc = prc / (morr_arr(i,j,k,MORRInd::qc3d) / morr_arr(i,j,k,MORRInd::nc3d));

                  // hm bug fix 3/20/12
                  nprc = std::min(nprc, morr_arr(i,j,k,MORRInd::nc3d) / dt);
                  nprc1 = std::min(nprc1, nprc);
                }
                // SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
                // THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW
                if (morr_arr(i,j,k,MORRInd::qni3d) >= 1.0e-8) {
                  nsagg = m_cons15 * morr_arr(i,j,k,MORRInd::asn) * std::pow(morr_arr(i,j,k,MORRInd::rho), ((2.0 + m_bs) / 3.0)) *
                    std::pow(morr_arr(i,j,k,MORRInd::qni3d), ((2.0 + m_bs) / 3.0)) *
                    std::pow((morr_arr(i,j,k,MORRInd::ns3d) * morr_arr(i,j,k,MORRInd::rho)), ((4.0 - m_bs) / 3.0)) / morr_arr(i,j,k,MORRInd::rho);
                }

                // ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
                // HERE USE CONTINUOUS COLLECTION EQUATION WITH
                // SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

                // SNOW
                if (morr_arr(i,j,k,MORRInd::qni3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  psacws = m_cons13 * morr_arr(i,j,k,MORRInd::asn) * morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::rho) *
                    morr_arr(i,j,k,MORRInd::n0s) / std::pow(morr_arr(i,j,k,MORRInd::lams), (m_bs + 3.0));

                  npsacws = m_cons13 * morr_arr(i,j,k,MORRInd::asn) * morr_arr(i,j,k,MORRInd::nc3d) * morr_arr(i,j,k,MORRInd::rho) *
                    morr_arr(i,j,k,MORRInd::n0s) / std::pow(morr_arr(i,j,k,MORRInd::lams), (m_bs + 3.0));
                }

                // COLLECTION OF CLOUD WATER BY GRAUPEL
                if (morr_arr(i,j,k,MORRInd::qg3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  psacwg = m_cons14 * morr_arr(i,j,k,MORRInd::agn) * morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::rho) *
                    morr_arr(i,j,k,MORRInd::n0g) / std::pow(morr_arr(i,j,k,MORRInd::lamg), (m_bg + 3.0));

                  npsacwg = m_cons14 * morr_arr(i,j,k,MORRInd::agn) * morr_arr(i,j,k,MORRInd::nc3d) * morr_arr(i,j,k,MORRInd::rho) *
                    morr_arr(i,j,k,MORRInd::n0g) / std::pow(morr_arr(i,j,k,MORRInd::lamg), (m_bg + 3.0));
                }
                // hm, add 12/13/06
                // CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
                // BEFORE RIMING CAN OCCUR
                // ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
                // TO HALLET-MOSSOP SPLINTERING
                if (morr_arr(i,j,k,MORRInd::qi3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  // PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
                  // FROM THOMPSON ET AL. 2004, MWR
                  if (1.0 / morr_arr(i,j,k,MORRInd::lami) >= 100.0e-6) {
                    psacwi = m_cons16 * morr_arr(i,j,k,MORRInd::ain) * morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::rho) *
                      morr_arr(i,j,k,MORRInd::n0i) / std::pow(morr_arr(i,j,k,MORRInd::lami), (m_bi + 3.0));

                    npsacwi = m_cons16 * morr_arr(i,j,k,MORRInd::ain) * morr_arr(i,j,k,MORRInd::nc3d) * morr_arr(i,j,k,MORRInd::rho) *
                      morr_arr(i,j,k,MORRInd::n0i) / std::pow(morr_arr(i,j,k,MORRInd::lami), (m_bi + 3.0));
                  }
                }

                // ACCRETION OF RAIN WATER BY SNOW
                // FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998
                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qni3d) >= 1.0e-8) {
                  amrex::Real ums_local = morr_arr(i,j,k,MORRInd::asn) * m_cons3 / std::pow(morr_arr(i,j,k,MORRInd::lams), m_bs);
                  amrex::Real umr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons4 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);
                  amrex::Real uns_local = morr_arr(i,j,k,MORRInd::asn) * m_cons5 / std::pow(morr_arr(i,j,k,MORRInd::lams), m_bs);
                  amrex::Real unr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons6 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);

                  // SET REASLISTIC LIMITS ON FALLSPEEDS
                  // bug fix, 10/08/09
                  dum = std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.54);
                  ums_local = std::min(ums_local, 1.2 * dum);
                  uns_local = std::min(uns_local, 1.2 * dum);
                  umr_local = std::min(umr_local, 9.1 * dum);
                  unr_local = std::min(unr_local, 9.1 * dum);

                  pracs = m_cons41 * (std::sqrt(std::pow(1.2 * umr_local - 0.95 * ums_local, 2) +
                                                0.08 * ums_local * umr_local) * morr_arr(i,j,k,MORRInd::rho) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0s) /
                                      std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * (5.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lams)) +
                                                           2.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lams), 2)) +
                                                           0.5 / (morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lams), 3))));

                  npracs = m_cons32 * morr_arr(i,j,k,MORRInd::rho) * std::sqrt(1.7 * std::pow(unr_local - uns_local, 2) +
                                                             0.3 * unr_local * uns_local) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0s) *
                    (1.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lams)) +
                     1.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lams), 2)) +
                     1.0 / (morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lams), 3)));

                  // MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
                  // AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
                  // RIME-SPLINTERING
                  pracs = std::min(pracs, morr_arr(i,j,k,MORRInd::qr3d) / dt);

                  // COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
                  // ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG
                  // hm modify for wrfv3.1
                  if (morr_arr(i,j,k,MORRInd::qni3d) >= 0.1e-3 && morr_arr(i,j,k,MORRInd::qr3d) >= 0.1e-3) {
                    psacr = m_cons31 * (std::sqrt(std::pow(1.2 * umr_local - 0.95 * ums_local, 2) +
                                                  0.08 * ums_local * umr_local) * morr_arr(i,j,k,MORRInd::rho) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0s) /
                                        std::pow(morr_arr(i,j,k,MORRInd::lams), 3) * (5.0 / (std::pow(morr_arr(i,j,k,MORRInd::lams), 3) * morr_arr(i,j,k,MORRInd::lamr)) +
                                                             2.0 / (std::pow(morr_arr(i,j,k,MORRInd::lams), 2) * std::pow(morr_arr(i,j,k,MORRInd::lamr), 2)) +
                                                             0.5 / (morr_arr(i,j,k,MORRInd::lams) * std::pow(morr_arr(i,j,k,MORRInd::lamr), 3))));
                  }
                }

                // COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990,
                // USED BY REISNER ET AL 1998
                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qg3d) >= 1.0e-8) {
                  amrex::Real umg_local = morr_arr(i,j,k,MORRInd::agn) * m_cons7 / std::pow(morr_arr(i,j,k,MORRInd::lamg), m_bg);
                  amrex::Real umr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons4 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);
                  amrex::Real ung_local = morr_arr(i,j,k,MORRInd::agn) * m_cons8 / std::pow(morr_arr(i,j,k,MORRInd::lamg), m_bg);
                  amrex::Real unr_local = morr_arr(i,j,k,MORRInd::arn) * m_cons6 / std::pow(morr_arr(i,j,k,MORRInd::lamr), m_br);

                  // SET REASLISTIC LIMITS ON FALLSPEEDS
                  // bug fix, 10/08/09
                  dum = std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.54);
                  umg_local = std::min(umg_local, 20.0 * dum);
                  ung_local = std::min(ung_local, 20.0 * dum);
                  umr_local = std::min(umr_local, 9.1 * dum);
                  unr_local = std::min(unr_local, 9.1 * dum);

                  pracg = m_cons41 * (std::sqrt(std::pow(1.2 * umr_local - 0.95 * umg_local, 2) +
                                                0.08 * umg_local * umr_local) * morr_arr(i,j,k,MORRInd::rho) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0g) /
                                      std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * (5.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lamg)) +
                                                           2.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 2)) +
                                                           0.5 / (morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 3))));

                  npracg = m_cons32 * morr_arr(i,j,k,MORRInd::rho) * std::sqrt(1.7 * std::pow(unr_local - ung_local, 2) +
                                                             0.3 * unr_local * ung_local) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::n0g) *
                    (1.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::lamg)) +
                     1.0 / (std::pow(morr_arr(i,j,k,MORRInd::lamr), 2) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 2)) +
                     1.0 / (morr_arr(i,j,k,MORRInd::lamr) * std::pow(morr_arr(i,j,k,MORRInd::lamg), 3)));

                  // MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
                  // AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
                  // RIME-SPLINTERING
                  pracg = std::min(pracg, morr_arr(i,j,k,MORRInd::qr3d) / dt);
                }

                // RIME-SPLINTERING - SNOW
                // HALLET-MOSSOP (1974)
                // NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER
                // hm add threshold snow and droplet mixing ratio for rime-splintering
                // to limit rime-splintering in stratiform clouds
                // these thresholds correspond with graupel thresholds in rh 1984
                //v1.4
                if (morr_arr(i,j,k,MORRInd::qni3d) >= 0.1e-3) {
                  if (morr_arr(i,j,k,MORRInd::qc3d) >= 0.5e-3 || morr_arr(i,j,k,MORRInd::qr3d) >= 0.1e-3) {
                    if (psacws > 0.0 || pracs > 0.0) {
                      if (morr_arr(i,j,k,MORRInd::t3d) < 270.16 && morr_arr(i,j,k,MORRInd::t3d) > 265.16) {
                        amrex::Real fmult = 0.0;

                        if (morr_arr(i,j,k,MORRInd::t3d) > 270.16) {
                          fmult = 0.0;
                        } else if (morr_arr(i,j,k,MORRInd::t3d) <= 270.16 && morr_arr(i,j,k,MORRInd::t3d) > 268.16) {
                          fmult = (270.16 - morr_arr(i,j,k,MORRInd::t3d)) / 2.0;
                        } else if (morr_arr(i,j,k,MORRInd::t3d) >= 265.16 && morr_arr(i,j,k,MORRInd::t3d) <= 268.16) {
                          fmult = (morr_arr(i,j,k,MORRInd::t3d) - 265.16) / 3.0;
                        } else if (morr_arr(i,j,k,MORRInd::t3d) < 265.16) {
                          fmult = 0.0;
                        }

                        // 1000 IS TO CONVERT FROM KG TO G
                        // SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW
                        if (psacws > 0.0) {
                          nmults = 35.0e4 * psacws * fmult * 1000.0;
                          qmults = nmults * m_mmult;

                          // CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
                          // THAN WAS RIMED ONTO SNOW
                          qmults = std::min(qmults, psacws);
                          psacws = psacws - qmults;
                        }

                        // RIMING AND SPLINTERING FROM ACCRETED RAINDROPS
                        if (pracs > 0.0) {
                          nmultr = 35.0e4 * pracs * fmult * 1000.0;
                          qmultr = nmultr * m_mmult;

                          // CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
                          // THAN WAS RIMED ONTO SNOW
                          qmultr = std::min(qmultr, pracs);
                          pracs = pracs - qmultr;
                        }
                      }
                    }
                  }
                }

                // RIME-SPLINTERING - GRAUPEL
                // HALLET-MOSSOP (1974)
                // NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER
                // hm add threshold snow mixing ratio for rime-splintering
                // to limit rime-splintering in stratiform clouds
                // v1.4
                if (morr_arr(i,j,k,MORRInd::qg3d) >= 0.1e-3) {
                  if (morr_arr(i,j,k,MORRInd::qc3d) >= 0.5e-3 || morr_arr(i,j,k,MORRInd::qr3d) >= 0.1e-3) {
                    if (psacwg > 0.0 || pracg > 0.0) {
                      if (morr_arr(i,j,k,MORRInd::t3d) < 270.16 && morr_arr(i,j,k,MORRInd::t3d) > 265.16) {
                        amrex::Real fmult = 0.0;

                        if (morr_arr(i,j,k,MORRInd::t3d) > 270.16) {
                          fmult = 0.0;
                        } else if (morr_arr(i,j,k,MORRInd::t3d) <= 270.16 && morr_arr(i,j,k,MORRInd::t3d) > 268.16) {
                          fmult = (270.16 - morr_arr(i,j,k,MORRInd::t3d)) / 2.0;
                        } else if (morr_arr(i,j,k,MORRInd::t3d) >= 265.16 && morr_arr(i,j,k,MORRInd::t3d) <= 268.16) {
                          fmult = (morr_arr(i,j,k,MORRInd::t3d) - 265.16) / 3.0;
                        } else if (morr_arr(i,j,k,MORRInd::t3d) < 265.16) {
                          fmult = 0.0;
                        }

                        // 1000 IS TO CONVERT FROM KG TO G
                        // SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL
                        if (psacwg > 0.0) {
                          nmultg = 35.0e4 * psacwg * fmult * 1000.0;
                          qmultg = nmultg * m_mmult;

                          // CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
                          // THAN WAS RIMED ONTO GRAUPEL
                          qmultg = std::min(qmultg, psacwg);
                          psacwg = psacwg - qmultg;
                        }

                        // RIMING AND SPLINTERING FROM ACCRETED RAINDROPS
                        if (pracg > 0.0) {
                          nmultrg = 35.0e4 * pracg * fmult * 1000.0;
                          qmultrg = nmultrg * m_mmult;

                          // CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
                          // THAN WAS RIMED ONTO GRAUPEL
                          qmultrg = std::min(qmultrg, pracg);
                          pracg = pracg - qmultrg;
                        }
                      }
                    }
                  }
                }

                // CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL
                if (psacws > 0.0) {
                  // ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
                  if (morr_arr(i,j,k,MORRInd::qni3d) >= 0.1e-3 && morr_arr(i,j,k,MORRInd::qc3d) >= 0.5e-3) {
                    // PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
                    pgsacw = std::min(psacws, m_cons17 * dt * morr_arr(i,j,k,MORRInd::n0s) * morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::qc3d) *
                                      morr_arr(i,j,k,MORRInd::asn) * morr_arr(i,j,k,MORRInd::asn) /
                                      (morr_arr(i,j,k,MORRInd::rho) * std::pow(morr_arr(i,j,k,MORRInd::lams), (2.0 * m_bs + 2.0))));

                    // MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
                    dum = std::max(m_rhosn / (m_rhog - m_rhosn) * pgsacw, 0.0);

                    // NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
                    nscng = dum / m_mg0 * morr_arr(i,j,k,MORRInd::rho);
                    // LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
                    nscng = std::min(nscng, morr_arr(i,j,k,MORRInd::ns3d) / dt);

                    // PORTION OF RIMING LEFT FOR SNOW
                    psacws = psacws - pgsacw;
                  }
                }

                // CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL
                if (pracs > 0.0) {
                  // ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
                  if (morr_arr(i,j,k,MORRInd::qni3d) >= 0.1e-3 && morr_arr(i,j,k,MORRInd::qr3d) >= 0.1e-3) {
                    // PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
                    dum = m_cons18 * std::pow(4.0 / morr_arr(i,j,k,MORRInd::lams), 3) * std::pow(4.0 / morr_arr(i,j,k,MORRInd::lams), 3) /
                      (m_cons18 * std::pow(4.0 / morr_arr(i,j,k,MORRInd::lams), 3) * std::pow(4.0 / morr_arr(i,j,k,MORRInd::lams), 3) +
                       m_cons19 * std::pow(4.0 / morr_arr(i,j,k,MORRInd::lamr), 3) * std::pow(4.0 / morr_arr(i,j,k,MORRInd::lamr), 3));
                    dum = std::min(dum, 1.0);
                    dum = std::max(dum, 0.0);

                    pgracs = (1.0 - dum) * pracs;
                    ngracs = (1.0 - dum) * npracs;

                    // LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
                    ngracs = std::min(ngracs, morr_arr(i,j,k,MORRInd::nr3d) / dt);
                    ngracs = std::min(ngracs, morr_arr(i,j,k,MORRInd::ns3d) / dt);

                    // AMOUNT LEFT FOR SNOW PRODUCTION
                    pracs = pracs - pgracs;
                    npracs = npracs - ngracs;

                    // CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
                    psacr = psacr * (1.0 - dum);
                  }
                }

                // FREEZING OF RAIN DROPS
                // FREEZING ALLOWED BELOW -4 C
                if (morr_arr(i,j,k,MORRInd::t3d) < 269.15 && morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  // IMMERSION FREEZING (BIGG 1953)
                  // hm fix 7/15/13 for consistency w/ original formula
                  mnuccr = m_cons20 * morr_arr(i,j,k,MORRInd::nr3d) * (std::exp(m_aimm * (273.15 - morr_arr(i,j,k,MORRInd::t3d))) - 1.0) /
                    std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) / std::pow(morr_arr(i,j,k,MORRInd::lamr), 3);

                  nnuccr = m_pi * morr_arr(i,j,k,MORRInd::nr3d) * m_bimm * (std::exp(m_aimm * (273.15 - morr_arr(i,j,k,MORRInd::t3d))) - 1.0) /
                    std::pow(morr_arr(i,j,k,MORRInd::lamr), 3);

                  // PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
                  nnuccr = std::min(nnuccr, morr_arr(i,j,k,MORRInd::nr3d) / dt);
                }

                // ACCRETION OF CLOUD LIQUID WATER BY RAIN
                // CONTINUOUS COLLECTION EQUATION WITH
                // GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED
                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qc3d) >= 1.0e-8) {
                  // 12/13/06 hm add, replace with newer formula from
                  // khairoutdinov and kogan 2000, mwr
                  dum = morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::qr3d);
                  pra = 67.0 * std::pow(dum, 1.15);
                  npra = pra / (morr_arr(i,j,k,MORRInd::qc3d) / morr_arr(i,j,k,MORRInd::nc3d));
                }

                // SELF-COLLECTION OF RAIN DROPS
                // FROM BEHENG(1994)
                // FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
                // AS DESCRINED ABOVE FOR AUTOCONVERSION
                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8) {
                  // include breakup add 10/09/09
                  dum1 = 300.0e-6;
                  if (1.0 / morr_arr(i,j,k,MORRInd::lamr) < dum1) {
                    dum = 1.0;
                  } else if (1.0 / morr_arr(i,j,k,MORRInd::lamr) >= dum1) {
                    dum = 2.0 - std::exp(2300.0 * (1.0 / morr_arr(i,j,k,MORRInd::lamr) - dum1));
                  }
                  nragg = -5.78 * dum * morr_arr(i,j,k,MORRInd::nr3d) * morr_arr(i,j,k,MORRInd::qr3d) * morr_arr(i,j,k,MORRInd::rho);
                }

                // AUTOCONVERSION OF CLOUD ICE TO SNOW
                // FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
                // HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
                // ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION
                if (morr_arr(i,j,k,MORRInd::qi3d) >= 1.0e-8 && qvqvsi >= 1.0) {
                  nprci = m_cons21 * (morr_arr(i,j,k,MORRInd::qv3d) - qvi) * morr_arr(i,j,k,MORRInd::rho) *
                    morr_arr(i,j,k,MORRInd::n0i) * std::exp(-morr_arr(i,j,k,MORRInd::lami) * m_dcs) * dv / abi;
                  prci = m_cons22 * nprci;
                  nprci = std::min(nprci, morr_arr(i,j,k,MORRInd::ni3d) / dt);
                }

                // ACCRETION OF CLOUD ICE BY SNOW
                // FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
                // AND DS >> DI FOR CONTINUOUS COLLECTION
                if (morr_arr(i,j,k,MORRInd::qni3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                  prai = m_cons23 * morr_arr(i,j,k,MORRInd::asn) * morr_arr(i,j,k,MORRInd::qi3d) * morr_arr(i,j,k,MORRInd::rho) * morr_arr(i,j,k,MORRInd::n0s) /
                    std::pow(morr_arr(i,j,k,MORRInd::lams), (m_bs + 3.0));
                  nprai = m_cons23 * morr_arr(i,j,k,MORRInd::asn) * morr_arr(i,j,k,MORRInd::ni3d) *
                    morr_arr(i,j,k,MORRInd::rho) * morr_arr(i,j,k,MORRInd::n0s) /
                    std::pow(morr_arr(i,j,k,MORRInd::lams), (m_bs + 3.0));
                  nprai = std::min(nprai, morr_arr(i,j,k,MORRInd::ni3d) / dt);
                }

                // hm, add 12/13/06, collision of rain and ice to produce snow or graupel
                // follows reisner et al. 1998
                // assumed fallspeed and size of ice crystal << than for rain
                if (morr_arr(i,j,k,MORRInd::qr3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::qi3d) >= 1.0e-8 && morr_arr(i,j,k,MORRInd::t3d) <= 273.15) {
                  // allow graupel formation from rain-ice collisions only if rain mixing ratio > 0.1 g/kg,
                  // otherwise add to snow
                  if (morr_arr(i,j,k,MORRInd::qr3d) >= 0.1e-3) {
                    niacr = m_cons24 * morr_arr(i,j,k,MORRInd::ni3d) * morr_arr(i,j,k,MORRInd::n0r)* morr_arr(i,j,k,MORRInd::arn) /
                      std::pow(morr_arr(i,j,k,MORRInd::lamr), (m_br + 3.0)) * morr_arr(i,j,k,MORRInd::rho);
                    piacr = m_cons25 * morr_arr(i,j,k,MORRInd::ni3d) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::arn) /
                      std::pow(morr_arr(i,j,k,MORRInd::lamr), (m_br + 3.0)) / std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::rho);
                    praci = m_cons24 * morr_arr(i,j,k,MORRInd::qi3d) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::arn) /
                      std::pow(morr_arr(i,j,k,MORRInd::lamr), (m_br + 3.0)) * morr_arr(i,j,k,MORRInd::rho);
                    niacr = std::min(niacr, morr_arr(i,j,k,MORRInd::nr3d) / dt);
                    niacr = std::min(niacr, morr_arr(i,j,k,MORRInd::ni3d) / dt);
                  } else {
                    niacrs = m_cons24 * morr_arr(i,j,k,MORRInd::ni3d) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::arn) /
                      std::pow(morr_arr(i,j,k,MORRInd::lamr), (m_br + 3.0)) * morr_arr(i,j,k,MORRInd::rho);
                    piacrs = m_cons25 * morr_arr(i,j,k,MORRInd::ni3d) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::arn) /
                      std::pow(morr_arr(i,j,k,MORRInd::lamr), (m_br + 3.0)) / std::pow(morr_arr(i,j,k,MORRInd::lamr), 3) * morr_arr(i,j,k,MORRInd::rho);
                    pracis = m_cons24 * morr_arr(i,j,k,MORRInd::qi3d) * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::arn) /
                      std::pow(morr_arr(i,j,k,MORRInd::lamr), (m_br + 3.0)) * morr_arr(i,j,k,MORRInd::rho);
                    niacrs = std::min(niacrs, morr_arr(i,j,k,MORRInd::nr3d) / dt);
                    niacrs = std::min(niacrs, morr_arr(i,j,k,MORRInd::ni3d) / dt);
                  }
                }

                // NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL
                if (m_inuc == 0) {
                  // ADD THRESHOLD ACCORDING TO GREG THOMSPON
                  if ((qvqvs >= 0.999 && morr_arr(i,j,k,MORRInd::t3d) <= 265.15) || qvqvsi >= 1.08) {
                    // hm, modify dec. 5, 2006, replace with cooper curve
                    kc2 = 0.005 * std::exp(0.304 * (273.15 - morr_arr(i,j,k,MORRInd::t3d))) * 1000.0; // CONVERT FROM L-1 TO M-3
                    // LIMIT TO 500 L-1
                    kc2 = std::min(kc2, 500.0e3);
                    kc2 = std::max(kc2 / morr_arr(i,j,k,MORRInd::rho), 0.0);  // CONVERT TO KG-1

                    if (kc2 > morr_arr(i,j,k,MORRInd::ni3d) + morr_arr(i,j,k,MORRInd::ns3d) + morr_arr(i,j,k,MORRInd::ng3d)) {
                      nnuccd = (kc2 - morr_arr(i,j,k,MORRInd::ni3d) - morr_arr(i,j,k,MORRInd::ns3d) - morr_arr(i,j,k,MORRInd::ng3d)) / dt;
                      mnuccd = nnuccd * m_mi0;
                    }
                  }
                } else if (m_inuc == 1) {
                  if (morr_arr(i,j,k,MORRInd::t3d) < 273.15 && qvqvsi > 1.0) {
                    kc2 = 0.16 * 1000.0 / morr_arr(i,j,k,MORRInd::rho);  // CONVERT FROM L-1 TO KG-1
                    if (kc2 > morr_arr(i,j,k,MORRInd::ni3d) + morr_arr(i,j,k,MORRInd::ns3d) + morr_arr(i,j,k,MORRInd::ng3d)) {
                      nnuccd = (kc2 - morr_arr(i,j,k,MORRInd::ni3d) - morr_arr(i,j,k,MORRInd::ns3d) - morr_arr(i,j,k,MORRInd::ng3d)) / dt;
                      mnuccd = nnuccd * m_mi0;
                    }
                  }
                }

                // CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR
                // NO VENTILATION FOR CLOUD ICE
                epsi = 0.0;
                if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                  epsi = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0i) * morr_arr(i,j,k,MORRInd::rho) * dv / (morr_arr(i,j,k,MORRInd::lami) * morr_arr(i,j,k,MORRInd::lami));
                }

                // VENTILATION FOR SNOW
                epss = 0.0;
                if (morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                  epss = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0s) * morr_arr(i,j,k,MORRInd::rho) * dv *
                    (m_f1s / (morr_arr(i,j,k,MORRInd::lams) * morr_arr(i,j,k,MORRInd::lams)) +
                     m_f2s * std::pow(morr_arr(i,j,k,MORRInd::asn) * morr_arr(i,j,k,MORRInd::rho) / morr_arr(i,j,k,MORRInd::mu), 0.5) *
                     std::pow(sc_schmidt, (1.0 / 3.0)) * m_cons10 /
                     std::pow(morr_arr(i,j,k,MORRInd::lams), m_cons35));
                }

                // Ventilation for graupel
                epsg = 0.0;
                if (morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                  epsg = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0g) * morr_arr(i,j,k,MORRInd::rho) * dv *
                    (m_f1s / (morr_arr(i,j,k,MORRInd::lamg) * morr_arr(i,j,k,MORRInd::lamg)) +
                     m_f2s * std::pow(morr_arr(i,j,k,MORRInd::agn) * morr_arr(i,j,k,MORRInd::rho) / morr_arr(i,j,k,MORRInd::mu), 0.5) *
                     std::pow(sc_schmidt, (1.0 / 3.0)) * m_cons11 /
                     std::pow(morr_arr(i,j,k,MORRInd::lamg), m_cons36));
                }

                // VENTILATION FOR RAIN
                epsr = 0.0;
                if (morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  epsr = 2.0 * m_pi * morr_arr(i,j,k,MORRInd::n0r) * morr_arr(i,j,k,MORRInd::rho) * dv *
                    (m_f1r / (morr_arr(i,j,k,MORRInd::lamr) * morr_arr(i,j,k,MORRInd::lamr)) +
                     m_f2r * std::pow(morr_arr(i,j,k,MORRInd::arn) * morr_arr(i,j,k,MORRInd::rho) / morr_arr(i,j,k,MORRInd::mu), 0.5) *
                     std::pow(sc_schmidt, (1.0 / 3.0)) * m_cons9 /
                     std::pow(morr_arr(i,j,k,MORRInd::lamr), m_cons34));
                }

                // ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
                // DUM IS FRACTION OF D*N(D) < DCS
                // LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
                if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                  dum = (1.0 - std::exp(-morr_arr(i,j,k,MORRInd::lami) * m_dcs) * (1.0 + morr_arr(i,j,k,MORRInd::lami) * m_dcs));
                  prd = epsi * (morr_arr(i,j,k,MORRInd::qv3d) - qvi) / abi * dum;
                } else {
                  dum = 0.0;
                }

                // ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
                if (morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                  prds = epss * (morr_arr(i,j,k,MORRInd::qv3d) - qvi) / abi +
                    epsi * (morr_arr(i,j,k,MORRInd::qv3d) - qvi) / abi * (1.0 - dum);
                } else {
                  // OTHERWISE ADD TO CLOUD ICE
                  prd = prd + epsi * (morr_arr(i,j,k,MORRInd::qv3d) - qvi) / abi * (1.0 - dum);
                }

                // VAPOR DPEOSITION ON GRAUPEL
                prdg = epsg * (morr_arr(i,j,k,MORRInd::qv3d) - qvi) / abi;

                // NO CONDENSATION ONTO RAIN, ONLY EVAP
                if (morr_arr(i,j,k,MORRInd::qv3d) < qvs) {
                  pre = epsr * (morr_arr(i,j,k,MORRInd::qv3d) - qvs) / ab;
                  pre = std::min(pre, 0.0);
                } else {
                  pre = 0.0;
                }

                // MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
                // FORMULA FROM REISNER 2 SCHEME
                dum = (morr_arr(i,j,k,MORRInd::qv3d) - qvi) / dt;

                fudgef = 0.9999;
                sum_dep = prd + prds + mnuccd + prdg;

                if ((dum > 0.0 && sum_dep > dum * fudgef) ||
                    (dum < 0.0 && sum_dep < dum * fudgef)) {
                  mnuccd = fudgef * mnuccd * dum / sum_dep;
                  prd = fudgef * prd * dum / sum_dep;
                  prds = fudgef * prds * dum / sum_dep;
                  prdg = fudgef * prdg * dum / sum_dep;
                }

                // IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES
                if (prd < 0.0) {
                  eprd = prd;
                  prd = 0.0;
                }
                if (prds < 0.0) {
                  eprds = prds;
                  prds = 0.0;
                }
                if (prdg < 0.0) {
                  eprdg = prdg;
                  prdg = 0.0;
                }
                // CONSERVATION OF WATER
                // THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
                // ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.

                // IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
                // THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

                // NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
                // BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
                // FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
                // TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
                // STEP

                // ****SENSITIVITY - NO ICE
                if (m_iliq == 1) {
                  mnuccc = 0.0;
                  nnuccc = 0.0;
                  mnuccr = 0.0;
                  nnuccr = 0.0;
                  mnuccd = 0.0;
                  nnuccd = 0.0;
                }

                // ****SENSITIVITY - NO GRAUPEL
                if (m_igraup == 1) {
                  pracg = 0.0;
                  psacr = 0.0;
                  psacwg = 0.0;
                  prdg = 0.0;
                  eprdg = 0.0;
                  evpmg = 0.0;
                  pgmlt = 0.0;
                  npracg = 0.0;
                  npsacwg = 0.0;
                  nscng = 0.0;
                  ngracs = 0.0;
                  nsubg = 0.0;
                  ngmltg = 0.0;
                  ngmltr = 0.0;

                  // fix 053011
                  piacrs = piacrs + piacr;
                  piacr = 0.0;

                  // fix 070713
                  pracis = pracis + praci;
                  praci = 0.0;
                  psacws = psacws + pgsacw;
                  pgsacw = 0.0;
                  pracs = pracs + pgracs;
                  pgracs = 0.0;
                }

                // CONSERVATION OF QC
                dum = (prc + pra + mnuccc + psacws + psacwi + qmults + psacwg + pgsacw + qmultg) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qc3d) && morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  ratio = morr_arr(i,j,k,MORRInd::qc3d) / dum;

                  prc = prc * ratio;
                  pra = pra * ratio;
                  mnuccc = mnuccc * ratio;
                  psacws = psacws * ratio;
                  psacwi = psacwi * ratio;
                  qmults = qmults * ratio;
                  qmultg = qmultg * ratio;
                  psacwg = psacwg * ratio;
                  pgsacw = pgsacw * ratio;
                }

                // CONSERVATION OF QI
                dum = (-prd - mnuccc + prci + prai - qmults - qmultg - qmultr - qmultrg
                       - mnuccd + praci + pracis - eprd - psacwi) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qi3d) && morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                  ratio = (morr_arr(i,j,k,MORRInd::qi3d) / dt + prd + mnuccc + qmults + qmultg + qmultr + qmultrg +
                                       mnuccd + psacwi) /
                    (prci + prai + praci + pracis - eprd);

                  prci = prci * ratio;
                  prai = prai * ratio;
                  praci = praci * ratio;
                  pracis = pracis * ratio;
                  eprd = eprd * ratio;
                }

                // CONSERVATION OF QR
                dum = ((pracs - pre) + (qmultr + qmultrg - prc) + (mnuccr - pra) +
                       piacr + piacrs + pgracs + pracg) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qr3d) && morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  ratio = (morr_arr(i,j,k,MORRInd::qr3d) / dt + prc + pra) /
                    (-pre + qmultr + qmultrg + pracs + mnuccr + piacr + piacrs + pgracs + pracg);

                  pre = pre * ratio;
                  pracs = pracs * ratio;
                  qmultr = qmultr * ratio;
                  qmultrg = qmultrg * ratio;
                  mnuccr = mnuccr * ratio;
                  piacr = piacr * ratio;
                  piacrs = piacrs * ratio;
                  pgracs = pgracs * ratio;
                  pracg = pracg * ratio;
                }

                // CONSERVATION OF QNI
                if (m_igraup == 0) {
                  dum = (-prds - psacws - prai - prci - pracs - eprds + psacr - piacrs - pracis) * dt;

                  if (dum > morr_arr(i,j,k,MORRInd::qni3d) && morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                    ratio = (morr_arr(i,j,k,MORRInd::qni3d) / dt + prds + psacws + prai + prci + pracs + piacrs + pracis) /
                      (-eprds + psacr);

                    eprds = eprds * ratio;
                    psacr = psacr * ratio;
                  }
                } else if (m_igraup == 1) {
                  // FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
                  dum = (-prds - psacws - prai - prci - pracs - eprds + psacr - piacrs - pracis - mnuccr) * dt;

                  if (dum > morr_arr(i,j,k,MORRInd::qni3d) && morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                    ratio = (morr_arr(i,j,k,MORRInd::qni3d) / dt + prds + psacws + prai + prci + pracs + piacrs + pracis + mnuccr) /
                      (-eprds + psacr);

                    eprds = eprds * ratio;
                    psacr = psacr * ratio;
                  }
                }

                // CONSERVATION OF QG
                dum = (-psacwg - pracg - pgsacw - pgracs - prdg - mnuccr - eprdg - piacr - praci - psacr) * dt;

                if (dum > morr_arr(i,j,k,MORRInd::qg3d) && morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                  ratio = (morr_arr(i,j,k,MORRInd::qg3d) / dt + psacwg + pracg + pgsacw + pgracs + prdg + mnuccr + psacr +
                                       piacr + praci) / (-eprdg);

                  eprdg = eprdg * ratio;
                }

                // TENDENCIES
                morr_arr(i,j,k,MORRInd::qv3dten) = morr_arr(i,j,k,MORRInd::qv3dten) + (-pre - prd - prds - mnuccd - eprd - eprds - prdg - eprdg);

                // bug fix hm, 3/1/11, include piacr and piacrs
                morr_arr(i,j,k,MORRInd::t3dten) = morr_arr(i,j,k,MORRInd::t3dten) + (pre * morr_arr(i,j,k,MORRInd::xxlv) +
                                                 (prd + prds + mnuccd + eprd + eprds + prdg + eprdg) * morr_arr(i,j,k,MORRInd::xxls) +
                                                 (psacws + psacwi + mnuccc + mnuccr + qmults + qmultg + qmultr + qmultrg + pracs +
                                                  psacwg + pracg + pgsacw + pgracs + piacr + piacrs) * morr_arr(i,j,k,MORRInd::xlf)) / morr_arr(i,j,k,MORRInd::cpm);

                morr_arr(i,j,k,MORRInd::qc3dten) = morr_arr(i,j,k,MORRInd::qc3dten) +
                  (-pra - prc - mnuccc + pcc -
                   psacws - psacwi - qmults - qmultg - psacwg - pgsacw);

                morr_arr(i,j,k,MORRInd::qi3dten) = morr_arr(i,j,k,MORRInd::qi3dten) +
                  (prd + eprd + psacwi + mnuccc - prci -
                   prai + qmults + qmultg + qmultr + qmultrg + mnuccd - praci - pracis);

                morr_arr(i,j,k,MORRInd::qr3dten) = morr_arr(i,j,k,MORRInd::qr3dten) +
                  (pre + pra + prc - pracs - mnuccr - qmultr - qmultrg -
                   piacr - piacrs - pracg - pgracs);
                if (m_igraup == 0) {
                  morr_arr(i,j,k,MORRInd::qni3dten) = morr_arr(i,j,k,MORRInd::qni3dten) +
                    (prai + psacws + prds + pracs + prci + eprds - psacr + piacrs + pracis);

                  morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) + (nsagg + nprci - nscng - ngracs + niacrs);

                  morr_arr(i,j,k,MORRInd::qg3dten) = morr_arr(i,j,k,MORRInd::qg3dten) + (pracg + psacwg + pgsacw + pgracs +
                                                     prdg + eprdg + mnuccr + piacr + praci + psacr);

                  morr_arr(i,j,k,MORRInd::ng3dten) = morr_arr(i,j,k,MORRInd::ng3dten) + (nscng + ngracs + nnuccr + niacr);
                } else if (m_igraup == 1) {
                  // FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
                  morr_arr(i,j,k,MORRInd::qni3dten) = morr_arr(i,j,k,MORRInd::qni3dten) +
                    (prai + psacws + prds + pracs + prci + eprds - psacr + piacrs + pracis + mnuccr);

                  morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) + (nsagg + nprci - nscng - ngracs + niacrs + nnuccr);
                }

                morr_arr(i,j,k,MORRInd::nc3dten) = morr_arr(i,j,k,MORRInd::nc3dten) + (-nnuccc - npsacws -
                                                   npra - nprc - npsacwi - npsacwg);

                morr_arr(i,j,k,MORRInd::ni3dten) = morr_arr(i,j,k,MORRInd::ni3dten) +
                  (nnuccc - nprci - nprai + nmults + nmultg + nmultr + nmultrg +
                   nnuccd - niacr - niacrs);

                morr_arr(i,j,k,MORRInd::nr3dten) = morr_arr(i,j,k,MORRInd::nr3dten) + (nprc1 - npracs - nnuccr +
                                                   nragg - niacr - niacrs - npracg - ngracs);

                // hm add, wrf-chem, add tendencies for c2prec
                c2prec = pra + prc + psacws + qmults + qmultg + psacwg +
                  pgsacw + mnuccc + psacwi;

                // CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
                // WATER SATURATION
                dumt = morr_arr(i,j,k,MORRInd::t3d) + dt * morr_arr(i,j,k,MORRInd::t3dten);
                dumqv = morr_arr(i,j,k,MORRInd::qv3d) + dt * morr_arr(i,j,k,MORRInd::qv3dten);

                // hm, add fix for low pressure, 5/12/10
                dum = std::min(0.99 * morr_arr(i,j,k,MORRInd::pres), calc_saturation_vapor_pressure(dumt, 0));
                dumqss = m_ep_2 * dum / (morr_arr(i,j,k,MORRInd::pres) - dum);

                dumqc = morr_arr(i,j,k,MORRInd::qc3d) + dt * morr_arr(i,j,k,MORRInd::qc3dten);
                dumqc = std::max(dumqc, 0.0);

                // SATURATION ADJUSTMENT FOR LIQUID
                dums = dumqv - dumqss;

                pcc = dums / (1.0 + std::pow(morr_arr(i,j,k,MORRInd::xxlv), 2) * dumqss / (morr_arr(i,j,k,MORRInd::cpm) * m_Rv * std::pow(dumt, 2))) / dt;

                if (pcc * dt + dumqc < 0.0) {
                  pcc = -dumqc / dt;
                }

                morr_arr(i,j,k,MORRInd::qv3dten) = morr_arr(i,j,k,MORRInd::qv3dten) - pcc;
                morr_arr(i,j,k,MORRInd::t3dten) = morr_arr(i,j,k,MORRInd::t3dten) + pcc * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm);
                morr_arr(i,j,k,MORRInd::qc3dten) = morr_arr(i,j,k,MORRInd::qc3dten) + pcc;
                // SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
                // THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
                // LOSS OF NUMBER CONCENTRATION
                if (eprd < 0.0) {
                  dum = eprd * dt / morr_arr(i,j,k,MORRInd::qi3d);
                  dum = std::max(-1.0, dum);
                  nsubi = dum * morr_arr(i,j,k,MORRInd::ni3d) / dt;
                }

                if (eprds < 0.0) {
                  dum = eprds * dt / morr_arr(i,j,k,MORRInd::qni3d);
                  dum = std::max(-1.0, dum);
                  nsubs = dum * morr_arr(i,j,k,MORRInd::ns3d) / dt;
                }

                if (pre < 0.0) {
                  dum = pre * dt / morr_arr(i,j,k,MORRInd::qr3d);
                  dum = std::max(-1.0, dum);
                  nsubr = dum * morr_arr(i,j,k,MORRInd::nr3d) / dt;
                }

                if (eprdg < 0.0) {
                  dum = eprdg * dt / morr_arr(i,j,k,MORRInd::qg3d);
                  dum = std::max(-1.0, dum);
                  nsubg = dum * morr_arr(i,j,k,MORRInd::ng3d) / dt;
                }

                // UPDATE TENDENCIES
                morr_arr(i,j,k,MORRInd::ni3dten) = morr_arr(i,j,k,MORRInd::ni3dten) + nsubi;
                morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) + nsubs;
                morr_arr(i,j,k,MORRInd::ng3dten) = morr_arr(i,j,k,MORRInd::ng3dten) + nsubg;
                morr_arr(i,j,k,MORRInd::nr3dten) = morr_arr(i,j,k,MORRInd::nr3dten) + nsubr;
            }
            ltrue = 1;
            }
            //            label_200:

           }
            for(int k=klo; k<=khi; k++) {
            // INITIALIZE PRECIP AND SNOW RATES
            morr_arr(i,j,k,MORRInd::precrt) = 0.0;
            morr_arr(i,j,k,MORRInd::snowrt) = 0.0;
        // hm added 7/13/13
            morr_arr(i,j,k,MORRInd::snowprt) = 0.0;
            morr_arr(i,j,k,MORRInd::grplprt) = 0.0;
            }
            nstep = 1;
            if(ltrue != 0) {
            //goto 400
            // CALCULATE SEDIMENTATION
            // THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
            // FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
            // STABILITY, I.E. COURANT# < 1
            // Loop from top to bottom (KTE to KTS)
            for(int k=khi; k>=klo; k--) {

              amrex::Real dum;                // DUM: General dummy variable

              amrex::Real di0;                // DI0: Characteristic diameter for ice
              amrex::Real ds0;                // DS0: Characteristic diameter for snow
              amrex::Real dg0;                // DG0: Characteristic diameter for graupel
              amrex::Real lammax;             // LAMMAX: Maximum value for slope parameter
              amrex::Real lammin;             // LAMMIN: Minimum value for slope parameter

              ds0 = 3.0;       // Size distribution parameter for snow
              di0 = 3.0;       // Size distribution parameter for cloud ice
              dg0 = 3.0;       // Size distribution parameter for graupel

              // Update prognostic variables with tendencies
              morr_arr(i,j,k,MORRInd::dumi) = morr_arr(i,j,k,MORRInd::qi3d) + morr_arr(i,j,k,MORRInd::qi3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumqs) = morr_arr(i,j,k,MORRInd::qni3d) + morr_arr(i,j,k,MORRInd::qni3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumr) = morr_arr(i,j,k,MORRInd::qr3d) + morr_arr(i,j,k,MORRInd::qr3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumfni) = morr_arr(i,j,k,MORRInd::ni3d) + morr_arr(i,j,k,MORRInd::ni3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumfns) = morr_arr(i,j,k,MORRInd::ns3d) + morr_arr(i,j,k,MORRInd::ns3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumfnr) = morr_arr(i,j,k,MORRInd::nr3d) + morr_arr(i,j,k,MORRInd::nr3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumc) = morr_arr(i,j,k,MORRInd::qc3d) + morr_arr(i,j,k,MORRInd::qc3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumfnc) = morr_arr(i,j,k,MORRInd::nc3d) + morr_arr(i,j,k,MORRInd::nc3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumg) = morr_arr(i,j,k,MORRInd::qg3d) + morr_arr(i,j,k,MORRInd::qg3dten) * dt;
              morr_arr(i,j,k,MORRInd::dumfng) = morr_arr(i,j,k,MORRInd::ng3d) + morr_arr(i,j,k,MORRInd::ng3dten) * dt;

              // SWITCH FOR CONSTANT DROPLET NUMBER
              if (iinum == 1) {
                morr_arr(i,j,k,MORRInd::dumfnc) = morr_arr(i,j,k,MORRInd::nc3d);
              }

              // MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
              morr_arr(i,j,k,MORRInd::dumfni) = amrex::max(0., morr_arr(i,j,k,MORRInd::dumfni));
              morr_arr(i,j,k,MORRInd::dumfns) = amrex::max(0., morr_arr(i,j,k,MORRInd::dumfns));
              morr_arr(i,j,k,MORRInd::dumfnc) = amrex::max(0., morr_arr(i,j,k,MORRInd::dumfnc));
              morr_arr(i,j,k,MORRInd::dumfnr) = amrex::max(0., morr_arr(i,j,k,MORRInd::dumfnr));
              morr_arr(i,j,k,MORRInd::dumfng) = amrex::max(0., morr_arr(i,j,k,MORRInd::dumfng));

              // CLOUD ICE
              if (morr_arr(i,j,k,MORRInd::dumi) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::dlami) = std::pow(m_cons12 * morr_arr(i,j,k,MORRInd::dumfni) / morr_arr(i,j,k,MORRInd::dumi), 1.0/di0);
                morr_arr(i,j,k,MORRInd::dlami) = amrex::max(morr_arr(i,j,k,MORRInd::dlami), m_lammini);
                morr_arr(i,j,k,MORRInd::dlami) = amrex::min(morr_arr(i,j,k,MORRInd::dlami), m_lammaxi);
              }

              // RAIN
              if (morr_arr(i,j,k,MORRInd::dumr) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::dlamr) = std::pow(m_pi * m_rhow * morr_arr(i,j,k,MORRInd::dumfnr) / morr_arr(i,j,k,MORRInd::dumr), 1.0/3.0);
                morr_arr(i,j,k,MORRInd::dlamr) = amrex::max(morr_arr(i,j,k,MORRInd::dlamr), m_lamminr);
                morr_arr(i,j,k,MORRInd::dlamr) = amrex::min(morr_arr(i,j,k,MORRInd::dlamr), m_lammaxr);
              }

              // CLOUD DROPLETS
              if (morr_arr(i,j,k,MORRInd::dumc) >= m_qsmall) {
                dum = morr_arr(i,j,k,MORRInd::pres) / (287.15 * morr_arr(i,j,k,MORRInd::t3d));
                morr_arr(i,j,k,MORRInd::pgam) = 0.0005714 * (morr_arr(i,j,k,MORRInd::nc3d) / 1.0e6 * dum) + 0.2714;
                morr_arr(i,j,k,MORRInd::pgam) = 1.0 / (morr_arr(i,j,k,MORRInd::pgam) * morr_arr(i,j,k,MORRInd::pgam)) - 1.0;
                morr_arr(i,j,k,MORRInd::pgam) = amrex::max(morr_arr(i,j,k,MORRInd::pgam), 2.0);
                morr_arr(i,j,k,MORRInd::pgam) = amrex::min(morr_arr(i,j,k,MORRInd::pgam), 10.0);

                morr_arr(i,j,k,MORRInd::dlamc) = std::pow(m_cons26 * morr_arr(i,j,k,MORRInd::dumfnc) * gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0) /
                                        (morr_arr(i,j,k,MORRInd::dumc) * gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.0)), 1.0/3.0);
                lammin = (morr_arr(i,j,k,MORRInd::pgam) + 1.0) / 60.0e-6;
                lammax = (morr_arr(i,j,k,MORRInd::pgam) + 1.0) / 1.0e-6;
                morr_arr(i,j,k,MORRInd::dlamc) = amrex::max(morr_arr(i,j,k,MORRInd::dlamc), lammin);
                morr_arr(i,j,k,MORRInd::dlamc) = amrex::min(morr_arr(i,j,k,MORRInd::dlamc), lammax);
              }

              // SNOW
              if (morr_arr(i,j,k,MORRInd::dumqs) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::dlams) = std::pow(m_cons1 * morr_arr(i,j,k,MORRInd::dumfns) / morr_arr(i,j,k,MORRInd::dumqs), 1.0/ds0);
                morr_arr(i,j,k,MORRInd::dlams) = amrex::max(morr_arr(i,j,k,MORRInd::dlams), m_lammins);
                morr_arr(i,j,k,MORRInd::dlams) = amrex::min(morr_arr(i,j,k,MORRInd::dlams), m_lammaxs);
              }

              // GRAUPEL
              if (morr_arr(i,j,k,MORRInd::dumg) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::dlamg) = std::pow(m_cons2 * morr_arr(i,j,k,MORRInd::dumfng) / morr_arr(i,j,k,MORRInd::dumg), 1.0/dg0);
                morr_arr(i,j,k,MORRInd::dlamg) = amrex::max(morr_arr(i,j,k,MORRInd::dlamg), m_lamming);
                morr_arr(i,j,k,MORRInd::dlamg) = amrex::min(morr_arr(i,j,k,MORRInd::dlamg), m_lammaxg);
              }

              // Calculate number-weighted and mass-weighted terminal fall speeds
              // CLOUD WATER
              if (morr_arr(i,j,k,MORRInd::dumc) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::unc) = morr_arr(i,j,k,MORRInd::acn) * gamma_function(1. + m_bc + morr_arr(i,j,k,MORRInd::pgam)) /
                             (std::pow(morr_arr(i,j,k,MORRInd::dlamc), m_bc) * gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.));
                morr_arr(i,j,k,MORRInd::umc) = morr_arr(i,j,k,MORRInd::acn) * gamma_function(4. + m_bc + morr_arr(i,j,k,MORRInd::pgam)) /
                             (std::pow(morr_arr(i,j,k,MORRInd::dlamc), m_bc) * gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.));
              } else {
                morr_arr(i,j,k,MORRInd::umc) = 0.;
                morr_arr(i,j,k,MORRInd::unc) = 0.;
              }

              // CLOUD ICE
              if (morr_arr(i,j,k,MORRInd::dumi) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::uni) = morr_arr(i,j,k,MORRInd::ain) * m_cons27 / std::pow(morr_arr(i,j,k,MORRInd::dlami), m_bi);
                morr_arr(i,j,k,MORRInd::umi) = morr_arr(i,j,k,MORRInd::ain) * m_cons28 / std::pow(morr_arr(i,j,k,MORRInd::dlami), m_bi);
              } else {
                morr_arr(i,j,k,MORRInd::umi) = 0.;
                morr_arr(i,j,k,MORRInd::uni) = 0.;
              }

              // RAIN
              if (morr_arr(i,j,k,MORRInd::dumr) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::unr) = morr_arr(i,j,k,MORRInd::arn) * m_cons6 / std::pow(morr_arr(i,j,k,MORRInd::dlamr), m_br);
                morr_arr(i,j,k,MORRInd::umr) = morr_arr(i,j,k,MORRInd::arn) * m_cons4 / std::pow(morr_arr(i,j,k,MORRInd::dlamr), m_br);
              } else {
                morr_arr(i,j,k,MORRInd::umr) = 0.;
                morr_arr(i,j,k,MORRInd::unr) = 0.;
              }

              // SNOW
              if (morr_arr(i,j,k,MORRInd::dumqs) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::ums) = morr_arr(i,j,k,MORRInd::asn) * m_cons3 / std::pow(morr_arr(i,j,k,MORRInd::dlams), m_bs);
                morr_arr(i,j,k,MORRInd::uns) = morr_arr(i,j,k,MORRInd::asn) * m_cons5 / std::pow(morr_arr(i,j,k,MORRInd::dlams), m_bs);
              } else {
                morr_arr(i,j,k,MORRInd::ums) = 0.;
                morr_arr(i,j,k,MORRInd::uns) = 0.;
              }

              // GRAUPEL
              if (morr_arr(i,j,k,MORRInd::dumg) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::umg) = morr_arr(i,j,k,MORRInd::agn) * m_cons7 / std::pow(morr_arr(i,j,k,MORRInd::dlamg), m_bg);
                morr_arr(i,j,k,MORRInd::ung) = morr_arr(i,j,k,MORRInd::agn) * m_cons8 / std::pow(morr_arr(i,j,k,MORRInd::dlamg), m_bg);
              } else {
                morr_arr(i,j,k,MORRInd::umg) = 0.;
                morr_arr(i,j,k,MORRInd::ung) = 0.;
              }

              // SET REALISTIC LIMITS ON FALLSPEED
              // Bug fix, 10/08/09
              dum = std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.54);
              morr_arr(i,j,k,MORRInd::ums) = std::min(morr_arr(i,j,k,MORRInd::ums), 1.2 * dum);
              morr_arr(i,j,k,MORRInd::uns) = std::min(morr_arr(i,j,k,MORRInd::uns), 1.2 * dum);

              // Fix 053011
              // Fix for correction by AA 4/6/11
              morr_arr(i,j,k,MORRInd::umi) = std::min(morr_arr(i,j,k,MORRInd::umi), 1.2 * std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.35));
              morr_arr(i,j,k,MORRInd::uni) = std::min(morr_arr(i,j,k,MORRInd::uni), 1.2 * std::pow(m_rhosu / morr_arr(i,j,k,MORRInd::rho), 0.35));
              morr_arr(i,j,k,MORRInd::umr) = std::min(morr_arr(i,j,k,MORRInd::umr), 9.1 * dum);
              morr_arr(i,j,k,MORRInd::unr) = std::min(morr_arr(i,j,k,MORRInd::unr), 9.1 * dum);
              morr_arr(i,j,k,MORRInd::umg) = std::min(morr_arr(i,j,k,MORRInd::umg), 20. * dum);
              morr_arr(i,j,k,MORRInd::ung) = std::min(morr_arr(i,j,k,MORRInd::ung), 20. * dum);

              // Set fall speed values
              morr_arr(i,j,k,MORRInd::fr) = morr_arr(i,j,k,MORRInd::umr);         // RAIN FALL SPEED
              morr_arr(i,j,k,MORRInd::fi) = morr_arr(i,j,k,MORRInd::umi);         // CLOUD ICE FALL SPEED
              morr_arr(i,j,k,MORRInd::fni) = morr_arr(i,j,k,MORRInd::uni);        // CLOUD ICE NUMBER FALL SPEED
              morr_arr(i,j,k,MORRInd::fs) = morr_arr(i,j,k,MORRInd::ums);         // SNOW FALL SPEED
              morr_arr(i,j,k,MORRInd::fns) = morr_arr(i,j,k,MORRInd::uns);        // SNOW NUMBER FALL SPEED
              morr_arr(i,j,k,MORRInd::fnr) = morr_arr(i,j,k,MORRInd::unr);        // RAIN NUMBER FALL SPEED
              morr_arr(i,j,k,MORRInd::fc) = morr_arr(i,j,k,MORRInd::umc);         // CLOUD WATER FALL SPEED
              morr_arr(i,j,k,MORRInd::fnc) = morr_arr(i,j,k,MORRInd::unc);        // CLOUD NUMBER FALL SPEED
              morr_arr(i,j,k,MORRInd::fg) = morr_arr(i,j,k,MORRInd::umg);         // GRAUPEL FALL SPEED
              morr_arr(i,j,k,MORRInd::fng) = morr_arr(i,j,k,MORRInd::ung);        // GRAUPEL NUMBER FALL SPEED

              // V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP
              if (morr_arr(i,j,k,MORRInd::fr) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fr) = morr_arr(i,j,k+1,MORRInd::fr);
              }
              if (morr_arr(i,j,k,MORRInd::fi) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fi) = morr_arr(i,j,k+1,MORRInd::fi);
              }
              if (morr_arr(i,j,k,MORRInd::fni) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fni) = morr_arr(i,j,k+1,MORRInd::fni);
              }
              if (morr_arr(i,j,k,MORRInd::fs) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fs) = morr_arr(i,j,k+1,MORRInd::fs);
              }
              if (morr_arr(i,j,k,MORRInd::fns) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fns) = morr_arr(i,j,k+1,MORRInd::fns);
              }
              if (morr_arr(i,j,k,MORRInd::fnr) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fnr) = morr_arr(i,j,k+1,MORRInd::fnr);
              }
              if (morr_arr(i,j,k,MORRInd::fc) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fc) = morr_arr(i,j,k+1,MORRInd::fc);
              }
              if (morr_arr(i,j,k,MORRInd::fnc) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fnc) = morr_arr(i,j,k+1,MORRInd::fnc);
              }
              if (morr_arr(i,j,k,MORRInd::fg) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fg) = morr_arr(i,j,k+1,MORRInd::fg);
              }
              if (morr_arr(i,j,k,MORRInd::fng) < 1.e-10) {
                morr_arr(i,j,k,MORRInd::fng) = morr_arr(i,j,k+1,MORRInd::fng);
              }

              // CALCULATE NUMBER OF SPLIT TIME STEPS
              // Find maximum fall speed at this point
              morr_arr(i,j,k,MORRInd::rgvm) = std::max({morr_arr(i,j,k,MORRInd::fr), morr_arr(i,j,k,MORRInd::fi), morr_arr(i,j,k,MORRInd::fs), morr_arr(i,j,k,MORRInd::fc),
                  morr_arr(i,j,k,MORRInd::fni), morr_arr(i,j,k,MORRInd::fnr), morr_arr(i,j,k,MORRInd::fns), morr_arr(i,j,k,MORRInd::fnc),
                  morr_arr(i,j,k,MORRInd::fg), morr_arr(i,j,k,MORRInd::fng)});

              // Calculate number of steps (dt and nstep would need to be defined elsewhere)
              nstep = std::max(static_cast<int>(morr_arr(i,j,k,MORRInd::rgvm) * dt / morr_arr(i,j,k,MORRInd::dzq) + 1.), nstep);
              // MULTIPLY VARIABLES BY RHO
              morr_arr(i,j,k,MORRInd::dumr) = morr_arr(i,j,k,MORRInd::dumr) * morr_arr(i,j,k,MORRInd::rho);       // Rain water content * density
              morr_arr(i,j,k,MORRInd::dumi) = morr_arr(i,j,k,MORRInd::dumi) * morr_arr(i,j,k,MORRInd::rho);       // Cloud ice content * density
              morr_arr(i,j,k,MORRInd::dumfni) = morr_arr(i,j,k,MORRInd::dumfni) * morr_arr(i,j,k,MORRInd::rho);   // Cloud ice number * density
              morr_arr(i,j,k,MORRInd::dumqs) = morr_arr(i,j,k,MORRInd::dumqs) * morr_arr(i,j,k,MORRInd::rho);     // Snow content * density
              morr_arr(i,j,k,MORRInd::dumfns) = morr_arr(i,j,k,MORRInd::dumfns) * morr_arr(i,j,k,MORRInd::rho);   // Snow number * density
              morr_arr(i,j,k,MORRInd::dumfnr) = morr_arr(i,j,k,MORRInd::dumfnr) * morr_arr(i,j,k,MORRInd::rho);   // Rain number * density
              morr_arr(i,j,k,MORRInd::dumc) = morr_arr(i,j,k,MORRInd::dumc) * morr_arr(i,j,k,MORRInd::rho);       // Cloud water content * density
              morr_arr(i,j,k,MORRInd::dumfnc) = morr_arr(i,j,k,MORRInd::dumfnc) * morr_arr(i,j,k,MORRInd::rho);   // Cloud droplet number * density
              morr_arr(i,j,k,MORRInd::dumg) = morr_arr(i,j,k,MORRInd::dumg) * morr_arr(i,j,k,MORRInd::rho);       // Graupel content * density
              morr_arr(i,j,k,MORRInd::dumfng) = morr_arr(i,j,k,MORRInd::dumfng) * morr_arr(i,j,k,MORRInd::rho);   // Graupel number * density
            }
            // Main time stepping loop for sedimentation
            for (int n = 1; n <= nstep; n++) {
              // Calculate initial fallout for each hydrometeor type for all levels
              for (int k = klo; k <= khi; k++) {
                morr_arr(i,j,k,MORRInd::faloutr) = morr_arr(i,j,k,MORRInd::fr) * morr_arr(i,j,k,MORRInd::dumr);
                morr_arr(i,j,k,MORRInd::falouti) = morr_arr(i,j,k,MORRInd::fi) * morr_arr(i,j,k,MORRInd::dumi);
                morr_arr(i,j,k,MORRInd::faloutni) = morr_arr(i,j,k,MORRInd::fni) * morr_arr(i,j,k,MORRInd::dumfni);
                morr_arr(i,j,k,MORRInd::falouts) = morr_arr(i,j,k,MORRInd::fs) * morr_arr(i,j,k,MORRInd::dumqs);
                morr_arr(i,j,k,MORRInd::faloutns) = morr_arr(i,j,k,MORRInd::fns) * morr_arr(i,j,k,MORRInd::dumfns);
                morr_arr(i,j,k,MORRInd::faloutnr) = morr_arr(i,j,k,MORRInd::fnr) * morr_arr(i,j,k,MORRInd::dumfnr);
                morr_arr(i,j,k,MORRInd::faloutc) = morr_arr(i,j,k,MORRInd::fc) * morr_arr(i,j,k,MORRInd::dumc);
                morr_arr(i,j,k,MORRInd::faloutnc) = morr_arr(i,j,k,MORRInd::fnc) * morr_arr(i,j,k,MORRInd::dumfnc);
                morr_arr(i,j,k,MORRInd::faloutg) = morr_arr(i,j,k,MORRInd::fg) * morr_arr(i,j,k,MORRInd::dumg);
                morr_arr(i,j,k,MORRInd::faloutng) = morr_arr(i,j,k,MORRInd::fng) * morr_arr(i,j,k,MORRInd::dumfng);
              }

              // Process top of model level
              int k = khi;

              // Calculate tendencies at top level
              morr_arr(i,j,k,MORRInd::faltndr) = morr_arr(i,j,k,MORRInd::faloutr) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndi) = morr_arr(i,j,k,MORRInd::falouti) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndni) = morr_arr(i,j,k,MORRInd::faloutni) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltnds) = morr_arr(i,j,k,MORRInd::falouts) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndns) = morr_arr(i,j,k,MORRInd::faloutns) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndnr) = morr_arr(i,j,k,MORRInd::faloutnr) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndc) = morr_arr(i,j,k,MORRInd::faloutc) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndnc) = morr_arr(i,j,k,MORRInd::faloutnc) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndg) = morr_arr(i,j,k,MORRInd::faloutg) / morr_arr(i,j,k,MORRInd::dzq);
              morr_arr(i,j,k,MORRInd::faltndng) = morr_arr(i,j,k,MORRInd::faloutng) / morr_arr(i,j,k,MORRInd::dzq);

              // Add fallout terms to Eulerian tendencies (scaled by time step and density)
              morr_arr(i,j,k,MORRInd::qrsten) = morr_arr(i,j,k,MORRInd::qrsten) - morr_arr(i,j,k,MORRInd::faltndr) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::qisten) = morr_arr(i,j,k,MORRInd::qisten) - morr_arr(i,j,k,MORRInd::faltndi) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::ni3dten) = morr_arr(i,j,k,MORRInd::ni3dten) - morr_arr(i,j,k,MORRInd::faltndni) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::qnisten) = morr_arr(i,j,k,MORRInd::qnisten) - morr_arr(i,j,k,MORRInd::faltnds) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) - morr_arr(i,j,k,MORRInd::faltndns) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::nr3dten) = morr_arr(i,j,k,MORRInd::nr3dten) - morr_arr(i,j,k,MORRInd::faltndnr) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::qcsten) = morr_arr(i,j,k,MORRInd::qcsten) - morr_arr(i,j,k,MORRInd::faltndc) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::nc3dten) = morr_arr(i,j,k,MORRInd::nc3dten) - morr_arr(i,j,k,MORRInd::faltndnc) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::qgsten) = morr_arr(i,j,k,MORRInd::qgsten) - morr_arr(i,j,k,MORRInd::faltndg) / nstep / morr_arr(i,j,k,MORRInd::rho);
              morr_arr(i,j,k,MORRInd::ng3dten) = morr_arr(i,j,k,MORRInd::ng3dten) - morr_arr(i,j,k,MORRInd::faltndng) / nstep / morr_arr(i,j,k,MORRInd::rho);

              // Update temporary working variables
              morr_arr(i,j,k,MORRInd::dumr) = morr_arr(i,j,k,MORRInd::dumr) - morr_arr(i,j,k,MORRInd::faltndr) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumi) = morr_arr(i,j,k,MORRInd::dumi) - morr_arr(i,j,k,MORRInd::faltndi) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumfni) = morr_arr(i,j,k,MORRInd::dumfni) - morr_arr(i,j,k,MORRInd::faltndni) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumqs) = morr_arr(i,j,k,MORRInd::dumqs) - morr_arr(i,j,k,MORRInd::faltnds) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumfns) = morr_arr(i,j,k,MORRInd::dumfns) - morr_arr(i,j,k,MORRInd::faltndns) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumfnr) = morr_arr(i,j,k,MORRInd::dumfnr) - morr_arr(i,j,k,MORRInd::faltndnr) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumc) = morr_arr(i,j,k,MORRInd::dumc) - morr_arr(i,j,k,MORRInd::faltndc) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumfnc) = morr_arr(i,j,k,MORRInd::dumfnc) - morr_arr(i,j,k,MORRInd::faltndnc) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumg) = morr_arr(i,j,k,MORRInd::dumg) - morr_arr(i,j,k,MORRInd::faltndg) * dt / nstep;
              morr_arr(i,j,k,MORRInd::dumfng) = morr_arr(i,j,k,MORRInd::dumfng) - morr_arr(i,j,k,MORRInd::faltndng) * dt / nstep;

              // Process remaining levels from top to bottom
              for (k = khi-1; k >= klo; k--) {
                // Calculate tendencies based on difference between levels
                morr_arr(i,j,k,MORRInd::faltndr) = (morr_arr(i,j,k+1,MORRInd::faloutr) - morr_arr(i,j,k,MORRInd::faloutr)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndi) = (morr_arr(i,j,k+1,MORRInd::falouti) - morr_arr(i,j,k,MORRInd::falouti)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndni) = (morr_arr(i,j,k+1,MORRInd::faloutni) - morr_arr(i,j,k,MORRInd::faloutni)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltnds) = (morr_arr(i,j,k+1,MORRInd::falouts) - morr_arr(i,j,k,MORRInd::falouts)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndns) = (morr_arr(i,j,k+1,MORRInd::faloutns) - morr_arr(i,j,k,MORRInd::faloutns)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndnr) = (morr_arr(i,j,k+1,MORRInd::faloutnr) - morr_arr(i,j,k,MORRInd::faloutnr)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndc) = (morr_arr(i,j,k+1,MORRInd::faloutc) - morr_arr(i,j,k,MORRInd::faloutc)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndnc) = (morr_arr(i,j,k+1,MORRInd::faloutnc) - morr_arr(i,j,k,MORRInd::faloutnc)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndg) = (morr_arr(i,j,k+1,MORRInd::faloutg) - morr_arr(i,j,k,MORRInd::faloutg)) / morr_arr(i,j,k,MORRInd::dzq);
                morr_arr(i,j,k,MORRInd::faltndng) = (morr_arr(i,j,k+1,MORRInd::faloutng) - morr_arr(i,j,k,MORRInd::faloutng)) / morr_arr(i,j,k,MORRInd::dzq);

                // Add fallout terms to Eulerian tendencies (positive here, as mass flows in from above)
                morr_arr(i,j,k,MORRInd::qrsten) = morr_arr(i,j,k,MORRInd::qrsten) + morr_arr(i,j,k,MORRInd::faltndr) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::qisten) = morr_arr(i,j,k,MORRInd::qisten) + morr_arr(i,j,k,MORRInd::faltndi) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::ni3dten) = morr_arr(i,j,k,MORRInd::ni3dten) + morr_arr(i,j,k,MORRInd::faltndni) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::qnisten) = morr_arr(i,j,k,MORRInd::qnisten) + morr_arr(i,j,k,MORRInd::faltnds) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) + morr_arr(i,j,k,MORRInd::faltndns) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::nr3dten) = morr_arr(i,j,k,MORRInd::nr3dten) + morr_arr(i,j,k,MORRInd::faltndnr) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::qcsten) = morr_arr(i,j,k,MORRInd::qcsten) + morr_arr(i,j,k,MORRInd::faltndc) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::nc3dten) = morr_arr(i,j,k,MORRInd::nc3dten) + morr_arr(i,j,k,MORRInd::faltndnc) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::qgsten) = morr_arr(i,j,k,MORRInd::qgsten) + morr_arr(i,j,k,MORRInd::faltndg) / nstep / morr_arr(i,j,k,MORRInd::rho);
                morr_arr(i,j,k,MORRInd::ng3dten) = morr_arr(i,j,k,MORRInd::ng3dten) + morr_arr(i,j,k,MORRInd::faltndng) / nstep / morr_arr(i,j,k,MORRInd::rho);
                // Update temporary working variables
                morr_arr(i,j,k,MORRInd::dumr) = morr_arr(i,j,k,MORRInd::dumr) + morr_arr(i,j,k,MORRInd::faltndr) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumi) = morr_arr(i,j,k,MORRInd::dumi) + morr_arr(i,j,k,MORRInd::faltndi) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumfni) = morr_arr(i,j,k,MORRInd::dumfni) + morr_arr(i,j,k,MORRInd::faltndni) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumqs) = morr_arr(i,j,k,MORRInd::dumqs) + morr_arr(i,j,k,MORRInd::faltnds) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumfns) = morr_arr(i,j,k,MORRInd::dumfns) + morr_arr(i,j,k,MORRInd::faltndns) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumfnr) = morr_arr(i,j,k,MORRInd::dumfnr) + morr_arr(i,j,k,MORRInd::faltndnr) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumc) = morr_arr(i,j,k,MORRInd::dumc) + morr_arr(i,j,k,MORRInd::faltndc) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumfnc) = morr_arr(i,j,k,MORRInd::dumfnc) + morr_arr(i,j,k,MORRInd::faltndnc) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumg) = morr_arr(i,j,k,MORRInd::dumg) + morr_arr(i,j,k,MORRInd::faltndg) * dt / nstep;
                morr_arr(i,j,k,MORRInd::dumfng) = morr_arr(i,j,k,MORRInd::dumfng) + morr_arr(i,j,k,MORRInd::faltndng) * dt / nstep;
              }
              // Get precipitation and snowfall accumulation during the time step
              // Factor of 1000 converts from m to mm, but division by density
              // of liquid water cancels this factor of 1000
              int kts=klo;
              morr_arr(i,j,klo,MORRInd::precrt) += (morr_arr(i,j,kts,MORRInd::faloutr) + morr_arr(i,j,kts,MORRInd::faloutc) + morr_arr(i,j,kts,MORRInd::falouts) +
                         morr_arr(i,j,kts,MORRInd::falouti) + morr_arr(i,j,kts,MORRInd::faloutg)) * dt / nstep;
              morr_arr(i,j,klo,MORRInd::snowrt) += (morr_arr(i,j,kts,MORRInd::falouts) + morr_arr(i,j,kts,MORRInd::falouti) + morr_arr(i,j,kts,MORRInd::faloutg)) * dt / nstep;

              // Added 7/13/13
              morr_arr(i,j,klo,MORRInd::snowprt) += (morr_arr(i,j,kts,MORRInd::falouti) + morr_arr(i,j,kts,MORRInd::falouts)) * dt / nstep;
              morr_arr(i,j,klo,MORRInd::grplprt) += morr_arr(i,j,kts,MORRInd::faloutg) * dt / nstep;
            }
            for(int k=klo; k<=khi; k++) {
              amrex::Real evs;                // EVS: Saturation vapor pressure
              amrex::Real eis;                // EIS: Ice saturation vapor pressure
              amrex::Real qvs;                // QVS: Saturation mixing ratio
              amrex::Real qvi;                // QVI: Ice saturation mixing ratio
              amrex::Real qvqvs;              // QVQVS: Saturation ratio
              amrex::Real qvqvsi;             // QVQVSI: Ice saturation ratio

              // ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES
              morr_arr(i,j,k,MORRInd::qr3dten) = morr_arr(i,j,k,MORRInd::qr3dten) + morr_arr(i,j,k,MORRInd::qrsten);
              morr_arr(i,j,k,MORRInd::qi3dten) = morr_arr(i,j,k,MORRInd::qi3dten) + morr_arr(i,j,k,MORRInd::qisten);
              morr_arr(i,j,k,MORRInd::qc3dten) = morr_arr(i,j,k,MORRInd::qc3dten) + morr_arr(i,j,k,MORRInd::qcsten);
              morr_arr(i,j,k,MORRInd::qg3dten) = morr_arr(i,j,k,MORRInd::qg3dten) + morr_arr(i,j,k,MORRInd::qgsten);
              morr_arr(i,j,k,MORRInd::qni3dten) = morr_arr(i,j,k,MORRInd::qni3dten) + morr_arr(i,j,k,MORRInd::qnisten);
              // PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs
              // bug fix
              if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall && morr_arr(i,j,k,MORRInd::t3d) < 273.15 && morr_arr(i,j,k,MORRInd::lami) >= 1.e-10) {
                if (1.0/morr_arr(i,j,k,MORRInd::lami) >= 2.0*m_dcs) {
                  morr_arr(i,j,k,MORRInd::qni3dten) = morr_arr(i,j,k,MORRInd::qni3dten) + morr_arr(i,j,k,MORRInd::qi3d)/dt + morr_arr(i,j,k,MORRInd::qi3dten);
                  morr_arr(i,j,k,MORRInd::ns3dten) = morr_arr(i,j,k,MORRInd::ns3dten) + morr_arr(i,j,k,MORRInd::ni3d)/dt + morr_arr(i,j,k,MORRInd::ni3dten);
                  morr_arr(i,j,k,MORRInd::qi3dten) = -morr_arr(i,j,k,MORRInd::qi3d)/dt;
                  morr_arr(i,j,k,MORRInd::ni3dten) = -morr_arr(i,j,k,MORRInd::ni3d)/dt;
                }
              }

              // Add tendencies to ensure consistency between mixing ratio and number concentration
              morr_arr(i,j,k,MORRInd::qc3d) = morr_arr(i,j,k,MORRInd::qc3d) + morr_arr(i,j,k,MORRInd::qc3dten)*dt;
              morr_arr(i,j,k,MORRInd::qi3d) = morr_arr(i,j,k,MORRInd::qi3d) + morr_arr(i,j,k,MORRInd::qi3dten)*dt;
              morr_arr(i,j,k,MORRInd::qni3d) = morr_arr(i,j,k,MORRInd::qni3d) + morr_arr(i,j,k,MORRInd::qni3dten)*dt;
              morr_arr(i,j,k,MORRInd::qr3d) = morr_arr(i,j,k,MORRInd::qr3d) + morr_arr(i,j,k,MORRInd::qr3dten)*dt;
              morr_arr(i,j,k,MORRInd::nc3d) = morr_arr(i,j,k,MORRInd::nc3d) + morr_arr(i,j,k,MORRInd::nc3dten)*dt;
              morr_arr(i,j,k,MORRInd::ni3d) = morr_arr(i,j,k,MORRInd::ni3d) + morr_arr(i,j,k,MORRInd::ni3dten)*dt;
              morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::ns3d) + morr_arr(i,j,k,MORRInd::ns3dten)*dt;
              morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::nr3d) + morr_arr(i,j,k,MORRInd::nr3dten)*dt;
              if (m_igraup == 0) {
                morr_arr(i,j,k,MORRInd::qg3d) = morr_arr(i,j,k,MORRInd::qg3d) + morr_arr(i,j,k,MORRInd::qg3dten)*dt;
                morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::ng3d) + morr_arr(i,j,k,MORRInd::ng3dten)*dt;
              }

              // ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
              morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) + morr_arr(i,j,k,MORRInd::t3dten)*dt;
              morr_arr(i,j,k,MORRInd::qv3d) = morr_arr(i,j,k,MORRInd::qv3d) + morr_arr(i,j,k,MORRInd::qv3dten)*dt;
              // SATURATION VAPOR PRESSURE AND MIXING RATIO
              // hm, add fix for low pressure, 5/12/10
              // Assuming POLYSVP is defined elsewhere
              evs = std::min(0.99 * morr_arr(i,j,k,MORRInd::pres), calc_saturation_vapor_pressure(morr_arr(i,j,k,MORRInd::t3d), 0));  // PA
              eis = std::min(0.99 * morr_arr(i,j,k,MORRInd::pres), calc_saturation_vapor_pressure(morr_arr(i,j,k,MORRInd::t3d), 1));  // PA

              // MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING
              if (eis > evs) {
                eis = evs; // temporary update: adjust ice saturation pressure
              }

              // SATURATION MIXING RATIOS
              qvs = m_ep_2 * evs / (morr_arr(i,j,k,MORRInd::pres) - evs); // budget equation: calculate water saturation mixing ratio
              qvi = m_ep_2 * eis / (morr_arr(i,j,k,MORRInd::pres) - eis); // budget equation: calculate ice saturation mixing ratio

              // SATURATION RATIOS
              qvqvs = morr_arr(i,j,k,MORRInd::qv3d) / qvs; // budget equation: calculate water saturation ratio
              qvqvsi = morr_arr(i,j,k,MORRInd::qv3d) / qvi; // budget equation: calculate ice saturation ratio
              // AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
              if (qvqvs < 0.9) {
                if (morr_arr(i,j,k,MORRInd::qr3d) < 1.0e-8) {
                  morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qr3d);
                  morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qr3d) * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm);
                  morr_arr(i,j,k,MORRInd::qr3d) = 0.0;
                }
                if (morr_arr(i,j,k,MORRInd::qc3d) < 1.0e-8) {
                  morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qc3d);
                  morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::xxlv) / morr_arr(i,j,k,MORRInd::cpm);
                  morr_arr(i,j,k,MORRInd::qc3d) = 0.0;
                }
              }
              if (qvqvsi < 0.9) {
                if (morr_arr(i,j,k,MORRInd::qi3d) < 1.0e-8) {
                  morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qi3d);
                  morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qi3d) * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm);
                  morr_arr(i,j,k,MORRInd::qi3d) = 0.0;
                }
                if (morr_arr(i,j,k,MORRInd::qni3d) < 1.0e-8) {
                  morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qni3d);
                  morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qni3d) * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm);
                  morr_arr(i,j,k,MORRInd::qni3d) = 0.0;
                }
                if (morr_arr(i,j,k,MORRInd::qg3d) < 1.0e-8) {
                  morr_arr(i,j,k,MORRInd::qv3d) += morr_arr(i,j,k,MORRInd::qg3d);
                  morr_arr(i,j,k,MORRInd::t3d) -= morr_arr(i,j,k,MORRInd::qg3d) * morr_arr(i,j,k,MORRInd::xxls) / morr_arr(i,j,k,MORRInd::cpm);
                  morr_arr(i,j,k,MORRInd::qg3d) = 0.0;
                }
              }
              // IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO
              if (morr_arr(i,j,k,MORRInd::qc3d) < m_qsmall) {
                morr_arr(i,j,k,MORRInd::qc3d) = 0.0;
                morr_arr(i,j,k,MORRInd::nc3d) = 0.0;
                morr_arr(i,j,k,MORRInd::effc) = 0.0;
              }
              if (morr_arr(i,j,k,MORRInd::qr3d) < m_qsmall) {
                morr_arr(i,j,k,MORRInd::qr3d) = 0.0;
                morr_arr(i,j,k,MORRInd::nr3d) = 0.0;
                morr_arr(i,j,k,MORRInd::effr) = 0.0;
              }
              if (morr_arr(i,j,k,MORRInd::qi3d) < m_qsmall) {
                morr_arr(i,j,k,MORRInd::qi3d) = 0.0;
                morr_arr(i,j,k,MORRInd::ni3d) = 0.0;
                morr_arr(i,j,k,MORRInd::effi) = 0.0;
              }
              if (morr_arr(i,j,k,MORRInd::qni3d) < m_qsmall) {
                morr_arr(i,j,k,MORRInd::qni3d) = 0.0;
                morr_arr(i,j,k,MORRInd::ns3d) = 0.0;
                morr_arr(i,j,k,MORRInd::effs) = 0.0;
              }
              if (morr_arr(i,j,k,MORRInd::qg3d) < m_qsmall) {
                morr_arr(i,j,k,MORRInd::qg3d) = 0.0;
                morr_arr(i,j,k,MORRInd::ng3d) = 0.0;
                morr_arr(i,j,k,MORRInd::effg) = 0.0;
              }
              /*
              // Skip calculations if there is no cloud/precipitation water
              if ((morr_arr(i,j,k,MORRInd::qc3d) < m_qsmall &&    // CLOUD WATER MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qi3d) < m_qsmall &&    // CLOUD ICE MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qni3d) < m_qsmall &&   // SNOW MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qr3d) < m_qsmall &&    // RAIN MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qg3d) < m_qsmall)) {    // GRAUPEL MIX RATIO (KG/KG)
                goto label_500;
              } else {*/
              if (!(morr_arr(i,j,k,MORRInd::qc3d) < m_qsmall &&    // CLOUD WATER MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qi3d) < m_qsmall &&    // CLOUD ICE MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qni3d) < m_qsmall &&   // SNOW MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qr3d) < m_qsmall &&    // RAIN MIXING RATIO (KG/KG)
                    morr_arr(i,j,k,MORRInd::qg3d) < m_qsmall)) {    // GRAUPEL MIX RATIO (KG/KG)
              // CALCULATE INSTANTANEOUS PROCESSES

              // ADD MELTING OF CLOUD ICE TO FORM RAIN
              if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall && morr_arr(i,j,k,MORRInd::t3d) >= 273.15) {
                morr_arr(i,j,k,MORRInd::qr3d) = morr_arr(i,j,k,MORRInd::qr3d) + morr_arr(i,j,k,MORRInd::qi3d);
                morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) - morr_arr(i,j,k,MORRInd::qi3d) * morr_arr(i,j,k,MORRInd::xlf) / morr_arr(i,j,k,MORRInd::cpm);
                morr_arr(i,j,k,MORRInd::qi3d) = 0.0;
                morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::nr3d) + morr_arr(i,j,k,MORRInd::ni3d);
                morr_arr(i,j,k,MORRInd::ni3d) = 0.0;
              }
              // ****SENSITIVITY - NO ICE
              if ((m_iliq != 1)) {

                // HOMOGENEOUS FREEZING OF CLOUD WATER
                if (morr_arr(i,j,k,MORRInd::t3d) <= 233.15 && morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  morr_arr(i,j,k,MORRInd::qi3d) = morr_arr(i,j,k,MORRInd::qi3d) + morr_arr(i,j,k,MORRInd::qc3d);
                  morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) + morr_arr(i,j,k,MORRInd::qc3d) * morr_arr(i,j,k,MORRInd::xlf) / morr_arr(i,j,k,MORRInd::cpm);
                  morr_arr(i,j,k,MORRInd::qc3d) = 0.0;
                  morr_arr(i,j,k,MORRInd::ni3d) = morr_arr(i,j,k,MORRInd::ni3d) + morr_arr(i,j,k,MORRInd::nc3d);
                  morr_arr(i,j,k,MORRInd::nc3d) = 0.0;
                }
                // HOMOGENEOUS FREEZING OF RAIN
                if (m_igraup == 0) {
                  if (morr_arr(i,j,k,MORRInd::t3d) <= 233.15 && morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                    morr_arr(i,j,k,MORRInd::qg3d) = morr_arr(i,j,k,MORRInd::qg3d) + morr_arr(i,j,k,MORRInd::qr3d);
                    morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) + morr_arr(i,j,k,MORRInd::qr3d) * morr_arr(i,j,k,MORRInd::xlf) / morr_arr(i,j,k,MORRInd::cpm);
                    morr_arr(i,j,k,MORRInd::qr3d) = 0.0;
                    morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::ng3d) + morr_arr(i,j,k,MORRInd::nr3d);
                    morr_arr(i,j,k,MORRInd::nr3d) = 0.0;
                  }
                } else if (m_igraup == 1) {
                  if (morr_arr(i,j,k,MORRInd::t3d) <= 233.15 && morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                    morr_arr(i,j,k,MORRInd::qni3d) = morr_arr(i,j,k,MORRInd::qni3d) + morr_arr(i,j,k,MORRInd::qr3d);
                    morr_arr(i,j,k,MORRInd::t3d) = morr_arr(i,j,k,MORRInd::t3d) + morr_arr(i,j,k,MORRInd::qr3d) * morr_arr(i,j,k,MORRInd::xlf) / morr_arr(i,j,k,MORRInd::cpm);
                    morr_arr(i,j,k,MORRInd::qr3d) = 0.0;
                    morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::ns3d) + morr_arr(i,j,k,MORRInd::nr3d);
                    morr_arr(i,j,k,MORRInd::nr3d) = 0.0;
                  }
                }

              }/* else {
                Real dontdoanything=m_iliq;//printf("m_iliq: %d\n",m_iliq);//                goto label_778;
              }*/

//            label_778:
                // MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE
                morr_arr(i,j,k,MORRInd::ni3d) = std::max(0.0, morr_arr(i,j,k,MORRInd::ni3d));
                morr_arr(i,j,k,MORRInd::ns3d) = std::max(0.0, morr_arr(i,j,k,MORRInd::ns3d));
                morr_arr(i,j,k,MORRInd::nc3d) = std::max(0.0, morr_arr(i,j,k,MORRInd::nc3d));
                morr_arr(i,j,k,MORRInd::nr3d) = std::max(0.0, morr_arr(i,j,k,MORRInd::nr3d));
                morr_arr(i,j,k,MORRInd::ng3d) = std::max(0.0, morr_arr(i,j,k,MORRInd::ng3d));

                // CLOUD ICE
                if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                  morr_arr(i,j,k,MORRInd::lami) = std::pow(m_cons12 * morr_arr(i,j,k,MORRInd::ni3d) / morr_arr(i,j,k,MORRInd::qi3d), 1.0/m_di);
                  // CHECK FOR SLOPE
                  // ADJUST VARS
                  if (morr_arr(i,j,k,MORRInd::lami) < m_lammini) {
                    morr_arr(i,j,k,MORRInd::lami) = m_lammini;
                    morr_arr(i,j,k,MORRInd::n0i) = std::pow(morr_arr(i,j,k,MORRInd::lami), 4) * morr_arr(i,j,k,MORRInd::qi3d) / m_cons12;
                    morr_arr(i,j,k,MORRInd::ni3d) = morr_arr(i,j,k,MORRInd::n0i) / morr_arr(i,j,k,MORRInd::lami);
                  } else if (morr_arr(i,j,k,MORRInd::lami) > m_lammaxi) {
                    morr_arr(i,j,k,MORRInd::lami) = m_lammaxi;
                    morr_arr(i,j,k,MORRInd::n0i) = std::pow(morr_arr(i,j,k,MORRInd::lami), 4) * morr_arr(i,j,k,MORRInd::qi3d) / m_cons12;
                    morr_arr(i,j,k,MORRInd::ni3d) = morr_arr(i,j,k,MORRInd::n0i) / morr_arr(i,j,k,MORRInd::lami);
                  }
                }

                // RAIN
                if (morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                  morr_arr(i,j,k,MORRInd::lamr) = std::pow(m_pi * m_rhow * morr_arr(i,j,k,MORRInd::nr3d) / morr_arr(i,j,k,MORRInd::qr3d), 1.0/3.0);

                  // CHECK FOR SLOPE
                  // ADJUST VARS
                  if (morr_arr(i,j,k,MORRInd::lamr) < m_lamminr) {
                    morr_arr(i,j,k,MORRInd::lamr) = m_lamminr;
                    morr_arr(i,j,k,MORRInd::n0r) = std::pow(morr_arr(i,j,k,MORRInd::lamr), 4) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                    morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::n0r) / morr_arr(i,j,k,MORRInd::lamr);
                  } else if (morr_arr(i,j,k,MORRInd::lamr) > m_lammaxr) {
                    morr_arr(i,j,k,MORRInd::lamr) = m_lammaxr;
                    morr_arr(i,j,k,MORRInd::n0r) = std::pow(morr_arr(i,j,k,MORRInd::lamr), 4) * morr_arr(i,j,k,MORRInd::qr3d) / (m_pi * m_rhow);
                    morr_arr(i,j,k,MORRInd::nr3d) = morr_arr(i,j,k,MORRInd::n0r) / morr_arr(i,j,k,MORRInd::lamr);
                  }
                }

                // CLOUD DROPLETS
                // MARTIN ET AL. (1994) FORMULA FOR PGAM
                if (morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                  amrex::Real dum = morr_arr(i,j,k,MORRInd::pres) / (287.15 * morr_arr(i,j,k,MORRInd::t3d));
                  morr_arr(i,j,k,MORRInd::pgam) = 0.0005714 * (morr_arr(i,j,k,MORRInd::nc3d) / 1.0e6 * dum) + 0.2714;
                  morr_arr(i,j,k,MORRInd::pgam) = 1.0/(std::pow(morr_arr(i,j,k,MORRInd::pgam), 2)) - 1.0;
                  morr_arr(i,j,k,MORRInd::pgam) = std::max(morr_arr(i,j,k,MORRInd::pgam), 2.0);
                  morr_arr(i,j,k,MORRInd::pgam) = std::min(morr_arr(i,j,k,MORRInd::pgam), 10.0);

                  // CALCULATE LAMC
                  morr_arr(i,j,k,MORRInd::lamc) = std::pow(m_cons26 * morr_arr(i,j,k,MORRInd::nc3d) * gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0) /
                                  (morr_arr(i,j,k,MORRInd::qc3d) * gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.0)), 1.0/3.0);

                  // LAMMIN, 60 MICRON DIAMETER
                  // LAMMAX, 1 MICRON
                  amrex::Real lammin = (morr_arr(i,j,k,MORRInd::pgam) + 1.0) / 60.0e-6;
                  amrex::Real lammax = (morr_arr(i,j,k,MORRInd::pgam) + 1.0) / 1.0e-6;

                  if (morr_arr(i,j,k,MORRInd::lamc) < lammin) {
                    morr_arr(i,j,k,MORRInd::lamc) = lammin;
                    morr_arr(i,j,k,MORRInd::nc3d) = std::exp(3.0 * std::log(morr_arr(i,j,k,MORRInd::lamc)) + std::log(morr_arr(i,j,k,MORRInd::qc3d)) +
                                           std::log(gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.0)) - std::log(gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0))) / m_cons26;
                  } else if (morr_arr(i,j,k,MORRInd::lamc) > lammax) {
                    morr_arr(i,j,k,MORRInd::lamc) = lammax;
                    morr_arr(i,j,k,MORRInd::nc3d) = std::exp(3.0 * std::log(morr_arr(i,j,k,MORRInd::lamc)) + std::log(morr_arr(i,j,k,MORRInd::qc3d)) +
                                           std::log(gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 1.0)) - std::log(gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0))) / m_cons26;
                  }
                }

                // SNOW
                if (morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                  morr_arr(i,j,k,MORRInd::lams) = std::pow(m_cons1 * morr_arr(i,j,k,MORRInd::ns3d) / morr_arr(i,j,k,MORRInd::qni3d), 1.0/m_ds);

                  // CHECK FOR SLOPE
                  // ADJUST VARS
                  if (morr_arr(i,j,k,MORRInd::lams) < m_lammins) {
                    morr_arr(i,j,k,MORRInd::lams) = m_lammins;
                    morr_arr(i,j,k,MORRInd::n0s) = std::pow(morr_arr(i,j,k,MORRInd::lams), 4) * morr_arr(i,j,k,MORRInd::qni3d) / m_cons1;
                    morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::n0s) / morr_arr(i,j,k,MORRInd::lams);
                  } else if (morr_arr(i,j,k,MORRInd::lams) > m_lammaxs) {
                    morr_arr(i,j,k,MORRInd::lams) = m_lammaxs;
                    morr_arr(i,j,k,MORRInd::n0s) = std::pow(morr_arr(i,j,k,MORRInd::lams), 4) * morr_arr(i,j,k,MORRInd::qni3d) / m_cons1;
                    morr_arr(i,j,k,MORRInd::ns3d) = morr_arr(i,j,k,MORRInd::n0s) / morr_arr(i,j,k,MORRInd::lams);
                  }
                }

                // GRAUPEL
                if (morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                  morr_arr(i,j,k,MORRInd::lamg) = std::pow(m_cons2 * morr_arr(i,j,k,MORRInd::ng3d) / morr_arr(i,j,k,MORRInd::qg3d), 1.0/m_dg);

                  // CHECK FOR SLOPE
                  // ADJUST VARS
                  if (morr_arr(i,j,k,MORRInd::lamg) < m_lamming) {
                    morr_arr(i,j,k,MORRInd::lamg) = m_lamming;
                    morr_arr(i,j,k,MORRInd::n0g) = std::pow(morr_arr(i,j,k,MORRInd::lamg), 4) * morr_arr(i,j,k,MORRInd::qg3d) / m_cons2;
                    morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::n0g) / morr_arr(i,j,k,MORRInd::lamg);
                  } else if (morr_arr(i,j,k,MORRInd::lamg) > m_lammaxg) {
                    morr_arr(i,j,k,MORRInd::lamg) = m_lammaxg;
                    morr_arr(i,j,k,MORRInd::n0g) = std::pow(morr_arr(i,j,k,MORRInd::lamg), 4) * morr_arr(i,j,k,MORRInd::qg3d) / m_cons2;
                    morr_arr(i,j,k,MORRInd::ng3d) = morr_arr(i,j,k,MORRInd::n0g) / morr_arr(i,j,k,MORRInd::lamg);
                  }
                }
              }

//            label_500:
              // CALCULATE EFFECTIVE RADIUS
              if (morr_arr(i,j,k,MORRInd::qi3d) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::effi) = 3.0 / morr_arr(i,j,k,MORRInd::lami) / 2.0 * 1.0e6;
              } else {
                morr_arr(i,j,k,MORRInd::effi) = 25.0;
              }

              if (morr_arr(i,j,k,MORRInd::qni3d) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::effs) = 3.0 / morr_arr(i,j,k,MORRInd::lams) / 2.0 * 1.0e6;
              } else {
                morr_arr(i,j,k,MORRInd::effs) = 25.0;
              }

              if (morr_arr(i,j,k,MORRInd::qr3d) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::effr) = 3.0 / morr_arr(i,j,k,MORRInd::lamr) / 2.0 * 1.0e6;
              } else {
                morr_arr(i,j,k,MORRInd::effr) = 25.0;
              }

              if (morr_arr(i,j,k,MORRInd::qc3d) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::effc) = gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 4.0) / gamma_function(morr_arr(i,j,k,MORRInd::pgam) + 3.0) / morr_arr(i,j,k,MORRInd::lamc) / 2.0 * 1.0e6;
              } else {
                morr_arr(i,j,k,MORRInd::effc) = 25.0;
              }

              if (morr_arr(i,j,k,MORRInd::qg3d) >= m_qsmall) {
                morr_arr(i,j,k,MORRInd::effg) = 3.0 / morr_arr(i,j,k,MORRInd::lamg) / 2.0 * 1.0e6;
              } else {
                morr_arr(i,j,k,MORRInd::effg) = 25.0;
              }

              // HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
              // TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
              // OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
              // HM, 12/28/12, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
              // OF EXCESSIVE AND PERSISTENT ANVIL
              // NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
              morr_arr(i,j,k,MORRInd::ni3d) = std::min(morr_arr(i,j,k,MORRInd::ni3d), 0.3e6 / morr_arr(i,j,k,MORRInd::rho));

              // ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
              if (iinum == 0 && m_iact == 2) {
                morr_arr(i,j,k,MORRInd::nc3d) = std::min(morr_arr(i,j,k,MORRInd::nc3d), (m_nanew1 + m_nanew2) / morr_arr(i,j,k,MORRInd::rho));
              }

              // SWITCH FOR CONSTANT DROPLET NUMBER
              if (iinum == 1) {
                // CHANGE NDCNST FROM CM-3 TO KG-1
                morr_arr(i,j,k,MORRInd::nc3d) = m_ndcnst * 1.0e6 / morr_arr(i,j,k,MORRInd::rho);
              }
            }

            }/* else {
              goto label_400;
            }
         label_400:*/
            //End of _micro
            if(use_morr_cpp_answer) {
              for(int k=klo; k<=khi; k++) {

            // Transfer 1D variables back to 3D arrays
            qcl_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qc3d);
            qci_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qi3d);
            qps_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qni3d);
            qpr_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qr3d);
            ni_arr(i,j,k) = morr_arr(i,j,k,MORRInd::ni3d);
            ns_arr(i,j,k) = morr_arr(i,j,k,MORRInd::ns3d);
            nr_arr(i,j,k) = morr_arr(i,j,k,MORRInd::nr3d);
            qpg_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qg3d);
            ng_arr(i,j,k) = morr_arr(i,j,k,MORRInd::ng3d);

            // Temperature and potential temperature conversion
            theta_arr(i,j,k) = morr_arr(i,j,k,MORRInd::t3d) / pii_arr(i,j,k); // Convert temp back to potential temp
            qv_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qv3d);

            //Deleted wrf-check, effc, and precr type data as not used by ERF
            /*
            // NEED gpu-compatabile summation for rain_accum, check SAM or Kessler for better example
            rain_accum_arr(i,j,k) = rain_accum_arr(i,j,k) + morr_arr(i,j,k,MORRInd::precrt);
            snow_accum_arr(i,j,k) = snow_accum_arr(i,j,k) + morr_arr(i,j,k,MORRInd::snowprt);
            graup_accum_arr(i,j,k) = graup_accum_arr(i,j,k) + morr_arr(i,j,k,MORRInd::grplprt);*/
            rainncv_arr(i,j,0) = morr_arr(i,j,klo,MORRInd::precrt);
            snowncv_arr(i,j,0) = morr_arr(i,j,klo,MORRInd::snowprt);
            graupelncv_arr(i,j,0) = morr_arr(i,j,klo,MORRInd::grplprt);
            sr_arr(i,j,0) = morr_arr(i,j,klo,MORRInd::snowrt) / (morr_arr(i,j,klo,MORRInd::precrt) + 1.e-12);
              }
            // Update precipitation accumulation variables
            // These are outside the k-loop in the original code
            rain_accum_arr(i,j,klo) = rain_accum_arr(i,j,klo) + morr_arr(i,j,klo,MORRInd::precrt);
            snow_accum_arr(i,j,klo) = snow_accum_arr(i,j,klo) + morr_arr(i,j,klo,MORRInd::snowprt);
            graup_accum_arr(i,j,klo) = graup_accum_arr(i,j,klo) + morr_arr(i,j,klo,MORRInd::grplprt);
            }
         });
          //          amrex::Print()<<FArrayBox(qv_arr)<<std::endl;
          }

          // amrex::Print()<<"fortran should run "<<run_morr_fort<<std::endl;

          if (run_morr_fort) {
#ifdef ERF_USE_MORR_FORT
          mp_morr_two_moment_c
          (
              1,  // ITIMESTEP - Use 1 for simplicity

              // 3D arrays in Fortran expected order (assume column-major for Fortran)
              theta_arr.dataPtr(),      // TH
              qv_arr.dataPtr(),         // QV
              qcl_arr.dataPtr(),        // QC
              qpr_arr.dataPtr(),        // QR
              qci_arr.dataPtr(),        // QI
              qps_arr.dataPtr(),        // QS
              qpg_arr.dataPtr(),        // QG
              ni_arr.dataPtr(),         // NI
              ns_arr.dataPtr(),         // NS
              nr_arr.dataPtr(),         // NR
              ng_arr.dataPtr(),         // NG

              rho_arr.dataPtr(),        // RHO
              pii_arr.dataPtr(),        // PII (Exner function)
              pres_arr.dataPtr(),       // P (in hPa, convert if needed)
              dt,                       // DT_IN
              dz_arr.dataPtr(),         // DZ
              w_arr.dataPtr(),          // W (vertical velocity)

              // 2D arrays for precipitation accounting
              rain_accum_arr.dataPtr(), // RAINNC
              rainncv_arr.dataPtr(),    // RAINNCV
              sr_arr.dataPtr(),         // SR
              snow_accum_arr.dataPtr(), // SNOWNC
              snowncv_arr.dataPtr(),    // SNOWNCV
              graup_accum_arr.dataPtr(),// GRAUPELNC
              graupelncv_arr.dataPtr(), // GRAUPELNCV

              // Radar reflectivity
              dummy_reflectivity_ptr,  // refl_10cm
              true,                     // diagflag
              0,   // do_radar_ref

              // Cumulus tendencies
              qrcuten_arr.dataPtr(),    // qrcuten
              qscuten_arr.dataPtr(),    // qscuten
              qicuten_arr.dataPtr(),    // qicuten

              // WRF-Chem flags
              flag_qndrop,              // F_QNDROP
              nullptr,                  // qndrop (not used here)
              ht_arr.dataPtr(),         // HT (terrain height - not used)

              // Domain dimensions
              ilo, ihi, jlo, jhi, klo, khi,  // IDS,IDE,JDS,JDE,KDS,KDE
              ilom, ihim, jlom, jhim, klom, khim,  // IMS,IME,JMS,JME,KMS,KME
              ilo, ihi, jlo, jhi, klo, khi,  // ITS,ITE,JTS,JTE,KTS,KTE

              // Optional WRF-Chem outputs
              false,                    // wetscav_on
              rainprod_arr.dataPtr(),   // rainprod
              evapprod_arr.dataPtr(),   // evapprod
              qlsink_arr.dataPtr(),     // QLSINK
              precr_arr.dataPtr(),      // PRECR
              preci_arr.dataPtr(),      // PRECI
              precs_arr.dataPtr(),      // PRECS
              precg_arr.dataPtr()       // PRECG
          );
#else
          amrex::Abort("Trying to run fortran without compiling with USE_MORR_FORT=TRUE");
#endif
        }
          //          amrex::Print()<<FArrayBox(qv_arr)<<std::endl;
          // After the call, all fields are updated
          // We don't need to copy results back since we passed direct pointers
          // to our class member arrays
        }
    }
