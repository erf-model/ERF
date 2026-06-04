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

constexpr Real xxx = Real(0.9189385332046727417803297);
/*
!------------------------------------------------------------------------------

      REAL(C_DOUBLE) FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL(C_DOUBLE) ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE one
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE two
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
!          one+EPS .GT. one
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        Real(966.961)
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        Real(177.803)
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        Real(35.040)
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        Real(171.624)
! IBM 3033       (D.P.)       16           63        Real(57.574)
! VAX D-FORMAT   (D.P.)        2          127        Real(34.844)
! VAX G-FORMAT   (D.P.)        2         1023        Real(171.489)
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   Real(5.45E+2465)   Real(7.11E-15)    Real(1.84E-2466)
! CYBER 180/855
!   UNDER NOS    (S.P.)   Real(1.26E+322)    Real(3.55E-15)    Real(3.14E-294)
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   Real(3.40E+38)     Real(1.19E-7)     Real(1.18E-38)
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
!              (ED.), SPRINGER VERLAG, BERLIN, Real(1976.)
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, Real(1968.)
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
    amrex::Real fact, res, sum, twelve, xbig, xden, xinf, xminin;
    amrex::Real xnum, y, y1, ysq, z;
    amrex::Real c[7];
    amrex::Real p[8];
    amrex::Real q[8];

    // Mathematical constants
    twelve = Real(12.0);

    // Machine dependent parameters
    xbig = Real(35.040);
    xminin = Real(1.18e-38);
    amrex::Real eps = Real(1.19e-7);
    xinf = Real(3.4e38);

    // Numerator and denominator coefficients for rational minimax approximation over (1,2)
    p[0] = -Real(1.71618513886549492533811e+0);
    p[1] =  Real(2.47656508055759199108314e+1);
    p[2] = -Real(3.79804256470945635097577e+2);
    p[3] =  Real(6.29331155312818442661052e+2);
    p[4] =  Real(8.66966202790413211295064e+2);
    p[5] = -Real(3.14512729688483675254357e+4);
    p[6] = -Real(3.61444134186911729807069e+4);
    p[7] =  Real(6.64561438202405440627855e+4);

    q[0] = -Real(3.08402300119738975254353e+1);
    q[1] =  Real(3.15350626979604161529144e+2);
    q[2] = -Real(1.01515636749021914166146e+3);
    q[3] = -Real(3.10777167157231109440444e+3);
    q[4] =  Real(2.25381184209801510330112e+4);
    q[5] =  Real(4.75584627752788110767815e+3);
    q[6] = -Real(1.34659959864969306392456e+5);
    q[7] = -Real(1.15132259675553483497211e+5);

    // Coefficients for minimax approximation over (12, inf)
    c[0] = -Real(1.910444077728e-03);
    c[1] = Real(8.4171387781295e-04);
    c[2] = -Real(5.952379913043012e-04);
    c[3] = Real(7.93650793500350248e-04);
    c[4] = -Real(2.777777777777681622553e-03);
    c[5] = Real(8.333333333333333331554247e-02);
    c[6] = Real(5.7083835261e-03);

    // Initialize variables
    parity = false;
    fact = one;
    n = 0;
    y = x;

//    printf("wrf_gamma: Initial y = %g\n", y);

    if (y <= Real(0)) {
        // Argument is negative
//        printf("wrf_gamma: Handling negative argument\n");
        y = -x;
        y1 = std::floor(y);
        res = y - y1;
        if (res != Real(0)) {
            if (y1 != std::floor(y1 * myhalf) * Real(2))
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
            // Real(0) < argument < one
//            printf("wrf_gamma: Sub-branch: 0 < y < 1\n");
            z = y;
            y = y + one;
        }
        else {
            // one < argument < Real(12.0), reduce argument if necessary
            n = static_cast<int>(y) - 1;
//            printf("wrf_gamma: Sub-branch: 1 <= y < 12, n = %d\n", n);
            y = y - static_cast<amrex::Real>(n);
            z = y - one;
        }

        // Evaluate approximation
//        printf("wrf_gamma: Before approximation: z = %g, y = %g\n", z, y);
        xnum = Real(0);
        xden = one;
        for (i = 0; i < 8; i++) {
            xnum = (xnum + p[i]) * z;
            xden = xden * z + q[i];
        }
        res = xnum / xden + one;
//        printf("wrf_gamma: After approximation: res = %g\n", res);

        if (y1 < y) {
            // Adjust result for case Real(0) < argument < one
            res = res / y1;
//            printf("wrf_gamma: Adjusted for y < 1: res = %g\n", res);
        }
        else if (y1 > y) {
            // Adjust for two < argument < Real(12.0)
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
            sum = sum + (y - myhalf) * std::log(y);
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

