#include <AMReX_GpuContainers.H>
#include "ERF_Morrison.H"
#include "ERF_IndexDefines.H"
#include "ERF_PlaneAverage.H"
#include "ERF_EOS.H"
#include "ERF_TileNoZ.H"

using namespace amrex;

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
    printf("wrf_gamma: Input value x = %g\n", x);

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

    printf("wrf_gamma: Initial y = %g\n", y);

    if (y <= zero) {
        // Argument is negative
        printf("wrf_gamma: Handling negative argument\n");
        y = -x;
        y1 = std::floor(y);
        res = y - y1;
        if (res != zero) {
            if (y1 != std::floor(y1 * half) * two)
                parity = true;
            Real pi=amrex::Math::pi<Real>();
            fact = -pi / std::sin(pi * res);
            y = y + one;
            printf("wrf_gamma: After reflection formula: y = %g, fact = %g, parity = %d\n", 
                   y, fact, parity);
        }
        else {
            printf("wrf_gamma: Singularity detected, returning xinf = %g\n", xinf);
            res = xinf;
            return res;
        }
    }

    // Argument is positive
    if (y < eps) {
        // Argument < eps
        printf("wrf_gamma: Small argument branch (y < eps)\n");
        if (y >= xminin) {
            res = one / y;
            printf("wrf_gamma: Small argument result: res = %g\n", res);
        }
        else {
            printf("wrf_gamma: Argument too small, returning xinf = %g\n", xinf);
            res = xinf;
            return res;
        }
    }
    else if (y < twelve) {
        // Medium range argument
        printf("wrf_gamma: Medium range branch (eps <= y < 12)\n");
        y1 = y;
        if (y < one) {
            // 0.0 < argument < 1.0
            printf("wrf_gamma: Sub-branch: 0 < y < 1\n");
            z = y;
            y = y + one;
        }
        else {
            // 1.0 < argument < 12.0, reduce argument if necessary
            n = static_cast<int>(y) - 1;
            printf("wrf_gamma: Sub-branch: 1 <= y < 12, n = %d\n", n);
            y = y - static_cast<amrex::Real>(n);
            z = y - one;
        }

        // Evaluate approximation
        printf("wrf_gamma: Before approximation: z = %g, y = %g\n", z, y);
        xnum = zero;
        xden = one;
        for (i = 0; i < 8; i++) {
            xnum = (xnum + p[i]) * z;
            xden = xden * z + q[i];
        }
        res = xnum / xden + one;
        printf("wrf_gamma: After approximation: res = %g\n", res);

        if (y1 < y) {
            // Adjust result for case 0.0 < argument < 1.0
            res = res / y1;
            printf("wrf_gamma: Adjusted for y < 1: res = %g\n", res);
        }
        else if (y1 > y) {
            // Adjust for 2.0 < argument < 12.0
            printf("wrf_gamma: Adjusting for y > 2 with %d multiplications\n", n);
            for (i = 0; i < n; i++) {
                res = res * y;
                y = y + one;
                printf("wrf_gamma: Multiplication %d: res = %g, y = %g\n", i+1, res, y);
            }
        }
    }
    else {
        // Large argument
        printf("wrf_gamma: Large argument branch (y >= 12)\n");
        if (y <= xbig) {
            ysq = y * y;
            sum = c[6];
            for (i = 0; i < 6; i++) {
                sum = sum / ysq + c[i];
                printf("wrf_gamma: Sum step %d: sum = %g\n", i+1, sum);
            }
            sum = sum / y - y + xxx;
            sum = sum + (y - half) * std::log(y);
            printf("wrf_gamma: Before exp: sum = %g\n", sum);
            res = std::exp(sum);
            printf("wrf_gamma: After exp: res = %g\n", res);
        }
        else {
            printf("wrf_gamma: Argument too large, returning xinf = %g\n", xinf);
            res = xinf;
            return res;
        }
    }

    // Final adjustments
    if (parity) {
        res = -res;
        printf("wrf_gamma: Applied parity adjustment: res = %g\n", res);
    }
    if (fact != one) {
        res = fact / res;
        printf("wrf_gamma: Applied reflection adjustment: res = %g\n", res);
    }
    
    printf("wrf_gamma: Final result = %g\n", res);
    return res;
}

/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 * @param[in] qc_in Cloud variables input
 * @param[in,out] qv_in Vapor variables input
 * @param[in] qi_in Ice variables input
 * @param[in] grids The boxes on which we will evolve the solution
 * @param[in] geom Geometry associated with these MultiFabs and grids
 * @param[in] dt_advance Timestep for the advance
 */
void
Morrison::Init (const MultiFab& cons_in,
                const BoxArray& grids,
                const Geometry& geom,
                const Real& dt_advance,
                std::unique_ptr<MultiFab>& z_phys_nd,
                std::unique_ptr<MultiFab>& detJ_cc)
{
    dt     = dt_advance;
    m_geom = geom;
    m_gtoe = grids;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc   = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_Morr::rain_accum, MicVar_Morr::snow_accum, MicVar_Morr::graup_accum};

    // initialize microphysics variables
    for (auto ivar = 0; ivar < MicVar_Morr::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.);
    }

    // Set class data members
    for ( MFIter mfi(cons_in, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        const auto& lo = lbound(box3d);
        const auto& hi = ubound(box3d);

        nlev = box3d.length(2);
        zlo  = lo.z;
        zhi  = hi.z;

        // parameters
        accrrc.resize({zlo},  {zhi});
        accrsi.resize({zlo},  {zhi});
        accrsc.resize({zlo},  {zhi});
        coefice.resize({zlo}, {zhi});
        evaps1.resize({zlo},  {zhi});
        evaps2.resize({zlo},  {zhi});
        accrgi.resize({zlo},  {zhi});
        accrgc.resize({zlo},  {zhi});
        evapg1.resize({zlo},  {zhi});
        evapg2.resize({zlo},  {zhi});
        evapr1.resize({zlo},  {zhi});
        evapr2.resize({zlo},  {zhi});

        // data (input)
        rho1d.resize({zlo}, {zhi});
        pres1d.resize({zlo}, {zhi});
        tabs1d.resize({zlo}, {zhi});
        gamaz.resize({zlo}, {zhi});
        zmid.resize({zlo}, {zhi});
    }

    int morr_rimed_ice = 0; // This is used to set something called "ihail"
    Print().SetPrecision(18)<<"wrf_gamma: "<<wrf_gamma(4+4.829501842840712377835644e+00)<<std::endl;
    morr_two_moment_init_c(morr_rimed_ice);
}


/**
 * Initializes the Microphysics module.
 *
 * @param[in] cons_in Conserved variables input
 */
void
Morrison::Copy_State_to_Micro (const MultiFab& cons_in)
{
    // Get the temperature, density, theta, qt and qp from input
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.growntilebox();

        auto states_array = cons_in.array(mfi);

        // Non-precipitating
        auto qv_array    = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
        auto qc_array    = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
        auto qi_array    = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
        auto qn_array    = mic_fab_vars[MicVar_Morr::qn]->array(mfi);
        auto qt_array    = mic_fab_vars[MicVar_Morr::qt]->array(mfi);

        // Precipitating
        auto qpr_array   = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
        auto qps_array   = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
        auto qpg_array   = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);
        auto qp_array    = mic_fab_vars[MicVar_Morr::qp]->array(mfi);

        auto rho_array   = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
        auto theta_array = mic_fab_vars[MicVar_Morr::theta]->array(mfi);
        auto tabs_array  = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
        auto pres_array  = mic_fab_vars[MicVar_Morr::pres]->array(mfi);

        // Get pressure, theta, temperature, density, and qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            theta_array(i,j,k) = states_array(i,j,k,RhoTheta_comp)/states_array(i,j,k,Rho_comp);

            qv_array(i,j,k)    = std::max(0.0,states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp));
            qc_array(i,j,k)    = std::max(0.0,states_array(i,j,k,RhoQ2_comp)/states_array(i,j,k,Rho_comp));
            qi_array(i,j,k)    = std::max(0.0,states_array(i,j,k,RhoQ3_comp)/states_array(i,j,k,Rho_comp));
            qn_array(i,j,k)    = qc_array(i,j,k) + qi_array(i,j,k);
            qt_array(i,j,k)    = qv_array(i,j,k) + qn_array(i,j,k);

            qpr_array(i,j,k)   = std::max(0.0,states_array(i,j,k,RhoQ4_comp)/states_array(i,j,k,Rho_comp));
            qps_array(i,j,k)   = std::max(0.0,states_array(i,j,k,RhoQ5_comp)/states_array(i,j,k,Rho_comp));
            qpg_array(i,j,k)   = std::max(0.0,states_array(i,j,k,RhoQ6_comp)/states_array(i,j,k,Rho_comp));
             qp_array(i,j,k)   = qpr_array(i,j,k) + qps_array(i,j,k) + qpg_array(i,j,k);

            tabs_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv_array(i,j,k));

            // NOTE: the Morrison Fortran version uses Pa not hPa so we don't divideby 100!
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv_array(i,j,k)); //  * 0.01;
        });
    }
}


void Morrison::Compute_Coefficients ()
{
    auto dz   = m_geom.CellSize(2);
    auto lowz = m_geom.ProbLo(2);

    auto accrrc_t  = accrrc.table();
    auto accrsi_t  = accrsi.table();
    auto accrsc_t  = accrsc.table();
    auto coefice_t = coefice.table();
    auto evaps1_t  = evaps1.table();
    auto evaps2_t  = evaps2.table();
    auto accrgi_t  = accrgi.table();
    auto accrgc_t  = accrgc.table();
    auto evapg1_t  = evapg1.table();
    auto evapg2_t  = evapg2.table();
    auto evapr1_t  = evapr1.table();
    auto evapr2_t  = evapr2.table();

    auto rho1d_t  = rho1d.table();
    auto pres1d_t = pres1d.table();
    auto tabs1d_t = tabs1d.table();

    auto gamaz_t  = gamaz.table();
    auto zmid_t   = zmid.table();

    Real gam3  = erf_gammafff(3.0             );
    Real gamr1 = erf_gammafff(3.0+b_rain      );
    Real gamr2 = erf_gammafff((5.0+b_rain)/2.0);
    Real gams1 = erf_gammafff(3.0+b_snow      );
    Real gams2 = erf_gammafff((5.0+b_snow)/2.0);
    Real gamg1 = erf_gammafff(3.0+b_grau      );
    Real gamg2 = erf_gammafff((5.0+b_grau)/2.0);

    // calculate the plane average variables
    PlaneAverage rho_ave(mic_fab_vars[MicVar_Morr::rho].get(), m_geom, m_axis);
    PlaneAverage theta_ave(mic_fab_vars[MicVar_Morr::theta].get(), m_geom, m_axis);
    PlaneAverage qv_ave(mic_fab_vars[MicVar_Morr::qv].get(), m_geom, m_axis);
    rho_ave.compute_averages(ZDir(), rho_ave.field());
    theta_ave.compute_averages(ZDir(), theta_ave.field());
    qv_ave.compute_averages(ZDir(), qv_ave.field());

    // get host variable rho, and rhotheta
    int ncell = rho_ave.ncell_line();

    Gpu::HostVector<Real> rho_h(ncell), theta_h(ncell), qv_h(ncell);
    rho_ave.line_average(0, rho_h);
    theta_ave.line_average(0, theta_h);
    qv_ave.line_average(0, qv_h);

    // copy data to device
    Gpu::DeviceVector<Real> rho_d(ncell), theta_d(ncell), qv_d(ncell);
    Gpu::copyAsync(Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    Gpu::streamSynchronize();

    Real* rho_dptr   = rho_d.data();
    Real* theta_dptr = theta_d.data();
    Real* qv_dptr    = qv_d.data();

    Real gOcp = m_gOcp;

    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real RhoTheta = rho_dptr[k]*theta_dptr[k];
        Real pressure = getPgivenRTh(RhoTheta, qv_dptr[k]);
        rho1d_t(k)    = rho_dptr[k];
        pres1d_t(k)   = pressure*0.01;
        // NOTE: Limit the temperature to the melting point of ice to avoid a divide by
        //       0 condition when computing the cold evaporation coefficients. This should
        //       not affect results since evporation requires snow/graupel to be present
        //       and thus T<273.16
        tabs1d_t(k)   = std::min(getTgivenRandRTh(rho_dptr[k], RhoTheta, qv_dptr[k]),273.16);
        zmid_t(k)     = lowz + (k+0.5)*dz;
        gamaz_t(k)    = gOcp*zmid_t(k);
    });

    if(round(gam3) != 2) {
        std::cout << "cannot compute gamma-function in Microphysics::Init" << std::endl;
        std::exit(-1);
    }

    // Populate all the coefficients
    ParallelFor(nlev, [=] AMREX_GPU_DEVICE (int k) noexcept
    {
        Real Prefactor;
        Real pratio = sqrt(1.29 / rho1d_t(k));
        //Real rrr1   = 393.0/(tabs1d_t(k)+120.0)*std::pow((tabs1d_t(k)/273.0),1.5);
        //Real rrr2   = std::pow((tabs1d_t(k)/273.0),1.94)*(1000.0/pres1d_t(k));
        Real estw   = 100.0*erf_esatw(tabs1d_t(k));
        Real esti   = 100.0*erf_esati(tabs1d_t(k));

        // accretion by snow:
        Real coef1   = 0.25 * PI * nzeros * a_snow * gams1 * pratio/pow((PI * rhos * nzeros/rho1d_t(k) ) , ((3.0+b_snow)/4.0));
        Real coef2   = exp(0.025*(tabs1d_t(k) - 273.15));
        accrsi_t(k)  =  coef1 * coef2 * esicoef;
        accrsc_t(k)  =  coef1 * esccoef;
        coefice_t(k) =  coef2;

        // evaporation of snow:
        coef1 = (lsub/(tabs1d_t(k)*R_v)-1.0)*lsub/(therco*tabs1d_t(k));
        coef2 = R_v * R_d / (diffelq * esti);
        Prefactor = 2.0 * PI * nzeros / (rho1d_t(k) * (coef1 + coef2));
        Prefactor *= (2.0/PI); // Shape factor snow
        evaps1_t(k) = Prefactor * 0.65 * sqrt(rho1d_t(k) / (PI * rhos * nzeros));
        evaps2_t(k) = Prefactor * 0.44 * sqrt(a_snow * rho1d_t(k) / muelq) * gams2
                    * sqrt(pratio) * pow(rho1d_t(k) / (PI * rhos * nzeros) , ((5.0+b_snow)/8.0));

        // accretion by graupel:
        coef1 = 0.25*PI*nzerog*a_grau*gamg1*pratio/pow((PI*rhog*nzerog/rho1d_t(k)) , ((3.0+b_grau)/4.0));
        coef2 = exp(0.025*(tabs1d_t(k) - 273.15));
        accrgi_t(k) = coef1 * coef2 * egicoef;
        accrgc_t(k) = coef1 * egccoef;

        // evaporation of graupel:
        coef1 = (lsub/(tabs1d_t(k)*R_v)-1.0)*lsub/(therco*tabs1d_t(k));
        coef2 = R_v * R_d / (diffelq * esti);
        Prefactor = 2.0 * PI * nzerog / (rho1d_t(k) * (coef1 + coef2)); // Shape factor for graupel is 1
        evapg1_t(k) = Prefactor * 0.78 * sqrt(rho1d_t(k) / (PI * rhog * nzerog));
        evapg2_t(k) = Prefactor * 0.31 * sqrt(a_grau * rho1d_t(k) / muelq) * gamg2
                    * sqrt(pratio) * pow(rho1d_t(k) / (PI * rhog * nzerog) , ((5.0+b_grau)/8.0));

        // accretion by rain:
        accrrc_t(k) = 0.25 * PI * nzeror * a_rain * gamr1 * pratio/pow((PI * rhor * nzeror / rho1d_t(k)) , ((3.0+b_rain)/4.))* erccoef;

        // evaporation of rain:
        coef1 = (lcond/(tabs1d_t(k)*R_v)-1.0)*lcond/(therco*tabs1d_t(k));
        coef2 = R_v * R_d / (diffelq * estw);
        Prefactor = 2.0 * PI * nzeror / (rho1d_t(k) * (coef1 + coef2)); // Shape factor for rain is 1
        evapr1_t(k) = Prefactor * 0.78 * sqrt(rho1d_t(k) / (PI * rhor * nzeror));
        evapr2_t(k) = Prefactor * 0.31 * sqrt(a_rain * rho1d_t(k) / muelq) * gamr2
                    * sqrt(pratio) * pow(rho1d_t(k) / (PI * rhor * nzeror) , ((5.0+b_rain)/8.0));
    });
}
