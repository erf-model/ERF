#ifndef ERF_ORB_COS_ZENITH_H
#define ERF_ORB_COS_ZENITH_H

#include <vector>
#include <cmath>

#include <ERF_Constants.H>

typedef double real;


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
orbital_decl (real& calday,
              real& eccen,
              real& mvelpp,
              real& lambm0,
              real& obliqr,
              real& delta,
              real& eccf)
{

    real lambm;  // Lambda m, mean long of perihelion (rad)
    real lmm;    // Intermediate argument involving lambm
    real lamb;   // Lambda, the earths long of perihelion
    real invrho; // Inverse normalized sun/earth distance
    real sinl;   // Sine of lmm
    static constexpr real dayspy = real(amrex::Real(365.0)); // Day per year
    static constexpr real ve     = real( amrex::Real(80.5)); // Calday of vernal equinox

    /*
     Compute eccentricity factor and solar declination using
     day value where a round day (such as amrex::Real(213.0)) refers to 0z at
     Greenwich longitude.

     Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
     Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
     35:2362-amrex::Real(2367.)

     To get the earths true longitude (position in orbit; lambda in Berger
     1978) which is necessary to find the eccentricity factor and declination,
     must first calculate the mean longitude (lambda m in Berger 1978) at
     the present day.  This is done by adding to lambm0 (the mean longitude
     at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
     an increment (delta lambda m in Berger 1978) that is the number of
     days past or before (a negative increment) the vernal equinox divided by
     the days in a model year times the 2*pi radians in a complete orbit.
    */

    lambm = lambm0 + (calday - ve)*two*PI/dayspy;
    lmm   = lambm  - mvelpp;

    // The earths true longitude, in radians, is then found from
    // the formula in Berger 1978:

    sinl  = std::sin(lmm);
    lamb  = lambm  + eccen*( two*sinl + eccen*( amrex::Real(1.25)*sin(two*lmm)
                                              + eccen*((amrex::Real(13.0)/amrex::Real(12.0))*std::sin(three*lmm) - fourth*sinl)));

    /*
     Using the obliquity, eccentricity, moving vernal equinox longitude of
     perihelion (plus), and earths true longitude, the declination (delta)
     and the normalized earth/sun distance (rho in Berger 1978; actually inverse
     rho will be used), and thus the eccentricity factor (eccf), can be
     calculated from formulas given in Berger amrex::Real(1978.)
    */

    invrho = (one + eccen*std::cos(lamb - mvelpp)) / (one - eccen*eccen);

    // Set solar declination and eccentricity factor

    delta  = std::asin(std::sin(obliqr)*std::sin(lamb));
    eccf   = invrho*invrho;
}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
orbital_params (int& iyear_AD,
                real& eccen,
                real& obliq,
                real& mvelp,
                real& obliqr,
                real& lambm0,
                real& mvelpp)
{
    /*
    !-------------------------------------------------------------------------------
    !
    ! Calculate earths orbital parameters using Dave Threshers formula which
    ! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term
    ! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy
    ! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
    !
    !-------------------------------------------------------------------------------
    */

    int poblen  = 47; // # of elements in series wrt obliquity
    int pecclen = 19; // # of elements in series wrt eccentricity
    int pmvelen = 78; // # of elements in series wrt vernal equinox
    static constexpr real psecdeg = real(one)/real(amrex::Real(3600.0)); // arc sec to deg conversion
    static constexpr real degrad = PI/real(amrex::Real(180.)); // degree to radian conversion factor

    // Cosine series data for computation of obliquity: amplitude (arc seconds),
    // rate (arc seconds/year), phase (degrees).
    // amplitudes for obliquity cos series
    std::vector<real> obamp =  {-amrex::Real(2462.2214466), -amrex::Real(857.3232075), -amrex::Real(629.3231835),
                                 -amrex::Real(414.2804924), -amrex::Real(311.7632587),  amrex::Real(308.9408604),
                                 -amrex::Real(162.5533601), -amrex::Real(116.1077911),  amrex::Real(101.1189923),
                                  -amrex::Real(67.6856209),   amrex::Real(24.9079067),   amrex::Real(22.5811241),
                                  -amrex::Real(21.1648355),  -amrex::Real(15.6549876),   amrex::Real(15.3936813),
                                   amrex::Real(14.6660938),  -amrex::Real(11.7273029),   amrex::Real(10.2742696),
                                    amrex::Real(6.4914588),    amrex::Real(5.8539148),   -amrex::Real(5.4872205),
                                   -amrex::Real(5.4290191),    amrex::Real(5.1609570),    amrex::Real(5.0786314),
                                   -amrex::Real(4.0735782),    amrex::Real(3.7227167),    amrex::Real(3.3971932),
                                   -amrex::Real(2.8347004),   -amrex::Real(2.6550721),   -amrex::Real(2.5717867),
                                   -amrex::Real(2.4712188),    amrex::Real(2.4625410),    amrex::Real(2.2464112),
                                   -amrex::Real(2.0755511),   -amrex::Real(1.9713669),   -amrex::Real(1.8813061),
                                   -amrex::Real(1.8468785),    amrex::Real(1.8186742),    amrex::Real(1.7601888),
                                   -amrex::Real(1.5428851),    amrex::Real(1.4738838),   -amrex::Real(1.4593669),
                                    amrex::Real(1.4192259),   -amrex::Real(1.1818980),    amrex::Real(1.1756474),
                                   -amrex::Real(1.1316126),    amrex::Real(1.0896928)};
    // rates for obliquity cosine series
    std::vector<real> obrate = {amrex::Real(31.609974), amrex::Real(32.620504), amrex::Real(24.172203),
                                amrex::Real(31.983787), amrex::Real(44.828336), amrex::Real(30.973257),
                                amrex::Real(43.668246), amrex::Real(32.246691), amrex::Real(30.599444),
                                amrex::Real(42.681324), amrex::Real(43.836462), amrex::Real(47.439436),
                                amrex::Real(63.219948), amrex::Real(64.230478),  amrex::Real(1.010530),
                                 amrex::Real(7.437771), amrex::Real(55.782177),  amrex::Real(0.373813),
                                amrex::Real(13.218362), amrex::Real(62.583231), amrex::Real(63.593761),
                                amrex::Real(76.438310), amrex::Real(45.815258),  amrex::Real(8.448301),
                                amrex::Real(56.792707), amrex::Real(49.747842), amrex::Real(12.058272),
                                amrex::Real(75.278220), amrex::Real(65.241008), amrex::Real(64.604291),
                                 amrex::Real(1.647247),  amrex::Real(7.811584), amrex::Real(12.207832),
                                amrex::Real(63.856665), amrex::Real(56.155990), amrex::Real(77.448840),
                                 amrex::Real(6.801054), amrex::Real(62.209418), amrex::Real(20.656133),
                                amrex::Real(48.344406), amrex::Real(55.145460), amrex::Real(69.000539),
                                amrex::Real(11.071350), amrex::Real(74.291298), amrex::Real(11.047742),
                                 amrex::Real(0.636717), amrex::Real(12.844549)};

    // phases for obliquity cosine series
    std::vector<real> obphas = {amrex::Real(251.9025), amrex::Real(280.8325), amrex::Real(128.3057),
                                amrex::Real(292.7252),  amrex::Real(15.3747), amrex::Real(263.7951),
                                amrex::Real(308.4258), amrex::Real(240.0099), amrex::Real(222.9725),
                                amrex::Real(268.7809), amrex::Real(316.7998), amrex::Real(319.6024),
                                amrex::Real(143.8050), amrex::Real(172.7351),  amrex::Real(28.9300),
                                amrex::Real(123.5968),  amrex::Real(20.2082),  amrex::Real(40.8226),
                                amrex::Real(123.4722), amrex::Real(155.6977), amrex::Real(184.6277),
                                amrex::Real(267.2772),  amrex::Real(55.0196), amrex::Real(152.5268),
                                 amrex::Real(49.1382), amrex::Real(204.6609),  amrex::Real(56.5233),
                                amrex::Real(200.3284), amrex::Real(201.6651), amrex::Real(213.5577),
                                 amrex::Real(17.0374), amrex::Real(164.4194),  amrex::Real(94.5422),
                                amrex::Real(131.9124),  amrex::Real(61.0309), amrex::Real(296.2073),
                                amrex::Real(135.4894), amrex::Real(114.8750), amrex::Real(247.0691),
                                amrex::Real(256.6114),  amrex::Real(32.1008), amrex::Real(143.6804),
                                 amrex::Real(16.8784), amrex::Real(160.6835),  amrex::Real(27.5932),
                                amrex::Real(348.1074),  amrex::Real(82.6496)};

    // Cosine/sine series data for computation of eccentricity and fixed vernal
    // equinox longitude of perihelion (fvelp): amplitude,
    // rate (arc seconds/year), phase (degrees).

    // ampl for eccen/fvelp cos/sin series
    std::vector<real> ecamp = { amrex::Real(0.01860798),  amrex::Real(0.01627522), -amrex::Real(0.01300660),
                                amrex::Real(0.00988829), -amrex::Real(0.00336700),  amrex::Real(0.00333077),
                               -amrex::Real(0.00235400),  amrex::Real(0.00140015),  amrex::Real(0.00100700),
                                amrex::Real(0.00085700),  amrex::Real(0.00064990),  amrex::Real(0.00059900),
                                amrex::Real(0.00037800), -amrex::Real(0.00033700),  amrex::Real(0.00027600),
                                amrex::Real(0.00018200), -amrex::Real(0.00017400), -amrex::Real(0.00012400),
                                amrex::Real(0.00001250)};

    // rates for eccen/fvelp cos/sin series
    std::vector<real> ecrate = { amrex::Real(4.2072050),  amrex::Real(7.3460910), amrex::Real(17.8572630),
                                amrex::Real(17.2205460), amrex::Real(16.8467330),  amrex::Real(5.1990790),
                                amrex::Real(18.2310760), amrex::Real(26.2167580),  amrex::Real(6.3591690),
                                amrex::Real(16.2100160),  amrex::Real(3.0651810), amrex::Real(16.5838290),
                                amrex::Real(18.4939800),  amrex::Real(6.1909530), amrex::Real(18.8677930),
                                amrex::Real(17.4255670),  amrex::Real(6.1860010), amrex::Real(18.4174410),
                                 amrex::Real(0.6678630)};

    // phases for eccen/fvelp cos/sin series
    std::vector<real> ecphas = { amrex::Real(28.620089), amrex::Real(193.788772), amrex::Real(308.307024),
                                amrex::Real(320.199637), amrex::Real(279.376984),  amrex::Real(87.195000),
                                amrex::Real(349.129677), amrex::Real(128.443387), amrex::Real(154.143880),
                                amrex::Real(291.269597), amrex::Real(114.860583), amrex::Real(332.092251),
                                amrex::Real(296.414411), amrex::Real(145.769910), amrex::Real(337.237063),
                                amrex::Real(152.092288), amrex::Real(126.839891), amrex::Real(210.667199),
                                 amrex::Real(72.108838)};

    // Sine series data for computation of moving vernal equinox longitude of
    // perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

    // amplitudes for mvelp sine series
    std::vector<real> mvamp = { amrex::Real(7391.0225890), amrex::Real(2555.1526947), amrex::Real(2022.7629188),
                               -amrex::Real(1973.6517951), amrex::Real(1240.2321818),  amrex::Real(953.8679112),
                                -amrex::Real(931.7537108),  amrex::Real(872.3795383),  amrex::Real(606.3544732),
                                -amrex::Real(496.0274038),  amrex::Real(456.9608039),  amrex::Real(346.9462320),
                                -amrex::Real(305.8412902),  amrex::Real(249.6173246), -amrex::Real(199.1027200),
                                 amrex::Real(191.0560889), -amrex::Real(175.2936572),  amrex::Real(165.9068833),
                                 amrex::Real(161.1285917),  amrex::Real(139.7878093), -amrex::Real(133.5228399),
                                 amrex::Real(117.0673811),  amrex::Real(104.6907281),   amrex::Real(95.3227476),
                                  amrex::Real(86.7824524),   amrex::Real(86.0857729),   amrex::Real(70.5893698),
                                 -amrex::Real(69.9719343),  -amrex::Real(62.5817473),   amrex::Real(61.5450059),
                                 -amrex::Real(57.9364011),   amrex::Real(57.1899832),  -amrex::Real(57.0236109),
                                 -amrex::Real(54.2119253),   amrex::Real(53.2834147),   amrex::Real(52.1223575),
                                 -amrex::Real(49.0059908),  -amrex::Real(48.3118757),  -amrex::Real(45.4191685),
                                 -amrex::Real(42.2357920),  -amrex::Real(34.7971099),   amrex::Real(34.4623613),
                                 -amrex::Real(33.8356643),   amrex::Real(33.6689362),  -amrex::Real(31.2521586),
                                 -amrex::Real(30.8798701),   amrex::Real(28.4640769),  -amrex::Real(27.1960802),
                                  amrex::Real(27.0860736),  -amrex::Real(26.3437456),   amrex::Real(24.7253740),
                                  amrex::Real(24.6732126),   amrex::Real(24.4272733),   amrex::Real(24.0127327),
                                  amrex::Real(21.7150294),  -amrex::Real(21.5375347),   amrex::Real(18.1148363),
                                 -amrex::Real(16.9603104),  -amrex::Real(16.1765215),   amrex::Real(15.5567653),
                                  amrex::Real(15.4846529),   amrex::Real(15.2150632),   amrex::Real(14.5047426),
                                 -amrex::Real(14.3873316),   amrex::Real(13.1351419),   amrex::Real(12.8776311),
                                  amrex::Real(11.9867234),   amrex::Real(11.9385578),   amrex::Real(11.7030822),
                                  amrex::Real(11.6018181),  -amrex::Real(11.2617293),  -amrex::Real(10.4664199),
                                  amrex::Real(10.4333970),  -amrex::Real(10.2377466),   amrex::Real(10.1934446),
                                 -amrex::Real(10.1280191),   amrex::Real(10.0289441),  -amrex::Real(10.0034259)};

    // rates for mvelp sine series
    std::vector<real> mvrate = {amrex::Real(31.609974), amrex::Real(32.620504), amrex::Real(24.172203),
                                 amrex::Real(0.636717), amrex::Real(31.983787),  amrex::Real(3.138886),
                                amrex::Real(30.973257), amrex::Real(44.828336),  amrex::Real(0.991874),
                                 amrex::Real(0.373813), amrex::Real(43.668246), amrex::Real(32.246691),
                                amrex::Real(30.599444),  amrex::Real(2.147012), amrex::Real(10.511172),
                                amrex::Real(42.681324), amrex::Real(13.650058),  amrex::Real(0.986922),
                                 amrex::Real(9.874455), amrex::Real(13.013341),  amrex::Real(0.262904),
                                 amrex::Real(0.004952),  amrex::Real(1.142024), amrex::Real(63.219948),
                                 amrex::Real(0.205021),  amrex::Real(2.151964), amrex::Real(64.230478),
                                amrex::Real(43.836462), amrex::Real(47.439436),  amrex::Real(1.384343),
                                 amrex::Real(7.437771), amrex::Real(18.829299),  amrex::Real(9.500642),
                                 amrex::Real(0.431696),  amrex::Real(1.160090), amrex::Real(55.782177),
                                amrex::Real(12.639528),  amrex::Real(1.155138),  amrex::Real(0.168216),
                                 amrex::Real(1.647247), amrex::Real(10.884985),  amrex::Real(5.610937),
                                amrex::Real(12.658184),  amrex::Real(1.010530),  amrex::Real(1.983748),
                                amrex::Real(14.023871),  amrex::Real(0.560178),  amrex::Real(1.273434),
                                amrex::Real(12.021467), amrex::Real(62.583231), amrex::Real(63.593761),
                                amrex::Real(76.438310),  amrex::Real(4.280910), amrex::Real(13.218362),
                                amrex::Real(17.818769),  amrex::Real(8.359495), amrex::Real(56.792707),
                                 amrex::Real(8.448301),  amrex::Real(1.978796),  amrex::Real(8.863925),
                                 amrex::Real(0.186365),  amrex::Real(8.996212),  amrex::Real(6.771027),
                                amrex::Real(45.815258), amrex::Real(12.002811), amrex::Real(75.278220),
                                amrex::Real(65.241008), amrex::Real(18.870667), amrex::Real(22.009553),
                                amrex::Real(64.604291), amrex::Real(11.498094),  amrex::Real(0.578834),
                                 amrex::Real(9.237738), amrex::Real(49.747842),  amrex::Real(2.147012),
                                 amrex::Real(1.196895),  amrex::Real(2.133898),  amrex::Real(0.173168)};

    // phases for mvelp sine series
    std::vector<real> mvphas = {amrex::Real(251.9025), amrex::Real(280.8325), amrex::Real(128.3057),
                                amrex::Real(348.1074), amrex::Real(292.7252), amrex::Real(165.1686),
                                amrex::Real(263.7951),  amrex::Real(15.3747),  amrex::Real(58.5749),
                                 amrex::Real(40.8226), amrex::Real(308.4258), amrex::Real(240.0099),
                                amrex::Real(222.9725), amrex::Real(106.5937), amrex::Real(114.5182),
                                amrex::Real(268.7809), amrex::Real(279.6869),  amrex::Real(39.6448),
                                amrex::Real(126.4108), amrex::Real(291.5795), amrex::Real(307.2848),
                                 amrex::Real(18.9300), amrex::Real(273.7596), amrex::Real(143.8050),
                                amrex::Real(191.8927), amrex::Real(125.5237), amrex::Real(172.7351),
                                amrex::Real(316.7998), amrex::Real(319.6024),  amrex::Real(69.7526),
                                amrex::Real(123.5968), amrex::Real(217.6432),  amrex::Real(85.5882),
                                amrex::Real(156.2147),  amrex::Real(66.9489),  amrex::Real(20.2082),
                                amrex::Real(250.7568),  amrex::Real(48.0188),   amrex::Real(8.3739),
                                 amrex::Real(17.0374), amrex::Real(155.3409),  amrex::Real(94.1709),
                                amrex::Real(221.1120),  amrex::Real(28.9300), amrex::Real(117.1498),
                                amrex::Real(320.5095), amrex::Real(262.3602), amrex::Real(336.2148),
                                amrex::Real(233.0046), amrex::Real(155.6977), amrex::Real(184.6277),
                                amrex::Real(267.2772),  amrex::Real(78.9281), amrex::Real(123.4722),
                                amrex::Real(188.7132), amrex::Real(180.1364),  amrex::Real(49.1382),
                                amrex::Real(152.5268),  amrex::Real(98.2198),  amrex::Real(97.4808),
                                amrex::Real(221.5376), amrex::Real(168.2438), amrex::Real(161.1199),
                                 amrex::Real(55.0196), amrex::Real(262.6495), amrex::Real(200.3284),
                                amrex::Real(201.6651), amrex::Real(294.6547),  amrex::Real(99.8233),
                                amrex::Real(213.5577), amrex::Real(154.1631), amrex::Real(232.7153),
                                amrex::Real(138.3034), amrex::Real(204.6609), amrex::Real(106.5938),
                                amrex::Real(250.4676), amrex::Real(332.3345),  amrex::Real(27.3039)};

    //---------------------------Local variables----------------------------------
    real obsum;   // Obliquity series summation
    real cossum;  // Cos series summation for eccentricity/fvelp
    real sinsum;  // Sin series summation for eccentricity/fvelp
    real fvelp;   // Fixed vernal equinox long of perihelion
    real mvsum;   // mvelp series summation
    real beta;    // Intermediate argument for lambm0
    real years;   // Years to time of interest ( pos <=> future)
    real eccen2;  // eccentricity squared
    real eccen3;  // eccentricity cubed
    real yb4_1950AD; // number of years before 1950 AD

    // Check for flag to use input orbit parameters
    if ( iyear_AD == ORB_UNDEF_INT ) {

        // The AMIP II settings (for a 1995 orbit) are adopted
        obliq  = amrex::Real(23.4441);
        eccen  = amrex::Real(0.016715);
        mvelp  = amrex::Real(102.7);
        eccen2 = eccen*eccen;
        eccen3 = eccen2*eccen;

    } else {  // Otherwise calculate based on years before present

        /*
          The following calculates the earths obliquity, orbital eccentricity
          (and various powers of it) and vernal equinox mean longitude of
          perihelion for years in the past (future = negative of years past),
          using constants (see parameter section) given in the program of:

          Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
          of Daily Insolation.  Contribution 18, Institute of Astronomy and
          Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.

          and formulas given in the paper (where less precise constants are also
          given):

          Berger, Andre.  amrex::Real(1978.)  Long-Term Variations of Daily Insolation and
          Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367

          The algorithm is valid only to 1,000,000 years past or hence.
          For a solution valid to 5-10 million years past see the above author.
          Algorithm below is better for years closer to present than is the
          5-10 million year solution.

          Years to time of interest must be negative of years before present
          (1950) in formulas that follow.
        */

        yb4_1950AD = real(amrex::Real(1950.0)) - real(iyear_AD);
        years = - yb4_1950AD;

        /*
          In the summations below, cosine or sine arguments, which end up in
          degrees, must be converted to radians via multiplication by degrad.

          Summation of cosine series for obliquity (epsilon in Berger 1978) in
          degrees. Convert the amplitudes and rates, which are in arc secs, into
          degrees via multiplication by psecdeg (arc seconds to degrees conversion
          factor).  For obliq, first term is Berger 1978 epsilon star; second
          term is series summation in degrees.
        */

        obsum = zero;
        for (int i(0); i<poblen; ++i) {
            obsum = obsum + obamp[i]*psecdeg*std::cos( (obrate[i]*psecdeg*years + obphas[i]) * degrad );
        }
        obliq = real(amrex::Real(23.320556)) + obsum;

        /*
          Summation of cosine and sine series for computation of eccentricity
          (eccen; e in Berger 1978) and fixed vernal equinox longitude of
          perihelion (fvelp; pi in Berger 1978), which is used for computation
          of moving vernal equinox longitude of perihelion.  Convert the rates,
          which are in arc seconds, into degrees via multiplication by psecdeg.
        */

        cossum = zero;
        for (int i(0); i<pecclen; ++i) {
            cossum = cossum + ecamp[i]*std::cos( (ecrate[i]*psecdeg*years+ecphas[i]) * degrad );
        }

        sinsum = zero;
        for (int i(0); i<pecclen; ++i) {
            sinsum = sinsum + ecamp[i]*std::sin( (ecrate[i]*psecdeg*years+ecphas[i]) * degrad );
        }

        // Use summations to calculate eccentricity

        eccen2 = cossum*cossum + sinsum*sinsum;
        eccen  = std::sqrt(eccen2);
        eccen3 = eccen2*eccen;

        // A series of cases for fvelp, which is in radians.
        if (std::fabs(cossum) <= amrex::Real(1.0e-8)) {
            if (sinsum == zero) {
                fvelp = zero;
            } else if (sinsum < zero) {
                fvelp = amrex::Real(1.5)*PI;
            } else if (sinsum > zero) {
                fvelp = myhalf*PI;
            }
        } else if (cossum < zero) {
            fvelp = std::atan(sinsum/cossum) + PI;
        } else {
            if (sinsum < zero) {
                fvelp = std::atan(sinsum/cossum) + two*PI;
            } else {
                fvelp = std::atan(sinsum/cossum);
            }
        }

        /*
          Summation of sin series for computation of moving vernal equinox long
          of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
          first term is fvelp in degrees; second term is Berger 1978 psi bar
          times years and in degrees; third term is Berger 1978 zeta; fourth
          term is series summation in degrees.  Convert the amplitudes and rates,
          which are in arc seconds, into degrees via multiplication by psecdeg.
          Series summation plus second and third terms constitute Berger 1978
          psi, which is the general precession.
        */
        mvsum = zero;
        for (int i(0); i<pmvelen; ++i) {
            mvsum = mvsum + mvamp[i]*psecdeg*std::sin( (mvrate[i]*psecdeg*years + mvphas[i]) * degrad);
        }
        mvelp = fvelp/degrad + real(amrex::Real(50.439273))*psecdeg*years + real(amrex::Real(3.392506)) + mvsum;

        // Cases to make sure mvelp is between 0 and amrex::Real(360.)
        if (mvelp < zero) {
            do {
                mvelp += amrex::Real(360.0);
            } while (mvelp < zero);
        }
        if (mvelp >= amrex::Real(360.0)) {
            do {
                mvelp -= amrex::Real(360.0);
            } while (mvelp >= amrex::Real(360.0));
        }

    }  // end of test on whether to calculate or use input orbital params

    // Orbit needs the obliquity in radians

    obliqr = obliq*degrad;

    /*
      180 degrees must be added to mvelp since observations are made from the
      earth and the sun is considered (wrongly for the algorithm) to go around
      the earth. For a more graphic explanation see Appendix B in:

      A. Berger, M. Loutre and C. Tricot. amrex::Real(1993.)  Insolation and Earth Orbital
      Periods.  J. of Geophysical Research 98:10,341-10,amrex::Real(362.)

      Additionally, orbit will need this value in radians. So mvelp becomes
      mvelpp (mvelp plus pi)
    */

    mvelpp = (mvelp + amrex::Real(180.0))*degrad;

    // Set up an argument used several times in lambm0 calculation ahead.

    beta = std::sqrt(one - eccen2);

    /*
      The mean longitude at the vernal equinox (lambda m nought in Berger
      1978; in radians) is calculated from the following formula given in
      Berger amrex::Real(1978.)  At the vernal equinox the true longitude (lambda in Berger
      1978) is zero
    */

    lambm0 = two*( (amrex::Real(.5)*eccen + amrex::Real(.125)*eccen3)*(one + beta)*std::sin(mvelpp)
                   - amrex::Real(.250)*eccen2*(amrex::Real(.5)    + beta)*std::sin(two*mvelpp)
                   + amrex::Real(.125)*eccen3*(one/three + beta)*std::sin(three*mvelpp) );
}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
real
orbital_avg_cos_zenith (real& jday,
                        real& lat,
                        real& lon,
                        real& declin,
                        real& dt_avg)
{
    // adjust latitude so that its tangent will be defined
    real del;
    if (lat == PIoTwo) {
        del = lat - real(amrex::Real(1.0e-05));
    } else if (lat == -PIoTwo) {
        del = lat + real(amrex::Real(1.0e-05));
    } else {
        del = lat;
    }

    // adjust declination so that its tangent will be defined
    real phi;
    if (declin == PIoTwo) {
        phi = declin - real(amrex::Real(1.0e-05));
    } else if (declin == -PIoTwo) {
        phi = declin + real(amrex::Real(1.0e-05));
    } else {
        phi = declin;
    }

    // define the cosine of the myhalf-day length
    // adjust for cases of all daylight or all night
    real h;
    real cos_h = - std::tan(del) * std::tan(phi);
    if (cos_h <= -one) {
        h = PI;
    } else if (cos_h >= one) {
        h = zero;
    } else {
        h = std::acos(cos_h);
    }

    // Define Local Time t and t + dt
    // adjust t to be between -pi and pi
    real t1 = (jday - int(jday)) * two*PI + lon - PI;
    if (t1 >= PI) {
        t1 = t1 - two*PI;
    } else if (t1 < -PI) {
        t1 = t1 + two*PI;
    }

    real dt = dt_avg / real(amrex::Real(86400.0)) * two*PI;
    real t2 = t1 + dt;

    // Compute Cosine Solar Zenith angle
    // define terms needed in the cosine zenith angle equation
    real aa = std::sin(lat) * std::sin(declin);
    real bb = std::cos(lat) * std::cos(declin);

    // define the hour angle
    // force it to be between -h and h
    // consider the situation when the night period is too short
    real tt1,tt2,tt3,tt4;
    if ( (t2 >= PI) && (t1 <= PI) && ((PI - h) <= dt) ) {
        tt2 = h;
        tt1 = std::min(std::max(t1,         -h),          h);
        tt4 = std::min(std::max(t2, two*PI - h), two*PI + h);
        tt3 = two*PI - h;
    } else if ( (t2 >= -PI) && (t1 <= -PI) && ((PI - h) <= dt) ) {
        tt2 = - two*PI + h;
        tt1 = std::min(std::max(t1, -two*PI - h), -two*PI + h);
        tt4 = std::min(std::max(t2,          -h),           h);
        tt3 = -h;
    } else {
        if (t2 > PI) {
            tt2 = std::min(std::max(t2 - two*PI, -h), h);
        } else if (t2 < -PI) {
            tt2 = std::min(std::max(t2 + two*PI, -h), h);
        } else {
            tt2 = std::min(std::max(t2         , -h), h);
        }

        if (t1 > PI) {
            tt1 = std::min(std::max(t1 - two*PI, -h), h);
        } else if (t1 < -PI) {
            tt1 = std::min(std::max(t1 + two*PI, -h), h);
        } else {
            tt1 = std::min(std::max(t1         , -h), h);
        }
        tt4 = zero;
        tt3 = zero;
    }

    // perform a time integration to obtain cosz if desired
    // output is valid over the period from t to t + dt
    if ( (tt2 > tt1) || (tt4 > tt3) ) {
       return (aa * (tt2 - tt1) + bb * (sin(tt2) - sin(tt1))) / dt +
              (aa * (tt4 - tt3) + bb * (sin(tt4) - sin(tt3))) / dt;
    } else {
        return zero;
    }
}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
real
orbital_cos_zenith (real& jday,
                    real& lat,
                    real& lon,
                    real& declin,
                    real dt_avg = -one,
                    real uniform_angle = -one,
                    real constant_zenith_angle_deg = -one)
{
    // Constant zenith angle is true
    if ( constant_zenith_angle_deg >= zero ) {
        return std::cos( constant_zenith_angle_deg * PI/amrex::Real(180.) );
    }

    // Uniform angle is true
    if ( uniform_angle >= zero) {
        return std::cos(uniform_angle);
    }

    // Perform the calculation of Orbital_Cos_Zenith
    bool use_dt_avg = false;
    if (dt_avg > zero) {
        use_dt_avg = true;
    }
    // If dt for the average cosz is specified, then call the shr_orb_avg_cosz
    if (use_dt_avg) {
        return orbital_avg_cos_zenith(jday, lat, lon, declin, dt_avg);
    } else {
        return  std::sin(lat)*std::sin(declin) - std::cos(lat)*std::cos(declin) *
                std::cos((jday-floor(jday))*real(two)*PI + lon);
    }
}
#endif
