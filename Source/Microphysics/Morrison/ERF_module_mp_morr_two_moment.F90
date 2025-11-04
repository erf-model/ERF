!WRF:MODEL_LAYER:PHYSICS
!

! THIS MODULE CONTAINS THE TWO-MOMENT MICROPHYSICS CODE DESCRIBED BY
!     MORRISON ET AL. (2009, MWR)

! CHANGES FOR V3.2, RELATIVE TO MOST RECENT (BUG-FIX) CODE FOR V3.1

! 1) ADDED ACCELERATED MELTING OF GRAUPEL/SNOW DUE TO COLLISION WITH RAIN, FOLLOWING LIN ET AL. (1983)
! 2) INCREASED MINIMUM LAMBDA FOR RAIN, AND ADDED RAIN DROP BREAKUP FOLLOWING MODIFIED VERSION
!     OF VERLINDE AND COTTON (1993)
! 3) CHANGE MINIMUM ALLOWED MIXING RATIOS IN DRY CONDITIONS (RH < 90%), THIS IMPROVES RADAR REFLECTIIVITY
!     IN LOW REFLECTIVITY REGIONS
! 4) BUG FIX TO MAXIMUM ALLOWED PARTICLE FALLSPEEDS AS A FUNCTION OF AIR DENSITY
! 5) BUG FIX TO CALCULATION OF LIQUID WATER SATURATION VAPOR PRESSURE (CHANGE IS VERY MINOR)
! 6) INCLUDE WRF CONSTANTS PER SUGGESTION OF JIMY

! bug fix, 5/12/10
! 7) bug fix for saturation vapor pressure in low pressure, to avoid division by zero
! 8) include 'EP2' WRF constant for saturation mixing ratio calculation, instead of hardwire constant

! CHANGES FOR V3.3
! 1) MODIFICATION FOR COUPLING WITH WRF-CHEM (PREDICTED DROPLET NUMBER CONCENTRATION) AS AN OPTION
! 2) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 3) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 4) CLEAN UP OF COMMENTS IN THE CODE

! additional minor bug fixes and small changes, 5/30/2011
! minor revisions by A. Ackerman April 2011:
! 1) replaced kinematic with dynamic viscosity
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! additional fixes:
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NOT SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SINCE
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS
! 7) fix for switch for running w/o graupel/hail (cloud ice and snow only)

! hm bug fix 3/16/12

! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted

! WRFV3.5
! hm/A. Ackerman bug fix 11/08/12

! 1) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass
!    collected by rain
! 2) minor changes to some comments
! 3) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!    cm-3. This was done to address the problem of excessive and persistent
!    anvil cirrus produced by the scheme.

! CHANGES FOR WRFV3.5.1
! 1) added output for snow+cloud ice and graupel time step and accumulated
!    surface precipitation
! 2) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!    Hailong Wang (PNNL)
! 3) very minor fix to immersion freezing rate formulation (negligible impact)
! 4) clarifications to code comments
! 5) minor change to shedding of rain, remove limit so that the number of
!    collected drops can smaller than number of shed drops
! 6) change of specific heat of liquid water from 4218 to 4187 J/kg/K

! CHANGES FOR WRFV3.6.1
! 1) minor bug fix to melting of snow and graupel, an extra factor of air density (RHO) was removed
!    from the calculation of PSMLT and PGMLT
! 2) redundant initialization of PSMLT (non answer-changing)

! CHANGES FOR WRFV3.8.1
! 1) changes and cleanup of code comments
! 2) correction to universal gas constant (very small change)

! CHANGES FOR WRFV4.3
! 1) fix to saturation vapor pressure polysvp to work at T < -80 C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

MODULE MODULE_MP_MORR_TWO_MOMENT
!   USE     module_wrf_error
!   USE module_mp_radar

    USE ISO_C_BINDING

!  USE WRF PHYSICS CONSTANTS
  use module_model_constants, ONLY: CP, G, R => r_d, RV => r_v, EP_2

!  USE module_state_description

   IMPLICIT NONE

   REAL(C_DOUBLE), PARAMETER :: PI = 3.1415926535897932384626434
   REAL(C_DOUBLE), PARAMETER :: xxx = 0.9189385332046727417803297

   PUBLIC  ::  MP_MORR_TWO_MOMENT
   PUBLIC  ::  POLYSVP
   PUBLIC  ::  SET_MORRISON_NDCNST

   ! PRIVATE :: GAMMA, DERF1
   PRIVATE :: GAMMA
   PRIVATE :: PI, xxx
   PRIVATE :: MORR_TWO_MOMENT_MICRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SWITCHES FOR MICROPHYSICS SCHEME
! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! IACT = 3, ACTIVATION CALCULATED IN MODULE_MIXACTIVATE

     INTEGER, PRIVATE ::  IACT

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION
! !!!NOTE: PREDICTED DROPLET CONCENTRATION NOT AVAILABLE IN THIS VERSION
! CONTACT HUGH MORRISON (morrison@ucar.edu) FOR FURTHER INFORMATION

     INTEGER, PRIVATE ::  INUM

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
     REAL(C_DOUBLE), PRIVATE ::      NDCNST

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

     INTEGER, PRIVATE ::  ILIQ

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS

     INTEGER, PRIVATE ::  INUC

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING
!             NON-EQULIBRIUM SUPERSATURATION,
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION, BASED ON THE
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

     INTEGER, PRIVATE ::  IBASE

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

     INTEGER, PRIVATE ::  ISUB

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

     INTEGER, PRIVATE ::  IGRAUP

! HM ADDED NEW OPTION FOR HAIL
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

     INTEGER, PRIVATE ::  IHAIL

! CLOUD MICROPHYSICS CONSTANTS

     REAL(C_DOUBLE), PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL(C_DOUBLE), PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
!     REAL(C_DOUBLE), PRIVATE ::      R           ! GAS CONSTANT FOR AIR
!     REAL(C_DOUBLE), PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
!     REAL(C_DOUBLE), PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
     REAL(C_DOUBLE), PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
     REAL(C_DOUBLE), PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
     REAL(C_DOUBLE), PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
     REAL(C_DOUBLE), PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
     REAL(C_DOUBLE), PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
     REAL(C_DOUBLE), PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL(C_DOUBLE), PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL(C_DOUBLE), PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
     REAL(C_DOUBLE), PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
     REAL(C_DOUBLE), PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
     REAL(C_DOUBLE), PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
     REAL(C_DOUBLE), PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
     REAL(C_DOUBLE), PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
     REAL(C_DOUBLE), PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
     REAL(C_DOUBLE), PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
!     REAL(C_DOUBLE), PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
     REAL(C_DOUBLE), PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
     REAL(C_DOUBLE), PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
     REAL(C_DOUBLE), PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
     REAL(C_DOUBLE), PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
     REAL(C_DOUBLE), PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
! hm, add for V3.2
     REAL(C_DOUBLE), PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

! CCN SPECTRA FOR IACT = 1

     REAL(C_DOUBLE), PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
     REAL(C_DOUBLE), PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

! AEROSOL PARAMETERS FOR IACT = 2

     REAL(C_DOUBLE), PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
     REAL(C_DOUBLE), PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
     REAL(C_DOUBLE), PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
     REAL(C_DOUBLE), PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
     REAL(C_DOUBLE), PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
     REAL(C_DOUBLE), PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
     REAL(C_DOUBLE), PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
     REAL(C_DOUBLE), PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
     REAL(C_DOUBLE), PRIVATE ::      BACT    ! ACTIVATION PARAMETER
     REAL(C_DOUBLE), PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
     REAL(C_DOUBLE), PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
     REAL(C_DOUBLE), PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
     REAL(C_DOUBLE), PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
     REAL(C_DOUBLE), PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
     REAL(C_DOUBLE), PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
     REAL(C_DOUBLE), PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL(C_DOUBLE), PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL(C_DOUBLE), PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL(C_DOUBLE), PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL(C_DOUBLE), PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
     REAL(C_DOUBLE), PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING

! CONSTANTS TO IMPROVE EFFICIENCY

     REAL(C_DOUBLE), PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
     REAL(C_DOUBLE), PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
     REAL(C_DOUBLE), PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
     REAL(C_DOUBLE), PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
     REAL(C_DOUBLE), PRIVATE :: CONS41


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MORR_TWO_MOMENT_INIT(morr_rimed_ice, morr_noice) ! RAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      INTEGER, INTENT(IN):: morr_rimed_ice ! RAS
      INTEGER, INTENT(IN):: morr_noice

      integer n,i

      print *,'IN MORR_TWO_MOMENT_INIT ',"with noice = ",morr_noice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
! SET PRIOR TO CODE COMPILATION

! INUM IS AUTOMATICALLY SET TO 0 FOR WRF-CHEM BELOW,
! ALLOWING PREDICTION OF DROPLET CONCENTRATION
! THUS, THIS PARAMETER SHOULD NOT BE CHANGED HERE
! AND SHOULD BE LEFT TO 1

      INUM = 1

! SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)
! IF NO COUPLING WITH WRF-CHEM

      NDCNST = 250.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE, THE FOLLOWING OPTIONS RELATED TO DROPLET ACTIVATION
! (IACT, IBASE, ISUB) ARE NOT AVAILABLE IN CURRENT VERSION
! FOR WRF-CHEM, DROPLET ACTIVATION IS PERFORMED
! IN 'MIX_ACTIVATE', NOT IN MICROPHYSICS SCHEME


! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

      IACT = 2

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING
!             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER,
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      IBASE = 2

! INCLUDE SUB-GRID VERTICAL VELOCITY (standard deviation of w) IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! currently, sub-grid w is constant of 0.5 m/s (not coupled with PBL/turbulence scheme)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      ISUB = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(morr_noice .eq. 1) then

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      ILIQ = 1

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

      INUC = 0

! SWITCH FOR GRAUPEL/HAIL NO GRAUPEL/HAIL
! IGRAUP = 0, INCLUDE GRAUPEL/HAIL
! IGRAUP = 1, NO GRAUPEL/HAIL

      IGRAUP = 1

! HM ADDED 11/7/07
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL
! NOTE ---> RECOMMEND IHAIL = 1 FOR CONTINENTAL DEEP CONVECTION

      !IHAIL = 0 !changed to namelist option (morr_rimed_ice) by RAS
      ! Check if namelist option is feasible, otherwise default to graupel - RAS
      IF (morr_rimed_ice .eq. 1) THEN
         IHAIL = 1
      ELSE
         IHAIL = 0
      ENDIF
else

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      ILIQ = 0

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

      INUC = 0

! SWITCH FOR GRAUPEL/HAIL NO GRAUPEL/HAIL
! IGRAUP = 0, INCLUDE GRAUPEL/HAIL
! IGRAUP = 1, NO GRAUPEL/HAIL

      IGRAUP = 0

! HM ADDED 11/7/07
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL
! NOTE ---> RECOMMEND IHAIL = 1 FOR CONTINENTAL DEEP CONVECTION

      !IHAIL = 0 !changed to namelist option (morr_rimed_ice) by RAS
      ! Check if namelist option is feasible, otherwise default to graupel - RAS
      IF (morr_rimed_ice .eq. 1) THEN
         IHAIL = 1
      ELSE
         IHAIL = 0
      ENDIF
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET PHYSICAL CONSTANTS

! FALLSPEED PARAMETERS (V=AD^B)
         AI = 700.
         AC = 3.E7
         AS = 11.72
         AR = 841.99667
         BI = 1.
         BC = 2.
         BS = 0.41
         BR = 0.8
         IF (IHAIL.EQ.0) THEN
             AG = 19.3
             BG = 0.37
         ELSE ! (MATSUN AND HUGGINS 1980)
             AG = 114.5
             BG = 0.5
         END IF

! CONSTANTS AND PARAMETERS
!        R = 287.15
!        RV = 461.5
!        CP = 1005.
         RHOSU = 85000./(287.15*273.15)
         RHOW = 997.
         RHOI = 500.
         RHOSN = 100.
         IF (IHAIL.EQ.0) THEN
             RHOG = 400.
         ELSE
             RHOG = 900.
         END IF
         AIMM = 0.66
         BIMM = 100.
         ECR = 1.
         DCS = 125.E-6
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3
         MG0 = 1.6E-10
         F1S = 0.86
         F2S = 0.28
         F1R = 0.78
!        F2R = 0.32
! fix 053011
         F2R = 0.308
!        G = 9.806

         QSMALL = 1.d-14
         print *,'SETTING SMALL TO ',QSMALL

         EII = 0.1
         ECI = 0.7
! HM, ADD FOR V3.2
! hm, 7/23/13
!         CPW = 4218.
         CPW = 4187.

! SIZE DISTRIBUTION PARAMETERS

         CI = RHOI*PI/6.
         DI = 3.
         CS = RHOSN*PI/6.
         DS = 3.
         CG = RHOG*PI/6.
         DG = 3.

! RADIUS OF CONTACT NUCLEI
         RIN = 0.1E-6

         MMULT = 4./3.*PI*RHOI*(5.E-6)**3

! SIZE LIMITS FOR LAMBDA

         LAMMAXI = 1./1.E-6
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./20.E-6
!         LAMMINR = 1./500.E-6
         LAMMINR = 1./2800.E-6

         LAMMAXS = 1./10.E-6
         LAMMINS = 1./2000.E-6
         LAMMAXG = 1./20.E-6
         LAMMING = 1./2000.E-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! note: these parameters only used by the non-wrf-chem version of the
!       scheme with predicted droplet number

! CCN SPECTRA FOR IACT = 1

! MARITIME
! MODIFIED FROM RASMUSSEN ET AL. 2002
! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

              K1 = 0.4
              C1 = 120.

! CONTINENTAL

!              K1 = 0.5
!              C1 = 1000.

! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

         MW = 0.018
         OSM = 1.
         VI = 3.
         EPSM = 0.7
         RHOA = 1777.
         MAP = 0.132
         MA = 0.0284
! hm fix 6/23/16
!         RR = 8.3187
         RR = 8.3145
         BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE
! (see morrison et al. 2007, JGR)
! MODE 1

         RM1 = 0.052E-6
         SIG1 = 2.04
         NANEW1 = 72.2E6
         F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
         F21 = 1.+0.25*LOG(SIG1)

! MODE 2

         RM2 = 1.3E-6
         SIG2 = 2.5
         NANEW2 = 1.8E6
         F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
         F22 = 1.+0.25*LOG(SIG2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CONSTANTS FOR EFFICIENCY

         CONS1=GAMMA(1.+DS)*CS
         CONS2=GAMMA(1.+DG)*CG
         CONS3=GAMMA(4.+BS)/6.
         CONS4=GAMMA(4.+BR)/6.
         CONS5=GAMMA(1.+BS)
         CONS6=GAMMA(1.+BR)
         CONS7=GAMMA(4.+BG)/6.
         CONS8=GAMMA(1.+BG)
         CONS9=GAMMA(5./2.+BR/2.)
         CONS10=GAMMA(5./2.+BS/2.)
         CONS11=GAMMA(5./2.+BG/2.)
         CONS12=GAMMA(1.+DI)*CI
         CONS13=GAMMA(BS+3.)*PI/4.*ECI
         CONS14=GAMMA(BG+3.)*PI/4.*ECI
         CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
         CONS16=GAMMA(BI+3.)*PI/4.*ECI
         CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
         CONS18=RHOSN*RHOSN
         CONS19=RHOW*RHOW
         CONS20=20.*PI*PI*RHOW*BIMM
         CONS21=4./(DCS*RHOI)
         CONS22=PI*RHOI*DCS**3/6.
         CONS23=PI/4.*EII*GAMMA(BS+3.)
         CONS24=PI/4.*ECR*GAMMA(BR+3.)
         CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
         CONS26=PI/6.*RHOW
         CONS27=GAMMA(1.+BI)
         CONS28=GAMMA(4.+BI)/6.
         CONS29=4./3.*PI*RHOW*(25.E-6)**3
         CONS30=4./3.*PI*RHOW
         CONS31=PI*PI*ECR*RHOSN
         CONS32=PI/2.*ECR
         CONS33=PI*PI*ECR*RHOG
         CONS34=5./2.+BR/2.
         CONS35=5./2.+BS/2.
         CONS36=5./2.+BG/2.
         CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
         CONS38=PI*PI/3.*RHOW
         CONS39=PI*PI/36.*RHOW*BIMM
         CONS40=PI/6.*BIMM
         CONS41=PI*PI*ECR*RHOW

END SUBROUTINE MORR_TWO_MOMENT_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SET_MORRISON_NDCNST(ndcnst_in) BIND(C, name='set_morrison_ndcnst_c')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE SETS THE CONSTANT DROPLET CONCENTRATION (NDCNST)
! ALLOWS RUNTIME CONFIGURATION FROM INPUT FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: ndcnst_in

      NDCNST = ndcnst_in

END SUBROUTINE SET_MORRISON_NDCNST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE MORR_TWO_MOMENT_MICRO)
! WHICH OPERATES ON 1D VERTICAL COLUMNS.
! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE.
! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK TO DRIVER MODEL.

! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND SLAVA TATARSKII (GEORGIA TECH).

! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MP_MORR_TWO_MOMENT(ITIMESTEP,                       &
                TH, QV, QC, QR, QI, QS, QG, NI, NS, NR, NG, &
                RHO, PII, P, DT_IN, DZ, HT, W,          &
                RAINNC, RAINNCV, SR,                    &
                SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,    & ! hm added 7/13/13
                refl_10cm, diagflag, do_radar_ref,      & ! GT added for reflectivity calcs
                qrcuten, qscuten, qicuten               & ! hm added
               ,F_QNDROP, qndrop                        & ! hm added, wrf-chem
               ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
               ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
               ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
               ,wetscav_on, rainprod, evapprod                      &
               ,QLSINK,PRECR,PRECI,PRECS,PRECG)

! QV - water vapor mixing ratio (kg/kg)
! QC - cloud water mixing ratio (kg/kg)
! QR - rain water mixing ratio (kg/kg)
! QI - cloud ice mixing ratio (kg/kg)
! QS - snow mixing ratio (kg/kg)
! QG - graupel mixing ratio (KG/KG)
! NI - cloud ice number concentration (1/kg)
! NS - Snow Number concentration (1/kg)
! NR - Rain Number concentration (1/kg)
! NG - Graupel number concentration (1/kg)
! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
! P - AIR PRESSURE (PA)
! W - VERTICAL AIR VELOCITY (M/S)
! TH - POTENTIAL TEMPERATURE (K)
! PII - exner function - used to convert potential temp to temp
! DZ - difference in height over interface (m)
! DT_IN - model time step (sec)
! ITIMESTEP - time step counter
! RAINNC - accumulated grid-scale precipitation (mm)
! RAINNCV - one time step grid scale precipitation (mm/time step)
! SNOWNC - accumulated grid-scale snow plus cloud ice (mm)
! SNOWNCV - one time step grid scale snow plus cloud ice (mm/time step)
! GRAUPELNC - accumulated grid-scale graupel (mm)
! GRAUPELNCV - one time step grid scale graupel (mm/time step)
! SR - one time step mass ratio of snow to total precip
! qrcuten, rain tendency from parameterized cumulus convection
! qscuten, snow tendency from parameterized cumulus convection
! qicuten, cloud ice tendency from parameterized cumulus convection

! variables below currently not in use, not coupled to PBL or radiation codes
! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! HM, ADDED FOR WRF-CHEM COUPLING
! QLSINK - TENDENCY OF CLOUD WATER TO RAIN, SNOW, GRAUPEL (KG/KG/S)
! CSED,ISED,SSED,GSED,RSED - SEDIMENTATION FLUXES (KG/M^2/S) FOR CLOUD WATER, ICE, SNOW, GRAUPEL, RAIN
! PRECI,PRECS,PRECG,PRECR - SEDIMENTATION FLUXES (KG/M^2/S) FOR ICE, SNOW, GRAUPEL, RAIN

! rainprod - total tendency of conversion of cloud water/ice and graupel to rain (kg kg-1 s-1)
! evapprod - tendency of evaporation of rain (kg kg-1 s-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! reflectivity currently not included!!!!
! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                                       ims, ime, jms, jme, kms, kme , &
                                       its, ite, jts, jte, kts, kte

   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme, kms:kme),          INTENT(INOUT) :: qv, qc, qr, qi, qs, qg, ni, ns, nr, TH, NG
   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme, kms:kme), optional,INTENT(INOUT) :: qndrop
   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme, kms:kme), optional,INTENT(INOUT) :: QLSINK, rainprod, evapprod, PRECI,PRECS,PRECG, &
     & PRECR
   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme, kms:kme),          INTENT(IN   ) :: pii, p, dz, rho, w

   REAL(C_DOUBLE), INTENT(IN):: dt_in
   INTEGER, INTENT(IN):: ITIMESTEP
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme, kms:kme) :: rainnc, snownc, graupelnc
    REAL(C_DOUBLE), INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: rainncv, sr,  snowncv, graupelncv
!   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme         ), INTENT(INOUT) :: RAINNC, RAINNCV, SR, SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV
   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:kme, kms:jme), INTENT(INOUT) :: refl_10cm
   REAL(C_DOUBLE), DIMENSION(ims:ime ,jms:jme         ), INTENT(IN   ) :: ht

   LOGICAL, optional, INTENT(IN) :: wetscav_on

   ! LOCAL VARIABLES

   REAL(C_DOUBLE), DIMENSION(its:ite, jts:jte, kts:kte) :: effi, effs, effr, EFFG
   REAL(C_DOUBLE), DIMENSION(its:ite, jts:jte, kts:kte) :: T, EFFC

   REAL(C_DOUBLE), DIMENSION(kts:kte) ::                                 &
                            QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, &
                            NI_TEND1D, NS_TEND1D, NR_TEND1D,             &
                            QC1D, QI1D, QR1D,NI1D, NS1D, NR1D, QS1D,     &
                            T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, W1D,     &
                            EFFC1D, EFFI1D, EFFS1D, EFFR1D,DZ1D,         &
                            QG_TEND1D, NG_TEND1D, QG1D, NG1D, EFFG1D,    &
! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
                            QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN,      &
! ADD CUMULUS TENDENCIES
                            QRCU1D, QSCU1D, QICU1D
! add cumulus tendencies

   REAL(C_DOUBLE), DIMENSION(ims:ime, jms:jme, kms:kme), INTENT(IN):: qrcuten, qscuten, qicuten

  LOGICAL, INTENT(IN), OPTIONAL :: F_QNDROP
  LOGICAL :: flag_qndrop  ! wrf-chem
  integer :: iinum ! wrf-chem

! wrf-chem
   REAL(C_DOUBLE), DIMENSION(kts:kte) :: nc1d, nc_tend1d,C2PREC,CSED,ISED,SSED,GSED,RSED
   REAL(C_DOUBLE), DIMENSION(kts:kte) :: rainprod1d, evapprod1d
! HM add reflectivity
   REAL(C_DOUBLE), DIMENSION(kts:kte) :: dBZ

   REAL(C_DOUBLE) PRECPRT1D, SNOWRT1D, SNOWPRT1D, GRPLPRT1D ! hm added 7/13/13

   INTEGER I,J,K

   REAL(C_DOUBLE) DT

   LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
   INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref

   LOGICAL :: has_wetscav

   ! return
   print *,'IN MORR_TWO_MOMENT MAIN ROUTINE '

! below for wrf-chem
   flag_qndrop = .false.
   IF ( PRESENT ( f_qndrop ) ) flag_qndrop = f_qndrop
!!!!!!!!!!!!!!!!!!!!!!

   IF( PRESENT( wetscav_on ) ) THEN
     has_wetscav = wetscav_on
   ELSE
     has_wetscav = .false.
   ENDIF

   ! Initialize tendencies (all set to 0) and transfer
   ! array to local variables
   DT = DT_IN

   ! Currently mixing of number concentrations also is neglected (not coupled with PBL schemes)

   print *,'IMS IME ',ims,ime
   print *,'JMS JME ',jms,jme
   print *,'KMS KME ',kms,kme

   print *,'ITS ITE ',its,ite
   print *,'JTS JTE ',jts,jte
   print *,'KTS KTE ',kts,kte

   DO I=ITS,ITE
   DO J=JTS,JTE
   DO K=KTS,KTE
        T(i,j,k) = TH(i,j,k) * PII(i,j,k)
   END DO
   END DO
   END DO

   do i=its,ite      ! i loop (east-west)
   do j=jts,jte      ! j loop (north-south)
   !
   ! Transfer 3D arrays into 1D for microphysical calculations
   !

! hm , initialize 1d tendency arrays to zero

      do k=kts,kte   ! k loop (vertical)

          QC_TEND1D(k)  = 0.
          QI_TEND1D(k)  = 0.
          QNI_TEND1D(k) = 0.
          QR_TEND1D(k)  = 0.
          NI_TEND1D(k)  = 0.
          NS_TEND1D(k)  = 0.
          NR_TEND1D(k)  = 0.
          T_TEND1D(k)   = 0.
          QV_TEND1D(k)  = 0.
          nc_tend1d(k) = 0. ! wrf-chem

          QC1D(k)       = QC(i,j,k)
          QI1D(k)       = QI(i,j,k)
          QS1D(k)       = QS(i,j,k)
          QR1D(k)       = QR(i,j,k)

          NI1D(k)       = NI(i,j,k)

          NS1D(k)       = NS(i,j,k)
          NR1D(k)       = NR(i,j,k)
! HM ADD GRAUPEL
          QG1D(K)       = QG(I,j,k)
          NG1D(K)       = NG(I,j,k)
          QG_TEND1D(K)  = 0.
          NG_TEND1D(K)  = 0.

          T1D(k)        = T(i,j,k)
          QV1D(k)       = QV(i,j,k)
          P1D(k)        = P(i,j,k)
          DZ1D(k)       = DZ(i,j,k)
          W1D(k)        = W(i,j,k)
! add cumulus tendencies, already decoupled
          qrcu1d(k)     = qrcuten(i,j,k)
          qscu1d(k)     = qscuten(i,j,k)
          qicu1d(k)     = qicuten(i,j,k)
      end do  !jdf added this

! below for wrf-chem
   IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
      iact = 3
      DO k = kts, kte
         nc1d(k)=qndrop(i,j,k)
         iinum=0
      ENDDO
   ELSE
      DO k = kts, kte
         nc1d(k)=0. ! temporary placeholder, set to constant in microphysics subroutine
         iinum=1
      ENDDO
   ENDIF

!jdf  end do

      call MORR_TWO_MOMENT_MICRO(i,j,ITIMESTEP,kts,kte,QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D,   &
                                 NI_TEND1D, NS_TEND1D, NR_TEND1D,               &
                                 QC1D, QI1D, QS1D, QR1D,NI1D, NS1D, NR1D,       &
                                 T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, DZ1D, W1D, &
                                 PRECPRT1D,SNOWRT1D,                            &
                                 SNOWPRT1D,GRPLPRT1D,                 & ! hm added 7/13/13
                                 EFFC1D,EFFI1D,EFFS1D,EFFR1D,DT,                &
                                 QG_TEND1D,NG_TEND1D,QG1D,NG1D,EFFG1D,          &
                                 qrcu1d, qscu1d, qicu1d,                        &
! ADD SEDIMENTATION TENDENCIES
                                 QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, &
                                 nc1d, nc_tend1d, iinum, C2PREC,CSED,ISED,SSED,GSED,RSED)

   !
   ! Transfer 1D arrays back into 3D arrays
   !
      do k=kts,kte

         ! hm, add tendencies to update global variables
         ! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
         ! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

          QC(i,j,k)        = QC1D(k)
          QI(i,j,k)        = QI1D(k)
          QS(i,j,k)        = QS1D(k)
          QR(i,j,k)        = QR1D(k)
          NI(i,j,k)        = NI1D(k)
          NS(i,j,k)        = NS1D(k)
          NR(i,j,k)        = NR1D(k)
          QG(I,j,k)        = QG1D(K)
          NG(I,j,k)        = NG1D(K)

          T(i,j,k)         = T1D(k)
          TH(I,j,k)        = T(i,j,k)/PII(i,j,k) ! CONVERT TEMP BACK TO POTENTIAL TEMP
          QV(i,j,k)        = QV1D(k)

          EFFC(i,j,k)      = EFFC1D(k)
          EFFI(i,j,k)      = EFFI1D(k)
          EFFS(i,j,k)      = EFFS1D(k)
          EFFR(i,j,k)      = EFFR1D(k)
          EFFG(I,j,k)      = EFFG1D(K)

! wrf-chem
          IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
             qndrop(i,j,k) = nc1d(k)
!jdf         CSED3D(i,j,k) = CSED(k)
          END IF
          IF ( PRESENT( QLSINK ) ) THEN
             if(qc(i,j,k)>1.e-10) then
                QLSINK(I,J,K)  = C2PREC(K)/QC(I,J,K)
             else
                QLSINK(I,J,K)  = 0.0
             endif
          END IF
          IF ( PRESENT( PRECR ) ) PRECR(I,J,K) = RSED(K)
          IF ( PRESENT( PRECI ) ) PRECI(I,J,K) = ISED(K)
          IF ( PRESENT( PRECS ) ) PRECS(I,J,K) = SSED(K)
          IF ( PRESENT( PRECG ) ) PRECG(I,J,K) = GSED(K)
! EFFECTIVE RADIUS FOR RADIATION CODE (currently not coupled)
! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
!          EFFCS(I,J,K)     = MIN(EFFC(I,J,K),50.)
!          EFFCS(I,J,K)     = MAX(EFFCS(I,J,K),1.)
!          EFFIS(I,J,K)     = MIN(EFFI(I,J,K),130.)
!          EFFIS(I,J,K)     = MAX(EFFIS(I,J,K),13.)

      end do

! hm modified so that m2005 precip variables correctly match wrf precip variables
      RAINNC(i,j,KTS)     = RAINNC(I,J,KTS)+PRECPRT1D
      RAINNCV(i,j)    = PRECPRT1D
      SNOWNC(i,j,KTS)     = SNOWNC(I,J,KTS)+SNOWPRT1D
      SNOWNCV(i,j)    = SNOWPRT1D
      GRAUPELNC(i,j,KTS)  = GRAUPELNC(I,J,KTS)+GRPLPRT1D
      GRAUPELNCV(i,j) = GRPLPRT1D
      SR(i,j) = SNOWRT1D/(PRECPRT1D+1.E-12)

   end do
   end do

END SUBROUTINE MP_MORR_TWO_MOMENT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MORR_TWO_MOMENT_MICRO(i,j,istep,kts,kte,                                                   &
                                       QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,                              &
                                       NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NI3D,NS3D,NR3D,   &
                                       T3DTEN,QV3DTEN,T3D,QV3D,PRES,DZQ,W3D,PRECRT,SNOWRT,            &
                                       SNOWPRT,GRPLPRT,                                               &
                                       EFFC,EFFI,EFFS,EFFR,DT,                                        &
                                       QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,qrcu1d,qscu1d, qicu1d,          &
                                       QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN,                           &
                                       nc3d,nc3dten,iinum,                                            & ! wrf-chem
                                       c2prec,CSED,ISED,SSED,GSED,RSED)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
! MORRISON ET AL. 2005 JAS AND MORRISON ET AL. 2009 MWR

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

! CODE STRUCTURE: MAIN SUBROUTINE IS 'MORR_TWO_MOMENT'. ALSO INCLUDED IN THIS FILE IS
! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
! 'FUNCTION GAMMA'.

! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! DECLARATIONS

      IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
! DEFINE ARRAY SIZES

! INPUT NUMBER OF GRID CELLS

! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
      INTEGER, INTENT( IN)  :: i,j,istep,kts,kte

      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  T3D                ! TEMPERATURE (K)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  DZQ                ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
! below for wrf-chem
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  nc3d
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  nc3dten
      integer, intent(in) :: iinum

! HM ADDED GRAUPEL VARIABLES
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)

! hm add cumulus tendencies for precip
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   qrcu1d
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   qscu1d
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   qicu1d

! OUTPUT VARIABLES

        REAL(C_DOUBLE) PRECRT                ! TOTAL PRECIP PER TIME STEP (mm)
        REAL(C_DOUBLE) SNOWRT                ! SNOW PER TIME STEP (mm)
! hm added 7/13/13
        REAL(C_DOUBLE) SNOWPRT      ! TOTAL CLOUD ICE PLUS SNOW PER TIME STEP (mm)
        REAL(C_DOUBLE) GRPLPRT      ! TOTAL GRAUPEL PER TIME STEP (mm)

        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
        REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)

! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

        REAL(C_DOUBLE) DT         ! MODEL TIME STEP (SEC)


!.....................................................................................................
! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
! REST OF THE MODEL.

! SIZE PARAMETER VARIABLES

     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: LAMC          ! SLOPE PARAMETER FOR DROPLETS (M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: LAMI          ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: LAMS          ! SLOPE PARAMETER FOR SNOW (M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: LAMR          ! SLOPE PARAMETER FOR RAIN (M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: LAMG          ! SLOPE PARAMETER FOR GRAUPEL (M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: CDIST1        ! PSD PARAMETER FOR DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: N0I           ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: N0S           ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: N0RR          ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: N0G           ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) :: PGAM          ! SPECTRAL SHAPE PARAMETER FOR DROPLETS

! MICROPHYSICAL PROCESSES

     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSUBC     ! LOSS OF NC DURING EVAP
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSUBI     ! LOSS OF NI DURING SUB.
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSUBS     ! LOSS OF NS DURING SUB.
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSUBR     ! LOSS OF NR DURING EVAP
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRD       ! DEP CLOUD ICE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRE       ! EVAP OF RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRDS      ! DEP SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRA       ! ACCRETION DROPLETS BY RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRC       ! AUTOCONVERSION DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PCC       ! COND/EVAP DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSAGG     ! SELF-COLLECTION OF SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRCI      ! CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSMLTS    ! CHANGE N MELTING SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
! HM ADDED 12/13/06
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  EPRD      ! SUBLIMATION CLOUD ICE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  EPRDS     ! SUBLIMATION SNOW
! HM ADDED GRAUPEL PROCESSES
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PRDG    ! DEP OF GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  EPRDG    ! SUB OF GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL

! TIME-VARYING ATMOSPHERIC PARAMETERS

     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EVS   ! SATURATION VAPOR PRESSURE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   EIS   ! ICE SATURATION VAPOR PRESSURE
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   QVS   ! SATURATION MIXING RATIO
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   QVI   ! ICE SATURATION MIXING RATIO
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   QVQVS ! SAUTRATION RATIO
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   QVQVSI! ICE SATURAION RATIO
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   MU    ! VISCOCITY OF AIR
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   SC    ! SCHMIDT NUMBER
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   XLF   ! LATENT HEAT OF FREEZING
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   RHO   ! AIR DENSITY
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   AB    ! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   ABI    ! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

! TIME-VARYING MICROPHYSICS PARAMETERS

     REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   DAP    ! DIFFUSIVITY OF AEROSOL
     REAL(C_DOUBLE)    NACNT                    ! NUMBER OF CONTACT IN
     REAL(C_DOUBLE)    FMULT                    ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
     REAL(C_DOUBLE)    COFFI                    ! ICE AUTOCONVERSION PARAMETER

! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
      REAL(C_DOUBLE) UNI, UMI,UMR
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::    FR, FI, FNI,FG,FNG
      REAL(C_DOUBLE) RGVM
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   FALOUTR,FALOUTI,FALOUTNI
      REAL(C_DOUBLE) FALTNDR,FALTNDI,FALTNDNI,RHO2
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   DUMQS,DUMFNS
      REAL(C_DOUBLE) UMS,UNS
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   FS,FNS, FALOUTS,FALOUTNS,FALOUTG,FALOUTNG
      REAL(C_DOUBLE) FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::    DUMC,DUMFNC
      REAL(C_DOUBLE) UNC,UMC,UNG,UMG
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   FC,FALOUTC,FALOUTNC
      REAL(C_DOUBLE) FALTNDC,FALTNDNC
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   FNC,DUMFNR,FALOUTNR
      REAL(C_DOUBLE) FALTNDNR
      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::   FNR

! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

      REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::    AIN,ARN,ASN,ACN,AGN

! EXTERNAL FUNCTION CALL RETURN VARIABLES

!      REAL(C_DOUBLE) GAMMA,      ! EULER GAMMA FUNCTION
!      REAL(C_DOUBLE) POLYSVP,    ! SAT. PRESSURE FUNCTION
!      REAL(C_DOUBLE) DERF1        ! ERROR FUNCTION

! DUMMY VARIABLES

     REAL(C_DOUBLE) DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS

! PROGNOSTIC SUPERSATURATION

     REAL(C_DOUBLE) DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
     REAL(C_DOUBLE) DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
     REAL(C_DOUBLE) EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
     REAL(C_DOUBLE) EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
     REAL(C_DOUBLE) EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
     REAL(C_DOUBLE) EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

! NEW DROPLET ACTIVATION VARIABLES
     REAL(C_DOUBLE) TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
     REAL(C_DOUBLE) TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
     REAL(C_DOUBLE) TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
     REAL(C_DOUBLE) TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
     REAL(C_DOUBLE) TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
     REAL(C_DOUBLE) DUMACT,DUM3

! COUNTING/INDEX VARIABLES

     INTEGER K,NSTEP,N ! ,I

! LTRUE IS ONLY USED TO SPEED UP THE CODE !!
! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN,
!               = 1, HYDROMETEORS IN COLUMN

      INTEGER LTRUE

! DROPLET ACTIVATION/FREEZING AEROSOL


     REAL(C_DOUBLE)    CT      ! DROPLET ACTIVATION PARAMETER
     REAL(C_DOUBLE)    TEMP1   ! DUMMY TEMPERATURE
     REAL(C_DOUBLE)    SAT1    ! DUMMY SATURATION
     REAL(C_DOUBLE)    SIGVL   ! SURFACE TENSION LIQ/VAPOR
     REAL(C_DOUBLE)    KEL     ! KELVIN PARAMETER
     REAL(C_DOUBLE)    KC2     ! TOTAL ICE NUCLEATION RATE

       REAL(C_DOUBLE) CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

! MORE WORKING/DUMMY VARIABLES

     REAL(C_DOUBLE) DUMQI,DUMNI,DC0,DS0,DG0
     REAL(C_DOUBLE) DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

! EFFECTIVE VERTICAL VELOCITY  (M/S)
     REAL(C_DOUBLE) WEF

! WORKING PARAMETERS FOR ICE NUCLEATION

      REAL(C_DOUBLE) ANUC,BNUC

! WORKING PARAMETERS FOR AEROSOL ACTIVATION

        REAL(C_DOUBLE) AACT,GAMM,GG,PSI,ETA1,ETA2,SM1,SM2,SMAX,UU1,UU2,ALPHA

! DUMMY SIZE DISTRIBUTION PARAMETERS

        REAL(C_DOUBLE) DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN

        INTEGER IDROP

! FOR WRF-CHEM
        REAL(C_DOUBLE), DIMENSION(KTS:KTE)::C2PREC,CSED,ISED,SSED,GSED,RSED
    REAL(C_DOUBLE), DIMENSION(KTS:KTE)                :: tqimelt ! melting of cloud ice (tendency)

! comment lines for wrf-chem since these are intent(in) in that case
!       REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
!       REAL(C_DOUBLE), DIMENSION(KTS:KTE) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        ! print *,'IN MICRO KTS:KTE ', i,j,kts, kte

! SET LTRUE INITIALLY TO 0

         LTRUE = 0

! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
         DO K = KTS,KTE

! NC3DTEN LOCAL ARRAY INITIALIZED
               NC3DTEN(K) = 0.
! INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO

                C2PREC(K)=0.
                CSED(K)=0.
                ISED(K)=0.
                SSED(K)=0.
                GSED(K)=0.
                RSED(K)=0.

! LATENT HEAT OF VAPORATION

            XXLV(K) = 3.1484E6-2370.*T3D(K)

! LATENT HEAT OF SUBLIMATION

            XXLS(K) = 3.15E6-2370.*T3D(K)+0.3337E6

            CPM(K) = CP*(1.+0.887*QV3D(K))

! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10

            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA
! Fortran version
! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AIR DENSITY

            RHO(K) = PRES(K)/(R*T3D(K))

! ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
! FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON

            IF (QRCU1D(K).GE.1.E-10) THEN
            DUM=1.8e5*(QRCU1D(K)*DT/(PI*RHOW*RHO(K)**3))**0.25
            NR3D(K)=NR3D(K)+DUM
            END IF
            IF (QSCU1D(K).GE.1.E-10) THEN
            DUM=3.e5*(QSCU1D(K)*DT/(CONS1*RHO(K)**3))**(1./(DS+1.))
            NS3D(K)=NS3D(K)+DUM
            END IF
            IF (QICU1D(K).GE.1.E-10) THEN
            DUM=QICU1D(K)*DT/(CI*(80.E-6)**DI)
            NI3D(K)=NI3D(K)+DUM
            END IF

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! hm modify 7/0/09 change limit to 1.e-8

             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
               END IF
             END IF
! Fortran version
             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
            END IF
            ! Fortran version
! HEAT OF FUSION

            XLF(K) = XXLS(K)-XXLV(K)

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
         QC3D(K) = 0.
         NC3D(K) = 0.
         EFFC(K) = 0.
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
         QR3D(K) = 0.
         NR3D(K) = 0.
         EFFR(K) = 0.
       END IF
       IF (QI3D(K).LT.QSMALL) THEN
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       END IF
       IF (QNI3D(K).LT.QSMALL) THEN
         QNI3D(K) = 0.
         NS3D(K) = 0.
         EFFS(K) = 0.
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
         QG3D(K) = 0.
         NG3D(K) = 0.
         EFFG(K) = 0.
       END IF
! Fortran version
! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

      QRSTEN(K) = 0.
      QISTEN(K) = 0.
      QNISTEN(K) = 0.
      QCSTEN(K) = 0.
      QGSTEN(K) = 0.

!..................................................................
! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

! fix 053011
            MU(K) = 1.496E-6*T3D(K)**1.5/(T3D(K)+120.)

! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

            DUM = (RHOSU/RHO(K))**0.54

! fix 053011
!            AIN(K) = DUM*AI
! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction
            AIN(K) = (RHOSU/RHO(K))**0.35*AI
            ARN(K) = DUM*AR
            ASN(K) = DUM*AS
!            ACN(K) = DUM*AC
! AA revision 4/1/11: temperature-dependent Stokes fall speed
            ACN(K) = G*RHOW/(18.*MU(K))
! HM ADD GRAUPEL 8/28/06
            AGN(K) = DUM*AG
!hm 4/7/09 bug fix, initialize lami to prevent later division by zero
            LAMI(K)=0.

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
! FOR THIS LEVEL

            IF ( QC3D(K).LT.QSMALL.AND. &
                 QI3D(K).LT.QSMALL.AND. &
                QNI3D(K).LT.QSMALL.AND. &
                 QR3D(K).LT.QSMALL.AND. &
                 QG3D(K).LT.QSMALL) THEN
            IF (T3D(K).LT.273.15.AND.QVQVSI(K).LT.0.999) then
               GOTO 200
            endif
            IF (T3D(K).GE.273.15.AND.QVQVS(K).LT.0.999) then
               GOTO 200
            endif
            END IF

! THERMAL CONDUCTIVITY FOR AIR

! fix 053011
            KAP(K) = 1.414E3*MU(K)

! DIFFUSIVITY OF WATER VAPOR

            DV(K) = 8.794E-5*T3D(K)**1.81/PRES(K)

! SCHMIT NUMBER

! fix 053011
            SC(K) = MU(K)/(RHO(K)*DV(K))

! PSYCHOMETIC CORRECTIONS

! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

            DUM = (RV*T3D(K)**2)

            DQSDT = XXLV(K)*QVS(K)/DUM
            DQSIDT =  XXLS(K)*QVI(K)/DUM

            ABI(K) = 1.+DQSIDT*XXLS(K)/CPM(K)
            AB(K) = 1.+DQSDT*XXLV(K)/CPM(K)
!
!.....................................................................
!.....................................................................
! CASE FOR TEMPERATURE ABOVE FREEZING

            IF (T3D(K).GE.273.15) THEN

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (iinum.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         END IF

! GET SIZE DISTRIBUTION PARAMETERS

! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
       IF (QNI3D(K).LT.1.E-6) THEN
          QR3D(K)=QR3D(K)+QNI3D(K)
          NR3D(K)=NR3D(K)+NS3D(K)
          T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
          QNI3D(K) = 0.
          NS3D(K) = 0.
       END IF
       IF (QG3D(K).LT.1.E-6) THEN
          QR3D(K)=QR3D(K)+QG3D(K)
          NR3D(K)=NR3D(K)+NG3D(K)
          T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
          QG3D(K) = 0.
          NG3D(K) = 0.
       END IF
       IF (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8) THEN
          GOTO 300
       ENDIF

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
      END IF
!.....................................................................
! ZERO OUT PROCESS RATES

            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            NSMLTS(K) = 0.
            NSMLTR(K) = 0.
            EVPMS(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            NSUBC(K) = 0.
            NSUBR(K) = 0.
            PRACG(K) = 0.
            NPRACG(K) = 0.
            PSMLT(K) = 0.
            PGMLT(K) = 0.
            EVPMG(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

         END IF

!.......................................................................
! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
! FORMULA FROM IKAWA AND SAITO (1991)

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! hm fix, 2/12/13
! for above freezing conditions to get accelerated melting of snow,
! we need collection of rain by snow (following Lin et al. 1983)
!            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
!                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
!                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
!                  (5./(LAMS(K)**3*LAMR(K))+                    &
!                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
!                  0.5/(LAMS(K)*LAMR(K)**3)))

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &
                  0.5/(LAMR(k)*LAMS(k)**3)))

! fix 053011, npracs no longer subtracted from snow
!            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
!                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
!                (1./(LAMR(K)**3*LAMS(K))+                      &
!                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
!                 1./(LAMR(K)*LAMS(K)**3))

         END IF

! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
! ASSUME SHED DROPS ARE 1 MM IN SIZE

         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+                              &
                                  0.5/(LAMR(k)*LAMG(k)**3)))

! ASSUME 1 MM DROPS ARE SHED, GET NUMBER SHED PER SEC

            DUM = PRACG(K)/5.2E-7

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! hm 7/15/13, remove limit so that the number of collected drops can smaller than
! number of shed drops
!            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
            NPRACG(K)=NPRACG(K)-DUM

            END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

         IF (QR3D(K).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
!            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
            NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF
! NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED

           IF (QV3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
           ELSE
              PRE(K) = 0.
           END IF
!.......................................................................
! MELTING OF SNOW

! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

          IF (QNI3D(K).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACS(K)

! hm fix 1/20/15
!             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
!                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
!                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS10/                   &
!                   (LAMS(K)**CONS35))+DUM
             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
                    XLF(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
                   (LAMS(K)**CONS35))+DUM

! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

      IF (QVQVS(K).LT.1.) THEN
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
! hm fix 8/4/08
        EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)
        EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
        PSMLT(K) = PSMLT(K)-EVPMS(K)
      END IF
      END IF

!.......................................................................
! MELTING OF GRAUPEL

! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

          IF (QG3D(K).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACG(K)

! hm fix 1/20/15
!             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/              &
!                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
!                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
!                    SC(K)**(1./3.)*CONS11/                   &
!                   (LAMG(K)**CONS36))+DUM
             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/               &
                    XLF(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
                   (LAMG(K)**CONS36))+DUM

! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

      IF (QVQVS(K).LT.1.) THEN
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))
! hm fix 8/4/08
        EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
        EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
        PGMLT(K) = PGMLT(K)-EVPMG(K)
      END IF
      END IF

! HM, V3.2
! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

      PRACG(K) = 0.
      PRACS(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
! CALCULATION

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K))*DT

      IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN

        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO

        END IF

! CONSERVATION OF SNOW

        DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

        IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

! NO SOURCE TERMS FOR SNOW AT T > FREEZING
        RATIO = QNI3D(K)/DUM

        PSMLT(K) = PSMLT(K)*RATIO
        EVPMS(K) = EVPMS(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO

        END IF

! CONSERVATION OF GRAUPEL

        DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT

        IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
        RATIO = QG3D(K)/DUM

        PGMLT(K) = PGMLT(K)*RATIO
        EVPMG(K) = EVPMG(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        END IF

! CONSERVATION OF QR
! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

        DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT

        IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

        RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                        (-PRE(K))
        PRE(K) = PRE(K)*RATIO
        END IF

!....................................

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

      T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)

      QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
      QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
      QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
! fix 053011
!      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
! HM, bug fix 5/12/08, npracg is subtracted from nr not ng
!      NG3DTEN(K) = NG3DTEN(K)
      NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
      NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))
! Fortran version
! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

        C2PREC(K) = PRA(K)+PRC(K)
      IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      END IF

        IF (EVPMS(K)+PSMLT(K).LT.0.) THEN
         DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSMLTS(K) = DUM*NS3D(K)/DT
        END IF
        IF (PSMLT(K).LT.0.) THEN
          DUM = PSMLT(K)*DT/QNI3D(K)
          DUM = MAX(-1.0,DUM)
          NSMLTR(K) = DUM*NS3D(K)/DT
        END IF
        IF (EVPMG(K)+PGMLT(K).LT.0.) THEN
         DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NGMLTG(K) = DUM*NG3D(K)/DT
        END IF
        IF (PGMLT(K).LT.0.) THEN
          DUM = PGMLT(K)*DT/QG3D(K)
          DUM = MAX(-1.0,DUM)
          NGMLTR(K) = DUM*NG3D(K)/DT
        END IF

         NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
         NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
         NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))

 300  CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K)+DT*QV3DTEN(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = EP_2*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)+DT*QC3DTEN(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
           PCC(K) = -DUMQC/DT
      END IF

      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
!.......................................................................
! ACTIVATION OF CLOUD DROPLETS
! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
! DROPLET CONCENTRATION IS SPECIFIED !!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

!.....................................................................
!.....................................................................
         ELSE  ! TEMPERATURE < 273.15

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (iinum.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         END IF

! CALCULATE SIZE DISTRIBUTION PARAMETERS
! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)
         N0I(K) = NI3D(K)*LAMI(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF
!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

! TO CALCULATE DROPLET FREEZING

        CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
   END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
      END IF
!.....................................................................
! ZERO OUT PROCESS RATES

            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
            PSACWG(K) = 0.
            PGSACW(K) = 0.
            PGRACS(K) = 0.
            PRDG(K) = 0.
            EPRDG(K) = 0.
            NPRACG(K) = 0.
            NPSACWG(K) = 0.
            NSCNG(K) = 0.
            NGRACS(K) = 0.
            NSUBG(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES
! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
!.......................................................................
! FREEZING OF CLOUD DROPLETS
! ONLY ALLOWED BELOW -4 C
        IF (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15) THEN

! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

! MEYERS CURVE

           NACNT = EXP(-2.80+0.262*(273.15-T3D(K)))*1000.

! COOPER CURVE
!        NACNT =  5.*EXP(0.304*(273.15-T3D(K)))

! FLECTHER
!     NACNT = 0.01*EXP(0.6*(273.15-T3D(K)))

! CONTACT FREEZING

! MEAN FREE PATH

            DUM = 7.37*T3D(K)/(288.*10.*PRES(K))/100.

! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
! BASED ON BROWNIAN DIFFUSION

            DAP(K) = CONS37*T3D(K)*(1.+DUM/RIN)/MU(K)

           MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                   LOG(GAMMA(PGAM(K)+5.))-4.*LOG(LAMC(K)))
           NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
                    GAMMA(PGAM(K)+2.)/                         &
                    LAMC(K)

! IMMERSION FREEZING (BIGG 1953)

!           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
!                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
!                   EXP(AIMM*(273.15-T3D(K)))

!           NNUCCC(K) = NNUCCC(K)+                                  &
!            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
!                *EXP(AIMM*(273.15-T3D(K)))

! hm 7/15/13 fix for consistency w/ original formula
           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                   (EXP(AIMM*(273.15-T3D(K)))-1.)

           NNUCCC(K) = NNUCCC(K)+                                  &
            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                *(EXP(AIMM*(273.15-T3D(K)))-1.)

! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
! MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC

           NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

        END IF

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

         END IF
!.......................................................................
! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

         IF (QNI3D(K).GE.1.E-8) THEN
             NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
            ((2.+BS)/3.)*QNI3D(K)**((2.+BS)/3.)*                  &
            (NS3D(K)*RHO(K))**((4.-BS)/3.)/                       &
            (RHO(K))
         END IF

!.......................................................................
! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
! HERE USE CONTINUOUS COLLECTION EQUATION WITH
! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

! SNOW

         IF (QNI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

           PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)
           NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)

         END IF

!............................................................................
! COLLECTION OF CLOUD WATER BY GRAUPEL

         IF (QG3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

           PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
           NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
            END IF
!.......................................................................
! HM, ADD 12/13/06
! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
! BEFORE RIMING CAN OCCUR
! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
! TO HALLET-MOSSOP SPLINTERING

         IF (QI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
! FROM THOMPSON ET AL. 2004, MWR

            IF (1./LAMI(K).GE.100.E-6) THEN

           PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           END IF
         END IF

!.......................................................................
! ACCRETION OF RAIN WATER BY SNOW
! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &
                  0.5/(LAMR(k)*LAMS(k)**3)))

            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                (1./(LAMR(K)**3*LAMS(K))+                      &
                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
                 1./(LAMR(K)*LAMS(K)**3))

! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

! HM MODIFY FOR WRFV3.1
!            IF (IHAIL.EQ.0) THEN
            IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
            PSACR(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                 N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                  (5./(LAMS(K)**3*LAMR(K))+                    &
                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                  0.5/(LAMS(K)*LAMR(K)**3)))
            END IF
!            END IF

         END IF

!.......................................................................

! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990,
! USED BY REISNER ET AL 1998
         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+                              &
                                  0.5/(LAMR(k)*LAMG(k)**3)))

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

            END IF

!.......................................................................
! RIME-SPLINTERING - SNOW
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

!v1.4
         IF (QNI3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWS(K).GT.0..OR.PRACS(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

               IF (PSACWS(K).GT.0.) THEN
                  NMULTS(K) = 35.E4*PSACWS(K)*FMULT*1000.
                  QMULTS(K) = NMULTS(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                  QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                  PSACWS(K) = PSACWS(K)-QMULTS(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACS(K).GT.0.) THEN
                   NMULTR(K) = 35.E4*PRACS(K)*FMULT*1000.
                   QMULTR(K) = NMULTR(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                   QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                   PRACS(K) = PRACS(K)-QMULTR(K)

               END IF

            END IF
         END IF
         END IF
         END IF

!.......................................................................
! RIME-SPLINTERING - GRAUPEL
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

!         IF (IHAIL.EQ.0) THEN
! v1.4
         IF (QG3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWG(K).GT.0..OR.PRACG(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

               IF (PSACWG(K).GT.0.) THEN
                  NMULTG(K) = 35.E4*PSACWG(K)*FMULT*1000.
                  QMULTG(K) = NMULTG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                  QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                  PSACWG(K) = PSACWG(K)-QMULTG(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACG(K).GT.0.) THEN
                   NMULTRG(K) = 35.E4*PRACG(K)*FMULT*1000.
                   QMULTRG(K) = NMULTRG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                   QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                   PRACG(K) = PRACG(K)-QMULTRG(K)

               END IF
               END IF
               END IF
            END IF
            END IF
!         END IF

!........................................................................
! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL

!           IF (IHAIL.EQ.0) THEN
           IF (PSACWS(K).GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QC3D(K).GE.0.5E-3) THEN

! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
             PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
                          ASN(K)*ASN(K)/ &
                           (RHO(K)*LAMS(K)**(2.*BS+2.)))

! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
             DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.)

! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
             NSCNG(K) = DUM/MG0*RHO(K)
! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
             NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

! PORTION OF RIMING LEFT FOR SNOW
             PSACWS(K) = PSACWS(K) - PGSACW(K)
             END IF
           END IF

! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

           IF (PRACS(K).GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
              DUM = CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3 &
                   /(CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3+ &
                   CONS19*(4./LAMR(K))**3*(4./LAMR(K))**3)
              DUM=MIN(DUM,1.)
              DUM=MAX(DUM,0.)
              PGRACS(K) = (1.-DUM)*PRACS(K)
            NGRACS(K) = (1.-DUM)*NPRACS(K)
! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
            NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
            NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

! AMOUNT LEFT FOR SNOW PRODUCTION
            PRACS(K) = PRACS(K) - PGRACS(K)
            NPRACS(K) = NPRACS(K) - NGRACS(K)
! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
            PSACR(K)=PSACR(K)*(1.-DUM)
            END IF
           END IF
!           END IF

!.......................................................................
! FREEZING OF RAIN DROPS
! FREEZING ALLOWED BELOW -4 C

         IF (T3D(K).LT.269.15.AND.QR3D(K).GE.QSMALL) THEN

! IMMERSION FREEZING (BIGG 1953)
!            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3 &
!                 /LAMR(K)**3

!            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3

! hm fix 7/15/13 for consistency w/ original formula
            MNUCCR(K) = CONS20*NR3D(K)*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3 &
                 /LAMR(K)**3

            NNUCCR(K) = PI*NR3D(K)*BIMM*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3

! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
            NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

         END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

         IF (QR3D(K).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
!            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
            NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
         END IF

!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

         IF (QI3D(K).GE.1.E-8 .AND.QVQVSI(K).GE.1.) THEN

!           COFFI = 2./LAMI(K)
!           IF (COFFI.GE.DCS) THEN
              NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
              PRCI(K) = CONS22*NPRCI(K)
              NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

!           END IF
         END IF

!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

         IF (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) THEN
            PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                     LAMS(K)**(BS+3.)
            NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                  RHO(K)*N0S(K)/                                 &
                  LAMS(K)**(BS+3.)
            NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
         END IF

!.......................................................................
! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
! FOLLOWS REISNER ET AL. 1998
! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

         IF (QR3D(K).GE.1.E-8.AND.QI3D(K).GE.1.E-8.AND.T3D(K).LE.273.15) THEN

! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
! OTHERWISE ADD TO SNOW

            IF (QR3D(K).GE.0.1E-3) THEN
            NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
            NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
            ELSE
            NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
            NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
            END IF
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

         IF (INUC.EQ.0) THEN

! add threshold according to Greg Thomspon

         if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
              QVQVSI(K).ge.1.08) then

! hm, modify dec. 5, 2006, replace with cooper curve
      kc2 = 0.005*exp(0.304*(273.15-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
      kc2 = min(kc2,500.e3)
      kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1

          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD(K) = NNUCCD(K)*MI0
          END IF

          END IF

          ELSE IF (INUC.EQ.1) THEN

          IF (T3D(K).LT.273.15.AND.QVQVSI(K).GT.1.) THEN

             KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD(K) = NNUCCD(K)*MI0
          END IF
          END IF

         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 101      CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

! NO VENTILATION FOR CLOUD ICE

        IF (QI3D(K).GE.QSMALL) THEN

         EPSI = 2.*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

      ELSE
         EPSI = 0.
      END IF

      IF (QNI3D(K).GE.QSMALL) THEN
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
      ELSE
      EPSS = 0.
      END IF

      IF (QG3D(K).GE.QSMALL) THEN
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))


      ELSE
      EPSG = 0.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF

! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
              IF (QI3D(K).GE.QSMALL) THEN
              DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
              PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
              ELSE
              DUM=0.
              END IF
! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
              IF (QNI3D(K).GE.QSMALL) THEN
              PRDS(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
                EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
! OTHERWISE ADD TO CLOUD ICE
              ELSE
              PRD(K) = PRD(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
              END IF
! VAPOR DPEOSITION ON GRAUPEL
              PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)

! NO CONDENSATION ONTO RAIN, ONLY EVAP

           IF (QV3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
           ELSE
              PRE(K) = 0.
           END IF

! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
! FORMULA FROM REISNER 2 SCHEME

           DUM = (QV3D(K)-QVI(K))/DT

           FUDGEF = 0.9999
           SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)

           IF( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
               (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
               MNUCCD(K) = FUDGEF*MNUCCD(K)*DUM/SUM_DEP
               PRD(K) = FUDGEF*PRD(K)*DUM/SUM_DEP
               PRDS(K) = FUDGEF*PRDS(K)*DUM/SUM_DEP
               PRDG(K) = FUDGEF*PRDG(K)*DUM/SUM_DEP
           ENDIF

! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

           IF (PRD(K).LT.0.) THEN
              EPRD(K)=PRD(K)
              PRD(K)=0.
           END IF
           IF (PRDS(K).LT.0.) THEN
              EPRDS(K)=PRDS(K)
              PRDS(K)=0.
           END IF
           IF (PRDG(K).LT.0.) THEN
              EPRDG(K)=PRDG(K)
              PRDG(K)=0.
           END IF
!.......................................................................
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CONSERVATION OF WATER
! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.

! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

! NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
! STEP

!****SENSITIVITY - NO ICE

      IF (ILIQ.EQ.1) THEN
      MNUCCC(K)=0.
      NNUCCC(K)=0.
      MNUCCR(K)=0.
      NNUCCR(K)=0.
      MNUCCD(K)=0.
      NNUCCD(K)=0.
      END IF

! ****SENSITIVITY - NO GRAUPEL
      IF (IGRAUP.EQ.1) THEN
            PRACG(K) = 0.
            PSACR(K) = 0.
            PSACWG(K) = 0.
            PRDG(K) = 0.
            EPRDG(K) = 0.
            EVPMG(K) = 0.
            PGMLT(K) = 0.
            NPRACG(K) = 0.
            NPSACWG(K) = 0.
            NSCNG(K) = 0.
            NGRACS(K) = 0.
            NSUBG(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.
! fix 053011
            PIACRS(K)=PIACRS(K)+PIACR(K)
            PIACR(K) = 0.
! fix 070713
            PRACIS(K)=PRACIS(K)+PRACI(K)
            PRACI(K) = 0.
            PSACWS(K)=PSACWS(K)+PGSACW(K)
            PGSACW(K) = 0.
            PRACS(K)=PRACS(K)+PGRACS(K)
            PGRACS(K) = 0.
       END IF

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT

      IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN
        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO
        MNUCCC(K) = MNUCCC(K)*RATIO
        PSACWS(K) = PSACWS(K)*RATIO
        PSACWI(K) = PSACWI(K)*RATIO
        QMULTS(K) = QMULTS(K)*RATIO
        QMULTG(K) = QMULTG(K)*RATIO
        PSACWG(K) = PSACWG(K)*RATIO
        PGSACW(K) = PGSACW(K)*RATIO
        END IF

! CONSERVATION OF QI

      DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
                -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT

      IF (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) THEN

        RATIO = (QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
                     MNUCCD(K)+PSACWI(K))/ &
                      (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))

        PRCI(K) = PRCI(K)*RATIO
        PRAI(K) = PRAI(K)*RATIO
        PRACI(K) = PRACI(K)*RATIO
        PRACIS(K) = PRACIS(K)*RATIO
        EPRD(K) = EPRD(K)*RATIO

        END IF

! CONSERVATION OF QR

      DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
             PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT

      IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

        RATIO = (QR3D(K)/DT+PRC(K)+PRA(K))/ &
             (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))

        PRE(K) = PRE(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO
        QMULTR(K) = QMULTR(K)*RATIO
        QMULTRG(K) = QMULTRG(K)*RATIO
        MNUCCR(K) = MNUCCR(K)*RATIO
        PIACR(K) = PIACR(K)*RATIO
        PIACRS(K) = PIACRS(K)*RATIO
        PGRACS(K) = PGRACS(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        END IF

! CONSERVATION OF QNI
! CONSERVATION FOR GRAUPEL SCHEME

        IF (IGRAUP.EQ.0) THEN

      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT

      IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

        RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K))

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       END IF

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       ELSE IF (IGRAUP.EQ.1) THEN

      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT

      IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

       RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K))

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       END IF

       END IF

! CONSERVATION OF QG

      DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT

      IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

        RATIO = (QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
                  PIACR(K)+PRACI(K))/(-EPRDG(K))

       EPRDG(K) = EPRDG(K)*RATIO

      END IF

! TENDENCIES

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))

! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
      T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
               *XXLV(K)+(PRD(K)+PRDS(K)+                            &
                MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
               (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
                QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
                +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*XLF(K))/CPM(K)

      QC3DTEN(K) = QC3DTEN(K)+                                      &
                 (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
                  PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
      QI3DTEN(K) = QI3DTEN(K)+                                      &
         (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
                  PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
      QR3DTEN(K) = QR3DTEN(K)+                                      &
                 (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
             -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))
      IF (IGRAUP.EQ.0) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
                    PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
      NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
      ELSE IF (IGRAUP.EQ.1) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))

      END IF

      NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
            -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))

      NI3DTEN(K) = NI3DTEN(K)+                                      &
       (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
               NNUCCD(K)-NIACR(K)-NIACRS(K))

      NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
                   +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))

! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

        C2PREC(K) = PRA(K)+PRC(K)+PSACWS(K)+QMULTS(K)+QMULTG(K)+PSACWG(K)+ &
       PGSACW(K)+MNUCCC(K)+PSACWI(K)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)+DT*T3DTEN(K)
      DUMQV = QV3D(K) + DT * QV3DTEN(K)

      ! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = EP_2*dum/(PRES(K)-dum)

      DUMQC = QC3D(K) + DT * QC3DTEN(K)

      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS

      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT

      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
           PCC(K) = -DUMQC/DT
      END IF

      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
! Fortran version
!.......................................................................
! ACTIVATION OF CLOUD DROPLETS
! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
! DROPLET CONCENTRATION IS SPECIFIED !!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

      IF (EPRD(K).LT.0.) THEN
         DUM = EPRD(K)*DT/QI3D(K)
            DUM = MAX(-1.,DUM)
         NSUBI(K) = DUM*NI3D(K)/DT
      END IF
      IF (EPRDS(K).LT.0.) THEN
         DUM = EPRDS(K)*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSUBS(K) = DUM*NS3D(K)/DT
      END IF
      IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      END IF
      IF (EPRDG(K).LT.0.) THEN
         DUM = EPRDG(K)*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NSUBG(K) = DUM*NG3D(K)/DT
      END IF

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
         NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
         NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
         NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
         NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)

         END IF !!!!!! TEMPERATURE

! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
         LTRUE = 1

 200     CONTINUE

        END DO

! INITIALIZE PRECIP AND SNOW RATES
      PRECRT = 0.
      SNOWRT = 0.
! hm added 7/13/13
      SNOWPRT = 0.
      GRPLPRT = 0.

! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE
      DO K = KTS,KTE
! Fortran version
     END DO
        IF (LTRUE.EQ.0) GOTO 400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!.......................................................................
! CALCULATE SEDIMENATION
! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
! STABILITY, I.E. COURANT# < 1

!.......................................................................

      NSTEP = 1

      DO K = KTE,KTS,-1

        DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
        DUMQS(K) = QNI3D(K)+QNI3DTEN(K)*DT
        DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
        DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
        DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
        DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
        DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
        DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
        DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
        DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT

! SWITCH FOR CONSTANT DROPLET NUMBER
        IF (iinum.EQ.1) THEN
        DUMFNC(K) = NC3D(K)
        END IF

! GET DUMMY LAMDA FOR SEDIMENTATION CALCULATIONS

! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
      DUMFNI(K) = MAX(0.,DUMFNI(K))
      DUMFNS(K) = MAX(0.,DUMFNS(K))
      DUMFNC(K) = MAX(0.,DUMFNC(K))
      DUMFNR(K) = MAX(0.,DUMFNR(K))
      DUMFNG(K) = MAX(0.,DUMFNG(K))

!......................................................................
! CLOUD ICE

      IF (DUMI(K).GE.QSMALL) THEN
        DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
        DLAMI=MAX(DLAMI,LAMMINI)
        DLAMI=MIN(DLAMI,LAMMAXI)
      END IF
!......................................................................
! RAIN

      IF (DUMR(K).GE.QSMALL) THEN
        DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1./3.)
        DLAMR=MAX(DLAMR,LAMMINR)
        DLAMR=MIN(DLAMR,LAMMAXR)
      END IF
!......................................................................
! CLOUD DROPLETS

      IF (DUMC(K).GE.QSMALL) THEN
         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

        DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.)/(DUMC(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
        LAMMIN = (PGAM(K)+1.)/60.E-6
        LAMMAX = (PGAM(K)+1.)/1.E-6
        DLAMC=MAX(DLAMC,LAMMIN)
        DLAMC=MIN(DLAMC,LAMMAX)
      END IF
!......................................................................
! SNOW

      IF (DUMQS(K).GE.QSMALL) THEN
        DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1./DS)
        DLAMS=MAX(DLAMS,LAMMINS)
        DLAMS=MIN(DLAMS,LAMMAXS)
      END IF
!......................................................................
! GRAUPEL

      IF (DUMG(K).GE.QSMALL) THEN
        DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1./DG)
        DLAMG=MAX(DLAMG,LAMMING)
        DLAMG=MIN(DLAMG,LAMMAXG)
      END IF
!......................................................................
! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

! CLOUD WATER

      IF (DUMC(K).GE.QSMALL) THEN
      UNC =  ACN(K)*GAMMA(1.+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.))
      UMC = ACN(K)*GAMMA(4.+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.))
      ELSE
      UMC = 0.
      UNC = 0.
      END IF

      IF (DUMI(K).GE.QSMALL) THEN
      UNI =  AIN(K)*CONS27/DLAMI**BI
      UMI = AIN(K)*CONS28/(DLAMI**BI)
      ELSE
      UMI = 0.
      UNI = 0.
      END IF

      IF (DUMR(K).GE.QSMALL) THEN
      UNR = ARN(K)*CONS6/DLAMR**BR
      UMR = ARN(K)*CONS4/(DLAMR**BR)
      ELSE
      UMR = 0.
      UNR = 0.
      END IF

      IF (DUMQS(K).GE.QSMALL) THEN
      UMS = ASN(K)*CONS3/(DLAMS**BS)
      UNS = ASN(K)*CONS5/DLAMS**BS
      ELSE
      UMS = 0.
      UNS = 0.
      END IF

      IF (DUMG(K).GE.QSMALL) THEN
      UMG = AGN(K)*CONS7/(DLAMG**BG)
      UNG = AGN(K)*CONS8/DLAMG**BG
      ELSE
      UMG = 0.
      UNG = 0.
      END IF

! SET REALISTIC LIMITS ON FALLSPEED

! bug fix, 10/08/09
        dum=(rhosu/rho(k))**0.54
        UMS=MIN(UMS,1.2*dum)
        UNS=MIN(UNS,1.2*dum)
! fix 053011
! fix for correction by AA 4/6/11
        UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)
        UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)
        UMR=MIN(UMR,9.1*dum)
        UNR=MIN(UNR,9.1*dum)
        UMG=MIN(UMG,20.*dum)
        UNG=MIN(UNG,20.*dum)

      FR(K) = UMR
      FI(K) = UMI
      FNI(K) = UNI
      FS(K) = UMS
      FNS(K) = UNS
      FNR(K) = UNR
      FC(K) = UMC
      FNC(K) = UNC
      FG(K) = UMG
      FNG(K) = UNG

! V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

        IF (K.LE.KTE-1) THEN
        IF (FR(K).LT.1.E-10) THEN
        FR(K)=FR(K+1)
        END IF
        IF (FI(K).LT.1.E-10) THEN
        FI(K)=FI(K+1)
        END IF
        IF (FNI(K).LT.1.E-10) THEN
        FNI(K)=FNI(K+1)
        END IF
        IF (FS(K).LT.1.E-10) THEN
        FS(K)=FS(K+1)
        END IF
        IF (FNS(K).LT.1.E-10) THEN
        FNS(K)=FNS(K+1)
        END IF
        IF (FNR(K).LT.1.E-10) THEN
        FNR(K)=FNR(K+1)
        END IF
        IF (FC(K).LT.1.E-10) THEN
        FC(K)=FC(K+1)
        END IF
        IF (FNC(K).LT.1.E-10) THEN
        FNC(K)=FNC(K+1)
        END IF
        IF (FG(K).LT.1.E-10) THEN
        FG(K)=FG(K+1)
        END IF
        IF (FNG(K).LT.1.E-10) THEN
        FNG(K)=FNG(K+1)
        END IF
        END IF ! K LE KTE-1

! CALCULATE NUMBER OF SPLIT TIME STEPS

      RGVM = MAX(FR(K),FI(K),FS(K),FC(K),FNI(K),FNR(K),FNS(K),FNC(K),FG(K),FNG(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEP = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEP)
! Fortran version
! MULTIPLY VARIABLES BY RHO
      DUMR(k) = DUMR(k)*RHO(K)
      DUMI(k) = DUMI(k)*RHO(K)
      DUMFNI(k) = DUMFNI(K)*RHO(K)
      DUMQS(k) = DUMQS(K)*RHO(K)
      DUMFNS(k) = DUMFNS(K)*RHO(K)
      DUMFNR(k) = DUMFNR(K)*RHO(K)
      DUMC(k) = DUMC(K)*RHO(K)
      DUMFNC(k) = DUMFNC(K)*RHO(K)
      DUMG(k) = DUMG(K)*RHO(K)
      DUMFNG(k) = DUMFNG(K)*RHO(K)

      END DO

      DO N = 1,NSTEP

      DO K = KTS,KTE
      FALOUTR(K) = FR(K)*DUMR(K)
      FALOUTI(K) = FI(K)*DUMI(K)
      FALOUTNI(K) = FNI(K)*DUMFNI(K)
      FALOUTS(K) = FS(K)*DUMQS(K)
      FALOUTNS(K) = FNS(K)*DUMFNS(K)
      FALOUTNR(K) = FNR(K)*DUMFNR(K)
      FALOUTC(K) = FC(K)*DUMC(K)
      FALOUTNC(K) = FNC(K)*DUMFNC(K)
      FALOUTG(K) = FG(K)*DUMG(K)
      FALOUTNG(K) = FNG(K)*DUMFNG(K)
      END DO

! TOP OF MODEL

      K = KTE
      FALTNDR = FALOUTR(K)/DZQ(k)
      FALTNDI = FALOUTI(K)/DZQ(k)
      FALTNDNI = FALOUTNI(K)/DZQ(k)
      FALTNDS = FALOUTS(K)/DZQ(k)
      FALTNDNS = FALOUTNS(K)/DZQ(k)
      FALTNDNR = FALOUTNR(K)/DZQ(k)
      FALTNDC = FALOUTC(K)/DZQ(k)
      FALTNDNC = FALOUTNC(K)/DZQ(k)
      FALTNDG = FALOUTG(K)/DZQ(k)
      FALTNDNG = FALOUTNG(K)/DZQ(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

      QRSTEN(K) = QRSTEN(K)-FALTNDR/NSTEP/RHO(k)
      QISTEN(K) = QISTEN(K)-FALTNDI/NSTEP/RHO(k)
      NI3DTEN(K) = NI3DTEN(K)-FALTNDNI/NSTEP/RHO(k)
      QNISTEN(K) = QNISTEN(K)-FALTNDS/NSTEP/RHO(k)
      NS3DTEN(K) = NS3DTEN(K)-FALTNDNS/NSTEP/RHO(k)
      NR3DTEN(K) = NR3DTEN(K)-FALTNDNR/NSTEP/RHO(k)
      QCSTEN(K) = QCSTEN(K)-FALTNDC/NSTEP/RHO(k)
      NC3DTEN(K) = NC3DTEN(K)-FALTNDNC/NSTEP/RHO(k)
      QGSTEN(K) = QGSTEN(K)-FALTNDG/NSTEP/RHO(k)
      NG3DTEN(K) = NG3DTEN(K)-FALTNDNG/NSTEP/RHO(k)

      DUMR(K) = DUMR(K)-FALTNDR*DT/NSTEP
      DUMI(K) = DUMI(K)-FALTNDI*DT/NSTEP
      DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT/NSTEP
      DUMQS(K) = DUMQS(K)-FALTNDS*DT/NSTEP
      DUMFNS(K) = DUMFNS(K)-FALTNDNS*DT/NSTEP
      DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT/NSTEP
      DUMC(K) = DUMC(K)-FALTNDC*DT/NSTEP
      DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT/NSTEP
      DUMG(K) = DUMG(K)-FALTNDG*DT/NSTEP
      DUMFNG(K) = DUMFNG(K)-FALTNDNG*DT/NSTEP

      DO K = KTE-1,KTS,-1
      FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)
      FALTNDI = (FALOUTI(K+1)-FALOUTI(K))/DZQ(K)
      FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))/DZQ(K)
      FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
      FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
      FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)
      FALTNDC = (FALOUTC(K+1)-FALOUTC(K))/DZQ(K)
      FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))/DZQ(K)
      FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
      FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)

! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

      QRSTEN(K) = QRSTEN(K)+FALTNDR/NSTEP/RHO(k)
      QISTEN(K) = QISTEN(K)+FALTNDI/NSTEP/RHO(k)
      NI3DTEN(K) = NI3DTEN(K)+FALTNDNI/NSTEP/RHO(k)
      QNISTEN(K) = QNISTEN(K)+FALTNDS/NSTEP/RHO(k)
      NS3DTEN(K) = NS3DTEN(K)+FALTNDNS/NSTEP/RHO(k)
      NR3DTEN(K) = NR3DTEN(K)+FALTNDNR/NSTEP/RHO(k)
      QCSTEN(K) = QCSTEN(K)+FALTNDC/NSTEP/RHO(k)
      NC3DTEN(K) = NC3DTEN(K)+FALTNDNC/NSTEP/RHO(k)
      QGSTEN(K) = QGSTEN(K)+FALTNDG/NSTEP/RHO(k)
      NG3DTEN(K) = NG3DTEN(K)+FALTNDNG/NSTEP/RHO(k)
! Fortran version
      DUMR(K) = DUMR(K)+FALTNDR*DT/NSTEP
      DUMI(K) = DUMI(K)+FALTNDI*DT/NSTEP
      DUMFNI(K) = DUMFNI(K)+FALTNDNI*DT/NSTEP
      DUMQS(K) = DUMQS(K)+FALTNDS*DT/NSTEP
      DUMFNS(K) = DUMFNS(K)+FALTNDNS*DT/NSTEP
      DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT/NSTEP
      DUMC(K) = DUMC(K)+FALTNDC*DT/NSTEP
      DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT/NSTEP
      DUMG(K) = DUMG(K)+FALTNDG*DT/NSTEP
      DUMFNG(K) = DUMFNG(K)+FALTNDNG*DT/NSTEP

! FOR WRF-CHEM, NEED PRECIP RATES (UNITS OF KG/M^2/S)
          CSED(K)=CSED(K)+FALOUTC(K)/NSTEP
          ISED(K)=ISED(K)+FALOUTI(K)/NSTEP
          SSED(K)=SSED(K)+FALOUTS(K)/NSTEP
          GSED(K)=GSED(K)+FALOUTG(K)/NSTEP
          RSED(K)=RSED(K)+FALOUTR(K)/NSTEP
      END DO

! GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
! FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
! OF LIQUID WATER CANCELS THIS FACTOR OF 1000

        PRECRT = PRECRT+(FALOUTR(KTS)+FALOUTC(KTS)+FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))  &
                     *DT/NSTEP
        SNOWRT = SNOWRT+(FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))*DT/NSTEP
! hm added 7/13/13
        SNOWPRT = SNOWPRT+(FALOUTI(KTS)+FALOUTS(KTS))*DT/NSTEP
        GRPLPRT = GRPLPRT+(FALOUTG(KTS))*DT/NSTEP
      END DO

        DO K=KTS,KTE
! Fortran version
! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

        QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
        QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
        QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
        QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
        QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)
! Fortran version
! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

!hm 4/7/09 bug fix
!        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15.AND.LAMI(K).GE.1.E-10) THEN
        IF (1./LAMI(K).GE.2.*DCS) THEN
           QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
           NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
           QI3DTEN(K) = -QI3D(K)/DT
           NI3DTEN(K) = -NI3D(K)/DT
        END IF
        END IF

! hm add tendencies here, then call sizeparameter
! to ensure consisitency between mixing ratio and number concentration

          QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
          QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
          QNI3D(k)        = QNI3D(k)+QNI3DTEN(k)*DT
          QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
          NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
          NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
          NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
          NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT
! Fortran version
          IF (IGRAUP.EQ.0) THEN
          QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
          NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
          END IF

! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
          T3D(K)         = T3D(K)+T3DTEN(k)*DT
          QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT
! Fortran version
! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)
! Fortran version
! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! hm 7/9/09 change limit to 1.e-8

             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
               END IF
             END IF
             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
             END IF
! Fortran version
!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
         QC3D(K) = 0.
         NC3D(K) = 0.
         EFFC(K) = 0.
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
         QR3D(K) = 0.
         NR3D(K) = 0.
         EFFR(K) = 0.
       END IF
       IF (QI3D(K).LT.QSMALL) THEN
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       END IF
       IF (QNI3D(K).LT.QSMALL) THEN
         QNI3D(K) = 0.
         NS3D(K) = 0.
         EFFS(K) = 0.
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
         QG3D(K) = 0.
         NG3D(K) = 0.
         EFFG(K) = 0.
       END IF
! Fortran version
!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

            IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) GOTO 500

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE INSTANTANEOUS PROCESSES

! ADD MELTING OF CLOUD ICE TO FORM RAIN

        IF (QI3D(K).GE.QSMALL.AND.T3D(K).GE.273.15) THEN
           QR3D(K) = QR3D(K)+QI3D(K)
           T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
           QI3D(K) = 0.
           NR3D(K) = NR3D(K)+NI3D(K)
           NI3D(K) = 0.
        END IF
! Fortran version
! ****SENSITIVITY - NO ICE
        IF (ILIQ.EQ.1) GOTO 778

! HOMOGENEOUS FREEZING OF CLOUD WATER

        IF (T3D(K).LE.233.15.AND.QC3D(K).GE.QSMALL) THEN
           QI3D(K)=QI3D(K)+QC3D(K)
           T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
           QC3D(K)=0.
           NI3D(K)=NI3D(K)+NC3D(K)
           NC3D(K)=0.
        END IF
! Fortran version
! HOMOGENEOUS FREEZING OF RAIN

        IF (IGRAUP.EQ.0) THEN

        IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
           QG3D(K) = QG3D(K)+QR3D(K)
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
           QR3D(K) = 0.
           NG3D(K) = NG3D(K)+ NR3D(K)
           NR3D(K) = 0.
        END IF

        ELSE IF (IGRAUP.EQ.1) THEN

        IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
           QNI3D(K) = QNI3D(K)+QR3D(K)
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
           QR3D(K) = 0.
           NS3D(K) = NS3D(K)+NR3D(K)
           NR3D(K) = 0.
        END IF

        END IF
! Fortran version
 778    CONTINUE

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF

      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         DUM = PRES(K)/(287.15*T3D(K))
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1
      NS3D(K) = N0S(K)/LAMS(K)
      END IF

      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF

      END IF

 500  CONTINUE

! CALCULATE EFFECTIVE RADIUS

      IF (QI3D(K).GE.QSMALL) THEN
         EFFI(K) = 3./LAMI(K)/2.*1.E6
      ELSE
         EFFI(K) = 25.
      END IF

      IF (QNI3D(K).GE.QSMALL) THEN
         EFFS(K) = 3./LAMS(K)/2.*1.E6
      ELSE
         EFFS(K) = 25.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
         EFFR(K) = 3./LAMR(K)/2.*1.E6
      ELSE
         EFFR(K) = 25.
      END IF

      IF (QC3D(K).GE.QSMALL) THEN
      EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
             GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
      ELSE
      EFFC(K) = 25.
      END IF

      IF (QG3D(K).GE.QSMALL) THEN
         EFFG(K) = 3./LAMG(K)/2.*1.E6
      ELSE
         EFFG(K) = 25.
      END IF

! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
!          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
! HM, 12/28/12, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
! OF EXCESSIVE AND PERSISTENT ANVIL
! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
          NI3D(K) = MIN(NI3D(K),0.3E6/RHO(K))

! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
          IF (iinum.EQ.0.AND.IACT.EQ.2) THEN
          NC3D(K) = MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K))
          END IF
! SWITCH FOR CONSTANT DROPLET NUMBER
          IF (iinum.EQ.1) THEN
! CHANGE NDCNST FROM CM-3 TO KG-1
             NC3D(K) = NDCNST*1.E6/RHO(K)
          END IF

      END DO !!! K LOOP

 400         CONTINUE

! ALL DONE !!!!!!!!!!!
      RETURN
      END SUBROUTINE MORR_TWO_MOMENT_MICRO

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL(C_DOUBLE) FUNCTION POLYSVP (T,TYPE)

!-------------------------------------------

!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

      IMPLICIT NONE

      REAL(C_DOUBLE) DUM
      REAL(C_DOUBLE) T
      INTEGER TYPE
! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE.EQ.1) THEN

!         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+                                              &
!          LOG10(6.1071))*100.
! hm 11/16/20, use Goff-Gratch for T < 195.8 K and Flatau et al. equal or above 195.8 K
      if (t.ge.195.8) then
         dt=t-273.15
         polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
         polysvp = polysvp*100.
      else

         polysvp = 10.**(-9.09718*(273.16/t-1.)-3.56654* &
          dlog10(273.16d0/t)+0.876793*(1.-t/273.16)+ &
          dlog10(6.1071d0))*100.

      end if

      END IF

! LIQUID

      IF (TYPE.EQ.0) THEN

!         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
!             5.02808*LOG10(373.16/T)-                                                                  &
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+                                &
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+                              &
!             LOG10(1013.246))*100.
! hm 11/16/20, use Goff-Gratch for T < 202.0 K and Flatau et al. equal or above 202.0 K
      if (t.ge.202.0) then
         dt = t-273.15
         polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
         polysvp = polysvp*100.
      else

! note: uncertain below -70 C, but produces physical values (non-negative) unlike flatau
         polysvp = 10.**(-7.90298*(373.16/t-1.)+ &
             5.02808*dlog10(373.16d0/t)- &
             1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+ &
             8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+ &
             dlog10(1013.246d0))*100.
      end if

      END IF


      END FUNCTION POLYSVP

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
      implicit none
      INTEGER I,N
      LOGICAL PARITY
      REAL(C_DOUBLE)                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL(C_DOUBLE), DIMENSION(7) :: C
      REAL(C_DOUBLE), DIMENSION(8) :: P
      REAL(C_DOUBLE), DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,                                 &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,    &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+xxx
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA
!------------------------------------------------------------------------------
END MODULE module_mp_morr_two_moment
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
