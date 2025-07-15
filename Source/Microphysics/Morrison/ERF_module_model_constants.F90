!WRF:MODEL_LAYER:CONSTANTS
!

 MODULE module_model_constants

   USE ISO_C_BINDING

   !  2. Following are constants for use in defining real number bounds.

   !  A really small number.

   REAL(C_DOUBLE)    , PARAMETER :: epsilon         = 1.E-15

   !  4. Following is information related to the physical constants.

   !  These are the physical constants used within the model.

! JM NOTE -- can we name this grav instead?
   REAL(C_DOUBLE)    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)

   REAL(C_DOUBLE)    , PARAMETER :: r_d          = 287.       ! gas constant of dry air (J deg^-1 kg^-1)
   REAL(C_DOUBLE)    , PARAMETER :: cp           = 7.*r_d/2.  !

   REAL(C_DOUBLE)    , PARAMETER :: r_v          = 461.6      ! gas constant for water vapor (J deg^-1 kg^-1)
   REAL(C_DOUBLE)    , PARAMETER :: cv           = cp-r_d     ! Specific heat of air at contant volume (J deg^-1 kg^-1)
   REAL(C_DOUBLE)    , PARAMETER :: cpv          = 4.*r_v
   REAL(C_DOUBLE)    , PARAMETER :: cvv          = cpv-r_v    !
   REAL(C_DOUBLE)    , PARAMETER :: cvpm         = -cv/cp
   REAL(C_DOUBLE)    , PARAMETER :: cliq         = 4190.      ! specific heat of liquid water at 0^oC
   REAL(C_DOUBLE)    , PARAMETER :: cice         = 2106.      ! specific heat of ice at 0^oC
   REAL(C_DOUBLE)    , PARAMETER :: psat         = 610.78
   REAL(C_DOUBLE)    , PARAMETER :: rcv          = r_d/cv     !
   REAL(C_DOUBLE)    , PARAMETER :: rcp          = r_d/cp
   REAL(C_DOUBLE)    , PARAMETER :: rovg         = r_d/g
   REAL(C_DOUBLE)    , PARAMETER :: c2           = cp * rcv
   real    , parameter :: mwdry        = 28.966     ! molecular weight of dry air (g/mole)

   REAL(C_DOUBLE)    , PARAMETER :: p1000mb      = 100000.    ! pressure at 1000 hPa (pa)
   REAL(C_DOUBLE)    , PARAMETER :: t0           = 300.       ! base state tempertaure (K)
   REAL(C_DOUBLE)    , PARAMETER :: p0           = p1000mb    ! base state surface pressure (pa)
   REAL(C_DOUBLE)    , PARAMETER :: cpovcv       = cp/(cp-r_d)
   REAL(C_DOUBLE)    , PARAMETER :: cvovcp       = 1./cpovcv
   REAL(C_DOUBLE)    , PARAMETER :: rvovrd       = r_v/r_d

   REAL(C_DOUBLE)    , PARAMETER :: reradius     = 1./6370.0e03  ! reciprocal of earth radius (m^-1)

   REAL(C_DOUBLE)    , PARAMETER :: asselin      = .025
!   REAL(C_DOUBLE)    , PARAMETER :: asselin      = .0
   REAL(C_DOUBLE)    , PARAMETER :: cb           = 25.

   REAL(C_DOUBLE)    , PARAMETER :: XLV0         = 3.15E6       !  constant defined for calculation of latent heating
   REAL(C_DOUBLE)    , PARAMETER :: XLV1         = 2370.        !  constant defined for calculation of latent heating
   REAL(C_DOUBLE)    , PARAMETER :: XLS0         = 2.905E6      !  constant defined for calculation of latent heating
   REAL(C_DOUBLE)    , PARAMETER :: XLS1         = 259.532      !  constant defined for calculation of latent heating

   REAL(C_DOUBLE)    , PARAMETER :: XLS          = 2.85E6      ! latent heat of sublimation of water at 0^oC (J kg^-1)
   REAL(C_DOUBLE)    , PARAMETER :: XLV          = 2.5E6       ! latent heat of vaporization of water at 0^oC (J kg^-1)
   REAL(C_DOUBLE)    , PARAMETER :: XLF          = 3.50E5      ! latent heat of fusion of water at 0^oC (J kg^-1)

   REAL(C_DOUBLE)    , PARAMETER :: rhowater     = 1000.       ! density of liquid water at 0^oC (kg m^-3)
   REAL(C_DOUBLE)    , PARAMETER :: rhosnow      = 100.        ! density of snow (kg m^-3)
   REAL(C_DOUBLE)    , PARAMETER :: rhoair0      = 1.28        ! density of dry air at 0^oC and 1000mb pressure (kg m^-3)

   REAL(C_DOUBLE)    , PARAMETER :: RE_QC_BG     = 2.49E-6     ! effective radius of cloud for background (m)
   REAL(C_DOUBLE)    , PARAMETER :: RE_QI_BG     = 4.99E-6     ! effective radius of ice for background (m)
   REAL(C_DOUBLE)    , PARAMETER :: RE_QS_BG     = 9.99E-6     ! effective radius of snow for background (m)
   REAL(C_DOUBLE)    , PARAMETER :: RE_QC_MAX    =  50.E-6     ! max effective radius of cloud allowed
   REAL(C_DOUBLE)    , PARAMETER :: RE_QI_MAX    = 125.E-6     ! max effective radius of ice allowed
   REAL(C_DOUBLE)    , PARAMETER :: RE_QS_MAX    = 999.E-6     ! max effective radius of snow allowed
!
! Now namelist-specified parameter: ccn_conc - RAS
!   REAL(C_DOUBLE)    , PARAMETER :: n_ccn0       = 1.0E8
!
   REAL(C_DOUBLE)    , PARAMETER :: piconst      = 3.1415926535897932384626433    ! constant of PI
   REAL(C_DOUBLE)    , PARAMETER :: DEGRAD       = piconst/180.                   ! radians per degree
   REAL(C_DOUBLE)    , PARAMETER :: DPD          = 360./365.

   REAL(C_DOUBLE)    , PARAMETER ::  SVP1=0.6112      ! constant for saturation vapor pressure calculation (dimensionless)
   REAL(C_DOUBLE)    , PARAMETER ::  SVP2=17.67       ! constant for saturation vapor pressure calculation (dimensionless)
   REAL(C_DOUBLE)    , PARAMETER ::  SVP3=29.65       ! constant for saturation vapor pressure calculation (K)
   REAL(C_DOUBLE)    , PARAMETER ::  SVPT0=273.15     ! constant for saturation vapor pressure calculation (K)
   REAL(C_DOUBLE)    , PARAMETER ::  EP_1=R_v/R_d-1.  !  constant for virtual temperature (r_v/r_d - 1) (dimensionless)
   REAL(C_DOUBLE)    , PARAMETER ::  EP_2=R_d/R_v     ! constant for specific humidity calculation (dimensionless)
   REAL(C_DOUBLE)    , PARAMETER ::  KARMAN=0.4               ! von Karman constant
   REAL(C_DOUBLE)    , PARAMETER ::  EOMEG=7.2921E-5          ! angular velocity of rotation (rad^-1)
   REAL(C_DOUBLE)    , PARAMETER ::  STBOLT=5.67051E-8        ! Stefan-Boltzmann constant (W m^-2 deg^-4)

   REAL(C_DOUBLE)    , PARAMETER ::  prandtl = 1./3.0   ! prandtl's mixing length (m)
                                              ! constants for w-damping option
   REAL(C_DOUBLE)    , PARAMETER ::  w_alpha = 0.3      ! strength m/s/s
   REAL(C_DOUBLE)    , PARAMETER ::  w_beta  = 1.0      ! activation cfl number

       REAL(C_DOUBLE) , PARAMETER ::  pq0=379.90516     !
       REAL(C_DOUBLE) , PARAMETER ::  epsq2=0.2         ! initial TKE for camuw PBL scheme (m2 s^-2)
       REAL(C_DOUBLE) , PARAMETER ::  a2=17.2693882
       REAL(C_DOUBLE) , PARAMETER ::  a3=273.16
       REAL(C_DOUBLE) , PARAMETER ::  a4=35.86
       REAL(C_DOUBLE) , PARAMETER ::  epsq=1.e-12      ! threshold specified for SPECIFIC HUMIDITY calculation in BMJ cumulus scheme (kg kg^-1)
       REAL(C_DOUBLE) , PARAMETER ::  p608=rvovrd-1.
       REAL(C_DOUBLE) , PARAMETER ::  climit=1.e-20
       REAL(C_DOUBLE) , PARAMETER ::  cm1=2937.4
       REAL(C_DOUBLE) , PARAMETER ::  cm2=4.9283
       REAL(C_DOUBLE) , PARAMETER ::  cm3=23.5518
!       REAL(C_DOUBLE) , PARAMETER ::  defc=8.0
!       REAL(C_DOUBLE) , PARAMETER ::  defm=32.0
       REAL(C_DOUBLE) , PARAMETER ::  defc=0.0
       REAL(C_DOUBLE) , PARAMETER ::  defm=99999.0
       REAL(C_DOUBLE) , PARAMETER ::  epsfc=1./1.05
       REAL(C_DOUBLE) , PARAMETER ::  epswet=0.0
       REAL(C_DOUBLE) , PARAMETER ::  fcdif=1./3.
       REAL(C_DOUBLE) , PARAMETER ::  fcm=0.00003
       REAL(C_DOUBLE) , PARAMETER ::  gma=-r_d*(1.-rcp)*0.5
       REAL(C_DOUBLE) , PARAMETER ::  p400=40000.0
       REAL(C_DOUBLE) , PARAMETER ::  phitp=15000.0
       REAL(C_DOUBLE) , PARAMETER ::  pi2=2.*3.1415926, pi1=3.1415926
       REAL(C_DOUBLE) , PARAMETER ::  plbtm=105000.0
       REAL(C_DOUBLE) , PARAMETER ::  plomd=64200.0
       REAL(C_DOUBLE) , PARAMETER ::  pmdhi=35000.0
       REAL(C_DOUBLE) , PARAMETER ::  q2ini=0.50
       REAL(C_DOUBLE) , PARAMETER ::  rfcp=0.25/cp
       REAL(C_DOUBLE) , PARAMETER ::  rhcrit_land=0.75
       REAL(C_DOUBLE) , PARAMETER ::  rhcrit_sea=0.80
       REAL(C_DOUBLE) , PARAMETER ::  rlag=14.8125
       REAL(C_DOUBLE) , PARAMETER ::  rlx=0.90
       REAL(C_DOUBLE) , PARAMETER ::  scq2=50.0
       REAL(C_DOUBLE) , PARAMETER ::  slopht=0.001
       REAL(C_DOUBLE) , PARAMETER ::  tlc=2.*0.703972477
       REAL(C_DOUBLE) , PARAMETER ::  wa=0.15
       REAL(C_DOUBLE) , PARAMETER ::  wght=0.35
       REAL(C_DOUBLE) , PARAMETER ::  wpc=0.075
       REAL(C_DOUBLE) , PARAMETER ::  z0land=0.10    ! surface roughness length over land (m)
       REAL(C_DOUBLE) , PARAMETER ::  z0max=0.008    !  maximum roughness length (m)
       REAL(C_DOUBLE) , PARAMETER ::  z0sea=0.001   ! roughness length over ocean (m)


   !  Earth

   !  The value for P2SI *must* be set to 1.0 for Earth
   !  Although, now we may not need this declaration here (see above)
   !REAL(C_DOUBLE)    , PARAMETER :: P2SI         = 1.0

   !  Orbital constants:

   INTEGER , PARAMETER :: PLANET_YEAR = 365   ! number of days in a calendar year
   REAL(C_DOUBLE) , PARAMETER :: OBLIQUITY = 23.5       ! solar obliquity (degree)
   REAL(C_DOUBLE) , PARAMETER :: ECCENTRICITY = 0.014   ! Orbital eccentricity
   REAL(C_DOUBLE) , PARAMETER :: SEMIMAJORAXIS = 1.0    ! Ratio of semi-major axis of planet / semi-major axis of earth
   REAL(C_DOUBLE) , PARAMETER :: zero_date = 0.0        ! Time of perihelion passage
   REAL(C_DOUBLE) , PARAMETER :: EQUINOX_FRACTION= 0.0  ! Fraction into the year (from perhelion) of the occurrence of the Northern Spring Equinox

! 2012103
#if (EM_CORE == 1)
! for calls to set_tiles
   INTEGER, PARAMETER :: ZONE_SOLVE_EM = 1
   INTEGER, PARAMETER :: ZONE_SFS = 2
#endif

 CONTAINS
   SUBROUTINE init_module_model_constants
   END SUBROUTINE init_module_model_constants
 END MODULE module_model_constants
