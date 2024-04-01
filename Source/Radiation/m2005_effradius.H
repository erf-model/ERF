/*-------------------------------------------------------------------------------------------
 This subroutine is used to calculate droplet and ice crystal effective radius, which will be
 used in the CAM radiation code. The method to calculate effective radius is taken out of the
 Morrision two moment scheme from M2005MICRO_GRAUPEL. It is also very similar to the subroutine
 effradius in the module of cldwat2m in the CAM source codes.
 Adopted by Minghuai Wang (Minghuai.Wang@pnl.gov).
 Calculate effective radius for radiation code
 If no cloud water, default value is:
   10 micron for droplets,
   25 micron for cloud ice.
 Be careful of the unit of effective radius : [micro meter]

 NOTE: this code is modified from E3SM
-------------------------------------------------------------------------------------------*/
#include "ERF_Constants.H"
#include "Rad_constants.H"
#include <cmath>

using namespace amrex;
using yakl::intrinsics::size;
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

YAKL_INLINE
void m2005_effradius (const real2d& ql, const real2d& nl, const real2d& qi,
                      const real2d& ni, const real2d& qs, const real2d& ns,
                      const real2d& cld, const real2d& pres, const real2d& tk,
                      const real2d& effl, const real2d& effi,
                      const real2d& deffi, const real2d& lamcrad,
                      const real2d& pgamrad, const real2d& des)
{
    // Main computation
    const real pi = 3.1415926535897932384626434;
    const real qsmall = 1.0e-14; // in the SAM source code (module_mp_graupel)
    const real rhow = 997.;      // in module_mp_graupel, SAM
    const real rhoi = 500.;      // in both CAM and SAM
    const real dcs = 125.e-6;   // in module_mp_graupel, SAM
    const real ci = rhoi*pi/6.;
    const real di = 3.;

    // for snow water
    const real rhos = 100.;
    const real cs = rhos*pi/6.;
    const real ds = 3.;
    const real mincld  = 0.0001;

    int ncol = pres.extent(0);
    int nlev = pres.extent(1);
    // Effective diameters of snow crystals
    parallel_for (SimpleBounds<2> (ncol, nlev), YAKL_LAMBDA (int i, int j)
    {
        auto rho  = pres(i,j)/(287.15*tk(i,j));   // air density [kg/m3]
        auto cldm = std::max(cld(i,j), mincld);
        auto qlic = std::min(5.e-3, std::max(0., ql(i,j)/cldm));
        auto qiic = std::min(5.e-3, std::max(0., qi(i,j)/cldm));
        auto nlic = std::max(nl(i,j), 0.)/cldm;
        auto niic = std::max(ni(i,j), 0.)/cldm;

        real res;
        if(qs(i,j) > 1.0e-7) {
            auto lammaxs=1./10.e-6;
            auto lammins=1./2000.e-6;
            auto lams = std::pow((std::tgamma(1.+ds)*cs*ns(i,j)/qs(i,j)), (1./ds));
            lams = std::min(lammaxs, std::max(lams,lammins));
            res = 1.5/lams*1.0e6;
        }
        else {
            res = 500.;
        }

        // from Hugh Morrision: rhos/917 accouts for assumptions about
        // ice density in the Mitchell optics.
        des(i,j) = res*rhos/917.*2.;

        // Effective radius of cloud ice droplet
        if( qiic >= qsmall ) {
            niic   = std::min(niic, qiic*1.e20);
            auto lammax = 1./1.e-6;         // in module_mp_graupel, SAM
            auto lammin = 1./(2.*dcs+100.e-6);    // in module_mp_graupel, SAM
            auto lami   = std::pow((std::tgamma(1.+di)*ci*niic/qiic), (1./di));
            lami        = std::min(lammax, std::max(lami, lammin));
            effi(i,j)   = 1.5/lami*1.e6;
        }
        else {
            effi(i,j)   = 25.;
        }

        // hm ice effective radius for david mitchell's optics
        // ac morrison indicates that this is effective diameter
        // ac morrison indicates 917 (for the density of pure ice..)
        deffi(i,j)  = effi(i,j)*rhoi/917.*2.;

        // Effective radius of cloud liquid droplet
        if( qlic >= qsmall ) {
            // Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM).
            // See also Morrison and Grabowski (2007, JAS, Eq. (2))
            nlic   = std::min(nlic, qlic*1.e20);

            // set the minimum droplet number as 20/cm3.
            // nlic   = max(nlic,20.e6_r8/rho) ! sghan minimum in #/cm3
            auto tempnc = nlic/rho/1.0e6;    // #/kg --> #/cm3

            // Should be the in-cloud dropelt number calculated as nlic*rho/1.0e6_r8 ????!!!! +++mhwang
            auto pgam   = 0.0005714*(nlic*rho/1.e6) + 0.2714;
            pgam   = 1./std::pow(pgam, 2.0)-1.;
            pgam   = std::min(10., std::max(pgam,2.));
            auto laml = std::pow((pi/6.*rhow*nlic*std::tgamma(pgam+4.)/(qlic*std::tgamma(pgam+1.))), (1./3.));
            auto lammin = (pgam+1.)/50.e-6;    // in cldwat2m, CAM
            auto lammax = (pgam+1.)/2.e-6;     // if lammax is too large, this will lead to crash in
            // src/physics/rrtmg/cloud_rad_props.F90 because
            // klambda-1 can be zero in gam_liquid_lw and gam_liquid_sw
            //  and g_lambda(kmu,klambda-1) will not be defined.
            laml         = std::min(std::max(laml, lammin),lammax);
            effl(i,j)    = std::tgamma(pgam+4.)/std::tgamma(pgam+3.)/laml/2.*1.e6;
            lamcrad(i,j) = laml;
            pgamrad(i,j) = pgam;
        }
        else {
            // we chose 10. over 25, since 10 is a more reasonable value for liquid droplet. +++mhwang
            effl(i,j)    = 10.;     // in cldwat2m, CAM
            lamcrad(i,j) = 0.0;
            pgamrad(i,j) = 0.0;
        }
    });
 } // m2005_effradius

