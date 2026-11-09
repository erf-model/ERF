#include <string>
#include <memory>
#include <cmath>

#include <AMReX_Math.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF.H"
#include "ERF_NumericalConstants.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_Morrison.H"
#include "ERF_MorrisonUtils.H"
#ifdef ERF_USE_MORR_FORT
#include <ERF_Morrison_Fortran_Interface.H>
#endif

using namespace amrex;

// Components of the per-box scratch FAB holding the Morrison column variables
namespace MORRInd {
    enum  {
           lamc = 0,    // Slope parameters and PSD variables
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
           qc3dten,     // CLOUD WATER MIXING RATIO TENDENCY
           qi3dten,     // CLOUD ICE MIXING RATIO TENDENCY
           qni3dten,    // SNOW MIXING RATIO TENDENCY
           qr3dten,     // RAIN MIXING RATIO TENDENCY
           ni3dten,     // CLOUD ICE NUMBER CONCENTRATION
           ns3dten,     // SNOW NUMBER CONCENTRATION
           nr3dten,     // RAIN NUMBER CONCENTRATION
           qc3d,        // CLOUD WATER MIXING RATIO
           qi3d,        // CLOUD ICE MIXING RATIO
           qni3d,       // SNOW MIXING RATIO
           qr3d,        // RAIN MIXING RATIO
           ni3d,        // CLOUD ICE NUMBER CONCENTRATION
           ns3d,        // SNOW NUMBER CONCENTRATION
           nr3d,        // RAIN NUMBER CONCENTRATION
           t3dten,      // TEMPERATURE TENDENCY
           qv3dten,     // WATER VAPOR MIXING RATIO TENDENCY
           t3d,         // TEMPERATURE
           qv3d,        // WATER VAPOR MIXING RATIO
           pres,        // ATMOSPHERIC PRESSURE
           dzq,         // DIFFERENCE IN HEIGHT ACROSS LEVEL
           w3d,         // GRID-SCALE VERTICAL VELOCITY
           nc3d,        // CLOUD DROPLET NUMBER CONCENTRATION
           nc3dten,     // CLOUD DROPLET NUMBER CONCENTRATION TENDENCY
           qg3dten,     // GRAUPEL MIX RATIO TENDENCY
           ng3dten,     // GRAUPEL NUMB CONC TENDENCY
           qg3d,        // GRAUPEL MIX RATIO
           ng3d,        // GRAUPEL NUMBER CONC
           qgsten,      // GRAUPEL SED TEND
           qrsten,      // RAIN SED TEND
           qisten,      // CLOUD ICE SED TEND
           qnisten,     // SNOW SED TEND
           qcsten,      // CLOUD WAT SED TEND
           qrcu1d,      // RAIN FROM CUMULUS PARAMETERIZATION
           qscu1d,      // SNOW FROM CUMULUS PARAMETERIZATION
           qicu1d,      // ICE FROM CUMULUS PARAMETERIZATION
           precrt,      // TOTAL PRECIP PER TIME STEP
           snowrt,      // SNOW PER TIME STEP
           snowprt,     // TOTAL CLOUD ICE PLUS SNOW PER TIME STEP
           grplprt,     // TOTAL GRAUPEL PER TIME STEP
           rainonlyprt, // RAIN ALONE PER TIME STEP
           snowonlyprt, // SNOW ALONE (NO CLOUD ICE) PER TIME STEP
           effc,        // DROPLET EFFECTIVE RADIUS
           effi,        // CLOUD ICE EFFECTIVE RADIUS
           effs,        // SNOW EFFECTIVE RADIUS
           effr,        // RAIN EFFECTIVE RADIUS
           effg,        // GRAUPEL EFFECTIVE RADIUS
           rho,         // AIR DENSITY P/(R T)
           mu,          // AIR VISCOSITY
           ain,         // Air-density corrected fall speed coefficients
           arn,
           asn,
           acn,
           agn,
           dumi,        // Sedimenting moments, multiplied by RHO
           dumr,
           dumfni,
           dumg,
           dumfng,
           dumqs,
           dumfns,
           dumc,
           dumfnc,
           dumfnr,
           fr,          // Fall speeds of the sedimenting moments
           fi,
           fni,
           fg,
           fng,
           fs,
           fns,
           fc,
           fnc,
           fnr,
           faloutr,     // Fluxes out of the bottom of each cell
           falouti,
           faloutni,
           falouts,
           faloutns,
           faloutg,
           faloutng,
           faloutc,
           faloutnc,
           faloutnr,
           xxls,        // Latent heat of sublimation
           xxlv,        // Latent heat of vaporization
           cpm,         // Specific heat at constant pressure for moist air
           xlf,         // Latent heat of freezing
           NumInds
    };
}

namespace {

// Properties of the air at one level that the process rates need
struct MorrisonCellEnv
{
    MorrisonSaturation sat;
    Real kap;  // KAP: thermal conductivity of air
    Real dv;   // DV: diffusivity of water vapor in air
    Real sc;   // SC: Schmidt number
    Real ab;   // AB: correction to condensation rate due to latent heating
    Real abi;  // ABI: correction to deposition rate due to latent heating
};

/**
 * Process rates and tendencies at T >= 273.15 K (rain, and the melting of snow and
 * graupel), followed by the saturation adjustment of cloud water.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void morrison_warm_processes (Array4<Real> const& m, int i, int j, int k,
                              const Real rho_erf, const MorrisonCellEnv& env,
                              const MorrisonParams& mp, const Real dt, const bool do_cond)
{
    Real& qc3d  = m(i,j,k,MORRInd::qc3d);
    Real& qr3d  = m(i,j,k,MORRInd::qr3d);
    Real& qni3d = m(i,j,k,MORRInd::qni3d);
    Real& qg3d  = m(i,j,k,MORRInd::qg3d);
    Real& nc3d  = m(i,j,k,MORRInd::nc3d);
    Real& nr3d  = m(i,j,k,MORRInd::nr3d);
    Real& ns3d  = m(i,j,k,MORRInd::ns3d);
    Real& ng3d  = m(i,j,k,MORRInd::ng3d);
    Real& t3d   = m(i,j,k,MORRInd::t3d);
    const Real qv3d = m(i,j,k,MORRInd::qv3d);
    const Real pres = m(i,j,k,MORRInd::pres);
    const Real rho  = m(i,j,k,MORRInd::rho);
    const Real mu   = m(i,j,k,MORRInd::mu);
    const Real arn  = m(i,j,k,MORRInd::arn);
    const Real asn  = m(i,j,k,MORRInd::asn);
    const Real agn  = m(i,j,k,MORRInd::agn);
    const Real xxlv = m(i,j,k,MORRInd::xxlv);
    const Real xxls = m(i,j,k,MORRInd::xxls);
    const Real xlf  = m(i,j,k,MORRInd::xlf);
    const Real cpm  = m(i,j,k,MORRInd::cpm);
    Real& lamr = m(i,j,k,MORRInd::lamr);
    Real& lams = m(i,j,k,MORRInd::lams);
    Real& lamg = m(i,j,k,MORRInd::lamg);
    Real& n0r  = m(i,j,k,MORRInd::n0r);
    Real& n0s  = m(i,j,k,MORRInd::n0s);
    Real& n0g  = m(i,j,k,MORRInd::n0g);

    const MorrisonSaturation& sat = env.sat;

    // Constant droplet number (INUM = 1): convert NDCNST from cm^-3 to kg^-1
    if (mp.inum == 1) {
        nc3d = morrison_ndcnst_to_number_mixing_ratio(mp.ndcnst, rho_erf);
    }

    // Melt very small snow and graupel mixing ratios, add to rain
    if (qni3d < Real(1.0e-6)) {
        morrison_convert_all(qni3d, ns3d, qr3d, nr3d, t3d, -xlf, cpm);
    }
    if (qg3d < Real(1.0e-6)) {
        morrison_convert_all(qg3d, ng3d, qr3d, nr3d, t3d, -xlf, cpm);
    }

    if (!(qc3d < mp.qsmall && qni3d < Real(1.0e-8) && qr3d < mp.qsmall && qg3d < Real(1.0e-8)))
    {
        ns3d = amrex::max(Real(0),ns3d);
        nc3d = amrex::max(Real(0),nc3d);
        nr3d = amrex::max(Real(0),nr3d);
        ng3d = amrex::max(Real(0),ng3d);

        // Size distribution parameters
        if (qr3d >= mp.qsmall) {
            morrison_exponential_psd(qr3d, nr3d, mp.pi * mp.rhow, one/three,
                                           mp.lamminr, mp.lammaxr, lamr, n0r);
        }
        if (qc3d >= mp.qsmall) {
            morrison_cloud_droplet_psd(qc3d, nc3d, pres, t3d, mp.cons26,
                                       m(i,j,k,MORRInd::pgam), m(i,j,k,MORRInd::lamc),
                                       m(i,j,k,MORRInd::cdist1));
        }
        if (qni3d >= mp.qsmall) {
            morrison_exponential_psd(qni3d, ns3d, mp.cons1, one/mp.ds,
                                           mp.lammins, mp.lammaxs, lams, n0s);
        }
        if (qg3d >= mp.qsmall) {
            morrison_exponential_psd(qg3d, ng3d, mp.cons2, one/mp.dg,
                                           mp.lamming, mp.lammaxg, lamg, n0g);
        }

        // Autoconversion of cloud liquid water to rain
        const MorrisonAutoconversion auto_c = morrison_cloud_autoconversion(qc3d, nc3d, rho, mp.cons29, dt);
        Real prc         = auto_c.prc;
        const Real nprc  = auto_c.nprc;
        const Real nprc1 = auto_c.nprc1;

        // Collection of snow by rain above freezing (Ikawa and Saito 1991)
        // hm fix, 2/12/13: for above freezing conditions to get accelerated melting of
        // snow, we need collection of rain by snow (following Lin et al. 1983)
        Real pracs = Real(0);
        if (qr3d >= Real(1.0e-8) && qni3d >= Real(1.0e-8)) {
            const Real rhofac = morrison_fallspeed_density_factor(mp.rhosu, rho);
            const MorrisonFallSpeed vs = morrison_snow_fallspeed(asn, lams, rhofac, mp);
            const MorrisonFallSpeed vr = morrison_rain_fallspeed(arn, lamr, rhofac, mp);
            pracs = morrison_rain_collection_mass(mp.cons41, vr.mass, vs.mass, rho, n0r, n0s, lamr, lams);
        }

        // Collection of graupel by rain above freezing; all rain collected by graupel
        // above freezing is assumed shed as 1 mm drops
        Real pracg  = Real(0);
        Real npracg = Real(0);
        if (qr3d >= Real(1.0e-8) && qg3d >= Real(1.0e-8)) {
            const Real rhofac = morrison_fallspeed_density_factor(mp.rhosu, rho);
            const MorrisonFallSpeed vg = morrison_graupel_fallspeed(agn, lamg, rhofac, mp);
            const MorrisonFallSpeed vr = morrison_rain_fallspeed(arn, lamr, rhofac, mp);
            pracg = morrison_rain_collection_mass(mp.cons41, vr.mass, vg.mass, rho, n0r, n0g, lamr, lamg);

            // Number of 1 mm drops shed per second
            const Real dum = pracg/Real(5.2e-7);

            npracg = morrison_rain_collection_number(mp.cons32, vr.number, vg.number, rho, n0r, n0g, lamr, lamg);
            // hm 7/15/13, remove limit so that the number of collected drops can be
            // smaller than the number of shed drops
            npracg = npracg - dum;
        }

        // Accretion of cloud liquid water by rain, and self-collection of rain
        const MorrisonRates acc = morrison_cloud_accretion_by_rain(qc3d, qr3d, nc3d);
        Real pra         = acc.mass;
        const Real npra  = acc.number;
        const Real nragg = morrison_rain_self_collection(qr3d, nr3d, lamr, rho);

        // Evaporation of rain (Rutledge and Hobbs 1983)
        Real epsr = Real(0);
        if (qr3d >= mp.qsmall) {
            epsr = morrison_phase_relaxation(mp.pi, n0r, rho, env.dv,
                       morrison_ventilation(mp.f1r, mp.f2r, arn, rho, mu, env.sc, mp.cons9, lamr, mp.cons34));
        }
        Real pre = morrison_rain_evaporation(qv3d, sat.qvs, epsr, env.ab);

        // Melting of snow and graupel (Rutledge and Hobbs 1984), which may persist above
        // freezing; in water supersaturation the melt forms rain, otherwise it also evaporates
        Real psmlt = Real(0);
        Real evpms = Real(0);
        morrison_melting_and_evaporation(qni3d, n0s, lams, asn, mp.cons10, mp.cons35, pracs,
                                         t3d, qv3d, rho, mu, xlf, env.kap, env.dv, env.sc, env.ab,
                                         sat, mp, psmlt, evpms);
        Real pgmlt = Real(0);
        Real evpmg = Real(0);
        morrison_melting_and_evaporation(qg3d, n0g, lamg, agn, mp.cons11, mp.cons36, pracg,
                                         t3d, qv3d, rho, mu, xlf, env.kap, env.dv, env.sc, env.ab,
                                         sat, mp, pgmlt, evpmg);

        // HM, V3.2: reset PRACG and PRACS to zero, this is done because there is no
        // transfer of mass from snow and graupel to rain directly from collection above
        // freezing, it is only used for enhancement of melting and shedding
        pracg = Real(0);
        pracs = Real(0);

        // Conservation of qc
        Real dum = (prc + pra) * dt;
        if (dum > qc3d && qc3d >= mp.qsmall) {
            const Real ratio = qc3d / dum;
            prc = prc * ratio;
            pra = pra * ratio;
        }

        // Conservation of snow; no source terms for snow at T > freezing
        dum = (-psmlt - evpms + pracs) * dt;
        if (dum > qni3d && qni3d >= mp.qsmall) {
            const Real ratio = qni3d / dum;
            psmlt = psmlt * ratio;
            evpms = evpms * ratio;
            pracs = pracs * ratio;
        }

        // Conservation of graupel; no source term for graupel above freezing
        dum = (-pgmlt - evpmg + pracg) * dt;
        if (dum > qg3d && qg3d >= mp.qsmall) {
            const Real ratio = qg3d / dum;
            pgmlt = pgmlt * ratio;
            evpmg = evpmg * ratio;
            pracg = pracg * ratio;
        }

        // Conservation of qr; HM 12/13/06, added conservation of rain since PRE is negative
        dum = (-pracs - pracg - pre - pra - prc + psmlt + pgmlt) * dt;
        if (dum > qr3d && qr3d >= mp.qsmall) {
            const Real ratio = (qr3d/dt + pracs + pracg + pra + prc - psmlt - pgmlt) / (-pre);
            pre = pre * ratio;
        }

        // Tendencies
        m(i,j,k,MORRInd::qv3dten) = m(i,j,k,MORRInd::qv3dten) + (-pre - evpms - evpmg);

        m(i,j,k,MORRInd::t3dten) = m(i,j,k,MORRInd::t3dten) + (pre * xxlv +
                                    (evpms + evpmg) * xxls +
                                    (psmlt + pgmlt - pracs - pracg) * xlf) / cpm;

        m(i,j,k,MORRInd::qc3dten)  = m(i,j,k,MORRInd::qc3dten)  + (-pra - prc);
        m(i,j,k,MORRInd::qr3dten)  = m(i,j,k,MORRInd::qr3dten)  + (pre + pra + prc - psmlt - pgmlt + pracs + pracg);
        m(i,j,k,MORRInd::qni3dten) = m(i,j,k,MORRInd::qni3dten) + (psmlt + evpms - pracs);
        m(i,j,k,MORRInd::qg3dten)  = m(i,j,k,MORRInd::qg3dten)  + (pgmlt + evpmg - pracg);

        // HM, bug fix 5/12/08, npracg is subtracted from nr not ng
        m(i,j,k,MORRInd::nc3dten) = m(i,j,k,MORRInd::nc3dten) + (-npra - nprc);
        m(i,j,k,MORRInd::nr3dten) = m(i,j,k,MORRInd::nr3dten) + (nprc1 + nragg - npracg);

        // Number losses that go with evaporation and melting
        const Real nsubr  = (pre < Real(0))           ? morrison_number_loss(pre, qr3d, nr3d, dt) : Real(0);
        const Real nsmlts = (evpms + psmlt < Real(0)) ? morrison_number_loss(evpms + psmlt, qni3d, ns3d, dt) : Real(0);
        const Real nsmltr = (psmlt < Real(0))         ? morrison_number_loss(psmlt, qni3d, ns3d, dt) : Real(0);
        const Real ngmltg = (evpmg + pgmlt < Real(0)) ? morrison_number_loss(evpmg + pgmlt, qg3d, ng3d, dt) : Real(0);
        const Real ngmltr = (pgmlt < Real(0))         ? morrison_number_loss(pgmlt, qg3d, ng3d, dt) : Real(0);

        m(i,j,k,MORRInd::ns3dten) = m(i,j,k,MORRInd::ns3dten) + nsmlts;
        m(i,j,k,MORRInd::ng3dten) = m(i,j,k,MORRInd::ng3dten) + ngmltg;
        m(i,j,k,MORRInd::nr3dten) = m(i,j,k,MORRInd::nr3dten) + (nsubr - nsmltr - ngmltr);
    }

    // Saturation adjustment to condense extra vapor above water saturation
    if (do_cond) {
        const Real pcc = morrison_condensation_rate(t3d, qv3d, qc3d,
                                                    m(i,j,k,MORRInd::t3dten), m(i,j,k,MORRInd::qv3dten),
                                                    m(i,j,k,MORRInd::qc3dten), pres, xxlv, cpm, mp, dt);
        m(i,j,k,MORRInd::qv3dten) -= pcc;
        m(i,j,k,MORRInd::t3dten)  += pcc * xxlv / cpm;
        m(i,j,k,MORRInd::qc3dten) += pcc;
    }
}

/**
 * Process rates and tendencies at T < 273.15 K, followed by the saturation
 * adjustment of cloud water and the number losses of sublimation and evaporation.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void morrison_cold_processes (Array4<Real> const& m, int i, int j, int k,
                              const Real rho_erf, const MorrisonCellEnv& env,
                              const MorrisonParams& mp, const Real dt, const bool do_cond)
{
    using amrex::Math::powi;

    Real& qc3d  = m(i,j,k,MORRInd::qc3d);
    Real& qi3d  = m(i,j,k,MORRInd::qi3d);
    Real& qr3d  = m(i,j,k,MORRInd::qr3d);
    Real& qni3d = m(i,j,k,MORRInd::qni3d);
    Real& qg3d  = m(i,j,k,MORRInd::qg3d);
    Real& nc3d  = m(i,j,k,MORRInd::nc3d);
    Real& ni3d  = m(i,j,k,MORRInd::ni3d);
    Real& nr3d  = m(i,j,k,MORRInd::nr3d);
    Real& ns3d  = m(i,j,k,MORRInd::ns3d);
    Real& ng3d  = m(i,j,k,MORRInd::ng3d);
    const Real t3d  = m(i,j,k,MORRInd::t3d);
    const Real qv3d = m(i,j,k,MORRInd::qv3d);
    const Real pres = m(i,j,k,MORRInd::pres);
    const Real rho  = m(i,j,k,MORRInd::rho);
    const Real mu   = m(i,j,k,MORRInd::mu);
    const Real ain  = m(i,j,k,MORRInd::ain);
    const Real arn  = m(i,j,k,MORRInd::arn);
    const Real asn  = m(i,j,k,MORRInd::asn);
    const Real agn  = m(i,j,k,MORRInd::agn);
    const Real xxlv = m(i,j,k,MORRInd::xxlv);
    const Real xxls = m(i,j,k,MORRInd::xxls);
    const Real xlf  = m(i,j,k,MORRInd::xlf);
    const Real cpm  = m(i,j,k,MORRInd::cpm);
    Real& lamr = m(i,j,k,MORRInd::lamr);
    Real& lams = m(i,j,k,MORRInd::lams);
    Real& lami = m(i,j,k,MORRInd::lami);
    Real& lamg = m(i,j,k,MORRInd::lamg);
    Real& n0r  = m(i,j,k,MORRInd::n0r);
    Real& n0s  = m(i,j,k,MORRInd::n0s);
    Real& n0i  = m(i,j,k,MORRInd::n0i);
    Real& n0g  = m(i,j,k,MORRInd::n0g);

    const MorrisonSaturation& sat = env.sat;
    const Real dv = env.dv;
    const Real sc = env.sc;

    // Constant droplet number (INUM = 1): convert NDCNST from cm^-3 to kg^-1
    if (mp.inum == 1) {
        nc3d = morrison_ndcnst_to_number_mixing_ratio(mp.ndcnst, rho_erf);
    }

    ni3d = amrex::max(Real(0),ni3d);
    ns3d = amrex::max(Real(0),ns3d);
    nc3d = amrex::max(Real(0),nc3d);
    nr3d = amrex::max(Real(0),nr3d);
    ng3d = amrex::max(Real(0),ng3d);

    // Size distribution parameters
    if (qr3d >= mp.qsmall) {
        morrison_exponential_psd(qr3d, nr3d, mp.pi * mp.rhow, one/three,
                                 mp.lamminr, mp.lammaxr, lamr, n0r);
    }
    if (qc3d >= mp.qsmall) {
        morrison_cloud_droplet_psd(qc3d, nc3d, pres, t3d, mp.cons26,
                                   m(i,j,k,MORRInd::pgam), m(i,j,k,MORRInd::lamc),
                                   m(i,j,k,MORRInd::cdist1));
    }
    if (qni3d >= mp.qsmall) {
        morrison_exponential_psd(qni3d, ns3d, mp.cons1, one/mp.ds,
                                       mp.lammins, mp.lammaxs, lams, n0s);
    }
    if (qi3d >= mp.qsmall) {
        morrison_exponential_psd(qi3d, ni3d, mp.cons12, one/three,
                                       mp.lammini, mp.lammaxi, lami, n0i);
    }
    if (qg3d >= mp.qsmall) {
        morrison_exponential_psd(qg3d, ng3d, mp.cons2, one/mp.dg,
                                       mp.lamming, mp.lammaxg, lamg, n0g);
    }

    // Freezing of cloud droplets, only allowed below -4C
    const MorrisonRates frz_c = morrison_cloud_droplet_freezing(
        qc3d, nc3d, t3d, pres, mu, m(i,j,k,MORRInd::cdist1), m(i,j,k,MORRInd::pgam),
        m(i,j,k,MORRInd::lamc), mp, dt);
    Real mnuccc = frz_c.mass;
    Real nnuccc = frz_c.number;

    // Autoconversion of cloud liquid water to rain
    const MorrisonAutoconversion auto_c = morrison_cloud_autoconversion(qc3d, nc3d, rho, mp.cons29, dt);
    Real prc         = auto_c.prc;
    const Real nprc  = auto_c.nprc;
    const Real nprc1 = auto_c.nprc1;

    // Snow aggregation
    const Real nsagg = morrison_snow_aggregation(qni3d, ns3d, asn, rho, mp);

    // Accretion of cloud droplets onto snow, graupel and cloud ice (continuous
    // collection with a simple gravitational kernel)
    Real psacws  = Real(0);
    Real npsacws = Real(0);
    if (qni3d >= Real(1.0e-8) && qc3d >= mp.qsmall) {
        psacws  = morrison_continuous_collection(mp.cons13, asn, qc3d, rho, n0s, lams, mp.bs);
        npsacws = morrison_continuous_collection(mp.cons13, asn, nc3d, rho, n0s, lams, mp.bs);
    }
    Real psacwg  = Real(0);
    Real npsacwg = Real(0);
    if (qg3d >= Real(1.0e-8) && qc3d >= mp.qsmall) {
        psacwg  = morrison_continuous_collection(mp.cons14, agn, qc3d, rho, n0g, lamg, mp.bg);
        npsacwg = morrison_continuous_collection(mp.cons14, agn, nc3d, rho, n0g, lamg, mp.bg);
    }
    // hm, add 12/13/06: cloud ice collecting droplets, assume that cloud ice mean diam
    // > 100 micron before riming can occur; assume that rime collected on cloud ice
    // does not lead to Hallet-Mossop splintering
    Real psacwi  = Real(0);
    Real npsacwi = Real(0);
    if (qi3d >= Real(1.0e-8) && qc3d >= mp.qsmall) {
        // Size dependent collection efficiency based on Stokes law (Thompson et al. 2004)
        if (one / lami >= Real(100.0e-6)) {
            psacwi  = morrison_continuous_collection(mp.cons16, ain, qc3d, rho, n0i, lami, mp.bi);
            npsacwi = morrison_continuous_collection(mp.cons16, ain, nc3d, rho, n0i, lami, mp.bi);
        }
    }

    // Accretion of rain water by snow (Ikawa and Saito 1991, used by Reisner et al. 1998),
    // and the collection of snow by rain needed for the graupel conversion
    Real pracs  = Real(0);
    Real npracs = Real(0);
    Real psacr  = Real(0);
    if (qr3d >= Real(1.0e-8) && qni3d >= Real(1.0e-8)) {
        const Real rhofac = morrison_fallspeed_density_factor(mp.rhosu, rho);
        const MorrisonFallSpeed vs = morrison_snow_fallspeed(asn, lams, rhofac, mp);
        const MorrisonFallSpeed vr = morrison_rain_fallspeed(arn, lamr, rhofac, mp);

        pracs  = morrison_rain_collection_mass(mp.cons41, vr.mass, vs.mass, rho, n0r, n0s, lamr, lams);
        npracs = morrison_rain_collection_number(mp.cons32, vr.number, vs.number, rho, n0r, n0s, lamr, lams);

        // Make sure PRACS doesn't exceed total rain mixing ratio, as this may otherwise
        // result in too much transfer of water during rime-splintering
        pracs = std::min(pracs, qr3d / dt);

        // Only calculate if snow and rain mixing ratios exceed 0.1 g/kg
        if (qni3d >= Real(0.1e-3) && qr3d >= Real(0.1e-3)) {
            psacr = morrison_rain_collection_mass(mp.cons31, vr.mass, vs.mass, rho, n0r, n0s, lams, lamr);
        }
    }

    // Collection of rainwater by graupel (Ikawa and Saito 1990, used by Reisner et al 1998)
    Real pracg  = Real(0);
    Real npracg = Real(0);
    if (qr3d >= Real(1.0e-8) && qg3d >= Real(1.0e-8)) {
        const Real rhofac = morrison_fallspeed_density_factor(mp.rhosu, rho);
        const MorrisonFallSpeed vg = morrison_graupel_fallspeed(agn, lamg, rhofac, mp);
        const MorrisonFallSpeed vr = morrison_rain_fallspeed(arn, lamr, rhofac, mp);

        pracg  = morrison_rain_collection_mass(mp.cons41, vr.mass, vg.mass, rho, n0r, n0g, lamr, lamg);
        npracg = morrison_rain_collection_number(mp.cons32, vr.number, vg.number, rho, n0r, n0g, lamr, lamg);

        // Make sure PRACG doesn't exceed total rain mixing ratio
        pracg = std::min(pracg, qr3d / dt);
    }

    // Rime-splintering on snow and on graupel
    Real nmults = Real(0), qmults = Real(0), nmultr = Real(0), qmultr = Real(0);
    morrison_rime_splintering(qni3d, qc3d, qr3d, t3d, mp.mmult, psacws, pracs,
                              nmults, qmults, nmultr, qmultr);
    Real nmultg = Real(0), qmultg = Real(0), nmultrg = Real(0), qmultrg = Real(0);
    morrison_rime_splintering(qg3d, qc3d, qr3d, t3d, mp.mmult, psacwg, pracg,
                              nmultg, qmultg, nmultrg, qmultrg);

    // Conversion of rimed cloud water onto snow to graupel/hail
    Real pgsacw = Real(0);
    Real nscng  = Real(0);
    if (psacws > Real(0)) {
        // Only allow conversion if QNI > 0.1 and QC > 0.5 g/kg following Rutledge and Hobbs (1984)
        if (qni3d >= Real(0.1e-3) && qc3d >= Real(0.5e-3)) {
            // Portion of riming converted to graupel (Reisner et al. 1998, originally IS1991)
            pgsacw = std::min(psacws, mp.cons17 * dt * n0s * qc3d * qc3d * asn * asn /
                                      (rho * std::pow(lams, (Real(2) * mp.bs + Real(2)))));

            // Mix rat converted into graupel as embryo (Reisner et al. 1998, orig M1990)
            const Real dum = std::max(mp.rhosn / (mp.rhog - mp.rhosn) * pgsacw, Real(0));

            // Number concentration of embryo graupel from riming of snow, limited to the snow number
            nscng = dum / mp.mg0 * rho;
            nscng = std::min(nscng, ns3d / dt);

            // Portion of riming left for snow
            psacws = psacws - pgsacw;
        }
    }

    // Conversion of rimed rainwater onto snow converted to graupel
    Real pgracs = Real(0);
    Real ngracs = Real(0);
    if (pracs > Real(0)) {
        // Only allow conversion if QNI > 0.1 and QR > 0.1 g/kg following Rutledge and Hobbs (1984)
        if (qni3d >= Real(0.1e-3) && qr3d >= Real(0.1e-3)) {
            // Portion of collected rainwater converted to graupel (Reisner et al. 1998)
            Real dum = mp.cons18 * powi<3>(Real(4.0) / lams) * powi<3>(Real(4.0) / lams) /
                (mp.cons18 * powi<3>(Real(4.0) / lams) * powi<3>(Real(4.0) / lams) +
                 mp.cons19 * powi<3>(Real(4.0) / lamr) * powi<3>(Real(4.0) / lamr));
            dum = std::min(dum, Real(1));
            dum = std::max(dum, Real(0));

            pgracs = (one - dum) * pracs;
            ngracs = (one - dum) * npracs;

            // Limit max number converted to min of either rain or snow number concentration
            ngracs = std::min(ngracs, nr3d / dt);
            ngracs = std::min(ngracs, ns3d / dt);

            // Amount left for snow production
            pracs = pracs - pgracs;
            npracs = npracs - ngracs;

            // Conversion to graupel due to collection of snow by rain
            psacr = psacr * (one - dum);
        }
    }

    // Freezing of rain drops, allowed below -4C
    const MorrisonRates frz_r = morrison_rain_freezing(qr3d, nr3d, lamr, t3d, mp, dt);
    Real mnuccr = frz_r.mass;
    Real nnuccr = frz_r.number;

    // Accretion of cloud liquid water by rain, and self-collection of rain
    const MorrisonRates acc = morrison_cloud_accretion_by_rain(qc3d, qr3d, nc3d);
    Real pra         = acc.mass;
    const Real npra  = acc.number;
    const Real nragg = morrison_rain_self_collection(qr3d, nr3d, lamr, rho);

    // Autoconversion of cloud ice to snow
    const MorrisonRates auto_i = morrison_ice_autoconversion(qi3d, ni3d, n0i, lami, qv3d, rho,
                                                                  dv, env.abi, sat, mp, dt);
    Real prci        = auto_i.mass;
    const Real nprci = auto_i.number;

    // Accretion of cloud ice by snow, assuming VS >> VI and DS >> DI for continuous collection
    Real prai  = Real(0);
    Real nprai = Real(0);
    if (qni3d >= Real(1.0e-8) && qi3d >= mp.qsmall) {
        prai  = morrison_continuous_collection(mp.cons23, asn, qi3d, rho, n0s, lams, mp.bs);
        nprai = morrison_continuous_collection(mp.cons23, asn, ni3d, rho, n0s, lams, mp.bs);
        nprai = std::min(nprai, ni3d / dt);
    }

    // hm, add 12/13/06, collision of rain and ice to produce snow or graupel; allow graupel
    // formation from rain-ice collisions only if rain mixing ratio > 0.1 g/kg, otherwise add to snow
    Real niacr = Real(0), piacr = Real(0), praci = Real(0);
    Real niacrs = Real(0), piacrs = Real(0), pracis = Real(0);
    if (qr3d >= Real(1.0e-8) && qi3d >= Real(1.0e-8) && t3d <= Real(273.15)) {
        const MorrisonRainIceCollision c = morrison_rain_ice_collision(qi3d, ni3d, nr3d, n0r, lamr, arn,
                                                                       rho, mp, dt);
        if (qr3d >= Real(0.1e-3)) {
            niacr = c.niacr;
            piacr = c.piacr;
            praci = c.praci;
        } else {
            niacrs = c.niacr;
            piacrs = c.piacr;
            pracis = c.praci;
        }
    }

    // Nucleation of cloud ice from homogeneous and heterogeneous freezing on aerosol
    const MorrisonRates nuc = morrison_ice_nucleation(t3d, rho, ni3d, ns3d, ng3d, sat, mp, dt);
    Real mnuccd = nuc.mass;
    Real nnuccd = nuc.number;

    // Evaporation/sublimation/deposition terms for qi, qni, qg and qr; no ventilation for cloud ice
    Real epsi = Real(0);
    if (qi3d >= mp.qsmall) {
        epsi = Real(2) * mp.pi * n0i * rho * dv / (lami * lami);
    }
    // Ventilation for snow, graupel and rain
    Real epss = Real(0);
    if (qni3d >= mp.qsmall) {
        epss = morrison_phase_relaxation(mp.pi, n0s, rho, dv,
                   morrison_ventilation(mp.f1s, mp.f2s, asn, rho, mu, sc, mp.cons10, lams, mp.cons35));
    }
    Real epsg = Real(0);
    if (qg3d >= mp.qsmall) {
        epsg = morrison_phase_relaxation(mp.pi, n0g, rho, dv,
                   morrison_ventilation(mp.f1s, mp.f2s, agn, rho, mu, sc, mp.cons11, lamg, mp.cons36));
    }
    Real epsr = Real(0);
    if (qr3d >= mp.qsmall) {
        epsr = morrison_phase_relaxation(mp.pi, n0r, rho, dv,
                   morrison_ventilation(mp.f1r, mp.f2r, arn, rho, mu, sc, mp.cons9, lamr, mp.cons34));
    }

    // Only include the region of the ice size distribution < DCS; DUM is the fraction
    // of D*N(D) < DCS (Harrington et al. 1995)
    Real prd = Real(0);
    Real dum = Real(0);
    if (qi3d >= mp.qsmall) {
        dum = (one - std::exp(-lami * mp.dcs) * (one + lami * mp.dcs));
        prd = epsi * (qv3d - sat.qvi) / env.abi * dum;
    }

    // Add deposition in the tail of the ice size distribution to snow if snow is present,
    // otherwise to cloud ice
    Real prds = Real(0);
    if (qni3d >= mp.qsmall) {
        prds = epss * (qv3d - sat.qvi) / env.abi +
            epsi * (qv3d - sat.qvi) / env.abi * (one - dum);
    } else {
        prd = prd + epsi * (qv3d - sat.qvi) / env.abi * (one - dum);
    }

    // Vapor deposition on graupel
    Real prdg = epsg * (qv3d - sat.qvi) / env.abi;

    // No condensation onto rain, only evaporation
    Real pre = morrison_rain_evaporation(qv3d, sat.qvs, epsr, env.ab);

    // Make sure not pushed into ice supersaturation/subsaturation (Reisner 2 scheme)
    dum = (qv3d - sat.qvi) / dt;

    const Real fudgef = Real(0.9999);
    const Real sum_dep = prd + prds + mnuccd + prdg;

    if ((dum > Real(0) && sum_dep > dum * fudgef) ||
        (dum < Real(0) && sum_dep < dum * fudgef)) {
        mnuccd = fudgef * mnuccd * dum / sum_dep;
        prd = fudgef * prd * dum / sum_dep;
        prds = fudgef * prds * dum / sum_dep;
        prdg = fudgef * prdg * dum / sum_dep;
    }

    // If cloud ice/snow/graupel vapor deposition is negative, assign it to sublimation
    Real eprd = Real(0);
    if (prd < Real(0)) {
        eprd = prd;
        prd = Real(0);
    }
    Real eprds = Real(0);
    if (prds < Real(0)) {
        eprds = prds;
        prds = Real(0);
    }
    Real eprdg = Real(0);
    if (prdg < Real(0)) {
        eprdg = prdg;
        prdg = Real(0);
    }

    // Conservation of water. This is adopted loosely from the Reisner code, but here only
    // processes that are negative are adjusted, rather than all processes. Below QSMALL
    // there is no depletion of water through microphysical processes, so conservation is
    // skipped. The check is not applied to number; each process that is a source of number
    // is limited so it can't exceed the number of the depleted species within the time step.

    // ****SENSITIVITY - NO ICE
    if (mp.iliq == 1) {
        mnuccc = Real(0);
        nnuccc = Real(0);
        mnuccr = Real(0);
        nnuccr = Real(0);
        mnuccd = Real(0);
        nnuccd = Real(0);
    }

    // ****SENSITIVITY - NO GRAUPEL
    Real nsubg = Real(0);
    if (mp.igraup == 1) {
        pracg   = Real(0);
        psacr   = Real(0);
        psacwg  = Real(0);
        prdg    = Real(0);
        eprdg   = Real(0);
        npracg  = Real(0);
        npsacwg = Real(0);
        nscng   = Real(0);
        ngracs  = Real(0);

        // fix 053011
        piacrs = piacrs + piacr;
        piacr = Real(0);

        // fix 070713
        pracis = pracis + praci;
        praci = Real(0);
        psacws = psacws + pgsacw;
        pgsacw = Real(0);
        pracs = pracs + pgracs;
        pgracs = Real(0);
    }

    // Conservation of qc
    dum = (prc + pra + mnuccc + psacws + psacwi + qmults + psacwg + pgsacw + qmultg) * dt;
    if (dum > qc3d && qc3d >= mp.qsmall) {
        const Real ratio = qc3d / dum;
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

    // Conservation of qi
    dum = (-prd - mnuccc + prci + prai - qmults - qmultg - qmultr - qmultrg
           - mnuccd + praci + pracis - eprd - psacwi) * dt;
    if (dum > qi3d && qi3d >= mp.qsmall) {
        const Real ratio = (qi3d / dt + prd + mnuccc + qmults + qmultg + qmultr + qmultrg +
                            mnuccd + psacwi) /
            (prci + prai + praci + pracis - eprd);
        prci = prci * ratio;
        prai = prai * ratio;
        praci = praci * ratio;
        pracis = pracis * ratio;
        eprd = eprd * ratio;
    }

    // Conservation of qr
    dum = ((pracs - pre) + (qmultr + qmultrg - prc) + (mnuccr - pra) +
           piacr + piacrs + pgracs + pracg) * dt;
    if (dum > qr3d && qr3d >= mp.qsmall) {
        const Real ratio = (qr3d / dt + prc + pra) /
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

    // Conservation of qni; for no graupel, the freezing of rain feeds snow
    if (mp.igraup == 0) {
        dum = (-prds - psacws - prai - prci - pracs - eprds + psacr - piacrs - pracis) * dt;
        if (dum > qni3d && qni3d >= mp.qsmall) {
            const Real ratio = (qni3d / dt + prds + psacws + prai + prci + pracs + piacrs + pracis) /
                (-eprds + psacr);
            eprds = eprds * ratio;
            psacr = psacr * ratio;
        }
    } else if (mp.igraup == 1) {
        dum = (-prds - psacws - prai - prci - pracs - eprds + psacr - piacrs - pracis - mnuccr) * dt;
        if (dum > qni3d && qni3d >= mp.qsmall) {
            const Real ratio = (qni3d / dt + prds + psacws + prai + prci + pracs + piacrs + pracis + mnuccr) /
                (-eprds + psacr);
            eprds = eprds * ratio;
            psacr = psacr * ratio;
        }
    }

    // Conservation of qg
    dum = (-psacwg - pracg - pgsacw - pgracs - prdg - mnuccr - eprdg - piacr - praci - psacr) * dt;
    if (dum > qg3d && qg3d >= mp.qsmall) {
        const Real ratio = (qg3d / dt + psacwg + pracg + pgsacw + pgracs + prdg + mnuccr + psacr +
                            piacr + praci) / (-eprdg);
        eprdg = eprdg * ratio;
    }

    // Tendencies; PCC is zero until the saturation adjustment below
    const Real pcc = Real(0);
    m(i,j,k,MORRInd::qv3dten) = m(i,j,k,MORRInd::qv3dten) + (-pre - prd - prds - mnuccd - eprd - eprds - prdg - eprdg);

    // bug fix hm, 3/1/11, include piacr and piacrs
    m(i,j,k,MORRInd::t3dten) = m(i,j,k,MORRInd::t3dten) + (pre * xxlv +
                                (prd + prds + mnuccd + eprd + eprds + prdg + eprdg) * xxls +
                                (psacws + psacwi + mnuccc + mnuccr + qmults + qmultg + qmultr + qmultrg + pracs +
                                 psacwg + pracg + pgsacw + pgracs + piacr + piacrs) * xlf) / cpm;

    m(i,j,k,MORRInd::qc3dten) = m(i,j,k,MORRInd::qc3dten) +
        (-pra - prc - mnuccc + pcc -
         psacws - psacwi - qmults - qmultg - psacwg - pgsacw);

    m(i,j,k,MORRInd::qi3dten) = m(i,j,k,MORRInd::qi3dten) +
        (prd + eprd + psacwi + mnuccc - prci -
         prai + qmults + qmultg + qmultr + qmultrg + mnuccd - praci - pracis);

    m(i,j,k,MORRInd::qr3dten) = m(i,j,k,MORRInd::qr3dten) +
        (pre + pra + prc - pracs - mnuccr - qmultr - qmultrg -
         piacr - piacrs - pracg - pgracs);

    if (mp.igraup == 0) {
        m(i,j,k,MORRInd::qni3dten) = m(i,j,k,MORRInd::qni3dten) +
            (prai + psacws + prds + pracs + prci + eprds - psacr + piacrs + pracis);

        m(i,j,k,MORRInd::ns3dten) = m(i,j,k,MORRInd::ns3dten) + (nsagg + nprci - nscng - ngracs + niacrs);

        m(i,j,k,MORRInd::qg3dten) = m(i,j,k,MORRInd::qg3dten) + (pracg + psacwg + pgsacw + pgracs +
                                     prdg + eprdg + mnuccr + piacr + praci + psacr);

        m(i,j,k,MORRInd::ng3dten) = m(i,j,k,MORRInd::ng3dten) + (nscng + ngracs + nnuccr + niacr);
    } else if (mp.igraup == 1) {
        // For no graupel, need to include freezing of rain for snow
        m(i,j,k,MORRInd::qni3dten) = m(i,j,k,MORRInd::qni3dten) +
            (prai + psacws + prds + pracs + prci + eprds - psacr + piacrs + pracis + mnuccr);

        m(i,j,k,MORRInd::ns3dten) = m(i,j,k,MORRInd::ns3dten) + (nsagg + nprci - nscng - ngracs + niacrs + nnuccr);
    }

    m(i,j,k,MORRInd::nc3dten) = m(i,j,k,MORRInd::nc3dten) + (-nnuccc - npsacws -
                                 npra - nprc - npsacwi - npsacwg);

    m(i,j,k,MORRInd::ni3dten) = m(i,j,k,MORRInd::ni3dten) +
        (nnuccc - nprci - nprai + nmults + nmultg + nmultr + nmultrg +
         nnuccd - niacr - niacrs);

    m(i,j,k,MORRInd::nr3dten) = m(i,j,k,MORRInd::nr3dten) + (nprc1 - npracs - nnuccr +
                                 nragg - niacr - niacrs - npracg - ngracs);

    // Saturation adjustment to condense extra vapor above water saturation
    if (do_cond) {
        const Real pcc_adj = morrison_condensation_rate(t3d, qv3d, qc3d,
                                                        m(i,j,k,MORRInd::t3dten), m(i,j,k,MORRInd::qv3dten),
                                                        m(i,j,k,MORRInd::qc3dten), pres, xxlv, cpm, mp, dt);
        m(i,j,k,MORRInd::qv3dten) = m(i,j,k,MORRInd::qv3dten) - pcc_adj;
        m(i,j,k,MORRInd::t3dten)  = m(i,j,k,MORRInd::t3dten) + pcc_adj * xxlv / cpm;
        m(i,j,k,MORRInd::qc3dten) = m(i,j,k,MORRInd::qc3dten) + pcc_adj;
    }

    // Sublimate, melt or evaporate number concentration, assuming a 1:1 ratio between
    // the loss of mass and the loss of number
    const Real nsubi = (eprd < Real(0))  ? morrison_number_loss(eprd, qi3d, ni3d, dt) : Real(0);
    const Real nsubs = (eprds < Real(0)) ? morrison_number_loss(eprds, qni3d, ns3d, dt) : Real(0);
    const Real nsubr = (pre < Real(0))   ? morrison_number_loss(pre, qr3d, nr3d, dt) : Real(0);
    if (eprdg < Real(0)) {
        nsubg = morrison_number_loss(eprdg, qg3d, ng3d, dt);
    }

    m(i,j,k,MORRInd::ni3dten) = m(i,j,k,MORRInd::ni3dten) + nsubi;
    m(i,j,k,MORRInd::ns3dten) = m(i,j,k,MORRInd::ns3dten) + nsubs;
    m(i,j,k,MORRInd::ng3dten) = m(i,j,k,MORRInd::ng3dten) + nsubg;
    m(i,j,k,MORRInd::nr3dten) = m(i,j,k,MORRInd::nr3dten) + nsubr;
}

/**
 * Microphysical process rates of one level (the first K loop of MORR_TWO_MOMENT_MICRO).
 * Returns false when the level has no hydrometeors and is subsaturated, in which case
 * the processes are skipped.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool morrison_cell_microphysics (Array4<Real> const& m, int i, int j, int k,
                                 const Real rho_erf, const MorrisonParams& mp,
                                 const Real dt, const bool do_cond)
{
    using amrex::Math::powi;

    Real& t3d   = m(i,j,k,MORRInd::t3d);
    Real& qv3d  = m(i,j,k,MORRInd::qv3d);
    Real& qc3d  = m(i,j,k,MORRInd::qc3d);
    Real& qi3d  = m(i,j,k,MORRInd::qi3d);
    Real& qr3d  = m(i,j,k,MORRInd::qr3d);
    Real& qni3d = m(i,j,k,MORRInd::qni3d);
    Real& qg3d  = m(i,j,k,MORRInd::qg3d);
    Real& rho   = m(i,j,k,MORRInd::rho);
    const Real pres = m(i,j,k,MORRInd::pres);

    m(i,j,k,MORRInd::nc3dten) = Real(0);

    const Real xxlv = morrison_latent_heat_vaporization(t3d);
    const Real xxls = morrison_latent_heat_sublimation(t3d);
    const Real cpm  = morrison_cpm(qv3d);
    m(i,j,k,MORRInd::xxlv) = xxlv;
    m(i,j,k,MORRInd::xxls) = xxls;
    m(i,j,k,MORRInd::cpm)  = cpm;

    MorrisonCellEnv env;
    env.sat = morrison_saturation(t3d, pres, qv3d, mp.ep_2);
    const MorrisonSaturation& sat = env.sat;

    rho = pres / (mp.R * t3d);

    // Add number concentration due to cumulus tendency. Assume N0 associated with cumulus
    // parameterization rain is 10^7 m^-4 and with snow 2 x 10^7 m^-4; for detrained cloud
    // ice, assume a mean volume diameter of 80 microns.
    if (m(i,j,k,MORRInd::qrcu1d) >= Real(1.0e-10)) {
        m(i,j,k,MORRInd::nr3d) += Real(1.8e5) * std::pow(m(i,j,k,MORRInd::qrcu1d) * dt /
                                     (mp.pi * mp.rhow * powi<3>(rho)), fourth);
    }
    if (m(i,j,k,MORRInd::qscu1d) >= Real(1.0e-10)) {
        m(i,j,k,MORRInd::ns3d) += Real(3.e5) * std::pow(m(i,j,k,MORRInd::qscu1d) * dt /
                                     (mp.cons1 * powi<3>(rho)), one / (mp.ds + one));
    }
    if (m(i,j,k,MORRInd::qicu1d) >= Real(1.0e-10)) {
        m(i,j,k,MORRInd::ni3d) += m(i,j,k,MORRInd::qicu1d) * dt / (mp.ci * std::pow(Real(80.0e-6), mp.di));
    }

    morrison_evaporate_small_at_subsaturation(sat, xxlv, xxls, cpm, qv3d, t3d,
                                              qc3d, qr3d, qi3d, qni3d, qg3d);

    // Heat of fusion
    const Real xlf = xxls - xxlv;
    m(i,j,k,MORRInd::xlf) = xlf;

    morrison_zero_below_qsmall(qc3d,  m(i,j,k,MORRInd::nc3d), m(i,j,k,MORRInd::effc), mp.qsmall);
    morrison_zero_below_qsmall(qr3d,  m(i,j,k,MORRInd::nr3d), m(i,j,k,MORRInd::effr), mp.qsmall);
    morrison_zero_below_qsmall(qi3d,  m(i,j,k,MORRInd::ni3d), m(i,j,k,MORRInd::effi), mp.qsmall);
    morrison_zero_below_qsmall(qni3d, m(i,j,k,MORRInd::ns3d), m(i,j,k,MORRInd::effs), mp.qsmall);
    morrison_zero_below_qsmall(qg3d,  m(i,j,k,MORRInd::ng3d), m(i,j,k,MORRInd::effg), mp.qsmall);

    // Initialize sedimentation tendencies for mixing ratio
    m(i,j,k,MORRInd::qrsten)  = Real(0);
    m(i,j,k,MORRInd::qisten)  = Real(0);
    m(i,j,k,MORRInd::qnisten) = Real(0);
    m(i,j,k,MORRInd::qcsten)  = Real(0);
    m(i,j,k,MORRInd::qgsten)  = Real(0);

    // Microphysics parameters varying in time/height
    const Real mu = morrison_air_viscosity(t3d);
    m(i,j,k,MORRInd::mu) = mu;

    const MorrisonFallSpeedCoefficients fc = morrison_fallspeed_coefficients(rho, mu, mp);
    m(i,j,k,MORRInd::ain) = fc.ain;
    m(i,j,k,MORRInd::arn) = fc.arn;
    m(i,j,k,MORRInd::asn) = fc.asn;
    m(i,j,k,MORRInd::acn) = fc.acn;
    m(i,j,k,MORRInd::agn) = fc.agn;

    // hm 4/7/09 bug fix, initialize LAMI to prevent later division by zero
    m(i,j,k,MORRInd::lami) = Real(0);

    // If there is no cloud/precip water, and if subsaturated, then skip microphysics for this level
    if (qc3d < mp.qsmall && qi3d < mp.qsmall && qni3d < mp.qsmall && qr3d < mp.qsmall && qg3d < mp.qsmall) {
        if ((t3d < Real(273.15) && sat.qvqvsi < Real(0.999)) || (t3d >= Real(273.15) && sat.qvqvs < Real(0.999))) {
            return false;
        }
    }

    env.kap = morrison_thermal_conductivity(mu);
    env.dv  = morrison_vapor_diffusivity(t3d, pres);
    env.sc  = mu / (rho * env.dv);

    // Psychrometric corrections: rate of change of saturation mixing ratio with temperature
    const Real dum = (mp.Rv * powi<2>(t3d));
    const Real dqsdt  = xxlv * sat.qvs / dum;
    const Real dqsidt = xxls * sat.qvi / dum;
    env.abi = one + dqsidt * xxls / cpm;
    env.ab  = one + dqsdt * xxlv / cpm;

    if (t3d >= Real(273.15)) {
        morrison_warm_processes(m, i, j, k, rho_erf, env, mp, dt, do_cond);
    } else {
        morrison_cold_processes(m, i, j, k, rho_erf, env, mp, dt, do_cond);
    }
    return true;
}

/**
 * Sedimentation of column (i,j) over [klo,khi] (Reisner et al. 1998). The fallout is
 * computed on split time steps so that the Courant number stays below one; the flux
 * out of the bottom of each cell is added to the face fluxes sed.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void morrison_column_sedimentation (Array4<Real> const& m, Array4<Real> const& sed,
                                    int i, int j, int klo, int khi,
                                    const MorrisonParams& mp, const Real dt)
{
    // The sedimenting moments: their working value, fall speed, flux out of the cell
    // bottom, the tendency they feed and their component of the face fluxes
    constexpr int nmom = 10;
    const int dum_c[nmom]    = {MORRInd::dumr,    MORRInd::dumi,    MORRInd::dumfni,   MORRInd::dumqs,
                                MORRInd::dumfns,  MORRInd::dumfnr,  MORRInd::dumc,     MORRInd::dumfnc,
                                MORRInd::dumg,    MORRInd::dumfng};
    const int f_c[nmom]      = {MORRInd::fr,      MORRInd::fi,      MORRInd::fni,      MORRInd::fs,
                                MORRInd::fns,     MORRInd::fnr,     MORRInd::fc,       MORRInd::fnc,
                                MORRInd::fg,      MORRInd::fng};
    const int falout_c[nmom] = {MORRInd::faloutr, MORRInd::falouti, MORRInd::faloutni, MORRInd::falouts,
                                MORRInd::faloutns,MORRInd::faloutnr,MORRInd::faloutc,  MORRInd::faloutnc,
                                MORRInd::faloutg, MORRInd::faloutng};
    const int ten_c[nmom]    = {MORRInd::qrsten,  MORRInd::qisten,  MORRInd::ni3dten,  MORRInd::qnisten,
                                MORRInd::ns3dten, MORRInd::nr3dten, MORRInd::qcsten,   MORRInd::nc3dten,
                                MORRInd::qgsten,  MORRInd::ng3dten};
    const int sed_c[nmom]    = {MorrSed::qr,      MorrSed::qi,      MorrSed::ni,       MorrSed::qs,
                                MorrSed::ns,      MorrSed::nr,      MorrSed::qc,       MorrSed::nc,
                                MorrSed::qg,      MorrSed::ng};
    const int q_c[nmom]      = {MORRInd::qr3d,    MORRInd::qi3d,    MORRInd::ni3d,     MORRInd::qni3d,
                                MORRInd::ns3d,    MORRInd::nr3d,    MORRInd::qc3d,     MORRInd::nc3d,
                                MORRInd::qg3d,    MORRInd::ng3d};
    const int qten_c[nmom]   = {MORRInd::qr3dten, MORRInd::qi3dten, MORRInd::ni3dten,  MORRInd::qni3dten,
                                MORRInd::ns3dten, MORRInd::nr3dten, MORRInd::qc3dten,  MORRInd::nc3dten,
                                MORRInd::qg3dten, MORRInd::ng3dten};

    int nstep = 1;

    // Fall speeds, top to bottom
    for (int k = khi; k >= klo; k--) {
        // Update prognostic variables with tendencies
        for (int n = 0; n < nmom; ++n) {
            m(i,j,k,dum_c[n]) = m(i,j,k,q_c[n]) + m(i,j,k,qten_c[n]) * dt;
        }

        // Switch for constant droplet number
        if (mp.inum == 1) {
            m(i,j,k,MORRInd::dumfnc) = m(i,j,k,MORRInd::nc3d);
        }

        // Make sure number concentrations are positive
        Real& dumfni = m(i,j,k,MORRInd::dumfni);
        Real& dumfns = m(i,j,k,MORRInd::dumfns);
        Real& dumfnc = m(i,j,k,MORRInd::dumfnc);
        Real& dumfnr = m(i,j,k,MORRInd::dumfnr);
        Real& dumfng = m(i,j,k,MORRInd::dumfng);
        dumfni = amrex::max(Real(0), dumfni);
        dumfns = amrex::max(Real(0), dumfns);
        dumfnc = amrex::max(Real(0), dumfnc);
        dumfnr = amrex::max(Real(0), dumfnr);
        dumfng = amrex::max(Real(0), dumfng);

        const Real dumi  = m(i,j,k,MORRInd::dumi);
        const Real dumr  = m(i,j,k,MORRInd::dumr);
        const Real dumc  = m(i,j,k,MORRInd::dumc);
        const Real dumqs = m(i,j,k,MORRInd::dumqs);
        const Real dumg  = m(i,j,k,MORRInd::dumg);
        const Real rho   = m(i,j,k,MORRInd::rho);

        const Real rhofac = morrison_fallspeed_density_factor(mp.rhosu, rho);

        // Mass- and number-weighted fall speeds from the slopes of the updated moments
        MorrisonFallSpeed vc{Real(0), Real(0)};
        if (dumc >= mp.qsmall) {
            Real& pgam = m(i,j,k,MORRInd::pgam);
            pgam = morrison_cloud_droplet_pgam(m(i,j,k,MORRInd::nc3d), m(i,j,k,MORRInd::pres), m(i,j,k,MORRInd::t3d));
            Real dlamc = morrison_cloud_droplet_lambda(dumc, dumfnc, mp.cons26,
                                                       gamma_function(pgam + one), gamma_function(pgam + Real(4.0)));
            dlamc = amrex::max(dlamc, morrison_cloud_droplet_lammin(pgam));
            dlamc = amrex::min(dlamc, morrison_cloud_droplet_lammax(pgam));
            vc = morrison_cloud_droplet_fallspeed(m(i,j,k,MORRInd::acn), dlamc, pgam, mp);
        }
        MorrisonFallSpeed vi{Real(0), Real(0)};
        if (dumi >= mp.qsmall) {
            const Real dlami = morrison_limited_slope(dumi, dumfni, mp.cons12, one/mp.di, mp.lammini, mp.lammaxi);
            vi = morrison_ice_fallspeed(m(i,j,k,MORRInd::ain), dlami, rho, mp);
        }
        MorrisonFallSpeed vr{Real(0), Real(0)};
        if (dumr >= mp.qsmall) {
            const Real dlamr = morrison_limited_slope(dumr, dumfnr, mp.pi * mp.rhow, one/three, mp.lamminr, mp.lammaxr);
            vr = morrison_rain_fallspeed(m(i,j,k,MORRInd::arn), dlamr, rhofac, mp);
        }
        MorrisonFallSpeed vs{Real(0), Real(0)};
        if (dumqs >= mp.qsmall) {
            const Real dlams = morrison_limited_slope(dumqs, dumfns, mp.cons1, one/mp.ds, mp.lammins, mp.lammaxs);
            vs = morrison_snow_fallspeed(m(i,j,k,MORRInd::asn), dlams, rhofac, mp);
        }
        MorrisonFallSpeed vg{Real(0), Real(0)};
        if (dumg >= mp.qsmall) {
            const Real dlamg = morrison_limited_slope(dumg, dumfng, mp.cons2, one/mp.dg, mp.lamming, mp.lammaxg);
            vg = morrison_graupel_fallspeed(m(i,j,k,MORRInd::agn), dlamg, rhofac, mp);
        }

        m(i,j,k,MORRInd::fr)  = vr.mass;
        m(i,j,k,MORRInd::fi)  = vi.mass;
        m(i,j,k,MORRInd::fni) = vi.number;
        m(i,j,k,MORRInd::fs)  = vs.mass;
        m(i,j,k,MORRInd::fns) = vs.number;
        m(i,j,k,MORRInd::fnr) = vr.number;
        m(i,j,k,MORRInd::fc)  = vc.mass;
        m(i,j,k,MORRInd::fnc) = vc.number;
        m(i,j,k,MORRInd::fg)  = vg.mass;
        m(i,j,k,MORRInd::fng) = vg.number;

        // V3.3 modify fallspeed below level of precip
        Real rgvm = Real(0);
        for (int n = 0; n < nmom; ++n) {
            if (m(i,j,k,f_c[n]) < Real(1.e-10)) {
                m(i,j,k,f_c[n]) = m(i,j,k+1,f_c[n]);
            }
            rgvm = std::max(rgvm, m(i,j,k,f_c[n]));
        }

        // Number of split time steps
        nstep = std::max(static_cast<int>(rgvm * dt / m(i,j,k,MORRInd::dzq) + one), nstep);

        // Multiply variables by RHO
        for (int n = 0; n < nmom; ++n) {
            m(i,j,k,dum_c[n]) = m(i,j,k,dum_c[n]) * rho;
        }
    }

    const Real dt_sub = dt / nstep;

    for (int n = 1; n <= nstep; n++) {
        // Fallout of each moment at all levels; falout is the flux out of the bottom of
        // cell k (face k) for this substep
        for (int k = klo; k <= khi; k++) {
            for (int c = 0; c < nmom; ++c) {
                m(i,j,k,falout_c[c]) = m(i,j,k,f_c[c]) * m(i,j,k,dum_c[c]);
                sed(i,j,k,sed_c[c]) += m(i,j,k,falout_c[c]) * dt_sub;
            }
        }

        // Top of the column, and the levels below it, which also get what falls from above
        for (int k = khi; k >= klo; k--) {
            const Real dzq = m(i,j,k,MORRInd::dzq);
            const Real rho = m(i,j,k,MORRInd::rho);
            for (int c = 0; c < nmom; ++c) {
                if (k == khi) {
                    const Real faltnd = m(i,j,k,falout_c[c]) / dzq;
                    m(i,j,k,ten_c[c]) = m(i,j,k,ten_c[c]) - faltnd / nstep / rho;
                    m(i,j,k,dum_c[c]) = m(i,j,k,dum_c[c]) - faltnd * dt / nstep;
                } else {
                    const Real faltnd = (m(i,j,k+1,falout_c[c]) - m(i,j,k,falout_c[c])) / dzq;
                    m(i,j,k,ten_c[c]) = m(i,j,k,ten_c[c]) + faltnd / nstep / rho;
                    m(i,j,k,dum_c[c]) = m(i,j,k,dum_c[c]) + faltnd * dt / nstep;
                }
            }
        }

        // Precipitation and snowfall accumulation during the time step. The factor of
        // 1000 that converts from m to mm cancels with the density of liquid water.
        const Real faloutr = m(i,j,klo,MORRInd::faloutr);
        const Real faloutc = m(i,j,klo,MORRInd::faloutc);
        const Real falouts = m(i,j,klo,MORRInd::falouts);
        const Real falouti = m(i,j,klo,MORRInd::falouti);
        const Real faloutg = m(i,j,klo,MORRInd::faloutg);
        m(i,j,klo,MORRInd::precrt)      += (faloutr + faloutc + falouts + falouti + faloutg) * dt / nstep;
        m(i,j,klo,MORRInd::snowrt)      += (falouts + falouti + faloutg) * dt / nstep;
        m(i,j,klo,MORRInd::snowprt)     += (falouti + falouts) * dt / nstep;
        m(i,j,klo,MORRInd::grplprt)     += faloutg * dt / nstep;
        m(i,j,klo,MORRInd::rainonlyprt) += faloutr * dt / nstep;
        m(i,j,klo,MORRInd::snowonlyprt) += falouts * dt / nstep;
    }
}

/**
 * Apply the tendencies of one level, then the instantaneous processes (melting of cloud
 * ice, homogeneous freezing), the final size distributions and the effective radii.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void morrison_cell_finalize (Array4<Real> const& m, int i, int j, int k,
                             const Real rho_erf, const MorrisonParams& mp, const Real dt)
{
    Real& qc3d  = m(i,j,k,MORRInd::qc3d);
    Real& qi3d  = m(i,j,k,MORRInd::qi3d);
    Real& qni3d = m(i,j,k,MORRInd::qni3d);
    Real& qr3d  = m(i,j,k,MORRInd::qr3d);
    Real& qg3d  = m(i,j,k,MORRInd::qg3d);
    Real& nc3d  = m(i,j,k,MORRInd::nc3d);
    Real& ni3d  = m(i,j,k,MORRInd::ni3d);
    Real& ns3d  = m(i,j,k,MORRInd::ns3d);
    Real& nr3d  = m(i,j,k,MORRInd::nr3d);
    Real& ng3d  = m(i,j,k,MORRInd::ng3d);
    Real& t3d   = m(i,j,k,MORRInd::t3d);
    Real& qv3d  = m(i,j,k,MORRInd::qv3d);
    Real& qc3dten  = m(i,j,k,MORRInd::qc3dten);
    Real& qi3dten  = m(i,j,k,MORRInd::qi3dten);
    Real& qni3dten = m(i,j,k,MORRInd::qni3dten);
    Real& qr3dten  = m(i,j,k,MORRInd::qr3dten);
    Real& qg3dten  = m(i,j,k,MORRInd::qg3dten);
    Real& ni3dten  = m(i,j,k,MORRInd::ni3dten);
    Real& ns3dten  = m(i,j,k,MORRInd::ns3dten);
    const Real pres = m(i,j,k,MORRInd::pres);
    const Real rho  = m(i,j,k,MORRInd::rho);
    const Real xxlv = m(i,j,k,MORRInd::xxlv);
    const Real xxls = m(i,j,k,MORRInd::xxls);
    const Real xlf  = m(i,j,k,MORRInd::xlf);
    const Real cpm  = m(i,j,k,MORRInd::cpm);
    Real& lami = m(i,j,k,MORRInd::lami);
    Real& lamr = m(i,j,k,MORRInd::lamr);
    Real& lams = m(i,j,k,MORRInd::lams);
    Real& lamg = m(i,j,k,MORRInd::lamg);
    Real& lamc = m(i,j,k,MORRInd::lamc);
    Real& pgam = m(i,j,k,MORRInd::pgam);

    // Add on sedimentation tendencies for mixing ratio to rest of tendencies
    qr3dten  = qr3dten  + m(i,j,k,MORRInd::qrsten);
    qi3dten  = qi3dten  + m(i,j,k,MORRInd::qisten);
    qc3dten  = qc3dten  + m(i,j,k,MORRInd::qcsten);
    qg3dten  = qg3dten  + m(i,j,k,MORRInd::qgsten);
    qni3dten = qni3dten + m(i,j,k,MORRInd::qnisten);

    // Put all cloud ice in snow category if mean diameter exceeds 2 * DCS
    if (qi3d >= mp.qsmall && t3d < Real(273.15) && lami >= Real(1.e-10)) {
        if (one/lami >= Real(2)*mp.dcs) {
            qni3dten = qni3dten + qi3d/dt + qi3dten;
            ns3dten  = ns3dten + ni3d/dt + ni3dten;
            qi3dten  = -qi3d/dt;
            ni3dten  = -ni3d/dt;
        }
    }

    // Add tendencies
    qc3d  = qc3d  + qc3dten*dt;
    qi3d  = qi3d  + qi3dten*dt;
    qni3d = qni3d + qni3dten*dt;
    qr3d  = qr3d  + qr3dten*dt;
    nc3d  = nc3d  + m(i,j,k,MORRInd::nc3dten)*dt;
    ni3d  = ni3d  + ni3dten*dt;
    ns3d  = ns3d  + ns3dten*dt;
    nr3d  = nr3d  + m(i,j,k,MORRInd::nr3dten)*dt;
    if (mp.igraup == 0) {
        qg3d = qg3d + qg3dten*dt;
        ng3d = ng3d + m(i,j,k,MORRInd::ng3dten)*dt;
    }

    // Add temperature and water vapor tendencies from microphysics
    t3d  = t3d  + m(i,j,k,MORRInd::t3dten)*dt;
    qv3d = qv3d + m(i,j,k,MORRInd::qv3dten)*dt;

    const MorrisonSaturation sat = morrison_saturation(t3d, pres, qv3d, mp.ep_2);
    morrison_evaporate_small_at_subsaturation(sat, xxlv, xxls, cpm, qv3d, t3d,
                                              qc3d, qr3d, qi3d, qni3d, qg3d);

    morrison_zero_below_qsmall(qc3d,  nc3d, m(i,j,k,MORRInd::effc), mp.qsmall);
    morrison_zero_below_qsmall(qr3d,  nr3d, m(i,j,k,MORRInd::effr), mp.qsmall);
    morrison_zero_below_qsmall(qi3d,  ni3d, m(i,j,k,MORRInd::effi), mp.qsmall);
    morrison_zero_below_qsmall(qni3d, ns3d, m(i,j,k,MORRInd::effs), mp.qsmall);
    morrison_zero_below_qsmall(qg3d,  ng3d, m(i,j,k,MORRInd::effg), mp.qsmall);

    if (!(qc3d < mp.qsmall && qi3d < mp.qsmall && qni3d < mp.qsmall && qr3d < mp.qsmall && qg3d < mp.qsmall))
    {
        // Instantaneous processes

        // Melting of cloud ice to form rain
        if (qi3d >= mp.qsmall && t3d >= Real(273.15)) {
            morrison_convert_all(qi3d, ni3d, qr3d, nr3d, t3d, -xlf, cpm);
        }

        // ****SENSITIVITY - NO ICE
        if (mp.iliq != 1) {
            // Homogeneous freezing of cloud water
            if (t3d <= Real(233.15) && qc3d >= mp.qsmall) {
                morrison_convert_all(qc3d, nc3d, qi3d, ni3d, t3d, xlf, cpm);
            }
            // Homogeneous freezing of rain, to graupel or (no graupel) to snow
            if (t3d <= Real(233.15) && qr3d >= mp.qsmall) {
                if (mp.igraup == 0) {
                    morrison_convert_all(qr3d, nr3d, qg3d, ng3d, t3d, xlf, cpm);
                } else if (mp.igraup == 1) {
                    morrison_convert_all(qr3d, nr3d, qni3d, ns3d, t3d, xlf, cpm);
                }
            }
        }

        // Make sure number concentrations aren't negative
        ni3d = std::max(Real(0), ni3d);
        ns3d = std::max(Real(0), ns3d);
        nc3d = std::max(Real(0), nc3d);
        nr3d = std::max(Real(0), nr3d);
        ng3d = std::max(Real(0), ng3d);

        // Size distribution parameters
        if (qi3d >= mp.qsmall) {
            morrison_exponential_psd(qi3d, ni3d, mp.cons12, one/mp.di,
                                     mp.lammini, mp.lammaxi, lami, m(i,j,k,MORRInd::n0i));
        }
        if (qr3d >= mp.qsmall) {
            morrison_exponential_psd(qr3d, nr3d, mp.pi * mp.rhow, one/three,
                                     mp.lamminr, mp.lammaxr, lamr, m(i,j,k,MORRInd::n0r));
        }
        if (qc3d >= mp.qsmall) {
            morrison_cloud_droplet_psd(qc3d, nc3d, pres, t3d, mp.cons26,
                                       pgam, lamc, m(i,j,k,MORRInd::cdist1));
        }
        if (qni3d >= mp.qsmall) {
            morrison_exponential_psd(qni3d, ns3d, mp.cons1, one/mp.ds,
                                     mp.lammins, mp.lammaxs, lams, m(i,j,k,MORRInd::n0s));
        }
        if (qg3d >= mp.qsmall) {
            morrison_exponential_psd(qg3d, ng3d, mp.cons2, one/mp.dg,
                                     mp.lamming, mp.lammaxg, lamg, m(i,j,k,MORRInd::n0g));
        }
    }

    // Effective radii (25 microns when the species is absent)
    m(i,j,k,MORRInd::effi) = (qi3d  >= mp.qsmall) ? morrison_effective_radius(lami) : Real(25.0);
    m(i,j,k,MORRInd::effs) = (qni3d >= mp.qsmall) ? morrison_effective_radius(lams) : Real(25.0);
    m(i,j,k,MORRInd::effr) = (qr3d  >= mp.qsmall) ? morrison_effective_radius(lamr) : Real(25.0);
    m(i,j,k,MORRInd::effc) = (qc3d  >= mp.qsmall) ? morrison_cloud_droplet_effective_radius(pgam, lamc) : Real(25.0);
    m(i,j,k,MORRInd::effg) = (qg3d  >= mp.qsmall) ? morrison_effective_radius(lamg) : Real(25.0);

    // HM, 12/28/12: upper bound on ice number, needed to prevent very large ice number due
    // to homogeneous freezing of droplets, especially when INUM = 1; lowered to address the
    // problem of excessive and persistent anvil. This may reduce sensitivity to aerosol/CCN.
    ni3d = std::min(ni3d, Real(0.3e6) / rho);

    // Droplet number cannot exceed the aerosol concentration
    if (mp.inum == 0 && mp.iact == 2) {
        nc3d = std::min(nc3d, (mp.nanew1 + mp.nanew2) / rho);
    }

    // Switch for constant droplet number: change NDCNST from cm^-3 to kg^-1
    if (mp.inum == 1) {
        nc3d = morrison_ndcnst_to_number_mixing_ratio(mp.ndcnst, rho_erf);
    }
}

} // namespace

    // wrapper to do all the updating
    void
    Morrison::Advance (const Real& dt_advance,
                       const SolverChoice& sc)
    {
        m_sed.reset(dt_advance);

        if (!m_columns.needs_gather()) {
            AdvanceColumns(dt_advance, sc);
            return;
        }

        // Some column is split between boxes stacked in z. Sedimentation treats each
        // box's k range as the whole column, so run on a copy of the fields in which
        // every column lies in one box, then copy the result back.
        amrex::Array<FabPtr, MicVar_Morr::NumVars> saved;
        m_columns.gather(mic_fab_vars, saved);
        MultiFab* z_phys_nd_saved = m_z_phys_nd;
        auto z_phys_nd_col = m_columns.gather_nodal(m_z_phys_nd);
        m_z_phys_nd = z_phys_nd_col.get();
        auto sed_saved = m_sed.down;
        m_sed.down = m_columns.gather_one(*sed_saved);
        const MultiFab* inflow_saved = m_sed_inflow;
        auto inflow_col = (m_sed_inflow) ? m_columns.gather_one(*m_sed_inflow) : nullptr;
        m_sed_inflow = inflow_col.get();

        AdvanceColumns(dt_advance, sc);

        m_sed_inflow = inflow_saved;
        sed_saved->ParallelCopy(*m_sed.down);
        m_sed.down = sed_saved;
        m_z_phys_nd = z_phys_nd_saved;
        m_columns.scatter(mic_fab_vars, saved);
    }

    void
    Morrison::AdvanceColumns (const Real& dt_advance,
                              const SolverChoice& sc)
    {
        // Expose for GPU
        bool do_cond = m_do_cond;

        // Store timestep
        Real dt = dt_advance;

        // Check if CPP or FORT answer is used
        ParmParse pp("erf");
        bool use_morr_cpp_answer = true;
        pp.queryAdd("use_morr_cpp_answer", use_morr_cpp_answer);
#ifndef ERF_USE_MORR_FORT
        // Without the Fortran bridge, false would skip the microphysics update entirely
        if (!use_morr_cpp_answer) {
            Abort("erf.use_morr_cpp_answer = false requires the Morrison Fortran bridge "
                  "(USE_MORR_FORT=TRUE or ERF_ENABLE_MORR_FORT=ON)");
        }
#endif

        // Ensure that only one of these is true
        bool run_morr_cpp  =  use_morr_cpp_answer;
        bool run_morr_fort = !run_morr_cpp;

        // The Fortran bridge takes each box's k range as the whole column (the box top
        // as the model top, its bottom as the ground), so it needs boxes that do span it.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!run_morr_fort || m_columns.boxes_span_domain(),
            "Morrison: erf.use_morr_cpp_answer = false requires every box to span the domain in z "
            "(no grids split in z, and no refined levels short of the domain top or bottom)");

        // Only a column standing on the domain bottom has a ground to accumulate precipitation on
        const int dom_klo = m_geom.Domain().smallEnd(2);

        // Constant droplet concentration (if INUM = 1), cm^-3, which the inputs may override
        Real ndcnst = Real(250.0);
        pp.queryAdd("morrison_ndcnst", ndcnst);

        // Switches and constants of the scheme
        const MorrisonParams mp = make_morrison_params(sc.moisture_type == MoistureType::Morrison_NoIce,
                                                       ndcnst);

        // Loop through the grids
        //
        // The Fortran bridge is a whole-FAB interface: it declares every array over the fab
        // box and would stage a full-FAB copy of ~14 scratch arrays per tile. Iterate one
        // box at a time on that path.
        //
        const IntVect morr_tile_size = run_morr_fort ? IntVect::TheZeroVector() : TileNoZ();
        for (MFIter mfi(*mic_fab_vars[MicVar_Morr::qcl],morr_tile_size); mfi.isValid(); ++mfi)
        {
          auto box = mfi.tilebox();

          if (!box.ok()) { // Avoid going farther if the box is inverted (i.e., ilo > ihi or jlo > jhi).
              continue;
          }

          // Get array data from class member variables
          auto const& theta_arr = mic_fab_vars[MicVar_Morr::theta]->array(mfi);

          auto const& qv_arr  = mic_fab_vars[MicVar_Morr::qv]->array(mfi);

          auto const& qcl_arr = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
          auto const& qpr_arr = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
          auto const& qci_arr = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
          auto const& qps_arr = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
          auto const& qpg_arr = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);

          auto const& nc_arr = mic_fab_vars[MicVar_Morr::nc]->array(mfi);
          auto const& ni_arr = mic_fab_vars[MicVar_Morr::ni]->array(mfi);
          auto const& nr_arr = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
          auto const& ns_arr = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
          auto const& ng_arr = mic_fab_vars[MicVar_Morr::ng]->array(mfi);

          auto const& rho_arr         = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
          auto const& pres_arr        = mic_fab_vars[MicVar_Morr::pres]->array(mfi);
          auto const& rain_accum_arr  = mic_fab_vars[MicVar_Morr::rain_accum]->array(mfi);
          auto const& snow_accum_arr  = mic_fab_vars[MicVar_Morr::snow_accum]->array(mfi);
          auto const& graupel_accum_arr = mic_fab_vars[MicVar_Morr::graupel_accum]->array(mfi);
          auto const& precip_accum_arr  = mic_fab_vars[MicVar_Morr::precip_accum]->array(mfi);
          auto const& w_arr           = mic_fab_vars[MicVar_Morr::omega]->array(mfi);
          auto const& sed_arr         = m_sed.down->array(mfi);
          const auto inflow_arr = (m_sed_inflow) ? m_sed_inflow->const_array(mfi)
                                                 : Array4<const Real>{};

          // Extract box dimensions (the x and y extents are for the Fortran bridge)
          [[maybe_unused]] const int ilo = box.loVect()[0];
          [[maybe_unused]] const int ihi = box.hiVect()[0];
          [[maybe_unused]] const int jlo = box.loVect()[1];
          [[maybe_unused]] const int jhi = box.hiVect()[1];
          const int klo = box.loVect()[2];
          const int khi = box.hiVect()[2];

          Box grown_box(box); grown_box.grow(3);

          // Every array handed to the Fortran bridge -- both the mic_fab_vars and the staging
          // FABs below -- is declared there as DIMENSION(ims:ime,jms:jme,kms:kme), and those
          // bounds are taken from mfi.fabbox(). So the staging FABs must be allocated on the
          // fab box too; sizing them from a tilebox grown by a hard-coded 3 mis-declares the
          // strides whenever ngrow_state != 3 or the MFIter tiles the box.
          // Note fab_box always contains grown_box (ComputeGhostCells() >= 2, so the state
          // carries at least 3 ghost cells), so this only ever widens the staged region.
          const Box& fab_box = mfi.fabbox();

#if defined(ERF_USE_MORR_FORT) && defined(AMREX_USE_GPU)
          Arena* Arena_Used = The_Pinned_Arena();
#else
          Arena* Arena_Used = The_Async_Arena();
#endif

          // Calculate Exner function (PII) to convert potential temperature to temperature
          // PII = (P/P0)^(R/cp)
          FArrayBox pii_fab(fab_box, 1, Arena_Used);
          auto const& pii_arr = pii_fab.array();

          const Real p0 = Real(100000.0); // Reference pressure (Pa)

          const Real rdcp = m_rdOcp; // R/cp ratio

          // Calculate Exner function
          ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // NOTE: the Morrison Fortran version uses Pa not hPa so we didn't divide p by 100
            //       so we don't need to multiply by 100 here
            pii_arr(i,j,k) = std::pow((pres_arr(i,j,k)) / p0, rdcp);
          });

          // Create arrays for height differences (dz)
          FArrayBox dz_fab(fab_box, 1, Arena_Used);
          auto const& dz_arr = dz_fab.array();

          // Calculate height differences
          const Real dz_val = m_geom.CellSize(2);
          const Array4<const Real> z_arr = (m_z_phys_nd) ? m_z_phys_nd->const_array(mfi) : Array4<const Real> {};
          // z_phys_nd carries ComputeGhostCells()+2 ghost nodes, i.e. one more than the state,
          // so the (i+1,j+1,k+1) reads stay in bounds over the whole fab box.
          ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            dz_arr(i,j,k) = (z_arr) ? Real(0.25) * ( (z_arr(i  ,j  ,k+1) - z_arr(i  ,j  ,k))
                                                   + (z_arr(i+1,j  ,k+1) - z_arr(i+1,j  ,k))
                                                   + (z_arr(i  ,j+1,k+1) - z_arr(i  ,j+1,k))
                                                   + (z_arr(i+1,j+1,k+1) - z_arr(i+1,j+1,k)) ) : dz_val;
          });

          // The 2D Fortran arrays are declared DIMENSION(ims:ime,jms:jme), so this slab must
          // carry the fab box extents in x and y.
          Box fab_boxD(fab_box); fab_boxD.makeSlab(2,0);

          // Arrays to store precipitation rates
          FArrayBox    rainncv_fab(fab_boxD, 1, Arena_Used);
          FArrayBox         sr_fab(fab_boxD, 1, Arena_Used);     // Ratio of snow to total precipitation
          FArrayBox    snowncv_fab(fab_boxD, 1, Arena_Used);
          FArrayBox graupelncv_fab(fab_boxD, 1, Arena_Used);

          auto const& rainncv_arr = rainncv_fab.array();
          auto const& sr_arr      = sr_fab.array();
          auto const& snowncv_arr = snowncv_fab.array();
          auto const& graupelncv_arr = graupelncv_fab.array();

          // Initialize precipitation rate arrays to Real(0)
          ParallelFor(fab_boxD, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainncv_arr(i,j,k)    = Real(0);
            sr_arr(i,j,k)         = Real(0);
            snowncv_arr(i,j,k)    = Real(0);
            graupelncv_arr(i,j,k) = Real(0);
          });

          Box boxD(box); boxD.makeSlab(2,0);

          if (run_morr_cpp) {

            // One FAB to rule them all
            FArrayBox morr_fab(grown_box, MORRInd::NumInds, Arena_Used);
            morr_fab.template setVal<RunOn::Device>(0);
            auto const& morr_arr = morr_fab.array();

            // Copy the state into the column variables
            ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
              morr_arr(i,j,k,MORRInd::qc3d)  = qcl_arr(i,j,k);   // CLOUD WATER MIXING RATIO
              morr_arr(i,j,k,MORRInd::qi3d)  = qci_arr(i,j,k);   // CLOUD ICE   MIXING RATIO
              morr_arr(i,j,k,MORRInd::qni3d) = qps_arr(i,j,k);   // SNOW        MIXING RATIO
              morr_arr(i,j,k,MORRInd::qr3d)  = qpr_arr(i,j,k);   // RAIN        MIXING RATIO
              morr_arr(i,j,k,MORRInd::qg3d)  = qpg_arr(i,j,k);   // GRAUPEL     MIXING RATIO

              morr_arr(i,j,k,MORRInd::nc3d)  = nc_arr(i,j,k);    // CLOUD WATER NUMBER CONCENTRATION
              morr_arr(i,j,k,MORRInd::ni3d)  = ni_arr(i,j,k);    // CLOUD ICE   NUMBER CONCENTRATION
              morr_arr(i,j,k,MORRInd::ns3d)  = ns_arr(i,j,k);    // SNOW        NUMBER CONCENTRATION
              morr_arr(i,j,k,MORRInd::nr3d)  = nr_arr(i,j,k);    // RAIN        NUMBER CONCENTRATION
              morr_arr(i,j,k,MORRInd::ng3d)  = ng_arr(i,j,k);    // GRAUPEL     NUMBER CONCENTRATION

              morr_arr(i,j,k,MORRInd::t3d)  = theta_arr(i,j,k) * pii_arr(i,j,k);  // TEMPERATURE
              morr_arr(i,j,k,MORRInd::qv3d) = qv_arr(i,j,k);                     // WATER VAPOR MIXING RATIO
              morr_arr(i,j,k,MORRInd::pres) = pres_arr(i,j,k);                   // ATMOSPHERIC PRESSURE
              morr_arr(i,j,k,MORRInd::dzq)  = dz_arr(i,j,k);                     // DIFFERENCE IN HEIGHT ACROSS LEVEL
              morr_arr(i,j,k,MORRInd::w3d)  = w_arr(i,j,k);                      // GRID-SCALE VERTICAL VELOCITY

              // There are no cumulus tendencies passed to Morrison (the Fortran zeros these too)
              morr_arr(i,j,k,MORRInd::qrcu1d) = Real(0); // RAIN FROM CUMULUS PARAMETERIZATION
              morr_arr(i,j,k,MORRInd::qscu1d) = Real(0); // SNOW FROM CUMULUS PARAMETERIZATION
              morr_arr(i,j,k,MORRInd::qicu1d) = Real(0); // ICE FROM CUMULUS PARAMETERIZATION
            });

            // MORR_TWO_MOMENT_MICRO, one column at a time
            ParallelFor(boxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
              // LTRUE: are there hydrometeors in the column
              bool ltrue = false;

              for (int k = klo; k <= khi; k++) {
                  if (morrison_cell_microphysics(morr_arr, i, j, k, rho_arr(i,j,k), mp, dt, do_cond)) {
                      ltrue = true;
                  }
              }

              // At the top of a fine patch, rain, snow, ... come in from the coarse level.
              // Add what enters during this step to the top cell as a tendency, before the
              // fall speeds are computed, so it falls with what is already there; record it
              // on the top face.
              if (inflow_arr && inflow_arr(i,j,khi+1,MorrSed::NumComps) > Real(0.5)) {
                  const Real fac = one / (morr_arr(i,j,khi,MORRInd::rho) * morr_arr(i,j,khi,MORRInd::dzq) * dt);
                  morr_arr(i,j,khi,MORRInd::qc3dten)  += inflow_arr(i,j,khi+1,MorrSed::qc) * fac;
                  morr_arr(i,j,khi,MORRInd::qi3dten)  += inflow_arr(i,j,khi+1,MorrSed::qi) * fac;
                  morr_arr(i,j,khi,MORRInd::qr3dten)  += inflow_arr(i,j,khi+1,MorrSed::qr) * fac;
                  morr_arr(i,j,khi,MORRInd::qni3dten) += inflow_arr(i,j,khi+1,MorrSed::qs) * fac;
                  morr_arr(i,j,khi,MORRInd::qg3dten)  += inflow_arr(i,j,khi+1,MorrSed::qg) * fac;
                  morr_arr(i,j,khi,MORRInd::nc3dten)  += inflow_arr(i,j,khi+1,MorrSed::nc) * fac;
                  morr_arr(i,j,khi,MORRInd::ni3dten)  += inflow_arr(i,j,khi+1,MorrSed::ni) * fac;
                  morr_arr(i,j,khi,MORRInd::nr3dten)  += inflow_arr(i,j,khi+1,MorrSed::nr) * fac;
                  morr_arr(i,j,khi,MORRInd::ns3dten)  += inflow_arr(i,j,khi+1,MorrSed::ns) * fac;
                  morr_arr(i,j,khi,MORRInd::ng3dten)  += inflow_arr(i,j,khi+1,MorrSed::ng) * fac;
                  for (int n = 0; n < MorrSed::NumComps; ++n) {
                      sed_arr(i,j,khi+1,n) += inflow_arr(i,j,khi+1,n);
                  }
                  ltrue = true;
              }

              if (ltrue) {
                  morrison_column_sedimentation(morr_arr, sed_arr, i, j, klo, khi, mp, dt);

                  for (int k = klo; k <= khi; k++) {
                      morrison_cell_finalize(morr_arr, i, j, k, rho_arr(i,j,k), mp, dt);
                  }
              }

              // Copy the column variables back to the state
              for (int k = klo; k <= khi; k++) {
                  qcl_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qc3d);
                  qci_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qi3d);
                  qps_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qni3d);
                  qpr_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qr3d);
                  nc_arr(i,j,k)  = morr_arr(i,j,k,MORRInd::nc3d);
                  ni_arr(i,j,k)  = morr_arr(i,j,k,MORRInd::ni3d);
                  ns_arr(i,j,k)  = morr_arr(i,j,k,MORRInd::ns3d);
                  nr_arr(i,j,k)  = morr_arr(i,j,k,MORRInd::nr3d);
                  qpg_arr(i,j,k) = morr_arr(i,j,k,MORRInd::qg3d);
                  ng_arr(i,j,k)  = morr_arr(i,j,k,MORRInd::ng3d);

                  theta_arr(i,j,k) = morr_arr(i,j,k,MORRInd::t3d) / pii_arr(i,j,k); // Convert temp back to potential temp
                  qv_arr(i,j,k)    = morr_arr(i,j,k,MORRInd::qv3d);
              }

              rainncv_arr(i,j,0)    = morr_arr(i,j,klo,MORRInd::precrt);
              snowncv_arr(i,j,0)    = morr_arr(i,j,klo,MORRInd::snowprt);
              graupelncv_arr(i,j,0) = morr_arr(i,j,klo,MORRInd::grplprt);
              sr_arr(i,j,0)         = morr_arr(i,j,klo,MORRInd::snowrt) / (morr_arr(i,j,klo,MORRInd::precrt) + Real(1.e-12));

              // Each species alone, and the total, which also holds the cloud water
              // and cloud ice that fell out
              if (klo == dom_klo) {
                  rain_accum_arr(i,j,klo)    += morr_arr(i,j,klo,MORRInd::rainonlyprt);
                  snow_accum_arr(i,j,klo)    += morr_arr(i,j,klo,MORRInd::snowonlyprt);
                  graupel_accum_arr(i,j,klo) += morr_arr(i,j,klo,MORRInd::grplprt);
                  precip_accum_arr(i,j,klo)  += morr_arr(i,j,klo,MORRInd::precrt);
              }
            });
          }

          if (run_morr_fort) {
#ifdef ERF_USE_MORR_FORT
#include  "ERF_Morrison_Advance_F.H"
#endif
          } // run_morr_fort

        }
    }
