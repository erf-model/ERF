#include "ERF_WSM6.H"

#include <AMReX_Gpu.H>
#include "ERF_EOS.H"

using namespace amrex;

void
WSM6::Init(const MultiFab& cons_in,
           const BoxArray&,
           const Geometry& geom,
           const Real& dt_advance,
           std::unique_ptr<MultiFab>& z_phys_nd,
           std::unique_ptr<MultiFab>& detJ_cc)
{
    dt = dt_advance;
    m_geom = geom;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_WSM6::rain_accum, MicVar_WSM6::snow_accum, MicVar_WSM6::graup_accum};

    for (int ivar = 0; ivar < MicVar_WSM6::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.0);
    }

    nlev = m_geom.Domain().length(2);
    zlo = m_geom.Domain().smallEnd(2);
    zhi = m_geom.Domain().bigEnd(2);

    initialize_coeffs();
}

void
WSM6::Copy_State_to_Micro(const MultiFab& cons_in)
{
    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        // Match Morrison behavior: refresh microphysics ghost zones from state.
        // WSM6 Fortran reads the full (ims:ime, jms:jme, kms:kme) slab.
        const auto& box3d = mfi.growntilebox();
        auto states = cons_in.array(mfi);

        auto rho = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto theta = mic_fab_vars[MicVar_WSM6::theta]->array(mfi);
        auto tabs = mic_fab_vars[MicVar_WSM6::tabs]->array(mfi);
        auto pres = mic_fab_vars[MicVar_WSM6::pres]->array(mfi);

        auto qv = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto qc = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto qi = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto qr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto qs = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto qg = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            rho(i,j,k) = states(i,j,k,Rho_comp);
            theta(i,j,k) = states(i,j,k,RhoTheta_comp) / states(i,j,k,Rho_comp);

            qv(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ1_comp) / states(i,j,k,Rho_comp));
            qc(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ2_comp) / states(i,j,k,Rho_comp));
            qi(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ3_comp) / states(i,j,k,Rho_comp));
            qr(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ4_comp) / states(i,j,k,Rho_comp));
            qs(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ5_comp) / states(i,j,k,Rho_comp));
            qg(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ6_comp) / states(i,j,k,Rho_comp));

            tabs(i,j,k) = getTgivenRandRTh(states(i,j,k,Rho_comp),
                                           states(i,j,k,RhoTheta_comp),
                                           qv(i,j,k));
            pres(i,j,k) = getPgivenRTh(states(i,j,k,RhoTheta_comp), qv(i,j,k));
        });
    }
}

void
WSM6::initialize_coeffs()
{
    using amrex::Real;

    // Exact port of Fortran rgmma() — ERF_module_mp_wsm6.F90 lines 1472-1495
    // Weierstrass infinite product form for Gamma(x).
    // Special case x==1 returns 0.0 matching Fortran exactly.
    // Never triggered in practice (bvtr1=1.8, bvts1=1.41, etc.)
    // but preserved for bit-compatible validation against Fortran reference.
    // CPU-only: called from initialize_coeffs(), not from GPU kernels.
    auto rgmma = [](Real x) -> Real {
        if (x == Real(1.0)) return Real(0.0);
        constexpr Real euler = Real(0.577215664901532);
        Real result = x * std::exp(euler * x);
        for (int i = 1; i <= 10000; ++i) {
            Real y = Real(i);
            result = result * (Real(1.0) + x/y) * std::exp(-x/y);
        }
        return Real(1.0) / result;
    };

    // Physical constants matching mp_wsm6_init argument list
    // den0: reference air density at 850mb (kg/m^3)
    // denr: liquid water density — rhoh2o from ERF_Constants.H
    // dens_arg: snow density — dens_snow constexpr = 100.0
    // cl: specific heat liquid water — Cp_l from ERF_Constants.H
    // cpv: specific heat water vapor — Cp_v from ERF_Constants.H
    const Real den0     = Real(1.28);
    const Real denr     = Real(rhoh2o);
    const Real dens_arg = dens_snow;
    const Real cl       = Real(Cp_l);
    const Real cpv_loc  = Real(Cp_v);

    // hail_opt branch: 5 regime-dependent coefficients
    if (m_hail_opt) {
        m_n0g       = Real(4.0e4);
        m_deng      = Real(700.0);
        m_avtg      = Real(285.0);
        m_bvtg      = Real(0.8);
        m_lamdagmax = Real(2.0e4);
    } else {
        m_n0g       = Real(4.0e6);
        m_deng      = Real(500.0);
        m_avtg      = Real(330.0);
        m_bvtg      = Real(0.8);
        m_lamdagmax = Real(6.0e4);
    }

    m_pi_wsm6 = Real(4.0) * std::atan(Real(1.0));
    m_xlv1    = cl - cpv_loc;

    m_qc0  = Real(4.0)/Real(3.0) * m_pi_wsm6 * denr
              * std::pow(r0, Real(3.0)) * xncr / den0;
    m_qck1 = Real(0.104) * Real(9.8) * peaut
              / std::pow(xncr * denr, Real(1.0)/Real(3.0))
              / xmyu * std::pow(den0, Real(4.0)/Real(3.0));
    m_pidnc = m_pi_wsm6 * denr / Real(6.0);

    // Rain coefficients
    m_bvtr1   = Real(1.0) + bvtr;
    m_bvtr2   = Real(2.5) + Real(0.5) * bvtr;
    m_bvtr3   = Real(3.0) + bvtr;
    m_bvtr4   = Real(4.0) + bvtr;
    m_bvtr6   = Real(6.0) + bvtr;
    m_g1pbr   = rgmma(m_bvtr1);
    m_g3pbr   = rgmma(m_bvtr3);
    m_g4pbr   = rgmma(m_bvtr4);
    m_g6pbr   = rgmma(m_bvtr6);
    m_g5pbro2 = rgmma(m_bvtr2);
    m_pvtr    = avtr * m_g4pbr / Real(6.0);
    m_eacrr   = Real(1.0);
    m_pacrr   = m_pi_wsm6 * n0r * avtr * m_g3pbr * Real(0.25) * m_eacrr;
    m_precr1  = Real(2.0) * m_pi_wsm6 * n0r * Real(0.78);
    m_precr2  = Real(2.0) * m_pi_wsm6 * n0r * Real(0.31)
                * std::pow(avtr, Real(0.5)) * m_g5pbro2;
    m_roqimax = Real(2.08e22) * std::pow(dimax, Real(8.0));

    // Snow coefficients
    m_bvts1   = Real(1.0) + bvts;
    m_bvts2   = Real(2.5) + Real(0.5) * bvts;
    m_bvts3   = Real(3.0) + bvts;
    m_bvts4   = Real(4.0) + bvts;
    m_g1pbs   = rgmma(m_bvts1);
    m_g3pbs   = rgmma(m_bvts3);
    m_g4pbs   = rgmma(m_bvts4);
    m_g5pbso2 = rgmma(m_bvts2);
    m_pvts    = avts * m_g4pbs / Real(6.0);
    m_pacrs   = m_pi_wsm6 * n0s * avts * m_g3pbs * Real(0.25);
    m_precs1  = Real(4.0) * n0s * Real(0.65);
    m_precs2  = Real(4.0) * n0s * Real(0.44)
                * std::pow(avts, Real(0.5)) * m_g5pbso2;
    m_pidn0r  = m_pi_wsm6 * denr * n0r;
    m_pidn0s  = m_pi_wsm6 * dens_arg * n0s;
    m_pacrc   = m_pi_wsm6 * n0s * avts * m_g3pbs * Real(0.25) * eacrc;

    // Graupel/hail coefficients
    m_bvtg1   = Real(1.0) + m_bvtg;
    m_bvtg2   = Real(2.5) + Real(0.5) * m_bvtg;
    m_bvtg3   = Real(3.0) + m_bvtg;
    m_bvtg4   = Real(4.0) + m_bvtg;
    m_g1pbg   = rgmma(m_bvtg1);
    m_g3pbg   = rgmma(m_bvtg3);
    m_g4pbg   = rgmma(m_bvtg4);
    m_pacrg   = m_pi_wsm6 * m_n0g * m_avtg * m_g3pbg * Real(0.25);
    m_g5pbgo2 = rgmma(m_bvtg2);
    m_pvtg    = m_avtg * m_g4pbg / Real(6.0);
    m_precg1  = Real(2.0) * m_pi_wsm6 * m_n0g * Real(0.78);
    m_precg2  = Real(2.0) * m_pi_wsm6 * m_n0g * Real(0.31)
                * std::pow(m_avtg, Real(0.5)) * m_g5pbgo2;
    m_pidn0g  = m_pi_wsm6 * m_deng * m_n0g;

    // Slope parameter limits
    m_rslopermax  = Real(1.0) / lamdarmax;
    m_rslopesmax  = Real(1.0) / lamdasmax;
    m_rslopegmax  = Real(1.0) / m_lamdagmax;
    m_rsloperbmax = std::pow(m_rslopermax, bvtr);
    m_rslopesbmax = std::pow(m_rslopesmax, bvts);
    m_rslopegbmax = std::pow(m_rslopegmax, m_bvtg);
    m_rsloper2max = m_rslopermax * m_rslopermax;
    m_rslopes2max = m_rslopesmax * m_rslopesmax;
    m_rslopeg2max = m_rslopegmax * m_rslopegmax;
    m_rsloper3max = m_rsloper2max * m_rslopermax;
    m_rslopes3max = m_rslopes2max * m_rslopesmax;
    m_rslopeg3max = m_rslopeg2max * m_rslopegmax;
}
