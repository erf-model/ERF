#include "ERF_WDM6.H"

#include <AMReX_Gpu.H>
#include "ERF_EOS.H"

using namespace amrex;

void
WDM6::Init(const MultiFab& cons_in,
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

    // Read CCN concentration from input
    amrex::ParmParse pp("wdm6");
    pp.query("ccn0", m_ccn0);  // default 100.0e6 m^-3

    // Graupel/hail regime selector. The Fortran gates on `hail_opt .eq. 1`
    // (ERF_module_mp_wdm6.F90:3259), so every value other than 1 selects the
    // graupel branch; the bool below reproduces that test exactly rather than
    // inventing a stricter 0/1 validation the Fortran does not perform.
    // Must be resolved before initialize_coeffs(), which branches on it.
    int hail_opt_in = 0;
    pp.query("hail_opt", hail_opt_in);  // default 0 = graupel
    m_hail_opt = (hail_opt_in == 1);

    amrex::Print() << "WDM6 Initialization: CCN0 = " << m_ccn0 << " #/m^3"
                   << ", hail_opt = " << (m_hail_opt ? 1 : 0)
                   << (m_hail_opt ? " (hail)" : " (graupel)") << "\n";

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_WDM6::rain_accum, MicVar_WDM6::snow_accum, MicVar_WDM6::graup_accum};

    // Select appropriate Arena based on execution mode
    // - Fortran bridge + GPU: Use managed memory for CPU (Fortran) ↔ GPU data transfer
    // - Fortran bridge CPU-only: Use pinned memory for CPU efficiency
    // - C++ GPU kernels: Use standard device memory
#if defined(ERF_USE_WDM6_FORT) && defined(AMREX_USE_GPU)
    Arena* Arena_Used = The_Managed_Arena();
#elif defined(ERF_USE_WDM6_FORT)
    Arena* Arena_Used = The_Pinned_Arena();
#else
    Arena* Arena_Used = The_Arena();
#endif

    for (int ivar = 0; ivar < MicVar_WDM6::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect(),
                                                        MFInfo().SetArena(Arena_Used));
        mic_fab_vars[ivar]->setVal(0.0);
    }

    // Initialize nn to ccn0 / rho (matches WRF itimestep==1 behavior)
    // This must be done here so nn is available before the first microphysics call
    // NOTE: Can't write to cons_in (it's const), so only initialize mic_fab_vars.
    // Copy_Micro_to_State will write it back to state after the first microphysics call.
    // IMPORTANT: Use growntilebox to include ghost zones, since Copy_State_to_Micro will
    // skip reading nn and expect it to be initialized everywhere.
    // nn is NOT seeded here. The Fortran seeds it inside the per-step driver
    // wdm6 (:219), immediately before the kernel, not in wdm6init, so the
    // faithful ERF home is Copy_State_to_Micro, which is the per-step pack.
    // Seeding here instead was also wrong on restart: Init runs again on a
    // restart, so it discarded the checkpoint's evolved RhoQ8 and the first
    // Copy_State_to_Micro kept the uniform background. See the state-based
    // seeding in Copy_State_to_Micro below for the condition and its rationale.

    nlev = m_geom.Domain().length(2);
    zlo = m_geom.Domain().smallEnd(2);
    zhi = m_geom.Domain().bigEnd(2);

    initialize_coeffs();
}

void
WDM6::Copy_State_to_Micro(const MultiFab& cons_in)
{
    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        // Match Morrison behavior: refresh microphysics ghost zones from state.
        const auto& box3d = mfi.growntilebox();
        auto states = cons_in.array(mfi);

        auto rho = mic_fab_vars[MicVar_WDM6::rho]->array(mfi);
        auto theta = mic_fab_vars[MicVar_WDM6::theta]->array(mfi);
        auto tabs = mic_fab_vars[MicVar_WDM6::tabs]->array(mfi);
        auto pres = mic_fab_vars[MicVar_WDM6::pres]->array(mfi);

        auto qv = mic_fab_vars[MicVar_WDM6::qv]->array(mfi);
        auto qc = mic_fab_vars[MicVar_WDM6::qc]->array(mfi);
        auto qi = mic_fab_vars[MicVar_WDM6::qi]->array(mfi);
        auto qr = mic_fab_vars[MicVar_WDM6::qr]->array(mfi);
        auto qs = mic_fab_vars[MicVar_WDM6::qs]->array(mfi);
        auto qg = mic_fab_vars[MicVar_WDM6::qg]->array(mfi);

        auto nn = mic_fab_vars[MicVar_WDM6::nn]->array(mfi);
        auto nc = mic_fab_vars[MicVar_WDM6::nc]->array(mfi);
        auto nr = mic_fab_vars[MicVar_WDM6::nr]->array(mfi);

        const Real ccn0_local = m_ccn0;  // CCN concentration in #/m³

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            rho(i,j,k) = states(i,j,k,Rho_comp);
            theta(i,j,k) = states(i,j,k,RhoTheta_comp) / states(i,j,k,Rho_comp);

            qv(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ1_comp) / states(i,j,k,Rho_comp));
            qc(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ2_comp) / states(i,j,k,Rho_comp));
            qi(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ3_comp) / states(i,j,k,Rho_comp));
            qr(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ4_comp) / states(i,j,k,Rho_comp));
            qs(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ5_comp) / states(i,j,k,Rho_comp));
            qg(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ6_comp) / states(i,j,k,Rho_comp));

            // Number concentrations
            nc(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ7_comp) / states(i,j,k,Rho_comp));
            nr(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ9_comp) / states(i,j,k,Rho_comp));

            // nn: seed the CCN reservoir where state carries nothing yet,
            // otherwise read it like every other state variable.
            //
            // This mirrors the Fortran, which does its seeding in the per-step
            // driver wdm6 (:219) rather than in wdm6init, so this is the
            // faithful location. The Fortran gates on itimestep == 1; RhoQ8 == 0
            // is equivalent for every reachable case and needs no step plumbed
            // through the Init interface. Cold start and restart from a step-0
            // checkpoint both have RhoQ8 == 0 and seed, matching itimestep == 1;
            // every later restart has RhoQ8 > 0 and restores. A genuine zero
            // cannot occur after initialization because Advance bounds nn well
            // above zero (:602). Note the Fortran block is in fact unreachable
            // from ERF: the bridge enters wdm62D directly via mp_wdm6_run and
            // skips the wdm6 wrapper, so BOTH legs take their nn from here.
            const Real nn_state = states(i,j,k,RhoQ8_comp) / states(i,j,k,Rho_comp);
            nn(i,j,k) = (nn_state > Real(0.0)) ? nn_state : ccn0_local;

            // WDM6: DO NOT enforce nc/nr minimums here!
            // WRF starts with nc=0, nr=0 and lets CCN activation build nc naturally during
            // the first microphysics call. Enforcing nc=10, nr=0.01 here prevents proper
            // activation from nn. Minimums are enforced in Advance() right before physics,
            // not during state copying.

            tabs(i,j,k) = getTgivenRandRTh(states(i,j,k,Rho_comp),
                                           states(i,j,k,RhoTheta_comp),
                                           qv(i,j,k));
            pres(i,j,k) = getPgivenRTh(states(i,j,k,RhoTheta_comp), qv(i,j,k));
        });
    }

    // After first Copy_State_to_Micro, nn has been preserved from Init().
    // DON'T clear the flag yet - wait until after Copy_Micro_to_State writes nn to state!
}

void
WDM6::initialize_coeffs()
{
    using amrex::Real;

    // Exact port of Fortran rgmma() function from WSM6
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

    // Physical constants
    const Real den0     = Real(1.28);
    const Real denr     = Real(rhoh2o);
    const Real dens_arg = dens_snow;
    const Real cl       = Real(Cp_l);
    const Real cpv_loc  = Real(Cp_v);

    // hail_opt branch
    if (m_hail_opt) {
        m_n0g       = Real(4.0e4);
        m_deng      = Real(700.0);
        m_avtg      = Real(285.0);
        m_bvtg      = wdm6_literal(0.8);
        m_lamdagmax = Real(2.0e4);
    } else {
        m_n0g       = Real(4.0e6);
        m_deng      = Real(500.0);
        m_avtg      = Real(330.0);
        m_bvtg      = wdm6_literal(0.8);
        m_lamdagmax = Real(6.0e4);
    }

    // LITERAL PRECISION CONTRACT, second hazard class. The Fortran writes
    //     pi = 4.*atan(1.)
    // at ERF_module_mp_wdm6.F90:3111. Neither literal carries a kind suffix, so
    // `1.` is a default REAL and the generic ATAN resolves to its SINGLE
    // precision specific; the whole expression is evaluated in float32 and only
    // then widened to kind_phys. Measured with gfortran:
    //     4.*atan(1.)         = 3.14159274101257324219
    //     4._kp*atan(1._kp)   = 3.14159265358979311600
    //     relative difference = 2.7828e-08
    // This is distinct from the unsuffixed-parameter hazard the contract already
    // covers: nothing here is a rounded constant, it is an intrinsic whose result
    // kind is chosen by its argument kind. pi is load-bearing -- it feeds pidn0s,
    // pidnc, pidnr, pidn0r, pidn0s, pidn0g, pacrs, pacrc and the pi*pi in G10c's
    // pfrzdtc -- so the error propagates broadly.
    // Routed through wdm6_literal like every other constant: it is constexpr but
    // not restricted to constant expressions, so it applies the same rounding to
    // a runtime value. Do not re-expand the toggle inline here; the contract is
    // that ERF_WDM6_F32_LITERALS is honoured in exactly one place.
    m_pi_wdm6 = wdm6_literal(Real(4.0) * std::atan(Real(1.0)));
    m_xlv1    = cl - cpv_loc;

    // Cloud droplet parameters. The Fortran, ERF_module_mp_wdm6.F90:3114-3116:
    //     qc0  = 4./3.*pi*denr*r0**3.*xncr0/den0
    //     qc1  = 4./3.*pi*denr*r0**3.*xncr1/den0
    //     qck1 = .104*9.8*peaut/(denr)**(1./3.)/xmyu*den0**(4./3.)
    // qc0 takes xncr0 and qc1 takes xncr1, unconditionally: the maritime versus
    // continental choice is made downstream by the slmsk branch that sets
    // qcr = qc0 or qc1 (:628-632), not here. The previous code used xncr for qc0
    // and a ccn0-dependent xncr_use for qc1, so with the default ccn0 = 1.0e8 it
    // selected xncr0 for BOTH.
    //
    // qck1 is the consequential one. WDM6 divides by denr**(1./3.) alone; the
    // previous code divided by (xncr*denr)**(1./3.), which is the WSM6 form at
    // ERF_module_mp_wsm6.F90:105. WSM6 is single-moment and folds its fixed
    // droplet number xncr into qck1; WDM6 is double-moment and already carries
    // number explicitly as ncr(:,:,2)**(-1./3.) in praut, so the WSM6 form
    // double-counts it. The error is exactly xncr**(1./3.) = (3.0e8)**(1/3) =
    // 669.4, and praut and nraut were both measured low by that factor at
    // (199,3,43) from bitwise-equal inputs.
    m_qc0  = wdm6_literal(4.0/3.0) * m_pi_wdm6 * denr
              * std::pow(r0, Real(3.0)) * xncr0 / den0;
    m_qc1  = wdm6_literal(4.0/3.0) * m_pi_wdm6 * denr
              * std::pow(r0, Real(3.0)) * xncr1 / den0;
    // A FOURTH HAZARD CLASS under the literal precision contract: the float
    // assumptions apply to the whole literal PRODUCT, not merely to each
    // literal. In the Fortran, .104 and 9.8 are both default REAL, so
    // `.104*9.8` is itself a default-REAL expression and the PRODUCT is
    // evaluated in single precision before anything widens it:
    //     float32(float32(.104)*float32(9.8)) = 1.01920008659362792969
    // Applying the float assumption to each operand separately and then
    // multiplying performs the product in double:
    //     float32(.104)*float32(9.8)          = 1.01920004086494486728
    // a 4.486723e-08 relative difference. It lands on qck1 as a constant
    // multiplicative factor, so praut carries exactly that relative error at
    // every k, and nraut = 3.5e9*den*praut inherits it unchanged.
    //
    // Proven by reconstruction at clean SHA ac1e445fa from bitwise-equal
    // PRE_G13A inputs, no Tier 2 instrumentation: of 61 live G13a cells across
    // columns (108,3) and (199,3) over two steps, 59 reproduce BOTH legs'
    // printed praut bitwise to 20 digits from the two qck1 variants; the other
    // 2 are cells where the min(praut, qc/dtcld) clamp binds and the legs
    // already agree bitwise.
    //
    // The nesting is load-bearing. wdm6_literal(0.104 * 9.8) is NOT equivalent
    // and is worse at 1.169636e-07: it rounds the double product, whereas the
    // Fortran rounds the product of already-rounded operands. Round the
    // operands first, then round their product. Under Option B the helper is
    // the identity, so this collapses to 0.104*9.8 in double, which is what
    // -fdefault-real-8 gives on the Fortran side.
    m_qck1 = wdm6_literal(wdm6_literal(0.104) * wdm6_literal(9.8)) * peaut
              / std::pow(denr, wdm6_literal(1.0/3.0))
              / xmyu * std::pow(den0, wdm6_literal(4.0/3.0));
    m_pidnc = m_pi_wdm6 * denr / Real(6.0);

    // Rain coefficients
    m_bvtr1   = Real(1.0) + bvtr;
    m_bvtr2   = Real(2.0) + bvtr;
    m_bvtr3   = Real(3.0) + bvtr;
    m_bvtr4   = Real(4.0) + bvtr;
    m_bvtr5   = Real(5.0) + bvtr;
    m_bvtr6   = Real(6.0) + bvtr;
    m_bvtr7   = Real(7.0) + bvtr;
    m_bvtr2o5 = Real(2.5) + Real(0.5) * bvtr;
    m_bvtr3o5 = Real(3.5) + Real(0.5) * bvtr;

    m_g1pbr   = rgmma(m_bvtr1);
    m_g2pbr   = rgmma(m_bvtr2);
    m_g3pbr   = rgmma(m_bvtr3);
    m_g4pbr   = rgmma(m_bvtr4);
    m_g5pbr   = rgmma(m_bvtr5);
    m_g6pbr   = rgmma(m_bvtr6);
    m_g7pbr   = rgmma(m_bvtr7);
    m_g5pbro2 = rgmma(m_bvtr2o5);
    m_g7pbro2 = rgmma(m_bvtr3o5);

    m_pvtr    = avtr * m_g5pbr / Real(24.0);
    m_pvtrn   = avtr * m_g2pbr;
    m_eacrr   = Real(1.0);
    m_pacrr   = m_pi_wdm6 * n0r * avtr * m_g3pbr * Real(0.25) * m_eacrr;
    // Fortran WDM6 :3200-3201
    //     precr1 = 2.*pi*1.56
    //     precr2 = 2.*pi*.31*avtr**.5*g7pbro2
    // NEITHER carries n0r. WSM6 :121-122 does -- precr1 = 2.*pi*n0r*.78 and
    // precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2 -- for the same reason qck1
    // differs between the schemes: WSM6 is single-moment and folds the fixed
    // intercept n0r into its coefficients, while WDM6 is double-moment and
    // multiplies by ncr(:,:,3) explicitly at the use site in prevp. precr1 had
    // already been ported from the WDM6 form; precr2 kept WSM6's n0r, making it
    // 8.0e6 times too large. That drove the raw prevp far enough negative to
    // hit its satdt/2 clamp on the native leg only: at (103,3,34) on step 3 the
    // bridge held prevp at -8.3286e-09 while the native leg sat exactly on the
    // clamp at -2.5824185419052320E-06, 310x larger.
    m_precr1  = wdm6_literal(2.0) * m_pi_wdm6 * wdm6_literal(1.56);
    m_precr2  = wdm6_literal(2.0) * m_pi_wdm6 * wdm6_literal(0.31)
                * std::pow(avtr, wdm6_literal(0.5)) * m_g7pbro2;
    // Fortran: roqimax = 2.08e22*dimax**8
    // Two mismatches, both confirmed against the Tier 2 psaut trace:
    //   - 2.08e22 is unsuffixed there, so it is float32-rounded before being
    //     promoted. Route it through wdm6_literal so it follows the same
    //     switch as the ERF_WDM6.H constants.
    //   - dimax**8 has an INTEGER exponent, which gfortran expands to repeated
    //     multiplication, not a libm pow call. Written as explicit squaring
    //     here so the two sides agree structurally rather than by luck. For
    //     this value the two routes happen to agree bitwise, but that is not
    //     guaranteed in general and should not be relied on.
    {
        const Real d2 = dimax * dimax;
        const Real d4 = d2 * d2;
        m_roqimax = wdm6_literal(2.08e22) * (d4 * d4);
    }
    // Fortran: xmmax = (dimax/dicon)**2, again an integer exponent.
    {
        const Real r = dimax / dicon;
        m_xmmax = r * r;
    }

    m_pidn0r  = m_pi_wdm6 * denr * n0r;
    m_pidnr   = Real(4.0) * m_pi_wdm6 * denr;

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
    m_pacrs   = m_pi_wdm6 * n0s * avts * m_g3pbs * Real(0.25);
    // Fortran :3248-3249
    //     precs1 = 4.*n0s*.65
    //     precs2 = 4.*n0s*.44*avts**.5*g5pbso2
    // .65 and .44 are unsuffixed so they carry the float assumption; 4. and
    // avts**.5 are binary-exact and need no routing. Proven at Tier 2 on the
    // G13E psdep decomposition at (199,3), k=85..98: coeres, rh-1 and work1
    // are all bitwise, precs1 diverged exactly 3.667978e-08 and precs2 exactly
    // 5.418604e-09, matching float32(.65) and float32(.44), and psdep's
    // 2.523449e-08 is the term-weighted mix of the two.
    m_precs1  = Real(4.0) * n0s * wdm6_literal(0.65);
    m_precs2  = Real(4.0) * n0s * wdm6_literal(0.44)
                * std::pow(avts, Real(0.5)) * m_g5pbso2;
    m_pidn0s  = m_pi_wdm6 * dens_arg * n0s;
    m_pacrc   = m_pi_wdm6 * n0s * avts * m_g3pbs * Real(0.25) * eacrc;

    // Graupel/hail coefficients
    m_bvtg1   = Real(1.0) + m_bvtg;
    m_bvtg2   = Real(2.5) + Real(0.5) * m_bvtg;
    m_bvtg3   = Real(3.0) + m_bvtg;
    m_bvtg4   = Real(4.0) + m_bvtg;
    m_g1pbg   = rgmma(m_bvtg1);
    m_g3pbg   = rgmma(m_bvtg3);
    m_g4pbg   = rgmma(m_bvtg4);
    m_pacrg   = m_pi_wdm6 * m_n0g * m_avtg * m_g3pbg * Real(0.25);
    m_g5pbgo2 = rgmma(m_bvtg2);
    m_pvtg    = m_avtg * m_g4pbg / Real(6.0);
    // Fortran :3264-3265, the graupel analogues of precs1/precs2 with the same
    // unsuffixed .78 and .31. Routed on the same evidence class as precs1 and
    // precs2 rather than on separate measurements: pgdep reads bitwise in the
    // current column only because its branch is inactive there.
    m_precg1  = Real(2.0) * m_pi_wdm6 * m_n0g * wdm6_literal(0.78);
    m_precg2  = Real(2.0) * m_pi_wdm6 * m_n0g * wdm6_literal(0.31)
                * std::pow(m_avtg, Real(0.5)) * m_g5pbgo2;
    m_pidn0g  = m_pi_wdm6 * m_deng * m_n0g;

    // Slope parameter limits
    m_rslopecmax  = Real(1.0) / lamdacmax;
    m_rslopec2max = m_rslopecmax * m_rslopecmax;
    m_rslopec3max = m_rslopec2max * m_rslopecmax;

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
