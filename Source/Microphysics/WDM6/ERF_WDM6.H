#ifndef ERF_WDM6_H
#define ERF_WDM6_H

#include <memory>
#include <string>
#include <vector>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Constants.H"
#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_NullMoist.H"
#include "ERF_SurfacePrecipitation.H"

namespace MicVar_WDM6 {
   enum {
      rho = 0,
      theta,
      tabs,
      pres,
      qv,
      qc,
      qi,
      qr,
      qs,
      qg,
      nn,          // total aerosol number concentration
      nc,          // cloud droplet number concentration
      nr,          // rain drop number concentration
      rain_accum,
      snow_accum,
      graup_accum,
      NumVars
   };
}

// Selects how the WDM6 parameter literals are interpreted. See the LITERAL
// PRECISION CONTRACT comment on the constant block below for the full
// rationale. Namespace scope, not a class member, because a constexpr member
// function cannot be used in a constant expression within its own class body.
constexpr amrex::Real wdm6_literal (double d)
{
#ifdef ERF_WDM6_F32_LITERALS
    // Match the Fortran as compiled today: its literals carry no kind suffix,
    // so they are evaluated in single precision and then widened.
    return amrex::Real(static_cast<float>(d));
#else
    // True doubles. Correct when the Fortran is built -fdefault-real-8.
    return amrex::Real(d);
#endif
}

class WDM6 : public NullMoist {
    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    WDM6() {}
    virtual ~WDM6() = default;

    void Define(SolverChoice& sc) override
    {
        m_moisture_type = sc.moisture_type;
        m_axis = sc.ave_plane;
        m_do_cond = (!sc.uses_shoc_family());
    }

    // WDM6 module parameter constants (from WRF module_mp_wdm6.F)
    //
    // LITERAL PRECISION CONTRACT
    // --------------------------
    // Every constant below is initialised through wdm6_literal() rather than
    // written directly, because the two implementations do not agree on what
    // these literals mean.
    //
    // In ERF_module_mp_wdm6.F90 the parameters are declared
    //     real(kind=kind_phys), parameter :: dimax = 500.e-6
    // with NO kind suffix on the literal. Fortran therefore evaluates the
    // literal as a default (single-precision) real and only then widens it to
    // kind_phys, so the stored value carries float32 rounding permanently.
    // All 22 parameters in that module are written this way; none is suffixed.
    // dimax for instance is +5.00000023748725652695E-04 in the Fortran, which
    // is bit-for-bit float32(5e-4) widened, against a true double
    // +5.00000000000000010408E-04 here.
    //
    // That single 2.37e-11 difference is the confirmed root cause of the first
    // materially failing group in the tag-frontier retreat: it clamps
    // diameter, which drives pidep in G13e, which amplifies into psaut in G13f
    // and on into the Milestone A and B plotfile failures.
    //
    // ERF_WDM6_F32_LITERALS selects which side to match:
    //   defined    reproduce the Fortran-as-compiled-today values bit-for-bit
    //   undefined  use true doubles, which is correct when the Fortran side is
    //              built with -fdefault-real-8 so its literals are doubles too
    //
    // Four constants (actk, bvtr, bvts, satmax) previously carried a hand-added
    // 'f' suffix. That was an undocumented partial fix covering 4 of 22 and it
    // would have silently fought a -fdefault-real-8 build. Those suffixes are
    // absorbed here so the choice is made in exactly one place.
    //
    // Note r0, peaut and alpha_wdm6 have no matching Fortran parameter of the
    // same name and value, so their correspondence is unverified; they are
    // routed through the same helper for uniformity, not because a match is
    // established.
    //
    // wdm6_literal is defined at namespace scope above rather than as a member:
    // a constexpr member function cannot be called in a constant expression
    // inside its own class body, because the class is incomplete until the
    // closing brace.
    static constexpr amrex::Real dtcldcr    = wdm6_literal(120.0);
    static constexpr amrex::Real n0r        = wdm6_literal(8.0e6);
    static constexpr amrex::Real avtr       = wdm6_literal(841.9);
    static constexpr amrex::Real bvtr       = wdm6_literal(0.8);
    static constexpr amrex::Real r0         = wdm6_literal(0.8e-5);
    static constexpr amrex::Real peaut      = wdm6_literal(0.55);
    static constexpr amrex::Real xncr       = wdm6_literal(3.0e8);
    static constexpr amrex::Real xncr0      = wdm6_literal(5.0e7);
    static constexpr amrex::Real xncr1      = wdm6_literal(5.0e8);
    static constexpr amrex::Real xmyu       = wdm6_literal(1.718e-5);
    static constexpr amrex::Real avts       = wdm6_literal(11.72);
    static constexpr amrex::Real bvts       = wdm6_literal(0.41);
    static constexpr amrex::Real lamdarmax  = wdm6_literal(5.0e4);
    static constexpr amrex::Real lamdarmin  = wdm6_literal(2.0e3);
    static constexpr amrex::Real lamdasmax  = wdm6_literal(1.0e5);
    static constexpr amrex::Real lamdacmax  = wdm6_literal(5.0e5);
    static constexpr amrex::Real lamdacmin  = wdm6_literal(2.0e4);
    static constexpr amrex::Real dicon      = wdm6_literal(11.9);
    static constexpr amrex::Real dimax      = wdm6_literal(500.0e-6);
    static constexpr amrex::Real pfrz1      = wdm6_literal(100.0);
    static constexpr amrex::Real pfrz2      = wdm6_literal(0.66);
    static constexpr amrex::Real qcrmin     = wdm6_literal(1.0e-9);
    static constexpr amrex::Real ncmin      = wdm6_literal(1.0e1);
    static constexpr amrex::Real nrmin      = wdm6_literal(1.0e-2);
    static constexpr amrex::Real eacrc      = wdm6_literal(1.0);
    static constexpr amrex::Real dens_snow  = wdm6_literal(100.0);
    static constexpr amrex::Real qs0        = wdm6_literal(6.0e-4);
    static constexpr amrex::Real n0smax     = wdm6_literal(1.0e11);
    static constexpr amrex::Real n0s        = wdm6_literal(2.0e6);
    static constexpr amrex::Real alpha_wdm6 = wdm6_literal(0.12);
    static constexpr amrex::Real satmax     = wdm6_literal(1.0048);
    static constexpr amrex::Real actk       = wdm6_literal(0.6);
    static constexpr amrex::Real actr       = wdm6_literal(1.5);
    static constexpr amrex::Real ncrk1      = wdm6_literal(3.03e3);
    static constexpr amrex::Real ncrk2      = wdm6_literal(2.59e15);
    static constexpr amrex::Real di100      = wdm6_literal(1.0e-4);
    static constexpr amrex::Real di600      = wdm6_literal(6.0e-4);
    static constexpr amrex::Real di2000     = wdm6_literal(2000.0e-6);
    static constexpr amrex::Real di82       = wdm6_literal(82.0e-6);
    static constexpr amrex::Real di15       = wdm6_literal(15.0e-6);

    void Init(const amrex::MultiFab& cons_in,
              const amrex::BoxArray& grids,
              const amrex::Geometry& geom,
              const amrex::Real& dt_advance,
              std::unique_ptr<amrex::MultiFab>& z_phys_nd,
              std::unique_ptr<amrex::MultiFab>& detJ_cc) override;

    void Set_dzmin(const amrex::Real dz_min) override { m_dzmin = dz_min; }

    // Held by pointer so a retag of the mask after this call is visible. See
    // the xland fill and the qcr selection in ERF_AdvanceWDM6.cpp for the
    // ERF-lmask-to-WRF-xland encoding difference.
    void Set_Lmask(amrex::iMultiFab* lmask) override { m_lmask = lmask; }

    void Copy_State_to_Micro(const amrex::MultiFab& cons_in) override;
    void Copy_Micro_to_State(amrex::MultiFab& cons_in) override;

    void Update_Micro_Vars(amrex::MultiFab& cons_in) override
    {
        Copy_State_to_Micro(cons_in);
    }

    void Update_State_Vars(amrex::MultiFab& cons_in,
                           const amrex::MultiFab& /*z_phys_nd*/) override
    {
        Copy_Micro_to_State(cons_in);
    }

    void Advance(const amrex::Real& dt_advance,
                 const SolverChoice& solverChoice) override;

    amrex::MultiFab* Qmoist_Ptr(const int& varIdx) override
    {
        AMREX_ALWAYS_ASSERT(varIdx < m_qmoist_size);
        return mic_fab_vars[MicVarMap[varIdx]].get();
    }

    int Qmoist_Size() override { return m_qmoist_size; }
    int Qstate_Moist_Size() override { return n_qstate_moist_size; }
    int Qstate_Moist_NumConc_Size () override { return n_qstate_moist_numconc_size; }

    void Qmoist_Restart_Vars(const SolverChoice&,
                             std::vector<int>& a_idx,
                             std::vector<std::string>& a_names) const override
    {
        a_idx = {0, 1, 2};
        a_names = {"RainAccum", "SnowAccum", "GraupAccum"};
    }

    SurfacePrecipAccumulationSources
    Get_Surface_Precip_Accumulation_Ptrs (const int&) const override
    {
        SurfacePrecipAccumulationSources sources;
        sources.total   = {mic_fab_vars[MicVar_WDM6::rain_accum].get(), rhoh2o / amrex::Real(1000.0)};
        sources.snow    = {mic_fab_vars[MicVar_WDM6::snow_accum].get(), rhoh2o / amrex::Real(1000.0)};
        sources.graupel = {mic_fab_vars[MicVar_WDM6::graup_accum].get(), rhoh2o / amrex::Real(1000.0)};
        return sources;
    }

private:
    int m_qmoist_size = 3;  // rain_accum, snow_accum, graup_accum
    int n_qstate_moist_size = 9;  // 6 mixing ratios + 3 number concentrations
    int n_qstate_moist_numconc_size = 3;  // nn, nc, nr

    amrex::Vector<int> MicVarMap;

    amrex::Geometry m_geom;
    amrex::Real dt{0.0};
    amrex::Real m_dzmin{0.0};
    amrex::Real m_ccn0{100.0e6};  // background CCN concentration (m^-3)
    int nlev{0}, zlo{0}, zhi{0};
    int m_axis{2};
    bool m_do_cond{true};
    MoistureType m_moisture_type{MoistureType::None};

    amrex::MultiFab* m_z_phys_nd{nullptr};
    amrex::MultiFab* m_detJ_cc{nullptr};
    amrex::iMultiFab* m_lmask{nullptr};  // land/water mask from ERF::lmask_lev

    amrex::Array<FabPtr, MicVar_WDM6::NumVars> mic_fab_vars;

    // No nn-initialized flag: Copy_State_to_Micro decides per cell from RhoQ8,
    // which is restart-correct and needs no cross-call bookkeeping.

    // WDM6 init coefficients (hail_opt-dependent)
    bool m_hail_opt{false};
    amrex::Real m_n0g{0}, m_deng{0}, m_avtg{0}, m_bvtg{0}, m_lamdagmax{0};

    void initialize_coeffs();

    // Derived coefficients (computed in initialize_coeffs())
    amrex::Real m_pi_wdm6{0}, m_xlv1{0};
    amrex::Real m_qc0{0}, m_qc1{0}, m_qck1{0}, m_pidnc{0};
    amrex::Real m_bvtr1{0}, m_bvtr2{0}, m_bvtr3{0}, m_bvtr4{0}, m_bvtr5{0}, m_bvtr6{0}, m_bvtr7{0};
    amrex::Real m_bvtr2o5{0}, m_bvtr3o5{0};
    amrex::Real m_g1pbr{0}, m_g2pbr{0}, m_g3pbr{0}, m_g4pbr{0}, m_g5pbr{0}, m_g6pbr{0}, m_g7pbr{0};
    amrex::Real m_g5pbro2{0}, m_g7pbro2{0};
    amrex::Real m_pvtr{0}, m_pvtrn{0}, m_eacrr{0}, m_pacrr{0};
    amrex::Real m_precr1{0}, m_precr2{0}, m_xmmax{0}, m_roqimax{0};
    amrex::Real m_pidn0r{0}, m_pidnr{0};
    amrex::Real m_bvts1{0}, m_bvts2{0}, m_bvts3{0}, m_bvts4{0};
    amrex::Real m_g1pbs{0}, m_g3pbs{0}, m_g4pbs{0}, m_g5pbso2{0};
    amrex::Real m_pvts{0}, m_pacrs{0}, m_precs1{0}, m_precs2{0};
    amrex::Real m_pidn0s{0}, m_pacrc{0};
    amrex::Real m_bvtg1{0}, m_bvtg2{0}, m_bvtg3{0}, m_bvtg4{0};
    amrex::Real m_g1pbg{0}, m_g3pbg{0}, m_g4pbg{0}, m_g5pbgo2{0};
    amrex::Real m_pvtg{0}, m_pacrg{0}, m_precg1{0}, m_precg2{0}, m_pidn0g{0};
    amrex::Real m_rslopecmax{0}, m_rslopec2max{0}, m_rslopec3max{0};
    amrex::Real m_rslopermax{0}, m_rslopesmax{0}, m_rslopegmax{0};
    amrex::Real m_rsloperbmax{0}, m_rslopesbmax{0}, m_rslopegbmax{0};
    amrex::Real m_rsloper2max{0}, m_rslopes2max{0}, m_rslopeg2max{0};
    amrex::Real m_rsloper3max{0}, m_rslopes3max{0}, m_rslopeg3max{0};
};

#endif
