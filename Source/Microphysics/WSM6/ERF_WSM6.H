#ifndef ERF_WSM6_H
#define ERF_WSM6_H

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

namespace MicVar_WSM6 {
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
      rain_accum,
      snow_accum,
      graup_accum,
      NumVars
   };
}

class WSM6 : public NullMoist {
    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    WSM6() {}
    virtual ~WSM6() = default;

    void Define(SolverChoice& sc) override
    {
        m_do_cond = (!sc.use_shoc);
    }

    // WSM6 module parameter constants (Fortran: parameter, private)
    static constexpr amrex::Real dtcldcr    = amrex::Real(120.0);
    static constexpr amrex::Real n0r        = amrex::Real(8.0e6);
    static constexpr amrex::Real avtr       = amrex::Real(841.9);
    static constexpr amrex::Real bvtr       = amrex::Real(0.8);
    static constexpr amrex::Real r0         = amrex::Real(0.8e-5);
    static constexpr amrex::Real peaut      = amrex::Real(0.55);
    static constexpr amrex::Real xncr       = amrex::Real(3.0e8);
    static constexpr amrex::Real xmyu       = amrex::Real(1.718e-5);
    static constexpr amrex::Real avts       = amrex::Real(11.72);
    static constexpr amrex::Real bvts       = amrex::Real(0.41);
    static constexpr amrex::Real lamdarmax  = amrex::Real(8.0e4);
    static constexpr amrex::Real lamdasmax  = amrex::Real(1.0e5);
    static constexpr amrex::Real dicon      = amrex::Real(11.9);
    static constexpr amrex::Real dimax      = amrex::Real(500.0e-6);
    static constexpr amrex::Real pfrz1      = amrex::Real(100.0);
    static constexpr amrex::Real pfrz2      = amrex::Real(0.66);
    static constexpr amrex::Real qcrmin     = amrex::Real(1.0e-9);
    static constexpr amrex::Real eacrc      = amrex::Real(1.0);
    static constexpr amrex::Real dens_snow  = amrex::Real(100.0);
    static constexpr amrex::Real qs0        = amrex::Real(6.0e-4);
    static constexpr amrex::Real n0smax     = amrex::Real(1.0e11);
    static constexpr amrex::Real n0s        = amrex::Real(2.0e6);
    static constexpr amrex::Real alpha_wsm6 = amrex::Real(0.12);

    void Init(const amrex::MultiFab& cons_in,
              const amrex::BoxArray& grids,
              const amrex::Geometry& geom,
              const amrex::Real& dt_advance,
              std::unique_ptr<amrex::MultiFab>& z_phys_nd,
              std::unique_ptr<amrex::MultiFab>& detJ_cc) override;

    void Set_dzmin(const amrex::Real dz_min) override { m_dzmin = dz_min; }

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

private:
    int m_qmoist_size = 3;
    int n_qstate_moist_size = 6;
    int n_qstate_moist_numconc_size = 0;

    amrex::Vector<int> MicVarMap;

    amrex::Geometry m_geom;
    amrex::Real dt{0.0};
    amrex::Real m_dzmin{0.0};
    int nlev{0}, zlo{0}, zhi{0};
    bool m_do_cond{true};

    amrex::MultiFab* m_z_phys_nd{nullptr};
    amrex::MultiFab* m_detJ_cc{nullptr};

    amrex::Array<FabPtr, MicVar_WSM6::NumVars> mic_fab_vars;

    // ---- WSM6 init coefficients (from mp_wsm6_init save variables) ----
    // hail_opt-dependent (set in Initialize() based on m_hail_opt)
    bool m_hail_opt{false};
    amrex::Real m_n0g{0}, m_deng{0}, m_avtg{0}, m_bvtg{0}, m_lamdagmax{0};

    void initialize_coeffs();

    // Derived coefficients (computed in initialize_coeffs(), called from Init())
    amrex::Real m_pi_wsm6{0}, m_xlv1{0};
    amrex::Real m_qc0{0}, m_qck1{0}, m_pidnc{0};
    amrex::Real m_bvtr1{0}, m_bvtr2{0}, m_bvtr3{0}, m_bvtr4{0}, m_bvtr6{0};
    amrex::Real m_g1pbr{0}, m_g3pbr{0}, m_g4pbr{0}, m_g5pbro2{0}, m_g6pbr{0};
    amrex::Real m_pvtr{0}, m_eacrr{0}, m_pacrr{0};
    amrex::Real m_precr1{0}, m_precr2{0}, m_roqimax{0};
    amrex::Real m_bvts1{0}, m_bvts2{0}, m_bvts3{0}, m_bvts4{0};
    amrex::Real m_g1pbs{0}, m_g3pbs{0}, m_g4pbs{0}, m_g5pbso2{0};
    amrex::Real m_pvts{0}, m_pacrs{0}, m_precs1{0}, m_precs2{0};
    amrex::Real m_pidn0r{0}, m_pidn0s{0}, m_pacrc{0};
    amrex::Real m_bvtg1{0}, m_bvtg2{0}, m_bvtg3{0}, m_bvtg4{0};
    amrex::Real m_g1pbg{0}, m_g3pbg{0}, m_g4pbg{0}, m_g5pbgo2{0};
    amrex::Real m_pvtg{0}, m_pacrg{0}, m_precg1{0}, m_precg2{0}, m_pidn0g{0};
    amrex::Real m_rslopermax{0}, m_rslopesmax{0}, m_rslopegmax{0};
    amrex::Real m_rsloperbmax{0}, m_rslopesbmax{0}, m_rslopegbmax{0};
    amrex::Real m_rsloper2max{0}, m_rslopes2max{0}, m_rslopeg2max{0};
    amrex::Real m_rsloper3max{0}, m_rslopes3max{0}, m_rslopeg3max{0};
};

namespace WSM6Ind {
    enum {
        // --- inout state (mirrors mp_wsm6_run arguments) ---
        t = 0,
        q,
        qc,
        qi,
        qr,
        qs,
        qg,
        den,
        p,
        delz,

        // --- sedimentation working arrays ---
        fallc,
        falkc,
        work1c,
        work2c,
        workr,
        worka,
        den_tmp,
        delz_tmp,

        // --- microphysics process rates ---
        pigen,
        pidep,
        pcond,
        prevp,
        psevp,
        pgevp,
        psdep,
        pgdep,
        praut,
        psaut,
        pgaut,
        piacr,
        pracw,
        praci,
        pracs,
        psacw,
        psaci,
        psacr,
        pgacw,
        pgaci,
        pgacr,
        pgacs,
        paacw,
        psmlt,
        pgmlt,
        pseml,
        pgeml,

        // --- derived thermodynamic fields ---
        qsum,
        xl,
        cpm,
        work2,
        denfac,
        xni,
        denqrs1,
        denqrs2,
        denqrs3,
        denqci,
        n0sfac,

        // --- multi-component arrays expanded from Fortran dimension(:,:,N) ---
        // rh/qsat: N=2, w=over liquid(1), i=over ice(2) — matches Morrison qsatw/qsati
        // slope/fall/qrs/falk/work1: N=3, _r=rain(1) _s=snow(2) _g=graupel(3)
        // workdiff: w=liquid diffusion work1(1), i=ice diffusion work1(2) — separate from work1_r/s/g
        rhw,  rhi,       // relative humidity: w=over liquid, i=over ice (Fortran rh(:,:,1/2))
        qsatw, qsati,    // saturation mixing ratio: w=over liquid, i=over ice (Fortran qsat(:,:,1/2))
        rslope_r,  rslope_s,  rslope_g,
        rslope2_r, rslope2_s, rslope2_g,
        rslope3_r, rslope3_s, rslope3_g,
        rslopeb_r, rslopeb_s, rslopeb_g,
        qrs_tmp_r, qrs_tmp_s, qrs_tmp_g,
        falk_r,  falk_s,  falk_g,
        fall_r,  fall_s,  fall_g,
        work1_r, work1_s, work1_g,  // fall speeds post slope_wsm6: 1=rain 2=snow 3=graupel
        workdiffw, workdiffi,        // diffusion denominators post slope_wsm6 call 2: w=liq i=ice

        NumInds
    };
}

constexpr int WSM6_MAX_LEVELS = 256;

#endif
