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
        m_axis = sc.ave_plane;
        m_do_cond = (!sc.use_shoc);
    }

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

    amrex::Vector<int> MicVarMap;

    amrex::Geometry m_geom;
    amrex::Real dt{0.0};
    amrex::Real m_dzmin{0.0};
    int nlev{0}, zlo{0}, zhi{0};
    int m_axis{2};
    bool m_do_cond{true};

    amrex::MultiFab* m_z_phys_nd{nullptr};
    amrex::MultiFab* m_detJ_cc{nullptr};

    amrex::Array<FabPtr, MicVar_WSM6::NumVars> mic_fab_vars;
};

#endif
