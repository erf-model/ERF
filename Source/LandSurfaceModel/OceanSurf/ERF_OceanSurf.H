#ifndef ERF_OCEANSURF_H
#define ERF_OCEANSURF_H

#include <memory>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <ERF_NullSurf.H>

namespace LsmVar_OceanSurf {
    enum {
        t_surf = 0,
        NumVars
    };
}

class OceanSurf : public NullSurf {
public:
    OceanSurf () = default;
    ~OceanSurf () override = default;

    void
    Define (SolverChoice& /*sc*/) override
    {
        amrex::ParmParse pp("erf");
        m_default_tsurf = amrex::Real(300.0);
        pp.query("most.surf_temp", m_default_tsurf);
    }

    void
    Init (const int& /*lev*/,
          const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Real& dt) override
    {
        m_dt = dt;
        m_lsm_geom = geom;

        amrex::Box domain = geom.Domain();
        khi_lsm    = domain.smallEnd(2);

        LsmVarMap.resize(m_lsm_size);
        LsmVarMap = {LsmVar_OceanSurf::t_surf};

        LsmVarName.resize(m_lsm_size);
        LsmVarName = {"t_surf"};

        amrex::IntVect ng(1,1,0);
        amrex::BoxArray ba = cons_in.boxArray();
        amrex::DistributionMapping dm = cons_in.DistributionMap();
        amrex::BoxList bl_lsm = ba.boxList();
        for (auto& b : bl_lsm) {
            b.setRange(2,0);
        }
        amrex::BoxArray ba_lsm(std::move(bl_lsm));

        const amrex::RealBox& dom_rb = geom.ProbDomain();
        amrex::RealBox lsm_rb = dom_rb;
        lsm_rb.setHi(2, dom_rb.hi(2));
        lsm_rb.setLo(2, dom_rb.lo(2));
        m_lsm_geom.define(ba_lsm.minimalBox(), lsm_rb, geom.Coord(), geom.isPeriodic());

        lsm_fab_data.resize(m_lsm_size);
        lsm_fab_data[0] = std::make_shared<amrex::MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_data[0]->setVal(m_default_tsurf);
    }

    void
    Advance (const amrex::Real& dt) override
    {
        // OceanSurf currently does not evolve internal state here; keep the
        // latest timestep for bookkeeping so this stays aligned with other LSMs.
        m_dt = dt;
    }

    amrex::MultiFab*
    Lsm_Data_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < OceanSurf::m_lsm_size && lsmIdx>=0);
        return lsm_fab_data[lsmIdx].get();
    }

    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& /*varIdx*/) override
    {
        return nullptr;
    }

    amrex::Geometry
    Lsm_Geom () override { return m_lsm_geom; }

    int
    Lsm_Data_Size () override { return OceanSurf::m_lsm_size; }

    int
    Lsm_Flux_Size () override { return 0; }

    std::string
    Lsm_DataName (const int& varIdx) override
    {
        int lsmIdx = LsmVarMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < OceanSurf::m_lsm_size && lsmIdx>=0);
        return LsmVarName[lsmIdx];
    }

    std::string
    Lsm_FluxName (const int& /*varIdx*/) override
    {
        return "";
    }

    int
    Lsm_DataIndex (std::string varname) override
    {
        std::string lc_varname = amrex::toLower(varname);
        if (lc_varname == "t_surf" || lc_varname == "theta") {
            return 0;
        }
        return -1;
    }

    int
    Lsm_FluxIndex (std::string /*varname*/) override
    {
        return -1;
    }

    void
    SetOceanState (const amrex::Vector<amrex::MultiFab*>* ocean_state)
    {
        // This pointer is currently packed and filled by ERF_Coupling.cpp.
        m_ocean_state = ocean_state;
    }

    [[nodiscard]]
    const amrex::Vector<amrex::MultiFab*>* OceanState () const
    {
        return m_ocean_state;
    }

private:
    using FabPtr = std::shared_ptr<amrex::MultiFab>;

    int m_lsm_size = 1;
    amrex::Vector<int> LsmVarMap;
    amrex::Vector<std::string> LsmVarName;
    amrex::Vector<FabPtr> lsm_fab_data;
    amrex::Geometry m_lsm_geom;
    amrex::Real m_dt{0.0};
    int khi_lsm{0};
    amrex::Real m_default_tsurf{amrex::Real(300.0)};
    const amrex::Vector<amrex::MultiFab*>* m_ocean_state = nullptr;
};
#endif
