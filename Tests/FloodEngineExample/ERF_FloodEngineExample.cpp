/**
 * @file ERF_FloodEngineExample.cpp
 * @brief The smallest complete flood engine directory, built from outside ERF
 *
 * A bucket for tests: each flood cell keeps the water that has fallen on it,
 * with no runoff, infiltration or evaporation. Each step it adds the step's
 * surface precipitation (liquid and frozen, as liquid water) to its depth,
 *
 *     h += (P_liquid + P_frozen) / rho_w,   rho_w = 1000 kg/m3,
 *
 * with P in kg/m2, so the depth in m is the accumulated precipitation in
 * kg/m2 divided by 1000. It writes the depth to the checkpoint and reads it
 * back on a restart, writes it to its own plotfile and offers it to the ERF
 * plotfile as flood_depth. A run checks the plumbing: the out-of-tree build,
 * the registry, the advance with the precipitation, the plotfiles and the
 * checkpoint. The CTest compares flood_depth with the microphysics' own
 * accumulation in the same plotfile. It builds through
 * ERF_FLOOD_ENGINE_DIRS (CMake) and FLOOD_HOME (GNU make); ERF's CI does
 * not build it.
 */
#include <ERF_FloodEngineRegistry.H>
#include <ERF_FloodGrid.H>

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <memory>
#include <string>
#include <vector>

using namespace amrex;

class FloodEngineExample : public FloodEngine
{
public:
    const char* name () const override { return "example"; }

    void verify_prerequisites (const ERF& /*erf*/, int /*lev*/) override {}

    void initialize (const ERF& /*erf*/, int lev, const FloodGrid& fg,
                     const MultiFab* /*z_phys_nd_atm*/) override
    {
        m_lev = lev;
        m_fg  = fg;
        m_depth = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 0);
        m_depth->setVal(0.0_rt);
        m_rate  = std::make_unique<MultiFab>(fg.ba, fg.dm, 2, 0);
    }

    void advance (const FloodForcing& forcing) override
    {
        AMREX_ALWAYS_ASSERT(forcing.precip_k0 != nullptr);
        // The step's precipitation of each column on every flood cell under it
        fill_flood_from_atm_k0(*m_rate, *forcing.precip_k0, m_fg, 0, 0);
        fill_flood_from_atm_k0(*m_rate, *forcing.precip_k0, m_fg, 1, 1);
        constexpr Real inv_rho_w = Real(1.0) / Real(1000.0);
        for (MFIter mfi(*m_depth); mfi.isValid(); ++mfi) {
            const auto h = m_depth->array(mfi);
            const auto p = m_rate->const_array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                h(i, j, k) += (p(i, j, k, 0) + p(i, j, k, 1)) * inv_rho_w;
            });
        }
    }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        VisMF::Write(*m_depth, MultiFabFileFullPrefix(m_lev, checkpointname, "Level_", "FloodExampleDepth"));
    }

    void read_checkpoint (const std::string& restart_chkfile, int /*step*/) override
    {
        VisMF::Read(*m_depth, MultiFabFileFullPrefix(m_lev, restart_chkfile, "Level_", "FloodExampleDepth"));
    }

    void write_plotfile (const std::string& plotfile_prefix, Real time, int step) const override
    {
        WriteSingleLevelPlotfile(Concatenate(plotfile_prefix, step, 5), *m_depth,
                                 {"flood_depth"}, m_fg.geom, time, step);
    }

    std::vector<std::string> plot_var_names () const override { return {"flood_depth"}; }

    const MultiFab* plot_var (const std::string& name) const override
    {
        return (name == "flood_depth") ? m_depth.get() : nullptr;
    }

private:
    int m_lev = 0;
    FloodGrid m_fg;
    std::unique_ptr<MultiFab> m_depth;   ///< water depth [m]
    std::unique_ptr<MultiFab> m_rate;    ///< scratch: the step's liquid (0) and frozen (1) precipitation [kg/m2]
};

ERF_REGISTER_FLOOD_ENGINE(FloodEngineExample, "example")
