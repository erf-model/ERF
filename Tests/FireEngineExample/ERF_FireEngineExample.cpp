/**
 * @file ERF_FireEngineExample.cpp
 * @brief The smallest complete fire engine directory, built from outside ERF
 *
 * A prescribed fire for tests: a disk centred at (fire.example.x0,
 * fire.example.y0) whose radius grows from fire.example.r0 at the constant
 * rate fire.example.ros, releasing fire.example.heat_flux [W/m2] inside it, a
 * fraction fire.example.latent_fraction of it as latent heat. It reads no
 * atmospheric field, so a run checks only the plumbing: the out-of-tree
 * build, the registry, the flux contract, the ERF plotfile fields and the
 * checkpoint. ERF's CI builds it through ERF_FIRE_ENGINE_DIRS (CMake) and
 * FIRE_HOME (GNU make).
 *
 * The flux it reports is the exact mean over the interval of the last
 * advance(): the engine takes fire.example.nsub substeps per advance and
 * averages the flux of each, as an engine that subcycles must.
 */
#include <ERF_FireEngineRegistry.H>
#include <ERF_FireGrid.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>

#include <AMReX_Math.H>

#include <fstream>
#include <sstream>
#include <memory>

using namespace amrex;

class FireEngineExample : public FireEngine
{
public:
    FireEngineExample ()
    {
        ParmParse pp("fire.example");
        pp.query("x0", m_x0);
        pp.query("y0", m_y0);
        pp.query("r0", m_r0);
        pp.query("ros", m_ros);
        pp.query("heat_flux", m_q0);
        pp.query("latent_fraction", m_latent_fraction);
        pp.query("nsub", m_nsub);
    }

    const char* name () const override { return "example"; }

    void verify_prerequisites (const ERF& /*erf*/, int /*lev*/) override
    {
        if (m_r0 < 0.0_rt || m_ros < 0.0_rt || m_q0 < 0.0_rt || m_nsub < 1 ||
            m_latent_fraction < 0.0_rt || m_latent_fraction > 1.0_rt) {
            Abort("[FIRE] fire.example: r0, ros and heat_flux must be >= 0, nsub >= 1"
                  " and latent_fraction in [0, 1]");
        }
    }

    void initialize (const ERF& /*erf*/, int /*lev*/, const FireGrid& fg,
                     const SurfaceLayer* /*surface_layer*/, const MultiFab* /*z_phys_nd_atm*/) override
    {
        m_geom = fg.geom;
        m_phi      = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 0);
        m_sensible = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 0);
        m_latent   = std::make_unique<MultiFab>(fg.ba, fg.dm, 1, 0);
        m_sensible->setVal(0.0_rt);
        m_latent->setVal(0.0_rt);
        fill_level_set();
    }

    void advance (const FireForcing& forcing) override
    {
        // Mean over [time, time + dt] of the flux of each substep
        m_sensible->setVal(0.0_rt);
        const Real dt_sub = forcing.dt / Real(m_nsub);
        for (int n = 0; n < m_nsub; ++n) {
            add_disk_flux(*m_sensible, forcing.time + (Real(n) + 0.5_rt) * dt_sub, 1.0_rt / Real(m_nsub));
        }
        MultiFab::Copy(*m_latent, *m_sensible, 0, 0, 1, 0);
        m_sensible->mult(1.0_rt - m_latent_fraction);
        m_latent->mult(m_latent_fraction);
        m_time = forcing.time + forcing.dt;
        fill_level_set();
    }

    void fill_fluxes (FireFluxes& fluxes) override
    {
        fluxes.sensible = m_sensible.get();
        fluxes.latent   = (m_latent_fraction > 0.0_rt) ? m_latent.get() : nullptr;
    }

    Real smoke_heat_per_kg () const override { return 1.8e7_rt; }

    // A prescribed fire reads nothing from the surface layer
    bool needs_surface_layer () const override { return false; }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream os(checkpointname + "/FireEngineExample");
            os.precision(17);
            os << m_time << "\n";
        }
    }

    void read_checkpoint (const std::string& restart_chkfile, int /*max_level*/, int /*step*/) override
    {
        Vector<char> buf;
        ParallelDescriptor::ReadAndBcastFile(restart_chkfile + "/FireEngineExample", buf);
        std::istringstream is(std::string(buf.dataPtr()));
        is >> m_time;
        fill_level_set();
    }

    std::vector<std::string> plot_var_names () const override
    {
        return {"fire_phi", "fire_heat_flux"};
    }

    const MultiFab* plot_var (const std::string& name) const override
    {
        if (name == "fire_phi")       { return m_phi.get(); }
        if (name == "fire_heat_flux") { return m_sensible.get(); }
        return nullptr;
    }

    const MultiFab* level_set () const override { return m_phi.get(); }

    Real burned_area () const override
    {
        const Real r = radius(m_time);
        return amrex::Math::pi<Real>() * r * r;
    }

private:
    Real radius (Real t) const { return m_r0 + m_ros * t; }

    void fill_level_set ()
    {
        const auto dx  = m_geom.CellSizeArray();
        const auto plo = m_geom.ProbLoArray();
        const Real x0 = m_x0, y0 = m_y0, r = radius(m_time);
        for (MFIter mfi(*m_phi); mfi.isValid(); ++mfi) {
            const auto phi = m_phi->array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real x = plo[0] + (i + 0.5_rt) * dx[0];
                const Real y = plo[1] + (j + 0.5_rt) * dx[1];
                phi(i, j, k) = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0)) - r;
            });
        }
    }

    void add_disk_flux (MultiFab& q, Real t, Real weight) const
    {
        const auto dx  = m_geom.CellSizeArray();
        const auto plo = m_geom.ProbLoArray();
        const Real x0 = m_x0, y0 = m_y0, r = radius(t), q0 = m_q0 * weight;
        for (MFIter mfi(q); mfi.isValid(); ++mfi) {
            const auto qa = q.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real x = plo[0] + (i + 0.5_rt) * dx[0];
                const Real y = plo[1] + (j + 0.5_rt) * dx[1];
                if ((x - x0) * (x - x0) + (y - y0) * (y - y0) <= r * r) { qa(i, j, k) += q0; }
            });
        }
    }

    Real m_x0 = 0.0_rt, m_y0 = 0.0_rt, m_r0 = 50.0_rt, m_ros = 0.5_rt;
    Real m_q0 = 5.0e4_rt, m_latent_fraction = 0.0_rt;
    int  m_nsub = 1;
    Real m_time = 0.0_rt;
    Geometry m_geom;
    std::unique_ptr<MultiFab> m_phi, m_sensible, m_latent;
};

ERF_REGISTER_FIRE_ENGINE(FireEngineExample, "example")
