/**
 * @file ERF_StructureModelExample.cpp
 * @brief The smallest complete structure model directory, built from outside ERF
 *
 * A fixed body for tests: one point at structure.example.position that pushes on
 * the air with the constant force structure.example.force [N], spread onto the
 * momentum sources with the actuator core's Gaussian kernel of width
 * structure.example.epsilon [m]. It also samples the wind at the point with the
 * core's sampler and keeps the impulse it has given the air, which it writes to the
 * checkpoint and reads back on a restart. A run checks the plumbing: the
 * out-of-tree build, the registry, the advance and momentum-source calls, the
 * ERF plotfile fields and the checkpoint. In a periodic box with no friction the
 * air's momentum grows by exactly force times time, which the CTest checks.
 * When the run has a fire (StructureStep::fire), it prints each step what the
 * fire view holds (the level and time of the fire's fields, the heat released
 * and the area burned), which the CTests with the example fire engine check.
 * ERF's CI builds it through ERF_STRUCTURE_MODEL_DIRS (CMake) and
 * STRUCTURE_HOME (GNU make).
 */
#include <ERF_StructureModelRegistry.H>
#include <ERF_ActuatorSampling.H>
#include <ERF_ActuatorSpreading.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <cmath>
#include <fstream>
#include <array>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace amrex;

namespace {

// What the fire view holds this step: the level and time of the fire's fields, the heat
// released over the fire grid [W] and the area burned [m2] (-1 where the engine has no
// such field). Every rank calls it: the sum and burned_area() reduce over the ranks.
void report_fire (const StructureStep& step)
{
    const StructureFire& f = *step.fire;
    const Real dA = f.geom.CellSize(0) * f.geom.CellSize(1);
    const Real heat = f.heat_release ? f.heat_release->sum(0) * dA : Real(-1.0);
    const Real burned = f.burned_area();
    Print() << std::setprecision(17) << "structure example fire: t = " << step.time << " fire_lev = " << f.lev
            << " fire_t = " << f.time << std::setprecision(10) << " heat = " << heat << " burned = " << burned << "\n";
}

} // namespace

class StructureModelExample : public StructureModel
{
public:
    explicit StructureModelExample (int max_level)
    {
        ParmParse pp("structure.example");
        pp.getarr("position", m_pos, 0, 3);
        pp.getarr("force", m_force, 0, 3);
        pp.get("epsilon", m_eps);
        pp.query("level", m_level);
        pp.query("test_support", m_test_support);
        if (m_test_support != "reach" && m_test_support != "nan" && m_test_support != "inverted") {
            Abort("structure.example.test_support must be reach, nan or inverted, got " + m_test_support);
        }
        for (int d = 0; d < 3; ++d) {
            if (!std::isfinite(m_pos[d]) || !std::isfinite(m_force[d])) {
                Abort("structure.example.position and structure.example.force must be finite");
            }
        }
        if (!(m_eps > Real(0.0)) || !std::isfinite(m_eps)) {
            Abort("structure.example.epsilon must be > 0 and finite, got " + std::to_string(m_eps));
        }
        if (m_level < -1 || m_level > max_level) {
            Abort("structure.example.level must be -1 (the finest level) or 0 to amr.max_level = "
                  + std::to_string(max_level) + ", got " + std::to_string(m_level));
        }
    }

    const char* name () const override { return "example"; }

    int level (int finest_level) const override
    {
        if (m_level > finest_level) {
            Abort("structure.example.level = " + std::to_string(m_level)
                  + " is above the finest level of this run (" + std::to_string(finest_level) + ")");
        }
        return (m_level < 0) ? finest_level : m_level;
    }

    void initialize (const StructureInit& init) override
    {
        if (!init.unit_map_factors) {
            Abort("structure model \"example\": the actuator core needs map factors of 1");
        }
        const Real* lo = init.geom.ProbLo();
        const Real* hi = init.geom.ProbHi();
        for (int d = 0; d < 3; ++d) {
            if (m_pos[d] < lo[d] || m_pos[d] > hi[d]) {
                Abort("structure.example.position lies outside the domain");
            }
        }
        // What a model that drives a solver of its own needs before the first step. Values
        // that mean "not given" are printed as words, so the CTests do not depend on how a
        // stream formats the largest double.
        std::ostringstream os;
        os << std::setprecision(10) << "structure example: level " << init.lev;
        if (init.fixed_dt > 0.0) { os << " steps with " << init.fixed_dt << " s"; }
        else                     { os << " has no fixed step"; }
        if (init.stop_time < std::numeric_limits<double>::max()) { os << ", the run stops at t = " << init.stop_time << " s"; }
        else                                                     { os << ", the run has no stop time"; }
        if (init.max_step < std::numeric_limits<int>::max()) { os << " or after level-0 step " << init.max_step; }
        else                                                 { os << " and no step limit"; }
        os << ", fully compressible: " << (init.fully_compressible ? "yes" : "no") << "\n";
        Print() << os.str();
        // A run that takes no step (a restart at its last step) ends where it starts
        m_end_time = init.time;
        if (!init.restart_chkfile.empty()) {
            // The impulse given so far, written by write_checkpoint
            const std::string file = init.restart_chkfile + "/example/impulse";
            if (!FileExists(file)) {
                Abort("erf.structures.model = example: the checkpoint " + init.restart_chkfile
                      + " holds no state of this model (no " + file + "); it was written without it");
            }
            Vector<char> buffer;
            ParallelDescriptor::ReadAndBcastFile(file, buffer);
            std::istringstream is(std::string(buffer.dataPtr()));
            is >> m_impulse[0] >> m_impulse[1] >> m_impulse[2];
            if (!is) { Abort("structure model \"example\": cannot read " + init.restart_chkfile + "/example/impulse"); }
        }
    }

    void advance (const StructureStep& step) override
    {
        if (!step.unit_map_factors) {
            Abort("structure model \"example\": the actuator core needs map factors of 1");
        }
        // The point and its kernel (3 epsilon) must lie inside the level's boxes, also after a regrid
        std::string outside;
        if (!erf_actuator::points_covered_by(step.xvel.boxArray(), step.geom, m_pos, Real(3.0) * m_eps, outside,
                                             step.z_phys_nd ? erf_actuator::CoverZ::Footprint
                                                            : erf_actuator::CoverZ::Reach)) {
            Abort("structure.example.position with 3 structure.example.epsilon around it is not inside the boxes of"
                  " level " + std::to_string(step.lev) + " (structure.example.level): " + outside);
        }
        m_lev = step.lev;

        // The wind at the point, from each component's own staggered grid
        erf_actuator::sample_velocity(step.xvel, step.yvel, step.zvel, step.z_phys_nd, step.geom,
                                      m_pos, m_wind);

        // The force on the air, spread once per step and held over its stages
        const BoxArray ba = amrex::convert(step.xvel.boxArray(), IntVect(0));
        const DistributionMapping& dm = step.xvel.DistributionMap();
        if (m_src_x.boxArray() != amrex::convert(ba, IntVect(1,0,0)) || m_src_x.DistributionMap() != dm) {
            m_src_x.define(amrex::convert(ba, IntVect(1,0,0)), dm, 1, 0);
            m_src_y.define(amrex::convert(ba, IntVect(0,1,0)), dm, 1, 0);
            m_src_z.define(amrex::convert(ba, IntVect(0,0,1)), dm, 1, 0);
        }
        erf_actuator::spread_forces(m_pos, m_force, m_eps, step.z_phys_nd, step.detJ_cc, step.geom,
                                    m_src_x, m_src_y, m_src_z);
        if (step.fire) { report_fire(step); }
        for (int d = 0; d < 3; ++d) { m_impulse[d] += static_cast<double>(m_force[d]) * step.dt; }
        ++m_steps;
        m_end_time = step.time + step.dt;

        Print() << std::setprecision(10) << "structure example: t = " << step.time + step.dt
                << " impulse_x = " << m_impulse[0] << " wind_x = " << m_wind[0] << "\n";
    }

    std::vector<RealBox> support () const override
    {
        // The point with its kernel's reach (3 epsilon), as advance() checks the coverage
        const Real r = Real(3.0) * m_eps;
        // For ERF's tests of its check of these boxes: a box no model should return
        if (m_test_support == "nan") {
            const Real nan = std::numeric_limits<Real>::quiet_NaN();
            return {RealBox(m_pos[0] - r, nan, m_pos[2] - r, m_pos[0] + r, m_pos[1] + r, m_pos[2] + r)};
        }
        if (m_test_support == "inverted") {
            return {RealBox(m_pos[0] + r, m_pos[1] - r, m_pos[2] - r, m_pos[0] - r, m_pos[1] + r, m_pos[2] + r)};
        }
        return {RealBox(m_pos[0] - r, m_pos[1] - r, m_pos[2] - r, m_pos[0] + r, m_pos[1] + r, m_pos[2] + r)};
    }

    void add_momentum_sources (int /*lev*/, MultiFab& xmom_src, MultiFab& ymom_src, MultiFab& zmom_src) const override
    {
        MultiFab::Add(xmom_src, m_src_x, 0, 0, 1, 0);
        MultiFab::Add(ymom_src, m_src_y, 0, 0, 1, 0);
        MultiFab::Add(zmom_src, m_src_z, 0, 0, 1, 0);
    }

    void finalize () override
    {
        Print() << std::setprecision(10) << "structure example: finished at t = " << m_end_time << " after "
                << m_steps << " steps of this run\n";
    }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        const std::string dir = checkpointname + "/example";
        if (ParallelDescriptor::IOProcessor()) {
            if (!UtilCreateDirectory(dir, 0755)) { CreateDirectoryFailed(dir); }
            std::ofstream os(dir + "/impulse");
            os << std::setprecision(std::numeric_limits<double>::max_digits10)
               << m_impulse[0] << " " << m_impulse[1] << " " << m_impulse[2] << "\n";
        }
        ParallelDescriptor::Barrier();
    }

    bool supports_restart () const override { return true; }

    std::vector<std::string> plot_var_names () const override
    {
        return {"structure_example_fx", "structure_example_fy", "structure_example_fz"};
    }

    void fill_plot_var (const std::string& name, int lev, MultiFab& dst, int comp) const override
    {
        // The face source averaged to the cell centres; zero before the first step
        // and on every level but the one the source was spread on
        const int dir = (name == "structure_example_fx") ? 0 : (name == "structure_example_fy") ? 1 : 2;
        const MultiFab& src = (dir == 0) ? m_src_x : (dir == 1) ? m_src_y : m_src_z;
        dst.setVal(Real(0.0), comp, 1, 0);
        if (lev != m_lev || !src.isDefined() || amrex::convert(src.boxArray(), IntVect(0)) != dst.boxArray()
            || src.DistributionMap() != dst.DistributionMap()) { return; }
        for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const auto s = src.const_array(mfi);
            const auto d = dst.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const IntVect iv(i, j, k);
                const IntVect ivp = iv + IntVect::TheDimensionVector(dir);
                d(i, j, k, comp) = Real(0.5) * (s(iv) + s(ivp));
            });
        }
    }

private:
    std::vector<Real> m_pos;
    std::vector<Real> m_force;
    Real m_eps = Real(0.0);
    int m_level = -1;
    std::string m_test_support = "reach";   // what support() returns: the kernel's reach, or a box ERF must refuse
    int m_lev = -1;   // the level of the last advance()
    std::array<double,3> m_impulse{{0.0, 0.0, 0.0}};
    int m_steps = 0;            // steps advanced in this run (not checkpointed)
    double m_end_time = 0.0;    // end of the last step advanced, or the start time [s]
    std::vector<Real> m_wind{Real(0.0), Real(0.0), Real(0.0)};
    MultiFab m_src_x, m_src_y, m_src_z;
};

ERF_REGISTER_STRUCTURE_MODEL(StructureModelExample, "example")
