#ifndef ERF_MASS_SOURCE_H_
#define ERF_MASS_SOURCE_H_

#include <functional>
#include <cmath>
#include <memory>
#include <set>
#include <sstream>
#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include "ERF_IndexDefines.H"

/**
 * A localized source of gas inside the domain (erf.mass_source.*), e.g. the release
 * from a ruptured pipeline.
 *
 * The gas is injected at a total rate mdot(t) [kg/s], spread uniformly over the cells
 * whose centres lie within erf.mass_source.radius of (center_x, center_y) horizontally
 * and less than erf.mass_source.height above the local ground.  The same source is
 * added to
 *   - rho, so the mass enters the continuity equation (FullyCompressible only: the
 *     anelastic projection would remove it),
 *   - rho*theta, at erf.mass_source.theta if given and otherwise at the cell's own theta,
 *     so the injection by itself does not change the potential temperature, and
 *   - the first passive scalar (rho*s, plotted as rhoadv_0) when erf.mass_source.tracer
 *     is true, so that s is the mass fraction of the released gas.
 * No momentum is added: the gas enters at rest.  The gas has the properties of dry air
 * (in a moist run it carries no water, so it dilutes the mixing ratios of the cells it
 * enters).  The cell volumes include the terrain metric and the map factors.  In a
 * periodic direction the region wraps around the domain.  The rate table is in the time
 * since the start of the run.  Unknown erf.mass_source keys stop the run.
 *
 * With mesh refinement every level adds the same source density S = mdot / V to its own
 * cells in the region, where V is the volume of the region on the composite grid: the
 * cells of each level that the next finer level does not cover.  The composite rate is
 * then mdot however the levels cut the region, and a level without source cells adds
 * nothing.  V, and the height of the ground under each column, are recomputed only when
 * the grids change.
 *
 * Inputs:
 *   erf.mass_source.center    = x y          [m]
 *   erf.mass_source.radius    = r            [m]
 *   erf.mass_source.height    = h            [m above the local ground]
 *   erf.mass_source.rate      = mdot         [kg/s >= 0]  constant, or
 *   erf.mass_source.rate_file = file         two columns "t mdot" [s, kg/s >= 0], '#' comments;
 *                                            linear in between, zero outside the table
 *   erf.mass_source.theta     = theta        [K > 0]  optional
 *   erf.mass_source.tracer    = true         optional
 */
namespace erf_mass_source {

/// Is the centre of cell (i,j,k) inside the source region?  z is the nodal height z_phys_nd
/// and zg the nodal height of the ground under the column (2D, at k = kg).  period holds the
/// domain length in x and y where the domain is periodic and 0 where it is not: there the
/// distance is to the nearest periodic image of the centre, so the region wraps around
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
in_source (int i, int j, int k,
           amrex::Array4<amrex::Real const> const& z,
           amrex::Array4<amrex::Real const> const& zg, int kg,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           amrex::GpuArray<amrex::Real,2> const& period,
           amrex::Real x0, amrex::Real y0, amrex::Real r2, amrex::Real hmax) noexcept
{
    amrex::Real ddx = problo[0] + (i + amrex::Real(0.5)) * dx[0] - x0;
    amrex::Real ddy = problo[1] + (j + amrex::Real(0.5)) * dx[1] - y0;
    if (period[0] > 0) { ddx -= period[0] * std::floor(ddx / period[0] + amrex::Real(0.5)); }
    if (period[1] > 0) { ddy -= period[1] * std::floor(ddy / period[1] + amrex::Real(0.5)); }
    if (ddx * ddx + ddy * ddy > r2) { return false; }
    const amrex::Real g  = amrex::Real(0.25) * (zg(i,j,kg) + zg(i+1,j,kg) + zg(i,j+1,kg) + zg(i+1,j+1,kg));
    const amrex::Real zc = amrex::Real(0.125) * (z(i,j,k  ) + z(i+1,j,k  ) + z(i,j+1,k  ) + z(i+1,j+1,k  ) +
                                                 z(i,j,k+1) + z(i+1,j,k+1) + z(i,j+1,k+1) + z(i+1,j+1,k+1));
    return zc - g < hmax;
}

/// The domain length in x and y where geom is periodic, 0 where it is not
inline amrex::GpuArray<amrex::Real,2>
periodic_lengths (const amrex::Geometry& geom)
{
    return {geom.isPeriodic(0) ? geom.ProbLength(0) : amrex::Real(0),
            geom.isPeriodic(1) ? geom.ProbLength(1) : amrex::Real(0)};
}

} // namespace erf_mass_source

class MassSource
{
public:
    /// The source the inputs ask for, or nullptr when no erf.mass_source input is given.
    /// Stops on an unknown erf.mass_source key (a misspelt one would otherwise be ignored,
    /// and a misspelt rate would turn the source off) and on inputs without a rate.
    static std::unique_ptr<MassSource> from_inputs ()
    {
        const std::string prefix = "erf.mass_source.";
        const std::set<std::string> known = {"rate", "rate_file", "center", "radius",
                                             "height", "theta", "tracer"};
        bool any = false;
        for (const auto& name : amrex::ParmParse::getEntries("erf.mass_source")) {
            if (name.rfind(prefix, 0) != 0) { continue; }
            any = true;
            if (known.count(name.substr(prefix.size())) == 0) {
                amrex::Abort("erf.mass_source: unknown input " + name + "; the inputs are rate, "
                             "rate_file, center, radius, height, theta and tracer");
            }
        }
        amrex::ParmParse pp("erf.mass_source");
        if (!pp.contains("rate") && !pp.contains("rate_file")) {
            if (any) {
                amrex::Abort("erf.mass_source: inputs are given but neither erf.mass_source.rate "
                             "nor erf.mass_source.rate_file");
            }
            return nullptr;
        }
        if (pp.contains("rate") && pp.contains("rate_file")) {
            amrex::Abort("erf.mass_source: give rate or rate_file, not both");
        }

        auto ms = std::make_unique<MassSource>();
        amrex::Vector<amrex::Real> c;
        pp.getarr("center", c);
        if (c.size() != 2) { amrex::Abort("erf.mass_source.center must be two values: x y"); }
        ms->m_x0 = c[0];
        ms->m_y0 = c[1];
        pp.get("radius", ms->m_radius);
        pp.get("height", ms->m_height);
        if (!(ms->m_radius > 0) || !(ms->m_height > 0)) {
            amrex::Abort("erf.mass_source.radius and height must be positive");
        }
        ms->m_fixed_theta = pp.query("theta", ms->m_theta);
        if (ms->m_fixed_theta && !(ms->m_theta > 0)) {
            amrex::Abort("erf.mass_source.theta must be positive");
        }
        pp.query("tracer", ms->m_tracer);

        if (pp.contains("rate")) {
            amrex::Real r;
            pp.get("rate", r);
            if (!(r >= 0)) { amrex::Abort("erf.mass_source.rate must be zero or positive"); }
            ms->m_t    = {amrex::Real(0)};
            ms->m_mdot = {r};
            ms->m_constant = true;
        } else {
            std::string fname;
            pp.get("rate_file", fname);
            ms->read_rate_file(fname);
        }

        amrex::Print() << "Mass source: centre (" << ms->m_x0 << ", " << ms->m_y0 << ") m, radius "
                       << ms->m_radius << " m, up to " << ms->m_height << " m above the ground; ";
        if (ms->m_constant) {
            amrex::Print() << "constant rate " << ms->m_mdot[0] << " kg/s";
        } else {
            amrex::Print() << ms->m_t.size() << "-point rate table, t = " << ms->m_t.front()
                           << " to " << ms->m_t.back() << " s, total "
                           << ms->total_mass() << " kg";
        }
        if (ms->m_fixed_theta) {
            amrex::Print() << "; theta " << ms->m_theta << " K";
        } else {
            amrex::Print() << "; at the theta of each cell";
        }
        amrex::Print() << (ms->m_tracer ? "; tracer in rhoadv_0" : "") << std::endl;
        return ms;
    }

    [[nodiscard]] amrex::Real x0 () const { return m_x0; }
    [[nodiscard]] amrex::Real y0 () const { return m_y0; }

    /// Mass flow rate [kg/s] at time t (since the start of the run)
    [[nodiscard]] amrex::Real rate (amrex::Real t) const
    {
        if (m_constant) { return m_mdot[0]; }
        if (t < m_t.front() || t > m_t.back()) { return amrex::Real(0); }
        amrex::Long n = 1;
        while (n < m_t.size() - 1 && m_t[n] < t) { ++n; }
        const amrex::Real w = (t - m_t[n-1]) / (m_t[n] - m_t[n-1]);
        return (amrex::Real(1) - w) * m_mdot[n-1] + w * m_mdot[n];
    }

    /// Total mass in the rate table [kg] (trapezoidal)
    [[nodiscard]] amrex::Real total_mass () const
    {
        amrex::Real m = 0;
        for (amrex::Long n = 1; n < m_t.size(); ++n) {
            m += amrex::Real(0.5) * (m_mdot[n] + m_mdot[n-1]) * (m_t[n] - m_t[n-1]);
        }
        return m;
    }

    /**
     * Bring the ground heights and the composite source volume up to date with the grids:
     * a no-op unless the BoxArray or DistributionMapping of a level, or the finest level,
     * has changed since the last call.  Aborts if no level has a cell centre in the region.
     * z_phys_nd and detJ must be allocated on every level.  fill_ground(lev, slab)
     * must fill a one-node-thick nodal slab at the bottom of the domain, on the footprint of
     * every box of the level, with the ground heights (ERF::fill_surface_nodes, which
     * interpolates them from coarser levels where the level does not reach the ground).
     */
    void update_geometry (int finest_level,
                          const amrex::Vector<amrex::Geometry>& geom,
                          const amrex::Vector<amrex::BoxArray>& grids,
                          const amrex::Vector<amrex::DistributionMapping>& dmap,
                          const amrex::Vector<amrex::IntVect>& ref_ratio,
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd,
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& detJ,
                          const amrex::Vector<const amrex::MultiFab*>& mf_mx,
                          const amrex::Vector<const amrex::MultiFab*>& mf_my,
                          const std::function<void(int, amrex::MultiFab&)>& fill_ground)
    {
        bool same = (finest_level + 1 == static_cast<int>(m_lev.size()));
        for (int lev = 0; same && lev <= finest_level; ++lev) {
            same = (m_lev[lev].ba == grids[lev]) && (m_lev[lev].dm == dmap[lev]);
        }
        if (same) { return; }

        m_lev.clear();
        m_lev.resize(finest_level + 1);
        m_volume = 0;
        for (int lev = 0; lev <= finest_level; ++lev) {
            LevelGeometry& lg = m_lev[lev];
            lg.ba = grids[lev];
            lg.dm = dmap[lev];

            // The nodal height of the ground under every column of every box, on the
            // box's own (i,j) range, so that a box that does not reach the bottom (a patch
            // aloft, or a grid split in z) can tell how high its cells are
            const int kg = geom[lev].Domain().smallEnd(2);
            AMREX_ALWAYS_ASSERT(z_phys_nd[lev] && detJ[lev]);
            amrex::BoxList bl = z_phys_nd[lev]->boxArray().boxList();
            for (auto& b : bl) { b.setRange(2, kg); }
            lg.ground = std::make_unique<amrex::MultiFab>(amrex::BoxArray(std::move(bl)),
                                                          z_phys_nd[lev]->DistributionMap(), 1, 0);
            fill_ground(lev, *lg.ground);

            // Cells not covered by the next finer level count towards the composite volume
            std::unique_ptr<amrex::iMultiFab> uncovered;
            if (lev < finest_level) {
                uncovered = std::make_unique<amrex::iMultiFab>(
                    amrex::makeFineMask(grids[lev], dmap[lev], grids[lev+1], ref_ratio[lev], 1, 0));
            }

            const auto dx     = geom[lev].CellSizeArray();
            const auto problo = geom[lev].ProbLoArray();
            const auto period = erf_mass_source::periodic_lengths(geom[lev]);
            const amrex::Real x0 = m_x0, y0 = m_y0, r2 = m_radius * m_radius, hmax = m_height;
            const amrex::Real dV0 = dx[0] * dx[1] * dx[2];

            amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
            amrex::ReduceData<amrex::Real, amrex::Real> reduce_data(reduce_op);
            for (amrex::MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.validbox();
                auto z  = z_phys_nd[lev]->const_array(mfi);
                auto zg = lg.ground->const_array(mfi);
                auto dj = detJ[lev]->const_array(mfi);
                auto mx = mf_mx[lev]->const_array(mfi);
                auto my = mf_my[lev]->const_array(mfi);
                const bool has_mask = (uncovered != nullptr);
                auto msk = has_mask ? uncovered->const_array(mfi) : amrex::Array4<int const>{};
                reduce_op.eval(bx, reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real, amrex::Real>
                    {
                        if (!erf_mass_source::in_source(i,j,k,z,zg,kg,problo,dx,period,x0,y0,r2,hmax)) {
                            return {amrex::Real(0), amrex::Real(0)};
                        }
                        const amrex::Real dV = dV0 * dj(i,j,k) / (mx(i,j,0) * my(i,j,0));
                        const amrex::Real u  = (!has_mask || msk(i,j,k) == 1) ? dV : amrex::Real(0);
                        return {dV, u};
                    });
            }
            auto hv = reduce_data.value(reduce_op);
            amrex::Real v[2] = {amrex::get<0>(hv), amrex::get<1>(hv)};
            amrex::ParallelDescriptor::ReduceRealSum(v, 2);
            lg.volume = v[0];
            m_volume += v[1];
        }

        if (!(m_volume > 0)) {
            amrex::Abort("erf.mass_source: no cell centre lies inside the source region on any level; "
                         "make radius or height larger, or check that center is inside the domain");
        }
        amrex::Print() << "Mass source: region volume " << m_volume << " m^3 on the composite grid";
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::Print() << (lev == 0 ? " (" : ", ") << "level " << lev << ": " << m_lev[lev].volume << " m^3";
        }
        amrex::Print() << ")" << std::endl;
    }

    /**
     * Add the source at time `time` to cc_src of `level` (per unit physical volume, as the
     * other cell-centred sources).  cons is the current conserved state and z_phys_nd the
     * nodal heights of the level.  update_geometry must have been called for the current grids.
     */
    void add_sources (int level,
                      amrex::MultiFab& cc_src,
                      const amrex::MultiFab& cons,
                      const amrex::MultiFab* z_phys_nd,
                      const amrex::Geometry& geom,
                      amrex::Real time) const
    {
        AMREX_ALWAYS_ASSERT(level < static_cast<int>(m_lev.size()) && m_volume > 0);
        const LevelGeometry& lg = m_lev[level];
        const amrex::Real mdot = rate(time);
        if (mdot == amrex::Real(0) || lg.volume == amrex::Real(0)) { return; }

        const auto dx     = geom.CellSizeArray();
        const auto problo = geom.ProbLoArray();
        const int  kg     = geom.Domain().smallEnd(2);
        const auto period = erf_mass_source::periodic_lengths(geom);
        const amrex::Real x0 = m_x0, y0 = m_y0, r2 = m_radius * m_radius, hmax = m_height;
        AMREX_ALWAYS_ASSERT(z_phys_nd != nullptr);

        const amrex::Real S = mdot / m_volume;  // [kg/m^3/s]
        const bool fixed_theta = m_fixed_theta;
        const amrex::Real theta_src = m_theta;
        const bool tracer = m_tracer;

        for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto z   = z_phys_nd->const_array(mfi);
            auto zg  = lg.ground->const_array(mfi);
            auto c   = cons.const_array(mfi);
            auto src = cc_src.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (!erf_mass_source::in_source(i,j,k,z,zg,kg,problo,dx,period,x0,y0,r2,hmax)) { return; }
                const amrex::Real th = fixed_theta ? theta_src : c(i,j,k,RhoTheta_comp) / c(i,j,k,Rho_comp);
                src(i,j,k,Rho_comp)      += S;
                src(i,j,k,RhoTheta_comp) += S * th;
                if (tracer) { src(i,j,k,RhoScalar_comp) += S; }
            });
        }
    }

private:
    void read_rate_file (const std::string& fname)
    {
        amrex::Vector<char> chars;
        amrex::ParallelDescriptor::ReadAndBcastFile(fname, chars);
        std::istringstream is(std::string(chars.dataPtr()));
        std::string line;
        while (std::getline(is, line)) {
            std::string data = line.substr(0, line.find('#'));
            if (data.find_first_not_of(" \t\r") == std::string::npos) { continue; }  // blank or comment
            std::istringstream ls(data);
            amrex::Real t, m;
            std::string extra;
            if (!(ls >> t >> m) || (ls >> extra)) {
                amrex::Abort("erf.mass_source.rate_file " + fname + ": a line is not two numbers, "
                             "time and rate: " + line);
            }
            if (!m_t.empty() && t <= m_t.back()) {
                amrex::Abort("erf.mass_source.rate_file " + fname + ": times must increase strictly");
            }
            if (!(m >= 0)) {
                amrex::Abort("erf.mass_source.rate_file " + fname + ": rates must be zero or positive");
            }
            m_t.push_back(t);
            m_mdot.push_back(m);
        }
        if (m_t.size() < 2) {
            amrex::Abort("erf.mass_source.rate_file " + fname + ": needs at least two rows");
        }
    }

    struct LevelGeometry {
        amrex::BoxArray ba;
        amrex::DistributionMapping dm;
        std::unique_ptr<amrex::MultiFab> ground;  // nodal ground height, one k layer per box
        amrex::Real volume = 0;                   // source volume of this level's own cells
    };

    amrex::Real m_x0 = 0, m_y0 = 0, m_radius = 0, m_height = 0;
    amrex::Real m_theta = 0;
    bool m_fixed_theta = false;
    bool m_tracer = true;
    bool m_constant = false;
    amrex::Vector<amrex::Real> m_t, m_mdot;

    amrex::Vector<LevelGeometry> m_lev;
    amrex::Real m_volume = 0;  // source volume on the composite grid
};

#endif
