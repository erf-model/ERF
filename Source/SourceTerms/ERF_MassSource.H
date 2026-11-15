#ifndef ERF_MASS_SOURCE_H_
#define ERF_MASS_SOURCE_H_

#include <memory>
#include <sstream>
#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
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
 *     so the injection by itself does not change the temperature, and
 *   - the first passive scalar (rho*s, plotted as rhoadv_0) when erf.mass_source.tracer
 *     is true, so that s is the mass fraction of the released gas.
 * No momentum is added: the gas enters at rest.  The gas has the properties of dry air
 * (in a moist run it carries no water, so it dilutes the mixing ratios of the cells it
 * enters).  The cell volumes include the terrain metric and the map factors.
 *
 * Each level normalizes the source by the volume of its own cells in the region, so the
 * composite rate is mdot when the region lies entirely inside the finest level that covers
 * any of it (or outside every refined level).
 *
 * Inputs:
 *   erf.mass_source.center    = x y          [m]
 *   erf.mass_source.radius    = r            [m]
 *   erf.mass_source.height    = h            [m above the local ground]
 *   erf.mass_source.rate      = mdot         [kg/s]  constant, or
 *   erf.mass_source.rate_file = file         two columns "t mdot" [s, kg/s], '#' comments;
 *                                            linear in between, zero outside the table
 *   erf.mass_source.theta     = theta        [K]     optional
 *   erf.mass_source.tracer    = true         optional
 */
namespace erf_mass_source {

/// Is the centre of cell (i,j,k) inside the source region?  z is the nodal height
/// z_phys_nd, or empty on a mesh without terrain (then the height of the centre above
/// the bottom is (k+1/2) dz)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
in_source (int i, int j, int k,
           amrex::Array4<amrex::Real const> const& z, bool has_z,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           amrex::Real x0, amrex::Real y0, amrex::Real r2, amrex::Real hmax) noexcept
{
    const amrex::Real xc = problo[0] + (i + amrex::Real(0.5)) * dx[0];
    const amrex::Real yc = problo[1] + (j + amrex::Real(0.5)) * dx[1];
    if ((xc - x0) * (xc - x0) + (yc - y0) * (yc - y0) > r2) { return false; }
    amrex::Real hag;
    if (has_z) {
        const amrex::Real zg = amrex::Real(0.25) * (z(i,j,0) + z(i+1,j,0) + z(i,j+1,0) + z(i+1,j+1,0));
        const amrex::Real zc = amrex::Real(0.125) * (z(i,j,k  ) + z(i+1,j,k  ) + z(i,j+1,k  ) + z(i+1,j+1,k  ) +
                                                     z(i,j,k+1) + z(i+1,j,k+1) + z(i,j+1,k+1) + z(i+1,j+1,k+1));
        hag = zc - zg;
    } else {
        hag = (k + amrex::Real(0.5)) * dx[2];
    }
    return hag < hmax;
}

} // namespace erf_mass_source

class MassSource
{
public:
    /// The source the inputs ask for, or nullptr when neither rate nor rate_file is given
    static std::unique_ptr<MassSource> from_inputs ()
    {
        amrex::ParmParse pp("erf.mass_source");
        if (!pp.contains("rate") && !pp.contains("rate_file")) { return nullptr; }
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
        if (ms->m_radius <= 0 || ms->m_height <= 0) {
            amrex::Abort("erf.mass_source.radius and height must be positive");
        }
        ms->m_fixed_theta = pp.query("theta", ms->m_theta);
        pp.query("tracer", ms->m_tracer);

        if (pp.contains("rate")) {
            amrex::Real r;
            pp.get("rate", r);
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
        amrex::Print() << (ms->m_tracer ? "; tracer in rhoadv_0" : "") << std::endl;
        return ms;
    }

    /// Mass flow rate [kg/s] at time t
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
     * Add the source at time `time` to cc_src (per unit physical volume, as the other
     * cell-centred sources).  cons is the current conserved state; z_phys_nd and detJ
     * may be null on a mesh without terrain.  mf_mx and mf_my are the cell-centred map
     * factors.
     */
    void add_sources (amrex::MultiFab& cc_src,
                      const amrex::MultiFab& cons,
                      const amrex::MultiFab* z_phys_nd,
                      const amrex::MultiFab* detJ,
                      const amrex::MultiFab& mf_mx,
                      const amrex::MultiFab& mf_my,
                      const amrex::Geometry& geom,
                      amrex::Real time) const
    {
        const amrex::Real mdot = rate(time);
        if (mdot == amrex::Real(0)) { return; }

        const auto dx     = geom.CellSizeArray();
        const auto problo = geom.ProbLoArray();
        const amrex::Real x0 = m_x0, y0 = m_y0, r2 = m_radius * m_radius, hmax = m_height;
        const bool has_z = (z_phys_nd != nullptr);
        const bool has_J = (detJ != nullptr);

        // Volume of the source region
        const amrex::Real dV0 = dx[0] * dx[1] * dx[2];
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto z  = has_z ? z_phys_nd->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            auto dj = has_J ? detJ->const_array(mfi)      : amrex::Array4<amrex::Real const>{};
            auto mx = mf_mx.const_array(mfi);
            auto my = mf_my.const_array(mfi);
            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real>
                {
                    if (!erf_mass_source::in_source(i,j,k,z,has_z,problo,dx,x0,y0,r2,hmax)) {
                        return {amrex::Real(0)};
                    }
                    const amrex::Real J = has_J ? dj(i,j,k) : amrex::Real(1);
                    return { dV0 * J / (mx(i,j,0) * my(i,j,0)) };
                });
        }
        amrex::Real vol = amrex::get<0>(reduce_data.value());
        amrex::ParallelDescriptor::ReduceRealSum(vol);
        if (vol <= amrex::Real(0)) {
            amrex::Abort("erf.mass_source: no cell centre lies inside the source region on this level; "
                         "increase radius or height");
        }

        const amrex::Real S = mdot / vol;  // [kg/m^3/s]
        const bool fixed_theta = m_fixed_theta;
        const amrex::Real theta_src = m_theta;
        const bool tracer = m_tracer;

        for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            auto z   = has_z ? z_phys_nd->const_array(mfi) : amrex::Array4<amrex::Real const>{};
            auto c   = cons.const_array(mfi);
            auto src = cc_src.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (!erf_mass_source::in_source(i,j,k,z,has_z,problo,dx,x0,y0,r2,hmax)) { return; }
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
            const auto p = line.find('#');
            if (p != std::string::npos) { line.erase(p); }
            std::istringstream ls(line);
            amrex::Real t, m;
            if (!(ls >> t)) { continue; }
            if (!(ls >> m)) { amrex::Abort("erf.mass_source.rate_file " + fname + ": line without a rate: " + line); }
            if (!m_t.empty() && t <= m_t.back()) {
                amrex::Abort("erf.mass_source.rate_file " + fname + ": times must increase strictly");
            }
            m_t.push_back(t);
            m_mdot.push_back(m);
        }
        if (m_t.size() < 2) {
            amrex::Abort("erf.mass_source.rate_file " + fname + ": needs at least two rows");
        }
    }

    amrex::Real m_x0 = 0, m_y0 = 0, m_radius = 0, m_height = 0;
    amrex::Real m_theta = 0;
    bool m_fixed_theta = false;
    bool m_tracer = true;
    bool m_constant = false;
    amrex::Vector<amrex::Real> m_t, m_mdot;
};

#endif
