// WakeLines: sampling points behind a rotor, the running average, its files and checkpoint.

#include "ERF_WakeLines.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "ERF_DiagnosticsLog.H"

using namespace amrex;

namespace erf_actuator {

WakeLines::WakeLines (std::string name, std::string output_root,
                      const std::array<Real,3>& hub, const std::array<Real,3>& axis,
                      Real diameter, const std::vector<Real>& lines_xD,
                      Real half_width, int num_points, Real z_min)
    : WakeLines(std::move(name), std::move(output_root), hub, axis, diameter, lines_xD, half_width, num_points,
                std::vector<Real>(lines_xD.size(), z_min))
{}

WakeLines::WakeLines (std::string name, std::string output_root,
                      const std::array<Real,3>& hub, const std::array<Real,3>& axis,
                      Real diameter, const std::vector<Real>& lines_xD,
                      Real half_width, int num_points, const std::vector<Real>& z_min)
    : m_name(std::move(name)), m_output_root(std::move(output_root)),
      m_hub(hub), m_axis(axis), m_diameter(diameter), m_xD(lines_xD),
      m_half_width(half_width), m_npts(num_points)
{
    if (!(m_diameter > 0.0)) { Abort("WakeLines " + m_name + ": the diameter " + std::to_string(m_diameter) + " m must be positive"); }
    if (m_npts < 2) { Abort("WakeLines " + m_name + ": num_points = " + std::to_string(m_npts) + "; at least two points per line"); }
    if (!(m_half_width > 0.0)) { Abort("WakeLines " + m_name + ": half_width = " + std::to_string(m_half_width) + " diameters must be positive"); }
    if (m_xD.empty()) { Abort("WakeLines " + m_name + ": no downstream distances (x/D) are given"); }
    if (z_min.size() != m_xD.size()) {
        Abort("WakeLines " + m_name + ": " + std::to_string(z_min.size()) + " ground heights for " +
              std::to_string(m_xD.size()) + " downstream distances; one per distance");
    }
    for (const Real zm : z_min) {
        if (!std::isfinite(zm)) { Abort("WakeLines " + m_name + ": the ground height z_min must be finite (m)"); }
    }
    // the downstream direction is the rotor axis projected on the horizontal plane: a shaft
    // tilt of a few degrees would otherwise put the far lines into the ground or the sky
    m_axis[2] = 0.0;
    const Real la = std::sqrt(m_axis[0]*m_axis[0] + m_axis[1]*m_axis[1]);
    if (!(la > 1.0e-6)) { Abort("WakeLines " + m_name + ": the axis (" + std::to_string(axis[0]) + ", " + std::to_string(axis[1]) + ", " + std::to_string(axis[2]) + ") must not be vertical"); }
    for (auto& a : m_axis) { a /= la; }
    // the lateral direction: horizontal and normal to the downstream direction
    m_lateral = {{-m_axis[1], m_axis[0], 0.0}};

    const std::array<Real,3> vertical{{0.0, 0.0, 1.0}};
    for (std::size_t n = 0; n < m_xD.size(); ++n) {
        const Real xD = m_xD[n];
        const Real zmin_n = z_min[n];
        std::array<Real,3> c;
        for (int d = 0; d < 3; ++d) { c[d] = m_hub[d] + xD * m_diameter * m_axis[d]; }
        for (int line = 0; line < 2; ++line) {
            const auto& e = (line == 0) ? m_lateral : vertical;
            // the vertical line stops at the ground; the lateral one is symmetric
            Real s_lo = -m_half_width;
            if (line == 1) { s_lo = std::max(s_lo, (zmin_n - c[2]) / m_diameter); }
            const Real s_hi = m_half_width;
            if (!(s_lo < s_hi)) { Abort("WakeLines " + m_name + ": the line centre at z = " + std::to_string(c[2]) + " m sits at or below the ground z_min = " + std::to_string(zmin_n) + " m"); }
            for (int i = 0; i < m_npts; ++i) {
                const Real s = s_lo + (s_hi - s_lo) * i / (m_npts - 1);
                m_s.push_back(s);
                for (int d = 0; d < 3; ++d) { m_pos.push_back(c[d] + s * m_diameter * e[d]); }
            }
        }
    }
    m_sum.assign(m_pos.size(), 0.0);
}

void
WakeLines::accumulate (const std::vector<Real>& vel)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.size() == m_pos.size(), "WakeLines::accumulate: one velocity per point");
    for (std::size_t i = 0; i < vel.size(); ++i) {
        if (!std::isfinite(vel[i])) {
            Abort("WakeLines " + m_name + ": the sampled velocity at point " + std::to_string(i / 3) + " (0-based) is not finite");
        }
    }
    for (std::size_t i = 0; i < vel.size(); ++i) { m_sum[i] += static_cast<double>(vel[i]); }
    ++m_count;
}

std::vector<Real>
WakeLines::average () const
{
    std::vector<Real> avg(m_pos.size(), 0.0);
    if (m_count > 0) {
        for (std::size_t i = 0; i < avg.size(); ++i) { avg[i] = static_cast<Real>(m_sum[i] / m_count); }
    }
    return avg;
}

// one row per point: [prefix,]xD,line,s,x,y,z,u,v,w
void
WakeLines::write_rows (std::ofstream& out, const std::vector<Real>& vel, const std::string& prefix) const
{
    int p = 0;
    for (const Real xD : m_xD) {
        for (const char* line : {"lateral", "vertical"}) {
            for (int i = 0; i < m_npts; ++i, ++p) {
                out << prefix << xD << "," << line << "," << m_s[p] << ","
                    << m_pos[3*p] << "," << m_pos[3*p+1] << "," << m_pos[3*p+2] << ","
                    << vel[3*p] << "," << vel[3*p+1] << "," << vel[3*p+2] << "\n";
            }
        }
    }
}

void
WakeLines::write_instantaneous (double time, const std::vector<Real>& vel, bool truncate) const
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vel.size() == m_pos.size(), "WakeLines::write_instantaneous: one velocity per point");
    std::ofstream out;
    if (open_log(out, m_output_root + "_wake.csv", truncate)) {
        out << "time,xD,line,s,x,y,z,u,v,w\n";
    }
    out << std::setprecision(10);
    std::ostringstream prefix;
    prefix << std::setprecision(10) << time << ",";
    write_rows(out, vel, prefix.str());
}

void
WakeLines::write_average () const
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream out;
    open_log(out, m_output_root + "_wake_avg.csv", true);
    out << "samples,xD,line,s,x,y,z,u,v,w\n" << std::setprecision(10);
    write_rows(out, average(), std::to_string(m_count) + ",");
}

void
WakeLines::write_state (const std::string& dir) const
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream out(dir + "/" + m_name + "_wake_avg.dat", std::ios::trunc);
    if (!out) { Abort("cannot write the wake-average checkpoint for " + m_name + " in '" + dir + "'"); }
    out << std::setprecision(17) << "count = " << m_count << "\nsize = " << m_sum.size() << "\n"
        << "hub = " << m_hub[0] << " " << m_hub[1] << " " << m_hub[2] << "\n"
        << "axis = " << m_axis[0] << " " << m_axis[1] << " " << m_axis[2] << "\n"
        << "diameter = " << m_diameter << "\n";
    for (const Real s : m_s) { out << s << "\n"; }
    for (const Real x : m_pos) { out << x << "\n"; }
    for (const double s : m_sum) { out << s << "\n"; }
}

bool
WakeLines::read_state (const std::string& dir)
{
    const std::string fname = dir + "/" + m_name + "_wake_avg.dat";
    // decided on the I/O rank for every rank: the read below is collective
    if (!file_exists_everywhere(fname)) {
        Print() << "WakeLines " << m_name << ": no " << fname << " in the checkpoint; the wake average starts afresh\n";
        return false;
    }
    Vector<char> chars;
    ParallelDescriptor::ReadAndBcastFile(fname, chars);
    std::istringstream in(std::string(chars.dataPtr(), chars.size()));
    std::string key, eq;
    std::size_t size = 0;
    int count = 0;
    std::array<Real,3> hub, axis;
    Real diameter = 0.0;
    if (!(in >> key >> eq >> count) || key != "count" || !(in >> key >> eq >> size) || key != "size" ||
        !(in >> key >> eq >> hub[0] >> hub[1] >> hub[2]) || key != "hub" ||
        !(in >> key >> eq >> axis[0] >> axis[1] >> axis[2]) || key != "axis" ||
        !(in >> key >> eq >> diameter) || key != "diameter") {
        Abort("malformed wake-average checkpoint '" + fname + "'");
    }
    if (size != m_sum.size()) {
        Abort("the wake-average checkpoint '" + fname + "' holds " + std::to_string(size) + " values but the lines of " +
              m_name + " have " + std::to_string(m_sum.size()) + "; the wake inputs must match the run being restarted");
    }
    const Real naxis = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    const bool finite = std::isfinite(hub[0]) && std::isfinite(hub[1]) && std::isfinite(hub[2]) && std::isfinite(naxis);
    if (count < 0 || !finite || !(std::isfinite(diameter) && diameter > 0.0) ||
        !(std::abs(naxis - Real(1.0)) < Real(1.0e-4) && std::abs(axis[2]) < Real(1.0e-4))) {
        Abort("malformed wake-average checkpoint '" + fname + "': the count must be >= 0, the diameter positive, "
              "the axis a horizontal unit vector and every value finite");
    }
    // the geometry is taken over as it was: a body model may rebuild the hub and diameter from
    // lower-precision data and differ at the last digits from one build of the lines to another
    std::vector<Real> s(m_s.size()), pos(size);
    std::vector<double> sum(size);
    auto read_all = [&] (auto& v) {
        for (auto& x : v) {
            if (!(in >> x)) { Abort("truncated wake-average checkpoint '" + fname + "'"); }
            if (!std::isfinite(x)) { Abort("malformed wake-average checkpoint '" + fname + "': a value is not finite"); }
        }
    };
    read_all(s);
    read_all(pos);
    read_all(sum);
    m_s = s;
    m_pos = pos;
    m_sum = sum;
    m_hub = hub;
    m_axis = axis;
    m_diameter = diameter;
    m_count = count;
    return true;
}

} // namespace erf_actuator
