#include <ERF_InflowProfile.H>
#include <ERF_Constants.H>

#include <AMReX.H>
#include <AMReX_Math.H>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>

using namespace amrex;

void
InflowProfile::copy_to_device ()
{
    auto copy = [] (const Vector<Real>& h, Gpu::DeviceVector<Real>& d) {
        d.resize(h.size());
        Gpu::copy(Gpu::hostToDevice, h.begin(), h.end(), d.begin());
    };
    copy(z, z_d);
    copy(u, u_d);
    copy(v, v_d);
    copy(theta, theta_d);
    copy(tke, tke_d);
}

namespace {

enum InflowColumn { col_z, col_u, col_v, col_w, col_theta, col_tke };

std::string
at_line (const int lineno)
{
    return "line " + std::to_string(lineno) + ": ";
}

} // namespace

/**
 * The layout of kynema-sgf's TabulatedProfile: one row per height, heights
 * strictly increasing. A comment line whose first word is "z" names the
 * columns (z u v w T tke, in any order; theta and k are accepted for T and
 * tke); other comment lines are ignored. Without such a header, four columns
 * are z u v T and five are z u v T tke. A w column is accepted and not used:
 * the vertical velocity on the face keeps its boundary value. Heights are
 * above the local ground.
 */
std::string
parse_inflow_profile (std::istream& is, InflowProfile& prof)
{
    std::vector<int> cols;
    std::vector<std::vector<Real>> rows;
    std::vector<int> row_line;

    std::string line;
    int lineno = 0;
    while (std::getline(is, line)) {
        ++lineno;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos) { continue; }

        if (line[first] == '#') {
            std::istringstream ws(line.substr(first + 1));
            std::vector<std::string> words;
            for (std::string w; ws >> w; ) {
                std::transform(w.begin(), w.end(), w.begin(),
                               [] (unsigned char c) { return static_cast<char>(std::tolower(c)); });
                words.push_back(w);
            }
            if (words.empty() || words[0] != "z") { continue; }
            if (!cols.empty() || !rows.empty()) {
                return at_line(lineno) + "the column header must come once, before the data";
            }
            for (const auto& w : words) {
                int c;
                if      (w == "z")                  { c = col_z; }
                else if (w == "u")                  { c = col_u; }
                else if (w == "v")                  { c = col_v; }
                else if (w == "w")                  { c = col_w; }
                else if (w == "t" || w == "theta")  { c = col_theta; }
                else if (w == "tke" || w == "k")    { c = col_tke; }
                else {
                    return at_line(lineno) + "unknown column '" + w + "' (accepted: z u v w T tke)";
                }
                if (std::find(cols.begin(), cols.end(), c) != cols.end()) {
                    return at_line(lineno) + "column '" + w + "' appears twice";
                }
                cols.push_back(c);
            }
            continue;
        }

        std::istringstream ls(line);
        std::vector<Real> vals;
        for (std::string tok; ls >> tok; ) {
            char* end = nullptr;
            const double d = std::strtod(tok.c_str(), &end);
            if (end == tok.c_str() || *end != '\0' || !std::isfinite(d)) {
                return at_line(lineno) + "'" + tok + "' is not a number";
            }
            vals.push_back(static_cast<Real>(d));
        }
        rows.push_back(vals);
        row_line.push_back(lineno);
    }

    if (rows.empty()) { return "no data rows"; }
    if (cols.empty()) {
        const auto nc = rows[0].size();
        if (nc == 4) {
            cols = {col_z, col_u, col_v, col_theta};
        } else if (nc == 5) {
            cols = {col_z, col_u, col_v, col_theta, col_tke};
        } else {
            return "a table without a '# z ...' header needs 4 columns (z u v T) or 5 (z u v T tke), not "
                   + std::to_string(nc);
        }
    }
    auto has = [&cols] (int c) { return std::find(cols.begin(), cols.end(), c) != cols.end(); };
    if (!has(col_z) || !has(col_u) || !has(col_v)) {
        return "the columns must include z, u and v";
    }
    if (rows.size() < 2) { return "at least two heights are needed"; }

    Vector<Real> z, u, v, theta, tke;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        if (rows[r].size() != cols.size()) {
            return at_line(row_line[r]) + "expected " + std::to_string(cols.size()) +
                   " columns, found " + std::to_string(rows[r].size());
        }
        for (std::size_t c = 0; c < cols.size(); ++c) {
            const Real x = rows[r][c];
            switch (cols[c]) {
                case col_z:     z.push_back(x);     break;
                case col_u:     u.push_back(x);     break;
                case col_v:     v.push_back(x);     break;
                case col_theta: theta.push_back(x); break;
                case col_tke:   tke.push_back(x);   break;
                default:                            break;
            }
        }
    }
    for (Long r = 1; r < z.size(); ++r) {
        if (!(z[r] > z[r-1])) {
            return at_line(row_line[r]) + "heights must increase strictly";
        }
    }
    for (Long r = 0; r < theta.size(); ++r) {
        if (!(theta[r] > Real(0.0))) { return at_line(row_line[r]) + "T must be positive"; }
    }
    for (Long r = 0; r < tke.size(); ++r) {
        if (tke[r] < Real(0.0)) { return at_line(row_line[r]) + "tke must not be negative"; }
    }

    prof.z = z;
    prof.u = u;
    prof.v = v;
    prof.theta = theta;
    prof.tke = tke;
    prof.has_theta = !theta.empty();
    prof.has_tke   = !tke.empty();
    prof.active    = true;
    return "";
}

void
read_inflow_profile_file (const std::string& fname, InflowProfile& prof)
{
    std::ifstream is(fname);
    if (!is.is_open()) {
        Abort("inflow_profile_file '" + fname + "' could not be opened");
    }
    const std::string err = parse_inflow_profile(is, prof);
    if (!err.empty()) {
        Abort("inflow_profile_file '" + fname + "': " + err);
    }
    prof.source = fname;
}

/**
 * u = (u_star / kappa) ln((z + z0)/z0), so the wind is zero at the ground and
 * `speed` at `height`, optionally capped at `max_speed`, blowing from
 * `direction`. tke = u_star^2 / Cmu0^2 at the ground, the value erf.dirichlet_k holds
 * in the first cell, tapering linearly to zero at u_star tke_zscale with the floor
 * of SurfaceLayer::init_tke_from_ustar. The table has the ground and 399
 * heights spaced geometrically from 1 cm to ztop.
 */
std::string
make_log_law_inflow_profile (const InflowLogLaw& p, const Real ztop, InflowProfile& prof)
{
    if (!(p.speed > Real(0.0)))      { return "inflow_log_law.speed must be positive"; }
    if (!(p.height > Real(0.0)))     { return "inflow_log_law.height must be positive"; }
    if (!(p.z0 > Real(0.0)))         { return "the roughness length (inflow_log_law.z0 or erf.most.z0) must be positive"; }
    if (!(p.tke_zscale > Real(0.0))) { return "inflow_log_law.tke_zscale must be positive"; }
    if (!(p.Cmu0 > Real(0.0)))       { return "Cmu0 must be positive"; }
    if (!(ztop > Real(0.0)))         { return "the domain height must be positive"; }

    const Real ustar = KAPPA * p.speed / std::log((p.height + p.z0) / p.z0);
    const Real dir   = p.direction * Math::pi<Real>() / Real(180.0);
    const Real ex    = -std::sin(dir);
    const Real ey    = -std::cos(dir);
    const Real tke0  = ustar * ustar / (p.Cmu0 * p.Cmu0);
    const Real small = Real(0.01);

    const int  n     = 400;
    const Real zmin  = std::min(Real(0.01), Real(0.5) * ztop);
    const Real ratio = std::pow(ztop / zmin, Real(1.0) / static_cast<Real>(n - 2));

    Vector<Real> z(n), u(n), v(n), tke(n);
    for (int m = 0; m < n; ++m) {
        const Real zm = (m == 0) ? Real(0.0) : zmin * std::pow(ratio, static_cast<Real>(m - 1));
        Real s = ustar / KAPPA * std::log((zm + p.z0) / p.z0);
        if (p.max_speed > Real(0.0)) { s = std::min(s, p.max_speed); }
        z[m]   = zm;
        u[m]   = s * ex;
        v[m]   = s * ey;
        tke[m] = tke0 * std::max((ustar * p.tke_zscale - zm) / (std::max(ustar, small) * p.tke_zscale), small);
    }

    prof.z = z;
    prof.u = u;
    prof.v = v;
    prof.theta.clear();
    prof.tke = tke;
    prof.has_theta = false;
    prof.has_tke   = true;
    prof.active    = true;
    prof.source    = "log_law";
    return "";
}
