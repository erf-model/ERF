// The strength checks of a lattice tower's members (ASCE 10-15) and the steel's reduction with
// temperature (EN 1993-1-2), with the member design file.

#include "ERF_MemberChecks.H"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

#include <AMReX.H>
#include <AMReX_BLassert.H>

namespace erf_towers {

namespace {

constexpr double pi = 3.14159265358979323846;

/** sqrt(29000): ASCE 10-15 writes its angle limits for Fy in ksi with E = 29000 ksi. */
const double sqrt_e_ksi = std::sqrt(29000.0);

std::string lower (std::string s)
{
    for (auto& c : s) { c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }
    return s;
}

bool to_number (const std::string& s, double& v)
{
    if (s.empty()) { return false; }
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

} // namespace

AngleProperties angle_properties (double b, double t)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(t > 0.0 && t < b, "angle_properties: an angle needs 0 < t < b");
    // the heel at the origin, one leg along x ([0, b] x [0, t]), the other along y ([0, t] x [t, b])
    const double a1 = b * t, a2 = t * (b - t);
    const double x1 = 0.5 * b, y1 = 0.5 * t;
    const double x2 = 0.5 * t, y2 = 0.5 * (b + t);
    AngleProperties p;
    p.A = a1 + a2;
    const double cx = (a1 * x1 + a2 * x2) / p.A;
    const double cy = (a1 * y1 + a2 * y2) / p.A;
    const double ixx = b * t * t * t / 12.0 + a1 * (y1 - cy) * (y1 - cy) + t * std::pow(b - t, 3) / 12.0 + a2 * (y2 - cy) * (y2 - cy);
    const double iyy = t * b * b * b / 12.0 + a1 * (x1 - cx) * (x1 - cx) + (b - t) * t * t * t / 12.0 + a2 * (x2 - cx) * (x2 - cx);
    const double ixy = a1 * (x1 - cx) * (y1 - cy) + a2 * (x2 - cx) * (y2 - cy);
    const double mid = 0.5 * (ixx + iyy);
    const double rad = std::sqrt(0.25 * (ixx - iyy) * (ixx - iyy) + ixy * ixy);
    p.Ig = ixx;
    p.Iu = mid + rad;
    p.Iv = mid - rad;
    p.J = (2.0 * b - t) * t * t * t / 3.0;
    p.rv = std::sqrt(p.Iv / p.A);
    p.centroid = cy;
    return p;
}

double effective_slenderness (MemberRole role, EndLoading ends, EndRestraint restraint, double s)
{
    if (role == MemberRole::Leg) { return s; }
    if (s <= 120.0) {
        switch (ends) {
            case EndLoading::Concentric: return s;
            case EndLoading::OneEccentric: return 30.0 + 0.75 * s;
            case EndLoading::BothEccentric: return 60.0 + 0.5 * s;
        }
    }
    switch (restraint) {
        case EndRestraint::None: return s;
        case EndRestraint::OneEnd: return 28.6 + 0.762 * s;
        case EndRestraint::BothEnds: return 46.2 + 0.615 * s;
    }
    return s;
}

double slenderness_limit (MemberRole role)
{
    switch (role) {
        case MemberRole::Leg: return 150.0;
        case MemberRole::Bracing: return 200.0;
        case MemberRole::Redundant: return 250.0;
    }
    return 200.0;
}

double local_buckling_stress (double w_t, double yield, double E)
{
    const double root = std::sqrt(E / yield);
    const double lim1 = 80.0 * root / sqrt_e_ksi;
    const double lim2 = 144.0 * root / sqrt_e_ksi;
    if (w_t <= lim1) { return yield; }
    if (w_t <= lim2) { return (1.677 - 0.677 * w_t / lim1) * yield; }
    return 0.0332 * pi * pi * E / (w_t * w_t);
}

double compression_stress (double effective, double yield, double E)
{
    const double cc = pi * std::sqrt(2.0 * E / yield);
    if (effective <= cc) {
        const double xi = effective / cc;
        return (1.0 - 0.5 * xi * xi) * yield;
    }
    return pi * pi * E / (effective * effective);
}

void steel_reduction (double theta, double& k_y, double& k_E)
{
    // EN 1993-1-2 Table 3.1: temperature (C), k_y (effective yield strength), k_E (slope of the linear elastic range)
    static const double rows[][3] = {{20.0, 1.0, 1.0},       {100.0, 1.0, 1.0},     {200.0, 1.0, 0.9},     {300.0, 1.0, 0.8},
                                     {400.0, 1.0, 0.7},      {500.0, 0.78, 0.6},    {600.0, 0.47, 0.31},   {700.0, 0.23, 0.13},
                                     {800.0, 0.11, 0.09},    {900.0, 0.06, 0.0675}, {1000.0, 0.04, 0.045}, {1100.0, 0.02, 0.0225},
                                     {1200.0, 0.0, 0.0}};
    constexpr int n = static_cast<int>(sizeof(rows) / sizeof(rows[0]));
    if (!(theta > rows[0][0])) { k_y = 1.0; k_E = 1.0; return; }
    if (theta >= rows[n - 1][0]) { k_y = 0.0; k_E = 0.0; return; }
    int i = 0;
    while (theta > rows[i + 1][0]) { ++i; }
    const double s = (theta - rows[i][0]) / (rows[i + 1][0] - rows[i][0]);
    k_y = rows[i][1] + s * (rows[i + 1][1] - rows[i][1]);
    k_E = rows[i][2] + s * (rows[i + 1][2] - rows[i][2]);
}

MemberCheck check_member (const MemberDesign& d, double length, double A, double r_section, double E, double theta,
                          double tension, double compression)
{
    MemberCheck c;
    c.member = d.member;
    c.role = d.role;
    c.length = length;
    double k_y = 1.0, k_E = 1.0;
    steel_reduction(theta, k_y, k_E);
    c.yield = k_y * d.yield;
    c.E = k_E * E;
    double area = A;
    double fy = c.yield;
    if (d.b > 0.0) {
        const AngleProperties ang = angle_properties(d.b, d.t);
        area = ang.A;
        c.r = ang.rv;
        c.w_t = (d.b - d.t) / d.t;
        c.thin = (c.w_t > 25.0);
        if (c.yield > 0.0 && c.E > 0.0) { fy = local_buckling_stress(c.w_t, c.yield, c.E); }
    } else {
        c.r = r_section;
    }
    c.slenderness = length / c.r;
    c.effective = effective_slenderness(d.role, d.ends, d.restraint, c.slenderness);
    c.limit = slenderness_limit(d.role);
    c.slender = ((d.role == MemberRole::Leg) ? c.slenderness : c.effective) > c.limit;
    c.compression_stress = (fy > 0.0 && c.E > 0.0) ? compression_stress(c.effective, fy, c.E) : 0.0;
    const double ft = (d.b > 0.0 && d.one_leg) ? 0.9 * c.yield : c.yield;
    c.tension_capacity = ft * d.net_area * area;
    c.compression_capacity = c.compression_stress * area;
    c.tension = tension;
    c.compression = compression;
    auto ratio = [] (double load, double capacity) {
        if (load <= 0.0) { return 0.0; }
        return capacity > 0.0 ? load / capacity : std::numeric_limits<double>::infinity();
    };
    c.utilisation = std::max(ratio(tension, c.tension_capacity), ratio(compression, c.compression_capacity));
    return c;
}

std::vector<MemberCheck> check_members (const Frame& frame, const std::vector<MemberDesign>& designs,
                                        const std::vector<std::array<double,12>>& f, const std::vector<double>& theta)
{
    const FrameInputs& in = frame.inputs();
    const auto& elems = frame.elements();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(designs.size() == in.members.size(), "check_members: one design per member is needed");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(f.size() == elems.size(), "check_members: one set of end forces per element is needed");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(theta.empty() || theta.size() == in.members.size(), "check_members: one temperature per member");
    const std::size_t nm = in.members.size();
    std::vector<double> length(nm, 0.0), tension(nm, 0.0), compression(nm, 0.0);
    std::vector<const FrameElement*> first(nm, nullptr);
    for (std::size_t e = 0; e < elems.size(); ++e) {
        const std::size_t m = elems[e].member;
        length[m] += elems[e].length;
        if (first[m] == nullptr) { first[m] = &elems[e]; }
        // the axial force at each end, tension positive: the force on the element along its local z at node_b, and against it at node_a
        for (const double n : {-f[e][2], f[e][8]}) {
            tension[m] = std::max(tension[m], n);
            compression[m] = std::max(compression[m], -n);
        }
    }
    std::vector<MemberCheck> out;
    out.reserve(nm);
    for (std::size_t m = 0; m < nm; ++m) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(designs[m].member == in.members[m].id, "check_members: the designs are not in the members' order");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(first[m] != nullptr, "check_members: a member has no element");
        const BeamProperties& p = first[m]->prop;
        const double r = std::sqrt(std::min(p.Ixx, p.Iyy) / p.A);
        // the section's E at 20 C: the element's own is already reduced for the member's temperature
        const FrameSection* s = in.section(in.members[m].section, in.members[m].shape);
        const double t = theta.empty() ? 20.0 : theta[m];
        out.push_back(check_member(designs[m], length[m], p.A, r, s->E, t, tension[m], compression[m]));
    }
    return out;
}

std::string match_designs (const FrameInputs& in, std::vector<MemberDesign>& designs, const std::string& where)
{
    std::map<int, std::size_t> row;
    for (std::size_t i = 0; i < designs.size(); ++i) {
        const MemberDesign& d = designs[i];
        const std::string key = where + ": member " + std::to_string(d.member) + " ";
        if (!row.emplace(d.member, i).second) { return key + "has two rows"; }
        bool found = false;
        for (const auto& m : in.members) { if (m.id == d.member) { found = true; break; } }
        if (!found) { return key + "is not a member of the frame " + in.file; }
        if (!(std::isfinite(d.yield) && d.yield > 0.0)) { return key + "needs a yield strength Fy > 0 (Pa)"; }
        if (!(std::isfinite(d.b) && std::isfinite(d.t) && d.b >= 0.0 && d.t >= 0.0)) { return key + "needs b, t >= 0 (m)"; }
        if ((d.b == 0.0) != (d.t == 0.0) || (d.b > 0.0 && !(d.t < d.b))) {
            return key + "needs an angle with 0 < t < b (m), or b = t = 0 for a member that is not an angle";
        }
        if (!(std::isfinite(d.net_area) && d.net_area > 0.0 && d.net_area <= 1.0)) { return key + "needs a net area ratio in (0, 1]"; }
    }
    std::vector<MemberDesign> ordered;
    ordered.reserve(in.members.size());
    for (const auto& m : in.members) {
        const auto it = row.find(m.id);
        if (it == row.end()) { return where + ": member " + std::to_string(m.id) + " of " + in.file + " has no row"; }
        const MemberDesign& d = designs[it->second];
        if (d.b > 0.0) {
            const FrameSection* s = in.section(m.section, m.shape);
            const double a_section = beam_properties(*s, in.theory).A;
            const double a_angle = angle_properties(d.b, d.t).A;
            if (std::abs(a_angle - a_section) > 0.1 * a_section) {
                std::ostringstream msg;
                msg << where << ": member " << m.id << " is an angle " << d.b << " x " << d.t << " m of area " << a_angle
                    << " m^2, but its property set " << m.section << " in " << in.file << " has area " << a_section
                    << " m^2 (more than 10 % apart)";
                return msg.str();
            }
        }
        ordered.push_back(d);
    }
    designs = ordered;
    return std::string();
}

const char* role_name (MemberRole role)
{
    switch (role) {
        case MemberRole::Leg: return "leg";
        case MemberRole::Bracing: return "bracing";
        case MemberRole::Redundant: return "redundant";
    }
    return "bracing";
}

std::string read_member_designs (const std::string& path, std::vector<MemberDesign>& designs)
{
    std::string text;
    if (!read_text_file(path, text)) { return "cannot read the member design file '" + path + "'"; }
    std::istringstream f(text);
    designs.clear();
    std::string line;
    std::size_t ln = 0;
    while (std::getline(f, line)) {
        ++ln;
        if (!line.empty() && line.back() == '\r') { line.pop_back(); }
        const auto p = line.find_first_not_of(" \t");
        if (p == std::string::npos || line[p] == '#' || line[p] == '!') { continue; }
        std::istringstream is(line);
        std::vector<std::string> t;
        std::string w;
        while (is >> w) { t.push_back(w); }
        const std::string at = path + " line " + std::to_string(ln) + ": ";
        if (t.size() != 9) {
            return at + "a row has 9 values: MemberID Role Fy(Pa) b(m) t(m) NetArea Bolted Ends Restraint";
        }
        MemberDesign d;
        double id = 0.0;
        if (!to_number(t[0], id) || id != std::floor(id) || std::abs(id) > 2.0e9) { return at + "the member id must be an integer"; }
        d.member = static_cast<int>(id);
        const std::string role = lower(t[1]);
        if (role == "leg") { d.role = MemberRole::Leg; }
        else if (role == "bracing") { d.role = MemberRole::Bracing; }
        else if (role == "redundant") { d.role = MemberRole::Redundant; }
        else { return at + "the role must be leg, bracing or redundant, not '" + t[1] + "'"; }
        if (!to_number(t[2], d.yield) || !to_number(t[3], d.b) || !to_number(t[4], d.t) || !to_number(t[5], d.net_area)) {
            return at + "Fy, b, t and NetArea must be numbers";
        }
        const std::string bolted = lower(t[6]), ends = lower(t[7]), restraint = lower(t[8]);
        if (bolted == "one") { d.one_leg = true; }
        else if (bolted == "both") { d.one_leg = false; }
        else { return at + "Bolted must be one or both (the legs of the angle bolted at its ends), not '" + t[6] + "'"; }
        if (ends == "concentric") { d.ends = EndLoading::Concentric; }
        else if (ends == "one") { d.ends = EndLoading::OneEccentric; }
        else if (ends == "both") { d.ends = EndLoading::BothEccentric; }
        else { return at + "Ends must be concentric, one or both (the ends with a normal framing eccentricity), not '" + t[7] + "'"; }
        if (restraint == "none") { d.restraint = EndRestraint::None; }
        else if (restraint == "one") { d.restraint = EndRestraint::OneEnd; }
        else if (restraint == "both") { d.restraint = EndRestraint::BothEnds; }
        else { return at + "Restraint must be none, one or both (the ends partially restrained against rotation), not '" + t[8] + "'"; }
        designs.push_back(d);
    }
    if (designs.empty()) { return path + ": no member rows"; }
    return std::string();
}

bool write_member_designs (const std::string& path, const std::vector<MemberDesign>& designs, const std::string& title)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out) { return false; }
    out << "# " << title << "\n"
        << "# MemberID  Role  Fy(Pa)  b(m)  t(m)  NetArea(-)  Bolted  Ends  Restraint\n";
    out << std::setprecision(10);
    for (const auto& d : designs) {
        const char* ends = (d.ends == EndLoading::Concentric) ? "concentric" : (d.ends == EndLoading::OneEccentric) ? "one" : "both";
        const char* restraint = (d.restraint == EndRestraint::None) ? "none" : (d.restraint == EndRestraint::OneEnd) ? "one" : "both";
        out << d.member << " " << role_name(d.role) << " " << d.yield << " " << d.b << " " << d.t << " " << d.net_area << " "
            << (d.one_leg ? "one" : "both") << " " << ends << " " << restraint << "\n";
    }
    return static_cast<bool>(out);
}

} // namespace erf_towers
