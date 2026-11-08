// MoorDynSystem: MoorDyn-C v2 C API calls with ERF's error handling.

#include "ERF_MoorDynSystem.H"

#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>

#include <AMReX.H>

namespace erf_moordyn {

namespace {
std::string lower (std::string s)
{
    for (auto& c : s) { c = static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }
    return s;
}

// whether every value is finite; the first offending index in `bad`
bool all_finite (const std::vector<double>& v, std::size_t& bad)
{
    for (std::size_t i = 0; i < v.size(); ++i) { if (!std::isfinite(v[i])) { bad = i; return false; } }
    return true;
}
} // namespace

std::string library_version ()
{
#ifdef ERF_MOORDYN_USE_STUB
    return "stub";
#else
    return ERF_MOORDYN_VERSION;
#endif
}

bool is_stub ()
{
#ifdef ERF_MOORDYN_USE_STUB
    return true;
#else
    return false;
#endif
}

std::string error_name (int code)
{
    switch (code) {
    case MOORDYN_SUCCESS:             return "MOORDYN_SUCCESS";
    case MOORDYN_INVALID_INPUT_FILE:  return "MOORDYN_INVALID_INPUT_FILE";
    case MOORDYN_INVALID_OUTPUT_FILE: return "MOORDYN_INVALID_OUTPUT_FILE";
    case MOORDYN_INVALID_INPUT:       return "MOORDYN_INVALID_INPUT";
    case MOORDYN_NAN_ERROR:           return "MOORDYN_NAN_ERROR";
    case MOORDYN_MEM_ERROR:           return "MOORDYN_MEM_ERROR";
    case MOORDYN_INVALID_VALUE:       return "MOORDYN_INVALID_VALUE";
    case MOORDYN_NON_IMPLEMENTED:     return "MOORDYN_NON_IMPLEMENTED";
    case MOORDYN_UNHANDLED_ERROR:     return "MOORDYN_UNHANDLED_ERROR";
    default:                          return "MoorDyn error " + std::to_string(code);
    }
}

FPTrapMask::FPTrapMask ()
{
    std::fegetexceptflag(&m_flags, FE_ALL_EXCEPT);
    m_traps = amrex::disableFPExcept(amrex::FPExcept::all);
}

FPTrapMask::~FPTrapMask ()
{
    // a flag MoorDyn raised under the mask would be reported by the next instruction once its trap is back on
    std::feclearexcept(FE_ALL_EXCEPT);
    std::fesetexceptflag(&m_flags, FE_ALL_EXCEPT);
    amrex::setFPExcept(m_traps);
}

std::string check_wave_kinematics_option (const std::string& input_file)
{
    std::ifstream in(input_file);
    if (!in) { return "cannot read the MoorDyn input file '" + input_file + "'"; }
    // the OPTIONS block runs from its "--- OPTIONS ---" header to the next "---" line; each option
    // line is "value name [description]"
    std::string row;
    bool options = false;
    while (std::getline(in, row)) {
        if (row.find("---") != std::string::npos) {
            options = (lower(row).find("option") != std::string::npos);
            continue;
        }
        if (!options) { continue; }
        std::istringstream tok(row);
        std::string value, name;
        if (!(tok >> value >> name) || lower(name) != "wavekin") { continue; }
        if (value == "1") { return std::string(); }
        return "the MoorDyn input file '" + input_file + "' sets WaveKin = " + value +
               "; ERF passes the wind through MoorDyn's external kinematics, which need WaveKin = 1 in its OPTIONS";
    }
    return "the MoorDyn input file '" + input_file + "' does not set WaveKin in its OPTIONS; ERF passes the wind through "
           "MoorDyn's external kinematics, which need WaveKin = 1";
}

std::unique_ptr<MoorDynSystem>
MoorDynSystem::create (const std::string& input_file, const std::string& log_file, int log_level, std::string& err)
{
    err.clear();
    MoorDyn sys = MoorDyn_Create(input_file.c_str());
    if (sys == nullptr) {
        err = "MoorDyn (" + library_version() + ") could not create a system from '" + input_file +
              "': the file is missing or malformed (MoorDyn's message is on stderr)";
        return nullptr;
    }
    std::unique_ptr<MoorDynSystem> s(new MoorDynSystem(sys, input_file));
    s->check(MoorDyn_SetVerbosity(sys, log_level), "MoorDyn_SetVerbosity");
    if (!log_file.empty()) {
        s->check(MoorDyn_SetLogFile(sys, log_file.c_str()), "MoorDyn_SetLogFile");
        s->check(MoorDyn_SetLogLevel(sys, log_level), "MoorDyn_SetLogLevel");
    }
    // MoorDyn-C reports kinematics points whatever WaveKin is, so the option is read from the file
    err = check_wave_kinematics_option(input_file);
    if (!err.empty()) { return nullptr; }
    return s;
}

MoorDynSystem::MoorDynSystem (MoorDyn sys, std::string file)
    : m_sys(sys), m_file(std::move(file)) {}

MoorDynSystem::~MoorDynSystem ()
{
    if (m_sys != nullptr) { MoorDyn_Close(m_sys); }
}

void MoorDynSystem::check (int rc, const std::string& what) const
{
    if (rc != MOORDYN_SUCCESS) {
        const std::string hint = (rc == MOORDYN_NAN_ERROR)
            ? ": the line integration diverged; reduce MoorDyn's internal step (erf.conductors.moordyn_cfl or "
              "erf.conductors.moordyn_dt for a conductor line; CFL and dtM in the input file)" : "";
        amrex::Abort(what + " failed with " + error_name(rc) + " for the MoorDyn system from '" + m_file + "'" + hint);
    }
}

unsigned MoorDynSystem::num_coupled_dof () const
{
    unsigned n = 0;
    check(MoorDyn_NCoupledDOF(m_sys, &n), "MoorDyn_NCoupledDOF");
    return n;
}

std::string MoorDynSystem::init (const std::vector<double>& x, const std::vector<double>& xd, bool compute_ic)
{
    const unsigned ndof = num_coupled_dof();
    if (x.size() != ndof || xd.size() != ndof) {
        return "MoorDyn system from '" + m_file + "' has " + std::to_string(ndof) + " coupled degrees of freedom but " +
               std::to_string(x.size()) + " positions and " + std::to_string(xd.size()) + " velocities were given";
    }
    std::size_t bad = 0;
    if (!all_finite(x, bad) || !all_finite(xd, bad)) {
        return "MoorDyn system from '" + m_file + "': coupled position or velocity component " + std::to_string(bad) +
               " (0-based) is not finite";
    }
    const double* px = ndof > 0 ? x.data() : nullptr;
    const double* pv = ndof > 0 ? xd.data() : nullptr;
    int rc = MOORDYN_SUCCESS;
    {
        const FPTrapMask mask;
        rc = compute_ic ? MoorDyn_Init(m_sys, px, pv) : MoorDyn_Init_NoIC(m_sys, px, pv);
    }
    if (rc != MOORDYN_SUCCESS) {
        return std::string(compute_ic ? "MoorDyn_Init" : "MoorDyn_Init_NoIC") + " failed with " + error_name(rc) +
               " for the system from '" + m_file + "'";
    }
    return std::string();
}

unsigned MoorDynSystem::external_kinematics_init (std::string& err)
{
    err.clear();
    unsigned n = 0;
    const int rc = MoorDyn_ExternalWaveKinInit(m_sys, &n);
    if (rc != MOORDYN_SUCCESS) {
        err = "MoorDyn_ExternalWaveKinInit failed with " + error_name(rc) + " for the system from '" + m_file + "'";
        return 0;
    }
    if (n == 0) {
        err = "the MoorDyn system from '" + m_file + "' takes no external fluid kinematics: its OPTIONS must set WaveKin = 1";
        return 0;
    }
    m_nkin = n;
    return n;
}

std::vector<double> MoorDynSystem::kinematics_points () const
{
    std::vector<double> r(3 * static_cast<std::size_t>(m_nkin), 0.0);
    if (m_nkin > 0) { check(MoorDyn_ExternalWaveKinGetCoordinates(m_sys, r.data()), "MoorDyn_ExternalWaveKinGetCoordinates"); }
    return r;
}

void MoorDynSystem::set_kinematics (const std::vector<double>& U, const std::vector<double>& Ud, double t)
{
    if (U.size() != 3 * static_cast<std::size_t>(m_nkin) || Ud.size() != U.size()) {
        amrex::Abort("MoorDynSystem::set_kinematics: " + std::to_string(3 * m_nkin) + " velocity and acceleration components are needed, " +
                     std::to_string(U.size()) + " and " + std::to_string(Ud.size()) + " were given (system from '" + m_file + "')");
    }
    std::size_t bad = 0;
    if (!all_finite(U, bad) || !all_finite(Ud, bad)) {
        amrex::Abort("MoorDynSystem::set_kinematics: a non-finite fluid velocity or acceleration at kinematics point " +
                     std::to_string(bad / 3) + " (0-based) for the system from '" + m_file + "'");
    }
    check(MoorDyn_ExternalWaveKinSet(m_sys, U.data(), Ud.data(), t), "MoorDyn_ExternalWaveKinSet");
}

void MoorDynSystem::step (const std::vector<double>& x, const std::vector<double>& xd, std::vector<double>& f, double& t, double dt)
{
    // MoorDyn returns the forces without stepping when dt <= 0
    if (!(std::isfinite(dt) && dt > 0.0)) {
        amrex::Abort("MoorDynSystem::step: the step must be finite and positive (s), " + std::to_string(dt) +
                     " given, for the system from '" + m_file + "'");
    }
    const unsigned ndof = num_coupled_dof();
    if (x.size() != ndof || xd.size() != ndof) {
        amrex::Abort("MoorDynSystem::step: " + std::to_string(ndof) + " coupled positions and velocities are needed, " +
                     std::to_string(x.size()) + " and " + std::to_string(xd.size()) + " were given (system from '" + m_file + "')");
    }
    f.assign(ndof, 0.0);
    const double* px = ndof > 0 ? x.data() : nullptr;
    const double* pv = ndof > 0 ? xd.data() : nullptr;
    double* pf = ndof > 0 ? f.data() : nullptr;
    double dt_in = dt;
    int rc = MOORDYN_SUCCESS;
    {
        const FPTrapMask mask;
        rc = MoorDyn_Step(m_sys, px, pv, pf, &t, &dt_in);
    }
    check(rc, "MoorDyn_Step");
    std::size_t bad = 0;
    if (!all_finite(f, bad)) {
        amrex::Abort("MoorDyn_Step returned a non-finite coupled force (component " + std::to_string(bad) +
                     ", 0-based) for the system from '" + m_file + "': reduce MoorDyn's internal step "
                     "(erf.conductors.moordyn_cfl or erf.conductors.moordyn_dt for a conductor line)");
    }
}

double MoorDynSystem::dt () const
{
    double v = 0.0;
    check(MoorDyn_GetDt(m_sys, &v), "MoorDyn_GetDt");
    return v;
}

void MoorDynSystem::set_dt (double v)
{
    if (!(std::isfinite(v) && v > 0.0)) {
        amrex::Abort("MoorDynSystem::set_dt: MoorDyn's internal step must be finite and positive (s), " + std::to_string(v) +
                     " given, for the system from '" + m_file + "'");
    }
    check(MoorDyn_SetDt(m_sys, v), "MoorDyn_SetDt");
}

unsigned MoorDynSystem::num_lines () const
{
    unsigned n = 0;
    check(MoorDyn_GetNumberLines(m_sys, &n), "MoorDyn_GetNumberLines");
    return n;
}

MoorDynLine MoorDynSystem::line (unsigned l) const
{
    MoorDynLine h = MoorDyn_GetLine(m_sys, l);
    if (h == nullptr) { amrex::Abort("MoorDyn_GetLine: no line " + std::to_string(l) + " in the system from '" + m_file + "'"); }
    return h;
}

MoorDynPoint MoorDynSystem::point (unsigned p) const
{
    MoorDynPoint h = MoorDyn_GetPoint(m_sys, p);
    if (h == nullptr) { amrex::Abort("MoorDyn_GetPoint: no point " + std::to_string(p) + " in the system from '" + m_file + "'"); }
    return h;
}

unsigned MoorDynSystem::line_num_nodes (unsigned l) const
{
    unsigned n = 0;
    check(MoorDyn_GetLineNumberNodes(line(l), &n), "MoorDyn_GetLineNumberNodes");
    return n;
}

double MoorDynSystem::line_unstretched_length (unsigned l) const
{
    double v = 0.0;
    check(MoorDyn_GetLineUnstretchedLength(line(l), &v), "MoorDyn_GetLineUnstretchedLength");
    return v;
}

std::array<double,3> MoorDynSystem::line_node_position (unsigned l, unsigned node) const
{
    std::array<double,3> r{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetLineNodePos(line(l), node, r.data()), "MoorDyn_GetLineNodePos");
    return r;
}

std::array<double,3> MoorDynSystem::line_node_velocity (unsigned l, unsigned node) const
{
    std::array<double,3> v{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetLineNodeVel(line(l), node, v.data()), "MoorDyn_GetLineNodeVel");
    return v;
}

std::array<double,3> MoorDynSystem::line_node_tension (unsigned l, unsigned node) const
{
    std::array<double,3> t{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetLineNodeTen(line(l), node, t.data()), "MoorDyn_GetLineNodeTen");
    return t;
}

std::array<double,3> MoorDynSystem::line_node_drag (unsigned l, unsigned node) const
{
    std::array<double,3> f{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetLineNodeDrag(line(l), node, f.data()), "MoorDyn_GetLineNodeDrag");
    return f;
}

std::array<double,3> MoorDynSystem::line_node_force (unsigned l, unsigned node) const
{
    std::array<double,3> f{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetLineNodeForce(line(l), node, f.data()), "MoorDyn_GetLineNodeForce");
    return f;
}

double MoorDynSystem::line_end_tension (unsigned l) const
{
    double v = 0.0;
    check(MoorDyn_GetLineFairTen(line(l), &v), "MoorDyn_GetLineFairTen");
    return v;
}

double MoorDynSystem::line_max_tension (unsigned l) const
{
    double v = 0.0;
    check(MoorDyn_GetLineMaxTen(line(l), &v), "MoorDyn_GetLineMaxTen");
    return v;
}

unsigned MoorDynSystem::num_points () const
{
    unsigned n = 0;
    check(MoorDyn_GetNumberPoints(m_sys, &n), "MoorDyn_GetNumberPoints");
    return n;
}

int MoorDynSystem::point_type (unsigned p) const
{
    int t = 0;
    check(MoorDyn_GetPointType(point(p), &t), "MoorDyn_GetPointType");
    return t;
}

std::array<double,3> MoorDynSystem::point_position (unsigned p) const
{
    std::array<double,3> r{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetPointPos(point(p), r.data()), "MoorDyn_GetPointPos");
    return r;
}

std::array<double,3> MoorDynSystem::point_force (unsigned p) const
{
    std::array<double,3> f{{0.0, 0.0, 0.0}};
    check(MoorDyn_GetPointForce(point(p), f.data()), "MoorDyn_GetPointForce");
    return f;
}

void MoorDynSystem::save (const std::string& path) { check(MoorDyn_Save(m_sys, path.c_str()), "MoorDyn_Save"); }

void MoorDynSystem::load (const std::string& path) { check(MoorDyn_Load(m_sys, path.c_str()), "MoorDyn_Load"); }

std::vector<std::uint64_t> MoorDynSystem::serialize () const
{
    std::size_t bytes = 0;
    check(MoorDyn_Serialize(m_sys, &bytes, nullptr), "MoorDyn_Serialize");
    std::vector<std::uint64_t> data((bytes + sizeof(std::uint64_t) - 1) / sizeof(std::uint64_t));
    check(MoorDyn_Serialize(m_sys, nullptr, data.data()), "MoorDyn_Serialize");
    m_serial_words = data.size();
    return data;
}

void MoorDynSystem::deserialize (const std::vector<std::uint64_t>& data)
{
    // MoorDyn_Deserialize reads the buffer without a size: a buffer from another system is undefined behaviour
    if (data.empty() || (m_serial_words > 0 && data.size() != m_serial_words)) {
        amrex::Abort("MoorDynSystem::deserialize: " + std::to_string(data.size()) + " words given, " + std::to_string(m_serial_words) +
                     " expected from serialize() for the system from '" + m_file + "'");
    }
    check(MoorDyn_Deserialize(m_sys, data.data()), "MoorDyn_Deserialize");
}

} // namespace erf_moordyn
