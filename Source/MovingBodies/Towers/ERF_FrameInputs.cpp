// Reading a SubDyn input file into FrameInputs, and the checks of a frame's inputs.

#include "ERF_FrameInputs.H"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <utility>

#include <AMReX_ParallelDescriptor.H>

namespace erf_towers {

namespace {

std::string upper (std::string s)
{
    for (auto& c : s) { c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
    return s;
}

std::string unquote (std::string s)
{
    if (s.size() >= 2 && (s.front() == '"' || s.front() == '\'') && s.back() == s.front()) { s = s.substr(1, s.size() - 2); }
    return s;
}

/** The whitespace- or comma-separated tokens of a line. */
std::vector<std::string> tokens (const std::string& line)
{
    std::string s = line;
    std::replace(s.begin(), s.end(), ',', ' ');
    std::istringstream is(s);
    std::vector<std::string> t;
    std::string w;
    while (is >> w) { t.push_back(w); }
    return t;
}

/** A number from a token, Fortran's D exponent accepted; false when the token is not one number. */
bool to_number (std::string s, double& v)
{
    for (auto& c : s) { if (c == 'd' || c == 'D') { c = 'e'; } }
    if (s.empty()) { return false; }
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

bool to_int (const std::string& s, int& v)
{
    double d = 0.0;
    if (!to_number(s, d) || d != std::floor(d) || std::abs(d) > 2.0e9) { return false; }
    v = static_cast<int>(d);
    return true;
}

bool dashed (const std::string& line)
{
    const auto p = line.find_first_not_of(" \t");
    return p != std::string::npos && line.compare(p, 3, "---") == 0;
}

/** Line-by-line reading of a SubDyn file, in the order SubDyn reads it. */
class Reader
{
public:
    Reader (std::vector<std::string> lines, std::string file) : m_lines(std::move(lines)), m_file(std::move(file)) {}

    /** Move past the next dashed section line that holds one of the keys (upper case); false when there is none. */
    bool section (const std::vector<std::string>& keys, std::string& title)
    {
        for (; m_pos < m_lines.size(); ++m_pos) {
            if (!dashed(m_lines[m_pos])) { continue; }
            const std::string u = upper(m_lines[m_pos]);
            for (const auto& k : keys) {
                if (u.find(k) != std::string::npos) { title = keys.front(); ++m_pos; return true; }
            }
        }
        return false;
    }

    /** The next line's tokens; false at the end of the file. */
    bool next (std::vector<std::string>& t, std::size_t& line_number)
    {
        if (m_pos >= m_lines.size()) { return false; }
        line_number = m_pos + 1;
        t = tokens(m_lines[m_pos++]);
        return true;
    }

    /** The lines from here to the next dashed line (not consumed). */
    std::vector<std::vector<std::string>> block ()
    {
        std::vector<std::vector<std::string>> b;
        while (m_pos < m_lines.size() && !dashed(m_lines[m_pos])) { b.push_back(tokens(m_lines[m_pos++])); }
        return b;
    }

    /**
     * A table as SubDyn writes it: a count line ("n Name"), two header lines, then n rows of at least
     * min_columns tokens. Returns an empty string, or the problem.
     */
    std::string table (const std::string& title, std::size_t min_columns, std::vector<std::vector<std::string>>& rows,
                       std::vector<std::size_t>& row_lines)
    {
        std::vector<std::string> t;
        std::size_t ln = 0;
        int n = 0;
        if (!next(t, ln) || t.empty() || !to_int(t[0], n) || n < 0) {
            return m_file + ": the " + title + " section needs its count on the line after its title";
        }
        for (int h = 0; h < 2; ++h) {
            if (!next(t, ln)) { return m_file + ": the " + title + " section ends before its two header lines"; }
        }
        rows.clear();
        row_lines.clear();
        for (int r = 0; r < n; ++r) {
            if (!next(t, ln)) { return m_file + ": the " + title + " section lists " + std::to_string(n) + " rows but the file ends"; }
            if (t.size() < min_columns) {
                return m_file + " line " + std::to_string(ln) + ": a row of the " + title + " section needs at least " +
                       std::to_string(min_columns) + " values";
            }
            rows.push_back(t);
            row_lines.push_back(ln);
        }
        return std::string();
    }

    const std::string& file () const { return m_file; }

private:
    std::vector<std::string> m_lines;
    std::string m_file;
    std::size_t m_pos{0};
};

std::string at_line (const std::string& file, std::size_t ln) { return file + " line " + std::to_string(ln) + ": "; }

/** Numbers from tokens[first, first+count) into out; false when one is not a number. */
bool numbers (const std::vector<std::string>& t, std::size_t first, std::size_t count, double* out)
{
    for (std::size_t i = 0; i < count; ++i) {
        if (first + i >= t.size() || !to_number(t[first + i], out[i])) { return false; }
    }
    return true;
}

std::string read_ssi (const std::string& path, FrameSupport& s)
{
    static const char* knames[21] = {"KXX", "KXY", "KYY", "KXZ", "KYZ", "KZZ", "KXTX", "KYTX", "KZTX", "KTXTX", "KXTY",
                                     "KYTY", "KZTY", "KTXTY", "KTYTY", "KXTZ", "KYTZ", "KZTZ", "KTXTZ", "KTYTZ", "KTZTZ"};
    std::string text;
    if (!read_text_file(path, text)) { return "cannot read the SSI file '" + path + "'"; }
    std::istringstream f(text);
    std::string line;
    int line_no = 0, entries = 0;
    while (std::getline(f, line)) {
        ++line_no;
        const auto t = tokens(line);
        if (t.empty()) { continue; }
        double v = 0.0;
        // comment and header lines start with text; an entry is a value then its name
        if (!to_number(t[0], v)) { continue; }
        std::string name = (t.size() >= 2) ? upper(t[1]) : std::string();
        const bool is_mass = !name.empty() && name[0] == 'M';
        if (is_mass) { name[0] = 'K'; }
        bool known = false;
        for (int i = 0; i < 21; ++i) {
            if (name == knames[i]) { (is_mass ? s.mass : s.stiffness)[static_cast<std::size_t>(i)] = v; known = true; }
        }
        if (!known) {
            return "the SSI file '" + path + "', line " + std::to_string(line_no) + ": '" + (t.size() >= 2 ? t[1] : std::string()) +
                   "' is not one of Kxx Kxy Kyy Kxz Kyz Kzz Kxtx Kytx Kztx Ktxtx Kxty Kyty Kzty Ktxty Ktyty Kxtz Kytz Kztz Ktxtz "
                   "Ktytz Ktztz or their M names";
        }
        ++entries;
    }
    if (entries == 0) { return "the SSI file '" + path + "' holds no stiffness or mass entry (a value then its name per line)"; }
    return std::string();
}

} // namespace

int FrameInputs::joint_index (int id) const
{
    for (std::size_t i = 0; i < joints.size(); ++i) { if (joints[i].id == id) { return static_cast<int>(i); } }
    return -1;
}

const FrameSection* FrameInputs::section (int id, SectionShape shape) const
{
    for (const auto& s : sections) { if (s.id == id && s.shape == shape) { return &s; } }
    return nullptr;
}

bool read_text_file (const std::string& path, std::string& text)
{
    long n = -1;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream f(path, std::ios::binary);
        if (f) {
            std::ostringstream all;
            all << f.rdbuf();
            text = all.str();
            n = static_cast<long>(text.size());
        }
    }
    amrex::ParallelDescriptor::Bcast(&n, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (n < 0) { return false; }
    text.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(&text[0], static_cast<std::size_t>(n), amrex::ParallelDescriptor::IOProcessorNumber()); }
    return true;
}

std::string read_subdyn (const std::string& path, FrameInputs& in)
{
    std::string text;
    if (!read_text_file(path, text)) { return "cannot read the SubDyn file '" + path + "'"; }
    std::istringstream f(text);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(f, line)) { if (!line.empty() && line.back() == '\r') { line.pop_back(); } lines.push_back(line); }
    in = FrameInputs();
    in.file = path;
    const std::string dir = (path.find('/') == std::string::npos) ? std::string() : path.substr(0, path.rfind('/') + 1);
    Reader r(std::move(lines), path);
    std::string title, err;
    std::vector<std::vector<std::string>> rows;
    std::vector<std::size_t> row_lines;

    // FEA and Craig-Bampton parameters: only FEMMod and NDiv matter for a static frame
    if (!r.section({"FEA"}, title)) { return path + ": no FEA and CRAIG-BAMPTON PARAMETERS section"; }
    int femmod = -1;
    bool have_ndiv = false;
    for (const auto& t : r.block()) {
        if (t.size() < 2) { continue; }
        const std::string name = upper(t[1]);
        if (name == "FEMMOD" && !to_int(t[0], femmod)) { return path + ": FEMMod is not an integer"; }
        if (name == "NDIV") {
            if (!to_int(t[0], in.divisions)) { return path + ": NDiv is not an integer"; }
            have_ndiv = true;
        }
    }
    if (femmod == 1) { in.theory = BeamTheory::EulerBernoulli; }
    else if (femmod == 3) { in.theory = BeamTheory::Timoshenko; }
    else { return path + ": FEMMod = " + std::to_string(femmod) + "; the frame model has Euler-Bernoulli (1) and Timoshenko (3) beams"; }
    if (!have_ndiv) { return path + ": no NDiv in the FEA section"; }

    // joints
    if (!r.section({"STRUCTURE JOINTS"}, title)) { return path + ": no STRUCTURE JOINTS section"; }
    if (!(err = r.table("STRUCTURE JOINTS", 4, rows, row_lines)).empty()) { return err; }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        FrameJoint j;
        if (!to_int(rows[i][0], j.id) || !numbers(rows[i], 1, 3, j.x.data())) {
            return at_line(path, row_lines[i]) + "a joint needs an integer id and x, y, z (m)";
        }
        int type = 1;
        if (rows[i].size() >= 5 && (!to_int(rows[i][4], type) || type != 1)) {
            return at_line(path, row_lines[i]) + "joint " + std::to_string(j.id) + " has JointType " + rows[i][4] +
                   "; the frame model joins members rigidly (JointType 1, cantilever) only";
        }
        in.joints.push_back(j);
    }

    // base reactions: 1 fixed, 0 free; an SSI file adds a spring on the free DOFs
    if (!r.section({"BASE REACTION"}, title)) { return path + ": no BASE REACTION JOINTS section"; }
    if (!(err = r.table("BASE REACTION JOINTS", 1, rows, row_lines)).empty()) { return err; }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        FrameSupport s;
        if (!to_int(rows[i][0], s.joint)) { return at_line(path, row_lines[i]) + "a reaction row needs an integer joint id"; }
        if (rows[i].size() >= 7) {
            for (int d = 0; d < 6; ++d) {
                int flag = -1;
                if (!to_int(rows[i][static_cast<std::size_t>(d) + 1], flag) || (flag != 0 && flag != 1)) {
                    return at_line(path, row_lines[i]) + "the six reaction flags of joint " + std::to_string(s.joint) +
                           " must be 1 (fixed) or 0 (free)";
                }
                s.fixed[static_cast<std::size_t>(d)] = (flag == 1);
            }
            if (rows[i].size() >= 8) {
                s.ssi_file = unquote(rows[i][7]);
                if (!s.ssi_file.empty() && s.ssi_file[0] != '/') { s.ssi_file = dir + s.ssi_file; }
                if (!(err = read_ssi(s.ssi_file, s)).empty()) { return at_line(path, row_lines[i]) + err; }
            }
        } else if (rows[i].size() != 1) {
            return at_line(path, row_lines[i]) + "a reaction row has a joint id and six flags (and an optional SSI file)";
        }
        in.supports.push_back(s);
    }

    if (!r.section({"INTERFACE JOINTS"}, title)) { return path + ": no INTERFACE JOINTS section"; }
    if (!(err = r.table("INTERFACE JOINTS", 1, rows, row_lines)).empty()) { return err; }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        int id = 0;
        if (!to_int(rows[i][0], id)) { return at_line(path, row_lines[i]) + "an interface row needs an integer joint id"; }
        in.interface_joints.push_back(id);
    }

    // members: id, joints, property sets at both ends, type, spin (deg) or COSM id
    if (!r.section({"MEMBERS"}, title)) { return path + ": no MEMBERS section"; }
    if (!(err = r.table("MEMBERS", 7, rows, row_lines)).empty()) { return err; }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& t = rows[i];
        FrameMember m;
        int p2 = 0;
        if (!to_int(t[0], m.id) || !to_int(t[1], m.joint_a) || !to_int(t[2], m.joint_b) || !to_int(t[3], m.section) || !to_int(t[4], p2)) {
            return at_line(path, row_lines[i]) +
                   "a member row starts with five integers: MemberID MJointID1 MJointID2 MPropSetID1 MPropSetID2";
        }
        const std::string type = upper(t[5]);
        if (type == "1" || type == "1C") { m.shape = SectionShape::Circular; }
        else if (type == "1R") { m.shape = SectionShape::Rectangular; }
        else if (type == "4") { m.shape = SectionShape::Arbitrary; }
        else {
            return at_line(path, row_lines[i]) + "member " + std::to_string(m.id) + " has MType " + t[5] +
                   "; the frame model has beam members only (1c circular, 1r rectangular, 4 arbitrary), "
                   "not cables (2), rigid links (3) or springs (5)";
        }
        if (p2 != m.section) {
            return at_line(path, row_lines[i]) + "member " + std::to_string(m.id) + " has different property sets at its ends (" +
                   std::to_string(m.section) + ", " + std::to_string(p2) + "); the frame model has uniform members only";
        }
        double spin_deg = 0.0;
        if (!to_number(t[6], spin_deg)) {
            return at_line(path, row_lines[i]) + "member " + std::to_string(m.id) + ": MSpin must be a number (deg)";
        }
        m.spin = spin_deg * 3.14159265358979323846 / 180.0;
        in.members.push_back(m);
    }

    // the three cross-section tables
    const std::pair<std::vector<std::string>, SectionShape> tables[] = {
        {{"CIRCULAR", "1/3"}, SectionShape::Circular}, {{"RECTANGULAR", "2/3"}, SectionShape::Rectangular},
        {{"ARBITRARY", "3/3"}, SectionShape::Arbitrary}};
    // each table's title, as the messages name it
    auto table_title = [] (const std::string& key) { return key + " BEAM CROSS-SECTION PROPERTIES"; };
    for (const auto& tb : tables) {
        if (!r.section(tb.first, title)) { return path + ": no " + table_title(tb.first.front()) + " section"; }
        const std::size_t ncol = (tb.second == SectionShape::Circular) ? 6 : (tb.second == SectionShape::Rectangular) ? 7 : 11;
        if (!(err = r.table(table_title(tb.first.front()), ncol, rows, row_lines)).empty()) { return err; }
        for (std::size_t i = 0; i < rows.size(); ++i) {
            FrameSection s;
            s.shape = tb.second;
            double v[10] = {};
            if (!to_int(rows[i][0], s.id) || !numbers(rows[i], 1, ncol - 1, v)) {
                return at_line(path, row_lines[i]) + "a row of " + table_title(tb.first.front()) + " needs an integer id and " +
                       std::to_string(ncol - 1) + " numbers";
            }
            s.E = v[0]; s.G = v[1]; s.rho = v[2];
            if (tb.second == SectionShape::Circular) { s.D = v[3]; s.t = v[4]; }
            else if (tb.second == SectionShape::Rectangular) { s.Sa = v[3]; s.Sb = v[4]; s.t = v[5]; }
            else { s.A = v[3]; s.Asx = v[4]; s.Asy = v[5]; s.Ixx = v[6]; s.Iyy = v[7]; s.J0 = v[8]; s.Jt = v[9]; }
            in.sections.push_back(s);
        }
    }

    // cable, rigid-link and spring properties and the cosine matrices: no member of the model uses them
    const std::pair<std::vector<std::string>, std::string> skipped[] = {
        {{"CABLE"}, "CABLE PROPERTIES"}, {{"RIGID"}, "RIGID LINK PROPERTIES"}, {{"SPRING"}, "SPRING ELEMENT PROPERTIES"},
        {{"COSINE", "COSM"}, "MEMBER COSINE MATRICES"}};
    for (const auto& sk : skipped) {
        if (!r.section(sk.first, title)) { return path + ": no " + sk.second + " section"; }
        if (!(err = r.table(sk.second, 1, rows, row_lines)).empty()) { return err; }
    }

    if (!r.section({"CONCENTRATED MASS"}, title)) { return path + ": no JOINT ADDITIONAL CONCENTRATED MASSES section"; }
    if (!(err = r.table("JOINT ADDITIONAL CONCENTRATED MASSES", 5, rows, row_lines)).empty()) { return err; }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        FrameMass m;
        double v[10] = {};
        // SubDyn's 5 columns (joint, mass, three inertias) or 11 (the products of inertia and the centre's offset too)
        const std::size_t nv = (rows[i].size() >= 11) ? 10 : 4;
        if (rows[i].size() != 5 && rows[i].size() < 11) {
            return at_line(path, row_lines[i]) + "a concentrated mass row has " + std::to_string(rows[i].size()) +
                   " values; give 5 (joint id, JMass, JMXX JMYY JMZZ) or 11 (also JMXY JMXZ JMYZ MCGX MCGY MCGZ)";
        }
        if (!to_int(rows[i][0], m.joint) || !numbers(rows[i], 1, nv, v)) {
            return at_line(path, row_lines[i]) +
                   "a concentrated mass row needs a joint id, JMass and JMXX JMYY JMZZ (and optionally JMXY JMXZ JMYZ MCGX MCGY MCGZ)";
        }
        m.mass = v[0];
        m.inertia = {{v[1], v[2], v[3], v[4], v[5], v[6]}};
        m.offset = {{v[7], v[8], v[9]}};
        in.masses.push_back(m);
    }
    return std::string();
}

std::string FrameInputs::validate () const
{
    const std::string where = file.empty() ? std::string("frame") : file;
    if (divisions < 1) { return where + ": NDiv must be >= 1"; }
    if (joints.size() < 2) { return where + ": a frame needs at least two joints"; }
    std::set<int> ids;
    for (const auto& j : joints) {
        if (!ids.insert(j.id).second) { return where + ": joint id " + std::to_string(j.id) + " appears twice"; }
        if (!(std::isfinite(j.x[0]) && std::isfinite(j.x[1]) && std::isfinite(j.x[2]))) {
            return where + ": joint " + std::to_string(j.id) + " has a non-finite position";
        }
    }
    std::set<std::pair<int,int>> sec_ids;
    for (const auto& s : sections) {
        const std::string key = where + ": property set " + std::to_string(s.id) + " (" +
            (s.shape == SectionShape::Circular ? "circular" : s.shape == SectionShape::Rectangular ? "rectangular" : "arbitrary") + ") ";
        if (!sec_ids.insert({s.id, static_cast<int>(s.shape)}).second) { return key + "appears twice in its table"; }
        if (!(std::isfinite(s.E) && s.E > 0.0)) { return key + "needs YoungE > 0 (Pa)"; }
        if (!(std::isfinite(s.G) && s.G > 0.0)) { return key + "needs ShearG > 0 (Pa)"; }
        if (!(std::isfinite(s.rho) && s.rho >= 0.0)) { return key + "needs MatDens >= 0 (kg/m^3)"; }
        if (s.shape == SectionShape::Circular) {
            if (!(std::isfinite(s.D) && s.D > 0.0)) { return key + "needs XsecD > 0 (m)"; }
            if (!(std::isfinite(s.t) && s.t >= 0.0 && s.t < 0.5 * s.D)) { return key + "needs 0 <= XsecT < XsecD/2 (m; 0: solid)"; }
        } else if (s.shape == SectionShape::Rectangular) {
            if (!(std::isfinite(s.Sa) && s.Sa > 0.0 && std::isfinite(s.Sb) && s.Sb > 0.0)) {
                return key + "needs XsecSa > 0 and XsecSb > 0 (m)";
            }
            if (!(std::isfinite(s.t) && s.t >= 0.0 && 2.0 * s.t < std::min(s.Sa, s.Sb))) {
                return key + "needs 0 <= XsecT < min(XsecSa, XsecSb)/2 (m; 0: solid)";
            }
        } else {
            const std::pair<const char*, double> pos[] = {{"XsecA", s.A}, {"XsecJxx", s.Ixx}, {"XsecJyy", s.Iyy}, {"XsecJt", s.Jt}};
            for (const auto& kv : pos) {
                if (!(std::isfinite(kv.second) && kv.second > 0.0)) { return key + "needs " + kv.first + " > 0"; }
            }
            if (!(std::isfinite(s.J0) && s.J0 >= 0.0)) { return key + "needs XsecJ0 >= 0 (m^4)"; }
            if (theory == BeamTheory::Timoshenko && !(std::isfinite(s.Asx) && s.Asx > 0.0 && std::isfinite(s.Asy) && s.Asy > 0.0)) {
                return key + "needs XsecAsx > 0 and XsecAsy > 0 (m^2) for Timoshenko beams (FEMMod 3)";
            }
        }
    }
    if (members.empty()) { return where + ": a frame needs at least one member"; }
    std::set<int> mids;
    std::vector<bool> used(joints.size(), false);
    for (const auto& m : members) {
        const std::string key = where + ": member " + std::to_string(m.id) + " ";
        if (!mids.insert(m.id).second) { return key + "appears twice"; }
        const int a = joint_index(m.joint_a), b = joint_index(m.joint_b);
        if (a < 0 || b < 0) {
            return key + "joins joint " + std::to_string(a < 0 ? m.joint_a : m.joint_b) + ", which is not in STRUCTURE JOINTS";
        }
        const auto& xa = joints[static_cast<std::size_t>(a)].x;
        const auto& xb = joints[static_cast<std::size_t>(b)].x;
        if (xa == xb) {
            return key + "has zero length: joints " + std::to_string(m.joint_a) + " and " + std::to_string(m.joint_b) + " coincide";
        }
        if (section(m.section, m.shape) == nullptr) {
            return key + "uses property set " + std::to_string(m.section) + ", which is not in the table of its member type";
        }
        if (!std::isfinite(m.spin)) { return key + "has a non-finite MSpin"; }
        used[static_cast<std::size_t>(a)] = true;
        used[static_cast<std::size_t>(b)] = true;
    }
    for (std::size_t i = 0; i < joints.size(); ++i) {
        if (!used[i]) { return where + ": joint " + std::to_string(joints[i].id) + " is on no member"; }
    }
    if (supports.empty()) { return where + ": a frame needs at least one base reaction joint"; }
    std::set<int> sj;
    for (const auto& s : supports) {
        const std::string key = where + ": base reaction joint " + std::to_string(s.joint) + " ";
        if (joint_index(s.joint) < 0) { return key + "is not in STRUCTURE JOINTS"; }
        if (!sj.insert(s.joint).second) { return key + "appears twice"; }
        for (const double k : s.stiffness) { if (!std::isfinite(k)) { return key + "has a non-finite SSI stiffness"; } }
        for (const double k : s.mass) { if (!std::isfinite(k)) { return key + "has a non-finite SSI mass"; } }
        if (std::none_of(s.fixed.begin(), s.fixed.end(), [] (bool b) { return b; }) &&
            std::all_of(s.stiffness.begin(), s.stiffness.end(), [] (double k) { return k == 0.0; })) {
            return key + "restrains nothing: every DOF is free (flag 0) and it has no SSI spring";
        }
        if (!s.ssi_file.empty() && std::all_of(s.fixed.begin(), s.fixed.end(), [] (bool b) { return b; })) {
            return key + "has an SSI file but no free DOF: the spring acts on free DOFs only (flag 0)";
        }
    }
    for (const int id : interface_joints) {
        if (joint_index(id) < 0) { return where + ": interface joint " + std::to_string(id) + " is not in STRUCTURE JOINTS"; }
    }
    if (!temperature.empty() && temperature.size() != members.size()) {
        return where + ": " + std::to_string(temperature.size()) + " member temperatures for " + std::to_string(members.size()) +
               " members";
    }
    for (std::size_t m = 0; m < temperature.size(); ++m) {
        if (!(std::isfinite(temperature[m]) && temperature[m] < 1200.0)) {
            return where + ": member " + std::to_string(members[m].id) + " has a temperature of " + std::to_string(temperature[m]) +
                   " C; the steel keeps no stiffness at 1200 C (EN 1993-1-2), so it must be below 1200 C";
        }
    }
    for (const auto& m : masses) {
        if (joint_index(m.joint) < 0) {
            return where + ": concentrated mass joint " + std::to_string(m.joint) + " is not in STRUCTURE JOINTS";
        }
        if (!(std::isfinite(m.mass) && m.mass >= 0.0)) {
            return where + ": the concentrated mass at joint " + std::to_string(m.joint) + " must be >= 0 (kg)";
        }
        for (const double v : m.offset) {
            if (!std::isfinite(v)) {
                return where + ": the concentrated mass at joint " + std::to_string(m.joint) + " has a non-finite offset";
            }
        }
        for (const double v : m.inertia) {
            if (!std::isfinite(v)) {
                return where + ": the concentrated mass at joint " + std::to_string(m.joint) + " has a non-finite inertia";
            }
        }
    }
    return std::string();
}

namespace {

/** A value at 17 significant digits, which reads back to the same double. */
std::string exact (double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << std::scientific << v;
    return os.str();
}

/** The upper-case names of SubDyn's 21 SSI entries, the upper triangle column by column. */
const char* const ssi_names[21] = {"Kxx", "Kxy", "Kyy", "Kxz", "Kyz", "Kzz", "Kxtx", "Kytx", "Kztx", "Ktxtx", "Kxty",
                                   "Kyty", "Kzty", "Ktxty", "Ktyty", "Kxtz", "Kytz", "Kztz", "Ktxtz", "Ktytz", "Ktztz"};

} // namespace

std::string write_subdyn (const FrameInputs& in, const std::string& path, const std::string& title)
{
    const std::string err = in.validate();
    if (!err.empty()) { return err; }
    if (in.interface_joints.empty()) { return path + ": SubDyn needs at least one interface joint"; }
    const std::string stem = (path.size() > 4 && path.compare(path.size() - 4, 4, ".dat") == 0) ? path.substr(0, path.size() - 4) : path;
    const std::string base = stem.substr(stem.rfind('/') == std::string::npos ? 0 : stem.rfind('/') + 1);
    std::ostringstream o;
    o << "----------- SubDyn MultiMember Support Structure Input File ---------------------------\n" << title << "\n"
      << "-------------------------- SIMULATION CONTROL -----------------------------------------\n"
      << "False            Echo        - Echo input data to \"<rootname>.SD.ech\" (flag)\n"
      << "\"DEFAULT\"        SDdeltaT    - Local Integration Step. If \"default\", the glue-code integration step will be used.\n"
      << "             3   IntMethod   - Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2].\n"
      << "False            SttcSolve   - Solve dynamics about static equilibrium point\n"
      << "-------------------- FEA and CRAIG-BAMPTON PARAMETERS ---------------------------------\n"
      << "             " << (in.theory == BeamTheory::EulerBernoulli ? 1 : 3) << "   FEMMod      - FEM switch: element model in the FEM. "
      << "[1= Euler-Bernoulli(E-B);  2=Tapered E-B (unavailable);  3= 2-node Timoshenko;  4= 2-node tapered Timoshenko (unavailable)]\n"
      << "             " << in.divisions << "   NDiv        - Number of sub-elements per member\n"
      << "             0   Nmodes      - Number of internal modes to retain. If Nmodes=0 --> Guyan Reduction. "
         "If Nmodes<0 --> retain all modes.\n"
      << "             0   JDampings   - Damping Ratios for each retained mode (% of critical)\n"
      << "             0   GuyanDampMod - Guyan damping {0=none, 1=Rayleigh Damping, 2=user specified 6x6 matrix}\n"
      << "  0.000, 0.000   RayleighDamp - Mass and stiffness proportional damping coefficients (Rayleigh Damping) "
         "[only if GuyanDampMod=1]\n"
      << "             6   GuyanDampSize - Guyan damping matrix (6x6) [only if GuyanDampMod=2]\n";
    for (int r = 0; r < 6; ++r) { o << "   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00\n"; }
    o << "------- INITIAL RIGID-BODY POSITION [used only for floating structure with more than one transition pieces] -------\n"
      << "RBSurge    RBSway     RBHeave    RBRoll     RBPitch    RBYaw\n"
      << "  (m)        (m)        (m)      (deg)      (deg)      (deg)\n"
      << "  0.0        0.0        0.0       0.0        0.0        0.0\n"
      << "---- STRUCTURE JOINTS: joints connect structure members (~Hydrodyn Input File) --------\n"
      << "   " << in.joints.size() << "   NJoints     - Number of joints (-)\n"
      << "JointID   JointXss   JointYss   JointZss   JointType   JointDirX   JointDirY   JointDirZ   JointStiff\n"
      << "  (-)       (m)        (m)        (m)         (-)         (-)         (-)         (-)       (Nm/rad)\n";
    for (const auto& j : in.joints) {
        o << "  " << j.id << "   " << exact(j.x[0]) << "   " << exact(j.x[1]) << "   " << exact(j.x[2]) << "   1   0.0   0.0   0.0   0.0\n";
    }
    o << "------------------- BASE REACTION JOINTS: 1/0 for Locked/Free DOF @ each Reaction Node ---------------------\n"
      << "   " << in.supports.size() << "   NReact      - Number of Joints with reaction forces\n"
      << "RJointID   RctTDXss    RctTDYss    RctTDZss    RctRDXss    RctRDYss    RctRDZss     SSIfile\n"
      << "  (-)       (flag)      (flag)      (flag)      (flag)      (flag)      (flag)      (string)\n";
    for (const auto& s : in.supports) {
        o << "   " << s.joint;
        for (const bool f : s.fixed) { o << "   " << (f ? 1 : 0); }
        const bool spring = std::any_of(s.stiffness.begin(), s.stiffness.end(), [] (double v) { return v != 0.0; }) ||
                            std::any_of(s.mass.begin(), s.mass.end(), [] (double v) { return v != 0.0; });
        if (spring) {
            const std::string ssi = stem + "_ssi_" + std::to_string(s.joint) + ".dat";
            std::ofstream f(ssi, std::ios::trunc);
            if (!f) { return "cannot write the SSI file '" + ssi + "'"; }
            f << "! SSI stiffness and mass of base reaction joint " << s.joint << ": value, then name\n";
            for (int k = 0; k < 21; ++k) { f << exact(s.stiffness[static_cast<std::size_t>(k)]) << "   " << ssi_names[k] << "\n"; }
            for (int k = 0; k < 21; ++k) {
                std::string name = ssi_names[k];
                name[0] = 'M';
                f << exact(s.mass[static_cast<std::size_t>(k)]) << "   " << name << "\n";
            }
            if (!f) { return "cannot write the SSI file '" + ssi + "'"; }
            o << "   \"" << base << "_ssi_" << s.joint << ".dat\"";
        }
        o << "\n";
    }
    o << "------- INTERFACE JOINTS: 1/0 for Locked (to the TP)/Free DOF @each Interface Joint "
         "(only Locked-to-TP implemented thus far (=rigid TP)) ---------\n"
      << "   " << in.interface_joints.size() << "   NInterf     - Number of interface joints locked to the Transition Piece (TP)\n"
      << "IJointID   TPID   ItfTDXss    ItfTDYss    ItfTDZss    ItfRDXss    ItfRDYss    ItfRDZss\n"
      << "  (-)      (-)     (flag)      (flag)      (flag)      (flag)      (flag)      (flag)\n";
    for (const int id : in.interface_joints) { o << "   " << id << "   1   1   1   1   1   1   1\n"; }
    o << "----------------------------------- MEMBERS -------------------------------------------\n"
      << "   " << in.members.size() << "   NMembers    - Number of members (-)\n"
      << "MemberID   MJointID1   MJointID2   MPropSetID1   MPropSetID2   MType   COSMID/MSpin\n"
      << "  (-)         (-)         (-)          (-)           (-)        (-)    (-)/(deg)\n";
    for (const auto& m : in.members) {
        const char* type = (m.shape == SectionShape::Circular) ? "1c" : (m.shape == SectionShape::Rectangular) ? "1r" : "4";
        o << "   " << m.id << "   " << m.joint_a << "   " << m.joint_b << "   " << m.section << "   " << m.section << "   " << type
          << "   " << exact(m.spin * 180.0 / 3.14159265358979323846) << "\n";
    }
    auto count = [&] (SectionShape shape) {
        return std::count_if(in.sections.begin(), in.sections.end(), [shape] (const FrameSection& s) { return s.shape == shape; });
    };
    auto material = [] (const FrameSection& s) { return exact(s.E) + "   " + exact(s.G) + "   " + exact(s.rho); };
    o << "------------------ CIRCULAR BEAM CROSS-SECTION PROPERTIES -----------------------------\n"
      << "   " << count(SectionShape::Circular) << "   NPropSetsCyl - Number of structurally unique circular cross-sections\n"
      << "PropSetID     YoungE          ShearG          MatDens          XsecD           XsecT\n"
      << "  (-)         (N/m2)          (N/m2)          (kg/m3)           (m)             (m)\n";
    for (const auto& s : in.sections) {
        if (s.shape == SectionShape::Circular) {
            o << "   " << s.id << "   " << material(s) << "   " << exact(s.D) << "   " << exact(s.t) << "\n";
        }
    }
    o << "----------------- RECTANGULAR BEAM CROSS-SECTION PROPERTIES ---------------------------\n"
      << "   " << count(SectionShape::Rectangular) << "   NPropSetsRec - Number of structurally unique rectangular cross-sections\n"
      << "PropSetID     YoungE          ShearG          MatDens          XsecSa         XsecSb          XsecT\n"
      << "  (-)         (N/m2)          (N/m2)          (kg/m3)           (m)            (m)             (m)\n";
    for (const auto& s : in.sections) {
        if (s.shape == SectionShape::Rectangular) {
            o << "   " << s.id << "   " << material(s) << "   " << exact(s.Sa) << "   " << exact(s.Sb) << "   " << exact(s.t) << "\n";
        }
    }
    o << "----------------- ARBITRARY BEAM CROSS-SECTION PROPERTIES -----------------------------\n"
      << "   " << count(SectionShape::Arbitrary) << "   NXPropSets   - Number of structurally unique arbitrary cross-sections\n"
      << "PropSetID     YoungE          ShearG          MatDens          XsecA          XsecAsx       XsecAsy"
         "       XsecJxx       XsecJyy        XsecJ0    XsecJt\n"
      << "  (-)         (N/m2)          (N/m2)          (kg/m3)          (m2)            (m2)          (m2)"
         "          (m4)          (m4)          (m4)       (m4)\n";
    for (const auto& s : in.sections) {
        if (s.shape == SectionShape::Arbitrary) {
            o << "   " << s.id << "   " << material(s) << "   " << exact(s.A) << "   " << exact(s.Asx) << "   " << exact(s.Asy) << "   "
              << exact(s.Ixx) << "   " << exact(s.Iyy) << "   " << exact(s.J0) << "   " << exact(s.Jt) << "\n";
        }
    }
    o << "-------------------------- CABLE PROPERTIES -------------------------------------------\n"
      << "             0   NCablePropSets   - Number of cable cable properties\n"
      << "PropSetID     EA          MatDens        T0         CtrlChannel\n"
      << "  (-)         (N)         (kg/m)        (N)             (-)\n"
      << "----------------------- RIGID LINK PROPERTIES -----------------------------------------\n"
      << "             0   NRigidPropSets - Number of rigid link properties\n"
      << "PropSetID   MatDens\n"
      << "  (-)       (kg/m)\n"
      << "----------------------- SPRING ELEMENT PROPERTIES -------------------------------------\n"
      << "             0   NSpringPropSets - Number of spring properties\n"
      << "PropSetID   k11     k12     k13     k14     k15     k16     k22     k23     k24     k25     k26     k33     k34     k35     k36"
         "     k44      k45      k46      k55      k56      k66\n"
      << "  (-)      (N/m)   (N/m)   (N/m)  (N/rad) (N/rad) (N/rad)  (N/m)   (N/m)  (N/rad) (N/rad) (N/rad)  (N/m)  (N/rad) (N/rad) (N/rad)"
         " (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad)\n"
      << "---------------------- MEMBER COSINE MATRICES COSM(i,j) -------------------------------\n"
      << "             0   NCOSMs      - Number of unique cosine matrices\n"
      << "COSMID    COSM11    COSM12    COSM13    COSM21    COSM22    COSM23    COSM31    COSM32    COSM33\n"
      << " (-)       (-)       (-)       (-)       (-)       (-)       (-)       (-)       (-)       (-)\n"
      << "------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------\n"
      << "   " << in.masses.size() << "   NCmass      - Number of joints with concentrated masses; Global Coordinate System\n"
      << "CMJointID       JMass            JMXX             JMYY             JMZZ          JMXY        JMXZ"
         "         JMYZ        MCGX      MCGY        MCGZ\n"
      << "  (-)            (kg)          (kg*m^2)         (kg*m^2)         (kg*m^2)      (kg*m^2)    (kg*m^2)"
         "     (kg*m^2)       (m)      (m)          (m)\n";
    for (const auto& c : in.masses) {
        o << "   " << c.joint << "   " << exact(c.mass);
        for (const double v : c.inertia) { o << "   " << exact(v); }
        for (const double v : c.offset) { o << "   " << exact(v); }
        o << "\n";
    }
    o << "---------------------------- OUTPUT: SUMMARY & OUTFILE --------------------------------\n"
      << "True             SumPrint    - Output a Summary File (flag)\n"
      << "0                OutCBModes  - Output Guyan and Craig-Bampton modes {0: No output, 1: JSON output}, (flag)\n"
      << "0                OutFEMModes - Output first 30 FEM modes {0: No output, 1: JSON output} (flag)\n"
      << "False            OutCOSM     - Output cosine matrices with the selected output member forces (flag)\n"
      << "False            OutAll      - [T/F] Output all members' end forces\n"
      << "             1   OutSwtch    - [1/2/3] Output requested channels to: 1=<rootname>.SD.out;  "
         "2=<rootname>.out (generated by FAST);  3=both files.\n"
      << "True             TabDelim    - Generate a tab-delimited output in the <rootname>.SD.out file\n"
      << "             1   OutDec      - Decimation of output in the <rootname>.SD.out file\n"
      << "\"ES11.4e2\"       OutFmt      - Output format for numerical results in the <rootname>.SD.out file\n"
      << "\"A11\"            OutSFmt     - Output format for header strings in the <rootname>.SD.out file\n"
      << "------------------------- MEMBER OUTPUT LIST ------------------------------------------\n"
      << "             0   NMOutputs   - Number of members whose forces/displacements/velocities/accelerations "
         "will be output (-) [Must be <= 99].\n"
      << "MemberID   NOutCnt    NodeCnt\n"
      << "  (-)        (-)        (-)\n"
      << "------------------------- SSOutList: The next line(s) contains a list of output parameters. ------\n"
      << "END of output channels and end of file. (the word \"END\" must appear in the first 3 columns of this line)\n";
    std::ofstream f(path, std::ios::trunc);
    if (!f) { return "cannot write the SubDyn file '" + path + "'"; }
    f << o.str();
    if (!f) { return "cannot write the SubDyn file '" + path + "'"; }
    return std::string();
}

} // namespace erf_towers
