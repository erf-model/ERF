#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Station time-series checker.  Reads the series ERF writes for a station
// (Output_Stations/<name>.dat: comment lines starting with '#', then one row
// per sample, the time first) and checks columns of it.
//
//   equal a=F:C b=G:D tol=E
//       Column C of F and column D of G must agree row by row to within tol,
//       at the same times, and must not be constant, so that two series that
//       should sample the same point are compared on something that moves.
//
//   analytic file=F col=C phi0=P a=A b=B tau=T tol=E
//       Every row must match the solution of d(phi)/dt = -(phi - (a + b t))/tau
//       with phi(0) = phi0,
//           phi(t) = a + b t - b tau + (phi0 - a + b tau) exp(-t / tau),
//       to within tol, and the series must move away from phi0 by more than
//       100 tol, so that a run in which nothing is nudged cannot pass.
//
//   approach on=F off=G col=C target=X factor=R
//       At the last row, the run with nudging (on) must be closer to target
//       than factor times the distance of the run without it (off), and the
//       two runs must differ, so that the check cannot pass on a no-op.
//
// Columns count from 1, as the header of the file does, and column 1 is time.

namespace {

int fail (const std::string& message)
{
    std::cerr << "StationSeriesCheck: FAIL: " << message << "\n";
    return 1;
}

bool read_series (const std::string& file, std::vector<std::vector<double>>& rows, std::string& error)
{
    std::ifstream is(file);
    if (!is) { error = "cannot open " + file; return false; }
    std::string line;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream ls(line);
        std::vector<double> row;
        double v;
        while (ls >> v) { row.push_back(v); }
        if (!row.empty()) { rows.push_back(row); }
    }
    if (rows.empty()) { error = file + " has no data rows"; return false; }
    return true;
}

bool get (const std::map<std::string, std::string>& args, const std::string& key, std::string& value)
{
    auto it = args.find(key);
    if (it == args.end()) { return false; }
    value = it->second;
    return true;
}

bool get (const std::map<std::string, std::string>& args, const std::string& key, double& value)
{
    std::string s;
    if (!get(args, key, s)) { return false; }
    char* end = nullptr;
    value = std::strtod(s.c_str(), &end);
    return end != s.c_str() && *end == '\0';
}

bool split_file_col (const std::string& spec, std::string& file, std::size_t& col)
{
    const auto colon = spec.rfind(':');
    if (colon == std::string::npos) { return false; }
    file = spec.substr(0, colon);
    col  = static_cast<std::size_t>(std::strtol(spec.substr(colon+1).c_str(), nullptr, 10));
    return col >= 2;
}

int check_equal (const std::map<std::string, std::string>& args)
{
    std::string a, b, fa, fb;
    double tol = 0;
    std::size_t ca = 0, cb = 0;
    if (!get(args, "a", a) || !get(args, "b", b) || !get(args, "tol", tol) ||
        !split_file_col(a, fa, ca) || !split_file_col(b, fb, cb)) {
        return fail("equal needs a=FILE:COL b=FILE:COL tol= with COL >= 2");
    }
    std::vector<std::vector<double>> ra, rb;
    std::string error;
    if (!read_series(fa, ra, error)) { return fail(error); }
    if (!read_series(fb, rb, error)) { return fail(error); }
    if (ra.size() != rb.size()) { return fail("the two series have different numbers of rows"); }

    double max_diff = 0.0, lo = 1.0e300, hi = -1.0e300;
    for (std::size_t n = 0; n < ra.size(); ++n) {
        if (ra[n].size() < ca || rb[n].size() < cb) { return fail("a row is missing the column"); }
        if (ra[n][0] != rb[n][0]) { return fail("the two series are at different times"); }
        const double va = ra[n][ca-1], vb = rb[n][cb-1];
        max_diff = std::max(max_diff, std::abs(va - vb));
        lo = std::min(lo, va); hi = std::max(hi, va);
    }
    std::cout << "StationSeriesCheck: equal " << a << " vs " << b << ": " << ra.size()
              << " rows, max |difference| = " << max_diff << ", range of the first " << hi - lo << "\n";
    if (hi - lo <= tol) { return fail("the series is constant, so the comparison is vacuous"); }
    if (max_diff > tol) { return fail("the series differ by more than tol"); }
    std::cout << "StationSeriesCheck: PASS\n";
    return 0;
}

int check_analytic (const std::map<std::string, std::string>& args)
{
    std::string file;
    double col = 0, phi0 = 0, a = 0, b = 0, tau = 0, tol = 0;
    if (!get(args, "file", file) || !get(args, "col", col) || !get(args, "phi0", phi0) ||
        !get(args, "a", a) || !get(args, "b", b) || !get(args, "tau", tau) || !get(args, "tol", tol)) {
        return fail("analytic needs file= col= phi0= a= b= tau= tol=");
    }
    std::vector<std::vector<double>> rows;
    std::string error;
    if (!read_series(file, rows, error)) { return fail(error); }

    const auto c = static_cast<std::size_t>(col) - 1;
    double max_err = 0.0, max_departure = 0.0;
    for (const auto& row : rows) {
        if (row.size() <= c) { return fail(file + " has no column " + std::to_string(c+1)); }
        const double t = row[0];
        const double exact = a + b*t - b*tau + (phi0 - a + b*tau) * std::exp(-t/tau);
        max_err = std::max(max_err, std::abs(row[c] - exact));
        max_departure = std::max(max_departure, std::abs(row[c] - phi0));
    }
    std::cout << "StationSeriesCheck: analytic " << file << " column " << c+1 << ": " << rows.size()
              << " rows, last t = " << rows.back()[0] << ", max |error| = " << max_err
              << ", max departure from phi0 = " << max_departure << "\n";
    if (max_departure <= 100.0 * tol) {
        return fail("the series hardly moves from its initial value, so the comparison is vacuous");
    }
    if (max_err > tol) {
        return fail("max |error| " + std::to_string(max_err) + " exceeds tol " + std::to_string(tol));
    }
    std::cout << "StationSeriesCheck: PASS\n";
    return 0;
}

int check_approach (const std::map<std::string, std::string>& args)
{
    std::string on, off;
    double col = 0, target = 0, factor = 0;
    if (!get(args, "on", on) || !get(args, "off", off) || !get(args, "col", col) ||
        !get(args, "target", target) || !get(args, "factor", factor)) {
        return fail("approach needs on= off= col= target= factor=");
    }
    std::vector<std::vector<double>> rows_on, rows_off;
    std::string error;
    if (!read_series(on, rows_on, error))   { return fail(error); }
    if (!read_series(off, rows_off, error)) { return fail(error); }

    const auto c = static_cast<std::size_t>(col) - 1;
    if (rows_on.back().size() <= c || rows_off.back().size() <= c) {
        return fail("no column " + std::to_string(c+1));
    }
    if (rows_on.back()[0] != rows_off.back()[0]) {
        return fail("the two runs end at different times");
    }
    const double v_on  = rows_on.back()[c];
    const double v_off = rows_off.back()[c];
    const double d_on  = std::abs(v_on - target);
    const double d_off = std::abs(v_off - target);
    std::cout << "StationSeriesCheck: approach column " << c+1 << " at t = " << rows_on.back()[0]
              << ": target " << target << ", nudged " << v_on << " (|d| " << d_on << "), free "
              << v_off << " (|d| " << d_off << "), required ratio " << factor << "\n";
    if (v_on == v_off) {
        return fail("the nudged and free runs agree exactly, so nothing was nudged");
    }
    if (!(d_on < factor * d_off)) {
        return fail("the nudged run is not closer to the target by the required factor");
    }
    std::cout << "StationSeriesCheck: PASS\n";
    return 0;
}

} // namespace

int main (int argc, char* argv[])
{
    if (argc < 2) { return fail("usage: StationSeriesCheck equal|analytic|approach key=value ..."); }
    const std::string mode = argv[1];
    std::map<std::string, std::string> args;
    for (int i = 2; i < argc; ++i) {
        const std::string s = argv[i];
        const auto eq = s.find('=');
        if (eq == std::string::npos) { return fail("argument '" + s + "' is not key=value"); }
        args[s.substr(0, eq)] = s.substr(eq+1);
    }
    if (mode == "equal")    { return check_equal(args); }
    if (mode == "analytic") { return check_analytic(args); }
    if (mode == "approach") { return check_approach(args); }
    return fail("unknown mode '" + mode + "'");
}
