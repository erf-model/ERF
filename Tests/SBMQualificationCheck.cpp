#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace {

bool read_real(const std::map<std::string, std::string>& values,
               const std::string& key, double& result)
{
    const auto it = values.find(key);
    if (it == values.end()) return false;
    try {
        std::size_t consumed = 0;
        result = std::stod(it->second, &consumed);
        return consumed == it->second.size() && std::isfinite(result);
    } catch (...) {
        return false;
    }
}

bool read_integer(const std::map<std::string, std::string>& values,
                  const std::string& key, int& result)
{
    const auto it = values.find(key);
    if (it == values.end()) return false;
    try {
        std::size_t consumed = 0;
        result = std::stoi(it->second, &consumed);
        return consumed == it->second.size();
    } catch (...) {
        return false;
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: erf_sbm_qualification_check DIAGNOSTIC\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) {
        std::cerr << "cannot open SBM diagnostic: " << argv[1] << '\n';
        return 2;
    }

    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(input, line)) {
        const auto separator = line.find('=');
        if (separator == std::string::npos || separator == 0) {
            std::cerr << "malformed diagnostic line: " << line << '\n';
            return 1;
        }
        values[line.substr(0, separator)] = line.substr(separator + 1);
    }

    const auto format = values.find("format");
    const auto method = values.find("method");
    if (format == values.end() || format->second != "erf-sbm-p1-diagnostic-v1" ||
        method == values.end() ||
        (method->second != "compressible" && method->second != "anelastic")) {
        std::cerr << "invalid SBM diagnostic identity\n";
        return 1;
    }

    int nbins = 0;
    int passed = 0;
    if (!read_integer(values, "nbins", nbins) ||
        (nbins != 4 && nbins != 16 && nbins != 64) ||
        !read_integer(values, "passed", passed) || passed != 1) {
        std::cerr << "SBM diagnostic failed its discrete checks\n";
        return 1;
    }

    double initial_mass = 0.0;
    double final_mass = 0.0;
    double mass_error = 0.0;
    double mass_tolerance = 0.0;
    double compact_mass_error = 0.0;
    double initial_variation = 0.0;
    double transport_change = 0.0;
    double projection_error = 0.0;
    double projection_tolerance = 0.0;
    double face_projection_error = 0.0;
    double face_tolerance = 0.0;
    const bool finite =
        read_real(values, "initial_mass", initial_mass) &&
        read_real(values, "final_mass", final_mass) &&
        read_real(values, "mass_error", mass_error) &&
        read_real(values, "mass_tolerance", mass_tolerance) &&
        read_real(values, "compact_mass_error", compact_mass_error) &&
        read_real(values, "initial_variation", initial_variation) &&
        read_real(values, "transport_change", transport_change) &&
        read_real(values, "projection_error", projection_error) &&
        read_real(values, "projection_tolerance", projection_tolerance) &&
        read_real(values, "face_projection_error", face_projection_error) &&
        read_real(values, "face_tolerance", face_tolerance);
    if (!finite || !(mass_tolerance > 0.0) || !(projection_tolerance > 0.0) ||
        !(face_tolerance > 0.0) || !(initial_variation > 0.0) ||
        !(transport_change > 0.0) || mass_error > mass_tolerance ||
        compact_mass_error > mass_tolerance ||
        projection_error > projection_tolerance ||
        face_projection_error > face_tolerance) {
        std::cerr << "SBM numerical invariant check failed\n";
        return 1;
    }

    std::cout << "SBM P1 qualification passed: method=" << method->second
              << " nbins=" << nbins << " mass_error=" << mass_error
              << " projection_error=" << projection_error
              << " face_projection_error=" << face_projection_error << '\n';
    return 0;
}
