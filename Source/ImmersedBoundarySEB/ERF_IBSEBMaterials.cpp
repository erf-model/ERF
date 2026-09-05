/**
 * \file ERF_IBSEBMaterials.cpp
 * \brief CSV reader of the face material library (rank-0 read, broadcast).
 */
#include "ERF_IBSEBMaterials.H"
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <sstream>
#include <set>
#include <cstring>
#include <algorithm>

using namespace amrex;

namespace {
std::string strip (std::string s)
{
    s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }), s.end());
    return s;
}
}

Vector<IBSEBMaterial>
read_ibseb_materials (const std::string& file)
{
    static_assert(std::is_trivially_copyable<IBSEBMaterial>::value, "IBSEBMaterial must be broadcastable");
    Vector<IBSEBMaterial> table;
    int n = 0;
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream in(file);
        if (!in.good()) { Abort("erf.ibseb.material_file: cannot open '" + file + "'"); }
        std::string line;
        if (!std::getline(in, line)) { Abort("erf.ibseb.material_file: empty file '" + file + "'"); }
        if (line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF) { line.erase(0, 3); }
        const std::string expected = "mat_id,name,albedo,emissivity,k_therm_W_per_mK,rho_cp_J_per_m3K,thickness_m,description";
        if (strip(line) != expected) {
            Abort("erf.ibseb.material_file: header must be '" + expected + "', got '" + line + "'");
        }
        std::set<int> seen;
        while (std::getline(in, line)) {
            if (strip(line).empty() || line[0] == '#') { continue; }
            std::stringstream ss(line);
            std::string f;
            IBSEBMaterial m;
            std::getline(ss, f, ','); m.mat_id = std::stoi(strip(f));
            std::getline(ss, f, ','); f = strip(f); std::strncpy(m.name, f.c_str(), sizeof(m.name) - 1); m.name[sizeof(m.name) - 1] = '\0';
            std::getline(ss, f, ','); m.albedo     = std::stod(strip(f));
            std::getline(ss, f, ','); m.emissivity = std::stod(strip(f));
            std::getline(ss, f, ','); m.k_therm    = std::stod(strip(f));
            std::getline(ss, f, ','); m.rho_cp     = std::stod(strip(f));
            std::getline(ss, f, ','); m.thickness  = std::stod(strip(f));
            if (!seen.insert(m.mat_id).second) { Abort("erf.ibseb.material_file: duplicate mat_id " + std::to_string(m.mat_id)); }
            if (m.albedo < 0.0 || m.albedo > 1.0 || m.emissivity <= 0.0 || m.emissivity > 1.0 ||
                m.k_therm <= 0.0 || m.rho_cp <= 0.0 || m.thickness <= 0.0) {
                Abort("erf.ibseb.material_file: material " + std::to_string(m.mat_id) + " has a value outside its range");
            }
            table.push_back(m);
        }
        n = static_cast<int>(table.size());
        if (n == 0) { Abort("erf.ibseb.material_file: no materials in '" + file + "'"); }
    }
    ParallelDescriptor::Bcast(&n, 1, ParallelDescriptor::IOProcessorNumber());
    table.resize(n);
    ParallelDescriptor::Bcast(reinterpret_cast<char*>(table.data()), n * sizeof(IBSEBMaterial),
                              ParallelDescriptor::IOProcessorNumber());
    return table;
}
