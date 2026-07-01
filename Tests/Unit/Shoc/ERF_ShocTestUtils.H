#ifndef ERF_SHOC_TEST_UTILS_H_
#define ERF_SHOC_TEST_UTILS_H_

#include "ERF_ShocColumnData.H"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace shoc_test
{
inline amrex::Arena*
test_arena ()
{
#ifdef AMREX_USE_GPU
    return amrex::The_Managed_Arena();
#else
    return amrex::The_Arena();
#endif
}

inline void
sync ()
{
    shoc::sync_if_needed();
}

// Run a SHOC test action and synchronize before host-side EXPECT/ASSERT reads.
template <typename F>
inline void
run_and_sync (F&& f)
{
    std::forward<F>(f)();
    sync();
}

inline std::string
fixture_dir ()
{
#ifdef ERF_SHOC_UNIT_FIXTURE_DIR
    return ERF_SHOC_UNIT_FIXTURE_DIR;
#else
    throw std::runtime_error("ERF_SHOC_UNIT_FIXTURE_DIR is not defined");
#endif
}

inline amrex::Vector<amrex::Real>
read_fixture_vector (const std::string& relpath)
{
    std::ifstream ifs(fixture_dir() + "/" + relpath);
    if (!ifs.is_open()) {
        throw std::runtime_error("Could not open SHOC fixture: " + relpath);
    }

    amrex::Vector<amrex::Real> values;
    std::string line;
    while (std::getline(ifs, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) {
            line.erase(hash);
        }
        std::istringstream iss(line);
        amrex::Real value = 0.0;
        while (iss >> value) {
            values.push_back(value);
        }
    }

    if (values.empty()) {
        throw std::runtime_error("Empty SHOC fixture: " + relpath);
    }
    return values;
}

inline std::map<std::string, amrex::Vector<amrex::Real>>
read_named_fixture_vectors (const std::string& relpath)
{
    std::ifstream ifs(fixture_dir() + "/" + relpath);
    if (!ifs.is_open()) {
        throw std::runtime_error("Could not open SHOC fixture: " + relpath);
    }

    std::map<std::string, amrex::Vector<amrex::Real>> values;
    std::string line;
    while (std::getline(ifs, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) {
            line.erase(hash);
        }

        const auto colon = line.find(':');
        if (colon == std::string::npos) {
            continue;
        }

        std::string key = line.substr(0, colon);
        key.erase(std::remove_if(key.begin(), key.end(), [](unsigned char ch) {
            return std::isspace(ch);
        }), key.end());
        if (key.empty()) {
            continue;
        }

        std::istringstream iss(line.substr(colon + 1));
        amrex::Real value = 0.0;
        while (iss >> value) {
            values[key].push_back(value);
        }
    }

    if (values.empty()) {
        throw std::runtime_error("Empty named SHOC fixture: " + relpath);
    }
    return values;
}

inline const amrex::Vector<amrex::Real>&
fixture_values (const std::map<std::string, amrex::Vector<amrex::Real>>& fixture,
                const std::string& key)
{
    const auto it = fixture.find(key);
    if (it == fixture.end()) {
        throw std::runtime_error("Missing SHOC fixture key: " + key);
    }
    return it->second;
}

inline ShocColumnData
make_column (int nlev)
{
    ShocColumnData col;
    ShocColumnLayout layout;
    layout.nx = 1;
    layout.ny = 1;
    layout.ncell = 1;
    layout.nlev = nlev;
    define_shoc_column_data(col, layout, test_arena(), shoc::InitRunOn::Host);
    sync();

    auto zt = col.zt.array();
    auto zi = col.zi.array();
    auto dz = col.dz.array();
    auto p_mid = col.p_mid.array();
    auto p_int = col.p_int.array();
    auto rho = col.rho.array();
    auto theta = col.theta.array();
    auto exner = col.exner.array();
    auto theta_v = col.theta_v.array();
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qc = col.qc.array();
    auto qi = col.qi.array();
    auto qw = col.qw.array();
    auto tabs = col.tabs.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();
    auto w = col.w.array();
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();

    for (int k = 0; k <= nlev; ++k) {
        zi(0,k,0) = 100.0 * k;
        p_int(0,k,0) = 100000.0 - 1000.0 * k;
    }

    for (int k = 0; k < nlev; ++k) {
        zt(0,k,0) = 50.0 + 100.0 * k;
        dz(0,k,0) = 100.0;
        p_mid(0,k,0) = 99500.0 - 1000.0 * k;
        rho(0,k,0) = 1.15 - 0.05 * k;
        theta(0,k,0) = 300.0 + 0.5 * k;
        exner(0,k,0) = 1.0;
        theta_v(0,k,0) = 301.0 + 0.8 * k;
        thetal(0,k,0) = 299.0 + 0.4 * k;
        qv(0,k,0) = 0.010 - 0.001 * k;
        qc(0,k,0) = 0.0;
        qi(0,k,0) = 0.0;
        qw(0,k,0) = qv(0,k,0);
        tabs(0,k,0) = 295.0 - 1.0 * k;
        tke(0,k,0) = 0.3 + 0.05 * k;
        u(0,k,0) = 2.0 + 0.7 * k;
        v(0,k,0) = 1.0 + 0.25 * k;
        w(0,k,0) = 0.0;
        tk(0,k,0) = 0.0;
        tkh(0,k,0) = 0.0;
    }

    shoc::set_fab_val(col.surf_sens_flux, 0.02, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, 1.0e-4, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, 0.0, shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, 0.0, shoc::InitRunOn::Host);

    return col;
}

inline ShocColumnData
make_randomized_column (int nlev, std::mt19937& gen)
{
    std::uniform_real_distribution<amrex::Real> theta_dist(285.0, 315.0);
    std::uniform_real_distribution<amrex::Real> qv_dist(1.0e-4, 0.018);
    std::uniform_real_distribution<amrex::Real> tke_dist(4.0e-4, 2.0);
    std::uniform_real_distribution<amrex::Real> vel_dist(-4.0, 8.0);
    std::uniform_real_distribution<amrex::Real> flux_dist(-0.03, 0.05);

    auto col = make_column(nlev);
    auto theta = col.theta.array();
    auto theta_v = col.theta_v.array();
    auto thetal = col.thetal.array();
    auto qv = col.qv.array();
    auto qw = col.qw.array();
    auto qi = col.qi.array();
    auto tabs = col.tabs.array();
    auto tke = col.tke.array();
    auto u = col.u.array();
    auto v = col.v.array();

    amrex::Real running_theta = theta_dist(gen);
    for (int k = 0; k < nlev; ++k) {
        running_theta += std::uniform_real_distribution<amrex::Real>(-0.2, 1.5)(gen);
        theta(0,k,0) = running_theta;
        theta_v(0,k,0) = running_theta + std::uniform_real_distribution<amrex::Real>(0.2, 2.0)(gen);
        thetal(0,k,0) = running_theta - std::uniform_real_distribution<amrex::Real>(0.0, 1.0)(gen);
        qv(0,k,0) = qv_dist(gen);
        qi(0,k,0) = 0.0;
        qw(0,k,0) = qv(0,k,0);
        tabs(0,k,0) = running_theta - std::uniform_real_distribution<amrex::Real>(2.0, 8.0)(gen);
        tke(0,k,0) = tke_dist(gen);
        u(0,k,0) = vel_dist(gen);
        v(0,k,0) = vel_dist(gen);
    }

    shoc::set_fab_val(col.surf_sens_flux, flux_dist(gen), shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_lat_flux, std::uniform_real_distribution<amrex::Real>(-5.0e-4, 5.0e-4)(gen), shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_u, std::uniform_real_distribution<amrex::Real>(-0.1, 0.1)(gen), shoc::InitRunOn::Host);
    shoc::set_fab_val(col.surf_tau_v, std::uniform_real_distribution<amrex::Real>(-0.1, 0.1)(gen), shoc::InitRunOn::Host);

    return col;
}
}

#endif
