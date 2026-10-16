#include <ERF_FireCouplingParams.H>
#include <ERF_FireEngineFactory.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <vector>

FireCouplingParams
FireCouplingParams::read ()
{
    FireCouplingParams c;
    amrex::ParmParse pp("erf.fire");
    // Every message names the key and the value
    auto fail = [](const std::string& msg) { amrex::Abort("[FIRE] " + msg); };
    auto one_of = [&](const char* key, const std::string& v, const std::vector<std::string>& allowed) {
        for (const auto& a : allowed) { if (v == a) { return; } }
        std::string list;
        for (const auto& a : allowed) { list += (list.empty() ? "" : ", ") + std::string("\"") + a + "\""; }
        fail(std::string("erf.fire.") + key + " = \"" + v + "\" is not one of " + list);
    };
    auto positive = [&](const char* key, amrex::Real v) {
        if (!(v > 0.0)) { fail(std::string("erf.fire.") + key + " must be > 0, got " + std::to_string(v)); }
    };
    auto nonneg = [&](const char* key, amrex::Real v) {
        if (!(v >= 0.0)) { fail(std::string("erf.fire.") + key + " must be >= 0, got " + std::to_string(v)); }
    };

    pp.query("enable", c.enable);
    pp.query("fire_debug", c.fire_debug);
    pp.query("grid_ratio", c.grid_ratio);
    if (c.grid_ratio < 1) {
        fail("erf.fire.grid_ratio must be >= 1, got " + std::to_string(c.grid_ratio));
    }
    pp.query("anchor_level", c.anchor_level);
    if (c.anchor_level < -1) {
        fail("erf.fire.anchor_level must be -1 (the finest level) or a level >= 0, got "
             + std::to_string(c.anchor_level));
    }

    // The engine: one built into this ERF, the first by default
    const std::vector<std::string> engines = fire_engine_names();
    if (!engines.empty()) { c.engine = engines.front(); }
    const bool engine_given = pp.query("engine", c.engine);
    if (engine_given || c.enable) {
        if (engines.empty()) {
            fail(std::string(c.enable ? "erf.fire.enable = true" : "erf.fire.engine is set")
                 + ", but this ERF has no fire engine built in: Source/FireCoupling holds only the "
                   "coupling to the atmosphere. Source/FireCoupling/ERF_FireEngineTemplate.H lists what "
                   "an engine provides and how to add one.");
        }
        one_of("engine", c.engine, engines);
    }

    pp.query("heat_flux_alfg", c.heat_flux_alfg);
    pp.query("fire_atm_feedback", c.fire_atm_feedback);
    pp.query("coupling_type", c.coupling_type);
    one_of("coupling_type", c.coupling_type, {"passive", "lagged", "synchronous"});
    positive("heat_flux_alfg", c.heat_flux_alfg);
    nonneg("fire_atm_feedback", c.fire_atm_feedback);   // 0 = one-way, the documented switch
    pp.query("inject_latent", c.inject_latent);
    c.source_mode_given = pp.query("source_mode", c.source_mode);
    one_of("source_mode", c.source_mode, {"overwrite", "add"});
    pp.query("heat_open_fraction", c.heat_open_fraction);
    pp.query("heat_tendency_density", c.heat_tendency_density);
    pp.query("heat_tendency_exner", c.heat_tendency_exner);
    // The keys of the heat injection, given with passive coupling, which
    // injects no heat: warn once per key
    if (c.is_passive()) {
        for (const char* k : {"heat_flux_alfg", "fire_atm_feedback", "inject_latent", "source_mode",
                              "heat_open_fraction", "heat_tendency_density", "heat_tendency_exner"}) {
            if (pp.contains(k)) {
                amrex::Print() << "[FIRE] WARNING: erf.fire." << k << " is ignored with erf.fire.coupling_type"
                               << " = passive (no heat is injected)\n";
            }
        }
    }
    pp.query("smoke_enable", c.smoke_enable);
    pp.query("smoke_emission_factor", c.smoke_emission_factor);
    nonneg("smoke_emission_factor", c.smoke_emission_factor);
    // A smoke key without the smoke switch: a deck that sets the smoke switch,
    // or the fire switch, to false explicitly is a variant switching smoke off,
    // so warn; a deck that set neither forgot the smoke switch, so stop
    if (!c.smoke_enable && pp.contains("smoke_emission_factor")) {
        if (pp.contains("smoke_enable")) {
            amrex::Print() << "[FIRE] WARNING: erf.fire.smoke_emission_factor is ignored with"
                           << " erf.fire.smoke_enable = false\n";
        } else if (pp.contains("enable") && !c.enable) {
            amrex::Print() << "[FIRE] WARNING: erf.fire.smoke_emission_factor is ignored with"
                           << " erf.fire.enable = false\n";
        } else {
            fail("erf.fire.smoke_emission_factor needs erf.fire.smoke_enable = true");
        }
    }
    return c;
}
