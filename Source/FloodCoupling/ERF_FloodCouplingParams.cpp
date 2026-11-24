#include <ERF_FloodCouplingParams.H>
#include <ERF_FloodEngineRegistry.H>
#include <AMReX_ParmParse.H>

#include <string>
#include <vector>

FloodCouplingParams
FloodCouplingParams::read ()
{
    FloodCouplingParams c;
    amrex::ParmParse pp("erf.flood");
    // Every message names the key and the value
    auto fail = [](const std::string& msg) { amrex::Abort("[FLOOD] " + msg); };
    auto quoted_list = [](const std::vector<std::string>& v) {
        std::string list;
        for (const auto& a : v) { list += (list.empty() ? "" : ", ") + std::string("\"") + a + "\""; }
        return list;
    };

    pp.query("enable", c.enable);
    pp.query("grid_ratio", c.grid_ratio);
    if (c.grid_ratio < 1) {
        fail("erf.flood.grid_ratio must be >= 1, got " + std::to_string(c.grid_ratio));
    }
    pp.query("anchor_level", c.anchor_level);
    if (c.anchor_level < -1) {
        fail("erf.flood.anchor_level must be -1 (the finest level) or a level >= 0, got "
             + std::to_string(c.anchor_level));
    }

    // The engine: one compiled into this executable (ERF_REGISTER_FLOOD_ENGINE).
    // With a single engine it is the default; with several the deck names one.
    const std::vector<std::string> engines = flood_engine_names();
    if (engines.size() == 1) { c.engine = engines.front(); }
    const bool engine_given = pp.query("engine", c.engine);
    if (engine_given || c.enable) {
        if (engines.empty()) {
            fail(std::string(c.enable ? "erf.flood.enable = true" : "erf.flood.engine is set")
                 + ", but no flood engine is compiled into this executable: Source/FloodCoupling holds "
                   "only the coupling to the atmosphere. Build with an engine directory (GNU make "
                   "FLOOD_HOME, CMake ERF_FLOOD_ENGINE_DIRS); Source/FloodCoupling/ERF_FloodEngineTemplate.H "
                   "lists what an engine provides.");
        }
        if (!engine_given && engines.size() > 1) {
            fail("several flood engines are compiled in (" + quoted_list(engines)
                 + "): set erf.flood.engine to one of them");
        }
        bool known = false;
        for (const auto& a : engines) { known = known || (a == c.engine); }
        if (!known) {
            fail("erf.flood.engine = \"" + c.engine + "\" is not one of " + quoted_list(engines));
        }
    }
    return c;
}
