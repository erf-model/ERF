#include <ERF_SurfaceCouplingParams.H>
#include <AMReX.H>

std::string
surface_quoted_list (const std::vector<std::string>& v)
{
    std::string list;
    for (const auto& a : v) { list += (list.empty() ? "" : ", ") + std::string("\"") + a + "\""; }
    return list;
}

void
read_surface_coupling_common (amrex::ParmParse& pp,
                              const SurfaceCouplingNames& names,
                              const std::vector<std::string>& engines,
                              bool& enable,
                              int& grid_ratio,
                              int& anchor_level,
                              std::string& engine)
{
    // Every message names the key and the value
    auto fail = [&](const std::string& msg) { amrex::Abort(names.tag() + " " + msg); };

    pp.query("enable", enable);
    pp.query("grid_ratio", grid_ratio);
    if (grid_ratio < 1) {
        fail(names.key("grid_ratio") + " must be >= 1, got " + std::to_string(grid_ratio));
    }
    pp.query("anchor_level", anchor_level);
    if (anchor_level < -1) {
        fail(names.key("anchor_level") + " must be -1 (the finest level) or a level >= 0, got "
             + std::to_string(anchor_level));
    }

    // The engine: one compiled into this executable (ERF_REGISTER_<X>_ENGINE).
    // With a single engine it is the default; with several the deck names one.
    if (engines.size() == 1) { engine = engines.front(); }
    const bool engine_given = pp.query("engine", engine);
    if (engine_given || enable) {
        const std::string dir = std::string("Source/") + names.title + "Coupling";
        if (engines.empty()) {
            fail((enable ? names.key("enable") + " = true" : names.key("engine") + " is set")
                 + ", but no " + names.lower + " engine is compiled into this executable: " + dir
                 + " holds only the coupling to the atmosphere. Build with an engine directory (GNU make "
                 + names.upper + "_HOME, CMake ERF_" + names.upper + "_ENGINE_DIRS); " + dir + "/ERF_"
                 + names.title + "EngineTemplate.H lists what an engine provides.");
        }
        if (!engine_given && engines.size() > 1) {
            fail(std::string("several ") + names.lower + " engines are compiled in ("
                 + surface_quoted_list(engines) + "): set " + names.key("engine") + " to one of them");
        }
        bool known = false;
        for (const auto& a : engines) { known = known || (a == engine); }
        if (!known) {
            fail(names.key("engine") + " = \"" + engine + "\" is not one of " + surface_quoted_list(engines));
        }
    }
}
