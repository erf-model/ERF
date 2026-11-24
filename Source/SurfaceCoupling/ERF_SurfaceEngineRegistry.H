#ifndef ERF_SURFACE_ENGINE_REGISTRY_H
#define ERF_SURFACE_ENGINE_REGISTRY_H

#include <ERF_SurfaceCouplingNames.H>
#include <AMReX.H>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

/**
 * @file ERF_SurfaceEngineRegistry.H
 * @brief The registry behind ERF_REGISTER_FIRE_ENGINE and ERF_REGISTER_FLOOD_ENGINE
 *
 * One registry per engine base class: SurfaceEngineRegistry<FireEngine,
 * FireCouplingParams> holds the fire engines, SurfaceEngineRegistry<FloodEngine,
 * FloodCouplingParams> the flood engines. Each coupling instantiates its
 * registry in one .cpp (ERF_FireEngineRegistry.cpp, ERF_FloodEngineRegistry.cpp)
 * behind its own functions (register_fire_engine, fire_engine_names,
 * make_fire_engine, and the flood ones), so the map lives in that
 * translation unit and the registration macros are unchanged.
 */
template <class Engine, class Params>
class SurfaceEngineRegistry
{
public:
    /// Creates an engine from the coupling inputs.
    using Creator = std::function<std::unique_ptr<Engine>(const Params&)>;

    /// Add an engine; aborts on an empty name or a name registered twice.
    static bool add (const SurfaceCouplingNames& names, const std::string& name, Creator creator)
    {
        auto& registry = map();
        const std::string macro = std::string("ERF_REGISTER_") + names.upper + "_ENGINE";
        if (name.empty()) {
            amrex::Abort(names.tag() + " " + macro + ": an engine registered an empty name");
        }
        if (registry.count(name) > 0) {
            amrex::Abort(names.tag() + " " + macro + ": two engines registered the name \"" + name
                         + "\"; each engine directory in the build must use its own name");
        }
        registry.emplace(name, std::move(creator));
        return true;
    }

    /// The registered names, sorted.
    static std::vector<std::string> names ()
    {
        std::vector<std::string> out;
        for (const auto& kv : map()) { out.push_back(kv.first); }
        return out;   // std::map keeps them sorted
    }

    /// The engine registered under name, or nullptr.
    static std::unique_ptr<Engine> make (const std::string& name, const Params& params)
    {
        const auto& registry = map();
        const auto it = registry.find(name);
        if (it == registry.end()) { return nullptr; }
        return it->second(params);
    }

private:
    // A function-local static, so it exists before the first registration
    // whatever the order in which the engines' translation units are initialized
    static std::map<std::string, Creator>& map ()
    {
        static std::map<std::string, Creator> registry;
        return registry;
    }
};

#endif
