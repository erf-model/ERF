#include <ERF_FireEngineRegistry.H>
#include <AMReX.H>

#include <algorithm>
#include <map>

namespace {
// A function-local static, so it exists before the first registration whatever
// the order in which the engines' translation units are initialised
std::map<std::string, FireEngineCreator>&
fire_engine_registry ()
{
    static std::map<std::string, FireEngineCreator> registry;
    return registry;
}
}

bool
register_fire_engine (const std::string& name, FireEngineCreator creator)
{
    auto& registry = fire_engine_registry();
    if (name.empty()) {
        amrex::Abort("[FIRE] ERF_REGISTER_FIRE_ENGINE: an engine registered an empty name");
    }
    if (registry.count(name) > 0) {
        amrex::Abort("[FIRE] ERF_REGISTER_FIRE_ENGINE: two engines registered the name \"" + name
                     + "\"; each engine directory in the build must use its own name");
    }
    registry.emplace(name, std::move(creator));
    return true;
}

std::vector<std::string>
fire_engine_names ()
{
    std::vector<std::string> names;
    for (const auto& kv : fire_engine_registry()) { names.push_back(kv.first); }
    return names;   // std::map keeps them sorted
}

std::unique_ptr<FireEngine>
make_fire_engine (const std::string& name, const FireCouplingParams& params)
{
    const auto& registry = fire_engine_registry();
    const auto it = registry.find(name);
    if (it == registry.end()) { return nullptr; }
    return it->second(params);
}
