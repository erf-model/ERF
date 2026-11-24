#include <ERF_FloodEngineRegistry.H>
#include <AMReX.H>

#include <map>

namespace {
// A function-local static, so it exists before the first registration whatever
// the order in which the engines' translation units are initialized
std::map<std::string, FloodEngineCreator>&
flood_engine_registry ()
{
    static std::map<std::string, FloodEngineCreator> registry;
    return registry;
}
}

bool
register_flood_engine (const std::string& name, FloodEngineCreator creator)
{
    auto& registry = flood_engine_registry();
    if (name.empty()) {
        amrex::Abort("[FLOOD] ERF_REGISTER_FLOOD_ENGINE: an engine registered an empty name");
    }
    if (registry.count(name) > 0) {
        amrex::Abort("[FLOOD] ERF_REGISTER_FLOOD_ENGINE: two engines registered the name \"" + name
                     + "\"; each engine directory in the build must use its own name");
    }
    registry.emplace(name, std::move(creator));
    return true;
}

std::vector<std::string>
flood_engine_names ()
{
    std::vector<std::string> names;
    for (const auto& kv : flood_engine_registry()) { names.push_back(kv.first); }
    return names;   // std::map keeps them sorted
}

std::unique_ptr<FloodEngine>
make_flood_engine (const std::string& name, const FloodCouplingParams& params)
{
    const auto& registry = flood_engine_registry();
    const auto it = registry.find(name);
    if (it == registry.end()) { return nullptr; }
    return it->second(params);
}
