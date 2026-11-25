#include <ERF_FloodEngineRegistry.H>
#include <ERF_FloodGrid.H>
#include <ERF_SurfaceEngineRegistry.H>

// The flood engines live in the surface-coupling registry for FloodEngine
using FloodRegistry = SurfaceEngineRegistry<FloodEngine, FloodCouplingParams>;

bool
register_flood_engine (const std::string& name, FloodEngineCreator creator)
{
    return FloodRegistry::add(flood_coupling_names(), name, std::move(creator));
}

std::vector<std::string>
flood_engine_names ()
{
    return FloodRegistry::names();
}

std::unique_ptr<FloodEngine>
make_flood_engine (const std::string& name, const FloodCouplingParams& params)
{
    return FloodRegistry::make(name, params);
}
