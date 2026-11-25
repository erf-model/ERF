#include <ERF_FireEngineRegistry.H>
#include <ERF_FireGrid.H>
#include <ERF_SurfaceEngineRegistry.H>

// The fire engines live in the surface-coupling registry for FireEngine
using FireRegistry = SurfaceEngineRegistry<FireEngine, FireCouplingParams>;

bool
register_fire_engine (const std::string& name, FireEngineCreator creator)
{
    return FireRegistry::add(fire_coupling_names(), name, std::move(creator));
}

std::vector<std::string>
fire_engine_names ()
{
    return FireRegistry::names();
}

std::unique_ptr<FireEngine>
make_fire_engine (const std::string& name, const FireCouplingParams& params)
{
    return FireRegistry::make(name, params);
}
