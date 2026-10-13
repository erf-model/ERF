#include <ERF_FireEngineFactory.H>

// No fire engine is built into this ERF: Source/FireCoupling holds the coupling
// to the atmosphere only. A code that adds an engine lists its name here and
// creates it below.

std::vector<std::string>
fire_engine_names ()
{
    return {};
}

std::unique_ptr<FireEngine>
make_fire_engine (const std::string& /*name*/, const FireCouplingParams& /*params*/)
{
    return nullptr;
}
