#ifndef ERF_FIRE_ENGINE_FACTORY_H
#define ERF_FIRE_ENGINE_FACTORY_H

#include <ERF_FireCouplingParams.H>
#include <ERF_FireEngine.H>

#include <memory>
#include <string>
#include <vector>

/**
 * @file ERF_FireEngineFactory.H
 * @brief The fire engines built into this ERF
 *
 * The one place an engine is registered: a code that adds an engine (see
 * ERF_FireEngineTemplate.H) adds its name to fire_engine_names() and a branch
 * creating it to make_fire_engine(), both in ERF_FireEngineFactory.cpp.
 */

/// Names erf.fire.engine accepts, the default first; empty when no engine is built in.
std::vector<std::string> fire_engine_names ();

/// The engine erf.fire.engine names, or nullptr when no engine has that name.
std::unique_ptr<FireEngine> make_fire_engine (const std::string& name, const FireCouplingParams& params);

#endif
