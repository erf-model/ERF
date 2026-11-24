#ifndef ERF_FLOOD_ENGINE_REGISTRY_H
#define ERF_FLOOD_ENGINE_REGISTRY_H

#include <ERF_FloodCouplingParams.H>
#include <ERF_FloodEngine.H>
#include <AMReX.H>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

/**
 * @file ERF_FloodEngineRegistry.H
 * @brief The flood engines compiled into this executable
 *
 * ERF lists no engine itself. Each engine registers itself from its own .cpp
 * with
 *
 *     ERF_REGISTER_FLOOD_ENGINE(MyEngine, "myengine")
 *
 * at namespace scope, and erf.flood.engine then selects it by that name. The
 * engine's directory (FLOOD_HOME in GNU make, ERF_FLOOD_ENGINE_DIRS in CMake)
 * adds that .cpp to the build, so adding an engine changes no file of ERF.
 *
 * The registration runs during static initialization, so the .cpp that holds
 * it must be linked as an object file: inside a static archive (.a) that
 * nothing else references, the linker drops it and the engine is missing.
 * ERF's CMake build compiles engines into its OBJECT library, and GNU make
 * links every object, so both keep it; an engine that ships a prebuilt
 * library keeps its registering .cpp as source.
 *
 * The registry follows the fire-engine registry (ERF_FireEngineRegistry.H)
 * and is independent of it: a build may hold fire engines, flood engines,
 * both or neither.
 */

/// Creates an engine from the coupling inputs.
using FloodEngineCreator = std::function<std::unique_ptr<FloodEngine>(const FloodCouplingParams&)>;

/**
 * @brief Add an engine to the registry; aborts on a name registered twice
 *
 * Called by ERF_REGISTER_FLOOD_ENGINE; returns true so the macro can hold the
 * result in a static variable.
 */
bool register_flood_engine (const std::string& name, FloodEngineCreator creator);

/// Names erf.flood.engine accepts, sorted; empty when no engine is compiled in.
std::vector<std::string> flood_engine_names ();

/// The engine erf.flood.engine names, or nullptr when no engine has that name.
std::unique_ptr<FloodEngine> make_flood_engine (const std::string& name, const FloodCouplingParams& params);

#define ERF_FLOOD_ENGINE_CONCAT_(a, b) a##b
#define ERF_FLOOD_ENGINE_CONCAT(a, b) ERF_FLOOD_ENGINE_CONCAT_(a, b)

/**
 * @brief Register engine class Class under name: erf.flood.engine = name
 *        selects it. Class is constructed from the coupling inputs,
 *        Class(const FloodCouplingParams&), or default-constructed when it has
 *        no such constructor.
 */
#define ERF_REGISTER_FLOOD_ENGINE(Class, name)                                        \
    namespace {                                                                       \
    [[maybe_unused]] const bool ERF_FLOOD_ENGINE_CONCAT(erf_flood_engine_registered_, __LINE__) = \
        register_flood_engine((name), [](const FloodCouplingParams& p) -> std::unique_ptr<FloodEngine> { \
            if constexpr (std::is_constructible_v<Class, const FloodCouplingParams&>) { \
                return std::make_unique<Class>(p);                                     \
            } else {                                                                  \
                amrex::ignore_unused(p);                                              \
                return std::make_unique<Class>();                                     \
            }                                                                         \
        });                                                                           \
    }

#endif
