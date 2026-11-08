#ifndef ERF_FIRE_ENGINE_REGISTRY_H
#define ERF_FIRE_ENGINE_REGISTRY_H

#include <ERF_FireCouplingParams.H>
#include <ERF_FireEngine.H>
#include <AMReX.H>

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

/**
 * @file ERF_FireEngineRegistry.H
 * @brief The fire engines compiled into this executable
 *
 * ERF lists no engine itself. Each engine registers itself from its own .cpp
 * with
 *
 *     ERF_REGISTER_FIRE_ENGINE(MyEngine, "myengine")
 *
 * at namespace scope, and erf.fire.engine then selects it by that name. The
 * engine's directory (FIRE_HOME in GNU make, ERF_FIRE_ENGINE_DIRS in CMake)
 * adds that .cpp to the build, so adding an engine changes no file of ERF.
 *
 * The registration runs during static initialisation, so the .cpp that holds
 * it must be linked as an object file: inside a static archive (.a) that
 * nothing else references, the linker drops it and the engine is missing.
 * ERF's CMake build compiles engines into its OBJECT library, and GNU make
 * links every object, so both keep it; an engine that ships a prebuilt
 * library keeps its registering .cpp as source (LSFIRE's wrapper does).
 */

/// Creates an engine from the coupling inputs.
using FireEngineCreator = std::function<std::unique_ptr<FireEngine>(const FireCouplingParams&)>;

/**
 * @brief Add an engine to the registry; aborts on a name registered twice
 *
 * Called by ERF_REGISTER_FIRE_ENGINE; returns true so the macro can hold the
 * result in a static variable.
 */
bool register_fire_engine (const std::string& name, FireEngineCreator creator);

/// Names erf.fire.engine accepts, sorted; empty when no engine is compiled in.
std::vector<std::string> fire_engine_names ();

/// The engine erf.fire.engine names, or nullptr when no engine has that name.
std::unique_ptr<FireEngine> make_fire_engine (const std::string& name, const FireCouplingParams& params);

#define ERF_FIRE_ENGINE_CONCAT_(a, b) a##b
#define ERF_FIRE_ENGINE_CONCAT(a, b) ERF_FIRE_ENGINE_CONCAT_(a, b)

/**
 * @brief Register engine class Class under name: erf.fire.engine = name
 *        selects it. Class is constructed from the coupling inputs,
 *        Class(const FireCouplingParams&), or default-constructed when it has
 *        no such constructor.
 */
#define ERF_REGISTER_FIRE_ENGINE(Class, name)                                         \
    namespace {                                                                       \
    [[maybe_unused]] const bool ERF_FIRE_ENGINE_CONCAT(erf_fire_engine_registered_, __LINE__) = \
        register_fire_engine((name), [](const FireCouplingParams& p) -> std::unique_ptr<FireEngine> { \
            if constexpr (std::is_constructible_v<Class, const FireCouplingParams&>) {  \
                return std::make_unique<Class>(p);                                     \
            } else {                                                                  \
                amrex::ignore_unused(p);                                              \
                return std::make_unique<Class>();                                     \
            }                                                                         \
        });                                                                           \
    }

#endif
