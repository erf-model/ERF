#ifndef ERF_DUST_MODEL_REGISTRY_H
#define ERF_DUST_MODEL_REGISTRY_H

#include <ERF_DustModel.H>

#include <functional>
#include <memory>
#include <string>
#include <vector>

/**
 * @file ERF_DustModelRegistry.H
 * @brief The dust models compiled into this executable
 *
 * ERF lists no model itself. Each model registers itself from its own .cpp with
 *
 *     ERF_REGISTER_DUST_MODEL(MyModel, "mymodel")
 *
 * at namespace scope, and erf.dust.model then selects it by that name. The
 * model's directory (DUST_HOME in GNU make, ERF_DUST_MODEL_DIRS in CMake) adds
 * that .cpp to the build, so adding a model changes no file of ERF.
 *
 * The registration runs during static initialisation, so the .cpp that holds it
 * must be linked as an object file (ERF's CMake build compiles model directories
 * into its OBJECT library, and GNU make links every object). This is the scheme of
 * the fire engines (ERF_FireEngineRegistry.H) and the structure models.
 */

/// Creates a model; the model reads its inputs in initialize().
using DustModelCreator = std::function<std::unique_ptr<DustModel>()>;

/**
 * @brief Add a model to the registry; aborts on an empty name or a name registered twice
 *
 * Called by ERF_REGISTER_DUST_MODEL; returns true so the macro can hold the
 * result in a static variable.
 */
bool register_dust_model (const std::string& name, DustModelCreator creator);

/// Names erf.dust.model accepts, sorted; empty when no model is compiled in.
std::vector<std::string> dust_model_names ();

/// The model registered as name, or nullptr when no model has that name.
std::unique_ptr<DustModel> make_dust_model (const std::string& name);

/**
 * @brief The model the inputs select, or nullptr when erf.dust.enable is not true
 *
 * erf.dust.model names the model; with one model compiled in it may be left
 * out. Aborts when dust is enabled and no model, or no model of that name, is
 * compiled in, listing the names this executable accepts, and when the model
 * fills more dust scalars than the build has (NSCALARS_DUST).
 */
std::unique_ptr<DustModel> make_dust_model_from_inputs ();

#define ERF_DUST_MODEL_CONCAT_(a, b) a##b
#define ERF_DUST_MODEL_CONCAT(a, b) ERF_DUST_MODEL_CONCAT_(a, b)

/**
 * @brief Register model class Class under name: erf.dust.model = name
 *        selects it. Class is default-constructed.
 */
#define ERF_REGISTER_DUST_MODEL(Class, name)                                                  \
    namespace {                                                                               \
    [[maybe_unused]] const bool ERF_DUST_MODEL_CONCAT(erf_dust_model_registered_, __LINE__) = \
        register_dust_model((name), []() -> std::unique_ptr<DustModel> {                      \
            return std::make_unique<Class>();                                                 \
        });                                                                                   \
    }

#endif
