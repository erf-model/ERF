#ifndef ERF_STRUCTURE_MODEL_REGISTRY_H
#define ERF_STRUCTURE_MODEL_REGISTRY_H

#include <ERF_StructureModel.H>

#include <functional>
#include <memory>
#include <string>
#include <vector>

/**
 * @file ERF_StructureModelRegistry.H
 * @brief The structure models compiled into this executable
 *
 * ERF lists no model itself. Each model registers itself from its own .cpp with
 *
 *     ERF_REGISTER_STRUCTURE_MODEL(MyModel, "mymodel")
 *
 * at namespace scope, and erf.structures.model then selects it by that name. The
 * model's directory (STRUCTURE_HOME in GNU make, ERF_STRUCTURE_MODEL_DIRS in
 * CMake) adds that .cpp to the build, so adding a model changes no file of ERF.
 *
 * The registration runs during static initialisation, so the .cpp that holds it
 * must be linked as an object file (ERF's CMake build compiles model directories
 * into its OBJECT library, and GNU make links every object). This is the scheme of
 * the fire engines (ERF_FireEngineRegistry.H).
 */

/// Creates a model; the model reads its own inputs in its constructor.
using StructureModelCreator = std::function<std::unique_ptr<StructureModel>(int max_level)>;

/**
 * @brief Add a model to the registry; aborts on an empty name or a name registered twice
 *
 * Called by ERF_REGISTER_STRUCTURE_MODEL; returns true so the macro can hold the
 * result in a static variable.
 */
bool register_structure_model (const std::string& name, StructureModelCreator creator);

/// Names erf.structures.model accepts, sorted; empty when no model is compiled in.
std::vector<std::string> structure_model_names ();

/// The model registered as name, or nullptr when no model has that name.
std::unique_ptr<StructureModel> make_structure_model (const std::string& name, int max_level);

/**
 * @brief The model erf.structures.model selects, or nullptr when the key is not set
 *
 * Aborts when the key names no registered model, listing the names this
 * executable accepts.
 */
std::unique_ptr<StructureModel> make_structure_model_from_inputs (int max_level);

/**
 * @brief Stop the run when the model's plot_var_names() cannot all be plotted
 *
 * Aborts, naming the model and the name, on an empty name, a name the model lists
 * twice, or a name in taken (the names ERF and the other models already provide).
 * ERF's plotfile assembly would otherwise drop the name from the list and still
 * fill it, and fail with a less helpful message.
 */
void check_structure_plot_names (const StructureModel& model, const std::vector<std::string>& taken);

#define ERF_STRUCTURE_MODEL_CONCAT_(a, b) a##b
#define ERF_STRUCTURE_MODEL_CONCAT(a, b) ERF_STRUCTURE_MODEL_CONCAT_(a, b)

/**
 * @brief Register model class Class under name: erf.structures.model = name
 *        selects it. Class is constructed as Class(int max_level).
 */
#define ERF_REGISTER_STRUCTURE_MODEL(Class, name)                                             \
    namespace {                                                                               \
    [[maybe_unused]] const bool ERF_STRUCTURE_MODEL_CONCAT(erf_structure_model_registered_, __LINE__) = \
        register_structure_model((name), [](int max_level) -> std::unique_ptr<StructureModel> { \
            return std::make_unique<Class>(max_level);                                        \
        });                                                                                   \
    }

#endif
