#include <ERF_StructureModelRegistry.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <map>

namespace {
// A function-local static, so it exists before the first registration whatever
// the order in which the models' translation units are initialised
std::map<std::string, StructureModelCreator>&
structure_model_registry ()
{
    static std::map<std::string, StructureModelCreator> registry;
    return registry;
}
}

bool
register_structure_model (const std::string& name, StructureModelCreator creator)
{
    auto& registry = structure_model_registry();
    if (name.empty()) {
        amrex::Abort("ERF_REGISTER_STRUCTURE_MODEL: a structure model registered an empty name");
    }
    if (registry.count(name) > 0) {
        amrex::Abort("ERF_REGISTER_STRUCTURE_MODEL: two structure models registered the name \"" + name
                     + "\"; each model directory in the build must use its own name");
    }
    registry.emplace(name, std::move(creator));
    return true;
}

std::vector<std::string>
structure_model_names ()
{
    std::vector<std::string> names;
    for (const auto& kv : structure_model_registry()) { names.push_back(kv.first); }
    return names;   // std::map keeps them sorted
}

std::unique_ptr<StructureModel>
make_structure_model (const std::string& name, int max_level)
{
    const auto& registry = structure_model_registry();
    const auto it = registry.find(name);
    if (it == registry.end()) { return nullptr; }
    return it->second(max_level);
}

std::unique_ptr<StructureModel>
make_structure_model_from_inputs (int max_level)
{
    amrex::ParmParse pp("erf.structures");
    std::string name;
    if (!pp.query("model", name)) { return nullptr; }

    std::unique_ptr<StructureModel> model = make_structure_model(name, max_level);
    if (!model) {
        std::string known;
        for (const auto& n : structure_model_names()) { known += (known.empty() ? "" : ", ") + ("\"" + n + "\""); }
        amrex::Abort("erf.structures.model = \"" + name + "\" is not one of the structure models compiled in ("
                     + (known.empty() ? std::string("none; add a model directory with ERF_STRUCTURE_MODEL_DIRS"
                                                    " in CMake or STRUCTURE_HOME in GNU make")
                                      : known) + ")");
    }
    return model;
}
