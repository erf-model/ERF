#include <ERF_DustModelRegistry.H>
#include <ERF_IndexDefines.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <map>

namespace {
// A function-local static, so it exists before the first registration whatever
// the order in which the models' translation units are initialised
std::map<std::string, DustModelCreator>&
dust_model_registry ()
{
    static std::map<std::string, DustModelCreator> registry;
    return registry;
}

std::string
known_dust_models ()
{
    std::string known;
    for (const auto& n : dust_model_names()) { known += (known.empty() ? "" : ", ") + ("\"" + n + "\""); }
    return known.empty() ? std::string("none; add a model directory with ERF_DUST_MODEL_DIRS"
                                       " in CMake or DUST_HOME in GNU make")
                         : known;
}
}

bool
register_dust_model (const std::string& name, DustModelCreator creator)
{
    auto& registry = dust_model_registry();
    if (name.empty()) {
        amrex::Abort("ERF_REGISTER_DUST_MODEL: a dust model registered an empty name");
    }
    if (registry.count(name) > 0) {
        amrex::Abort("ERF_REGISTER_DUST_MODEL: two dust models registered the name \"" + name
                     + "\"; each model directory in the build must use its own name");
    }
    registry.emplace(name, std::move(creator));
    return true;
}

std::vector<std::string>
dust_model_names ()
{
    std::vector<std::string> names;
    for (const auto& kv : dust_model_registry()) { names.push_back(kv.first); }
    return names;   // std::map keeps them sorted
}

std::unique_ptr<DustModel>
make_dust_model (const std::string& name)
{
    const auto& registry = dust_model_registry();
    const auto it = registry.find(name);
    if (it == registry.end()) { return nullptr; }
    return it->second();
}

std::unique_ptr<DustModel>
make_dust_model_from_inputs ()
{
    amrex::ParmParse pp("erf.dust");
    bool enable = false;
    pp.query("enable", enable);
    if (!enable) { return nullptr; }

    std::string name;
    if (!pp.query("model", name)) {
        const std::vector<std::string> names = dust_model_names();
        if (names.size() != 1) {
            amrex::Abort("erf.dust.enable = true needs erf.dust.model to name one of the dust models compiled in ("
                         + known_dust_models() + ")");
        }
        name = names.front();
    }

    std::unique_ptr<DustModel> model = make_dust_model(name);
    if (!model) {
        amrex::Abort("erf.dust.model = \"" + name + "\" is not one of the dust models compiled in ("
                     + known_dust_models() + ")");
    }
    if (model->num_scalars() < 1 || model->num_scalars() > NSCALARS_DUST) {
        amrex::Abort("dust model \"" + name + "\" fills " + std::to_string(model->num_scalars())
                     + " dust scalar(s), but this build has " + std::to_string(NSCALARS_DUST)
                     + "; rebuild with ERF_NUM_DUST_SCALARS (CMake) or NUM_DUST_SCALARS (GNU make) = "
                     + std::to_string(model->num_scalars()));
    }
    return model;
}
