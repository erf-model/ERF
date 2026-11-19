#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include "../Source/Prob/ERF_CloudChamber.H"

int main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    MoistureType model = MoistureType::None;
    amrex::ParmParse pp_erf("erf");
    pp_erf.query_enum_case_insensitive("moisture_model", model);
    DiffChoice dc;
    dc.init_params("erf");
    amrex::IntVect ncells;
    amrex::ParmParse("amr").getarr("n_cell", ncells);
    const amrex::Geometry geom(amrex::Box(amrex::IntVect(0), ncells - amrex::IntVect(1)));
    const auto config = erf_cloud_chamber::parse_config(geom, model, dc);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        config.has_moisture == (model != MoistureType::None),
        "Cloud Chamber: parsed moisture flag disagrees with the selected model");
    amrex::Finalize();
    return 0;
}
