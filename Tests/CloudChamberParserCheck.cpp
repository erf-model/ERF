#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include "../Source/Prob/ERF_CloudChamber.H"

int main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    MoistureType model = MoistureType::None;
    amrex::ParmParse pp_erf("erf");
    pp_erf.query_enum_case_insensitive("moisture_model", model);
    const amrex::Real problo[AMREX_SPACEDIM] = {0.0, 0.0, 0.0};
    const amrex::Real probhi[AMREX_SPACEDIM] = {2.0, 2.0, 1.0};
    const auto config = erf_cloud_chamber::parse_config(problo, probhi, model);
    amrex::ignore_unused(config);
    amrex::Finalize();
    return 0;
}
