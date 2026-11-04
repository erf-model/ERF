#include "ERF_Synthetic.H"

void read_synthetic (SyntheticParams& sp)
{
    ParmParse pp("erf");
    pp.query("cfl", sp.cfl);                       // a real with a default from the header
    pp.get("nsteps", sp.nsteps);                   // required: a get
    MoistureType moisture_type = MoistureType::None;
    pp.query_enum_case_insensitive("moisture_model", moisture_type);
#ifdef ERF_USE_WIDGET
    bool widget_on = false;
    pp.query("widget_on", widget_on);              // read only when ERF_USE_WIDGET is defined
#endif
    amrex::Vector<amrex::Real> Cs;
    query_one_or_per_level(pp, "Cs", Cs, 0, 1);    // one value or one per level

    ParmParse pa("amr");
    int max_level = 0;
    pa.query("max_level", max_level);              // also read by AMReX: a shared key
    pp.query("eye_latitude", sp.eye_latitude);
    pp.query("eye_longitude", sp.eye_longitude);
    amrex::Vector<MoistureType> species;
    pp.queryarr("moisture_species", species);     // a list of enum values
}
