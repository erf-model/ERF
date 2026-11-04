void check_removed ()
{
    ParmParse pp("erf");
    if (pp.contains("old_knob")) {
        amrex::Abort("erf.old_knob has been removed; use erf.new_knob instead");
    }
    if (pp.contains("probe_only")) {
        amrex::Print() << "probe_only is set\n";
    }
}

void check_retired_more ()
{
    ParmParse pp("erf");
    bool use_terrain = false;
    if (pp.query("use_terrain", use_terrain) > 0) {
        amrex::Error("erf.use_terrain is deprecated; set erf.terrain_type instead");
    }
    if (pp.countval("no_substepping") > 0) {
        amrex::Error("no_substepping flag is deprecated");
    }
    int also_typed = 0;
    if (pp.contains("also_typed")) {
        pp.get("also_typed", also_typed);
    }
}
