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
