void more_shapes ()
{
    ParmParse pp("erf");
    int nvals_old = pp.countval("old_flag");       // ERF_DataStruct.H's no_substepping shape
    if (nvals_old > 0) {
        amrex::Error("The old_flag flag is deprecated -- set new_flag instead");
    }
    if (pp.contains("legacy")) {                    // AMReX_Amr.cpp's nosub shape: read, mapped, still accepted
        int legacy;
        pp.query("legacy", legacy);
        if (legacy > 0) {
            use_legacy = true;
        } else {
            amrex::Error("legacy <= 0 not allowed");
        }
    }
    int phv(0), chv(1);                            // AMReX_Amr.cpp's constructor-style declarators
    pp.query("chv", chv);
    std::vector<amrex::Real> box_lo(3), box_hi(3);  // ERF_WriteBndryPlanes.cpp's
    pp.queryarr("box_hi", box_hi);
}

void zlevels ()
{
    ParmParse pp("erf");
    int n_zlevels = pp.countval("z_levels");      // ERF_InitZLevels.cpp's shape: a size check, then the read
    if (n_zlevels > 0)
    {
        if (n_zlevels != 5) {
            amrex::Abort("You must specify a z_level for every value of k");
        }
        amrex::Vector<amrex::Real> zlev;
        pp.queryarr("z_levels", zlev);
    }
}
