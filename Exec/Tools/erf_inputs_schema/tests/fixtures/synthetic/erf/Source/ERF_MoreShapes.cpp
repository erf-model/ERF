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
