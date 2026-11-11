void read_kinds ()
{
    ParmParse pp("xlo");
    {
        std::vector<amrex::Real> v;
        pp.queryarr("velocity", v, 0, 3);          // ERF_InitBCs.cpp's block-local list
    }
    amrex::Real z, u, v, w, th;                    // ... and its scalar of the same name
    pp.query("speed", v);
    {
        int k = 1;
        pp.query("k_int", k);
    }
    {
        amrex::Real k = 0.5;
        pp.query("k_real", k);                     // the nearest declaration before the read decides
    }
}
