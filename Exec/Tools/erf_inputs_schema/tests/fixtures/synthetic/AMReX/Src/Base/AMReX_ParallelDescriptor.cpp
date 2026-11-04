// A tiny AMReX-like source for the schema generator's tests.
void init_synthetic ()
{
    ParmParse pp("amrex");
    bool use_gpu_aware_mpi = false;
    pp.query("use_gpu_aware_mpi", use_gpu_aware_mpi);
    ParmParse pa("amr");
    int max_level = 0;
    pa.query("max_level", max_level);
}
