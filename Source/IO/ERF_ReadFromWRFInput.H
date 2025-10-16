#ifdef ERF_USE_NETCDF
void
read_from_wrfinput (int lev,
                    const amrex::Box& subdomain,
                    const std::string& fname,
                    amrex::FArrayBox& NC_fab,
                    const std::string& NC_name,
                    amrex::Geometry& geom,
                    int& use_theta_m,
                    int& success);
#endif
