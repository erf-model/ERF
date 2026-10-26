/**
 * \file ERF_ReadFromWRFInput.H
 */
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

void
check_wrfinput_geometry (int lev,
                         const std::string& fname,
                         const amrex::Geometry& geom,
                         int& use_theta_m);

int
wrfinput_var_box (const amrex::Box& subdomain_to_read,
                  const std::string& fname,
                  const std::string& NC_name,
                  amrex::Box& var_box);

void
read_region_from_wrfinput (const std::string& fname,
                           const std::string& NC_name,
                           const amrex::Box& var_box,
                           const amrex::Box& region,
                           amrex::FArrayBox& fab,
                           int max_readers);
#endif
