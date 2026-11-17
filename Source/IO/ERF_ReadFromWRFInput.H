/**
 * \file ERF_ReadFromWRFInput.H
 */
#ifndef ERF_READ_FROM_WRFINPUT_H_
#define ERF_READ_FROM_WRFINPUT_H_

#include <AMReX_Box.H>
#include <AMReX_IntVect.H>

namespace erf_wrfinput {

/**
 * Shift horizontal WRF file coordinates into the ERF index space used by
 * init_from_wrfinput.  Keep k unchanged: WRF nests carry the parent's eta
 * column, while their horizontal origins may be cropped or offset.
 */
inline amrex::IntVect
file_to_erf_shift (const amrex::Box& erf_file_box,
                   const amrex::Box& file_subdomain)
{
    amrex::IntVect shift(0);
    shift[0] = erf_file_box.smallEnd(0) - file_subdomain.smallEnd(0);
    shift[1] = erf_file_box.smallEnd(1) - file_subdomain.smallEnd(1);
    return shift;
}

/** Map one file-coordinate box into the ERF index space used by init_from_wrfinput. */
inline amrex::Box
file_box_to_erf (amrex::Box file_box, const amrex::IntVect& shift)
{
    file_box.shift(shift);
    return file_box;
}

} // namespace erf_wrfinput

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

#endif
