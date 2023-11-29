
#include <rrtmgp_const.h>
#include "Albedo.H"

void set_albedo(const real1d& coszrs, real2d& albedo_dir, real2d& albedo_dif)
{
  // Albedos for land type I (Briegleb)
  auto nswbands = albedo_dir.extent(0);
  auto ncol     = albedo_dif.extent(1);
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  parallel_for(SimpleBounds<2>(nswbands,ncol), YAKL_LAMBDA (int ibnd, int icol) {
    albedo_dir(ibnd, icol) = 1.4 * 0.24 / ( 1. + 0.8 * coszrs(icol));
    albedo_dif(ibnd, icol) = 1.2 * 0.24;
  });
}

