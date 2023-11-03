/*------------------------------------------------------------
   Computes surface albedos
   For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
   Approximation for Solar Radiation in the NCAR Community Climate Model,
   Journal of Geophysical Research, Vol 97, D7, pp7603-7612).

   NOTE: this is the simplest formula for albedo calculation, we have to
         update to the latest implementation from microphysics output.
*/
#ifndef ERF_ALBEDO_H_
#define ERF_ALBEDO_H_
void set_albedo(const real1d& coszrs, real2d& albedo_dir, real2d& albedo_dif)
{
  // Albedos for land type I (Briegleb)
  auto nswbands = albedo_dir.extent(0);
  auto ncol     = albedo_dif.extent(1);
  yakl::c::parallel_for(yakl::c::Bounds<2>(nswbands,ncol), YAKL_LAMBDA (int ibnd, int icol) {  
    albedo_dir(ibnd, icol) = 1.4 * 0.24 / ( 1. + 0.8 * coszrs(icol));
    albedo_dif(ibnd, icol) = 1.2 * 0.24;
  });
}
#endif
