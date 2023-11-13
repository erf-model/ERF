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
void set_albedo(const real1d& coszrs, real2d& albedo_dir, real2d& albedo_dif);
#endif
