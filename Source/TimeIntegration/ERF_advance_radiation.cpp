#include <ERF.H>

using namespace amrex;

#if defined(ERF_USE_RRTMGP)
void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const Real& dt_advance)
{
   bool do_sw_rad {true};
   bool do_lw_rad {true};
   bool do_aero_rad {true};
   bool do_snow_opt {true};
   bool is_cmip6_volcano {false};

    rad.initialize(cons,
                   qheating_rates[lev].get(),
                   lat_m[lev].get(),
                   lon_m[lev].get(),
                   qmoist[lev],
                   grids[lev],
                   Geom(lev),
                   dt_advance,
                   do_sw_rad,
                   do_lw_rad,
                   do_aero_rad,
                   do_snow_opt,
                   is_cmip6_volcano);
    rad.run();
    rad.on_complete();
}
#endif
