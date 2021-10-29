#include <AMReX_ParmParse.H>

#include "ERF.H"
#include "Tagging.H"

TaggingParm tparm;

void
ERF::read_tagging_params()
{
  amrex::ParmParse pp("tagging");

  pp.query("denerr", tparm.denerr);
  pp.query("max_denerr_lev", tparm.max_denerr_lev);
  pp.query("dengrad", tparm.dengrad);
  pp.query("max_dengrad_lev", tparm.max_dengrad_lev);

  pp.query("presserr", tparm.presserr);
  pp.query("max_presserr_lev", tparm.max_presserr_lev);
  pp.query("pressgrad", tparm.pressgrad);
  pp.query("max_pressgrad_lev", tparm.max_pressgrad_lev);

  pp.query("velerr", tparm.velerr);
  pp.query("max_velerr_lev", tparm.max_velerr_lev);
  pp.query("velgrad", tparm.velgrad);
  pp.query("max_velgrad_lev", tparm.max_velgrad_lev);

  pp.query("vorterr", tparm.vorterr);
  pp.query("max_vorterr_lev", tparm.max_vorterr_lev);

  pp.query("temperr", tparm.temperr);
  pp.query("max_temperr_lev", tparm.max_temperr_lev);
  pp.query("tempgrad", tparm.tempgrad);
  pp.query("max_tempgrad_lev", tparm.max_tempgrad_lev);

  pp.query("ftracerr", tparm.ftracerr);
  pp.query("max_ftracerr_lev", tparm.max_ftracerr_lev);
  pp.query("ftracgrad", tparm.ftracgrad);
  pp.query("max_ftracgrad_lev", tparm.max_ftracgrad_lev);

  pp.query("vfracerr", tparm.vfracerr);
  pp.query("max_vfracerr_lev", tparm.max_vfracerr_lev);

  pp.query("tag_region", tparm.tag_region);
  if (tparm.tag_region)
  {
      tparm.region_lo.resize(3);
      tparm.region_hi.resize(3);
      pp.getarr("region_lo", tparm.region_lo);
      pp.getarr("region_hi", tparm.region_hi);
  }
}
