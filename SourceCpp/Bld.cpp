#include <AMReX_LevelBld.H>

#include "ERF.H"

class ERFBld : public amrex::LevelBld
{
  virtual void variableSetUp();
  virtual void variableCleanUp();
  virtual amrex::AmrLevel* operator()();
  virtual amrex::AmrLevel* operator()(
    amrex::Amr& papa,
    int lev,
    const amrex::Geometry& level_geom,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    amrex::Real time);
};

ERFBld ERF_bld;

amrex::LevelBld*
getLevelBld()
{
  return &ERF_bld;
}

void
ERFBld::variableSetUp()
{
  ERF::variableSetUp();
}

void
ERFBld::variableCleanUp()
{
  ERF::variableCleanUp();
}

amrex::AmrLevel*
ERFBld::operator()()
{
  return new ERF;
}

amrex::AmrLevel*
ERFBld::operator()(
  amrex::Amr& papa,
  int lev,
  const amrex::Geometry& level_geom,
  const amrex::BoxArray& ba,
  const amrex::DistributionMapping& dm,
  amrex::Real time)
{
  return new ERF(papa, lev, level_geom, ba, dm, time);
}
