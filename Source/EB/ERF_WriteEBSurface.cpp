
#include <AMReX_WriteEBSurface.H>
#include <ERF.H>

void ERF::WriteMyEBSurface ()
{
  using namespace amrex;

  amrex::Print() << "Writing the geometry to a vtp file.\n" << std::endl;

  // Only write at the finest level!
  int lev = finest_level;

  BoxArray & ba            = grids[lev];
  DistributionMapping & dm = dmap[lev];

  const EBFArrayBoxFactory* ebfact = &EBFactory(lev);

  WriteEBSurface(ba,dm,Geom(lev),ebfact);
}
