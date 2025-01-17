#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB_utils.H>

#include <ERF_EB.H>

using namespace amrex;

/************************************************************************
 *                                                                      *
 * Function to create a simple box geometry:                            *
 *        -> box.{Lo,Hi} vector storing box lo/hi                       *
 *        -> box.offset  vector storing box offset                      *
 *                                                                      *
 * NOTE: walls are placed _outside_ domain for periodic directions.     *
 *                                                                      *
 ***********************************************************************/
void
eb_::
make_box ( Geometry const& a_geom)
{
  ParmParse pp_box("eb.box");

  bool inside = true;

  Vector<Real> boxLo(AMREX_SPACEDIM), boxHi(AMREX_SPACEDIM);

  Real const* dx = a_geom.CellSize();

  if ( !almostEqual(dx[0],dx[1]) || !almostEqual(dx[1],dx[2])) {
    amrex::Error(" EB Error: Mesh spacing must be uniform!.\n");
  }

  Real offset = 0.01*dx[0];

  for (int i = 0; i < AMREX_SPACEDIM; i++) {
      boxLo[i] = a_geom.ProbLo(i);
      boxHi[i] = a_geom.ProbHi(i);
  }

  pp_box.queryarr("lo", boxLo,  0, AMREX_SPACEDIM);
  pp_box.queryarr("hi", boxHi,  0, AMREX_SPACEDIM);

  pp_box.query("internal_flow", inside);
  pp_box.query("offset", offset);

  Real xlo = boxLo[0] + offset;
  Real xhi = boxHi[0] - offset;

  // This ensures that the walls won't even touch the ghost cells. By
  // putting them one domain width away
  if (a_geom.isPeriodic(0))
  {
      xlo = 2.0*a_geom.ProbLo(0) - a_geom.ProbHi(0);
      xhi = 2.0*a_geom.ProbHi(0) - a_geom.ProbLo(0);
  }


  Real ylo = boxLo[1] + offset;
  Real yhi = boxHi[1] - offset;

  // This ensures that the walls won't even touch the ghost cells. By
  // putting them one domain width away
  if (a_geom.isPeriodic(1))
  {
      ylo = 2.0*a_geom.ProbLo(1) - a_geom.ProbHi(1);
      yhi = 2.0*a_geom.ProbHi(1) - a_geom.ProbLo(1);
  }

  Real zlo = boxLo[2] + offset;
  Real zhi = boxHi[2] - offset;

  // This ensures that the walls won't even touch the ghost cells. By
  // putting them one domain width away
  if (a_geom.isPeriodic(2))
  {
      zlo = 2.0*a_geom.ProbLo(2) - a_geom.ProbHi(2);
      zhi = 2.0*a_geom.ProbHi(2) - a_geom.ProbLo(2);
  }

  Array<Real,3>  point_lox{ xlo, 0.0, 0.0};
  Array<Real,3> normal_lox{-1.0, 0.0, 0.0};
  Array<Real,3>  point_hix{ xhi, 0.0, 0.0};
  Array<Real,3> normal_hix{ 1.0, 0.0, 0.0};

  Array<Real,3>  point_loy{0.0, ylo, 0.0};
  Array<Real,3> normal_loy{0.0,-1.0, 0.0};
  Array<Real,3>  point_hiy{0.0, yhi, 0.0};
  Array<Real,3> normal_hiy{0.0, 1.0, 0.0};

  Array<Real,3>  point_loz{0.0, 0.0, zlo};
  Array<Real,3> normal_loz{0.0, 0.0,-1.0};
  Array<Real,3>  point_hiz{0.0, 0.0, zhi};
  Array<Real,3> normal_hiz{0.0, 0.0, 1.0};

  EB2::PlaneIF plane_lox(point_lox,normal_lox);
  EB2::PlaneIF plane_hix(point_hix,normal_hix);

  EB2::PlaneIF plane_loy(point_loy,normal_loy);
  EB2::PlaneIF plane_hiy(point_hiy,normal_hiy);

  EB2::PlaneIF plane_loz(point_loz,normal_loz);
  EB2::PlaneIF plane_hiz(point_hiz,normal_hiz);

  auto box = EB2::makeUnion(plane_lox, plane_hix, plane_loy,
                            plane_hiy, plane_loz, plane_hiz );


  if( inside ) {

    auto gshop = EB2::makeShop(box);
    build_level(a_geom, gshop);

  } else {

    auto xob = EB2::makeComplement(box);
    auto gshop = EB2::makeShop(xob);

    build_level(a_geom, gshop);
  }

}

