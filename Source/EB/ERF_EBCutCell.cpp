#include <AMReX_Print.H>
#include <ERF_EBCutCell.H>

using namespace amrex;

#ifndef AMREX_USE_GPU
void
eb_cut_cell_::
debug ( int const a_face ) {

  amrex::Print() << "\n\nDEBUG THIS "
    << "--------------------------------------"
    << "\nisCovered?       " << isCovered()
    << "\nisRegular?       " << isRegular()
    << "\nisSingleValued?  " << isSingleValued()
    << "\n";

  if ( isCovered() || isRegular() ) { return; }

  amrex::RealVect v0(m_rbox.lo(0), m_rbox.lo(1), m_rbox.lo(2));
  amrex::RealVect v7(m_rbox.hi(0), m_rbox.hi(1), m_rbox.hi(2));

  amrex::Print() << "\n"
          << "lo:     " << v0 << '\n'
          << "hi:     " << v7 << '\n'
          << "p:      " << m_eb_point << '\n'
          << "n:      " << m_eb_normal << '\n'
          << "invert? " << m_invert << "\n\n";

  amrex::Print() << "Edge intersections:\n";
  calc_edge_intersections();
  amrex::Print() << '\n';

  if ( a_face == -1 || a_face == 1 ) { m_F1.report(1, v0); }
  if ( a_face == -1 || a_face == 2 ) { m_F2.report(2, v0); }
  if ( a_face == -1 || a_face == 3 ) { m_F3.report(3, v0); }
  if ( a_face == -1 || a_face == 4 ) { m_F4.report(4, v0); }
  if ( a_face == -1 || a_face == 5 ) { m_F5.report(5, v0); }
  if ( a_face == -1 || a_face == 6 ) { m_F6.report(6, v0); }
  if ( a_face == -1 || a_face == 7 ) { m_F7.report(7, v0); }
}
#endif
