#ifndef ERF_EB_CUT_CELL_H_
#define ERF_EB_CUT_CELL_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_RealBox.H>
#include <AMReX_RealVect.H>
#include <AMReX_Algorithm.H>
#include <AMReX_EBCellFlag.H>

#include "ERF_EBPolygon.H"
#include "ERF_EBUtils.H"

class eb_cut_cell_ {

  public:

    AMREX_GPU_HOST_DEVICE
    eb_cut_cell_ ( amrex::EBCellFlag const& a_flag,
                   amrex::RealBox    const& a_rbox,
                   amrex::RealVect   const& a_point,
                   amrex::RealVect   const& a_normal );

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isCovered () const noexcept { return m_flag.isCovered(); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isRegular () const noexcept { return m_flag.isRegular(); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isSingleValued () const noexcept { return m_flag.isSingleValued(); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_covered () {
      m_flag.setCovered();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_regular () {

      m_flag.setRegular();

      m_F1.set_area( m_rbox.length(0)*m_rbox.length(1) );
      m_F3.set_area( m_rbox.length(0)*m_rbox.length(1) );

      m_F2.set_area( m_rbox.length(1)*m_rbox.length(2) );
      m_F4.set_area( m_rbox.length(1)*m_rbox.length(2) );

      m_F5.set_area( m_rbox.length(0)*m_rbox.length(2) );
      m_F6.set_area( m_rbox.length(0)*m_rbox.length(2) );

      m_F7.set_area(0.);

    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real volume () const {

      if (m_flag.isCovered() ) { return 0.; }
      if (m_flag.isRegular() ) { return m_rbox.volume(); }

      amrex::Real volume(0.);

      amrex::Real const* lo = m_rbox.lo();
      amrex::RealVect v0(lo[0], lo[1], lo[2]);

      if (m_F2.ok() ) { volume += m_F2.area() * m_F2.distance(v0); }
      if (m_F3.ok() ) { volume += m_F3.area() * m_F3.distance(v0); }
      if (m_F6.ok() ) { volume += m_F6.area() * m_F6.distance(v0); }
      if (m_F7.ok() ) { volume += m_F7.area() * m_F7.distance(v0); }

      volume /= 3.;

      return m_invert*volume + (1.-m_invert)*(m_rbox.volume()-volume);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real areaLo ( int const idir ) const noexcept {
      AMREX_ASSERT( idir >=0 && idir < AMREX_SPACEDIM );
      if (m_flag.isCovered() ) { return 0.; }
      if (m_flag.isRegular() ) { return m_lo_faces[idir]->area(); }
      amrex::Real const area(m_lo_faces[idir]->area());
      return m_invert*area + (1.-m_invert)*(m_rbox_area[idir] - area);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real areaHi ( int const idir ) const noexcept {
      AMREX_ASSERT( idir >= 0 && idir < AMREX_SPACEDIM );
      if (m_flag.isCovered() ) { return 0.; }
      if (m_flag.isRegular() ) { return m_hi_faces[idir]->area(); }
      amrex::Real const area(m_hi_faces[idir]->area());
      return m_invert*area + (1.-m_invert)*(m_rbox_area[idir] - area);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real areaBoun () const noexcept {
      return m_F7.area();
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect centLo ( int const idir ) const noexcept {
      AMREX_ASSERT( idir >=0 && idir < AMREX_SPACEDIM );
      amrex::RealVect cent = m_cellface_cent[ m_lo_faces_id[idir] ];
      if (m_flag.isCovered() || m_flag.isRegular() ) {
        // Default cent
      } else {
        amrex::Real const area_R = m_lo_faces[idir]->area();
        if (amrex::almostEqual(area_R,0.0) || amrex::almostEqual(area_R, m_rbox_area[idir]) ){
          // Default cent
        } else {
          amrex::RealVect   cent_O = cent;
          amrex::RealVect   cent_R = m_lo_faces[idir]->get_centroid();
          amrex::Real const area_C = m_rbox_area[idir] - area_R;
          cent = m_invert * cent_R + (1.-m_invert) * ((1.+area_R/area_C)*cent_O - area_R/area_C*cent_R);
        }
      }
      return cent;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect centHi ( int const idir ) const noexcept {
      AMREX_ASSERT( idir >=0 && idir < AMREX_SPACEDIM );
      amrex::RealVect cent = m_cellface_cent[ m_hi_faces_id[idir] ];
      if (m_flag.isCovered() || m_flag.isRegular() ) {
        // Default cent
      } else {
        amrex::Real const area_R = m_hi_faces[idir]->area();
        if (amrex::almostEqual(area_R,0.0) || amrex::almostEqual(area_R, m_rbox_area[idir]) ){
          // Default cent
        } else {
          amrex::RealVect   cent_O = cent;
          amrex::RealVect   cent_R = m_hi_faces[idir]->get_centroid();
          amrex::Real const area_C = m_rbox_area[idir] - area_R;
          cent = m_invert * cent_R + (1.-m_invert) * ((1.+area_R/area_C)*cent_O - area_R/area_C*cent_R);
        }
      }
      return cent;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect centBoun () const noexcept {
      amrex::RealVect cent{0.,0.,0.};
      if (m_flag.isSingleValued()) {
        cent = m_F7.get_centroid();
      }
      return cent;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect normBoun () const noexcept {
      amrex::RealVect normal{0.,0.,0.};
      if (m_flag.isSingleValued()) {
        normal = m_eb_normal;
      }
      return normal;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect centVol () const noexcept {
      amrex::RealVect vcent = m_cell_cent;
      if (m_flag.isSingleValued()) {
        amrex::Real xm = m_rbox.lo(0);
        amrex::Real xp = m_rbox.hi(0);
        amrex::Real ym = m_rbox.lo(1);
        amrex::Real yp = m_rbox.hi(1);
        amrex::Real zm = m_rbox.lo(2);
        amrex::Real zp = m_rbox.hi(2);

        amrex::Real axm = areaLo(0);
        amrex::Real axp = areaHi(0);
        amrex::Real aym = areaLo(1);
        amrex::Real ayp = areaHi(1);
        amrex::Real azm = areaLo(2);
        amrex::Real azp = areaHi(2);

        amrex::Real barea = m_F7.area();
        amrex::RealVect bcent = m_F7.get_centroid();

        amrex::Real vol = volume();

        vcent[0] = 0.5 * ( - xm * xm * axm + xp * xp * axp + bcent[0] * bcent[0] * m_eb_normal[0] * barea ) / vol;
        vcent[1] = 0.5 * ( - ym * ym * aym + yp * yp * ayp + bcent[1] * bcent[1] * m_eb_normal[1] * barea ) / vol;
        vcent[2] = 0.5 * ( - zm * zm * azm + zp * zp * azp + bcent[2] * bcent[2] * m_eb_normal[2] * barea ) / vol;
      }
      return vcent;
    }

    void debug ( int const a_face = -1 );

  private:

    amrex::RealBox  const m_rbox;
    amrex::RealVect const m_eb_point;
    amrex::RealVect const m_eb_normal;

    amrex::Real m_invert;

    amrex::RealVect m_rbox_area;

    amrex::EBCellFlag m_flag;

    // Cell faces

    polygon_ m_F1;
    polygon_ m_F2;
    polygon_ m_F3;
    polygon_ m_F4;
    polygon_ m_F5;
    polygon_ m_F6;

    static int constexpr m_max_faces = 6;
    amrex::Array<amrex::RealVect,m_max_faces> m_cellface_cent; // Centroids of cell faces
    amrex::RealVect m_cell_cent; // Centroid of cell
    amrex::RealVect m_dx;

    // cut face
    polygon_ m_F7;

    amrex::Array<polygon_ const* const,3> m_lo_faces {&m_F4, &m_F5, &m_F1};
    amrex::Array<polygon_ const* const,3> m_hi_faces {&m_F2, &m_F6, &m_F3};

    amrex::Array<int const,3> m_lo_faces_id {3,4,0};
    amrex::Array<int const,3> m_hi_faces_id {1,5,2};

    AMREX_GPU_HOST_DEVICE
    void calc_edge_intersections ();

    AMREX_GPU_HOST_DEVICE
    void set_covered_regular_cell_vertices ();
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
eb_cut_cell_::
eb_cut_cell_ ( amrex::EBCellFlag const& a_flag,
               amrex::RealBox    const& a_rbox,
               amrex::RealVect   const& a_point,
               amrex::RealVect   const& a_normal )
  : m_rbox(a_rbox)
  , m_eb_point(a_point)
  , m_eb_normal(a_normal)
  , m_invert(0.0)
  , m_F1(a_point, a_normal)
  , m_F2(a_point, a_normal)
  , m_F3(a_point, a_normal)
  , m_F4(a_point, a_normal)
  , m_F5(a_point, a_normal)
  , m_F6(a_point, a_normal)
  , m_cellface_cent({amrex::RealVect(0.), amrex::RealVect(0.), amrex::RealVect(0.),
                 amrex::RealVect(0.), amrex::RealVect(0.), amrex::RealVect(0.)})
{
  using namespace amrex;

  m_rbox_area[0] = m_rbox.length(1)*m_rbox.length(2);
  m_rbox_area[1] = m_rbox.length(0)*m_rbox.length(2);
  m_rbox_area[2] = m_rbox.length(0)*m_rbox.length(1);

  RealVect v0(m_rbox.lo(0), m_rbox.lo(1), m_rbox.lo(2));
  RealVect v1(m_rbox.hi(0), m_rbox.lo(1), m_rbox.lo(2));
  RealVect v2(m_rbox.lo(0), m_rbox.hi(1), m_rbox.lo(2));
  RealVect v3(m_rbox.lo(0), m_rbox.lo(1), m_rbox.hi(2));
  RealVect v4(m_rbox.hi(0), m_rbox.lo(1), m_rbox.hi(2));
  RealVect v5(m_rbox.hi(0), m_rbox.hi(1), m_rbox.lo(2));
  RealVect v6(m_rbox.lo(0), m_rbox.hi(1), m_rbox.hi(2));
  RealVect v7(m_rbox.hi(0), m_rbox.hi(1), m_rbox.hi(2));

  // Centroids of cell faces

  m_cellface_cent[0] = 0.25 * ( v0 + v2 + v5 + v1 ); // F1
  m_cellface_cent[1] = 0.25 * ( v1 + v5 + v7 + v4 ); // F2
  m_cellface_cent[2] = 0.25 * ( v3 + v4 + v7 + v6 ); // F3
  m_cellface_cent[3] = 0.25 * ( v0 + v3 + v6 + v2 ); // F4
  m_cellface_cent[4] = 0.25 * ( v0 + v1 + v4 + v3 ); // F5
  m_cellface_cent[5] = 0.25 * ( v2 + v5 + v7 + v6 ); // F6

  // Cell centroid

  m_cell_cent[0] = 0.5 * ( m_rbox.lo(0) + m_rbox.hi(0) );
  m_cell_cent[1] = 0.5 * ( m_rbox.lo(1) + m_rbox.hi(1) );
  m_cell_cent[2] = 0.5 * ( m_rbox.lo(2) + m_rbox.hi(2) );

  // Cell size

  m_dx[0] = m_rbox.hi(0) - m_rbox.lo(0);
  m_dx[1] = m_rbox.hi(1) - m_rbox.lo(1);
  m_dx[2] = m_rbox.hi(2) - m_rbox.lo(2);

  if (a_flag.isCovered() ) {

    set_covered();

  } else if (a_flag.isRegular() ) {

    set_regular();

  } else { // Check that the box and plane intersect.

    RealVect c = 0.5*(v0 + v7);
    RealVect e = v7 - c;

    Real r = e[0]*amrex::Math::abs(a_normal[0]) +
             e[1]*amrex::Math::abs(a_normal[1]) +
             e[2]*amrex::Math::abs(a_normal[2]);

    Real s = amrex::Math::abs(c.dotProduct(a_normal)
               - a_point.dotProduct(a_normal));

    if (s > r) {
      if ((a_normal.dotProduct(v0) - a_normal.dotProduct(a_point)) > 0.)
      { set_covered(); } else { set_regular(); }
    } else { m_flag.setSingleValued(); }
  }

  if ( m_flag.isSingleValued() ) {

    m_invert = ((m_eb_normal.dotProduct(v0) - m_eb_normal.dotProduct(m_eb_point)) > 0.) ? 0.0 : 1.0;

    calc_edge_intersections();

  } // end singleValued

  m_F1.define();
  m_F2.define();
  m_F3.define();
  m_F4.define();
  m_F5.define();
  m_F6.define();
  m_F7.define();

  // For covered and regular cut cells, add vertices to utilize e.g., get_centroid.
  if ( !m_flag.isSingleValued() ) {
    set_covered_regular_cell_vertices();
  }

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
eb_cut_cell_::
calc_edge_intersections ()
{
  using namespace amrex;

  RealVect v0(m_rbox.lo(0), m_rbox.lo(1), m_rbox.lo(2));
  RealVect v1(m_rbox.hi(0), m_rbox.lo(1), m_rbox.lo(2));
  RealVect v2(m_rbox.lo(0), m_rbox.hi(1), m_rbox.lo(2));
  RealVect v3(m_rbox.lo(0), m_rbox.lo(1), m_rbox.hi(2));
  RealVect v4(m_rbox.hi(0), m_rbox.lo(1), m_rbox.hi(2));
  RealVect v5(m_rbox.hi(0), m_rbox.hi(1), m_rbox.lo(2));
  RealVect v6(m_rbox.lo(0), m_rbox.hi(1), m_rbox.hi(2));
  RealVect v7(m_rbox.hi(0), m_rbox.hi(1), m_rbox.hi(2));

#ifndef AMREX_USE_GPU
  constexpr bool print_F1 = false;
  constexpr bool print_F2 = false;
  constexpr bool print_F3 = false;
  constexpr bool print_F4 = false;
  constexpr bool print_F5 = false;
  constexpr bool print_F6 = false;
  constexpr bool print_F7 = false;
#endif

#ifdef AMREX_USE_GPU
  m_F1.add_vertex(v0);
  m_F4.add_vertex(v0);
  m_F5.add_vertex(v0);
#else
  m_F1.add_vertex(v0,"Initial: add v0 -> F1",print_F1);
  m_F4.add_vertex(v0,"Initial: add v0 -> F4",print_F4);
  m_F5.add_vertex(v0,"Initial: add v0 -> F5",print_F5);
#endif

  int add_v7(1);

  RealVect vIP;
  Real distIP;
  constexpr Real tol = 1.e-15;

  amrex::Array<bool, 8> vertex_intersected{}; // true if vertex has been intersected within a path.

  // Path 1
  {
    int cuts(0);
    for (int i = 0; i < 8; ++i) vertex_intersected[i] = false;

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v0, v1, vIP, distIP)) {

      amrex::Real length = (v0-v1).vectorLength();

      if ( Math::abs(distIP) < tol ) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v0);
        m_F4.add_vertex(v0);
        m_F5.add_vertex(v0);
        m_F7.add_vertex(v0);
#else
        m_F1.add_vertex(v0,"Path 1: v0--v1: add v0 -> F1",print_F1);
        m_F4.add_vertex(v0,"Path 1: v0--v1: add v0 -> F4",print_F4);
        m_F5.add_vertex(v0,"Path 1: v0--v1: add v0 -> F5",print_F5);
        m_F7.add_vertex(v0,"Path 1: v0--v1: add v0 -> F7",print_F7);
#endif
        if (!vertex_intersected[0]) {
          vertex_intersected[0] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol ) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v1);
        m_F2.add_vertex(v1);
        m_F5.add_vertex(v1);
        m_F7.add_vertex(v1);
#else
        m_F1.add_vertex(v1,"Path 1: v0--v1: add v1 -> F1",print_F1);
        m_F2.add_vertex(v1,"Path 1: v0--v1: add v1 -> F2",print_F2);
        m_F5.add_vertex(v1,"Path 1: v0--v1: add v1 -> F5",print_F5);
        m_F7.add_vertex(v1,"Path 1: v0--v1: add v1 -> F7",print_F7);
#endif
        if (!vertex_intersected[1]) {
          vertex_intersected[1] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(vIP);
        m_F5.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F1.add_vertex(vIP,"Path 1: v0--v1: add vIP -> F1",print_F1);
        m_F5.add_vertex(vIP,"Path 1: v0--v1: add vIP -> F5",print_F5);
        m_F7.add_vertex(vIP,"Path 1: v0--v1: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }

    } // Path 1: v0--v1

    if (cuts%2 == 0) {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(v1);
      m_F2.add_vertex(v1);
      m_F5.add_vertex(v1);
#else
      m_F1.add_vertex(v1,"Path 1: after v0--v1, cuts%2 == 0: add v1 -> F1",print_F1);
      m_F2.add_vertex(v1,"Path 1: after v0--v1, cuts%2 == 0: add v1 -> F2",print_F2);
      m_F5.add_vertex(v1,"Path 1: after v0--v1, cuts%2 == 0: add v1 -> F5",print_F5);
#endif
    }

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v1, v4, vIP, distIP)) {

      amrex::Real length = (v1-v4).vectorLength();

      if ( Math::abs(distIP) < tol ) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v1);
        m_F2.add_vertex(v1);
        m_F5.add_vertex(v1);
        m_F7.add_vertex(v1);
#else
        m_F1.add_vertex(v1,"Path 1: v1--v4: add v1 -> F1",print_F1);
        m_F2.add_vertex(v1,"Path 1: v1--v4: add v1 -> F2",print_F2);
        m_F5.add_vertex(v1,"Path 1: v1--v4: add v1 -> F5",print_F5);
        m_F7.add_vertex(v1,"Path 1: v1--v4: add v1 -> F7",print_F7);
#endif
        if (!vertex_intersected[1]) {
          vertex_intersected[1] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol ) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v4);
        m_F3.add_vertex(v4);
        m_F5.add_vertex(v4);
        m_F7.add_vertex(v4);
#else
        m_F2.add_vertex(v4,"Path 1: v1--v4: add v4 -> F2",print_F2);
        m_F3.add_vertex(v4,"Path 1: v1--v4: add v4 -> F3",print_F3);
        m_F5.add_vertex(v4,"Path 1: v1--v4: add v4 -> F5",print_F5);
        m_F7.add_vertex(v4,"Path 1: v1--v4: add v4 -> F7",print_F7);
#endif
        if (!vertex_intersected[4]) {
          vertex_intersected[4] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(vIP);
        m_F5.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F2.add_vertex(vIP,"Path 1: v1--v4: add vIP -> F2",print_F2);
        m_F5.add_vertex(vIP,"Path 1: v1--v4: add vIP -> F5",print_F5);
        m_F7.add_vertex(vIP,"Path 1: v1--v4: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }

    } // P1: v1--v4

    if (cuts%2 == 0) {
#ifdef AMREX_USE_GPU
      m_F2.add_vertex(v4);
      m_F3.add_vertex(v4);
      m_F5.add_vertex(v4);
#else
      m_F2.add_vertex(v4,"Path 1: after v1--v4, cuts%2 == 0: add v4 -> F2",print_F2);
      m_F3.add_vertex(v4,"Path 1: after v1--v4, cuts%2 == 0: add v4 -> F3",print_F3);
      m_F5.add_vertex(v4,"Path 1: after v1--v4, cuts%2 == 0: add v4 -> F5",print_F5);
#endif
    }

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v4, v7, vIP, distIP)) {

      amrex::Real length = (v4-v7).vectorLength();

      if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v4);
        m_F3.add_vertex(v4);
        m_F5.add_vertex(v4);
        m_F7.add_vertex(v4);
#else
        m_F2.add_vertex(v4,"Path 1: v4--v7: add v4 -> F2",print_F2);
        m_F3.add_vertex(v4,"Path 1: v4--v7: add v4 -> F3",print_F3);
        m_F5.add_vertex(v4,"Path 1: v4--v7: add v4 -> F5",print_F5);
        m_F7.add_vertex(v4,"Path 1: v4--v7: add v4 -> F7",print_F7);
#endif
        if (!vertex_intersected[4]) {
          vertex_intersected[4] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v7);
        m_F3.add_vertex(v7);
        m_F6.add_vertex(v7);
        m_F7.add_vertex(v7);
#else
        m_F2.add_vertex(v7,"Path 1: v4--v7: add v7 -> F2",print_F2);
        m_F3.add_vertex(v7,"Path 1: v4--v7: add v7 -> F3",print_F3);
        m_F6.add_vertex(v7,"Path 1: v4--v7: add v7 -> F6",print_F6);
        m_F7.add_vertex(v7,"Path 1: v4--v7: add v7 -> F7",print_F7);
#endif
        if (!vertex_intersected[7]) {
          vertex_intersected[7] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(vIP);
        m_F3.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F2.add_vertex(vIP,"Path 1: v4--v7: add vIP -> F2",print_F2);
        m_F3.add_vertex(vIP,"Path 1: v4--v7: add vIP -> F3",print_F3);
        m_F7.add_vertex(vIP,"Path 1: v4--v7: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }

    } // P1: v4--v7

    if (cuts == 2 && add_v7) {
      if (!utils::intersect_plane_edge(m_eb_point, m_eb_normal, v0, v7, vIP, distIP)) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v7);
        m_F3.add_vertex(v7);
        m_F6.add_vertex(v7);
#else
        m_F2.add_vertex(v7,"Path 1: after v4--v7, cuts == 2 && add_v7: add v7 -> F2",print_F2);
        m_F3.add_vertex(v7,"Path 1: after v4--v7, cuts == 2 && add_v7: add v7 -> F3",print_F3);
        m_F6.add_vertex(v7,"Path 1: after v4--v7, cuts == 2 && add_v7: add v7 -> F6",print_F6);
#endif
        add_v7 = 0;
      }
    }
  } // end Path 1

  // Path 4
  if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v1, v5, vIP, distIP)) {

    amrex::Real length = (v1-v5).vectorLength();

    if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(v1);
      m_F2.add_vertex(v1);
      m_F5.add_vertex(v1);
      m_F7.add_vertex(v1);
#else
      m_F1.add_vertex(v1,"Path 4: v1--v5: add v1 -> F1",print_F1);
      m_F2.add_vertex(v1,"Path 4: v1--v5: add v1 -> F2",print_F2);
      m_F5.add_vertex(v1,"Path 4: v1--v5: add v1 -> F5",print_F5);
      m_F7.add_vertex(v1,"Path 4: v1--v5: add v1 -> F7",print_F7);
#endif
    } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(v5);
      m_F2.add_vertex(v5);
      m_F6.add_vertex(v5);
      m_F7.add_vertex(v5);
#else
      m_F1.add_vertex(v5,"Path 4: v1--v5: add v5 -> F1",print_F1);
      m_F2.add_vertex(v5,"Path 4: v1--v5: add v5 -> F2",print_F2);
      m_F6.add_vertex(v5,"Path 4: v1--v5: add v5 -> F6",print_F6);
      m_F7.add_vertex(v5,"Path 4: v1--v5: add v5 -> F7",print_F7);
#endif
    } else {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(vIP);
      m_F2.add_vertex(vIP);
      m_F7.add_vertex(vIP);
#else
      m_F1.add_vertex(vIP,"Path 4: v1--v5: add vIP -> F1",print_F1);
      m_F2.add_vertex(vIP,"Path 4: v1--v5: add vIP -> F2",print_F2);
      m_F7.add_vertex(vIP,"Path 4: v1--v5: add vIP -> F7",print_F7);
#endif
    }
  }

  // Path 2
  { int cuts(0);
    for (int i = 0; i < 8; ++i) vertex_intersected[i] = false;

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v0, v2, vIP, distIP)) {

      amrex::Real length = (v0-v2).vectorLength();

      if ( Math::abs(distIP) < tol ) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v0);
        m_F4.add_vertex(v0);
        m_F5.add_vertex(v0);
        m_F7.add_vertex(v0);
#else
        m_F1.add_vertex(v0,"Path 2: v0--v2: add v0 -> F1",print_F1);
        m_F4.add_vertex(v0,"Path 2: v0--v2: add v0 -> F4",print_F4);
        m_F5.add_vertex(v0,"Path 2: v0--v2: add v0 -> F5",print_F5);
        m_F7.add_vertex(v0,"Path 2: v0--v2: add v0 -> F7",print_F7);
#endif
        if (!vertex_intersected[0]) {
          vertex_intersected[0] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol ) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v2);
        m_F4.add_vertex(v2);
        m_F6.add_vertex(v2);
        m_F7.add_vertex(v2);
#else
        m_F1.add_vertex(v2,"Path 2: v0--v2: add v2 -> F1",print_F1);
        m_F4.add_vertex(v2,"Path 2: v0--v2: add v2 -> F4",print_F4);
        m_F6.add_vertex(v2,"Path 2: v0--v2: add v2 -> F6",print_F6);
        m_F7.add_vertex(v2,"Path 2: v0--v2: add v2 -> F7",print_F7);
#endif
        if (!vertex_intersected[2]) {
          vertex_intersected[2] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(vIP);
        m_F4.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F1.add_vertex(vIP,"Path 2: v0--v2: add vIP -> F1",print_F1);
        m_F4.add_vertex(vIP,"Path 2: v0--v2: add vIP -> F4",print_F4);
        m_F7.add_vertex(vIP,"Path 2: v0--v2: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }

    } // P2: v0--v2

    if (cuts%2 == 0) {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(v2);
      m_F4.add_vertex(v2);
      m_F6.add_vertex(v2);
#else
      m_F1.add_vertex(v2,"Path 2: after v0--v2, cuts%2 == 0: add v2 -> F1",print_F1);
      m_F4.add_vertex(v2,"Path 2: after v0--v2, cuts%2 == 0: add v2 -> F4",print_F4);
      m_F6.add_vertex(v2,"Path 2: after v0--v2, cuts%2 == 0: add v2 -> F6",print_F6);
#endif
    }

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v2, v5, vIP, distIP)) {

      amrex::Real length = (v2-v5).vectorLength();

      if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v2);
        m_F4.add_vertex(v2);
        m_F6.add_vertex(v2);
        m_F7.add_vertex(v2);
#else
        m_F1.add_vertex(v2,"Path 2: v2--v5: add v2 -> F1",print_F1);
        m_F4.add_vertex(v2,"Path 2: v2--v5: add v2 -> F4",print_F4);
        m_F6.add_vertex(v2,"Path 2: v2--v5: add v2 -> F6",print_F6);
        m_F7.add_vertex(v2,"Path 2: v2--v5: add v2 -> F7",print_F7);
#endif
        if (!vertex_intersected[2]) {
          vertex_intersected[2] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v5);
        m_F2.add_vertex(v5);
        m_F6.add_vertex(v5);
        m_F7.add_vertex(v5);
#else
        m_F1.add_vertex(v5,"Path 2: v2--v5: add v5 -> F1",print_F1);
        m_F2.add_vertex(v5,"Path 2: v2--v5: add v5 -> F2",print_F2);
        m_F6.add_vertex(v5,"Path 2: v2--v5: add v5 -> F6",print_F6);
        m_F7.add_vertex(v5,"Path 2: v2--v5: add v5 -> F7",print_F7);
#endif
        if (!vertex_intersected[5]) {
          vertex_intersected[5] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(vIP);
        m_F6.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F1.add_vertex(vIP,"Path 2: v2--v5: add vIP -> F1",print_F1);
        m_F6.add_vertex(vIP,"Path 2: v2--v5: add vIP -> F6",print_F6);
        m_F7.add_vertex(vIP,"Path 2: v2--v5: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }
    }

    if (cuts%2 == 0) {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(v5);
      m_F2.add_vertex(v5);
      m_F6.add_vertex(v5);
#else
      m_F1.add_vertex(v5,"Path 2: after v2--v5, cuts%2 == 0: add v5 -> F1",print_F1);
      m_F2.add_vertex(v5,"Path 2: after v2--v5, cuts%2 == 0: add v5 -> F2",print_F2);
      m_F6.add_vertex(v5,"Path 2: after v2--v5, cuts%2 == 0: add v5 -> F6",print_F6);
#endif
    }

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v5, v7, vIP, distIP)) {

      amrex::Real length = (v5-v7).vectorLength();

      if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v5);
        m_F2.add_vertex(v5);
        m_F6.add_vertex(v5);
        m_F7.add_vertex(v5);
#else
        m_F1.add_vertex(v5,"Path 2: v5--v7: add v5 -> F1",print_F1);
        m_F2.add_vertex(v5,"Path 2: v5--v7: add v5 -> F2",print_F2);
        m_F6.add_vertex(v5,"Path 2: v5--v7: add v5 -> F6",print_F6);
        m_F7.add_vertex(v5,"Path 2: v5--v7: add v5 -> F7",print_F7);
#endif
        if (!vertex_intersected[5]) {
          vertex_intersected[5] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v7);
        m_F3.add_vertex(v7);
        m_F6.add_vertex(v7);
        m_F7.add_vertex(v7);
#else
        m_F2.add_vertex(v7,"Path 2: v5--v7: add v7 -> F2",print_F2);
        m_F3.add_vertex(v7,"Path 2: v5--v7: add v7 -> F3",print_F3);
        m_F6.add_vertex(v7,"Path 2: v5--v7: add v7 -> F6",print_F6);
        m_F7.add_vertex(v7,"Path 2: v5--v7: add v7 -> F7",print_F7);
#endif
        if (!vertex_intersected[7]) {
          vertex_intersected[7] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(vIP);
        m_F6.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F2.add_vertex(vIP,"Path 2: v5--v7: add vIP -> F2",print_F2);
        m_F6.add_vertex(vIP,"Path 2: v5--v7: add vIP -> F6",print_F6);
        m_F7.add_vertex(vIP,"Path 2: v5--v7: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }
    } // Path 2: v5--v7

    if (cuts == 2 && add_v7) {
      if (!utils::intersect_plane_edge(m_eb_point, m_eb_normal, v0, v7, vIP, distIP)) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v7);
        m_F3.add_vertex(v7);
        m_F6.add_vertex(v7);
#else
        m_F2.add_vertex(v7,"Path 2: after v5--v7, cuts == 2 && add_v7: add v7 -> F2",print_F2);
        m_F3.add_vertex(v7,"Path 2: after v5--v7, cuts == 2 && add_v7: add v7 -> F3",print_F3);
        m_F6.add_vertex(v7,"Path 2: after v5--v7, cuts == 2 && add_v7: add v7 -> F6",print_F6);
#endif
        add_v7 = 0;
      }
    }

  } // end Path 2

  // Path 5
  if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v2, v6, vIP, distIP)) {

    amrex::Real length = (v2-v6).vectorLength();

    if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
      m_F1.add_vertex(v2);
      m_F4.add_vertex(v2);
      m_F6.add_vertex(v2);
      m_F7.add_vertex(v2);
#else
      m_F1.add_vertex(v2,"Path 5: v2--v6: add v2 -> F1",print_F1);
      m_F4.add_vertex(v2,"Path 5: v2--v6: add v2 -> F4",print_F4);
      m_F6.add_vertex(v2,"Path 5: v2--v6: add v2 -> F6",print_F6);
      m_F7.add_vertex(v2,"Path 5: v2--v6: add v2 -> F7",print_F7);
#endif
    } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
      m_F3.add_vertex(v6);
      m_F4.add_vertex(v6);
      m_F6.add_vertex(v6);
      m_F7.add_vertex(v6);
#else
      m_F3.add_vertex(v6,"Path 5: v2--v6: add v6 -> F3",print_F3);
      m_F4.add_vertex(v6,"Path 5: v2--v6: add v6 -> F4",print_F4);
      m_F6.add_vertex(v6,"Path 5: v2--v6: add v6 -> F6",print_F6);
      m_F7.add_vertex(v6,"Path 5: v2--v6: add v6 -> F7",print_F7);
#endif
    } else {
#ifdef AMREX_USE_GPU
      m_F4.add_vertex(vIP);
      m_F6.add_vertex(vIP);
      m_F7.add_vertex(vIP);
#else
      m_F4.add_vertex(vIP,"Path 5: v2--v6: add vIP -> F4",print_F4);
      m_F6.add_vertex(vIP,"Path 5: v2--v6: add vIP -> F6",print_F6);
      m_F7.add_vertex(vIP,"Path 5: v2--v6: add vIP -> F7",print_F7);
#endif
    }
  } // end Path 5


  // Path 3
  { int cuts(0);
    for (int i = 0; i < 8; ++i) vertex_intersected[i] = false;

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v0, v3, vIP, distIP)) {

      amrex::Real length = (v0-v3).vectorLength();

      if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
        m_F1.add_vertex(v0);
        m_F4.add_vertex(v0);
        m_F5.add_vertex(v0);
        m_F7.add_vertex(v0);
#else
        m_F1.add_vertex(v0,"Path 3: v0--v3: add v0 -> F1",print_F1);
        m_F4.add_vertex(v0,"Path 3: v0--v3: add v0 -> F4",print_F4);
        m_F5.add_vertex(v0,"Path 3: v0--v3: add v0 -> F5",print_F5);
        m_F7.add_vertex(v0,"Path 3: v0--v3: add v0 -> F7",print_F7);
#endif
        if (!vertex_intersected[0]) {
          vertex_intersected[0] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
        m_F3.add_vertex(v3);
        m_F4.add_vertex(v3);
        m_F5.add_vertex(v3);
        m_F7.add_vertex(v3);
#else
        m_F3.add_vertex(v3,"Path 3: v0--v3: add v3 -> F3",print_F3);
        m_F4.add_vertex(v3,"Path 3: v0--v3: add v3 -> F4",print_F4);
        m_F5.add_vertex(v3,"Path 3: v0--v3: add v3 -> F5",print_F5);
        m_F7.add_vertex(v3,"Path 3: v0--v3: add v3 -> F7",print_F7);
#endif
        if (!vertex_intersected[3]) {
          vertex_intersected[3] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F4.add_vertex(vIP);
        m_F5.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F4.add_vertex(vIP,"Path 3: v0--v3: add vIP -> F4",print_F4);
        m_F5.add_vertex(vIP,"Path 3: v0--v3: add vIP -> F5",print_F5);
        m_F7.add_vertex(vIP,"Path 3: v0--v3: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }

    } // P3: v0--v3

    if (cuts%2 == 0) {
#ifdef AMREX_USE_GPU
      m_F3.add_vertex(v3);
      m_F4.add_vertex(v3);
      m_F5.add_vertex(v3);
#else
      m_F3.add_vertex(v3,"Path 3: after v0--v3, cuts%2 == 0: add v3 -> F3",print_F3);
      m_F4.add_vertex(v3,"Path 3: after v0--v3, cuts%2 == 0: add v3 -> F4",print_F4);
      m_F5.add_vertex(v3,"Path 3: after v0--v3, cuts%2 == 0: add v3 -> F5",print_F5);
#endif
    }

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v3, v6, vIP, distIP)) {

      amrex::Real length = (v3-v6).vectorLength();

      if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
        m_F3.add_vertex(v3);
        m_F4.add_vertex(v3);
        m_F5.add_vertex(v3);
        m_F7.add_vertex(v3);
#else
        m_F3.add_vertex(v3,"Path 3: v3--v6: add v3 -> F3",print_F3);
        m_F4.add_vertex(v3,"Path 3: v3--v6: add v3 -> F4",print_F4);
        m_F5.add_vertex(v3,"Path 3: v3--v6: add v3 -> F5",print_F5);
        m_F7.add_vertex(v3,"Path 3: v3--v6: add v3 -> F7",print_F7);
#endif
        if (!vertex_intersected[3]) {
          vertex_intersected[3] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
        m_F3.add_vertex(v6);
        m_F4.add_vertex(v6);
        m_F6.add_vertex(v6);
        m_F7.add_vertex(v6);
#else
        m_F3.add_vertex(v6,"Path 3: v3--v6: add v6 -> F3",print_F3);
        m_F4.add_vertex(v6,"Path 3: v3--v6: add v6 -> F4",print_F4);
        m_F6.add_vertex(v6,"Path 3: v3--v6: add v6 -> F6",print_F6);
        m_F7.add_vertex(v6,"Path 3: v3--v6: add v6 -> F7",print_F7);
#endif
        if (!vertex_intersected[6]) {
          vertex_intersected[6] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F3.add_vertex(vIP);
        m_F4.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F3.add_vertex(vIP,"Path 3: v3--v6: add vIP -> F3",print_F3);
        m_F4.add_vertex(vIP,"Path 3: v3--v6: add vIP -> F4",print_F4);
        m_F7.add_vertex(vIP,"Path 3: v3--v6: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }

    } // P3: v3--v6

    if (cuts%2 == 0) {
#ifdef AMREX_USE_GPU
      m_F3.add_vertex(v6);
      m_F4.add_vertex(v6);
      m_F6.add_vertex(v6);
#else
      m_F3.add_vertex(v6,"Path 3: after v3--v6, cuts%2 == 0: add v6 -> F3",print_F3);
      m_F4.add_vertex(v6,"Path 3: after v3--v6, cuts%2 == 0: add v6 -> F4",print_F4);
      m_F6.add_vertex(v6,"Path 3: after v3--v6, cuts%2 == 0: add v6 -> F6",print_F6);
#endif
    }

    if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v6, v7, vIP, distIP)) {

      amrex::Real length = (v6-v7).vectorLength();

      if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
        m_F3.add_vertex(v6);
        m_F4.add_vertex(v6);
        m_F6.add_vertex(v6);
        m_F7.add_vertex(v6);
#else
        m_F3.add_vertex(v6,"Path 3: v6--v7: add v6 -> F3",print_F3);
        m_F4.add_vertex(v6,"Path 3: v6--v7: add v6 -> F4",print_F4);
        m_F6.add_vertex(v6,"Path 3: v6--v7: add v6 -> F6",print_F6);
        m_F7.add_vertex(v6,"Path 3: v6--v7: add v6 -> F7",print_F7);
#endif
        if (!vertex_intersected[6]) {
          vertex_intersected[6] = true;
          ++cuts;
        }
      } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v7);
        m_F3.add_vertex(v7);
        m_F6.add_vertex(v7);
        m_F7.add_vertex(v7);
#else
        m_F2.add_vertex(v7,"Path 3: v6--v7: add v7 -> F2",print_F2);
        m_F3.add_vertex(v7,"Path 3: v6--v7: add v7 -> F3",print_F3);
        m_F6.add_vertex(v7,"Path 3: v6--v7: add v7 -> F6",print_F6);
        m_F7.add_vertex(v7,"Path 3: v6--v7: add v7 -> F7",print_F7);
#endif
        if (!vertex_intersected[7]) {
          vertex_intersected[7] = true;
          ++cuts;
        }
      } else {
#ifdef AMREX_USE_GPU
        m_F3.add_vertex(vIP);
        m_F6.add_vertex(vIP);
        m_F7.add_vertex(vIP);
#else
        m_F3.add_vertex(vIP,"Path 3: v6--v7: add vIP -> F3",print_F3);
        m_F6.add_vertex(vIP,"Path 3: v6--v7: add vIP -> F6",print_F6);
        m_F7.add_vertex(vIP,"Path 3: v6--v7: add vIP -> F7",print_F7);
#endif
        ++cuts;
      }
    } // P3: v6--v7

    if (cuts == 2 && add_v7) {
      if (!utils::intersect_plane_edge(m_eb_point, m_eb_normal, v0, v7, vIP, distIP)) {
#ifdef AMREX_USE_GPU
        m_F2.add_vertex(v7);
        m_F3.add_vertex(v7);
        m_F6.add_vertex(v7);
#else
        m_F2.add_vertex(v7,"Path 3: after v6--v7, cuts == 2 && add_v7: add v7 -> F2",print_F2);
        m_F3.add_vertex(v7,"Path 3: after v6--v7, cuts == 2 && add_v7: add v7 -> F3",print_F3);
        m_F6.add_vertex(v7,"Path 3: after v6--v7, cuts == 2 && add_v7: add v7 -> F6",print_F6);
#endif
        add_v7 = 0;
      }
    }

  } // end Path 3

  // Path 6
  if (utils::intersect_plane_edge(m_eb_point, m_eb_normal, v3, v4, vIP, distIP)) {

    amrex::Real length = (v3-v4).vectorLength();

    if ( Math::abs(distIP) < tol) {
#ifdef AMREX_USE_GPU
      m_F3.add_vertex(v3);
      m_F4.add_vertex(v3);
      m_F5.add_vertex(v3);
      m_F7.add_vertex(v3);
#else
      m_F3.add_vertex(v3,"Path 6: v3--v4: add v3 -> F3",print_F3);
      m_F4.add_vertex(v3,"Path 6: v3--v4: add v3 -> F4",print_F4);
      m_F5.add_vertex(v3,"Path 6: v3--v4: add v3 -> F5",print_F5);
      m_F7.add_vertex(v3,"Path 6: v3--v4: add v3 -> F7",print_F7);
#endif
    } else if ( Math::abs(distIP - length) < tol) {
#ifdef AMREX_USE_GPU
      m_F2.add_vertex(v4);
      m_F3.add_vertex(v4);
      m_F5.add_vertex(v4);
      m_F7.add_vertex(v4);
#else
      m_F2.add_vertex(v4,"Path 6: v3--v4: add v4 -> F2",print_F2);
      m_F3.add_vertex(v4,"Path 6: v3--v4: add v4 -> F3",print_F3);
      m_F5.add_vertex(v4,"Path 6: v3--v4: add v4 -> F5",print_F5);
      m_F7.add_vertex(v4,"Path 6: v3--v4: add v4 -> F7",print_F7);
#endif
    } else {
#ifdef AMREX_USE_GPU
      m_F3.add_vertex(vIP);
      m_F5.add_vertex(vIP);
      m_F7.add_vertex(vIP);
#else
      m_F3.add_vertex(vIP,"Path 6: v3--v4: add vIP -> F3",print_F3);
      m_F5.add_vertex(vIP,"Path 6: v3--v4: add vIP -> F5",print_F5);
      m_F7.add_vertex(vIP,"Path 6: v3--v4: add vIP -> F7",print_F7);
#endif
    }
  } // end Path 6

#ifndef AMREX_USE_GPU
  if (print_F1 || print_F2 || print_F3 || print_F4 || print_F5 || print_F6 || print_F7){
    amrex::Print()<<"\n";
  }
#endif

}



AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
eb_cut_cell_::
set_covered_regular_cell_vertices ()
{
  using namespace amrex;

  RealVect v0(m_rbox.lo(0), m_rbox.lo(1), m_rbox.lo(2));
  RealVect v1(m_rbox.hi(0), m_rbox.lo(1), m_rbox.lo(2));
  RealVect v2(m_rbox.lo(0), m_rbox.hi(1), m_rbox.lo(2));
  RealVect v3(m_rbox.lo(0), m_rbox.lo(1), m_rbox.hi(2));
  RealVect v4(m_rbox.hi(0), m_rbox.lo(1), m_rbox.hi(2));
  RealVect v5(m_rbox.hi(0), m_rbox.hi(1), m_rbox.lo(2));
  RealVect v6(m_rbox.lo(0), m_rbox.hi(1), m_rbox.hi(2));
  RealVect v7(m_rbox.hi(0), m_rbox.hi(1), m_rbox.hi(2));

  // Add vertices in the order of outward normal vector

  m_F1.add_vertex(v0);
  m_F1.add_vertex(v2);
  m_F1.add_vertex(v5);
  m_F1.add_vertex(v1);

  m_F2.add_vertex(v1);
  m_F2.add_vertex(v5);
  m_F2.add_vertex(v7);
  m_F2.add_vertex(v4);

  m_F3.add_vertex(v3);
  m_F3.add_vertex(v4);
  m_F3.add_vertex(v7);
  m_F3.add_vertex(v6);

  m_F4.add_vertex(v0);
  m_F4.add_vertex(v3);
  m_F4.add_vertex(v6);
  m_F4.add_vertex(v2);

  m_F5.add_vertex(v0);
  m_F5.add_vertex(v1);
  m_F5.add_vertex(v4);
  m_F5.add_vertex(v3);

  m_F6.add_vertex(v2);
  m_F6.add_vertex(v5);
  m_F6.add_vertex(v7);
  m_F6.add_vertex(v6);

}

#endif
