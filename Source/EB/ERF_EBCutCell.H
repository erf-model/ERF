#ifndef ERF_EB_CUT_CELL_H_
#define ERF_EB_CUT_CELL_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_RealBox.H>
#include <AMReX_RealVect.H>
#include <AMReX_Algorithm.H>
#include <AMReX_EBCellFlag.H>

#include "ERF_EBPolygon.H"

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
intersect_plane_edge ( amrex::RealVect const& a_plane_point,
                       amrex::RealVect const& a_plane_normal,
                       amrex::RealVect const& a_edge_point0,
                       amrex::RealVect const& a_edge_point1,
                       amrex::RealVect& a_intersection_point,
                       amrex::Real& a_intersection_dist )
{
  amrex::RealVect const edge(a_edge_point1 - a_edge_point0);
  amrex::Real const edge_length = edge.vectorLength();

  AMREX_ALWAYS_ASSERT(edge_length > 0.);

  amrex::RealVect edge_normal = edge / edge_length;

  amrex::Real np_dot_ne = a_plane_normal.dotProduct(edge_normal);

  // if ( amrex::Math::abs(np_dot_ne) < std::numeric_limits<amrex::Real>::min() )
  if ( amrex::Math::abs(np_dot_ne) < 10.0 * std::numeric_limits<amrex::Real>::epsilon() )
  { return 0; }

  a_intersection_dist = a_plane_normal.dotProduct(a_plane_point)
                      - a_plane_normal.dotProduct(a_edge_point0);
  a_intersection_dist /= np_dot_ne;

  a_intersection_point = a_edge_point0 + a_intersection_dist*edge_normal;

  if (0. <= a_intersection_dist && a_intersection_dist <= edge_length) { return 1;}

  return 0;
}

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

    struct path_data;

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

  amrex::RealVect v0(m_rbox.lo(0), m_rbox.lo(1), m_rbox.lo(2));
  amrex::RealVect v1(m_rbox.hi(0), m_rbox.lo(1), m_rbox.lo(2));
  amrex::RealVect v2(m_rbox.lo(0), m_rbox.hi(1), m_rbox.lo(2));
  amrex::RealVect v3(m_rbox.lo(0), m_rbox.lo(1), m_rbox.hi(2));
  amrex::RealVect v4(m_rbox.hi(0), m_rbox.lo(1), m_rbox.hi(2));
  amrex::RealVect v5(m_rbox.hi(0), m_rbox.hi(1), m_rbox.lo(2));
  amrex::RealVect v6(m_rbox.lo(0), m_rbox.hi(1), m_rbox.hi(2));
  amrex::RealVect v7(m_rbox.hi(0), m_rbox.hi(1), m_rbox.hi(2));

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

    amrex::RealVect c = 0.5*(v0 + v7);
    amrex::RealVect e = v7 - c;

    amrex::Real r = e[0]*amrex::Math::abs(a_normal[0]) +
             e[1]*amrex::Math::abs(a_normal[1]) +
             e[2]*amrex::Math::abs(a_normal[2]);

    amrex::Real s = amrex::Math::abs(c.dotProduct(a_normal)
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

struct eb_cut_cell_::path_data
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  path_data() = default;

  bool intersected{};
  bool intersected_start{};
  bool intersected_end{};
  amrex::RealVect vIP{};
  amrex::Real    distIP{};
  amrex::Real    edge_length{};

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void set(const amrex::RealVect& eb_point,
          const amrex::RealVect& eb_normal,
          const amrex::RealVect& v_start,
          const amrex::RealVect& v_end)
  {
    constexpr amrex::Real tol = 1.e-15;
    amrex::RealVect v_intersect;
    amrex::Real    dist_intersect;
    intersected = intersect_plane_edge(
      eb_point, eb_normal, v_start, v_end, v_intersect, dist_intersect
    );
    vIP = v_intersect;
    distIP = dist_intersect;
    edge_length = (v_start - v_end).vectorLength();
    // check if intersection is at vertices
    intersected_start = (amrex::Math::abs(dist_intersect) < tol);
    intersected_end   = (amrex::Math::abs(edge_length - dist_intersect) < tol);
  }
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
eb_cut_cell_::
calc_edge_intersections ()
{
  using namespace amrex;

  amrex::RealVect v0(m_rbox.lo(0), m_rbox.lo(1), m_rbox.lo(2));
  amrex::RealVect v1(m_rbox.hi(0), m_rbox.lo(1), m_rbox.lo(2));
  amrex::RealVect v2(m_rbox.lo(0), m_rbox.hi(1), m_rbox.lo(2));
  amrex::RealVect v3(m_rbox.lo(0), m_rbox.lo(1), m_rbox.hi(2));
  amrex::RealVect v4(m_rbox.hi(0), m_rbox.lo(1), m_rbox.hi(2));
  amrex::RealVect v5(m_rbox.hi(0), m_rbox.hi(1), m_rbox.lo(2));
  amrex::RealVect v6(m_rbox.lo(0), m_rbox.hi(1), m_rbox.hi(2));
  amrex::RealVect v7(m_rbox.hi(0), m_rbox.hi(1), m_rbox.hi(2));

// #ifndef AMREX_USE_GPU

  bool print_initial = ( Math::abs(m_eb_point[0]-1.666667e-01)<1.e-4 &&
                         Math::abs(m_eb_point[1]+4.194018e-01)<1.e-4 &&
                         Math::abs(m_eb_point[2]+1.666667e-01)<1.e-4 );
  const bool print_F1 = print_initial && false;
  const bool print_F2 = print_initial && false;
  const bool print_F3 = print_initial && false;
  const bool print_F4 = print_initial && false;
  const bool print_F5 = print_initial && false;
  const bool print_F6 = print_initial && false;
  const bool print_F7 = print_initial && false;
// #endif

  m_F1.add_vertex(v0); if(print_F1) printf("Initial: add v0 -> F1\n");
  m_F4.add_vertex(v0); if(print_F4) printf("Initial: add v0 -> F4\n");
  m_F5.add_vertex(v0); if(print_F5) printf("Initial: add v0 -> F5\n");

  int add_v7(1);

  amrex::RealVect vIP;
  amrex::Real distIP;

  amrex::Array<bool,8> vertex_intersected{}; // true if vertex has been intersected within a path.

  //------------------------------------------------------------
  // STEP 1: Predefine edge-plane intersections along 6 paths.
  // Check how the edges intersect the plane and whether the configuration is valid.
  // If the configuration is not valid, adjust it.
  //------------------------------------------------------------

  amrex::Array<path_data,3> p1, p2, p3;
  path_data p4, p5, p6;

  // Path 1: v0->v1->v4->v7
  p1[0].set(m_eb_point, m_eb_normal, v0, v1);
  p1[1].set(m_eb_point, m_eb_normal, v1, v4);
  p1[2].set(m_eb_point, m_eb_normal, v4, v7);

  // Path 2: v0->v2->v5->v7
  p2[0].set(m_eb_point, m_eb_normal, v0, v2);
  p2[1].set(m_eb_point, m_eb_normal, v2, v5);
  p2[2].set(m_eb_point, m_eb_normal, v5, v7);

  // Path 3: v0->v3->v6->v7
  p3[0].set(m_eb_point, m_eb_normal, v0, v3);
  p3[1].set(m_eb_point, m_eb_normal, v3, v6);
  p3[2].set(m_eb_point, m_eb_normal, v6, v7);

  // Path 4: v1->v5
  p4.set(m_eb_point, m_eb_normal, v1, v5);

  // Path 5: v2->v6
  p5.set(m_eb_point, m_eb_normal, v2, v6);

  // Path 6: v3->v4
  p6.set(m_eb_point, m_eb_normal, v3, v4);

  //------------------------------------------------------------
  // STEP 2: Add vertices to faces.
  //------------------------------------------------------------

  // Path 1
  {
    int cuts(0);
    for (int i = 0; i < 8; ++i) vertex_intersected[i] = false;

    if (p1[0].intersected) {
      if (p1[0].intersected_start) {
        m_F1.add_vertex(v0); if(print_F1) printf("Path 1: v0--v1: add v0 -> F1\n");
        m_F4.add_vertex(v0); if(print_F4) printf("Path 1: v0--v1: add v0 -> F4\n");
        m_F5.add_vertex(v0); if(print_F5) printf("Path 1: v0--v1: add v0 -> F5\n");
        m_F7.add_vertex(v0); if(print_F7) printf("Path 1: v0--v1: add v0 -> F7\n");
        if (!vertex_intersected[0]) {
          vertex_intersected[0] = true;
          ++cuts;
        }
      } else if (p1[0].intersected_end) {
        m_F1.add_vertex(v1); if(print_F1) printf("Path 1: v0--v1: add v1 -> F1\n");
        m_F2.add_vertex(v1); if(print_F2) printf("Path 1: v0--v1: add v1 -> F2\n");
        m_F5.add_vertex(v1); if(print_F5) printf("Path 1: v0--v1: add v1 -> F5\n");
        m_F7.add_vertex(v1); if(print_F7) printf("Path 1: v0--v1: add v1 -> F7\n");
        if (!vertex_intersected[1]) {
          vertex_intersected[1] = true;
          ++cuts;
        }
      } else {
        m_F1.add_vertex(p1[0].vIP); if(print_F1) printf("Path 1: v0--v1: add vIP -> F1\n");
        m_F5.add_vertex(p1[0].vIP); if(print_F5) printf("Path 1: v0--v1: add vIP -> F5\n");
        m_F7.add_vertex(p1[0].vIP); if(print_F7) printf("Path 1: v0--v1: add vIP -> F7\n");
        ++cuts;
      }

    } // Path 1: v0--v1

    if (cuts%2 == 0) {
      m_F1.add_vertex(v1); if(print_F1) printf("Path 1: after v0--v1, cuts-mod-2==0: add v1 -> F1\n");
      m_F2.add_vertex(v1); if(print_F2) printf("Path 1: after v0--v1, cuts-mod-2==0: add v1 -> F2\n");
      m_F5.add_vertex(v1); if(print_F5) printf("Path 1: after v0--v1, cuts-mod-2==0: add v1 -> F5\n");
    }

    if (p1[1].intersected) {
      if (p1[1].intersected_start) {
        m_F1.add_vertex(v1); if(print_F1) printf("Path 1: v1--v4: add v1 -> F1\n");
        m_F2.add_vertex(v1); if(print_F2) printf("Path 1: v1--v4: add v1 -> F2\n");
        m_F5.add_vertex(v1); if(print_F5) printf("Path 1: v1--v4: add v1 -> F5\n");
        m_F7.add_vertex(v1); if(print_F7) printf("Path 1: v1--v4: add v1 -> F7\n");
        if (!vertex_intersected[1]) {
          vertex_intersected[1] = true;
          ++cuts;
        }
      } else if (p1[1].intersected_end) {
        m_F2.add_vertex(v4); if(print_F2) printf("Path 1: v1--v4: add v4 -> F2\n");
        m_F3.add_vertex(v4); if(print_F3) printf("Path 1: v1--v4: add v4 -> F3\n");
        m_F5.add_vertex(v4); if(print_F5) printf("Path 1: v1--v4: add v4 -> F5\n");
        m_F7.add_vertex(v4); if(print_F7) printf("Path 1: v1--v4: add v4 -> F7\n");
        if (!vertex_intersected[4]) {
          vertex_intersected[4] = true;
          ++cuts;
        }
      } else {
        m_F2.add_vertex(p1[1].vIP); if(print_F2) printf("Path 1: v1--v4: add vIP -> F2\n");
        m_F5.add_vertex(p1[1].vIP); if(print_F5) printf("Path 1: v1--v4: add vIP -> F5\n");
        m_F7.add_vertex(p1[1].vIP); if(print_F7) printf("Path 1: v1--v4: add vIP -> F7\n");
        ++cuts;
      }

    } // P1: v1--v4

    if (cuts%2 == 0) {
      m_F2.add_vertex(v4); if(print_F2) printf("Path 1: after v1--v4, cuts-mod-2==0: add v4 -> F2\n");
      m_F3.add_vertex(v4); if(print_F3) printf("Path 1: after v1--v4, cuts-mod-2==0: add v4 -> F3\n");
      m_F5.add_vertex(v4); if(print_F5) printf("Path 1: after v1--v4, cuts-mod-2==0: add v4 -> F5\n");
    }

    if (p1[2].intersected) {
      if (p1[2].intersected_start) {
        m_F2.add_vertex(v4); if(print_F2) printf("Path 1: v4--v7: add v4 -> F2\n");
        m_F3.add_vertex(v4); if(print_F3) printf("Path 1: v4--v7: add v4 -> F3\n");
        m_F5.add_vertex(v4); if(print_F5) printf("Path 1: v4--v7: add v4 -> F5\n");
        m_F7.add_vertex(v4); if(print_F7) printf("Path 1: v4--v7: add v4 -> F7\n");
        if (!vertex_intersected[4]) {
          vertex_intersected[4] = true;
          ++cuts;
        }
      } else if (p1[2].intersected_end) {
        m_F2.add_vertex(v7); if(print_F2) printf("Path 1: v4--v7: add v7 -> F2\n");
        m_F3.add_vertex(v7); if(print_F3) printf("Path 1: v4--v7: add v7 -> F3\n");
        m_F6.add_vertex(v7); if(print_F6) printf("Path 1: v4--v7: add v7 -> F6\n");
        m_F7.add_vertex(v7); if(print_F7) printf("Path 1: v4--v7: add v7 -> F7\n");
        if (!vertex_intersected[7]) {
          vertex_intersected[7] = true;
          ++cuts;
        }
      } else {
        m_F2.add_vertex(p1[2].vIP); if(print_F2) printf("Path 1: v4--v7: add vIP -> F2\n");
        m_F3.add_vertex(p1[2].vIP); if(print_F3) printf("Path 1: v4--v7: add vIP -> F3\n");
        m_F7.add_vertex(p1[2].vIP); if(print_F7) printf("Path 1: v4--v7: add vIP -> F7\n");
        ++cuts;
      }

    } // P1: v4--v7

    if (cuts == 2 && add_v7) {
      if (!intersect_plane_edge(m_eb_point, m_eb_normal, v0, v7, vIP, distIP)) {
        m_F2.add_vertex(v7); if(print_F2) printf("Path 1: after v4--v7, cuts == 2 && add_v7: add v7 -> F2\n");
        m_F3.add_vertex(v7); if(print_F3) printf("Path 1: after v4--v7, cuts == 2 && add_v7: add v7 -> F3\n");
        m_F6.add_vertex(v7); if(print_F6) printf("Path 1: after v4--v7, cuts == 2 && add_v7: add v7 -> F6\n");
        add_v7 = 0;
      }
    }
  } // end Path 1

  // Path 4
  if (p4.intersected) {
    if (p4.intersected_start) {
      m_F1.add_vertex(v1); if(print_F1) printf("Path 4: v1--v5: add v1 -> F1\n");
      m_F2.add_vertex(v1); if(print_F2) printf("Path 4: v1--v5: add v1 -> F2\n");
      m_F5.add_vertex(v1); if(print_F5) printf("Path 4: v1--v5: add v1 -> F5\n");
      m_F7.add_vertex(v1); if(print_F7) printf("Path 4: v1--v5: add v1 -> F7\n");
    } else if (p4.intersected_end) {
      m_F1.add_vertex(v5); if(print_F1) printf("Path 4: v1--v5: add v5 -> F1\n");
      m_F2.add_vertex(v5); if(print_F2) printf("Path 4: v1--v5: add v5 -> F2\n");
      m_F6.add_vertex(v5); if(print_F6) printf("Path 4: v1--v5: add v5 -> F6\n");
      m_F7.add_vertex(v5); if(print_F7) printf("Path 4: v1--v5: add v5 -> F7\n");
    } else {
      m_F1.add_vertex(p4.vIP); if(print_F1) printf("Path 4: v1--v5: add vIP -> F1\n");
      m_F2.add_vertex(p4.vIP); if(print_F2) printf("Path 4: v1--v5: add vIP -> F2\n");
      m_F7.add_vertex(p4.vIP); if(print_F7) printf("Path 4: v1--v5: add vIP -> F7\n");
    }
  }

  // Path 2
  { int cuts(0);
    for (int i = 0; i < 8; ++i) vertex_intersected[i] = false;

    if (p2[0].intersected) {
      if (p2[0].intersected_start) {
        m_F1.add_vertex(v0); if(print_F1) printf("Path 2: v0--v2: add v0 -> F1\n");
        m_F4.add_vertex(v0); if(print_F4) printf("Path 2: v0--v2: add v0 -> F4\n");
        m_F5.add_vertex(v0); if(print_F5) printf("Path 2: v0--v2: add v0 -> F5\n");
        m_F7.add_vertex(v0); if(print_F7) printf("Path 2: v0--v2: add v0 -> F7\n");
        if (!vertex_intersected[0]) {
          vertex_intersected[0] = true;
          ++cuts;
        }
      } else if (p2[0].intersected_end) {
        m_F1.add_vertex(v2); if(print_F1) printf("Path 2: v0--v2: add v2 -> F1\n");
        m_F4.add_vertex(v2); if(print_F4) printf("Path 2: v0--v2: add v2 -> F4\n");
        m_F6.add_vertex(v2); if(print_F6) printf("Path 2: v0--v2: add v2 -> F6\n");
        m_F7.add_vertex(v2); if(print_F7) printf("Path 2: v0--v2: add v2 -> F7\n");
        if (!vertex_intersected[2]) {
          vertex_intersected[2] = true;
          ++cuts;
        }
      } else {
        m_F1.add_vertex(p2[0].vIP); if(print_F1) printf("Path 2: v0--v2: add vIP -> F1\n");
        m_F4.add_vertex(p2[0].vIP); if(print_F4) printf("Path 2: v0--v2: add vIP -> F4\n");
        m_F7.add_vertex(p2[0].vIP); if(print_F7) printf("Path 2: v0--v2: add vIP -> F7\n");
        ++cuts;
      }

    } // P2: v0--v2

    if (cuts%2 == 0) {
      m_F1.add_vertex(v2); if(print_F1) printf("Path 2: after v0--v2, cuts-mod-2==0: add v2 -> F1\n");
      m_F4.add_vertex(v2); if(print_F4) printf("Path 2: after v0--v2, cuts-mod-2==0: add v2 -> F4\n");
      m_F6.add_vertex(v2); if(print_F6) printf("Path 2: after v0--v2, cuts-mod-2==0: add v2 -> F6\n");
    }

    if (p2[1].intersected) {
      if (p2[1].intersected_start) {
        m_F1.add_vertex(v2); if(print_F1) printf("Path 2: v2--v5: add v2 -> F1\n");
        m_F4.add_vertex(v2); if(print_F4) printf("Path 2: v2--v5: add v2 -> F4\n");
        m_F6.add_vertex(v2); if(print_F6) printf("Path 2: v2--v5: add v2 -> F6\n");
        m_F7.add_vertex(v2); if(print_F7) printf("Path 2: v2--v5: add v2 -> F7\n");
        if (!vertex_intersected[2]) {
          vertex_intersected[2] = true;
          ++cuts;
        }
      } else if (p2[1].intersected_end) {
        m_F1.add_vertex(v5); if(print_F1) printf("Path 2: v2--v5: add v5 -> F1\n");
        m_F2.add_vertex(v5); if(print_F2) printf("Path 2: v2--v5: add v5 -> F2\n");
        m_F6.add_vertex(v5); if(print_F6) printf("Path 2: v2--v5: add v5 -> F6\n");
        m_F7.add_vertex(v5); if(print_F7) printf("Path 2: v2--v5: add v5 -> F7\n");
        if (!vertex_intersected[5]) {
          vertex_intersected[5] = true;
          ++cuts;
        }
      } else {
        m_F1.add_vertex(p2[1].vIP); if(print_F1) printf("Path 2: v2--v5: add vIP -> F1\n");
        m_F6.add_vertex(p2[1].vIP); if(print_F6) printf("Path 2: v2--v5: add vIP -> F6\n");
        m_F7.add_vertex(p2[1].vIP); if(print_F7) printf("Path 2: v2--v5: add vIP -> F7\n");
        ++cuts;
      }
    }

    if (cuts%2 == 0) {
      m_F1.add_vertex(v5); if(print_F1) printf("Path 2: after v2--v5, cuts-mod-2==0: add v5 -> F1\n");
      m_F2.add_vertex(v5); if(print_F2) printf("Path 2: after v2--v5, cuts-mod-2==0: add v5 -> F2\n");
      m_F6.add_vertex(v5); if(print_F6) printf("Path 2: after v2--v5, cuts-mod-2==0: add v5 -> F6\n");
    }

    if (p2[2].intersected) {
      if (p2[2].intersected_start) {
        m_F1.add_vertex(v5); if(print_F1) printf("Path 2: v5--v7: add v5 -> F1\n");
        m_F2.add_vertex(v5); if(print_F2) printf("Path 2: v5--v7: add v5 -> F2\n");
        m_F6.add_vertex(v5); if(print_F6) printf("Path 2: v5--v7: add v5 -> F6\n");
        m_F7.add_vertex(v5); if(print_F7) printf("Path 2: v5--v7: add v5 -> F7\n");
        if (!vertex_intersected[5]) {
          vertex_intersected[5] = true;
          ++cuts;
        }
      } else if (p2[2].intersected_end) {
        m_F2.add_vertex(v7); if(print_F2) printf("Path 2: v5--v7: add v7 -> F2\n");
        m_F3.add_vertex(v7); if(print_F3) printf("Path 2: v5--v7: add v7 -> F3\n");
        m_F6.add_vertex(v7); if(print_F6) printf("Path 2: v5--v7: add v7 -> F6\n");
        m_F7.add_vertex(v7); if(print_F7) printf("Path 2: v5--v7: add v7 -> F7\n");
        if (!vertex_intersected[7]) {
          vertex_intersected[7] = true;
          ++cuts;
        }
      } else {
        m_F2.add_vertex(p2[2].vIP); if(print_F2) printf("Path 2: v5--v7: add vIP -> F2\n");
        m_F6.add_vertex(p2[2].vIP); if(print_F6) printf("Path 2: v5--v7: add vIP -> F6\n");
        m_F7.add_vertex(p2[2].vIP); if(print_F7) printf("Path 2: v5--v7: add vIP -> F7\n");
        ++cuts;
      }
    } // Path 2: v5--v7

    if (cuts == 2 && add_v7) {
      if (!intersect_plane_edge(m_eb_point, m_eb_normal, v0, v7, vIP, distIP)) {
        m_F2.add_vertex(v7); if(print_F2) printf("Path 2: after v5--v7, cuts == 2 && add_v7: add v7 -> F2\n");
        m_F3.add_vertex(v7); if(print_F3) printf("Path 2: after v5--v7, cuts == 2 && add_v7: add v7 -> F3\n");
        m_F6.add_vertex(v7); if(print_F6) printf("Path 2: after v5--v7, cuts == 2 && add_v7: add v7 -> F6\n");
        add_v7 = 0;
      }
    }

  } // end Path 2

  // Path 5
  if (p5.intersected) {
    if (p5.intersected_start) {
      m_F1.add_vertex(v2); if(print_F1) printf("Path 5: v2--v6: add v2 -> F1\n");
      m_F4.add_vertex(v2); if(print_F4) printf("Path 5: v2--v6: add v2 -> F4\n");
      m_F6.add_vertex(v2); if(print_F6) printf("Path 5: v2--v6: add v2 -> F6\n");
      m_F7.add_vertex(v2); if(print_F7) printf("Path 5: v2--v6: add v2 -> F7\n");
    } else if (p5.intersected_end) {
      m_F3.add_vertex(v6); if(print_F3) printf("Path 5: v2--v6: add v6 -> F3\n");
      m_F4.add_vertex(v6); if(print_F4) printf("Path 5: v2--v6: add v6 -> F4\n");
      m_F6.add_vertex(v6); if(print_F6) printf("Path 5: v2--v6: add v6 -> F6\n");
      m_F7.add_vertex(v6); if(print_F7) printf("Path 5: v2--v6: add v6 -> F7\n");
    } else {
      m_F4.add_vertex(p5.vIP); if(print_F4) printf("Path 5: v2--v6: add vIP -> F4\n");
      m_F6.add_vertex(p5.vIP); if(print_F6) printf("Path 5: v2--v6: add vIP -> F6\n");
      m_F7.add_vertex(p5.vIP); if(print_F7) printf("Path 5: v2--v6: add vIP -> F7\n");
    }
  } // end Path 5


  // Path 3
  { int cuts(0);
    for (int i = 0; i < 8; ++i) vertex_intersected[i] = false;

    if (p3[0].intersected) {
      if (p3[0].intersected_start) {
        m_F1.add_vertex(v0); if(print_F1) printf("Path 3: v0--v3: add v0 -> F1\n");
        m_F4.add_vertex(v0); if(print_F4) printf("Path 3: v0--v3: add v0 -> F4\n");
        m_F5.add_vertex(v0); if(print_F5) printf("Path 3: v0--v3: add v0 -> F5\n");
        m_F7.add_vertex(v0); if(print_F7) printf("Path 3: v0--v3: add v0 -> F7\n");
        if (!vertex_intersected[0]) {
          vertex_intersected[0] = true;
          ++cuts;
        }
      } else if (p3[0].intersected_end) {
        m_F3.add_vertex(v3); if(print_F3) printf("Path 3: v0--v3: add v3 -> F3\n");
        m_F4.add_vertex(v3); if(print_F4) printf("Path 3: v0--v3: add v3 -> F4\n");
        m_F5.add_vertex(v3); if(print_F5) printf("Path 3: v0--v3: add v3 -> F5\n");
        m_F7.add_vertex(v3); if(print_F7) printf("Path 3: v0--v3: add v3 -> F7\n");
        if (!vertex_intersected[3]) {
          vertex_intersected[3] = true;
          ++cuts;
        }
      } else {
        m_F4.add_vertex(p3[0].vIP); if(print_F4) printf("Path 3: v0--v3: add vIP -> F4\n");
        m_F5.add_vertex(p3[0].vIP); if(print_F5) printf("Path 3: v0--v3: add vIP -> F5\n");
        m_F7.add_vertex(p3[0].vIP); if(print_F7) printf("Path 3: v0--v3: add vIP -> F7\n");
        ++cuts;
      }

    } // P3: v0--v3

    if (cuts%2 == 0) {
      m_F3.add_vertex(v3); if(print_F3) printf("Path 3: after v0--v3, cuts-mod-2 == 0: add v3 -> F3\n");
      m_F4.add_vertex(v3); if(print_F4) printf("Path 3: after v0--v3, cuts-mod-2 == 0: add v3 -> F4\n");
      m_F5.add_vertex(v3); if(print_F5) printf("Path 3: after v0--v3, cuts-mod-2 == 0: add v3 -> F5\n");
    }

    if (p3[1].intersected) {
      if (p3[1].intersected_start) {
        m_F3.add_vertex(v3); if(print_F3) printf("Path 3: v3--v6: add v3 -> F3\n");
        m_F4.add_vertex(v3); if(print_F4) printf("Path 3: v3--v6: add v3 -> F4\n");
        m_F5.add_vertex(v3); if(print_F5) printf("Path 3: v3--v6: add v3 -> F5\n");
        m_F7.add_vertex(v3); if(print_F7) printf("Path 3: v3--v6: add v3 -> F7\n");
        if (!vertex_intersected[3]) {
          vertex_intersected[3] = true;
          ++cuts;
        }
      } else if (p3[1].intersected_end) {
        m_F3.add_vertex(v6); if(print_F3) printf("Path 3: v3--v6: add v6 -> F3\n");
        m_F4.add_vertex(v6); if(print_F4) printf("Path 3: v3--v6: add v6 -> F4\n");
        m_F6.add_vertex(v6); if(print_F6) printf("Path 3: v3--v6: add v6 -> F6\n");
        m_F7.add_vertex(v6); if(print_F7) printf("Path 3: v3--v6: add v6 -> F7\n");
        if (!vertex_intersected[6]) {
          vertex_intersected[6] = true;
          ++cuts;
        }
      } else {
        m_F3.add_vertex(p3[1].vIP); if(print_F3) printf("Path 3: v3--v6: add vIP -> F3\n");
        m_F4.add_vertex(p3[1].vIP); if(print_F4) printf("Path 3: v3--v6: add vIP -> F4\n");
        m_F7.add_vertex(p3[1].vIP); if(print_F7) printf("Path 3: v3--v6: add vIP -> F7\n");
        ++cuts;
      }

    } // P3: v3--v6

    if (cuts%2 == 0) {
      m_F3.add_vertex(v6); if(print_F3) printf("Path 3: after v3--v6, cuts-mod-2 == 0: add v6 -> F3\n");
      m_F4.add_vertex(v6); if(print_F4) printf("Path 3: after v3--v6, cuts-mod-2 == 0: add v6 -> F4\n");
      m_F6.add_vertex(v6); if(print_F6) printf("Path 3: after v3--v6, cuts-mod-2 == 0: add v6 -> F6\n");
    }

    if (p3[2].intersected) {
      if (p3[2].intersected_start) {
        m_F3.add_vertex(v6); if(print_F3) printf("Path 3: v6--v7: add v6 -> F3\n");
        m_F4.add_vertex(v6); if(print_F4) printf("Path 3: v6--v7: add v6 -> F4\n");
        m_F6.add_vertex(v6); if(print_F6) printf("Path 3: v6--v7: add v6 -> F6\n");
        m_F7.add_vertex(v6); if(print_F7) printf("Path 3: v6--v7: add v6 -> F7\n");
        if (!vertex_intersected[6]) {
          vertex_intersected[6] = true;
          ++cuts;
        }
      } else if (p3[2].intersected_end) {
        m_F2.add_vertex(v7); if(print_F2) printf("Path 3: v6--v7: add v7 -> F2\n");
        m_F3.add_vertex(v7); if(print_F3) printf("Path 3: v6--v7: add v7 -> F3\n");
        m_F6.add_vertex(v7); if(print_F6) printf("Path 3: v6--v7: add v7 -> F6\n");
        m_F7.add_vertex(v7); if(print_F7) printf("Path 3: v6--v7: add v7 -> F7\n");
        if (!vertex_intersected[7]) {
          vertex_intersected[7] = true;
          ++cuts;
        }
      } else {
        m_F3.add_vertex(p3[2].vIP); if(print_F3) printf("Path 3: v6--v7: add vIP -> F3\n");
        m_F6.add_vertex(p3[2].vIP); if(print_F6) printf("Path 3: v6--v7: add vIP -> F6\n");
        m_F7.add_vertex(p3[2].vIP); if(print_F7) printf("Path 3: v6--v7: add vIP -> F7\n");
        ++cuts;
      }
    } // P3: v6--v7

    if (cuts == 2 && add_v7) {
      if (!intersect_plane_edge(m_eb_point, m_eb_normal, v0, v7, vIP, distIP)) {
        m_F2.add_vertex(v7); if(print_F2) printf("Path 3: after v6--v7, cuts == 2 && add_v7: add v7 -> F2\n");
        m_F3.add_vertex(v7); if(print_F3) printf("Path 3: after v6--v7, cuts == 2 && add_v7: add v7 -> F3\n");
        m_F6.add_vertex(v7); if(print_F6) printf("Path 3: after v6--v7, cuts == 2 && add_v7: add v7 -> F6\n");
        add_v7 = 0;
      }
    }

  } // end Path 3

  // Path 6
  if (p6.intersected) {
    if (p6.intersected_start) {
      m_F3.add_vertex(v3); if(print_F3) printf("Path 6: v3--v4: add v3 -> F3\n");
      m_F4.add_vertex(v3); if(print_F4) printf("Path 6: v3--v4: add v3 -> F4\n");
      m_F5.add_vertex(v3); if(print_F5) printf("Path 6: v3--v4: add v3 -> F5\n");
      m_F7.add_vertex(v3); if(print_F7) printf("Path 6: v3--v4: add v3 -> F7\n");
    } else if (p6.intersected_end) {
      m_F2.add_vertex(v4); if(print_F2) printf("Path 6: v3--v4: add v4 -> F2\n");
      m_F3.add_vertex(v4); if(print_F3) printf("Path 6: v3--v4: add v4 -> F3\n");
      m_F5.add_vertex(v4); if(print_F5) printf("Path 6: v3--v4: add v4 -> F5\n");
      m_F7.add_vertex(v4); if(print_F7) printf("Path 6: v3--v4: add v4 -> F7\n");
    } else {
      m_F3.add_vertex(p6.vIP); if(print_F3) printf("Path 6: v3--v4: add vIP -> F3\n");
      m_F5.add_vertex(p6.vIP); if(print_F5) printf("Path 6: v3--v4: add vIP -> F5\n");
      m_F7.add_vertex(p6.vIP); if(print_F7) printf("Path 6: v3--v4: add vIP -> F7\n");
    }
  } // end Path 6

  if (print_F1 || print_F2 || print_F3 || print_F4 || print_F5 || print_F6 || print_F7){
    printf("\n");
  }

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
