#ifndef ERF_EB_POLYGON_H_
#define ERF_EB_POLYGON_H_

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

#include "ERF_Constants.H"

class polygon_ {

  public:

    AMREX_GPU_HOST_DEVICE
    polygon_ ( amrex::RealVect a_point,
               amrex::RealVect a_normal )
      : m_cell_face(1)
      , m_eb_point(a_point)
      , m_eb_normal(a_normal)
      , m_defined(0)
      , m_num_vertices(0)
      , m_area(0.)
      , m_sorted(0)
      , m_vertices({amrex::RealVect(0.), amrex::RealVect(0.), amrex::RealVect(0.),
                  amrex::RealVect(0.), amrex::RealVect(0.), amrex::RealVect(0.)})
      , m_zdir(0.)
    {}

    AMREX_GPU_HOST_DEVICE
    polygon_ ( )
      : m_cell_face(0)
      , m_eb_point(0.)
      , m_eb_normal(0.)
      , m_defined(0)
      , m_num_vertices(0)
      , m_area(0.)
      , m_sorted(0)
      , m_vertices({amrex::RealVect(0.), amrex::RealVect(0.), amrex::RealVect(0.),
                  amrex::RealVect(0.), amrex::RealVect(0.), amrex::RealVect(0.)})
      , m_zdir(0.)
    {}

    AMREX_GPU_HOST_DEVICE
    void add_vertex ( amrex::RealVect const& a_v ) {
      for ( int i(0); i<m_num_vertices; ++i) {
        if ( amrex::almostEqual(m_vertices[i][0],a_v[0]) &&
             amrex::almostEqual(m_vertices[i][1],a_v[1]) &&
             amrex::almostEqual(m_vertices[i][2],a_v[2]) ) {
          return;
        }
      }
      AMREX_ASSERT( m_num_vertices < m_max_vertices );
      m_vertices[m_num_vertices] = a_v;
      m_num_vertices++;
    }

    AMREX_GPU_HOST_DEVICE
    int get_num_vertices ( ) {
      return m_num_vertices;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set_area ( amrex::Real const& a_area ) { m_area = a_area; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void define () {

      AMREX_ALWAYS_ASSERT( m_defined == 0 );  // TODO ---------------------------- remove ALWAYS

      m_defined = 1;

      // We need at least 3 vertices.
      if (m_num_vertices < 3) { return; }

      // Check to see if the vertices of the face are inside or outside of the plane.
      // If they are outside, this face doesn't belong to the volume.
      if ( m_cell_face ) {
        for ( int i(0); i<m_num_vertices; ++i) {
          if ((m_eb_normal.dotProduct(m_vertices[i] - m_eb_point)) >= 0.) {
          } else {
          }
        }
      }

      // Calculate the centroid of the polygon.
      amrex::RealVect centroid = get_centroid();

      // Shift the vertices relative to the centroid
      amrex::Array<amrex::RealVect,m_max_vertices> vertex_cent;
      for ( int i(0); i<m_num_vertices; ++i) {
        vertex_cent[i] = m_vertices[i] - centroid;
      }

      // Compute the normal vector m_zdir by cross product of two vectors in vertex_cent.
      // Choose the vectors with the largest cross-product magnitude.
      amrex::RealVect v_normal;
      amrex::Real max_norm2 = -1.0;
      const amrex::RealVect& vertex_cent_0 = vertex_cent[0];

      for (int i = 1; i < m_num_vertices; ++i)
      {
        amrex::RealVect vi = vertex_cent[i];
        amrex::RealVect v_cross = vertex_cent_0.crossProduct(vi);
        amrex::Real n2 = v_cross.radSquared();
        if (n2 > max_norm2)
        {
          max_norm2 = n2;
          v_normal = v_cross;
        }
      }
      if (!amrex::almostEqual(max_norm2, 0.0))
      {
          v_normal /= std::sqrt(max_norm2);
      }
      m_zdir = v_normal;

      //

      m_theta[0] = 0.0;
      for ( int i(1); i<m_num_vertices; ++i) {

        m_theta[i] = std::atan2(m_zdir.dotProduct( vertex_cent[0].crossProduct(vertex_cent[i]) ),
                              vertex_cent[0].dotProduct(vertex_cent[i]));

        m_theta[i] += ((m_theta[i] >= 0.) ? 0.0 : 2.0*PI);
      }

      // Sort counter clockwise based on theta.
      for (int i(0); i<m_num_vertices; ++i) {
        for (int j(0); j < m_num_vertices-i-1; ++j) {
          if ( m_theta[j] > m_theta[j+1] ) {
            amrex::Swap(m_theta[j], m_theta[j+1]);
            amrex::Swap(m_vertices[j], m_vertices[j+1]);
            amrex::Swap(vertex_cent[j], vertex_cent[j+1]);
          }
        } // j-loop
      } // i-loop
      m_sorted = 1;

      // Compute areas of triangles

      for (int i(0); i<m_num_vertices; ++i) {
        int const j( (i+1 == m_num_vertices) ? 0 : i+1 );
        amrex::RealVect vi_cross_vj = vertex_cent[i].crossProduct(vertex_cent[j]);
        m_area += 0.5*vi_cross_vj.vectorLength();
      }
      AMREX_ALWAYS_ASSERT( m_area > 0. ); // <------------------------------- TODO remove ALWAYS
    } // void define

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int ok ( ) const noexcept
    { return ((m_area > 0. || (m_area == 0. && m_defined == 1)) ? 1 : 0); }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real area ( ) const noexcept {
      AMREX_ALWAYS_ASSERT( ok() ); // <-------------------------------------- TODO remove ALWAYS
      return m_area;
    }

    // Distance from polygon to a_point
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real distance ( amrex::RealVect const& a_point ) const noexcept {
      AMREX_ALWAYS_ASSERT( m_defined == 1 );  // <--------------------------- TODO remove ALWAYS
      amrex::RealVect x0 = a_point - m_vertices[0];
      return amrex::Math::abs(x0.dotProduct(m_zdir));
    }

    // Centroid of the polygon based on the sub-triangulation
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect get_centroid ( ) const noexcept {
      amrex::RealVect cent(0.);
      if (m_num_vertices==3) {
        cent = (m_vertices[0] + m_vertices[1] + m_vertices[2]) / 3.0;
      } else {
        // Loop over sub-triangles
        amrex::Real area(0.);
        for ( int i(0); i<m_num_vertices-2; ++i) {
          amrex::RealVect const v0 (m_vertices[i+1] - m_vertices[0]);
          amrex::RealVect const v1 (m_vertices[i+2] - m_vertices[0]);
          amrex::RealVect v0_cross_v1 = v0.crossProduct(v1);
          amrex::Real area_tri = 0.5 * v0_cross_v1.vectorLength();
          amrex::RealVect cent_tri = (m_vertices[0] + m_vertices[i+1] + m_vertices[i+2]) / 3.0;
          area += area_tri;
          cent += area_tri * cent_tri;
        }
        cent = cent / area;
      }
      return cent;
    }

    // Unit normal vector
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect normal () const noexcept {
      return m_zdir;
    }

#ifndef AMREX_USE_GPU
    void report ( int const a_id, amrex::RealVect a_v0 )
#else
    void report ( int const /*a_id*/, amrex::RealVect /*a_v0*/ )
#endif
    {

#ifndef AMREX_USE_GPU
      amrex::RealVect centroid = get_centroid();
      amrex::Print() << "Face " << a_id
        << " -------------------------------------------\n"
        << "\nok? " << (ok() ? "yes" : "no") << "\n\n";
      for (int i(0); i<m_num_vertices; ++i) {
        amrex::Print() << "v" << i << ":  " << m_vertices[i] << '\n';
      }
      amrex::Print() << "\nvc: " << "  " << centroid << "\n";

      amrex::Real const dist = distance( a_v0 );
      // amrex::Real const vol  = ok() ? m_area * dist : 0.;

      amrex::Print() << "\narea:     " << m_area
                     << "\ndistance: " << dist
                     << "\nvolume:   " << dist*m_area
      << "\n==================================================\n\n";
#endif
    }

    void debug( int const a_id ) {
      amrex::Print() << "EBPolygon: id = " << a_id << ", m_num_vertices = " << m_num_vertices << '\n';
      for (int i(0); i<m_num_vertices; ++i) {
        amrex::Print() << "EBPolygon: v" << i << ":  " << m_vertices[i] << '\n';
      }
    }
    void debug() {
      debug(-1);
    }

  private:

    static int constexpr m_max_vertices = 6;

    int const m_cell_face;

    amrex::RealVect const m_eb_point;
    amrex::RealVect const m_eb_normal;

    int m_defined;

    int m_num_vertices;

    amrex::Real m_area;

    int m_sorted;

    amrex::Array<amrex::RealVect,m_max_vertices> m_vertices;

    amrex::GpuArray<amrex::Real,m_max_vertices> m_theta;

    amrex::RealVect m_zdir; // normal to polygon
};


#endif
