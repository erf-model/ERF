#ifndef ERF_EB_UTIL_H_
#define ERF_EB_UTIL_H_

#include <AMReX_RealVect.H>

namespace utils {

// Harmonic average
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real avg_h ( amrex::Real XXXm, amrex::Real XXXp) {
  if ( XXXm + XXXp < std::numeric_limits<amrex::Real>::min()) { return 0.; }
  else { return ((XXXm * XXXp) / (0.5*(XXXm + XXXp))); }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real trD (int i, int j, int k,
       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
       amrex::Array4<amrex::Real const> const& vel_x,
       amrex::Array4<amrex::Real const> const& vel_y,
       amrex::Array4<amrex::Real const> const& vel_z)
{
  return (vel_x(i+1,j,k)-vel_x(i,j,k))/dx[0] +
         (vel_y(i,j+1,k)-vel_y(i,j,k))/dx[1] +
         (vel_z(i,j,k+1)-vel_z(i,j,k))/dx[2];
}


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

  if ( amrex::Math::abs(np_dot_ne) < std::numeric_limits<amrex::Real>::min() )
  { return 0; }

  a_intersection_dist = a_plane_normal.dotProduct(a_plane_point - a_edge_point0);
  a_intersection_dist /= np_dot_ne;

  a_intersection_point = a_edge_point0 + a_intersection_dist*edge_normal;

  if (0. <= a_intersection_dist && a_intersection_dist <= edge_length) { return 1;}

  return 0;
}

}
#endif
