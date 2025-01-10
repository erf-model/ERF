#ifndef ERF_TERRAIN_IF_H_
#define ERF_TERRAIN_IF_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <cmath>
#include <algorithm>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class TerrainIF
    : amrex::GPUable
{
public:

    TerrainIF (amrex::FArrayBox& a_z_terrain, amrex::Geometry& a_geom)
        : terr_arr(a_z_terrain.const_array()),
          dx(a_geom.CellSize(0)),
          dy(a_geom.CellSize(1))
        {}

    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
        const noexcept
    {
        const int i = static_cast<int>(x / dx);
        const int j = static_cast<int>(y / dy);

        return -(z - terr_arr(i,j,0));
    }

    AMREX_GPU_HOST_DEVICE
    inline amrex::Real operator() (const amrex::RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:
    amrex::Array4<amrex::Real const> terr_arr;
    amrex::Real dx, dy;
};

#endif
