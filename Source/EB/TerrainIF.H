#ifndef TERRAIN_IF_H_
#define TERRAIN_IF_H_

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
        : m_terr(a_z_terrain),
          m_geom(a_geom)
        {
            amrex::Print() << " EB type = Terrain " << std::endl;
        }

    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
        const noexcept
    {
        amrex::Real dx = m_geom.CellSizeArray()[0];
        amrex::Real dy = m_geom.CellSizeArray()[1];

        const int i = static_cast<int>(x / dx);
        const int j = static_cast<int>(y / dy);

        amrex::Array4<amrex::Real const> const& terr_arr = m_terr.const_array();

        if (i == 128 and j == 3) amrex::Print() <<" TERRARR " << amrex::IntVect(i,j,0) << " " << z << " " << terr_arr(i,j,0) << std::endl;

        return -(z - terr_arr(i,j,0));
    }

    inline amrex::Real operator() (const amrex::RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:
    amrex::FArrayBox&  m_terr;
    amrex::Geometry&   m_geom;
};

#endif
