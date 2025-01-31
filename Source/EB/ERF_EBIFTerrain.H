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

    TerrainIF (amrex::FArrayBox const& a_z_terrain, amrex::Geometry const& a_geom,
               amrex::Gpu::DeviceVector<amrex::Real>& a_dz_stretched)
        : terr_arr(a_z_terrain.const_array()),
          dz_s(a_dz_stretched.data()),
          dx(a_geom.CellSize(0)),
          dy(a_geom.CellSize(1)),
          dz(a_geom.CellSize(2)),
          prob_lo_x(a_geom.ProbLo(0)),
          prob_lo_y(a_geom.ProbLo(1)),
          prob_lo_z(a_geom.ProbLo(2))
        {}

    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
        const noexcept
    {
        const int i1 = static_cast<int>(std::floor((x-prob_lo_x) / dx));
        const int j1 = static_cast<int>(std::floor((y-prob_lo_y) / dy));
        const int k1 = static_cast<int>(std::floor((z-prob_lo_z) / dz));

        const int i2 = i1+1;
        const int j2 = j1+1;

        amrex::Real x1 = prob_lo_x + i1*dx;
        amrex::Real y1 = prob_lo_y + j1*dy;
        amrex::Real z1 = prob_lo_z + k1*dz;

        amrex::Real x2 = prob_lo_x + i2*dx;
        amrex::Real y2 = prob_lo_y + j2*dy;

        amrex::Real remainder_z = (z - z1)/dz;

        amrex::Real z_stretched = prob_lo_z;

        for (int kk = 0; kk < k1; ++kk) {
            z_stretched += dz_s[kk];
        }
        z_stretched += remainder_z * dz_s[k1];

        // Do bilinear interpolation of the terrain surface
        amrex::Real denom = (x2-x1)*(y2-y1);
        amrex::Real  w_11 = (x2-x)*(y2-y)/denom;
        amrex::Real  w_12 = (x2-x)*(y-y1)/denom;
        amrex::Real  w_21 = (x-x1)*(y2-y)/denom;
        amrex::Real  w_22 = (x-x1)*(y-y1)/denom;

        amrex::Real terr_z = w_11*terr_arr(i1,j1,0) + w_12*terr_arr(i1,j2,0) + w_21*terr_arr(i2,j1,0) + w_22*terr_arr(i2,j2,0);

        return -(z_stretched - terr_z);
    }

    AMREX_GPU_HOST_DEVICE
    inline amrex::Real operator() (const amrex::RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:
    amrex::Array4<amrex::Real const> terr_arr;
    amrex::Real const* dz_s;
    amrex::Real dx, dy, dz;
    amrex::Real prob_lo_x, prob_lo_y, prob_lo_z;
};

#endif
