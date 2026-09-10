/**
 * \file ERF_EBIFTerrain.H
 * \brief Defines the implicit function used to build EB terrain geometry.
 */
#ifndef ERF_TERRAIN_IF_H_
#define ERF_TERRAIN_IF_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <cmath>
#include <algorithm>

//! AMReX EB implicit-function sign convention: >0 is body, 0 is boundary, <0 is fluid.

/**
 * \brief Terrain implicit function backed by gridded terrain heights.
 *
 * The operator interpolates terrain height in x-y, maps the query z coordinate
 * through stretched vertical spacing, and returns the AMReX EB signed value.
 */
class TerrainIF
    : amrex::GPUable
{
public:

    /**
     * \brief Construct a terrain implicit function from height and geometry data.
     * \param a_z_terrain Terrain elevation values on the horizontal mesh.
     * \param a_geom Geometry used for grid spacing and problem bounds.
     * \param a_dz_stretched Device vector of stretched vertical cell sizes.
     */
    TerrainIF (amrex::FArrayBox const& a_z_terrain, amrex::Geometry const& a_geom,
               amrex::Gpu::DeviceVector<amrex::Real>& a_dz_stretched)
        : terr_arr(a_z_terrain.const_array()),
          dz_s(a_dz_stretched.data()),
          dx(a_geom.CellSize(0)),
          dy(a_geom.CellSize(1)),
          dz(a_geom.CellSize(2)),
          prob_lo_x(a_geom.ProbLo(0)),
          prob_lo_y(a_geom.ProbLo(1)),
          prob_lo_z(a_geom.ProbLo(2)),
          prob_hi_x(a_geom.ProbHi(0)),
          prob_hi_y(a_geom.ProbHi(1)),
          i_hi(static_cast<int>(std::round((a_geom.ProbHi(0)-a_geom.ProbLo(0))/a_geom.CellSize(0)))),
          j_hi(static_cast<int>(std::round((a_geom.ProbHi(1)-a_geom.ProbLo(1))/a_geom.CellSize(1)))),
          k_hi(static_cast<int>(std::round((a_geom.ProbHi(2)-a_geom.ProbLo(2))/a_geom.CellSize(2))))
        {}

    /**
     * \brief Evaluate the signed implicit-function value at physical coordinates.
     *
     * The coordinates are passed through AMREX_D_DECL for spatial-dimension
     * portability.
     *
     * \return Positive inside terrain, zero at the terrain surface, negative in fluid.
     */
    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
        const noexcept
    {
        int i1{};
        int j1{};
        int i2{};
        int j2{};

        amrex::Real x1{};
        amrex::Real x2{};
        amrex::Real y1{};
        amrex::Real y2{};
        amrex::Real w1{};
        amrex::Real w2{};
        amrex::Real terr_z{amrex::Real(0.0)};

        // z_stretched
        const int k1 = amrex::Clamp(static_cast<int>(std::floor((z-prob_lo_z) / dz)), 0, k_hi - 1);
        const amrex::Real z1 = prob_lo_z + k1*dz;
        const amrex::Real remainder_z = (z - z1)/dz;
        amrex::Real z_stretched = prob_lo_z;
        for (int kk = 0; kk < k1; ++kk) {
            z_stretched += dz_s[kk];
        }
        z_stretched += remainder_z * dz_s[k1];

        // Interpolation (nine subregions of x-y plane)
        if        (x <= prob_lo_x && y <= prob_lo_y) {
            terr_z = terr_arr(0   ,0   ,0);
        } else if (x >= prob_hi_x && y <= prob_lo_y) {
            terr_z = terr_arr(i_hi,0   ,0);
        } else if (x >= prob_hi_x && y >= prob_hi_y) {
            terr_z = terr_arr(i_hi,j_hi,0);
        } else if (x <= prob_lo_x && y >= prob_hi_y) {
            terr_z = terr_arr(0   ,j_hi,0);
        } else if (x > prob_lo_x && x < prob_hi_x && (y <= prob_lo_y || y >= prob_hi_y) ) {
            i1 = static_cast<int>(std::floor((x-prob_lo_x) / dx));
            i2 = i1+1;
            x1 = prob_lo_x + i1*dx;
            x2 = x1 + dx;
            w1 = (x2-x)/dx;
            w2 = (x-x1)/dx;
            if        (y <= prob_lo_y){
                terr_z = w1*terr_arr(i1,0   ,0) + w2*terr_arr(i2,0   ,0);
            } else if (y >= prob_hi_y) {
                terr_z = w1*terr_arr(i1,j_hi,0) + w2*terr_arr(i2,j_hi,0);
            }
        } else if (y > prob_lo_y && y < prob_hi_y && (x <= prob_lo_x || x >= prob_hi_x) ) {
            j1 = static_cast<int>(std::floor((y-prob_lo_y) / dy));
            j2 = j1+1;
            y1 = prob_lo_y + j1*dy;
            y2 = y1 + dy;
            w1 = (y2-y)/dy;
            w2 = (y-y1)/dy;
            if        (x <= prob_lo_x){
                terr_z = w1*terr_arr(0   ,j1,0) + w2*terr_arr(0   ,j2,0);
            } else if (x >= prob_hi_x) {
                terr_z = w1*terr_arr(i_hi,j1,0) + w2*terr_arr(i_hi,j2,0);
            }
        } else {
            // Do bilinear interpolation of the terrain surface
            i1 = static_cast<int>(std::floor((x-prob_lo_x) / dx));
            i2 = i1+1;
            j1 = static_cast<int>(std::floor((y-prob_lo_y) / dy));
            j2 = j1+1;
            x1 = prob_lo_x + i1*dx;
            x2 = x1 + dx;
            y1 = prob_lo_y + j1*dy;
            y2 = y1 + dy;

            const amrex::Real denom = dx*dy;
            const amrex::Real w11 = (x2-x)*(y2-y)/denom;
            const amrex::Real w12 = (x2-x)*(y-y1)/denom;
            const amrex::Real w21 = (x-x1)*(y2-y)/denom;
            const amrex::Real w22 = (x-x1)*(y-y1)/denom;
            terr_z = w11*terr_arr(i1,j1,0) + w12*terr_arr(i1,j2,0) + w21*terr_arr(i2,j1,0) + w22*terr_arr(i2,j2,0);
        }
        return -(z_stretched - terr_z);
    }

    /**
     * \brief Evaluate the signed implicit-function value at a point array.
     * \param p Physical coordinate array.
     * \return Positive inside terrain, zero at the terrain surface, negative in fluid.
     */
    AMREX_GPU_HOST_DEVICE
    inline amrex::Real operator() (const amrex::RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:
    //! Terrain height array sampled on the horizontal mesh.
    amrex::Array4<amrex::Real const> terr_arr;
    //! Device pointer to stretched vertical cell sizes.
    amrex::Real const* dz_s;
    //! Uniform computational grid spacing in each direction.
    amrex::Real dx, dy, dz;
    //! Low problem bounds used for interpolation.
    amrex::Real prob_lo_x, prob_lo_y, prob_lo_z;
    //! High problem bounds used for interpolation.
    amrex::Real prob_hi_x, prob_hi_y;
    //! Highest terrain-array index in the x and y directions.
    int i_hi, j_hi;
    //! Number of vertical cells (one past the highest valid dz_s index).
    int k_hi;
    /**
     * \var TerrainIF::dy
     * \brief Uniform computational grid spacing in y.
     */
    /**
     * \var TerrainIF::dz
     * \brief Uniform computational grid spacing in z.
     */
    /**
     * \var TerrainIF::prob_lo_y
     * \brief Low problem bound in y.
     */
    /**
     * \var TerrainIF::prob_lo_z
     * \brief Low problem bound in z.
     */
    /**
     * \var TerrainIF::prob_hi_y
     * \brief High problem bound in y.
     */
    /**
     * \var TerrainIF::j_hi
     * \brief Highest terrain-array index in y.
     */
};

#endif
