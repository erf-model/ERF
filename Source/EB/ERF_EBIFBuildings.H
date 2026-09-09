/**
 * \file ERF_EBIFBuildings.H
 * \brief Defines the implicit function wrapper for STL-based building geometry.
 *
 * This file provides a BuildingsIF class that wraps AMReX's STL infrastructure
 * to enable 3D building geometry representation in ERF. Unlike the 2D height-based
 * TerrainIF approach, this provides accurate representation of vertical walls and
 * complex building structures.
 */

#ifndef ERF_BUILDINGS_IF_H_
#define ERF_BUILDINGS_IF_H_

#include <AMReX_EB_STL_utils.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>

#include <string>
#include <memory>

/**
 * \brief Buildings implicit function backed by STL triangle mesh.
 *
 * This class provides an implicit function interface for building geometry
 * loaded from STL files. It uses AMReX's STL infrastructure with BVH
 * (Bounding Volume Hierarchy) acceleration for efficient signed distance queries.
 *
 * The operator evaluates signed distance at any 3D point, following the
 * AMReX EB convention: positive = inside solid (building), zero = boundary,
 * negative = fluid.
 *
 * Note: This class pre-samples the STL geometry onto a reference grid and
 * interpolates at query points. For very large or complex geometries, consider
 * using AMReX's IndexSpaceSTL directly instead.
 */
class BuildingsIF
{
public:
    /**
     * \brief Construct a buildings implicit function from an STL file.
     *
     * \param stl_file Path to the STL mesh file (binary or ASCII format).
     * \param scale Uniform scaling factor applied to STL coordinates (typically 1.0).
     * \param center Translation vector to recenter the mesh [m].
     * \param reverse_normal Non-zero to flip triangle normals if needed.
     * \param geom Domain geometry for grid setup and BVH optimization.
     *
     * The constructor:
     * 1. Loads the STL file using AMReX::STLtools
     * 2. Builds a BVH acceleration structure for fast distance queries
     * 3. Pre-samples signed distances onto a reference grid
     * 4. Stores interpolation data for operator() calls
     */
    BuildingsIF(std::string const& stl_file,
                amrex::Real scale,
                amrex::Array<amrex::Real, 3> const& center,
                int reverse_normal,
                amrex::Geometry const& geom)
        : m_geom(geom)
    {
        using namespace amrex;

        // Load STL file with BVH optimization
        // This broadcasts triangle data to all ranks internally
        m_stl_tools.setBVHOptimization(true);
        m_stl_tools.read_stl_file(stl_file, scale, center, reverse_normal);

        // Pre-sample signed distance field onto a reference grid
        // This grid covers the entire domain with sufficient resolution
        // to capture building features (match computational grid)
        Box const& domain = geom.Domain();
        IntVect ngrow_vec(2); // Ghost cells for interpolation

        Box sampled_box = domain;
        sampled_box.grow(ngrow_vec);

        // Create distributed BoxArray and DistributionMapping for parallel computation
        // Chunk the domain into reasonable-sized boxes for load balancing
        BoxArray ba(sampled_box);
        ba.maxSize(32);
        DistributionMapping dm(ba);

        // Create distributed MultiFab to hold signed distance values
        // Each rank will compute its portion in parallel
        MultiFab mf_distributed(ba, dm, 1, ngrow_vec);

        // Compute signed distances in parallel across all ranks
        // fillSignedDistance uses ParallelFor internally, distributing work
        // across the MultiFab's boxes according to the DistributionMapping
        m_stl_tools.fillSignedDistance(mf_distributed, ngrow_vec, geom);

        // Now copy the distributed result to a single-box MultiFab
        // that contains the full domain on all ranks for operator() access
        BoxArray ba_single(sampled_box);
        DistributionMapping dm_single(ba_single);  // Single box on all ranks

        m_levelset_mf = std::make_shared<MultiFab>(ba_single, dm_single, 1, 0);

        // ParallelCopy from distributed to replicated
        // Signature: ParallelCopy(src, src_comp, dst_comp, num_comp, src_nghost, dst_nghost)
        IntVect dst_nghost(0);  // Destination has no ghost cells
        m_levelset_mf->ParallelCopy(mf_distributed, 0, 0, 1, ngrow_vec, dst_nghost);

        // Store bounds for interpolation
        m_data_box = sampled_box;

        // Set up local array view for operator() access
        m_local_array = (*m_levelset_mf)[0].const_array();
    }

    /**
     * \brief Evaluate signed distance at a point.
     *
     * \param x Physical x-coordinate [m]
     * \param y Physical y-coordinate [m]
     * \param z Physical z-coordinate [m]
     * \return Signed distance value: >0 inside building, <0 in fluid, =0 at boundary.
     *
     * This operator performs trilinear interpolation of the pre-sampled signed
     * distance field. For points outside the sampled domain, it returns a
     * negative value (fluid region).
     *
     * Note: The signed distance field is computed in parallel across ranks
     * during construction, then replicated for fast operator() access.
     */
    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator()(AMREX_D_DECL(amrex::Real x,
                                        amrex::Real y,
                                        amrex::Real z)) const noexcept
    {
        using namespace amrex;

        // Get domain bounds and cell size
        const Real* problo = m_geom.ProbLo();
        const Real* dx_arr = m_geom.CellSize();

        // Convert physical coordinates to cell indices (floating point)
        Real i_real = (x - problo[0]) / dx_arr[0];
        Real j_real = (y - problo[1]) / dx_arr[1];
        Real k_real = (z - problo[2]) / dx_arr[2];

        // Get integer cell indices
        int i = static_cast<int>(std::floor(i_real));
        int j = static_cast<int>(std::floor(j_real));
        int k = static_cast<int>(std::floor(k_real));

        // Check if point is within data bounds
        if (!m_data_box.contains(IntVect(AMREX_D_DECL(i, j, k))) ||
            !m_data_box.contains(IntVect(AMREX_D_DECL(i + 1, j + 1, k + 1)))) {
            // Outside sampled region - assume fluid (negative distance)
            return -1.0;
        }

        // Trilinear interpolation weights
        Real wx = i_real - static_cast<Real>(i);
        Real wy = j_real - static_cast<Real>(j);
        Real wz = k_real - static_cast<Real>(k);

        // Access replicated levelset data
        auto const& arr = m_local_array;

        // Trilinear interpolation
        Real phi_000 = arr(i, j, k);
        Real phi_100 = arr(i + 1, j, k);
        Real phi_010 = arr(i, j + 1, k);
        Real phi_110 = arr(i + 1, j + 1, k);
        Real phi_001 = arr(i, j, k + 1);
        Real phi_101 = arr(i + 1, j, k + 1);
        Real phi_011 = arr(i, j + 1, k + 1);
        Real phi_111 = arr(i + 1, j + 1, k + 1);

        Real phi_00 = phi_000 * (1.0 - wx) + phi_100 * wx;
        Real phi_01 = phi_001 * (1.0 - wx) + phi_101 * wx;
        Real phi_10 = phi_010 * (1.0 - wx) + phi_110 * wx;
        Real phi_11 = phi_011 * (1.0 - wx) + phi_111 * wx;

        Real phi_0 = phi_00 * (1.0 - wy) + phi_10 * wy;
        Real phi_1 = phi_01 * (1.0 - wy) + phi_11 * wy;

        Real phi = phi_0 * (1.0 - wz) + phi_1 * wz;

        return phi;
    }

    /**
     * \brief Evaluate signed distance at a point array.
     *
     * \param p Physical coordinate array [x, y, z] in meters.
     * \return Signed distance value.
     */
    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator()(const amrex::RealArray& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:
    //! Domain geometry (for coordinate transformations)
    amrex::Geometry m_geom;

    //! STL tools object (for loading and processing STL files)
    amrex::STLtools m_stl_tools;

    //! Pre-sampled signed distance field (distributed MultiFab for parallel computation)
    std::shared_ptr<amrex::MultiFab> m_levelset_mf;

    //! Bounding box of sampled data
    amrex::Box m_data_box;

    //! Local array view for GPU-compatible operator() access
    mutable amrex::Array4<const amrex::Real> m_local_array;
};

#endif // ERF_BUILDINGS_IF_H_
