/**
 * \file ERF_SolverUtils.H
 */
#ifndef ERF_SOLVER_UTILS_H_
#define ERF_SOLVER_UTILS_H_

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include "AMReX_LO_BCTYPES.H"
#include <string>
#include <utility>

#ifdef ERF_USE_FFT
#include <AMReX_FFT.H>
#endif

using namespace amrex;

/**
 * Define the domain boundary conditions for the (optional) Poisson solve
 * if we want to enforce that the initial conditions satisfy the constraint
 *
 * @param lev_geom Geometry used to identify periodic low-side boundaries
 * @param l_domain_bc_type Domain boundary-condition names
 * @return Low-side linear-operator boundary condition in each direction
 */
inline Array<LinOpBCType,AMREX_SPACEDIM>
get_lo_projection_bc (Geometry const& lev_geom,
                      Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type)
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if ( lev_geom.isPeriodic(dir) )
        {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = l_domain_bc_type[Orientation(dir,Orientation::low)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                r[dir] = LinOpBCType::Dirichlet;
            } else {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    return r;
}

/**
 * Define high-side domain boundary conditions for the projection solve.
 *
 * @param lev_geom Geometry used to identify periodic high-side boundaries
 * @param l_domain_bc_type Domain boundary-condition names
 * @return High-side linear-operator boundary condition in each direction
 */
inline Array<LinOpBCType,AMREX_SPACEDIM>
get_hi_projection_bc (Geometry const& lev_geom,
                      Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type)
{
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if ( lev_geom.isPeriodic(dir) )
        {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc_type = l_domain_bc_type[Orientation(dir,Orientation::high)];
            if (bc_type == "Outflow" || bc_type == "Open") {
                r[dir] = LinOpBCType::Dirichlet;
            } else {
                r[dir] = LinOpBCType::Neumann;
            }
        }
    }
    return r;
}

/**
 * Boundary type of the terrain Poisson operator on one face of the solve box:
 * periodic, even (homogeneous Neumann, ghost = interior) or odd (homogeneous
 * Dirichlet, ghost = -interior).  The FFT preconditioner and the multigrid
 * operator both derive their boundary conditions from this.
 */
enum struct TerrainBC : int { periodic = 0, even = 1, odd = 2 };

/**
 * Define the terrain Poisson boundary types on the specified solve box.
 *
 * A side of the solve box that is not a domain face is even (the solve is
 * closed by a zero-flux condition there); on a domain face Outflow and Open
 * give odd, except laterally when real boundary conditions are in use.
 *
 * @param lev_geom Geometry for the full level domain
 * @param l_domain_bc_type Domain boundary-condition names
 * @param bounding_box Box covered by the solve
 * @param use_real_bcs Whether real lateral boundary conditions are active
 * @return Low/high boundary-type pair in each direction
 */
inline Array<std::pair<TerrainBC,TerrainBC>,AMREX_SPACEDIM>
get_terrain_bc (Geometry const& lev_geom,
                Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
                Box const& bounding_box, bool use_real_bcs) noexcept
{
    Array<std::pair<TerrainBC,TerrainBC>,AMREX_SPACEDIM> r;

    TerrainBC first, second;

    for (int dir = 0; dir <= 2; dir++)
    {
        auto bc_type_lo = l_domain_bc_type[Orientation(dir,Orientation::low)];
        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) )
        {
            first = TerrainBC::periodic;
        }

        else if ( (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) &&
                  (bc_type_lo == "Outflow" || bc_type_lo == "Open") && (dir == 2 || !use_real_bcs) )
        {
            first = TerrainBC::odd;
        }

        else
        {
            first = TerrainBC::even;
        }

        auto bc_type_hi = l_domain_bc_type[Orientation(dir,Orientation::high)];
        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().bigEnd(dir) == bounding_box.bigEnd(dir)) )
        {
            second = TerrainBC::periodic;
        }

        else if ( (lev_geom.Domain().bigEnd(dir) == bounding_box.bigEnd(dir)) &&
                  (bc_type_hi == "Outflow" || bc_type_hi == "Open") && (dir == 2 || !use_real_bcs) )
        {
            second = TerrainBC::odd;
        }

        else
        {
            second = TerrainBC::even;
        }

        r[dir] = std::make_pair(first,second);

    } // dir

    return r;
}

/**
 * Linear-operator boundary conditions for the multigrid form of the terrain
 * Poisson operator, evaluated on the domain faces with the same rules as
 * get_terrain_bc: periodic where the level is periodic, Dirichlet for odd,
 * Neumann for even.  A box face inside the domain is a coarse/fine face and is
 * handled by MLLinOp::setCoarseFineBC.
 *
 * @param lev_geom Geometry for the full level domain
 * @param l_domain_bc_type Domain boundary-condition names
 * @param use_real_bcs Whether real lateral boundary conditions are active
 * @param bclo Low-side boundary condition in each direction (output)
 * @param bchi High-side boundary condition in each direction (output)
 */
inline void
get_terrain_projection_bc (Geometry const& lev_geom,
                           Array<std::string,2*AMREX_SPACEDIM> const& l_domain_bc_type,
                           bool use_real_bcs,
                           Array<LinOpBCType,AMREX_SPACEDIM>& bclo,
                           Array<LinOpBCType,AMREX_SPACEDIM>& bchi) noexcept
{
    auto tbc = get_terrain_bc(lev_geom, l_domain_bc_type, lev_geom.Domain(), use_real_bcs);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (lev_geom.isPeriodic(dir)) {
            bclo[dir] = LinOpBCType::Periodic;
            bchi[dir] = LinOpBCType::Periodic;
        } else {
            bclo[dir] = (tbc[dir].first  == TerrainBC::odd) ? LinOpBCType::Dirichlet : LinOpBCType::Neumann;
            bchi[dir] = (tbc[dir].second == TerrainBC::odd) ? LinOpBCType::Dirichlet : LinOpBCType::Neumann;
        }
    }
}

#ifdef ERF_USE_FFT

/**
 * Define FFT boundary conditions on the specified solve box (see get_terrain_bc).
 *
 * @param lev_geom Geometry for the full level domain
 * @param l_domain_bc_type Domain boundary-condition names
 * @param bounding_box Box covered by the FFT solve
 * @param use_real_bcs Whether real lateral boundary conditions are active
 * @return Low/high FFT boundary-condition pair in each direction
 */
inline Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
get_fft_bc (Geometry const& lev_geom,
            Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
            Box const& bounding_box, bool use_real_bcs) noexcept
{
    auto to_fft = [] (TerrainBC b) {
        return (b == TerrainBC::periodic) ? FFT::Boundary::periodic
             : (b == TerrainBC::odd)      ? FFT::Boundary::odd
                                          : FFT::Boundary::even;
    };
    auto tbc = get_terrain_bc(lev_geom, l_domain_bc_type, bounding_box, use_real_bcs);
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        r[dir] = std::make_pair(to_fft(tbc[dir].first), to_fft(tbc[dir].second));
    }
    return r;
}
#endif
#endif
