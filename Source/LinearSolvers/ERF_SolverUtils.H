/**
 * \file ERF_SolverUtils.H
 */
#ifndef ERF_SOLVER_UTILS_H_
#define ERF_SOLVER_UTILS_H_

#include "AMReX_LO_BCTYPES.H"

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

#ifdef ERF_USE_FFT

/**
 * Define FFT boundary conditions on the specified solve box.
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
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;

    FFT::Boundary first, second;

    for (int dir = 0; dir <= 2; dir++)
    {
        const bool on_lo = (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir));
        const bool on_hi = (lev_geom.Domain().bigEnd(dir)   == bounding_box.bigEnd(dir));

        // A direction is periodic only if the region spans it.  A region that reaches one end
        // of a periodic direction only has coarse cells across the other end, so both of its
        // faces are Neumann faces (AMReX's FFT solvers need both ends periodic or neither).
        if (lev_geom.isPeriodic(dir) && on_lo && on_hi)
        {
            r[dir] = std::make_pair(FFT::Boundary::periodic, FFT::Boundary::periodic);
            continue;
        }

        auto bc_type_lo = l_domain_bc_type[Orientation(dir,Orientation::low)];
        if ( on_lo && !lev_geom.isPeriodic(dir) &&
             (bc_type_lo == "Outflow" || bc_type_lo == "Open") && (dir == 2 || !use_real_bcs) )
        {
            first = FFT::Boundary::odd;
        }
        else
        {
            first = FFT::Boundary::even;
        }

        auto bc_type_hi = l_domain_bc_type[Orientation(dir,Orientation::high)];
        if ( on_hi && !lev_geom.isPeriodic(dir) &&
             (bc_type_hi == "Outflow" || bc_type_hi == "Open") && (dir == 2 || !use_real_bcs) )
        {
            second = FFT::Boundary::odd;
        }
        else
        {
            second = FFT::Boundary::even;
        }

        r[dir] = std::make_pair(first,second);

    } // dir

    return r;
}
#endif
#endif
