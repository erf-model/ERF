/** \addtogroup Utilities
 * @{
 */

#ifndef UTILS_H
#define UTILS_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

/**
 * Collection of utility functions
 *
 */

void
create_umac_grown (int lev, int nGrow, amrex::BoxArray& fine_grids,
                   const amrex::Geometry& fine_geom,
                   const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> u_mac_crse,
                   const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> u_mac_fine,
                   const amrex::IntVect& crse_ratio);

// A placeholder for filling dirichlet boundary conditions on u_mac.
struct FCFillExtDirDummy
{
    AMREX_GPU_HOST
    constexpr FCFillExtDirDummy( ) {}

    AMREX_GPU_DEVICE
    void operator()(
       const amrex::IntVect& /*iv*/,
       amrex::Array4<amrex::Real> const& /*dummy*/,
       const int /*dcomp*/,
       const int numcomp,
       amrex::GeometryData const& /*geom*/,
       const amrex::Real /*time*/,
       const amrex::BCRec* bcr,
       const int bcomp,
       const int /*orig_comp*/) const
    {
        // Abort if an calling this function was supposed to fill an ext_dir BC. 
        for (int n = bcomp; n < bcomp+numcomp; ++n) {
            const amrex::BCRec& bc = bcr[n];
            if ( AMREX_D_TERM(   bc.lo(0) == amrex::BCType::ext_dir || bc.hi(0) == amrex::BCType::ext_dir,
                              || bc.lo(1) == amrex::BCType::ext_dir || bc.hi(1) == amrex::BCType::ext_dir,
                              || bc.lo(2) == amrex::BCType::ext_dir || bc.hi(2) == amrex::BCType::ext_dir ) ) { 
               amrex::Abort("FCFillExtDirDummy is a dummy placeholder and should not be used in combination with BCType::ext_dir");
            }
        }
    }
};
#endif
/** @}*/
