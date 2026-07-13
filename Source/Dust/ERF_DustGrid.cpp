/**
 * @file ERF_DustGrid.cpp
 * @brief Implementation of the 2D dust emission grid factory function.
 */

#include <ERF_DustGrid.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>

DustGrid create_dust_grid(const amrex::BoxArray&          ba_atm,
                          const amrex::DistributionMapping& dm_atm,
                          const amrex::Geometry&          geom_atm,
                          int                             grid_ratio)
{
    // Step 1: Extract k=0 slab from each box in ba_atm
    amrex::BoxArray ba_2d = ba_atm;
    ba_2d.setSmall(2, 0);
    ba_2d.setBig(2, 0);

    // Step 2: Refine the 2D BoxArray by IntVect{grid_ratio, grid_ratio, 1}
    ba_2d.refine(amrex::IntVect{grid_ratio, grid_ratio, 1});

    // Step 3: Reuse dm_atm directly (refinement does not change the number of boxes)
    amrex::DistributionMapping dm_2d = dm_atm;

    // Step 4: Construct a refined 2D Geometry
    // Get physical domain from atmospheric geometry
    amrex::Real xlo = geom_atm.ProbLo(0);
    amrex::Real ylo = geom_atm.ProbLo(1);
    amrex::Real xhi = geom_atm.ProbHi(0);
    amrex::Real yhi = geom_atm.ProbHi(1);
    // z extent set to 1.0 m; the dust grid is a 2D slab and z is unused
    amrex::Real zlo = geom_atm.ProbLo(2);
    amrex::Real zhi = zlo + 1.0;

    amrex::RealBox prob_domain(xlo, ylo, zlo, xhi, yhi, zhi);
    amrex::Geometry geom_2d(ba_2d.getDomain(), prob_domain);

    DustGrid dg;
    dg.ba = ba_2d;
    dg.dm = dm_2d;
    dg.geom = geom_2d;
    dg.grid_ratio = grid_ratio;

    return dg;
}
