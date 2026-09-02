/**
 * @file ERF_DustGrid.cpp
 * @brief Implementation of the 2D dust emission grid factory function.
 */

#include <ERF_DustGrid.H>
#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_IntVect.H>

using namespace amrex;

DustGrid
create_dust_grid(const BoxArray& ba_atm,
                 const DistributionMapping& dm_atm,
                 const Geometry& geom_atm,
                 int grid_ratio)
{
    DustGrid dg;
    dg.grid_ratio = grid_ratio;

    // Step 1: Extract k=0 2D slice
    // Create a 2D BoxArray by taking the k=0 slice of each box
    Box domain_2d = geom_atm.Domain();
    domain_2d.setSmall(2, 0);
    domain_2d.setBig(2, 0);

    Vector<Box> box_list_2d;
    for (int i = 0; i < ba_atm.size(); ++i) {
        Box b = ba_atm[i];
        // Extract k=0 slice
        b.setSmall(2, 0);
        b.setBig(2, 0);
        box_list_2d.push_back(b);
    }
    BoxList bl_2d(std::move(box_list_2d));   // rvalue — matches Vector<Box>&&
    BoxArray ba_2d(bl_2d);


    // Step 2: Refine the 2D BoxArray by {grid_ratio, grid_ratio, 1}
    IntVect ref_ratio(grid_ratio, grid_ratio, 1);
    BoxArray ba_dust = amrex::refine(ba_2d, ref_ratio);

    // Step 3: DistributionMapping is unchanged — refine() on a BoxArray does not
    // change the number of boxes, only their size, so dm_atm maps identically
    // to ba_dust (box i in ba_dust is owned by the same rank as box i in ba_atm).
    DistributionMapping dm_dust = dm_atm;

    // Step 4: Create 2D Geometry with REFINED index-space domain
    // The index-space domain is scaled by grid_ratio so that cell size = 
    // physical_size / (grid_ratio * n_cells) = (physical_size / n_cells) / grid_ratio
    // Keep x-y physical domain unchanged — same extent as atmospheric grid.
    Box atm_domain   = geom_atm.Domain();
    Box atm_2d_full  = makeSlab(atm_domain, 2, 0);  // full atmospheric 2D domain
    // Scale hi-end: new hi = old_hi * grid_ratio + (grid_ratio-1) (zero-based indexing)
    Box dust_domain(atm_2d_full.smallEnd(),
                    IntVect(atm_2d_full.bigEnd(0) * grid_ratio + (grid_ratio - 1),
                            atm_2d_full.bigEnd(1) * grid_ratio + (grid_ratio - 1),
                            0));

    RealBox prob_domain_2d = geom_atm.ProbDomain();
    prob_domain_2d.setHi(2, prob_domain_2d.lo(2) + 1.0); // z extent = 1 m (dummy)

    // Inherit the atmospheric periodicity in x and y. Hard-coding this
    // non-periodic silently disables every FillBoundary on the dust grid: the
    // dust grid is one box spanning the domain, so all of its ghost cells are
    // domain-boundary ghosts, and a non-periodic Periodicity() leaves them
    // untouched. Anything that then reads a dust ghost cell picks up whatever
    // the allocator supplied. Same defect, and same fix, as ERF_FireGrid.cpp.
    // z is a dummy 1 m slab and stays non-periodic.
    Geometry geom_dust_2d(dust_domain, prob_domain_2d,
                          CoordSys::cartesian,
                          {geom_atm.isPeriodic(0), geom_atm.isPeriodic(1), 0});

    dg.ba = ba_dust;
    dg.dm = dm_dust;
    dg.geom = geom_dust_2d;

    return dg;
}
