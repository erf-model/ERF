#include <ERF_FireGrid.H>
#include <AMReX_IntVect.H>

using namespace amrex;

FireGrid
create_fire_grid(const BoxArray& ba_atm,
                 const DistributionMapping& dm_atm,
                 const Geometry& geom_atm,
                 int C)
{
    FireGrid fg;
    fg.C = C;

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
    BoxArray ba_2d(box_list_2d);

    // Step 2: Refine the 2D BoxArray by {C, C, 1}
    IntVect ref_ratio(C, C, 1);
    BoxArray ba_fire = amrex::refine(ba_2d, ref_ratio);

    // Step 3: Refine DistributionMapping
    // Each MPI rank that owns box i in ba_atm now owns the corresponding
    // refined boxes in ba_fire
    DistributionMapping dm_fire = amrex::refine(dm_atm, ref_ratio);

    // Step 4: Create 2D Geometry
    // Keep x-y physical domain, drop z
    RealBox prob_domain_2d = geom_atm.ProbDomain();
    prob_domain_2d.setHi(2, prob_domain_2d.lo(2) + 1.0); // z extent = 1 m (dummy)

    Geometry geom_fire_2d(domain_2d, prob_domain_2d,
                          CoordSys::cartesian,
                          {false, false, false});

    fg.ba = ba_fire;
    fg.dm = dm_fire;
    fg.geom = geom_fire_2d;

    return fg;
}
