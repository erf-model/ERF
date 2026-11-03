#include <gtest/gtest.h>
#include <limits>

#include <ERFPCParticleToMesh.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParGDB.H>

using namespace amrex;

TEST(SpeciesProjection, CompositeCICConservesParticleMassAtWallsAndLevelEdges)
{
    const Box coarse_domain(IntVect(0), IntVect(3));
    const RealBox rb(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4));
    const Vector<Geometry> geom{
        Geometry(coarse_domain, rb, CoordSys::cartesian, {1, 0, 0}),
        Geometry(refine(coarse_domain, 2), rb, CoordSys::cartesian, {1, 0, 0})};
    Box fine_box = geom[1].Domain();
    fine_box.setBig(2, 3);
    const Vector<BoxArray> ba{BoxArray(coarse_domain), BoxArray(fine_box)};
    const Vector<DistributionMapping> dm{DistributionMapping(ba[0]), DistributionMapping(ba[1])};
    const Vector<IntVect> rr{IntVect(2)};
    ParGDB gdb(geom, dm, ba, rr);

    // Fine wall, fine-side interface, coarse-side interface; x also crosses
    // the periodic seam. No advection, microphysics, or plotting can repair this deposit.
    for (Real zpos : {Real(0.1), Real(1.9), Real(2.1)}) {
        ERFPC pc(&gdb, "species_projection_test", ERFPCOptions{false, false});
        pc.InitRandom(1, 42, {});
        for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
            for (ERFPC::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
                auto* particles = pti.GetArrayOfStructs().dataPtr();
                ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept {
                    particles[i].pos(0) = ParticleReal(0.1);
                    particles[i].pos(1) = ParticleReal(1.5);
                    particles[i].pos(2) = ParticleReal(zpos);
                });
            }
        }
        pc.Redistribute();
        Real mass = Real(0);
        for (int lev = 0; lev < 2; ++lev) {
            MultiFab deposit(ba[lev], dm[lev], 1, 1);
            BoxArray nodes(ba[lev]);
            nodes.surroundingNodes();
            MultiFab height(nodes, dm[lev], 1, 2);
            const Real dz = geom[lev].CellSize(2);
            for (MFIter mfi(height); mfi.isValid(); ++mfi) {
                const auto h = height.array(mfi);
                ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    h(i,j,k) = Real(k) * dz;
                });
            }
            pc.ERFPCParticleToMesh(deposit, height, lev, 0,
                [=] AMREX_GPU_DEVICE (const ERFPC::ParticleTileType::ConstParticleTileDataType&, int) {
                    return ParticleReal(1);
                }, true);
            if (lev == 0) {
                const iMultiFab mask = makeFineMask(deposit, ba[1], rr[0], 1, 0);
                for (MFIter mfi(deposit); mfi.isValid(); ++mfi) {
                    const auto d = deposit.array(mfi);
                    const auto m = mask.const_array(mfi);
                    ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (!m(i,j,k)) { d(i,j,k) = Real(0); }
                    });
                }
            }
            const auto dx = geom[lev].CellSizeArray();
            mass += deposit.sum(0) * dx[0] * dx[1] * dx[2];
        }
        EXPECT_NEAR(mass, Real(1), Real(64) * std::numeric_limits<Real>::epsilon())
            << "particle z=" << zpos;
    }
}
