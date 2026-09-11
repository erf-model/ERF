#include <cmath>

#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <gtest/gtest.h>

#include "ERF_TerrainMetrics.H"

// Motivation: with erf.terrain_smoothing = 0 (BTF) a fine level's mesh must not depend on
// how its BoxArray is split in z.  Each box maps its columns from its own lowest node, so a
// box stacked on another box of the level has to continue the column of the box below; it
// used to start again from its lowest node, which holds the mesh interpolated from the
// coarse level, and the fine mesh above a split followed the coarse terrain.

using namespace amrex;

namespace {

constexpr int nx = 16;
constexpr int ny = 16;
constexpr int nz = 12;

// Placeholder for nodes that the reference mesh does not cover
constexpr Real not_covered = Real(-1.0e30);

Geometry
make_geom (bool periodic)
{
    const Box domain(IntVect(0,0,0), IntVect(nx-1,ny-1,nz-1));
    const RealBox rb(Real(0.), Real(0.), Real(0.), Real(160.), Real(160.), Real(120.));
    const int per = periodic ? 1 : 0;
    return Geometry(domain, rb, CoordSys::cartesian, {per, per, 0});
}

// Stretched z levels, so that no level is a round number
Vector<Real>
make_z_levels ()
{
    Vector<Real> z_levels(nz+1);
    for (int k = 0; k <= nz; ++k) {
        z_levels[k] = Real(10.)*k + Real(0.5)*k*k;
    }
    return z_levels;
}

//
// The state a fine level is in when the BTF branch runs: the mesh interpolated from the
// coarse level everywhere, i.e. BTF over a smooth coarse terrain, and the fine terrain in
// the k = 0 slab.  Every node gets a value that depends only on its index, as it does in
// ERF, so any box that holds a node holds the same value there.
//
void
fill_interpolated_mesh (MultiFab& z_phys_nd, Vector<Real> const& z_levels)
{
    Gpu::DeviceVector<Real> z_levels_d(z_levels.size());
    Gpu::copy(Gpu::hostToDevice, z_levels.begin(), z_levels.end(), z_levels_d.begin());
    const Real* z_lev = z_levels_d.data();
    const Real z_top = z_levels[nz];

    for (MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        auto const z = z_phys_nd.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k == 0) {
                z(i,j,k) = Real(25.) + Real(12.)*std::sin(Real(0.31)*i + Real(0.1))*std::cos(Real(0.23)*j);
            } else {
                const int kk = amrex::max(0, amrex::min(k, nz));
                const Real h_coarse = Real(20.) + Real(10.)*std::sin(Real(0.3)*i)*std::cos(Real(0.2)*j);
                z(i,j,k) = z_lev[kk] + h_coarse * (z_top - z_lev[kk]) / z_top;
            }
        });
    }
    Gpu::streamSynchronize();
}

// The fine-level mesh on the cells of ba_cc, built by make_terrain_fitted_coords or, for
// per_box = true, by init_which_terrain_grid on each box on its own
MultiFab
build_fine_mesh (const BoxArray& ba_cc, const Geometry& geom, bool per_box = false)
{
    BoxArray ba_nd(ba_cc);
    ba_nd.surroundingNodes();
    DistributionMapping dm(ba_nd);
    MultiFab z_phys_nd(ba_nd, dm, 1, IntVect(3));

    const Vector<Real> z_levels = make_z_levels();
    fill_interpolated_mesh(z_phys_nd, z_levels);

    if (per_box) {
        z_phys_nd.setDomainBndry(bogus_large_value, 0, 1, geom);
        init_which_terrain_grid(1, geom, z_phys_nd, z_levels);
    } else {
        GpuArray<ERF_BC, AMREX_SPACEDIM*2> phys_bc_type;
        for (auto& bc : phys_bc_type) { bc = ERF_BC::slip_wall; }
        make_terrain_fitted_coords(1, geom, z_phys_nd, z_levels, phys_bc_type);
    }
    return z_phys_nd;
}

struct Mismatch
{
    Long valid = 0;   // valid nodes of the tested mesh that differ from the reference
    Long ghost = 0;   // ghost nodes that differ, where the reference covers them
    Real max_abs = Real(0.);
};

// Compare every node of test with the reference mesh, bit for bit
Mismatch
compare_meshes (MultiFab const& test, MultiFab const& ref)
{
    MultiFab ref_on_test(test.boxArray(), test.DistributionMap(), 1, test.nGrowVect());
    ref_on_test.setVal(not_covered);
    ref_on_test.ParallelCopy(ref, 0, 0, 1, ref.nGrowVect(), test.nGrowVect());
    // Where the reference holds a node as valid, compare with that value
    ref_on_test.ParallelCopy(ref, 0, 0, 1, IntVect(0), test.nGrowVect());

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpMax> reduce_op;
    ReduceData<Long, Long, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(test); mfi.isValid(); ++mfi) {
        const Box vbx = mfi.validbox();
        auto const a = test.const_array(mfi);
        auto const b = ref_on_test.const_array(mfi);
        reduce_op.eval(mfi.fabbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            if (b(i,j,k) == not_covered) { return {Long(0), Long(0), Real(0.)}; }
            const bool differs = (a(i,j,k) != b(i,j,k));
            const bool valid = vbx.contains(IntVect(i,j,k));
            return {Long(differs && valid), Long(differs && !valid), std::abs(a(i,j,k) - b(i,j,k))};
        });
    }

    auto const result = reduce_data.value(reduce_op);
    Mismatch m;
    m.valid   = amrex::get<0>(result);
    m.ghost   = amrex::get<1>(result);
    m.max_abs = amrex::get<2>(result);
    ParallelDescriptor::ReduceLongSum(m.valid);
    ParallelDescriptor::ReduceLongSum(m.ghost);
    ParallelDescriptor::ReduceRealMax(m.max_abs);
    return m;
}

void
expect_same_mesh (const BoxArray& ba_test, const BoxArray& ba_ref, bool periodic)
{
    const Geometry geom = make_geom(periodic);
    const MultiFab ref  = build_fine_mesh(ba_ref,  geom);
    const MultiFab test = build_fine_mesh(ba_test, geom);
    const Mismatch m = compare_meshes(test, ref);
    EXPECT_EQ(m.valid, 0) << "max |diff| " << m.max_abs << (periodic ? " (periodic)" : "");
    EXPECT_EQ(m.ghost, 0) << "max |diff| " << m.max_abs << (periodic ? " (periodic)" : "");
}

// A fine region touching the lateral domain boundary in x, over the full height
const Box full_region(IntVect(0,4,0), IntVect(7,11,nz-1));

} // namespace

TEST(TerrainFineColumns, SplitInZMatchesWholeColumns)
{
    BoxArray ba_ref(full_region);
    ba_ref.maxSize(IntVect(4,4,nz));

    BoxArray ba_split(full_region);
    ba_split.maxSize(IntVect(4,4,4));

    for (bool periodic : {false, true}) {
        expect_same_mesh(ba_split, ba_ref, periodic);
    }
}

TEST(TerrainFineColumns, UnevenSplitInZMatchesWholeColumns)
{
    BoxArray ba_ref(full_region);
    ba_ref.maxSize(IntVect(4,4,nz));

    BoxList bl;
    Box lower(full_region); lower.setBig(2, 4);
    Box upper(full_region); upper.setSmall(2, 5);
    bl.push_back(lower);
    bl.push_back(upper);
    BoxArray ba_split(std::move(bl));
    ba_split.maxSize(IntVect(4,4,nz));

    for (bool periodic : {false, true}) {
        expect_same_mesh(ba_split, ba_ref, periodic);
    }
}

// Boxes stacked with different lateral extents: the lower boxes are cut in x, the upper in y
TEST(TerrainFineColumns, MisalignedStackMatchesWholeColumns)
{
    BoxArray ba_ref(full_region);

    Box lower(full_region); lower.setBig(2, 5);
    Box upper(full_region); upper.setSmall(2, 6);
    BoxList bl_lower(lower); bl_lower.maxSize(IntVect(4,8,nz));
    BoxList bl_upper(upper); bl_upper.maxSize(IntVect(8,4,nz));
    bl_lower.join(bl_upper);
    const BoxArray ba_split(std::move(bl_lower));

    for (bool periodic : {false, true}) {
        expect_same_mesh(ba_split, ba_ref, periodic);
    }
}

// A patch aloft still starts from its own lowest node; split in z, it must give the same mesh
TEST(TerrainFineColumns, SplitPatchAloftMatchesWholePatch)
{
    Box aloft(full_region);
    aloft.setSmall(2, 3);

    const BoxArray ba_ref(aloft);

    BoxArray ba_split(aloft);
    ba_split.maxSize(IntVect(8,8,3));

    for (bool periodic : {false, true}) {
        expect_same_mesh(ba_split, ba_ref, periodic);
    }
}

// The per-box build is what the fine level used to get; it does depend on the split, which
// is what the tests above would catch
TEST(TerrainFineColumns, PerBoxBuildDependsOnSplit)
{
    const Geometry geom = make_geom(true);

    BoxArray ba_ref(full_region);
    ba_ref.maxSize(IntVect(4,4,nz));

    BoxArray ba_split(full_region);
    ba_split.maxSize(IntVect(4,4,4));

    const MultiFab ref  = build_fine_mesh(ba_ref,   geom);
    const MultiFab test = build_fine_mesh(ba_split, geom, true);
    EXPECT_GT(compare_meshes(test, ref).valid, 0);
}
