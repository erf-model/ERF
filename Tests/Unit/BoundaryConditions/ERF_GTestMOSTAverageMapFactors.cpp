#include <AMReX_FabArrayBase.H>
#include <AMReX_Gpu.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include "ERF_GTestMOSTAverageCommon.H"

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_MOSTAverage.H>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>

using namespace amrex;

using erf_most_average_test::tolerance;
using erf_most_average_test::test_theta;
using erf_most_average_test::test_u;
using erf_most_average_test::test_v;
using erf_most_average_test::test_w;

// With a uniform map factor m, a mesh of computational spacing m*Delta has the physical
// spacing Delta of an m = 1 mesh. On the same physical terrain, the MOST average must find
// the same reference cell, the same reference point in cells, the same rotated velocities
// and the same subgrid velocity scale on both meshes (the invariance of erf-model/ERF#5366).

namespace {

constexpr int  ncell   = 16;
constexpr Real mfac    = Real(2.0);   // uniform map factor of the coarse computational mesh
constexpr Real delta   = Real(10.0);  // physical horizontal spacing of both meshes
constexpr Real dz      = Real(10.0);
constexpr Real slope_x = Real(0.5);   // physical terrain slopes dh/dx and dh/dy
constexpr Real slope_y = Real(0.25);
constexpr Real zref    = Real(40.0);

// Exposes the protected index, position, rotation and subgrid-scale data
class MOSTAverageProbe : public MOSTAverage
{
public:
    using MOSTAverage::MOSTAverage;

    [[nodiscard]] const iMultiFab& i_indices () const { return *m_i_indx[0]; }
    [[nodiscard]] const iMultiFab& j_indices () const { return *m_j_indx[0]; }
    [[nodiscard]] const iMultiFab& k_indices () const { return *m_k_indx[0]; }
    [[nodiscard]] const MultiFab& x_positions () const { return *m_x_pos[0]; }
    [[nodiscard]] const MultiFab& y_positions () const { return *m_y_pos[0]; }
    [[nodiscard]] const MultiFab& z_positions () const { return *m_z_pos[0]; }
    [[nodiscard]] const MultiFab& rotated (const int comp) const { return *m_rot_fields[0][comp]; }
    [[nodiscard]] Real subgrid_velocity_scale () const { return m_Vsg[0]; }
};

class ScopedMFIterTileSize
{
public:
    explicit ScopedMFIterTileSize (const IntVect& tile_size)
        : m_saved(FabArrayBase::mfiter_tile_size)
    {
        FabArrayBase::mfiter_tile_size = tile_size;
    }

    ~ScopedMFIterTileSize ()
    {
        FabArrayBase::mfiter_tile_size = m_saved;
    }

private:
    IntVect m_saved;
};

class ScopedMapFacMOSTParams
{
public:
    explicit ScopedMapFacMOSTParams (const char* prefix)
        : m_pp(prefix)
    {
        m_pp.add("most.average_policy", 0);
        m_pp.add("most.radius", 0);
    }

    void add (const char* name, const int value) { m_pp.add(name, value); }
    void add (const char* name, const Real value) { m_pp.add(name, value); }

    ~ScopedMapFacMOSTParams ()
    {
        for (const auto* name : m_names) {
            m_pp.remove(name);
        }
    }

private:
    ParmParse m_pp;
    const std::array<const char*, 7> m_names{{
        "most.average_policy",
        "most.radius",
        "most.zref",
        "most.terrain_rotate",
        "most.use_interpolation",
        "most.use_normal_vector",
        "most.include_subgrid_vel"
    }};
};

// A level with computational spacing h_dx in x and y, uniform map factor m (none if m <= 0),
// and the planar terrain z = slope_x x + slope_y y in physical coordinates. The velocities
// carry ng_vel ghost cells and the cell-centered fields and heights ng_vel+1, as in ERF
struct MapFacLevel
{
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab cons;
    MultiFab xvel;
    MultiFab yvel;
    MultiFab zvel;
    std::unique_ptr<MultiFab> theta;
    std::unique_ptr<MultiFab> qv;
    std::unique_ptr<MultiFab> qr;
    std::unique_ptr<MultiFab> z_phys_nd;
    Vector<std::unique_ptr<MultiFab>> mapfac;
    Vector<MultiFab*> vars_old;

    MapFacLevel (const Real h_dx, const Real m, const int ng_vel = 3)
        : geom(make_geometry(h_dx)),
          ba(geom.Domain()),
          dm(ba),
          cons(ba, dm, 1, ng_vel+1),
          xvel(convert(ba, IntVect(1,0,0)), dm, 1, ng_vel),
          yvel(convert(ba, IntVect(0,1,0)), dm, 1, ng_vel),
          zvel(convert(ba, IntVect(0,0,1)), dm, 1, ng_vel),
          theta(std::make_unique<MultiFab>(ba, dm, 1, ng_vel+1)),
          qv(std::make_unique<MultiFab>(ba, dm, 1, ng_vel+1)),
          z_phys_nd(std::make_unique<MultiFab>(convert(ba, IntVect(1)), dm, 1, ng_vel+1))
    {
        cons.setVal(one);
        xvel.setVal(test_u);
        yvel.setVal(test_v);
        zvel.setVal(test_w);
        theta->setVal(test_theta);
        qv->setVal(Real(0.01));
        vars_old = {&cons, &xvel, &yvel, &zvel};

        // Node (i,j) is at the physical position (i*delta, j*delta) on both meshes
        for (MFIter mfi(*z_phys_nd); mfi.isValid(); ++mfi) {
            const Box gbx = mfi.fabbox();
            auto z_arr = z_phys_nd->array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                z_arr(i,j,k) = slope_x * static_cast<Real>(i) * delta
                             + slope_y * static_cast<Real>(j) * delta
                             + static_cast<Real>(k) * dz;
            });
        }
        Gpu::streamSynchronize();

        if (m > zero) {
            BoxList bl2d = ba.boxList();
            for (auto& b : bl2d) { b.setRange(2,0); }
            const BoxArray ba2d(std::move(bl2d));
            mapfac.resize(MapFacType::num);
            mapfac[MapFacType::m_x] = std::make_unique<MultiFab>(ba2d, dm, 1, IntVect(3,3,0));
            mapfac[MapFacType::u_x] = std::make_unique<MultiFab>(convert(ba2d,IntVect(1,0,0)), dm, 1, IntVect(3,3,0));
            mapfac[MapFacType::v_x] = std::make_unique<MultiFab>(convert(ba2d,IntVect(0,1,0)), dm, 1, IntVect(3,3,0));
            for (auto& mf : mapfac) { mf->setVal(m); }
        }
    }

    static Geometry make_geometry (const Real h_dx)
    {
        const Box domain(IntVect(0), IntVect(ncell-1));
        const RealBox real_box({AMREX_D_DECL(zero, zero, zero)},
                               {AMREX_D_DECL(ncell*h_dx, ncell*h_dx, ncell*dz)});
        const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
        return Geometry(domain, &real_box, 0, is_periodic.data());
    }

    void make (MOSTAverageProbe& averages)
    {
        averages.make_MOSTAverage_at_level(0, vars_old, theta, qv, qr, z_phys_nd,
                                           mapfac.empty() ? nullptr : &mapfac);
        Gpu::streamSynchronize();
    }
};

// Largest |ind(i,j) - (idx_dir + offset)| over the valid cells
int
max_index_error (const iMultiFab& ind, const int dir, const int offset)
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    for (MFIter mfi(ind); mfi.isValid(); ++mfi) {
        const auto arr = ind.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<int>
            {
                const int idx = (dir == 0) ? i : j;
                const int err = arr(i,j,k) - (idx + offset);
                return {(err < 0) ? -err : err};
            });
    }
    return get<0>(reduce_data.value(reduce_op));
}

// Largest |pos(i,j)/h_dx - (idx_dir + 1/2) - offset| over the valid cells
Real
max_position_error (const MultiFab& pos, const int dir, const Real h_dx, const Real offset)
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    for (MFIter mfi(pos); mfi.isValid(); ++mfi) {
        const auto arr = pos.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real>
            {
                const int idx = (dir == 0) ? i : j;
                return {std::abs(arr(i,j,k)/h_dx - (static_cast<Real>(idx) + myhalf) - offset)};
            });
    }
    return get<0>(reduce_data.value(reduce_op));
}

// Largest |a - b| over the valid cells; a and b share a layout
Real
max_difference (const MultiFab& a, const MultiFab& b)
{
    MultiFab diff(a.boxArray(), a.DistributionMap(), 1, 0);
    MultiFab::Copy(diff, a, 0, 0, 1, 0);
    MultiFab::Subtract(diff, b, 0, 0, 1, 0);
    return diff.norm0(0, 0);
}

int
max_difference (const iMultiFab& a, const iMultiFab& b)
{
    iMultiFab diff(a.boxArray(), a.DistributionMap(), 1, 0);
    iMultiFab::Copy(diff, a, 0, 0, 1, 0);
    iMultiFab::Subtract(diff, b, 0, 0, 1, 0);
    return std::max(diff.max(0), -diff.min(0));
}

// Number of reference points, over the columns the averages read, whose interpolation
// stencil leaves the data of a field that reads them
int
unheld_points (const MOSTAverageProbe& averages, const MapFacLevel& level)
{
    const auto plo = level.geom.ProbLoArray();
    const auto dxi = level.geom.InvCellSizeArray();
    const IntVect ng_read = averages.get_ng_fill(0);
    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<int> reduce_data(reduce_op);
    for (MFIter mfi(*level.theta); mfi.isValid(); ++mfi) {
        const auto bounds = averages.interp_stencil_bounds(0, mfi);
        Box bx = amrex::grow(amrex::convert(mfi.validbox(), IntVect(1,1,0)), ng_read);
        bx.makeSlab(2, 0);
        const auto x_arr = averages.x_positions().const_array(mfi);
        const auto y_arr = averages.y_positions().const_array(mfi);
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<int>
            {
                Real xp = x_arr(i,j,k);
                Real yp = y_arr(i,j,k);
                bounds.clamp_to_held(i, j, xp, yp, plo, dxi);
                return {(xp != x_arr(i,j,k) || yp != y_arr(i,j,k)) ? 1 : 0};
            });
    }
    return get<0>(reduce_data.value(reduce_op));
}

// The unit normal of the terrain scaled by z_r, in physical lengths
std::array<Real, 3>
normal_displacement (const Real z_r = zref)
{
    const Real imag = one / std::sqrt(one + slope_x*slope_x + slope_y*slope_y);
    return {{-slope_x * imag * z_r, -slope_y * imag * z_r, imag * z_r}};
}

const Orientation zlo(Direction::z, Orientation::low);

} // namespace

// Motivation: without interpolation, the normal vector picks the column delta/(dx/m) cells
// away, where delta is the physical displacement along the physical normal
TEST(MOSTAverageMapFactors, NormalIndicesMatchFinerMesh)
{
    const char* prefix = "unit_most_mapfac_norm_indices";
    ScopedMapFacMOSTParams params(prefix);
    params.add("most.use_normal_vector", 1);
    params.add("most.zref", zref);

    MapFacLevel coarse(mfac*delta, mfac);
    MapFacLevel fine(delta, Real(-1.0));

    MOSTAverageProbe avg_coarse(zlo, Vector<Geometry>{coarse.geom}, true, prefix,
                                MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                                Vector<Vector<Real>>{});
    MOSTAverageProbe avg_fine(zlo, Vector<Geometry>{fine.geom}, true, prefix,
                              MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                              Vector<Vector<Real>>{});
    coarse.make(avg_coarse);
    fine.make(avg_fine);

    const auto d = normal_displacement();
    const int di = static_cast<int>(std::round(d[0] / delta));
    const int dj = static_cast<int>(std::round(d[1] / delta));
    ASSERT_NE(di, 0);
    ASSERT_NE(dj, 0);

    EXPECT_EQ(max_index_error(avg_fine.i_indices(), 0, di), 0);
    EXPECT_EQ(max_index_error(avg_fine.j_indices(), 1, dj), 0);
    EXPECT_EQ(max_index_error(avg_coarse.i_indices(), 0, di), 0);
    EXPECT_EQ(max_index_error(avg_coarse.j_indices(), 1, dj), 0);

    EXPECT_EQ(max_difference(avg_coarse.k_indices(), avg_fine.k_indices()), 0);
    EXPECT_LE(max_difference(*avg_coarse.get_zref(0), *avg_fine.get_zref(0)), tolerance(zref));
}

// Motivation: with interpolation, the reference point is x0 + m delta_x in computational
// coordinates, delta_x/(dx/m) cells from the column on both meshes
TEST(MOSTAverageMapFactors, NormalPositionsMatchFinerMesh)
{
    const char* prefix = "unit_most_mapfac_norm_positions";
    ScopedMapFacMOSTParams params(prefix);
    params.add("most.use_normal_vector", 1);
    params.add("most.use_interpolation", 1);
    params.add("most.zref", zref);

    MapFacLevel coarse(mfac*delta, mfac);
    MapFacLevel fine(delta, Real(-1.0));

    MOSTAverageProbe avg_coarse(zlo, Vector<Geometry>{coarse.geom}, true, prefix,
                                MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                                Vector<Vector<Real>>{});
    MOSTAverageProbe avg_fine(zlo, Vector<Geometry>{fine.geom}, true, prefix,
                              MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                              Vector<Vector<Real>>{});
    coarse.make(avg_coarse);
    fine.make(avg_fine);

    const auto d = normal_displacement();
    const Real di = d[0] / delta;
    const Real dj = d[1] / delta;
    const Real tol = tolerance(Real(ncell));

    EXPECT_LE(max_position_error(avg_fine.x_positions(), 0, delta, di), tol);
    EXPECT_LE(max_position_error(avg_fine.y_positions(), 1, delta, dj), tol);
    EXPECT_LE(max_position_error(avg_coarse.x_positions(), 0, mfac*delta, di), tol);
    EXPECT_LE(max_position_error(avg_coarse.y_positions(), 1, mfac*delta, dj), tol);

    const Real ztol = tolerance(ncell*dz);
    EXPECT_LE(max_difference(avg_coarse.z_positions(), avg_fine.z_positions()), ztol);
    EXPECT_LE(max_difference(*avg_coarse.get_zref(0), *avg_fine.get_zref(0)), ztol);
    EXPECT_NEAR(avg_coarse.get_zref(0)->min(0), d[2], ztol);
    EXPECT_NEAR(avg_coarse.get_zref(0)->max(0), d[2], ztol);
}

// Motivation: the tangent vectors use the physical slopes m_x h_xi and m_y h_eta
TEST(MOSTAverageMapFactors, RotatedVelocitiesMatchFinerMesh)
{
    const char* prefix = "unit_most_mapfac_rotation";
    ScopedMapFacMOSTParams params(prefix);
    params.add("most.terrain_rotate", 1);

    MapFacLevel coarse(mfac*delta, mfac);
    MapFacLevel fine(delta, Real(-1.0));

    MOSTAverageProbe avg_coarse(zlo, Vector<Geometry>{coarse.geom}, true, prefix,
                                MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                                Vector<Vector<Real>>{});
    MOSTAverageProbe avg_fine(zlo, Vector<Geometry>{fine.geom}, true, prefix,
                              MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                              Vector<Vector<Real>>{});
    coarse.make(avg_coarse);
    fine.make(avg_fine);
    avg_coarse.set_rotated_fields(0);
    avg_fine.set_rotated_fields(0);
    Gpu::streamSynchronize();

    const Real u_rot = (test_u + slope_x*test_w) / std::sqrt(one + slope_x*slope_x);
    const Real v_rot = (test_v + slope_y*test_w) / std::sqrt(one + slope_y*slope_y);

    for (auto* avg : {&avg_coarse, &avg_fine}) {
        EXPECT_NEAR(avg->rotated(0).min(0), u_rot, tolerance(u_rot));
        EXPECT_NEAR(avg->rotated(0).max(0), u_rot, tolerance(u_rot));
        EXPECT_NEAR(avg->rotated(1).min(0), v_rot, tolerance(v_rot));
        EXPECT_NEAR(avg->rotated(1).max(0), v_rot, tolerance(v_rot));
    }
    EXPECT_LE(max_difference(avg_coarse.rotated(0), avg_fine.rotated(0)), tolerance(u_rot));
    EXPECT_LE(max_difference(avg_coarse.rotated(1), avg_fine.rotated(1)), tolerance(v_rot));
}

// Motivation: the subgrid velocity scale takes the physical cell size dx/m, and all-ones
// map factors change nothing
TEST(MOSTAverageMapFactors, SubgridVelocityScaleUsesPhysicalCellSize)
{
    const char* prefix = "unit_most_mapfac_subgrid_vel";
    ScopedMapFacMOSTParams params(prefix);
    params.add("most.include_subgrid_vel", 1);

    // A physical spacing of 10 km, above the 5 km where the scale starts
    const Real big = Real(10000.0);
    MapFacLevel coarse(mfac*big, mfac);
    MapFacLevel fine(big, Real(-1.0));
    MapFacLevel fine_ones(big, one);

    MOSTAverageProbe avg_coarse(zlo, Vector<Geometry>{coarse.geom}, false, prefix,
                                MeshType::ConstantDz, TerrainType::None,
                                Vector<Vector<Real>>{});
    MOSTAverageProbe avg_fine(zlo, Vector<Geometry>{fine.geom}, false, prefix,
                              MeshType::ConstantDz, TerrainType::None,
                              Vector<Vector<Real>>{});
    MOSTAverageProbe avg_fine_ones(zlo, Vector<Geometry>{fine_ones.geom}, false, prefix,
                                   MeshType::ConstantDz, TerrainType::None,
                                   Vector<Vector<Real>>{});
    coarse.make(avg_coarse);
    fine.make(avg_fine);
    fine_ones.make(avg_fine_ones);

    const Real vsg = Real(0.32) * std::pow(big/Real(5000.) - one, Real(0.33));
    EXPECT_NEAR(avg_fine.subgrid_velocity_scale(), vsg, tolerance(vsg));
    EXPECT_NEAR(avg_coarse.subgrid_velocity_scale(), vsg, tolerance(vsg));
    EXPECT_NEAR(avg_fine_ones.subgrid_velocity_scale(), vsg, tolerance(vsg));
}

// Motivation: the reference points were required to lie in the grown tile, so on a CPU build
// a point displaced more than half a cell across an internal tile edge aborted, although the
// interpolation reads the box's own data there. With ERF's fewest ghost cells (2 for the
// velocities, 3 for the cell-centered fields) and 4-cell tiles, a displacement of 0.65 cells
// must pass and give the same reference points as one tile per box.
TEST(MOSTAverageInterpStencil, TileEdgesDoNotLimitTheReferencePoints)
{
    const char* prefix = "unit_most_interp_stencil_tiles";
    const Real z_r = Real(15.0);
    ScopedMapFacMOSTParams params(prefix);
    params.add("most.use_normal_vector", 1);
    params.add("most.use_interpolation", 1);
    params.add("most.zref", z_r);

    ScopedMFIterTileSize tile_size(IntVect(4));
    MapFacLevel level(delta, Real(-1.0), 2);
    MOSTAverageProbe averages(zlo, Vector<Geometry>{level.geom}, true, prefix,
                              MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                              Vector<Vector<Real>>{});
    level.make(averages);

    const auto d = normal_displacement(z_r);
    ASSERT_GT(std::abs(d[0] / delta), myhalf);
    const Real tol = tolerance(Real(ncell));
    EXPECT_LE(max_position_error(averages.x_positions(), 0, delta, d[0] / delta), tol);
    EXPECT_LE(max_position_error(averages.y_positions(), 1, delta, d[1] / delta), tol);
    EXPECT_EQ(unheld_points(averages, level), 0);
}

// Motivation: the averages read the reference points of ghost columns up to one cell from
// the edge of the fields' data, so a point displaced more than a cell outward read past the
// ghost cells, whatever their number. Here the x-velocity average reads the column one ghost
// cell out, displaced 1.75 cells: its stencil needs the third ghost cell of the velocity,
// which holds two. Such points must be moved back into the data, and the valid cells keep
// their exact points.
TEST(MOSTAverageInterpStencil, PointsPastTheGhostCellsAreClamped)
{
    for (const int ng_vel : {2, 3}) {
        const std::string prefix = "unit_most_interp_stencil_clamp_" + std::to_string(ng_vel);
        ScopedMapFacMOSTParams params(prefix.c_str());
        params.add("most.use_normal_vector", 1);
        params.add("most.use_interpolation", 1);
        params.add("most.zref", zref);

        MapFacLevel level(delta, Real(-1.0), ng_vel);
        MOSTAverageProbe averages(zlo, Vector<Geometry>{level.geom}, true, prefix,
                                  MeshType::ConstantDz, TerrainType::StaticFittedMesh,
                                  Vector<Vector<Real>>{});
        level.make(averages);

        EXPECT_EQ(unheld_points(averages, level), 0) << "ng_vel = " << ng_vel;

        const auto d = normal_displacement();
        const Real tol = tolerance(Real(ncell));
        EXPECT_LE(max_position_error(averages.x_positions(), 0, delta, d[0] / delta), tol)
            << "ng_vel = " << ng_vel;
        EXPECT_LE(max_position_error(averages.y_positions(), 1, delta, d[1] / delta), tol)
            << "ng_vel = " << ng_vel;
    }
}
