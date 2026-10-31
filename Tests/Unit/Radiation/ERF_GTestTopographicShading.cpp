#include <gtest/gtest.h>
#include <ERF_TopographicShading.H>
#include <ERF_TopographicShadingKernels.H>
#include <ERF_OrbCosZenith.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;
using T = TopographicShading;

namespace {
// Split the ridge and its shadow into separate boxes, distributed round-robin
// when built into erf_parallel_tests. The target's local terrain alone is flat.
struct TerrainCase {
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    MultiFab z, map;
    TopographicShading topo;
    TopographicShadingChoice choice;

    TerrainCase (bool periodic = false, Real dx = 1, Real dy = 1) {
        const Box dom(IntVect(0,0,0), IntVect(15,3,3));
        RealBox rb({0,0,0}, {16*dx,4*dy,20});
        int per[] = {int(periodic),0,0};
        geom.define(dom, &rb, 0, per);
        ba = BoxArray(makeSlab(dom,2,0));
        ba.maxSize(IntVect(4,4,1));
        Vector<int> owners(ba.size());
        for (int b = 0; b < int(ba.size()); ++b) { owners[b] = b % ParallelDescriptor::NProcs(); }
        dm = DistributionMapping(owners);
        BoxList bl = ba.boxList();
        for (auto& b : bl) { b.setRange(2,0,4); }
        z.define(convert(BoxArray(std::move(bl)),IntVect(1)),dm,1,0);
        map.define(ba,dm,1,0);
        map.setVal(1);
        choice.topo_shading = true;
        choice.slope_rad = true;
        choice.fixed_mu = std::sqrt(Real(0.5));
        choice.fixed_azimuth = 90;
        choice.shadlen = 20;
        topo.define_level(0,ba,dm,geom);
    }
    void terrain (int kind) {
        for (MFIter mfi(z); mfi.isValid(); ++mfi) {
            auto a = z.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i,int j,int k) {
                // kind 1: narrow ridge of height 10 at node 10 (mass points 9/10 are 5).
                // kind 2: planar slope dh/d(grid x) = 1/4.
                a(i,j,k) = kind == 1 ? (i == 10 ? Real(10) : Real(0)) :
                             (kind == 2 ? Real(0.25)*i : Real(0));
            });
        }
    }
    void initialize () {
        topo.initialize_geometry(0,choice,z,map,map,nullptr,nullptr);
        const Real flux[] = {choice.fixed_mu,400,100,200,30,70,350};
        for (int n = 0; n < T::NRad; ++n) { topo.raw(0)->setVal(flux[n],n,1); }
        topo.set_sun(0,1.5,0);
    }
};

template <typename F>
void check_cells (const MultiFab& mf, F check)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        FArrayBox host(mfi.validbox(),mf.nComp(),The_Pinned_Arena());
        Gpu::dtoh_memcpy(host.dataPtr(), mf[mfi].dataPtr(), host.size()*sizeof(Real));
        const auto a = host.const_array();
        const Box b = mfi.validbox();
        for (int j=b.smallEnd(1); j<=b.bigEnd(1); ++j) {
            for (int i=b.smallEnd(0); i<=b.bigEnd(0); ++i) { check(a,i,j); }
        }
    }
}
constexpr Real tol = sizeof(Real) == 8 ? Real(1.e-10) : Real(1.e-3);
}

TEST(TopographicShading, SolarVectorMatchesOrbitalConvention)
{
    for (double day : {80.3,80.5,80.7,250.9}) {
        for (Real lat : {Real(-35),Real(0),Real(40)}) {
            const auto s = erf_topography::sun_direction(day,0.2,lat,-100);
            const double expected = orbital_cos_zenith_instant(day,lat*PI/180.,-100.*PI/180.,0.2);
            EXPECT_NEAR(s.z,expected,tol);
            EXPECT_NEAR(s.x*s.x+s.y*s.y+s.z*s.z,1,tol);
        }
    }
    EXPECT_GT(erf_topography::sun_direction(80.3,0,0,0).x,0);
    EXPECT_LT(erf_topography::sun_direction(80.7,0,0,0).x,0);
}

TEST(TopographicShading, FlatGroundAndRepeatedApplicationPreserveForcing)
{
    TerrainCase c;
    c.terrain(0); c.initialize();
    MultiFab out(c.ba,c.dm,T::NRad,0);
    Vector<std::unique_ptr<MultiFab>> aliases;
    Vector<MultiFab*> fields;
    for (int n=0;n<T::NRad;++n) {
        aliases.push_back(std::make_unique<MultiFab>(out,make_alias,n,1));
        fields.push_back(aliases.back().get());
    }
    c.topo.apply(0,c.choice,fields);
    c.topo.apply(0,c.choice,fields);
    MultiFab::Subtract(out,*c.topo.raw(0),0,0,T::NRad,0);
    for (int n=0;n<T::NRad;++n) { EXPECT_NEAR(out.norm0(n),0,tol); }
    EXPECT_EQ(c.topo.diagnostics(0)->max(T::Shadow),0);
    EXPECT_EQ(c.topo.diagnostics(0)->min(T::DirectFactor),1);
}

TEST(TopographicShading, RidgeShadowsAcrossBoxesAndPreservesDiffuse)
{
    TerrainCase c;
    c.choice.slope_rad = false;
    c.terrain(1); c.initialize();
    c.topo.apply(0,c.choice,{});
    check_cells(*c.topo.diagnostics(0), [](auto a,int i,int j) {
        // At a 45-degree elevation, the 5 m ridge at mass point 9
        // blocks points 5..8, but does not block the exactly tangent point 4.
        const bool shadow = i>=5 && i<=8;
        EXPECT_EQ(a(i,j,0,T::Shadow),shadow ? 1 : 0) << "i=" << i;
        EXPECT_NEAR(a(i,j,0,T::TerrainSW),shadow ? 100 : 400,tol);
        EXPECT_EQ(a(i,j,0,T::UnshadedSW),400);
    });
}

TEST(TopographicShading, SlopeUsesPhysicalMapFactors)
{
    TerrainCase c(false,2,3);
    c.choice.topo_shading = false;
    c.terrain(2); c.map.setVal(2); c.initialize();
    MultiFab sw(c.ba,c.dm,1,0), dir(c.ba,c.dm,1,0), dif(c.ba,c.dm,1,0), lw(c.ba,c.dm,1,0), mu(c.ba,c.dm,1,0);
    const Real factor = Real(0.75)/std::sqrt(Real(1.0625));
    c.topo.apply(0,c.choice,{&mu,&sw,&dir,nullptr,&dif,nullptr,&lw});
    EXPECT_NEAR(sw.min(0),100+300*factor,tol);
    EXPECT_NEAR(sw.max(0),100+300*factor,tol);
    EXPECT_NEAR(dir.min(0),100*factor,tol);
    EXPECT_EQ(dif.min(0),30);
    EXPECT_EQ(lw.min(0),350);
    EXPECT_EQ(mu.min(0),c.choice.fixed_mu);
    c.topo.apply(0,c.choice,{&mu,&sw,&dir,nullptr,&dif,nullptr,&lw});
    EXPECT_NEAR(sw.min(0),100+300*factor,tol);
}

TEST(TopographicShading, SearchDistanceAndPeriodicBoundary)
{
    TerrainCase c(true);
    c.terrain(1); c.initialize();
    c.choice.fixed_azimuth = 270;
    c.choice.fixed_mu = Real(0.2);
    c.topo.apply(0,c.choice,{});
    // Sun to the west: point 0 sees the previous periodic copy of the ridge.
    check_cells(*c.topo.diagnostics(0), [](auto a,int i,int j) {
        if (i==0) { EXPECT_EQ(a(i,j,0,T::Shadow),1); }
    });
    c.choice.shadlen = Real(0.5);
    c.topo.apply(0,c.choice,{});
    EXPECT_EQ(c.topo.diagnostics(0)->max(T::Shadow),0);
}

TEST(TopographicShading, NightAndZenithHaveNoTerrainShadow)
{
    TerrainCase c;
    c.terrain(1); c.initialize();
    c.choice.fixed_mu = 1;
    c.topo.apply(0,c.choice,{});
    EXPECT_EQ(c.topo.diagnostics(0)->max(T::Shadow),0);
    c.choice.fixed_mu = -1;
    c.choice.latitude = 0; c.choice.longitude = 0;
    c.topo.initialize_geometry(0,c.choice,c.z,c.map,c.map,nullptr,nullptr);
    c.topo.set_sun(0,1.0,0);
    c.topo.apply(0,c.choice,{});
    EXPECT_EQ(c.topo.diagnostics(0)->max(T::Shadow),0);
    EXPECT_EQ(c.topo.diagnostics(0)->max(T::Incidence),0);
}

TEST(TopographicShading, FinePatchSeesAncestorRidgeAndUsesRawFlux)
{
    TerrainCase c;
    c.choice.slope_rad = false;
    c.terrain(1); c.initialize(); c.topo.apply(0,c.choice,{});
    Geometry fine = c.geom;
    fine.refine(IntVect(2,2,1));
    BoxArray ba(Box(IntVect(10,0,0),IntVect(15,7,0)));
    DistributionMapping dm(ba);
    MultiFab map(ba,dm,1,0); map.setVal(1);
    BoxArray nodes = ba;
    BoxList bl = nodes.boxList();
    for (auto& b:bl) { b.setRange(2,0,4); }
    MultiFab z(convert(BoxArray(std::move(bl)),IntVect(1)),dm,1,0); z.setVal(0);
    c.topo.define_level(1,ba,dm,fine);
    c.topo.initialize_geometry(1,c.choice,z,map,map,nullptr,nullptr);
    c.topo.interpolate_from_coarse(1);
    c.topo.apply(1,c.choice,{});
    EXPECT_EQ(c.topo.raw(1)->min(T::Shortwave),400);
    EXPECT_EQ(c.topo.diagnostics(1)->min(T::Shadow),1);
    EXPECT_EQ(c.topo.diagnostics(1)->max(T::TerrainSW),100);
    // Rebuilding a parent invalidates descendant terrain and parent raw forcing.
    c.topo.define_level(0,c.ba,c.dm,c.geom);
    EXPECT_FALSE(c.topo.geometry_ready(1));
    EXPECT_FALSE(c.topo.forcing_ready(0));
}

TEST(TopographicShading, GeographicSunRotatesOntoProjectedGrid)
{
    TerrainCase c;
    c.terrain(2); c.initialize();
    MultiFab lat(c.ba,c.dm,1,IntVect(1,1,0)), lon(c.ba,c.dm,1,IntVect(1,1,0));
    for (MFIter mfi(lat); mfi.isValid(); ++mfi) {
        auto la = lat.array(mfi), lo = lon.array(mfi);
        ParallelFor(lat[mfi].box(), [=] AMREX_GPU_DEVICE (int i,int j,int k) {
            la(i,j,k) = Real(40)+Real(0.01)*i;
            lo(i,j,k) = Real(-100)-Real(0.01)*j;
        });
    }
    c.topo.initialize_geometry(0,c.choice,c.z,c.map,c.map,&lat,&lon);
    c.topo.apply(0,c.choice,{});
    // grid x is north: east sunlight is parallel to this ridge/slope.
    EXPECT_NEAR(c.topo.diagnostics(0)->min(T::DirectFactor),1/std::sqrt(Real(1.0625)),tol);
    EXPECT_EQ(c.topo.diagnostics(0)->max(T::Shadow),0);
}

TEST(TopographicShading, SpectralForcingWritesSlmSurfaceGhostPlaneOnly)
{
    TerrainCase c;
    c.terrain(1); c.initialize();
    BoxList bl = c.ba.boxList();
    for (auto& b:bl) { b.setRange(2,-2,2); }
    MultiFab soil(BoxArray(std::move(bl)),c.dm,1,1);
    soil.setVal(-999);
    c.topo.apply(0,c.choice,{nullptr,nullptr,&soil});
    EXPECT_EQ(soil.min(0),-999); // valid soil layers untouched
    EXPECT_EQ(soil.max(0),-999);
    MultiFab surface(c.ba,c.dm,1,0);
    for (MFIter mfi(surface); mfi.isValid(); ++mfi) {
        auto d = surface.array(mfi); auto s = soil.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i,int j,int k) {
            d(i,j,k) = s(i,j,0);
        });
    }
    EXPECT_EQ(surface.min(0),0); // direct VIS extinguished in shadow
    EXPECT_GE(surface.max(0),100);
}

TEST(TopographicShading, ObliqueRayInterpolatesAndUsesBothCellWidths)
{
    const Box b(IntVect(0,0,0),IntVect(7,7,0));
    FArrayBox terrain(b,1,The_Pinned_Arena());
    terrain.setVal(0);
    auto a = terrain.array();
    // Sun travels one x cell and two y cells per sample: dx=2, dy=1,
    // geographic x/y components equal. At (2,4), terrain crosses the ray.
    a(2,4,0) = 10;
    const Real h = std::sqrt(Real(0.48));
    EXPECT_TRUE(erf_topography::shadowed(0,0,terrain.const_array(),0,7,0,7,
        false,false,2,1,{h,h,Real(0.2)},10));
    EXPECT_FALSE(erf_topography::shadowed(0,0,terrain.const_array(),0,7,0,7,
        false,false,1,1,{h,h,Real(0.2)},10));
}
