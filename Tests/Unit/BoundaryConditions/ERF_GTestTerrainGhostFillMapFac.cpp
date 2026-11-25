#include <gtest/gtest.h>

#include <memory>

#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_PhysBCFunct.H>

/**
 * @file ERF_GTestTerrainGhostFillMapFac.cpp
 * @brief Zero-normal-gradient (foextrap) ghost fill below sloping terrain with map factors
 *
 * impose_vertical_{cons,xvel,yvel}_bcs correct the zero-gradient fill below a
 * terrain-fitted mesh so that the derivative along the physical surface normal
 * is zero. With map factors the physical slopes are m_x h_xi and m_y h_eta and
 * the physical horizontal derivatives are m_x d/dxi and m_y d/deta.
 *
 * The terrain is a plane, h = a xi + b eta, and the field is linear in the
 * computational horizontal coordinates and the physical height,
 *     phi = c1 xi + c2 eta + gz z,
 * so its physical gradient is (m_x c1, m_y c2, gz). With gz = m_x^2 a c1 + m_y^2 b c2
 * that gradient is orthogonal to the physical normal (-m_x a, -m_y b, 1), and
 * every one-sided and centered difference the fill takes is exact, so the ghost
 * cells must hold phi itself. Distinct m_x and m_y catch a missing or swapped
 * factor; with the slopes in computational coordinates the error is O(0.1).
 */

namespace {

using namespace amrex;

constexpr Real a  = Real( 0.3);   // terrain slope d h / d xi
constexpr Real b  = Real(-0.2);   // terrain slope d h / d eta
constexpr Real c1 = Real( 1.5);
constexpr Real c2 = Real( 0.7);

enum class Loc { cell, xface, yface };

Real gz_for (Real mx, Real my) { return mx*mx*a*c1 + my*my*b*c2; }

// Offsets of the variable's location from its (i,j) index, in units of the cell size
Real off_x (Loc loc) { return (loc == Loc::xface) ? Real(0) : Real(0.5); }
Real off_y (Loc loc) { return (loc == Loc::yface) ? Real(0) : Real(0.5); }

// phi at (i,j,k) for a variable at offsets (ox,oy) on a unit mesh
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_phi (int i, int j, int k, Real ox, Real oy, Real gz)
{
    Real xi  = Real(i) + ox;
    Real eta = Real(j) + oy;
    Real z   = a*xi + b*eta + Real(k) + Real(0.5);
    return c1*xi + c2*eta + gz*z;
}

/**
 * Fill the ghost cells below the bottom with the functor for loc and return the largest
 * difference from phi there.
 *
 * @param[in] loc         where the variable lives (cons, u, v)
 * @param[in] with_mapfac pass map factors m_x = 2, m_y = 1.25; otherwise pass none,
 *                        which the fill must treat as unit map factors
 */
Real ghost_fill_error (Loc loc, bool with_mapfac)
{
    const Box domain(IntVect(0), IntVect(7,7,3));
    const RealBox rb(Real(0), Real(0), Real(0), Real(8), Real(8), Real(4));
    const Geometry geom(domain, rb, CoordSys::cartesian, {0, 0, 0});
    const auto dxInv = geom.InvCellSizeArray();

    Vector<BCRec> bcs(BCVars::NumTypes);
    for (auto& bc : bcs) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bc.setLo(d, ERFBCType::foextrap);
            bc.setHi(d, ERFBCType::foextrap);
        }
    }
    Gpu::DeviceVector<BCRec> device_bcs(bcs.size());
    Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), device_bcs.begin());
    Array<Array<Real,6>, AMREX_SPACEDIM+NBCVAR_max> extdir{}, neumann{};

    // Nodal heights of the planar terrain on a unit mesh
    const BoxArray nba(convert(domain, IntVect(1)));
    const DistributionMapping dm(nba);
    auto z = std::make_unique<MultiFab>(nba, dm, 1, 2);
    for (MFIter mfi(*z); mfi.isValid(); ++mfi) {
        const auto zarr = z->array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            zarr(i,j,k) = a*Real(i) + b*Real(j) + Real(k);
        });
    }
    const Array4<const Real> z_nd = z->const_array(0);

    // The functor's own map-factor pointers are not used by impose_vertical_*_bcs
    Vector<std::unique_ptr<MultiFab>> mapfac(MapFacType::num);

    const Real mx = with_mapfac ? Real(2)    : Real(1);
    const Real my = with_mapfac ? Real(1.25) : Real(1);
    const Real gz = gz_for(mx, my);
    const Real ox = off_x(loc), oy = off_y(loc);

    // Fill box: the variable's valid box plus one ghost layer below the bottom
    Box vbx = domain;
    if (loc == Loc::xface) { vbx.surroundingNodes(0); }
    if (loc == Loc::yface) { vbx.surroundingNodes(1); }
    Box bx = vbx; bx.setSmall(2, -1);

    // 2D map factors at the variable's location
    Box mbx = grow(vbx, 2); mbx.setSmall(2, 0); mbx.setBig(2, 0);
    FArrayBox mfx_fab(mbx, 1), mfy_fab(mbx, 1);
    mfx_fab.setVal<RunOn::Device>(mx);
    mfy_fab.setVal<RunOn::Device>(my);
    const Array4<const Real> mf_x = with_mapfac ? mfx_fab.const_array() : Array4<const Real>{};
    const Array4<const Real> mf_y = with_mapfac ? mfy_fab.const_array() : Array4<const Real>{};

    const int ncomp = (loc == Loc::cell) ? NVAR_max : 1;
    const int comp  = (loc == Loc::cell) ? RhoTheta_comp : 0;
    FArrayBox fab(grow(vbx, 1), ncomp);
    {
        const auto arr = fab.array();
        ParallelFor(fab.box(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            // Ghost cells below the bottom start with a value the fill must overwrite
            arr(i,j,k,n) = (k < 0) ? Real(-999) : exact_phi(i,j,k,ox,oy,gz);
        });
    }

    if (loc == Loc::cell) {
        ERFPhysBCFunct_cons bc(0, geom, bcs, device_bcs, extdir, neumann,
                               TerrainType::StaticFittedMesh, mapfac, z, false, nullptr, {});
        bc.impose_vertical_cons_bcs(fab.array(), bx, domain, z_nd, mf_x, mf_y, dxInv,
                                    comp, 1, 0.0);
    } else if (loc == Loc::xface) {
        ERFPhysBCFunct_u bc(0, geom, bcs, device_bcs, extdir, neumann,
                            TerrainType::StaticFittedMesh, mapfac, z, false, nullptr);
        bc.impose_vertical_xvel_bcs(fab.array(), bx, domain, z_nd, mf_x, mf_y, dxInv,
                                    BCVars::xvel_bc, 0.0);
    } else {
        ERFPhysBCFunct_v bc(0, geom, bcs, device_bcs, extdir, neumann,
                            TerrainType::StaticFittedMesh, mapfac, z, false, nullptr);
        bc.impose_vertical_yvel_bcs(fab.array(), bx, domain, z_nd, mf_x, mf_y, dxInv,
                                    BCVars::yvel_bc, 0.0);
    }
    Gpu::streamSynchronize();

    FArrayBox host(fab.box(), ncomp, The_Pinned_Arena());
    Gpu::dtoh_memcpy(host.dataPtr(), fab.dataPtr(), fab.nBytes());
    const auto h = host.const_array();

    Real err = Real(0);
    const auto lo = lbound(vbx), hi = ubound(vbx);
    for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            err = std::max(err, std::abs(h(i,j,-1,comp) - exact_phi(i,j,-1,ox,oy,gz)));
        }
    }
    return err;
}

} // namespace

constexpr amrex::Real tol = amrex::Real(1.e-12);

TEST(TerrainGhostFillMapFac, CellCenteredFillHasZeroPhysicalNormalGradient)
{
    EXPECT_LT(ghost_fill_error(Loc::cell, true), tol);
}

TEST(TerrainGhostFillMapFac, XFaceFillHasZeroPhysicalNormalGradient)
{
    EXPECT_LT(ghost_fill_error(Loc::xface, true), tol);
}

TEST(TerrainGhostFillMapFac, YFaceFillHasZeroPhysicalNormalGradient)
{
    EXPECT_LT(ghost_fill_error(Loc::yface, true), tol);
}

// Without map-factor arrays the fill must use unit map factors, not read through null
TEST(TerrainGhostFillMapFac, MissingMapFactorsMeanUnitMapFactors)
{
    EXPECT_LT(ghost_fill_error(Loc::cell,  false), tol);
    EXPECT_LT(ghost_fill_error(Loc::xface, false), tol);
    EXPECT_LT(ghost_fill_error(Loc::yface, false), tol);
}
