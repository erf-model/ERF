#include <algorithm>
#include <cmath>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Utils.H"

// Motivation: compressible levels coupled two-way, with map factors other than one, averaged
// the face momenta down as plain face averages.  The flux through a face is (rho u) ax / m_y
// (x), (rho v) ay / m_x (y) or (rho w) / (m_x m_y) (z), so with map factors and face areas that
// vary along a coarse face the plain average is not the coarse momentum whose flux is the mean
// of the fine fluxes over it (erf-model/ERF#5359).  average_down_faces_mapfac averages the
// fluxes: each covered coarse face must carry the mean of the fine fluxes over it, the faces the
// fine level does not cover keep their value, and a uniform state with one map factor stays.

using namespace amrex;

namespace {

constexpr int  ncx = 8, ncy = 8, ncz = 4;
constexpr Real untouched = Real(-77.0);
constexpr Real tol_rel = (sizeof(Real) == 8) ? Real(1.e-13) : Real(1.e-5);

// Smooth fields of the physical position (x, y, z in coarse cell units), so that coarse and
// fine faces see the same functions at their own locations
Real mapfac_x (Real x, Real y) { return Real(1.1) + Real(0.15)*std::sin(Real(0.7)*x) * std::cos(Real(0.5)*y); }
Real mapfac_y (Real x, Real y) { return Real(0.9) + Real(0.12)*std::cos(Real(0.6)*x + Real(0.4)*y); }
Real area     (Real x, Real y, Real z) { return Real(1.0) + Real(0.2)*std::sin(Real(0.9)*x - Real(0.3)*y + Real(0.5)*z); }
Real momentum (Real x, Real y, Real z) { return Real(2.0) + std::sin(Real(1.3)*x) * std::cos(Real(0.8)*y) + Real(0.3)*z; }

// The position of index (i,j,k) of a MultiFab of type 'nodal' on a mesh refined by r
// (r = 1 coarse), in coarse cell units
RealVect position (int i, int j, int k, const IntVect& nodal, const IntVect& r)
{
    return RealVect((Real(i) + (nodal[0] ? Real(0) : Real(0.5))) / Real(r[0]),
                    (Real(j) + (nodal[1] ? Real(0) : Real(0.5))) / Real(r[1]),
                    (Real(k) + (nodal[2] ? Real(0) : Real(0.5))) / Real(r[2]));
}

template <typename F>
void fill_host (MultiFab& mf, const IntVect& r, F&& f)
{
    const IntVect nodal = mf.ixType().toIntVect();
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf[mfi].array();
        const Box& b = mfi.fabbox();
        amrex::Loop(b, [&] (int i, int j, int k) {
            const RealVect p = position(i, j, k, nodal, r);
            a(i,j,k) = f(p[0], p[1], p[2]);
        });
    }
}

// Two-dimensional (no extent in z) map-factor MultiFab with the given face type, as ERF keeps them
MultiFab make_mf2d (const BoxArray& ba, const DistributionMapping& dm, const IntVect& nodal2d)
{
    BoxList bl = ba.boxList();
    for (auto& b : bl) { b.setRange(2, 0); }
    BoxArray ba2d(std::move(bl));
    return MultiFab(convert(ba2d, nodal2d), dm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
}

struct Levels
{
    IntVect rr;
    BoxArray cba, fba;
    DistributionMapping cdm, fdm;
    Box patch;   // the coarse cells the fine level covers

    explicit Levels (const IntVect& r)
        : rr(r),
          cba(Box(IntVect(0), IntVect(ncx-1, ncy-1, ncz-1))),
          patch(IntVect(2, 3, 0), IntVect(5, 6, 1))
    {
        cba.maxSize(4);
        cdm.define(cba);
        // Several fine boxes, each coarsenable by rr, as ERF's blocking factor makes them
        fba = BoxArray(refine(patch, rr));
        fba.maxSize(IntVect(2*rr[0], 2*rr[1], 2*rr[2]));
        fdm.define(fba);
    }
};

// Run average_down_faces_mapfac on the faces of direction dir and check every coarse face
void check_direction (const Levels& L, int dir, bool uniform)
{
    const IntVect nodal = IntVect::TheDimensionVector(dir);
    const IntVect nodal2d(nodal[0], nodal[1], 0);
    const IntVect one(1);
    auto info = MFInfo().SetArena(The_Pinned_Arena());

    MultiFab fine(convert(L.fba, nodal), L.fdm, 1, 0, info);
    MultiFab crse(convert(L.cba, nodal), L.cdm, 1, 0, info);
    MultiFab mf_f  = make_mf2d(L.fba, L.fdm, nodal2d), mf_c  = make_mf2d(L.cba, L.cdm, nodal2d);
    MultiFab mf2_f = make_mf2d(L.fba, L.fdm, nodal2d), mf2_c = make_mf2d(L.cba, L.cdm, nodal2d);
    MultiFab ar_f(convert(L.fba, nodal), L.fdm, 1, 0, info);
    MultiFab ar_c(convert(L.cba, nodal), L.cdm, 1, 0, info);

    // x-faces divide by m_y, y-faces by m_x, z-faces by m_x m_y (as the product of two)
    auto first  = [&] (Real x, Real y, Real) { return uniform ? Real(1.3) : (dir == 0 ? mapfac_y(x,y) : mapfac_x(x,y)); };
    auto second = [&] (Real x, Real y, Real) { return uniform ? Real(1.3) : mapfac_y(x,y); };
    fill_host(mf_f, L.rr, first);   fill_host(mf_c, one, first);
    fill_host(mf2_f, L.rr, second); fill_host(mf2_c, one, second);
    // The z-faces have unit area on every mesh ERF builds
    auto ar = [&] (Real x, Real y, Real z) { return (uniform || dir == 2) ? Real(1.0) : area(x,y,z); };
    fill_host(ar_f, L.rr, ar); fill_host(ar_c, one, ar);
    fill_host(fine, L.rr, [&] (Real x, Real y, Real z) { return uniform ? Real(2.5) : momentum(x,y,z); });
    crse.setVal(untouched);

    const bool two  = (dir == 2);
    const bool area_on = (dir != 2);
    average_down_faces_mapfac(fine, crse, L.rr, Periodicity::NonPeriodic(), mf_f, mf_c,
                              two ? &mf2_f : nullptr, two ? &mf2_c : nullptr,
                              area_on ? &ar_f : nullptr, area_on ? &ar_c : nullptr);

    // The fine fluxes and the coarse layout, gathered on the host
    const Box fdom = convert(refine(L.patch, L.rr), nodal);
    FArrayBox ff(fdom, 1, The_Pinned_Arena());
    FArrayBox fw(fdom, 1, The_Pinned_Arena());   // the fine weight a / m
    fine.copyTo(ff);
    {
        MultiFab w(convert(L.fba, nodal), L.fdm, 1, 0, info);
        for (MFIter mfi(w); mfi.isValid(); ++mfi) {
            auto wa = w[mfi].array(); auto aa = ar_f[mfi].const_array();
            auto m1 = mf_f[mfi].const_array(); auto m2 = mf2_f[mfi].const_array();
            amrex::Loop(mfi.validbox(), [&] (int i, int j, int k) {
                wa(i,j,k) = aa(i,j,k) / (m1(i,j,0) * (two ? m2(i,j,0) : Real(1.0)));
            });
        }
        w.copyTo(fw);
    }
    const Box cdom = convert(Box(IntVect(0), IntVect(ncx-1, ncy-1, ncz-1)), nodal);
    FArrayBox cf(cdom, 1, The_Pinned_Arena());
    crse.copyTo(cf);
    const Box cpatch = convert(L.patch, nodal);

    auto fa = ff.const_array(); auto wa = fw.const_array(); auto ca = cf.const_array();
    Real max_plain_err = 0;
    amrex::Loop(cdom, [&] (int i, int j, int k) {
        if (!cpatch.contains(IntVect(i,j,k))) {
            EXPECT_EQ(ca(i,j,k), untouched) << "uncovered face (" << i << "," << j << "," << k << ")";
            return;
        }
        // The fine faces on coarse face (i,j,k): the same plane in dir, rr cells across the others
        Real flux_sum = 0, plain_sum = 0;
        int n = 0;
        IntVect lo(i*L.rr[0], j*L.rr[1], k*L.rr[2]);
        IntVect hi = lo + L.rr - 1; hi[dir] = lo[dir];
        amrex::Loop(Box(lo, hi), [&] (int ii, int jj, int kk) {
            flux_sum  += fa(ii,jj,kk) * wa(ii,jj,kk);
            plain_sum += fa(ii,jj,kk);
            ++n;
        });
        const RealVect p = position(i, j, k, nodal, one);
        const Real m_c = first(p[0], p[1], p[2]) * (two ? second(p[0], p[1], p[2]) : Real(1.0));
        const Real a_c = ar(p[0], p[1], p[2]);
        const Real expected = flux_sum / Real(n) * m_c / a_c;
        EXPECT_NEAR(ca(i,j,k), expected, tol_rel * std::abs(expected))
            << "dir " << dir << ", coarse face (" << i << "," << j << "," << k << ")";
        max_plain_err = std::max(max_plain_err, std::abs(plain_sum / Real(n) - expected) / std::abs(expected));
    });
    if (uniform) {
        // One map factor and unit areas: the mean of a uniform momentum is that momentum
        amrex::Loop(cpatch, [&] (int i, int j, int k) {
            EXPECT_NEAR(ca(i,j,k), Real(2.5), tol_rel * Real(2.5));
        });
    } else {
        // The plain face average, what the coarse faces got before, is off by much more than
        // the tolerance, so the check above tells the two apart
        EXPECT_GT(max_plain_err, Real(1.e3) * tol_rel) << "dir " << dir;
    }
}

} // namespace

TEST(AverageDownMapFac, CoarseFaceCarriesTheMeanFineMassFlux)
{
    for (const IntVect& rr : {IntVect(2,2,2), IntVect(2,2,1), IntVect(3,3,1)}) {
        Levels L(rr);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            SCOPED_TRACE("rr " + std::to_string(rr[0]) + std::to_string(rr[1]) + std::to_string(rr[2]));
            check_direction(L, dir, false);
        }
    }
}

TEST(AverageDownMapFac, UniformStateWithOneMapFactorStaysUniform)
{
    Levels L(IntVect(2,2,2));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        check_direction(L, dir, true);
    }
}
