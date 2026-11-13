// Contract of erf_structure::restrict_face_sources (the structure model's source on the coarser
// levels, ERF::make_structure_coarse_sources): the force the source puts on the air, the sum
// of the source times the faces' control volumes J_f dx dy dz, is the same on every coarser
// level as on the model's level, over a hill (where the Jacobian varies from column to column
// and with height), one and two levels down, with an anisotropic refinement ratio, with the
// fine level split into boxes, across a periodic boundary and at the ground; the force on the
// faces inside the domain (not on the ground and top, where the boundary condition sets w) is
// kept too. Without the Jacobian the same restriction misses the force on this hill by 3.8e-5
// (x), 1.6e-3 (y) and 3.4e-3 (z), so the check is sensitive to it. On a uniform-dz mesh (no
// Jacobian) the force is kept too. scale_by_face_jacobian takes the one cell inside the
// domain on a non-periodic domain face, so it never reads detJ below the ground.

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <gtest/gtest.h>

#include "ERF_StructureCoarseSource.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>

namespace {

using amrex::Real;

constexpr Real Lx = 1600.0, Ly = 1200.0, H = 800.0;

// A basic terrain-following hill: detJ of a column is (H - h) / H, at the column's centre
Real hill (Real x, Real y)
{
    constexpr Real pi = 3.14159265358979323846;
    return Real(300.0) * (Real(0.5) + Real(0.5) * std::cos(Real(2.0) * pi * x / Lx))
                       * (Real(0.5) + Real(0.5) * std::sin(Real(2.0) * pi * y / Ly));
}

// A field with no structure the restriction could cancel: a smooth part and a hashed part,
// a function of the face's position only, so faces shared by two boxes agree
Real source_value (int dir, int i, int j, int k)
{
    const unsigned h = static_cast<unsigned>(i * 73856093) ^ static_cast<unsigned>(j * 19349663)
                     ^ static_cast<unsigned>(k * 83492791) ^ static_cast<unsigned>(dir * 2654435761u);
    return Real(1.0) + Real(0.3) * std::sin(Real(0.7) * i + Real(1.3) * j + Real(0.4) * k)
         + Real(0.5) * Real(h % 1000u) / Real(1000.0);
}

struct Hierarchy {
    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<amrex::IntVect> rr;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> detj;

    // Level 0: 16 x 12 x 8 cells, periodic in x and y. Level 1 (ratio 2): an interior patch and a
    // patch across the periodic x boundary. Level 2 (ratio rr2): inside the interior patch.
    // Boxes of at most max_grid cells, so each level is split.
    Hierarchy (int nlev, const amrex::IntVect& rr2, bool terrain, int max_grid = 4)
    {
        const amrex::Box dom0(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 11, 7));
        const amrex::RealBox rb({AMREX_D_DECL(Real(0.0), Real(0.0), Real(0.0))}, {AMREX_D_DECL(Lx, Ly, H)});
        const std::array<int,3> periodic{{1, 1, 0}};
        rr = {amrex::IntVect(2, 2, 2), rr2};

        amrex::Box dom = dom0;
        for (int l = 0; l < nlev; ++l) {
            if (l > 0) { dom.refine(rr[l-1]); }
            geom.emplace_back(dom, &rb, 0, periodic.data());
        }

        amrex::BoxList bl0(dom0);
        grids.emplace_back(bl0);
        if (nlev > 1) {
            amrex::BoxList bl1;
            bl1.push_back(amrex::Box(amrex::IntVect(8, 4, 2), amrex::IntVect(23, 19, 11)));   // interior
            bl1.push_back(amrex::Box(amrex::IntVect(0, 6, 0), amrex::IntVect(3, 13, 7)));     // at x = 0
            bl1.push_back(amrex::Box(amrex::IntVect(28, 6, 0), amrex::IntVect(31, 13, 7)));   // at x = Lx
            grids.emplace_back(bl1);
        }
        if (nlev > 2) {
            amrex::Box b2(amrex::IntVect(10, 6, 4), amrex::IntVect(21, 17, 9));   // in level-1 cells
            b2.refine(rr2);
            grids.emplace_back(amrex::BoxList(b2));
        }
        for (int l = 0; l < nlev; ++l) {
            grids[l].maxSize(max_grid * (l + 1));
            dmap.emplace_back(grids[l]);
        }

        detj.resize(nlev);
        if (terrain) {
            for (int l = 0; l < nlev; ++l) {
                detj[l] = std::make_unique<amrex::MultiFab>(grids[l], dmap[l], 1, 1);
                // below the ground make_J leaves detJ unset: NaN there, which the restriction must not read
                detj[l]->setVal(std::numeric_limits<Real>::quiet_NaN());
                for (amrex::MFIter mfi(*detj[l], false); mfi.isValid(); ++mfi) {
                    auto ja = detj[l]->array(mfi);
                    amrex::Box gbx = amrex::grow(mfi.validbox(), 1);
                    gbx.setSmall(2, std::max(gbx.smallEnd(2), 0));   // nothing below the ground
                    amrex::LoopOnCpu(gbx, [&](int i, int j, int k) { ja(i,j,k) = cell_jacobian(l, i, j, k); });
                }
            }
        }
    }

    amrex::Vector<const amrex::MultiFab*> detj_ptrs (bool use) const
    {
        amrex::Vector<const amrex::MultiFab*> p(detj.size(), nullptr);
        for (int l = 0; l < detj.size(); ++l) { if (use) { p[l] = detj[l].get(); } }
        return p;
    }

    // The source on the faces of level l (no ghost faces)
    std::array<amrex::MultiFab,3> source (int l) const
    {
        std::array<amrex::MultiFab,3> s;
        for (int d = 0; d < 3; ++d) {
            s[d].define(amrex::convert(grids[l], amrex::IntVect::TheDimensionVector(d)), dmap[l], 1, 0);
            for (amrex::MFIter mfi(s[d], false); mfi.isValid(); ++mfi) {
                auto a = s[d].array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i,j,k) = source_value(d, i, j, k); });
            }
        }
        return s;
    }

    // The Jacobian of cell (i, j, k) of level l: (H - h) / H of its column times a factor that
    // varies with height, so the z faces of a column see it vary too
    Real cell_jacobian (int l, int i, int j, int k) const
    {
        constexpr Real pi = 3.14159265358979323846;
        const Real zc = (k + Real(0.5)) * geom[l].CellSize(2);
        return (H - hill((i + Real(0.5)) * geom[l].CellSize(0), (j + Real(0.5)) * geom[l].CellSize(1))) / H
             * (Real(1.0) + Real(0.3) * std::sin(pi * zc / H));
    }

    // J_f of face (i, j, k) normal to dir on level l, from the formula above, not from the code
    // under test: the mean of the two cells beside the face (x and y are periodic), the one cell
    // inside on the ground and top faces
    Real face_jacobian (int l, int dir, int i, int j, int k) const
    {
        const int nz = geom[l].Domain().bigEnd(2) + 1;
        if (dir == 2 && k == 0)  { return cell_jacobian(l, i, j, 0); }
        if (dir == 2 && k == nz) { return cell_jacobian(l, i, j, nz - 1); }
        return Real(0.5) * (cell_jacobian(l, i - (dir == 0), j - (dir == 1), k - (dir == 2)) + cell_jacobian(l, i, j, k));
    }

    // The force along dir: the source times J_f dx dy dz over the faces of level l, each face
    // once; with inside_only, without the ground and top faces
    Real force (const amrex::MultiFab& src, int dir, int l, bool terrain, bool inside_only = false) const
    {
        const int nz = geom[l].Domain().bigEnd(2) + 1;
        const auto mask = src.OwnerMask(geom[l].periodicity());
        Real sum = 0;
        for (amrex::MFIter mfi(src, false); mfi.isValid(); ++mfi) {
            const auto a = src.const_array(mfi);
            const auto m = mask->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (!m(i,j,k)) { return; }
                if (inside_only && dir == 2 && (k == 0 || k == nz)) { return; }
                sum += a(i,j,k) * (terrain ? face_jacobian(l, dir, i, j, k) : Real(1.0));
            });
        }
        amrex::ParallelAllReduce::Sum(sum, amrex::ParallelContext::CommunicatorSub());
        const Real dV = geom[l].CellSize(0) * geom[l].CellSize(1) * geom[l].CellSize(2);
        return sum * dV;
    }
};

constexpr Real tol = (sizeof(Real) == 8) ? Real(1.0e-13) : Real(1.0e-5);

// Restrict level L's source to level lev and compare the force on each component
std::array<Real,3> relative_force_error (const Hierarchy& h, int L, int lev, bool use_jacobian,
                                         bool inside_only = false)
{
    std::array<amrex::MultiFab,3> src = h.source(L);
    std::array<Real,3> fine{};
    for (int d = 0; d < 3; ++d) { fine[d] = h.force(src[d], d, L, h.detj[L] != nullptr, inside_only); }
    erf_structure::restrict_face_sources(src, L, lev, h.grids, h.dmap, h.geom, h.rr, h.detj_ptrs(use_jacobian));
    std::array<Real,3> err{};
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(src[d].boxArray(), amrex::convert(h.grids[lev], amrex::IntVect::TheDimensionVector(d)));
        const Real crse = h.force(src[d], d, lev, h.detj[lev] != nullptr, inside_only);
        err[d] = std::abs(crse - fine[d]) / std::abs(fine[d]);
    }
    return err;
}

} // namespace

TEST(StructureCoarseSource, KeepsTheForceOverAHill)
{
    for (int nlev : {2, 3}) {
        const Hierarchy h(nlev, amrex::IntVect(2, 2, 1), true);
        const int L = nlev - 1;
        for (int lev = L - 1; lev >= 0; --lev) {
            for (bool inside : {false, true}) {
                SCOPED_TRACE("levels " + std::to_string(nlev) + ", down to level " + std::to_string(lev)
                             + (inside ? ", faces inside the domain" : ", all faces"));
                const auto err = relative_force_error(h, L, lev, true, inside);
                for (int d = 0; d < 3; ++d) {
                    SCOPED_TRACE("component " + std::to_string(d));
                    EXPECT_LT(err[d], tol);
                }
            }
        }
    }
}

TEST(StructureCoarseSource, IsotropicRatiosOverAHill)
{
    const Hierarchy h(3, amrex::IntVect(2, 2, 2), true);
    const auto err = relative_force_error(h, 2, 0, true);
    for (int d = 0; d < 3; ++d) {
        SCOPED_TRACE("component " + std::to_string(d));
        EXPECT_LT(err[d], tol);
    }
}

TEST(StructureCoarseSource, WithoutTheJacobianTheHillChangesTheForce)
{
    // The check above is sensitive to the Jacobian: the index-space restriction alone (as
    // before the Jacobian was applied) misses the force on this hill by 3.8e-5 (x), 1.6e-3 (y)
    // and 3.4e-3 (z), above the tolerance of the checks above in single precision too
    const Hierarchy h(2, amrex::IntVect(2, 2, 1), true);
    const auto err = relative_force_error(h, 1, 0, false);
    for (int d = 0; d < 3; ++d) {
        SCOPED_TRACE("component " + std::to_string(d));
        EXPECT_GT(err[d], Real(2.0e-5));
    }
}

TEST(StructureCoarseSource, KeepsTheForceOnAUniformMesh)
{
    const Hierarchy h(3, amrex::IntVect(2, 2, 1), false);
    for (int lev : {1, 0}) {
        SCOPED_TRACE("down to level " + std::to_string(lev));
        const auto err = relative_force_error(h, 2, lev, false);
        for (int d = 0; d < 3; ++d) {
            SCOPED_TRACE("component " + std::to_string(d));
            EXPECT_LT(err[d], tol);
        }
    }
}

TEST(StructureCoarseSource, FaceJacobianOnTheDomainFaces)
{
    // On the ground and top faces J_f is the Jacobian of the one cell inside the domain (detJ
    // below the ground is NaN here); elsewhere the mean of the two cells beside the face, also
    // across the periodic seams.
    // Multiplying and dividing gives the source back.
    const Hierarchy h(1, amrex::IntVect(2, 2, 1), true);
    std::array<amrex::MultiFab,3> src = h.source(0);
    amrex::MultiFab z(src[2].boxArray(), src[2].DistributionMap(), 1, 0);
    z.setVal(Real(1.0));
    erf_structure::scale_by_face_jacobian(z, 2, *h.detj[0], h.geom[0], false);
    for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
        const auto a = z.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const Real want = h.face_jacobian(0, 2, i, j, k);
            EXPECT_NEAR(a(i,j,k), want, Real(10.0) * std::numeric_limits<Real>::epsilon()) << i << " " << j << " " << k;
        });
    }
    // x and y faces (periodic): the mean of the two cells beside the face, also across the seam
    for (int d = 0; d < 2; ++d) {
        amrex::MultiFab f(src[d].boxArray(), src[d].DistributionMap(), 1, 0);
        f.setVal(Real(1.0));
        erf_structure::scale_by_face_jacobian(f, d, *h.detj[0], h.geom[0], false);
        for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
            const auto a = f.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                const Real want = h.face_jacobian(0, d, i, j, k);
                EXPECT_NEAR(a(i,j,k), want, Real(10.0) * std::numeric_limits<Real>::epsilon())
                    << "component " << d << " at " << i << " " << j << " " << k;
            });
        }
    }
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab s(src[d].boxArray(), src[d].DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(s, src[d], 0, 0, 1, 0);
        erf_structure::scale_by_face_jacobian(s, d, *h.detj[0], h.geom[0], false);
        erf_structure::scale_by_face_jacobian(s, d, *h.detj[0], h.geom[0], true);
        amrex::MultiFab::Subtract(s, src[d], 0, 0, 1, 0);
        EXPECT_LT(s.norm0(), Real(1.0e3) * std::numeric_limits<Real>::epsilon()) << "component " << d;
    }
}

TEST(StructureCoarseSource, ForceDiagnostic)
{
    // face_source_force (printed with erf.v > 0) against the force computed here from the hill
    for (bool terrain : {true, false}) {
        SCOPED_TRACE(terrain ? "hill" : "uniform");
        const Hierarchy h(2, amrex::IntVect(2, 2, 1), terrain);
        for (int l = 0; l < 2; ++l) {
            const std::array<amrex::MultiFab,3> src = h.source(l);
            const auto f = erf_structure::face_source_force(src, h.geom[l], h.detj[l].get());
            for (int d = 0; d < 3; ++d) {
                const Real want = h.force(src[d], d, l, terrain);
                EXPECT_LT(std::abs(f[d] - want), tol * std::abs(want)) << "level " << l << " component " << d;
            }
        }
    }
}
