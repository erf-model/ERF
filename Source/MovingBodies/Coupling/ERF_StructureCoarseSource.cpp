/**
 * @file ERF_StructureCoarseSource.cpp
 * @brief The restriction of a structure model's momentum source onto the coarser levels
 */
#include "ERF_StructureCoarseSource.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>

using namespace amrex;

namespace {

/**
 * Restrict a face field from a level onto the faces normal to dir of the next coarser level:
 * coarse face I takes the fine faces r I + m, |m| < r, along dir with the weights
 * (r - |m|) / r^2, averaged over the fine faces across it. Applied to the source times J_f,
 * this is the force-space restriction of restrict_face_sources divided by the coarse
 * control volume (dx dy dz of the coarse level is r_dir r_t1 r_t2 times the fine one).
 * On a non-periodic domain face, where the boundary condition sets the momentum and a
 * source does nothing, a fine face inside the domain gives its whole force to the first
 * coarse face inside instead of sharing it with the domain face (weight r / r^2 there);
 * the fine domain face's own force stays on the coarse domain face.
 * fine must have r - 1 ghost faces along dir, filled where another box of its level holds
 * them and zero outside the level.
 */
void
restrict_face_field (const MultiFab& fine, MultiFab& crse, int dir, const IntVect& rr, const Geometry& cgeom)
{
    AMREX_ASSERT(fine.boxArray().coarsenable(rr));
    BoxArray cba = fine.boxArray();
    cba.coarsen(rr);
    MultiFab crse_on_fine(cba, fine.DistributionMap(), 1, 0);

    const int rn = rr[dir];
    const int t1 = (dir + 1) % AMREX_SPACEDIM;
    const int t2 = (dir + 2) % AMREX_SPACEDIM;
    const int r1 = rr[t1];
    const int r2 = rr[t2];
    const Real inv_w = Real(1.0) / Real(rn * rn * r1 * r2);

    // The coarse domain faces normal to dir, where the boundary condition sets the momentum
    const bool walls = !cgeom.isPeriodic(dir);
    const int clo = cgeom.Domain().smallEnd(dir);
    const int chi = cgeom.Domain().bigEnd(dir) + 1;
    AMREX_ALWAYS_ASSERT(!walls || chi - clo >= 2);

    for (MFIter mfi(crse_on_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto f = fine.const_array(mfi);
        const auto c = crse_on_fine.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const IntVect I(AMREX_D_DECL(i, j, k));
            const IntVect base = I * rr;
            const int Id = I[dir];
            Real s = Real(0.0);
            for (int m = 1 - rn; m < rn; ++m) {
                Real w = Real(rn - amrex::Math::abs(m));
                if (walls && m != 0) {
                    if ((Id == clo && m > 0) || (Id == chi && m < 0)) {
                        continue;                       // given to the face inside
                    }
                    if ((Id == clo + 1 && m < 0) || (Id == chi - 1 && m > 0)) {
                        w = Real(rn);                   // the whole force, the domain face's share too
                    }
                }
                for (int a = 0; a < r1; ++a) {
                    for (int b = 0; b < r2; ++b) {
                        IntVect iv = base;
                        iv[dir] += m;
                        iv[t1]  += a;
                        iv[t2]  += b;
                        s += w * f(iv);
                    }
                }
            }
            c(i, j, k) = s * inv_w;
        });
    }

    crse.setVal(Real(0.0));
    crse.ParallelCopy(crse_on_fine, 0, 0, 1, IntVect(0), IntVect(0), cgeom.periodicity());
}

} // namespace

namespace erf_structure {

void
scale_by_face_jacobian (MultiFab& mf, int dir, const MultiFab& detJ, const Geometry& geom, bool divide)
{
    AMREX_ALWAYS_ASSERT(detJ.nGrow() >= 1);
    // A face on a non-periodic domain face takes the Jacobian of its one cell inside the
    // domain: detJ_cc is not set below the ground (make_J)
    const bool periodic = geom.isPeriodic(dir);
    const int dlo = geom.Domain().smallEnd(dir);
    const int dhi = geom.Domain().bigEnd(dir) + 1;

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto a = mf.array(mfi);
        const auto J = detJ.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const IntVect hi(AMREX_D_DECL(i, j, k));   // the cell on the high side of the face
            IntVect lo = hi;                           // and the one on its low side
            lo[dir] -= 1;
            Real Jf;
            if (!periodic && hi[dir] == dlo) {
                Jf = J(hi);
            } else if (!periodic && hi[dir] == dhi) {
                Jf = J(lo);
            } else {
                Jf = Real(0.5) * (J(lo) + J(hi));
            }
            a(i, j, k) = divide ? a(i, j, k) / Jf : a(i, j, k) * Jf;
        });
    }
}

void
restrict_face_sources (std::array<MultiFab,AMREX_SPACEDIM>& src, int L, int lev,
                       const Vector<BoxArray>& grids,
                       const Vector<DistributionMapping>& dmap,
                       const Vector<Geometry>& geom,
                       const Vector<IntVect>& ref_ratio,
                       const Vector<const MultiFab*>& detJ)
{
    AMREX_ALWAYS_ASSERT(lev >= 0 && lev < L);

    // The source times J_f on level L: the force per nominal cell volume dx dy dz
    if (detJ[L]) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            scale_by_face_jacobian(src[d], d, *detJ[L], geom[L], false);
        }
    }

    // Down one level at a time to level lev, each with the ghost faces the next restriction reads
    for (int l = L - 1; l >= lev; --l) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            IntVect ng_fine(0);
            ng_fine[d] = ref_ratio[l][d] - 1;
            if (src[d].nGrowVect() != ng_fine) {
                MultiFab grown(src[d].boxArray(), src[d].DistributionMap(), 1, ng_fine);
                grown.setVal(Real(0.0));
                MultiFab::Copy(grown, src[d], 0, 0, 1, 0);
                src[d] = std::move(grown);
            }
            src[d].FillBoundary(geom[l+1].periodicity());
            MultiFab crse(convert(grids[l], IntVect::TheDimensionVector(d)), dmap[l], 1, 0);
            restrict_face_field(src[d], crse, d, ref_ratio[l], geom[l]);
            src[d] = std::move(crse);
        }
    }

    // Back from force per nominal volume to the source on level lev
    if (detJ[lev]) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            scale_by_face_jacobian(src[d], d, *detJ[lev], geom[lev], true);
        }
    }
}

std::array<Real,AMREX_SPACEDIM>
face_source_force (const std::array<MultiFab,AMREX_SPACEDIM>& src, const Geometry& geom, const MultiFab* detJ)
{
    std::array<Real,AMREX_SPACEDIM> force{};
    const Real dV = AMREX_D_TERM(geom.CellSize(0), * geom.CellSize(1), * geom.CellSize(2));
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        MultiFab f(src[d].boxArray(), src[d].DistributionMap(), 1, 0);
        MultiFab::Copy(f, src[d], 0, 0, 1, 0);
        if (detJ) { scale_by_face_jacobian(f, d, *detJ, geom, false); }
        force[d] = f.sum_unique(0, false, geom.periodicity()) * dV;
    }
    return force;
}

} // namespace erf_structure
