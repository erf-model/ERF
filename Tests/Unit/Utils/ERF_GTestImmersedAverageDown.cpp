#include <cmath>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_IndexDefines.H"
#include "ERF_Utils.H"

// Motivation: with immersed forcing every AMR level builds its own volume fraction and forces its
// own solid and partial cells, so near the surface the levels disagree about the geometry. The
// plain average down mixed solid values into coarse fluid cells (a half-solid coarse cell under a
// fine level held half the fine fluid value) and fluid values into coarse solid cells.
// if_average_down (erf.if_average_down = fluid_mass) weights the fine cells by their fluid mass and
// leaves coarse cells that are solid on their own level alone. The surface h = 1.25 sits inside a coarse cell (dz = 1) and,
// for rz = 2, inside a fine cell too. The density falls with height, as the anelastic base state
// does: a fluid-only weight without the density biased a uniform theta by -0.04 K in a
// two-level anelastic run, so the weight carries the density.

using namespace amrex;

namespace {

constexpr int  ncx = 8;
constexpr int  ncy = 2;
constexpr int  ncz = 6;
constexpr int  ncomp_state = 3;   // Rho_comp, then two conserved scalars (rho phi)
constexpr Real h_surf    = Real(1.25);
constexpr Real old_value = Real(-99.);
// tolerances: the sums carry a few ulps per fine cell
constexpr Real tol_abs = (sizeof(Real) == 8) ? Real(1.e-12) : Real(1.e-5);
constexpr Real tol_rel = (sizeof(Real) == 8) ? Real(1.e-10) : Real(1.e-4);

// Solid fraction of cell k of height dz under a flat surface at h (1 below the mesh)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real beta_of (int k, Real dz, Real h)
{
    return amrex::min(Real(1.), amrex::max(Real(0.), (h - k*dz) / dz));
}

// Density at the centre of cell k of height dz, falling with height
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rho_of (int k, Real dz)
{
    return Real(1.2) - Real(0.1) * (k + Real(0.5)) * dz;
}

// Fine cell volume weight, varying in x and z as detJ / (m_x m_y) does on a stretched or mapped
// grid (the code passes the real one; any constant factor cancels)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real vol_of (int i, int k)
{
    return Real(1.) + Real(0.05) * (k % 3) + Real(0.02) * (i % 2);
}

// Specific value of scalar n (0 or 1) in fine cell (i,k): n = 0 uniform (theta), n = 1 varies in
// x and z and is zero in solid cells, as the forcing leaves it
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real phi_of (int i, int k, int n, Real beta)
{
    if (n == 0) { return Real(300.); }
    return (beta >= Real(1.)) ? Real(0.) : Real(10.)*k + i;
}

struct Case
{
    IntVect ratio;
    int     scomp;           // Rho_comp: the density is averaged too (compressible); 1: not (anelastic)
    Real    h_c;             // surface height the coarse level sees (the levels may disagree)
    bool    uniform_vol;     // unit volume weight, or the varying one
};

void
fill_vol (MultiFab& vol, bool uniform)
{
    for (MFIter mfi(vol); mfi.isValid(); ++mfi) {
        auto const v = vol.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            v(i,j,k) = uniform ? Real(1.) : vol_of(i, k);
        });
    }
}

void
fill_beta (MultiFab& beta, Real dz, Real h)
{
    for (MFIter mfi(beta); mfi.isValid(); ++mfi) {
        auto const b = beta.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            b(i,j,k) = beta_of(k, dz, h);
        });
    }
}

// Fine state: rho, rho phi_0, rho phi_1
void
fill_fine_state (MultiFab& S, Real dz)
{
    for (MFIter mfi(S); mfi.isValid(); ++mfi) {
        auto const s = S.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real r = rho_of(k, dz);
            const Real b = beta_of(k, dz, h_surf);
            s(i,j,k,Rho_comp) = r;
            s(i,j,k,1) = r * phi_of(i, k, 0, b);
            s(i,j,k,2) = r * phi_of(i, k, 1, b);
        });
    }
}

// All of a coarse MultiFab on the host (every rank), on one box
struct HostData
{
    Box box;
    std::vector<Real> v;
    Real at (int i, int j, int k, int n) const
    {
        return v[static_cast<std::size_t>(box.index(IntVect(i,j,k)) + n*box.numPts())];
    }
};

HostData
to_host (const MultiFab& mf, const Box& domain)
{
    MultiFab one(BoxArray(domain), DistributionMapping(BoxArray(domain)), mf.nComp(), 0);
    one.ParallelCopy(mf, 0, 0, mf.nComp());
    HostData h{domain, std::vector<Real>(static_cast<std::size_t>(domain.numPts()*mf.nComp()), Real(-1.e30))};
    for (MFIter mfi(one); mfi.isValid(); ++mfi) {
        Gpu::copy(Gpu::deviceToHost, one[mfi].dataPtr(), one[mfi].dataPtr() + one[mfi].size(), h.v.begin());
    }
    Gpu::streamSynchronize();
    ParallelDescriptor::Bcast(h.v.data(), h.v.size(), one.DistributionMap()[0]);
    return h;
}

// A fine patch over coarse x cells 2-5 of 8, all of y and z, split into several boxes on its own
// distribution, so the sums cross ranks
void
check (const Case& c, Real h_f = h_surf)
{
    const Box domain_c(IntVect(0,0,0), IntVect(ncx-1,ncy-1,ncz-1));
    BoxArray ba_c(domain_c); ba_c.maxSize(IntVect(4,2,ncz));
    const Box patch(IntVect(2,0,0), IntVect(5,ncy-1,ncz-1));
    BoxArray ba_f(amrex::refine(patch, c.ratio)); ba_f.maxSize(IntVect(4,2,64));
    DistributionMapping dm_c(ba_c), dm_f(ba_f);
    const Real dzf = Real(1.) / c.ratio[2];

    MultiFab beta_c(ba_c, dm_c, 1, 1), beta_f(ba_f, dm_f, 1, 1);
    fill_beta(beta_c, Real(1.), c.h_c);
    fill_beta(beta_f, dzf, h_f);
    MultiFab S_f(ba_f, dm_f, ncomp_state, 0), S_c(ba_c, dm_c, ncomp_state, 0);
    MultiFab vol_f(ba_f, dm_f, 1, 0);
    fill_fine_state(S_f, dzf);
    fill_vol(vol_f, c.uniform_vol);
    S_c.setVal(old_value);
    // the coarse cell's own density (read when the density is not averaged)
    for (MFIter mfi(S_c); mfi.isValid(); ++mfi) {
        auto const s = S_c.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            s(i,j,k,Rho_comp) = rho_of(k, Real(1.));
        });
    }
    Gpu::streamSynchronize();

    const int ncomp = ncomp_state - c.scomp;
    if_average_down(S_f, S_c, beta_f, beta_c, vol_f, c.scomp, ncomp, c.ratio);
    const HostData h = to_host(S_c, domain_c);

    const int rx = c.ratio[0], ry = c.ratio[1], rz = c.ratio[2];
    for (int K = 0; K < ncz; ++K) {
    for (int J = 0; J < ncy; ++J) {
    for (int I = 0; I < ncx; ++I) {
        SCOPED_TRACE(testing::Message() << "ratio " << c.ratio << ", scomp " << c.scomp
                     << ", h_c " << c.h_c << ", uniform vol " << c.uniform_vol << ", cell " << IntVect(I,J,K));
        Real w = 0, wrho = 0, wphi[2] = {0, 0};
        for (int kk = K*rz; kk < (K+1)*rz; ++kk) {
        for (int jj = J*ry; jj < (J+1)*ry; ++jj) {
        for (int ii = I*rx; ii < (I+1)*rx; ++ii) {
            const Real bf = beta_of(kk, dzf, h_f);
            const Real r  = rho_of(kk, dzf);
            const Real vf = (Real(1.) - bf) * (c.uniform_vol ? Real(1.) : vol_of(ii, kk));
            w    += vf;
            wrho += vf * r;
            for (int n = 0; n < 2; ++n) { wphi[n] += vf * r * phi_of(ii, kk, n, bf); }
        }}}
        const bool covered = (I >= 2 && I <= 5);
        const bool apply   = covered && beta_of(K, Real(1.), c.h_c) < Real(1.) && w > 0;
        const Real rho_own = rho_of(K, Real(1.));
        // expected coarse density and conserved scalars
        Real rho_c = rho_own;
        if (apply && c.scomp == Rho_comp) { rho_c = wrho / w; }
        EXPECT_NEAR(h.at(I,J,K,Rho_comp), rho_c, tol_abs);
        for (int n = 0; n < 2; ++n) {
            Real expect = old_value;
            if (apply) { expect = rho_c * wphi[n] / wrho; }
            EXPECT_NEAR(h.at(I,J,K,1+n), expect, tol_rel * (Real(1.) + std::abs(expect)));
        }
        // the mass weighting keeps a uniform theta exactly uniform
        if (apply) {
            EXPECT_NEAR(h.at(I,J,K,1) / h.at(I,J,K,Rho_comp), Real(300.), tol_rel * Real(300.));
        }
    }}}
}

} // namespace

// Anelastic (the density is not averaged) and compressible (it is), isotropic ratio 2, an
// anisotropic 4 x 1 x 2, and 2 x 2 x 4 with the surface on a fine face inside the coarse cell
TEST(ImmersedAverageDown, CellsAreFluidMassWeighted)
{
    for (int scomp : {1, Rho_comp}) {
        check({IntVect(2,2,2), scomp, h_surf, true});
        check({IntVect(4,1,2), scomp, h_surf, true});
        check({IntVect(2,2,4), scomp, h_surf, true});
        check({IntVect(2,2,1), scomp, h_surf, true});   // no vertical refinement
    }
}

// The fine cell volume (detJ / (m_x m_y) on a stretched or mapped grid) enters the weight: with a
// volume varying inside the coarse cell a uniform theta must still come back uniform, and the
// scalar must be the volume-and-mass-weighted mean
TEST(ImmersedAverageDown, FineVolumeEntersTheWeight)
{
    for (int scomp : {1, Rho_comp}) {
        check({IntVect(2,2,2), scomp, h_surf, false});
        check({IntVect(4,1,2), scomp, h_surf, false});
    }
}

// A coarse fluid cell over fine cells that are all solid (the levels disagree the other way:
// the coarse level sees the surface at z = 1, its cell k = 1 is fluid; with h_f = 2 the fine
// cells 2 and 3 under it are solid) keeps its value: no fluid mass to average
TEST(ImmersedAverageDown, CoarseFluidOverSolidFineCellsKeepsItsValue)
{
    check({IntVect(2,2,2), 1, Real(1.), true}, Real(2.));
}

// The levels disagree about the geometry: the coarse level sees the surface at z = 2, so its cell
// k = 1 is solid while the fine cells over it still hold fluid. That coarse cell keeps its value.
TEST(ImmersedAverageDown, CoarseSolidKeepsItsValueWhenTheLevelsDisagree)
{
    check({IntVect(2,2,2), 1, Real(2.), true});
    check({IntVect(2,2,2), Rho_comp, Real(2.), true});
}
