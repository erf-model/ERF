/**
 * \file ERF_ImmersedPoisson.cpp
 */
#ifdef ERF_USE_FFT

#include "ERF_ImmersedPoisson.H"
#include "ERF_NumericalConstants.H"
#include "ERF_SolverUtils.H"

using namespace amrex;

ImmersedPoisson::ImmersedPoisson (Geometry const& geom, Geometry const& lev_geom,
                                  BoxArray const& ba, DistributionMapping const& dm,
                                  Array<std::string,2*AMREX_SPACEDIM> const& domain_bcs_type,
                                  Array<MultiFab const*,AMREX_SPACEDIM> const& sigma,
                                  bool use_real_bcs)
    : m_geom(geom), m_grids(ba), m_dmap(dm), m_sigma(sigma)
{
    m_bc_fft = get_fft_bc(lev_geom, domain_bcs_type, ba.minimalBox(), use_real_bcs);
    m_fft_precond = std::make_unique<FFT::Poisson<MultiFab>>(geom, m_bc_fft);
}

void ImmersedPoisson::apply_bcs (MultiFab& phi)
{
    phi.FillBoundary(m_geom.periodicity());

    const Box& domain = m_geom.Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        if (m_geom.isPeriodic(dir)) { continue; }
        for (int side = 0; side < 2; ++side)
        {
            const FFT::Boundary bc = (side == 0) ? m_bc_fft[dir].first : m_bc_fft[dir].second;
            if (bc == FFT::Boundary::periodic) { continue; }
            // even: zero normal gradient on the boundary face; odd: phi = 0 on the boundary face
            const Real sgn = (bc == FFT::Boundary::odd) ? -one : one;
            const int edge = (side == 0) ? domain.smallEnd(dir) : domain.bigEnd(dir);
            const IntVect shift = (side == 0) ? -IntVect::TheDimensionVector(dir)
                                              :  IntVect::TheDimensionVector(dir);
            for (MFIter mfi(phi); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                if (bx.smallEnd(dir) > edge || bx.bigEnd(dir) < edge) { continue; }
                const Array4<Real>& p = phi.array(mfi);
                ParallelFor(makeSlab(bx, dir, edge), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    p(i+shift[0], j+shift[1], k+shift[2]) = sgn * p(i, j, k);
                });
            }
        }
    }
}

void ImmersedPoisson::apply (MultiFab& lhs, MultiFab const& rhs)
{
    AMREX_ASSERT(rhs.nGrowVect().allGT(0));

    MultiFab& xx = const_cast<MultiFab&>(rhs);
    apply_bcs(xx);

    auto const& dxinv = m_geom.InvCellSizeArray();
    const Real idx2 = dxinv[0]*dxinv[0];
    const Real idy2 = dxinv[1]*dxinv[1];
    const Real idz2 = dxinv[2]*dxinv[2];

    auto const& y  = lhs.arrays();
    auto const& x  = xx.const_arrays();
    auto const& sx = m_sigma[0]->const_arrays();
    auto const& sy = m_sigma[1]->const_arrays();
    auto const& sz = m_sigma[2]->const_arrays();

    ParallelFor(lhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        auto const& p = x[b];
        const Real pc = p(i,j,k);
        y[b](i,j,k) = idx2 * ( sx[b](i+1,j,k) * (p(i+1,j,k) - pc) - sx[b](i,j,k) * (pc - p(i-1,j,k)) )
                    + idy2 * ( sy[b](i,j+1,k) * (p(i,j+1,k) - pc) - sy[b](i,j,k) * (pc - p(i,j-1,k)) )
                    + idz2 * ( sz[b](i,j,k+1) * (p(i,j,k+1) - pc) - sz[b](i,j,k) * (pc - p(i,j,k-1)) );
    });
    Gpu::streamSynchronize();
}

void ImmersedPoisson::getFluxes (MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM>& fluxes)
{
    apply_bcs(phi);

    auto const& dxinv = m_geom.InvCellSizeArray();
    const Real dx_inv = dxinv[0];
    const Real dy_inv = dxinv[1];
    const Real dz_inv = dxinv[2];

    auto const& x  = phi.const_arrays();
    auto const& fx = fluxes[0].arrays();
    auto const& fy = fluxes[1].arrays();
    auto const& fz = fluxes[2].arrays();
    ParallelFor(fluxes[0], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        fx[b](i,j,k) = -(x[b](i,j,k) - x[b](i-1,j,k)) * dx_inv;
    });
    ParallelFor(fluxes[1], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        fy[b](i,j,k) = -(x[b](i,j,k) - x[b](i,j-1,k)) * dy_inv;
    });
    ParallelFor(fluxes[2], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        fz[b](i,j,k) = -(x[b](i,j,k) - x[b](i,j,k-1)) * dz_inv;
    });
    Gpu::streamSynchronize();
}

void ImmersedPoisson::assign (MultiFab& lhs, MultiFab const& rhs)
{
    MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
}

void ImmersedPoisson::scale (MultiFab& lhs, Real fac)
{
    lhs.mult(fac);
}

Real ImmersedPoisson::dotProduct (MultiFab const& v1, MultiFab const& v2)
{
    return MultiFab::Dot(v1, 0, v2, 0, 1, 0);
}

void ImmersedPoisson::increment (MultiFab& lhs, MultiFab const& rhs, Real a)
{
    MultiFab::Saxpy(lhs, a, rhs, 0, 0, 1, 0);
}

void ImmersedPoisson::linComb (MultiFab& lhs, Real a, MultiFab const& rhs_a,
                               Real b, MultiFab const& rhs_b)
{
    MultiFab::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, 1, 0);
}

MultiFab ImmersedPoisson::makeVecRHS ()
{
    return MultiFab(m_grids, m_dmap, 1, 0);
}

MultiFab ImmersedPoisson::makeVecLHS ()
{
    return MultiFab(m_grids, m_dmap, 1, 1);
}

Real ImmersedPoisson::norm2 (MultiFab const& v)
{
    return v.norm2();
}

void ImmersedPoisson::precond (MultiFab& lhs, MultiFab const& rhs)
{
    if (m_use_precond) {
        MultiFab& rhs_tmp = const_cast<MultiFab&>(rhs);
        lhs.setVal(zero);
        m_fft_precond->solve(lhs, rhs_tmp);
    } else {
        MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
    }
}

void ImmersedPoisson::setToZero (MultiFab& v)
{
    v.setVal(zero);
}
#endif
