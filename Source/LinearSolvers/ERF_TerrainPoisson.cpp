#include "ERF_TerrainPoisson.H"

using namespace amrex;


TerrainPoisson::TerrainPoisson (amrex::Geometry const& geom, amrex::BoxArray const& ba,
                                amrex::DistributionMapping const& dm,
                                amrex::MultiFab const* z_phys_nd)
    : m_geom(geom),
      m_grids(ba),
      m_dmap(dm),
      m_zphys(z_phys_nd)
{
}

void TerrainPoisson::apply(amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
{
    AMREX_ASSERT(rhs.nGrowVect().allGT(0));

    auto domlo = lbound(m_geom.Domain());
    auto domhi = ubound(m_geom.Domain());

    MultiFab& xx = const_cast<MultiFab&>(rhs);

    auto const& dxinv = m_geom.InvCellSizeArray();

    auto const& y = lhs.arrays();
    auto const& x =  xx.arrays();
    auto const& zpa = m_zphys->const_arrays();

    // Impose periodic and internal boundary conditions
    xx.FillBoundary(m_geom.periodicity());

    // Impose top and bottom Neumann bcs
    amrex::ParallelFor(xx, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        if (k == domlo.z) {
            x[b](i,j,k-1) = x[b](i,j,k);
        } else if (k == domhi.z) {
            x[b](i,j,k+1) = x[b](i,j,k);
        }
    });

    auto const& xc = xx.const_arrays();
    amrex::ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        terrpoisson_adotx(i,j,k,y[b], xc[b], zpa[b], dxinv[0], dxinv[1], dxinv[2]);
    });
}

void TerrainPoisson::getFluxes(amrex::MultiFab const& phi,
                               amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fluxes)
{
    auto const& dxinv = m_geom.InvCellSizeArray();

    auto domlo = lbound(m_geom.Domain());
    auto domhi = ubound(m_geom.Domain());

    auto const& x   = phi.const_arrays();
    auto const& zpa = m_zphys->const_arrays();

    auto const& fx = fluxes[0].arrays();
    amrex::ParallelFor(fluxes[0], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        terrpoisson_flux_x(i,j,k,fx[b],x[b],zpa[b],dxinv[0],dxinv[2]);
    });

    auto const& fy = fluxes[1].arrays();
    amrex::ParallelFor(fluxes[1], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        terrpoisson_flux_y(i,j,k,fy[b],x[b],zpa[b],dxinv[1],dxinv[2]);
    });

    auto const& fz = fluxes[2].arrays();
    amrex::ParallelFor(fluxes[2], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        terrpoisson_flux_z(i,j,k,fz[b],x[b],zpa[b],dxinv[0],dxinv[1],dxinv[2]);

        if (k == domlo.z) {
            fz[b](i,j,k) = 0.0;
        } else if (k == domhi.z+1) {
            fz[b](i,j,k) = 0.0;
        }
    });
}

void TerrainPoisson::assign(amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
{
    MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
}

void TerrainPoisson::scale(amrex::MultiFab& lhs, amrex::Real fac)
{
    lhs.mult(fac);
}

Real TerrainPoisson::dotProduct(amrex::MultiFab const& v1, amrex::MultiFab const& v2)
{
    return MultiFab::Dot(v1, 0, v2, 0, 1, 0);
}

void TerrainPoisson::increment(amrex::MultiFab& lhs, amrex::MultiFab const& rhs, Real a)
{
    MultiFab::Saxpy(lhs, a, rhs, 0, 0, 1, 0);
}

void TerrainPoisson::linComb(amrex::MultiFab& lhs, Real a, amrex::MultiFab const& rhs_a,
                             Real b, amrex::MultiFab const& rhs_b)
{
    MultiFab::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, 1, 0);
}


MultiFab TerrainPoisson::makeVecRHS()
{
    return MultiFab(m_grids, m_dmap, 1, 0);
}

MultiFab TerrainPoisson::makeVecLHS()
{
    return MultiFab(m_grids, m_dmap, 1, 1);
}

Real TerrainPoisson::norm2(amrex::MultiFab const& v)
{
    return v.norm2();
}

void TerrainPoisson::precond(amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
{
    MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
}

void TerrainPoisson::setToZero(amrex::MultiFab& v)
{
    v.setVal(0);
}
