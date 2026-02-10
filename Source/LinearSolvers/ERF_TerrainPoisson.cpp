#ifdef ERF_USE_FFT

#include "ERF_TerrainPoisson.H"
#include "ERF_SolverUtils.H"

using namespace amrex;

TerrainPoisson::TerrainPoisson (Geometry const& geom, BoxArray const& ba,
                                DistributionMapping const& dm,
                                Array<std::string,2*AMREX_SPACEDIM>& domain_bcs_type,
                                Gpu::DeviceVector<Real>& stretched_dz_lev_d,
                                const MultiFab& ax, const MultiFab& ay,
                                const MultiFab& az, const MultiFab& dJ,
                                MultiFab const* z_phys_nd,
                                bool use_real_bcs)
    : m_geom(geom),
      m_grids(ba),
      m_dmap(dm),
      m_domain_bcs_type(domain_bcs_type),
      m_stretched_dz_d(stretched_dz_lev_d),
      m_ax(ax),
      m_ay(ay),
      m_az(az),
      m_dJ(dJ),
      m_zphys(z_phys_nd)
{
    if (!m_2D_fft_precond) {
        Box bounding_box = ba.minimalBox();
        bc_fft = get_fft_bc(geom,domain_bcs_type,bounding_box,use_real_bcs);
        m_2D_fft_precond = std::make_unique<FFT::PoissonHybrid<MultiFab>>(geom,bc_fft);
    }
}

void TerrainPoisson::usePrecond (bool use_precond_in)
{
    m_use_precond = use_precond_in;
}

void TerrainPoisson::apply (MultiFab& lhs, MultiFab const& rhs)
{
    AMREX_ASSERT(rhs.nGrowVect().allGT(0));

    MultiFab& xx = const_cast<MultiFab&>(rhs);

    auto const& dxinv = m_geom.InvCellSizeArray();

    auto const& y = lhs.arrays();
    auto const& zpa = m_zphys->const_arrays();
    auto const& axa = m_ax.const_arrays();
    auto const& aya = m_ay.const_arrays();
    auto const& aza = m_az.const_arrays();
    auto const& dJa = m_dJ.const_arrays();

    apply_bcs(xx);

    auto const& xc = xx.const_arrays();
    ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        terrpoisson_adotx(i, j, k, y[b], xc[b], axa[b], aya[b], aza[b], dJa[b], zpa[b], dxinv[0], dxinv[1], dxinv[2]);
    });
}

void TerrainPoisson::apply_bcs (MultiFab& phi)
{
    auto domlo = lbound(m_geom.Domain());
    auto domhi = ubound(m_geom.Domain());

    phi.FillBoundary(m_geom.periodicity());

    if (!m_geom.isPeriodic(0)) {
        for (MFIter mfi(phi,true); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            const Array4<Real>& phi_arr = phi.array(mfi);
            if (bx.smallEnd(0) <= domlo.x) {
                if (bc_fft[0].first == FFT::Boundary::even) {
                    ParallelFor(makeSlab(bx,0,domlo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i-1,j,k) =  phi_arr(i,j,k);
                    });
                } else if (bc_fft[0].first == FFT::Boundary::odd) {
                    ParallelFor(makeSlab(bx,0,domlo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i-1,j,k) =  -phi_arr(i,j,k);
                    });
                }
            } // lo x
            if (bx.bigEnd(0) >= domhi.x) {
                if (bc_fft[0].second == FFT::Boundary::even) {
                    ParallelFor(makeSlab(bx,0,domhi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i+1,j,k) =  phi_arr(i,j,k);
                    });
                } else if (bc_fft[0].second == FFT::Boundary::odd) {
                    ParallelFor(makeSlab(bx,0,domhi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i+1,j,k) =  -phi_arr(i,j,k);
                    });
                }
            } // hi x
        } // mfi
    } // not periodic in x

    if (!m_geom.isPeriodic(1)) {
        for (MFIter mfi(phi,true); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            Box bx2(bx); bx2.grow(0,1);
            const Array4<Real>& phi_arr = phi.array(mfi);
            if (bx.smallEnd(1) <= domlo.y) {
                if (bc_fft[1].first == FFT::Boundary::even) {
                    ParallelFor(makeSlab(bx2,1,domlo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j-1,k) =  phi_arr(i,j,k);
                    });
                } else if (bc_fft[1].first == FFT::Boundary::odd) {
                    ParallelFor(makeSlab(bx2,1,domlo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j-1,k) =  -phi_arr(i,j,k);
                    });
                }
            } // lo y
            if (bx.bigEnd(1) >= domhi.y) {
                if (bc_fft[1].second == FFT::Boundary::even) {
                    ParallelFor(makeSlab(bx2,1,domhi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j+1,k) =  phi_arr(i,j,k);
                    });
                } else if (bc_fft[1].second == FFT::Boundary::odd) {
                    ParallelFor(makeSlab(bx2,1,domhi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j+1,k) =  -phi_arr(i,j,k);
                    });
                }
            } // hi y

        } // mfi
    } // not periodic in y

    auto bc_type_lo = m_domain_bcs_type[Orientation(2,Orientation::low)];
    auto bc_type_hi = m_domain_bcs_type[Orientation(2,Orientation::high)];

    for (MFIter mfi(phi,true); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        Box gbx(bx); gbx.grow(0,1); gbx.grow(1,1);
        const Array4<Real>& phi_arr = phi.array(mfi);
        if (bx.smallEnd(2) <= domlo.z) {
            ParallelFor(makeSlab(gbx,2,domlo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                phi_arr(i,j,k-1) = phi_arr(i,j,k);
            });
        } // lo z
        if (bx.bigEnd(2) >= domhi.z) {
            if (bc_fft[2].second == FFT::Boundary::even) {
                ParallelFor(makeSlab(gbx,2,domhi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    phi_arr(i,j,k+1) =  phi_arr(i,j,k);
                });
            } else if (bc_fft[2].second == FFT::Boundary::odd) {
                ParallelFor(makeSlab(gbx,2,domhi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    phi_arr(i,j,k+1) =  -phi_arr(i,j,k);
                });
            }
        } // hi z
    } // mfi

    phi.FillBoundary(m_geom.periodicity());
}

void TerrainPoisson::getFluxes (MultiFab& phi,
                                Array<MultiFab,AMREX_SPACEDIM>& fluxes)
{
    auto const& dxinv = m_geom.InvCellSizeArray();

    auto const& x   = phi.const_arrays();
    auto const& zpa = m_zphys->const_arrays();

    apply_bcs(phi);

    auto const& fx = fluxes[0].arrays();
    ParallelFor(fluxes[0], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        fx[b](i,j,k) = terrpoisson_flux_x(i,j,k,x[b],zpa[b],dxinv[0]);
    });

    auto const& fy = fluxes[1].arrays();
    ParallelFor(fluxes[1], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        fy[b](i,j,k) = terrpoisson_flux_y(i,j,k,x[b],zpa[b],dxinv[1]);
    });

    auto const& fz = fluxes[2].arrays();
    ParallelFor(fluxes[2], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        fz[b](i,j,k) = terrpoisson_flux_z(i,j,k,x[b],zpa[b],dxinv[0],dxinv[1]);
    });
}

void TerrainPoisson::assign (MultiFab& lhs, MultiFab const& rhs)
{
    MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
}

void TerrainPoisson::scale (MultiFab& lhs, Real fac)
{
    lhs.mult(fac);
}

Real TerrainPoisson::dotProduct (MultiFab const& v1, MultiFab const& v2)
{
    return MultiFab::Dot(v1, 0, v2, 0, 1, 0);
}

void TerrainPoisson::increment (MultiFab& lhs, MultiFab const& rhs, Real a)
{
    MultiFab::Saxpy(lhs, a, rhs, 0, 0, 1, 0);
}

void TerrainPoisson::linComb (MultiFab& lhs, Real a, MultiFab const& rhs_a,
                              Real b, MultiFab const& rhs_b)
{
    MultiFab::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, 1, 0);
}


MultiFab TerrainPoisson::makeVecRHS ()
{
    return MultiFab(m_grids, m_dmap, 1, 0);
}

MultiFab TerrainPoisson::makeVecLHS ()
{
    return MultiFab(m_grids, m_dmap, 1, 1);
}

Real TerrainPoisson::norm2 (MultiFab const& v)
{
    return v.norm2();
}

void TerrainPoisson::precond (MultiFab& lhs, MultiFab const& rhs)
{
#ifdef ERF_USE_FFT
    if (m_use_precond)
    {
        // Make a version that isn't constant
        MultiFab& rhs_tmp = const_cast<MultiFab&>(rhs);

        lhs.setVal(0.);
        m_2D_fft_precond->solve(lhs, rhs_tmp, m_stretched_dz_d);

#if 0
        AMREX_ASSERT(m_zphys_fft.local_size() <= 1);
        FArrayBox const* zfab = nullptr;
        if (m_zphys_fft.local_size() == 1) {
            zfab = m_zphys_fft.fabPtr(m_zphys_fft.IndexArray()[0]);
        }
        auto za = zfab ? zfab->const_array() : Array4<Real const>{};
        auto dxinv = m_geom.InvCellSize(0);
        auto dyinv = m_geom.InvCellSize(1);
        m_2D_fft_precond->solve(lhs, rhs,
            [=] AMREX_GPU_DEVICE (int ii, int jj, int k) -> Real
            {
                int i = 0; int j = 0;
                Real hzeta_inv_on_cc = Real(4.0) / ( (za(i,j,k+1) + za(i+1,j,k+1) + za(i,j+1,k+1) + za(i+1,j+1,k+1))
                                                    -(za(i,j,k  ) + za(i+1,j,k  ) + za(i,j+1,k  ) + za(i+1,j+1,k  )) );
                eal hzeta_inv_on_zlo = Real(8.0) / ( (za(i,j,k+1) + za(i+1,j,k+1) + za(i,j+1,k+1) + za(i+1,j+1,k+1))
                                                    -(za(i,j,k-1) + za(i+1,j,k-1) + za(i,j+1,k-1) + za(i+1,j+1,k-1)) );
                Real h_xi_on_zlo  = Real(0.5) * (za(i+1,j+1,k  ) + za(i+1,j,k  ) - za(i,j+1,k  ) - za(i,j,k  )) * dxinv;
                Real h_eta_on_zlo = Real(0.5) * (za(i+1,j+1,k  ) + za(i,j+1,k  ) - za(i+1,j,k  ) - za(i,j,k  )) * dyinv;
                return hzeta_inv_on_cc * (Real(1.0) + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
            },
            [=] AMREX_GPU_DEVICE (int ii, int jj, int k) -> Real
            {
                Real hzeta_inv_on_cc = Real(4.0) / ( (za(i,j,k+1) + za(i+1,j,k+1) + za(i,j+1,k+1) + za(i+1,j+1,k+1))
                                                    -(za(i,j,k  ) + za(i+1,j,k  ) + za(i,j+1,k  ) + za(i+1,j+1,k  )) );
                Real hzeta_inv_on_zhi = Real(8.0) / ( (za(i,j,k+2) + za(i+1,j,k+2) + za(i,j+1,k+2) + za(i+1,j+1,k+2))
                                                     -(za(i,j,k  ) + za(i+1,j,k  ) + za(i,j+1,k  ) + za(i+1,j+1,k  )) );
                Real h_xi_on_zhi  = Real(0.5) * (za(i+1,j+1,k+1) + za(i+1,j,k+1) - za(i,j+1,k+1) - za(i,j,k+1)) * dxinv;
                Real h_eta_on_zhi = Real(0.5) * (za(i+1,j+1,k+1) + za(i,j+1,k+1) - za(i+1,j,k+1) - za(i,j,k+1)) * dyinv;
                return hzeta_inv_on_cc * (Real(1.0) + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
            });
#endif
    } else
#endif
    {
        MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
    }
}

void TerrainPoisson::setToZero(MultiFab& v)
{
    v.setVal(0);
}
#endif
