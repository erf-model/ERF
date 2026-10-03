/**
 * \file ERF_TerrainPoisson.cpp
 */
#include "ERF_TerrainPoisson.H"
#include "ERF_SolverUtils.H"

using namespace amrex;

/**
 * @brief Construct a terrain-following Poisson operator.
 *
 * @param[in] geom Geometry defining the domain.
 * @param[in] lev_geom Geometry of the level, used to determine the FFT boundary conditions.
 * @param[in] ba BoxArray for the grid hierarchy.
 * @param[in] dm Distribution mapping for the grid.
 * @param[in] domain_bcs_type Boundary condition types for the domain.
 * @param[in] stretched_dz_lev_d Device vector of stretched vertical grid spacings.
 * @param[in] ax Metric term ax.
 * @param[in] ay Metric term ay.
 * @param[in] az Metric term az.
 * @param[in] dJ Jacobian of the coordinate transformation.
 * @param[in] z_phys_nd Nodal physical height field.
 * @param[in] use_real_bcs Flag to use real boundary conditions.
 * @param[in] build_fft_precond Whether to set up the FFT preconditioner (FFT builds only).
 */
TerrainPoisson::TerrainPoisson (Geometry const& geom, Geometry const& lev_geom,
                                BoxArray const& ba,
                                DistributionMapping const& dm,
                                Array<std::string,2*AMREX_SPACEDIM>& domain_bcs_type,
                                Gpu::DeviceVector<Real>& stretched_dz_lev_d,
                                const MultiFab& ax, const MultiFab& ay,
                                const MultiFab& az, const MultiFab& dJ,
                                MultiFab const* z_phys_nd,
                                bool use_real_bcs,
                                bool build_fft_precond)
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
    Box bounding_box = ba.minimalBox();
    m_bc = get_terrain_bc(lev_geom,domain_bcs_type,bounding_box,use_real_bcs);

#ifdef ERF_USE_FFT
    if (build_fft_precond) {
        auto bc_fft = get_fft_bc(lev_geom,domain_bcs_type,bounding_box,use_real_bcs);
        m_2D_fft_precond = std::make_unique<FFT::PoissonHybrid<MultiFab>>(geom,bc_fft);
    }
#else
    amrex::ignore_unused(build_fft_precond);
#endif
}

void TerrainPoisson::usePrecond (bool use_precond_in)
{
    m_use_precond = use_precond_in;
}

void TerrainPoisson::setPrecondFunction (PrecondFn fn)
{
    m_precond_fn = std::move(fn);
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
                if (m_bc[0].first == TerrainBC::even) {
                    ParallelFor(makeSlab(bx,0,domlo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i-1,j,k) =  phi_arr(i,j,k);
                    });
                } else if (m_bc[0].first == TerrainBC::odd) {
                    ParallelFor(makeSlab(bx,0,domlo.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i-1,j,k) =  -phi_arr(i,j,k);
                    });
                }
            } // lo x
            if (bx.bigEnd(0) >= domhi.x) {
                if (m_bc[0].second == TerrainBC::even) {
                    ParallelFor(makeSlab(bx,0,domhi.x), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i+1,j,k) =  phi_arr(i,j,k);
                    });
                } else if (m_bc[0].second == TerrainBC::odd) {
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
                if (m_bc[1].first == TerrainBC::even) {
                    ParallelFor(makeSlab(bx2,1,domlo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j-1,k) =  phi_arr(i,j,k);
                    });
                } else if (m_bc[1].first == TerrainBC::odd) {
                    ParallelFor(makeSlab(bx2,1,domlo.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j-1,k) =  -phi_arr(i,j,k);
                    });
                }
            } // lo y
            if (bx.bigEnd(1) >= domhi.y) {
                if (m_bc[1].second == TerrainBC::even) {
                    ParallelFor(makeSlab(bx2,1,domhi.y), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j+1,k) =  phi_arr(i,j,k);
                    });
                } else if (m_bc[1].second == TerrainBC::odd) {
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
            if (m_bc[2].second == TerrainBC::even) {
                ParallelFor(makeSlab(gbx,2,domhi.z), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    phi_arr(i,j,k+1) =  phi_arr(i,j,k);
                });
            } else if (m_bc[2].second == TerrainBC::odd) {
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
    if (m_use_precond)
    {
        if (m_precond_fn) {
            m_precond_fn(lhs, rhs);
            return;
        }
#ifdef ERF_USE_FFT
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_2D_fft_precond != nullptr,
                                         "TerrainPoisson: the FFT preconditioner was not built");
        // Make a version that isn't constant
        MultiFab& rhs_tmp = const_cast<MultiFab&>(rhs);

        lhs.setVal(0.);
        m_2D_fft_precond->solve(lhs, rhs_tmp, m_stretched_dz_d);
#else
        amrex::Abort("TerrainPoisson: no preconditioner available; rebuild with FFT or set one with setPrecondFunction");
#endif
    } else
    {
        MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
    }
}

void TerrainPoisson::setToZero (MultiFab& v)
{
    v.setVal(0);
}
