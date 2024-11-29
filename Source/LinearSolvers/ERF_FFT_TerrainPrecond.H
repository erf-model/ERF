#ifndef ERF_FFT_TERRAIN_PRECOND_H_
#define ERF_FFT_TERRAIN_PRECOND_H_

#include <AMReX_FFT.H>
#include <AMReX_FFT_Poisson.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

namespace detail {
template <typename MF>
void fill_physbc (MF& mf, Geometry const& geom,
                  Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc);
}

/**
 * \brief 3D Preconditioner for terrain problems with periodic, Neumann or Dirichlet boundaries
 *  in the first two dimensions and Neumann in the last dimension.  The last dimension can have
 *  non-uniform mesh spacing.
 */
template <typename MF = MultiFab>
class PoissonTerrainPrecond
{
public:
    using T = typename MF::value_type;

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    PoissonTerrainPrecond (Geometry const& geom,
                           Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        : m_geom(geom), m_bc(bc)
    {
        bool periodic_xy = true;
        for (int idim = 0; idim < 2; ++idim) {
            periodic_xy = periodic_xy && (bc[idim].first == Boundary::periodic);
            AMREX_ALWAYS_ASSERT((bc[idim].first == Boundary::periodic &&
                                 bc[idim].second == Boundary::periodic) ||
                                (bc[idim].first != Boundary::periodic &&
                                 bc[idim].second != Boundary::periodic));
        }
        Info info{};
        info.setBatchMode(true);
        if (periodic_xy) {
            m_r2c = std::make_unique<R2C<typename MF::value_type>>(m_geom.Domain(), info);
        } else {
            m_r2x = std::make_unique<R2X<typename MF::value_type>> (m_geom.Domain(), m_bc, info);
        }
    }

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit PoissonTerrainPrecond (Geometry const& geom)
        : m_geom(geom),
          m_bc{AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::even,Boundary::even))},
          m_r2c(std::make_unique<R2C<typename MF::value_type>>
                (geom.Domain(), Info().setBatchMode(true)))
    {
        AMREX_ALWAYS_ASSERT(geom.isPeriodic(0) && geom.isPeriodic(1));
    }

    /*
     * \brief Solve del dot grad soln = rhs
     *
     * If soln has ghost cells, one layer of ghost cells will be filled
     * except for the corners of the physical domain where they are not used
     * in the cross stencil of the operator. The two MFs could be the same MF.
     */
    void solve (MF& soln, MF const& rhs, MF const& height);

    // This is public for cuda
    template <typename FA, typename DZ>
    void solve_z (FA& spmf, MF const& height, DZ const& dz);

private:

    template <typename DZ>
    void solve_doit (MF& soln, MF const& rhs, MF const& height, DZ const& dz); // has to be public for cuda

    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    std::unique_ptr<R2X<typename MF::value_type>> m_r2x;
    std::unique_ptr<R2C<typename MF::value_type>> m_r2c;
    MF m_spmf_r;
    FabArray<BaseFab<GpuComplex<T>>> m_spmf_c;
};

template <typename MF>
void PoissonTerrainPrecond<MF>::solve (MF& soln, MF const& rhs, MF const& height)
{
    auto delz = T(m_geom.CellSize(AMREX_SPACEDIM-1));
    solve_doit(soln, rhs, height, fft_poisson_detail::DZ<T>{delz});
}

template <typename MF>
template <typename DZ>
void PoissonTerrainPrecond<MF>::solve_doit (MF& soln, MF const& rhs, MF const& height, DZ const& dz)
{
    BL_PROFILE("FFT::PoissonTerrainPrecond::solve");

    AMREX_ASSERT(soln.is_cell_centered() && rhs.is_cell_centered());

    IntVect const& ng = amrex::elemwiseMin(soln.nGrowVect(), IntVect(1));

    if (m_r2c)
    {
        if (m_spmf_c.empty()) {
            Box cdomain = m_geom.Domain();
            cdomain.setBig(0,cdomain.length(0)/2);
            auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                        {AMREX_D_DECL(true,true,false)});
            DistributionMapping dm = detail::make_iota_distromap(cba.size());
            m_spmf_c.define(cba, dm, 1, 0);
        }

        m_r2c->forward(rhs, m_spmf_c);

        solve_z(m_spmf_c, height, dz);

        m_r2c->backward_doit(m_spmf_c, soln, ng, m_geom.periodicity());
    }
    else
    {
         if (m_r2x->m_cy.empty()) { // spectral data is real
            auto sba = amrex::decompose(m_geom.Domain(),ParallelContext::NProcsSub(),
                                        {AMREX_D_DECL(true,true,false)});
            DistributionMapping dm = detail::make_iota_distromap(sba.size());
            m_spmf_r.define(sba, dm, 1, 0);
            m_r2x->forward(rhs, m_spmf_r);
            solve_z(m_spmf_r, height, dz);
            m_r2x->backward(m_spmf_r, soln, ng, m_geom.periodicity());
        } else { // spectral data is complex. one of the first two dimensions is periodic.
            Box cdomain = m_geom.Domain();
            if (m_bc[0].first == Boundary::periodic) {
                cdomain.setBig(0,cdomain.length(0)/2);
            } else {
                cdomain.setBig(1,cdomain.length(1)/2);
            }
            auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                        {AMREX_D_DECL(true,true,false)});
            DistributionMapping dm = detail::make_iota_distromap(cba.size());
            m_spmf_c.define(cba, dm, 1, 0);
            m_r2x->forward(rhs, m_spmf_c);
            solve_z(m_spmf_c, height, dz);
            m_r2x->backward(m_spmf_c, soln, ng, m_geom.periodicity());
        }
    }

    detail::fill_physbc(soln, m_geom, m_bc);
} // solve

template <typename MF>
template <typename FA, typename DZ>
void PoissonTerrainPrecond<MF>::solve_z (FA& spmf, MF const& height, DZ const& dz)
{
    BL_PROFILE("PoissonTerrainPrecond::solve_z");

    auto facx = Math::pi<T>()/T(m_geom.Domain().length(0));
    auto facy = Math::pi<T>()/T(m_geom.Domain().length(1));

    if (m_bc[0].first == Boundary::periodic) { facx *= T(2); }
    if (m_bc[1].first == Boundary::periodic) { facy *= T(2); }

    auto dxfac = T(2)/T(m_geom.CellSize(0)*m_geom.CellSize(0));
    auto dyfac = T(2)/T(m_geom.CellSize(1)*m_geom.CellSize(1));
    auto scale = (m_r2x) ? m_r2x->scalingFactor() : m_r2c->scalingFactor();

    auto dxinv = T(m_geom.InvCellSize(0));
    auto dyinv = T(m_geom.InvCellSize(1));

    GpuArray<T,AMREX_SPACEDIM-1> offset{T(0),T(0)};
    for (int idim = 0; idim < AMREX_SPACEDIM-1; ++idim) {
        if (m_bc[idim].first == Boundary::odd &&
            m_bc[idim].second == Boundary::odd)
        {
            offset[idim] = T(1);
        }
        else if ((m_bc[idim].first == Boundary::odd &&
                  m_bc[idim].second == Boundary::even) ||
                 (m_bc[idim].first == Boundary::even &&
                  m_bc[idim].second == Boundary::odd))
        {
            offset[idim] = T(0.5);
        }
    }

    bool has_dirichlet = (offset[0] != T(0)) || (offset[1] != T(0));

    auto nz = m_geom.Domain().length(2);

    for (MFIter mfi(spmf); mfi.isValid(); ++mfi)
    {
        auto const& spectral = spmf.array(mfi);
        auto const& box = mfi.validbox();
        auto const& xybox = amrex::makeSlab(box, 2, 0);

        auto const zp = height.const_array(mfi);

#ifdef AMREX_USE_GPU
        FArrayBox tridiag_workspace(box,4);
        auto const& ald     = tridiag_workspace.array(0);
        auto const& bd      = tridiag_workspace.array(1);
        auto const& cud     = tridiag_workspace.array(2);
        auto const& scratch = tridiag_workspace.array(3);

        amrex::ParallelFor(xybox, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            T a = facx*(i+offset[0]);
            T b = facy*(j+offset[1]);
            T k2 = dxfac * (std::cos(a)-T(1))
                +  dyfac * (std::cos(b)-T(1));

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {
                Real hzeta_inv_on_cc = T(4.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                              -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                if(k==0) {

                    Real hzeta_inv_on_zhi = T(8.0) / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = T(0.5) * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = T(0.5) * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald(i,j,k) = T(0.0);
                    cud(i,j,k) = hzeta_inv_on_cc * (1.0 + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
                    bd(i,j,k) = k2 - ald(i,j,k) - cud(i,j,k);

                } else if (k == nz-1) {

                    Real hzeta_inv_on_zlo = T(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;
                    ald(i,j,k) = hzeta_inv_on_cc * (1.0 + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud(i,j,k) = 0.;
                    bd(i,j,k) = k2 - ald(i,j,k) - cud(i,j,k);
                    if (i == 0 && j == 0) {
                        bd(i,j,k) *= T(2.0);
                    }
                } else {
                    Real hzeta_inv_on_zlo = T(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;

                    Real hzeta_inv_on_zhi = T(8.0) / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = T(0.5) * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = T(0.5) * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald(i,j,k) = hzeta_inv_on_cc * (1.0 + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud(i,j,k) = hzeta_inv_on_cc * (1.0 + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
                    bd(i,j,k) = k2 - ald(i,j,k) - cud(i,j,k);

                }
            }

            scratch(i,j,0) = cud(i,j,0)/bd(i,j,0);
            spectral(i,j,0) = spectral(i,j,0)/bd(i,j,0);

            for (int k = 1; k < nz; k++) {
                if (k < nz-1) {
                    scratch(i,j,k) = cud(i,j,k) / (bd(i,j,k) - ald(i,j,k) * scratch(i,j,k-1));
                }
                spectral(i,j,k) = (spectral(i,j,k) - ald(i,j,k) * spectral(i,j,k - 1))
                    / (bd(i,j,k) - ald(i,j,k) * scratch(i,j,k-1));
            }

            for (int k = nz - 2; k >= 0; k--) {
                spectral(i,j,k) -= scratch(i,j,k) * spectral(i,j,k + 1);
            }

            for (int k = 0; k < nz; ++k) {
                spectral(i,j,k) *= scale;
            }
        });
        Gpu::streamSynchronize();

#else

        Gpu::DeviceVector<T> ald(nz);
        Gpu::DeviceVector<T> bd(nz);
        Gpu::DeviceVector<T> cud(nz);
        Gpu::DeviceVector<T> scratch(nz);

        amrex::LoopOnCpu(xybox, [&] (int i, int j, int)
        {
            T a = facx*(i+offset[0]);
            T b = facy*(j+offset[1]);
            T k2 = dxfac*(std::cos(a)-T(1))
                +  dyfac*(std::cos(b)-T(1));

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {

                Real hzeta_inv_on_cc = T(4.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                              -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );

                if(k==0) {

                    Real hzeta_inv_on_zhi = T(8.0) / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = T(0.5) * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = T(0.5) * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald[k] = 0.;
                    cud[k] = hzeta_inv_on_cc * (T(1.0) + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
                    bd[k]  = k2 -ald[k]-cud[k];

                } else if (k == nz-1) {

                    Real hzeta_inv_on_zlo = T(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = T(0.5) * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = T(0.5) * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;

                    ald[k] = hzeta_inv_on_cc * (T(1.0) + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud[k] = 0.;
                    bd[k]  = k2 -ald[k]-cud[k];

                    if (i == 0 && j == 0) {
                        bd[k] *= T(2.0);
                    }
                } else {

                    Real hzeta_inv_on_zlo = T(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = T(0.5) * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = T(0.5) * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;

                    Real hzeta_inv_on_zhi = T(8.0) / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = T(0.5) * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = T(0.5) * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald[k] = hzeta_inv_on_cc * (T(1.0) + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud[k] = hzeta_inv_on_cc * (T(1.0) + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
                    bd[k]  = k2 - ald[k] - cud[k];
                }
            }

            scratch[0] = cud[0]/bd[0];
            spectral(i,j,0) = spectral(i,j,0)/bd[0];

            for (int k = 1; k < nz; k++) {
                if (k < nz-1) {
                    scratch[k] = cud[k] / (bd[k] - ald[k] * scratch[k-1]);
                }
                spectral(i,j,k) = (spectral(i,j,k) - ald[k] * spectral(i,j,k - 1))
                    / (bd[k] - ald[k] * scratch[k-1]);
            }

            for (int k = nz - 2; k >= 0; k--) {
                spectral(i,j,k) -= scratch[k] * spectral(i,j,k + 1);
            }

            for (int k = 0; k < nz; ++k) {
                spectral(i,j,k) *= scale;
            }
        });
#endif
    } // mfi
} // solve_z

} // namespace
#endif
