#ifndef ERF_FFT_TERRAIN_PRECOND_H_
#define ERF_FFT_TERRAIN_PRECOND_H_

#include <AMReX_FFT.H>
#include <AMReX_FFT_Poisson.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

/**
 * \brief 3D Preconditioner for terrain problems with periodic boundaries in the first two
 * dimensions and Neumann in the last dimension.
 */
template <typename MF = MultiFab>
class PoissonTerrainPrecond
{
public:
    using T = typename MF::value_type;

    template <typename FA=MF, std::enable_if_t<IsFabArray_v<FA>,int> = 0>
    explicit PoissonTerrainPrecond (Geometry const& geom)
        : m_geom(geom), m_r2c(geom.Domain(), Info().setBatchMode(true))
    {
        AMREX_ALWAYS_ASSERT(geom.isPeriodic(0) && geom.isPeriodic(1));
    }

    void solve (MF& soln, MF const& rhs, MF const& height);

    template <typename DZ>
    void solve_doit (MF& soln, MF const& rhs, MF const& height, DZ const& dz); // has to be public for cuda

private:
    Geometry m_geom;
    R2C<typename MF::value_type, Direction::both> m_r2c;
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

#if (AMREX_SPACEDIM < 3)
    amrex::ignore_unused(soln, rhs, dz);
#else
    auto facx = T(2)*Math::pi<T>()/T(m_geom.ProbLength(0));
    auto facy = T(2)*Math::pi<T>()/T(m_geom.ProbLength(1));

    auto dx =T(m_geom.CellSize(0));
    auto dy = T(m_geom.CellSize(1));

    auto dxinv = T(m_geom.InvCellSize(0));
    auto dyinv = T(m_geom.InvCellSize(1));

    auto scale = T(1.0)/(T(m_geom.Domain().length(0)) *
                         T(m_geom.Domain().length(1)));
    auto ny = m_geom.Domain().length(1);
    auto nz = m_geom.Domain().length(2);

    Box cdomain = m_geom.Domain();
    cdomain.setBig(0,cdomain.length(0)/2);
    auto cba = amrex::decompose(cdomain, ParallelContext::NProcsSub(),
                                {AMREX_D_DECL(true,true,false)});
    DistributionMapping dm = detail::make_iota_distromap(cba.size());
    FabArray<BaseFab<GpuComplex<T> > > spmf(cba, dm, 1, 0);

    m_r2c.forward(rhs, spmf);

    for (MFIter mfi(spmf); mfi.isValid(); ++mfi)
    {
        auto const& spectral = spmf.array(mfi);
        auto const& box = mfi.validbox();
        auto const& xybox = amrex::makeSlab(box, 2, 0);

        auto const zp = height.const_array(mfi);

#ifdef AMREX_USE_GPU
        // xxxxx TODO: We need to explore how to optimize this
        // function. Maybe we can use cusparse. Maybe we should make
        // z-direction to be the unit stride direction.

        FArrayBox tridiag_workspace(box,4);
        auto const& ald = tridiag_workspace.array(0);
        auto const& bd = tridiag_workspace.array(1);
        auto const& cud = tridiag_workspace.array(2);
        auto const& scratch = tridiag_workspace.array(3);

        amrex::ParallelFor(xybox, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            T a = facx*i;
            T b = (j < ny/2) ? facy*j : facy*(ny-j);

            T k2 = T(2)*(std::cos(a*dx)-T(1))/(dx*dx)
                +  T(2)*(std::cos(b*dy)-T(1))/(dy*dy);

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {
                Real hzeta_inv_on_cc = 4.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                              -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                if(k==0) {

                    Real hzeta_inv_on_zhi = 8.0 / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = 0.5 * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = 0.5 * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald(i,j,k) = 0.;
                    cud(i,j,k) = hzeta_inv_on_cc * (1.0 + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
                    bd(i,j,k) = k2 - ald(i,j,k) - cud(i,j,k);

                } else if (k == nz-1) {

                    Real hzeta_inv_on_zlo = 8.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;
                    ald(i,j,k) = hzeta_inv_on_cc * (1.0 + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud(i,j,k) = 0.;
                    bd(i,j,k) = k2 - ald(i,j,k) - cud(i,j,k);
                    if (i == 0 && j == 0) {
                        bd(i,j,k) *= 2.0;
                    }
                } else {
                    Real hzeta_inv_on_zlo = 8.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;

                    Real hzeta_inv_on_zhi = 8.0 / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = 0.5 * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = 0.5 * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

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

        Gpu::DeviceVector<GpuComplex<Real>> ald(nz);
        Gpu::DeviceVector<GpuComplex<Real>> bd(nz);
        Gpu::DeviceVector<GpuComplex<Real>> cud(nz);
        Gpu::DeviceVector<GpuComplex<Real>> scratch(nz);

        amrex::LoopOnCpu(xybox, [&] (int i, int j, int)
        {
            T a = facx*i;
            T b = (j < ny/2) ? facy*j : facy*(ny-j);

            T k2 = T(2)*(std::cos(a*dx)-T(1))/(dx*dx)
                +  T(2)*(std::cos(b*dy)-T(1))/(dy*dy);

            // Tridiagonal solve with homogeneous Neumann
            for(int k=0; k < nz; k++) {

                Real hzeta_inv_on_cc = 4.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                              -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );

                if(k==0) {

                    Real hzeta_inv_on_zhi = 8.0 / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = 0.5 * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = 0.5 * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald[k] = 0.;
                    cud[k] = hzeta_inv_on_cc * (1.0 + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
                    bd[k]  = k2 -ald[k]-cud[k];

                } else if (k == nz-1) {

                    Real hzeta_inv_on_zlo = 8.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;

                    ald[k] = hzeta_inv_on_cc * (1.0 + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud[k] = 0.;
                    bd[k]  = k2 -ald[k]-cud[k];

                    if (i == 0 && j == 0) {
                        bd[k] *= 2.0;
                    }
                } else {

                    Real hzeta_inv_on_zlo = 8.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
                    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k  ) + zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i,j,k  )) * dxinv;
                    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k  ) + zp(i,j+1,k  ) - zp(i+1,j,k  ) - zp(i,j,k  )) * dyinv;

                    Real hzeta_inv_on_zhi = 8.0 / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
                    Real h_xi_on_zhi  = 0.5 * (zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1)) * dxinv;
                    Real h_eta_on_zhi = 0.5 * (zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1)) * dyinv;

                    ald[k] = hzeta_inv_on_cc * (1.0 + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo) * hzeta_inv_on_zlo;
                    cud[k] = hzeta_inv_on_cc * (1.0 + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi) * hzeta_inv_on_zhi;
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
    }

    m_r2c.backward(spmf, soln);
#endif
}

}

#endif
