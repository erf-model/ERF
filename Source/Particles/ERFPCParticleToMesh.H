/*! \file ERFPCParticleToMesh.H
 *  \brief Template implementation of ERFPC::ERFPCParticleToMesh
 *
 *  Separated from ERFPC.H to keep the header lightweight.
 *  Include this in .cpp files that call ERFPCParticleToMesh().
 */

#ifndef ERFPC_PARTICLE_TO_MESH_H_
#define ERFPC_PARTICLE_TO_MESH_H_

#include <ERFPC.H>
#include <AMReX_ParticleInterpolators.H>

template<typename ValueFunc>
void ERFPC::ERFPCParticleToMesh(amrex::MultiFab& a_mf, int a_lev, int a_comp,
                                 ValueFunc&& value_func) const
{
    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const amrex::Real* zlevels = m_zlevels_d.empty() ? nullptr : m_zlevels_d.data();
    const int nz_levels = static_cast<int>(m_zlevels_d.size());

    a_mf.setVal(0.0);

    amrex::ParticleToMesh(*this, a_mf, a_lev,
        [=] AMREX_GPU_DEVICE (const typename ERFPC::ParticleTileType::ConstParticleTileDataType& ptd,
                              int i, amrex::Array4<amrex::Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            auto particle_value = value_func(ptd, i);

            if (zlevels) {
                amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + amrex::Real(0.5);
                int ix = static_cast<int>(amrex::Math::floor(lx)) - 1;
                amrex::Real wx1 = lx - static_cast<amrex::Real>(ix + 1);
                amrex::Real wx0 = amrex::Real(1.0) - wx1;

                amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + amrex::Real(0.5);
                int iy = static_cast<int>(amrex::Math::floor(ly)) - 1;
                amrex::Real wy1 = ly - static_cast<amrex::Real>(iy + 1);
                amrex::Real wy0 = amrex::Real(1.0) - wy1;

                int kz = p.idata(ERFParticlesIntIdxAoS::k);
                int nz = nz_levels - 1;

                amrex::Real dz_k = zlevels[kz+1] - zlevels[kz];
                amrex::Real z_center_k = amrex::Real(0.5) * (zlevels[kz] + zlevels[kz+1]);

                int kz_lo, kz_hi;
                amrex::Real wz_lo, wz_hi;
                amrex::Real dz_lo, dz_hi;

                if (p.pos(2) >= z_center_k) {
                    kz_lo = kz;
                    kz_hi = kz + 1;
                    dz_lo = dz_k;
                    if (kz_hi < nz) {
                        dz_hi = zlevels[kz_hi+1] - zlevels[kz_hi];
                        amrex::Real z_center_hi = amrex::Real(0.5) * (zlevels[kz_hi] + zlevels[kz_hi+1]);
                        amrex::Real span = z_center_hi - z_center_k;
                        wz_hi = (p.pos(2) - z_center_k) / span;
                        wz_lo = amrex::Real(1.0) - wz_hi;
                    } else {
                        kz_hi = kz;
                        dz_hi = dz_k;
                        wz_lo = amrex::Real(1.0);
                        wz_hi = amrex::Real(0.0);
                    }
                } else {
                    kz_lo = kz - 1;
                    kz_hi = kz;
                    dz_hi = dz_k;
                    if (kz_lo >= 0) {
                        dz_lo = zlevels[kz_lo+1] - zlevels[kz_lo];
                        amrex::Real z_center_lo = amrex::Real(0.5) * (zlevels[kz_lo] + zlevels[kz_lo+1]);
                        amrex::Real span = z_center_k - z_center_lo;
                        wz_lo = (z_center_k - p.pos(2)) / span;
                        wz_hi = amrex::Real(1.0) - wz_lo;
                    } else {
                        kz_lo = kz;
                        dz_lo = dz_k;
                        wz_lo = amrex::Real(0.0);
                        wz_hi = amrex::Real(1.0);
                    }
                }

                amrex::Real inv_vol_lo = dxi[0] * dxi[1] / dz_lo;
                amrex::Real inv_vol_hi = dxi[0] * dxi[1] / dz_hi;

                amrex::Real pval_lo = particle_value * wz_lo * inv_vol_lo;
                amrex::Real pval_hi = particle_value * wz_hi * inv_vol_hi;

                amrex::Gpu::Atomic::AddNoRet(&rho(ix  , iy  , kz_lo, a_comp), wx0*wy0*pval_lo);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix+1, iy  , kz_lo, a_comp), wx1*wy0*pval_lo);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix  , iy+1, kz_lo, a_comp), wx0*wy1*pval_lo);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix+1, iy+1, kz_lo, a_comp), wx1*wy1*pval_lo);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix  , iy  , kz_hi, a_comp), wx0*wy0*pval_hi);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix+1, iy  , kz_hi, a_comp), wx1*wy0*pval_hi);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix  , iy+1, kz_hi, a_comp), wx0*wy1*pval_hi);
                amrex::Gpu::Atomic::AddNoRet(&rho(ix+1, iy+1, kz_hi, a_comp), wx1*wy1*pval_hi);
            } else {
                amrex::ParticleInterpolator::Linear interp(p, plo, dxi);
                amrex::ParticleReal inv_cell_volume = dxi[0] * dxi[1] * dxi[2];
                interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                    [=] AMREX_GPU_DEVICE (const ERFPC::ParticleType&, int) {
                        return particle_value * inv_cell_volume;
                    });
            }
        });
}

#endif
