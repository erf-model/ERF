/*! \file ERFPCParticleToMesh.H
 *  \brief Template implementation of ERFPC::ERFPCParticleToMesh
 *
 *  Separated from ERFPC.H to keep the header lightweight.
 *  Include this in .cpp files that call ERFPCParticleToMesh().
 */

#ifndef ERFPC_PARTICLE_TO_MESH_H_
#define ERFPC_PARTICLE_TO_MESH_H_

#include <ERFPC.H>
#include <ERF_TerrainConversion.H>

template<typename ValueFunc>
void ERFPC::ERFPCParticleToMesh(amrex::MultiFab& a_mf,
                                const amrex::MultiFab& a_z_phys_nd,
                                int a_lev, int a_comp,
                                ValueFunc&& value_func) const
{
    using namespace amrex;

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const Real dx  = geom.CellSize(0);
    const Real dy  = geom.CellSize(1);
    const Box  domain = geom.Domain();
    const int  k_max  = domain.bigEnd(AMREX_SPACEDIM-1);

    a_mf.setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParConstIterType pti(*this, a_lev); pti.isValid(); ++pti)
    {
        const int grid = pti.index();
        const auto& ptile = ParticlesAt(a_lev, pti);
        const auto& aos = ptile.GetArrayOfStructs();
        const int np = aos.numParticles();
        if (np == 0) { continue; }

        auto rho     = a_mf[grid].array();
        auto zheight = a_z_phys_nd[grid].array();

        const auto ptd = ptile.getConstParticleTileData();

        ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
        {
            const auto& p = ptd.m_aos[i];
            if (p.id() <= 0) { return; }

            const Real pval = static_cast<Real>(value_func(ptd, i));

            const Real lx = static_cast<Real>((p.pos(0) - plo[0]) * dxi[0]) + Real(0.5);
            const int  ix = static_cast<int>(Math::floor(lx)) - 1;
            const Real wx1 = lx - static_cast<Real>(ix + 1);
            const Real wx0 = Real(1.0) - wx1;

            const Real ly = static_cast<Real>((p.pos(1) - plo[1]) * dxi[1]) + Real(0.5);
            const int  iy = static_cast<int>(Math::floor(ly)) - 1;
            const Real wy1 = ly - static_cast<Real>(iy + 1);
            const Real wy0 = Real(1.0) - wy1;

            const int kz = static_cast<int>(amrex::Math::floor((p.pos(AMREX_SPACEDIM-1) - plo[AMREX_SPACEDIM-1])
                                                               * dxi[AMREX_SPACEDIM-1]));
            const int ix_p = static_cast<int>(Math::floor((p.pos(0) - plo[0]) * dxi[0]));
            const int iy_p = static_cast<int>(Math::floor((p.pos(1) - plo[1]) * dxi[1]));
            const Real fx  = static_cast<Real>((p.pos(0) - plo[0]) * dxi[0]) - static_cast<Real>(ix_p);
            const Real fy  = static_cast<Real>((p.pos(1) - plo[1]) * dxi[1]) - static_cast<Real>(iy_p);

            auto zface = [=] AMREX_GPU_DEVICE (int k_face) -> Real {
                return zheight(ix_p  , iy_p  , k_face) * (Real(1.0)-fx) * (Real(1.0)-fy)
                     + zheight(ix_p+1, iy_p  , k_face) *            fx  * (Real(1.0)-fy)
                     + zheight(ix_p  , iy_p+1, k_face) * (Real(1.0)-fx) *            fy
                     + zheight(ix_p+1, iy_p+1, k_face) *            fx  *            fy;
            };

            // Effective k bounds for safe zface access: must respect both the
            // geometry domain (k_max) and the local fab's z extent. With AMR
            // partial-z fabs (especially with many small boxes from large rank
            // counts), a particle at the top/bottom of its fab can have kz
            // such that kz+2 / kz-1 read outside the fab. Use the fallback
            // (one-sided) span formula in those cases.
            const int kz_hi_safe = amrex::min(k_max, zheight.end[2] - 2);
            const int kz_lo_safe = amrex::max(0,     zheight.begin[2] + 1);

            const Real z_lo_k = zface(kz);
            const Real z_hi_k = zface(kz + 1);
            const Real z_ctr_k = Real(0.5) * (z_lo_k + z_hi_k);

            // Particle pos(2) is computational zeta; convert to physical z so
            // the cell-center comparison and weight span are in the same units
            // as zface() (which is built from the height array).
            const Real p_z = static_cast<Real>(ERF::ParticlePos::z_from_zeta(
                p.pos(0), p.pos(1), p.pos(AMREX_SPACEDIM-1),
                plo, dxi, zheight));

            int  kz_lo, kz_hi;
            Real wz_lo, wz_hi;
            if (p_z >= z_ctr_k) {
                kz_lo = kz;
                kz_hi = kz + 1;
                Real span;
                if (kz < kz_hi_safe) {
                    const Real z_hi2    = zface(kz + 2);
                    const Real z_ctr_hi = Real(0.5) * (z_hi_k + z_hi2);
                    span = z_ctr_hi - z_ctr_k;
                } else {
                    span = z_hi_k - z_lo_k;
                }
                wz_hi = (p_z - z_ctr_k) / span;
                wz_lo = Real(1.0) - wz_hi;
            } else {
                kz_hi = kz;
                kz_lo = kz - 1;
                Real span;
                if (kz > kz_lo_safe) {
                    const Real z_lo2    = zface(kz - 1);
                    const Real z_ctr_lo = Real(0.5) * (z_lo2 + z_lo_k);
                    span = z_ctr_k - z_ctr_lo;
                } else {
                    span = z_hi_k - z_lo_k;
                }
                wz_lo = (z_ctr_k - p_z) / span;
                wz_hi = Real(1.0) - wz_lo;
            }

            // Per-corner node-averaged cell thickness for volume normalization.
            auto dz_cell = [=] AMREX_GPU_DEVICE (int ic, int jc, int kc) -> Real {
                const Real dz = Real(0.25) * (
                      zheight(ic  , jc  , kc+1) - zheight(ic  , jc  , kc)
                    + zheight(ic+1, jc  , kc+1) - zheight(ic+1, jc  , kc)
                    + zheight(ic  , jc+1, kc+1) - zheight(ic  , jc+1, kc)
                    + zheight(ic+1, jc+1, kc+1) - zheight(ic+1, jc+1, kc) );
                return dz;
            };

            const Real inv_dxdy = Real(1.0) / (dx * dy);

            auto deposit = [=] AMREX_GPU_DEVICE (int ic, int jc, int kc, Real w_horiz, Real w_z)
            {
                // Clamp kc only for dz lookup; still write to ghost rho cell.
                const int kc_dz = amrex::max(0, amrex::min(k_max, kc));
                const Real dz = dz_cell(ic, jc, kc_dz);
                const Real contrib = pval * w_horiz * w_z * inv_dxdy / dz;
                Gpu::Atomic::AddNoRet(&rho(ic, jc, kc, a_comp), contrib);
            };

            deposit(ix  , iy  , kz_lo, wx0*wy0, wz_lo);
            deposit(ix+1, iy  , kz_lo, wx1*wy0, wz_lo);
            deposit(ix  , iy+1, kz_lo, wx0*wy1, wz_lo);
            deposit(ix+1, iy+1, kz_lo, wx1*wy1, wz_lo);
            deposit(ix  , iy  , kz_hi, wx0*wy0, wz_hi);
            deposit(ix+1, iy  , kz_hi, wx1*wy0, wz_hi);
            deposit(ix  , iy+1, kz_hi, wx0*wy1, wz_hi);
            deposit(ix+1, iy+1, kz_hi, wx1*wy1, wz_hi);
        });
    }

    a_mf.SumBoundary(geom.periodicity());
}

#endif
