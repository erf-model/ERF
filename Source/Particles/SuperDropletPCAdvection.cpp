#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "IndexDefines.H"

using namespace amrex;

/*! Evolve particles for one time step */
void SuperDropletPC::AdvectParticles ( int                                        a_lev,
                                       Real                                       a_dt,
                                       Vector<Vector<MultiFab>>&                  a_flow_vars,
                                       const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd,
                                       const bool                                 a_advect,
                                       const bool                                 a_grav )
{
    BL_PROFILE("SuperDropletPC::AdvectParticles()");

    MultiFab* flow_vel( &a_flow_vars[a_lev][Vars::xvel] );
    const std::unique_ptr<MultiFab>& z_height = a_z_phys_nd[a_lev];

    AMREX_ASSERT(OK(a_lev, a_lev, flow_vel[0].nGrow()-1));
    AMREX_ASSERT(a_lev >= 0 && a_lev < GetParticles().size());
    AMREX_D_TERM(AMREX_ASSERT(flow_vel[0].nGrow() >= 1);,
                 AMREX_ASSERT(flow_vel[1].nGrow() >= 1);,
                 AMREX_ASSERT(flow_vel[2].nGrow() >= 1););

    AMREX_D_TERM(AMREX_ASSERT(!flow_vel[0].contains_nan());,
                 AMREX_ASSERT(!flow_vel[1].contains_nan());,
                 AMREX_ASSERT(!flow_vel[2].contains_nan()););

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    for (int ipass = 0; ipass < 2; ipass++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            auto& soa  = ptile.GetStructOfArrays();
            const int n = aos.numParticles();
            auto *p_pbox = aos().data();

            Array<ParticleReal*,AMREX_SPACEDIM> v_ptr;
            v_ptr[0] = soa.GetRealData(SuperDropletsRealIdxSoA::vx).data();
            v_ptr[1] = soa.GetRealData(SuperDropletsRealIdxSoA::vy).data();
            v_ptr[2] = soa.GetRealData(SuperDropletsRealIdxSoA::vz).data();

            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&(flow_vel[0][grid]),
                                                                  &(flow_vel[1][grid]),
                                                                  &(flow_vel[2][grid])) };
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
                const umacarr {{AMREX_D_DECL((*fab[0]).array(),
                                             (*fab[1]).array(),
                                             (*fab[2]).array() )}};

            bool use_terrain = (z_height != nullptr);
            auto zheight = use_terrain ? (*z_height)[grid].array() : Array4<Real>{};

            int rt_offset = SuperDropletsRealIdxSoA::ncomps;
            auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();

            MatPropFunctions::TerminalVelocity term_vel { m_vapour_mat->term_vel_a(),
                                                          m_vapour_mat->term_vel_b() };

            ParallelFor(n, [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }

                ParticleReal v[AMREX_SPACEDIM];
                if (a_advect) {
                    if (use_terrain) {
                        mac_interpolate_mapped_z(p, plo, dxi, umacarr, zheight, v);
                    } else {
                        mac_interpolate(p, plo, dxi, umacarr, v);
                    }
                } else {
                    v[0] = v[1] = v[2] = 0.0;
                }

                if (a_grav) {
                    ParticleReal terminal_vel = term_vel(radius_ptr[i]);
                    terminal_vel *= (1.0 - std::exp(-p.pos(2)*dxi[2]));
                    v[2] -= terminal_vel;
                }

                if (ipass == 0) {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                        v_ptr[dim][i] = p.pos(dim);
                        p.pos(dim) += static_cast<ParticleReal>(ParticleReal(0.5)*a_dt*v[dim]);
                    }
                    // Update z-coordinate carried by the particle
                    update_location_idata(p,plo,dxi,zheight);
                } else {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                        p.pos(dim) = v_ptr[dim][i] + static_cast<ParticleReal>(a_dt*v[dim]);
                        v_ptr[dim][i] = v[dim];
                    }
                    // Update z-coordinate carried by the particle
                    update_location_idata(p,plo,dxi,zheight);
                }
            });
        }
    }

}

#endif
