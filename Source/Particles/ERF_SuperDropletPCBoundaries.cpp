#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_InterpolationUtils.H"

using namespace amrex;
using namespace SDPCDefn;


/*! Handle the boundaries for the particles */
void SuperDropletPC::applyBoundaryTreatment ( int                   a_lev,
                                              const Vector<MFPtr>&  a_z_phys_nd,
                                              const BCTypeArr&      a_bctypes,
                                              const bool            a_recycle )
{
    BL_PROFILE("SuperDropletPC::applyBoundaryTreatment()");
    const MFPtr& z_height = a_z_phys_nd[a_lev];

    const auto ctx = buildProcessContext(a_lev);
    const auto save_inac = m_save_inactive;

    // number of super-droplets per cell
    int num_sd_per_cell = m_num_sd_per_cell;
    // number of physical particles per cell
    Real num_par_per_cell = zero;
    for (int i = 0; i < m_num_initializations; i++) {
        num_par_per_cell += m_initializations[i]->numParticlesPerCell(ctx.cell_volume);
    }
    auto multiplicity = (num_sd_per_cell > 0 ? num_par_per_cell / num_sd_per_cell : zero);

    Long num_deactivated_particles = 0;

    forEachParticleTile(a_lev, ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        auto zheight = (*z_height)[grid].array();

        Gpu::Buffer<Long> deactivated_particles({0});
        auto* deactivated_particles_ptr = deactivated_particles.data();

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            // check for ground impact
            {
                auto z_ground = ctx.plo[2];
                {
                    auto iv = getParticleCell(p, ctx.plo, ctx.dxi, ctx.domain);
                    z_ground = zheight(iv[0],iv[1],ctx.domain.smallEnd(2));
                }
                if (p.pos(2) < z_ground) {
                    p.pos(2) = z_ground + Real(0.01)*ctx.dx[2];
                    ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;
                    ptrs.active_ptr[i] = 0;
                    if ((!a_recycle) && (!save_inac)) { p.id() = -1; }
                    Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));
                    update_location_idata(p,ctx.plo,ctx.dxi,zheight);
                }
            }

            // check for top boundary exits
            {
                auto z_roof = ctx.phi[2];
                {
                    auto iv = getParticleCell(p, ctx.plo, ctx.dxi, ctx.domain);
                    z_roof = zheight(iv[0],iv[1],ctx.domain.bigEnd(2)+1);
                }
                if (p.pos(2) > z_roof) {
                    p.pos(2) = z_roof - ctx.dx[2];
                    ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;
                    ptrs.active_ptr[i] = 0;
                    if ((!a_recycle) && (!save_inac)) { p.id() = -1; }
                    Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));
                    update_location_idata(p,ctx.plo,ctx.dxi,zheight);
                }
            }

            // check if particles have exited the domain along x and y
            for (int d = 0; d < 2; d++) {

                // domain bounds
                auto x_min = ctx.plo[d];
                auto x_max = ctx.phi[d];

                auto bc_lo = a_bctypes[Orientation(d,Orientation::low)];
                auto bc_hi = a_bctypes[Orientation(d,Orientation::high)];

                if (p.pos(d) < x_min) {

                    if ((bc_lo == ERF_BC::slip_wall) || (bc_lo == ERF_BC::no_slip_wall)) {

                        p.pos(d) = x_min + Real(0.01)*ctx.dx[d];
                        ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;
                        ptrs.active_ptr[i] = 0;
                        if ((!a_recycle) && (!save_inac)) { p.id() = -1; }
                        Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));

                    } else {

                        auto delta = x_min - p.pos(d);
                        p.pos(d) = x_max - delta;
                        if (!ctx.is_periodic[d]) {
                            ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;

                            for (int ctr = 0; ctr < ctx.num_species; ctr++) {
                                ptrs.sp_mass_ptrs[ctr][i] = zero;
                            }

                            ptrs.mult_ptr[i] = multiplicity;

                            SuperDropletPC::updateParticleAttributes(
                                i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                                ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);
                        }

                    }

                } else if (p.pos(d) > x_max) {

                    if ((bc_hi == ERF_BC::slip_wall) || (bc_hi == ERF_BC::no_slip_wall)) {

                        p.pos(d) = x_max - Real(0.01)*ctx.dx[d];
                        ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;
                        ptrs.active_ptr[i] = 0;
                        if ((!a_recycle) && (!save_inac)) { p.id() = -1; }
                        Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));

                    } else {

                        auto delta = p.pos(d) - x_max;
                        p.pos(d) = x_min + delta;
                        if (!ctx.is_periodic[d]) {
                            ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;

                            for (int ctr = 0; ctr < ctx.num_species; ctr++) {
                                ptrs.sp_mass_ptrs[ctr][i] = zero;
                            }

                            ptrs.mult_ptr[i] = multiplicity;

                            SuperDropletPC::updateParticleAttributes(
                                i, ptrs.radius_ptr, ptrs.mass_ptr, ctx.idx_water, ctx.rho_water,
                                ctx.num_species, ctx.num_aerosols, ptrs.sp_sol_arr, ptrs.ae_sol_arr,
                                ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs, ptrs.sp_rho_arr, ptrs.ae_rho_arr);
                        }

                    }

                }

            }
        });
        Gpu::synchronize();
        num_deactivated_particles += *(deactivated_particles.copyToHost());
    }); // end forEachParticleTile

    ParallelDescriptor::ReduceLongSum( &num_deactivated_particles, 1 );
    Print() << "SuperDropletPC(" << m_name << "): "
            << "deactivated " << num_deactivated_particles << " super-droplets at boundaries.\n";
}

#endif
