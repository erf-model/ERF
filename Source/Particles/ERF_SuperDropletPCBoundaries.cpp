#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

#include <AMReX_TracerParticle_mod_K.H>
#include "ERF_InterpolationUtils.H"

using namespace amrex;
using namespace SDPCDefn;


/*! Handle the boundaries for the particles, then redistribute them */
void SuperDropletPC::applyBoundaryTreatment ( int                   a_lev,
                                              const Vector<MFPtr>&  a_z_phys_nd,
                                              const BCTypeArr&      a_bctypes,
                                              const bool            a_recycle )
{
    BL_PROFILE("SuperDropletPC::applyBoundaryTreatment()");
    const Long num_inflow = treatBoundaryCrossings(a_lev, a_bctypes, a_recycle);
    Redistribute();
    // super-droplets that entered through an open boundary, now on the grid of the cell they entered
    if (num_inflow > 0) { setInflowMultiplicity(a_z_phys_nd); }
}

/*! Deactivate, reflect or wrap the particles that crossed a domain boundary; returns the
 *  number that entered through an open boundary with their multiplicity still to be set
 *  (stretched or terrain-fitted mesh) */
Long SuperDropletPC::treatBoundaryCrossings ( int              a_lev,
                                              const BCTypeArr& a_bctypes,
                                              const bool       a_recycle )
{
    BL_PROFILE("SuperDropletPC::treatBoundaryCrossings()");

    const auto proc_ctx = buildProcessContext(a_lev);
    const auto save_inac = m_save_inactive;
    const bool variable_dz = m_variable_dz;

    // On a stretched or terrain-fitted mesh the multiplicity follows the volume of the
    // cell entered; a negative value marks it for setInflowMultiplicity after Redistribute
    Real multiplicity = -one;
    if (!variable_dz) {
        // number of super-droplets per cell
        int num_sd_per_cell = m_num_sd_per_cell;
        // number of physical particles per cell
        Real num_par_per_cell = zero;
        for (int i = 0; i < m_num_initializations; i++) {
            num_par_per_cell += m_initializations[i]->numParticlesPerCell(proc_ctx.cell_volume);
        }
        multiplicity = (num_sd_per_cell > 0 ? num_par_per_cell / num_sd_per_cell : zero);
    }

    Long num_deactivated_particles = 0;
    Long num_inflow_particles = 0;

    forEachParticleTile(a_lev, proc_ctx,
        [&](ParIterType& /*pti*/, int /*grid*/, ParticleType* p_pbox,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& ctx)
    {
        Gpu::Buffer<Long> deactivated_particles({0});
        auto* deactivated_particles_ptr = deactivated_particles.data();
        Gpu::Buffer<Long> inflow_particles({0});
        auto* inflow_particles_ptr = variable_dz ? inflow_particles.data() : nullptr;

        ParallelFor(ptrs.num_particles, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }
            if (ptrs.active_ptr[i] == 0) { return; }

            // check for ground impact; in computational zeta the ground is the
            // bottom of the domain for both flat and terrain-fitted meshes
            {
                const auto z_ground = ctx.plo[2];
                if (p.pos(2) < z_ground) {
                    p.pos(2) = z_ground + Real(0.01)*ctx.dx[2];
                    ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;
                    ptrs.active_ptr[i] = 0;
                    if ((!a_recycle) && (!save_inac)) { p.id() = -1; }
                    Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));
                }
            }

            // check for top boundary exits; in computational zeta the roof is the
            // top of the domain for both flat and terrain-fitted meshes
            {
                const auto z_roof = ctx.phi[2];
                if (p.pos(2) > z_roof) {
                    p.pos(2) = z_roof - ctx.dx[2];
                    ptrs.v_ptr[0][i] = ptrs.v_ptr[1][i] = ptrs.v_ptr[2][i] = ptrs.vterm_ptr[i] = zero;
                    ptrs.active_ptr[i] = 0;
                    if ((!a_recycle) && (!save_inac)) { p.id() = -1; }
                    Gpu::Atomic::Add(deactivated_particles_ptr, Long(1));
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

                            ptrs.mult_ptr[i]  = multiplicity;
                            if (inflow_particles_ptr) { Gpu::Atomic::Add(inflow_particles_ptr, Long(1)); }
                            ptrs.a_ptr[i]     = zero;
                            ptrs.c_ptr[i]     = zero;
                            ptrs.mrime_ptr[i] = zero;
                            ptrs.nmono_ptr[i] = zero;

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

                            ptrs.mult_ptr[i]  = multiplicity;
                            if (inflow_particles_ptr) { Gpu::Atomic::Add(inflow_particles_ptr, Long(1)); }
                            ptrs.a_ptr[i]     = zero;
                            ptrs.c_ptr[i]     = zero;
                            ptrs.mrime_ptr[i] = zero;
                            ptrs.nmono_ptr[i] = zero;

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
#ifdef AMREX_USE_OMP
#pragma omp atomic
#endif
        num_deactivated_particles += *(deactivated_particles.copyToHost());
        if (variable_dz) {
#ifdef AMREX_USE_OMP
#pragma omp atomic
#endif
            num_inflow_particles += *(inflow_particles.copyToHost());
        }
    }); // end forEachParticleTile

    ParallelDescriptor::ReduceLongSum( &num_deactivated_particles, 1 );
    if (variable_dz) { ParallelDescriptor::ReduceLongSum( &num_inflow_particles, 1 ); }
    Print() << "SuperDropletPC(" << m_name << "): "
            << "deactivated " << num_deactivated_particles << " super-droplets at boundaries.\n";
    return num_inflow_particles;
}

/*! Set the multiplicity of super-droplets that entered through an open lateral boundary
 *  (marked by a negative multiplicity) to the initial multiplicity of the cell they are in */
void SuperDropletPC::setInflowMultiplicity ( const Vector<MFPtr>& a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPC::setInflowMultiplicity()");

    const auto init_counts_d = initCellCounts();
    const SDCellCounts* init_counts = init_counts_d.data();
    const int num_init = static_cast<int>(init_counts_d.size());

    for (int lev = 0; lev <= finestLevel(); lev++) {
        const auto& geom = Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto dx = geom.CellSizeArray();
        const auto domain = geom.Domain();
        const MFPtr& z_height = a_z_phys_nd[lev];

        forEachParticleTile(lev,
            [&](ParIterType& /*pti*/, int grid, ParticleType* p_pbox,
                const SDProcess::ParticlePointers& ptrs, int np)
        {
            const SDCellVolume cell_vol{dx, (*z_height)[grid].const_array(), true};
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                const ParticleType& p = p_pbox[i];
                if (p.id() <= 0) { return; }
                if (ptrs.mult_ptr[i] >= zero) { return; }
                const auto iv = getParticleCell(p, plo, dxi, domain);
                ptrs.mult_ptr[i] = SD_init_multiplicity(init_counts, num_init, cell_vol(iv[0],iv[1],iv[2]));
            });
        });
    }
}

#endif
