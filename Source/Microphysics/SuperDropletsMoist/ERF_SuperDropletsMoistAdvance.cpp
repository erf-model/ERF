#include "ERF_IndexDefines.H"
#include "ERF_SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Advance the moisture model for a timestep
 *
 * Evolve the super-droplet particles for a timestep - this includes:
 * one Injection of new particles if configured
 * two Phase change (condensation/evaporation) if enabled
 * three Advection if enabled
 * Real(4.) Coalescence if enabled
 * Real(5.) Recycling of particles
 * Real(6.) Computing diagnostics at specified intervals
 *
 * \param[in] a_dt Timestep size
 * \param[in] a_iter Current iteration number
 * \param[in] a_time Current simulation time
 * \param[in,out] a_flow_vars Flow variables for all components
 * \param[in] a_z Terrain height information
 * \param[in] a_bc Boundary condition types
 */
void SuperDropletsMoist::Advance ( const Real& a_dt,
                                   const int& a_iter,
                                   const Real&  a_time,
                                   Vector<Vector<MultiFab>>& a_flow_vars,
                                   const Vector<MFPtr>& a_z,
                                   const BCTypeArr& a_bc )
{
    BL_PROFILE("SuperDropletsMoist::Advance()");
    AMREX_ALWAYS_ASSERT((m_current_lev >= 0) && (m_current_lev < static_cast<int>(m_mic_fab_vars.size())));

    // update dt
    m_dt = a_dt;

    // Current level (set by wrapper)
    const int current_lev = m_current_lev;
    const int lev = m_current_lev;

    // Inject particles on level 0 only
    if (current_lev == 0) {
        m_super_droplets->InjectParticles(a_time, a_z[0], m_dt);
    }

    // Print diagnostics once per iteration (level 0)
    if (current_lev == 0) {
        auto num_particles = m_super_droplets->TotalNumberOfParticles();
        auto num_SD = m_super_droplets->NumSuperDroplets();
        auto num_SD_inactive = m_super_droplets->NumSDDeactivated();

        amrex::Print() << "SuperDropletsMoist: iteration=" << a_iter+1
                       << ", dt=" << a_dt <<", evolving "
                       << num_SD
                       << " super-droplets representing "
                       << num_particles
                       << " particles.\n";

        amrex::Print() << "    Per-level super-droplet counts:";
        for (int l = 0; l <= m_super_droplets->finestLevel(); l++) {
            amrex::Print() << " level " << l << ": "
                           << m_super_droplets->NumberOfParticlesAtLevel(l);
        }
        amrex::Print() << "\n";

        amrex::Print() << "    Number of deactivated super-droplets: "
                       << num_SD_inactive
                       << " ("
                       << (num_SD > 0 ? amrex::Real(num_SD_inactive)/amrex::Real(num_SD)*100 : 0)
                       << "%).\n";
    }

    // Verify m_mic_fab_vars matches current level geometry
    int local_check = (m_mic_fab_vars[lev].empty() ||
                       m_mic_fab_vars[lev][MicVar_SD::rho]->boxArray() != a_flow_vars[current_lev][Vars::cons].boxArray()) ? 1 : 0;
    int global_check = local_check;
    ParallelDescriptor::ReduceIntMax(global_check);
    if (global_check > 0) {
        if (local_check > 0) {
            amrex::Warning("SuperDropletsMoist::Advance: m_mic_fab_vars doesn't match level geometry, skipping");
        }
        return;
    }

    // Compute mass/size change due to evaporation/condensation
    if (m_flag_phase_change) {
        phaseChange ( a_dt, a_z, current_lev );
    }

    // Advect particles
    if (m_flag_advection) {
        m_super_droplets->AdvectParticles ( current_lev,
                                            a_time,
                                            a_dt,
                                            &a_flow_vars[current_lev][Vars::xvel],
                                            *(m_mic_fab_vars[lev][MicVar_SD::rho]),
                                            *(m_mic_fab_vars[lev][MicVar_SD::pressure]),
                                            *(m_mic_fab_vars[lev][MicVar_SD::temperature]),
                                            a_z,
                                            a_bc,
                                            m_recycle_particles );
    }

    // Coalescence of super-droplets
    if (m_flag_coalescence) {
        m_super_droplets->Coalescence(  current_lev,
                                        a_dt,
                                        *m_mic_fab_vars[lev][MicVar_SD::pressure],
                                        *m_mic_fab_vars[lev][MicVar_SD::temperature] );
    }

    // Recycle super-droplets
    m_super_droplets->Recycle( current_lev, a_z, a_iter, a_dt, m_recycle_particles );

    // Diagnostics on level 0 only (avoids double-counting)
    if (current_lev == 0) {
        m_super_droplets->Diagnostics(a_iter, current_lev, a_time, (((a_iter+1)%m_diagnostics_iter==0) && (m_diagnostics_iter>0)));
    }
}

#endif
