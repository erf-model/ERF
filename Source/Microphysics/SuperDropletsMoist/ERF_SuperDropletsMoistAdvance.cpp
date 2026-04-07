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

    // update dt
    m_dt = a_dt;

    // inject particles
    m_super_droplets->InjectParticles(a_time, a_z[0], m_dt);

    auto num_particles = m_super_droplets->TotalNumberOfParticles();
    auto num_SD = m_super_droplets->NumSuperDroplets();
    auto num_SD_inactive = m_super_droplets->NumSDDeactivated();

    amrex::Print() << "SuperDropletsMoist: iteration=" << a_iter+1
                   << ", dt=" << a_dt <<", evolving "
                   << num_SD
                   << " super-droplets representing "
                   << num_particles
                   << " particles.\n";

    amrex::Print() << "    Number of deactivated super-droplets: "
                   << num_SD_inactive
                   << " ("
                   << (num_SD > 0 ? Real(num_SD_inactive)/Real(num_SD)*100 : 0)
                   << "%).\n";

    // Compute mass/size change due to evaporation/condensation
    if (m_flag_phase_change) {
        phaseChange ( a_dt, a_z, true );
    }

    // Advect particles
    if (m_flag_advection) {
        m_super_droplets->AdvectParticles ( 0,
                                            a_time,
                                            a_dt,
                                            &a_flow_vars[0][Vars::xvel],
                                            *(m_mic_fab_vars[MicVar_SD::rho]),
                                            *(m_mic_fab_vars[MicVar_SD::pressure]),
                                            *(m_mic_fab_vars[MicVar_SD::temperature]),
                                            a_z,
                                            a_bc,
                                            m_recycle_particles );
    }

    // Coalescence of super-droplets
    if (m_flag_coalescence) {
        m_super_droplets->Coalescence(  0,
                                        a_dt,
                                        *m_mic_fab_vars[MicVar_SD::pressure],
                                        *m_mic_fab_vars[MicVar_SD::temperature] );
    }

    // Recycle super-droplets
    m_super_droplets->Recycle( 0, a_z, a_iter, a_dt, m_recycle_particles );

    m_super_droplets->Diagnostics(a_iter, a_time, (((a_iter+1)%m_diagnostics_iter==0) && (m_diagnostics_iter>0)));
}

#endif
