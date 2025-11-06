#include "ERF_IndexDefines.H"
#include "ERF_SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

/*! Advance the moisture model for a timestep: evolve the super-droplet particles
    for a timestep - this includes nucleation, advection, and coalescence */
void SuperDropletsMoist::Advance ( const Real& a_dt, /*!< Timestep */
                                   const int& a_iter, /*!< Iteration number */
                                   const Real&  a_time, /*!< Simulation time */
                                   Vector<Vector<MultiFab>>& a_flow_vars, /*!< flow variables (*all*) */
                                   const Vector<MFPtr>& a_z, /*!< terrain */
                                   const BCTypeArr& a_bc /*! Boundary types */)
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
                   << (num_SD > 0 ? amrex::Real(num_SD_inactive)/amrex::Real(num_SD)*100 : 0)
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
