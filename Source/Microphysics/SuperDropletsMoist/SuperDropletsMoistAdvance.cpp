#include "SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

/*! Advance the moisture model for a timestep: evolve the super-droplet particles
    for a timestep - this includes nucleation, advection, and coalescence */
void SuperDropletsMoist::Advance ( const Real& a_dt, /*!< Timestep */
                                   const int& a_iter, /*!< Iteration number */
                                   const Real& /* a_time */, /*!< Simulation time */
                                   Vector<Vector<MultiFab>>& a_flow_vars, /*!< flow variables (*all*) */
                                   const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */)
{
    amrex::Print() << "SuperDropletsMoist: iteration=" << a_iter+1
                   << ", dt=" << a_dt <<", evolving "
                   << m_super_droplets->NumSuperDroplets()
                   << " super-droplets representing "
                   << m_super_droplets->TotalNumberOfParticles()
                   << " particles.\n";

    // Compute mass/size change due to evaporation/condensation
    if (m_flag_phase_change) {
        phaseChange ( a_dt, a_z, true );
    }

    // Advect particles
    if (m_flag_advection) {
        m_super_droplets->AdvectParticles (0, a_dt, a_flow_vars, a_z);
    }

    // Coalescence of super-droplets
    if (m_flag_coalescence) {
        m_super_droplets->Coalescence(0,a_dt,*m_mic_fab_vars[MicVar_SD::temperature]);
    }

    m_super_droplets->Diagnostics(a_iter, ((a_iter+1)%m_diagnostics_iter == 0));
}

#endif
