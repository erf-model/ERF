#include "SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

/*! Advance the moisture model for a timestep: evolve the super-droplet particles
    for a timestep - this includes nucleation, advection, and coalescence */
void SuperDropletsMoist::Advance ( const Real& a_dt, /*!< Timestep */
                                   Vector<Vector<MultiFab>>& a_flow_vars, /*!< flow variables (*all*) */
                                   const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */)
{
    // Compute mass/size change due to evaporation/condensation
    phaseChange ( a_dt, a_z );

    // Advect particles
    m_super_droplets->AdvectParticles (0, a_dt, a_flow_vars, a_z);
}

#endif
