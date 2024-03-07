#include <SuperDropletPC.H>

#ifdef ERF_USE_PARTICLES

#include <IndexDefines.H>
#include <ERF_Constants.H>

using namespace amrex;

/*! Evolve particles for one time step */
void SuperDropletPC::EvolveParticles ( int                                        a_lev,
                                       Real                                       a_dt_lev,
                                       Vector<Vector<MultiFab>>&                  a_flow_vars,
                                       const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPCPC::EvolveParticles()");
    AMREX_ASSERT( a_lev == m_lev );

    AdvectParticles( a_lev, a_dt_lev, a_flow_vars, a_z_phys_nd, m_advect_w_flow, m_advect_w_gravity );

    // TODO: coalescence();

    Redistribute();

    return;
}

/*! Evolve particles for one time step */
void SuperDropletPC::EvolveParticles ( int                                        a_lev,
                                       Real                                       a_dt_lev,
                                       Vector<Vector<MultiFab>>&                  a_flow_vars,
                                       MultiFab&                                  a_qv,
                                       const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    BL_PROFILE("SuperDropletPCPC::EvolveParticles()");
    AMREX_ASSERT( a_lev == m_lev );

    if (m_nucleate_particles) {
        //MultiFab* cons_vars( &a_flow_vars[a_lev][Vars::cons] );
        // TODO: nucleation( const_vars );
    }

    MassChange( a_lev, a_dt_lev, a_flow_vars, a_qv, a_z_phys_nd );
    AdvectParticles( a_lev, a_dt_lev, a_flow_vars, a_z_phys_nd, m_advect_w_flow, m_advect_w_gravity );

    // TODO: coalescence();

    Redistribute();

    return;
}

#endif
