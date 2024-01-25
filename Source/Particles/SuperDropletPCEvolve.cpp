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

    if (m_nucleate_particles) {
        MultiFab* cons_vars( &a_flow_vars[a_lev][Vars::cons] );
        // TODO: nucleation( const_vars );
    }

    {
        MultiFab* cons_vars( &a_flow_vars[a_lev][Vars::cons] );
        // TODO: massChange( const_vars );
    }

    if (m_advect_w_flow) {
        MultiFab* flow_vel( &a_flow_vars[a_lev][Vars::xvel] );
        AdvectWithFlow( flow_vel, a_lev, a_dt_lev, a_z_phys_nd[a_lev] );
    }

    if (m_advect_w_gravity) {
        AdvectWithGravity( a_lev, a_dt_lev, a_z_phys_nd[a_lev] );
    }

    // TODO: coalescence();

    return;
}

#endif
