#include <ERF.H>

using namespace amrex;

/**
 * Function that advances the solution at one level for a single time step --
 * this sets up the multirate time integrator and calls the integrator's advance function
 *
 * @param[in] lev level of refinement (coarsest level is 0)
 * @param[in] dt_lev time step for this time advance
 */

void ERF::advance_particles(int lev, Real dt_lev)
{
    // Update tracer particles on level 0
    if (lev == 0 && use_tracer_particles) {
        MultiFab* mac_vel = &vars_new[lev][Vars::xvel];
        tracer_particles->AdvectWithUmac(mac_vel, lev, dt_lev, *z_phys_nd[lev]);
    }
}
