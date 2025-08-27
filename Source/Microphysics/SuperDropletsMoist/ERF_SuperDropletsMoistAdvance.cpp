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
                   << amrex::Real(num_SD_inactive)/amrex::Real(num_SD)*100
                   << "%).\n";

    // update dt
    m_dt = a_dt;

    // Compute mass/size change due to evaporation/condensation
    if (m_flag_phase_change) {
        phaseChange ( a_dt, a_z );
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
                                            a_bc );
    }

    // Coalescence of super-droplets
    if (m_flag_coalescence) {
        // Compute moist density
        MultiFab mf_moist_density(  m_mic_fab_vars[MicVar_SD::rho]->boxArray(),
                                    m_mic_fab_vars[MicVar_SD::rho]->DistributionMap(),
                                    1,
                                    m_mic_fab_vars[MicVar_SD::rho]->nGrowVect() );
        for ( MFIter mfi(mf_moist_density); mfi.isValid(); ++mfi) {

            Box bx = mfi.tilebox();
            bx.grow(mf_moist_density.nGrowVect());

            auto qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
            auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->const_array(mfi);
            auto rhom_arr = mf_moist_density.array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         { rhom_arr(i,j,k) = (1.0+qv_arr(i,j,k)) * rho_arr(i,j,k); } );
        }
        mf_moist_density.FillBoundary();

        m_super_droplets->Coalescence(  0,
                                        a_dt,
                                        *m_mic_fab_vars[MicVar_SD::pressure],
                                        mf_moist_density,
                                        *m_mic_fab_vars[MicVar_SD::temperature],
                                        *m_mic_fab_vars[MicVar_SD::q_v],
                                        a_z );
    }

    // Recycle super-droplets
    if (m_recycle_particles) {
        m_super_droplets->Recycle( 0, a_z );
    }

    m_super_droplets->Diagnostics(a_iter, a_time, (((a_iter+1)%m_diagnostics_iter==0) && (m_diagnostics_iter>0)));
}

#endif
