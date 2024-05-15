#include "SuperDropletsMoist.H"
#include "EOS.H"

#ifdef ERF_USE_PARTICLES

/*! Compute phase change for a timestep: shrink/grow the super-droplet particles
    for a timestep, depending on the ambient flow conditions (saturation ratio,
    saturation pressure, and temperature). Update the Eulerial vapour and condensate
    mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange ( const Real& a_dt, /*!< Timestep */
                                       const Vector<std::unique_ptr<MultiFab>>& a_z, /*!< terrain */
                                       const bool a_update_qv,
                                       const int a_max_particle_substeps )
{
    // Get vapour material properties object
    auto& vapour_mat = m_super_droplets->getVapourMaterial();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[MicVar_SD::pressure]->nGrowVect() );
    vapour_mat.computeSaturationPressure( mf_sat_pressure,
                                          (*m_mic_fab_vars[MicVar_SD::temperature]) );
    mf_sat_pressure.FillBoundary();

    for (int substep = 0; substep < m_num_substeps_phase_change; substep++) {

        auto dt_s = a_dt / static_cast<Real>(m_num_substeps_phase_change);

        // Compute saturation ratio
        vapour_mat.computeSaturationVapFrac(    (*m_mic_fab_vars[MicVar_SD::rh]),
                                                (*m_mic_fab_vars[MicVar_SD::temperature]),
                                                (*m_mic_fab_vars[MicVar_SD::pressure]) );

        for (   MFIter mfi((*m_mic_fab_vars[MicVar_SD::rh]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[MicVar_SD::rh]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[MicVar_SD::rh]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });

        }

        (*m_mic_fab_vars[MicVar_SD::rh]).FillBoundary();

        // Compute total water content qt
        computeQt();

        // Compute super-droplets mass change
        m_super_droplets->MassChange (  0,
                                        dt_s,
                                        (*m_mic_fab_vars[MicVar_SD::temperature]),
                                        mf_sat_pressure,
                                        (*m_mic_fab_vars[MicVar_SD::rh]),
                                        a_z,
                                        a_max_particle_substeps );

        if (a_update_qv) {
            // Compute new condensate mixing ratio
            computeQc();
            // Update vapour mixing ratio
            for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_v]); mfi.isValid(); ++mfi) {

                Box bx = mfi.tilebox();
                bx.grow(m_mic_fab_vars[MicVar_SD::q_v]->nGrowVect());

                auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
                auto dqc_arr = m_mic_fab_vars[MicVar_SD::dqcdt]->array(mfi);
                auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->const_array(mfi);
                auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);

                ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    auto old_qv = q_v_arr(i,j,k);
                    q_v_arr(i,j,k) = std::max(0.0, q_t_arr(i,j,k) - q_c_arr(i,j,k));
                    dqc_arr(i,j,k) = - (q_v_arr(i,j,k) - old_qv) / dt_s;
                });

            }
        }
    }

}

#endif
