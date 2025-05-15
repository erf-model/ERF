#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"

#ifdef ERF_USE_PARTICLES

/*! Compute phase change for a timestep: shrink/grow the super-droplet particles
    for a timestep, depending on the ambient flow conditions (saturation ratio,
    saturation pressure, and temperature). Update the Eulerial vapour and condensate
    mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange ( const Real& a_dt, /*!< Timestep */
                                       const Vector<std::unique_ptr<MultiFab>>& a_z, /*!< terrain */
                                       const bool a_update_qv )
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
                                        (*m_mic_fab_vars[MicVar_SD::pressure]),
                                        mf_sat_pressure,
                                        (*m_mic_fab_vars[MicVar_SD::rh]),
                                        a_z );

        // Compute new condensate mixing ratio
        computeQcQr();

        if (a_update_qv) {

            // Update vapour mixing ratio
            for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_v]); mfi.isValid(); ++mfi) {

                Box bx = mfi.tilebox();
                bx.grow(m_mic_fab_vars[MicVar_SD::q_v]->nGrowVect());

                auto qt_arr = m_mic_fab_vars[MicVar_SD::q_t]->const_array(mfi);
                auto qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
                auto dqc_arr = m_mic_fab_vars[MicVar_SD::dqcdt]->array(mfi);
                auto qc_arr = m_mic_fab_vars[MicVar_SD::q_c]->array(mfi);
                auto qr_arr = m_mic_fab_vars[MicVar_SD::q_r]->array(mfi);
                auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);
                auto T_arr = m_mic_fab_vars[MicVar_SD::temperature]->array(mfi);

                auto fac_cond = m_fac_cond;
                ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    AMREX_ALWAYS_ASSERT(qr_arr(i,j,k) >= 0.0);
                    auto old_qv = qv_arr(i,j,k);
                    auto qw = qc_arr(i,j,k) + qr_arr(i,j,k);
                    if (qw > qt_arr(i,j,k)) {
                        qv_arr(i,j,k) = 0.0;
                        if (qr_arr(i,j,k) > qt_arr(i,j,k)) {
                            qc_arr(i,j,k) = 0.0;
                        } else {
                            qc_arr(i,j,k) = qt_arr(i,j,k) - qr_arr(i,j,k);
                        }
                    } else {
                        qv_arr(i,j,k) = qt_arr(i,j,k) - (qc_arr(i,j,k) + qr_arr(i,j,k));
                    }
                    AMREX_ALWAYS_ASSERT(qv_arr(i,j,k) >= 0.0);
                    AMREX_ALWAYS_ASSERT(qc_arr(i,j,k) >= 0.0);
                    dqc_arr(i,j,k) = - (qv_arr(i,j,k) - old_qv) / dt_s;

                    auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
                    theta_arr(i,j,k) += theta_over_T * fac_cond * (old_qv-qv_arr(i,j,k));
                });

            }

            const auto& gvec = (*m_mic_fab_vars[MicVar_SD::temperature]).nGrowVect();
            // Update pressure and temperature
            for (MFIter mfi(*m_mic_fab_vars[MicVar_SD::temperature], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                Box bx = mfi.tilebox();
                bx.grow(gvec);

                const Array4<Real>& t_arr = m_mic_fab_vars[MicVar_SD::temperature]->array(mfi);
                const Array4<Real>& p_arr = m_mic_fab_vars[MicVar_SD::pressure]->array(mfi);

                const Array4<Real const>& rho_arr = m_mic_fab_vars[MicVar_SD::rho]->const_array(mfi);
                const Array4<Real const>& theta_arr = m_mic_fab_vars[MicVar_SD::theta]->const_array(mfi);
                const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    t_arr(i,j,k,0) = getTgivenRandRTh(  rho_arr(i,j,k,0),
                                                        rho_arr(i,j,k,0)*theta_arr(i,j,k,0),
                                                        qv_arr(i,j,k,0));
                    p_arr(i,j,k,0) = getPgivenRTh(  rho_arr(i,j,k,0)*theta_arr(i,j,k,0),
                                                    qv_arr(i,j,k,0));
                });
            }

            // Update saturation ratio
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

        }
    }

}

#endif
