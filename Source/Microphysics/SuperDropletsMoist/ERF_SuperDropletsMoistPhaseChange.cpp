#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Compute phase change for a timestep
 *
 * Shrink/grow the super-droplet particles for a timestep, depending on the ambient
 * flow conditions (saturation ratio, saturation pressure, and temperature).
 * This function handles all species defined in the model and updates the
 * Eulerian vapor and condensate mixing ratios accordingly.
 *
 * For each species, this function:
 * one Computes saturation pressure based on temperature
 * two Computes saturation ratio
 * three Computes phase change for all particles
 * Real(4.) Updates condensate mixing ratio
 * Real(5.) Updates vapor mixing ratio if requested
 * Real(6.) Updates temperature and pressure fields to account for latent heat
 *
 * \param[in] a_dt Timestep size
 * \param[in] a_z Array containing terrain height information
 * \param[in] a_update_qv Flag to enable updating vapor mixing ratio
 */
void SuperDropletsMoist::phaseChange ( const Real& a_dt,
                                       const Vector<std::unique_ptr<MultiFab>>& a_z,
                                       const bool a_update_qv )
{
    BL_PROFILE("SuperDropletsMoist::phaseChange()");
    for (int is = 0; is < m_num_species; is++) {
        auto& species = m_species[is];
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(species);
        bool is_water = vapour_mat.m_is_water;
        const auto idx_w = m_idx_w;

        MultiFab* qv_ptr = nullptr;
        MultiFab* qc_ptr = nullptr;
        MultiFab* qt_ptr = nullptr;
        MultiFab* sr_ptr = nullptr;
        if (is == m_idx_w) {
            qv_ptr = m_mic_fab_vars[MicVar_SD::q_v].get();
            qc_ptr = m_mic_fab_vars[MicVar_SD::q_c].get();
            qt_ptr = m_mic_fab_vars[MicVar_SD::q_t].get();
            sr_ptr = m_mic_fab_vars[MicVar_SD::rh].get();
        } else {
            qv_ptr = m_mic_fab_vars[s_qv_idx(is)].get();
            qc_ptr = m_mic_fab_vars[s_qc_idx(is)].get();
            qt_ptr = m_mic_fab_vars[s_qt_idx(is)].get();
            sr_ptr = m_mic_fab_vars[s_sr_idx(is)].get();
        }

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
            vapour_mat.computeSaturationVapFrac(    (*sr_ptr),
                                                    (*m_mic_fab_vars[MicVar_SD::temperature]),
                                                    (*m_mic_fab_vars[MicVar_SD::pressure]) );

            for (   MFIter mfi((*sr_ptr),
                    TilingIfNotGPU()); mfi.isValid();
                    ++mfi ) {

                Box bx = mfi.tilebox();
                bx.grow( sr_ptr->nGrowVect() );

                const Array4<Real>& sr_arr = sr_ptr->array(mfi);
                const Array4<Real const>& qv_arr = qv_ptr->const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });

            }

            (*sr_ptr).FillBoundary();

            // Compute total water content qt
            computeQt(is);

            // Compute super-droplets mass change
            m_super_droplets->MassChange (  0,
                                            dt_s,
                                            m_species[is],
                                            (*m_mic_fab_vars[MicVar_SD::temperature]),
                                            (*m_mic_fab_vars[MicVar_SD::pressure]),
                                            mf_sat_pressure,
                                            (*sr_ptr),
                                            a_z,
                                            is_water );

            // Compute new condensate mixing ratio
            computeQc(is, *a_z[0]);

            if (a_update_qv) {

                // Update vapour mixing ratio
                for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_v]); mfi.isValid(); ++mfi) {

                    Box bx = mfi.tilebox();
                    bx.grow(m_mic_fab_vars[MicVar_SD::q_v]->nGrowVect());

                    auto qt_arr = qt_ptr->const_array(mfi);
                    auto qv_arr = qv_ptr->array(mfi);
                    auto qc_arr = qc_ptr->array(mfi);

                    auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);
                    auto T_arr = m_mic_fab_vars[MicVar_SD::temperature]->array(mfi);
                    auto dqc_arr = m_mic_fab_vars[MicVar_SD::dqcdt]->array(mfi);
                    auto qr_arr = m_mic_fab_vars[MicVar_SD::q_r]->array(mfi);

                    auto fac_cond = vapour_mat.m_lat_vap / m_Cp;

                    ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        auto old_qv = qv_arr(i,j,k);
                        if (is == idx_w) {
                            auto qw = qc_arr(i,j,k) + qr_arr(i,j,k);
                            if (qw > qt_arr(i,j,k)) {
                                qv_arr(i,j,k) = zero;
                                if (qr_arr(i,j,k) > qt_arr(i,j,k)) {
                                    qc_arr(i,j,k) = zero;
                                    qr_arr(i,j,k) = qt_arr(i,j,k);
                                } else {
                                    qc_arr(i,j,k) = qt_arr(i,j,k) - qr_arr(i,j,k);
                                }
                            } else {
                                qv_arr(i,j,k) = qt_arr(i,j,k) - qw;
                            }
                            AMREX_ALWAYS_ASSERT(qr_arr(i,j,k) >= zero);
                        } else {
                            if (qc_arr(i,j,k) > qt_arr(i,j,k)) {
                                qv_arr(i,j,k) = zero;
                                qc_arr(i,j,k) = qt_arr(i,j,k);
                            } else {
                                qv_arr(i,j,k) = qt_arr(i,j,k) - qc_arr(i,j,k);
                            }
                        }
                        AMREX_ALWAYS_ASSERT(qv_arr(i,j,k) >= zero);
                        AMREX_ALWAYS_ASSERT(qc_arr(i,j,k) >= zero);

                        if (is == idx_w) { dqc_arr(i,j,k) = - (qv_arr(i,j,k) - old_qv) / dt_s; }

                        auto theta_over_T = theta_arr(i,j,k)/T_arr(i,j,k);
                        theta_arr(i,j,k) += theta_over_T * fac_cond * (old_qv-qv_arr(i,j,k));
                    });

                }

                // Update pressure and temperature
                const auto& gvec = (*m_mic_fab_vars[MicVar_SD::temperature]).nGrowVect();
                auto* qvw_ptr = m_mic_fab_vars[MicVar_SD::q_v].get(); // qv for water
                for (MFIter mfi(*m_mic_fab_vars[MicVar_SD::temperature], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                    Box bx = mfi.tilebox();
                    bx.grow(gvec);

                    const Array4<Real>& t_arr = m_mic_fab_vars[MicVar_SD::temperature]->array(mfi);
                    const Array4<Real>& p_arr = m_mic_fab_vars[MicVar_SD::pressure]->array(mfi);

                    const Array4<Real const>& rho_arr = m_mic_fab_vars[MicVar_SD::rho]->const_array(mfi);
                    const Array4<Real const>& theta_arr = m_mic_fab_vars[MicVar_SD::theta]->const_array(mfi);
                    const Array4<Real const>& qv_arr = qvw_ptr->const_array(mfi);

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
                vapour_mat.computeSaturationVapFrac(    (*sr_ptr),
                                                        (*m_mic_fab_vars[MicVar_SD::temperature]),
                                                        (*m_mic_fab_vars[MicVar_SD::pressure]) );

                for (   MFIter mfi((*m_mic_fab_vars[MicVar_SD::rh]),
                        TilingIfNotGPU()); mfi.isValid();
                        ++mfi ) {

                    Box bx = mfi.tilebox();
                    bx.grow( m_mic_fab_vars[MicVar_SD::rh]->nGrowVect() );

                    const Array4<Real>& sr_arr = sr_ptr->array(mfi);
                    const Array4<Real const>& qv_arr = qv_ptr->const_array(mfi);

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });

                }

                (*sr_ptr).FillBoundary();

            }
        }
    }
}

#endif
