#include "SuperDropletsMoist.H"
#include "EOS.H"

#ifdef ERF_USE_PARTICLES

/*! Compute phase change for a timestep: shrink/grow the super-droplet particles
    for a timestep, depending on the ambient flow conditions (saturation ratio,
    saturation pressure, and temperature). Update the Eulerial vapour and condensate
    mixing ratios accordingly. */
void SuperDropletsMoist::phaseChange ( const Real& a_dt, /*!< Timestep */
                                       const Vector<std::unique_ptr<MultiFab>>& a_z /*!< terrain */)
{
    // Get vapour material properties object
    auto& vapour_mat = m_super_droplets->getVapourMaterial();

    // Compute saturation pressure
    MultiFab mf_sat_pressure(   m_mic_fab_vars[MicVar_SD::pressure]->boxArray(),
                                m_mic_fab_vars[MicVar_SD::pressure]->DistributionMap(),
                                1,
                                m_mic_fab_vars[MicVar_SD::pressure]->nGrowVect() );
    vapour_mat.computeSaturationPressure( mf_sat_pressure, (*m_mic_fab_vars[MicVar_SD::temperature]) );

    // Compute saturation ratio
    MultiFab mf_sat_ratio(  m_mic_fab_vars[MicVar_SD::pressure]->boxArray(),
                            m_mic_fab_vars[MicVar_SD::pressure]->DistributionMap(),
                            1,
                            m_mic_fab_vars[MicVar_SD::pressure]->nGrowVect() );
    vapour_mat.computeSaturationVapFrac(    mf_sat_ratio,
                                            (*m_mic_fab_vars[MicVar_SD::temperature]),
                                            (*m_mic_fab_vars[MicVar_SD::pressure]) );
    for (MFIter mfi(mf_sat_ratio, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& sr_arr = mf_sat_ratio.array(mfi);
        const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });
    }

    // Compute total water content qt
    computeQt();

    // Compute super-droplets mass change
    m_super_droplets->MassChange (  0,
                                    a_dt,
                                    (*m_mic_fab_vars[MicVar_SD::temperature]),
                                    mf_sat_pressure,
                                    mf_sat_ratio,
                                    a_z );

    // Compute new condensate mixing ratio
    computeQc();

    // Update vapour mixing ratio
    for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_v]); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
        auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->const_array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { q_v_arr(i,j,k) = std::max(0.0, q_t_arr(i,j,k) - q_c_arr(i,j,k)); });
    }

}

#endif
