#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"

#ifdef ERF_USE_PARTICLES

/*! Copy moisture model variables from the conserved state vector to the
    member multifabs in this object, compute and save pressure and temperature */
void SuperDropletsMoist::Copy_State_to_Micro (  const MultiFab& a_cons_vars /*!< Conserved variables */)
{
    BL_PROFILE("SuperDropletsMoist::Copy_State_to_Micro()");
    const auto& gvec = a_cons_vars.nGrowVect();

    // Copy density and vapour mixing ratio from state variables
    // Note: do *not* copy other q
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow(gvec);
        auto states_arr = a_cons_vars.const_array(mfi);

        // state variables
        {
            auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->array(mfi);
            auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rho_arr(i,j,k) = states_arr(i,j,k,Rho_comp);
                theta_arr(i,j,k) = states_arr(i,j,k,RhoTheta_comp)/states_arr(i,j,k,Rho_comp);
            });
        }

        // water/ice
        {
            auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);
            auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
            auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
            auto q_i_arr = m_mic_fab_vars[MicVar_SD::q_i]->const_array(mfi);
            auto q_s_arr = m_mic_fab_vars[MicVar_SD::q_s]->const_array(mfi);
            auto q_g_arr = m_mic_fab_vars[MicVar_SD::q_g]->const_array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                q_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp) / states_arr(i,j,k,Rho_comp);
                q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k)
                                 + q_r_arr(i,j,k)
                                 + q_i_arr(i,j,k)
                                 + q_s_arr(i,j,k)
                                 + q_g_arr(i,j,k);
            });
        }

        // other species
        for (int is = m_istart_sp; is < m_num_species; is++) {
            auto q_v_arr = m_mic_fab_vars[s_qv_idx(is,m_istart_sp)]->array(mfi);
            auto q_t_arr = m_mic_fab_vars[s_qt_idx(is,m_istart_sp)]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(is,m_istart_sp)]->const_array(mfi);
            auto qv_comp = q_qv_idx(is,m_istart_sp);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                q_v_arr(i,j,k) = states_arr(i,j,k,qv_comp) / states_arr(i,j,k,Rho_comp);
                q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k);
            });
        }
    }

    // Compute pressure and temperature
    for (MFIter mfi(*m_mic_fab_vars[MicVar_SD::temperature], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        const Array4<Real      >& t_arr  = m_mic_fab_vars[MicVar_SD::temperature]->array(mfi);
        const Array4<Real      >& p_arr  = m_mic_fab_vars[MicVar_SD::pressure]->array(mfi);
        const Array4<Real const>& S_arr  = a_cons_vars.const_array(mfi);
        const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            t_arr(i,j,k,0) = getTgivenRandRTh(S_arr(i,j,k,Rho_comp),S_arr(i,j,k,RhoTheta_comp),qv_arr(i,j,k,0));
            p_arr(i,j,k,0) = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp),qv_arr(i,j,k,0));
        });
    }

    AMREX_ASSERT( !m_mic_fab_vars[MicVar_SD::pressure]->contains_nan() );
    AMREX_ASSERT( !m_mic_fab_vars[MicVar_SD::temperature]->contains_nan() );

    // compute saturation ratio for water
    {
        // Get vapour material properties object for water
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::H2O);

        // Compute saturation ratio
        vapour_mat.computeSaturationVapFrac( (*m_mic_fab_vars[MicVar_SD::rh_w]),
                                             (*m_mic_fab_vars[MicVar_SD::temperature]),
                                             (*m_mic_fab_vars[MicVar_SD::pressure]) );

        for (   MFIter mfi((*m_mic_fab_vars[MicVar_SD::rh_w]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[MicVar_SD::rh_w]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[MicVar_SD::rh_w]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });

        }
    }

    // compute saturation ratio for ice
    if (m_with_ice) {
        // Get vapour material properties object for ice
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::ice);

        // Compute saturation ratio
        vapour_mat.computeSaturationVapFrac( (*m_mic_fab_vars[MicVar_SD::rh_i]),
                                             (*m_mic_fab_vars[MicVar_SD::temperature]),
                                             (*m_mic_fab_vars[MicVar_SD::pressure]) );

        for (   MFIter mfi((*m_mic_fab_vars[MicVar_SD::rh_i]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[MicVar_SD::rh_i]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[MicVar_SD::rh_i]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = (sr_arr(i,j,k,0) > 0.0 ? qv_arr(i,j,k,0) / sr_arr(i,j,k,0) : 0.0); });

        }
    }

    // other species
    for (int is = m_istart_sp; is < m_num_species; is++) {
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(m_species[is]);
        vapour_mat.computeSaturationVapFrac((*m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]),
                                            (*m_mic_fab_vars[MicVar_SD::temperature]),
                                            (*m_mic_fab_vars[MicVar_SD::pressure]) );
        for (   MFIter mfi((*m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[s_qv_idx(is,m_istart_sp)]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });
        }
    }

    for (auto i(0); i < MicVar_SD::NumVars; i++) {
        m_mic_fab_vars[i]->FillBoundary(m_geom.periodicity());
    }
}

/*! Copy moisture model variables to the conserved state vector from the
    member multifabs in this object */
void SuperDropletsMoist::Copy_Micro_to_State (  MultiFab& a_cons_vars /*!< Conserved variables */)
{
    BL_PROFILE("SuperDropletsMoist::Copy_Micro_to_state()");
    const auto& gvec = a_cons_vars.nGrowVect();

    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto states_arr = a_cons_vars.array(mfi);

        // state variables
        {
            auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->const_array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,RhoTheta_comp) = states_arr(i,j,k,Rho_comp)*theta_arr(i,j,k);
            });
        }

        // water
        {
            auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
            auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
            auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
            auto q_i_arr = m_mic_fab_vars[MicVar_SD::q_i]->const_array(mfi);
            auto q_s_arr = m_mic_fab_vars[MicVar_SD::q_s]->const_array(mfi);
            auto q_g_arr = m_mic_fab_vars[MicVar_SD::q_g]->const_array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,RhoQ1_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
                states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
                states_arr(i,j,k,RhoQ3_comp) = states_arr(i,j,k,Rho_comp)*q_i_arr(i,j,k);
                states_arr(i,j,k,RhoQ4_comp) = states_arr(i,j,k,Rho_comp)*q_r_arr(i,j,k);
                states_arr(i,j,k,RhoQ5_comp) = states_arr(i,j,k,Rho_comp)*q_s_arr(i,j,k);
                states_arr(i,j,k,RhoQ6_comp) = states_arr(i,j,k,Rho_comp)*q_g_arr(i,j,k);
            });
        }

        // other species
        for (int is = m_istart_sp; is < m_num_species; is++) {
            auto q_v_arr = m_mic_fab_vars[s_qv_idx(is,m_istart_sp)]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(is,m_istart_sp)]->array(mfi);
            auto qv_comp = q_qv_idx(is,m_istart_sp);
            auto qc_comp = q_qc_idx(is,m_istart_sp);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,qv_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
                states_arr(i,j,k,qc_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
            });
        }
    }

    a_cons_vars.FillBoundary(m_geom.periodicity());
}

/*! Update microphysics variables */
void SuperDropletsMoist::Update_Micro_Vars (MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Update_Micro_Vars()");
    Copy_State_to_Micro(a_cons_vars);
}

/*! Compute qc, qt, and rain accumulation,and update state variables
 *  from microphysics variables */
void SuperDropletsMoist::Update_State_Vars (MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Update_State_Vars()");
    computeQcQrWater();
    computeQiQgQsWater();
    computeQtWater();
    rainAccumulation();

    computeQcSpecies();
    computeQtSpecies();
    speciesAccumulation();

    if (!m_kinematic_mode) { Copy_Micro_to_State(a_cons_vars); }
}

/*! Convert a multifab containing density of something to its mixing ratio */
void SuperDropletsMoist::densityToRatio (  MultiFab& a_var, /*!< Multifab */
                                           const int a_comp /*!< Component */ )
{
    BL_PROFILE("SuperDropletsMoist::densityToRatio()");
    const auto& gvec = a_var.nGrowVect();

    for ( MFIter mfi(a_var); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->const_array(mfi);
        auto fab_arr = a_var.array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { fab_arr(i,j,k,a_comp) /= rho_arr(i,j,k); });
    }

    a_var.FillBoundary(m_geom.periodicity());
}

/*! Convert a multifab containing the mixing ratio of something to its density */
void SuperDropletsMoist::ratioToDensity (  MultiFab& a_var, /*!< Multifab */
                                           const int a_comp /*!< Component */ )
{
    BL_PROFILE("SuperDropletsMoist::ratioToDensity()");
    const auto& gvec = a_var.nGrowVect();

    for ( MFIter mfi(a_var); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->const_array(mfi);
        auto fab_arr = a_var.array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { fab_arr(i,j,k,a_comp) *= rho_arr(i,j,k); });
    }

    a_var.FillBoundary(m_geom.periodicity());
}

/*! compute cloud/rain mixing ratio for water */
void SuperDropletsMoist::computeQcQrWater ()
{
    BL_PROFILE("SuperDropletsMoist::computeQcQrWater()");
    m_super_droplets->cloudRainDensity( *(m_mic_fab_vars[MicVar_SD::q_c]),
                                        0,
                                        m_r_rain );
    m_super_droplets->cloudRainDensity( *(m_mic_fab_vars[MicVar_SD::q_r]),
                                        m_r_rain,
                                        1.0 );

    if (m_dimensionality == SDMSimulationDim::one_d_z) {
        for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_c]); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            int imin = bx.smallEnd(0);
            int jmin = bx.smallEnd(1);
            auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->array(mfi);
            auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                q_c_arr(i,j,k) = q_c_arr(imin,jmin,k);
                q_r_arr(i,j,k) = q_r_arr(imin,jmin,k);
            });
        }
    }

    densityToRatio(*(m_mic_fab_vars[MicVar_SD::q_c]));
    densityToRatio(*(m_mic_fab_vars[MicVar_SD::q_r]));
}

/*! compute ice/graup/snow mixing ratio for water */
void SuperDropletsMoist::computeQiQgQsWater ()
{
    BL_PROFILE("SuperDropletsMoist::computeQiQgQsWater()");
    if (m_with_ice) {
        m_super_droplets->iceSnowGraupelDensity ( *(m_mic_fab_vars[MicVar_SD::q_i]),
                                                  0.0,
                                                  m_rime_ratio );
        m_super_droplets->iceSnowGraupelDensity ( *(m_mic_fab_vars[MicVar_SD::q_g]),
                                                  m_rime_ratio,
                                                  1.0 );

        // TODO: What is snow?

        if (m_dimensionality == SDMSimulationDim::one_d_z) {
            for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_c]); mfi.isValid(); ++mfi) {
                Box bx = mfi.tilebox();
                int imin = bx.smallEnd(0);
                int jmin = bx.smallEnd(1);
                auto q_i_arr = m_mic_fab_vars[MicVar_SD::q_i]->array(mfi);
                auto q_g_arr = m_mic_fab_vars[MicVar_SD::q_g]->array(mfi);
                auto q_s_arr = m_mic_fab_vars[MicVar_SD::q_s]->array(mfi);

                ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    q_i_arr(i,j,k) = q_i_arr(imin,jmin,k);
                    q_g_arr(i,j,k) = q_g_arr(imin,jmin,k);
                    q_s_arr(i,j,k) = q_s_arr(imin,jmin,k);
                });
            }
        }

        densityToRatio(*(m_mic_fab_vars[MicVar_SD::q_i]));
        densityToRatio(*(m_mic_fab_vars[MicVar_SD::q_g]));
        densityToRatio(*(m_mic_fab_vars[MicVar_SD::q_s]));
    }
}

/*! compute qt (total) for water */
void SuperDropletsMoist::computeQtWater ()
{
    BL_PROFILE("SuperDropletsMoist::computeQtWater()");
    for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_t]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(m_mic_fab_vars[MicVar_SD::q_t]->nGrowVect());

        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
        auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
        auto q_i_arr = m_mic_fab_vars[MicVar_SD::q_i]->const_array(mfi);
        auto q_s_arr = m_mic_fab_vars[MicVar_SD::q_s]->const_array(mfi);
        auto q_g_arr = m_mic_fab_vars[MicVar_SD::q_g]->const_array(mfi);
        auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k)
                             + q_r_arr(i,j,k)
                             + q_i_arr(i,j,k)
                             + q_s_arr(i,j,k)
                             + q_g_arr(i,j,k);
        });
    }
}

/*! Compute rain accumulation */
void SuperDropletsMoist::rainAccumulation ()
{
    BL_PROFILE("SuperDropletsMoist::rainAccumulation()");
    auto domain = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::H2O);
    auto mat_density = vapour_mat.m_density;

    MultiFab mf_zflux( m_mic_fab_vars[MicVar_SD::rain_accum]->boxArray(),
                       m_mic_fab_vars[MicVar_SD::rain_accum]->DistributionMap(),
                       1,
                       m_mic_fab_vars[MicVar_SD::rain_accum]->nGrowVect() );
    m_super_droplets->speciesMassFlux(mf_zflux, m_idx_w, 2);

    for ( MFIter mfi((*m_mic_fab_vars[MicVar_SD::rain_accum]),TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box bx = mfi.tilebox();
        const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
        const Array4<Real>& rain_accum_arr = m_mic_fab_vars[MicVar_SD::rain_accum]->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (k == k_lo) {
                auto rain_accum = std::max(0.0, -zflux_arr(i,j,k)*dt/mat_density);
                rain_accum_arr(i,j,k) += (rain_accum * 1000.0 /* [m] -> [mm] */);
            }
        });
    }
}

/*! Compute snow accumulation */
void SuperDropletsMoist::snowAccumulation ()
{
    BL_PROFILE("SuperDropletsMoist::snowAccumulation()");
    if (m_with_ice) {
        auto domain = m_geom.Domain();
        int k_lo = domain.smallEnd(2);
        auto dt = m_dt;

        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::ice);
        auto mat_density = vapour_mat.m_density;

        MultiFab mf_zflux( m_mic_fab_vars[MicVar_SD::snow_accum]->boxArray(),
                           m_mic_fab_vars[MicVar_SD::snow_accum]->DistributionMap(),
                           1,
                           m_mic_fab_vars[MicVar_SD::snow_accum]->nGrowVect() );
        m_super_droplets->speciesMassFlux(mf_zflux, m_idx_i, 2);

        for ( MFIter mfi((*m_mic_fab_vars[MicVar_SD::snow_accum]),TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
            const Array4<Real>& snow_accum_arr = m_mic_fab_vars[MicVar_SD::snow_accum]->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (k == k_lo) {
                    auto snow_accum = std::max(0.0, -zflux_arr(i,j,k)*dt/mat_density);
                    snow_accum_arr(i,j,k) += (snow_accum * 1000.0 /* [m] -> [mm] */);
                }
            });
        }
    }
}

/*! compute condensate mixing ratio */
void SuperDropletsMoist::computeQcSpecies (const int a_i)
{
    BL_PROFILE("SuperDropletsMoist::computeQcSpecies()");
    m_super_droplets->speciesMassDensity( *(m_mic_fab_vars[s_qc_idx(a_i,m_istart_sp)]), a_i );
    if (m_dimensionality == SDMSimulationDim::one_d_z) {
        for ( MFIter mfi(*m_mic_fab_vars[s_qc_idx(a_i,m_istart_sp)]); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            int imin = bx.smallEnd(0);
            int jmin = bx.smallEnd(1);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(a_i,m_istart_sp)]->array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            { q_c_arr(i,j,k) = q_c_arr(imin,jmin,k); });
        }
    }
    densityToRatio(*(m_mic_fab_vars[s_qc_idx(a_i,m_istart_sp)]));
}

/*! compute qt (total) */
void SuperDropletsMoist::computeQtSpecies (const int a_i)
{
    BL_PROFILE("SuperDropletsMoist::computeQtSpecies()");
    for ( MFIter mfi(*m_mic_fab_vars[s_qt_idx(a_i,m_istart_sp)]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(m_mic_fab_vars[MicVar_SD::q_t]->nGrowVect());

        auto q_c_arr = m_mic_fab_vars[s_qc_idx(a_i,m_istart_sp)]->const_array(mfi);
        auto q_v_arr = m_mic_fab_vars[s_qv_idx(a_i,m_istart_sp)]->const_array(mfi);
        auto q_t_arr = m_mic_fab_vars[s_qt_idx(a_i,m_istart_sp)]->array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k); });
    }
}

/*! Compute ground accumulation for non-water species */
void SuperDropletsMoist::speciesAccumulation ()
{
    BL_PROFILE("SuperDropletsMoist::speciesAccumulation()");
    auto domain = m_geom.Domain();
    const auto dx = m_geom.CellSizeArray();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    for (int is = m_istart_sp; is < m_num_species; is++) {
        MultiFab mf_zflux( m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]->boxArray(),
                           m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]->DistributionMap(),
                           1,
                           m_mic_fab_vars[s_sr_idx(is,m_istart_sp)]->nGrowVect() );
        m_super_droplets->speciesMassFlux(mf_zflux, is, 2);

        for ( MFIter mfi((*m_mic_fab_vars[MicVar_SD::rain_accum]),TilingIfNotGPU());
              mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
            const Array4<Real>& accum_arr = m_mic_fab_vars[s_accum_idx(is,m_istart_sp)]->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (k == k_lo) {
                    auto accum = std::max(0.0, -zflux_arr(i,j,k)*dt*dx[0]*dx[1]);
                    accum_arr(i,j,k) += accum;
                }
            });
        }
    }
}

#endif

