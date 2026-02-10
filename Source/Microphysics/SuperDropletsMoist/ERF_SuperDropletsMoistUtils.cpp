#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! \brief Copy moisture model variables from conserved state to member MultiFabs
 *
 * This function extracts moisture-related variables from the conserved state vector
 * and copies them to the internal member MultiFabs. It also computes derived quantities
 * such as pressure, temperature, and saturation ratios for all defined species.
 *
 * The function performs:
 * 1. Copy of density and potential temperature from state variables
 * 2. Copy and computation of water-related variables (vapor, total water)
 * 3. Copy and computation of other species variables
 * 4. Computation of pressure and temperature fields
 * 5. Computation of saturation ratios for all species
 *
 * \param[in] a_cons_vars MultiFab containing the conserved state variables
 */
void SuperDropletsMoist::Copy_State_to_Micro (const MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Copy_State_to_Micro()");
    const auto& gvec = a_cons_vars.nGrowVect();

    // Copy density and vapour mixing ratio from state variables
    // Note: do *not* copy qc
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

        // water
        {
            auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);
            auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
            auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                q_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp) / states_arr(i,j,k,Rho_comp);
                q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k) + q_r_arr(i,j,k);
            });
        }

        // other species
        for (int is = 1; is < m_num_species; is++) {
            auto q_v_arr = m_mic_fab_vars[s_qv_idx(is)]->array(mfi);
            auto q_t_arr = m_mic_fab_vars[s_qt_idx(is)]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(is)]->const_array(mfi);
            auto qv_comp = q_qv_idx(is);
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

    // water
    {
        // Get vapour material properties object for water
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(Species::Name::H2O);

        // Compute saturation ratio
        vapour_mat.computeSaturationVapFrac( (*m_mic_fab_vars[MicVar_SD::rh]),
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
    }

    // other species
    for (int is = 1; is < m_num_species; is++) {
        auto& vapour_mat = m_super_droplets->getSpeciesMaterial(m_species[is]);
        vapour_mat.computeSaturationVapFrac((*m_mic_fab_vars[s_sr_idx(is)]),
                                            (*m_mic_fab_vars[MicVar_SD::temperature]),
                                            (*m_mic_fab_vars[MicVar_SD::pressure]) );
        for (   MFIter mfi((*m_mic_fab_vars[s_sr_idx(is)]),
                TilingIfNotGPU()); mfi.isValid();
                ++mfi ) {

            Box bx = mfi.tilebox();
            bx.grow( m_mic_fab_vars[s_sr_idx(is)]->nGrowVect() );

            const Array4<Real>& sr_arr = m_mic_fab_vars[s_sr_idx(is)]->array(mfi);
            const Array4<Real const>& qv_arr = m_mic_fab_vars[s_qv_idx(is)]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            { sr_arr(i,j,k,0) = qv_arr(i,j,k,0) / sr_arr(i,j,k,0); });
        }
    }

    for (auto i(0); i < MicVar_SD::NumVars; i++) {
        m_mic_fab_vars[i]->FillBoundary(m_geom.periodicity());
    }
}

/*! \brief Copy moisture model variables to the conserved state vector
 *
 * This function copies moisture-related variables from the internal member MultiFabs
 * to the conserved state vector. It updates:
 * 1. Potential temperature field in conserved variables
 * 2. Water-related fields (vapor, cloud, rain) in conserved variables
 * 3. Other species fields in conserved variables
 *
 * All mixing ratios are converted to density-weighted variables (rho*q) when
 * stored in the conserved state vector.
 *
 * \param[out] a_cons_vars MultiFab containing conserved state variables to be updated
 */
void SuperDropletsMoist::Copy_Micro_to_State (MultiFab& a_cons_vars)
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

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,RhoQ1_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
                states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
                states_arr(i,j,k,RhoQ3_comp) = states_arr(i,j,k,Rho_comp)*q_r_arr(i,j,k);
            });
        }

        // other species
        for (int is = 1; is < m_num_species; is++) {
            auto q_v_arr = m_mic_fab_vars[s_qv_idx(is)]->array(mfi);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(is)]->array(mfi);
            auto qv_comp = q_qv_idx(is);
            auto qc_comp = q_qc_idx(is);
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                states_arr(i,j,k,qv_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
                states_arr(i,j,k,qc_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
            });
        }
    }

    a_cons_vars.FillBoundary(m_geom.periodicity());
}

/*! \brief Update microphysics variables from conserved state
 *
 * This function updates the microphysics variables by copying data
 * from the conserved state variables. It's a convenience wrapper
 * around Copy_State_to_Micro.
 *
 * \param[in] a_cons_vars MultiFab containing conserved state variables
 */
void SuperDropletsMoist::Update_Micro_Vars (MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Update_Micro_Vars()");
    Copy_State_to_Micro(a_cons_vars);
}

/*! \brief Compute derived quantities and update state variables
 *
 * This function computes various derived quantities such as cloud water,
 * total water, and accumulation values, then updates the conserved state
 * variables if not in kinematic mode. It performs:
 * 1. Computation of cloud/rain water and total water for water
 * 2. Calculation of rain accumulation at the ground surface
 * 3. Computation of condensate and total for other species
 * 4. Calculation of species and aerosol accumulation at the ground
 * 5. Update of conserved state variables with computed values (if not in kinematic mode)
 *
 * \param[in,out] a_cons_vars MultiFab containing conserved state variables to be updated
 */
void SuperDropletsMoist::Update_State_Vars (MultiFab& a_cons_vars)
{
    BL_PROFILE("SuperDropletsMoist::Update_State_Vars()");
    computeQcQrWater();
    computeQtWater();
    rainAccumulation();

    computeQcSpecies();
    computeQtSpecies();
    speciesAccumulation();
    aerosolAccumulation();

    if (!m_kinematic_mode) { Copy_Micro_to_State(a_cons_vars); }
}

/*! \brief Convert a density field to a mixing ratio field
 *
 * This function converts a field containing density values (kg/m^3)
 * to mixing ratio values (kg/kg) by dividing by the air density.
 * The operation is performed in-place on the provided MultiFab.
 *
 * \param[in,out] a_var MultiFab containing the field to convert
 * \param[in] a_comp Component index within the MultiFab to convert
 */
void SuperDropletsMoist::densityToRatio (MultiFab& a_var,
                                         const int a_comp)
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

/*! \brief Convert a mixing ratio field to a density field
 *
 * This function converts a field containing mixing ratio values (kg/kg)
 * to density values (kg/m^3) by multiplying by the air density.
 * The operation is performed in-place on the provided MultiFab.
 *
 * \param[in,out] a_var MultiFab containing the field to convert
 * \param[in] a_comp Component index within the MultiFab to convert
 */
void SuperDropletsMoist::ratioToDensity (MultiFab& a_var,
                                         const int a_comp)
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

/*! \brief Compute cloud and rain water mixing ratios from superdroplets
 *
 * This function computes cloud and rain water mixing ratios based on the current
 * superdroplet population. It:
 * 1. Calculates cloud water density (droplets smaller than rain threshold)
 * 2. Calculates rain water density (droplets larger than rain threshold)
 * 3. Handles special dimensionality case for 1D simulations
 * 4. Converts density values to mixing ratios
 *
 * The distinction between cloud and rain water is based on the configured rain threshold
 * radius (m_r_rain).
 */
void SuperDropletsMoist::computeQcQrWater ()
{
    BL_PROFILE("SuperDropletsMoist::computeQcQrWater()");
    m_super_droplets->speciesMassDensity( *(m_mic_fab_vars[MicVar_SD::q_c]),
                                          m_idx_w,
                                          0,
                                          m_r_rain );
    m_super_droplets->speciesMassDensity( *(m_mic_fab_vars[MicVar_SD::q_r]),
                                          m_idx_w,
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

/*! compute qt (total) for water */
void SuperDropletsMoist::computeQtWater ()
{
    BL_PROFILE("SuperDropletsMoist::computeQtWater()");
    for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_t]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(m_mic_fab_vars[MicVar_SD::q_t]->nGrowVect());

        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
        auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
        auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k) + q_r_arr(i,j,k); });
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

    for ( MFIter mfi((*m_mic_fab_vars[MicVar_SD::rain_accum]),TilingIfNotGPU());
          mfi.isValid(); ++mfi ) {
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

/*! compute condensate mixing ratio */
void SuperDropletsMoist::computeQcSpecies (const int a_i)
{
    BL_PROFILE("SuperDropletsMoist::computeQcSpecies()");
    m_super_droplets->speciesMassDensity( *(m_mic_fab_vars[s_qc_idx(a_i)]), a_i );
    if (m_dimensionality == SDMSimulationDim::one_d_z) {
        for ( MFIter mfi(*m_mic_fab_vars[s_qc_idx(a_i)]); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            int imin = bx.smallEnd(0);
            int jmin = bx.smallEnd(1);
            auto q_c_arr = m_mic_fab_vars[s_qc_idx(a_i)]->array(mfi);

            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            { q_c_arr(i,j,k) = q_c_arr(imin,jmin,k); });
        }
    }
    densityToRatio(*(m_mic_fab_vars[s_qc_idx(a_i)]));
}

/*! compute qt (total) */
void SuperDropletsMoist::computeQtSpecies (const int a_i)
{
    BL_PROFILE("SuperDropletsMoist::computeQtSpecies()");
    for ( MFIter mfi(*m_mic_fab_vars[s_qt_idx(a_i)]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(m_mic_fab_vars[MicVar_SD::q_t]->nGrowVect());

        auto q_c_arr = m_mic_fab_vars[s_qc_idx(a_i)]->const_array(mfi);
        auto q_v_arr = m_mic_fab_vars[s_qv_idx(a_i)]->const_array(mfi);
        auto q_t_arr = m_mic_fab_vars[s_qt_idx(a_i)]->array(mfi);

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

    for (int is = 1; is < m_num_species; is++) {
        MultiFab mf_zflux( m_mic_fab_vars[s_sr_idx(is)]->boxArray(),
                           m_mic_fab_vars[s_sr_idx(is)]->DistributionMap(),
                           1,
                           m_mic_fab_vars[s_sr_idx(is)]->nGrowVect() );
        m_super_droplets->speciesMassFlux(mf_zflux, is, 2);

        for ( MFIter mfi((*m_mic_fab_vars[MicVar_SD::rain_accum]),TilingIfNotGPU());
              mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
            const Array4<Real>& accum_arr = m_mic_fab_vars[s_accum_idx(is)]->array(mfi);

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

/*! Compute ground accumulation for non-water species */
void SuperDropletsMoist::aerosolAccumulation ()
{
    BL_PROFILE("SuperDropletsMoist::aerosolAccumulation()");
    auto domain = m_geom.Domain();
    const auto dx = m_geom.CellSizeArray();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    for (int ia = 0; ia < m_num_aerosols; ia++) {
        MultiFab mf_zflux( m_mic_fab_vars[MicVar_SD::rain_accum]->boxArray(),
                           m_mic_fab_vars[MicVar_SD::rain_accum]->DistributionMap(),
                           1,
                           m_mic_fab_vars[MicVar_SD::rain_accum]->nGrowVect() );
        m_super_droplets->aerosolMassFlux(mf_zflux, ia, 2);

        for ( MFIter mfi((*m_mic_fab_vars[MicVar_SD::rain_accum]),TilingIfNotGPU());
              mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            const Array4<Real const>& zflux_arr = mf_zflux.const_array(mfi);
            const Array4<Real>& accum_arr = m_mic_fab_vars[a_accum_idx(m_num_species,ia)]->array(mfi);

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

