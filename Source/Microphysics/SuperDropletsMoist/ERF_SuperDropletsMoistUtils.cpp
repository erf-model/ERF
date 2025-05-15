#include "ERF_SuperDropletsMoist.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"

#ifdef ERF_USE_PARTICLES

/*! Copy moisture model variables from the conserved state vector to the
    member multifabs in this object, compute and save pressure and temperature */
void SuperDropletsMoist::Copy_State_to_Micro (  const MultiFab& a_cons_vars /*!< Conserved variables */)
{
    const auto& gvec = a_cons_vars.nGrowVect();

    // Copy density and vapour mixing ratio from state variables
    // Note: do *not* copy qc
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto states_arr = a_cons_vars.const_array(mfi);
        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->array(mfi);
        auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);
        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
        auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
        auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_arr(i,j,k) = states_arr(i,j,k,Rho_comp);
            theta_arr(i,j,k) = states_arr(i,j,k,RhoTheta_comp)/states_arr(i,j,k,Rho_comp);
        });

        if (m_update_qv) {
            ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            { q_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp) / states_arr(i,j,k,Rho_comp); });
        }

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k) + q_r_arr(i,j,k); });

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

    // Get vapour material properties object
    auto& vapour_mat = m_super_droplets->getVapourMaterial();

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

    for (auto i(0); i < MicVar_SD::NumVars; i++) {
        m_mic_fab_vars[i]->FillBoundary(m_geom.periodicity());
    }
}

/*! Copy moisture model variables to the conserved state vector from the
    member multifabs in this object */
void SuperDropletsMoist::Copy_Micro_to_State (  MultiFab& a_cons_vars /*!< Conserved variables */)
{
    const auto& gvec = a_cons_vars.nGrowVect();

    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto states_arr = a_cons_vars.array(mfi);
        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
        auto q_r_arr = m_mic_fab_vars[MicVar_SD::q_r]->const_array(mfi);
        auto theta_arr = m_mic_fab_vars[MicVar_SD::theta]->const_array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoTheta_comp) = states_arr(i,j,k,Rho_comp)*theta_arr(i,j,k);
            states_arr(i,j,k,RhoQ1_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
            states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
            states_arr(i,j,k,RhoQ3_comp) = states_arr(i,j,k,Rho_comp)*q_r_arr(i,j,k);
        });
    }

    a_cons_vars.FillBoundary(m_geom.periodicity());
}

/*! Update microphysics variables */
void SuperDropletsMoist::Update_Micro_Vars (MultiFab& a_cons_vars)
{
    Copy_State_to_Micro(a_cons_vars);
}

/*! Compute qc, qt, and rain accumulation,and update state variables
 *  from microphysics variables */
void SuperDropletsMoist::Update_State_Vars (MultiFab& a_cons_vars)
{
    computeQcQr();
    computeQt();
    rainAccumulation();
    if (!m_kinematic_mode) { Copy_Micro_to_State(a_cons_vars); }
}

/*! Convert a multifab containing density of something to its mixing ratio */
void SuperDropletsMoist::densityToRatio (  MultiFab& a_var, /*!< Multifab */
                                           const int a_comp /*!< Component */ )
{
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

/*! compute condensate mixing ratio */
void SuperDropletsMoist::computeQcQr ()
{
    m_super_droplets->massDensityCondensate(*(m_mic_fab_vars[MicVar_SD::q_c]), 0, m_r_rain);
    m_super_droplets->massDensityCondensate(*(m_mic_fab_vars[MicVar_SD::q_r]), m_r_rain, 1.0);

    if (m_dimensionality == SDMSimulationDim::one_d_z) {
        for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_t]); mfi.isValid(); ++mfi) {
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

/*! compute qt (total) */
void SuperDropletsMoist::computeQt ()
{
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
    auto domain = m_geom.Domain();
    int k_lo = domain.smallEnd(2);
    auto dt = m_dt;

    auto& vapour_mat = m_super_droplets->getVapourMaterial();
    auto mat_density = vapour_mat.m_density;

    MultiFab mf_zflux( m_mic_fab_vars[MicVar_SD::rain_accum]->boxArray(),
                       m_mic_fab_vars[MicVar_SD::rain_accum]->DistributionMap(),
                       1,
                       m_mic_fab_vars[MicVar_SD::rain_accum]->nGrowVect() );
    m_super_droplets->massFluxCondensate(mf_zflux, 2);

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

#endif

