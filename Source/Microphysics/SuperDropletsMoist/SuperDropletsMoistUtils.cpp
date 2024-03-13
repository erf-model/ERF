#include "SuperDropletsMoist.H"
#include "EOS.H"
#include "IndexDefines.H"

#ifdef ERF_USE_PARTICLES

/*! Copy moisture model variables from the conserved state vector to the
    member multifabs in this object, compute and save pressure and temperature */
void SuperDropletsMoist::Copy_State_to_Micro (  const MultiFab& a_cons_vars /*!< Conserved variables */)
{
    const auto& gvec = a_cons_vars.nGrowVect();

    // Copy density and vapour mixing ratio from state variables
    // Note: do *not* copy condensate mixing ratio
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        bx.grow(gvec);

        auto states_arr = a_cons_vars.array(mfi);

        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->array(mfi);
        auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);
        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);

        ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_arr(i,j,k) = states_arr(i,j,k,Rho_comp);
            q_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp) / states_arr(i,j,k,Rho_comp);
            q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k);
        });
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
            t_arr(i,j,k,0) = getTgivenRandRTh(S_arr(i,j,k,Rho_comp),S_arr(i,j,k,RhoTheta_comp));
            p_arr(i,j,k,0) = getPgivenRTh(S_arr(i,j,k,RhoTheta_comp),qv_arr(i,j,k,0));
        });
    }

    AMREX_ASSERT( !m_mic_fab_vars[MicVar_SD::pressure]->contains_nan() );
    AMREX_ASSERT( !m_mic_fab_vars[MicVar_SD::temperature]->contains_nan() );
}

/*! Copy moisture model variables to the conserved state vector from the
    member multifabs in this object */
void SuperDropletsMoist::Copy_Micro_to_State (  MultiFab& a_cons_vars /*!< Conserved variables */)
{
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto states_arr = a_cons_vars.array(mfi);

        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->array(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoQ1_comp) = states_arr(i,j,k,Rho_comp)*q_v_arr(i,j,k);
            states_arr(i,j,k,RhoQ2_comp) = states_arr(i,j,k,Rho_comp)*q_c_arr(i,j,k);
        });
    }
}

/*! Convert a multifab containing density of something to its mixing ratio */
void SuperDropletsMoist::densityToRatio (  MultiFab& a_var, /*!< Multifab */
                                           const int a_comp /*!< Component */ )
{
    for ( MFIter mfi(a_var); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();

        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->array(mfi);
        auto fab_arr = a_var.array(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            fab_arr(i,j,k,a_comp) /= rho_arr(i,j,k);
        });
    }
}

/*! Convert a multifab containing the mixing ratio of something to its density */
void SuperDropletsMoist::ratioToDensity (  MultiFab& a_var, /*!< Multifab */
                                           const int a_comp /*!< Component */ )
{
    for ( MFIter mfi(a_var); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();

        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->array(mfi);
        auto fab_arr = a_var.array(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            fab_arr(i,j,k,a_comp) *= rho_arr(i,j,k);
        });
    }
}

#endif

