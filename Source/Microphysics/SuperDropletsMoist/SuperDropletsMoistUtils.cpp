#include "SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

/*! Copy moisture model variables from the conserved state vector to the
    member multifabs in this object */
void SuperDropletsMoist::Copy_State_to_Micro (  const MultiFab& a_cons_vars /*!< Conserved variables */)
{
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto states_arr = a_cons_vars.array(mfi);

        auto rho_arr = m_mic_fab_vars[MicVar_SD::rho]->array(mfi);
        auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);
        auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->array(mfi);
        auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->array(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_arr(i,j,k) = states_arr(i,j,k,Rho_comp);
            q_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp) / states_arr(i,j,k,Rho_comp);
            q_c_arr(i,j,k) = states_arr(i,j,k,RhoQ2_comp) / states_arr(i,j,k,Rho_comp);
            q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k);
        });
    }
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

