#include "SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES

/*! Copy moisture model variables from the conserved state vector to the
    member multifabs in this object */
void SuperDropletsMoist::Copy_State_to_Micro (  const MultiFab& a_cons_vars /*!< Conserved variables */)
{
    for ( MFIter mfi(a_cons_vars); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto states_arr = a_cons_vars.array(mfi);

        auto rho_v_arr = m_mic_fab_vars[MicVar_SD::rho_v]->array(mfi);
        auto rho_c_arr = m_mic_fab_vars[MicVar_SD::rho_c]->array(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rho_v_arr(i,j,k) = states_arr(i,j,k,RhoQ1_comp);
            rho_c_arr(i,j,k) = states_arr(i,j,k,RhoQ2_comp);
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

        auto rho_v_arr = m_mic_fab_vars[MicVar_SD::rho_v]->array(mfi);
        auto rho_c_arr = m_mic_fab_vars[MicVar_SD::rho_c]->array(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            states_arr(i,j,k,RhoQ1_comp) = rho_v_arr(i,j,k);
            states_arr(i,j,k,RhoQ2_comp) = rho_c_arr(i,j,k);
        });
    }
}

#endif

