/* This follows the rhoAlpha calculations in ERF_DiffusionSrcForState_*.cpp
 * and averages rhoAlpha from cell centers to z faces */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void getRhoAlpha (int i, int j, int k,
                  amrex::Real& rhoAlpha_lo, amrex::Real& rhoAlpha_hi,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& mu_turb,
                  const amrex::Real* d_alpha_eff,
                  const int* d_eddy_diff_idz,
                  int prim_index, int prim_scal_index,
                  bool l_consA, bool l_turb)
{
    if (l_consA && l_turb) {
        rhoAlpha_lo = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) ) * d_alpha_eff[prim_scal_index]
                    + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_scal_index])
                            + mu_turb(i,j,k-1, d_eddy_diff_idz[prim_scal_index]) );
        rhoAlpha_hi = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k+1,Rho_comp) ) * d_alpha_eff[prim_scal_index]
                    + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_scal_index])
                            + mu_turb(i,j,k+1, d_eddy_diff_idz[prim_scal_index]) );
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = d_alpha_eff[prim_index]
                    + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_index])
                            + mu_turb(i,j,k-1, d_eddy_diff_idz[prim_index]) );
        rhoAlpha_hi =  d_alpha_eff[prim_index]
                    + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_index])
                            + mu_turb(i,j,k+1, d_eddy_diff_idz[prim_index]) );
    }
    else if (l_consA) // without an LES/PBL model
    {
        rhoAlpha_lo = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) ) * d_alpha_eff[prim_index];
        rhoAlpha_hi = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k+1,Rho_comp) ) * d_alpha_eff[prim_index];
    }
    else // with MolecDiffType::Constant or None - without an LES/PBL model
    {
        rhoAlpha_lo = d_alpha_eff[prim_index];
        rhoAlpha_hi = d_alpha_eff[prim_index];
    }
}
