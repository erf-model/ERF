/**
 * Average rhoAlpha from cell centers to the lower and upper z-faces.
 *
 * @param[in] i cell-centered i-index
 * @param[in] j cell-centered j-index
 * @param[in] k cell-centered k-index
 * @param[out] rhoAlpha_lo rhoAlpha averaged to the lower z-face
 * @param[out] rhoAlpha_hi rhoAlpha averaged to the upper z-face
 * @param[in] cell_data conserved cell-centered state
 * @param[in] mu_turb turbulent diffusivity
 * @param[in] d_alpha_eff effective molecular diffusivity by primitive component
 * @param[in] d_eddy_diff_idz vertical eddy-diffusivity component map
 * @param[in] prim_index primitive component index
 * @param[in] prim_scal_index scalar primitive component index
 * @param[in] l_consA flag for conservative constant-alpha diffusion
 * @param[in] l_turb flag for turbulent diffusion
 */
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
        rhoAlpha_lo = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) ) * d_alpha_eff[prim_scal_index]
                    + myhalf * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_scal_index])
                               + mu_turb(i,j,k-1, d_eddy_diff_idz[prim_scal_index]) );
        rhoAlpha_hi = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k+1,Rho_comp) ) * d_alpha_eff[prim_scal_index]
                    + myhalf * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_scal_index])
                               + mu_turb(i,j,k+1, d_eddy_diff_idz[prim_scal_index]) );
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = d_alpha_eff[prim_index]
                    + myhalf * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_index])
                               + mu_turb(i,j,k-1, d_eddy_diff_idz[prim_index]) );
        rhoAlpha_hi =  d_alpha_eff[prim_index]
                    + myhalf * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_index])
                               + mu_turb(i,j,k+1, d_eddy_diff_idz[prim_index]) );
    }
    else if (l_consA) // without an LES/PBL model
    {
        rhoAlpha_lo = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) ) * d_alpha_eff[prim_index];
        rhoAlpha_hi = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k+1,Rho_comp) ) * d_alpha_eff[prim_index];
    }
    else // with MolecDiffType::Constant or None - without an LES/PBL model
    {
        rhoAlpha_lo = d_alpha_eff[prim_index];
        rhoAlpha_hi = d_alpha_eff[prim_index];
    }
}
