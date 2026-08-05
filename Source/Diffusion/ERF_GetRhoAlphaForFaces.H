/**
 * Average rhoAlpha from cell centers to lower and upper z-edges of staggered faces.
 *
 * @param[in] i face i-index
 * @param[in] j face j-index
 * @param[in] k face k-index
 * @param[in] ioff x-direction staggering offset
 * @param[in] joff y-direction staggering offset
 * @param[out] rhoAlpha_lo rhoAlpha averaged to the lower z-edge
 * @param[out] rhoAlpha_hi rhoAlpha averaged to the upper z-edge
 * @param[in] cell_data conserved cell-centered state
 * @param[in] mu_turb turbulent viscosity
 * @param[in] mu_eff effective molecular viscosity. NOTE: unlike the mu_eff in
 *            ERF_MakeTauTerms.cpp and ComputeStress*, this includes the factor
 *            of 2 only for the w (diagonal) correction stress, because the u/v
 *            (off-diagonal) correction strains already carry a factor of 1/2.
 * @param[in] l_consA flag for conservative constant-alpha diffusion
 * @param[in] l_turb flag for turbulent diffusion
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void getRhoAlphaForFaces (int i, int j, int k, int ioff, int joff,
                  amrex::Real& rhoAlpha_lo, amrex::Real& rhoAlpha_hi,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& mu_turb,
                  const amrex::Real mu_eff,
                  bool l_consA, bool l_turb)
{
    if (l_consA && l_turb) {
        rhoAlpha_lo = fourth * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                               + cell_data(i,j,k-1,Rho_comp       ) + cell_data(i-ioff,j-joff,k-1,Rho_comp       ) ) * mu_eff
                    + fourth * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                               +   mu_turb(i,j,k-1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k-1,EddyDiff::Mom_v) );
        rhoAlpha_hi = fourth * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                               + cell_data(i,j,k+1,Rho_comp       ) + cell_data(i-ioff,j-joff,k+1,Rho_comp       ) ) * mu_eff
                    + fourth * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                               +   mu_turb(i,j,k+1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k+1,EddyDiff::Mom_v) );
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = mu_eff
                    + fourth * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                               +   mu_turb(i,j,k-1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k-1,EddyDiff::Mom_v) );
        rhoAlpha_hi = mu_eff
                    + fourth * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                               +   mu_turb(i,j,k+1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k+1,EddyDiff::Mom_v) );
    }
    else if (l_consA) // without an LES/PBL model
    {
        rhoAlpha_lo = fourth * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                               + cell_data(i,j,k-1,Rho_comp       ) + cell_data(i-ioff,j-joff,k-1,Rho_comp       ) ) * mu_eff;
        rhoAlpha_hi = fourth * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                               + cell_data(i,j,k+1,Rho_comp       ) + cell_data(i-ioff,j-joff,k+1,Rho_comp       ) ) * mu_eff;
    }
    else // with MolecDiffType::Constant or None - without an LES/PBL model
    {
        rhoAlpha_lo = mu_eff;
        rhoAlpha_hi = mu_eff;
    }
}

/**
 * Return rhoAlpha at cell centers below and above a z-staggered face.
 *
 * @param[in] i face i-index
 * @param[in] j face j-index
 * @param[in] k face k-index
 * @param[out] rhoAlpha_lo rhoAlpha below the face
 * @param[out] rhoAlpha_hi rhoAlpha above the face
 * @param[in] cell_data conserved cell-centered state
 * @param[in] mu_turb turbulent viscosity
 * @param[in] mu_eff effective molecular viscosity. NOTE: unlike the mu_eff in
 *            ERF_MakeTauTerms.cpp and ComputeStress*, this includes the factor
 *            of 2 only for the w (diagonal) correction stress, because the u/v
 *            (off-diagonal) correction strains already carry a factor of 1/2.
 * @param[in] l_consA flag for conservative constant-alpha diffusion
 * @param[in] l_turb flag for turbulent diffusion
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void getRhoAlphaForFaces (int i, int j, int k,
                  amrex::Real& rhoAlpha_lo, amrex::Real& rhoAlpha_hi,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& mu_turb,
                  const amrex::Real mu_eff,
                  bool l_consA, bool l_turb)
{
    if (l_consA && l_turb) {
        rhoAlpha_lo = cell_data(i,j,k-1,Rho_comp) * mu_eff
                    +   mu_turb(i,j,k-1,EddyDiff::Mom_v);
        rhoAlpha_hi = cell_data(i,j,k  ,Rho_comp) * mu_eff
                    +   mu_turb(i,j,k  ,EddyDiff::Mom_v);
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = mu_eff + mu_turb(i,j,k-1,EddyDiff::Mom_v);
        rhoAlpha_hi = mu_eff + mu_turb(i,j,k  ,EddyDiff::Mom_v);
    }
    else if (l_consA) // without an LES/PBL model
    {
        rhoAlpha_lo = cell_data(i,j,k-1,Rho_comp) * mu_eff;
        rhoAlpha_hi = cell_data(i,j,k  ,Rho_comp) * mu_eff;
    }
    else // with MolecDiffType::Constant or None - without an LES/PBL model
    {
        rhoAlpha_lo = mu_eff;
        rhoAlpha_hi = mu_eff;
    }
}
