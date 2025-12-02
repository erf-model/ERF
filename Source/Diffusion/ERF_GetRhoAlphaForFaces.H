/* This averages rhoAlpha from cell centers to lo/hi z edges of
 * faces with staggering in i (ioff=1) or j (joff=1) */
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
        rhoAlpha_lo = 0.25 * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                             + cell_data(i,j,k-1,Rho_comp       ) + cell_data(i-ioff,j-joff,k-1,Rho_comp       ) ) * mu_eff
                    + 0.25 * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                             +   mu_turb(i,j,k-1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k-1,EddyDiff::Mom_v) );
        rhoAlpha_hi = 0.25 * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                             + cell_data(i,j,k+1,Rho_comp       ) + cell_data(i-ioff,j-joff,k+1,Rho_comp       ) ) * mu_eff
                    + 0.25 * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                             +   mu_turb(i,j,k+1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k+1,EddyDiff::Mom_v) );
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = mu_eff
                    + 0.25 * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                             +   mu_turb(i,j,k-1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k-1,EddyDiff::Mom_v) );
        rhoAlpha_hi = mu_eff
                    + 0.25 * (   mu_turb(i,j,k  ,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k  ,EddyDiff::Mom_v)
                             +   mu_turb(i,j,k+1,EddyDiff::Mom_v) +   mu_turb(i-ioff,j-joff,k+1,EddyDiff::Mom_v) );
    }
    else if (l_consA) // without an LES/PBL model
    {
        rhoAlpha_lo = 0.25 * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                             + cell_data(i,j,k-1,Rho_comp       ) + cell_data(i-ioff,j-joff,k-1,Rho_comp       ) ) * mu_eff;
        rhoAlpha_hi = 0.25 * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                             + cell_data(i,j,k+1,Rho_comp       ) + cell_data(i-ioff,j-joff,k+1,Rho_comp       ) ) * mu_eff;
    }
    else // with MolecDiffType::Constant or None - without an LES/PBL model
    {
        rhoAlpha_lo = mu_eff;
        rhoAlpha_hi = mu_eff;
    }
}

/* This returns rhoAlpha at cell centers surrounding
 * faces with staggering in k */
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
