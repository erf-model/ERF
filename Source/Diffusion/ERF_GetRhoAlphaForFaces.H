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
 * @param[in] eddy_comp eddy-viscosity component to average (EddyDiff::Mom_v, or
 *            EddyDiff::Mom_h for the terrain-metric coefficient)
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void getRhoAlphaForFaces (int i, int j, int k, int ioff, int joff,
                  amrex::Real& rhoAlpha_lo, amrex::Real& rhoAlpha_hi,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& mu_turb,
                  const amrex::Real mu_eff,
                  bool l_consA, bool l_turb,
                  const int eddy_comp = EddyDiff::Mom_v)
{
    if (l_consA && l_turb) {
        rhoAlpha_lo = fourth * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                               + cell_data(i,j,k-1,Rho_comp       ) + cell_data(i-ioff,j-joff,k-1,Rho_comp       ) ) * mu_eff
                    + fourth * (   mu_turb(i,j,k  ,eddy_comp) +   mu_turb(i-ioff,j-joff,k  ,eddy_comp)
                               +   mu_turb(i,j,k-1,eddy_comp) +   mu_turb(i-ioff,j-joff,k-1,eddy_comp) );
        rhoAlpha_hi = fourth * ( cell_data(i,j,k  ,Rho_comp       ) + cell_data(i-ioff,j-joff,k  ,Rho_comp       )
                               + cell_data(i,j,k+1,Rho_comp       ) + cell_data(i-ioff,j-joff,k+1,Rho_comp       ) ) * mu_eff
                    + fourth * (   mu_turb(i,j,k  ,eddy_comp) +   mu_turb(i-ioff,j-joff,k  ,eddy_comp)
                               +   mu_turb(i,j,k+1,eddy_comp) +   mu_turb(i-ioff,j-joff,k+1,eddy_comp) );
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = mu_eff
                    + fourth * (   mu_turb(i,j,k  ,eddy_comp) +   mu_turb(i-ioff,j-joff,k  ,eddy_comp)
                               +   mu_turb(i,j,k-1,eddy_comp) +   mu_turb(i-ioff,j-joff,k-1,eddy_comp) );
        rhoAlpha_hi = mu_eff
                    + fourth * (   mu_turb(i,j,k  ,eddy_comp) +   mu_turb(i-ioff,j-joff,k  ,eddy_comp)
                               +   mu_turb(i,j,k+1,eddy_comp) +   mu_turb(i-ioff,j-joff,k+1,eddy_comp) );
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
 * Return rhoAlpha at the cell centers below (k-1) and above (k) a z-staggered face: the
 * coefficient of dw/dz in the correction stress tau33i (ComputeStress*), which carries
 * twice the eddy viscosity.
 *
 * @param[in] i face i-index
 * @param[in] j face j-index
 * @param[in] k face k-index
 * @param[out] rhoAlpha_lo rhoAlpha below the face
 * @param[out] rhoAlpha_hi rhoAlpha above the face
 * @param[in] cell_data conserved cell-centered state
 * @param[in] mu_turb turbulent viscosity
 * @param[in] mu_eff effective molecular viscosity, with the factor of 2 of the w
 *            (diagonal) correction stress
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
                    + two * mu_turb(i,j,k-1,EddyDiff::Mom_v);
        rhoAlpha_hi = cell_data(i,j,k  ,Rho_comp) * mu_eff
                    + two * mu_turb(i,j,k  ,EddyDiff::Mom_v);
    }
    else if (l_turb) // with MolecDiffType::Constant or None
    {
        rhoAlpha_lo = mu_eff + two * mu_turb(i,j,k-1,EddyDiff::Mom_v);
        rhoAlpha_hi = mu_eff + two * mu_turb(i,j,k  ,EddyDiff::Mom_v);
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

/**
 * rhoAlpha on the lower and upper z-faces of the control volume of row k of the implicit
 * momentum solve for the u (stagdir 0), v (stagdir 1) or w (stagdir 2) column at (i,j):
 * the z-edges at k and k+1 for u and v, and the cells k-1 and k for w.
 */
template <int stagdir>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void getRhoAlphaForMomRow (int i, int j, int k,
                  amrex::Real& rhoAlpha_lo, amrex::Real& rhoAlpha_hi,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& mu_turb,
                  const amrex::Real mu_eff,
                  bool l_consA, bool l_turb)
{
    if constexpr (stagdir == 2) {
        getRhoAlphaForFaces(i, j, k, rhoAlpha_lo, rhoAlpha_hi,
                            cell_data, mu_turb, mu_eff, l_consA, l_turb);
    } else {
        getRhoAlphaForFaces(i, j, k, (stagdir == 0) ? 1 : 0, (stagdir == 1) ? 1 : 0,
                            rhoAlpha_lo, rhoAlpha_hi,
                            cell_data, mu_turb, mu_eff, l_consA, l_turb);
    }
}

/**
 * Density on face k of the u (stagdir 0), v (stagdir 1) or w (stagdir 2) column at (i,j),
 * averaged from the two cells on either side of it, as FillCoarseFineGhostsForImplicit
 * forms the momenta.
 */
template <int stagdir>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real getRhoForMomFace (int i, int j, int k,
                              const amrex::Array4<const amrex::Real>& cell_data)
{
    constexpr int ioff = (stagdir == 0) ? 1 : 0;
    constexpr int joff = (stagdir == 1) ? 1 : 0;
    constexpr int koff = (stagdir == 2) ? 1 : 0;
    return myhalf * (cell_data(i,j,k,Rho_comp) + cell_data(i-ioff,j-joff,k-koff,Rho_comp));
}
