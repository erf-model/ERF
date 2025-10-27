            Real rhoAlpha_lo;
            Real rhoAlpha_hi;

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
