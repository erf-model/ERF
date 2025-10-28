        auto get_rhoAlpha = [&](int i, int j, int k, Real& rhoAlpha_lo, Real& rhoAlpha_hi) {
            DiffChoice dc = solverChoice.diffChoice;
            Real mu_eff = (l_consA) ? 2.0 * dc.dynamic_viscosity / dc.rho0_trans
                                    : 2.0 * dc.dynamic_viscosity;

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
        };
