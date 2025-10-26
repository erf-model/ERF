    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    DiffChoice diffChoice = solverChoice.diffChoice;
    TurbChoice turbChoice = solverChoice.turbChoice[level];

    bool l_consA  = (diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = turbChoice.use_kturb;

    // Theta, KE, Scalar
    Vector<Real> alpha_eff(NPRIMVAR_max, 0.0);
    if (l_consA) {
        for (int i = 0; i < NPRIMVAR_max; ++i) {
           switch (i) {
               case PrimTheta_comp:
                    alpha_eff[PrimTheta_comp] = diffChoice.alpha_T;
                    break;
               case PrimScalar_comp:
                    alpha_eff[PrimScalar_comp] = diffChoice.alpha_C;
                    break;
               case PrimQ1_comp:
                    alpha_eff[PrimQ1_comp] = diffChoice.alpha_C;
                    break;
               case PrimQ2_comp:
                    alpha_eff[PrimQ2_comp] = diffChoice.alpha_C;
                    break;
               case PrimQ3_comp:
                   alpha_eff[PrimQ3_comp] = diffChoice.alpha_C;
                   break;
               case PrimQ4_comp:
                   alpha_eff[PrimQ4_comp] = diffChoice.alpha_C;
                   break;
               case PrimQ5_comp:
                   alpha_eff[PrimQ5_comp] = diffChoice.alpha_C;
                   break;
               case PrimQ6_comp:
                   alpha_eff[PrimQ6_comp] = diffChoice.alpha_C;
                   break;
               default:
                    alpha_eff[i] = 0.0;
                    break;
          }
       }
    } else {
        for (int i = 0; i < NPRIMVAR_max; ++i) {
           switch (i) {
               case PrimTheta_comp:
                    alpha_eff[PrimTheta_comp] = diffChoice.rhoAlpha_T;
                    break;
               case PrimScalar_comp:
                    alpha_eff[PrimScalar_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ1_comp:
                    alpha_eff[PrimQ1_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ2_comp:
                    alpha_eff[PrimQ2_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ3_comp:
                    alpha_eff[PrimQ3_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ4_comp:
                   alpha_eff[PrimQ4_comp] = diffChoice.rhoAlpha_C;
                   break;
               case PrimQ5_comp:
                   alpha_eff[PrimQ5_comp] = diffChoice.rhoAlpha_C;
                   break;
               case PrimQ6_comp:
                   alpha_eff[PrimQ6_comp] = diffChoice.rhoAlpha_C;
                   break;
               default:
                    alpha_eff[i] = 0.0;
                    break;
          }
       }
    }

//  Vector<int> eddy_diff_idx{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::Scalar_h,
//                            EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h ,
//                            EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h };
//  Vector<int> eddy_diff_idy{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::Scalar_h,
//                            EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h ,
//                            EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h };
    Vector<int> eddy_diff_idz{EddyDiff::Theta_v, EddyDiff::KE_v, EddyDiff::Scalar_v,
                              EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v ,
                              EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v };

    // Device vectors
    Gpu::AsyncVector<Real> alpha_eff_d;
//  Gpu::AsyncVector<int>  eddy_diff_idx_d;
//  Gpu::AsyncVector<int>  eddy_diff_idy_d;
    Gpu::AsyncVector<int>  eddy_diff_idz_d;
    alpha_eff_d.resize(alpha_eff.size());
//  eddy_diff_idx_d.resize(eddy_diff_idx.size());
//  eddy_diff_idy_d.resize(eddy_diff_idy.size());
    eddy_diff_idz_d.resize(eddy_diff_idz.size());

    Gpu::copy(Gpu::hostToDevice, alpha_eff.begin()    , alpha_eff.end()    , alpha_eff_d.begin());
//  Gpu::copy(Gpu::hostToDevice, eddy_diff_idx.begin(), eddy_diff_idx.end(), eddy_diff_idx_d.begin());
//  Gpu::copy(Gpu::hostToDevice, eddy_diff_idy.begin(), eddy_diff_idy.end(), eddy_diff_idy_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idz.begin(), eddy_diff_idz.end(), eddy_diff_idz_d.begin());

    // Capture pointers for device code
    Real* d_alpha_eff     = alpha_eff_d.data();
//  int*  d_eddy_diff_idx = eddy_diff_idx_d.data();
//  int*  d_eddy_diff_idy = eddy_diff_idy_d.data();
    int*  d_eddy_diff_idz = eddy_diff_idz_d.data();
