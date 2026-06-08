    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);
    amrex::ignore_unused(dom_lo, dom_hi);

    DiffChoice diffChoice = solverChoice.diffChoice;
    TurbChoice turbChoice = solverChoice.turbChoice[level];

    bool l_consA  = (diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = turbChoice.use_kturb;

    // Theta, KE, Scalar
    Vector<Real> alpha_eff(NPRIMVAR_max, zero);
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
                    alpha_eff[i] = zero;
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
                    alpha_eff[i] = zero;
                    break;
          }
       }
    }

    // The number of quantities following Scalar_h must be equal to the number
    //   of quantities in the conserved state following RhoScalar_comp -- this is
    //   currently all the moist quantities, which tops out at 11 for the Morrison model
    Vector<int> eddy_diff_idz{EddyDiff::Theta_v, EddyDiff::KE_v, EddyDiff::Scalar_v,
                              EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v,
                              EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v,
                              EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v,
                              EddyDiff::Q_v    , EddyDiff::Q_v};

    // Device vectors
    Gpu::AsyncVector<Real> alpha_eff_d;
    Gpu::AsyncVector<int>  eddy_diff_idz_d;
    alpha_eff_d.resize(alpha_eff.size());
    eddy_diff_idz_d.resize(eddy_diff_idz.size());

    Gpu::copy(Gpu::hostToDevice, alpha_eff.begin()    , alpha_eff.end()    , alpha_eff_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idz.begin(), eddy_diff_idz.end(), eddy_diff_idz_d.begin());

    // Capture pointers for device code
    Real* d_alpha_eff     = alpha_eff_d.data();
    [[maybe_unused]] int*  d_eddy_diff_idz = eddy_diff_idz_d.data();
