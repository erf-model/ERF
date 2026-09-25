/** \cond */
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);
    amrex::ignore_unused(dom_lo, dom_hi);

    DiffChoice diffChoice = solverChoice.diffChoice;
    TurbChoice turbChoice = solverChoice.turbChoice[level];

    bool l_consA  = (diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = turbChoice.use_kturb;

    // Molecular diffusivity by primitive component.  Only theta, the passive
    //   scalar and the water species through Q6 have one defined; everything
    //   else -- KE, the number concentrations, and any non-moist species a
    //   scheme appends after its water species -- is left at zero.  This
    //   matches ResolveNativeScalarDiffusionPolicy() in
    //   ERF_NativeScalarDiffusion.H, which the explicit paths use instead of
    //   these tables.
    Vector<Real> alpha_eff(NPRIMVAR_max, zero);
    {
        const Real alpha_T_eff = l_consA ? diffChoice.alpha_T : diffChoice.rhoAlpha_T;
        const Real alpha_C_eff = l_consA ? diffChoice.alpha_C : diffChoice.rhoAlpha_C;

        alpha_eff[PrimTheta_comp]  = alpha_T_eff;
        alpha_eff[PrimScalar_comp] = alpha_C_eff;
        for (int i = PrimQ1_comp; i <= PrimQ6_comp; ++i) {
            alpha_eff[i] = alpha_C_eff;
        }
    }

    // Vertical eddy diffusivity component by primitive component.  Every
    //   component above the passive scalar takes the Q category, including the
    //   non-moist species a scheme may append; the table spans NPRIMVAR_max,
    //   which tracks NVAR_max, so it always covers the full primitive range the
    //   callers can request.
    Vector<int> eddy_diff_idz(NPRIMVAR_max, EddyDiff::Q_v);
    eddy_diff_idz[PrimTheta_comp]  = EddyDiff::Theta_v;
    eddy_diff_idz[PrimKE_comp]     = EddyDiff::KE_v;
    eddy_diff_idz[PrimScalar_comp] = EddyDiff::Scalar_v;

    // Device vectors
    Gpu::AsyncVector<Real> alpha_eff_d;
    Gpu::AsyncVector<int>  eddy_diff_idz_d;
    alpha_eff_d.resize(alpha_eff.size());
    eddy_diff_idz_d.resize(eddy_diff_idz.size());

    Gpu::copy(Gpu::hostToDevice, alpha_eff.begin()    , alpha_eff.end()    , alpha_eff_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idz.begin(), eddy_diff_idz.end(), eddy_diff_idz_d.begin());

    // Capture pointers for device code
    [[maybe_unused]] Real* d_alpha_eff     = alpha_eff_d.data();
    [[maybe_unused]] int*  d_eddy_diff_idz = eddy_diff_idz_d.data();
/** \endcond */
