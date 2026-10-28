/** \cond */
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);
    amrex::ignore_unused(dom_lo, dom_hi);

    DiffChoice diffChoice = solverChoice.diffChoice;
    TurbChoice turbChoice = solverChoice.turbChoice[level];

    bool l_consA  = (diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = turbChoice.use_kturb;

    // Molecular diffusivity by primitive component (theta, every passive
    //   scalar and the moisture species, number concentrations included), the
    //   coefficients the explicit paths take from ResolveNativeScalarDiffusionPolicy()
    Vector<Real> alpha_eff = VerticalMolecularAlphaByPrim(diffChoice);

    // Vertical eddy diffusivity component by primitive component.  Every
    //   component above the passive scalars takes the Q category, including the
    //   non-moist species a scheme may append; the table spans NPRIMVAR_max,
    //   which tracks NVAR_max, so it always covers the full primitive range the
    //   callers can request.  Every scalar slot (base and smoke) takes the
    //   Scalar category.
    Vector<int> eddy_diff_idz(NPRIMVAR_max, EddyDiff::Q_v);
    eddy_diff_idz[PrimTheta_comp] = EddyDiff::Theta_v;
    eddy_diff_idz[PrimKE_comp]    = EddyDiff::KE_v;
    for (int ns = 0; ns < NSCALARS; ++ns) {
        eddy_diff_idz[PrimScalar_comp + ns] = EddyDiff::Scalar_v;
    }

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
