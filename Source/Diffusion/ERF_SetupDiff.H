    // include common definitions used in the implicit diffusion solve
#include "ERF_SetupVertDiff.H"

    const Real dx_inv = cellSizeInv[0];
    const Real dy_inv = cellSizeInv[1];

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);

    const int end_comp   = start_comp + num_comp - 1;

    bool l_use_keqn = turbChoice.use_keqn;
    bool l_use_mynn = turbChoice.use_pbl_tke;

    // used by ERF_AddTKESources.H
    Real l_inv_theta0 = (turbChoice.theta_ref > 0) ? 1.0 / turbChoice.theta_ref : 1.0;

    // vertical diffusivities are defined in ERF_SetupVertDiff.H
    Vector<int> eddy_diff_idx{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::Scalar_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h ,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h };
    Vector<int> eddy_diff_idy{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::Scalar_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h ,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h };
//  Vector<int> eddy_diff_idz{EddyDiff::Theta_v, EddyDiff::KE_v, EddyDiff::Scalar_v,
//                            EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v ,
//                            EddyDiff::Q_v    , EddyDiff::Q_v,  EddyDiff::Q_v };

    // Device vectors
//  Gpu::AsyncVector<Real> alpha_eff_d;
    Gpu::AsyncVector<int>  eddy_diff_idx_d;
    Gpu::AsyncVector<int>  eddy_diff_idy_d;
//  Gpu::AsyncVector<int>  eddy_diff_idz_d;
//  alpha_eff_d.resize(alpha_eff.size());
    eddy_diff_idx_d.resize(eddy_diff_idx.size());
    eddy_diff_idy_d.resize(eddy_diff_idy.size());
//  eddy_diff_idz_d.resize(eddy_diff_idz.size());

//  Gpu::copy(Gpu::hostToDevice, alpha_eff.begin()    , alpha_eff.end()    , alpha_eff_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idx.begin(), eddy_diff_idx.end(), eddy_diff_idx_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idy.begin(), eddy_diff_idy.end(), eddy_diff_idy_d.begin());
//  Gpu::copy(Gpu::hostToDevice, eddy_diff_idz.begin(), eddy_diff_idz.end(), eddy_diff_idz_d.begin());

    // Capture pointers for device code
//  Real* d_alpha_eff     = alpha_eff_d.data();
    int*  d_eddy_diff_idx = eddy_diff_idx_d.data();
    int*  d_eddy_diff_idy = eddy_diff_idy_d.data();
//  int*  d_eddy_diff_idz = eddy_diff_idz_d.data();
