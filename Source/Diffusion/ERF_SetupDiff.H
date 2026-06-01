    // include common definitions used in the implicit diffusion solve
#include "ERF_SetupVertDiff.H"

    const Real dx_inv = cellSizeInv[0];
    const Real dy_inv = cellSizeInv[1];

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);

    [[maybe_unused]] const int end_comp   = start_comp + num_comp - 1;
    [[maybe_unused]] bool l_use_keqn = turbChoice.use_keqn;
    [[maybe_unused]] bool l_use_mynn = turbChoice.use_pbl_tke;

    // used by ERF_AddTKESources.H
    [[maybe_unused]] const bool use_ref_theta = (turbChoice.theta_ref > 0);
    [[maybe_unused]] const Real l_inv_theta0 = (use_ref_theta) ? one / turbChoice.theta_ref : one;

    // vertical diffusivities are defined in ERF_SetupVertDiff.H
    // The number of quantities following Scalar_h must be equal to the number
    //   of quantities in the conserved state following RhoScalar_comp -- this is
    //   currently all the moist quantities, which tops out at 11 for the Morrison model
    Vector<int> eddy_diff_idx{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::Scalar_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h};
    Vector<int> eddy_diff_idy{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::Scalar_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h,  EddyDiff::Q_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h};

    // Device vectors
    Gpu::AsyncVector<int>  eddy_diff_idx_d;
    Gpu::AsyncVector<int>  eddy_diff_idy_d;

    eddy_diff_idx_d.resize(eddy_diff_idx.size());
    eddy_diff_idy_d.resize(eddy_diff_idy.size());

    Gpu::copy(Gpu::hostToDevice, eddy_diff_idx.begin(), eddy_diff_idx.end(), eddy_diff_idx_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idy.begin(), eddy_diff_idy.end(), eddy_diff_idy_d.begin());

    // Capture pointers for device code
    int*  d_eddy_diff_idx = eddy_diff_idx_d.data();
    int*  d_eddy_diff_idy = eddy_diff_idy_d.data();
