#include <ERF_Diffusion.H>
#include "ERF_EddyViscosity.H"

using namespace amrex;

/**
 * Function for computing the strain rates for EB.
 *
 * @param[in] bxcc cell center box for tau_ii
 * @param[in] tbxxy nodal xy box for tau_12
 * @param[in] tbxxz nodal xz box for tau_13
 * @param[in] tbxyz nodal yz box for tau_23
 * @param[in] u x-direction velocity
 * @param[in] v y-direction velocity
 * @param[in] w z-direction velocity
 * @param[out] tau11 11 strain
 * @param[out] tau22 22 strain
 * @param[out] tau33 33 strain
 * @param[out] tau12 12 strain
 * @param[out] tau13 13 strain
 * @param[out] tau23 23 strain
 * @param[in] bc_ptr container with boundary condition types
 * @param[in] dxInv inverse cell size array
 * @param[in] tau13i contribution to strain from du/dz
 * @param[in] tau23i contribution to strain from dv/dz
 */
void
ComputeStrain_EB (const MFIter& mfi,
                 Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Box domain,
                 const Array4<const Real>& u,
                 const Array4<const Real>& v,
                 const Array4<const Real>& w,
                 Array4<Real>& tau11,
                 Array4<Real>& tau22,
                 Array4<Real>& tau33,
                 Array4<Real>& tau12,
                 Array4<Real>& tau13,
                 Array4<Real>& tau23,
                 const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                 const BCRec* bc_ptr,
                 const eb_& ebfact,
                 Array4<Real>& tau13i, Array4<Real>& tau23i)
{
    // Convert domain to each index type to test if we are on Dirichlet boundary
    Box domain_xy = convert(domain, tbxxy.ixType());
    Box domain_xz = convert(domain, tbxxz.ixType());
    Box domain_yz = convert(domain, tbxyz.ixType());

    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // EB
    // Array4<const EBCellFlag> cflag = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const EBCellFlag> u_cflag = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const EBCellFlag> v_cflag = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const EBCellFlag> w_cflag = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();

    // Dirichlet on left or right plane
    bool xl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
         xl_v_dir = ( xl_v_dir && (tbxxy.smallEnd(0) == domain_xy.smallEnd(0)) );

    bool xh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );
         xh_v_dir = ( xh_v_dir && (tbxxy.bigEnd(0) == domain_xy.bigEnd(0)) );

    bool xl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
         xl_w_dir = ( xl_w_dir && (tbxxz.smallEnd(0) == domain_xz.smallEnd(0)) );

    bool xh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );
         xh_w_dir = ( xh_w_dir && (tbxxz.bigEnd(0) == domain_xz.bigEnd(0)) );

    // Dirichlet on front or back plane
    bool yl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
         yl_u_dir = ( yl_u_dir && (tbxxy.smallEnd(1) == domain_xy.smallEnd(1)) );

    bool yh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );
         yh_u_dir = ( yh_u_dir && (tbxxy.bigEnd(1) == domain_xy.bigEnd(1)) );

    bool yl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
         yl_w_dir = ( yl_w_dir && (tbxyz.smallEnd(1) == domain_yz.smallEnd(1)) );

    bool yh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir_upwind)   ||
                      (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );
         yh_w_dir = ( yh_w_dir && (tbxyz.bigEnd(1) == domain_yz.bigEnd(1)) );

    // Dirichlet on top or bottom plane
    bool zl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
         zl_u_dir = ( zl_u_dir && (tbxxz.smallEnd(2) == domain_xz.smallEnd(2)) );

    bool zh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );
         zh_u_dir = ( zh_u_dir && (tbxxz.bigEnd(2) == domain_xz.bigEnd(2)) );

    bool zl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
         zl_v_dir = ( zl_v_dir && (tbxyz.smallEnd(2) == domain_yz.smallEnd(2)) );

    bool zh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );
         zh_v_dir = ( zh_v_dir && (tbxyz.bigEnd(2) == domain_yz.bigEnd(2)) );

    //***********************************************************************************
    // X-Dirichlet
    //***********************************************************************************
    if (xl_v_dir) {
        Box planexy = tbxxy; planexy.setBig(0, planexy.smallEnd(0) );
        tbxxy.growLo(0,-1);
        bool need_to_test = (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (!need_to_test || u(dom_lo.x,j,k) >= 0.) {
                tau12(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                     + (-(8./3.) * v(i-1,j,k) + 3. * v(i,j,k) - (1./3.) * v(i+1,j,k))*dxInv[0] );
            } else {
                tau12(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1] +
                                       (v(i, j, k) - v(i-1, j, k))*dxInv[0] );
            }
        });
    }
    if (xh_v_dir) {
        // note: tilebox xy should be nodal, so i|i-1|i-2 at the bigEnd is analogous to i-1|i|i+1 at the smallEnd
        Box planexy = tbxxy; planexy.setSmall(0, planexy.bigEnd(0) );
        tbxxy.growHi(0,-1);
        bool need_to_test = (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (!need_to_test || u(dom_hi.x+1,j,k) <= 0.) {
                tau12(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                     - (-(8./3.) * v(i,j,k) + 3. * v(i-1,j,k) - (1./3.) * v(i-2,j,k))*dxInv[0] );
            } else {
                tau12(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1] +
                                       (v(i, j, k) - v(i-1, j, k))*dxInv[0] );
            }
        });
    }

    if (xl_w_dir) {
        Box planexz = tbxxz; planexz.setBig(0, planexz.smallEnd(0) );
        tbxxz.growLo(0,-1);
        bool need_to_test = (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real du_dz = (u(i, j, k) - u(i, j, k-1))*dxInv[2];
            if (!need_to_test || u(dom_lo.x,j,k) >= 0.) {
                tau13(i,j,k) = 0.5 * ( du_dz
                                     + (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dxInv[0] );
            } else {
                tau13(i,j,k) = 0.5 * ( du_dz
                                     + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );
            }

            if (tau13i) tau13i(i,j,k) = 0.5 * du_dz;
        });
    }
    if (xh_w_dir) {
        // note: tilebox xz should be nodal, so i|i-1|i-2 at the bigEnd is analogous to i-1|i|i+1 at the smallEnd
        Box planexz = tbxxz; planexz.setSmall(0, planexz.bigEnd(0) );
        tbxxz.growHi(0,-1);
        bool need_to_test = (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real du_dz = (u(i, j, k) - u(i, j, k-1))*dxInv[2];
            if (!need_to_test || u(dom_hi.x+1,j,k) <= 0.) {
                tau13(i,j,k) = 0.5 * ( du_dz
                                     - (-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dxInv[0] );
            } else {
                tau13(i,j,k) = 0.5 * ( du_dz
                                     + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );
            }

            if (tau13i) tau13i(i,j,k) = 0.5 * du_dz;
        });
    }

    //***********************************************************************************
    // Y-Dirichlet
    //***********************************************************************************
    if (yl_u_dir) {
        Box planexy = tbxxy; planexy.setBig(1, planexy.smallEnd(1) );
        tbxxy.growLo(1,-1);
        bool need_to_test = (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (!need_to_test || v(i,dom_lo.y,k) >= 0.) {
                tau12(i,j,k) = 0.5 * ( (-(8./3.) * u(i,j-1,k) + 3. * u(i,j,k) - (1./3.) * u(i,j+1,k))*dxInv[1]
                                     + (v(i, j, k) - v(i-1, j, k))*dxInv[0] );
            } else {
                tau12(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                     + (v(i, j, k) - v(i-1, j, k))*dxInv[0] );
            }
        });
    }
    if (yh_u_dir) {
        // note: tilebox xy should be nodal, so j|j-1|j-2 at the bigEnd is analogous to j-1|j|j+1 at the smallEnd
        Box planexy = tbxxy; planexy.setSmall(1, planexy.bigEnd(1) );
        tbxxy.growHi(1,-1);
        bool need_to_test = (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (!need_to_test || v(i,dom_hi.y+1,k) <= 0.) {
                tau12(i,j,k) = 0.5 * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j-1,k) - (1./3.) * u(i,j-2,k))*dxInv[1]
                                      + (v(i, j, k) - v(i-1, j, k))*dxInv[0] );
            } else {
                tau12(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                     + (v(i, j, k) - v(i-1, j, k))*dxInv[0] );
            }
        });
    }

    if (yl_w_dir) {
        Box planeyz = tbxyz; planeyz.setBig(1, planeyz.smallEnd(1) );
        tbxyz.growLo(1,-1);
        bool need_to_test = (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real dv_dz = (v(i, j, k) - v(i, j, k-1))*dxInv[2];
            if (!need_to_test || v(i,dom_lo.y,k) >= 0.) {
                tau23(i,j,k) = 0.5 * ( dv_dz
                                     + (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dxInv[1] );
            } else {
                tau23(i,j,k) = 0.5 * ( dv_dz
                                     + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );
            }

            if (tau23i) tau23i(i,j,k) = 0.5 * dv_dz;
        });
    }
    if (yh_w_dir) {
        // note: tilebox yz should be nodal, so j|j-1|j-2 at the bigEnd is analogous to j-1|j|j+1 at the smallEnd
        Box planeyz = tbxyz; planeyz.setSmall(1, planeyz.bigEnd(1) );
        tbxyz.growHi(1,-1);
        bool need_to_test = (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir_upwind) ? true : false;

        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real dv_dz = (v(i, j, k) - v(i, j, k-1))*dxInv[2];
            if (!need_to_test || v(i,dom_hi.y+1,k) <= 0.) {
                tau23(i,j,k) = 0.5 * ( dv_dz
                                     - (-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dxInv[1] );
            } else {
                tau23(i,j,k) = 0.5 * ( dv_dz
                                     + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );
            }

            if (tau23i) tau23i(i,j,k) = 0.5 * dv_dz;
        });
    }

    //***********************************************************************************
    // Z-Dirichlet
    //***********************************************************************************
    if (zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);

        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real du_dz = (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2];
            tau13(i,j,k) = 0.5 * ( du_dz
                                 + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );

            if (tau13i) tau13i(i,j,k) = 0.5 * du_dz;
        });
    }
    if (zh_u_dir) {
        // note: tilebox xz should be nodal, so k|k-1|k-2 at the bigEnd is analogous to k-1|k|k+1 at the smallEnd
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);

        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real du_dz = -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dxInv[2];
            tau13(i,j,k) = 0.5 * ( du_dz
                                 +  (w(i, j, k) - w(i-1, j, k))*dxInv[0] );

            if (tau13i) tau13i(i,j,k) = 0.5 * du_dz;
        });
    }

    if (zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);

        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real dv_dz = (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dxInv[2];
            tau23(i,j,k) = 0.5 * ( dv_dz
                                 + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );

            if (tau23i) tau23i(i,j,k) = 0.5 * dv_dz;
        });
    }
    if (zh_v_dir) {
        // note: tilebox yz should be nodal, so k|k-1|k-2 at the bigEnd is analogous to k-1|k|k+1 at the smallEnd
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);

        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real dv_dz = -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dxInv[2];
            tau23(i,j,k) = 0.5 * ( dv_dz
                                 +  (w(i, j, k) - w(i, j-1, k))*dxInv[1] );

            if (tau23i) tau23i(i,j,k) = 0.5 * dv_dz;
        });
    }

    // Fill the remaining cells
    //***********************************************************************************
    // Cell centered strains
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

        Real du_dx{0.0};
        if (u_cflag(i+1,j,k).isCovered() && !u_cflag(i,j,k).isCovered()) {
            du_dx = ( 2.0*u(i, j, k) - 3.0*u(i-1, j, k) + u(i-2, j, k))*dxInv[0];
        } else if (u_cflag(i,j,k).isCovered() && !u_cflag(i+1,j,k).isCovered()) {
            du_dx = (- 2.0*u(i+1, j, k) + 3.0*u(i+2, j, k) - u(i+3, j, k))*dxInv[0];
        } else if (!u_cflag(i+1,j,k).isCovered() && !u_cflag(i,j,k).isCovered()) {
            du_dx = (u(i+1, j, k) - u(i, j, k))*dxInv[0];
        }

        Real dv_dy{0.0};
        if (v_cflag(i,j+1,k).isCovered() && !v_cflag(i,j,k).isCovered()) {
            dv_dy = ( 2.0*v(i, j, k) - 3.0*v(i, j-1, k) + v(i, j-2, k))*dxInv[1];
        } else if (v_cflag(i,j,k).isCovered() && !v_cflag(i,j+1,k).isCovered()) {
            dv_dy = (- 2.0*v(i, j+1, k) + 3.0*v(i, j+2, k) - v(i, j+3, k))*dxInv[1];
        } else if (!v_cflag(i,j+1,k).isCovered() && !v_cflag(i,j,k).isCovered()) {
            dv_dy = (v(i, j+1, k) - v(i, j, k))*dxInv[1];
        }

        Real dw_dz{0.0};
        if (w_cflag(i,j,k+1).isCovered() && !w_cflag(i,j,k).isCovered()) {
            dw_dz = ( 2.0*w(i, j, k) - 3.0*w(i, j, k-1) + w(i, j, k-2))*dxInv[2];
        } else if (w_cflag(i,j,k).isCovered() && !w_cflag(i,j,k+1).isCovered()) {
            dw_dz = (- 2.0*w(i, j, k+1) + 3.0*w(i, j, k+2) - w(i, j, k+3))*dxInv[2];
        } else if (!w_cflag(i,j,k+1).isCovered() && !w_cflag(i,j,k).isCovered()) {
            dw_dz = (w(i, j, k+1) - w(i, j, k))*dxInv[2];
        }

        tau11(i,j,k) = du_dx;
        tau22(i,j,k) = dv_dy;
        tau33(i,j,k) = dw_dz;
    });

    // Off-diagonal strains
    ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

        Real du_dy{0.0};
        if (u_cflag(i,j,k).isCovered() && !u_cflag(i,j-1,k).isCovered()) {
            du_dy = ( 2.0*u(i, j-1, k) - 3.0*u(i, j-2, k) + u(i, j-3, k))*dxInv[1];
        } else if (u_cflag(i,j-1,k).isCovered() && !u_cflag(i,j,k).isCovered()) {
            du_dy = (- 2.0*u(i, j, k) + 3.0*u(i, j+1, k) - u(i, j+2, k))*dxInv[1];
        } else if (!u_cflag(i,j,k).isCovered() && !u_cflag(i,j-1,k).isCovered()) {
            du_dy = (u(i, j, k) - u(i, j-1, k))*dxInv[1];
        }

        Real dv_dx{0.0};
        if (v_cflag(i,j,k).isCovered() && !v_cflag(i-1,j,k).isCovered()) {
            dv_dx = ( 2.0*v(i-1, j, k) - 3.0*v(i-2, j, k) + v(i-3, j, k))*dxInv[0];
        } else if (v_cflag(i-1,j,k).isCovered() && !v_cflag(i,j,k).isCovered()) {
            dv_dx = (- 2.0*v(i, j, k) + 3.0*v(i+1, j, k) - v(i+2, j, k))*dxInv[0];
        } else if (!v_cflag(i,j,k).isCovered() && !v_cflag(i-1,j,k).isCovered()) {
            dv_dx = (v(i, j, k) - v(i-1, j, k))*dxInv[0];
        }

        tau12(i,j,k) = 0.5 * ( du_dy + dv_dx );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

        Real du_dz{0.0};
        if (u_cflag(i,j,k).isCovered() && !u_cflag(i,j,k-1).isCovered()) {
            du_dz = ( 2.0*u(i, j, k-1) - 3.0*u(i, j, k-2) + u(i, j, k-3))*dxInv[2];
        } else if (u_cflag(i,j,k-1).isCovered() && !u_cflag(i,j,k).isCovered()) {
            du_dz = (- 2.0*u(i, j, k) + 3.0*u(i, j, k+1) - u(i, j, k+2))*dxInv[2];
        } else if (!u_cflag(i,j,k).isCovered() && !u_cflag(i,j,k-1).isCovered()) {
            du_dz = (u(i, j, k) - u(i, j, k-1))*dxInv[2];
        }

        Real dw_dx{0.0};
        if (w_cflag(i,j,k).isCovered() && !w_cflag(i-1,j,k).isCovered()) {
            dw_dx = ( 2.0*w(i-1, j, k) - 3.0*w(i-2, j, k) + w(i-3, j, k))*dxInv[0];
        } else if (w_cflag(i-1,j,k).isCovered() && !w_cflag(i,j,k).isCovered()) {
            dw_dx = (- 2.0*w(i, j, k) + 3.0*w(i+1, j, k) - w(i+2, j, k))*dxInv[0];
        } else if (!w_cflag(i,j,k).isCovered() && !w_cflag(i-1,j,k).isCovered()) {
            dw_dx = (w(i, j, k) - w(i-1, j, k))*dxInv[0];
        }

        tau13(i,j,k) = 0.5 * ( du_dz + dw_dx );

        if (tau13i) tau13i(i,j,k) = 0.5 * du_dz;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

        Real dv_dz{0.0};
        if (v_cflag(i,j,k).isCovered() && !v_cflag(i,j,k-1).isCovered()) {
            dv_dz = ( 2.0*v(i, j, k-1) - 3.0*v(i, j, k-2) + v(i, j, k-3))*dxInv[2];
        } else if (v_cflag(i,j,k-1).isCovered() && !v_cflag(i,j,k).isCovered()) {
            dv_dz = (- 2.0*v(i, j, k) + 3.0*v(i, j, k+1) - v(i, j, k+2))*dxInv[2];
        } else if (!v_cflag(i,j,k).isCovered() && !v_cflag(i,j,k-1).isCovered()) {
            dv_dz = (v(i, j, k) - v(i, j, k-1))*dxInv[2];
        }

        Real dw_dy{0.0};
        if (w_cflag(i,j,k).isCovered() && !w_cflag(i,j-1,k).isCovered()) {
            dw_dy = ( 2.0*w(i, j-1, k) - 3.0*w(i, j-2, k) + w(i, j-3, k))*dxInv[1];
        } else if (w_cflag(i,j-1,k).isCovered() && !w_cflag(i,j,k).isCovered()) {
            dw_dy = (- 2.0*w(i, j, k) + 3.0*w(i, j+1, k) - w(i, j+2, k))*dxInv[1];
        } else if (!w_cflag(i,j,k).isCovered() && !w_cflag(i,j-1,k).isCovered()) {
            dw_dy = (w(i, j, k) - w(i, j-1, k))*dxInv[1];
        }

        tau23(i,j,k) = 0.5 * ( dv_dz + dw_dy );

        if (tau23i) tau23i(i,j,k) = 0.5 * dv_dz;
    });
}
