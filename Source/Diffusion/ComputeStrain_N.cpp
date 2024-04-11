#include <Diffusion.H>

using namespace amrex;

/**
 * Function for computing the strain rates without terrain.
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
 * @param[in] mf_m map factor at cell center
 * @param[in] mf_u map factor at x-face
 * @param[in] mf_v map factor at y-face
 */
void
ComputeStrain_N (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Box domain,
                 const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                 Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                 Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                 const BCRec* bc_ptr, const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                 const Array4<const Real>& mf_m, const Array4<const Real>& mf_u, const Array4<const Real>& mf_v)
{
    // Conver domain to each index type to test if we are on dirichlet boundary
    Box domain_xy = convert(domain, tbxxy.ixType());
    Box domain_xz = convert(domain, tbxxz.ixType());
    Box domain_yz = convert(domain, tbxyz.ixType());

    // Dirichlet on left or right plane
    bool xl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
         xl_v_dir = ( xl_v_dir && (tbxxy.smallEnd(0) == domain_xy.smallEnd(0)) );

    bool xh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );
         xh_v_dir = ( xh_v_dir && (tbxxy.bigEnd(0) == domain_xy.bigEnd(0)) );

    bool xl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
         xl_w_dir = ( xl_w_dir && (tbxxz.smallEnd(0) == domain_xz.smallEnd(0)) );

    bool xh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );
         xh_w_dir = ( xh_w_dir && (tbxxz.bigEnd(0) == domain_xz.bigEnd(0)) );

    // Dirichlet on front or back plane
    bool yl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
         yl_u_dir = ( yl_u_dir && (tbxxy.smallEnd(1) == domain_xy.smallEnd(1)) );

    bool yh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );
         yh_u_dir = ( yh_u_dir && (tbxxy.bigEnd(1) == domain_xy.bigEnd(1)) );

    bool yl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
         yl_w_dir = ( yl_w_dir && (tbxyz.smallEnd(1) == domain_yz.smallEnd(1)) );

    bool yh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir)          ||
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

    // X-Dirichlet
    //***********************************************************************************
    if (xl_v_dir) {
        Box planexy = tbxxy; planexy.setBig(0, planexy.smallEnd(0) );
        tbxxy.growLo(0,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i, j-1, k)/mf_u(i,j-1,0))*dxInv[1] +
                                   (-(8./3.) * v(i-1,j,k)/mf_v(i-1,j,0) + 3. * v(i,j,k)/mf_v(i,j,0) - (1./3.) * v(i+1,j,k)/mf_v(i+1,j,0))*dxInv[0] ) * mf_u(i,j,0)*mf_u(i,j,0);
        });
    }
    if (xh_v_dir) {
        // note: tilebox xy should be nodal, so i|i-1|i-2 at the bigEnd is analogous to i-1|i|i+1 at the smallEnd
        Box planexy = tbxxy; planexy.setSmall(0, planexy.bigEnd(0) );
        tbxxy.growHi(0,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i, j-1, k)/mf_u(i,j-1,0))*dxInv[1] +
                                  -(-(8./3.) * v(i,j,k)/mf_v(i,j,0) + 3. * v(i-1,j,k)/mf_v(i-1,j,0) - (1./3.) * v(i-2,j,k)/mf_v(i-2,j,0))*dxInv[0] ) * mf_u(i,j,0)*mf_u(i,j,0);
        });
    }

    if (xl_w_dir) {
        Box planexz = tbxxz; planexz.setBig(0, planexz.smallEnd(0) );
        tbxxz.growLo(0,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2] +
                                   (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dxInv[0]*mf_u(i,j,0) );
        });
    }
    if (xh_w_dir) {
        // note: tilebox xz should be nodal, so i|i-1|i-2 at the bigEnd is analogous to i-1|i|i+1 at the smallEnd
        Box planexz = tbxxz; planexz.setSmall(0, planexz.bigEnd(0) );
        tbxxz.growHi(0,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2] +
                                  -(-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dxInv[0]*mf_u(i,j,0) );
        });
    }

    // Y-Dirichlet
    //***********************************************************************************
    if (yl_u_dir) {
        Box planexy = tbxxy; planexy.setBig(1, planexy.smallEnd(1) );
        tbxxy.growLo(1,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau12(i,j,k) = 0.5 * ( (-(8./3.) * u(i,j-1,k)/mf_u(i,j-1,0) + 3. * u(i,j,k)/mf_u(i,j,0) - (1./3.) * u(i,j+1,k)/mf_m(i,j+1,0))*dxInv[1] +
                                   (v(i, j, k)/mf_v(i,j,0) - v(i-1, j, k)/mf_v(i-1,j,0))*dxInv[0] ) * mf_u(i,j,0)*mf_u(i,j,0);
        });
    }
    if (yh_u_dir) {
        // note: tilebox xy should be nodal, so j|j-1|j-2 at the bigEnd is analogous to j-1|j|j+1 at the smallEnd
        Box planexy = tbxxy; planexy.setSmall(1, planexy.bigEnd(1) );
        tbxxy.growHi(1,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau12(i,j,k) = 0.5 * ( -(-(8./3.) * u(i,j,k)/mf_u(i,j,0) + 3. * u(i,j-1,k)/mf_u(i,j-1,0) - (1./3.) * u(i,j-2,k)/mf_u(i,j-2,0))*dxInv[1] +
                                   (v(i, j, k)/mf_v(i,j,0) - v(i-1, j, k)/mf_v(i-1,j,0))*dxInv[0] ) * mf_u(i,j,0)*mf_u(i,j,0);
        });
    }

    if (yl_w_dir) {
        Box planeyz = tbxyz; planeyz.setBig(1, planeyz.smallEnd(1) );
        tbxyz.growLo(1,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2] +
                                   (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dxInv[1]*mf_v(i,j,0)*mf_v(i,j,0) );
        });
    }
    if (yh_w_dir) {
        // note: tilebox yz should be nodal, so j|j-1|j-2 at the bigEnd is analogous to j-1|j|j+1 at the smallEnd
        Box planeyz = tbxyz; planeyz.setSmall(1, planeyz.bigEnd(1) );
        tbxyz.growHi(1,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2] +
                                   -(-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dxInv[1]*mf_v(i,j,0)*mf_v(i,j,0) );
        });
    }

    // Z-Dirichlet
    //***********************************************************************************
    if (zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau13(i,j,k) = 0.5 * ( (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2] +
                                   (w(i, j, k) - w(i-1, j, k))*dxInv[0]*mf_u(i,j,0) );
        });
    }
    if (zh_u_dir) {
        // note: tilebox xz should be nodal, so k|k-1|k-2 at the bigEnd is analogous to k-1|k|k+1 at the smallEnd
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau13(i,j,k) = 0.5 * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dxInv[2] +
                                   (w(i, j, k) - w(i-1, j, k))*dxInv[0]*mf_u(i,j,0) );
        });
    }

    if (zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau23(i,j,k) = 0.5 * ( (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dxInv[2] +
                                   (w(i, j, k) - w(i, j-1, k))*dxInv[1]*mf_v(i,j,0) );
        });
    }
    if (zh_v_dir) {
        // note: tilebox yz should be nodal, so k|k-1|k-2 at the bigEnd is analogous to k-1|k|k+1 at the smallEnd
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau23(i,j,k) = 0.5 * ( -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dxInv[2] +
                                   (w(i, j, k) - w(i, j-1, k))*dxInv[1]*mf_v(i,j,0) );
        });
    }

    // Fill the remaining cells
    //***********************************************************************************
    // Cell centered strains
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau11(i,j,k) = (u(i+1, j  , k  )/mf_u(i+1,j,0) - u(i, j, k)/mf_u(i,j,0))*dxInv[0]*mf_u(i,j,0)*mf_u(i,j,0);
        tau22(i,j,k) = (v(i  , j+1, k  )/mf_v(i,j+1,0) - v(i, j, k)/mf_v(i,j,0))*dxInv[1]*mf_v(i,j,0)*mf_v(i,j,0);
        tau33(i,j,k) = (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
    });

    // Off-diagonal strains
    ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i, j-1, k)/mf_u(i,j-1,0))*dxInv[1] +
                               (v(i, j, k)/mf_v(i,j,0) - v(i-1, j, k)/mf_v(i-1,j,0))*dxInv[0] ) * mf_u(i,j,0)*mf_u(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2] + (w(i, j, k) - w(i-1, j, k))*dxInv[0]*mf_u(i,j,0) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2] + (w(i, j, k) - w(i, j-1, k))*dxInv[1]*mf_v(i,j,0) );
    });

}
