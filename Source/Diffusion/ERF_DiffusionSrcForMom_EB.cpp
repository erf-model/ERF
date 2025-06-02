#include <AMReX.H>
#include <AMReX_EB_Slopes_K.H>
#include <ERF_EB.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>

using namespace amrex;

/**
 * Function for computing the momentum RHS for diffusion operator without terrain.
 *
 * @param[in] mfi MultiFab Iterator
 * @param[in] domain computational domain
 * @param[in]  bxx nodal x box for x-mom
 * @param[in]  bxy nodal y box for y-mom
 * @param[in]  bxz nodal z box for z-mom
 * @param[out] rho_u_rhs RHS for x-mom
 * @param[out] rho_v_rhs RHS for y-mom
 * @param[out] rho_w_rhs RHS for z-mom
 * @param[in]  tau11 11 stress
 * @param[in]  tau22 22 stress
 * @param[in]  tau33 33 stress
 * @param[in]  tau12 12 stress
 * @param[in]  tau13 13 stress
 * @param[in]  tau23 23 stress
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_m map factor at cell center
 * @param[in] ebfact EB factories for cell- and face-centered variables
 */
void
DiffusionSrcForMom_EB (const MFIter& mfi,
                    const Box& domain,
                    const Box& bxx, const Box& bxy , const Box& bxz,
                    const Array4<Real>& rho_u_rhs  ,
                    const Array4<Real>& rho_v_rhs  ,
                    const Array4<Real>& rho_w_rhs  ,
                    const Array4<const Real>& u    ,
                    const Array4<const Real>& v    ,
                    const Array4<const Real>& w    ,
                    const Array4<const Real>& tau11,
                    const Array4<const Real>& tau22,
                    const Array4<const Real>& tau33,
                    const Array4<const Real>& tau12,
                    const Array4<const Real>& tau13,
                    const Array4<const Real>& tau23,
                    const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                    const Array4<const Real>& mf_mx,
                    const Array4<const Real>& mf_ux,
                    const Array4<const Real>& mf_vx,
                    const Array4<const Real>& mf_my,
                    const Array4<const Real>& mf_uy,
                    const Array4<const Real>& mf_vy,
                    const eb_& ebfact,
                    const BCRec* d_bcrec_ptr)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_EB()",DiffusionSrcForMom_EB);

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];
    auto dx = 1.0/dxInv[0], dy = 1.0/dxInv[1], dz = 1.0/dxInv[2];
    Real vol = dx * dy * dz;

    bool l_simple = true;

    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

    int n = 0;

    bool extdir_ilo = (d_bcrec_ptr[n].lo(0)==BCType::ext_dir || d_bcrec_ptr[n].lo(0)==BCType::hoextrap);
    bool extdir_ihi = (d_bcrec_ptr[n].hi(0)==BCType::ext_dir || d_bcrec_ptr[n].hi(0)==BCType::hoextrap);
    bool extdir_jlo = (d_bcrec_ptr[n].lo(1)==BCType::ext_dir || d_bcrec_ptr[n].lo(1)==BCType::hoextrap);
    bool extdir_jhi = (d_bcrec_ptr[n].hi(1)==BCType::ext_dir || d_bcrec_ptr[n].hi(1)==BCType::hoextrap);
    bool extdir_klo = (d_bcrec_ptr[n].lo(2)==BCType::ext_dir || d_bcrec_ptr[n].lo(2)==BCType::hoextrap);
    bool extdir_khi = (d_bcrec_ptr[n].hi(2)==BCType::ext_dir || d_bcrec_ptr[n].hi(2)==BCType::hoextrap);

    // EB u-factory
    Array4<const EBCellFlag> u_cellflg = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > u_volfrac = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > u_volcent = (ebfact.get_u_const_factory())->getCentroid().const_array(mfi);
    Array4<const Real      > u_afrac_x = (ebfact.get_u_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > u_afrac_y = (ebfact.get_u_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > u_afrac_z = (ebfact.get_u_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > u_bcent = (ebfact.get_u_const_factory())->getBndryCent().const_array(mfi);
    Array4<const Real      > u_fcx   = (ebfact.get_u_const_factory())->getFaceCent()[0]->const_array(mfi);
    Array4<const Real      > u_fcy   = (ebfact.get_u_const_factory())->getFaceCent()[1]->const_array(mfi);
    Array4<const Real      > u_fcz   = (ebfact.get_u_const_factory())->getFaceCent()[2]->const_array(mfi);

    // EB v-factory
    Array4<const EBCellFlag> v_cellflg = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > v_volfrac = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > v_volcent = (ebfact.get_v_const_factory())->getCentroid().const_array(mfi);
    Array4<const Real      > v_afrac_x = (ebfact.get_v_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > v_afrac_y = (ebfact.get_v_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > v_afrac_z = (ebfact.get_v_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > v_bcent = (ebfact.get_v_const_factory())->getBndryCent().const_array(mfi);
    Array4<const Real      > v_fcx   = (ebfact.get_v_const_factory())->getFaceCent()[0]->const_array(mfi);
    Array4<const Real      > v_fcy   = (ebfact.get_v_const_factory())->getFaceCent()[1]->const_array(mfi);
    Array4<const Real      > v_fcz   = (ebfact.get_v_const_factory())->getFaceCent()[2]->const_array(mfi);

    // EB w-factory
    Array4<const EBCellFlag> w_cellflg = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > w_volfrac = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > w_volcent = (ebfact.get_w_const_factory())->getCentroid().const_array(mfi);
    Array4<const Real      > w_afrac_x = (ebfact.get_w_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > w_afrac_y = (ebfact.get_w_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > w_afrac_z = (ebfact.get_w_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > w_bcent = (ebfact.get_w_const_factory())->getBndryCent().const_array(mfi);
    Array4<const Real      > w_fcx   = (ebfact.get_w_const_factory())->getFaceCent()[0]->const_array(mfi);
    Array4<const Real      > w_fcy   = (ebfact.get_w_const_factory())->getFaceCent()[1]->const_array(mfi);
    Array4<const Real      > w_fcz   = (ebfact.get_w_const_factory())->getFaceCent()[2]->const_array(mfi);

    // Two methods are implemented:
    // (1) Simple method: approximate the gradient at EB
    // (2) Compute the gradient at EB using Least-Squares fitting

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (u_volfrac(i,j,k)>0.) {
            // Inv Jacobian
            Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);
            // Area corrections
            Real Imfy_hi = 1. / mf_my(i  ,j,0);
            Real Imfy_lo = 1. / mf_my(i-1,j,0);
            Real Imfx_hi = 1. / (0.5 * (mf_vx(i,j+1,0) + mf_vx(i-1,j+1,0)));
            Real Imfx_lo = 1. / (0.5 * (mf_vx(i,j  ,0) + mf_vx(i-1,j  ,0)));
            rho_u_rhs(i,j,k) -= ( (tau11(i  , j  , k  ) * Imfy_hi * u_afrac_x(i+1,j  ,k  ) 
                                -  tau11(i-1, j  , k  ) * Imfy_lo * u_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                + (tau12(i  , j+1, k  ) * Imfx_hi * u_afrac_y(i  ,j+1,k  ) 
                                -  tau12(i  , j  , k  ) * Imfx_lo * u_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                + (tau13(i  , j  , k+1) * u_afrac_z(i  ,j  ,k+1) 
                                -  tau13(i  , j  , k  ) * u_afrac_z(i  ,j  ,k  ) ) * dzinv ) / u_volfrac(i,j,k);
            // Boundary flux (simple version)
            if (u_cellflg(i,j,k).isSingleValued()) {

                if (l_simple) {

                    Real axm = u_afrac_x(i  ,j  ,k  );
                    Real axp = u_afrac_x(i+1,j  ,k  );
                    Real aym = u_afrac_y(i  ,j  ,k  );
                    Real ayp = u_afrac_y(i  ,j+1,k  );
                    Real azm = u_afrac_z(i  ,j  ,k  );
                    Real azp = u_afrac_z(i  ,j  ,k+1);
    
                    Real adx = (axm-axp) * dy * dz;
                    Real ady = (aym-ayp) * dx * dz;
                    Real adz = (azm-azp) * dx * dy;
    
                    Real barea = std::sqrt(adx*adx + ady*ady + adz*adz);
    
                    Real dist_x = u_bcent(i,j,k,0)*dx;
                    Real dist_y = u_bcent(i,j,k,1)*dy;
                    Real dist_z = u_bcent(i,j,k,2)*dz;
                    Real dn = std::sqrt( dist_x * dist_x + dist_y * dist_y + dist_z * dist_z );

                    rho_u_rhs(i,j,k) -= - barea / vol * u(i,j,k) / dn / u_volfrac(i,j,k);

                } else {

                    // tau11, tau12, tau13

                    // GpuArray<Real,AMREX_SPACEDIM> slopes_eb;

                    // slopes_eb = amrex_calc_slopes_extdir_eb(
                    //     i, j, k, n, u, u_volcent, u_volfrac,
                    //     AMREX_D_DECL(u_fcx, u_fcy, u_fcz), u_cellflg,
                    //     AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                    //     AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                    //     AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                    //     AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                    //     2);

                    

                }
            }
        }

    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (v_volfrac(i,j,k)>0.) {
            // Inv Jacobian
            Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);
            // Area corrections
            Real Imfy_hi = 1. / (0.5 * (mf_uy(i+1,j,0) + mf_uy(i+1,j-1,0)));
            Real Imfy_lo = 1. / (0.5 * (mf_uy(i  ,j,0) + mf_uy(i  ,j-1,0)));
            Real Imfx_hi = 1. / mf_mx(i  ,j,0);
            Real Imfx_lo = 1. / mf_mx(i-1,j,0);
            rho_v_rhs(i,j,k) -= ( (tau12(i+1, j  , k  ) * Imfy_hi * v_afrac_x(i+1,j  ,k  )
                                -  tau12(i  , j  , k  ) * Imfy_lo * v_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                + (tau22(i  , j  , k  ) * Imfx_hi * v_afrac_y(i  ,j+1,k  )
                                -  tau22(i  , j-1, k  ) * Imfx_lo * v_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                + (tau23(i  , j  , k+1) * v_afrac_z(i  ,j  ,k+1)
                                -  tau23(i  , j  , k  ) * v_afrac_z(i  ,j  ,k  ) ) * dzinv ) / v_volfrac(i,j,k);
            // Boundary flux (simple version)
            if (v_cellflg(i,j,k).isSingleValued()) {

                if (l_simple) {

                    Real axm = v_afrac_x(i  ,j  ,k  );
                    Real axp = v_afrac_x(i+1,j  ,k  );
                    Real aym = v_afrac_y(i  ,j  ,k  );
                    Real ayp = v_afrac_y(i  ,j+1,k  );
                    Real azm = v_afrac_z(i  ,j  ,k  );
                    Real azp = v_afrac_z(i  ,j  ,k+1);
    
                    Real adx = (axm-axp) * dy * dz;
                    Real ady = (aym-ayp) * dx * dz;
                    Real adz = (azm-azp) * dx * dy;
    
                    Real barea = std::sqrt(adx*adx + ady*ady + adz*adz);
    
                    Real dist_x = v_bcent(i,j,k,0)*dx;
                    Real dist_y = v_bcent(i,j,k,1)*dy;
                    Real dist_z = v_bcent(i,j,k,2)*dz;
                    Real dn = std::sqrt( dist_x * dist_x + dist_y * dist_y + dist_z * dist_z );

                    rho_v_rhs(i,j,k) -= - barea / vol * v(i,j,k) / dn / v_volfrac(i,j,k);

                } else {

                    // tau12
                    // tau22
                    // tau23


                }
            }
        }                           
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (w_volfrac(i,j,k)>0.) {
            // Inv Jacobian
            Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
            // Area corrections
            Real Imfy_hi = 1. / mf_uy(i+1,j  ,0);
            Real Imfy_lo = 1. / mf_uy(i  ,j  ,0);
            Real Imfx_hi = 1. / mf_vx(i  ,j+1,0);
            Real Imfx_lo = 1. / mf_vx(i  ,j  ,0);
            rho_w_rhs(i,j,k) -= ( (tau13(i+1, j  , k  ) * Imfy_hi * w_afrac_x(i+1,j  ,k  )
                                -  tau13(i  , j  , k  ) * Imfy_lo * w_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                + (tau23(i  , j+1, k  ) * Imfx_hi * w_afrac_y(i  ,j+1,k  )
                                -  tau23(i  , j  , k  ) * Imfx_lo * w_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                + (tau33(i  , j  , k  ) * w_afrac_z(i  ,j  ,k+1)        
                                -  tau33(i  , j  , k-1) * w_afrac_z(i  ,j  ,k  ) ) * dzinv ) / w_volfrac(i,j,k);
            // Boundary flux (simple version)
            if (w_cellflg(i,j,k).isSingleValued()) {

                Real axm = w_afrac_x(i  ,j  ,k  );
                Real axp = w_afrac_x(i+1,j  ,k  );
                Real aym = w_afrac_y(i  ,j  ,k  );
                Real ayp = w_afrac_y(i  ,j+1,k  );
                Real azm = w_afrac_z(i  ,j  ,k  );
                Real azp = w_afrac_z(i  ,j  ,k+1);
 
                Real adx = (axm-axp) * dy * dz;
                Real ady = (aym-ayp) * dx * dz;
                Real adz = (azm-azp) * dx * dy;
 
                Real barea = std::sqrt(adx*adx + ady*ady + adz*adz);
 
                Real dist_x = w_bcent(i,j,k,0)*dx;
                Real dist_y = w_bcent(i,j,k,1)*dy;
                Real dist_z = w_bcent(i,j,k,2)*dz;
                Real dn = std::sqrt( dist_x * dist_x + dist_y * dist_y + dist_z * dist_z );

                rho_w_rhs(i,j,k) -= - barea / vol * w(i,j,k) / dn / w_volfrac(i,j,k);
            }
        }
    });

}
