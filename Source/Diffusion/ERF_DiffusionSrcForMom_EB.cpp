#include <AMReX.H>
#include <AMReX_EB_Slopes_K.H>
#include <ERF_EB.H>
#include <ERF_Diffusion.H>
#include <ERF_IndexDefines.H>
#include <ERF_EBSlopes.H>
#include <ERF_DiffStruct.H>
#include <ERF_EBStruct.H>

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
                    [[maybe_unused]] const Box& domain,
                    const Box& bxx, const Box& bxy , const Box& bxz,
                    const Array4<Real>& rho_u_rhs  ,
                    const Array4<Real>& rho_v_rhs  ,
                    const Array4<Real>& rho_w_rhs  ,
                    const Array4<const Real>& u_arr,
                    const Array4<const Real>& v_arr,
                    const Array4<const Real>& w_arr,
                    const Array4<const Real>& tau11,
                    const Array4<const Real>& tau22,
                    const Array4<const Real>& tau33,
                    const Array4<const Real>& tau12,
                    const Array4<const Real>& tau13,
                    const Array4<const Real>& tau23,
                    const Real* dx_arr,
                    const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                    const Array4<const Real>& mf_mx,
                    const Array4<const Real>& mf_ux,
                    const Array4<const Real>& mf_vx,
                    const Array4<const Real>& mf_my,
                    const Array4<const Real>& mf_uy,
                    const Array4<const Real>& mf_vy,
                    const SolverChoice& solverChoice,
                    const eb_& ebfact,
                   [[maybe_unused]] const BCRec* d_bcrec_ptr)
{
    BL_PROFILE_VAR("DiffusionSrcForMom_EB()",DiffusionSrcForMom_EB);

    DiffChoice dc = solverChoice.diffChoice;
    const bool l_use_constAlpha = ( dc.molec_diff_type == MolecDiffType::ConstantAlpha );
    Real mu_eff = (l_use_constAlpha) ? 2.0 * dc.dynamic_viscosity / dc.rho0_trans
                                     : 2.0 * dc.dynamic_viscosity;

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];
    Real dx = dx_arr[0], dy = dx_arr[1], dz = dx_arr[2];
    Real vol = dx * dy * dz;

    EBChoice ebChoice = solverChoice.ebChoice;
    const bool l_no_slip = (ebChoice.eb_boundary_type == EBBoundaryType::NoSlipWall);

    const bool l_constraint_x = solverChoice.diffChoice.eb_diff_constraint_x;
    const bool l_constraint_y = solverChoice.diffChoice.eb_diff_constraint_y;
    const bool l_constraint_z = solverChoice.diffChoice.eb_diff_constraint_z;

    // int n = 0;

    // EB u-factory
    Array4<const EBCellFlag> u_cellflg = (ebfact.get_u_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > u_volfrac = (ebfact.get_u_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > u_volcent = (ebfact.get_u_const_factory())->getCentroid().const_array(mfi);
    Array4<const Real      > u_afrac_x = (ebfact.get_u_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > u_afrac_y = (ebfact.get_u_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > u_afrac_z = (ebfact.get_u_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > u_bcent = (ebfact.get_u_const_factory())->getBndryCent().const_array(mfi);
    Array4<const Real      > u_bnorm = (ebfact.get_u_const_factory())->getBndryNorm().const_array(mfi);

    // EB v-factory
    Array4<const EBCellFlag> v_cellflg = (ebfact.get_v_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > v_volfrac = (ebfact.get_v_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > v_volcent = (ebfact.get_v_const_factory())->getCentroid().const_array(mfi);
    Array4<const Real      > v_afrac_x = (ebfact.get_v_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > v_afrac_y = (ebfact.get_v_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > v_afrac_z = (ebfact.get_v_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > v_bcent = (ebfact.get_v_const_factory())->getBndryCent().const_array(mfi);
    Array4<const Real      > v_bnorm = (ebfact.get_v_const_factory())->getBndryNorm().const_array(mfi);

    // EB w-factory
    Array4<const EBCellFlag> w_cellflg = (ebfact.get_w_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > w_volfrac = (ebfact.get_w_const_factory())->getVolFrac().const_array(mfi);
    Array4<const Real      > w_volcent = (ebfact.get_w_const_factory())->getCentroid().const_array(mfi);
    Array4<const Real      > w_afrac_x = (ebfact.get_w_const_factory())->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > w_afrac_y = (ebfact.get_w_const_factory())->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > w_afrac_z = (ebfact.get_w_const_factory())->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > w_bcent = (ebfact.get_w_const_factory())->getBndryCent().const_array(mfi);
    Array4<const Real      > w_bnorm = (ebfact.get_w_const_factory())->getBndryNorm().const_array(mfi);

    ParallelFor(bxx, bxy, bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (u_volfrac(i,j,k)>0.) {

            // Inv Jacobian
            Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

            Real diffContrib  = ( (tau11(i  , j  , k  ) * u_afrac_x(i+1,j  ,k  )
                                -  tau11(i-1, j  , k  ) * u_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                + (tau12(i  , j+1, k  ) * u_afrac_y(i  ,j+1,k  )
                                -  tau12(i  , j  , k  ) * u_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                + (tau13(i  , j  , k+1) * u_afrac_z(i  ,j  ,k+1)
                                -  tau13(i  , j  , k  ) * u_afrac_z(i  ,j  ,k  )) * dzinv );
            diffContrib      /= u_volfrac(i,j,k);

            rho_u_rhs(i,j,k) -= diffContrib;

            if (!l_constraint_x && l_no_slip && u_cellflg(i,j,k).isSingleValued()) {

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

                const RealVect bcent_eb {u_bcent(i,j,k,0), u_bcent(i,j,k,1), u_bcent(i,j,k,2)};
                const Real Dirichlet_u {0.};
                const Real Dirichlet_v {0.};
                const Real Dirichlet_w {0.};

                GpuArray<Real,AMREX_SPACEDIM> slopes_u;
                GpuArray<Real,AMREX_SPACEDIM> slopes_v;
                GpuArray<Real,AMREX_SPACEDIM> slopes_w;

                slopes_u = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
                slopes_v = erf_calc_slopes_eb_Dirichlet_staggered( Vars::xvel, Vars::yvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
                slopes_w = erf_calc_slopes_eb_Dirichlet_staggered( Vars::xvel, Vars::zvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);

                Real dudx = slopes_u[0];
                Real dudy = slopes_u[1];
                Real dudz = slopes_u[2];
                Real dvdx = slopes_v[0];
                Real dvdy = slopes_v[1];
                Real dwdx = slopes_w[0];
                Real dwdz = slopes_w[2];

                Real tau11_eb = ( dudx - ( dudx + dvdy + dwdz ) / 3. );
                Real tau12_eb = 0.5 * (dudy + dvdx);
                Real tau13_eb = 0.5 * (dudz + dwdx);

                Real dudn = -(u_bnorm(i,j,k,0) * tau11_eb + u_bnorm(i,j,k,1) * tau12_eb + u_bnorm(i,j,k,2) * tau13_eb);

                rho_u_rhs(i,j,k) -= mu_eff * barea * dudn / (vol * u_volfrac(i,j,k));
            }
        }

    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (v_volfrac(i,j,k)>0.) {

            // Inv Jacobian
            Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

            Real diffContrib  = ( (tau12(i+1, j  , k  ) * v_afrac_x(i+1,j  ,k  )
                                -  tau12(i  , j  , k  ) * v_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                + (tau22(i  , j  , k  ) * v_afrac_y(i  ,j+1,k  )
                                -  tau22(i  , j-1, k  ) * v_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                + (tau23(i  , j  , k+1) * v_afrac_z(i  ,j  ,k+1)
                                -  tau23(i  , j  , k  ) * v_afrac_z(i  ,j  ,k  ) ) * dzinv );
            diffContrib      /= v_volfrac(i,j,k);

            rho_v_rhs(i,j,k) -= diffContrib;

            if (!l_constraint_y && l_no_slip && v_cellflg(i,j,k).isSingleValued()) {

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

                const RealVect bcent_eb {v_bcent(i,j,k,0), v_bcent(i,j,k,1), v_bcent(i,j,k,2)};
                const Real Dirichlet_u {0.};
                const Real Dirichlet_v {0.};
                const Real Dirichlet_w {0.};

                GpuArray<Real,AMREX_SPACEDIM> slopes_u;
                GpuArray<Real,AMREX_SPACEDIM> slopes_v;
                GpuArray<Real,AMREX_SPACEDIM> slopes_w;

                slopes_u = erf_calc_slopes_eb_Dirichlet_staggered( Vars::yvel, Vars::xvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
                slopes_v = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
                slopes_w = erf_calc_slopes_eb_Dirichlet_staggered( Vars::yvel, Vars::zvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);

                Real dudx = slopes_u[0];
                Real dudy = slopes_u[1];
                Real dvdx = slopes_v[0];
                Real dvdy = slopes_v[1];
                Real dvdz = slopes_v[2];
                Real dwdy = slopes_w[1];
                Real dwdz = slopes_w[2];

                Real tau22_eb = ( dvdy - ( dudx + dvdy + dwdz ) / 3. );
                Real tau12_eb = 0.5 * (dudy + dvdx);
                Real tau23_eb = 0.5 * (dvdz + dwdy);

                Real dvdn = -(v_bnorm(i,j,k,0) * tau12_eb + v_bnorm(i,j,k,1) * tau22_eb + v_bnorm(i,j,k,2) * tau23_eb);

                rho_v_rhs(i,j,k) -= mu_eff * barea * dvdn / (vol * v_volfrac(i,j,k));
            }
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (w_volfrac(i,j,k)>0.) {

            // Inv Jacobian
            Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

            Real diffContrib  = ( (tau13(i+1, j  , k  ) * w_afrac_x(i+1,j  ,k  )
                                -  tau13(i  , j  , k  ) * w_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                + (tau23(i  , j+1, k  ) * w_afrac_y(i  ,j+1,k  )
                                -  tau23(i  , j  , k  ) * w_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                + (tau33(i  , j  , k  ) * w_afrac_z(i  ,j  ,k+1)
                                -  tau33(i  , j  , k-1) * w_afrac_z(i  ,j  ,k  ) ) * dzinv );
            diffContrib      /= w_volfrac(i,j,k);
            rho_w_rhs(i,j,k) -= diffContrib;

            if (!l_constraint_z && l_no_slip && w_cellflg(i,j,k).isSingleValued()) {

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

                const RealVect bcent_eb {w_bcent(i,j,k,0), w_bcent(i,j,k,1), w_bcent(i,j,k,2)};
                const Real Dirichlet_u {0.};
                const Real Dirichlet_v {0.};
                const Real Dirichlet_w {0.};

                GpuArray<Real,AMREX_SPACEDIM> slopes_u;
                GpuArray<Real,AMREX_SPACEDIM> slopes_v;
                GpuArray<Real,AMREX_SPACEDIM> slopes_w;

                slopes_u = erf_calc_slopes_eb_Dirichlet_staggered( Vars::zvel, Vars::xvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
                slopes_v = erf_calc_slopes_eb_Dirichlet_staggered( Vars::zvel, Vars::yvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
                slopes_w = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);

                Real dudx = slopes_u[0];
                Real dudz = slopes_u[2];
                Real dvdy = slopes_v[1];
                Real dvdz = slopes_v[2];
                Real dwdx = slopes_w[0];
                Real dwdy = slopes_w[1];
                Real dwdz = slopes_w[2];

                Real tau33_eb = ( dwdz - ( dudx + dvdy + dwdz ) / 3. );
                Real tau13_eb = 0.5 * (dudz + dwdx);
                Real tau23_eb = 0.5 * (dvdz + dwdy);

                Real dwdn = -(w_bnorm(i,j,k,0) * tau13_eb + w_bnorm(i,j,k,1) * tau23_eb + w_bnorm(i,j,k,2) * tau33_eb);

                rho_w_rhs(i,j,k) -= mu_eff * barea * dwdn / (vol * w_volfrac(i,j,k));
            }
        }
    });

}
