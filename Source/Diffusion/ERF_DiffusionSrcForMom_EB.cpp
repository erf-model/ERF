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
 * Compute orthonormal tangent vectors at an EB boundary given the normal vector.
 * Uses Gram-Schmidt orthogonalization against standard basis vectors.
 *
 * @param[in]  nx, ny, nz  Components of the normal vector
 * @param[out] tbx_x, tbx_y, tbx_z  Components of first tangent vector (from e_x)
 * @param[out] tby_x, tby_y, tby_z  Components of second tangent vector (from e_y)
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_tangent_vectors (Real nx, Real ny, Real nz,
                               Real& tbx_x, Real& tbx_y, Real& tbx_z,
                               Real& tby_x, Real& tby_y, Real& tby_z)
{
    // x-tangential vector: t_bx = (e_x - (e_x · n)n) / ||e_x - (e_x · n)n||
    // e_x = (1,0,0), so e_x · n = nx
    tbx_x = one - nx * nx;
    tbx_y = - nx * ny;
    tbx_z = - nx * nz;
    Real tbx_norm = std::sqrt(tbx_x*tbx_x + tbx_y*tbx_y + tbx_z*tbx_z);
    tbx_x /= tbx_norm;
    tbx_y /= tbx_norm;
    tbx_z /= tbx_norm;

    // y-tangential vector: t_by = (e_y - (e_y · n)n) / ||e_y - (e_y · n)n||
    // e_y = (0,1,0), so e_y · n = ny
    tby_x = - ny * nx;
    tby_y = one - ny * ny;
    tby_z = - ny * nz;
    Real tby_norm = std::sqrt(tby_x*tby_x + tby_y*tby_y + tby_z*tby_z);
    tby_x /= tby_norm;
    tby_y /= tby_norm;
    tby_z /= tby_norm;
}

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
                    const Array4<const Real>& u_tau_eb13,
                    const Array4<const Real>& u_tau_eb23,
                    const Array4<const Real>& v_tau_eb13,
                    const Array4<const Real>& v_tau_eb23,
                    const Array4<const Real>& w_tau_eb13,
                    const Array4<const Real>& w_tau_eb23,
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
    Real mu_eff = (l_use_constAlpha) ? two * dc.dynamic_viscosity / dc.rho0_trans
                                     : two * dc.dynamic_viscosity;

    auto dxinv = dxInv[0], dyinv = dxInv[1], dzinv = dxInv[2];
    Real dx = dx_arr[0], dy = dx_arr[1], dz = dx_arr[2];
    Real vol = dx * dy * dz;

    EBChoice ebChoice = solverChoice.ebChoice;
    const bool l_no_slip = (ebChoice.eb_boundary_type == EBBoundaryType::NoSlipWall);
    const bool l_surface_layer = (ebChoice.eb_boundary_type == EBBoundaryType::SurfaceLayer);
    const bool l_constraint_x = solverChoice.diffChoice.eb_diff_constraint_x;
    const bool l_constraint_y = solverChoice.diffChoice.eb_diff_constraint_y;
    const bool l_constraint_z = solverChoice.diffChoice.eb_diff_constraint_z;

    // int n = 0;

    // EB u-factory
    const auto* u_factory = ebfact.get_u_const_factory();
    Array4<const EBCellFlag> u_cellflg = u_factory->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > u_volfrac = u_factory->getVolFrac().const_array(mfi);
    Array4<const Real      > u_volcent = u_factory->getCentroid().const_array(mfi);
    Array4<const Real      > u_afrac_x = u_factory->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > u_afrac_y = u_factory->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > u_afrac_z = u_factory->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > u_bcent = u_factory->getBndryCent().const_array(mfi);
    Array4<const Real      > u_bnorm = u_factory->getBndryNorm().const_array(mfi);

    // EB v-factory
    const auto* v_factory = ebfact.get_v_const_factory();
    Array4<const EBCellFlag> v_cellflg = v_factory->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > v_volfrac = v_factory->getVolFrac().const_array(mfi);
    Array4<const Real      > v_volcent = v_factory->getCentroid().const_array(mfi);
    Array4<const Real      > v_afrac_x = v_factory->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > v_afrac_y = v_factory->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > v_afrac_z = v_factory->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > v_bcent = v_factory->getBndryCent().const_array(mfi);
    Array4<const Real      > v_bnorm = v_factory->getBndryNorm().const_array(mfi);

    // EB w-factory
    const auto* w_factory = ebfact.get_w_const_factory();
    Array4<const EBCellFlag> w_cellflg = w_factory->getMultiEBCellFlagFab()[mfi].const_array();
    Array4<const Real      > w_volfrac = w_factory->getVolFrac().const_array(mfi);
    Array4<const Real      > w_volcent = w_factory->getCentroid().const_array(mfi);
    Array4<const Real      > w_afrac_x = w_factory->getAreaFrac()[0]->const_array(mfi);
    Array4<const Real      > w_afrac_y = w_factory->getAreaFrac()[1]->const_array(mfi);
    Array4<const Real      > w_afrac_z = w_factory->getAreaFrac()[2]->const_array(mfi);
    Array4<const Real      > w_bcent = w_factory->getBndryCent().const_array(mfi);
    Array4<const Real      > w_bnorm = w_factory->getBndryNorm().const_array(mfi);

    // x-momentum
    ParallelFor(bxx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (u_volfrac(i,j,k)>zero) {

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

            if (!l_constraint_x && u_cellflg(i,j,k).isSingleValued()) {

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

                Real dudn = zero;

                if (l_no_slip || l_surface_layer) {

                    RealVect bcent_eb {u_bcent(i,j,k,0), u_bcent(i,j,k,1), u_bcent(i,j,k,2)};

                    Real Dirichlet_u {zero};
                    Real Dirichlet_v {zero};
                    Real Dirichlet_w {zero};

                    Real nx = u_bnorm(i,j,k,0);
                    Real ny = u_bnorm(i,j,k,1);
                    Real nz = u_bnorm(i,j,k,2);

                    if (l_surface_layer) {

                        // Average v and w onto the x-face
                        Real velx = u_arr(i,j,k);
                        Real vely = (v_volfrac(i-1,j  ,k) * v_arr(i-1,j  ,k) + v_volfrac(i,j  ,k) * v_arr(i,j  ,k)
                                    + v_volfrac(i-1,j+1,k) * v_arr(i-1,j+1,k) + v_volfrac(i,j+1,k) * v_arr(i,j+1,k))
                                    / (v_volfrac(i-1,j,k) + v_volfrac(i,j,k) + v_volfrac(i-1,j+1,k) + v_volfrac(i,j+1,k));

                        Real velz = (w_volfrac(i-1,j,k  ) * w_arr(i-1,j,k  ) + w_volfrac(i,j,k  ) * w_arr(i,j,k  )
                                    + w_volfrac(i-1,j,k+1) * w_arr(i-1,j,k+1) + w_volfrac(i,j,k+1) * w_arr(i,j,k+1))
                                    / (w_volfrac(i-1,j,k) + w_volfrac(i,j,k) + w_volfrac(i-1,j,k+1) + w_volfrac(i,j,k+1));

                        // Impose tangential velocity as Dirichlet condition
                        Real v_dot_n = velx * nx + vely * ny + velz * nz;
                        Dirichlet_u = velx - v_dot_n * nx;
                        Dirichlet_v = vely - v_dot_n * ny;
                        Dirichlet_w = velz - v_dot_n * nz;
                    }

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
                    Real dvdz = slopes_v[2];
                    Real dwdx = slopes_w[0];
                    Real dwdy = slopes_w[1];
                    Real dwdz = slopes_w[2];

                    Real tau11_eb = ( dudx - ( dudx + dvdy + dwdz ) / three );
                    Real tau12_eb = myhalf * (dudy + dvdx);
                    Real tau13_eb = myhalf * (dudz + dwdx);

                    if (l_no_slip) {

                        dudn = - mu_eff * (nx * tau11_eb + ny * tau12_eb + nz * tau13_eb);

                    } else if (l_surface_layer) {

                        Real tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z;
                        compute_tangent_vectors(nx, ny, nz, tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z);

                        Real tau22_eb = ( dvdy - ( dudx + dvdy + dwdz ) / three );
                        Real tau33_eb = ( dwdz - ( dudx + dvdy + dwdz ) / three );
                        Real tau23_eb = myhalf * (dvdz + dwdy);

                        Real tauzz = mu_eff * ( nx*nx*tau11_eb + ny*ny*tau22_eb + nz*nz*tau33_eb
                                            + two * (nx*ny*tau12_eb + ny*nz*tau23_eb + nx*nz*tau13_eb ));

                        dudn = - tbx_x * u_tau_eb13(i,j,k) - tby_x * u_tau_eb23(i,j,k) - nx * tauzz;
                    }
                }

                rho_u_rhs(i,j,k) -= barea * dudn / (vol * u_volfrac(i,j,k));
            }
        }

    });

    // y-momentum
    ParallelFor(bxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (v_volfrac(i,j,k)>zero) {

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

            if (!l_constraint_y && v_cellflg(i,j,k).isSingleValued()) {

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

                Real dvdn = 0.0;

                if (l_no_slip || l_surface_layer) {

                    RealVect bcent_eb {v_bcent(i,j,k,0), v_bcent(i,j,k,1), v_bcent(i,j,k,2)};

                    Real Dirichlet_u {zero};
                    Real Dirichlet_v {zero};
                    Real Dirichlet_w {zero};

                    Real nx = v_bnorm(i,j,k,0);
                    Real ny = v_bnorm(i,j,k,1);
                    Real nz = v_bnorm(i,j,k,2);

                    if (l_surface_layer) {

                        // Average u and w onto the y-face
                        Real velx = (u_volfrac(i  ,j-1,k) * u_arr(i  ,j-1,k) + u_volfrac(i+1,j-1,k) * u_arr(i+1,j-1,k)
                                    + u_volfrac(i+1,j  ,k) * u_arr(i+1,j  ,k) + u_volfrac(i  ,j  ,k) * u_arr(i  ,j  ,k))
                                    / (u_volfrac(i,j-1,k) + u_volfrac(i+1,j-1,k) + u_volfrac(i+1,j,k) + u_volfrac(i,j,k));
                        Real vely = v_arr(i,j,k);
                        Real velz = (w_volfrac(i,j-1,k  ) * w_arr(i,j-1,k  ) + w_volfrac(i,j,k  ) * w_arr(i,j,k  )
                                    + w_volfrac(i,j  ,k+1) * w_arr(i,j  ,k+1) + w_volfrac(i,j-1,k+1) * w_arr(i,j-1,k+1))
                                    / (w_volfrac(i,j-1,k) + w_volfrac(i,j,k) + w_volfrac(i,j,k+1) + w_volfrac(i,j-1,k+1));

                        // Impose tangential velocity as Dirichlet condition
                        Real v_dot_n = velx * nx + vely * ny + velz * nz;
                        Dirichlet_u = velx - v_dot_n * nx;
                        Dirichlet_v = vely - v_dot_n * ny;
                        Dirichlet_w = velz - v_dot_n * nz;
                    }

                    GpuArray<Real,AMREX_SPACEDIM> slopes_u;
                    GpuArray<Real,AMREX_SPACEDIM> slopes_v;
                    GpuArray<Real,AMREX_SPACEDIM> slopes_w;

                    slopes_u = erf_calc_slopes_eb_Dirichlet_staggered( Vars::yvel, Vars::xvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
                    slopes_v = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
                    slopes_w = erf_calc_slopes_eb_Dirichlet_staggered( Vars::yvel, Vars::zvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);

                    Real dudx = slopes_u[0];
                    Real dudy = slopes_u[1];
                    Real dudz = slopes_u[2];
                    Real dvdx = slopes_v[0];
                    Real dvdy = slopes_v[1];
                    Real dvdz = slopes_v[2];
                    Real dwdx = slopes_w[0];
                    Real dwdy = slopes_w[1];
                    Real dwdz = slopes_w[2];

                    Real tau22_eb = ( dvdy - ( dudx + dvdy + dwdz ) / three );
                    Real tau12_eb = myhalf * (dudy + dvdx);
                    Real tau23_eb = myhalf * (dvdz + dwdy);

                    if (l_no_slip) {

                        dvdn = - mu_eff * (nx * tau12_eb + ny * tau22_eb + nz * tau23_eb);

                    } else if (l_surface_layer) {

                        Real tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z;
                        compute_tangent_vectors(nx, ny, nz, tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z);

                        Real tau11_eb = ( dudx - ( dudx + dvdy + dwdz ) / three );
                        Real tau33_eb = ( dwdz - ( dudx + dvdy + dwdz ) / three );
                        Real tau13_eb = myhalf * (dudz + dwdx);

                        Real tauzz = mu_eff * ( nx*nx*tau11_eb + ny*ny*tau22_eb + nz*nz*tau33_eb
                                            + two * (nx*ny*tau12_eb + ny*nz*tau23_eb + nx*nz*tau13_eb ));

                        dvdn = - tbx_y * v_tau_eb13(i,j,k) - tby_y * v_tau_eb23(i,j,k) - ny * tauzz;
                    }
                }

                rho_v_rhs(i,j,k) -= barea * dvdn / (vol * v_volfrac(i,j,k));
            }
        }
    });

    // z-momentum
    ParallelFor(bxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (w_volfrac(i,j,k)>zero) {

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

            if (!l_constraint_z && w_cellflg(i,j,k).isSingleValued()) {

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

                Real dwdn = zero;

                if (l_no_slip || l_surface_layer) {

                    const RealVect bcent_eb {w_bcent(i,j,k,0), w_bcent(i,j,k,1), w_bcent(i,j,k,2)};

                    Real Dirichlet_u {zero};
                    Real Dirichlet_v {zero};
                    Real Dirichlet_w {zero};

                    Real nx = w_bnorm(i,j,k,0);
                    Real ny = w_bnorm(i,j,k,1);
                    Real nz = w_bnorm(i,j,k,2);

                    if (l_surface_layer) {

                        // Average u and v onto the z-face
                        Real velx = (u_volfrac(i  ,j,k-1) * u_arr(i  ,j,k-1) + u_volfrac(i+1,j,k-1) * u_arr(i+1,j,k-1)
                                    + u_volfrac(i+1,j,k  ) * u_arr(i+1,j,k  ) + u_volfrac(i  ,j,k  ) * u_arr(i  ,j,k  ))
                                    / (u_volfrac(i,j,k-1) + u_volfrac(i+1,j,k-1) + u_volfrac(i+1,j,k) + u_volfrac(i,j,k));
                        Real vely = (v_volfrac(i,j  ,k-1) * v_arr(i,j  ,k-1) + v_volfrac(i,j+1,k-1) * v_arr(i,j+1,k-1)
                                    + v_volfrac(i,j+1,k  ) * v_arr(i,j+1,k  ) + v_volfrac(i,j  ,k  ) * v_arr(i,j  ,k  ))
                                    / (v_volfrac(i,j,k-1) + v_volfrac(i,j+1,k-1) + v_volfrac(i,j+1,k) + v_volfrac(i,j,k));
                        Real velz = w_arr(i,j,k);

                        // Impose tangential velocity as Dirichlet condition
                        Real v_dot_n = velx * nx + vely * ny + velz * nz;
                        Dirichlet_u = velx - v_dot_n * nx;
                        Dirichlet_v = vely - v_dot_n * ny;
                        Dirichlet_w = velz - v_dot_n * nz;
                    }

                    GpuArray<Real,AMREX_SPACEDIM> slopes_u;
                    GpuArray<Real,AMREX_SPACEDIM> slopes_v;
                    GpuArray<Real,AMREX_SPACEDIM> slopes_w;

                    slopes_u = erf_calc_slopes_eb_Dirichlet_staggered( Vars::zvel, Vars::xvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
                    slopes_v = erf_calc_slopes_eb_Dirichlet_staggered( Vars::zvel, Vars::yvel, dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
                    slopes_w = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);

                    Real dudx = slopes_u[0];
                    Real dudy = slopes_u[1];
                    Real dudz = slopes_u[2];
                    Real dvdx = slopes_v[0];
                    Real dvdy = slopes_v[1];
                    Real dvdz = slopes_v[2];
                    Real dwdx = slopes_w[0];
                    Real dwdy = slopes_w[1];
                    Real dwdz = slopes_w[2];

                    Real tau33_eb = ( dwdz - ( dudx + dvdy + dwdz ) / three );
                    Real tau13_eb = myhalf * (dudz + dwdx);
                    Real tau23_eb = myhalf * (dvdz + dwdy);

                    if (l_no_slip) {

                        dwdn = - mu_eff * (nx * tau13_eb + ny * tau23_eb + nz * tau33_eb);

                    } else if (l_surface_layer) {

                        Real tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z;
                        compute_tangent_vectors(nx, ny, nz, tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z);

                        Real tau11_eb = ( dudx - ( dudx + dvdy + dwdz ) / three );
                        Real tau22_eb = ( dvdy - ( dudx + dvdy + dwdz ) / three );
                        Real tau12_eb = myhalf * (dudy + dvdx);

                        Real tauzz = mu_eff * ( nx*nx*tau11_eb + ny*ny*tau22_eb + nz*nz*tau33_eb
                                            + two * (nx*ny*tau12_eb + ny*nz*tau23_eb + nx*nz*tau13_eb ));

                        dwdn = - tbx_z * w_tau_eb13(i,j,k) - tby_z * w_tau_eb23(i,j,k) - nz * tauzz;

                    }
                }

                rho_w_rhs(i,j,k) -= barea * dwdn / (vol * w_volfrac(i,j,k));
            }
        }
    });

}
