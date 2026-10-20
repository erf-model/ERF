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
 * Compute the tangent directions at an EB boundary given the unit normal vector.
 *
 * t_bx and t_by are the unit vectors along the projections of e_x and e_y onto the
 * tangent plane. They are the directions in which the stored EB surface stresses
 * tau_eb13 and tau_eb23 act, since those are built from the x- and y-components of
 * the tangential velocity (see ERF_EBMOSTStress.H). For a horizontal surface they are
 * exactly e_x and e_y. Note that they are individually unit but are not orthogonal to
 * each other unless nx*ny is zero -- t_bx . t_by = -nx*ny.
 *
 * e_x has no projection onto the tangent plane when n is parallel to e_x, and likewise
 * for e_y, so one of the two is undefined for an axis-aligned normal such as (1,0,0)
 * (issue #3533). That limit is genuinely two-sided -- approaching n = (1,0,0) from
 * nz > 0 and from nz < 0 gives t_bx of -e_z and +e_z -- so there is no continuous
 * choice, and the degenerate vector is returned as zero. That keeps the momentum RHS
 * finite, and it keeps the component that does have a limit continuous: the x-component
 * of t_bx is sqrt(1-nx*nx), which tends to zero here as well. It also avoids driving
 * w-momentum with tau_eb13, an x-associated stress, along an arbitrarily signed
 * direction. Only nx and ny can be degenerate, and never both at once.
 *
 * @param[in]  nx x-component of the normal vector
 * @param[in]  ny y-component of the normal vector
 * @param[in]  nz z-component of the normal vector
 * @param[out] tbx_x x-component of first tangent vector
 * @param[out] tbx_y y-component of first tangent vector
 * @param[out] tbx_z z-component of first tangent vector
 * @param[out] tby_x x-component of second tangent vector
 * @param[out] tby_y y-component of second tangent vector
 * @param[out] tby_z z-component of second tangent vector
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_tangent_vectors (Real nx, Real ny, Real nz,
                               Real& tbx_x, Real& tbx_y, Real& tbx_z,
                               Real& tby_x, Real& tby_y, Real& tby_z)
{
    // Below this the projection is all round-off and its direction is meaningless.
    // The squared norms are formed as ny^2+nz^2 and nx^2+nz^2 rather than as the
    // algebraically equal 1-nx^2 and 1-ny^2 to avoid cancellation for a near-axis normal
    // The squared norm of a unit normal carries an absolute roundoff of order
    // real_eps, so the floor is a (generous) multiple of it: 2.2e-12 in double,
    // 1.2e-03 in single, i.e. a normal within ~2 degrees of the axis.
    constexpr Real tol = Real(1.e4)*real_eps;

    // x-tangential vector: t_bx = (e_x - (e_x · n)n) / ||e_x - (e_x · n)n||
    // e_x = (1,0,0), so e_x · n = nx
    Real tbx_norm2 = ny*ny + nz*nz;
    if (tbx_norm2 > tol) {
        Real tbx_norm_inv = one / std::sqrt(tbx_norm2);
        tbx_x = ( one - nx * nx) * tbx_norm_inv;
        tbx_y = (     - nx * ny) * tbx_norm_inv;
        tbx_z = (     - nx * nz) * tbx_norm_inv;
    } else {
        tbx_x = zero;
        tbx_y = zero;
        tbx_z = zero;
    }

    // y-tangential vector: t_by = (e_y - (e_y · n)n) / ||e_y - (e_y · n)n||
    // e_y = (0,1,0), so e_y · n = ny
    Real tby_norm2 = nx*nx + nz*nz;
    if (tby_norm2 > tol) {
        Real tby_norm_inv = one / std::sqrt(tby_norm2);
        tby_x = (     - ny * nx) * tby_norm_inv;
        tby_y = ( one - ny * ny) * tby_norm_inv;
        tby_z = (     - ny * nz) * tby_norm_inv;
    } else {
        tby_x = zero;
        tby_y = zero;
        tby_z = zero;
    }
}

namespace {

/**
 * Compute the EB wall contribution for one staggered momentum component.
 *
 * The momentum direction and wall model are compile-time parameters so that a
 * correction kernel contains only the selected component and boundary model.
 * This routine is called only for single-valued cells.
 */
template <int MomDir, EBBoundaryType WallModel>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real eb_wall_correction (
    int i, int j, int k,
    const Array4<const Real>& u_arr,
    const Array4<const Real>& v_arr,
    const Array4<const Real>& w_arr,
    const Array4<const EBCellFlag>& u_cellflg,
    const Array4<const EBCellFlag>& v_cellflg,
    const Array4<const EBCellFlag>& w_cellflg,
    const Array4<const Real>& u_volcent,
    const Array4<const Real>& v_volcent,
    const Array4<const Real>& w_volcent,
    const Array4<const Real>& u_vfrac,
    const Array4<const Real>& v_vfrac,
    const Array4<const Real>& w_vfrac,
    const Array4<const Real>& mom_volfrac,
    const Array4<const Real>& mom_afrac_x,
    const Array4<const Real>& mom_afrac_y,
    const Array4<const Real>& mom_afrac_z,
    const Array4<const Real>& mom_bcent,
    const Array4<const Real>& mom_bnorm,
    const Array4<const Real>& mom_tau_eb13,
    const Array4<const Real>& mom_tau_eb23,
    Real dx, Real dy, Real dz, Real vol, Real mu_eff) noexcept
{
    static_assert(MomDir >= 0 && MomDir < AMREX_SPACEDIM);
    static_assert(WallModel == EBBoundaryType::NoSlipWall ||
                  WallModel == EBBoundaryType::SurfaceLayer);

    RealVect bcent_eb {mom_bcent(i,j,k,0), mom_bcent(i,j,k,1), mom_bcent(i,j,k,2)};

    Real nx = mom_bnorm(i,j,k,0);
    Real ny = mom_bnorm(i,j,k,1);
    Real nz = mom_bnorm(i,j,k,2);

    Real Dirichlet_u {zero};
    Real Dirichlet_v {zero};
    Real Dirichlet_w {zero};

    if constexpr (WallModel == EBBoundaryType::SurfaceLayer) {
        Real velx;
        Real vely;
        Real velz;

        if constexpr (MomDir == 0) {
            // Average v and w onto the x-face.
            velx = u_arr(i,j,k);
            vely = (v_vfrac(i-1,j  ,k) * v_arr(i-1,j  ,k) + v_vfrac(i,j  ,k) * v_arr(i,j  ,k)
                  + v_vfrac(i-1,j+1,k) * v_arr(i-1,j+1,k) + v_vfrac(i,j+1,k) * v_arr(i,j+1,k))
                  / (v_vfrac(i-1,j,k) + v_vfrac(i,j,k) + v_vfrac(i-1,j+1,k) + v_vfrac(i,j+1,k));
            velz = (w_vfrac(i-1,j,k  ) * w_arr(i-1,j,k  ) + w_vfrac(i,j,k  ) * w_arr(i,j,k  )
                  + w_vfrac(i-1,j,k+1) * w_arr(i-1,j,k+1) + w_vfrac(i,j,k+1) * w_arr(i,j,k+1))
                  / (w_vfrac(i-1,j,k) + w_vfrac(i,j,k) + w_vfrac(i-1,j,k+1) + w_vfrac(i,j,k+1));
        } else if constexpr (MomDir == 1) {
            // Average u and w onto the y-face.
            velx = (u_vfrac(i  ,j-1,k) * u_arr(i  ,j-1,k) + u_vfrac(i+1,j-1,k) * u_arr(i+1,j-1,k)
                  + u_vfrac(i+1,j  ,k) * u_arr(i+1,j  ,k) + u_vfrac(i  ,j  ,k) * u_arr(i  ,j  ,k))
                  / (u_vfrac(i,j-1,k) + u_vfrac(i+1,j-1,k) + u_vfrac(i+1,j,k) + u_vfrac(i,j,k));
            vely = v_arr(i,j,k);
            velz = (w_vfrac(i,j-1,k  ) * w_arr(i,j-1,k  ) + w_vfrac(i,j,k  ) * w_arr(i,j,k  )
                  + w_vfrac(i,j  ,k+1) * w_arr(i,j  ,k+1) + w_vfrac(i,j-1,k+1) * w_arr(i,j-1,k+1))
                  / (w_vfrac(i,j-1,k) + w_vfrac(i,j,k) + w_vfrac(i,j,k+1) + w_vfrac(i,j-1,k+1));
        } else {
            // Average u and v onto the z-face.
            velx = (u_vfrac(i  ,j,k-1) * u_arr(i  ,j,k-1) + u_vfrac(i+1,j,k-1) * u_arr(i+1,j,k-1)
                  + u_vfrac(i+1,j,k  ) * u_arr(i+1,j,k  ) + u_vfrac(i  ,j,k  ) * u_arr(i  ,j,k  ))
                  / (u_vfrac(i,j,k-1) + u_vfrac(i+1,j,k-1) + u_vfrac(i+1,j,k) + u_vfrac(i,j,k));
            vely = (v_vfrac(i,j  ,k-1) * v_arr(i,j  ,k-1) + v_vfrac(i,j+1,k-1) * v_arr(i,j+1,k-1)
                  + v_vfrac(i,j+1,k  ) * v_arr(i,j+1,k  ) + v_vfrac(i,j  ,k  ) * v_arr(i,j  ,k  ))
                  / (v_vfrac(i,j,k-1) + v_vfrac(i,j+1,k-1) + v_vfrac(i,j+1,k) + v_vfrac(i,j,k));
            velz = w_arr(i,j,k);
        }

        Real v_dot_n = velx * nx + vely * ny + velz * nz;
        Dirichlet_u = velx - v_dot_n * nx;
        Dirichlet_v = vely - v_dot_n * ny;
        Dirichlet_w = velz - v_dot_n * nz;
    }

    GpuArray<Real,AMREX_SPACEDIM> slopes_u;
    GpuArray<Real,AMREX_SPACEDIM> slopes_v;
    GpuArray<Real,AMREX_SPACEDIM> slopes_w;

    if constexpr (MomDir == 0) {
        slopes_u = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
        slopes_v = erf_calc_slopes_eb_Dirichlet_staggered<Vars::xvel,Vars::yvel>(dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
        slopes_w = erf_calc_slopes_eb_Dirichlet_staggered<Vars::xvel,Vars::zvel>(dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);
    } else if constexpr (MomDir == 1) {
        slopes_u = erf_calc_slopes_eb_Dirichlet_staggered<Vars::yvel,Vars::xvel>(dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
        slopes_v = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
        slopes_w = erf_calc_slopes_eb_Dirichlet_staggered<Vars::yvel,Vars::zvel>(dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);
    } else {
        slopes_u = erf_calc_slopes_eb_Dirichlet_staggered<Vars::zvel,Vars::xvel>(dx, dy, dz, i, j, k, bcent_eb, Dirichlet_u, u_arr, u_volcent, u_cellflg);
        slopes_v = erf_calc_slopes_eb_Dirichlet_staggered<Vars::zvel,Vars::yvel>(dx, dy, dz, i, j, k, bcent_eb, Dirichlet_v, v_arr, v_volcent, v_cellflg);
        slopes_w = erf_calc_slopes_eb_Dirichlet          (                         dx, dy, dz, i, j, k, bcent_eb, Dirichlet_w, w_arr, w_volcent, w_cellflg);
    }

    Real dudx = slopes_u[0];
    Real dudy = slopes_u[1];
    Real dudz = slopes_u[2];
    Real dvdx = slopes_v[0];
    Real dvdy = slopes_v[1];
    Real dvdz = slopes_v[2];
    Real dwdx = slopes_w[0];
    Real dwdy = slopes_w[1];
    Real dwdz = slopes_w[2];

    Real divergence = dudx + dvdy + dwdz;
    Real dmomdn;

    if constexpr (WallModel == EBBoundaryType::NoSlipWall) {
        if constexpr (MomDir == 0) {
            Real tau11_eb = dudx - divergence / three;
            Real tau12_eb = myhalf * (dudy + dvdx);
            Real tau13_eb = myhalf * (dudz + dwdx);
            dmomdn = -mu_eff * (nx * tau11_eb + ny * tau12_eb + nz * tau13_eb);
        } else if constexpr (MomDir == 1) {
            Real tau22_eb = dvdy - divergence / three;
            Real tau12_eb = myhalf * (dudy + dvdx);
            Real tau23_eb = myhalf * (dvdz + dwdy);
            dmomdn = -mu_eff * (nx * tau12_eb + ny * tau22_eb + nz * tau23_eb);
        } else {
            Real tau33_eb = dwdz - divergence / three;
            Real tau13_eb = myhalf * (dudz + dwdx);
            Real tau23_eb = myhalf * (dvdz + dwdy);
            dmomdn = -mu_eff * (nx * tau13_eb + ny * tau23_eb + nz * tau33_eb);
        }
    } else {
        Real tau11_eb = dudx - divergence / three;
        Real tau22_eb = dvdy - divergence / three;
        Real tau33_eb = dwdz - divergence / three;
        Real tau12_eb = myhalf * (dudy + dvdx);
        Real tau13_eb = myhalf * (dudz + dwdx);
        Real tau23_eb = myhalf * (dvdz + dwdy);

        Real tauzz = mu_eff * ( nx*nx*tau11_eb + ny*ny*tau22_eb + nz*nz*tau33_eb
                              + two * (nx*ny*tau12_eb + ny*nz*tau23_eb + nx*nz*tau13_eb) );

        Real tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z;
        compute_tangent_vectors(nx, ny, nz, tbx_x, tbx_y, tbx_z, tby_x, tby_y, tby_z);

        if constexpr (MomDir == 0) {
            dmomdn = -tbx_x * mom_tau_eb13(i,j,k) - tby_x * mom_tau_eb23(i,j,k) - nx * tauzz;
        } else if constexpr (MomDir == 1) {
            dmomdn = -tbx_y * mom_tau_eb13(i,j,k) - tby_y * mom_tau_eb23(i,j,k) - ny * tauzz;
        } else {
            dmomdn = -tbx_z * mom_tau_eb13(i,j,k) - tby_z * mom_tau_eb23(i,j,k) - nz * tauzz;
        }
    }

    Real axm = mom_afrac_x(i  ,j  ,k  );
    Real axp = mom_afrac_x(i+1,j  ,k  );
    Real aym = mom_afrac_y(i  ,j  ,k  );
    Real ayp = mom_afrac_y(i  ,j+1,k  );
    Real azm = mom_afrac_z(i  ,j  ,k  );
    Real azp = mom_afrac_z(i  ,j  ,k+1);

    Real adx = (axm-axp) * dy * dz;
    Real ady = (aym-ayp) * dx * dz;
    Real adz = (azm-azp) * dx * dy;
    Real barea = std::sqrt(adx*adx + ady*ady + adz*adz);

    return barea * dmomdn / (vol * mom_volfrac(i,j,k));
}

/**
 * Launch an EB wall-correction kernel over the cached cut-face list on GPU.
 * CPU execution retains the full-box kernel.
 */
template <int MomDir, EBBoundaryType WallModel>
void launch_eb_wall_correction (
    Box const& bx,
    [[maybe_unused]] eb_::CutCellList const& cut_cells,
    Array4<Real> const& mom_rhs,
    const Array4<const Real>& u_arr,
    const Array4<const Real>& v_arr,
    const Array4<const Real>& w_arr,
    const Array4<const EBCellFlag>& u_cellflg,
    const Array4<const EBCellFlag>& v_cellflg,
    const Array4<const EBCellFlag>& w_cellflg,
    const Array4<const Real>& u_volcent,
    const Array4<const Real>& v_volcent,
    const Array4<const Real>& w_volcent,
    const Array4<const Real>& u_vfrac,
    const Array4<const Real>& v_vfrac,
    const Array4<const Real>& w_vfrac,
    const Array4<const Real>& mom_volfrac,
    const Array4<const Real>& mom_afrac_x,
    const Array4<const Real>& mom_afrac_y,
    const Array4<const Real>& mom_afrac_z,
    const Array4<const Real>& mom_bcent,
    const Array4<const Real>& mom_bnorm,
    const Array4<const Real>& mom_tau_eb13,
    const Array4<const Real>& mom_tau_eb23,
    Real dx, Real dy, Real dz, Real vol, Real mu_eff) noexcept
{
    static_assert(MomDir >= 0 && MomDir < AMREX_SPACEDIM);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        Dim3 const* const cut_cells_ptr = cut_cells.data();
        int const num_cut_cells = static_cast<int>(cut_cells.size());
        ParallelFor(num_cut_cells, [=] AMREX_GPU_DEVICE (int cut_index) noexcept {
            Dim3 const cell = cut_cells_ptr[cut_index];
            if (!bx.contains(cell)) { return; }
            mom_rhs(cell.x,cell.y,cell.z) -= eb_wall_correction<MomDir,WallModel>(
                cell.x, cell.y, cell.z, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                u_vfrac, v_vfrac, w_vfrac, mom_volfrac,
                mom_afrac_x, mom_afrac_y, mom_afrac_z, mom_bcent, mom_bnorm,
                mom_tau_eb13, mom_tau_eb23, dx, dy, dz, vol, mu_eff);
        });
        return;
    }
#endif

    Array4<const EBCellFlag> mom_cellflg;
    if constexpr (MomDir == 0) {
        mom_cellflg = u_cellflg;
    } else if constexpr (MomDir == 1) {
        mom_cellflg = v_cellflg;
    } else {
        mom_cellflg = w_cellflg;
    }

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (mom_cellflg(i,j,k).isSingleValued()) {
            mom_rhs(i,j,k) -= eb_wall_correction<MomDir,WallModel>(
                i, j, k, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                u_vfrac, v_vfrac, w_vfrac, mom_volfrac,
                mom_afrac_x, mom_afrac_y, mom_afrac_z, mom_bcent, mom_bnorm,
                mom_tau_eb13, mom_tau_eb23, dx, dy, dz, vol, mu_eff);
        }
    });
}

} // namespace

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
 * @param[in]  u_arr x-direction velocity
 * @param[in]  v_arr y-direction velocity
 * @param[in]  w_arr z-direction velocity
 * @param[in]  tau11 11 stress
 * @param[in]  tau22 22 stress
 * @param[in]  tau33 33 stress
 * @param[in]  tau12 12 stress
 * @param[in]  tau13 13 stress
 * @param[in]  tau23 23 stress
 * @param[in]  u_tau_eb13 EB tangential stress for x-momentum on z-faces
 * @param[in]  u_tau_eb23 EB tangential stress for x-momentum on y-faces
 * @param[in]  v_tau_eb13 EB tangential stress for y-momentum on z-faces
 * @param[in]  v_tau_eb23 EB tangential stress for y-momentum on x-faces
 * @param[in]  w_tau_eb13 EB tangential stress for z-momentum on y-faces
 * @param[in]  w_tau_eb23 EB tangential stress for z-momentum on x-faces
 * @param[in]  dx_arr cell size array
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_mx x map factor at cell centers
 * @param[in]  mf_ux x map factor at x-faces
 * @param[in]  mf_vx x map factor at y-faces
 * @param[in]  mf_my y map factor at cell centers
 * @param[in]  mf_uy y map factor at x-faces
 * @param[in]  mf_vy y map factor at y-faces
 * @param[in] solverChoice container with diffusion parameters
 * @param[in] ebfact EB factories for cell- and face-centered variables
 * @param[in] d_bcrec_ptr boundary condition records on device
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

    const auto* u_factory = ebfact.get_u_const_factory();
    const auto* v_factory = ebfact.get_v_const_factory();
    const auto* w_factory = ebfact.get_w_const_factory();

    auto const& u_cut_cells = ebfact.get_u_cut_cells(mfi.index());
    auto const& v_cut_cells = ebfact.get_v_cut_cells(mfi.index());
    auto const& w_cut_cells = ebfact.get_w_cut_cells(mfi.index());

    const auto& u_flagfab = u_factory->getMultiEBCellFlagFab()[mfi];
    const auto& v_flagfab = v_factory->getMultiEBCellFlagFab()[mfi];
    const auto& w_flagfab = w_factory->getMultiEBCellFlagFab()[mfi];

    FabType u_type = u_flagfab.getType();
    FabType v_type = v_flagfab.getType();
    FabType w_type = w_flagfab.getType();

    const bool l_wall_model = l_no_slip || l_surface_layer;
    const bool l_u_correction = l_wall_model && !l_constraint_x &&
                                u_type == FabType::singlevalued;
    const bool l_v_correction = l_wall_model && !l_constraint_y &&
                                v_type == FabType::singlevalued;
    const bool l_w_correction = l_wall_model && !l_constraint_z &&
                                w_type == FabType::singlevalued;
    const bool l_any_correction = l_u_correction || l_v_correction || l_w_correction;

    Array4<const Real> u_volfrac = u_factory->getVolFrac().const_array(mfi);
    Array4<const Real> v_volfrac = v_factory->getVolFrac().const_array(mfi);
    Array4<const Real> w_volfrac = w_factory->getVolFrac().const_array(mfi);

    Array4<const EBCellFlag> u_cellflg = u_flagfab.const_array();
    Array4<const EBCellFlag> v_cellflg = v_flagfab.const_array();
    Array4<const EBCellFlag> w_cellflg = w_flagfab.const_array();
    Array4<const Real> u_volcent{};
    Array4<const Real> v_volcent{};
    Array4<const Real> w_volcent{};

    if (l_any_correction) {
        if (u_type == FabType::singlevalued) {
            u_volcent = u_factory->getCentroid().const_array(mfi);
        }
        if (v_type == FabType::singlevalued) {
            v_volcent = v_factory->getCentroid().const_array(mfi);
        }
        if (w_type == FabType::singlevalued) {
            w_volcent = w_factory->getCentroid().const_array(mfi);
        }
    }

    Array4<const Real> u_afrac_x{};
    Array4<const Real> u_afrac_y{};
    Array4<const Real> u_afrac_z{};
    Array4<const Real> v_afrac_x{};
    Array4<const Real> v_afrac_y{};
    Array4<const Real> v_afrac_z{};
    Array4<const Real> w_afrac_x{};
    Array4<const Real> w_afrac_y{};
    Array4<const Real> w_afrac_z{};

    if (u_type == FabType::singlevalued) {
        u_afrac_x = u_factory->getAreaFrac()[0]->const_array(mfi);
        u_afrac_y = u_factory->getAreaFrac()[1]->const_array(mfi);
        u_afrac_z = u_factory->getAreaFrac()[2]->const_array(mfi);
    }
    if (v_type == FabType::singlevalued) {
        v_afrac_x = v_factory->getAreaFrac()[0]->const_array(mfi);
        v_afrac_y = v_factory->getAreaFrac()[1]->const_array(mfi);
        v_afrac_z = v_factory->getAreaFrac()[2]->const_array(mfi);
    }
    if (w_type == FabType::singlevalued) {
        w_afrac_x = w_factory->getAreaFrac()[0]->const_array(mfi);
        w_afrac_y = w_factory->getAreaFrac()[1]->const_array(mfi);
        w_afrac_z = w_factory->getAreaFrac()[2]->const_array(mfi);
    }

    Array4<const Real> u_bcent{};
    Array4<const Real> u_bnorm{};
    Array4<const Real> v_bcent{};
    Array4<const Real> v_bnorm{};
    Array4<const Real> w_bcent{};
    Array4<const Real> w_bnorm{};

    if (l_u_correction) {
        u_bcent = u_factory->getBndryCent().const_array(mfi);
        u_bnorm = u_factory->getBndryNormal().const_array(mfi);
    }
    if (l_v_correction) {
        v_bcent = v_factory->getBndryCent().const_array(mfi);
        v_bnorm = v_factory->getBndryNormal().const_array(mfi);
    }
    if (l_w_correction) {
        w_bcent = w_factory->getBndryCent().const_array(mfi);
        w_bnorm = w_factory->getBndryNormal().const_array(mfi);
    }

    const Array4<const Real> unused_real{};

    // x-momentum
    if (u_type == FabType::regular) {

        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

            Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

            rho_u_rhs(i,j,k) -= ( (tau11(i  , j  , k  ) - tau11(i-1, j  , k  ) ) * dxinv * mfsq
                                + (tau12(i  , j+1, k  ) - tau12(i  , j  , k  ) ) * dyinv * mfsq
                                + (tau13(i  , j  , k+1) - tau13(i  , j  , k  ) ) * dzinv );
        });

    } else if (u_type == FabType::singlevalued) {

        ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            EBCellFlag const cell_flag = u_cellflg(i,j,k);

            if (cell_flag.isRegular()) {
                Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

                rho_u_rhs(i,j,k) -= ( (tau11(i  , j  , k  ) - tau11(i-1, j  , k  ) ) * dxinv * mfsq
                                    + (tau12(i  , j+1, k  ) - tau12(i  , j  , k  ) ) * dyinv * mfsq
                                    + (tau13(i  , j  , k+1) - tau13(i  , j  , k  ) ) * dzinv );
            } else if (cell_flag.isSingleValued()) {
                Real vfrac = u_volfrac(i,j,k);
                if (vfrac <= zero) { return; }

                // Inv Jacobian
                Real mfsq = mf_ux(i,j,0) * mf_uy(i,j,0);

                Real diffContrib  = ( (tau11(i  , j  , k  ) * u_afrac_x(i+1,j  ,k  )
                                    -  tau11(i-1, j  , k  ) * u_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                    + (tau12(i  , j+1, k  ) * u_afrac_y(i  ,j+1,k  )
                                    -  tau12(i  , j  , k  ) * u_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                    + (tau13(i  , j  , k+1) * u_afrac_z(i  ,j  ,k+1)
                                    -  tau13(i  , j  , k  ) * u_afrac_z(i  ,j  ,k  )) * dzinv );
                diffContrib      /= vfrac;

                rho_u_rhs(i,j,k) -= diffContrib;
            }
        });
    } // u_type

    if (l_u_correction) {
        if (l_no_slip) {
            launch_eb_wall_correction<0,EBBoundaryType::NoSlipWall>(
                bxx, u_cut_cells, rho_u_rhs, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                unused_real, unused_real, unused_real, u_volfrac,
                u_afrac_x, u_afrac_y, u_afrac_z, u_bcent, u_bnorm,
                unused_real, unused_real, dx, dy, dz, vol, mu_eff);
        } else {
            launch_eb_wall_correction<0,EBBoundaryType::SurfaceLayer>(
                bxx, u_cut_cells, rho_u_rhs, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                u_volfrac, v_volfrac, w_volfrac, u_volfrac,
                u_afrac_x, u_afrac_y, u_afrac_z, u_bcent, u_bnorm,
                u_tau_eb13, u_tau_eb23, dx, dy, dz, vol, mu_eff);
        }
    }

    // y-momentum
    if (v_type == FabType::regular) {

        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

            Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

            rho_v_rhs(i,j,k) -= ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  ) ) * dxinv * mfsq
                                + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  ) ) * dyinv * mfsq
                                + (tau23(i  , j  , k+1) - tau23(i  , j  , k  ) ) * dzinv );
        });
    } else if (v_type == FabType::singlevalued) {

        ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            EBCellFlag const cell_flag = v_cellflg(i,j,k);

            if (cell_flag.isRegular()) {
                Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

                rho_v_rhs(i,j,k) -= ( (tau12(i+1, j  , k  ) - tau12(i  , j  , k  ) ) * dxinv * mfsq
                                    + (tau22(i  , j  , k  ) - tau22(i  , j-1, k  ) ) * dyinv * mfsq
                                    + (tau23(i  , j  , k+1) - tau23(i  , j  , k  ) ) * dzinv );
            } else if (cell_flag.isSingleValued()) {
                Real vfrac = v_volfrac(i,j,k);
                if (vfrac <= zero) { return; }

                // Inv Jacobian
                Real mfsq = mf_vx(i,j,0) * mf_vy(i,j,0);

                Real diffContrib  = ( (tau12(i+1, j  , k  ) * v_afrac_x(i+1,j  ,k  )
                                    -  tau12(i  , j  , k  ) * v_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                    + (tau22(i  , j  , k  ) * v_afrac_y(i  ,j+1,k  )
                                    -  tau22(i  , j-1, k  ) * v_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                    + (tau23(i  , j  , k+1) * v_afrac_z(i  ,j  ,k+1)
                                    -  tau23(i  , j  , k  ) * v_afrac_z(i  ,j  ,k  ) ) * dzinv );
                diffContrib      /= vfrac;

                rho_v_rhs(i,j,k) -= diffContrib;
            }
        });
    } // v_type

    if (l_v_correction) {
        if (l_no_slip) {
            launch_eb_wall_correction<1,EBBoundaryType::NoSlipWall>(
                bxy, v_cut_cells, rho_v_rhs, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                unused_real, unused_real, unused_real, v_volfrac,
                v_afrac_x, v_afrac_y, v_afrac_z, v_bcent, v_bnorm,
                unused_real, unused_real, dx, dy, dz, vol, mu_eff);
        } else {
            launch_eb_wall_correction<1,EBBoundaryType::SurfaceLayer>(
                bxy, v_cut_cells, rho_v_rhs, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                u_volfrac, v_volfrac, w_volfrac, v_volfrac,
                v_afrac_x, v_afrac_y, v_afrac_z, v_bcent, v_bnorm,
                v_tau_eb13, v_tau_eb23, dx, dy, dz, vol, mu_eff);
        }
    }

    // z-momentum
    if (w_type == FabType::regular) {

        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

            Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

            rho_w_rhs(i,j,k) -= ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  ) ) * dxinv * mfsq
                                + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  ) ) * dyinv * mfsq
                                + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1) ) * dzinv );
        });

    } else if (w_type == FabType::singlevalued) {

        ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            EBCellFlag const cell_flag = w_cellflg(i,j,k);

            if (cell_flag.isRegular()) {
                Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

                rho_w_rhs(i,j,k) -= ( (tau13(i+1, j  , k  ) - tau13(i  , j  , k  ) ) * dxinv * mfsq
                                    + (tau23(i  , j+1, k  ) - tau23(i  , j  , k  ) ) * dyinv * mfsq
                                    + (tau33(i  , j  , k  ) - tau33(i  , j  , k-1) ) * dzinv );
            } else if (cell_flag.isSingleValued()) {
                Real vfrac = w_volfrac(i,j,k);
                if (vfrac <= zero) { return; }

                // Inv Jacobian
                Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);

                Real diffContrib  = ( (tau13(i+1, j  , k  ) * w_afrac_x(i+1,j  ,k  )
                                    -  tau13(i  , j  , k  ) * w_afrac_x(i  ,j  ,k  ) ) * dxinv * mfsq
                                    + (tau23(i  , j+1, k  ) * w_afrac_y(i  ,j+1,k  )
                                    -  tau23(i  , j  , k  ) * w_afrac_y(i  ,j  ,k  ) ) * dyinv * mfsq
                                    + (tau33(i  , j  , k  ) * w_afrac_z(i  ,j  ,k+1)
                                    -  tau33(i  , j  , k-1) * w_afrac_z(i  ,j  ,k  ) ) * dzinv );
                diffContrib      /= vfrac;

                rho_w_rhs(i,j,k) -= diffContrib;
            }
        });
    } // w_type

    if (l_w_correction) {
        if (l_no_slip) {
            launch_eb_wall_correction<2,EBBoundaryType::NoSlipWall>(
                bxz, w_cut_cells, rho_w_rhs, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                unused_real, unused_real, unused_real, w_volfrac,
                w_afrac_x, w_afrac_y, w_afrac_z, w_bcent, w_bnorm,
                unused_real, unused_real, dx, dy, dz, vol, mu_eff);
        } else {
            launch_eb_wall_correction<2,EBBoundaryType::SurfaceLayer>(
                bxz, w_cut_cells, rho_w_rhs, u_arr, v_arr, w_arr,
                u_cellflg, v_cellflg, w_cellflg, u_volcent, v_volcent, w_volcent,
                u_volfrac, v_volfrac, w_volfrac, w_volfrac,
                w_afrac_x, w_afrac_y, w_afrac_z, w_bcent, w_bnorm,
                w_tau_eb13, w_tau_eb23, dx, dy, dz, vol, mu_eff);
        }
    }
}
