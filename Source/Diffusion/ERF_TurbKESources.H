#ifndef ERF_TURB_KE_SOURCES_H_
#define ERF_TURB_KE_SOURCES_H_

#include <AMReX.H>
#include <AMReX_BCRec.H>

#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_NumericalConstants.H"
#include "ERF_SurfaceLayer.H"
#include "ERF_PBLModels.H"

/**
 * Add the production and dissipation sources of the subgrid TKE (Deardorff or
 * k-eqn RANS) or of the MYNN QKE to the RhoKE RHS.
 *
 * These sources consume the face fluxes written by the scalar diffusion of the
 * same RK stage, so they must be called after every DiffusionSrcForState_* call
 * of the stage has returned: the RhoTheta call in erf_slow_rhs_pre (hfx_z) and
 * the RhoQ1 call in erf_slow_rhs_post (qfx1_z).
 *
 * @param[in]    bx cell center box to loop over
 * @param[in]    domain box of the whole domain
 * @param[in]    u velocity in x-dir
 * @param[in]    v velocity in y-dir
 * @param[in]    cell_data conserved cell center vars
 * @param[in]    cell_prim primitive cell center vars
 * @param[inout] cell_rhs RHS for cell center vars
 * @param[in]    SmnSmn_a strain rate magnitude
 * @param[in]    hfx_z rho<w'theta'> on z-faces
 * @param[in]    qfx1_z rho<w'qv'> on z-faces (may be null)
 * @param[in]    diss dissipation of TKE
 * @param[in]    tau13 rho<u'w'>; holds the surface-layer stress on the zlo face
 * @param[in]    tau23 rho<v'w'>; holds the surface-layer stress on the zlo face
 * @param[in]    mu_turb turbulent viscosity
 * @param[in]    solverChoice container of solver and diffusion parameters
 * @param[in]    level AMR level
 * @param[in]    tm_arr theta mean array
 * @param[in]    grav_gpu gravity vector
 * @param[in]    bc_ptr container with boundary conditions
 * @param[in]    use_SurfLayer whether surface-layer fluxes are applied in the diffusion
 * @param[in]    SurfLayer surface layer objects, one per face
 * @param[in]    pbl_derivative_dz_inv inverse vertical spacings for the mesh type
 */
template <typename PBLDerivativeDzInvFunctor>
void
AddTurbKESources (const amrex::Box& bx, const amrex::Box& domain,
                  const amrex::Array4<const amrex::Real>& u,
                  const amrex::Array4<const amrex::Real>& v,
                  const amrex::Array4<const amrex::Real>& cell_data,
                  const amrex::Array4<const amrex::Real>& cell_prim,
                  const amrex::Array4<amrex::Real>& cell_rhs,
                  const amrex::Array4<const amrex::Real>& SmnSmn_a,
                  const amrex::Array4<const amrex::Real>& hfx_z,
                  const amrex::Array4<const amrex::Real>& qfx1_z,
                  const amrex::Array4<const amrex::Real>& diss,
                  const amrex::Array4<const amrex::Real>& tau13,
                  const amrex::Array4<const amrex::Real>& tau23,
                  const amrex::Array4<const amrex::Real>& mu_turb,
                  const SolverChoice& solverChoice,
                  const int level,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                  const amrex::BCRec* bc_ptr,
                  const bool use_SurfLayer,
                  const amrex::Vector<std::unique_ptr<SurfaceLayer>>& SurfLayer,
                  const PBLDerivativeDzInvFunctor& pbl_derivative_dz_inv)
{
    BL_PROFILE("AddTurbKESources()");

    const TurbChoice& turbChoice = solverChoice.turbChoice[level];

    const bool l_use_keqn      = turbChoice.use_keqn;
    const bool l_implicit_diss = turbChoice.use_keqn && turbChoice.implicit_tke_dissipation;
    const bool l_use_mynn      = turbChoice.use_pbl_tke;

    const bool SurfLayer_zlo = use_SurfLayer && SurfLayer[amrex::Orientation::zlo()] != nullptr;

    // Clamp the divisor to the smallest normal Real before the select: an optimised
    //   build if-converts the select into an unconditional 1/theta_ref (also when
    //   the select is on the divisor), which trips fpe_trap_zero with the default
    //   theta_ref = 0. The clamp is exact for any normal theta_ref > 0 and its
    //   reciprocal is finite in single and double precision.
    const bool use_ref_theta = (turbChoice.theta_ref > 0);
    const amrex::Real inv_theta_ref = one / amrex::max(turbChoice.theta_ref,
                                                       std::numeric_limits<amrex::Real>::min());
    const amrex::Real l_inv_theta0 = (use_ref_theta) ? inv_theta_ref : one;
    const amrex::Real l_abs_g      = std::abs(grav_gpu[2]);

    // Using Deardorff (see Sullivan et al 1994)
    //    or k-eqn RANS (see Axell & Liungman 2001)
    //
    // Note: At this point, the thermal diffusivity ("Khv" field in ERF) and the
    //       subgrid dissipation ("diss" here) have been updated by the closure
    //       at the beginning of each timestep.
    //       The strain rate magnitude is updated at the beginning of the first
    //       RK stage only, therefore the shear production term also does not
    //       change between RK stages.
    //       hfx_z (Hfx3) is z-nodal, (0,0,1): it lives at the w-faces, not at the
    //       cell centre, so a single hfx_z(i,j,k) is the flux at the LOWER face of
    //       the cell and using it as the cell-centred value is not formally correct.
    //       hfx_z holds the face fluxes of the theta diffusion of this stage
    //       (the full flux, independent of the implicit/explicit split); at the
    //       bottom face of the domain that is the surface-layer flux when a
    //       surface layer is used.  The buoyancy production therefore averages the
    //       two faces of every cell, the first cell included, for the cell-centred
    //       vertical heat flux.  The first cell is not special-cased to the surface
    //       flux alone: that would again place a w-face value at the cell centre
    //       (a MOST-profile gradient evaluated at the first cell centre would be
    //       the consistent alternative).
    //       On a terrain-fitted mesh hfx_z holds the flux along zeta only: the slope
    //       corrections and the Jacobian are applied where the fluxes enter the RHS, so
    //       over sloping terrain this source omits those terms (as does the hfx3 output).
    //       REQUIRED WRITER: the RhoTheta call of DiffusionSrcForState_{N,S,T}
    //       in erf_slow_rhs_pre, on the same RK stage, which writes every z-face
    //       of each valid box (zbx = surroundingNodes(bx,2), top face included).
    //       Running the RhoKE source without that same-stage theta pass would
    //       read stale faces.
    if (l_use_keqn) {
        int qty_index = RhoKE_comp;
        // The source reads hfx_z at k+1: the z-faces of bx must be in hfx_z
        AMREX_ASSERT(hfx_z.contains(bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)) &&
                     hfx_z.contains(bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)+1));
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Add Buoyancy Source
            // where the SGS buoyancy flux tau_{theta,i} = -KH * dtheta/dx_i,
            // such that for dtheta/dz < 0, there is a positive (upward) heat
            // flux; the TKE buoyancy production is then
            //   B = g/theta_0 * tau_{theta,w}
            // for a dry atmosphere.
            // TODO: To account for moisture, the Brunt-Vaisala frequency,
            //   N^2 = g[1/theta * dtheta/dz + ...]
            // **should** be a function of the water vapor and total water
            // mixing ratios, depending on whether conditions are saturated or
            // not (see the WRF model description, Skamarock et al 2019).
            // NOTE: reads hfx_z at k+1, one face above bx; Hfx3 is z-nodal and
            //       the theta diffusion writes every face of its valid box, top included.
            const amrex::Real hfx_cc = amrex::Real(0.5) * ( hfx_z(i,j,k) + hfx_z(i,j,k+1) );
            amrex::Real buoyancy_src = l_abs_g * l_inv_theta0 * hfx_cc;
            if (!use_ref_theta) {
                // l_inv_theta0 == 1, divide by actual theta
                buoyancy_src /= cell_prim(i,j,k,PrimTheta_comp);
            }
            cell_rhs(i,j,k,qty_index) += buoyancy_src;

            // TKE shear production
            //   P = -tau_ij * S_ij = 2 * mu_turb * S_ij * S_ij
            // Note: This assumes that the horizontal and vertical diffusivities
            // of momentum are equal
            cell_rhs(i,j,k,qty_index) += two*mu_turb(i,j,k,EddyDiff::Mom_v) * SmnSmn_a(i,j,k);

            // TKE dissipation (explicit); with erf.implicit_tke_dissipation the
            // sink is folded into the update in erf_slow_rhs_post instead
            if (!l_implicit_diss) {
                cell_rhs(i,j,k,qty_index) -= diss(i,j,k);
            }
        });
    }

    // Using PBL
    if (l_use_mynn) {
        int qty_index = RhoKE_comp;
        auto pbl_mynn_B1_l = turbChoice.pbl_mynn.B1;

        // At k = klo the resolved vertical gradients are contaminated by the
        // foextrap/hoextrap ghost cells (ERF #4037). When a SurfaceLayer owns the
        // zlo face, ComputeQKESourceTerms builds the production terms there from
        // the surface stresses and scalar fluxes that are actually applied to the
        // momentum, heat and moisture equations instead.
        const bool l_surface_layer_on_zlo = SurfLayer_zlo;

        amrex::ParallelFor(bx,[=,moist_indices=solverChoice.moisture_indices] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
            bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].hi(2) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir) );

            const auto pbl_dz_inv = pbl_derivative_dz_inv(i,j,k);

            // This computes shear production, buoyancy production, and dissipation terms only.
            cell_rhs(i, j, k, qty_index) += ComputeQKESourceTerms(i,j,k,u,v,cell_data,cell_prim,
                                                                  mu_turb,domain,pbl_dz_inv,
                                                                  pbl_mynn_B1_l,tm_arr(i,j,0),
                                                                  moist_indices,
                                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                                                  u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                                                  v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                                                  tau13, tau23, hfx_z, qfx1_z,
                                                                  l_surface_layer_on_zlo);
        });
    }
}

#endif
