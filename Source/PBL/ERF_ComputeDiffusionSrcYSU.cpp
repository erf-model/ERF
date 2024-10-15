#include "ERF_ABLMost.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"

using namespace amrex;

void
DiffusionSrcForStateYSU (const Box& bx, const Box& domain,
                         int start_comp, int num_comp,
                         const bool& exp_most,
                         const bool& rot_most,
                         const Array4<const Real>& u,
                         const Array4<const Real>& v,
                         const Array4<const Real>& cell_data,
                         const Array4<const Real>& cell_prim,
                         const Array4<Real>& cell_rhs,
                         const Array4<Real>& /*xflux*/,
                         const Array4<Real>& /*yflux*/,
                         const Array4<Real>& zflux,
                         const Array4<const Real>& z_nd,
                         const Array4<const Real>& /*ax*/,
                         const Array4<const Real>& /*ay*/,
                         const Array4<const Real>& az,
                         const Array4<const Real>& detJ,
                         const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                         const Array4<const Real>& SmnSmn_a,
                         const Array4<const Real>& mf_m,
                         const Array4<const Real>& mf_u,
                         const Array4<const Real>& mf_v,
                               Array4<      Real>& /*hfx_x*/,
                               Array4<      Real>& /*hfx_y*/,
                               Array4<      Real>& hfx_z,
                               Array4<      Real>& /*qfx1_x*/,
                               Array4<      Real>& /*qfx1_y*/,
                               Array4<      Real>& /*qfx1_z*/,
                               Array4<      Real>& /*qfx2_z*/,
                               Array4<      Real>& diss,
                         const Array4<const Real>& mu_turb,
                         const SolverChoice &solverChoice,
                         const int level,
                         const Array4<const Real>& tm_arr,
                         const GpuArray<Real,AMREX_SPACEDIM> grav_gpu,
                         const BCRec* bc_ptr,
                         const bool use_most)
{
    BL_PROFILE_VAR("DiffusionSrcForStateYSU()",DiffusionSrcForStateYSU);

    // YSU SRC: Diffusion only in the z-direction
    Print() << "Computing YSU diffusion SRC for state" << std::endl;

    // For now, only theta is diffused. TODO: Moisture
    int eddy_diff_idz[1] {EddyDiff::Theta_v};
    int entrainment_idz[1] {EddyDiff::Theta_ent_YSU};
    if (num_comp != 1 ) {
        Abort("DiffusionSrcForStateYSU(): num_comp must be 1");
    } else if ( start_comp != RhoTheta_comp) {
        Abort("DiffusionSrcForStateYSU(): start_comp must be RhoTheta_comp");
    }

    // Geometry preparations
    const Box zbx = surroundingNodes(bx,2);
    const auto& dom_hi = ubound(domain);
    const Real dz_inv = cellSizeInv[2];

    ParallelFor(zbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - 1;
            const int diff_idx  = eddy_diff_idz[prim_index];
            const int entrainment_idx =  eddy_diff_idz[prim_index];

            const amrex::Real rhoAlpha = 0.5 * ( mu_turb(i, j, k  , diff_idx)
                                               + mu_turb(i, j, k-1, diff_idx) );

            const amrex::Real entrainment = 0.5 * ( mu_turb(i, j, k  , entrainment_idx)
                                                  + mu_turb(i, j, k-1, entrainment_idx) );

            Real met_h_zeta = az(i,j,k);

            Real GradCz;
            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == 0);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool most_on_zlo    = ( use_most && exp_most &&
                                  (bc_ptr[bc_comp].lo(2) == ERFBCType::foextrap) && k == 0);

            if (ext_dir_on_zlo) {
                GradCz = dz_inv * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                        + 3. * cell_prim(i, j, k  , prim_index)
                                   - (1./3.) * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                GradCz = dz_inv * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                        - 3. * cell_prim(i, j, k-1, prim_index)
                                   + (1./3.) * cell_prim(i, j, k-2, prim_index) );
            } else {
                GradCz = dz_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            // Always True for now
            if (qty_index == RhoTheta_comp) {
                if (!most_on_zlo) {
                    zflux(i,j,k,qty_index) = -rhoAlpha * GradCz / met_h_zeta - entrainment;
                    hfx_z(i,j,k) = zflux(i,j,k,qty_index);
                } else {
                    zflux(i,j,k,qty_index) = hfx_z(i,j,0);
                }
            }
        });


    // Use fluxes to compute RHS
    //-----------------------------------------------------------------------------------
    for (int n(0); n < num_comp; ++n)
    {
        int qty_index = start_comp + n;
        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real stateContrib = (zflux(i  ,j  ,k+1,qty_index) - zflux(i, j, k, qty_index)) * dz_inv;  // Diffusive flux in z-dir

            stateContrib /= detJ(i,j,k);

            cell_rhs(i,j,k,qty_index) -= stateContrib;
        });
    }
}

void
DiffusionSrcForMomYSU (const Box& bxcc, const Box& bxxz , const Box& bxyz,
                       const Array4<const Real>& uvel ,
                       const Array4<const Real>& vvel ,
                       const Array4<const Real>& wvel ,
                       const Array4<Real>& tau33,
                       const Array4<Real>& tau13, const Array4<Real>& tau23,
                       const Array4<Real>& tau31, const Array4<Real>& tau32,
                       const Array4<const Real>& detJ ,
                       const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                       const Array4<const Real>& mf_m,
                       const Array4<const Real>& /*mf_u*/,
                       const Array4<const Real>& /*mf_v*/)
{
    BL_PROFILE_VAR("DiffusionSrcForMomYSU",DiffusionSrcForMomYSU);
    Print() << "Computing YSU diffusion SRC for mom" << std::endl;

    // YSU PBL: always only vertical diffusion
    // For now - only diffuse x and y velocity
}
