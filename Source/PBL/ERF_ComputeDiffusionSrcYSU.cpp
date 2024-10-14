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
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& ax,
                        const Array4<const Real>& ay,
                        const Array4<const Real>& az,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& SmnSmn_a,
                        const Array4<const Real>& mf_m,
                        const Array4<const Real>& mf_u,
                        const Array4<const Real>& mf_v,
                              Array4<      Real>& hfx_x,
                              Array4<      Real>& hfx_y,
                              Array4<      Real>& hfx_z,
                              Array4<      Real>& qfx1_x,
                              Array4<      Real>& qfx1_y,
                              Array4<      Real>& qfx1_z,
                              Array4<      Real>& qfx2_z,
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
    Print() << "Computing YSU diffusion SRC for state" << std::endl;
}

void
DiffusionSrcForMomYSU (const Box& bxx, const Box& bxy , const Box& bxz,
                      const Array4<Real>& rho_u_rhs  ,
                      const Array4<Real>& rho_v_rhs  ,
                      const Array4<Real>& rho_w_rhs  ,
                      const Array4<const Real>& tau11, const Array4<const Real>& tau22, const Array4<const Real>& tau33,
                      const Array4<const Real>& tau12, const Array4<const Real>& tau13,
                      const Array4<const Real>& tau21, const Array4<const Real>& tau23,
                      const Array4<const Real>& tau31, const Array4<const Real>& tau32,
                      const Array4<const Real>& detJ ,
                      const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                      const Array4<const Real>& mf_m,
                      const Array4<const Real>& /*mf_u*/,
                      const Array4<const Real>& /*mf_v*/)
{
    BL_PROFILE_VAR("DiffusionSrcForMomYSU",DiffusionSrcForMomYSU);
    Print() << "Computing YSU diffusion SRC for mom" << std::endl;
}
