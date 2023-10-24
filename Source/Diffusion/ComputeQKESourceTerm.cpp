#include <ComputeQKESourceTerm.H>

void
ComputeEquilibriumQKE (amrex::Box bx_QKE,
                       const amrex::Array4<const amrex::Real>& uvel,
                       const amrex::Array4<const amrex::Real>& vvel,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& tm_arr,
                       const amrex::Array4<const amrex::Real>& K_turb,
                       amrex::Array4<amrex::Real>& QKE_equil_arr,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Box& domain,
                       amrex::Real pbl_B1_l)
{
    // Store QKE equilibrium for Level 2 limiting of MYNN2.5 model
    int izmin = domain.smallEnd(2);
    amrex::Real dz_inv = cellSizeInv[2];
    amrex::ParallelFor(bx_QKE,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Gradients (ghost cells have been filled)
        int lk = amrex::max(k,izmin);
        amrex::Real dthetadz = 0.5*(cell_prim(i,j,lk+1,PrimTheta_comp) - cell_prim(i,j,lk-1,PrimTheta_comp))*dz_inv;
        amrex::Real dudz = 0.25*(uvel(i,j,lk+1) - uvel(i,j,lk-1) + uvel(i+1,j,lk+1) - uvel(i+1,j,lk-1))*dz_inv;
        amrex::Real dvdz = 0.25*(vvel(i,j,lk+1) - vvel(i,j,lk-1) + vvel(i,j+1,lk+1) - vvel(i,j+1,lk-1))*dz_inv;

        // Balance relation
        amrex::Real qke       = cell_prim(i,j,k,PrimQKE_comp);
        amrex::Real qvel      = std::sqrt(qke);
        amrex::Real l_comb    = K_turb(i,j,k,EddyDiff::PBL_lengthscale);
        amrex::Real shearProd = dudz*dudz + dvdz*dvdz;
        amrex::Real buoyProd  = -(CONST_GRAV/tm_arr(i,j,0)) * dthetadz;
        amrex::Real lSM       = K_turb(i,j,k,EddyDiff::Mom_v)   / (qvel + 1.0e-16);
        amrex::Real lSH       = K_turb(i,j,k,EddyDiff::Theta_v) / (qvel + 1.0e-16);
        QKE_equil_arr(i,j,k)  = amrex::max(pbl_B1_l * l_comb * ( lSM * shearProd + lSH * buoyProd ), 0.0);
    });
}
