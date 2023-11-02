#include <ComputeQKESourceTerm.H>

void
ComputeEquilibriumQKE (amrex::Box bx_QKE,
                       const amrex::Array4<const amrex::Real>& uvel,
                       const amrex::Array4<const amrex::Real>& vvel,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const amrex::Array4<const amrex::Real>& tm_arr,
                       const amrex::Array4<const amrex::Real>& K_turb,
                       amrex::Array4<amrex::Real>& QKE_equil_arr,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                       const amrex::Box& domain,
                       const TurbChoice& turbChoice)
{
    // Constants
    const int izmin = domain.smallEnd(2);
    const amrex::Real dz_inv = cellSizeInv[2];

    const amrex::Real A1 = turbChoice.pbl_A1;
    const amrex::Real A2 = turbChoice.pbl_A2;
    const amrex::Real B1 = turbChoice.pbl_B1;

    const amrex::Real G1  = turbChoice.pbl_G1;
    const amrex::Real G2  = turbChoice.pbl_G2;
    const amrex::Real F1  = turbChoice.pbl_F1;
    const amrex::Real F2  = turbChoice.pbl_F2;
    const amrex::Real Rf1 = turbChoice.pbl_Rf1;
    const amrex::Real Rf2 = turbChoice.pbl_Rf2;
    const amrex::Real Rfc = turbChoice.pbl_Rfc;
    const amrex::Real Ri1 = turbChoice.pbl_Ri1;
    const amrex::Real Ri2 = turbChoice.pbl_Ri2;
    const amrex::Real Ri3 = turbChoice.pbl_Ri3;

    amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();

    // Store QKE equilibrium for Level 2 limiting of MYNN2.5 model
    amrex::ParallelFor(bx_QKE,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Gradients (ghost cells have been filled)
        int lk = amrex::max(k,izmin);
        amrex::Real dthetadz = 0.5*(cell_prim(i,j,lk+1,PrimTheta_comp) - cell_prim(i,j,lk-1,PrimTheta_comp))*dz_inv;
        amrex::Real dudz = 0.25*(uvel(i,j,lk+1) - uvel(i,j,lk-1) + uvel(i+1,j,lk+1) - uvel(i+1,j,lk-1))*dz_inv;
        amrex::Real dvdz = 0.25*(vvel(i,j,lk+1) - vvel(i,j,lk-1) + vvel(i,j+1,lk+1) - vvel(i,j+1,lk-1))*dz_inv;

        // MYNN2 Limiting (Nakanishi & Niino 09 DOI:10.2151/jmsj.87.895; Appendix A)
        //=========================================================================
        amrex::Real l_comb    = K_turb(i,j,k,EddyDiff::PBL_lengthscale);
        amrex::Real shearProd = dudz*dudz + dvdz*dvdz;
        amrex::Real buoyProd  = -(CONST_GRAV/tm_arr(i,j,0)) * dthetadz;
        // Gradient Richardson number
        amrex::Real Ri = -buoyProd / (shearProd + eps);
        // Flux Richardson number
        amrex::Real Rf = Ri1 * (Ri + Ri2 - std::sqrt(Ri*Ri - Ri3*Ri + Ri2*Ri2));
        // Stability functions
        amrex::Real SH2 = 3.0*A2*(G1 + G2)*(Rfc - Rf)/(1.0 - Rf);
        amrex::Real SM2 = ( (A1*F1)/(A2*F2) ) * ( (Rf1 - Rf)/(Rf2 - Rf) ) * SH2;
        // Equilibrium QKE
        QKE_equil_arr(i,j,k) = amrex::max(B1*l_comb*l_comb*SM2*(1.0-Rf)*shearProd, 0.0);
    });
}
