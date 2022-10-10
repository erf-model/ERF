#include <Diffusion.H>
#include <EddyViscosity.H>
#include <ComputeQKESourceTerm.H>

using namespace amrex;

void
DiffusionSrcForState_N (const amrex::Box& bx, const amrex::Box& domain, int n_start, int n_end,
                        const amrex::Array4<const amrex::Real>& u,
                        const amrex::Array4<const amrex::Real>& v,
                        const amrex::Array4<const amrex::Real>& w,
                        const amrex::Array4<const amrex::Real>& cell_data,
                        const amrex::Array4<const amrex::Real>& cell_prim,
                        const amrex::Array4<const amrex::Real>& source_fab,
                        const amrex::Array4<amrex::Real>& cell_rhs,
                        const amrex::Array4<amrex::Real>& xflux,
                        const amrex::Array4<amrex::Real>& yflux,
                        const amrex::Array4<amrex::Real>& zflux,
                        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                        const amrex::Array4<const amrex::Real>& K_turb,
                        const SolverChoice &solverChoice,
                        const amrex::Real& theta_mean,
                        const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                        const amrex::BCRec* bc_ptr)
{
    BL_PROFILE_VAR("DiffusionSrcForState_N()",DiffusionSrcForState_N);

    const amrex::Real dx_inv = cellSizeInv[0];
    const amrex::Real dy_inv = cellSizeInv[1];
    const amrex::Real dz_inv = cellSizeInv[2];

    bool l_use_QKE       = solverChoice.use_QKE && solverChoice.advect_QKE;
    bool l_use_deardorff = (solverChoice.les_type == LESType::Deardorff);
    Real l_Delta         = std::pow(dx_inv * dy_inv * dz_inv,-1./3.);
    Real l_C_e           = solverChoice.Ce;

    bool l_consA  = (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = ( (solverChoice.les_type == LESType::Smagorinsky) ||
                      (solverChoice.les_type == LESType::Deardorff  ) ||
                      (solverChoice.pbl_type == PBLType::MYNN25     ) );

    int l_use_terrain = solverChoice.use_terrain;

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);

    const int ncomp      = n_end - n_start + 1;
    const int qty_offset = RhoTheta_comp;

    // Theta, KE, QKE, Scalar
    Vector<Real> alpha_eff;
    if (l_consA) {
        alpha_eff = {solverChoice.alpha_T   , 0., 0., solverChoice.alpha_C   };
    } else {
        alpha_eff = {solverChoice.rhoAlpha_T, 0., 0., solverChoice.rhoAlpha_C};
    }
    Vector<int> eddy_diff_idx{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::QKE_h, EddyDiff::Scalar_h};
    Vector<int> eddy_diff_idy{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::QKE_h, EddyDiff::Scalar_h};
    Vector<int> eddy_diff_idz{EddyDiff::Theta_v, EddyDiff::KE_v, EddyDiff::QKE_v, EddyDiff::Scalar_v};

    // Device vectors
    Gpu::AsyncVector<Real> alpha_eff_d;
    Gpu::AsyncVector<int>  eddy_diff_idx_d,eddy_diff_idy_d,eddy_diff_idz_d;
    alpha_eff_d.resize(alpha_eff.size());
    eddy_diff_idx_d.resize(eddy_diff_idx.size());
    eddy_diff_idy_d.resize(eddy_diff_idy.size());
    eddy_diff_idz_d.resize(eddy_diff_idz.size());

    Gpu::copy(Gpu::hostToDevice, alpha_eff.begin(), alpha_eff.end(), alpha_eff_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idx.begin(), eddy_diff_idx.end(), eddy_diff_idx_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idy.begin(), eddy_diff_idy.end(), eddy_diff_idy_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idz.begin(), eddy_diff_idz.end(), eddy_diff_idz_d.begin());

    // Capture pointers for device code
    Real* d_alpha_eff     = alpha_eff_d.data();
    int*  d_eddy_diff_idx = eddy_diff_idx_d.data();
    int*  d_eddy_diff_idy = eddy_diff_idy_d.data();
    int*  d_eddy_diff_idz = eddy_diff_idz_d.data();

    // Compute fluxes at each face
    if (l_consA && l_turb) {
        amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            amrex::Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( K_turb(i  , j, k, d_eddy_diff_idx[prim_index])
                              + K_turb(i-1, j, k, d_eddy_diff_idx[prim_index]) );

            xflux(i,j,k,qty_index) = rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i-1, j, k, prim_index)) * dx_inv;
        });
        amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            amrex::Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( K_turb(i, j  , k, d_eddy_diff_idy[prim_index])
                              + K_turb(i, j-1, k, d_eddy_diff_idy[prim_index]) );

            yflux(i,j,k,qty_index) = rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv;
        });
        amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            amrex::Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( K_turb(i, j, k  , d_eddy_diff_idz[prim_index])
                              + K_turb(i, j, k-1, d_eddy_diff_idz[prim_index]) );

            zflux(i,j,k,qty_index) = rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
        });
    } else if (l_turb) {
        amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real Alpha = d_alpha_eff[prim_index];
            Alpha += 0.5 * ( K_turb(i  , j, k, d_eddy_diff_idx[prim_index])
                           + K_turb(i-1, j, k, d_eddy_diff_idx[prim_index]) );

            xflux(i,j,k,qty_index) = Alpha * (cell_prim(i, j, k, prim_index) - cell_prim(i-1, j, k, prim_index)) * dx_inv;
        });
        amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real Alpha = d_alpha_eff[prim_index];
            Alpha += 0.5 * ( K_turb(i, j  , k, d_eddy_diff_idy[prim_index])
                           + K_turb(i, j-1, k, d_eddy_diff_idy[prim_index]) );

            yflux(i,j,k,qty_index) = Alpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv;
        });
        amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real Alpha = d_alpha_eff[prim_index];
            Alpha += 0.5 * ( K_turb(i, j, k  , d_eddy_diff_idz[prim_index])
                           + K_turb(i, j, k-1, d_eddy_diff_idz[prim_index]) );

            zflux(i,j,k,qty_index) = Alpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
        });
    } else if(l_consA) {
        amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            amrex::Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            xflux(i,j,k,qty_index) = rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i-1, j, k, prim_index)) * dx_inv;
        });
        amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            amrex::Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            yflux(i,j,k,qty_index) = rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv;
        });
        amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            amrex::Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            zflux(i,j,k,qty_index) = rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
        });
    } else {
        amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real Alpha = d_alpha_eff[prim_index];

            xflux(i,j,k,qty_index) = Alpha * (cell_prim(i, j, k, prim_index) - cell_prim(i-1, j, k, prim_index)) * dx_inv;
        });
        amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real Alpha = d_alpha_eff[prim_index];

            yflux(i,j,k,qty_index) = Alpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv;
        });
        amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = n_start + n;
            const int prim_index = qty_index - qty_offset;

            amrex::Real Alpha = d_alpha_eff[prim_index];

            zflux(i,j,k,qty_index) = Alpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
        });
    }


    // Use fluxes to compute RHS
    for (int qty_index = n_start; qty_index <= n_end; qty_index++)
    {
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cell_rhs(i,j,k,qty_index) += (xflux(i+1,j  ,k  ,qty_index) - xflux(i, j, k, qty_index)) * dx_inv   // Diffusive flux in x-dir
                                        +(yflux(i  ,j+1,k  ,qty_index) - yflux(i, j, k, qty_index)) * dy_inv   // Diffusive flux in y-dir
                                        +(zflux(i  ,j  ,k+1,qty_index) - zflux(i, j, k, qty_index)) * dz_inv;  // Diffusive flux in z-dir

            // Add source terms. TODO: Put this under an if condition when we implement source term
            cell_rhs(i,j,k,qty_index) += source_fab(i,j,k,qty_index);
        });
    }

    // Using Deardorff
    if (l_use_deardorff && n_end >= RhoKE_comp) {
        int qty_index = RhoKE_comp;
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Add Buoyancy Source
            Real theta     = cell_prim(i,j,k,PrimTheta_comp);
            Real dtheta_dz = 0.5*(cell_prim(i,j,k+1,PrimTheta_comp)-cell_prim(i,j,k-1,PrimTheta_comp))*dz_inv;
            Real E         = cell_prim(i,j,k,PrimKE_comp);
            Real length;
            if (dtheta_dz <= 0.) {
                length = l_Delta;
            } else {
                length = 0.76*std::sqrt(E)*(grav_gpu[2]/theta)*dtheta_dz;
            }
            Real KH   = 0.1 * (1.+2.*length/l_Delta) * std::sqrt(E);
            cell_rhs(i,j,k,qty_index) += cell_data(i,j,k,Rho_comp) * grav_gpu[2] * KH * dtheta_dz;

            // Add TKE production
            cell_rhs(i,j,k,qty_index) += ComputeTKEProduction(i,j,k,u,v,w,K_turb,cellSizeInv,domain,bc_ptr,l_use_terrain);

            // Add dissipation
            if (std::abs(E) > 0.) {
                cell_rhs(i,j,k,qty_index) += cell_data(i,j,k,Rho_comp) * l_C_e *
                    std::pow(E,1.5) / length;
            }
        });
    }

    // Using Deardorff
    if (l_use_QKE && n_end >= RhoQKE_comp) {
        int qty_index = RhoQKE_comp;
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cell_rhs(i, j, k, qty_index) += ComputeQKESourceTerms(i,j,k,u,v,cell_data,cell_prim,
                                                                  K_turb,cellSizeInv,domain,solverChoice,theta_mean);
        });
    }

}
