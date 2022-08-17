#include <DiffusionFluxForState.H>
#include <StressTerm.H>

using namespace amrex;

void
DiffusionSrcForState(const amrex::Box& bx, const amrex::Box& domain, int n_start, int n_end,
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
                     const amrex::Array4<amrex::Real>& K_turb,
                     const SolverChoice &solverChoice,
                     const amrex::Real& theta_mean,
                     const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> grav_gpu,
                     const amrex::BCRec* bc_ptr)
{
  if ( (solverChoice.molec_diff_type == MolecDiffType::None) &&
       (solverChoice.les_type        ==       LESType::None) &&
       (solverChoice.pbl_type        ==       PBLType::None) ) {
          return;
  } else {
      const amrex::Real dx_inv = cellSizeInv[0];
      const amrex::Real dy_inv = cellSizeInv[1];
      const amrex::Real dz_inv = cellSizeInv[2];

      bool l_use_QKE       = solverChoice.use_QKE && solverChoice.advect_QKE;
      bool l_use_deardorff = (solverChoice.les_type == LESType::Deardorff);
      Real l_Delta         = std::pow(dx_inv * dy_inv * dz_inv,-1./3.);
      Real l_C_e           = solverChoice.Ce;

      const int l_use_terrain = solverChoice.use_terrain;

      const Box xbx = surroundingNodes(bx,0);
      const Box ybx = surroundingNodes(bx,1);
      const Box zbx = surroundingNodes(bx,2);

      const int ncomp      = n_end - n_start + 1;
      const int qty_offset = RhoTheta_comp;

      amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
          const int  qty_index = n_start + n;
          const int prim_index = qty_index - qty_offset;
          xflux(i,j,k,qty_index) = ComputeDiffusionFluxForState(i, j, k, cell_data, cell_prim, prim_index,
                                                                dx_inv, K_turb, solverChoice, Coord::x);
      });
      amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
          const int  qty_index = n_start + n;
          const int prim_index = qty_index - qty_offset;
          yflux(i,j,k,qty_index) = ComputeDiffusionFluxForState(i, j, k, cell_data, cell_prim, prim_index,
                                                                dy_inv, K_turb, solverChoice, Coord::y);
      });
      amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
          const int  qty_index = n_start + n;
          const int prim_index = qty_index - qty_offset;
          zflux(i,j,k,qty_index) = ComputeDiffusionFluxForState(i, j, k, cell_data, cell_prim, prim_index,
                                                                dz_inv, K_turb, solverChoice, Coord::z);
      });

      for (int qty_index = n_start; qty_index <= n_end; qty_index++)
      {
          amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              cell_rhs(i, j, k, qty_index) +=
                   (xflux(i+1,j  ,k  ,qty_index) - xflux(i, j, k, qty_index)) * dx_inv   // Diffusive flux in x-dir
                  +(yflux(i  ,j+1,k  ,qty_index) - yflux(i, j, k, qty_index)) * dy_inv   // Diffusive flux in y-dir
                  +(zflux(i  ,j  ,k+1,qty_index) - zflux(i, j, k, qty_index)) * dz_inv;  // Diffusive flux in z-dir

              if (l_use_deardorff && qty_index == RhoKE_comp)
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
                  cell_rhs(i, j, k, qty_index) += cell_data(i,j,k,Rho_comp) * grav_gpu[2] * KH * dtheta_dz;

                  // Add TKE production
                  cell_rhs(i, j, k, qty_index) += ComputeTKEProduction(i,j,k,u,v,w,K_turb,cellSizeInv,domain,bc_ptr,l_use_terrain);

                  // Add dissipation
                  if (std::abs(E) > 0.) {
                      cell_rhs(i, j, k, qty_index) += cell_data(i,j,k,Rho_comp) * l_C_e *
                          std::pow(E,1.5) / length;
                  }
              }

              // QKE : similar terms to TKE
              if (l_use_QKE && qty_index == RhoQKE_comp) {
                  cell_rhs(i, j, k, qty_index) += ComputeQKESourceTerms(i,j,k,u,v,cell_data,cell_prim,
                                                                        K_turb,cellSizeInv,domain,solverChoice,theta_mean);
              }

              // Add source terms. TODO: Put this under an if condition when we implement source term
              cell_rhs(i, j, k, qty_index) += source_fab(i, j, k, qty_index);
          });
      } // qty_index
  }
}
