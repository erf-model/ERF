#ifndef _COMPUTE_EXPANSIONRATE_H_
#define _COMPUTE_EXPANSIONRATE_H_

#ifdef ERF_USE_TERRAIN

#include <TerrainMetrics.H>
#include <IndexDefines.H>

AMREX_GPU_DEVICE
inline amrex::Real
ComputeExpansionRate(const int &i, const int &j, const int &k,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     const enum MomentumEqn &momentumEqn,
                     const enum DiffusionDir &diffDir,
                     const amrex::Array4<const amrex::Real>& z_nd,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv)
{
    amrex::Real dx_inv = cellSizeInv[0];
    amrex::Real dy_inv = cellSizeInv[1];
    amrex::Real dz_inv = cellSizeInv[2];

    amrex::Real GradU,GradV,GradOmega;
    amrex::Real Omega_hi,Omega_lo;

    amrex::Real met_h_xi,met_h_eta,met_h_zeta;
     
    amrex::Real met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi;
    amrex::Real met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo;

    amrex::Real met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi;
    amrex::Real met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo;
   
    amrex::Real expansionRate = 0;

    switch (momentumEqn) {
        case MomentumEqn::x:
            switch (diffDir) {
                case DiffusionDir::x: // D11
                  ComputeMetricAtIface(i  ,j  ,k  ,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i-1,j  ,k  ,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);

                  ComputeMetricAtJface(i-1,j+1,k  ,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i-1,j  ,k  ,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);

                  Omega_hi = OmegaFromW(i-1,j  ,k+1, w(i-1,j  ,k+1),u,v,z_nd,cellSizeInv);
                  Omega_lo = OmegaFromW(i-1,j  ,k  , w(i-1,j  ,k  ),u,v,z_nd,cellSizeInv);
                  
                  // D11 (i+1/2)
                  expansionRate = (u(i  ,j  ,k  )*met_u_h_zeta_hi - u(i-1,j  ,k  )*met_u_h_zeta_lo)*dx_inv +
                                  (v(i-1,j+1,k  )*met_v_h_zeta_hi - v(i-1,j  ,k  )*met_v_h_zeta_lo)*dy_inv +
                                  (Omega_hi - Omega_lo)*dz_inv;
                  break;
                case DiffusionDir::y: // D12
                  expansionRate = 0.0;
                  break;
                case DiffusionDir::z: // D13
                  expansionRate = 0.0;
                  break;
                default:
                  amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        case MomentumEqn::y:
            switch (diffDir) {
                case DiffusionDir::x: // D21
                  expansionRate = 0.0;
                  break;
                case DiffusionDir::y: // D22
                  ComputeMetricAtIface(i+1,j-1,k  ,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i  ,j-1,k  ,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);

                  ComputeMetricAtJface(i  ,j  ,k  ,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i  ,j-1,k  ,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);

                  Omega_hi = OmegaFromW(i  ,j-1,k+1, w(i  ,j-1,k+1),u,v,z_nd,cellSizeInv);
                  Omega_lo = OmegaFromW(i  ,j-1,k  , w(i  ,j-1,k  ),u,v,z_nd,cellSizeInv);
                  
                  expansionRate = (u(i+1,j-1,k  )*met_u_h_zeta_hi - u(i  ,j-1,k  )*met_u_h_zeta_lo)*dx_inv +
                                  (v(i  ,j  ,k  )*met_v_h_zeta_hi - v(i  ,j-1,k  )*met_v_h_zeta_lo)*dy_inv +
                                  (Omega_hi - Omega_lo)*dz_inv;
                  break;
                case DiffusionDir::z: // D23
                  expansionRate = 0.0;
                  break;
                default:
                  amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        case MomentumEqn::z:
            switch (diffDir) {
                case DiffusionDir::x: // D31

                  ComputeMetricAtIface(i+1,j  ,k  ,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i-1,j  ,k  ,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradU = ( u(i+1,j  ,k  )*met_u_h_zeta_hi - u(i-1,j  ,k  )*met_u_h_zeta_lo )*0.5*dx_inv;
                  
                  ComputeMetricAtIface(i+1,j  ,k-1,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i-1,j  ,k-1,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradU += ( u(i+1,j  ,k-1)*met_u_h_zeta_hi - u(i-1,j  ,k-1)*met_u_h_zeta_lo )*0.5*dx_inv;
                  GradU *= 0.5;
                  //**************************************************************************************************
                  ComputeMetricAtJface(i  ,j+1,k  ,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i  ,j  ,k  ,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradV = ( v(i  ,j+1,k  )*met_v_h_zeta_hi - v(i  ,j  ,k  )*met_v_h_zeta_lo )*dy_inv;
                  
                  ComputeMetricAtJface(i  ,j+1,k-1,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i  ,j  ,k-1,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradV += ( v(i  ,j+1,k-1)*met_v_h_zeta_hi - v(i  ,j  ,k-1)*met_v_h_zeta_lo )*dy_inv;
                  
                  ComputeMetricAtJface(i-1,j+1,k  ,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i-1,j  ,k  ,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradV += ( v(i-1,j+1,k  )*met_v_h_zeta_hi - v(i-1,j  ,k  )*met_v_h_zeta_lo )*dy_inv;

                  ComputeMetricAtJface(i-1,j+1,k-1,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i-1,j  ,k-1,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradV += ( v(i-1,j+1,k-1)*met_v_h_zeta_hi - v(i-1,j  ,k-1)*met_v_h_zeta_lo )*dy_inv;
                  GradV *= 0.25;
                  //**************************************************************************************************
                  Omega_hi   = OmegaFromW(i  ,j  ,k+1, w(i  ,j  ,k+1),u,v,z_nd,cellSizeInv);
                  Omega_lo   = OmegaFromW(i  ,j  ,k-1, w(i  ,j  ,k-1),u,v,z_nd,cellSizeInv);
                  GradOmega  = (Omega_hi - Omega_lo)*0.5*dz_inv;

                  Omega_hi   = OmegaFromW(i-1,j  ,k+1, w(i-1,j  ,k+1),u,v,z_nd,cellSizeInv);
                  Omega_lo   = OmegaFromW(i-1,j  ,k-1, w(i-1,j  ,k-1),u,v,z_nd,cellSizeInv);
                  GradOmega += (Omega_hi - Omega_lo)*0.5*dz_inv;
                  GradOmega *= 0.5;
                  //**************************************************************************************************
                  ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,
                                            met_h_xi,met_h_eta,met_h_zeta,
                                            cellSizeInv,z_nd,TerrainMet::h_xi_zeta);
                  
                  expansionRate = -(met_h_xi/met_h_zeta)*(GradU + GradV + GradOmega);                  
                  break;
                case DiffusionDir::y: // D32

                  ComputeMetricAtIface(i+1,j  ,k  ,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i  ,j  ,k  ,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradU = ( u(i+1,j  ,k  )*met_u_h_zeta_hi - u(i  ,j  ,k  )*met_u_h_zeta_lo )*dx_inv;

                  ComputeMetricAtIface(i+1,j  ,k-1,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i  ,j  ,k-1,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradU += ( u(i+1,j  ,k-1)*met_u_h_zeta_hi - u(i  ,j  ,k-1)*met_u_h_zeta_lo )*dx_inv;

                  ComputeMetricAtIface(i+1,j-1,k  ,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i  ,j-1,k  ,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradU += ( u(i+1,j-1,k  )*met_u_h_zeta_hi - u(i  ,j-1,k  )*met_u_h_zeta_lo )*dx_inv;

                  ComputeMetricAtIface(i+1,j-1,k-1,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i  ,j-1,k-1,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradU += ( u(i+1,j-1,k-1)*met_u_h_zeta_hi - u(i  ,j-1,k-1)*met_u_h_zeta_lo )*dx_inv;
                  GradU *= 0.25;
                  //**************************************************************************************************
                  ComputeMetricAtJface(i  ,j+1,k  ,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i  ,j-1,k  ,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradV = ( v(i  ,j+1,k  )*met_v_h_zeta_hi - v(i  ,j-1,k  )*met_v_h_zeta_lo )*0.5*dy_inv;
                  
                  ComputeMetricAtJface(i  ,j+1,k-1,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i  ,j-1,k-1,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  GradV += ( v(i  ,j+1,k-1)*met_v_h_zeta_hi - v(i  ,j-1,k-1)*met_v_h_zeta_lo )*0.5*dy_inv;
                  GradV *= 0.5;
                  //**************************************************************************************************
                  Omega_hi   = OmegaFromW(i  ,j  ,k+1, w(i  ,j  ,k+1),u,v,z_nd,cellSizeInv);
                  Omega_lo   = OmegaFromW(i  ,j  ,k-1, w(i  ,j  ,k-1),u,v,z_nd,cellSizeInv);
                  GradOmega  = (Omega_hi - Omega_lo)*0.5*dz_inv;

                  Omega_hi   = OmegaFromW(i  ,j-1,k+1, w(i  ,j-1,k+1),u,v,z_nd,cellSizeInv);
                  Omega_lo   = OmegaFromW(i  ,j-1,k-1, w(i  ,j-1,k-1),u,v,z_nd,cellSizeInv);
                  GradOmega += (Omega_hi - Omega_lo)*0.5*dz_inv;
                  GradOmega *= 0.5;
                  //**************************************************************************************************
                  ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,
                                            met_h_xi,met_h_eta,met_h_zeta,
                                            cellSizeInv,z_nd,TerrainMet::h_eta_zeta);
                  
                  expansionRate = -(met_h_eta/met_h_zeta)*(GradU + GradV + GradOmega);
                  break;
                case DiffusionDir::z: // D33
                  ComputeMetricAtIface(i+1,j  ,k-1,
                                       met_u_h_xi_hi,met_u_h_eta_hi,met_u_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtIface(i  ,j  ,k-1,
                                       met_u_h_xi_lo,met_u_h_eta_lo,met_u_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  
                  ComputeMetricAtJface(i  ,j+1,k-1,
                                       met_v_h_xi_hi,met_v_h_eta_hi,met_v_h_zeta_hi,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  ComputeMetricAtJface(i  ,j  ,k-1,
                                       met_v_h_xi_lo,met_v_h_eta_lo,met_v_h_zeta_lo,
                                       cellSizeInv,z_nd,TerrainMet::h_zeta);
                  
                  Omega_hi = OmegaFromW(i  ,j  ,k  , w(i  ,j  ,k  ),u,v,z_nd,cellSizeInv);
                  Omega_lo = OmegaFromW(i  ,j  ,k-1, w(i  ,j  ,k-1),u,v,z_nd,cellSizeInv);
                  
                  expansionRate = (u(i+1,j  ,k-1)*met_u_h_zeta_hi - u(i  ,j  ,k-1)*met_u_h_zeta_lo)*dx_inv +
                                  (v(i  ,j+1,k-1)*met_v_h_zeta_hi - v(i  ,j  ,k-1)*met_v_h_zeta_lo)*dy_inv +
                                  (Omega_hi - Omega_lo)*dz_inv;
                  
                  ComputeMetricAtCellCenter(i  ,j  ,k-1,
                                            met_h_xi,met_h_eta,met_h_zeta,
                                            cellSizeInv,z_nd,TerrainMet::h_zeta);

                  expansionRate /= met_h_zeta;
                  
                  break;
                default:
                  amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        default:
          amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return (1.0/3.0) * expansionRate;
}

#else
AMREX_GPU_DEVICE
inline amrex::Real
ComputeExpansionRate(const int &i, const int &j, const int &k,
                     const amrex::Array4<const amrex::Real>& u,
                     const amrex::Array4<const amrex::Real>& v,
                     const amrex::Array4<const amrex::Real>& w,
                     const enum MomentumEqn &momentumEqn,
                     const enum DiffusionDir &diffDir,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv) {
    amrex::Real dx_inv = cellSizeInv[0];
    amrex::Real dy_inv = cellSizeInv[1];
    amrex::Real dz_inv = cellSizeInv[2];

    amrex::Real expansionRate = 0;

    switch (momentumEqn) {
        case MomentumEqn::x:
            switch (diffDir) {
                case DiffusionDir::x: // D11
                    // D11 (i+1/2)
                    expansionRate = (u(i  , j  , k  ) - u(i-1, j, k))*dx_inv +
                                    (v(i-1, j+1, k  ) - v(i-1, j, k))*dy_inv +
                                    (w(i-1, j  , k+1) - w(i-1, j, k))*dz_inv;
                    break;
                case DiffusionDir::y: // D12
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::z: // D13
                    expansionRate = 0.0;
                    break;
                default:
                    amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        case MomentumEqn::y:
            switch (diffDir) {
                case DiffusionDir::x: // D21
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::y: // D22
                    expansionRate = (u(i+1, j-1, k  ) - u(i, j-1, k))*dx_inv +
                                    (v(i  , j  , k  ) - v(i, j-1, k))*dy_inv +
                                    (w(i  , j-1, k+1) - w(i, j-1, k))*dz_inv;
                    break;
                case DiffusionDir::z: // D23
                    expansionRate = 0.0;
                    break;
                default:
                    amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        case MomentumEqn::z:
            switch (diffDir) {
                case DiffusionDir::x: // D31
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::y: // D32
                    expansionRate = 0.0;
                    break;
                case DiffusionDir::z: // D33
                    expansionRate = (u(i+1, j  , k-1) - u(i, j, k-1))*dx_inv +
                                    (v(i  , j+1, k-1) - v(i, j, k-1))*dy_inv +
                                    (w(i  , j  , k  ) - w(i, j, k-1))*dz_inv;
                    break;
                default:
                    amrex::Abort("Error: Diffusion direction is unrecognized");
            }
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return (1.0/3.0) * expansionRate;
}
#endif

#endif
