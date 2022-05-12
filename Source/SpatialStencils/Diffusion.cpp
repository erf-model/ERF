#include <SpatialStencils.H>
#include <StressTerm.H>

using namespace amrex;

#ifdef ERF_USE_TERRAIN

AMREX_GPU_DEVICE
Real
DiffusionSrcForMom(const int &i, const int &j, const int &k,
                   const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                   const Array4<const Real>& cons,
                   const enum MomentumEqn &momentumEqn,
                   const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                   const Array4<Real>& K_turb,
                   const SolverChoice &solverChoice,
                   const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                   const Box& domain, const amrex::BCRec* bc_ptr)
{
  auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
  Real diffContrib = 0.0;

  int l_spatial_order = solverChoice.spatial_order;

  switch (momentumEqn) {
  case MomentumEqn::x:
    Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
    
    Real met_h_xi_11_hi,met_h_eta_11_hi,met_h_zeta_11_hi;
    Real met_h_xi_11_lo,met_h_eta_11_lo,met_h_zeta_11_lo;

    Real met_h_xi_12_hi,met_h_eta_12_hi,met_h_zeta_12_hi;
    Real met_h_xi_12_lo,met_h_eta_12_lo,met_h_zeta_12_lo;

    Real met_h_xi_13_hi,met_h_eta_13_hi,met_h_zeta_13_hi;
    Real met_h_xi_13_lo,met_h_eta_13_lo,met_h_zeta_13_lo;
    
    tau11Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);

    ComputeMetricAtCellCenter(i  ,j,k,met_h_xi_11_hi,met_h_eta_11_hi,met_h_zeta_11_hi,
                              cellSizeInv,z_nd,TerrainMet::h_zeta);
    ComputeMetricAtCellCenter(i-1,j,k,met_h_xi_11_lo,met_h_eta_11_lo,met_h_zeta_11_lo,
                              cellSizeInv,z_nd,TerrainMet::h_zeta);

    ComputeMetricAtEdgeCenterK(i,j+1,k,met_h_xi_12_hi,met_h_eta_12_hi,met_h_zeta_12_hi,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);
    ComputeMetricAtEdgeCenterK(i,j  ,k,met_h_xi_12_lo,met_h_eta_12_lo,met_h_zeta_12_lo,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);

    ComputeMetricAtEdgeCenterJ(i,j+1,k,met_h_xi_13_hi,met_h_eta_13_hi,met_h_zeta_13_hi,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);
    ComputeMetricAtEdgeCenterJ(i,j  ,k,met_h_xi_13_lo,met_h_eta_13_lo,met_h_zeta_13_lo,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);

    diffContrib = (tau11Next*met_h_zeta_11_hi - tau11Prev*met_h_zeta_11_lo) * dxInv  // Contribution to x-mom eqn from diffusive flux in x-dir
                + (tau12Next*met_h_zeta_12_hi - tau12Prev*met_h_zeta_12_lo) * dyInv  // Contribution to x-mom eqn from diffusive flux in y-dir
                + (tau13Next*met_h_zeta_13_hi - tau13Prev*met_h_zeta_13_lo) * dzInv; // Contribution to x-mom eqn from diffusive flux in z-dir

    diffContrib /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k)); // Terrain grid stretching

    if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
    {
      diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, u(i,j,k), Coord::x, l_spatial_order) /
        solverChoice.rho0_trans;
    }
    break;
  case MomentumEqn::y:
    Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;

    Real met_h_xi_22_hi,met_h_eta_22_hi,met_h_zeta_22_hi;
    Real met_h_xi_22_lo,met_h_eta_22_lo,met_h_zeta_22_lo;

    Real met_h_xi_21_hi,met_h_eta_21_hi,met_h_zeta_21_hi;
    Real met_h_xi_21_lo,met_h_eta_21_lo,met_h_zeta_21_lo;

    Real met_h_xi_23_hi,met_h_eta_23_hi,met_h_zeta_23_hi;
    Real met_h_xi_23_lo,met_h_eta_23_lo,met_h_zeta_23_lo;
    
    tau21Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);

    ComputeMetricAtCellCenter(i,j  ,k,met_h_xi_22_hi,met_h_eta_22_hi,met_h_zeta_22_hi,
                              cellSizeInv,z_nd,TerrainMet::h_zeta);
    ComputeMetricAtCellCenter(i,j-1,k,met_h_xi_22_lo,met_h_eta_22_lo,met_h_zeta_22_lo,
                              cellSizeInv,z_nd,TerrainMet::h_zeta);
    
    ComputeMetricAtEdgeCenterK(i+1,j,k,met_h_xi_21_hi,met_h_eta_21_hi,met_h_zeta_21_hi,
                                 cellSizeInv,z_nd,TerrainMet::all);
    ComputeMetricAtEdgeCenterK(i  ,j,k,met_h_xi_21_lo,met_h_eta_21_lo,met_h_zeta_21_lo,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);

    ComputeMetricAtEdgeCenterI(i,j,k+1,met_h_xi_23_hi,met_h_eta_23_hi,met_h_zeta_23_hi,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);
    ComputeMetricAtEdgeCenterI(i,j,k  ,met_h_xi_23_lo,met_h_eta_23_lo,met_h_zeta_23_lo,
                                 cellSizeInv,z_nd,TerrainMet::h_zeta);

    diffContrib = (tau21Next*met_h_zeta_21_hi - tau21Prev*met_h_zeta_21_lo) * dxInv  // Contribution to y-mom eqn from diffusive flux in x-dir
                + (tau22Next*met_h_zeta_22_hi - tau22Prev*met_h_zeta_22_lo) * dyInv  // Contribution to y-mom eqn from diffusive flux in y-dir
                + (tau23Next*met_h_zeta_23_hi - tau23Prev*met_h_zeta_23_lo) * dzInv; // Contribution to y-mom eqn from diffusive flux in z-dir

    diffContrib /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k)); // Terrain grid stretching

    if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
    {
      diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, v(i,j,k), Coord::y, l_spatial_order) /
        solverChoice.rho0_trans;
    }
    break;
  case MomentumEqn::z:
    Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev, normv;

    Real met_h_xi_33_hi,met_h_eta_33_hi,met_h_zeta_33_hi;
    Real met_h_xi_33_lo,met_h_eta_33_lo,met_h_zeta_33_lo;

    Real met_h_xi_31_hi,met_h_eta_31_hi,met_h_zeta_31_hi;
    Real met_h_xi_31_lo,met_h_eta_31_lo,met_h_zeta_31_lo;

    Real met_h_xi_32_hi,met_h_eta_32_hi,met_h_zeta_32_hi;
    Real met_h_xi_32_lo,met_h_eta_32_lo,met_h_zeta_32_lo;

    Real met_h_xi,met_h_eta,met_h_zeta;
    Real tau11Bar,tau22Bar,tau12Bar,tau21Bar,tau13Bar,tau23Bar;

    // Average of tau11 to tau31 location
    tau11Bar =  ComputeStressTerm(i+1,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar += ComputeStressTerm(i+2,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar += ComputeStressTerm(i+1,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar += ComputeStressTerm(i+2,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar *= 0.25;

    // Average of tau21 to tau31 location
    tau21Bar =  ComputeStressTerm(i+1,j  ,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar += ComputeStressTerm(i+1,j+1,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar += ComputeStressTerm(i+1,j  ,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar += ComputeStressTerm(i+1,j+1,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar *= 0.25;

    // Metrics at tau31Next
    ComputeMetricAtEdgeCenterJ(i+1,j,k,met_h_xi,met_h_eta,met_h_zeta,
                               cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    // Compute tau31 from jT(S-D)
    tau31Next = -met_h_xi * tau11Bar - met_h_eta * tau21Bar 
               + ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                   DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                   z_nd, domain, bc_ptr);

    
    // Average of tau11 to tau31 location
    tau11Bar =  ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar += ComputeStressTerm(i+1,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar += ComputeStressTerm(i  ,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar += ComputeStressTerm(i+1,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau11Bar *= 0.25;

    // Average of tau21 to tau31 location
    tau21Bar =  ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar += ComputeStressTerm(i  ,j+1,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar += ComputeStressTerm(i  ,j  ,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar += ComputeStressTerm(i  ,j+1,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau21Bar *= 0.25;

    // Metrics at tau31Prev
    ComputeMetricAtEdgeCenterJ(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                               cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    tau31Prev = -met_h_xi * tau11Bar - met_h_eta * tau21Bar
               + ComputeStressTerm(i  ,j  ,k  , u, v, w, momentumEqn,
                                   DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                   z_nd, domain, bc_ptr);
    //************************************************************************************
    // Average of tau12 to tau32 location
    tau12Bar =  ComputeStressTerm(i  ,j+1,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar += ComputeStressTerm(i+1,j+1,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar += ComputeStressTerm(i  ,j+1,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar += ComputeStressTerm(i+1,j+1,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar *= 0.25;

    // Average of tau22 to tau32 location
    tau22Bar =  ComputeStressTerm(i  ,j+1,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar += ComputeStressTerm(i  ,j+2,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar += ComputeStressTerm(i  ,j+1,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar += ComputeStressTerm(i  ,j+2,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar *= 0.25;

    // Metrics at tau32Next
    ComputeMetricAtEdgeCenterI(i,j+1,k,met_h_xi,met_h_eta,met_h_zeta,
                               cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    tau32Next = -met_h_xi * tau12Bar - met_h_eta * tau22Bar
               + ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                   DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                   z_nd, domain, bc_ptr);


    // Average of tau12 to tau32 location
    tau12Bar =  ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar += ComputeStressTerm(i+1,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar += ComputeStressTerm(i  ,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar += ComputeStressTerm(i+1,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau12Bar *= 0.25;

    // Average of tau22 to tau32 location
    tau22Bar =  ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar += ComputeStressTerm(i  ,j+1,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar += ComputeStressTerm(i  ,j  ,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar += ComputeStressTerm(i  ,j+1,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau22Bar *= 0.25;

    // Metrics at tau32Prev
    ComputeMetricAtEdgeCenterI(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                               cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    tau32Prev = -met_h_xi * tau12Bar - met_h_eta * tau22Bar
               + ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                   DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                   z_nd, domain, bc_ptr);
    //************************************************************************************
    // Average of tau31 to tau33 location
    tau13Bar =  ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar += ComputeStressTerm(i+1,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar += ComputeStressTerm(i  ,j  ,k+1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar += ComputeStressTerm(i+1,j  ,k+1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar *= 0.25;

    // Average of tau23 to tau33 location
    tau23Bar =  ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar += ComputeStressTerm(i  ,j+1,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar += ComputeStressTerm(i  ,j  ,k+1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar += ComputeStressTerm(i  ,j+1,k+1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar *= 0.25;

    // Metrics at tau33Next
    ComputeMetricAtCellCenter(i,j,k,met_h_xi,met_h_eta,met_h_zeta,
                              cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    tau33Next = -met_h_xi * tau13Bar - met_h_eta * tau23Bar
               + ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                   DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                   z_nd, domain, bc_ptr);

    // Average of tau31 to tau33 location
    tau13Bar =  ComputeStressTerm(i  ,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar += ComputeStressTerm(i+1,j  ,k-1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar += ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar += ComputeStressTerm(i+1,j  ,k  , u, v, w, MomentumEqn::x,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau13Bar *= 0.25;

    // Average of tau23 to tau33 location
    tau23Bar =  ComputeStressTerm(i  ,j  ,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar += ComputeStressTerm(i  ,j+1,k-1, u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar += ComputeStressTerm(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar += ComputeStressTerm(i  ,j+1,k  , u, v, w, MomentumEqn::y,
                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                  z_nd, domain, bc_ptr);
    tau23Bar *= 0.25;

    // Metrics at tau33Prev
    ComputeMetricAtCellCenter(i,j,k-1,met_h_xi,met_h_eta,met_h_zeta,
                              cellSizeInv,z_nd,TerrainMet::h_xi_eta);
    
    tau33Prev = -met_h_xi * tau13Bar - met_h_eta * tau23Bar
               + ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                   DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                   z_nd, domain, bc_ptr);

    diffContrib = (tau31Next - tau31Prev) * dxInv  // Contribution to z-mom eqn from diffusive flux in x-dir
                + (tau32Next - tau32Prev) * dyInv  // Contribution to z-mom eqn from diffusive flux in y-dir
                + (tau33Next - tau33Prev) * dzInv; // Contribution to z-mom eqn from diffusive flux in z-dir

    normv = (k == 0) ? detJ(i,j,k) : 0.5*( detJ(i,j,k) + detJ(i,j,k-1) ); // Terrain grid stretching
    diffContrib /= normv;

    if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
    {
      diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, w(i,j,k), Coord::z, l_spatial_order) /
        solverChoice.rho0_trans;
    }
    break;
  default:
    amrex::Abort("Error: Momentum equation is unrecognized");
  }

  return diffContrib;
}

#else
AMREX_GPU_DEVICE
Real
DiffusionSrcForMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                            const Array4<const Real>& cons,
                            const enum MomentumEqn &momentumEqn,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const Array4<Real>& K_turb,
                            const SolverChoice &solverChoice,
                            const Box& domain, const amrex::BCRec* bc_ptr)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real diffContrib = 0.0;

    int l_spatial_order = solverChoice.spatial_order;

    switch (momentumEqn) {
        case MomentumEqn::x:
            Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
            tau11Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau11Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau12Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau12Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau13Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau13Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);

            diffContrib = (tau11Next - tau11Prev) * dxInv  // Contribution to x-mom eqn from diffusive flux in x-dir
                        + (tau12Next - tau12Prev) * dyInv  // Contribution to x-mom eqn from diffusive flux in y-dir
                        + (tau13Next - tau13Prev) * dzInv; // Contribution to x-mom eqn from diffusive flux in z-dir
            if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
            {
                diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, u(i,j,k), Coord::x, l_spatial_order) /
                               solverChoice.rho0_trans;
            }
            break;
        case MomentumEqn::y:
            Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
            tau21Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau21Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau22Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau22Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau23Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau23Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);

            diffContrib = (tau21Next - tau21Prev) * dxInv  // Contribution to y-mom eqn from diffusive flux in x-dir
                        + (tau22Next - tau22Prev) * dyInv  // Contribution to y-mom eqn from diffusive flux in y-dir
                        + (tau23Next - tau23Prev) * dzInv; // Contribution to y-mom eqn from diffusive flux in z-dir
            if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
            {
                diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, v(i,j,k), Coord::y, l_spatial_order) /
                               solverChoice.rho0_trans;
            }
            break;
        case MomentumEqn::z:
            Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;
            tau31Next = ComputeStressTerm(i+1, j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau31Prev = ComputeStressTerm(i  , j, k, u, v, w, momentumEqn,
                                          DiffusionDir::x, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau32Next = ComputeStressTerm(i, j+1, k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau32Prev = ComputeStressTerm(i, j  , k, u, v, w, momentumEqn,
                                          DiffusionDir::y, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau33Next = ComputeStressTerm(i, j, k+1, u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);
            tau33Prev = ComputeStressTerm(i, j, k  , u, v, w, momentumEqn,
                                          DiffusionDir::z, cellSizeInv, K_turb, solverChoice, domain, bc_ptr);

            diffContrib = (tau31Next - tau31Prev) * dxInv  // Contribution to z-mom eqn from diffusive flux in x-dir
                        + (tau32Next - tau32Prev) * dyInv  // Contribution to z-mom eqn from diffusive flux in y-dir
                        + (tau33Next - tau33Prev) * dzInv; // Contribution to z-mom eqn from diffusive flux in z-dir
            if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
            {
                diffContrib *= InterpolateFromCellOrFace(i, j, k, cons, Rho_comp, w(i,j,k), Coord::z, l_spatial_order) /
                               solverChoice.rho0_trans;
            }
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return diffContrib;
}
#endif

AMREX_GPU_DEVICE
amrex::Real ComputeDiffusionFluxForState(const int &i, const int &j, const int &k,
                     const Array4<const Real>& cell_data,
                     const Array4<const Real>& cell_prim, const int & prim_index,
                     const amrex::Real invCellWidth,
                     const Array4<Real>& K_turb,
                     const SolverChoice &solverChoice,
                     const enum Coord& coordDir)
{
  // Get indices of states to left and right of the face on which we want the flux
  const int il = i - (coordDir == Coord::x);
  const int ir = i;
  const int jl = j - (coordDir == Coord::y);
  const int jr = j;
  const int kl = k - (coordDir == Coord::z);
  const int kr = k;

  // Get diffusion coefficients
  amrex::Real rhoAlpha_molec;
  int eddy_diff_idx;

  amrex::Real rhoFace;
  if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
    rhoFace = (cell_data(il, jl, kl, Rho_comp) + cell_data(ir, jr, kr, Rho_comp)) * 0.5;
  }

  switch(prim_index) {
      case PrimTheta_comp: // Potential Temperature
          if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
              rhoAlpha_molec = rhoFace * solverChoice.alpha_T;
          } else {
              rhoAlpha_molec = solverChoice.rhoAlpha_T;
          }
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::Theta_v;
          } else {
              eddy_diff_idx = EddyDiff::Theta_h;
          }
          break;

      case PrimKE_comp: // Turbulent KE
          rhoAlpha_molec = 0.;
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::KE_v;
          } else {
              eddy_diff_idx = EddyDiff::KE_h;
          }
          break;

      case PrimQKE_comp: // Turbulent QKE
          rhoAlpha_molec = 0.;
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::QKE_v;
          } else {
              eddy_diff_idx = EddyDiff::QKE_h;
          }
          break;

      case PrimScalar_comp: // Scalar
          if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) {
              rhoAlpha_molec = rhoFace * solverChoice.alpha_C;
          } else {
              rhoAlpha_molec = solverChoice.rhoAlpha_C;
          }
          if (coordDir == Coord::z) {
              eddy_diff_idx = EddyDiff::Scalar_v;
          } else {
              eddy_diff_idx = EddyDiff::Scalar_h;
          }
          break;
      default:
          amrex::Abort("Error: Diffusion term for the data index isn't implemented");
  }

  amrex::Real rhoAlpha = 0.0;

  if ( (solverChoice.molec_diff_type == MolecDiffType::Constant) ||
       (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha) ) {
    rhoAlpha += rhoAlpha_molec;
  }

  if ( (solverChoice.les_type == LESType::Smagorinsky) ||
       (solverChoice.les_type == LESType::Deardorff  ) ||
       (solverChoice.pbl_type == PBLType::MYNN25     ) ) {
    rhoAlpha += 0.5*(K_turb(ir,jr,kr,eddy_diff_idx) + K_turb(il,jl,kl,eddy_diff_idx));
  }

  // Compute the flux
  amrex::Real diffusionFlux = rhoAlpha * invCellWidth *
      (cell_prim(ir, jr, kr, prim_index) - cell_prim(il, jl, kl, prim_index));

  return diffusionFlux;
}

AMREX_GPU_DEVICE
Real
DiffusionSrcForState(const int &i, const int &j, const int &k,
                              const Array4<const Real>& cell_data,
                              const Array4<const Real>& cell_prim, const int & qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const Array4<Real>& K_turb,
                              const SolverChoice &solverChoice)
{
  const amrex::Real dx_inv = cellSizeInv[0];
  const amrex::Real dy_inv = cellSizeInv[1];
  const amrex::Real dz_inv = cellSizeInv[2];

  const int prim_index = qty_index - RhoTheta_comp;

  // TODO : could be more efficient to compute and save all fluxes before taking divergence (now all fluxes are computed 2x);
  xflux(i+1,j,k,qty_index) = ComputeDiffusionFluxForState(i+1, j, k, cell_data, cell_prim, prim_index, dx_inv, K_turb, solverChoice, Coord::x);
  xflux(i  ,j,k,qty_index) = ComputeDiffusionFluxForState(i  , j, k, cell_data, cell_prim, prim_index, dx_inv, K_turb, solverChoice, Coord::x);

  yflux(i,j+1,k,qty_index) = ComputeDiffusionFluxForState(i, j+1, k, cell_data, cell_prim, prim_index, dy_inv, K_turb, solverChoice, Coord::y);
  yflux(i,j  ,k,qty_index) = ComputeDiffusionFluxForState(i, j  , k, cell_data, cell_prim, prim_index, dy_inv, K_turb, solverChoice, Coord::y);

  zflux(i,j,k+1,qty_index) = ComputeDiffusionFluxForState(i, j, k+1, cell_data, cell_prim, prim_index, dz_inv, K_turb, solverChoice, Coord::z);
  zflux(i,j,k  ,qty_index) = ComputeDiffusionFluxForState(i, j, k  , cell_data, cell_prim, prim_index, dz_inv, K_turb, solverChoice, Coord::z);

  Real diffusionSrc =
      (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dx_inv   // Diffusive flux in x-dir
     +(yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dy_inv   // Diffusive flux in y-dir
     +(zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dz_inv;  // Diffusive flux in z-dir

  return diffusionSrc;
}
